#include <cmath>
#include "../../../src/Accelerator/hybrid.h"
#include "../../../src/Constants/behavior.h"
#include "../../../src/Constants/scaling.h"
#include "../../../src/Constants/symbol_values.h"
#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/DataTypes/stormm_vector_types.h"
#include "../../../src/FileManagement/file_listing.h"
#include "../../../src/Math/rounding.h"
#include "../../../src/Math/series_ops.h"
#include "../../../src/Math/statistical_enumerators.h"
#include "../../../src/Math/summation.h"
#include "../../../src/MoleculeFormat/molecule_parsing.h"
#include "../../../src/Parsing/parse.h"
#include "../../../src/Reporting/error_format.h"
#include "../../../src/Structure/clash_detection.h"
#include "../../../src/Structure/isomerization.h"
#include "../../../src/Topology/atomgraph_abstracts.h"
#include "../../../src/Topology/topology_util.h"
#include "../../../src/Trajectory/coordinateframe.h"
#include "setup.h"

namespace conf_app {
namespace setup {

using stormm::card::HybridFormat;
using stormm::chemistry::ChiralOrientation;
using stormm::chemistry::ConformationEdit;
using stormm::chemistry::MapRotatableGroups;
using stormm::chemistry::permutationsAreLinked;
using stormm::constants::PrecisionModel;
using stormm::constants::warp_size_int;
#ifndef STORMM_USE_HPC
using stormm::data_types::double2;
#endif
using stormm::diskutil::getBaseName;
using stormm::errors::rtErr;
using stormm::errors::rtWarn;
using stormm::namelist::ConformerControls;
using stormm::namelist::SamplingIntensity;
using stormm::parse::char4ToString;
using stormm::stmath::DataOrder;
using stormm::stmath::incrementingSeries;
using stormm::stmath::prefixSumInPlace;
using stormm::stmath::PrefixSumType;
using stormm::stmath::roundUp;
using stormm::stmath::sum;
using stormm::stmath::loadScalarStateValues;
using stormm::structure::ClashReport;
using stormm::structure::detectClash;
using stormm::structure::flipChiralCenter;
using stormm::structure::maskFromSdfDataItem;
using stormm::structure::rotateAboutBond;
using stormm::synthesis::SynthesisMapReader;
using stormm::topology::ChemicalDetailsKit;
using stormm::topology::ImplicitSolventModel;
using stormm::topology::MobilitySetting;
using stormm::topology::isBonded;
using stormm::trajectory::CoordinateFrame;

//-------------------------------------------------------------------------------------------------
AtomMask getCoreMask(const ConformerControls &conf_input, const MdlMol &sdf_example,
                     const PhaseSpace &ps, const AtomGraph *ag, const ChemicalFeatures &chemfe,
                     const ExceptionResponse policy) {
  AtomMask result(ag);
  bool core_found = false;
  if (conf_input.getCoreDataItemName().size() > 0LLU && sdf_example.getAtomCount() > 0) {
    result.addAtoms(maskFromSdfDataItem(conf_input.getCoreDataItemName(), sdf_example, ag, chemfe,
                                        sdf_example.exportCoordinateFrame(), policy),
                    sdf_example.exportCoordinateFrame(), chemfe);
    core_found = (result.getMaskedAtomCount() > 0);
  }
  if (core_found == false && conf_input.getCoreAtomMask().size() > 0LLU) {

    // If present, a core atom mask common to all systems will fill in for situations in which a
    // list of core atoms is not defined in an SD file entry.
    result.addAtoms(conf_input.getCoreAtomMask(), CoordinateFrame(ps), chemfe);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void setGenerativeConditions(const UserSettings &ui, SystemCache *sc, StopWatch *tm) {

  // Establish a category for this function and stash time up to this point in the "miscellaneous"
  // category.
  const int cat_no = tm->addCategory("Set generative conditions");
  tm->assignTime(0);
  
  // Add implicit solvent conditions.
  std::vector<AtomGraph*> all_topologies = sc->getTopologyPointer();
  const size_t ntop = all_topologies.size();
  const ImplicitSolventModel igb = ui.getSolventNamelistInfo().getImplicitSolventModel();
  if (igb != ImplicitSolventModel::NONE) {
    for (size_t i = 0; i < ntop; i++) {
      all_topologies[i]->setImplicitSolventModel(igb);
    }
  }

  // Create the core masks for all topologies.  Build a list of chemical features for each unique
  // topology in the systems cache.
  const ConformerControls& conf_input = ui.getConformerNamelistInfo();  
  const int nsys = sc->getSystemCount();
  for (int i = 0; i < nsys; i++) {
    AtomGraph* iag_ptr = sc->getSystemTopologyPointer(i);
    const AtomMask imask = getCoreMask(conf_input, sc->getStructureDataEntry(i),
                                       sc->getCoordinates(i), iag_ptr, sc->getFeatures(i),
                                       ui.getExceptionBehavior());
    const std::vector<int> icore_atoms = imask.getMaskedAtomList();
    const size_t n_core_atoms = icore_atoms.size();
    std::vector<bool> icore_mobility(icore_atoms.size());
    for (size_t i = 0LLU; i < n_core_atoms; i++) {
      icore_mobility[i] = iag_ptr->getAtomMobility(icore_atoms[i]);
    }
    iag_ptr->modifyAtomMobility(icore_atoms, MobilitySetting::OFF);

    // With the topology's core atoms now frozen, compute the rotatable bond groups.
    sc->getFeaturesPointer(i)->findRotatableBondGroups();

    // Restore the core atoms' original mobility settings (henceforth, they will be restrained to
    // their original positions using harmonic restraints).
    for (size_t i = 0; i < n_core_atoms; i++) {
      const MobilitySetting imb = (icore_mobility[i]) ? MobilitySetting::ON : MobilitySetting::OFF;
      iag_ptr->modifyAtomMobility(icore_atoms[i], imb);
    }
  }

  // Record the time needed for this procedure.
  tm->assignTime(cat_no);
}
  
//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis buildSamplingWorkspace(const SystemCache &sc, const ConformerControls &confcon,
                                           const MinimizeControls &mincon,
                                           Xoshiro256ppGenerator *xrs, SynthesisCacheMap *scmap,
                                           SynthesisPermutor *syper, StopWatch *tm) {
  tm->assignTime(0);
  PhaseSpaceSynthesis seed_structures(sc.getCoordinates(), sc.getSystemTopologyPointer());
  std::vector<ChemicalFeatures*> sc_unique_features(sc.getTopologyCount());
  for (int i = 0; i < sc.getTopologyCount(); i++) {
    const int cache_example = sc.getSystemExampleIndex(i);
    sc_unique_features[i] = const_cast<ChemicalFeatures*>(sc.getFeaturesPointer(cache_example));
  }

  // The seed structures may need to have chirality set to particular values.  This will only
  // occur if chiral sampling is turned off.  The initial settings will be retained by any
  // downstream conformers, barring extreme clashes or strain that further perturbs the
  // orientations.
  const int n_chiral_settings = confcon.getChiralSettingCount();
  for (int i = 0; i < n_chiral_settings; i++) {
    const std::string& i_label = confcon.getChiralSettingLabel(i);
    const std::vector<int> i_systems = sc.getMatchingSystemIndices(i_label);
    const ChiralOrientation i_setting = confcon.getChiralSetting(i);
    const size_t i_nsys = i_systems.size();
    for (size_t j = 0; j < i_nsys; j++) {

      // The preliminary coordinate synthesis, seed_structures, aligns with the systems cache sc.
      // The chemical features of the jth system in the cache are valid for the jth structure in
      // the synthesis.
      const ChemicalFeatures& ijfeat = sc.getFeatures(i_systems[j]);
      const std::vector<ChiralInversionProtocol> ij_protocols = ijfeat.getChiralInversionMethods();
      const std::vector<IsomerPlan> ij_plan = ijfeat.getChiralInversionGroups();
      const std::vector<int> ij_centers = ijfeat.getChiralCenters();
      CoordinateFrame ij_cf = seed_structures.exportCoordinates(i_systems[j],
                                                                HybridFormat::HOST_ONLY);
      const AtomMask ij_mask(confcon.getChiralSettingMask(i),
                             seed_structures.getSystemTopologyPointer(i_systems[j]), ijfeat,
                             seed_structures.exportCoordinates(i_systems[j]));
      const int nij_atoms = ij_mask.getMaskedAtomCount();
      const std::vector<int> ij_atoms = ij_mask.getMaskedAtomList();
      for (int k = 0; k < nij_atoms; k++) {
        
        // The atom in question is known by its topological index, but its position in the list of
        // chiral centers must be known.
        const int ij_ccen_idx = locateValue(ij_centers, ij_atoms[k], DataOrder::ASCENDING);
        if (ij_ccen_idx >= ij_centers.size()) {
          const AtomGraph* ij_ag = ijfeat.getTopologyPointer();
          const int kres_idx = ij_ag->getResidueIndex(ij_atoms[k]);
          rtWarn("Chiral center " + char4ToString(ij_ag->getAtomName(ij_atoms[k])) +
                 " in residue " + char4ToString(ij_ag->getResidueName(kres_idx)) + "(" +
                 std::to_string(kres_idx) + ") was not located in the list of " +
                 std::to_string(ij_centers.size()) + " identified chiral centers.",
                 "buildSamplingWorkspace");
          continue;
        }
        const ChiralOrientation k_orient = ijfeat.getAtomChirality(ij_atoms[k]);
        if (k_orient != ChiralOrientation::NONE && k_orient != i_setting) {
          flipChiralCenter(&ij_cf, ij_ccen_idx, ij_centers, ij_protocols, ij_plan);
        }
      }
    }
  }
  
  SynthesisPermutor tmp_syper(sc_unique_features, seed_structures, confcon);
  ClashReport clrep(mincon.getAbsoluteClashDistance(), mincon.getVdwClashRatio());
  std::vector<int> correspondence;
  PhaseSpaceSynthesis result = tmp_syper.buildSynthesis(confcon.getSamplingIntensity(), xrs,
                                                        &correspondence, PrecisionModel::SINGLE,
                                                        confcon.getSamplingTrialLimit() *
                                                        tmp_syper.getSystemCount(), &clrep,
                                                        confcon.getMaxSeedingAttempts(),
                                                        confcon.getClashPairTolerance());
  
  // The correspondence traces structures in the result to structures in the synthesis linked to
  // tmp_syper, and the structures in that have a 1:1 correspondence with systems in the original
  // cache created from user input.
  scmap->setCache(correspondence, sc.getSelfPointer());
  scmap->setSynthesis(result);
  
  // The new synthesis can now be applied to the permutor object, and by an assignment operation
  // the permutor object developed herein can be passed back up to the main program.
  tmp_syper.applySynthesis(result);

  *syper = tmp_syper;
  tm->assignTime(tm_coordinate_expansion);
  return result;
}

//-------------------------------------------------------------------------------------------------
const std::vector<RestraintApparatus*>
buildReplicaRestraints(const SynthesisCacheMap &scmap) {
  const SynthesisMapReader scmapr = scmap.data();
  const SystemCache *sc = scmap.getCachePointer();
  std::vector<RestraintApparatus*> result(scmapr.nsynth);
  for (int i = 0; i < scmapr.nsynth; i++) {
    result[i] = const_cast<RestraintApparatus*>(sc->getRestraintPointer(scmapr.cache_origins[i]));
  }
  return result;  
}
  
} // namespace setup
} // namespace conf_app
