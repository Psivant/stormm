#include <map>
#include <string>
#include <vector>
#include "../../../src/copyright.h"
#include "../../../src/Accelerator/hpc_config.h"
#include "../../../src/Accelerator/hybrid.h"
#include "../../../src/Chemistry/atommask.h"
#include "../../../src/Chemistry/chemical_features.h"
#include "../../../src/Constants/behavior.h"
#include "../../../src/Constants/generalized_born.h"
#include "../../../src/DataTypes/stormm_vector_types.h"
#include "../../../src/FileManagement/file_listing.h"
#include "../../../src/Math/series_ops.h"
#include "../../../src/MoleculeFormat/pdb.h"
#include "../../../src/MolecularMechanics/mm_evaluation.h"
#include "../../../src/Namelists/command_line_parser.h"
#include "../../../src/Namelists/input_transcript.h"
#include "../../../src/Namelists/namelist_inventory.h"
#include "../../../src/Namelists/nml_emulate.h"
#include "../../../src/Namelists/nml_files.h"
#include "../../../src/Namelists/nml_minimize.h"
#include "../../../src/Namelists/nml_precision.h"
#include "../../../src/Namelists/nml_random.h"
#include "../../../src/Namelists/nml_solvent.h"
#include "../../../src/Namelists/user_settings.h"
#include "../../../src/Numerics/split_fixed_precision.h"
#include "../../../src/Parsing/parse.h"
#include "../../../src/Parsing/parsing_enumerators.h"
#include "../../../src/Potential/scorecard.h"
#include "../../../src/Potential/surface_area.h"
#include "../../../src/Reporting/error_format.h"
#include "../../../src/Reporting/help_messages.h"
#include "../../../src/Reporting/report_table.h"
#include "../../../src/Reporting/section_contents.h"
#include "../../../src/Synthesis/atomgraph_synthesis.h"
#include "../../../src/Synthesis/phasespace_synthesis.h"
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/Topology/atomgraph_enumerators.h"
#include "../../../src/Trajectory/coordinateframe.h"
#include "../../../src/Trajectory/coordinate_graft.h"
#include "../../../src/UnitTesting/stopwatch.h"
#include "emulate_analysis.h"
#include "emulate_enumerators.h"

using stormm::testing::StopWatch;
using namespace stormm::card;
using namespace stormm::chemistry;
using namespace stormm::constants;
using namespace stormm::data_types;
using namespace stormm::diskutil;
using namespace stormm::display;
using namespace stormm::energy;
using namespace stormm::errors;
using namespace stormm::generalized_born_defaults;
using namespace stormm::mm;
using namespace stormm::namelist;
using namespace stormm::numerics;
using namespace stormm::parse;
using namespace stormm::review;
using namespace stormm::stmath;
using namespace stormm::structure;
using namespace stormm::synthesis;
using namespace stormm::topology;
using namespace stormm::trajectory;
using namespace emulation;

//-------------------------------------------------------------------------------------------------
// Loop over all (masked) atoms in a structure and contribute the values of basis functions for
// interacting atoms to the appropriate columns of a matrix.
//
// Arguments:
//   
//-------------------------------------------------------------------------------------------------
void contributeToFit(size_t row_index, const size_t total_rows, PhaseSpaceWriter *psw,
                     const std::vector<int> &source_id, const std::vector<int> &parameter_map,
                     const EmulatorControls &emulcon, const std::vector<int> &masked_atoms,
                     const StaticExclusionMask &sem, Hybrid<double> *amat, const double mult,
                     const FittingContribution mode) {

  // Unpack the basis functions.
  const int nbss = emulcon.getBasisFunctionCount();
  std::vector<double> bss_fnc_width(nbss);
  std::vector<double> bss_fnc_start(nbss);
  double max_dist = 0.0;
  for (int i = 0; i < nbss; i++) {
    bss_fnc_width[i] = emulcon.getSupportWidth(i);
    bss_fnc_start[i] = emulcon.getSupportStart(i);
    max_dist = std::max(max_dist, bss_fnc_width[i] + bss_fnc_start[i]);
  }
  const double max_sq_dist = max_dist * max_dist;
  
  // Get the number of atom sources.
  const int nsrc = emulcon.getSourceCount();
  const int nsrc_sq = nsrc * nsrc;
  
  // Detect a complete mask and include all atoms.
  const bool mask_all = (masked_atoms[0] < 0); 
  const int ilim = (mask_all) ? psw->natom : masked_atoms.size();
  double* amat_ptr = amat->data();
  for (int i = 0; i < ilim; i++) {
    const size_t atomi = (mask_all) ? i : masked_atoms[i];
    if (atomi < 0 || source_id[atomi] < 0) {
      continue;
    }
    const int isrc_ofs = nsrc * source_id[atomi];
    for (int j = 0; j < i; j++) {
      const size_t atomj = (mask_all) ? j : masked_atoms[j];
      if (atomj < 0 || source_id[atomj] < 0) {
        continue;
      }
      if (! sem.testExclusion(atomi, atomj)) {

        // This non-excluded interaction may be relevant.  First, test the atom types to ensure
        // that parameters are available in the fit.
        if (parameter_map[isrc_ofs + source_id[atomj]] >= 0) {

          // Check whether the two atoms are within range of each other.
          const double dx   = psw->xcrd[atomj] - psw->xcrd[atomi];
          const double dy   = psw->ycrd[atomj] - psw->ycrd[atomi];
          const double dz   = psw->zcrd[atomj] - psw->zcrd[atomi];
          const double r_sq = (dx * dx) + (dy * dy) + (dz * dz);
          if (r_sq < max_sq_dist) {
            const double r = sqrt(r_sq);
            switch (mode) {
            case FittingContribution::ENERGY:
              for (int k = 0; k < nbss; k++) {
                if (r < bss_fnc_start[k]) {
                  const size_t col_idx = parameter_map[(k * nsrc_sq) + isrc_ofs +
                                                       source_id[atomj]];
                  amat_ptr[(col_idx * total_rows) + row_index] += 1.0;
                }
                else if (r < bss_fnc_start[k] + bss_fnc_width[k]) {
                  const size_t col_idx = parameter_map[(k * nsrc_sq) + isrc_ofs +
                                                       source_id[atomj]];
                  const double rp = (r - bss_fnc_start[k]) / bss_fnc_width[k];
                  const double bss_fnc_val = (((2.0 * rp) - 3.0) * rp * rp) + 1.0;
                  amat_ptr[(col_idx * total_rows) + row_index] += bss_fnc_val;
                }
              }
              break;
            case FittingContribution::FORCE:
              for (int k = 0; k < nbss; k++) {
                if (r > bss_fnc_start[k] && r < bss_fnc_start[k] + bss_fnc_width[k]) {
                  const size_t col_idx = parameter_map[(k * nsrc_sq) + isrc_ofs +
                                                       source_id[atomj]];
                  const double rp = (r - bss_fnc_start[k]) / bss_fnc_width[k];
                  amat_ptr[(col_idx * total_rows) + row_index] += ((6.0 * rp) - 6.0) * r;
                }
              }
              break;
            }
          }
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
// 
//
// Arguments:
//   ps:
//-------------------------------------------------------------------------------------------------
double energyContextForFit(PhaseSpace *ps, const AtomGraph *ag, const StaticExclusionMask &sem,
                           const std::vector<StateVariable> &mm_parts,
                           const NeckGeneralizedBornKit<double> &ngb_kit, ScoreCard *sc) {
  PhaseSpaceWriter psw = ps->data();
  const ValenceKit<double> vk = ag->getDoublePrecisionValenceKit();
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  const StaticExclusionMaskReader semr = sem.data();
  const size_t n_mm = mm_parts.size();
  double result = 0.0;
  for (size_t i = 0; i < n_mm; i++) {
    switch (mm_parts[i]) {
    case StateVariable::BOND:
      result += evaluateBondTerms(vk, psw, sc, EvaluateForce::NO);
      break;
    case StateVariable::ANGLE:
      result += evaluateAngleTerms(vk, psw, sc, EvaluateForce::NO);
      break;
    case StateVariable::PROPER_DIHEDRAL:
      {
        const double2 dihe_e = evaluateDihedralTerms(vk, psw, sc, EvaluateForce::NO);
        result += dihe_e.x;
      }
      break;
    case StateVariable::IMPROPER_DIHEDRAL:
      {
        const double2 impr_e = evaluateDihedralTerms(vk, psw, sc, EvaluateForce::NO);
        result += impr_e.y;
      }
      break;
    case StateVariable::UREY_BRADLEY:
      result += evaluateUreyBradleyTerms(vk, psw, sc, EvaluateForce::NO);
      break;
    case StateVariable::CHARMM_IMPROPER:
      result += evaluateCharmmImproperTerms(vk, psw, sc, EvaluateForce::NO);
      break;
    case StateVariable::CMAP:
      result += evaluateCmapTerms(vk, psw, sc, EvaluateForce::NO);
      break;
    case StateVariable::ELECTROSTATIC:
      {
        const double2 nb_e = evaluateNonbondedEnergy(nbk, semr, psw, sc, EvaluateForce::NO,
                                                     EvaluateForce::NO);
        const double2 nb_14_e = evaluateAttenuated14Terms(vk, nbk, psw, sc, EvaluateForce::NO,
                                                          EvaluateForce::NO);
        result += nb_e.x + nb_14_e.x;
      }
      break;
    case StateVariable::VDW:
      {
        const double2 nb_e = evaluateNonbondedEnergy(nbk, semr, psw, sc, EvaluateForce::NO,
                                                     EvaluateForce::NO);
        const double2 nb_14_e = evaluateAttenuated14Terms(vk, nbk, psw, sc, EvaluateForce::NO,
                                                          EvaluateForce::NO);
        result += nb_e.y + nb_14_e.y;
      }
      break;
    case StateVariable::GENERALIZED_BORN:
      {
        const ImplicitSolventKit<double> isk = ag->getDoublePrecisionImplicitSolventKit();
        result += evaluateGeneralizedBornEnergy(nbk, semr, isk, ngb_kit, psw, sc,
                                                EvaluateForce::NO);
      }
      break;
    case StateVariable::ELEC_ONE_FOUR:
    case StateVariable::VDW_ONE_FOUR:
    case StateVariable::RESTRAINT:
    case StateVariable::SURFACE_AREA:
    case StateVariable::KINETIC:
    case StateVariable::PRESSURE:
    case StateVariable::VIRIAL_11:
    case StateVariable::VIRIAL_12:
    case StateVariable::VIRIAL_22:
    case StateVariable::VIRIAL_13:
    case StateVariable::VIRIAL_23:
    case StateVariable::VIRIAL_33:
    case StateVariable::VOLUME:
    case StateVariable::TEMPERATURE_ALL:
    case StateVariable::TEMPERATURE_PROTEIN:
    case StateVariable::TEMPERATURE_LIGAND:
    case StateVariable::TEMPERATURE_SOLVENT:
    case StateVariable::DU_DLAMBDA:
    case StateVariable::POTENTIAL_ENERGY:
    case StateVariable::TOTAL_ENERGY:
    case StateVariable::ALL_STATES:
      break;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(int argc, const char* argv[]) {

  // Wall time tracking
  StopWatch master_timer("Master timings for mmgbsa.stormm");
  master_timer.addCategory("Input parsing");
  master_timer.addCategory("Work unit building");
#ifdef STORMM_USE_HPC
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  Hybrid<int> force_gpu_to_engage(1);
#else
  const GpuDetails gpu = null_gpu;
#endif
  // Parse the command line
  CommandLineParser clip("emulate.stormm", "A program for creating custom pairwise interaction "
                         "potentials for reproducing target energies of known structures.");
  clip.addStandardApplicationInputs({ "-i", "-O", "-o", "-except" });
  NamelistEmulator *t_nml = clip.getNamelistPointer();
  const std::vector<std::string> my_namelist_names = { "&files", "&emulator", "&minimize",
                                                       "&restraint", "&solvent", "&report" };
  clip.addControlBlocks(my_namelist_names);
  if (displayNamelistHelp(argc, argv, my_namelist_names) && clip.doesProgramExitOnHelp()) {
    return 0;
  }
  clip.parseUserInput(argc, argv); 
    
  // If testing is requested, perform it and exit.
#if 0
  if (t_nml->getBoolValue("-unittest")) {
    return runUnitTests();
  }
  if (t_nml->getBoolValue("-regtest")) {
    return runRegressionTests();
  }
#endif
  // Take in the user input from the input file.  Immediately test for valid inputs.
  const UserSettings ui(clip, { "-pe", "-ce", "-rg" });
  FilesControls ficon = ui.getFilesNamelistInfo();
  const ExceptionResponse policy = translateExceptionResponse(t_nml->getStringValue("-except"));
  const PrintSituation prnt_protocol = t_nml->getBoolValue("-O") ? PrintSituation::OVERWRITE :
                                                                   PrintSituation::OPEN_NEW;
  EmulatorControls emulcon = ui.getEmulatorNamelistInfo();
  const NeckGeneralizedBornTable ngb_tab;
  const NeckGeneralizedBornKit<double> ngb_kit = ngb_tab.dpData();
  const int nsrc_def  = emulcon.getSourceCount();
  const int n_targets = emulcon.getTargetCount();
  const int nbss_fnc  = emulcon.getBasisFunctionCount();
  
  // Load all of the systems
  SystemCache sysche(ficon, policy, MapRotatableGroups::YES, prnt_protocol, &master_timer);
  std::vector<bool> sources_present(nsrc_def, false);
  std::vector<std::vector<int>> source_id;
  source_id.reserve(sysche.getSystemCount());
  for (int i = 0; i < sysche.getSystemCount(); i++) {
    const AtomGraph *ag_i = sysche.getSystemTopologyPointer(i);
    const ChemicalFeatures& chemfe_i = sysche.getFeatures(i);
    const ChemicalFeaturesReader chemfe_ir = chemfe_i.data();
    const int natom = ag_i->getAtomCount();
    source_id.emplace_back(natom, -1);
    std::vector<bool> iarom_mask(natom, false);
    for (int j = 0; j < chemfe_ir.aromatic_group_count; j++) {
      for (int k = chemfe_ir.aromatic_group_bounds[j]; k < chemfe_ir.aromatic_group_bounds[j + 1];
           k++) {
        iarom_mask[k] = true;
      }
    }
    for (int j = 0; j < natom; j++) {
      bool unmatched = true;
      int k = 0;
      while (unmatched && k < nsrc_def) {
        const NBAtomSource& ksrc = emulcon.getSource(k);
        if (ksrc.atomMatches(j, chemfe_i, iarom_mask)) {
          source_id[i][j] = k;
          sources_present[k] = true;
          unmatched = false;
        }
        k++;
      }
    }
  }
  
  // Check the cache of systems against the various targets--are all targets valid?
  const int ncache_systems = sysche.getSystemCount();
  for (int i = 0; i < n_targets; i++) {
    const EmulationTarget& itrg = emulcon.getTarget(i);
    const std::string& fstr_label = itrg.getFirstStructureLabel();
    if (sysche.getLabelCacheIndex(fstr_label) == ncache_systems) {
      rtErr("Structure label " + fstr_label + " for target " + std::to_string(i) + " did not "
            "match any label in the cache enumerated in the &files control block.", "main");
    }
    if (itrg.getFirstStructureFrame() >= sysche.getSystemCountWithinLabel(fstr_label)) {
      rtErr("Structure label " + fstr_label + " for target " + std::to_string(i) + " did not "
            "match any label in the cache enumerated in the &files control block.", "main");
    }
  }
  
  // Count the columns of the matrix, beginning with the number of unique parameters and moving on
  // to the number of floating baselines for different systems.  The parameter_map array tracks the
  // column to which each parameter pair is assigned.
  std::vector<bool> sources_in_range(nbss_fnc * nsrc_def * nsrc_def, false);
  int n_parameters = 0;
  std::vector<int> parameter_map(nbss_fnc * nsrc_def * nsrc_def, -1);
  for (int bcon = 0; bcon < nbss_fnc; bcon++) {
    const double max_dist = emulcon.getSupportWidth(bcon) + emulcon.getSupportStart(bcon);
    const double max_sq_dist = max_dist * max_dist;
    for (int i = 0; i < sysche.getSystemCount(); i++) {
      const PhaseSpaceReader i_psr = sysche.getCoordinatePointer(i)->data();
      for (int j = 0; j < i_psr.natom; j++) {
        const int src_ij = source_id[i][j];
        if (src_ij >= 0) {
          for (int k = 0; k < j; k++) {
            const int src_ik = source_id[i][k];
            if (src_ik >= 0 &&
                parameter_map[(((bcon * nsrc_def) + src_ij) * nsrc_def) + src_ik] == -1) {
              const double dx = i_psr.xcrd[k] - i_psr.xcrd[j];
              const double dy = i_psr.ycrd[k] - i_psr.ycrd[j];
              const double dz = i_psr.zcrd[k] - i_psr.zcrd[j];
              const double r2 = (dx * dx) + (dy * dy) + (dz * dz);
              if (r2 < max_sq_dist) {
                parameter_map[(((bcon * nsrc_def) + src_ij) * nsrc_def) + src_ik] = n_parameters;
                parameter_map[(((bcon * nsrc_def) + src_ik) * nsrc_def) + src_ij] = n_parameters;
                n_parameters++;
              }
            }
          }
        }
      }
    }
  }
  
  // Count the number of floating baseline energies, extra parameters that augment the fitting
  // matrix but do not carry over to the resulting force field.  
  const std::vector<AtomGraph*> topology_image = sysche.getSystemTopologyPointer();
  const int unique_topl_count = topology_image.size();
  std::vector<std::vector<std::string>> mask_image(topology_image.size());
  std::vector<std::vector<int>> floating_baseline_counts(topology_image.size());
  std::vector<int> baseline_needed(sysche.getTopologyCount());
  std::vector<AtomGraph*> target_top1(n_targets, nullptr);
  std::vector<AtomGraph*> target_top2(n_targets, nullptr);
  std::vector<AtomGraph*> target_top3(n_targets, nullptr);
  std::vector<std::vector<std::vector<int>>> baseline_usage(topology_image.size(),
                                                            std::vector<std::vector<int>>());
  for (int i = 0; i < n_targets; i++) {
    const EmulationTarget& itrg = emulcon.getTarget(i);
    if (itrg.useFloatingBaseline()) {

      // Identify the system and subset of atoms that has a floating baseline
      const int chc_idx = sysche.getCacheIndex(itrg.getFirstStructureLabel(),
                                               itrg.getFirstStructureFrame());
      const AtomGraph* ag_ptr = sysche.getSystemTopologyPointer(chc_idx);
      const std::string &mask_ref = itrg.getFirstStructureMask();
      const int topl_idx = sysche.getSystemTopologyIndex(chc_idx);
      const int mask_pos = findStringInVector(mask_image[topl_idx], mask_ref);
      if (mask_pos == mask_image[topl_idx].size()) {
        mask_image[topl_idx].push_back(mask_ref);
        floating_baseline_counts[topl_idx].push_back(1);
        baseline_usage[topl_idx].push_back(std::vector<int>(1, i));
      }
      else {
        floating_baseline_counts[topl_idx][mask_pos] += 1;
        baseline_usage[topl_idx][mask_pos].push_back(i);
      }
    }
  }
  std::vector<int> baseline_map(n_targets, -1);
  std::vector<std::string> baseline_desc;
  int baseline_counter = 0;
  for (int i = 0; i < unique_topl_count; i++) {
    const int ni_unique_mask = mask_image[i].size();
    for (int j = 0; j < ni_unique_mask; j++) {
      const int nj_systems_sharing = baseline_usage[i][j].size();
      for (int k = 0; k < nj_systems_sharing; k++) {
        baseline_map[baseline_usage[i][j][k]] = baseline_counter;
      }
      baseline_desc.push_back(topology_image[i]->getFileName() + " " + mask_image[i][j]);
      baseline_counter++;
    }
  }
  const int n_baseline = baseline_counter;

  // Lay out the column descriptions
  std::vector<std::string> column_descriptors(n_parameters + n_baseline);
  for (int i = 0; i < nsrc_def; i++) {
    const NBAtomSource& i_src = emulcon.getSource(i);
    for (int j = 0; j < nsrc_def; j++) {
      const NBAtomSource& j_src = emulcon.getSource(j);
      for (int k = 0; k < nbss_fnc; k++) {
        const int parm_idx = (((k * nsrc_def) + i) * nsrc_def) + j;
        if (parameter_map[parm_idx] >= 0) {
          column_descriptors[parameter_map[parm_idx]] = "Scaling factor, " + i_src.getLabel() +
                                                        " :: " + j_src.getLabel() + ", basis " +
                                                        std::to_string(k + 1);
        }
      }
    }
  }
  for (int i = 0; i < n_baseline; i++) {
    column_descriptors[i + n_parameters] = "Adjustment, " + baseline_desc[i];
  }
  const size_t n_columns = n_parameters + n_baseline;

  // Map the rows of the matrix.
  const size_t n_rows = n_targets;

  // CHECK
  printf("Target = [\n");
  // END CHECK
  
  // Create the fitting matrix and vestor of target energies.
  Hybrid<double> amat(n_rows * n_columns, "A_matrix");
  Hybrid<double> bvec(n_rows, "b_vector");
  Hybrid<double> xvec(n_columns, "x_vector");
  ScoreCard sc(1, 1, 32);
  for (int i = 0; i < n_targets; i++) {
    const EmulationTarget& itrg = emulcon.getTarget(i);

    // Contribute the complex, ligand, and receptor to the fitting matrix as appropriate.
    std::vector<int> masked_atoms;
    double mm_energy = 0.0;
    const std::vector<StateVariable>& mm_parts = emulcon.getMMContext();
    const size_t n_mm = mm_parts.size();
    if (itrg.getSecondStructureLabel().size() > 0) {
      const int chc_idx = sysche.getCacheIndex(itrg.getFirstStructureLabel(),
                                               itrg.getFirstStructureFrame());
      PhaseSpace ps_i = sysche.getCoordinates(chc_idx);
      PhaseSpaceWriter ps_iw = ps_i.data();
      if (itrg.getFirstStructureMask() == std::string(default_emul_subset_mask)) {
        masked_atoms.resize(1, -1);
      }
      else {
        const CoordinateFrameReader cf_ir(ps_i);
        const AtomMask mask(itrg.getFirstStructureMask(), sysche.getSystemTopologyPointer(chc_idx),
                            sysche.getFeatures(chc_idx), cf_ir);
        masked_atoms = mask.getMaskedAtomList();
      }
      contributeToFit(i, n_rows, &ps_iw, source_id[chc_idx], parameter_map, emulcon, masked_atoms,
                      sysche.getSystemStaticMask(chc_idx), &amat, itrg.getWeight(),
                      FittingContribution::ENERGY);

      // The topology is needed to compute contextual molecular mechanics energies.
      const AtomGraph* ag_i = sysche.getSystemTopologyPointer(chc_idx);
      mm_energy += energyContextForFit(&ps_i, ag_i, sysche.getSystemStaticMask(chc_idx), mm_parts,
                                       ngb_kit, &sc);
    }
    if (itrg.getSecondStructureLabel().size() > 0) {
      const int chc_idx = sysche.getCacheIndex(itrg.getSecondStructureLabel(),
                                               itrg.getSecondStructureFrame());
      PhaseSpace ps_i = sysche.getCoordinates(chc_idx);
      PhaseSpaceWriter ps_iw = ps_i.data();
      if (itrg.getSecondStructureMask() == std::string(default_emul_subset_mask)) {
        masked_atoms.resize(1, -1);
      }
      else {
        const CoordinateFrameReader cf_ir(sysche.getCoordinatePointer(chc_idx));
        const AtomMask mask(itrg.getSecondStructureMask(),
                            sysche.getSystemTopologyPointer(chc_idx), sysche.getFeatures(chc_idx),
                            cf_ir);
        masked_atoms = mask.getMaskedAtomList();
      }
      contributeToFit(i, n_rows, &ps_iw, source_id[chc_idx], parameter_map, emulcon, masked_atoms,
                      sysche.getSystemStaticMask(chc_idx), &amat, -1.0 * itrg.getWeight(),
                      FittingContribution::ENERGY);

      // The topology is needed to compute contextual molecular mechanics energies.
      const AtomGraph* ag_i = sysche.getSystemTopologyPointer(chc_idx);
      sc.initialize();
      mm_energy -= energyContextForFit(&ps_i, ag_i, sysche.getSystemStaticMask(chc_idx), mm_parts,
                                       ngb_kit, &sc);
    }
    if (itrg.getThirdStructureLabel().size() > 0) {
      const int chc_idx = sysche.getCacheIndex(itrg.getThirdStructureLabel(),
                                               itrg.getThirdStructureFrame());
      PhaseSpace ps_i = sysche.getCoordinates(chc_idx);
      PhaseSpaceWriter ps_iw = ps_i.data();
      if (itrg.getThirdStructureMask() == std::string(default_emul_subset_mask)) {
        masked_atoms.resize(1, -1);
      }
      else {
        const CoordinateFrameReader cf_ir(sysche.getCoordinatePointer(chc_idx));
        const AtomMask mask(itrg.getThirdStructureMask(), sysche.getSystemTopologyPointer(chc_idx),
                            sysche.getFeatures(chc_idx), cf_ir);
        masked_atoms = mask.getMaskedAtomList();
      }
      contributeToFit(i, n_rows, &ps_iw, source_id[chc_idx], parameter_map, emulcon, masked_atoms,
                      sysche.getSystemStaticMask(chc_idx), &amat, -1.0 * itrg.getWeight(),
                      FittingContribution::ENERGY);

      // The topology is needed to compute contextual molecular mechanics energies.
      const AtomGraph* ag_i = sysche.getSystemTopologyPointer(chc_idx);
      sc.initialize();
      mm_energy -= energyContextForFit(&ps_i, ag_i, sysche.getSystemStaticMask(chc_idx), mm_parts,
                                       ngb_kit, &sc);
    }

    // CHECK
    printf("  %12.4lf %12.4lf  %s\n", itrg.getTargetEnergy(), mm_energy,
           itrg.getFirstStructureLabel().c_str());
    // END CHECK
    
    bvec.putHost(itrg.getTargetEnergy() - mm_energy, i);
  }

  // CHECK
  printf("];\n");
  printf("There are %zu parameters and %zu targets.\n", n_columns, n_rows);
  const std::vector<double> amat_cpy = amat.readHost();
  const std::vector<double> bvec_cpy = bvec.readHost();
  // END CHECK
  
  // Perform the fitting operation.
#ifdef STORMM_USE_HPC
#else
  // CHECK
#if 0
  printf("A x b = [\n");
  for (int i = 0; i < n_rows; i++) {
    for (int j = 0; j < n_columns; j++) {
      printf(" %12.7lf", amat.data()[(j * n_rows) + i]);
    }
    printf(" %12.7lf\n", bvec.readHost(i));
  }
  printf("];\n");
  printf("\n");
#endif
  // END CHECK
  
  qrSolver(amat.data(), xvec.data(), bvec.data(), n_rows, n_columns);
#endif
  
  // Test the result of the fit.

  // CHECK
#if 0
  printf("Acpy = [\n");
  for (int i = 0; i < n_rows; i++) {
    for (int j = 0; j < n_columns; j++) {
      printf(" %12.7lf", amat_cpy.data()[(j * n_rows) + i]);
    }
    printf("\n");
  }
  printf("];\n");
#endif
  std::vector<double> b_solved(n_rows);
  matrixVectorMultiply(amat_cpy.data(), xvec.data(), b_solved.data(), n_rows, n_columns);
#if 0
  printf("b_solved = [\n");
  for (int i = 0; i < n_rows; i++) {
    printf(" %12.7lf\n", b_solved[i]);
  }
  printf("];\n");
#endif
  printf("Correlation: %9.4lf\n", pearson(bvec_cpy.data(), b_solved.data(), n_rows));
  printf("RMSD:        %9.4lf\n", rmsError(bvec_cpy.data(), b_solved.data(), n_rows));
  // END CHECK
  
  return 0;
}
