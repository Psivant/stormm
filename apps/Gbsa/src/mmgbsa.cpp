#include <string>
#include <vector>
#include "../../../src/copyright.h"
#include "../../../src/Accelerator/hpc_config.h"
#include "../../../src/Constants/behavior.h"
#include "../../../src/Constants/generalized_born.h"
#include "../../../src/DataTypes/stormm_vector_types.h"
#include "../../../src/FileManagement/file_listing.h"
#include "../../../src/Math/series_ops.h"
#include "../../../src/MoleculeFormat/pdb.h"
#ifdef STORMM_USE_HPC
#  include "../../../src/MolecularMechanics/hpc_minimization.h"
#endif
#include "../../../src/MolecularMechanics/minimization.h"
#include "../../../src/Namelists/command_line_parser.h"
#include "../../../src/Namelists/input_transcript.h"
#include "../../../src/Namelists/namelist_inventory.h"
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
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/Topology/atomgraph_enumerators.h"
#include "../../../src/Trajectory/coordinateframe.h"
#include "../../../src/Trajectory/coordinate_graft.h"
#include "../../../src/UnitTesting/stopwatch.h"
#include "mmgbsa_analysis.h"
#include "mmgbsa_carveout.h"
#include "mmgbsa_problem_set.h"
#include "mmgbsa_testing.h"
#include "mmgbsa_enumerators.h"
#include "nml_mmgbsa.h"

using stormm::testing::StopWatch;
using namespace stormm::card;
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
using namespace stormm::topology;
using namespace stormm::trajectory;
using namespace stormm::synthesis;
using namespace mmgbsa;

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
  CommandLineParser clip("mmgbsa.stormm", "A program for calculating MMGBSA energies of ligands "
                         "based on structures of the individual structures and the complex.");
  clip.addStandardApplicationInputs({ "-i", "-O", "-o", "-except" });
  NamelistEmulator *t_nml = clip.getNamelistPointer();
  t_nml->addKeyword("-nsph", NamelistType::INTEGER, std::to_string(default_sasa_point_count));
  t_nml->addHelp("-nsph", "Set the number of sphere points that will be used to estimate surface "
                 "area");
  t_nml->addKeyword("-rec_p", NamelistType::STRING, std::string(""));
  t_nml->addHelp("-rec_p", "File path to the receptor topology file.  If this structure is not "
                 "provided on the command line, the information will be sought from a -sys "
                 "keyword will be sought from the &files namelist of the input file with the "
                 "label \"receptor\".");
  t_nml->addKeyword("-rec_c", NamelistType::STRING, std::string(""));
  t_nml->addHelp("-rec_c", "File path to the receptor input coordinates file.  If this structure "
                 "is not provided on the command line, the information will be sought from a -sys "
                 "keyword will be sought from the &files namelist of the input file with the "
                 "label \"receptor\".");
  t_nml->addKeyword("-lig_p", NamelistType::STRING, std::string(""));
  t_nml->addHelp("-lig_p", "File path to the ligand topology file");
  t_nml->addKeyword("-lig_c", NamelistType::STRING, std::string(""));
  t_nml->addHelp("-lig_c", "File path to the ligand input coordinates file.  This may contain "
                 "many frames with different poses of the ligand, but must exist in the same "
                 "'coordinate frame as the receptor.");
  t_nml->addKeyword("-unittest", NamelistType::BOOLEAN);
  t_nml->addHelp("-unittest", "Run a series of unit tests to ensure the program's exclusive code "
                 "functions properly.");
  t_nml->addKeyword("-regtest", NamelistType::BOOLEAN);
  t_nml->addHelp("-regtest", "Run a series of regression tests to ensure the program's exclusive "
                 "code functions properly.");
  const std::vector<std::string> my_namelist_names = { "&files", "&mmgbsa", "&minimize",
                                                       "&restraint", "&solvent", "&report",
                                                       "&precision" };
  clip.addControlBlocks(my_namelist_names);
  const std::vector<NamelistToken> mmgbsa_specific_namelists = {
    NamelistToken(std::string("&mmgbsa"), mmgbsaInput)
  };
  clip.addCustomNamelists(mmgbsa_specific_namelists);
  if (displayNamelistHelp(argc, argv, my_namelist_names, mmgbsa_specific_namelists) &&
      clip.doesProgramExitOnHelp()) {
    return 0;
  }
  clip.parseUserInput(argc, argv); 
  
  // If testing is requested, perform it and exit.
  if (t_nml->getBoolValue("-unittest")) {
    return runUnitTests();
  }
  if (t_nml->getBoolValue("-regtest")) {
    return runRegressionTests();
  }
  
  // Take in the user input from the input file.  Immediately test for valid inputs.
  const UserSettings ui(clip, { "-pe", "-ce", "-rg" });
  FilesControls ficon = ui.getFilesNamelistInfo();
  const ExceptionResponse policy = translateExceptionResponse(t_nml->getStringValue("-except"));
  const PrintSituation prnt_protocol = t_nml->getBoolValue("-O") ? PrintSituation::OVERWRITE :
                                                                   PrintSituation::OPEN_NEW;
  int start_line = 0;
  bool gbsa_found = false;
  MMGBSAControls gbsacon(ui.getInputFileName(), &start_line, &gbsa_found,
                         ui.getExceptionBehavior(), WrapTextSearch::YES);
  
  // Search for a receptor among the various structures.  This will make a copy of the &files
  // namelist in the UserSettings class object, but the copy will be modifiable in the context of
  // this program.  Check for input errors.
  const bool cli_top = (t_nml->getKeywordStatus("-rec_p") == InputStatus::USER_SPECIFIED);
  const bool cli_crd = (t_nml->getKeywordStatus("-rec_c") == InputStatus::USER_SPECIFIED);
  if (cli_top && cli_crd) {
    MoleculeSystem recmol(t_nml->getStringValue("-rec_p"), t_nml->getStringValue("-rec_c"), "", "",
                          "receptor", 0, 1, 1, CoordinateFileKind::UNKNOWN,
                          CoordinateFileKind::UNKNOWN, CoordinateFileKind::UNKNOWN);
    ficon.addSystem(recmol);
  }
  else if ((cli_top && cli_crd == false) || (cli_top == false && cli_crd)) {
    rtErr("A receptor " +
          ((cli_top) ? std::string("topology") : std::string("input coordinates")) + " file was "
          "provided but the corresponding " +
          ((cli_top) ? std::string("input_coordinates") : std::string("topology")) + " was not.",
          "mmgbsa");
  }

  // Load all of the systems
  SystemCache sysche(ficon, policy, MapRotatableGroups::YES, prnt_protocol, &master_timer);
  if (sysche.getFirstMatchingSystemIndex("receptor", ExceptionResponse::SILENT) ==
      sysche.getSystemCount()) {
    rtErr("No receptor was specified.  The receptor is designated by adding \"-label receptor\" "
          "to its -sys { ... } keyword / value group instance.", "mmgbsa");
  }
  if (sysche.getFirstMatchingSystemIndex("ligand", ExceptionResponse::SILENT) ==
      sysche.getSystemCount()) {
    rtErr("No ligands were specified.  A ligand is designated by adding \"-label ligand\" to its "
          "-sys { ... } keyword / value group instance.", "mmgbsa");
  }
  
  // Find the receptor among all of the different structures in the &files namelist (including
  // command line edits).
  const std::vector<int> rec_cache_idx = sysche.getMatchingSystemIndices("receptor");
  const std::vector<int> lig_cache_idx = sysche.getMatchingSystemIndices("ligand");
  
  // Combine the topologies of each ligand with that of the receptor.  The goal is to make a list
  // of the supplied ligand poses with an index into the system cache's topology list for the
  // ligand's own topology (the "x" member of the tuple), and index into the system cache's
  // topology list for the corresponding receptor topology (the "y" member of the tuple), and the
  SolventControls solvcon = ui.getSolventNamelistInfo();
  PrecisionControls preccon = ui.getPrecisionNamelistInfo();
  MMGBSAProblemSet sandbox(sysche.getSelfPointer());
  const int gpos_bits = preccon.getGlobalPosScalingBits();
  const int vel_bits = preccon.getVelocityScalingBits();
  const int frc_bits = preccon.getForceScalingBits();
#ifdef STORMM_USE_HPC
  AtomGraphSynthesis poly_ag = sandbox.exportTopologySynthesis(gpu, policy);
  PhaseSpaceSynthesis poly_ps = sandbox.exportCoordinateSynthesis(gpu, HybridFormat::EXPEDITED,
                                                                  gpos_bits, vel_bits, frc_bits);
#else
  AtomGraphSynthesis poly_ag = sandbox.exportTopologySynthesis();
  PhaseSpaceSynthesis poly_ps = sandbox.exportCoordinateSynthesis(null_gpu,
                                                                  HybridFormat::HOST_ONLY,
                                                                  gpos_bits, vel_bits, frc_bits);
  const std::vector<AtomGraph*>& sandbox_ag_v = sandbox.getUniqueTopologies();
  std::vector<StaticExclusionMask> sem_v;
  for (size_t i = 0; i < sandbox_ag_v.size(); i++) {
    sem_v.emplace_back(sandbox_ag_v[i]);
  }
#endif
  StaticExclusionMaskSynthesis poly_sems(sandbox.getUniqueTopologies(),
                                         sandbox.getMasterTopologyIndices());
  const NeckGeneralizedBornTable ngb_tab;
  poly_ag.setImplicitSolventModel(solvcon.getImplicitSolventModel(), ngb_tab);
  InitializationTask init_task;
  switch (poly_ag.getImplicitSolventModel()) {
  case ImplicitSolventModel::NONE:
    init_task = InitializationTask::GENERAL_MINIMIZATION;
    break;
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    init_task = InitializationTask::GB_MINIMIZATION;
    break;
  }
  poly_ag.loadNonbondedWorkUnits(poly_sems, init_task, 0, gpu);
#ifdef STORMM_USE_HPC
  poly_ag.upload();
  poly_sems.upload();
  poly_ps.upload();
#endif

  // Issue a call to the API, to perform energy minimization
  const MinimizeControls mincon = ui.getMinimizeNamelistInfo();
  const PrecisionModel general_prec = (preccon.getValenceMethod() == PrecisionModel::DOUBLE ||
                                       preccon.getNonbondedMethod() == PrecisionModel::DOUBLE) ?
                                      PrecisionModel::DOUBLE : PrecisionModel::SINGLE;
#ifdef STORMM_USE_HPC
  ScoreCard sc = launchMinimization(poly_ag, poly_sems, &poly_ps, mincon, gpu, general_prec, 32,
                                    &master_timer, nullptr,
                                    "Receptor, ligand, and complex minimization");
  poly_ps.download();
  sc.download();
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    PhaseSpace ps_i = poly_ps.exportSystem(i);
    const AtomGraph* ag_i = poly_ag.getSystemTopologyPointer(i);
    const std::vector<int> all_mols = incrementingSeries<int>(0, ag_i->getMoleculeCount());
    const double sa_ri = surfaceArea(&ps_i, ag_i, all_mols, 272, default_sasa_probe_radius,
                                     default_sasa_energy, EvaluateForce::NO,
                                     SasaReference::LENNARD_JONES_SIGMA, PrecisionModel::DOUBLE);
    const llint sa_fpi = sa_ri * sc.getEnergyScalingFactor<double>();
    sc.contribute(StateVariable::SURFACE_AREA, sa_fpi, i);
  }
#else
  ScoreCard sc(poly_ps.getSystemCount(),
               (mincon.getTotalCycles() + mincon.getDiagnosticPrintFrequency() - 1) /
               mincon.getDiagnosticPrintFrequency());
  for (int i = 0; i < poly_ps.getSystemCount(); i++) {
    PhaseSpace ps_i = poly_ps.exportSystem(i);
    const AtomGraph* ag_i = poly_ag.getSystemTopologyPointer(i);
    const int problem_topl_idx = sandbox.getMasterTopologyIndex(i);
    ScoreCard sc_i = minimize(&ps_i, ag_i, sem_v[problem_topl_idx], mincon,
                              sc.getEnergyScaleBits());
    const std::vector<int> all_mols = incrementingSeries<int>(0, ag_i->getMoleculeCount());
    const double sa_ri = surfaceArea(&ps_i, ag_i, all_mols, 272, default_sasa_probe_radius,
                                     default_sasa_energy, EvaluateForce::NO,
                                     SasaReference::LENNARD_JONES_SIGMA, PrecisionModel::DOUBLE);
    const llint sa_fpi = sa_ri * sc.getEnergyScalingFactor<double>();
    sc_i.contribute(StateVariable::SURFACE_AREA, sa_fpi, i);
    sc.importCard(sc_i, i, 0);
  }
#endif
  sc.computeTotalEnergy();
  const std::vector<double> sc_totals = sc.reportInstantaneousStates(StateVariable::TOTAL_ENERGY);

  // Create the report file.
  ReportControls repcon = ui.getReportNamelistInfo();
  std::vector<SectionContents> outp;
  addLigandTopologySummary(sandbox, ficon, repcon, &outp);
  outp.emplace_back("MM-GBSA General Results", ficon.getReportFile());
  outp.back().addNarration("Total ligand poses: " + std::to_string(sandbox.getComplexCount()) +
                           "\nTotal unique ligands: " +
                           std::to_string((poly_ps.getUniqueTopologyCount() - 1) / 2));
  const int n_complex = sandbox.getComplexCount();
  const int n_receptor = sandbox.getReceptorStates();

  // Print the bound structures, if requested.
  if (gbsacon.printFinalStructures()) {
    
    // Determine the receptor carveout
    const AtomGraph *receptor_ag = sysche.getTopologiesMatchingLabel("receptor")[0];
    const std::vector<int> binding_site_atoms = receptorCarveOut(gbsacon, sysche);
    const int binding_site_natom = binding_site_atoms.size();
    const int receptor_natom = receptor_ag->getAtomCount();
    AtomGraph binding_site_ag(*receptor_ag, binding_site_atoms);
    binding_site_ag.setImplicitSolventModel(receptor_ag->getImplicitSolventModel(),
                                            receptor_ag->getDielectricConstant(),
                                            receptor_ag->getSaltConcentration(),
                                            translateAtomicRadiusSet(receptor_ag->getPBRadiiSet(),
                                                                     ui.getExceptionBehavior()),
                                            ui.getExceptionBehavior());
    const PsSynthesisReader poly_psr = poly_ps.data();
    switch (gbsacon.getEnergyReportDepth()) {
    case EnergyReport::COMPLETE:
      for (int i = n_receptor + n_complex; i < poly_ps.getSystemCount(); i++) {

        // Export the coordinates as a new[]-allocated CoordinateFrame (avoid the overhead of
        // extra GPU memory allocations when they will not be needed).  Then, create a reduced
        // set of coordinates from the frame along with a topology to match.  The mapping that
        // follows assumes the receptor will lie before the ligand in any topology.
        const CoordinateFrame cf_complex = poly_ps.exportCoordinates(i, HybridFormat::HOST_ONLY);
        const int cmpx_natom = cf_complex.getAtomCount();
        const int lgnd_natom = poly_psr.atom_counts[i - n_complex];
        const int reduced_natom = binding_site_natom + lgnd_natom;
        CoordinateFrame cf_reduced(reduced_natom);
        std::vector<int2> reduction_mapping(reduced_natom);
        for (int j = 0; j < binding_site_natom; j++) {
          reduction_mapping[j] = { j, binding_site_atoms[j] };
        }
        for (int j = 0; j < lgnd_natom; j++) {
          const int j_offset = binding_site_natom + j;
          reduction_mapping[j_offset] = { j_offset, receptor_natom + j };
        }
        const CoordinateFrameReader cfr_complex = cf_complex.data();
        CoordinateFrameWriter cfw_reduced = cf_reduced.data();
        coordGraft(&cfw_reduced, cfr_complex, reduction_mapping.data(), reduced_natom);
        const AtomGraph *lig_topl = poly_ps.getSystemTopologyPointer(i - n_complex);
        AtomGraph ag_reduced(binding_site_ag, *lig_topl, MoleculeOrdering::RETAIN_ORDER);
        const ChemicalDetailsKit cdk_red = ag_reduced.getChemicalDetailsKit();
        switch (gbsacon.getPrintedStructureFormat()) {
        case PrintedPoseFormat::AMBER:
          {
            const std::string fname_base(gbsacon.getOutputStructureBase() + "_" +
                                         getBaseName(lig_topl->getFileName()) + "_" +
                                         std::to_string(i - n_complex - n_receptor));
            ag_reduced.printToFile(fname_base + std::string(".parm7"), TopologyKind::AMBER,
                                   ui.getPrintingPolicy());
            cf_reduced.exportToFile(fname_base + std::string(".inpcrd"),
                                    CoordinateFileKind::AMBER_INPCRD, ui.getPrintingPolicy());
          }
          break;
        case PrintedPoseFormat::PDB:
          {
            Pdb rcsb_i(ag_reduced, cf_reduced);
            rcsb_i.writeToFile(gbsacon.getOutputStructureBase() + "_" +
                               getBaseName(lig_topl->getFileName()) + "_" +
                               std::to_string(i - n_complex - n_receptor) + std::string(".pdb"),
                               ui.getPrintingPolicy());
          }
          break;
        }
      }
      break;
    case EnergyReport::COMPOUND_AVERAGES:
      break;
    case EnergyReport::BEST_RESULTS:
      {
        std::vector<bool> coverage(poly_ps.getSystemCount(), false);
        for (int i = n_receptor; i < n_receptor + n_complex; i++) {

          // Find the ligand topology index, then all examples of the ligand within the data set.
          // Select the bext pose for the ligand to display.
          if (coverage[i]) {
            continue;
          }
          const AtomGraph *lig_topl = poly_ps.getSystemTopologyPointer(i);
          const int topl_idx = poly_ps.getUniqueTopologyIndex(i);
          const std::vector<int> all_exi = poly_ps.getSystemIndicesByTopology(topl_idx);
          const int n_alternates = all_exi.size();
          double best_nrg;
          int best_idx;
          for (int j = 0; j < n_alternates; j++) {

            // Compute the binding energy for the complex--this may be repeated in a subsequent
            // reporting step but is trivial.  The receptor energy can be taken as any of the
            // values for different receptor poses--the ranking will not change.
            const double cmplx_nrg = sc_totals[all_exi[j] + n_complex];
            const double rcptr_nrg = sc_totals[0];
            const double lgnd_nrg  = sc_totals[all_exi[j]];
            if (j == 0 || cmplx_nrg - rcptr_nrg - lgnd_nrg <= best_nrg) {
              best_nrg = cmplx_nrg - rcptr_nrg - lgnd_nrg;
              best_idx = all_exi[j];
            }
          }
          const CoordinateFrame cf_complex = poly_ps.exportCoordinates(best_idx + n_complex,
                                                                       HybridFormat::HOST_ONLY);
          const int lgnd_natom = poly_psr.atom_counts[best_idx];
          const int reduced_natom = binding_site_natom + lgnd_natom;
          CoordinateFrame cf_reduced(reduced_natom);
          std::vector<int2> reduction_mapping(reduced_natom);
          for (int j = 0; j < binding_site_natom; j++) {
            reduction_mapping[j] = { j, binding_site_atoms[j] };
          }
          for (int j = 0; j < lgnd_natom; j++) {
            const int j_offset = binding_site_natom + j;
            reduction_mapping[j_offset] = { j_offset, receptor_natom + j };
          }
          const CoordinateFrameReader cfr_complex = cf_complex.data();
          CoordinateFrameWriter cfw_reduced = cf_reduced.data();
          coordGraft(&cfw_reduced, cfr_complex, reduction_mapping.data(), reduced_natom);
          AtomGraph ag_reduced(binding_site_ag, *lig_topl, MoleculeOrdering::RETAIN_ORDER);
          const ChemicalDetailsKit cdk_red = ag_reduced.getChemicalDetailsKit();
          switch (gbsacon.getPrintedStructureFormat()) {
          case PrintedPoseFormat::AMBER:
            {
              const std::string fname_base(gbsacon.getOutputStructureBase() + "_" +
                                           getBaseName(lig_topl->getFileName()) + "_" +
                                           std::to_string(i - n_receptor));
              ag_reduced.printToFile(fname_base + std::string(".parm7"), TopologyKind::AMBER,
                                     ui.getPrintingPolicy());
              cf_reduced.exportToFile(fname_base + std::string(".inpcrd"),
                                      CoordinateFileKind::AMBER_INPCRD, ui.getPrintingPolicy());
            }
            break;
          case PrintedPoseFormat::PDB:
            {
              Pdb rcsb_i(ag_reduced, cf_reduced);
              rcsb_i.writeToFile(gbsacon.getOutputStructureBase() + "_" +
                                 getBaseName(lig_topl->getFileName()) + "_" +
                                 std::to_string(i - n_receptor) + std::string(".pdb"),
                                 PrintSituation::OVERWRITE);
            }
            break;
          }
        }
      }
      break;
    }
  }
  
  // Compute the mean receptor energy.  This is only valid if the receptor's apo state is
  // energy-minimized to calculate an energy.  If receptor structures are taken from the context
  // of a complex, then the receptor state corresponding to each ligand provides dg_receptor.
  double dg_receptor;
  switch (gbsacon.getWeightingScheme()) {
  case WeightingScheme::FLAT:
    dg_receptor = mean<double>(sc_totals.data(), sandbox.getReceptorStates());
    break;
  case WeightingScheme::BOLTZMANN:
    dg_receptor = boltzmannWeightedMean<double>(sc_totals.data(), n_receptor,
                                                gbsacon.getTemperature());
    break;
  }

  // Prepare to recalculate the energies of receptors taken from the context of each complex.
  ScoreCard sc_rcpt_from_cmpx(sandbox.getComplexCount(), 1, 32);
  if (gbsacon.takeReceptorFromComplex()) {
    const PsSynthesisReader poly_psr = poly_ps.data();

    // One by one, take the receptor structure from each complex, recompute the energies, and
    // replace the energy values in the ScoreCard.
    const int natom_rcpt = poly_ps.getAtomCount(0);
    Hybrid<int2> rcpt_map(natom_rcpt, 0);
    int2* rcpt_map_ptr = rcpt_map.data();
    for (int i = 0; i < sandbox.getComplexCount(); i++) {
      PhaseSpace tmp_receptor(natom_rcpt, HybridFormat::HOST_ONLY);
      CoordinateFrameWriter tmp_receptor_cfr(&tmp_receptor);
      const int synth_idx = i + sandbox.getReceptorStates() + sandbox.getComplexCount();
      const int cmpx_start = poly_ps.getAtomOffset(synth_idx);
      for (int j = 0; j < natom_rcpt; j++) {
        rcpt_map_ptr[j] = { j, cmpx_start + j };
      }
      const int rcpt_cache_idx = sandbox.getReceptorTopologyCacheIndex();
      coordGraft(&tmp_receptor_cfr, poly_psr, rcpt_map_ptr, natom_rcpt);
      const AtomGraph *ag_i = sysche.getSystemTopologyPointer(rcpt_cache_idx);
      evalRestrainedMMGB(&tmp_receptor, &sc_rcpt_from_cmpx, ag_i, ngb_tab,
                         sysche.getSystemStaticMask(rcpt_cache_idx),
                         sysche.getRestraintPointer(rcpt_cache_idx), EvaluateForce::NO, i, 0);
      const std::vector<int> all_mols = incrementingSeries<int>(0, ag_i->getMoleculeCount());
      const double sa_ri = surfaceArea(&tmp_receptor, ag_i, all_mols, 272,
                                       default_sasa_probe_radius, default_sasa_energy,
                                       EvaluateForce::NO, SasaReference::LENNARD_JONES_SIGMA,
                                       PrecisionModel::DOUBLE);
      const llint sa_fpi = sa_ri * sc.getEnergyScalingFactor<double>();
      sc_rcpt_from_cmpx.contribute(StateVariable::SURFACE_AREA, sa_fpi, i);
    }
  }
  
  // Switch the output mode to "complete" if there is only one example of each ligand in the data
  // set.
  if (poly_ps.getUniqueTopologyCount() == poly_ps.getSystemCount()) {
    gbsacon.setEnergyReportDepth(EnergyReport::COMPLETE);
  }
  const std::vector<AtomGraph*>& unique_ag = poly_ps.getUniqueTopologies();
  const std::vector<int> unique_ag_examples = poly_ps.getUniqueTopologyExampleIndices();
  const int n_unique_topol = unique_ag.size();
  switch (gbsacon.getEnergyReportDepth()) {
  case EnergyReport::COMPLETE:
    {
      // Print the complex, ligand, receptor, and difference in total binding energies
      std::vector<double> dg_vector(4 * n_complex);
      for (int i = 0; i < n_complex; i++) {
        dg_vector[i] = sc_totals[i + n_receptor + n_complex];
        dg_vector[i + n_complex] = sc_totals[i + n_receptor];
        dg_vector[i + (3 * n_complex)] = dg_vector[i] - dg_vector[i + n_complex] - dg_receptor;
      }
      if (gbsacon.takeReceptorFromComplex()) {
        sc_rcpt_from_cmpx.computeTotalEnergy();
        for (int i = 0; i < n_complex; i++) {
          dg_vector[i + (2 * n_complex)] =
            sc_rcpt_from_cmpx.reportInstantaneousStates(StateVariable::TOTAL_ENERGY, i);
        }
      }
      else {
        for (int i = 0; i < n_complex; i++) {
          dg_vector[i + (2 * n_complex)] = dg_receptor;
        }
      }
      std::vector<std::string> dgs_vector(5 * n_complex);
      for (int i = 0; i < 4 * n_complex; i++) {
        dgs_vector[i] = realToString(dg_vector[i], 2, NumberFormat::STANDARD_REAL);
      }
      std::vector<std::string> base_topl_names(n_complex);
      for (int i = 0; i < n_complex; i++) {
        const AtomGraph* ag_i = poly_ps.getSystemTopologyPointer(i + 1);
        base_topl_names[i] = ag_i->getFileName();
      }
      const std::vector<std::string> common_paths = extractCommonPaths(&base_topl_names, 8, 2);
      const char comment_symbol = commentSymbol(repcon.getOutputSyntax());
      for (int i = 4 * n_complex; i < 5 * n_complex; i++) {
        dgs_vector[i] = std::string(1, comment_symbol) + " " +
                        base_topl_names[i - (4 * n_complex)];
      }
      const std::vector<std::string> headers = { std::string("Complex"), std::string("Ligand"),
                                                 std::string("Receptor"),
                                                 std::string("ddG Binding"),
                                                 std::string("Structure") };
      const std::vector<JustifyText> justifications = { JustifyText::RIGHT, JustifyText::RIGHT,
                                                        JustifyText::RIGHT, JustifyText::RIGHT,
                                                        JustifyText::LEFT };
      ReportTable ddg_tab(dgs_vector, headers, "ddG", repcon.getReportFileWidth(), justifications);
      ddg_tab.unprotectContent();
      outp.back().addNarration("The energy of the relaxed ligand and receptor structures can be "
                               "subtracted from the relaxed energy of the complex in order to "
                               "obtain an estimate of the binding energy.  Results below are "
                               "presented for all systems in kcal/mol.");
      if (common_paths.size() > 0) {
        outp.back().addNarration("In the table that follows, file names are abbreviated.  The "
                                 "following key illustrates the paths:\n" +
                                 listCommonPaths(common_paths));
      }
      outp.back().addTable(ddg_tab);
      
      // Print the decomposition
      if (gbsacon.printEnergyDecomposition() ||
          repcon.getTranscript().getKeywordStatus("energy") != InputStatus::MISSING) {
        outp.emplace_back("Molecular Mechanics Energy Decomposition", ficon.getReportFile());
        outp.back().addNarration("The molecular mechanics energy can be split into a handful of "
                                 "intuitive quantities (albeit some without a physical "
                                 "foundation).  The following tables present the energy "
                                 "components of each receptor-ligand interaction and their "
                                 "changes upon binding.");
        std::vector<StateVariable> nrg_dcmp;
        if (repcon.getTranscript().getKeywordStatus("energy") == InputStatus::MISSING) {
          nrg_dcmp = { StateVariable::BOND, StateVariable::ANGLE, StateVariable::PROPER_DIHEDRAL,
                       StateVariable::IMPROPER_DIHEDRAL, StateVariable::UREY_BRADLEY,
                       StateVariable::CHARMM_IMPROPER, StateVariable::CMAP, StateVariable::VDW,
                       StateVariable::VDW_ONE_FOUR, StateVariable::ELECTROSTATIC,
                       StateVariable::ELEC_ONE_FOUR, StateVariable::GENERALIZED_BORN,
                       StateVariable::RESTRAINT, StateVariable::SURFACE_AREA };
        }
        else {
          nrg_dcmp = repcon.getReportedQuantities();
        }
        std::vector<std::string> eparts_vector(n_complex * 5);
        for (size_t i = 0; i < nrg_dcmp.size(); i++) {
          const std::vector<double> estt = sc.reportInstantaneousStates(nrg_dcmp[i]);
          for (int j = 0; j < n_complex; j++) {
            eparts_vector[j                  ] = realToString(estt[1 + n_complex + j], 2,
                                                              NumberFormat::STANDARD_REAL);
            eparts_vector[j +       n_complex] = realToString(estt[1 + j], 2,
                                                              NumberFormat::STANDARD_REAL);
            const double j_ddg = estt[1 + n_complex + j] - estt[1 + j] - estt[0];
            eparts_vector[j + (3 * n_complex)] = realToString(j_ddg, 2,
                                                              NumberFormat::STANDARD_REAL);
            eparts_vector[j + (4 * n_complex)] = std::string(1, comment_symbol) + " " +
                                                 base_topl_names[j];
          }
          if (gbsacon.takeReceptorFromComplex()) {
            const std::vector<double> estt_rec =
              sc_rcpt_from_cmpx.reportInstantaneousStates(nrg_dcmp[i]);
            for (int j = 0; j < n_complex; j++) {
              eparts_vector[j + (2 * n_complex)] = realToString(estt_rec[j], 2,
                                                                NumberFormat::STANDARD_REAL);
            }
          }
          else {
            for (int j = 0; j < n_complex; j++) {
              eparts_vector[j + (2 * n_complex)] = realToString(estt[0], 2,
                                                                NumberFormat::STANDARD_REAL);
            }
          }
          ReportTable edcmp_tab(eparts_vector, headers,
                                std::string("ddg_") + lowercase(getEnumerationName(nrg_dcmp[i])),
                                repcon.getReportFileWidth(), justifications);
          edcmp_tab.unprotectContent();
          bool print_table = false;
          switch (nrg_dcmp[i]) {
          case StateVariable::BOND:
          case StateVariable::ANGLE:
          case StateVariable::PROPER_DIHEDRAL:
          case StateVariable::IMPROPER_DIHEDRAL:
          case StateVariable::VDW:
          case StateVariable::VDW_ONE_FOUR:
          case StateVariable::ELECTROSTATIC:
          case StateVariable::ELEC_ONE_FOUR:
          case StateVariable::GENERALIZED_BORN:            
          case StateVariable::SURFACE_AREA:
            print_table = true;
            break;
          case StateVariable::RESTRAINT:
            print_table = true;
            break;
          case StateVariable::UREY_BRADLEY:
          case StateVariable::CHARMM_IMPROPER:
          case StateVariable::CMAP:
            for (int i = 0; i < n_complex; i++) {
              const AtomGraph* ag_i = poly_ps.getSystemTopologyPointer(i + n_complex + 1);
              if (ag_i->getUreyBradleyTermCount() > 0 || ag_i->getCharmmImprTermCount() > 0 ||
                  ag_i->getCmapTermCount() > 0) {
                print_table = true;
              }
            }
            break;
          case StateVariable::KINETIC:
          case StateVariable::PRESSURE:
          case StateVariable::VIRIAL_11:
          case StateVariable::VIRIAL_12:
          case StateVariable::VIRIAL_13:
          case StateVariable::VIRIAL_22:
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
          if (print_table) {
            outp.back().addTable(edcmp_tab);
          }
        }
      }
    }
    break;
  case EnergyReport::COMPOUND_AVERAGES:
    {
      const int n_ligand = (poly_ps.getUniqueTopologyCount() - n_receptor) / 2;
      std::vector<bool> coverage(poly_ps.getSystemCount(), false);
      std::vector<double> dg_vector(4 * n_ligand);
      int lidx = 0;
      std::vector<int> snapshot_counts(n_ligand);
      for (int i = n_receptor; i < n_receptor + n_complex; i++) {
        if (coverage[i]) {
          continue;
        }

        // Get the unique topology index of the ligand, then all examples of the ligand within the
        // synthesis.  Loop over those examples and assemble averages for the ligand in isolation
        // and the corresponding complexes.
        const int topl_idx = poly_ps.getUniqueTopologyIndex(i);
        const std::vector<int> all_exi = poly_ps.getSystemIndicesByTopology(topl_idx);
        const int n_exi = all_exi.size();
        std::vector<double> iligand_nrg(n_exi), icomplex_nrg(n_exi);
        for (int j = 0; j < n_exi; j++) {
          iligand_nrg[j] = sc_totals[all_exi[j]];
          icomplex_nrg[j] = sc_totals[all_exi[j] + n_complex];
          coverage[all_exi[j]] = true;
          coverage[all_exi[j] + n_complex] = true;
        }
        double dg_ligand, dg_complex;
        switch (gbsacon.getWeightingScheme()) {
        case WeightingScheme::FLAT:
          dg_ligand  = mean(iligand_nrg);
          dg_complex = mean(icomplex_nrg);
          break;
        case WeightingScheme::BOLTZMANN:
          dg_ligand  = boltzmannWeightedMean(iligand_nrg,  gbsacon.getTemperature());
          dg_complex = boltzmannWeightedMean(icomplex_nrg, gbsacon.getTemperature());
          break;
        }
        dg_vector[lidx] = dg_complex;
        dg_vector[lidx + n_ligand] = dg_ligand;
        dg_vector[lidx + (2 * n_ligand)] = dg_receptor;
        dg_vector[lidx + (3 * n_ligand)] = dg_complex - dg_ligand - dg_receptor;
        snapshot_counts[lidx] = n_exi;
        lidx++;
      }
      std::vector<std::string> dgs_vector(6 * n_ligand);
      for (int j = 0; j < 4 * n_ligand; j++) {
        dgs_vector[j] = realToString(dg_vector[j], 2, NumberFormat::STANDARD_REAL);
      }
      for (int j = 4 * n_ligand; j < 5 * n_ligand; j++) {
        dgs_vector[j] = std::to_string(snapshot_counts[j - (4 * n_ligand)]);
      }
      std::vector<std::string> base_topl_names(n_ligand);
      for (int j = 0; j < n_ligand; j++) {
        const AtomGraph* ag_j = poly_ps.getSystemTopologyPointer(unique_ag_examples[j] + 1);
        base_topl_names[j] = ag_j->getFileName();
      }
      const std::vector<std::string> common_paths = extractCommonPaths(&base_topl_names, 8, 2);
      const char comment_symbol = commentSymbol(repcon.getOutputSyntax());
      for (int j = 5 * n_ligand; j < 6 * n_ligand; j++) {
        dgs_vector[j] = std::string(1, comment_symbol) + " " +
                        base_topl_names[j - (5 * n_ligand)];
      }
      ReportTable ddg_tab(dgs_vector, { std::string("Complex"), std::string("Ligand"),
                                        std::string("Receptor"), std::string("ddG Binding"),
                                        std::string("Snapshots"), std::string("Structure") },
                          "ddG", repcon.getReportFileWidth(),
                          { JustifyText::RIGHT, JustifyText::RIGHT, JustifyText::RIGHT,
                            JustifyText::RIGHT, JustifyText::RIGHT, JustifyText::LEFT });
      ddg_tab.unprotectContent();
      outp.back().addNarration("The energy of the relaxed ligand and receptor structures can be "
                               "subtracted from the relaxed energy of the complex in order to "
                               "obtain an estimate of the binding energy.  Results below are "
                               "presented for all systems in kcal/mol, after averaging using " +
                               getEnumerationName(gbsacon.getWeightingScheme()) + " weights.");
      if (common_paths.size() > 0) {
        outp.back().addNarration("In the table that follows, file names are abbreviated.  The "
                                 "following key illustrates the paths:\n" +
                                 listCommonPaths(common_paths));
      }
      outp.back().addTable(ddg_tab);
    }
    break;
  case EnergyReport::BEST_RESULTS:
    break;
  }
  printAllSections(ficon.getReportFile(), ui.getPrintingPolicy(), outp, repcon.getOutputSyntax(),
                   ListEnumeration::NUMBERED, ListEnumeration::NUMBERED,
                   repcon.getReportFileWidth());

  return 0;
}
