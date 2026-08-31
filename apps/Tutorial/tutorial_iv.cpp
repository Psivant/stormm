#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Chemistry/chemical_features.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/statistical_enumerators.h"
#include "../../src/Math/summation.h"
#include "../../src/Namelists/command_line_parser.h"
#include "../../src/Namelists/namelist_emulator.h"
#include "../../src/Namelists/namelist_enumerators.h"
#include "../../src/Numerics/split_fixed_precision.h"
#include "../../src/Potential/local_exclusionmask.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Structure/local_arrangement.h"
#include "../../src/Structure/structure_enumerators.h"
#include "../../src/Synthesis/systemcache.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Topology/atomgraph_enumerators.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/UnitTesting/unit_test_enumerators.h"

using stormm::chemistry::ChemicalFeatures;
using stormm::constants::CaseSensitivity;
using stormm::data_types::int95_t;
#ifndef STORMM_USE_HPC
using stormm::data_types::double2;
#endif
using stormm::diskutil::getBaseName;
using stormm::energy::LocalExclusionMask;
using stormm::energy::LocalExclusionMaskReader;
using stormm::errors::rtErr;
using stormm::namelist::CommandLineParser;
using stormm::namelist::DefaultIsObligatory;
using stormm::namelist::InputRepeats;
using stormm::namelist::NamelistEmulator;
using stormm::namelist::NamelistType;
using stormm::numerics::hostInt95Subtract;
using stormm::numerics::hostInt95ToDouble;
using stormm::parse::strcmpCased;
using stormm::review::stormmWatermark;
using stormm::stmath::sum;
using stormm::stmath::DataOrder;
using stormm::structure::imageCoordinates;
using stormm::structure::ImagingMethod;
using stormm::symbols::charmm_gromacs_bioq;
using stormm::synthesis::AtomGraphSynthesis;
using stormm::synthesis::PhaseSpaceSynthesis;
using stormm::synthesis::PsSynthesisReader;
using stormm::synthesis::SyNonbondedKit;
using stormm::testing::check;
using stormm::testing::RelationalOperator;
using stormm::testing::section;
using stormm::testing::StopWatch;
using stormm::testing::TestEnvironment;
using stormm::testing::TestPriority;
using stormm::testing::TestVerbosity;
using stormm::topology::AtomGraph;
using stormm::topology::NonbondedKit;
using stormm::topology::UnitCellType;
using stormm::trajectory::PhaseSpace;
using stormm::trajectory::PhaseSpaceReader;

//-------------------------------------------------------------------------------------------------
// Update the exclusion lists for interactions already eliminated from any given atom.  Return TRUE
// if the exclusion was already known, FALSE otherwise.
//
// Arguments:
//   excl_counted:  List of lists, one list for each atom containing the known exclusions of that
//                  atom
//   iatom:         The index of the ith atom, expected to be greater than jatom
//   jatom:         Index of the jth atom, expected to be less than iatom
//-------------------------------------------------------------------------------------------------
bool exclusionKnown(std::vector<std::vector<int>> *excl_counted, const int iatom,
                    const int jatom) {
  const size_t ni = excl_counted->at(iatom).size();
  if (ni > 0 && locateValue(excl_counted->at(iatom), jatom, DataOrder::ASCENDING) < ni) {
    return true;
  }
  else {
    excl_counted->at(iatom).push_back(jatom);
    if (ni > 0) {
      std::sort(excl_counted->at(iatom).begin(), excl_counted->at(iatom).end(),
                [](const int a, const int b) { return a < b; });
    }
    excl_counted->at(jatom).push_back(iatom);
    if (excl_counted->at(jatom).size() > 1) {
      std::sort(excl_counted->at(jatom).begin(), excl_counted->at(jatom).end(),
                [](const int a, const int b) { return a < b; });
    }
    return false;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(int argc, const char* argv[]) {

  // Prepare a stopwatch to time some non-bonded calculations
  StopWatch the_clock("STORMM Tutorial IV");
  const int file_read_tm  = the_clock.addCategory("File Reading");
  const int chem_work_tm  = the_clock.addCategory("Chemical Features Detection");
  const int excl_work_tm  = the_clock.addCategory("Non-bonded Exclusion Detection");
  const int basic_nonb_tm = the_clock.addCategory("Basic Non-bonded Evaulation");
  const int clean_nonb_tm = the_clock.addCategory("Non-bonded Evaulation Cleanup");
  const int excl_nonb_tm  = the_clock.addCategory("Excluded Non-bonded Evaulation");
  
  // Prepare for command line input
  CommandLineParser clip("tutorial_iv", "A demonstration of topology and coordinate manipulation "
                         "in STORMM");
  NamelistEmulator *t_nml = clip.getNamelistPointer();
  t_nml->addKeyword("-topol", NamelistType::STRING, std::string(""), DefaultIsObligatory::NO,
                    InputRepeats::YES);
  t_nml->addHelp("-topol", "Select a topology file.  This input may be repeated.\n");
  t_nml->addKeyword("-coord", NamelistType::STRING, std::string(""), DefaultIsObligatory::NO,
                    InputRepeats::YES);
  t_nml->addHelp("-topol", "Select an input coordinate file.  This input may be repeated.\n");
  t_nml->addKeyword("energy_loop", NamelistType::INTEGER, std::to_string(1));
  t_nml->addHelp("energy_loop", "Specify the number of times to loop over the electrostatic "
                 "energy calculation for each system, in order to get a good reading of the "
                 "average time to perform the nested loop over all atoms.");

  // Initialize the testing environment
  TestEnvironment tenv(argc, argv, &clip);
  
  // Parse command line input
  clip.parseUserInput(argc, argv);
  
  // Read in command-line variables
  const int ntop = t_nml->getKeywordEntries("-topol");
  const int ncrd = t_nml->getKeywordEntries("-coord");
  std::vector<std::string> user_defined_topologies;
  user_defined_topologies.reserve(ntop);
  std::vector<std::string> user_defined_coordinates;
  user_defined_coordinates.reserve(ncrd);
  for (int i = 0; i < ntop; i++) {
    user_defined_topologies.push_back(t_nml->getStringValue("-topol", i));
  }
  for (int i = 0; i < ncrd; i++) {
    user_defined_coordinates.push_back(t_nml->getStringValue("-coord", i));
  }
  the_clock.assignTime();

  // Read each topology in series
  std::vector<AtomGraph> agv;
  if (ntop == 0) {
    rtErr("At least one topology must be provided.");
  }
  for (int i = 0; i < ntop; i++) {
    agv.emplace_back(user_defined_topologies[i]);
  }
  std::vector<PhaseSpace> psv;
  for (int i = 0; i < ncrd; i++) {
    psv.emplace_back(user_defined_coordinates[i]);
  }
  check(ntop, RelationalOperator::EQUAL, ncrd, "The number of topologies (" +
        std::to_string(ntop) + ") and coordinate files (" + std::to_string(ncrd) + ") supplied "
        "by the user are inconsistent.", TestPriority::NON_CRITICAL);
  the_clock.assignTime(file_read_tm);

  // Check the number of atoms in each topology, under the order in which they were entered.
  const int failures_pt_a = gbl_test_results.getOverallFailureCount();
  for (int i = 0; i < std::min(ncrd, ntop); i++) {
    check(agv[i].getAtomCount(), RelationalOperator::EQUAL, psv[i].getAtomCount(), "The number of "
          "atoms in supplied topology " + std::to_string(i + 1) + " (" +
          std::to_string(agv[i].getAtomCount()) + ") does not match the number of atoms in the "
          "corresponding coordinate set (" + std::to_string(psv[i].getAtomCount()) + ").",
          TestPriority::NON_CRITICAL);
  }
  const bool all_match = (failures_pt_a == gbl_test_results.getOverallFailureCount());
  the_clock.assignTime();

  // Create ChemicalFeatures objects for each topology.  For a complete analysis of the chirality,
  // all topologies must have matching structures.
  std::vector<ChemicalFeatures> ftv;
  for (int i = 0; i < ntop; i++) {
    if (all_match) {
      ftv.emplace_back(agv[i], psv[i]);
    }
    else {
      ftv.emplace_back(agv[i]);
    }
  }
  the_clock.assignTime(chem_work_tm);

  // Create LocalExclusionMask objects for each topology.  This is the preferred means for
  // detecting whether two atoms share a non-bonded exclusion.
  std::vector<LocalExclusionMask> exv;
  for (int i = 0; i < ntop; i++) {
    exv.emplace_back(agv[i]);
  }
  the_clock.assignTime(excl_work_tm);

  // Run some diagnostics on the various topologies
  printf("Properties of each system:\n\n");
  printf("                               Net      Molecular  Chiral  Aromatic Rotatable\n");
  printf("    System       Atoms Bonds  Charge     Weight    Centers  Groups    Bonds  \n");
  printf("--------------   ----- ----- -------- ------------ ------- -------- ---------\n");
  for (int i = 0; i < ntop; i++) {
    printf("%-14.14s : %5d %5d %7.4lf %12.4lf %7d %8d %9d\n",
           getBaseName(agv[i].getFileName()).c_str(), agv[i].getAtomCount(),
           agv[i].getBondTermCount(), sum<double>(agv[i].getPartialCharge<double>()),
           agv[i].getTotalMass(), ftv[i].getChiralCenterCount(),
           ftv[i].getAromaticGroupCount(), ftv[i].getRotatableBondCount());
  }
  printf("\n");

  // Do some operations with individual atoms
  if (all_match == false || ntop != ncrd) {
    the_clock.printResults();
    rtErr("Exiting due to prior errors.");
  }
  the_clock.assignTime();
  const int nloop = t_nml->getIntValue("energy_loop");
  for (int top_idx = 0; top_idx < ntop; top_idx++) {
    const bool is_periodic = (psv[top_idx].getUnitCellType() != UnitCellType::NONE);
    const NonbondedKit<double> inbk = agv[top_idx].getDoublePrecisionNonbondedKit();
    const PhaseSpaceReader ipsr = psv[top_idx].data();
    const LocalExclusionMaskReader ilmr = exv[top_idx].data();
    for (int loop = 0; loop < nloop; loop++) {

      // Compute the electrostatic interaction of all particles to one another.  Sum the energy.
      double basic_elec_nrg = 0.0;
      for (int i = 1; i < ipsr.natom; i++) {
        const double posx = ipsr.xcrd[i];
        const double posy = ipsr.ycrd[i];
        const double posz = ipsr.zcrd[i];
        const double qi = inbk.charge[i];
        for (int j = 0; j < i; j++) {
          const double qij = inbk.charge[j] * qi;
          double dx = ipsr.xcrd[j] - posx;
          double dy = ipsr.ycrd[j] - posy;
          double dz = ipsr.zcrd[j] - posz;

          // Imaging is only needed if there is a unit cell to consider.  While imageCoordinates()
          // would detect a non-periodic system and return, pre-evaluating the condition allows the
          // inner loop to skip the switch evaluation and, if not-inlined, the function call
          // itself.
          if (is_periodic) {
            imageCoordinates<double, double>(&dx, &dy, &dz, ipsr.umat, ipsr.invu, ipsr.unit_cell,
                                             ImagingMethod::MINIMUM_IMAGE);
          }
          const double dr = sqrt((dx * dx) + (dy * dy) + (dz * dz));
          basic_elec_nrg += qij / dr;
        }
      }
      the_clock.assignTime(basic_nonb_tm);

      // Erase excluded interactions from the basic result.
      std::vector<std::vector<int>> excl_counted(ipsr.natom);
      for (int i = 0; i < ipsr.natom; i++) {
        excl_counted[i].reserve(16);
      }
      the_clock.assignTime();
      for (int i = 0; i < ipsr.natom; i++) {
        for (int j = inbk.nb11_bounds[i]; j < inbk.nb11_bounds[i + 1]; j++) {
          const size_t jatom = inbk.nb11x[j];
          if (exclusionKnown(&excl_counted, i, jatom) == false) {
            const double dx = ipsr.xcrd[jatom] - ipsr.xcrd[i];
            const double dy = ipsr.ycrd[jatom] - ipsr.ycrd[i];
            const double dz = ipsr.zcrd[jatom] - ipsr.zcrd[i];
            const double dr = sqrt((dx * dx) + (dy * dy) + (dz * dz));
            const double qij = inbk.charge[i] * inbk.charge[jatom];
            basic_elec_nrg -= qij / dr;
          }
        }
        for (int j = inbk.nb12_bounds[i]; j < inbk.nb12_bounds[i + 1]; j++) {
          const size_t jatom = inbk.nb12x[j];
          if (exclusionKnown(&excl_counted, i, jatom) == false) {
            const double dx = ipsr.xcrd[jatom] - ipsr.xcrd[i];
            const double dy = ipsr.ycrd[jatom] - ipsr.ycrd[i];
            const double dz = ipsr.zcrd[jatom] - ipsr.zcrd[i];
            const double dr = sqrt((dx * dx) + (dy * dy) + (dz * dz));
            const double qij = inbk.charge[i] * inbk.charge[jatom];
            basic_elec_nrg -= qij / dr;
          }
        }
        for (int j = inbk.nb13_bounds[i]; j < inbk.nb13_bounds[i + 1]; j++) {
          const size_t jatom = inbk.nb13x[j];
          if (exclusionKnown(&excl_counted, i, jatom) == false) {
            const double dx = ipsr.xcrd[jatom] - ipsr.xcrd[i];
            const double dy = ipsr.ycrd[jatom] - ipsr.ycrd[i];
            const double dz = ipsr.zcrd[jatom] - ipsr.zcrd[i];
            const double dr = sqrt((dx * dx) + (dy * dy) + (dz * dz));
            const double qij = inbk.charge[i] * inbk.charge[jatom];
            basic_elec_nrg -= qij / dr;
          }
        }
        for (int j = inbk.nb14_bounds[i]; j < inbk.nb14_bounds[i + 1]; j++) {
          const size_t jatom = inbk.nb14x[j];
          if (exclusionKnown(&excl_counted, i, jatom) == false) {
            const double dx = ipsr.xcrd[jatom] - ipsr.xcrd[i];
            const double dy = ipsr.ycrd[jatom] - ipsr.ycrd[i];
            const double dz = ipsr.zcrd[jatom] - ipsr.zcrd[i];
            const double dr = sqrt((dx * dx) + (dy * dy) + (dz * dz));
            const double qij = inbk.charge[i] * inbk.charge[jatom];
            basic_elec_nrg -= qij / dr;
          }
        }
      }
      basic_elec_nrg *= charmm_gromacs_bioq;
      the_clock.assignTime(clean_nonb_tm);

      // Compute the electrostatic interaction of all particles to one another, omitting all
      // excluded interactions by reference to the LocalExclusionMask.
      double masked_elec_nrg = 0.0;
      for (int i = 1; i < ipsr.natom; i++) {
        const double posx = ipsr.xcrd[i];
        const double posy = ipsr.ycrd[i];
        const double posz = ipsr.zcrd[i];
        const double qi = inbk.charge[i];
        for (int j = 0; j < i; j++) {
          if (! testExclusion(ilmr, j, i)) {
            const double qij = inbk.charge[j] * qi;
            double dx = ipsr.xcrd[j] - posx;
            double dy = ipsr.ycrd[j] - posy;
            double dz = ipsr.zcrd[j] - posz;

            // Imaging is only needed if there is a unit cell to consider.  While
            // imageCoordinates() would detect a non-periodic system and return, pre-evaluating
            // the condition allows the inner loop to skip the switch evaluation and, if
            // not-inlined, the function call itself.
            if (is_periodic) {
              imageCoordinates<double, double>(&dx, &dy, &dz, ipsr.umat, ipsr.invu, ipsr.unit_cell,
                                               ImagingMethod::MINIMUM_IMAGE);
            }
            const double dr = sqrt((dx * dx) + (dy * dy) + (dz * dz));
            masked_elec_nrg += qij / dr;
          }
        }
      }
      masked_elec_nrg *= charmm_gromacs_bioq;
      the_clock.assignTime(excl_nonb_tm);

      // Print the evaluated energy
      if (loop == 0) {
        printf("Electrostatic energy (%14.14s / %14.14s) :\n"
               "    %12.4lf kcal/mol (exclusions erased) vs.\n"
               "    %12.4lf kcal/mol (exclusions omitted)\n",
               getBaseName(agv[top_idx].getFileName()).c_str(),
               getBaseName(psv[top_idx].getFileName()).c_str(), basic_elec_nrg, masked_elec_nrg);
      }
    }
  }
  printf("\n");
  
  // Create a synthesis of the various systems
  std::vector<AtomGraph*> agv_ptr;
  std::vector<int> synth_list;
  for (int i = 0; i < ntop; i++) {
    agv_ptr.push_back(&agv[i]);
    for (int j = 0; j < 2; j++) {
      synth_list.push_back(i);
    }
  }
  PhaseSpaceSynthesis poly_ps(psv, agv_ptr, synth_list);
  AtomGraphSynthesis poly_ag(agv_ptr, synth_list);
  const LocalExclusionMask poly_lem(poly_ag);
  const LocalExclusionMaskReader poly_lemr = poly_lem.data();
  const SyNonbondedKit<double, double2> poly_nbk = poly_ag.getDoublePrecisionNonbondedKit();
  const PsSynthesisReader poly_psr = poly_ps.data();
  std::vector<double> synth_elec_nrg(poly_psr.system_count);
  for (int sys_idx = 1; sys_idx < poly_psr.system_count; sys_idx += 2) {
    const int llim = poly_psr.atom_starts[sys_idx];
    const int hlim = llim + poly_psr.atom_counts[sys_idx];
    const double* umat_ptr = &poly_psr.umat[32 * sys_idx];
    const double* invu_ptr = &poly_psr.invu[32 * sys_idx];
    const bool is_periodic = (poly_psr.unit_cell != UnitCellType::NONE);
    double masked_elec_nrg = 0.0;
    for (int i = llim + 1; i < hlim; i++) {
      for (int j = llim; j < i; j++) {
        if (! testExclusion(poly_lemr, i, j)) {
          const int95_t idx = hostInt95Subtract(poly_psr.xcrd[j], poly_psr.xcrd_ovrf[j],
                                                poly_psr.xcrd[i], poly_psr.xcrd_ovrf[i]);
          const int95_t idy = hostInt95Subtract(poly_psr.ycrd[j], poly_psr.ycrd_ovrf[j],
                                                poly_psr.ycrd[i], poly_psr.ycrd_ovrf[i]);
          const int95_t idz = hostInt95Subtract(poly_psr.zcrd[j], poly_psr.zcrd_ovrf[j],
                                                poly_psr.zcrd[i], poly_psr.zcrd_ovrf[i]);
          double dx = hostInt95ToDouble(idx) * poly_psr.inv_gpos_scale;
          double dy = hostInt95ToDouble(idy) * poly_psr.inv_gpos_scale;
          double dz = hostInt95ToDouble(idz) * poly_psr.inv_gpos_scale;
          if (is_periodic) {
            imageCoordinates<double, double>(&dx, &dy, &dz, umat_ptr, invu_ptr, poly_psr.unit_cell,
                                             ImagingMethod::MINIMUM_IMAGE);
          }
          const double qij = poly_nbk.charge[i] * poly_nbk.charge[j];
          const double dr = sqrt((dx * dx) + (dy * dy) + (dz * dz));
          masked_elec_nrg += qij / dr;
        }
      }
    }
    synth_elec_nrg[sys_idx] = masked_elec_nrg * charmm_gromacs_bioq;
    printf("Electrostatic Energy, Synthesis %14.14s : %12.4lf\n\n",
           getBaseName(poly_ps.getSystemTopologyPointer(sys_idx)->getFileName()).c_str(),
           synth_elec_nrg[sys_idx]);
  }

  // Display results and return success
  the_clock.printResults();
  printTestSummary(tenv.getVerbosity());
  if (tenv.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return 0;
}
