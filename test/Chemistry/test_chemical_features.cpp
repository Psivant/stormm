#include "copyright.h"
#include "../../src/Chemistry/atommask.h"
#include "../../src/Chemistry/atom_equivalence.h"
#include "../../src/Chemistry/chemical_features.h"
#include "../../src/Chemistry/chemistry_enumerators.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Parsing/polynumeric.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Topology/atomgraph_enumerators.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/Trajectory/trajectory_enumerators.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::constants::ExceptionResponse;
using stormm::symbols::amber_ancient_bioq;
#ifndef STORMM_USE_HPC
using stormm::data_types::int2;
using stormm::data_types::char4;
#endif
using stormm::data_types::ullint;
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getBaseName;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::osSeparator;
using stormm::errors::rtWarn;
using stormm::stmath::maxValue;
using stormm::parse::polyNumericVector;
using stormm::parse::NumberFormat;
using stormm::data_types::operator!=;
using stormm::data_types::operator==;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using stormm::topology::AtomGraph;
using stormm::topology::ChemicalDetailsKit;
using stormm::topology::MobilitySetting;
using stormm::trajectory::CoordinateFileKind;
using stormm::trajectory::PhaseSpace;

using namespace stormm::chemistry;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Lay out the testing environment
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
  StopWatch timer;

  // Section 1: chemical bonding structures
  section("Chemical bonding patterns");

  // Section 2: formal charges and bond orders
  section("Resonance structures");

  // Section 3: object manipulation
  section("A first-class C++ object");

  // Section 4: test bond pattern matching (a free function critical to the AtomEquivalence oject)
  section("Bond pattern matching");
  
  // Test the existence of topology and coordinate files
  const char osc = osSeparator();
  const std::string base_chem_name = oe.getStormmSourcePath() + osc + "test" + osc + "Chemistry";
  const std::string base_nml_name  = oe.getStormmSourcePath() + osc + "test" + osc + "Namelists";
  const std::string base_nmlp_name = base_nml_name + osc + "topol";
  const std::string base_nmlc_name = base_nml_name + osc + "coord";
  const std::string base_top_name  = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_name  = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";

  // Ligands with an eight-membered ring that will avoid anti-aromaticity by puckering
  const std::string mol1_top_name  = base_chem_name + osc + "lig1_c8h8.top";
  const std::string mol1_crd_name  = base_chem_name + osc + "lig1_c8h8.inpcrd";
  const std::string mol2_top_name  = base_chem_name + osc + "lig2_c8h8.top";
  const std::string mol2_crd_name  = base_chem_name + osc + "lig2_c8h8.inpcrd";

  // Additional ligands, with and without water
  const std::string mol3_top_name    = base_chem_name + osc + "lig3_1fsj.top";
  const std::string mol3_crd_name    = base_chem_name + osc + "lig3_1fsj.inpcrd";
  const std::string drug_top_name    = base_top_name + osc + "drug_example.top";
  const std::string drug_crd_name    = base_crd_name + osc + "drug_example.inpcrd";
  const std::string drug_vs_top_name = base_top_name + osc + "drug_example_vs.top";
  const std::string drug_vs_crd_name = base_crd_name + osc + "drug_example_vs.inpcrd";

  // Small peptides and proteins
  const std::string ala_top_name     = base_nmlp_name + osc + "ala.top";
  const std::string ala_crd_name     = base_nmlc_name + osc + "ala.inpcrd";
  const std::string gly_top_name     = base_nmlp_name + osc + "gly_gly.top";
  const std::string gly_crd_name     = base_nmlc_name + osc + "gly_gly.inpcrd";
  const std::string phe_top_name     = base_nmlp_name + osc + "phe.top";
  const std::string phe_crd_name     = base_nmlc_name + osc + "phe.inpcrd";
  const std::string trpcage_top_name = base_top_name + osc + "trpcage_in_water.top";
  const std::string trpcage_crd_name = base_crd_name + osc + "trpcage_in_water.inpcrd";
  const std::string ubiquit_top_name = base_top_name + osc + "ubiquitin.top";
  const std::string ubiquit_crd_name = base_crd_name + osc + "ubiquitin.inpcrd";

  // Water boxes
  const std::string tip5p_top_name = base_top_name + osc + "tip5p.top";
  const std::string tip5p_crd_name = base_crd_name + osc + "tip5p.rst";
  
  // Check the existence of all files
  const std::vector<std::string> top_files = { mol1_top_name, mol2_top_name, mol3_top_name,
                                               drug_top_name, drug_vs_top_name, ala_top_name,
                                               gly_top_name, phe_top_name, trpcage_top_name,
                                               ubiquit_top_name, tip5p_top_name };
  const std::vector<std::string> crd_files = { mol1_crd_name, mol2_crd_name, mol3_crd_name,
                                               drug_crd_name, drug_vs_crd_name, ala_crd_name,
                                               gly_crd_name, phe_crd_name, trpcage_crd_name,
                                               ubiquit_crd_name, tip5p_crd_name };
  bool files_exist = true;
  const size_t nsys = top_files.size();
  for (size_t i = 0; i < nsys; i++) {
    files_exist = (files_exist && getDrivePathType(top_files[i]) == DrivePathType::FILE);
    files_exist = (files_exist && getDrivePathType(crd_files[i]) == DrivePathType::FILE);
  }
  const TestPriority do_tests = (files_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (files_exist == false) {
    rtWarn("Some required files were not found.  Check the ${STORMM_SOURCE} environment variable, "
           "currently set to " + oe.getStormmSourcePath() + ", for validity.  Subsequent tests "
           "will be skipped.", "test_chemical_features");
  }
  std::vector<AtomGraph> sys_ag(nsys);
  std::vector<PhaseSpace> sys_ps(nsys);
  if (files_exist) {
    for (size_t i = 0; i < nsys; i++) {
      sys_ag[i].buildFromPrmtop(top_files[i], ExceptionResponse::SILENT, amber_ancient_bioq,
                                1.2, 2.0, 0.01);
      sys_ps[i].buildFromFile(crd_files[i], CoordinateFileKind::UNKNOWN);
    }
  }
  std::vector<int> first_mol_size(nsys);
  std::vector<ChemicalFeatures> sys_chem;
  const MapRotatableGroups mapg_yes = MapRotatableGroups::YES;
  const MapRotatableGroups mapg_no  = MapRotatableGroups::NO;
  if (files_exist) {
    sys_chem.reserve(nsys);
    for (size_t i = 0; i < nsys; i++) {
      const int2 first_mol_lims = sys_ag[i].getMoleculeLimits(0);
      first_mol_size[i] = first_mol_lims.y - first_mol_lims.x;
      sys_chem.emplace_back(&sys_ag[i], CoordinateFrameReader(sys_ps[i]),
                            (first_mol_size[i] < 120) ? mapg_yes : mapg_no, 300.0, &timer);
    }
  }
  else {
    sys_chem.resize(nsys);  
  }
  std::vector<int> ring_counts(nsys);
  std::vector<int> fused_ring_counts(nsys);
  std::vector<int> mutable_ring_counts(nsys);
  std::vector<int> aromatic_group_counts(nsys);  
  std::vector<int> polar_h_counts(nsys);
  std::vector<int> hbond_donor_counts(nsys);
  std::vector<int> hbond_acceptor_counts(nsys);
  std::vector<int> chiral_center_counts(nsys);
  std::vector<int> rotatable_bond_counts(nsys);
  std::vector<int> cis_trans_bond_counts(nsys);
  for (size_t i = 0; i < nsys; i++) {
    ring_counts[i] = sys_chem[i].getRingCount();
    fused_ring_counts[i] = sys_chem[i].getFusedRingCount();
    mutable_ring_counts[i] = sys_chem[i].getRingCount();
    aromatic_group_counts[i] = sys_chem[i].getAromaticGroupCount();
    polar_h_counts[i] = sys_chem[i].getPolarHydrogenCount();
    hbond_donor_counts[i] = sys_chem[i].getHydrogenBondDonorCount();
    hbond_acceptor_counts[i] = sys_chem[i].getHydrogenBondAcceptorCount();
    chiral_center_counts[i] = sys_chem[i].getChiralCenterCount();
    rotatable_bond_counts[i] = sys_chem[i].getRotatableBondCount();
    cis_trans_bond_counts[i] = sys_chem[i].getCisTransBondCount();
  }
  std::vector<int> ring_cnt_ans           = {    1,    1,    4, 1226, 1226,    0,    0,    1,
                                              1567,  486,  216 };
  std::vector<int> fused_ring_cnt_ans     = {    0,    0,    0,    0,    0,    0,    0,    0,
                                                 1,    0,    0 };
  std::vector<int> mutable_ring_cnt_ans   = {    1,    1,    4, 1226, 1226,    0,    0,    1,
                                              1567,  486,  216 };
  std::vector<int> aromatic_group_cnt_ans = {    1,    1,    3,    1,    1,    0,    0,    1,
                                                 2,    3,    0 };
  std::vector<int> polar_h_cnt_ans        = {    0,    0,    6, 2450, 2450,    2,    3,    2,
                                              3155, 1088,  432 };
  std::vector<int> hbond_donor_cnt_ans    = {    0,    0,    4, 1226, 1226,    2,    3,    2,
                                              1587,  580,  216 };
  std::vector<int> hbond_acceptor_cnt_ans = {    0,    1,    7, 1228, 1228,    4,    6,    4,
                                              1614,  681,  216 };
  std::vector<int> chiral_center_cnt_ans  = {    0,    0,    0,    1,    1,    1,    0,    1,
                                                18,   82,    0 };
  std::vector<int> rotatable_bond_cnt_ans = {    1,    3,    9,    8,    8,    4,    7,    6,
                                                 0,    0,    0 };
  std::vector<int> cis_trans_bond_cnt_ans = {    0,    0,    0,    0,    0,    0,    0,    0,
                                                 0,    0,    0 };
  std::vector<int> trp_cage_lchir = sys_chem[8].getChiralCenters(ChiralOrientation::SINISTER);
  std::vector<int> trp_cage_lchir_ans;
  const ChemicalDetailsKit trp_cdk = sys_ag[8].getChemicalDetailsKit();
  for (int i = 0; i < trp_cdk.natom; i++) {
    if ((trp_cdk.atom_names[i] == char4({'C', 'A', ' ', ' '}) &&
         trp_cdk.res_names[sys_ag[8].getResidueIndex(i)] != char4({'G', 'L', 'Y', ' '})) ||
        (trp_cdk.atom_names[i] == char4({'C', 'B', ' ', ' '}) &&
         trp_cdk.res_names[sys_ag[8].getResidueIndex(i)] == char4({'I', 'L', 'E', ' '}))) {
      trp_cage_lchir_ans.push_back(i);
    }
  }
  section(1);
  check(ring_counts, RelationalOperator::EQUAL, ring_cnt_ans, "Overall counts of ring systems do "
        "not meet expectations.", do_tests);
  check(fused_ring_counts, RelationalOperator::EQUAL, fused_ring_cnt_ans, "Counts of fused ring "
        "systems do not meet expectations.", do_tests);
  check(mutable_ring_counts, RelationalOperator::EQUAL, mutable_ring_cnt_ans, "Counts of mutable "
        "ring systems do not meet expectations.", do_tests);
  check(aromatic_group_counts, RelationalOperator::EQUAL, aromatic_group_cnt_ans, "Counts of "
        "aromatic ring systems do not meet expectations.", do_tests);
  check(polar_h_counts, RelationalOperator::EQUAL, polar_h_cnt_ans, "Counts of polar hydrogens "
        "do not meet expectations.", do_tests);
  check(hbond_donor_counts, RelationalOperator::EQUAL, hbond_donor_cnt_ans, "Counts of hydrogen "
        "bond donors do not meet expectations.", do_tests);
  check(hbond_acceptor_counts, RelationalOperator::EQUAL, hbond_acceptor_cnt_ans, "Counts of "
        "hydrogen bond acceptors do not meet expectations.", do_tests);
  check(chiral_center_counts, RelationalOperator::EQUAL, chiral_center_cnt_ans, "Counts of chiral "
        "centers in various systems do not agree.", do_tests);
  check(rotatable_bond_counts, RelationalOperator::EQUAL, rotatable_bond_cnt_ans, "Counts of "
        "rotatable bonds do not meet expectations.", do_tests);
  check(cis_trans_bond_counts, RelationalOperator::EQUAL, cis_trans_bond_cnt_ans, "Counts of "
        "cis-trans bonds do not meet expectations.", do_tests);
  check(trp_cage_lchir_ans, RelationalOperator::EQUAL, trp_cage_lchir, "Chiral center indices for "
        "L-chiral centers (should be amino acid CA atoms, excluding glycine, plus the isoleucine "
        "CB atom) do not meet expectations.", do_tests);
  const std::string fc_name = base_chem_name + osc + "formal_charges.m";
  const std::string bo_name = base_chem_name + osc + "bond_orders.m";
  const std::string ro_name = base_chem_name + osc + "rotating_groups.m";
  const std::string ch_name = base_chem_name + osc + "chiral_atoms.m";
  const bool snps_exist = (getDrivePathType(fc_name) == DrivePathType::FILE &&
                           getDrivePathType(bo_name) == DrivePathType::FILE &&
                           getDrivePathType(ro_name) == DrivePathType::FILE &&
                           getDrivePathType(ch_name) == DrivePathType::FILE);
  const TestPriority do_snps = (snps_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (snps_exist == false) {
    rtWarn("Snapshot files " + fc_name + ", " + bo_name + ", and " + ro_name + " must be "
           "accessible in order to check formal charge and bond order calculations, "
           "respectively.  Check the ${STORMM_SOURCE} environment variable for validity.  "
           "Subsequent tests will be skipped.", "test_chemical_features");
  }
  bool ch_unwritten = true;
  section(2);
  const TestPriority check_lewis_struct = (do_snps == TestPriority::CRITICAL) ?
                                          TestPriority::NON_CRITICAL : TestPriority::ABORT;
  for (size_t i = 0; i < nsys; i++) {
    snapshot(fc_name, polyNumericVector(sys_chem[i].getFormalCharges()), std::string("fc_") +
             std::to_string(i), 1.0e-6, "Formal charges computed for the system described by " +
             sys_ag[i].getFileName() + " do not meet expectations.", oe.takeSnapshot(), 1.0e-8,
             NumberFormat::STANDARD_REAL,
             (i == 0LLU) ? PrintSituation::OVERWRITE : PrintSituation::APPEND, do_snps);
    snapshot(fc_name, polyNumericVector(sys_chem[i].getZeroKelvinFormalCharges()),
             std::string("fczk_") + std::to_string(i), 1.0e-6, "Formal charges computed for a "
             "representative state of the system described by " + sys_ag[i].getFileName() +
             " do not meet expectations.", oe.takeSnapshot(), 1.0e-8, NumberFormat::INTEGER,
             PrintSituation::APPEND, check_lewis_struct);
    check(sum<double>(sys_chem[i].getZeroKelvinFormalCharges()), RelationalOperator::EQUAL,
          sum<double>(sys_chem[i].getFormalCharges()), "The sum of formal charges computed for a "
          "representative Lewis structure and the resonance structure of " +
          getBaseName(sys_ag[i].getFileName()) + " do not agree.", do_tests);
    snapshot(bo_name, polyNumericVector(sys_chem[i].getBondOrders()), std::string("bo_") +
             std::to_string(i), 1.0e-6, "Bond orders computed for the system described by " +
             sys_ag[i].getFileName() + " do not meet expectations.", oe.takeSnapshot(), 1.0e-8,
             NumberFormat::STANDARD_REAL,
             (i == 0LLU) ? PrintSituation::OVERWRITE : PrintSituation::APPEND, do_snps);
    snapshot(bo_name, polyNumericVector(sys_chem[i].getZeroKelvinBondOrders()),
             std::string("bozk_") + std::to_string(i), 1.0e-6, "Bond orders computed for a "
             "representative state of the system described by " + sys_ag[i].getFileName() +
             " do not meet expectations.", oe.takeSnapshot(), 1.0e-8, NumberFormat::INTEGER,
             PrintSituation::APPEND, check_lewis_struct);
    check(sum<double>(sys_chem[i].getZeroKelvinBondOrders()), RelationalOperator::EQUAL,
          sum<double>(sys_chem[i].getBondOrders()), "The sum of all bond orders computed for a "
          "representative Lewis structure and the resonance structure of " +
          getBaseName(sys_ag[i].getFileName()) + " do not agree.", do_tests);
    if (sys_chem[i].getChiralCenterCount() > 0) {
      snapshot(ch_name, polyNumericVector(sys_chem[i].getChiralCenters()), std::string("chcen_") +
               std::to_string(i), 1.0e-6, "Chiral centers detected for the system described by " +
               sys_ag[i].getFileName() + " do not meet expectations.", oe.takeSnapshot(), 1.0e-8,
               NumberFormat::INTEGER,
               (ch_unwritten) ? PrintSituation::OVERWRITE : PrintSituation::APPEND, do_snps);
      ch_unwritten = false;
    }
    check(sum<double>(sys_chem[i].getFormalCharges()), RelationalOperator::EQUAL,
          Approx(sum<double>(sys_ag[i].getPartialCharge<double>())).margin(1.0e-4), "The sum of "
          "formal charges computed for " + sys_ag[i].getFileName() + " does not match the sum of "
          "partial charges given in the topology.", do_tests);
  }
  
  // Check the rotatable bond groups, and the inversion groups generated for smaller structures
  bool ro_unwritten = true;
  for (size_t i = 0; i < nsys; i++) {
    const int nrotor = sys_chem[i].getRotatableBondCount();
    if (nrotor > 0) {
      const std::vector<IsomerPlan> all_rotors = sys_chem[i].getRotatableBondGroups();
      std::vector<int> rotor_ids;
      for (int j = 0; j < nrotor; j++) {
        rotor_ids.push_back(all_rotors[j].getRootAtom());
        rotor_ids.push_back(all_rotors[j].getPivotAtom());
        const int nr_atom = all_rotors[j].getMovingAtomCount();
        for (int k = 0; k < nr_atom; k++) {
          rotor_ids.push_back(all_rotors[j].getMovingAtom(k));
        }
      }
      snapshot(ro_name, polyNumericVector(rotor_ids), std::string("rotators_") + std::to_string(i),
               1.0e-6, "Rotatable atom IDs obtained for the system described by " +
               sys_ag[i].getFileName() + " do not meet expectations.", oe.takeSnapshot(), 1.0e-8,
               NumberFormat::INTEGER,
               (ro_unwritten) ? PrintSituation::OVERWRITE : PrintSituation::APPEND, do_snps);
      ro_unwritten = false;
    }
    const int nchiral = sys_chem[i].getChiralCenterCount();
    if (nchiral > 0) {
      const std::vector<IsomerPlan> all_invertors = sys_chem[i].getChiralInversionGroups();
      std::vector<int> invertor_ids;
      for (int j = 0; j < nchiral; j++) {
        invertor_ids.push_back(all_invertors[j].getRootAtom());
        invertor_ids.push_back(all_invertors[j].getPivotAtom());
        const int nr_atom = all_invertors[j].getMovingAtomCount();
        for (size_t k = 0; k < nr_atom; k++) {
          invertor_ids.push_back(all_invertors[j].getMovingAtom(k));
        }
      }
      snapshot(ro_name, polyNumericVector(invertor_ids), std::string("invertors_") +
               std::to_string(i), 1.0e-6, "Chiral inversion group atom IDs obtained for the "
               "system described by " + sys_ag[i].getFileName() + " do not meet expectations.",
               oe.takeSnapshot(), 1.0e-8, NumberFormat::INTEGER,
               (ro_unwritten) ? PrintSituation::OVERWRITE : PrintSituation::APPEND, do_snps);
      ro_unwritten = false;
    }    
  }

  // Ring detection in some polycyclic systems
  const std::vector<std::string> ring_top_names = {
    base_chem_name + osc + "anthracene_like.top",
    base_chem_name + osc + "morphine_like.top",
    base_chem_name + osc + "pyrene_like.top"
  };
  const std::vector<std::string> ring_crd_names = {
    base_chem_name + osc + "anthracene_like.inpcrd",
    base_chem_name + osc + "morphine_like.inpcrd",
    base_chem_name + osc + "pyrene_like.inpcrd"
  };
  const size_t nring_mols = ring_top_names.size();
  bool rings_exist = true;
  for (size_t i = 0; i < nring_mols; i++) {
    rings_exist = (rings_exist && getDrivePathType(ring_top_names[i]) == DrivePathType::FILE &&
                   getDrivePathType(ring_top_names[i]) == DrivePathType::FILE);
  }
  if (rings_exist == false) {
    rtWarn("Topologies and coordinates for polycyclic test molecules, i.e. " + ring_top_names[0] +
           " and " + ring_crd_names[0] + ", were not found.  Check the ${STORMM_SOURCE} "
           "environment variable for validity.  Subsequent tests will be skipped.",
           "test_chemical_features");
  }
  std::vector<AtomGraph> ring_ag(nring_mols);
  std::vector<PhaseSpace> ring_ps(nring_mols);
  for (size_t i = 0; i < nring_mols; i++) {
    ring_ag[i].buildFromPrmtop(ring_top_names[i], ExceptionResponse::SILENT, amber_ancient_bioq,
                               1.2, 2.0, 0.01);
    ring_ps[i].buildFromFile(ring_crd_names[i], CoordinateFileKind::UNKNOWN);    
  }
  std::vector<ChemicalFeatures> ring_chem;
  if (rings_exist) {
    ring_chem.reserve(nring_mols);
    for (size_t i = 0; i < nring_mols; i++) {
      const int2 first_mol_lims = ring_ag[i].getMoleculeLimits(0);
      first_mol_size[i] = first_mol_lims.y - first_mol_lims.x;
      ring_chem.emplace_back(&ring_ag[i], CoordinateFrameReader(ring_ps[i]), mapg_yes, 300.0,
                             &timer);
    }
  }
  else {
    ring_chem.resize(nring_mols);
  }
  const TestPriority do_rings = (rings_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  std::vector<int> polycyclic_ring_counts(nring_mols);
  std::vector<int> polycyclic_rotator_counts(nring_mols);
  std::vector<int> polycyclic_cistrans_counts(nring_mols);
  for (size_t i = 0; i < nring_mols; i++) {
    polycyclic_ring_counts[i] = ring_chem[i].getRingCount();
    polycyclic_rotator_counts[i] = ring_chem[i].getRotatableBondCount();
    polycyclic_cistrans_counts[i] = ring_chem[i].getCisTransBondCount();
  }
  const std::vector<int> polycyclic_ring_counts_answer = { 4, 5, 8 };
  const std::vector<int> polycyclic_rotator_counts_answer = { 0, 0, 0 };
  const std::vector<int> polycyclic_cistrans_counts_answer = { 0, 0, 0 };
  check(polycyclic_ring_counts, RelationalOperator::EQUAL, polycyclic_ring_counts_answer,
        "The number of rings detected in polycyclic systems does not meet expectations.",
        do_rings);
  check(polycyclic_rotator_counts, RelationalOperator::EQUAL, polycyclic_rotator_counts_answer,
        "The number of rotatable bonds detected in polycyclic systems does not meet "
        "expectations. (There should be no rotatable bonds.)", do_rings);
  check(polycyclic_cistrans_counts, RelationalOperator::EQUAL, polycyclic_cistrans_counts_answer,
        "The number of cis-trans isomerization sites detected in polycyclic systems does not meet "
        "expectations. (There should be no bonds capable of cis-trans isomerization.)", do_rings);

  // Try creating some chemical features objects from just the topologies (no coordinates)
  std::vector<ChemicalFeatures> sys_chem_agonly;
  sys_chem_agonly.reserve(nsys);
  for (int i = 0; i < nsys; i++) {
    sys_chem_agonly.emplace_back(&sys_ag[i], mapg_no, 300.0, &timer);
  }

  // Try creating a chemical features object, then mapping rotatable groups post-facto
  std::vector<ChemicalFeatures> sys_chem_pf;
  sys_chem_pf.reserve(nsys);
  for (int i = 0; i < nsys; i++) {
    sys_chem_pf.emplace_back(&sys_ag[i], CoordinateFrameReader(sys_ps[i]), mapg_no, 300.0,
                             &timer);
  }
  int j = 0;
  std::vector<int> rotatable_bond_counts_pf(nsys, 0), cis_trans_bond_counts_pf(nsys, 0);
  for (int i = 0; i < nsys; i++) {
    if (first_mol_size[i] < 120) {
      sys_chem_pf[i].findRotatableBondGroups(&timer);
    }
    rotatable_bond_counts_pf[i] = sys_chem_pf[i].getRotatableBondCount();
    cis_trans_bond_counts_pf[i] = sys_chem_pf[i].getCisTransBondCount();
  }
  section(3);
  check(rotatable_bond_counts_pf, RelationalOperator::EQUAL, rotatable_bond_cnt_ans, "Counts of "
        "rotatable bonds do not meet expectations.", do_tests);
  check(cis_trans_bond_counts_pf, RelationalOperator::EQUAL, cis_trans_bond_cnt_ans, "Counts of "
        "cis-trans bonds do not meet expectations.", do_tests);

  // Immobilize some atoms in the Phe topology at the side chain, then inspect whether the
  // rest of the molecule becomes mobile as a result.
  AtomMask phe_side_chain(":PHE & @CB,CG,CD1,CD2,CE1,CE2,CZ", &sys_ag[7], sys_chem_pf[7],
                          sys_ps[7]);
  const std::vector<int> phe_side_chain_heavy = phe_side_chain.getMaskedAtomList();
  sys_ag[7].modifyAtomMobility(phe_side_chain_heavy, MobilitySetting::OFF);
  sys_chem_pf[7].findRotatableBondGroups(&timer);
  check(sys_chem_pf[7].getRotatableBondCount(), RelationalOperator::EQUAL,
        rotatable_bond_cnt_ans[7], "The number of rotatable bonds was reduced as a consequence "
        "of immobilizing some atoms, despite the fact that rotation is still possible.", do_tests);
  const std::vector<IsomerPlan> free_phe_rotors = sys_chem[7].getRotatableBondGroups();
  const std::vector<IsomerPlan> cnst_phe_rotors = sys_chem_pf[7].getRotatableBondGroups();
  const int nrotor = free_phe_rotors.size();
  std::vector<int> free_phe_moving_atoms(nrotor);
  std::vector<int> cnst_phe_moving_atoms(nrotor);
  for (int i = 0; i < nrotor; i++) {
    free_phe_moving_atoms[i] = free_phe_rotors[i].getMovingAtomCount();
    cnst_phe_moving_atoms[i] = cnst_phe_rotors[i].getMovingAtomCount();
  }
  const std::vector<int> free_phe_moving_atoms_ans = { 5, 5, 10, 13, 7, 7 };
  const std::vector<int> cnst_phe_moving_atoms_ans = { 5, 5, 20, 17, 7, 7 };
  check(free_phe_moving_atoms, RelationalOperator::EQUAL, free_phe_moving_atoms_ans, "The number "
        "of moving atoms in each Ace-Phe-Nme rotating group does not meet expectations.",
        do_tests);
  check(cnst_phe_moving_atoms, RelationalOperator::EQUAL, cnst_phe_moving_atoms_ans, "The number "
        "of moving atoms in each constrained Ace-Phe-Nme rotating group does not meet "
        "expectations.", do_tests);
  phe_side_chain.addAtoms(":PHE & @CA", CoordinateFrame(sys_ps[7]));
  const std::vector<int> phe_ca_sc_heavy = phe_side_chain.getMaskedAtomList();  
  sys_ag[7].modifyAtomMobility(phe_ca_sc_heavy, MobilitySetting::OFF);
  sys_chem_pf[7].findRotatableBondGroups(&timer);
  check(sys_chem_pf[7].getRotatableBondCount(), RelationalOperator::EQUAL,
        rotatable_bond_cnt_ans[7] - 1, "The number of rotatable bonds was not reduced as a "
        "consequence of immobilizing even the CA atom.", do_tests);

  // Test bond pattern matching
  TestSystemManager bmp_tsm(base_top_name, "top", { "symmetry_L1", "stereo_L1", "trpcage",
                                                    "symmetry_C1", "symmetry_C2", "symmetry_C3",
                                                    "symmetry_C4", "symmetry_C5", "symmetry_C6" },
                            base_crd_name, "inpcrd",
                            { "symmetry_L1", "stereo_L1", "trpcage", "symmetry_C1", "symmetry_C2",
                              "symmetry_C3", "symmetry_C4", "symmetry_C5", "symmetry_C6" });
  const AtomGraph& symm_ag = bmp_tsm.getTopologyReference(0);
  const CoordinateFrame symm_cf = bmp_tsm.exportCoordinateFrame(0);
  const ChemicalFeatures symm_chemfe(symm_ag, symm_cf, MapRotatableGroups::YES, 300.0, &timer);
  const std::vector<double> symm_formal_charges = symm_chemfe.getFormalCharges();
  const std::vector<double> symm_free_electrons = symm_chemfe.getFreeElectrons();
  const std::vector<ChiralOrientation> symm_chiralities = symm_chemfe.getAtomChirality();
  const std::vector<ullint> symm_rings = symm_chemfe.getRingInclusion();
  const std::vector<int> nitro_atoms = { 6, 9, 14, 16, 23, 25 };
  std::vector<int> nitro_equiv_results(5);
  const std::vector<int> random_atoms = { 8, 2, 9, 18, 7, 24, 0, 1, 3, 20 };
  std::vector<int> decoy_equiv_results(5);
  const std::vector<int> decoy_equiv_ans = { 0, 0, 1, 0, 1 };
  for (int i = 0; i < 5; i++) {
    nitro_equiv_results[i] = matchBondingPattern(symm_ag, symm_formal_charges, symm_free_electrons,
                                                 symm_rings, symm_chiralities, nitro_atoms[0],
                                                 nitro_atoms[i + 1]);
    decoy_equiv_results[i] = matchBondingPattern(symm_ag, symm_formal_charges, symm_free_electrons,
                                                 symm_rings, symm_chiralities, random_atoms[2 * i],
                                                 random_atoms[(2 * i) + 1]);
  }
  check(nitro_equiv_results, RelationalOperator::EQUAL, std::vector<int>(5, 1), "Equivalent "
        "bonding patterns between atoms N1 and N3 of a molecule with three-fold and two-fold "
        "compounded symmetry were not detected as they should have been.",
        bmp_tsm.getTestingStatus());
  check(decoy_equiv_results, RelationalOperator::EQUAL, decoy_equiv_ans, "Equivalent bonding "
        "patterns between decoy atoms of molecule with three-fold and two-fold compounded "
        "symmetry were not detected as they should have been.", bmp_tsm.getTestingStatus());
  const int equiv_timings = timer.addCategory("Symmetry equivalence");

  // Test atom equivalence determination
  const AtomEquivalence symm_eq(symm_ag, symm_formal_charges, symm_free_electrons, symm_rings,
                                symm_chemfe.getAtomChirality(), &timer);
  check(symm_eq.getGroupCount(), RelationalOperator::EQUAL, 4, "The number of symmetric atom "
        "groups with interchangeable atoms in a highly symmetric molecule was not determined as "
        "expected.", bmp_tsm.getTestingStatus());
  check(static_cast<int>(symm_eq.getAsymmetricAtoms().size()), RelationalOperator::EQUAL, 1,
        "The asymmetric component of the small molecule should include only one atom, N2 at the "
        "center.", bmp_tsm.getTestingStatus());
  const AtomGraph &trpi_ag = bmp_tsm.getTopologyReference(2);
  const CoordinateFrame trpi_cf = bmp_tsm.exportCoordinateFrame(2);
  const ChemicalFeatures trpi_chemfe (trpi_ag, trpi_cf, MapRotatableGroups::NO, 300.0, &timer);
  const std::vector<double> trpi_formal_charges = trpi_chemfe.getFormalCharges();
  const std::vector<double> trpi_free_electrons = trpi_chemfe.getFreeElectrons();
  const std::vector<ChiralOrientation> trpi_chiralities = trpi_chemfe.getAtomChirality();
  const std::vector<ullint> trpi_rings = trpi_chemfe.getRingInclusion();
  const AtomEquivalence trpi_eq(trpi_ag, trpi_formal_charges, trpi_free_electrons, trpi_rings,
                                trpi_chiralities, &timer);
  const int n_trpi_groups = trpi_eq.getGroupCount();
  check(n_trpi_groups, RelationalOperator::EQUAL, 52, "The number of symmetric atom groups found "
        "for Trp-cage miniprotein does not meet expectations.", bmp_tsm.getTestingStatus());
  std::vector<int> trpi_group_sizes(n_trpi_groups), trpi_group_orders(n_trpi_groups);
  bool trpi_arg_group_found = false;
  const ChemicalDetailsKit trpi_cdk = trpi_ag.getChemicalDetailsKit();
  int trpi_symm_atom_count = 0;
  for (int i = 0; i < n_trpi_groups; i++) {
    trpi_group_sizes[i] = trpi_eq.getGroupSize(i);
    trpi_group_orders[i] = trpi_eq.getGroupOrder(i);
    trpi_symm_atom_count += trpi_group_sizes[i] * trpi_group_orders[i];
    if (trpi_group_sizes[i] == 3 && trpi_group_orders[i] == 2) {
      const std::vector<int> trpi_arg_group = trpi_eq.getGroup(i);
      const int narg_group = trpi_arg_group.size();
      const char4 arg_res = { 'A', 'R', 'G', ' ' };
      const char4 nh1_atm = { 'N', 'H', '1', ' ' };
      if (trpi_cdk.res_names[trpi_ag.getResidueIndex(trpi_arg_group[0])] == arg_res) {
        for (int j = 0; j < narg_group; j++) {
          const int j_atom = trpi_arg_group[j];
          if (fabs(trpi_cdk.masses[j_atom] - 14.01) < 1.0e-4 &&
              trpi_cdk.atom_names[j_atom] == nh1_atm) {
            trpi_arg_group_found = true;
          }
        }
      }
    }
  }
  check(maxValue(trpi_group_sizes), RelationalOperator::EQUAL, 4, "The maximum size of a group "
        "in the Trp-cage miniprotein does not meet expectations.", bmp_tsm.getTestingStatus());
  check(maxValue(trpi_group_orders), RelationalOperator::EQUAL, 3, "The maximum order of a group "
        "in the Trp-cage miniprotein does not meet expectations.", bmp_tsm.getTestingStatus());
  check(trpi_arg_group_found, "The ARG guanidino head group did not appear in the list of "
        "Trp-cage symmetry groups.", bmp_tsm.getTestingStatus());
  const std::vector<int> trpi_asymm_atoms = trpi_eq.getAsymmetricAtoms();
  check(static_cast<int>(trpi_asymm_atoms.size()), RelationalOperator::EQUAL, 186, "The number of "
        "asymmetric atoms in Trp-cage miniprotein does not account for atoms not in symmetry "
        "groups.  This may indicate that nested symmetry groups (dependencies) were found in "
        "Trp-cage miniprotein.", bmp_tsm.getTestingStatus());
  check(trpi_symm_atom_count, RelationalOperator::EQUAL, 134, "The number of atoms in symmetry "
        "groups in Trp-cage miniprotein does not meet expectations.  Various nested groups should "
        "be present, raising the total.", bmp_tsm.getTestingStatus());
  bool no_asymm_overlap = true;
  const int trpi_asymm_atom_count = trpi_asymm_atoms.size();
  bool asymm_distinct = true;
  for (int i = 0; i < n_trpi_groups; i++) {
    const int jlim = trpi_eq.getGroupSize(i) * trpi_eq.getGroupOrder(i);
    const std::vector<int> jgroup = trpi_eq.getGroup(i);
    for (int j = 0; j < jlim; j++) {
      bool atom_found = false;
      for (int k = 0; k < trpi_asymm_atom_count; k++) {
        atom_found = (atom_found || (jgroup[j] == trpi_asymm_atoms[k]));
      }
      asymm_distinct = (asymm_distinct && atom_found == false);
    }
  }
  check(asymm_distinct, "Overlap was found between atoms in the asymmetric regions of Trp-cage "
        "miniprotein and one or more of its symmetry-equivalent groups.",
        bmp_tsm.getTestingStatus());
  for (int i = 0; i < 9; i++) {
    const CoordinateFrame ifrm(bmp_tsm.exportCoordinateFrame(i));
    const AtomGraph& iag_ref = bmp_tsm.getTopologyReference(i);
    const ChemicalFeatures ichemfe(iag_ref, ifrm);
    const AtomRank irnks(iag_ref, ichemfe);
    const std::vector<int>& atom_ranks = irnks.getRanks();
    const int natom = iag_ref.getAtomCount();
    std::vector<bool> bonds_same(natom * natom);
    for (int j = 0; j < natom; j++) {
      for (int k = j; k < natom; k++) {
        bonds_same[(j * natom) + k] = matchBondingPattern(iag_ref, ichemfe, j, k);
        bonds_same[(k * natom) + j] = bonds_same[(j * natom) + k];
      }
    }
    int false_positives = 0;
    int false_negatives = 0;
    for (int j = 0; j < natom; j++) {
      for (int k = j + 1; k < natom; k++) {
        if (atom_ranks[k] == atom_ranks[j]) {
          if (bonds_same[(j * natom) + k] == false) {
            false_positives++;
          }
        }
        else {
          if (bonds_same[(j * natom) + k]) {
            false_negatives++;
          }
        }
      }
    }
    check(false_positives, RelationalOperator::EQUAL, 0, "Atom ranks incorrectly identify atoms "
          "as equivalent when a direct comparison of their bonding patterns indicates they are "
          "not.", bmp_tsm.getTestingStatus());
    check(false_negatives, RelationalOperator::EQUAL, 0, "Atom ranks fail to identify atoms as "
          "equivalent when a direct comparison of their bonding patterns indicates they are.",
          bmp_tsm.getTestingStatus());
  }
  
  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
