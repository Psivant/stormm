#include "copyright.h"
#include "../../src/Chemistry/atommask.h"
#include "../../src/Chemistry/chemical_features.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::card::HybridTargetLevel;
using stormm::data_types::uint;
#ifndef STORMM_USE_HPC
using stormm::data_types::char4;
#endif
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::osSeparator;
using stormm::errors::rtWarn;
using stormm::parse::NumberFormat;
using stormm::parse::polyNumericVector;
using stormm::data_types::operator==;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using stormm::topology::AtomGraph;
using stormm::topology::ChemicalDetailsKit;
using stormm::trajectory::CoordinateFileKind;
using stormm::trajectory::CoordinateFrame;
using stormm::trajectory::PhaseSpace;

using namespace stormm::chemistry;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Section 1
  section("Masks based on atom or residue numbers");
    
  // Section 2
  section("Masks based on atom or residue names");

  // Section 3
  section("Masks with wildcards in names");

  // Section 4
  section("Basic AND, OR, and NOT operators");
  
  // Section 5
  section("Masks with multiple scopes");

  // Section 6
  section("Masks with ranged operators");

  // Section 7
  section("Mask accessibility features");

  // Test the most basic atom masks
  section(1);
  AtomGraph trpcage;
  CoordinateFrame trpcage_crd;
  const char osc = osSeparator();
  const std::string trpcage_top_name = oe.getStormmSourcePath() + osc + "test" +  osc +
                                       "Topology" + osc + "trpcage_in_water.top";
  const std::string trpcage_crd_name = oe.getStormmSourcePath() + osc + "test" +  osc +
                                       "Trajectory" + osc + "trpcage_in_water.inpcrd";
  const bool trpcage_exists = (getDrivePathType(trpcage_top_name) == DrivePathType::FILE &&
                               getDrivePathType(trpcage_crd_name) == DrivePathType::FILE);
  const TestPriority do_trpcage = (trpcage_exists) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (trpcage_exists) {
    trpcage.buildFromPrmtop(trpcage_top_name);
    trpcage_crd.buildFromFile(trpcage_crd_name, CoordinateFileKind::AMBER_INPCRD);
  }
  else {
    rtWarn("Atom masks refer to a specific topology, but either or both of the files " +
           trpcage_top_name + " and " + trpcage_crd_name + " were not found.  Check the "
           "$STORMM_SOURCE environment variable to ensure that it is set to the STORMM root "
           "directory where src/ and test/ subdirectories can be found.  Some subsequent tests "
           "will be skipped.", "test_atommask");
  }
  const ChemicalFeatures trpcage_chem = (trpcage_exists) ?
                                        ChemicalFeatures(&trpcage, trpcage_crd) :
                                        ChemicalFeatures();
  const std::string mask_a_input(":1-5");
  const AtomMask mask_a = (trpcage_exists) ? AtomMask(mask_a_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "First five residues") : AtomMask();
  std::vector<int> mask_a_answer(trpcage.getAtomCount(), 0);
  for (int i = 0; i < trpcage.getResidueLimits(4).y; i++) {
    mask_a_answer[i] = 1;
  }
  const std::vector<uint> mavec = mask_a.getRawMask();
  const std::string trpcage_mask_snp = oe.getStormmSourcePath() + osc + "test" + osc +
                                       "Chemistry" + osc + "trpcage_mask_image.m";
  const bool trp_snap_exist = (getDrivePathType(trpcage_mask_snp) == DrivePathType::FILE);
  if (trp_snap_exist == false && oe.takeSnapshot() != SnapshotOperation::SNAPSHOT) {
    rtWarn("The snapshot file " + trpcage_mask_snp + ", needed for snapshot tests on some mask "
           "results, could not be found.  If this occurred without also triggering a warning on "
           "files for the Trp-cage system itself, the installation may be corrupted.  Check the "
           "${STORMM_SOURCE} environment variable to ensure that it is set to the STORMM root "
           "directory where src/ and test/ subdirectories can be found.  Some tests will be "
           "skipped.", "test_atommask");
  }
  const TestPriority do_trpsnap = (trp_snap_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  const std::vector<bool> amask_result = mask_a.getMask();
  check(std::vector<int>(amask_result.begin(), amask_result.end()), RelationalOperator::EQUAL,
        mask_a_answer, "Atom mask \"" + mask_a_input + "\" does not generate the expected "
        "pattern.", do_trpcage);
  section(7);
  snapshot(trpcage_mask_snp, polyNumericVector(mavec), "rawmask_a", NumberFormat::UNSIGNED_INTEGER,
           "Raw mask data for \"" + mask_a_input + "\" acting on AtomGraph built from " +
           trpcage_top_name + " does not meet expectations.", oe.takeSnapshot(), 1.0e-1, 1.0e-8,
           PrintSituation::OVERWRITE, do_trpcage);
  const std::string mask_b_input(":1-5, 15-19, 237");
  const AtomMask mask_b = (trpcage_exists) ? AtomMask(mask_b_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "Scattered residues") : AtomMask();
  const std::vector<uint> mbvec = mask_b.getRawMask();
  std::vector<int> mask_b_answer = mask_a_answer;
  for (int i = trpcage.getResidueLimits(14).x; i < trpcage.getResidueLimits(18).y; i++) {
    mask_b_answer[i] = 1;
  }
  for (int i = trpcage.getResidueLimits(236).x; i < trpcage.getResidueLimits(236).y; i++) {
    mask_b_answer[i] = 1;
  }
  const std::vector<bool> bmask_result = mask_b.getMask();
  section(1);
  check(std::vector<int>(bmask_result.begin(), bmask_result.end()), RelationalOperator::EQUAL,
        mask_b_answer, "Atom mask \"" + mask_b_input + "\" does not generate the expected "
        "pattern.", do_trpcage);
  snapshot(trpcage_mask_snp, polyNumericVector(mbvec), "rawmask_b", NumberFormat::UNSIGNED_INTEGER,
           "Raw mask data for \"" + mask_b_input + "\" acting on AtomGraph built from " +
           trpcage_top_name + " does not meet expectations.", oe.takeSnapshot(), 1.0e-1, 1.0e-8,
           PrintSituation::APPEND, do_trpcage);

  // Switch to a ChemicalDetailsKit for simpler access to the required information.  There are
  // many more checks to perform.  Try atom selections.
  ChemicalDetailsKit trp_cdk = trpcage.getChemicalDetailsKit(HybridTargetLevel::HOST);
  const std::string mask_c_input("@2-5,151-155,298,1014-1207,45");
  const AtomMask mask_c = (trpcage_exists) ? AtomMask(mask_c_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "Scattered atoms") : AtomMask();
  std::vector<int> mask_c_answer(trp_cdk.natom, 0);
  for (int i = 0; i < trp_cdk.natom; i++) {
    const int ip1 = i + 1;
    mask_c_answer[i] = ((ip1 >= 2 && ip1 <= 5) || (ip1 >= 151 && ip1 <= 155) || ip1 == 298 ||
                        (ip1 >= 1014 && ip1 <= 1207) || ip1 == 45);
  }
  const std::vector<bool> cmask_result = mask_c.getMask();
  check(std::vector<int>(cmask_result.begin(), cmask_result.end()), RelationalOperator::EQUAL,
        mask_c_answer, "Atom mask \"" + mask_c_input + "\" does not generate the expected "
        "pattern.", do_trpcage);

  // Try using the OR operator to combine atom selections, rather than just commas
  section(4);
  const std::string mask_d_input("@2-5,45,298 | @151-155 | @1014-1207");
  const AtomMask mask_d = (trpcage_exists) ? AtomMask(mask_d_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "Scattered atoms linked by OR operators") :
                                             AtomMask();
  const std::vector<bool> dmask_result = mask_d.getMask();
  check(std::vector<int>(dmask_result.begin(), dmask_result.end()), RelationalOperator::EQUAL,
        mask_c_answer, "Atom mask " + mask_d_input + " does not generate the expected pattern.",
        do_trpcage);
  check(mask_d.getRawMask(), RelationalOperator::EQUAL, mask_c.getRawMask(), "Atom masks \"" +
        mask_c_input + "\" and \"" + mask_d_input + "\" do not generate identical raw bitstreams "
        "as they should.", do_trpcage);

  // Test the AND operator
  const std::string mask_e_input("@2-5,45,298 & :1-18");
  const AtomMask mask_e = (trpcage_exists) ? AtomMask(mask_e_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "Scattered atoms limited by an "
                                                      "AND operator") : AtomMask();
  const std::vector<bool> emask_result = mask_e.getMask();
  std::vector<int> mask_e_answer(trp_cdk.natom, 0);
  for (int i = 0; i < trp_cdk.natom; i++) {
    const int ip1 = i + 1;
    mask_e_answer[i] = ((ip1 >= 2 && ip1 <= 5) || ip1 == 45);
  }
  check(std::vector<int>(emask_result.begin(), emask_result.end()), RelationalOperator::EQUAL,
        mask_e_answer, "Atom mask \"" + mask_e_input + "\" does not generate the expected "
        "pattern.", do_trpcage);

  // Try a combination of OR and AND operators
  const std::string mask_f_input("@2-5,45,298 | @151-155 | @1014-1207 & :1-18");
  const AtomMask mask_f = (trpcage_exists) ? AtomMask(mask_f_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "Scattered atoms linked by OR operators and "
                                                      "limited by an AND operator") : AtomMask();
  const std::vector<bool> fmask_result = mask_f.getMask();
  std::vector<int> mask_f_answer(trp_cdk.natom, 0);
  for (int i = 0; i < trp_cdk.natom; i++) {
    const int ip1 = i + 1;
    mask_f_answer[i] = (ip1 >= 2 && ip1 <= 5) || ip1 == 45 || ip1 == 298 ||
                       (ip1 >= 151 && ip1 <= 155);
  }
  check(std::vector<int>(fmask_result.begin(), fmask_result.end()), RelationalOperator::EQUAL,
        mask_f_answer, "Atom mask \"" + mask_f_input + "\" does not generate the expected "
        "pattern.", do_trpcage);

  // Throw in a NOT operator to modify one of the atom selections
  const std::string mask_g_input("! @2-5,45,298 | @151-155 | @1014-1207 & :1-18");
  const AtomMask mask_g = (trpcage_exists) ? AtomMask(mask_g_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "Scattered atoms linked by OR operators, "
                                                      "limited by an AND operator, partially "
                                                      "inverted by a NOT operator.") : AtomMask();
  const std::vector<bool> gmask_result = mask_g.getMask();
  std::vector<int> mask_g_answer(trp_cdk.natom, 0);
  for (int i = 0; i < trp_cdk.natom; i++) {
    const int ip1 = i + 1;
    mask_g_answer[i] = !((ip1 >= 2 && ip1 <= 5) || ip1 == 45 || ip1 == 298);
  }
  check(std::vector<int>(gmask_result.begin(), gmask_result.end()), RelationalOperator::EQUAL,
        mask_g_answer, "Atom mask \"" + mask_g_input + "\" does not generate the expected "
        "pattern.", do_trpcage);

  // Try making the spacing dumb--no one would write a mask like this or purpose but it's just
  // ugly, not broken.
  const std::string mask_h_input("!@2 -5, 45,298|@151- 155| @1014-1207&:1- 18  ");
  const AtomMask mask_h = (trpcage_exists) ? AtomMask(mask_h_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "Scattered atoms linked by OR operators, "
                                                      "limited by an AND operator, inverted by a "
                                                      "NOT operator, with stupid spacing") :
                                             AtomMask();
  const std::vector<bool> hmask_result = mask_h.getMask();
  check(std::vector<int>(hmask_result.begin(), hmask_result.end()), RelationalOperator::EQUAL,
        mask_g_answer, "Atom mask \"" + mask_h_input + "\" does not generate the expected "
        "pattern.", do_trpcage);

  // Try selecting atoms and residues by name, not number
  section(2);
  const std::string mask_i_input("@CA | @HA,CD,O & !:WAT"); 
  const AtomMask mask_i = (trpcage_exists) ? AtomMask(mask_i_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "A mask based on atom names, excluding "
                                                      "water residues") : AtomMask();
  std::vector<int> mask_i_answer(trp_cdk.natom, 0);
  const char4 ca_c4 = {'C', 'A', ' ', ' '};
  const char4 ha_c4 = {'H', 'A', ' ', ' '};
  const char4 cd_c4 = {'C', 'D', ' ', ' '};
  const char4  o_c4 = {'O', ' ', ' ', ' '};
  for (int i = 0; i < trp_cdk.natom; i++) {
    if (trp_cdk.atom_names[i] == ca_c4 || trp_cdk.atom_names[i] == ha_c4 ||
        trp_cdk.atom_names[i] == cd_c4 || trp_cdk.atom_names[i] ==  o_c4) {
      mask_i_answer[i] = 1;
    }
  }
  const char4 wat_c4 = {'W', 'A', 'T', ' '};
  for (int i = 0; i < trp_cdk.nres; i++) {
    if (trp_cdk.res_names[i] == wat_c4) {
      for (int j = trp_cdk.res_limits[i]; j < trp_cdk.res_limits[i + 1]; j++) {
        mask_i_answer[j] = 0;
      }
    }
  }
  const std::vector<bool> imask_result = mask_i.getMask();
  check(std::vector<int>(imask_result.begin(), imask_result.end()), RelationalOperator::EQUAL,
        mask_i_answer, "Atom mask \"" + mask_i_input + "\" does not generate the expected "
        "pattern.", do_trpcage);
  const std::string mask_j_input(":TRP | :ARG,ASN,PRO & !:WAT"); 
  const AtomMask mask_j = (trpcage_exists) ? AtomMask(mask_j_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "A mask based on residue names") :
                                             AtomMask();
  std::vector<int> mask_j_answer(trp_cdk.natom, 0);
  const char4 ala_c4 = {'A', 'L', 'A', ' '};
  const char4 arg_c4 = {'A', 'R', 'G', ' '};
  const char4 asn_c4 = {'A', 'S', 'N', ' '};
  const char4 pro_c4 = {'P', 'R', 'O', ' '};
  const char4 thr_c4 = {'T', 'H', 'R', ' '};
  const char4 trp_c4 = {'T', 'R', 'P', ' '};
  const char4 tyr_c4 = {'T', 'Y', 'R', ' '};
  for (int i = 0; i < trp_cdk.nres; i++) {
    if (trp_cdk.res_names[i] == trp_c4 || trp_cdk.res_names[i] == arg_c4 ||
        trp_cdk.res_names[i] == asn_c4 || trp_cdk.res_names[i] == pro_c4) {
      for (int j = trp_cdk.res_limits[i]; j < trp_cdk.res_limits[i + 1]; j++) {
        mask_j_answer[j] = 1;
      }
    }
  }
  const std::vector<bool> jmask_result = mask_j.getMask();
  check(std::vector<int>(jmask_result.begin(), jmask_result.end()), RelationalOperator::EQUAL,
        mask_j_answer, "Atom mask \"" + mask_j_input + "\" does not generate the expected "
        "pattern.", do_trpcage);

  // Try selecting atoms and residues by names with wildcards
  section(3);
  const std::string mask_k_input(":T?? | :A??,PR? & !:TYR,ASP");  
  const AtomMask mask_k = (trpcage_exists) ? AtomMask(mask_k_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "A mask based on residue names, with "
                                                      "wildcard characters") : AtomMask();
  std::vector<int> mask_k_answer(trp_cdk.natom, 0);
  for (int i = 0; i < trp_cdk.nres; i++) {
    if (trp_cdk.res_names[i] == ala_c4 || trp_cdk.res_names[i] == arg_c4 ||
        trp_cdk.res_names[i] == asn_c4 || trp_cdk.res_names[i] == pro_c4 ||
        trp_cdk.res_names[i] == thr_c4 || trp_cdk.res_names[i] == trp_c4 ||
        trp_cdk.res_names[i] == tyr_c4) {
      for (int j = trp_cdk.res_limits[i]; j < trp_cdk.res_limits[i + 1]; j++) {
        mask_k_answer[j] = 1;
      }
    }
  }
  const std::vector<bool> kmask_result = mask_k.getMask();
  check(std::vector<int>(kmask_result.begin(), kmask_result.end()), RelationalOperator::EQUAL,
        mask_k_answer, "Atom mask \"" + mask_k_input + "\" does not generate the expected "
        "pattern.", do_trpcage);
  const std::string mask_l_input(":T* | :A=,P* & !:?Y=,A*P,=T");
  const AtomMask mask_l = (trpcage_exists) ? AtomMask(mask_l_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "A mask based on residue names, with "
                                                      "wildcard stretches") : AtomMask();
  const std::vector<bool> lmask_result = mask_l.getMask();
  check(std::vector<int>(lmask_result.begin(), lmask_result.end()), RelationalOperator::EQUAL,
        mask_k_answer, "Atom mask \"" + mask_l_input + "\" does not generate the expected "
        "pattern.", do_trpcage);
  check(mask_l.getMaskedAtomCount(), RelationalOperator::EQUAL, 141, "Atom mask \"" +
        mask_l_input + "\" does not generate the expected number of masked atoms.\n");
  const std::string mask_m_input("@?A | @*D,O & !:W*??");
  const AtomMask mask_m = (trpcage_exists) ? AtomMask(mask_m_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "A mask based on atom names, with wildcards "
                                                      "of all kinds") : AtomMask();
  const std::vector<bool> mmask_result = mask_m.getMask();
  check(std::vector<int>(mmask_result.begin(), mmask_result.end()), RelationalOperator::EQUAL,
        mask_i_answer, "Atom mask " + mask_m_input + " does not generate the expected pattern.",
        do_trpcage);
  
  // Apply parentheses to get multiple scopes
  section(5);
  const std::string mask_n_input(":7-8 | (:6-12 & !:7-8)"); 
  const AtomMask mask_n = (trpcage_exists) ? AtomMask(mask_n_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "A continuous stretch of residues, where "
                                                      "two residues are eliminated in the context "
                                                      "of a protected scope and added back by the "
                                                      "preceding OR operator") : AtomMask();
  std::vector<int> mask_n_answer(trp_cdk.natom, 0);
  for (int i = 0; i < trp_cdk.natom; i++) {
    const int atom_rn = trp_cdk.res_numbers[i];
    mask_n_answer[i] = (atom_rn >= 6 && atom_rn <= 12);
  }
  const std::vector<bool> nmask_result = mask_n.getMask();
  check(std::vector<int>(nmask_result.begin(), nmask_result.end()), RelationalOperator::EQUAL,
        mask_n_answer, "Atom mask \"" + mask_n_input + "\" does not generate the expected "
        "pattern.", do_trpcage);
  const std::string mask_o_input("[{(:7-8 | ([:6-12] & !{:7-8}))}]"); 
  const AtomMask mask_o = (trpcage_exists) ? AtomMask(mask_o_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "A continuous stretch of residues, where "
                                                      "two residues are eliminated in the context "
                                                      "of a protected scope and added back by the "
                                                      "preceding OR operator") : AtomMask();
  const std::vector<bool> omask_result = mask_o.getMask();
  check(std::vector<int>(omask_result.begin(), omask_result.end()), RelationalOperator::EQUAL,
        mask_n_answer, "Atom mask \"" + mask_o_input + "\" does not generate the expected "
        "pattern.", do_trpcage);

  // Try selections using elements (@/)
  section(1);
  const std::string mask_p_input("@/6,8"); 
  const AtomMask mask_p = (trpcage_exists) ? AtomMask(mask_p_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "Carbon and oxygen, selected by atomic "
                                                      "number") : AtomMask();
  std::vector<int> mask_p_answer(trp_cdk.natom, 0);
  for (int i = 0; i < trp_cdk.natom; i++) {
    const int znum = trp_cdk.z_numbers[i];
    mask_p_answer[i] = (znum == 6 || znum == 8);
  }
  const std::vector<bool> pmask_result = mask_p.getMask();
  check(std::vector<int>(pmask_result.begin(), pmask_result.end()), RelationalOperator::EQUAL,
        mask_p_answer, "Atom mask \"" + mask_p_input + "\" does not generate the expected "
        "pattern.", do_trpcage);
  const std::string mask_q_input("@/6-8"); 
  const AtomMask mask_q = (trpcage_exists) ? AtomMask(mask_q_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "Carbon, nitrogen, and oxygen, selected by "
                                                      "atomic number") :
                                             AtomMask();
  std::vector<int> mask_q_answer(trp_cdk.natom, 0);
  for (int i = 0; i < trp_cdk.natom; i++) {
    const int znum = trp_cdk.z_numbers[i];
    mask_q_answer[i] = (znum >= 6 && znum <= 8);
  }
  const std::vector<bool> qmask_result = mask_q.getMask();
  check(std::vector<int>(qmask_result.begin(), qmask_result.end()), RelationalOperator::EQUAL,
        mask_q_answer, "Atom mask \"" + mask_q_input + "\" does not generate the expected "
        "pattern.", do_trpcage);
  section(2);
  const std::string mask_r_input("@/C?,O="); 
  const AtomMask mask_r = (trpcage_exists) ? AtomMask(mask_r_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "Carbon and oxygen, selected by atomic "
                                                      "number (wildcards were added that could "
                                                      "increase the number of detected atoms, but "
                                                      "not for this simple peptide in water "
                                                      "system)") : AtomMask();
  const std::vector<bool> rmask_result = mask_r.getMask();
  check(std::vector<int>(rmask_result.begin(), rmask_result.end()), RelationalOperator::EQUAL,
        mask_p_answer, "Atom mask \"" + mask_r_input + "\" does not generate the expected "
        "pattern.", do_trpcage);
  const std::string mask_s_input("@/ C, N, O"); 
  const AtomMask mask_s = (trpcage_exists) ? AtomMask(mask_s_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "Carbon, nitrogen, and oxygen, selected by "
                                                      "atomic number (wildcards were added that "
                                                      "could increase the number of detected "
                                                      "atoms, but not for this simple peptide in "
                                                      "water system)") :
                                             AtomMask();
  const std::vector<bool> smask_result = mask_s.getMask();
  check(std::vector<int>(smask_result.begin(), smask_result.end()), RelationalOperator::EQUAL,
        mask_q_answer, "Atom mask \"" + mask_s_input + "\" does not generate the expected "
        "pattern.", do_trpcage);

  // Select atoms by atom type (name)
  const std::string mask_t_input("@%C, N, OW"); 
  AtomMask mask_t = (trpcage_exists) ? AtomMask(mask_t_input, &trpcage, trpcage_chem, trpcage_crd,
                                                MaskInputMode::AMBMASK, "Backbone sp2 carbon, "
                                                "nitrogen, and TIP3P water oxygen, selected by "
                                                "atomic number (wildcards were added that could "
                                                "increase the number of detected atoms, but not "
                                                "for this simple peptide in water system)") :
                                       AtomMask();
  std::vector<int> mask_t_answer(trp_cdk.natom);
  const char4 atc_c4 = {'C', ' ', ' ', ' '};
  const char4 atn_c4 = {'N', ' ', ' ', ' '};
  const char4 atow_c4 = {'O', 'W', ' ', ' '};
  for (int i = 0; i < trp_cdk.natom; i++) {
    const char4 tmp_at = trp_cdk.atom_types[i];
    mask_t_answer[i] = (tmp_at == atc_c4 || tmp_at == atn_c4 || tmp_at == atow_c4);
  }
  const std::vector<bool> tmask_result = mask_t.getMask();
  check(std::vector<int>(tmask_result.begin(), tmask_result.end()), RelationalOperator::EQUAL,
        mask_t_answer, "Atom mask \"" + mask_t_input + "\" does not generate the expected "
        "pattern.", do_trpcage);

  // Additional checks on mask getter functions
  section(7);
  check(mask_a.getMaskedAtomCount(), RelationalOperator::EQUAL, 92, "Atom mask \"" + mask_a_input +
        "\" does not report the correct number of masked atoms.", do_trpcage);
  const std::vector<int> mask_e_list_answer = { 1, 2, 3, 4, 44 };
  check(mask_e.getMaskedAtomList(), RelationalOperator::EQUAL, mask_e_list_answer, "Atom mask \"" +
        mask_e_input + "\" does not report the correct list of masked atoms.", do_trpcage);
  std::string mask_x_input("");
  const AtomMask mask_x = (trpcage_exists) ? AtomMask(mask_x_input, &trpcage, trpcage_chem,
                                                      trpcage_crd, MaskInputMode::AMBMASK,
                                                      "No selection") :
                                             AtomMask();
  check(mask_x.getMaskedAtomCount(), RelationalOperator::EQUAL, 0, "A blank atom mask creates a "
        "significant group of atoms.", do_trpcage);
  check(strncmp(mask_t.getDescription().c_str(), "Backbone sp2 c", 14) == 0, "The description "
        "returned by a mask does not meet expectations.", do_trpcage);
  check(mask_t.getInputText(), RelationalOperator::EQUAL, std::string("@%C, N, OW"), "A mask "
        "input string was not properly transcribed.", do_trpcage);

  // Check mask post-processing operations: adding atoms
  const std::vector<int> atom_adds = { 0, 1, 2, 3, 4, 5, 6, 7, 16, 17, 18, 19, 24, 25, 32 };
  mask_t.addAtoms(atom_adds);
  check(mask_t.getMaskedAtomCount(), RelationalOperator::EQUAL, 1611, "The number of atoms in a "
        "mask after adding a new list does not meet expectations.", do_trpcage);
  check(mask_t.getMaskedAtomList().size(), RelationalOperator::EQUAL, 1611, "The list of masked "
        "atoms is no longer consistent with the reported size after new additions.", do_trpcage);
  CHECK_THROWS_SOFT(mask_t.addAtoms(std::vector<int>(2, 5700000), ExceptionResponse::DIE),
                    "Atoms with bogus indices were incorrectly added to a mask.", do_trpcage);
  
  // Try ranged operations
  const std::string mask_u_input(":1-20 @< 5.0");
  
  // Probe the chemical features
  AtomGraph trpi_ag, drug_ag;
  PhaseSpace trpi_ps, drug_ps;
  const std::string trpi_top_name = oe.getStormmSourcePath() + osc + "test" +  osc + "Topology" +
                                    osc + "trpcage.top";
  const std::string trpi_crd_name = oe.getStormmSourcePath() + osc + "test" +  osc +
                                    "Trajectory" + osc + "trpcage.inpcrd";
  const bool trpi_exists = (getDrivePathType(trpi_top_name) == DrivePathType::FILE &&
                            getDrivePathType(trpi_crd_name) == DrivePathType::FILE);
  const TestPriority do_trpi = (trpi_exists) ? TestPriority::CRITICAL : TestPriority::ABORT;
  const std::string drug_top_name = oe.getStormmSourcePath() + osc + "test" +  osc + "Topology" +
                                    osc + "drug_example.top";
  const std::string drug_crd_name = oe.getStormmSourcePath() + osc + "test" +  osc +
                                    "Trajectory" + osc + "drug_example.inpcrd";
  const bool drug_exists = (getDrivePathType(drug_top_name) == DrivePathType::FILE &&
                            getDrivePathType(drug_crd_name) == DrivePathType::FILE);
  const TestPriority do_drug = (drug_exists) ? TestPriority::CRITICAL : TestPriority::ABORT;
  bool missing_files = false;
  if (trpi_exists) {
    trpi_ag.buildFromPrmtop(trpi_top_name);
    trpi_ps.buildFromFile(trpi_crd_name, CoordinateFileKind::AMBER_INPCRD);
  }
  else {
    missing_files = true;
  }
  if (drug_exists) {
    drug_ag.buildFromPrmtop(drug_top_name);
    drug_ps.buildFromFile(drug_crd_name, CoordinateFileKind::AMBER_INPCRD);
  }
  else {
    missing_files = true;
  }
  if (missing_files) {
    rtWarn("Additional topology and coordinate files needed by other tests were not found.  Check "
           "the $STORMM_SOURCE environment variable.  Some subsequent tests will be skipped.",
           "test_atommask");
  }

  const ChemicalFeatures trpi_chem(&trpi_ag, trpi_ps);
  const ChemicalFeatures drug_chem(&drug_ag, drug_ps);
  
  // Summary evaluation
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
