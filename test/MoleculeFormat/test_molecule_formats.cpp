#include <vector>
#include "copyright.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/scaling.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/MoleculeFormat/mdlmol.h"
#include "../../src/MoleculeFormat/molecule_format_enumerators.h"
#include "../../src/MoleculeFormat/molecule_parsing.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Parsing/textfile.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/render_molecule.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::constants;
using namespace stormm::data_types;
using namespace stormm::diskutil;
using namespace stormm::errors;
using namespace stormm::parse;
using namespace stormm::review;
using namespace stormm::structure;
using namespace stormm::topology;
using namespace stormm::trajectory;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Obtain environment variables or command-line input, if available
  TestEnvironment oe(argc, argv, TmpdirStatus::REQUIRED);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
  const char osc = osSeparator();

  // Section 1
  section("Test the MDL MOL entry");

  // Open an MDL MOL file with a single entry in V2000 format
  const std::string base_name = oe.getStormmSourcePath() + osc + "test" + osc + "MoleculeFormat";
  const std::string mdl_name      = base_name + osc + "sulfonamide.mol";
  const std::string sdf_name      = base_name + osc + "sulfonamide_rots.sdf";
  const std::string chemaxon_name = base_name + osc + "sdf_chemaxon.sdf";
  const std::string tether_name = base_name + osc + "tethered_atoms.sdf";
  const bool files_ready = (getDrivePathType(mdl_name) == DrivePathType::FILE &&
                            getDrivePathType(sdf_name) == DrivePathType::FILE &&
                            getDrivePathType(chemaxon_name) == DrivePathType::FILE &&
                            getDrivePathType(tether_name) == DrivePathType::FILE);
  const TestPriority do_tests = (files_ready) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (files_ready == false) {
    rtWarn("Files for MDL MOL format molecules were not found.  Check the STORMM source path, "
           "currently " + oe.getStormmSourcePath() + ", to ensure that it contains a valid "
           "subdirectory test/MoleculeFormat.  Subsequent tests will be skipped.");
  }
  TextFile mdl_a(mdl_name);
  MdlMol sulf_a(mdl_a);
  check(sulf_a.getAtomCount(), RelationalOperator::EQUAL, 16, "The atom count from an MDL MOL "
        "file is incorrect.", do_tests);
  check(sulf_a.getBondCount(), RelationalOperator::EQUAL, 15, "The bond count from an MDL MOL "
        "file is incorrect.", do_tests);
  const std::vector<double> sulf_a_ycrds = { -1.0845,  0.1506,  0.2311,  1.6123, -0.3787, -0.8405,
                                             -0.0173, -1.8467, -0.8392, -1.4937,  0.6911, -1.5531,
                                             -1.3843,  0.5307, -0.6815,  0.7010 };
  check(sulf_a.getCoordinates(CartesianDimension::Y), RelationalOperator::EQUAL,
        Approx(sulf_a_ycrds).margin(verytiny), "Coordinates for an MDL MOL entry were not taken "
        "up as expected.  The Y coordinate values do not match.", do_tests);
  const std::vector<int> sulf_a_znumbers = { 6, 7, 16, 8, 8, 6, 6, 1, 1, 1, 1, 1, 1, 1, 1, 1 };
  check(sulf_a.getAtomicNumbers(), RelationalOperator::EQUAL, sulf_a_znumbers, "Atomic numbers "
        "inferred for an MDL MOL entry do notmeet expectations.", do_tests);
  std::vector<MdlMol> sulf_rtm = readStructureDataFile(sdf_name);
  const std::vector<double> sulf_rtm_c_xcrds = {  2.2071,  1.4262, -0.0626, -0.5187,  0.1066,
                                                 -1.1940, -1.9484,  1.8501,  3.2455,  2.1797,
                                                  1.4358, -1.8772, -0.5739, -2.5472, -2.6168,
                                                 -1.2508 };
  check(sulf_rtm[2].getCoordinates(CartesianDimension::X), RelationalOperator::EQUAL,
        Approx(sulf_rtm_c_xcrds).margin(verytiny), "Coordinates for an MDL MOL entry within an "
        "SD file were not taken up as expected.  The X coordinate values do not match.", do_tests);

  // Read an SD file with V2000 MDL MOL entries containing some properties
  std::vector<MdlMol> chemaxon_mols = readStructureDataFile(chemaxon_name);
  check(chemaxon_mols[1].getPropertiesCount(), RelationalOperator::EQUAL, 1, "An incorrect number "
        "of properties were found in the second V2000 MDL MOL entry of " + chemaxon_name + ".",
        do_tests);
  const std::vector<int> fc_result = { chemaxon_mols[1].getFormalCharge(22),
                                       chemaxon_mols[1].getFormalCharge(22) };
  check(fc_result, RelationalOperator::EQUAL, std::vector<int>(2, 1), "Formal charges were not "
        "properly interpreted from an MDL MOL V2000 format property.", do_tests);

  // Read an SD file containing specific data items from which information can be extracted.
  std::vector<MdlMol> tether_mols = readStructureDataFile(tether_name);
  const std::vector<double> nrg_rep = realFromSdfDataItem("Internal_energy_repulsive",
                                                          tether_mols[0]);
  check(nrg_rep, RelationalOperator::EQUAL,
        Approx(std::vector<double>(1, 13.451076)).margin(1.0e-7), "A real value was not read "
        "correctly from one of the data items of an SDF file.", do_tests);
  const std::vector<int> hbond_donors = intFromSdfDataItem("HBond_Donors", tether_mols[0]);
  check(hbond_donors, RelationalOperator::EQUAL, std::vector<int>(1, 2), "An integer value was "
        "not read correctly from one of the data items of an SDF file.", do_tests);
  const std::vector<int> teth_atom_list = intFromSdfDataItem("TETHERED ATOMS", tether_mols[0]);
  const std::vector<int> teth_atom_ans = { 16, 11, 12, 13, 14, 15, 26, 24, 18, 17, 10,  9,  8,  7,
                                            6, 33, 32,  5, 25 };
  check(teth_atom_list, RelationalOperator::EQUAL, teth_atom_ans, "A list of integers pulled from "
        "an SD file data item does not meet expectations.", do_tests);
  const std::vector<double> real_storage_list = realFromSdfDataItem("REAL_NUMBERS",
                                                                    tether_mols[0]);
  const std::vector<double> real_storage_ans = { 0.1, 814.232, 97.354, -3.245 };
  check(real_storage_list, RelationalOperator::EQUAL, real_storage_ans, "A list of real numbers "
        "pulled from an SD file data item does not meet expectations.", do_tests);

  // Create an SD file de novo
  const std::string base_top_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::vector<std::string> sys_names = { "trpcage", "symmetry_L1", "symmetry_L1_vs",
                                               "bromobenzene_vs_iso", "symmetry_C3_in_water" };
  TestSystemManager tsm(base_top_name, "top", sys_names, base_crd_name, "inpcrd", sys_names);
  std::vector<int> recon_atom_counts(tsm.getSystemCount());
  std::vector<int> recon_bond_counts(tsm.getSystemCount());
  const std::vector<int>  recon_atom_counts_ans = { 304, 28, 28, 12, 25 };
  const std::vector<int>  recon_bond_counts_ans = { 310, 30, 30, 12, 26 };
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    CoordinateFrame sysi_cf = tsm.exportCoordinateFrame(i);
    MdlMol sysi_mdl(tsm.getTopologyPointer(i), sysi_cf);
    TextFile sysi_tf(sysi_mdl.writeMdl(MdlMolVersion::V2000), TextOrigin::RAM);
    MdlMol sysi_mdl_reconstructed(sysi_tf);
    recon_atom_counts[i] = sysi_mdl_reconstructed.getAtomCount();
    recon_bond_counts[i] = sysi_mdl_reconstructed.getBondCount();

    // Check the Trp-cage formal charges, which have been placed in a property.
    if (sys_names[i] == "trpcage") {
      check(sysi_mdl_reconstructed.getFormalCharge(302) == -1 ||
            sysi_mdl_reconstructed.getFormalCharge(303) == -1, "One of the terminal carboxylate "
            "oxygens of Trp-cage should have a formal charge of -1.", tsm.getTestingStatus());
      check(sysi_mdl_reconstructed.getFormalCharge(151) == 1, "The lysine head group charge on a "
            "residue in Trp-cage should have a formal charge of +1.", tsm.getTestingStatus());
      CHECK_THROWS_SOFT(sysi_mdl_reconstructed.getFormalCharge(10000), "A bogus atom index was "
                        "permitted into a formal charge query.", tsm.getTestingStatus());
    }
  }
  check(recon_atom_counts, RelationalOperator::EQUAL, recon_atom_counts_ans, "Atom counts "
        "recovered from MDL MOL objects built from topology and CoordinateFrame objects do not "
        "meet expectations.", tsm.getTestingStatus());
  check(recon_bond_counts, RelationalOperator::EQUAL, recon_bond_counts_ans, "Bond counts "
        "recovered from MDL MOL objects built from topology and CoordinateFrame objects do not "
        "meet expectations.", tsm.getTestingStatus());

  // Process a series of structures and test rendering methods
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    const NonbondedKit<double> nbk = tsm.getTopologyPointer(i)->getDoublePrecisionNonbondedKit();
    const ChemicalDetailsKit cdk = tsm.getTopologyPointer(i)->getChemicalDetailsKit();
    const int nw_atom = (nbk.natom + (sizeof(uint) * 8) - 1) / (sizeof(uint) * 8);
    const int nw_bond = (nbk.nb12_bounds[nbk.natom] + (sizeof(uint) * 8) - 1) / (sizeof(uint) * 8);
    std::vector<uint> atom_coverage(nw_atom, 0U), bond_coverage(nw_bond, 0U);
    for (int j = 0; j < cdk.nmol; j++) {
      const std::vector<int> atom_trace = traceBondLines(nbk, cdk, j, &atom_coverage,
                                                         &bond_coverage);
    }
    std::string missing_atoms;
    std::string missing_bonds;
    int n_missing_atoms = 0;
    int n_missing_bonds = 0;
    for (int j = 0; j < nbk.natom; j++) {
      if (readBitFromMask(atom_coverage, j) == 0) {
        if (n_missing_atoms < 8) {
          if (n_missing_atoms > 0) {
            missing_atoms += ", ";
          }
          missing_atoms += char4ToString(cdk.atom_names[j]);
        }
        n_missing_atoms++;
      }
      for (int k = nbk.nb12_bounds[j]; k < nbk.nb12_bounds[j + 1]; k++) {
        if (readBitFromMask(bond_coverage, k) == 0) {
          const size_t katom = nbk.nb12x[k];
          for (int m = nbk.nb12_bounds[katom]; m < nbk.nb12_bounds[katom + 1]; m++) {
            if (nbk.nb12x[m] == k && readBitFromMask(bond_coverage, m) == 0) {
              if (n_missing_bonds < 8) {
                if (n_missing_bonds > 0) {
                  missing_bonds += ", ";
                }
                missing_bonds += char4ToString(cdk.atom_names[j]) + " -- " +
                                 char4ToString(cdk.atom_names[katom]);
              }
              n_missing_bonds++;
            }
          }
        }
      }
    }
    check(n_missing_atoms, RelationalOperator::EQUAL, 0, "Some atoms were not covered in the "
          "bond tracing for the molecular system originating in topology " +
          getBaseName(tsm.getTopologyPointer(i)->getFileName()) + ".  Examples of missing "
          "atoms: " + missing_atoms + ".", tsm.getTestingStatus(i));
    check(n_missing_bonds, RelationalOperator::EQUAL, 0, "Some bonds were not covered in the "
          "bond tracing for the molecular system originating in topology " +
          getBaseName(tsm.getTopologyPointer(i)->getFileName()) + ".  Examples of missing "
          "bonds: " + missing_bonds + ".", tsm.getTestingStatus(i));
  }
  
  // Print results
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
