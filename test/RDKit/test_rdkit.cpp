#include "copyright.h"

#include "../../src/Reporting/summary_file.h"
#include "../../src/UnitTesting/unit_test.h"

#include <GraphMol/RDKitBase.h>
#include <GraphMol/GraphMol.h>
#include <GraphMol/SmilesParse/SmilesParse.h>

using stormm::errors::rtWarn;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Create an RDKit molecule and run basic tests.
//-------------------------------------------------------------------------------------------------
void testCreateRDKitMol() {
  std::shared_ptr<RDKit::ROMol> implicit_mol(RDKit::SmilesToMol("C"));
  auto molecule_was_created = (implicit_mol != nullptr);
  check(molecule_was_created, "RDKit ROMol could not be created from SMILES");
  check(implicit_mol->getNumAtoms(), RelationalOperator::EQUAL, 1, "RDKit ROMol has wrong number "
        "of implicit atoms");
  std::shared_ptr<RDKit::ROMol> explicit_mol( RDKit::MolOps::addHs(*implicit_mol) );
  check(explicit_mol->getNumAtoms(), RelationalOperator::EQUAL, 5, "RDKit ROMol has wrong number "
        "of explicit atoms");
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Create an RDMit molecule
  section(1);
  testCreateRDKitMol();

  // Print results
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
