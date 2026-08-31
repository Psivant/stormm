#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Chemistry/smiles.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/Math/summation.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/unit_test.h"
#include "../../src/Reporting/summary_file.h"

#ifndef STORMM_USE_HPC
using stormm::data_types::char2;
#endif
using stormm::errors::rtErr;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using stormm::stmath::sum;

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
  StopWatch timer;

  // Section 1
  section("Test atom interpretation");

  // Section 2
  section("Test atom chains");
  
  // Begin with major atoms, which may contain appended hydrogens
  section(1);
  const SmilesAtom radioactive_ethyl_carbon("13C", { 2 });
  check(radioactive_ethyl_carbon.getIsotopicNumber(), RelationalOperator::EQUAL, 13,
        "The isotopic number reported for what should be a radioactive ethyl carbon does not meet "
        "expectations.");
  check(radioactive_ethyl_carbon.getAtomicNumber(), RelationalOperator::EQUAL, 6, "The atomic "
        "number of a carbon atom is misinterpreted.");
  check(radioactive_ethyl_carbon.getAttachedHydrogenCount(), RelationalOperator::EQUAL, 2,
        "The number of hydrogens attached to an ethyl carbon is misassigned.");
  check(radioactive_ethyl_carbon.getFormalCharge(), RelationalOperator::EQUAL, 0, "The formal "
        "charge of a neutral carbon atom was not interpreted correctly.");
  const SmilesAtom oxide_oxygen("17O-", { 1 });
  check(oxide_oxygen.getIsotopicNumber(), RelationalOperator::EQUAL, 17,
        "The isotopic number reported for a rare isotope of oxygen does not meet expectations.");
  check(oxide_oxygen.getAtomicNumber(), RelationalOperator::EQUAL, 8, "The atomic number of an "
        "oxygen atom is misinterpreted.");
  check(oxide_oxygen.getAttachedHydrogenCount(), RelationalOperator::EQUAL, 0, "The number of "
        "hydrogens attached to a deprotonated oxygen atom is misassigned.");
  check(oxide_oxygen.getFormalCharge(), RelationalOperator::EQUAL, -1, "The formal charge "
        "of an ionized oxygen atom was not interpreted correctly.");
  const SmilesAtom organic_chlorine("37Cl", { 1 });
  check(organic_chlorine.getIsotopicNumber(), RelationalOperator::EQUAL, 37, "The isotopic number "
        "reported for a rare isotope of chlorine does not meet expectations.");
  check(organic_chlorine.getAtomicNumber(), RelationalOperator::EQUAL, 17, "The atomic number of "
        "a chlorine atom is misinterpreted.");
  check(organic_chlorine.getAttachedHydrogenCount(), RelationalOperator::EQUAL, 0, "A nonzero "
        "number of hydrogens was attached to an organic chlorine atom.");
  const SmilesAtom iron_sulfur_cluster("[56Fe++]", { 1, 1, 1, 1 });
  check(iron_sulfur_cluster.getIsotopicNumber(), RelationalOperator::EQUAL, 56, "The isotopic "
        "number reported for common iron does not meet expectations.");
  check(iron_sulfur_cluster.getAtomicNumber(), RelationalOperator::EQUAL, 26, "The atomic number "
        "of iron is misinterpreted.");
  check(iron_sulfur_cluster.getAttachedHydrogenCount(), RelationalOperator::EQUAL, 0, "A nonzero "
        "number of hydrogens was attached to an iron as might be at the center of an iron-sulfur "
        "cluster.");
  check(iron_sulfur_cluster.getFormalCharge(), RelationalOperator::EQUAL, 2, "The formal charge "
        "of an iron atom, as would be found in an iron-sulfur cluster, was not interpreted "
        "correctly.");
  const SmilesAtom chiral_carbon("C@TH1", { 1, 1, 1 });
  check(chiral_carbon.getChirality() == ChiralOrientation::RECTUS, "A chiral carbon did not "
        "receive the correct chiral designation.");
  const SmilesAtom chiral_carbon_ii("C@@", { 1, 1, 1 });
  check(chiral_carbon_ii.getChirality() == ChiralOrientation::SINISTER, "A chiral carbon did not "
        "receive the correct chiral designation when designated " +
        getEnumerationName(ChiralOrientation::SINISTER) + " in shorthand.");
  CHECK_THROWS(SmilesAtom chiral_carbon_iii("C@@", { 1, 1 }), "A chiral carbon was created "
               "without enough bonds to unique atoms.");
  SmilesAtom chiral_nitrogen("15N@@+", { 1, 1, 1 });
  check(chiral_nitrogen.getFormalCharge(), RelationalOperator::EQUAL, 1, "A chiral nitrogen's "
        "formal charge was not interpreted correctly.");
  CHECK_THROWS(SmilesAtom chiral_nitrogen_ii("15N@@+", { 1, 1 }), "A nitrogen with formal charge, "
               "which nonetheless could not be chiral, was given a chiral designation.");
  CHECK_THROWS(SmilesAtom chiral_nitrogen_iii("[15N@@+]", { 1, 1 }), "A nitrogen with formal "
               "charge, which nonetheless could not be chiral, was given a chiral designation "
               "when [ ] brackets were included in the token.");
  SmilesAtom chiral_nitrogen_iv("[15N@TH2+]", { 1, 1, 1 });
  check(chiral_nitrogen_iv.getFormalCharge(), RelationalOperator::EQUAL, 1, "A chiral nitrogen's "
        "formal charge was not interpreted correctly when its chirality was specified in long "
        "form.");
  check(chiral_nitrogen_iv.getChirality() == ChiralOrientation::SINISTER, "A chiral nitrogen's "
        "chirality was mis-interpreted when specified in long form.");

  // Non-branching molecules test the ability to create chains of major atoms with proper
  // protonation and bonding patterns.
  SmilesString ethane("CC");
  check(ethane.getAtomCount(), RelationalOperator::EQUAL, 8, "The total number of atoms found in "
        "ethane (SMILES \"" + ethane.getBasis() + "\") is incorrect.");
  SmilesString propanol("CCCO");
  check(propanol.getAtomCount(), RelationalOperator::EQUAL, 12, "The total number of atoms found "
        "in propanol (SMILES \"" + propanol.getBasis() + "\") is incorrect.");
  SmilesString tert_butanol("C(C)(C)CO");
  check(tert_butanol.getAtomCount(), RelationalOperator::EQUAL, 15, "The total number of atoms "
        "found in tert-butanol (SMILES \"" + tert_butanol.getBasis() + "\") is incorrect.");
  SmilesString tert_butanoate("C(C)(C)C([O-])=O");
  check(tert_butanoate.getAtomCount(), RelationalOperator::EQUAL, 13, "The total number of atoms "
        "found in tert-butanoate (SMILES \"" + tert_butanoate.getBasis() + "\") is incorrect.");
  SmilesString tert_butanoate_ii("C(C)(C)C(=O)[O-]");
  check(tert_butanoate_ii.getAtomCount(), RelationalOperator::EQUAL, 13, "The total number of "
        "atoms found in a second tert-butanoate representation (SMILES \"" +
        tert_butanoate_ii.getBasis() + "\") is incorrect.");
  SmilesString cyclohexane("C1CCCCC1");
  check(cyclohexane.getAtomCount(), RelationalOperator::EQUAL, 18, "The total number of atoms "
        "found in cyclohexane (SMILES \"" + cyclohexane.getBasis() + "\") is incorrect.");
  SmilesString cyclohexene("C1CCCCC=1");
  check(cyclohexene.getAtomCount(), RelationalOperator::EQUAL, 16, "The total number of atoms "
        "found in cyclohexene (SMILES \"" + cyclohexene.getBasis() + "\") is incorrect.");
  SmilesString cyclohexene_ii("C1CCCCC=1");
  check(cyclohexene_ii.getAtomCount(), RelationalOperator::EQUAL, 16, "The total number of atoms "
        "found in cyclohexene (SMILES \"" + cyclohexene_ii.getBasis() + "\") is incorrect.");
  SmilesString benzene("c1ccccc1");
  check(benzene.getAtomCount(), RelationalOperator::EQUAL, 12, "The total number of atoms "
        "found in benzene (SMILES \"" + benzene.getBasis() + "\") is incorrect.");
  check(benzene.getIsotopicNumber(8), RelationalOperator::EQUAL, 12, "The isotopic number "
        "(assigned from natural abundances) of a carbon in SMILES \"" + benzene.getBasis() +
        "\" is incorrect.");
  check(benzene.getBondCount(), RelationalOperator::EQUAL, 12, "The total number of bonds "
        "found in benzene (SMILES \"" + benzene.getBasis() + "\") is incorrect.");
  check(benzene.getBondCount(5), RelationalOperator::EQUAL, 1, "The total number of bonds "
        "involving one of the hydrogen atoms in benzene (SMILES \"" + benzene.getBasis() +
        "\") is incorrect.");
  check(benzene.getBondCount(6), RelationalOperator::EQUAL, 3, "The total number of bonds "
        "involving one of the carbon atoms in benzene (SMILES \"" + benzene.getBasis() +
        "\") is incorrect.");
  SmilesString benzene_ii("c1c=cc=cc=1");
  check(benzene_ii.getAtomCount(), RelationalOperator::EQUAL, 12, "The total number of atoms "
        "found in the kekule representation of benzene (SMILES \"" + benzene_ii.getBasis() +
        "\") is incorrect.");
  check(benzene_ii.getBondOrder(0, 1), RelationalOperator::EQUAL, 1, "The order of what should "
        "be a bond between an aromatic carbon and a hydrogen is incorrect (SMILES \"" +
        benzene_ii.getBasis() + "\").");
  check(benzene_ii.getBondOrder(0, 2), RelationalOperator::EQUAL, -1, "The order of what should "
        "be a bond in an aromatic ring is incorrect (SMILES \"" + benzene_ii.getBasis() + "\").");
  check(benzene_ii.getBondOrder(0, 10), RelationalOperator::EQUAL, -1, "The order of what should "
        "be a bond in an aromatic ring is incorrect (SMILES \"" + benzene_ii.getBasis() + "\").");
  check(benzene_ii.getAtomicNumber(8), RelationalOperator::EQUAL, 6, "The atomic number of an "
        "atom in SMILES \"" + benzene_ii.getBasis() + "\" is incorrect.");
  SmilesString two_methanes("C.C");
  check(two_methanes.getAtomCount(), RelationalOperator::EQUAL, 10, "The total number of atoms in "
        "a pair of methanes, specified by a non-bonded connection, is incorrect.");
  
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
