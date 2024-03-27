#include "copyright.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/energy_enumerators.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Topology/atomgraph_refinement.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"
#include "test_amber_prmtop.h"

using stormm::card::Hybrid;
using stormm::constants::tiny;
using stormm::data_types::ulint;
#ifndef STORMM_USE_HPC
using stormm::data_types::char4;
#endif
using stormm::data_types::char4ToUint;
using stormm::data_types::uintToChar4;
using stormm::diskutil::DrivePathType;
using stormm::diskutil::getBaseName;
using stormm::diskutil::getDrivePathType;
using stormm::diskutil::osSeparator;
using stormm::energy::getEnumerationName;
using stormm::energy::VdwCombiningRule;
using stormm::errors::rtWarn;
using stormm::stmath::sum;
using stormm::parse::polyNumericVector;
using stormm::parse::stringToChar4;
using stormm::parse::TextFile;
using stormm::data_types::operator==;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Enumerated integer analyses of every topology in a list
//-------------------------------------------------------------------------------------------------
enum class IntInfoCode {
  MOLECULE_COUNT,             // Record the numbers of molecules
  LARGEST_RESIDUE_SIZE,       // Record the largest residue
  FIRST_MOLECULE_SIZE,        // Record the number of atoms in the first molecule
  LAST_SOLUTE_RESIDUE,        // Record the last solute residue
  LAST_SOLUTE_ATOM,           // Record the last solute atoms
  FIRST_SOLVENT_MOLECULE,     // Record the last solute atoms
  FIRST_SOLVENT_ATOM,         // Record the last solute atoms
  FIRST_RESIDUE_UPPER_LIMIT,  // Record the upper bound of atoms in the first residue
  SUM_ATOMIC_NUMBERS,         // Compute the sum of all atoms' atomic numbers
  ATOM_TYPE_COUNT,            // Record the number of total atom types
  MOST_COMMON_ATOM_TYPE,      // Determine the atom type that is most common in the topology
  TYPICAL_CMAP_DIMENSION,     // Determine the most common dimension for all CMAP terms
  SUM_BOND_INDICES,           // Compute the sum of up to 100 bond atom I and J indices
  SUM_ANGL_INDICES,           // Compute the sum of up to 100 bond angle atom I, J, and K indices
  SUM_DIHE_INDICES,           // Compute the sum of up to 100 dihedral atom I, J, K, and L indices
  SUM_UBRD_INDICES,           // Compute the sum of up to 100 Urey-Bradley atom I and K indices
  SUM_CIMP_INDICES,           // Compute the sum of up to 100 CHARMM improper atom indices
  SUM_CMAP_INDICES,           // Compute the sum of up to 100 CMAP atom I, J, K, L, and M indices
  CHARGE_TYPE_COUNT,          // Show the number of unique charges detected
  VIRTUAL_SITE_FRAME_SUM,     // Show the sum of virtual site frame atom indices
  CONSTRAINT_GROUP_COUNT,     // Show the number of constraint group
  CONSTRAINT_GROUP_SUM        // Show a sum invovling atom indices of constraint groups
};

//-------------------------------------------------------------------------------------------------
// Enumerated char4 name and type label analyses of every topology in a list
//-------------------------------------------------------------------------------------------------
enum class CharInfoCode {
  FIFTH_ATOM_NAME,        // Every topology has at least five atoms.  Get the fifth atom's name.
  FIRST_RESIDUE_NAME,     // Get the first residue's name
  MOST_COMMON_TYPE_NAME,  // Produce the name of the most common atom type
};

//-------------------------------------------------------------------------------------------------
// Create integer vectors for a particular property based on inputs from a series of topologies.
//
// Arguments:
//   topols:         The list of topologies
//   all_top_exist:  Flag to indicate that the topologies do indeed exist and should therefore be
//                   analyzed
//   analysis_code:  Codification of the analysis to perform (see the IntInfoCode enumerator)
//-------------------------------------------------------------------------------------------------
std::vector<int> integerAGProp(const std::vector<AtomGraph*> &topols, const bool all_top_exist,
                               const IntInfoCode analysis_code) {

  // Return a vector of blank data if no topologies have been filled out
  const int ntop = topols.size();
  if (all_top_exist == false) {
    return std::vector<int>(ntop, 0);
  }
  
  // Loop over all topologies and apply the requested analysis
  std::vector<int> result(ntop);
  for (int i = 0; i < ntop; i++) {
    switch (analysis_code) {
    case IntInfoCode::MOLECULE_COUNT:
      result[i] = topols[i]->getMoleculeCount();
      break;
    case IntInfoCode::LARGEST_RESIDUE_SIZE:
      result[i] = topols[i]->getLargestResidueSize();
      break;
    case IntInfoCode::FIRST_MOLECULE_SIZE:
      result[i] = topols[i]->getMoleculeLimits()[1];
      break;
    case IntInfoCode::LAST_SOLUTE_RESIDUE:
      result[i] = topols[i]->getLastSoluteResidue();
      break;
    case IntInfoCode::LAST_SOLUTE_ATOM:
      result[i] = topols[i]->getLastSoluteAtom();
      break;
    case IntInfoCode::FIRST_SOLVENT_MOLECULE:
      result[i] = topols[i]->getFirstSolventMolecule();
      break;
    case IntInfoCode::FIRST_SOLVENT_ATOM:
      result[i] = topols[i]->getFirstSolventAtom();
      break;
    case IntInfoCode::FIRST_RESIDUE_UPPER_LIMIT:
      result[i] = topols[i]->getResidueLimits().readHost(1);
      break;
    case IntInfoCode::SUM_ATOMIC_NUMBERS:
      result[i] = sum<int>(topols[i]->getAtomicNumber());
      break;
    case IntInfoCode::ATOM_TYPE_COUNT:
      {
        const std::vector<int> ati = topols[i]->getLennardJonesIndex();
        const size_t natom = topols[i]->getAtomCount();
        int maxtype = 0;
        for (size_t j = 0; j < natom; j++) {
          maxtype = std::max(ati[j], maxtype);
        }
        if (maxtype + 1 != topols[i]->getLJTypeCount()) {
          rtWarn("Topology " + topols[i]->getFileName() + " states that it contains " +
                 std::to_string(topols[i]->getLJTypeCount()) + "  atom types, but its "
                 "Lennard-Jones arrays only index " + std::to_string(maxtype + 1) + ".",
                 "test_amber_prmtop", "integerAGProp");
        }
        result[i] = maxtype + 1;
      }
      break;
    case IntInfoCode::MOST_COMMON_ATOM_TYPE:
      {
        const std::vector<int> ati = topols[i]->getLennardJonesIndex();
        const size_t natom = topols[i]->getAtomCount();
        const int n_lj_type = topols[i]->getLJTypeCount();
        std::vector<int> atype_populations(n_lj_type, 0);
        for (size_t j = 0; j < natom; j++) {
          if (ati[j] >= 0 && ati[j] < n_lj_type) {
            atype_populations[ati[j]] += 1;
          }
          else {
            rtWarn("Topology " + topols[i]->getFileName() + " contains a Lennard-Jones index for "
                   "atom " + std::to_string(j + 1) + " that is outside the number of unique atom "
                   "types it claims to have.", "test_amber_prmtop", "integerAGProp");
          }
        }
        int max_pop = 0;
        int most_common_type = -1;
        for (int i = 0; i < n_lj_type; i++) {
          if (i == 0 || max_pop < atype_populations[i]) {
            most_common_type = i;
            max_pop = atype_populations[i];
          }
        }
        result[i] = most_common_type;
      }
      break;
    case IntInfoCode::TYPICAL_CMAP_DIMENSION:
      {
        const int nsurf = topols[i]->getCmapSurfaceCount();
        const int ncmap = topols[i]->getCmapTermCount();

        // Test if there are CMAPs in the topology--otherwise leave the result as zero
        if (nsurf > 0 && ncmap > 0) {
          std::vector<int> surface_utilization(nsurf, 0);
          for (int j = 0; j < ncmap; j++) {
            const CmapTerm<double> tcm = topols[i]->getCmapTerm<double>(j);
            surface_utilization[tcm.surf_idx] += 1;
          }
          int max_util, max_loc;
          for (int j = 0; j < nsurf; j++) {
            if (j == 0 || max_util < surface_utilization[j]) {
              max_util = surface_utilization[j];
              max_loc = j;
            }
          }
          result[i] = topols[i]->getCmapDimension(max_loc);
        }
      }
      break;
    case IntInfoCode::SUM_BOND_INDICES:
      {
        const int nbond = topols[i]->getBondTermCount();
        int bsum = 0;
        for (int j = 0; j < nbond; j++) {
          const BondTerm<double> tbond = topols[i]->getBondTerm<double>(j);
          bsum += tbond.j_atom - tbond.i_atom;
        }
        result[i] = bsum;
      }
      break;
    case IntInfoCode::SUM_ANGL_INDICES:
      {
        const int nangl = topols[i]->getAngleTermCount();
        int asum = 0;
        for (int j = 0; j < nangl; j++) {
          const AngleTerm<double> tangl = topols[i]->getAngleTerm<double>(j);
          asum += (2 * tangl.j_atom) - tangl.i_atom - tangl.k_atom;
        }
        result[i] = asum;
      }
      break;
    case IntInfoCode::SUM_DIHE_INDICES:
      {
        const int ndihe = topols[i]->getDihedralTermCount();
        int hsum = 0;
        for (int j = 0; j < ndihe; j++) {
          const DihedralTerm<double> tdihe = topols[i]->getDihedralTerm<double>(j);
          hsum += (3 * tdihe.k_atom) - tdihe.i_atom - tdihe.j_atom - tdihe.l_atom;
        }
        result[i] = hsum;
      }
      break;
    case IntInfoCode::SUM_UBRD_INDICES:
      {
        const int nubrd = topols[i]->getUreyBradleyTermCount();
        int usum = 0;
        for (int j = 0; j < nubrd; j++) {
          const UreyBradleyTerm<double> tubrd = topols[i]->getUreyBradleyTerm<double>(j);
          usum += tubrd.k_atom - tubrd.i_atom;
        }
        result[i] = usum;
      }
      break;
    case IntInfoCode::SUM_CIMP_INDICES:
      {
        const int ncimp = topols[i]->getCharmmImprTermCount();
        int csum = 0;
        for (int j = 0; j < ncimp; j++) {
          const CharmmImprTerm<double> tcimp = topols[i]->getCharmmImprTerm<double>(j);
          csum += (3 * tcimp.k_atom) - tcimp.i_atom - tcimp.j_atom - tcimp.l_atom;
        }
        result[i] = csum;
      }
      break;
    case IntInfoCode::SUM_CMAP_INDICES:
      {
        const int ncmap = topols[i]->getCmapTermCount();
        int msum = 0;
        for (int j = 0; j < ncmap; j++) {
          const CmapTerm<double> tcmap = topols[i]->getCmapTerm<double>(j);
          msum += (4 * tcmap.k_atom) - tcmap.i_atom - tcmap.j_atom - tcmap.l_atom - tcmap.m_atom;
        }
        result[i] = msum;
      }
      break;
    case IntInfoCode::CHARGE_TYPE_COUNT:
      result[i] = topols[i]->getChargeTypeCount();
      break;
    case IntInfoCode::VIRTUAL_SITE_FRAME_SUM:
      {
        const int nvs = topols[i]->getVirtualSiteCount();
        int fsum = 0;
        for (int j = 0; j < nvs; j++) {
          const int parent_atom_idx = topols[i]->getVirtualSiteFrameAtom(j, 1);
          for (int k = 2; k <= 4; k++) {
            if (topols[i]->getVirtualSiteFrameAtom(j, k) >= 0) {
              fsum += topols[i]->getVirtualSiteFrameAtom(j, k) - parent_atom_idx;
            }
          }
        }
        result[i] = fsum;
      }
      break;
    case IntInfoCode::CONSTRAINT_GROUP_COUNT:
      result[i] = topols[i]->getConstraintGroupCount();
      break;
    case IntInfoCode::CONSTRAINT_GROUP_SUM:
      {
        const int ngrp = topols[i]->getConstraintGroupCount();
        int fsum = 0;
        for (int j = 0; j < ngrp; j++) {
          fsum += sum<int>(topols[i]->getConstraintGroupAtoms(j)) * (1 - (2 * (fsum >= 0)));
        }
        result[i] = fsum;
      }
      break;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Create char4 vectors for a particular property based on inputs from a series of topologies.
//
// Arguments:
//   topols:         The list of topologies
//   all_top_exist:  Flag to indicate that the topologies do indeed exist and should therefore be
//                   analyzed
//   analysis_code:  Codification of the analysis to perform (see the IntInfoCode enumerator)
//-------------------------------------------------------------------------------------------------
std::vector<char4> charAGProp(const std::vector<AtomGraph*> &topols, const bool all_top_exist,
                              const CharInfoCode analysis_code) {

  // Return a vector of blank data if no topologies have been filled out
  const int ntop = topols.size();
  if (all_top_exist == false) {
    char4 blank = { 'N', 'O', 'N', 'E' };
    return std::vector<char4>(ntop, blank);
  }

  // Loop over all topologies and apply the requested analysis
  std::vector<char4> result(ntop);
  for (int i = 0; i < ntop; i++) {
    switch (analysis_code) {
    case CharInfoCode::FIFTH_ATOM_NAME:
      if (topols[i]->getAtomCount() >= 5) {
        result[i] = topols[i]->getAtomName(4);
      }
      break;
    case CharInfoCode::FIRST_RESIDUE_NAME:
      result[i] = topols[i]->getResidueName(0);
      break;
    case CharInfoCode::MOST_COMMON_TYPE_NAME:
      {
        const int natom = topols[i]->getAtomCount();
        const Hybrid<char4>& iag_at = topols[i]->getAtomType();
        std::vector<char4> all_type_names = iag_at.readHost();
        std::vector<bool> found(natom, false);
        int max_repeats = 0;
        for (int j = 0; j < natom; j++) {
          if (found[j] == true) {
            continue;
          }
          int n_repeats = 0;
          for (int k = j; k < natom; k++) {
            if (all_type_names[k] == all_type_names[j]) {
              found[k] = true;
              n_repeats++;
            }
          }
          if (n_repeats > max_repeats) {
            result[i] = all_type_names[j];
            max_repeats = n_repeats;
          }
        }
      }
      break;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
// Streamline the process of printing a list of results for all topologies.
//
// Arguments:
//   topols:      The list of topologies
//   method:      The (integer, real, or string) analysis method to use
//-------------------------------------------------------------------------------------------------
void showAnswer(const std::vector<AtomGraph*> &topols, const IntInfoCode method,
                const char* hint) {
  const std::vector<int> answer = integerAGProp(topols, true, method);
  const int ntop = topols.size();
  printf("%s = [\n", hint);
  for (int i = 0; i < ntop; i++) {
    printf("%-60.60s  %8d\n", topols[i]->getFileName().c_str(), answer[i]);
  }
  printf("];\n");
  printf("%s = [\n", hint);
  int j = 0;
  for (int i = 0; i < ntop; i++) {
    printf("%8d", answer[i]);
    if (i < ntop - 1) {
      printf(", ");
    }
    j++;
    if (j == 11) {
      printf("\n");
      j = 0;
    }
  }
  if (j > 0) {
    printf("\n");
  }
  printf("];\n");
}

template <typename T>
void showAnswer(const std::vector<AtomGraph*> &topols, const RealInfoCode method,
                const char* hint) {
  const std::vector<double> answer = realAGProp<T>(topols, true, method);
  const int ntop = topols.size();
  printf("%s = [\n", hint);
  for (int i = 0; i < ntop; i++) {
    printf("%-60.60s  %12.8lf\n", topols[i]->getFileName().c_str(), answer[i]);
  }
  printf("];\n");
  printf("%s = [\n", hint);
  int j = 0;
  for (int i = 0; i < ntop; i++) {
    printf("%12.8lf", answer[i]);
    if (i < ntop - 1) {
      printf(", ");
    }
    j++;
    if (j == 7) {
      printf("\n");
      j = 0;
    }
  }
  if (j > 0) {
    printf("\n");
  }
  printf("];\n");
}

void showAnswer(const std::vector<AtomGraph*> &topols, const CharInfoCode method,
                const char* hint) {
  const std::vector<char4> answer = charAGProp(topols, true, method);
  const int ntop = topols.size();
  printf("%s = [\n", hint);
  for (int i = 0; i < ntop; i++) {
    printf("%-60.60s  %c%c%c%c\n", topols[i]->getFileName().c_str(), answer[i].x, answer[i].y,
           answer[i].z, answer[i].w);
  }
  printf("];\n");
  printf("%s = [\n", hint);
  int j = 0;
  for (int i = 0; i < ntop; i++) {
    printf("%c%c%c%c", answer[i].x, answer[i].y, answer[i].z, answer[i].w);
    if (i < ntop - 1) {
      printf(", ");
    }
    j++;
    if (j == 16) {
      printf("\n");
      j = 0;
    }
  }
  if (j > 0) {
    printf("\n");
  }
  printf("];\n");
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // The temporary directory will not be needed.  Results will be compared against snapshots.
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
  
  // Section 1
  section("Atomic details");

  //Section 2
  section("Valence terms");

  // Section 3
  section("Non-bonded information");

  // Section 4
  section("Analysis on complete topologies");

  // Section 5
  section("Error capture");

  // Section 6
  section("C++ manipulation of topologies");
  
  // Read many topologies.  Test one of them as an indicator whether they are all available.
  const char osc = osSeparator();
  const std::string base_top_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string tip3p_top_name = base_top_name + osc + "tip3p.top";
  const std::string tip4p_top_name = base_top_name + osc + "tip4p.top";
  const std::string tip4p_error_top_name = base_top_name + osc + "tip4p_error.top";
  const std::string tip5p_top_name = base_top_name + osc + "tip5p.top";
  const std::string trpcage_top_name = base_top_name + osc + "trpcage.top";
  const std::string trpcage_in_water_top_name = base_top_name + osc + "trpcage_in_water.top";
  const std::string trpcage_noz_name = base_top_name + osc + "trpcage_no_z.top";
  const std::string brbz_top_name = base_top_name + osc + "bromobenzene.top";
  const std::string brbz_vs_top_name = base_top_name + osc + "bromobenzene_vs.top";
  const std::string dhfr_top_name = base_top_name + osc + "dhfr_cmap.top";
  const std::string camp_top_name = base_top_name + osc + "camp_in_water.top";
  const std::string dna_top_name = base_top_name + osc + "dna.top";
  const std::string rna_top_name = base_top_name + osc + "rna.top";
  const std::string ubiquitin_top_name = base_top_name + osc + "ubiquitin.top";
  const std::string drug_top_name = base_top_name + osc + "drug_example.top";
  const std::string ala_dipeptide_top_name = base_top_name + osc + "ala_dipeptide.top";
  const bool all_top_exist =
    (getDrivePathType(tip3p_top_name) == DrivePathType::FILE &&
     getDrivePathType(tip4p_top_name) == DrivePathType::FILE &&
     getDrivePathType(tip4p_error_top_name) == DrivePathType::FILE &&
     getDrivePathType(tip5p_top_name) == DrivePathType::FILE &&
     getDrivePathType(trpcage_top_name) == DrivePathType::FILE &&
     getDrivePathType(trpcage_in_water_top_name) == DrivePathType::FILE &&
     getDrivePathType(trpcage_noz_name) == DrivePathType::FILE &&
     getDrivePathType(brbz_top_name) == DrivePathType::FILE &&
     getDrivePathType(brbz_vs_top_name) == DrivePathType::FILE &&
     getDrivePathType(dhfr_top_name) == DrivePathType::FILE &&
     getDrivePathType(camp_top_name) == DrivePathType::FILE &&
     getDrivePathType(dna_top_name) == DrivePathType::FILE &&
     getDrivePathType(rna_top_name) == DrivePathType::FILE &&
     getDrivePathType(ubiquitin_top_name) == DrivePathType::FILE &&
     getDrivePathType(drug_top_name) == DrivePathType::FILE &&
     getDrivePathType(ala_dipeptide_top_name) == DrivePathType::FILE);
  const TestPriority top_check = (all_top_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  AtomGraph tip3p, tip4p, tip4p_error, tip5p, trpcage, trpcage_water, trpcage_noz, brbz, brbz_vs,
    dhfr, camp, dna, rna, ubiquitin, drug, ala_dipeptide;
  if (all_top_exist) {
    ExceptionResponse policy = ExceptionResponse::SILENT;
    tip3p.buildFromPrmtop(tip3p_top_name, policy);
    tip4p.buildFromPrmtop(tip4p_top_name, policy);
    tip4p_error.buildFromPrmtop(tip4p_error_top_name, policy);
    tip5p.buildFromPrmtop(tip5p_top_name, policy);
    trpcage.buildFromPrmtop(trpcage_top_name, policy);
    trpcage_water.buildFromPrmtop(trpcage_in_water_top_name, policy);
    trpcage_noz.buildFromPrmtop(trpcage_noz_name, policy);
    brbz.buildFromPrmtop(brbz_top_name, policy);
    brbz_vs.buildFromPrmtop(brbz_vs_top_name, policy);
    dhfr.buildFromPrmtop(dhfr_top_name, policy);
    camp.buildFromPrmtop(camp_top_name, policy);
    dna.buildFromPrmtop(dna_top_name, policy);
    rna.buildFromPrmtop(rna_top_name, policy);
    ubiquitin.buildFromPrmtop(ubiquitin_top_name, policy);
    drug.buildFromPrmtop(drug_top_name, policy);
    ala_dipeptide.buildFromPrmtop(ala_dipeptide_top_name, policy);
  }
  else {
    const std::string stormmsrc_base_name(std::string("${STORMM_SOURCE}") + osc + "test" + osc +
                                        "Topology");
    rtWarn("One or more topologies required by subsequent tests were not found.  Check the "
           "$STORMM_SOURCE environment variable to make sure it indicates the location of the "
           "STORMM source tree, with src/ and test/ as subdirectories.  This test program will "
           "skip nearly all tests until files such as " + getBaseName(tip3p_top_name) + ", " +
           getBaseName(tip4p_top_name) + ", and " + getBaseName(brbz_vs_top_name) +
           "can be found in directory " + stormmsrc_base_name + ".", "test_amber_prmtop");
  }
  const bool snapshots_exist =
    (getDrivePathType(base_top_name + osc + "tip3p_atoms.m") == DrivePathType::FILE &&
     getDrivePathType(base_top_name + osc + "trpcage_details.m") == DrivePathType::FILE &&
     getDrivePathType(base_top_name + osc + "ala_charges.m") == DrivePathType::FILE);
  const TestPriority snap_check = (snapshots_exist) ? TestPriority::CRITICAL : TestPriority::ABORT;
  if (snapshots_exist == false && oe.takeSnapshot() != SnapshotOperation::SNAPSHOT) {
    const std::string stormmsrc_base_name(std::string("${STORMM_SOURCE}") + osc + "test" + osc +
                                        "Topology");
    rtWarn("Snapshot required by subsequent tests were not found.  Check the $STORMM_SOURCE "
           "environment variable.  Some subsequent tests will be skipped until the files " +
           base_top_name + osc + "tip3p_atoms.m," + base_top_name + osc + "trpcage_details.m," +
           " and " + base_top_name + osc + "ala_charges.m can be found in directory " +
           stormmsrc_base_name + ".", "test_amber_prmtop");
  }

  // Check various atomic details
  section(1);
  const int nt_atom = (all_top_exist) ? tip3p.getDescriptor(TopologyDescriptor::ATOM_COUNT) : 0; 
  check(nt_atom, RelationalOperator::EQUAL, 768, "Atom count for TIP3P water topology is "
        "incorrect.", top_check);
  PolyNumeric one_pn;
  one_pn.i = 0;
  const std::vector<PolyNumeric> tip3p_zs = (all_top_exist) ?
                                            polyNumericVector(tip3p.getAtomicNumber()) :
                                            std::vector<PolyNumeric>(1, one_pn);
  snapshot(base_top_name + osc + "tip3p_atoms.m", tip3p_zs, "z_numbers", NumberFormat::INTEGER,
           "Atomic numbers for the tip3p topology do not meet expectations.", oe.takeSnapshot(),
           1.0e-4, 1.0e-8, PrintSituation::OVERWRITE, snap_check);
  check(static_cast<ulint>(TopologyDescriptor::N_VALUES), RelationalOperator::EQUAL,
        static_cast<ulint>(SanderDescriptor::N_VALUES), "The number of enumerations in "
        "TopologyDescriptor and SanderDescriptor do not match.");
  const int nbrbz_atom = (all_top_exist) ? brbz.getDescriptor(SanderDescriptor::NATOM) : 0;
  check(nbrbz_atom, RelationalOperator::EQUAL, 12, "An incorrect number of atoms is recorded in "
        "Amber prmtop file " + brbz.getFileName() + ".", top_check);
  const int nbrbz_vs = (all_top_exist) ? brbz_vs.getDescriptor(SanderDescriptor::NUMEXTRA) : 0;
  check(nbrbz_vs, RelationalOperator::EQUAL, 13, "An incorrect number of virtual sites is found "
        "in Amber prmtop file " + brbz.getFileName() + ".", top_check);
  check(tip3p.getResidueIndex(56), RelationalOperator::EQUAL, 18, "The residue index of a "
        "particular atom was not computed correctly.");
  const int ntrp_atom = trpcage.getAtomCount();
  std::vector<int> ridx_result(ntrp_atom);
  for (int i = 0; i < ntrp_atom; i++) {
    ridx_result[i] = trpcage.getResidueIndex(i);
  }
  check(trpcage.getResidueIndex(), RelationalOperator::EQUAL, ridx_result, "A scan of all atoms "
        "in the Trp-cage system (isolated boundary conditions) did not produce complete "
        "agreement when computing residue indices for individual atoms versus taking the "
        "limit-based result en masse.");
  trpcage.setWaterResidueName("WAT");
  check(trpcage.getWaterResidueName(), RelationalOperator::EQUAL, std::string("WAT "),
        "A topology's water residue name was not updated as expected.");

  // Check valence term details
  section(2);
  const int nbrbz_phi = (all_top_exist) ? brbz.getDescriptor(SanderDescriptor::NPHIA) +
                                          brbz.getDescriptor(SanderDescriptor::NPHIH) : 0;
  check(nbrbz_phi, RelationalOperator::EQUAL, 30, "An incorrect number of total dihedrals is "
        "recorded in Amber prmtop file " + brbz.getFileName() + ".", top_check);
  const int ntrp_ang = (all_top_exist) ? trpcage_water.getDescriptor(SanderDescriptor::NTHETH) +
                                         trpcage_water.getDescriptor(SanderDescriptor::NTHETA) : 0;
  check(ntrp_ang, RelationalOperator::EQUAL, 565, "An incorrect total number of angles is "
        "recorded in Amber prmtop file " + trpcage_water.getFileName() + ".", top_check);
  const int nt_bondh = (all_top_exist) ? tip3p.getDescriptor(SanderDescriptor::NBONH) : 0;
  check(nt_bondh, RelationalOperator::EQUAL, 768, "An incorrect number of bonds involving "
        "hydrogen atoms is recorded in Amber prmtop file " + tip3p.getFileName(), top_check);
  const int nt_bonda = (all_top_exist) ? tip3p.getDescriptor(SanderDescriptor::NBONA) : -1;
  check(nt_bonda, RelationalOperator::EQUAL, 0, "An incorrect number of bonds involving heavy "
        "atoms is recorded in Amber prmtop file " + tip3p.getFileName(), top_check);
  const std::vector<double> cmap_stencil_answer = {
    -0.000000000000,  0.803847577293, -0.215390309173,  0.057713659398, -0.015464328418,
     0.004143654273, -0.001110288675,  0.000297500427, -0.000079713033,  0.000021351705,
    -0.000005693788,  0.000001423447, -0.000000000000, -0.000001423447,  0.000005693788,
    -0.000021351705,  0.000079713033, -0.000297500427,  0.001110288675, -0.004143654273,
     0.015464328418, -0.057713659398,  0.215390309173, -0.803847577293 };
  check(cubicSplineDerivativeStencil(24), RelationalOperator::EQUAL,
        Approx(cmap_stencil_answer).margin(tiny), "The cubic spline derivative stencil (needed "
        "for solving CMAP derivatives) was not computed correctly.");

  // Check non-bonded properties
  section(3);
  const std::vector<PolyNumeric> trp_q = (all_top_exist) ?
                                         polyNumericVector(trpcage.getPartialCharge<double>()) :
                                         std::vector<PolyNumeric>(1, one_pn);
  snapshot(base_top_name + osc + "trpcage_details.m", trp_q, "charges",
           NumberFormat::STANDARD_REAL, "Smoothed charges for the trpcage topology do not meet "
           "expectations.", oe.takeSnapshot(), 1.0e-8, 1.0e-12, PrintSituation::OVERWRITE,
           snap_check);
  const std::vector<PolyNumeric> trp_znum = (all_top_exist) ?
                                            polyNumericVector(trpcage_noz.getAtomicNumber()) :
                                            std::vector<PolyNumeric>(1, one_pn);
  snapshot(base_top_name + osc + "trpcage_details.m", trp_znum, "inferred_atomic_numbers",
           NumberFormat::INTEGER, "Inferred atomic numbers for the trpcage (mass "
           "repartitioned, no Z numbers) topology do not meet expectations.", oe.takeSnapshot(),
           1.0e-4, 1.0e-8, PrintSituation::APPEND, snap_check);
  const std::vector<PolyNumeric> alad_q =
     (all_top_exist) ? polyNumericVector(ala_dipeptide.getPartialCharge<double>()) :
                       std::vector<PolyNumeric>(1, one_pn);
  snapshot(base_top_name + osc + "ala_charges.m", alad_q, "charges", NumberFormat::STANDARD_REAL,
           "Smoothed charges for the alanine dipeptide topology do not meet "
           "expectations.", oe.takeSnapshot(), 1.0e-8, 1.0e-12, PrintSituation::OVERWRITE,
           snap_check);

  // Copy and move constructor testing
  section(6);
  std::vector<AtomGraph> dry_topologies_pb;
  dry_topologies_pb.reserve(6);
  dry_topologies_pb.push_back(brbz);
  if (all_top_exist) {
    dry_topologies_pb.push_back(AtomGraph(brbz_vs_top_name, ExceptionResponse::SILENT));
  }
  else {
    dry_topologies_pb.push_back(AtomGraph());
  }
  if (all_top_exist) {
    dry_topologies_pb.push_back(AtomGraph(trpcage_top_name, ExceptionResponse::SILENT));
  }
  else {
    dry_topologies_pb.push_back(AtomGraph());
  }
  dry_topologies_pb.push_back(trpcage_noz);
  if (all_top_exist) {
    dry_topologies_pb.push_back(AtomGraph(dhfr_top_name, ExceptionResponse::SILENT));
  }
  else {
    dry_topologies_pb.push_back(AtomGraph());
  }
  dry_topologies_pb.push_back(dna);
  const std::vector<WaterModel> dry_water_answers = { WaterModel::NONE, WaterModel::NONE,
                                                      WaterModel::NONE, WaterModel::NONE,
                                                      WaterModel::NONE, WaterModel::NONE };
  bool dryness_check = true;
  for (size_t i = 0; i < dry_topologies_pb.size(); i++) {
    if (all_top_exist) {
      dryness_check = (dryness_check &&
                       identifyWaterModel(dry_topologies_pb[i]) == WaterModel::NONE);
    }
  }
  check(dryness_check, "Water models were identified in a Standard Template Library vector of dry "
        "topologies created by the push_back() method.", top_check);
  bool blank_topology_vector_growth_success;
  try {
    std::vector<AtomGraph> agvn;
    for (int i = 0; i < 16; i++) {
      agvn.emplace_back();
    }
    blank_topology_vector_growth_success = true;
  }
  catch (std::runtime_error) {
    blank_topology_vector_growth_success = false;
  }
  check(blank_topology_vector_growth_success, "A Standard Template Library vector of blank "
        "AtomGraph objects does not grow properly.  Most likely, one of the constructors or "
        "assignment operators is having trouble with some POINTER-kind Hybrid object not yet "
        "pointing to anything as there is no content.");
  
  // Test the move assignment operator by erasing elements from a list of topologies.
  std::vector<AtomGraph> clones = dry_topologies_pb;
  clones.erase(clones.begin() + 1);
  clones.erase(clones.begin() + 1);
  check(clones.size(), RelationalOperator::EQUAL, 4, "The Standard Template Library erase() "
        "method does not work as intended.  This may indicate a problem with the AtomGraph move "
        "assignment operator, or another move assignment operator at a deeper level.", top_check);

  // Test the self pointers and their resilience against various copy and move operations.
  const AtomGraph *tip3p_self_ptr = tip3p.getSelfPointer();
  const WaterModel should_be_tip3p = identifyWaterModel(*tip3p_self_ptr);
  const AtomGraph& tip4p_ref = tip4p;
  check(tip4p_ref.getSelfPointer() == tip4p.getSelfPointer(), "A reference to a topology does not "
        "return the same self pointer as the original object.", top_check);
  check(should_be_tip3p == WaterModel::TIP3P, "The self-pointer of a TIP3P water system's "
        "topology does not correctly identify the water model.", top_check);
  const std::vector<int> clone_system_sizes_ans = { 12, 304, 2489, 632 };
  std::vector<int> clone_system_sizes(clones.size());
  for (size_t i = 0; i < clones.size(); i++) {
    const AtomGraph *iclone_ptr = clones[i].getSelfPointer();
    clone_system_sizes[i] = iclone_ptr->getAtomCount();
  }
  check(clone_system_sizes, RelationalOperator::EQUAL, clone_system_sizes_ans, "AtomGraph "
        "self-pointers are no longer reporting accurate atom counts after going through an STL "
        "vector manipulation that invokes the move assignment operator.", top_check);
  
  // Create a vector of many topologies for subsequent analyses.  This will implicitly test copy
  // and move constructors, as the vector of all topologies will be created based on a product of
  // many such operations.
  std::vector<AtomGraph> all_top_stack;
  if (all_top_exist) {
    all_top_stack.push_back(tip3p);
    all_top_stack.push_back(AtomGraph(tip4p_top_name, ExceptionResponse::SILENT));
    all_top_stack.push_back(AtomGraph(tip4p_error_top_name, ExceptionResponse::SILENT));
    all_top_stack.push_back(tip5p);
    all_top_stack.push_back(AtomGraph(trpcage_top_name, ExceptionResponse::SILENT));
    all_top_stack.push_back(AtomGraph(trpcage_in_water_top_name, ExceptionResponse::SILENT));
    all_top_stack.push_back(trpcage_noz);
    all_top_stack.push_back(AtomGraph(brbz_top_name, ExceptionResponse::SILENT));
    all_top_stack.push_back(AtomGraph(brbz_vs_top_name, ExceptionResponse::SILENT));
    all_top_stack.push_back(dhfr);
    all_top_stack.push_back(AtomGraph(camp_top_name, ExceptionResponse::SILENT));
    all_top_stack.push_back(AtomGraph(dna_top_name, ExceptionResponse::SILENT));
    all_top_stack.push_back(rna);
    all_top_stack.push_back(AtomGraph(ubiquitin_top_name, ExceptionResponse::SILENT));
    all_top_stack.push_back(AtomGraph(drug_top_name, ExceptionResponse::SILENT));
  }
  else {
    for (int i = 0; i < 15; i++) {
      all_top_stack.push_back(AtomGraph());
    }
  }
  section(4);
  const std::vector<AtomGraph*> all_topologies = { &all_top_stack[ 0], &all_top_stack[ 1],
                                                   &all_top_stack[ 2], &all_top_stack[ 3],
                                                   &all_top_stack[ 4], &all_top_stack[ 5],
                                                   &all_top_stack[ 6], &all_top_stack[ 7],
                                                   &all_top_stack[ 8], &all_top_stack[ 9],
                                                   &all_top_stack[10], &all_top_stack[11],
                                                   &all_top_stack[12], &all_top_stack[13],
                                                   &all_top_stack[14] };
  const std::vector<int> molecule_count_answer = { 256, 256, 256, 216, 1, 1561, 1, 1, 1, 1, 1173,
                                                   2, 1185, 481, 1225 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::MOLECULE_COUNT),
        RelationalOperator::EQUAL, molecule_count_answer, "Molecule counts for topologies are "
        "incorrect.", top_check);
  const std::vector<int> largest_residue_answer = { 3, 4, 4, 5, 24, 24, 24, 12, 25, 25, 33, 33, 35,
                                                    24, 53 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::LARGEST_RESIDUE_SIZE),
        RelationalOperator::EQUAL, largest_residue_answer, "Largest residue sizes for topologies "
        "are incorrect.", top_check);
  const std::vector<int> first_mol_size_answer = { 3, 4, 4, 5, 304, 304, 304, 12, 25, 2489, 33,
                                                   316, 320, 1185, 53 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::FIRST_MOLECULE_SIZE),
        RelationalOperator::EQUAL, first_mol_size_answer, "First molecule sizes for topologies "
        "are incorrect.", top_check);
  const std::vector<int> last_solute_res_answer = { -1, -1, -1, -1, 19, 19, 19, 0, 0, 158, 0, 19,
                                                    37, 72, 0 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::LAST_SOLUTE_RESIDUE),
        RelationalOperator::EQUAL, last_solute_res_answer, "Last solute residues for topologies "
        "are incorrect.", top_check);
  const std::vector<int> last_solute_atom_answer = { -1, -1, -1, -1, 303, 303, 303, 11, 24, 2488,
                                                     32, 631, 657, 1184, 52 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::LAST_SOLUTE_ATOM),
        RelationalOperator::EQUAL, last_solute_atom_answer, "Last solute atoms for topologies "
        "are incorrect.", top_check);
  const std::vector<int> first_solvent_mol_answer = { 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 20, 1,
                                                      1 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::FIRST_SOLVENT_MOLECULE),
        RelationalOperator::EQUAL, first_solvent_mol_answer, "First solvent molecules for "
        "topologies are incorrect.", top_check);
  const std::vector<int> first_solvent_atom_answer = { 0, 0, 0, 0, -1, 304, -1, -1, -1, -1, 33,
                                                       -1, 658, 1185, 53 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::FIRST_SOLVENT_ATOM),
        RelationalOperator::EQUAL, first_solvent_atom_answer, "First solvent atoms for topologies "
        "are incorrect.", top_check);
  const std::vector<int> first_res_uplim_answer = { 3, 4, 4, 5, 16, 16, 16, 12, 25, 19, 33, 31,
                                                    29, 19, 53 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::FIRST_RESIDUE_UPPER_LIMIT),
        RelationalOperator::EQUAL, first_res_uplim_answer, "First residue upper atom limits for "
        "topologies are incorrect.", top_check);
  const std::vector<int> sum_znum_answer = { 2560, 2560, 2560, 2160, 1159, 16759, 1159, 76, 76,
                                             9595, 11889, 3122, 15098, 9199, 12414 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::SUM_ATOMIC_NUMBERS),
        RelationalOperator::EQUAL, sum_znum_answer, "Sums of atomic numbers from topologies do "
        "not meet expectations.", top_check);
  const std::vector<int> atom_type_count_answer = { 2, 3, 3, 2, 13, 14, 12, 3, 4, 33, 14, 15, 16,
                                                    16, 13 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::ATOM_TYPE_COUNT),
        RelationalOperator::EQUAL, atom_type_count_answer, "The number of Lennard-Jones atom "
        "types in topologies do not meet expectations.  This check include as implicit test of "
        "whether the number of Lennard-Jones types matches the stated atom type count for the "
        "topology.", top_check);
  const std::vector<int> common_atom_type_answer = { 1, 1, 1, 1, 4, 10, 4, 0, 3, 2, 13, 2, 0, 10,
                                                     12 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::MOST_COMMON_ATOM_TYPE),
        RelationalOperator::EQUAL, common_atom_type_answer, "The most common Lennard-Jones "
        "indices found in each topology do not meet expectations.", top_check);
  const std::vector<int> common_cmap_answer = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 24, 0, 0, 0, 24, 0 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::TYPICAL_CMAP_DIMENSION),
        RelationalOperator::EQUAL, common_cmap_answer, "The most common CMAP surfaces "
        "found in each topology do not meet expectations.", top_check);
  section(2);
  const std::vector<int> bond_index_sum_answer = { -1016, -1778, -1778, 1080, 756, -5484, 756, 26,
                                                   26, -1092, -4588, 2094, -2598, -403, -4637 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::SUM_BOND_INDICES),
        RelationalOperator::EQUAL, bond_index_sum_answer, "The bond index sums did not meet "
        "expectations for all topologies.", top_check);
  const std::vector<int> angl_index_sum_answer = { 0, 0, 0, 0, 888, 888, 888, 30, 30, 7090, 149,
                                                   3340, 3284, 3259, 231 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::SUM_ANGL_INDICES),
        RelationalOperator::EQUAL, angl_index_sum_answer, "The bond angle index sums did not meet "
        "expectations for all topologies.", top_check);
  const std::vector<int> dihe_index_sum_answer = { 0, 0, 0, 0, -7642, -7642, -5194, -74, -74,
                                                   -39529, -463, -13161, -14817, -34988, -1932 };

  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::SUM_DIHE_INDICES),
        RelationalOperator::EQUAL, dihe_index_sum_answer, "The dihedral index sums did not meet "
        "expectations for all topologies.", top_check);
  const std::vector<int> ubrd_index_sum_answer = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 7106, 0, 0, 0, 0,
                                                   0 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::SUM_UBRD_INDICES),
        RelationalOperator::EQUAL, ubrd_index_sum_answer, "The Urey-Bradley index sums did not "
        "meet expectations for all topologies.", top_check);
  const std::vector<int> cimp_index_sum_answer = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 4130, 0, 0, 0, 0,
                                                   0 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::SUM_CIMP_INDICES),
        RelationalOperator::EQUAL, cimp_index_sum_answer, "The CHARMM improper dihedral index "
        "sums did not meet expectations for all topologies.", top_check);
  const std::vector<int> cmap_index_sum_answer = {0, 0, 0, 0, 0, 0, 0, 0, 0, -2926, 0, 0, 0, -1344,
                                                   0 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::SUM_CMAP_INDICES),
        RelationalOperator::EQUAL, cmap_index_sum_answer, "The CMAP index sums did not meet "
        "expectations for all topologies.", top_check);
  const std::vector<int> charge_count_answer = { 2, 3, 5, 3, 148, 151, 141, 9, 9, 56, 33, 76, 79,
                                                 175, 40 };
  section(1);
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::CHARGE_TYPE_COUNT),
        RelationalOperator::EQUAL, charge_count_answer, "Charge parameter counts were not "
        "computed as expected.", top_check);
  const std::vector<int> vs_frame_sum_answer = { 0, 768, 768, 1296, 0, 0, 0, 0, -1, 0, 0, 0, 0,
                                                 1440, 0 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::VIRTUAL_SITE_FRAME_SUM),
        RelationalOperator::EQUAL, vs_frame_sum_answer, "Virtual site frame sum was not computed "
        "as expected.", top_check);
  const std::vector<int> cnst_grp_count_answer = { 0, 0, 0, 0, 96, 96, 96, 5, 5, 790, 9, 164,
                                                   184, 366, 18 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::CONSTRAINT_GROUP_COUNT),
        RelationalOperator::EQUAL, cnst_grp_count_answer, "Constraint group counts did not match "
        "the expected results.", top_check);
  const std::vector<int> cnst_grp_sum_answer = { 0, 0, 0, 0, 165, 165, 165, -9, -9, 1538, -5,
                                                 -366, -894, 445, -103 };
  check(integerAGProp(all_topologies, all_top_exist, IntInfoCode::CONSTRAINT_GROUP_SUM),
        RelationalOperator::EQUAL, cnst_grp_sum_answer, "Constraint group sums did not match the "
        "expected results.", top_check);
  section(4);
  const std::vector<double> first_mol_mass = {    18.016,    18.016,    18.016,    18.016,
                                                2170.450,  2170.450,  2170.450,   157.000,
                                                 157.000, 17989.315,   328.208,  3018.992,
                                                3150.940,  8198.542,   318.450 };
  const std::vector<double> first_mol_massf = {    18.01600003,    18.01600003,    18.01600003,
                                                   18.01600003,  2170.45003104,  2170.45003104,
                                                 2170.45002127,   157.00000298,   157.00000298,
                                                17989.31466782,   328.20800292,  3018.99202657,
                                                 3150.94002604,  8198.54211676,   318.45000529 };
  check(realAGProp<double>(all_topologies, all_top_exist, RealInfoCode::FIRST_MOLECULE_MASS),
        RelationalOperator::EQUAL, Approx(first_mol_mass).margin(1.0e-7), "The mass of the first "
        "molecule from each topology was not computed as expected, in double precision.",
        top_check);
  check(realAGProp<float>(all_topologies, all_top_exist, RealInfoCode::FIRST_MOLECULE_MASS),
        RelationalOperator::EQUAL, Approx(first_mol_massf).margin(1.0e-7), "The mass of the first "
        "molecule from each topology was not computed as expected, in single precision.",
        top_check);
  const std::vector<double> net_sys_charge = {   0.00,   0.00,   0.00,   0.00,   1.00,   1.00,
                                                 1.00,   0.00,  -0.25, -11.00,  -1.00, -18.00,
                                                 0.00,  -3.00,   0.00 };
  section(3);
  check(realAGProp<float>(all_topologies, all_top_exist, RealInfoCode::TOTAL_CHARGE),
        RelationalOperator::EQUAL, Approx(net_sys_charge).margin(1.0e-7), "Total charges of each "
        "system do not meet expectations when computed in single precision.", top_check);
  check(realAGProp<double>(all_topologies, all_top_exist, RealInfoCode::TOTAL_CHARGE),
        RelationalOperator::EQUAL, Approx(net_sys_charge).margin(1.0e-7), "Total charges of each "
        "system do not meet expectations when computed in single precision.", top_check);
  const std::vector<double> max_abs_q = { 0.83400011, 1.04843998, 1.04843998, 0.24100006,
                                          0.93461990, 0.93462014, 0.94070017, 0.20714998,
                                          0.25000000, 0.80002332, 1.27849996, 1.16590011,
                                          1.16620004, 1.35828412, 0.83400023 };
  check(realAGProp<float>(all_topologies, all_top_exist, RealInfoCode::MAX_ABS_CHARGE),
        RelationalOperator::EQUAL, Approx(max_abs_q).margin(1.0e-7), "Maximum absolute charges of "
        "each system do not meet expectations when computed in single precision.", top_check);
  check(realAGProp<double>(all_topologies, all_top_exist, RealInfoCode::MAX_ABS_CHARGE),
        RelationalOperator::EQUAL, Approx(max_abs_q).margin(1.0e-7), "Maximum absolute charges of "
        "each system do not meet expectations when computed in double precision.", top_check);
  const std::vector<double> mean_ljsig = { 3.15075241, 3.16434930, 3.16434930, 3.11992731,
                                           2.81266724, 3.09622254, 2.78547205, 3.05260506,
                                           3.08295501, 2.75337292, 3.14475719, 2.89930334,
                                           3.06696950, 2.88887753, 3.13747153 };
  const std::vector<double> mean_ljsigf = { 3.15075254, 3.16434932, 3.16434932, 3.11992741,
                                            2.81266724, 3.09622266, 2.78547205, 3.05260503,
                                            3.08295510, 2.75337290, 3.14475733, 2.89930335,
                                            3.06696951, 2.88887752, 3.13747166 };
  check(realAGProp<double>(all_topologies, all_top_exist, RealInfoCode::AVERAGE_LJ_SIGMA),
        RelationalOperator::EQUAL, Approx(mean_ljsig).margin(1.0e-7), "Average (nonzero) "
        "Lennard-Jones sigma radii of particles in each system do not meet expectations when "
        "computed in double precision.", top_check);
  check(realAGProp<float>(all_topologies, all_top_exist, RealInfoCode::AVERAGE_LJ_SIGMA),
        RelationalOperator::EQUAL, Approx(mean_ljsigf).margin(3.5e-7), "Average (nonzero) "
        "Lennard-Jones sigma radii of particles in each system do not meet expectations when "
        "computed in single precision.", top_check);
  const std::vector<double> mean_ljeps = { 0.15200000, 0.16275000, 0.16275000, 0.16000000,
                                           0.07599833, 0.13974167, 0.07599833, 0.08887500,
                                           0.08425000, 0.06467212, 0.15072010, 0.09753981,
                                           0.13999725, 0.11436276, 0.14804542 };
  const std::vector<double> mean_ljepsf = { 0.15200000, 0.16275001, 0.16275001, 0.16000000,
                                            0.07599833, 0.13974166, 0.07599833, 0.08887500,
                                            0.08424999, 0.06467212, 0.15072009, 0.09753981,
                                            0.13999724, 0.11436276, 0.14804541 };
  check(realAGProp<double>(all_topologies, all_top_exist, RealInfoCode::AVERAGE_LJ_EPS),
        RelationalOperator::EQUAL, Approx(mean_ljeps).margin(1.0e-7), "Average (nonzero) "
        "Lennard-Jones epsilon depths of particles in each system do not meet expectations when "
        "computed in double precision.", top_check);
  check(realAGProp<float>(all_topologies, all_top_exist, RealInfoCode::AVERAGE_LJ_EPS),
        RelationalOperator::EQUAL, Approx(mean_ljepsf).margin(1.0e-7), "Average (nonzero) "
        "Lennard-Jones epsilon depths of particles in each system do not meet expectations when "
        "computed in single precision.", top_check);
  section(2);
  const std::vector<double> mean_bond_k = { 553.00000000, 553.00000000, 553.00000000, 553.00000000,
                                            384.80000000, 542.55070140, 384.80000000, 367.69166667,
                                            396.52500000, 342.82078874, 551.29265203, 384.82941176,
                                            526.22495816, 486.71125402, 550.07812399 };
  const std::vector<double> mean_bond_kf = { 553.00000000, 553.00000000, 553.00000000,
                                             553.00000000, 384.80000000, 542.55070140,
                                             384.80000000, 367.69166667, 396.52500000,
                                             342.82078874, 551.29265203, 384.82941176,
                                             526.22495816, 486.71125402, 550.07812399 };
  check(realAGProp<double>(all_topologies, all_top_exist, RealInfoCode::AVERAGE_BOND_STIFFNESS),
        RelationalOperator::EQUAL, Approx(mean_bond_k).margin(1.0e-7), "Average bond stiffness "
        "of each system does not meet expectations when computed in double precision.", top_check);
  check(realAGProp<double>(all_topologies, all_top_exist, RealInfoCode::AVERAGE_BOND_STIFFNESS),
        RelationalOperator::EQUAL, Approx(mean_bond_kf).margin(1.0e-7), "Average bond stiffness "
        "of each system does not meet expectations when computed in single precision.", top_check);
  const std::vector<double> mean_angl_t = { 0.00000000, 0.00000000, 0.00000000, 0.00000000,
                                            1.97325960, 1.97325960, 1.97411651, 2.09199132,
                                            2.09199132, 1.97582047, 1.98536469, 1.97810821,
                                            1.97804200, 1.96204695, 1.95809620 };
  const std::vector<double> mean_angl_tf = { 0.00000000, 0.00000000, 0.00000000, 0.00000000,
                                             1.97325962, 1.97325962, 1.97411653, 2.09199135,
                                             2.09199135, 1.97582048, 1.98536469, 1.97810823,
                                             1.97804202, 1.96204697, 1.95809620 };
  check(realAGProp<double>(all_topologies, all_top_exist, RealInfoCode::AVERAGE_ANGL_THETA),
        RelationalOperator::EQUAL, Approx(mean_angl_t).margin(1.0e-7), "Average angle equilibria "
        "of each system do not meet expectations when computed in double precision.", top_check);
  check(realAGProp<double>(all_topologies, all_top_exist, RealInfoCode::AVERAGE_ANGL_THETA),
        RelationalOperator::EQUAL, Approx(mean_angl_tf).margin(1.0e-7), "Average angle equilibria "
        "of each system do not meet expectations when computed in single precision.", top_check);
  const std::vector<double> mean_dihe_comp = {  0.00000000,  0.00000000,  0.00000000,  0.00000000,
                                                3.10533885,  3.10533885,  4.81443331, 19.60354656,
                                               19.60354656,  3.53448926,  5.44725081,  6.62202276,
                                                5.85881812,  2.73704401,  4.04773160 };
  const std::vector<double> mean_dihe_compf = {  0.00000000,  0.00000000,  0.00000000,  0.00000000,
                                                 3.10533880,  3.10533880,  4.81443322, 19.60354614,
                                                19.60354614,  3.53448930,  5.44725085,  6.62202267,
                                                 5.85881807,  2.73704395,  4.04773151 };
  check(realAGProp<double>(all_topologies, all_top_exist, RealInfoCode::AVERAGE_DIHE_COMP),
        RelationalOperator::EQUAL, Approx(mean_dihe_comp).margin(1.0e-7), "Average dihedral angle "
        "parameter composites of each system do not meet expectations when computed in double "
        "precision.", top_check);
  check(realAGProp<float>(all_topologies, all_top_exist, RealInfoCode::AVERAGE_DIHE_COMP),
        RelationalOperator::EQUAL, Approx(mean_dihe_compf).margin(1.0e-7), "Average dihedral "
        "angle parameter composites of each system do not meet expectations when computed in "
        "single precision.", top_check);
  const std::vector<double> mean_ubrd_comp = {  0.00000000,  0.00000000,  0.00000000,  0.00000000,
                                                0.00000000,  0.00000000,  0.00000000,  0.00000000,
                                                0.00000000, 45.86209669,  0.00000000,  0.00000000,
                                                0.00000000,  0.00000000,  0.00000000 };
  const std::vector<double> mean_ubrd_compf = {  0.00000000,  0.00000000,  0.00000000,  0.00000000,
                                                 0.00000000,  0.00000000,  0.00000000,  0.00000000,
                                                 0.00000000, 45.86209630,  0.00000000,  0.00000000,
                                                 0.00000000,  0.00000000,  0.00000000 };
  const std::vector<double> mean_cimp_comp = {  0.00000000,  0.00000000,  0.00000000,  0.00000000,
                                                0.00000000,  0.00000000,  0.00000000,  0.00000000,
                                                0.00000000, 61.64712919,  0.00000000,  0.00000000,
                                                0.00000000,  0.00000000,  0.00000000 };
  const std::vector<double> mean_cmap_val = {  0.00000000,  0.00000000,  0.00000000,  0.00000000,
                                               0.00000000,  0.00000000,  0.00000000,  0.00000000,
                                               0.00000000, -1.11202291,  0.00000000,  0.00000000,
                                               0.00000000,  1.49289928,  0.00000000 };
  const std::vector<double> mean_cmap_valf = {  0.00000000,  0.00000000,  0.00000000,  0.00000000,
                                                0.00000000,  0.00000000,  0.00000000,  0.00000000,
                                                0.00000000, -1.11202292,  0.00000000,  0.00000000,
                                                0.00000000,  1.49289928,  0.00000000 };
  check(realAGProp<double>(all_topologies, all_top_exist, RealInfoCode::AVERAGE_UBRD_COMP),
        RelationalOperator::EQUAL, Approx(mean_ubrd_comp).margin(1.0e-7), "Average Urey-Bradley "
        "parameter composites of each system do not meet expectations when computed in double "
        "precision.", top_check);
  check(realAGProp<float>(all_topologies, all_top_exist, RealInfoCode::AVERAGE_UBRD_COMP),
        RelationalOperator::EQUAL, Approx(mean_ubrd_compf).margin(1.0e-7), "Average Urey-Bradley "
        "parameter composites of each system do not meet expectations when computed in single "
        "precision.", top_check);
  check(realAGProp<double>(all_topologies, all_top_exist, RealInfoCode::AVERAGE_CIMP_COMP),
        RelationalOperator::EQUAL, Approx(mean_cimp_comp).margin(1.0e-7), "Average CHARMM "
        "improper dihedral parameter composites of each system do not meet expectations when "
        "computed in double precision.", top_check);
  check(realAGProp<float>(all_topologies, all_top_exist, RealInfoCode::AVERAGE_CIMP_COMP),
        RelationalOperator::EQUAL, Approx(mean_cimp_comp).margin(1.0e-7), "Average CHARMM "
        "improper dihedral parameter composites of each system do not meet expectations when "
        "computed in single precision.", top_check);
  check(realAGProp<double>(all_topologies, all_top_exist, RealInfoCode::AVERAGE_CMAP_VALUE),
        RelationalOperator::EQUAL, Approx(mean_cmap_val).margin(1.0e-7), "Average CMAP surface "
        "values of each system do not meet expectations when computed in double precision.",
        top_check);
  check(realAGProp<float>(all_topologies, all_top_exist, RealInfoCode::AVERAGE_CMAP_VALUE),
        RelationalOperator::EQUAL, Approx(mean_cmap_val).margin(1.0e-7), "Average CMAP surface "
        "values of each system do not meet expectations when computed in single precision.",
        top_check);

  // Check the ability to detect Lennard-Jones combining rules
  check(inferCombiningRule(tip3p) == VdwCombiningRule::LORENTZ_BERTHELOT, "TIP3P water should "
        "register as having a " + getEnumerationName(VdwCombiningRule::LORENTZ_BERTHELOT) +
        " Lennard-Jones combining rule.", top_check);
  check(inferCombiningRule(tip4p_error) == VdwCombiningRule::LORENTZ_BERTHELOT, "TIP4P water, "
        "even when improperly created, should register as having a " +
        getEnumerationName(VdwCombiningRule::LORENTZ_BERTHELOT) + " Lennard-Jones combining rule.",
        top_check);
  check(inferCombiningRule(trpcage_noz) == VdwCombiningRule::LORENTZ_BERTHELOT, "The ancient "
        "Amber Trp-cage topology should register as having a " +
        getEnumerationName(VdwCombiningRule::LORENTZ_BERTHELOT) + " Lennard-Jones combining rule.",
        top_check);
  const std::string base_crd_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const TestSystemManager tsm(base_top_name, "top", { "tamavidin" }, base_crd_name, "inpcrd",
                              { "tamavidin" });
  const VdwCombiningRule tama_rule = (tsm.getTestingStatus() == TestPriority::CRITICAL) ?
    inferCombiningRule(tsm.getTopologyPointer(0)) : VdwCombiningRule::NBFIX;
  check(tama_rule == VdwCombiningRule::LORENTZ_BERTHELOT, "The tamavidin topology should register "
        "as having a " + getEnumerationName(VdwCombiningRule::LORENTZ_BERTHELOT) +
        " Lennard-Jones combining rule.", tsm.getTestingStatus());

  // Create a matrix that should register as having a geometric Lennard-Jones combining rule.
  const std::vector<double> mock_sigma   = { 3.15, 2.90, 1.50, 3.00 };
  const std::vector<double> mock_epsilon = { 0.10, 0.20, 0.05, 0.30 };
  std::vector<double> mock_lja(16), mock_ljb(16);
  for (int i = 0; i < 4; i++) {
    for (int j = 0; j <= i; j++) {
      const double sig_ij = sqrt(mock_sigma[i] * mock_sigma[j]);
      const double eps_ij = sqrt(mock_epsilon[i] * mock_epsilon[j]);
      const double sig6_ij = sig_ij * sig_ij * sig_ij * sig_ij * sig_ij * sig_ij;
      const size_t ij_idx = (4 * j) + i;
      const size_t ji_idx = (4 * i) + j;
      mock_ljb[ij_idx] = 4.0 * eps_ij * sig6_ij;
      mock_lja[ij_idx] = mock_ljb[ij_idx] * sig6_ij;
      if (i != j) {
        mock_lja[ji_idx] = mock_lja[ij_idx];
        mock_ljb[ji_idx] = mock_ljb[ij_idx];
      }
    }
  }
  const VdwCombiningRule geom_test_rule = inferCombiningRule(mock_lja, mock_ljb);
  check(geom_test_rule == VdwCombiningRule::GEOMETRIC, "A matrix that should register as having "
        "a " + getEnumerationName(VdwCombiningRule::GEOMETRIC) + " Lennard-Jones combining rule "
        "instead has a " + getEnumerationName(geom_test_rule) + ".");

  // Modify the test to have NBFix characteristics.
  mock_lja[(4 * 2) + 3] += 8.1;
  mock_lja[(4 * 3) + 2] += 8.1;
  const VdwCombiningRule nbfx_test_rule = inferCombiningRule(mock_lja, mock_ljb);
  check(nbfx_test_rule == VdwCombiningRule::NBFIX, "A matrix that should register as having a " +
        getEnumerationName(VdwCombiningRule::NBFIX) + " Lennard-Jones combining rule instead has "
        "a " + getEnumerationName(nbfx_test_rule) + ".");

  // Check the correspondence in atom type names and Lennard-Jones indices
  bool type_dictionary_passes = false;
  if (tsm.getSystemCount() == 1) {
    type_dictionary_passes = true;
    const std::vector<std::vector<char4>> atyp_dict =
      tsm.getTopologyPointer(0)->getAtomTypeNameTable();
  }
  
  // Resume checking atomic details
  section(1);
  const std::vector<char4> fifth_atom_name_answer = {
    stringToChar4("H1  "), stringToChar4("O   "), stringToChar4("O   "), stringToChar4("EP2 "),
    stringToChar4("CA  "), stringToChar4("CA  "), stringToChar4("CA  "), stringToChar4("C1  "),
    stringToChar4("C1  "), stringToChar4("CA  "), stringToChar4("H1  "), stringToChar4("H5'2"),
    stringToChar4("H5''"), stringToChar4("CA  "), stringToChar4("H2  ") };
  const std::vector<char4> first_residue_name_answer = {
    stringToChar4("WAT "), stringToChar4("WAT "), stringToChar4("WAT "), stringToChar4("WAT "),
    stringToChar4("ASN "), stringToChar4("ASN "), stringToChar4("ASN "), stringToChar4("LIG "),
    stringToChar4("LIG "), stringToChar4("MET "), stringToChar4("MOL "), stringToChar4("DG5 "),
    stringToChar4("C5  "), stringToChar4("MET "), stringToChar4("LIG ") };
  const std::vector<char4> most_common_type_answer = {
    stringToChar4("HW  "), stringToChar4("HW  "), stringToChar4("HW  "), stringToChar4("HW  "),
    stringToChar4("HC  "), stringToChar4("HW  "), stringToChar4("HC  "), stringToChar4("ca  "),
    stringToChar4("EP  "), stringToChar4("HA  "), stringToChar4("HW  "), stringToChar4("CT  "),
    stringToChar4("HW  "), stringToChar4("HW  "), stringToChar4("HW  ") };
  check(charAGProp(all_topologies, all_top_exist, CharInfoCode::FIFTH_ATOM_NAME),
        RelationalOperator::EQUAL, fifth_atom_name_answer, "Atom names of topologies were not "
        "transcribed as expected.", top_check);
  check(charAGProp(all_topologies, all_top_exist, CharInfoCode::FIRST_RESIDUE_NAME),
        RelationalOperator::EQUAL, first_residue_name_answer, "Residue names of topologies were "
        "not transcribed as expected.", top_check);
  check(charAGProp(all_topologies, all_top_exist, CharInfoCode::MOST_COMMON_TYPE_NAME),
        RelationalOperator::EQUAL, most_common_type_answer, "Atom type name of topologies were "
        "not tallied as expected.", top_check);

  // Perform additional analyses on each topology
  section(4);
  const std::vector<AtomGraph*> dry_topologies = { &brbz, &brbz_vs, &trpcage, &trpcage_noz, &dhfr,
                                                   &dna };
  for (size_t i = 0; i < dry_topologies.size(); i++) {
    const bool wcmp = (all_top_exist) ?
                      (identifyWaterModel(*dry_topologies[i]) == WaterModel::NONE) : false;
    check(wcmp, "Some sort of water (" +
          getEnumerationName(identifyWaterModel(*dry_topologies[i])) + ") was identified in " +
          dry_topologies[i]->getFileName() + ", which should have no water.", top_check);
  }
  const std::vector<AtomGraph*> wet_topologies = { &tip3p, &tip4p, &tip4p_error, &tip5p,
                                                   &trpcage_water, &camp, &rna, &ubiquitin,
                                                   &drug };
  const std::vector<WaterModel> water_types = { WaterModel::TIP3P, WaterModel::TIP4P_EW,
                                                WaterModel::MULTI_CHIMERA, WaterModel::TIP5P,
                                                WaterModel::TIP3P, WaterModel::TIP3P,
                                                WaterModel::SPC_E, WaterModel::OPC,
                                                WaterModel::TIP3P };
  for (size_t i = 0; i < wet_topologies.size(); i++) {
    const bool wcmp = (all_top_exist) ?
                      (identifyWaterModel(*wet_topologies[i]) == water_types[i]) : false;
    check(wcmp, "Topology " + wet_topologies[i]->getFileName() + " was found to have " +
          getEnumerationName(identifyWaterModel(*wet_topologies[i])) + " water, but it contains " +
          getEnumerationName(water_types[i]) + ".", top_check);
  }

  // Test traps for bad input
  section(5);
  AtomGraph ntrpcage;
  CHECK_THROWS_SOFT(ntrpcage.buildFromPrmtop(trpcage_top_name, ExceptionResponse::DIE, 332.063711,
                                             1.2, 2.0, 0.001, 0.002), "A very high charge "
                    "discretization increment was permitted during charge smoothing.", top_check);
  CHECK_THROWS_SOFT(std::vector<int> bad_znums = trpcage.getAtomicNumber(-1, ntrp_atom),
                    "A nonsensical index was processed in creating an atomic numbers list.",
                    top_check);
  CHECK_THROWS_SOFT(std::vector<bool> bad_mobility = trpcage.getAtomMobility(ntrp_atom, 0),
                    "A nonsensical index was processed in creating an atomic numbers list.",
                    top_check);
  CHECK_THROWS_SOFT(std::vector<int> bad_content = trpcage.getMoleculeContents(2), "A nonsensical "
                    "index was processed in creating a molecule's atom contents list.", top_check);
  CHECK_THROWS_SOFT(float bad_q = trpcage.getPartialCharge<float>(-1), "A nonsensical index was "
                    "processed in reading an atom's partial charge.", top_check);
  CHECK_THROWS_SOFT(std::vector<double> bad_qs = trpcage.getPartialCharge<double>(-1, ntrp_atom),
                    "A nonsensical index was processed in reading atomic partial charges.",
                    top_check);
  CHECK_THROWS_SOFT(float bad_ljsig = trpcage.getLennardJonesSigma<float>(-1), "A nonsensical "
                    "index was processed in reading a Lenard-Jones type's sigma parameter.",
                    top_check);
  CHECK_THROWS_SOFT(AngleTerm<double> bad_angl = trpcage.getAngleTerm<double>(-1), "A nonsensical "
                    "index was processed in reading an angle term.", top_check);
  CHECK_THROWS_SOFT(DihedralTerm<double> bad_dihe = tip3p.getDihedralTerm<double>(2), "An attempt "
                    "was made to extract a non-existent dihedral term.", top_check);
  CHECK_THROWS_SOFT(double bad_mass = trpcage.getAtomicMass<double>(-5, MassForm::INVERSE),
                    "An invalid index was given to the atomic mass getter function.", top_check);

  // Summary evaluation
  printTestSummary(oe.getVerbosity());  
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
