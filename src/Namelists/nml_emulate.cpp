#include <algorithm>
#include "copyright.h"
#include "Chemistry/periodic_table.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "Reporting/error_format.h"
#include "namelist_element.h"
#include "namelist_emulator.h"
#include "namelist_enumerators.h"
#include "nml_emulate.h"

namespace stormm {
namespace namelist {

using chemistry::elemental_symbols;
using constants::CaseSensitivity;
using energy::default_emul_support_width;
using energy::default_emul_support_start;
using energy::default_emul_basis_zero_hold;
using energy::default_emul_restraint_weight;
using energy::default_emul_target_weight;
using energy::default_force_zero_hold;
using energy::NBAtomSource;
using energy::NBAtomSource;
using errors::rtErr;
using errors::rtWarn;
using parse::NumberFormat;
using parse::stringToChar2;
using parse::realToString;
using parse::TextOrigin;
using stmath::locateValue;
  
//-------------------------------------------------------------------------------------------------
EmulatorControls::EmulatorControls(const ExceptionResponse policy_in) :
    policy{policy_in},
    source_list{}, target_list{}, balance_list{}, mm_context{}, support_widths{}, support_starts{},
    bf_zero_holds{}, nml_transcript{"mmgbsa"}
{
  // Load in a blank namelist so that certain keywords will be present, as if this were the means
  // by which the data was loaded.
  std::string tfs("&mmgbsa\n&end\n");
  TextFile tf(tfs, TextOrigin::RAM);
  int start_line = 0;
  bool found;
  nml_transcript = emulatorInput(tf, &start_line, &found, ExceptionResponse::SILENT);
}
  
//-------------------------------------------------------------------------------------------------
EmulatorControls::EmulatorControls(const TextFile &tf, int *start_line, bool *found_nml,
                                   const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    EmulatorControls(policy_in)
{
  const NamelistEmulator t_nml = emulatorInput(tf, start_line, found_nml, policy, wrap);
  nml_transcript = t_nml;

  // Take in the list of atom sources
  const int nsrc = t_nml.getKeywordEntries("source");
  if (nsrc == 0) {
    rtErr("At least one atom source must be specfied.", "EmulatorControls");
  }
  const size_t n_elements = elemental_symbols.size();
  std::vector<bool> source_covered(nsrc, false);
  for (int i = 0; i < nsrc; i++) {
    if (source_covered[i]) {
      continue;
    }
    source_list.push_back(extractNBAtomSourceInput(t_nml, i));

    // Check subsequent source declarations for a similar label and expand the current source if
    // warranted.
    const std::string i_label = t_nml.getStringValue("source", "label", i);
    for (int j = i + 1; j < nsrc; j++) {
      if (source_covered[j]) {
        continue;
      }
      if (t_nml.getStringValue("source", "label", j) == i_label) {
        source_list.back().merge(extractNBAtomSourceInput(t_nml, j));
        source_covered[j] = true;
      }
    }
  }

  // Re-order the sources to put the most particular first and broader categories last.  Ensure
  // that, if one class of sources is completely subsumed within another, it gets precedence in
  // the list.  This bubble sort is needed because, while it may be true that A' is within A and
  // B' is within B, if B is not within A then B could be placed after A.  But, if A is compared
  // to B' first, A would not be within B' and therefore placed after B', which is within B and
  // therefore belongs before B.  In such a case, A would come after B, B' would come after A, and
  // B would come after B', which is impossible.  Only by making sequential comparisons of all
  // possible pairs can the ordering be resolved in a sensible manner.
  const size_t nrefined_src = source_list.size();
  for (int i = nrefined_src - 1; i > 0; i--) {
    for (int j = i - 1; j >= 0; j--) {
      if (source_list[i].subsumedBy(source_list[j])) {
        const NBAtomSource tmp = source_list[i];
        source_list[i] = source_list[j];
        source_list[j] = tmp;
      }
    }
  }

  // Take in directives about the basis functions for each pair potential.
  generic_zero_hold = t_nml.getRealValue("regularization");
  const int nbss = t_nml.getKeywordEntries("basis");
  std::vector<double3> tmp_basis(nbss);
  for (int i = 0; i < nbss; i++) {
    tmp_basis[i].x = t_nml.getRealValue("basis", "support", i);
    tmp_basis[i].y = t_nml.getRealValue("basis", "start", i);
    if (t_nml.getKeywordStatus("basis", "regularization", i) != InputStatus::MISSING) {
      tmp_basis[i].z = t_nml.getRealValue("basis", "regularization", i);
    }
    else {
      tmp_basis[i].z = generic_zero_hold;
    }
  }

  // Prune duplicate basis functions before committing them to the official lists.
  std::sort(tmp_basis.begin(), tmp_basis.end(),
            [](double3 a, double3 b) { if (a.x > b.x) {
                                         return true;
                                       }
                                       else if (a.x < b.x) {
                                         return false;
                                       }
                                       else {
                                         return (a.y > b.y);
                                       }
                                     });
  for (int i = 0; i < nbss; i++) {
    support_widths.push_back(tmp_basis[i].x);
    support_starts.push_back(tmp_basis[i].y);
    bf_zero_holds.push_back(tmp_basis[i].z);
    int j = i + 1;
    while (j < nbss &&
           fabs(tmp_basis[j].x - tmp_basis[i].x) < 1.0e-4 &&
           fabs(tmp_basis[j].y - tmp_basis[i].y) < 1.0e-4) {
      j++;
      i++;
    }
  }

  // Read all targets from the input.
  const int ntrg = t_nml.getKeywordEntries("target");
  for (int i = 0; i < ntrg; i++) {
    const std::string trg_label = t_nml.getStringValue("target", "label", i);
    const double trg_val = t_nml.getRealValue("target", "energy", i);
    const double trg_weight = t_nml.getRealValue("target", "weight", i);
    
    // If a structure is specified, an alternate structure may or may not be.
    if (t_nml.getKeywordStatus("target", "structure", i) == InputStatus::USER_SPECIFIED) {
      const int strc_idx = t_nml.getIntValue("target", "frame", i);
      if (t_nml.getKeywordStatus("target", "alternate", i) == InputStatus::USER_SPECIFIED) {
        const int altr_idx = t_nml.getIntValue("target", "altr_frame", i);
        target_list.emplace_back(trg_label, trg_val, trg_weight,
                                 t_nml.getStringValue("target", "structure", i),
                                 t_nml.getStringValue("target", "alternate", i), strc_idx,
                                 altr_idx);                                 
      }
      else {
        target_list.emplace_back(trg_label, trg_val, trg_weight,
                                 t_nml.getStringValue("target", "structure", i), strc_idx);
      }
    }
    else {
      if (t_nml.getKeywordStatus("target", "complex", i) != InputStatus::USER_SPECIFIED ||
          t_nml.getKeywordStatus("target", "ligand", i) != InputStatus::USER_SPECIFIED ||
          t_nml.getKeywordStatus("target", "receptor", i) != InputStatus::USER_SPECIFIED) {
        std::string lacking_msg;
        if (t_nml.getKeywordStatus("target", "complex", i) != InputStatus::USER_SPECIFIED) {
          lacking_msg += "complex";
        }
        if (t_nml.getKeywordStatus("target", "ligand", i) != InputStatus::USER_SPECIFIED) {
          if (lacking_msg.size() > 0) {
            lacking_msg += ", ";
          }
          lacking_msg += "ligand";
        }
        if (t_nml.getKeywordStatus("target", "receptor", i) != InputStatus::USER_SPECIFIED) {
          if (lacking_msg.size() > 0) {
            lacking_msg += ", ";
          }
          lacking_msg += "receptor";
        }
        rtErr("A structure, pair of structures, or a complex / receptor / ligand combination "
              "must be provided. (Missing " + lacking_msg + ".)", "EmulatorControls");
      }
      const int cmpx_idx = t_nml.getIntValue("target", "cmpx_frame", i);
      const int lgnd_idx = t_nml.getIntValue("target", "lgnd_frame", i);
      const int rcpt_idx = t_nml.getIntValue("target", "rcpt_frame", i);
      target_list.emplace_back(trg_label, trg_val, trg_weight,
                               t_nml.getStringValue("target", "complex", i),
                               t_nml.getStringValue("target", "ligand", i),
                               t_nml.getStringValue("target", "receptor", i),
                               cmpx_idx, lgnd_idx, rcpt_idx,
                               t_nml.getStringValue("target", "cmpx_mask", i),
                               t_nml.getStringValue("target", "lgnd_mask", i),
                               t_nml.getStringValue("target", "rcpt_mask", i));
    }
  }

  // Read information on the molecular mechanics context.  A default setting is implemented here,
  // to overcome the fact that boolean keywords default to FALSE.
  if (t_nml.getKeywordStatus("context") == InputStatus::MISSING) {
    mm_context = { StateVariable::BOND, StateVariable::ANGLE, StateVariable::PROPER_DIHEDRAL,
                   StateVariable::IMPROPER_DIHEDRAL, StateVariable::UREY_BRADLEY,
                   StateVariable::CHARMM_IMPROPER, StateVariable::CMAP, StateVariable::VDW,
                   StateVariable::VDW_ONE_FOUR };
  }
  else {
    if (t_nml.getBoolValue("context", "bond")) {
      mm_context.push_back(StateVariable::BOND);
    }
    if (t_nml.getBoolValue("context", "angle")) {
      mm_context.push_back(StateVariable::ANGLE);
    }
    if (t_nml.getBoolValue("context", "dihedral")) {
      mm_context.push_back(StateVariable::PROPER_DIHEDRAL);
      mm_context.push_back(StateVariable::IMPROPER_DIHEDRAL);
    }
    if (t_nml.getBoolValue("context", "urey_bradley")) {
      mm_context.push_back(StateVariable::UREY_BRADLEY);
    }
    if (t_nml.getBoolValue("context", "charmm_improper")) {
      mm_context.push_back(StateVariable::CHARMM_IMPROPER);
    }
    if (t_nml.getBoolValue("context", "cmap")) {
      mm_context.push_back(StateVariable::CMAP);
    }
    if (t_nml.getBoolValue("context", "vdw")) {
      mm_context.push_back(StateVariable::VDW);
      mm_context.push_back(StateVariable::VDW_ONE_FOUR);
    }
    if (t_nml.getBoolValue("context", "electrostatic")) {
      mm_context.push_back(StateVariable::ELECTROSTATIC);
      mm_context.push_back(StateVariable::ELEC_ONE_FOUR);
    }
    if (t_nml.getBoolValue("context", "gb")) {
      mm_context.push_back(StateVariable::GENERALIZED_BORN);
    }
  }

  // Read the force balancing directives.
  if (t_nml.getKeywordStatus("balance") != InputStatus::MISSING) {
    const int nbal = t_nml.getKeywordEntries("balance");
    for (int i = 0; i < nbal; i++) {
      balance_list.emplace_back(t_nml.getStringValue("balance", "structure", i),
                                t_nml.getIntValue("balance", "index", i),
                                t_nml.getStringValue("balance", "mask", i),
                                t_nml.getRealValue("balance", "regularization", i));
    }
  }
}

//-------------------------------------------------------------------------------------------------
int EmulatorControls::getSourceCount() const {
  return source_list.size();
}
  
//-------------------------------------------------------------------------------------------------
const NBAtomSource& EmulatorControls::getSource(const int index) const {
  if (index < 0 || index >= source_list.size()) {
    rtErr("Index " + std::to_string(index) + " is invalid for a list of " +
          std::to_string(source_list.size()) + " sources.", "EmulatorControls", "getSource");
  }
  return source_list[index];
}

//-------------------------------------------------------------------------------------------------
int EmulatorControls::getTargetCount() const {
  return target_list.size();
}

//-------------------------------------------------------------------------------------------------
const EmulationTarget& EmulatorControls::getTarget(int t_index) const {
  if (t_index < 0 || t_index >= target_list.size()) {
    rtErr("Target index " + std::to_string(t_index) + " is invalid for a list of " +
          std::to_string(target_list.size()) + " targets.", "EmulatorControls", "getTarget");
  }
  return target_list[t_index];
}

//-------------------------------------------------------------------------------------------------
int EmulatorControls::getBalanceCount() const {
  return balance_list.size();
}

//-------------------------------------------------------------------------------------------------
const ForceZero& EmulatorControls::getBalance(int b_index) const {
  if (b_index < 0 || b_index >= balance_list.size()) {
    rtErr("Atomic force balancing instruction index " + std::to_string(b_index) + " is invalid "
          "for a list of " + std::to_string(balance_list.size()) + " instructions.",
          "EmulatorControls", "getBalance");
  }
  return balance_list[b_index];
}

//-------------------------------------------------------------------------------------------------
int EmulatorControls::getBasisFunctionCount() const {
  return support_widths.size();
}

//-------------------------------------------------------------------------------------------------
double EmulatorControls::getSupportWidth(const int index) const {
  if (index < 0 || index >= support_widths.size()) {
    rtErr("Index " + std::to_string(index) + " is invalid for a series of " +
          std::to_string(support_widths.size()) + " basis functions.", "EmulatorControls",
          "getSupportWidth");
  }
  return support_widths[index];
}

//-------------------------------------------------------------------------------------------------
double EmulatorControls::getSupportStart(const int index) const {
  if (index < 0 || index >= support_starts.size()) {
    rtErr("Index " + std::to_string(index) + " is invalid for a series of " +
          std::to_string(support_starts.size()) + " basis functions.", "EmulatorControls",
          "getSupportStart");
  }
  return support_starts[index];
}

//-------------------------------------------------------------------------------------------------
double EmulatorControls::getRegularization() const {
  return generic_zero_hold;
}

//-------------------------------------------------------------------------------------------------
double EmulatorControls::getRegularization(const int index) const {
  if (index < 0 || index >= bf_zero_holds.size()) {
    rtErr("Index " + std::to_string(index) + " is invalid for a series of " +
          std::to_string(bf_zero_holds.size()) + " basis functions.", "EmulatorControls",
          "getRegulation");
  }
  return bf_zero_holds[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<StateVariable>& EmulatorControls::getMMContext() const {
  return mm_context;
}

//-------------------------------------------------------------------------------------------------
const NamelistEmulator EmulatorControls::getTranscript() const {
  return nml_transcript;
}

//-------------------------------------------------------------------------------------------------
NBAtomSource extractNBAtomSourceInput(const NamelistEmulator &t_nml, const int src_idx,
                                      const std::string &key_name) {
  const size_t n_elements = elemental_symbols.size();
  const char2 iupac_symbol = stringToChar2(t_nml.getStringValue(key_name, "atom", src_idx));
  int atomic_number = -1;
  for (size_t i = 0; i < n_elements; i++) {
    if (iupac_symbol == elemental_symbols[i]) {
      atomic_number = i;
    }
  }
  if (atomic_number == -1) {
    rtErr("Element " + t_nml.getStringValue(key_name, "atom", src_idx) + " (source index " +
          std::to_string(src_idx) + ") was not found in the IUPAC periodic table.",
          "extractNBAtomSourceInput");
  }
  const int n_bonds = t_nml.getIntValue(key_name, "bonds", src_idx);
  const bool inc_aromatic = t_nml.getBoolValue(key_name, "aromatic", src_idx);
  const bool inc_nonaromatic = t_nml.getBoolValue(key_name, "non-aromatic", src_idx);
  const bool inc_charged = t_nml.getBoolValue(key_name, "charged", src_idx);
  const bool inc_noncharged = t_nml.getBoolValue(key_name, "non-charged", src_idx);
  std::vector<bool> in_aromatic_group, has_formal_charge;
  if ((inc_aromatic && inc_nonaromatic) || (! inc_aromatic && ! inc_nonaromatic)) {
    in_aromatic_group = { true, false };
  }
  else if (inc_aromatic) {
    in_aromatic_group = { true };
  }
  else {
    in_aromatic_group = { false };
  }
  if ((inc_charged && inc_noncharged) || (! inc_charged && ! inc_noncharged)) {
    has_formal_charge = { true, false };
  }
  else if (inc_charged) {
    has_formal_charge = { true };
  }
  else {
    has_formal_charge = { false };
  }
  std::vector<int> bond_options = (n_bonds == -1) ? std::vector<int>() :
                                                    std::vector<int>(1, n_bonds);
  NBAtomSource result(t_nml.getStringValue(key_name, "label", src_idx),
                      std::vector<int>(1, atomic_number), bond_options, in_aromatic_group,
                      has_formal_charge);
  return result;
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator emulatorInput(const TextFile &tf, int *start_line, bool *found,
                               const ExceptionResponse policy, const WrapTextSearch wrap) {
  NamelistEmulator t_nml("emulator", CaseSensitivity::AUTOMATIC, policy, "Collects instructions "
                         "for fitting non-bonded pair potentials to reproduce target energies.");

  // Sources are atoms, defined by chemical characteristics
  const std::string source_help("Define a source for non-bonded energy emulation potentials.  All "
                                "sources will interact with all others in a matrix of "
                                "potentials.");
  const std::vector<std::string> source_keys_help = {
    "IUPAC atomic element symbol of an atom in the source type",
    "Indicate that the sources must include aromatic atoms.  Omitting both this and "
    "\"non-aromatic\" from the input will allow atoms to be included regardless of their "
    "participation in an aromatic ring system.  Including both terms will also allow all atoms "
    "matching other descriptors to be included.",
    "Indicate that sources must include atoms outside of aromatic ring systems.",
    "The number of bonds to other atoms that atoms in the source must have.",
    "Indicate that atoms carrying a net formal charge (including fractional formal charges, after "
    "considering resonance states of the molecule) should be counted among the sources.  To "
    "specify both \"charged\" and \"non-charged\" will include all atoms matching other "
    "descriptors, regardless of their formal charge state.  To omit both inputs will likewise "
    "remove charge states from the selection criteria.",
    "Indicate that atoms carrying zero formal charge should be counted among the sources.",
    "A label to name this source group.  Identical labels in repeated instances of the \"source\" "
    "keyword can be used to place multiple atomic elements with specific charge or bonding "
    "qualities under the same group of sources." };
  t_nml.addKeyword("source",
                   { "atom", "bonds", "aromatic", "non-aromatic", "charged", "non-charged",
                     "label" },
                   { NamelistType::STRING, NamelistType::INTEGER, NamelistType::BOOLEAN,
                     NamelistType::BOOLEAN, NamelistType::BOOLEAN, NamelistType::BOOLEAN,
                     NamelistType::STRING },
                   { std::string("C"), std::to_string(-1), std::string("false"),
                     std::string("false"), std::string("false"), std::string("false"),
                     std::string("src0") },
                   DefaultIsObligatory::NO, InputRepeats::YES, source_help, source_keys_help,
                   { KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                     KeyRequirement::OPTIONAL });

  // Basis functions are the elements of each pair potential.  The matrix of all pairs of atom
  // types has one of each basis function for every unique pair combination.
  const std::string basis_help("Define a basis function for the fitted potentials.  All are "
                               "piecewise cubic splines, set to one until the beginning of the "
                               "support, a cubic spline with zero derivatives at the end points "
                               "over the width of the support, and zero beyond the support.");
  const std::vector<std::string> basis_keys_help = {
    "Support width of the basis function, defining the range over which its value is changing",
    "The lower end point of the support region",
    "A restraint towards zero applicable to instances of this basis function in the fitted "
    "potential for any atom source pair.  If left unspecified, the general restraint scaling "
    "factor will be used."
  };
  t_nml.addKeyword("basis", { "support", "start", "regularization" },
                   { NamelistType::REAL, NamelistType::REAL, NamelistType::REAL },
                   { std::to_string(default_emul_support_width),
                     std::to_string(default_emul_support_start),
                     std::to_string(default_emul_basis_zero_hold) }, DefaultIsObligatory::NO,
                   InputRepeats::YES, basis_help, basis_keys_help,
                   { KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                     KeyRequirement::OPTIONAL });
  t_nml.addKeyword("regularization", NamelistType::REAL,
                   std::to_string(default_emul_basis_zero_hold), DefaultIsObligatory::NO,
                   InputRepeats::NO, "Strength of harmonic penalty used to tether the fitted "
                   "force constants to a particular target value.");

  // Energy targets are the values that a structure or combination of structures should yield
  // when calculated in the context of all applicable molecular mechanics terms and the fitted
  // pairwise potentials.
  const std::string target_help("Specify a target energy for a particular structure, or the "
                                "difference in energy between the structures of a complex, "
                                "ligand, and receptor.");
  const std::vector<std::string> target_keys_help = {
    "The label of this target, for bookkeeping and troubleshooting purposes.  If unassigned, a "
    "unique label will be generated.",
    "The label assigned to the system in the &files namelist.  This can only apply to target "
    "energies for individual structures.  For complexes, use the \"complex\", \"ligand\", and "
    "\"receptor\" keywords.",
    "Index of the structure within the system's label group.  Assigning all systems a different "
    "label in the &files namelist will operate with the default of zero.",
    "Label assigned to the alternate structure in the &files namelist",
    "Index of the alternate structure within its label group",
    "The label assigned to the complex structure in the &files namelist.  If this subkey is "
    "invoked, the \"ligand\" and \"receptor\" subkeys must also be included.",
    "Index of the complex structure within its label group",
    "Label assigned to the ligand structure in the &files namelist",
    "Index of the ligand structure within its label group",
    "Label assigned to the receptor structure in the &files namelist",
    "Index of the receptor structure within its label group",
    "Mask of atoms in the complex which, including all terms subsumed by the selection, should "
    "produce the taret relative energy.  If specified, this selection should correspond to masks "
    "in the ligand and receptor.",
    "Mask of atoms in the ligand which, including all terms subsumed by the selection, should "
    "produce the target relative energy",
    "Mask of atoms in the receptor which, including all terms subsumed by the selection, should "
    "produce the target relative energy",
    "Target energy for the system, or the target energy of the complex minus ligand minus "
    "receptor.",
    "The weight assigned to the target within the overall fitting problem.  Use this to emphasize "
    "some structure:energy correspondences over others.",
    "Assign an extra parameter to this target such that the fitted parameters express a relative "
    "energy in the end, relative to some arbitrary baseline for this system.  This can only be "
    "used with \"structure\" labels for a single system, not a pair of systems or the complex of "
    "ligand and receptor (which are already relative energies).  To use this when only a single "
    "structure and target uses the floating baseline will give the parameters nothing to fit."
  };
  t_nml.addKeyword("target", { "label", "structure", "frame", "alternate", "altr_index", "complex",
                               "cmpx_frame", "ligand", "lgnd_frame", "receptor", "rcpt_frame",
                               "cmpx_mask", "lgnd_mask", "rcpt_mask", "energy", "weight",
                               "float" },
                   { NamelistType::STRING, NamelistType::STRING, NamelistType::INTEGER,
                     NamelistType::STRING, NamelistType::INTEGER, NamelistType::STRING,
                     NamelistType::INTEGER, NamelistType::STRING, NamelistType::INTEGER,
                     NamelistType::STRING, NamelistType::INTEGER, NamelistType::STRING,
                     NamelistType::STRING, NamelistType::STRING, NamelistType::REAL,
                     NamelistType::REAL, NamelistType::BOOLEAN },
                   { std::string(""), std::string(""), std::to_string(0), std::string(""),
                     std::to_string(0), std::string(""), std::to_string(0), std::string(""),
                     std::to_string(0), std::string(""), std::to_string(0), std::string(":*"),
                     std::string(":*"), std::string(":*"), std::string(""),
                     std::to_string(default_emul_target_weight), std::string("") },
                   DefaultIsObligatory::NO, InputRepeats::YES, target_help, target_keys_help,
                   { KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL, KeyRequirement::REQUIRED,
                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL });

  // The molecular mechanics context is essential for the fitted pairwise potentials.  It is
  // against this backdrop that the energies emerging from the pair potentials have physical
  // meaning.  The pair potentials are a correction to the molecular mechanics context.
  const std::string mm_context_help("The molecular mechanics components of the energy function "
                                    "provide essential context for applying the fitted pairwise "
                                    "potentials.  The potentials can be thought of as a "
                                    "correction to the energy returned by the molecular mechanics "
                                    "model.  If no terms are specified, all valence interactions "
                                    "plus van-der Waals (Lennard-Jones) interactions will be "
                                    "presumed as the context.");
  const std::vector<std::string> mm_context_keys_help = {
    "Include harmonic bonds in the context",
    "Include harmonic angle terms in the context",
    "Include cosine-based dihedral terms in the context",
    "Include Urey-Bradley (harmonic stretching) angle terms",
    "Include harmonic planarity enforcement terms in the context",
    "Include CHARMM correction maps in the context",
    "Include non-bonded Lennard-Jones terms in the context",
    "Include electrostatic terms in the context",
    "Include Generalized Born energy in the context"
  };
  t_nml.addKeyword("context", { "bond", "angle", "dihedral", "urey_bradley",
                                "charmm_improper", "cmap", "vdw", "electrostatic", "gb" },
                   { NamelistType::BOOLEAN, NamelistType::BOOLEAN, NamelistType::BOOLEAN,
                     NamelistType::BOOLEAN, NamelistType::BOOLEAN, NamelistType::BOOLEAN,
                     NamelistType::BOOLEAN, NamelistType::BOOLEAN, NamelistType::BOOLEAN },
                   { std::string(""), std::string(""), std::string(""), std::string(""),
                     std::string(""), std::string(""), std::string(""), std::string(""),
                     std::string("") }, DefaultIsObligatory::NO, InputRepeats::YES,
                   mm_context_help, mm_context_keys_help,
                   { KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                     KeyRequirement::OPTIONAL });

  // The configurations of docked complexes or known biomolecular structures should remain intact
  // under the emualted potentials.  Provide a means for the user to specify that selected groups
  // of atoms be under little to no force.
  const std::string balance_help("Known binding modes and biomolecular structures should remain "
                                 "intact under the new potentials.  This keyword provides a means "
                                 "to specify groups of atoms that should remain \"in balance,\" "
                                 "with as little force as possible under the new potentials.  In "
                                 "practice, this means inserting an extra fitting equation for "
                                 "each atom which applies a harmonic restraint such that the sum "
                                 "of contributions from all fitted potentials and the broader "
                                 "molecular mechanics context is harmonically restrained towards "
                                 "zero.");
  const std::vector<std::string> balance_keys_help = {
    "Label of the structure containing atoms to balance, referencing one system from the &files "
    "namelist input",
    "Optional index of the structure within its label group.  The default negative value implies "
    "that all structures in the named label group should have balance applied to atoms in the "
    "selection mask.",
    "A mask string selecting atoms to balance",
    "The strength of regularization.  Analogous to the way in which regularizations on specific "
    "basis functions for each pair of sources scale with the number of instances of each pair, "
    "the weight scales with the number of structures in the label group which contribute to the "
    "data fitting.",
    "Relative weight of this target, used to scale the importance of keeping forces on atoms in "
    "the named structure small as opposed to forces on atoms in another structure"
  };
  t_nml.addKeyword("balance", { "structure", "index", "mask", "regularization", "weight" },
                   { NamelistType::STRING, NamelistType::INTEGER, NamelistType::STRING,
                     NamelistType::REAL, NamelistType::REAL, },
                   { std::string(""), std::to_string(0), std::string(":*"),
                     std::to_string(default_force_zero_hold),
                     std::to_string(default_emul_target_weight),
                   }, DefaultIsObligatory::NO,
                   InputRepeats::YES, balance_help, balance_keys_help,
                   { KeyRequirement::REQUIRED, KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL });

  // Read the namelist from the transcribed input
  *start_line = readNamelist(tf, &t_nml, *start_line, wrap, tf.getLineCount(), found);
  return t_nml;
}

} // namespace namelist
} // namespace stormm
