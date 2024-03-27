#include "copyright.h"
#include "Parsing/parse.h"
#include "namelist_element.h"
#include "nml_ffmorph.h"

namespace stormm {
namespace namelist {

using modeling::getEnumerationName;
using parse::stringToChar4;
using parse::char4ToString;

//-------------------------------------------------------------------------------------------------
FFMorphControls::FFMorphControls(const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    policy{policy_in},
    harmonic_bonds{}, harmonic_angles{}, cosine_dihedrals{}, urey_bradley_angles{},
    charmm_impropers{}, cmap_surfaces{}, attn14_scalings{}, charge_properties{},
    van_der_waals_properties{}, virtual_sites{},
    nml_transcript{"ffmorph"}
{}

//-------------------------------------------------------------------------------------------------
FFMorphControls::FFMorphControls(const TextFile &tf, int *start_line, bool *found_nml,
                                 const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    FFMorphControls(policy_in)
{
  NamelistEmulator t_nml = ffmorphInput(tf, start_line, found_nml, policy, wrap);
  nml_transcript = t_nml;
  
  // Load each kind of parameter edit
  const int nbond    = t_nml.getKeywordEntries("bond");
  const int nangl    = t_nml.getKeywordEntries("angle");
  const int ndihe    = t_nml.getKeywordEntries("dihedral");
  const int nubrd    = t_nml.getKeywordEntries("urey_bradley");
  const int ncimp    = t_nml.getKeywordEntries("charmm_improper");
  const int ncmap    = t_nml.getKeywordEntries("cmap");
  const int nattn_14 = t_nml.getKeywordEntries("attenuation");
  const int ncharge  = t_nml.getKeywordEntries("charge");
  const int nljparm  = t_nml.getKeywordEntries("vdw");
  const int nvsite   = t_nml.getKeywordEntries("virtual_site");
  harmonic_bonds.reserve(nbond);
  harmonic_angles.reserve(nangl);
  cosine_dihedrals.reserve(ndihe);
  urey_bradley_angles.reserve(nubrd);
  charmm_impropers.reserve(ncimp);
  cmap_surfaces.reserve(ncmap);
  attn14_scalings.reserve(nattn_14);
  charge_properties.reserve(ncharge);
  van_der_waals_properties.reserve(nljparm);
  virtual_sites.reserve(nvsite);
  const InputStatus stt_missing = InputStatus::MISSING;
  int bond_con = 0;
  for (int i = 0; i < nbond; i++) {
    const char4 ti = stringToChar4(t_nml.getStringValue("bond", "-ti", i));
    const char4 tj = stringToChar4(t_nml.getStringValue("bond", "-tj", i));
    const bool keq_provided = (t_nml.getKeywordStatus("bond", "-k", i) != stt_missing);
    const bool leq_provided = (t_nml.getKeywordStatus("bond", "-l0", i) != stt_missing);
    const double keq = (keq_provided) ? t_nml.getRealValue("bond", "-k", i) : 0.0;
    const double leq = (leq_provided) ? t_nml.getRealValue("bond", "-l0", i) : 0.0;
    if (keq_provided == false && leq_provided == false) {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Neither the equilibrium constant nor the stiffness constant are slated for "
              "modification in a harmonic bond for atom types " + char4ToString(ti) + " and " +
              char4ToString(tj) + ".", "FFMorphControls");
      case ExceptionResponse::WARN:
        rtWarn("Neither the equilibrium constant nor the stiffness constant are slated for "
               "modification in a harmonic bond for atom types " + char4ToString(ti) + " and " +
               char4ToString(tj) + ".  This modification will be skipped.", "FFMorphControls");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
    else {
      harmonic_bonds.emplace_back(ParameterKind::BOND, ti, tj);
      if (keq_provided) {
        harmonic_bonds[bond_con].setStiffness(keq);
      }
      if (leq_provided) {
        harmonic_bonds[bond_con].setEquilibrium(leq);
      }
      bond_con++;
    }
  }
  int angl_con = 0;
  for (int i = 0; i < nangl; i++) {
    const char4 ti = stringToChar4(t_nml.getStringValue("angle", "-ti", i));
    const char4 tj = stringToChar4(t_nml.getStringValue("angle", "-tj", i));
    const char4 tk = stringToChar4(t_nml.getStringValue("angle", "-tk", i));
    const bool keq_provided = (t_nml.getKeywordStatus("angle", "-k", i) != stt_missing);
    const bool teq_provided = (t_nml.getKeywordStatus("angle", "-theta0", i) != stt_missing);
    const double keq = (keq_provided) ? t_nml.getRealValue("angle", "-k", i) : 0.0;
    const double teq = (teq_provided) ? t_nml.getRealValue("angle", "-theta0", i) : 0.0;
    if (keq_provided == false && teq_provided == false) {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Neither the equilibrium constant nor the stiffness constant are slated for "
              "modification in a harmonic angle for atom types " + char4ToString(ti) + ", " +
              char4ToString(tj) + ", and " + char4ToString(tk) + ".", "FFMorphControls");
      case ExceptionResponse::WARN:
        rtWarn("Neither the equilibrium constant nor the stiffness constant are slated for "
               "modification in a harmonic angle for atom types " + char4ToString(ti) + ", " +
               char4ToString(tj) + ", and " + char4ToString(tk) + ".  This modification will be "
               "skipped.", "FFMorphControls");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
    else {
      harmonic_angles.emplace_back(ParameterKind::ANGLE, ti, tj, tk);
      if (keq_provided) {
        harmonic_angles[angl_con].setStiffness(keq);
      }
      if (teq_provided) {
        harmonic_angles[angl_con].setEquilibrium(teq);
      }
      angl_con++;
    }
  }
  int dihe_con = 0;
  for (int i = 0; i < ndihe; i++) {
    const char4 ti = stringToChar4(t_nml.getStringValue("dihedral", "-ti", i));
    const char4 tj = stringToChar4(t_nml.getStringValue("dihedral", "-tj", i));
    const char4 tk = stringToChar4(t_nml.getStringValue("dihedral", "-tk", i));
    const char4 tl = stringToChar4(t_nml.getStringValue("dihedral", "-tl", i));
    const int periodicity = t_nml.getIntValue("dihedral", "-n", i);
    const bool amp_provided = (t_nml.getKeywordStatus("dihedral", "-amp", i) != stt_missing);
    const bool phi_provided = (t_nml.getKeywordStatus("dihedral", "-phi", i) != stt_missing);
    const double amp = (amp_provided) ? t_nml.getRealValue("dihedral", "-amp", i) : 0.0;
    const double phi = (phi_provided) ? t_nml.getRealValue("dihedral", "-phi", i) : 0.0;
    if (amp_provided == false && phi_provided == false) {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Neither the amplitude nor the phase angle are slated for modification in a "
              "cosine-based dihedral angle for atom types " + char4ToString(ti) + ", " +
              char4ToString(tj) + ", " + char4ToString(tk) + ", and " + char4ToString(tl) +
              " with periodicity " + std::to_string(periodicity) + ".", "FFMorphControls");
      case ExceptionResponse::WARN:
        rtWarn("Neither the amplitude nor the phase angle are slated for modification in a "
               "cosine-based dihedral angle for atom types " + char4ToString(ti) + ", " +
               char4ToString(tj) + ", " + char4ToString(tk) + ", and " + char4ToString(tl) +
               " with periodicity " + std::to_string(periodicity) + ".  This modification will be "
               "skipped.", "FFMorphControls");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
    else {
      cosine_dihedrals.emplace_back(ParameterKind::DIHEDRAL, ti, tj, tk, tl);
      if (amp_provided) {
        cosine_dihedrals[dihe_con].setAmplitude(amp);
      }
      if (phi_provided) {
        cosine_dihedrals[dihe_con].setPhaseAngle(phi);
      }
      dihe_con++;
    }    
  }
  int ubrd_con = 0;
  for (int i = 0; i < nubrd; i++) {
    const char4 ti = stringToChar4(t_nml.getStringValue("urey_bradley", "-ti", i));
    const char4 tj = stringToChar4(t_nml.getStringValue("urey_bradley", "-tj", i));
    const char4 tk = stringToChar4(t_nml.getStringValue("urey_bradley", "-tk", i));
    const bool keq_provided = (t_nml.getKeywordStatus("urey_bradley", "-k", i) != stt_missing);
    const bool leq_provided = (t_nml.getKeywordStatus("urey_bradley", "-l0", i) != stt_missing);
    const double keq = (keq_provided) ? t_nml.getRealValue("urey_bradley", "-k", i) : 0.0;
    const double leq = (leq_provided) ? t_nml.getRealValue("urey_bradley", "-l0", i) : 0.0;
    if (keq_provided == false && leq_provided == false) {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Neither the stiffness constant nor the equilibrium separation were specified in a "
              "Urey-Bradley modification involving atom types " + char4ToString(ti) + ", " +
              char4ToString(tj) + ", and " + char4ToString(tk) + ".", "FFMorphControls");
      case ExceptionResponse::WARN:
        rtWarn("Neither the stiffness constant nor the equilibrium separation were specified in a "
               "Urey-Bradley modification involving atom types " + char4ToString(ti) + ", " +
               char4ToString(tj) + ", and " + char4ToString(tk) + ".  This modification will be "
               "skipped.", "FFMorphControls");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
    else {
      urey_bradley_angles.emplace_back(ParameterKind::UREY_BRADLEY, ti, tj, tk);
      if (keq_provided) {
        urey_bradley_angles[ubrd_con].setStiffness(keq);
      }
      if (leq_provided) {
        urey_bradley_angles[ubrd_con].setEquilibrium(leq);
      }
      ubrd_con++;
    }
  }
  int cimp_con = 0;
  for (int i = 0; i < ncimp; i++) {
    const char4 ti = stringToChar4(t_nml.getStringValue("charmm_improper", "-ti", i));
    const char4 tj = stringToChar4(t_nml.getStringValue("charmm_improper", "-tj", i));
    const char4 tk = stringToChar4(t_nml.getStringValue("charmm_improper", "-tk", i));
    const char4 tl = stringToChar4(t_nml.getStringValue("charmm_improper", "-tl", i));
    const bool keq_provided = (t_nml.getKeywordStatus("charmm_improper", "-k", i) != stt_missing);
    const bool phi_provided = (t_nml.getKeywordStatus("charmm_improper", "-phi", i) !=
                               stt_missing);
    const double keq = (keq_provided) ? t_nml.getRealValue("charmm_improper", "-k", i) : 0.0;
    const double phi = (phi_provided) ? t_nml.getRealValue("charmm_improper", "-phi", i) : 0.0;
    if (keq_provided == false && phi_provided == false) {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Neither the stiffness constant nor the phase angle were specified in a "
              "CHARMM improper modification involving atom types " + char4ToString(ti) + ", " +
              char4ToString(tj) + ", " + char4ToString(tk) + ", and " + char4ToString(tl) + ".",
              "FFMorphControls");
      case ExceptionResponse::WARN:
        rtWarn("Neither the stiffness constant nor the phase angle were specified in a "
               "CHARMM improper modification involving atom types " + char4ToString(ti) + ", " +
               char4ToString(tj) + ", " + char4ToString(tk) + ", and " + char4ToString(tl) + ".  "
               "This modification will be skipped.", "FFMorphControls");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
    else {
      charmm_impropers.emplace_back(ParameterKind::CHARMM_IMPROPER, ti, tj, tk, tl);
      if (keq_provided) {
        charmm_impropers[cimp_con].setStiffness(keq);
      }
      if (phi_provided) {
        charmm_impropers[cimp_con].setPhaseAngle(phi);
      }
      cimp_con++;
    }
  }
  int cmap_con = 0;
  for (int i = 0; i < ncmap; i++) {
    const char4 atomi = stringToChar4(t_nml.getStringValue("cmap", "-atom_i", i));
    const char4 atomj = stringToChar4(t_nml.getStringValue("cmap", "-atom_j", i));
    const char4 atomk = stringToChar4(t_nml.getStringValue("cmap", "-atom_k", i));
    const char4 atoml = stringToChar4(t_nml.getStringValue("cmap", "-atom_l", i));
    const char4 atomm = stringToChar4(t_nml.getStringValue("cmap", "-atom_m", i));
    const char4 resii = stringToChar4(t_nml.getStringValue("cmap", "-atom_i", i));
    const char4 resij = stringToChar4(t_nml.getStringValue("cmap", "-atom_j", i));
    const char4 resik = stringToChar4(t_nml.getStringValue("cmap", "-atom_k", i));
    const char4 resil = stringToChar4(t_nml.getStringValue("cmap", "-atom_l", i));
    const char4 resim = stringToChar4(t_nml.getStringValue("cmap", "-atom_m", i));
    const std::vector<std::string> roman_numerals = { "i", "ii", "iii", "iv" };
    std::vector<double> potential_values;
    std::vector<int2> locations;
    for (size_t j = 0; j < roman_numerals.size(); j++) {
      if (t_nml.getKeywordStatus("cmap", "-u_" + roman_numerals[j], i) != stt_missing) {
        if (t_nml.getKeywordStatus("cmap", "-phi_" + roman_numerals[j], i) == stt_missing ||
            t_nml.getKeywordStatus("cmap", "-psi_" + roman_numerals[j], i) == stt_missing) {
          switch (policy) {
          case ExceptionResponse::DIE:
            rtErr("A potential value was provided for point " + roman_numerals[j] + " but the "
                  "grid indexing is not present in a CMAP modification for atoms " +
                  char4ToString(resii) + ":" + char4ToString(atomi) + ", " +
                  char4ToString(resij) + ":" + char4ToString(atomj) + ", " +
                  char4ToString(resik) + ":" + char4ToString(atomk) + ", " +
                  char4ToString(resil) + ":" + char4ToString(atoml) + ", and " +
                  char4ToString(resim) + ":" + char4ToString(atomm) + ".", "FFMorphControls");
          case ExceptionResponse::WARN:
            rtErr("A potential value was provided for point " + roman_numerals[j] + " but the "
                  "grid indexing is not present in a CMAP modification for atoms " +
                  char4ToString(resii) + ":" + char4ToString(atomi) + ", " +
                  char4ToString(resij) + ":" + char4ToString(atomj) + ", " +
                  char4ToString(resik) + ":" + char4ToString(atomk) + ", " +
                  char4ToString(resil) + ":" + char4ToString(atoml) + ", and " +
                  char4ToString(resim) + ":" + char4ToString(atomm) + ".  This modification will "
                  "be skipped.", "FFMorphControls");
            break;
          case ExceptionResponse::SILENT:
            break;
          }
        }
        else {
          potential_values.push_back(t_nml.getRealValue("cmap", "-u_" + roman_numerals[j]));
          
          locations.push_back({ t_nml.getIntValue("cmap", "-phi_" + roman_numerals[j]),
                                t_nml.getIntValue("cmap", "-psi_" + roman_numerals[j]) });
        }
      }
    }
    if (potential_values.size() == 0LLU) {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("No potential surface modifications were specified in a CMAP modification involving "
              "atoms " + char4ToString(resii) + ":" + char4ToString(atomi) + ", " +
              char4ToString(resij) + ":" + char4ToString(atomj) + ", " + char4ToString(resik) +
              ":" + char4ToString(atomk) + ", " + char4ToString(resil) + ":" +
              char4ToString(atoml) + ", and " + char4ToString(resim) + char4ToString(atomm) + ".",
              "FFMorphControls");
      case ExceptionResponse::WARN:
        rtErr("No potential surface modifications were specified in a CMAP modification involving "
              "atoms " + char4ToString(resii) + ":" + char4ToString(atomi) + ", " +
              char4ToString(resij) + ":" + char4ToString(atomj) + ", " + char4ToString(resik) +
              ":" + char4ToString(atomk) + ", " + char4ToString(resil) + ":" +
              char4ToString(atoml) + ", and " + char4ToString(resim) + char4ToString(atomm) + ".  "
              "This modification will be skipped.", "FFMorphControls");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
    else {
      cmap_surfaces.emplace_back(ParameterKind::CMAP, atomi, atomj, atomk, atoml, atomm, resii,
                                 resij, resik, resil, resim, potential_values, locations);
      cmap_con++;
    }
  }
  int attn_con = 0;
  for (int i = 0; i < nattn_14; i++) {
    const char4 ti = stringToChar4(t_nml.getStringValue("attenuation", "-ti", i));
    const char4 tj = stringToChar4(t_nml.getStringValue("attenuation", "-tj", i));
    const char4 tk = stringToChar4(t_nml.getStringValue("attenuation", "-tk", i));
    const char4 tl = stringToChar4(t_nml.getStringValue("attenuation", "-tl", i));
    const bool qq_provided = (t_nml.getKeywordStatus("attenuation", "-qq", i) != stt_missing);
    const bool vdw_provided = (t_nml.getKeywordStatus("attenuation", "-vdw", i) != stt_missing);
    const double qq = (qq_provided) ? t_nml.getRealValue("attenuation", "-qq", i) : 0.0;
    const double vdw = (vdw_provided) ? t_nml.getRealValue("attenuation", "-vdw", i) : 0.0;
    if (qq_provided == false && vdw_provided == false) {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Neither the electrostatic nor van-der Waals scaling factors were specified in an "
              "attenuated non-bonded interaction modification involving atom types " +
              char4ToString(ti) + ", " + char4ToString(tj) + ", " + char4ToString(tk) + ", and " +
              char4ToString(tl) + ".", "FFMorphControls");
      case ExceptionResponse::WARN:
        rtWarn("Neither the electrostatic nor van-der Waals scaling factors were specified in an "
               "attenuated non-bonded interaction modification involving atom types " +
               char4ToString(ti) + ", " + char4ToString(tj) + ", " + char4ToString(tk) + ", and " +
               char4ToString(tl) + ".  This modification will be skipped.", "FFMorphControls");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
    else {
      attn14_scalings.emplace_back(ParameterKind::ATTN_14_SCALE, ti, tj, tk, tl);
      if (qq_provided) {
        attn14_scalings[attn_con].setChargeScaling(qq);
      }
      if (vdw_provided) {
        attn14_scalings[attn_con].setVanDerWaalsScaling(vdw);
      }
      attn_con++;
    }
  }
  int chrg_con = 0;
  for (int i = 0; i < ncharge; i++) {
    const char4 atomi = stringToChar4(t_nml.getStringValue("attenuation", "-atom", i));
    const char4 resii = stringToChar4(t_nml.getStringValue("attenuation", "-resi", i));
    const bool q_provided = (t_nml.getKeywordStatus("charge", "-q", i) != stt_missing);
    
  }
}

//-------------------------------------------------------------------------------------------------
int FFMorphControls::getEditCount(const ParameterKind kind) const {
  switch (kind) {
  case ParameterKind::BOND:
    return harmonic_bonds.size();
  case ParameterKind::ANGLE:
    return harmonic_angles.size();
  case ParameterKind::DIHEDRAL:
    return cosine_dihedrals.size();
  case ParameterKind::UREY_BRADLEY:
    return urey_bradley_angles.size();
  case ParameterKind::CHARMM_IMPROPER:
    return charmm_impropers.size();
  case ParameterKind::CMAP:
    return cmap_surfaces.size();
  case ParameterKind::ATTN_14_SCALE:
    return attn14_scalings.size();
  case ParameterKind::CHARGE:
    return charge_properties.size();
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
    return van_der_waals_properties.size();
  case ParameterKind::VIRTUAL_SITE_FRAME:
    return virtual_sites.size();
  case ParameterKind::NONE:
    rtErr("A valid parameter kind must be specified.", "FFMorphControls", "getEditCount");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement FFMorphControls::getModelEdit(const ParameterKind kind, const int index) const {
  const int nedits = getEditCount(kind);
  if (index >= nedits) {
    rtErr("Index " + std::to_string(index) + " was requested from an array of " +
          std::to_string(nedits) + " stated edits for " + getEnumerationName(kind) +
          " parameters.", "FFMorphControls", "getModelEdit");
  }
  switch (kind) {
  case ParameterKind::BOND:
    return harmonic_bonds[index];
  case ParameterKind::ANGLE:
    return harmonic_angles[index];
  case ParameterKind::DIHEDRAL:
    return cosine_dihedrals[index];
  case ParameterKind::UREY_BRADLEY:
    return urey_bradley_angles[index];
  case ParameterKind::CHARMM_IMPROPER:
    return charmm_impropers[index];
  case ParameterKind::CMAP:
    return cmap_surfaces[index];
  case ParameterKind::ATTN_14_SCALE:
    return attn14_scalings[index];
  case ParameterKind::CHARGE:
    return charge_properties[index];
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
    return van_der_waals_properties[index];
  case ParameterKind::VIRTUAL_SITE_FRAME:
    return virtual_sites[index];
  case ParameterKind::NONE:

    // This case is trapped in the call to getEditCount() above.
    break;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
const NamelistEmulator& FFMorphControls::getTranscript() const {
  return nml_transcript;
}
  
//-------------------------------------------------------------------------------------------------
NamelistEmulator ffmorphInput(const TextFile &tf, int *start_line, bool *found,
                              const ExceptionResponse policy, const WrapTextSearch wrap) {
  NamelistEmulator t_nml("ffmorph", CaseSensitivity::AUTOMATIC, policy, "Permits user control of "
                         "specific parameters within the topologies at hand.  The control extends "
                         "as far as changing individual parameters' values, but not the atoms to "
                         "which the parameters apply or the insertion of new potential terms.  "
                         "This namelist drives a metamorphosis of the available topologies, but "
                         "not a change in the non-bonded pair list or valence parameter array "
                         "sizes.  In all contextes, unspecified parameters will leave any "
                         "existing settings unchanged.  The modified topologies are not written "
                         "back to disk.");
  const std::string bond_help("Modify the parameters of a bond between two atom types.");
  const std::vector<std::string> bond_keys_help = {
    "Atom type for the I atom in the bond (required)",
    "Atom type for the J atom in the bond (required)",
    "Bond stiffness constant, units of kcal/mol-Angstrom",
    "Bond equilibrium constant, units of Angstroms" };
  t_nml.addKeyword(NamelistElement("bond", { "-ti", "-tj", "-k", "-l0" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::REAL, NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string("") }, DefaultIsObligatory::NO, InputRepeats::YES,
                                   bond_help, bond_keys_help,
                                   { KeyRequirement::REQUIRED, KeyRequirement::REQUIRED,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL }));
  const std::string angle_help("Modify the parameters of an angle between three atom types.");
  const std::vector<std::string> angle_keys_help = {
    "Atom type for the I atom in the angle (required)",
    "Atom type for the J atom in the angle (required)",
    "Atom type for the K atom in the angle (required)",
    "Angle stiffness constant, units of kcal/mol-radian",
    "Angle equilibrium constant, units of degrees" };
  t_nml.addKeyword(NamelistElement("angle", { "-ti", "-tj", "-tk", "-k", "-theta0" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::REAL,
                                     NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string("") }, DefaultIsObligatory::NO,
                                   InputRepeats::YES, angle_help, angle_keys_help,
                                   { KeyRequirement::REQUIRED, KeyRequirement::REQUIRED,
                                     KeyRequirement::REQUIRED, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL }));
  const std::string dihedral_help("Modify the parameters of a dihedral interaction between four "
                                  "atom types.");
  const std::vector<std::string> dihedral_keys_help = {
    "Atom type for the I atom in the dihedral angle (required)",
    "Atom type for the J atom in the dihedral angle (required)",
    "Atom type for the K atom in the dihedral angle (required)",
    "Atom type for the L atom in the dihedral angle (required)",
    "Cosine function periodicity", "Dihedral cosine function amplitude (kcal/mol)",
    "Dihedral phase angle (degrees)", "Indicator of whether the parameter set pertains to a "
    "proper or improper torsion" };
  t_nml.addKeyword(NamelistElement("dihedral",
                                   { "-ti", "-tj", "-tk", "-tl", "-n", "-amp", "-phi", "-kind" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::INTEGER, NamelistType::REAL,
                                     NamelistType::REAL, NamelistType::STRING },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string("PROPER") },
                                   DefaultIsObligatory::NO, InputRepeats::YES, dihedral_help,
                                   dihedral_keys_help,
                                   { KeyRequirement::REQUIRED, KeyRequirement::REQUIRED,
                                     KeyRequirement::REQUIRED, KeyRequirement::REQUIRED,
                                     KeyRequirement::REQUIRED, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL, KeyRequirement::REQUIRED }));
  const std::string ubrd_help("Modify the parameters of a Urey-Bradley interaction between three "
                              "atom types (the spring constant applies only between the first and "
                              "third atom types, but the central atom type is essential for "
                              "determining where the parameter should be applied).");
  const std::vector<std::string> ubrd_keys_help = {
    "Atom type for the I atom in the Urey-Bradley interaction (required)",
    "Atom type for the J atom in the Urey-Bradley interaction (required)",
    "Atom type for the K atom in the Urey-Bradley interaction (required)",
    "Urey-Bradley stretching constant (kcal/mol-Angstrom)",
    "Urey-Bradley equilibrium constant (Angstrom)" };
  t_nml.addKeyword(NamelistElement("urey_bradley", { "-ti", "-tj", "-tk", "-k", "-l0" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::REAL,
                                     NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string("") }, DefaultIsObligatory::NO,
                                   InputRepeats::YES, ubrd_help, ubrd_keys_help,
                                   { KeyRequirement::REQUIRED, KeyRequirement::REQUIRED,
                                     KeyRequirement::REQUIRED, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL }));
  const std::string cimp_help("Modify the parameters of a CHARMM improper dihedral between four "
                              "atom types.");
  const std::vector<std::string> cimp_keys_help = {
    "Atom type for the I atom in the CHARMM improper dihedral (required)",
    "Atom type for the J atom in the CHARMM improper dihedral (required)",
    "Atom type for the K atom in the CHARMM improper dihedral (required)",
    "Atom type for the L atom in the CHARMM improper dihedral (required)",
    "Harmonic stiffness of the restoring force (kcal/mol-radian)",
    "Phase angle for the equilibrium of the angle between two planes (degrees)" };
  t_nml.addKeyword(NamelistElement("charmm_improper", { "-ti", "-tj", "-tk", "-tl", "-k", "-phi" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::REAL, NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string("") },
                                   DefaultIsObligatory::NO, InputRepeats::YES, cimp_help,
                                   cimp_keys_help,
                                   { KeyRequirement::REQUIRED, KeyRequirement::REQUIRED,
                                     KeyRequirement::REQUIRED, KeyRequirement::REQUIRED,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL }));
  const std::string cmap_help("Modify selected elements of a CMAP surface.");
  const std::vector<std::string> cmap_keys_help = {
    "Name of the I atom in the CMAP interaction (required)",
    "Name of the J atom in the CMAP interaction (required)",
    "Name of the K atom in the CMAP interaction (required)",
    "Name of the L atom in the CMAP interaction (required)",
    "Name of the M atom in the CMAP interaction (required)",
    "Grid phi index (first dimension) of modified point I",
    "Grid phi index (first dimension) of modified point II",
    "Grid phi index (first dimension) of modified point III",
    "Grid phi index (first dimension) of modified point IV",
    "Grid psi index (second dimension) of modified point I",
    "Grid psi index (second dimension) of modified point II",
    "Grid psi index (second dimension) of modified point III",
    "Grid psi index (second dimension) of modified point IV",
    "Modified CMAP value at point I", "Modified CMAP value at point II",
    "Modified CMAP value at point III", "Modified CMAP value at point IV" };
  t_nml.addKeyword(NamelistElement("cmap", { "-atom_i", "-atom_j", "-atom_k", "-atom_l", "-atom_m",
                                             "-resi_i", "-resi_j", "-resi_k", "-resi_l", "-resi_m",
                                             "-phi_i", "-phi_ii", "-phi_iii", "-phi_iv", "-psi_i",
                                             "-psi_ii", "-psi_iii", "-psi_iv", "-u_i", "-u_ii",
                                             "-u_iii", "-u_iv" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::INTEGER, NamelistType::INTEGER,
                                     NamelistType::INTEGER, NamelistType::INTEGER,
                                     NamelistType::INTEGER, NamelistType::INTEGER,
                                     NamelistType::INTEGER, NamelistType::INTEGER,
                                     NamelistType::REAL, NamelistType::REAL,
                                     NamelistType::REAL, NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string(""),
                                     std::string("") }, DefaultIsObligatory::NO, InputRepeats::YES,
                                   cmap_help, cmap_keys_help,
                                   { KeyRequirement::REQUIRED, KeyRequirement::REQUIRED,
                                     KeyRequirement::REQUIRED, KeyRequirement::REQUIRED,
                                     KeyRequirement::REQUIRED, KeyRequirement::REQUIRED,
                                     KeyRequirement::REQUIRED, KeyRequirement::REQUIRED,
                                     KeyRequirement::REQUIRED, KeyRequirement::REQUIRED,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL }));
  const std::string attn14_help("Modify the non-bonded scaling parameters of an attenuated 1:4 "
                                "interaction.");
  const std::vector<std::string> attn14_keys_help = {
    "Atom type for the I atom in the attenuated 1:4 interaction (required)",
    "Atom type for the J atom in the attenuated 1:4 interaction (required)",
    "Atom type for the K atom in the attenuated 1:4 interaction (required)",
    "Atom type for the L atom in the attenuated 1:4 interaction (required)",
    "Electrostatic non-bonded inverse scaling factor (a factor of two sets electrostatic "
    "interactions to occur at half strength", "Inverse scaling factor for van-der Waals "
    "interactions" };
  t_nml.addKeyword(NamelistElement("attenuation", { "-ti", "-tj", "-tk", "-tl", "-qq", "-vdw" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::REAL, NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string("") },
                                   DefaultIsObligatory::NO, InputRepeats::YES, attn14_help,
                                   attn14_keys_help,
                                   { KeyRequirement::REQUIRED, KeyRequirement::REQUIRED,
                                     KeyRequirement::REQUIRED, KeyRequirement::REQUIRED,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL }));
  const std::string charge_help("Modify the properties of a charged atom.");
  const std::vector<std::string> charge_keys_help = {
    "Name of the atom to alter (this is not an atom type, and is required)", "Residue name of the "
    "atom to alter (required)", "Charge value, atomic units (e, the proton charge)" };
  t_nml.addKeyword(NamelistElement("charge", { "-atom", "-resi", "-q" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string("") },
                                   DefaultIsObligatory::NO, InputRepeats::YES, charge_help,
                                   charge_keys_help,
                                   { KeyRequirement::REQUIRED, KeyRequirement::REQUIRED,
                                     KeyRequirement::OPTIONAL }));
  const std::string vdw_help("Modify the van-der Waals properties of an atom.");
  const std::vector<std::string> vdw_keys_help = {
    "Atom type to alter (the van-der Waals atom types are one and the same with those used in "
    "valence parameters)", "Sigma (Lennard-Jones radius) value (Angstrom)",
    "Epsilon (Lennard-Jones well depth) value, units of kcal/mol", "Rho (tertiary parameter for "
    "some Lennard-Jones models or Buckingham potentials--the definition will depend on the "
    "context of each topology)"};
  t_nml.addKeyword(NamelistElement("vdw", { "-t", "-sig", "-eps", "-rho" },
                                   { NamelistType::STRING, NamelistType::REAL,
                                     NamelistType::REAL, NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string("") },
                                   DefaultIsObligatory::NO, InputRepeats::YES, vdw_help,
                                   vdw_keys_help,
                                   { KeyRequirement::REQUIRED, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL }));
  const std::string vsite_help("Modify the dimensions of a virtual site frame (the non-bonded "
                               "properties of the virtual site can be altered with the 'charge' "
                               "and 'vdw' keywords--this keyword pertains to the geometry of the "
                               "virtual site amidst its frame atoms)");
  const std::vector<std::string> vsite_keys_help = {
    "Virtual site atom name", "Parent atom name.  This is the real atom to which the virtual site "
    "likely transfers most of the forces it accumulates, although in principle the distinction "
    "between the parent atom and other frame atoms is the way they factor into various chain "
    "rules.  This and other atom names are case sensitive.", "Frame atom 2 name", "Frame atom 3 "
    "name (if applicable to the frame type)", "Frame atom 4 name (if applicable to the frame "
    "type)", "The virtual site frame type, i.e. Flex-2, case-insensitive)", "First frame "
    "dimension (the meaning and units of this and other dimensions depend on the context of the "
    "frame type)", "Second frame dimension", "Third frame dimension" };
  t_nml.addKeyword(NamelistElement("virtual_site", { "-vsatom", "-parent", "-frame2", "-frame3",
                                                     "-frame4", "-ft", "-dim1", "-dim2", "-dim3" },
                                   { NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::STRING, NamelistType::STRING,
                                     NamelistType::REAL, NamelistType::REAL,
                                     NamelistType::REAL },
                                   { std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string(""),
                                     std::string(""), std::string(""), std::string("") },
                                   DefaultIsObligatory::NO, InputRepeats::YES, vsite_help,
                                   vsite_keys_help,
                                   { KeyRequirement::REQUIRED, KeyRequirement::REQUIRED,
                                     KeyRequirement::REQUIRED, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL, KeyRequirement::REQUIRED,
                                     KeyRequirement::OPTIONAL, KeyRequirement::OPTIONAL,
                                     KeyRequirement::OPTIONAL }));

  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, wrap, tf.getLineCount(), found);
  return t_nml;
}
  
} // namespace namelist 
} // namespace stormm
