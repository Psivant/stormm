#include <cmath>
#include "copyright.h"
#include "Constants/symbol_values.h"
#include "Parsing/parse.h"
#include "forcefield_element.h"

namespace stormm {
namespace modeling {

using parse::uppercase;
  
//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in) :
    kind{kind_in}, atom_name_i{' ', ' ', ' ', ' '}, atom_name_j{' ', ' ', ' ', ' '},
    atom_name_k{' ', ' ', ' ', ' '}, atom_name_l{' ', ' ', ' ', ' '},
    atom_name_m{' ', ' ', ' ', ' '}, residue_name_i{' ', ' ', ' ', ' '},
    residue_name_j{' ', ' ', ' ', ' '}, residue_name_k{' ', ' ', ' ', ' '},
    residue_name_l{' ', ' ', ' ', ' '}, residue_name_m{' ', ' ', ' ', ' '}, property_a{0.0},
    property_b{0.0}, property_c{0.0}, activate_a{false}, activate_b{false}, activate_c{false},
    surface{}, locations{}, frame_type{VirtualSiteKind::NONE}
{}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in, const char4 atom_i_in) :
    ForceFieldElement(kind_in)
{
  switch (kind_in) {
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::NONE:
    break;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::VIRTUAL_SITE_FRAME:
    rtErr("Construction with only one atom type is not acceptable for a parameter of type \"" +
          getEnumerationName(kind_in) + "\".", "ForceFieldElement");
  }
  atom_name_i = atom_i_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in, const char4 atom_i_in,
                                     const char4 atom_j_in) :
    ForceFieldElement(kind_in)
{
  switch (kind_in) {
  case ParameterKind::BOND:
    break;
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("Construction with two atom types is not acceptable for a parameter of type \"" +
          getEnumerationName(kind_in) + "\".", "ForceFieldElement");
  }
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const char4 atom_k_in) :
    ForceFieldElement(kind_in)
{
  switch (kind_in) {
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
    break;
  case ParameterKind::BOND:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("Construction with three atom types is not acceptable for a parameter of type \"" +
          getEnumerationName(kind_in) + "\".", "ForceFieldElement");
  }
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  atom_name_k = atom_k_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const char4 atom_k_in,
                                     const char4 atom_l_in, const TorsionKind tkind_in) :
    ForceFieldElement(kind_in)
{
  switch (kind_in) {
  case ParameterKind::DIHEDRAL:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::ATTN_14_SCALE:
    break;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CMAP:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("Construction with four atom types is not acceptable for a parameter of type \"" +
          getEnumerationName(kind_in) + "\".", "ForceFieldElement");
  }
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  atom_name_k = atom_k_in;
  atom_name_l = atom_l_in;
  torsion_kind = tkind_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in,
                                     const VirtualSiteKind frame_type_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const char4 atom_k_in,
                                     const char4 residue_i_in, const char4 residue_j_in,
                                     const char4 residue_k_in) :
    ForceFieldElement(kind_in)
{
  switch (kind_in) {
  case ParameterKind::VIRTUAL_SITE_FRAME:
    switch (frame_type_in) {
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
      break;
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
    case VirtualSiteKind::FIXED_4:
    case VirtualSiteKind::NONE:
      rtErr("Construction with three atom and residue names does not properly specify a virtual "
            "site of frame type " + getEnumerationName(frame_type_in)+ ".",
            "ForceFieldElement");
    }
    break;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::NONE:
    rtErr("Construction with three atom and residue names is not acceptable for a parameter of "
          "type \"" + getEnumerationName(kind_in) + "\".", "ForceFieldElement");
  }
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  atom_name_k = atom_k_in;
  residue_name_i = residue_i_in;
  residue_name_j = residue_j_in;
  residue_name_k = residue_k_in;
  frame_type = frame_type_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in,
                                     const VirtualSiteKind frame_type_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const char4 atom_k_in,
                                     const char4 atom_l_in, const char4 residue_i_in,
                                     const char4 residue_j_in, const char4 residue_k_in,
                                     const char4 residue_l_in) :
  ForceFieldElement(kind_in)
{
  switch (kind_in) {
  case ParameterKind::VIRTUAL_SITE_FRAME:
    switch (frame_type_in) {
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      break;
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
    case VirtualSiteKind::FIXED_4:
    case VirtualSiteKind::NONE:
      rtErr("Construction with three atom and residue names does not properly specify a virtual "
            "site of frame type " + getEnumerationName(frame_type_in)+ ".",
            "ForceFieldElement");
    }
    break;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::NONE:
    rtErr("Construction with four atom and residue names is not acceptable for a parameter of "
          "type \"" + getEnumerationName(kind_in) + "\".", "ForceFieldElement");
  }
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  atom_name_k = atom_k_in;
  atom_name_l = atom_l_in;
  residue_name_i = residue_i_in;
  residue_name_j = residue_j_in;
  residue_name_k = residue_k_in;
  residue_name_l = residue_l_in;
  frame_type = frame_type_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in,
                                     const VirtualSiteKind frame_type_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const char4 atom_k_in,
                                     const char4 atom_l_in, const char4 atom_m_in,
                                     const char4 residue_i_in, const char4 residue_j_in,
                                     const char4 residue_k_in, const char4 residue_l_in,
                                     const char4 residue_m_in) :
  ForceFieldElement(kind_in)
{
  switch (kind_in) {
  case ParameterKind::VIRTUAL_SITE_FRAME:
    switch (frame_type_in) {
    case VirtualSiteKind::FLEX_3:
    case VirtualSiteKind::FIXED_3:
    case VirtualSiteKind::FAD_3:
    case VirtualSiteKind::OUT_3:
      break;
    case VirtualSiteKind::FLEX_2:
    case VirtualSiteKind::FIXED_2:
    case VirtualSiteKind::FIXED_4:
    case VirtualSiteKind::NONE:
      rtErr("Construction with five atom and residue names does not properly specify a virtual "
            "site of frame type " + getEnumerationName(frame_type_in)+ ".",
            "ForceFieldElement");
    }
    break;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::NONE:
    rtErr("Construction with five atom and residue names plus a virtual site frame type is not "
          "acceptable for a parameter of type \"" + getEnumerationName(kind_in) + "\".",
          "ForceFieldElement");
  }
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  atom_name_k = atom_k_in;
  atom_name_l = atom_l_in;
  atom_name_m = atom_m_in;
  residue_name_i = residue_i_in;
  residue_name_j = residue_j_in;
  residue_name_k = residue_k_in;
  residue_name_l = residue_l_in;
  residue_name_m = residue_m_in;
  frame_type = frame_type_in;
}

//-------------------------------------------------------------------------------------------------
ForceFieldElement::ForceFieldElement(const ParameterKind kind_in, const char4 atom_i_in,
                                     const char4 atom_j_in, const char4 atom_k_in,
                                     const char4 atom_l_in, const char4 atom_m_in,
                                     const char4 residue_i_in, const char4 residue_j_in,
                                     const char4 residue_k_in, const char4 residue_l_in,
                                     const char4 residue_m_in,
                                     const std::vector<double> &surface_in,
                                     const std::vector<int2> &locations_in) :
  ForceFieldElement(kind_in)
{
  switch (kind_in) {
  case ParameterKind::CMAP:
    break;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("Construction with five atom and residue names is not acceptable for a parameter of "
          "type \"" + getEnumerationName(kind_in) + "\".", "ForceFieldElement");
  }
  atom_name_i = atom_i_in;
  atom_name_j = atom_j_in;
  atom_name_k = atom_k_in;
  atom_name_l = atom_l_in;
  atom_name_m = atom_m_in;
  residue_name_i = residue_i_in;
  residue_name_j = residue_j_in;
  residue_name_k = residue_k_in;
  residue_name_l = residue_l_in;
  residue_name_m = residue_m_in;
  const size_t nsurf = surface_in.size();
  surface.resize(nsurf);
  locations.resize(nsurf);
  if (surface_in.size() != locations_in.size()) {
    rtErr("The provided surface of " + std::to_string(nsurf) + " elements does not have phi / psi "
          "indices for all points.", "ForceFieldElement");
  }
  for (size_t i = 0; i < nsurf; i++) {
    surface[i] = surface_in[i];
    locations[i] = locations_in[i];
  }
}

//-------------------------------------------------------------------------------------------------
ParameterKind ForceFieldElement::getKind() const {
  return kind;
}

//-------------------------------------------------------------------------------------------------
char4 ForceFieldElement::getNameOfAtom(const char atom_rank) const {
  switch (kind) {
  case ParameterKind::CMAP:
    switch (uppercase(atom_rank)) {
    case 'I':
      return atom_name_i;
    case 'J':
      return atom_name_j;
    case 'K':
      return atom_name_k;
    case 'L':
      return atom_name_l;
    case 'M':
      return atom_name_m;
    }
    break;
  case ParameterKind::CHARGE:
    switch(uppercase(atom_rank)) {
    case 'I':
      return atom_name_i;
    default:
      rtErr("Invalid atom rank " + std::to_string(atom_rank) + " for a charge parameter.  To get "
            "partial charge information, call getNameOfAtom() with no argument.",
            "ForceFieldElement", "getNameOfAtom");
    }
    break;
  case ParameterKind::VIRTUAL_SITE_FRAME:
    int number_rank;
    switch (uppercase(atom_rank)) {
    case 'I':
      number_rank = 0;
    case 'J':
      number_rank = 1;
      break;
    case 'K':
      number_rank = 2;
      break;
    case 'L':
      number_rank = 3;
      break;
    case 'M':
      number_rank = 4;
      break;
    default:
      rtErr("Invalid atom rank " + std::to_string(atom_rank) + " for a virtual site frame.",
            "ForceFieldElement", "getNameOfAtom");
    }
    if (number_rank == 0) {
      rtWarn("Virtual site frames contain valid information that can be accessed as atom " +
             std::to_string(atom_rank) + ", but is more sensibly accessed by "
             "getVirtualSiteAtom().", "ForceFieldElement", "getNameOfAtom");
    }
    else {
      rtWarn("Virtual site frames contain valid information that can be accessed as atom " +
             std::to_string(atom_rank) + ", but is more sensibly accessed by getFrameAtom(" +
             std::to_string(number_rank) + ").", "ForceFieldElement", "getNameOfAtom");
    }
    switch (number_rank) {
    case 0:
      return atom_name_i;
    case 1:
      return atom_name_j;
    case 2:
      return atom_name_k;
    case 3:
      return atom_name_l;
    case 4:
      return atom_name_m;
    }
    break;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::NONE:
    rtErr("There is no valid atom " + std::to_string(atom_rank) + " named in a parameter of "
          "type \"" + getEnumerationName(kind) + "\".", "ForceFieldElement", "getNameOfAtom");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
char4 ForceFieldElement::getTypeOfAtom(const char atom_rank) const {
  bool problem = false;
  switch (uppercase(atom_rank)) {
  case 'I':
    switch (kind) {
    case ParameterKind::BOND:
    case ParameterKind::ANGLE:
    case ParameterKind::DIHEDRAL:
    case ParameterKind::UREY_BRADLEY:
    case ParameterKind::CHARMM_IMPROPER:
    case ParameterKind::ATTN_14_SCALE:
    case ParameterKind::LENNARD_JONES:
    case ParameterKind::BUCKINGHAM:
      return atom_name_i;
    case ParameterKind::CMAP:
    case ParameterKind::CHARGE:
    case ParameterKind::VIRTUAL_SITE_FRAME:
    case ParameterKind::NONE:
      problem = true;
      break;
    }
    break;
  case 'J':
    switch (kind) {
    case ParameterKind::BOND:
    case ParameterKind::ANGLE:
    case ParameterKind::DIHEDRAL:
    case ParameterKind::UREY_BRADLEY:
    case ParameterKind::CHARMM_IMPROPER:
    case ParameterKind::ATTN_14_SCALE:
      return atom_name_j;
    case ParameterKind::LENNARD_JONES:
    case ParameterKind::BUCKINGHAM:
    case ParameterKind::CMAP:
    case ParameterKind::CHARGE:
    case ParameterKind::VIRTUAL_SITE_FRAME:
    case ParameterKind::NONE:
      problem = true;
      break;
    }
    break;
  case 'K':
    switch (kind) {
    case ParameterKind::ANGLE:
    case ParameterKind::DIHEDRAL:
    case ParameterKind::UREY_BRADLEY:
    case ParameterKind::CHARMM_IMPROPER:
    case ParameterKind::ATTN_14_SCALE:
      return atom_name_k;
    case ParameterKind::BOND:
    case ParameterKind::LENNARD_JONES:
    case ParameterKind::BUCKINGHAM:
    case ParameterKind::CMAP:
    case ParameterKind::CHARGE:
    case ParameterKind::VIRTUAL_SITE_FRAME:
    case ParameterKind::NONE:
      problem = true;
      break;
    }
    break;
  case 'L':
    switch (kind) {
    case ParameterKind::DIHEDRAL:
    case ParameterKind::CHARMM_IMPROPER:
    case ParameterKind::ATTN_14_SCALE:
      return atom_name_l;
    case ParameterKind::BOND:
    case ParameterKind::ANGLE:
    case ParameterKind::UREY_BRADLEY:
    case ParameterKind::LENNARD_JONES:
    case ParameterKind::BUCKINGHAM:
    case ParameterKind::CMAP:
    case ParameterKind::CHARGE:
    case ParameterKind::VIRTUAL_SITE_FRAME:
    case ParameterKind::NONE:
      problem = true;
      break;
    }
    break;
  default:
    problem = true;
    break;
  }
  if (problem) {
    rtErr("There is no valid atom " + std::to_string(atom_rank) + " named in a parameter of "
          "type \"" + getEnumerationName(kind) + "\".", "ForceFieldElement", "getTypeOfAtom");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
char4 ForceFieldElement::getNameOfResidue(const char atom_rank) const {
  switch (kind) {
  case ParameterKind::CMAP:
    switch (uppercase(atom_rank)) {
    case 'I':
      return residue_name_i;
    case 'J':
      return residue_name_j;
    case 'K':
      return residue_name_k;
    case 'L':
      return residue_name_l;
    case 'M':
      return residue_name_m;
    }
    break;
  case ParameterKind::CHARGE:
    switch(uppercase(atom_rank)) {
    case 'I':
      return residue_name_i;
    default:
      rtErr("Invalid atom rank " + std::to_string(atom_rank) + " for a charge parameter.  To get "
            "partial charge information, call getNameOfAtom() with no argument.",
            "ForceFieldElement", "getNameOfAtom");
    }
    break;
  case ParameterKind::VIRTUAL_SITE_FRAME:
    int number_rank;
    switch (uppercase(atom_rank)) {
    case 'I':
      number_rank = 0;
    case 'J':
      number_rank = 1;
      break;
    case 'K':
      number_rank = 2;
      break;
    case 'L':
      number_rank = 3;
      break;
    case 'M':
      number_rank = 4;
      break;
    default:
      rtErr("Invalid atom rank " + std::to_string(atom_rank) + " for a virtual site frame.",
            "ForceFieldElement", "getNameOfResidue");
    }
    if (number_rank == 0) {
      rtWarn("Virtual site frames contain valid information that can be accessed as atom " +
             std::to_string(atom_rank) + ", but is more sensibly accessed by "
             "getVirtualSiteResidue().", "ForceFieldElement", "getNameOfResidue");
    }
    else {
      rtWarn("Virtual site frames contain valid information that can be accessed as atom " +
             std::to_string(atom_rank) + ", but is more sensibly accessed by getFrameAtom(" +
             std::to_string(number_rank) + ").", "ForceFieldElement", "getNameOfResidue");
    }
    switch (number_rank) {
    case 0:
      return residue_name_i;
    case 1:
      return residue_name_j;
    case 2:
      return residue_name_k;
    case 3:
      return residue_name_l;
    case 4:
      return residue_name_m;
    }
    break;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::NONE:
    rtErr("There is no valid atom " + std::to_string(atom_rank) + " named in a parameter of "
          "type \"" + getEnumerationName(kind) + "\".", "ForceFieldElement", "getNameOfResidue");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getCharge() const {
  switch (kind) {
  case ParameterKind::CHARGE:
    return property_a;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no charge property.",
          "ForceFieldElement", "getCharge");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getSigma() const {
  switch (kind) {
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
    return property_a;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no sigma property.",
          "ForceFieldElement", "getSigma");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getEpsilon() const {
  switch (kind) {
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
    return property_b;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no epsilon property.",
          "ForceFieldElement", "getEpsilon");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getRho() const {
  switch (kind) {
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
    return property_c;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no rho property.",
          "ForceFieldElement", "getRho");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getStiffnessConstant() const {
  switch (kind) {
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
    return property_a;
  case ParameterKind::DIHEDRAL:
    rtWarn("A " + getEnumerationName(kind) + "\" has an amplitude, which is more correctly "
           "accessed with getAmplitude(), but the amplitude is equivalent to a stiffness.",
           "ForceFieldElement", "getStiffnessConstant");
    return property_a;
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no stiffness constant.",
          "ForceFieldElement", "getStiffnessConstant");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getEquilibriumConstant() const {
  switch (kind) {
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
    return property_b;
  case ParameterKind::DIHEDRAL:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no equilibrium "
          "constant.", "ForceFieldElement", "getEquilibriumConstant");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getAmplitude() const {
  switch (kind) {
  case ParameterKind::DIHEDRAL:
    return property_a;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no amplitude.",
          "ForceFieldElement", "getAmplitude");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getPhaseAngle() const {
  switch (kind) {
  case ParameterKind::DIHEDRAL:
  case ParameterKind::CHARMM_IMPROPER:
    return property_b;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no phase angle.",
          "ForceFieldElement", "getPhaseAngle");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getPeriodicity() const {
  switch (kind) {
  case ParameterKind::DIHEDRAL:
    return property_c;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no periodicity.",
          "ForceFieldElement", "getPeriodicity");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getElectrostaticScaling() const {
  switch (kind) {
  case ParameterKind::ATTN_14_SCALE:
    return property_a;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no 1:4 scaling factor.",
          "ForceFieldElement", "getElectrostaticScaling");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double ForceFieldElement::getVanDerWaalsScaling() const {
  switch (kind) {
  case ParameterKind::ATTN_14_SCALE:
    return property_b;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no 1:4 scaling factor.",
          "ForceFieldElement", "getVanDerWaalsScaling");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
TorsionKind ForceFieldElement::getTorsionKind() const {
  switch (kind) {
  case ParameterKind::DIHEDRAL:
    return torsion_kind;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" does not have a "
          "torsion kind associated with it.", "ForceFieldElement", "getTorsionKind");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<double> ForceFieldElement::getSurfaceValues() const {
  switch (kind) {
  case ParameterKind::CMAP:
    return surface;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" does not have any "
          "two-dimensional potential surface elements.", "ForceFieldElement", "getSurfaceValues");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<int2> ForceFieldElement::getSurfaceIndices() const {
  switch (kind) {
  case ParameterKind::CMAP:
    return locations;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" does not have any "
          "two-dimensional potential surface elements.", "ForceFieldElement", "getSurface");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
VirtualSiteKind ForceFieldElement::getVirtualSiteFrameType() const {
  switch (kind) {
  case ParameterKind::VIRTUAL_SITE_FRAME:
    return frame_type;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" does not have a "
          "virtual site frame type.", "ForceFieldElement", "getVirtualSiteFrameType");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool ForceFieldElement::testSigmaModification() const {
  switch (kind) {
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
    return activate_a;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no sigma property.",
          "ForceFieldElement", "testSigmaModification");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool ForceFieldElement::testEpsilonModification() const {
  switch (kind) {
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
    return activate_b;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no epsilon property.",
          "ForceFieldElement", "testEpsilonModification");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool ForceFieldElement::testRhoModification() const {
  switch (kind) {
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
    return activate_c;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no rho property.",
          "ForceFieldElement", "testRhoModification");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool ForceFieldElement::testStiffnessModification() const {
  switch (kind) {
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
    return activate_a;
  case ParameterKind::DIHEDRAL:
    rtWarn("A " + getEnumerationName(kind) + "\" has an amplitude, which is more correctly "
           "accessed with getAmplitude(), but the amplitude is equivalent to a stiffness.",
           "ForceFieldElement", "testStiffnessModification");
    return activate_a;
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no stiffness constant.",
          "ForceFieldElement", "testStiffnessModification");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool ForceFieldElement::testEquilibriumModification() const {
  switch (kind) {
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
    return activate_b;
  case ParameterKind::DIHEDRAL:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no equilibrium "
          "constant.", "ForceFieldElement", "testEquilibriumModification");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool ForceFieldElement::testAmplitudeModification() const {
  switch (kind) {
  case ParameterKind::DIHEDRAL:
    return property_a;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no amplitude.",
          "ForceFieldElement", "testAmplitudeModification");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool ForceFieldElement::testPhaseAngleModification() const {
  switch (kind) {
  case ParameterKind::DIHEDRAL:
  case ParameterKind::CHARMM_IMPROPER:
    return property_b;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no phase angle.",
          "ForceFieldElement", "testPhaseAngleModification");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
bool ForceFieldElement::testPeriodicityModification() const {
  switch (kind) {
  case ParameterKind::DIHEDRAL:
    return property_c;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no periodicity.",
          "ForceFieldElement", "testPeriodicityModification");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void ForceFieldElement::setStiffness(const double stiffness_in) {
  switch (kind) {
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
    property_a = stiffness_in;
    activate_a = true;
    break;
  case ParameterKind::DIHEDRAL:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no stiffness property.",
          "ForceFieldElement", "setStiffness");
  }
}

//-------------------------------------------------------------------------------------------------
void ForceFieldElement::setEquilibrium(const double equilibrium_in) {
  switch (kind) {
  case ParameterKind::BOND:
  case ParameterKind::UREY_BRADLEY:
    property_b = equilibrium_in;
    activate_b = true;
    break;
  case ParameterKind::ANGLE:
    property_b = equilibrium_in * symbols::pi / 180.0;
    activate_b = true;
    break;
  case ParameterKind::DIHEDRAL:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no equilibrium "
          "property.", "ForceFieldElement", "setEquilibrium");
  }
}

//-------------------------------------------------------------------------------------------------
void ForceFieldElement::setPhaseAngle(const double phase_angle_in) {
  switch (kind) {
  case ParameterKind::DIHEDRAL:
  case ParameterKind::CHARMM_IMPROPER:
    property_b = phase_angle_in * symbols::pi / 180.0;
    activate_b = true;
    break;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no phase angle "
          "property.", "ForceFieldElement", "setPhaseAngle");
  }
}

//-------------------------------------------------------------------------------------------------
void ForceFieldElement::setAmplitude(const double amplitude_in) {
  switch (kind) {
  case ParameterKind::DIHEDRAL:
    property_a = amplitude_in;
    activate_a = true;
    break;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no amplitude property.",
          "ForceFieldElement", "setAmplitude");
  }
}

//-------------------------------------------------------------------------------------------------
void ForceFieldElement::setPeriodicity(const double periodicity_in) {
  switch (kind) {
  case ParameterKind::DIHEDRAL:
    property_c = periodicity_in;
    activate_c = true;
    break;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no periodicity.",
          "ForceFieldElement", "setPeriodicity");
  }
}

//-------------------------------------------------------------------------------------------------
void ForceFieldElement::setChargeScaling(const double scaling_in) {
  switch (kind) {
  case ParameterKind::ATTN_14_SCALE:
    property_a = scaling_in;
    activate_a = true;
    break;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no scaling factors.",
          "ForceFieldElement", "setChargeScaling");
  }
}

//-------------------------------------------------------------------------------------------------
void ForceFieldElement::setVanDerWaalsScaling(const double scaling_in) {
  switch (kind) {
  case ParameterKind::ATTN_14_SCALE:
    property_b = scaling_in;
    activate_b = true;
    break;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no scaling factors.",
          "ForceFieldElement", "setVanDerWaalsScaling");
  }
}

//-------------------------------------------------------------------------------------------------
void ForceFieldElement::setCharge(const double charge_in) {
  switch (kind) {
  case ParameterKind::CHARGE:
    property_a = charge_in;
    activate_a = true;
    break;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no charge property.",
          "ForceFieldElement", "setCharge");
  }
}

//-------------------------------------------------------------------------------------------------
void ForceFieldElement::setSigma(const double sigma_in) {
  switch (kind) {
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
    property_a = sigma_in;
    activate_a = true;
    break;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no Sigma property.",
          "ForceFieldElement", "setSigma");
  }
}

//-------------------------------------------------------------------------------------------------
void ForceFieldElement::setEpsilon(const double epsilon_in) {
  switch (kind) {
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
    property_b = epsilon_in;
    activate_b = true;
    break;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no Epsilon property.",
          "ForceFieldElement", "setEpsilon");
  }
}

//-------------------------------------------------------------------------------------------------
void ForceFieldElement::setRho(const double rho_in) {
  switch (kind) {
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
    property_c = rho_in;
    activate_c = true;
    break;
  case ParameterKind::BOND:
  case ParameterKind::ANGLE:
  case ParameterKind::DIHEDRAL:
  case ParameterKind::UREY_BRADLEY:
  case ParameterKind::CHARMM_IMPROPER:
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    rtErr("A parameter of type \"" + getEnumerationName(kind) + "\" has no Rho property.",
          "ForceFieldElement", "setRho");
  }
}

//-------------------------------------------------------------------------------------------------
void ForceFieldElement::apply(AtomGraph *ag, const ExceptionResponse policy) const {
  switch (kind) {
  case ParameterKind::BOND:
    ag->setBondParameters(property_a, property_b, activate_a, activate_b, atom_name_i,
                          atom_name_j, policy);
  case ParameterKind::ANGLE:
    ag->setAngleParameters(property_a, property_b, activate_a, activate_b, atom_name_i,
                           atom_name_j, atom_name_k, policy);
  case ParameterKind::DIHEDRAL:
    ag->setDihedralParameters(property_a, property_b, activate_a, activate_b, atom_name_i,
                              atom_name_j, atom_name_k, atom_name_l, property_c, policy);
  case ParameterKind::UREY_BRADLEY:
    ag->setUreyBradleyParameters(property_a, property_b, activate_a, activate_b, atom_name_i,
                                 atom_name_j, atom_name_k, policy);
  case ParameterKind::CHARMM_IMPROPER:
    ag->setCharmmImprParameters(property_a, property_b, activate_a, activate_b, atom_name_i,
                                atom_name_j, atom_name_k, atom_name_l, policy);
  case ParameterKind::CMAP:
  case ParameterKind::ATTN_14_SCALE:
  case ParameterKind::CHARGE:
  case ParameterKind::LENNARD_JONES:
  case ParameterKind::BUCKINGHAM:
  case ParameterKind::VIRTUAL_SITE_FRAME:
  case ParameterKind::NONE:
    break;
  }
}

} // namespace modeling
} // namespace stormm
