#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "atommask.h"
#include "znumber.h"

namespace stormm {
namespace chemistry {

using card::HybridKind;
using card::HybridTargetLevel;
using stmath::accumulateBitmask;
using stmath::readBitFromMask;
using parse::char4ToString;
using parse::NumberFormat;
using parse::TextGuard;
using parse::realToString;
using parse::resolveScopes;
using parse::separateText;
using parse::stringToChar4;
using parse::strcmpWildCard;
using data_types::operator==;
using topology::ChemicalDetailsKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
  
//-------------------------------------------------------------------------------------------------
SelectionItem::SelectionItem(const SelectionItemKind kind_in, const char4 name_in,
                             const std::vector<WildCardKind> &wildcards_in) :
    kind{kind_in},
    name{name_in},
    wildcards{wildcards_in},
    scan_begin{0},
    scan_end{0}
{}

//-------------------------------------------------------------------------------------------------
SelectionItem::SelectionItem(const SelectionItemKind kind_in, const int begin_in,
                             const int end_in) :
    kind{kind_in},
    name{'N', 'O', 'N', 'E'},
    wildcards{},
    scan_begin{begin_in},
    scan_end{(end_in >= begin_in) ? end_in : begin_in + 1}
{}

//-------------------------------------------------------------------------------------------------
MaskComponent::MaskComponent(const MaskOperator op_in, const std::string &basis_in,
                             const int start_idx, const int end_idx) :
    kind{MaskComponentKind::OPERATOR},
    op{op_in},
    range{0.0},
    primitive_mask{},
    text_basis{basis_in},
    text_limits{start_idx, end_idx}
{}

//-------------------------------------------------------------------------------------------------
MaskComponent::MaskComponent(const MaskOperator op_in, const double range_in,
                             const std::string &basis_in, const int start_idx, const int end_idx) :
    kind{MaskComponentKind::OPERATOR},
    op{op_in},
    range{range_in},
    primitive_mask{},
    text_basis{basis_in},
    text_limits{start_idx, end_idx}
{}

//-------------------------------------------------------------------------------------------------
MaskComponent::MaskComponent(const std::vector<SelectionItem> &parts_in, const AtomGraph *ag,
                             const ChemicalFeatures &chemfe, const std::string &basis_in,
                             const int start_idx, const int end_idx) :
    kind{MaskComponentKind::MASK},
    op{MaskOperator::NONE},
    range{0.0},
    primitive_mask{std::vector<uint>((ag->getAtomCount() + (sizeof(uint) * 8) - 1) /
                                     (sizeof(uint) * 8), 0U)},
    text_basis{basis_in},
    text_limits{start_idx, end_idx}
{
  // Fetch an abstract into the topology's atom names for rapid access
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit(HybridTargetLevel::HOST);

  // Pre-compute some factors for the bitmask conversion (also needed for deconstructing the
  // bitmasks produced by chemical features masks)
  const int n_prim = primitive_mask.size();
  const int n_bits = sizeof(uint) * 8;
  const int n_uint = (cdk.natom + n_bits - 1) / n_bits;
  
  // This is the only mask that requires any additional calculation.  Loop over all parts.
  const int n_parts = parts_in.size();
  std::vector<bool> active_mask(cdk.natom, false);

  // Loop over all parts of the mask
  for (int i = 0; i < n_parts; i++) {
    switch (parts_in[i].kind) {
    case SelectionItemKind::ATOM_NAME:
      {
        const std::string query = char4ToString(parts_in[i].name);
        std::string target(4, ' ');
        for (int j = 0; j < cdk.natom; j++) {
          target[0] = cdk.atom_names[j].x;
          target[1] = cdk.atom_names[j].y;
          target[2] = cdk.atom_names[j].z;
          target[3] = cdk.atom_names[j].w;

          // The more elaborate evaluation, using strcmpWildCard(), is needed to ensure that
          // leading whitespace does not interfere with the comparison and that wildcard characters
          // will be properly handled.
          active_mask[j] = (active_mask[j] ||
                            strcmpWildCard(target, query, parts_in[i].wildcards));
        }
      }
      break;
    case SelectionItemKind::ATOM_NUMBER:
      for (int j = 0; j < cdk.natom; j++) {
        active_mask[j] = (active_mask[j] || (cdk.atom_numbers[j] >= parts_in[i].scan_begin &&
                                             cdk.atom_numbers[j] <= parts_in[i].scan_end));
      }
      break;
    case SelectionItemKind::ELEMENT_NAME:
      {
        const std::string query = char4ToString(parts_in[i].name);
        std::string target(4, ' ');
        for (int j = 0; j < cdk.natom; j++) {
          const char2 tmp_zn = zNumberToSymbol(cdk.z_numbers[j]);
          target[0] = tmp_zn.x;
          target[1] = tmp_zn.y;
          active_mask[j] = (active_mask[j] ||
                            strcmpWildCard(target, query, parts_in[i].wildcards));
        }
      }
      break;
    case SelectionItemKind::ELEMENT_NUMBER:
      for (int j = 0; j < cdk.natom; j++) {
        active_mask[j] = (active_mask[j] || (cdk.z_numbers[j] >= parts_in[i].scan_begin &&
                                             cdk.z_numbers[j] <= parts_in[i].scan_end));
      }
      break;
    case SelectionItemKind::ATOM_TYPE:
      {
        const std::string query = char4ToString(parts_in[i].name);
        std::string target(4, ' ');
        for (int j = 0; j < cdk.natom; j++) {
          target[0] = cdk.atom_types[j].x;
          target[1] = cdk.atom_types[j].y;
          target[2] = cdk.atom_types[j].z;
          target[3] = cdk.atom_types[j].w;
          active_mask[j] = (active_mask[j] ||
                            strcmpWildCard(target, query, parts_in[i].wildcards));
        }
      }
      break;
    case SelectionItemKind::RESIDUE_NAME:
      {
        const std::string query = char4ToString(parts_in[i].name);
        std::string target(4, ' ');
        for (int j = 0; j < cdk.nres; j++) {
          target[0] = cdk.res_names[j].x;
          target[1] = cdk.res_names[j].y;
          target[2] = cdk.res_names[j].z;
          target[3] = cdk.res_names[j].w;
          if (strcmpWildCard(target, query, parts_in[i].wildcards)) {
            for (int k = cdk.res_limits[j]; k < cdk.res_limits[j + 1]; k++) {
              active_mask[k] = true;
            }
          }
        }
      }
      break;
    case SelectionItemKind::RESIDUE_NUMBER:
      for (int j = 0; j < cdk.natom; j++) {
        active_mask[j] = (active_mask[j] || (cdk.res_numbers[j] >= parts_in[i].scan_begin &&
                                             cdk.res_numbers[j] <= parts_in[i].scan_end));
      }
      break;
#if 0
    case SelectionItemKind::RING_SIZE:
      {
        const std::vector<uint> ring_mask = chemfe.getRingMask(parts_in[i].scan_begin,
                                                               parts_in[i].scan_end);
        for (int j = 0; j < cdk.natom; j++) {
          const int jidx = j / n_bits;
          const int jbit = j - (jidx * n_bits);
          active_mask[j] = (active_mask[j] || ((ring_mask[jidx] >> jbit) & 0x1));
        }
      }
      break;
    case SelectionItemKind::AROMATICITY:
      {
        const std::vector<uint> aromatic_mask = chemfe.getAromaticMask(parts_in[i].scan_begin,
                                                                       parts_in[i].scan_end);
        for (int j = 0; j < cdk.natom; j++) {
          const int jidx = j / n_bits;
          const int jbit = j - (jidx * n_bits);
          active_mask[j] = (active_mask[j] || ((aromatic_mask[jidx] >> jbit) & 0x1));
        }
      }
      break;
    case SelectionItemKind::CHIRALITY:
      {
        const std::vector<uint> chiral_mask =
          chemfe.getChiralityMask(static_cast<ChiralOrientation>(parts_in[i].begin));
        for (int j = 0; j < cdk.natom; j++) {
          const int jidx = j / n_bits;
          const int jbit = j - (jidx * n_bits);
          active_mask[j] = (active_mask[j] || ((chiral_mask[jidx] >> jbit) & 0x1));
        }
      }
      break;
#endif
    case SelectionItemKind::NONE:
      break;
    }
  }

  // Convert the active mask into a primitive mask
  int k = 0;
  for (int i = 0; i < n_prim; i++) {
    ulint pbuffer = 0U;
    for (int j = 0; j < n_bits; j++) {
      if (active_mask[k]) {
        pbuffer |= (0x1 << j);
      }
      k++;
    }
    primitive_mask[i] = pbuffer;
  }
}

//-------------------------------------------------------------------------------------------------
MaskComponent::MaskComponent(const std::vector<uint> &primitive_mask_in,
                             const std::string &basis_in, const int start_idx, const int end_idx) :
    kind{MaskComponentKind::MASK},
    op{MaskOperator::NONE},
    range{0.0},
    primitive_mask{primitive_mask_in},
    text_basis{basis_in},
    text_limits{start_idx, end_idx}
{}
  
//-------------------------------------------------------------------------------------------------
MaskComponentKind MaskComponent::getKind() const {
  return kind;
}

//-------------------------------------------------------------------------------------------------
MaskOperator MaskComponent::getOperator() const {
  if (kind == MaskComponentKind::MASK) {
    rtErr("This mask component contains a primitive mask.", "MaskComponent", "getOperator");
  }
  return op;
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> MaskComponent::getMask() const {
  if (kind == MaskComponentKind::OPERATOR) {
    rtErr("This mask component contains an operator.", "MaskComponent", "getMask");
  }
  return primitive_mask;
}

//-------------------------------------------------------------------------------------------------
std::string MaskComponent::getText() const {
  return text_basis;
}

//-------------------------------------------------------------------------------------------------
int2 MaskComponent::getTextLimits() const {
  return text_limits;
}
  
//-------------------------------------------------------------------------------------------------
void MaskComponent::applyNotOperator() {
  for (size_t i = 0; i < primitive_mask.size(); i++) {
    primitive_mask[i] = ~(primitive_mask[i]);
  }
}

//-------------------------------------------------------------------------------------------------
void MaskComponent::applyAndOperator(const std::vector<uint> &other) {
  for (size_t i = 0; i < primitive_mask.size(); i++) {
    primitive_mask[i] &= other[i];
  }
}

//-------------------------------------------------------------------------------------------------
void MaskComponent::applyOrOperator(const std::vector<uint> &other) {
  for (size_t i = 0; i < primitive_mask.size(); i++) {
    primitive_mask[i] |= other[i];
  }
}

//-------------------------------------------------------------------------------------------------
MaskComponent MaskComponent::applyRangeOperator(const std::vector<uint> &other,
                                                const AtomGraph *ag,
                                                const CoordinateFrameReader &cfr) {

  // This and other range-based operators will benefit from being able to stride through the
  // primitive mask by looking at each unsigned integer as single number.  If it is zero, then
  // there are no bits representing particles in the mask in those indices.  For most machines,
  // this means 32x faster scanning.
  const int uint_bits = sizeof(uint) * 8;
  std::vector<size_t> nodes;
  for (int i = 0; i < primitive_mask.size(); i++) {
    uint numbermask = primitive_mask[i];
    uint multiplier = 0x1;
    int j = 0;
    while (numbermask > 0) {
      if ((numbermask >> j) & 0x1) {
        numbermask -= multiplier;
        nodes.push_back((uint_bits * i) + j);
      }
      j++;
      multiplier *= 0x2;
    }
  }
  const double cut_squared = range * range;

  // Create a temporary mask
  bool ibval;
  switch (op) {
  case MaskOperator::ATOM_LT:
  case MaskOperator::ATOM_LE:
  case MaskOperator::RESIDUE_LT:
  case MaskOperator::RESIDUE_LE:
    ibval = false;
    break;
  case MaskOperator::ATOM_GT:
  case MaskOperator::ATOM_GE:
  case MaskOperator::RESIDUE_GT:
  case MaskOperator::RESIDUE_GE:
    ibval = true;
    break;
  case MaskOperator::AND:
  case MaskOperator::OR:
  case MaskOperator::NOT:
  case MaskOperator::NONE:
    break;
  }
  std::vector<bool> tmp_mask(cfr.natom, ibval);
  const size_t nnode = nodes.size();
  for (size_t i = 0; i < cfr.natom; i++) {
    const double xpos = cfr.xcrd[i];
    const double ypos = cfr.ycrd[i];
    const double zpos = cfr.zcrd[i];
    switch (op) {
    case MaskOperator::ATOM_LT:
    case MaskOperator::RESIDUE_LT:
      {
        size_t j = 0;
        while (j < nnode && !tmp_mask[i]) {
          const double dx = cfr.xcrd[nodes[j]] - xpos;
          const double dy = cfr.ycrd[nodes[j]] - ypos;
          const double dz = cfr.zcrd[nodes[j]] - zpos;
          tmp_mask[i] = ((dx * dx) + (dy * dy) + (dz * dz) < cut_squared);
          j++;
        }
      }
      break;
    case MaskOperator::ATOM_GT:
    case MaskOperator::RESIDUE_GT:
      {
        size_t j = 0;
        while (j < nnode && tmp_mask[i]) {
          const double dx = cfr.xcrd[nodes[j]] - xpos;
          const double dy = cfr.ycrd[nodes[j]] - ypos;
          const double dz = cfr.zcrd[nodes[j]] - zpos;
          tmp_mask[i] = ((dx * dx) + (dy * dy) + (dz * dz) > cut_squared);
          j++;
        }
      }
      break;
    case MaskOperator::ATOM_LE:
    case MaskOperator::RESIDUE_LE:
      {
        size_t j = 0;
        while (j < nnode && !tmp_mask[i]) {
          const double dx = cfr.xcrd[nodes[j]] - xpos;
          const double dy = cfr.ycrd[nodes[j]] - ypos;
          const double dz = cfr.zcrd[nodes[j]] - zpos;
          tmp_mask[i] = ((dx * dx) + (dy * dy) + (dz * dz) <= cut_squared);
          j++;
        }
      }
      break;
    case MaskOperator::ATOM_GE:
    case MaskOperator::RESIDUE_GE:
      {
        size_t j = 0;
        while (j < nnode && tmp_mask[i]) {
          const double dx = cfr.xcrd[nodes[j]] - xpos;
          const double dy = cfr.ycrd[nodes[j]] - ypos;
          const double dz = cfr.zcrd[nodes[j]] - zpos;
          tmp_mask[i] = ((dx * dx) + (dy * dy) + (dz * dz) >= cut_squared);
          j++;
        }
      }
      break;
    case MaskOperator::AND:
    case MaskOperator::OR:
    case MaskOperator::NOT:
    case MaskOperator::NONE:
      break;
    }
  }

  // For residue-based range masks, step through the new atom mask and determine whether each
  // residue, as a whole, meets the criteria.  Only perform this stage, and only bother to get
  // the residue limits out of the topology, if the operatior is residue-based.
  if (op == MaskOperator::RESIDUE_LT || op == MaskOperator::RESIDUE_LE ||
      op == MaskOperator::RESIDUE_GT || op == MaskOperator::RESIDUE_GE) {
    const int nres = ag->getResidueCount();
    const int* tmp_residue_limits = ag->getResidueLimits().data();
    for (int i = 0; i < nres; i++) {
      bool res_works;
      if (op == MaskOperator::RESIDUE_LT || op == MaskOperator::RESIDUE_GT) {
        res_works = true;
        for (int j = tmp_residue_limits[i]; j < tmp_residue_limits[i + 1]; j++) {
          res_works = (res_works && tmp_mask[j]);
        }
      }
      else {
        res_works = false;
        for (int j = tmp_residue_limits[i]; j < tmp_residue_limits[i + 1]; j++) {
          res_works = (res_works || tmp_mask[j]);
        }
      }
      if (res_works) {
        for (int j = tmp_residue_limits[i]; j < tmp_residue_limits[i + 1]; j++) {
          tmp_mask[j] = true;
        }
      }
    }
  }

  // Convert the output mask into a bit string, similar to the input mask
  int j = 0;
  uint out_number = 0x0;
  std::vector<uint> outmask((cfr.natom + (uint_bits - 1)) / uint_bits, 0x0);
  size_t pos = 0;
  for (int i = 0; i < cfr.natom; i++) {
    ulint test = tmp_mask[i];
    out_number |= (test << j);
    j++;
    if (j == uint_bits) {
      outmask[pos] = out_number;
      j = 0;
      out_number = 0x0;
      pos++;
    }
  }
  if (j > 0) {
    outmask[pos] = out_number;
  }
  return MaskComponent(outmask, text_basis, text_limits.x, text_limits.y);
}
  
//-------------------------------------------------------------------------------------------------
AtomMask::AtomMask(const AtomGraph *ag_in) :
    recommended_scan{MaskTraversalMode::COMPLETE},
    style{MaskInputMode::AMBMASK},
    masked_atom_count{0},
    raw_mask{},
    segments{},
    input_text{std::string("")},
    description{std::string("")},
    ag_pointer{ag_in}
{
  if (ag_pointer != nullptr) {
    raw_mask = std::vector<uint>((ag_pointer->getAtomCount() + (sizeof(uint) * 8) - 1) /
                                 (sizeof(uint) * 8), 0U);
  }
}

//-------------------------------------------------------------------------------------------------
AtomMask::AtomMask(const std::string &input_text_in, const AtomGraph *ag_in,
                   const ChemicalFeatures &chemfe, const CoordinateFrameReader &cfr,
                   const MaskInputMode mode, const std::string &description_in) :
    recommended_scan{MaskTraversalMode::COMPLETE},
    style{mode},
    raw_mask{std::vector<uint>((ag_in->getAtomCount() + (sizeof(uint) * 8) - 1) /
                               (sizeof(uint) * 8), 0U)},
    segments{},
    input_text{input_text_in},
    description{description_in},
    ag_pointer{ag_in}
{
  // Parse the mask text into its scopes.  Aside from the outermost scope (level 0), scopes are
  // denoted by '(' ... ')', '{' ... '}', and '[' ... ']'.
  const std::vector<TextGuard> scope_bounds = { TextGuard("(", ")"), TextGuard("{", "}"),
                                                TextGuard("[", "]") };
  const std::vector<int> scope_levels = resolveScopes(input_text_in, scope_bounds);

  // Parse the input test using a recursive strategy.  The recursive function stores a sequence of
  // operators ('!', '&', '|', '>:', '<:', '>@', and '<@') and masks obtained from evaluating
  // lists of atoms or residues.  At this point, all scope delimiters count as white space and
  // scope navigation proceeds by the scope levels array.
  int eval_pos = 0;
  switch (mode) {
  case MaskInputMode::AMBMASK:
    raw_mask = parseMask(scope_levels, &eval_pos, cfr, chemfe);
    break;
  case MaskInputMode::VMD:
    break;
  }

  // Get a count of the atoms that have been masked out
  const int natom = (ag_pointer == nullptr) ? 0 : ag_pointer->getAtomCount();
  const int uint_bits = sizeof(uint) * 8;
  int tmp_ma_count = 0;
  for (int i = 0; i < natom; i++) {
    const int pull_index = i / uint_bits;
    tmp_ma_count += ((raw_mask[pull_index] >> (i - (pull_index * uint_bits))) & 0x1);
  }
  masked_atom_count = tmp_ma_count;  
}

//-------------------------------------------------------------------------------------------------
AtomMask::AtomMask(const std::string &input_text_in, const AtomGraph *ag_in,
                   const ChemicalFeatures &chemfe, const CoordinateFrame &cf,
                   const MaskInputMode mode, const std::string &description_in) :
    AtomMask(input_text_in, ag_in, chemfe, cf.data(), mode, description_in)
{}
  
//-------------------------------------------------------------------------------------------------
AtomMask::AtomMask(const std::string &input_text_in, const AtomGraph *ag_in,
                   const ChemicalFeatures &chemfe, const PhaseSpace &ps,
                   const MaskInputMode mode, const std::string &description_in) :
    AtomMask(input_text_in, ag_in, chemfe, CoordinateFrameReader(ps), mode, description_in)
{}
  
//-------------------------------------------------------------------------------------------------
MaskTraversalMode AtomMask::getRecommendation() const {
  return recommended_scan;
}

//-------------------------------------------------------------------------------------------------
const std::vector<uint>& AtomMask::getRawMask() const {
  return raw_mask;
}

//-------------------------------------------------------------------------------------------------
int AtomMask::getMaskedAtomCount() const {
  return masked_atom_count;
}

//-------------------------------------------------------------------------------------------------
std::vector<bool> AtomMask::getMask() const {
  const int natom = (ag_pointer == nullptr) ? 0 : ag_pointer->getAtomCount();
  std::vector<bool> result(natom, false);
  const int uint_bits = sizeof(uint) * 8;
  for (int i = 0; i < natom; i++) {
    const int pull_index = i / uint_bits;
    result[i] = ((raw_mask[pull_index] >> (i - (pull_index * uint_bits))) & 0x1);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> AtomMask::getMaskedAtomList() const {
  const int natom = (ag_pointer == nullptr) ? 0 : ag_pointer->getAtomCount();
  std::vector<int> result;
  const int uint_bits = sizeof(uint) * 8;
  for (int i = 0; i < natom; i++) {
    const int pull_index = i / uint_bits;
    if ((raw_mask[pull_index] >> (i - (pull_index * uint_bits))) & 0x1) {
      result.push_back(i);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int2> AtomMask::getSegments() const {
  return segments;
}

//-------------------------------------------------------------------------------------------------
std::string AtomMask::getInputText() const {
  return input_text;
}

//-------------------------------------------------------------------------------------------------
MaskInputMode AtomMask::getInputKind() const {
  return style;
}

//-------------------------------------------------------------------------------------------------
std::string AtomMask::getDescription() const {
  return description;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* AtomMask::getTopologyPointer() const {
  return ag_pointer;
}

//-------------------------------------------------------------------------------------------------
bool AtomMask::isAtomInMask(const char4 atom_name) const {
  const Hybrid<char4>& ag_atom_names = ag_pointer->getResidueName();
  const char4* ag_atom_names_ptr = ag_atom_names.data();
  const size_t natom = ag_atom_names.size();
  for (size_t i = 0; i < natom; i++) {
    if (ag_atom_names_ptr[i] == atom_name) {
      return true;
    }
  }
  return false;
}

//-------------------------------------------------------------------------------------------------
bool AtomMask::isAtomInMask(const char4 res_name, const char4 atom_name) const {
  const Hybrid<int>& ag_res_lims = ag_pointer->getResidueLimits();
  const Hybrid<char4>& ag_res_names = ag_pointer->getResidueName();
  const Hybrid<char4>& ag_atom_names = ag_pointer->getAtomName();
  const int* ag_res_lims_ptr = ag_res_lims.data();
  const char4* ag_res_names_ptr = ag_res_names.data();
  const char4* ag_atom_names_ptr = ag_atom_names.data();
  const size_t nres = ag_res_lims.size();
  for (size_t i = 0; i < nres; i++) {
    if (ag_res_names_ptr[i] == res_name) {
      for (size_t j = ag_res_lims_ptr[i]; j < ag_res_lims_ptr[i + 1]; j++) {
        if (ag_atom_names_ptr[j] == atom_name) {
          return true;
        }
      }
    }
  }
  return false;
}

//-------------------------------------------------------------------------------------------------
bool AtomMask::isAtomInMask(int atom_index) const {
  return readBitFromMask(raw_mask, atom_index);
}

//-------------------------------------------------------------------------------------------------
void AtomMask::addAtoms(const std::vector<int> &new_indices, const ExceptionResponse policy) {
  const int n_adds = new_indices.size();
  const int n_atoms = ag_pointer->getAtomCount();
  std::string mask_addendum;
  bool addendum_started = false;
  for (int i = 0; i < n_adds; i++) {
    if (new_indices[i] < 0 || new_indices[i] >= n_atoms) {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Atom " + std::to_string(new_indices[i] + 1) + " is invalid for a topology with " +
              std::to_string(n_atoms) + " atoms.", "AtomMask", "addAtoms");
      case ExceptionResponse::WARN:
        rtWarn("Atom " + std::to_string(new_indices[i] + 1) + " is invalid for a topology with " +
               std::to_string(n_atoms) + " atoms.  This entry will not be added.", "AtomMask",
               "addAtoms");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
    else if (readBitFromMask(raw_mask, new_indices[i]) == 0) {
      accumulateBitmask(&raw_mask, new_indices[i]);
      segments.push_back({new_indices[i], new_indices[i] + 1});
      masked_atom_count += 1;
      if (addendum_started) {
        mask_addendum.append("," + std::to_string(new_indices[i]));
      }
      else {
        mask_addendum.append("@" + std::to_string(new_indices[i]));
        addendum_started = true;
      }
    }
  }

  // Update the mask string with the added atoms
  if (addendum_started) {
    input_text = "(" + input_text + ") | (" + mask_addendum + ")";
  }
}

//-------------------------------------------------------------------------------------------------
void AtomMask::addAtoms(const std::vector<char4> &new_names) {
  const int n_names = new_names.size();
  const int n_atoms = ag_pointer->getAtomCount();
  const char4* topology_names = ag_pointer->getAtomName().data();
  std::string mask_addendum;
  bool addendum_started = false;
  for (int i = 0; i < n_names; i++) {
    const char4 tname = new_names[i];
    for (int j = 0; j < n_atoms; j++) {
      if (readBitFromMask(raw_mask, j) == 0 && tname == topology_names[j]) {
        accumulateBitmask(&raw_mask, j);
        segments.push_back({j, j + 1});
        masked_atom_count += 1;
        if (addendum_started) {
          mask_addendum.append("," + char4ToString(new_names[i]));
        }
        else {
          mask_addendum.append("@" + char4ToString(new_names[i]));
          addendum_started = true;
        }
      }
    }
  }

  // Update the mask string with the added atom names
  if (addendum_started) {
    input_text = "(" + input_text + ") | (" + mask_addendum + ")";
  }
}

//-------------------------------------------------------------------------------------------------
void AtomMask::addAtoms(const AtomMask &new_mask, const CoordinateFrame &cf,
                        const ChemicalFeatures &chemfe) {

  // Determine if the two masks point to the same topology.  If so, the addition can be a simple
  // union of sets.  Otherwise, take the input text from the new mask and apply it to the current
  // mask's topology.
  if (new_mask.getTopologyPointer() == ag_pointer) {
    const std::vector<uint>& new_raw_mask = new_mask.getRawMask();
    const size_t nr_val = raw_mask.size();
    for (size_t i = 0LLU; i < nr_val; i++) {
      raw_mask[i] |= new_raw_mask[i];
    }
    const int natom = ag_pointer->getAtomCount();
    int update_ac = 0;
    for (int i = 0; i < natom; i++) {
      update_ac += readBitFromMask(raw_mask, i);
    }
    if (update_ac > masked_atom_count) {
      input_text = "(" + input_text + ") | (" + new_mask.getInputText() + ")";
    }
    masked_atom_count = update_ac;
  }
  else {
    addAtoms(new_mask.getInputText(), cf, chemfe);
  }
}

//-------------------------------------------------------------------------------------------------
void AtomMask::addAtoms(const std::string &new_mask, const CoordinateFrame &cf,
                        const ChemicalFeatures &chemfe) {
  AtomMask applied_criteria(new_mask, ag_pointer, chemfe, cf, style);
  addAtoms(applied_criteria, cf, chemfe);
}

//-------------------------------------------------------------------------------------------------
void AtomMask::addAtoms(const AtomMask &new_mask, const CoordinateFrame &cf) {
  const ChemicalFeatures chemfe(ag_pointer, cf.data());
  addAtoms(new_mask, cf, chemfe);
}

//-------------------------------------------------------------------------------------------------
void AtomMask::addAtoms(const std::string &new_mask, const CoordinateFrame &cf) {
  const ChemicalFeatures chemfe(ag_pointer, cf.data());
  addAtoms(new_mask, cf, chemfe);
}

//-------------------------------------------------------------------------------------------------
double AtomMask::extractRangeValue(const int start_char, const int final_char, int *position) {

  // Determine the last character that could be part of the number.
  int limit = start_char;
  while (limit < final_char &&
         input_text[limit] != '(' && input_text[limit] == ')' && input_text[limit] == '[' &&
         input_text[limit] != ']' && input_text[limit] == '{' && input_text[limit] == '}' &&
         input_text[limit] != '&' && input_text[limit] == '|' && input_text[limit] == '!') {
    limit++;
  }
  if (limit == start_char) {
    rtErr("No range could be found after character " + std::to_string(start_char) +
          " in atom mask \"" + input_text + "\".", "AtomMask", "extractRangeValue");
  }
  const int rlen = limit = start_char;
  std::string buffer(rlen, ' ');
  for (int i = start_char; i < limit; i++) {
    buffer[i - start_char] = input_text[i];
  }
  if (verifyNumberFormat(input_text.data(), NumberFormat::STANDARD_REAL, start_char, rlen) ||
      verifyNumberFormat(input_text.data(), NumberFormat::SCIENTIFIC, start_char, rlen)) {

    // Advance the counter so that the number does not become part of subsequent interpretations
    *position = limit;
    return stod(buffer);
  }
  else {
    rtErr("No range could be found in segment \"" + buffer + "\" of atom mask \"" + input_text +
          "\".", "AtomMask", "extractRangeValue");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<SelectionItem> AtomMask::evaluateInclusions(const std::string &inclusion_list) {
  const int incl_len = inclusion_list.size();
  
  // Determine the kind of selection items.  There is some ambiguity with ELEMENT_NAME, ATOM_NAME,
  // and RESIDUE_NAME.  In what follows, the selections could tain numbered elements, atoms, and
  // residues, respectively.  Only ATOM_TYPE is unambiguous at this stage.  ??_NAME will be
  // taken to mean ??_NAME or ??_NUMBER until the items can be parsed further.
  SelectionItemKind init_kind, final_kind;
  int start_idx;
  if (inclusion_list[0] == '@') {
    if (incl_len > 2 && inclusion_list[1] == '%') {
      init_kind = SelectionItemKind::ATOM_TYPE;
      start_idx = 2;
    }
    else if (incl_len > 2 && inclusion_list[1] == '/') {
      init_kind = SelectionItemKind::ELEMENT_NAME;
      start_idx = 2;
    }
    else {
      init_kind = SelectionItemKind::ATOM_NAME;
      start_idx = 1;
    }
  }
  else if (inclusion_list[0] == ':') {
    init_kind = SelectionItemKind::RESIDUE_NAME;
    start_idx = 1;
  }
  
  // Separate the text based on commas.  Remove all whitespace from the initial string.
  std::string packed_list = inclusion_list;
  int j = 0;
  for (int i = 0; i < incl_len; i++) {
    packed_list[j] = packed_list[i];
    j += (inclusion_list[i] != ' ');
  }
  const int packed_incl_len = j;
  packed_list.resize(packed_incl_len);
  const std::vector<std::string> delimiters = { std::string(",") };
  const std::vector<std::string> dashes = { std::string(",") };
  std::vector<std::string> raw_items = separateText(packed_list.substr(start_idx, incl_len),
                                                    delimiters);
  const int n_items = raw_items.size();
  std::vector<SelectionItem> result;
  for (int i = 0; i < n_items; i++) {

    // Is this a number or a number range?
    const char tc = raw_items[i][0];
    const int ri_len = raw_items[i].size(); 
    if ((tc >= '0' && tc <= '9') || tc == '-') {

      // Get the final SelectionItem kind
      switch (init_kind) {
      case SelectionItemKind::ATOM_TYPE:
        rtErr("Inclusion list " + inclusion_list + " of atom mask " + input_text + " leads to an "
              "atom type name beginning with a digit in entry " + std::to_string(i + 1) + ".",
              "AtomMask", "evaluateInclusions");
      case SelectionItemKind::ATOM_NAME:
        final_kind = SelectionItemKind::ATOM_NUMBER;
        break;
      case SelectionItemKind::ELEMENT_NAME:
        final_kind = SelectionItemKind::ELEMENT_NUMBER;
        break;
      case SelectionItemKind::RESIDUE_NAME:
        final_kind = SelectionItemKind::RESIDUE_NUMBER;
        break;
      case SelectionItemKind::ATOM_NUMBER:
      case SelectionItemKind::ELEMENT_NUMBER:
      case SelectionItemKind::RESIDUE_NUMBER:
      case SelectionItemKind::NONE:
        break;
      }

      // Parse the number or number series.  The first character may be a dash, as the residue
      // and atom number search will run based on atom and residue indices from some input file
      // (encoded in special arrays in the topology) rather than the C-array indexing.
      int value_a, value_b;
      int break_char = 0;
      for (int j = 1; j < ri_len; j++) {
        if (raw_items[i][j] == '-' && raw_items[i][j - 1] >= '0' && raw_items[i][j - 1] <= '9') {
          break_char = j;
          break;
        }
      }
      bool problem = false;
      if (break_char == 0) {
        
        // There is only one number.  Read it.
        if (verifyNumberFormat(raw_items[i].c_str(), NumberFormat::INTEGER)) {
          value_a = stol(raw_items[i]);
          value_b = value_a;
        }
        else {
          problem = true;
        }
      }
      else if (break_char == ri_len - 1) {
        problem = true;
      }
      else {

        // There are two numbers.  Read them.
        if (verifyNumberFormat(raw_items[i].c_str(), NumberFormat::INTEGER, 0, break_char)) {
          std::string tmp_number;
          for (int j = 0; j < break_char; j++) {
            tmp_number += raw_items[i][j];
          }
          value_a = stol(tmp_number);
        }
        else {
          problem = true;
        }
        if (verifyNumberFormat(raw_items[i].c_str(), NumberFormat::INTEGER, break_char + 1,
                               ri_len - break_char - 1)) {
          std::string tmp_number;
          for (int j = break_char + 1; j < ri_len; j++) {
            tmp_number += raw_items[i][j];
          }
          value_b = stol(tmp_number);
        }
        else {
          problem = true;
        }
      }
      if (problem) {
        rtErr("Invalid numerical expression \"" + raw_items[i] + "\" in item " +
              std::to_string(i + 1) + " of list " + inclusion_list + " in atom mask " +
              input_text + ".", "AtomMask", "evaluateInclusions");
      }
      result.push_back(SelectionItem(final_kind, value_a, value_b));
    }
    else {

      // Even if the name is longer than four characters, it may contain escape sequences that
      // protect certain characters from being interpreted as wildcards.  The only accepted
      // escape sequence for an atom, atom type, or residue name is "\".  Determine which
      // characters are wildcards.  Note that only names of four characters or less may contain
      // wildcards.
      std::vector<WildCardKind> wildcards(std::max(ri_len, 4), WildCardKind::NONE);
      for (size_t j = 0; j < raw_items[i].size(); j++) {
        const char tc = raw_items[i][j];
        if (tc == '\\') {
          raw_items[i].erase(j);
        }
        else {
          if (tc == '*' || tc == '=') {
            wildcards[j] = WildCardKind::FREE_STRETCH;
          }
          else if (tc == '?') {
            wildcards[j] = WildCardKind::FREE_CHARACTER;
          }
        }
      }

      // The given name may still be longer than four characters, and the topology may contain
      // codifications of names that are longer than four characters.  Address this possibility.
      final_kind = init_kind;
      if (raw_items[i].size() > 4) {
        std::vector<char4> name_codes;
        switch (final_kind) {
        case SelectionItemKind::ATOM_NAME:
          name_codes = ag_pointer->matchOverflowAtomName(raw_items[i]);
          break;
        case SelectionItemKind::ELEMENT_NAME:
          rtErr("Invalid element name " + raw_items[i] + ".  No element can have a name longer "
                "than two characters.", "AtomMask", "evaluateInclusions");
        case SelectionItemKind::ATOM_TYPE:
          name_codes = ag_pointer->matchOverflowAtomType(raw_items[i]);
          break;
        case SelectionItemKind::RESIDUE_NAME:
          name_codes = ag_pointer->matchOverflowResidueName(raw_items[i]);
          break;
        case SelectionItemKind::ATOM_NUMBER:
        case SelectionItemKind::ELEMENT_NUMBER:
        case SelectionItemKind::RESIDUE_NUMBER:
        case SelectionItemKind::NONE:
          break;
        }

        // No wildcards in a codified name, but there may be many codified names matching the
        // one query string due to wildcards
        const int n_codes = name_codes.size();
        for (int j = 0; j < 4; j++) {
          wildcards[j] = WildCardKind::NONE;
        }
        for (int j = 0; j < n_codes; j++) {
          result.push_back(SelectionItem(final_kind, name_codes[j], wildcards));
        }
      }
      else {
        result.push_back(SelectionItem(final_kind, stringToChar4(raw_items[i]), wildcards));
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> AtomMask::parseMask(const std::vector<int> &scope_levels, int *position,
                                      const CoordinateFrameReader &cfr,
                                      const ChemicalFeatures &chemfe) {

  // Bail out immediately if the mask contains no information
  if (scope_levels.size() == 0LLU) {
    const int uint_bits = sizeof(uint) * 8;
    return std::vector<uint>((cfr.natom + uint_bits - 1) / uint_bits, 0U);
  }

  // For a given starting position, determine the scope level.  This function will make recursive
  // calls to itself if it sees the scope level increasing, and return its accumulated result if
  // it sees the scope level decrease.
  const int start_position = *position;
  const int operating_level = scope_levels[start_position];
  const int n_char = input_text.size();

  // Find the extent of the current scope
  int i = start_position;
  while (i < n_char) {
    if (scope_levels[i] < operating_level) {
      break;
    }
    i++;
  }
  const int end_position = i;
  
  // Accumulate a series of masks and operations within the current scope
  std::vector<MaskComponent> scope_cmp;
  i = start_position;
  while (i < end_position) {
    
    // Identify operators
    if ((input_text[i] == '<' || input_text[i] == '>') && i < end_position - 1 &&
        (input_text[i + 1] == '@' || input_text[i + 1] == ':')) {

      // Assign a "<" or ">" ranged operator
      const int op_init_pos = i;
      if (input_text[i] == '<') {
        if (input_text[i + 1] == '@') {
          const double range = extractRangeValue(i + 2, end_position, &i);
          scope_cmp.push_back(MaskComponent(MaskOperator::ATOM_LT, range, "<@ (atoms within " +
                                            realToString(range) + " Angstrom)", op_init_pos, i));
        }
        else {
          const double range = extractRangeValue(i + 2, end_position, &i);
          scope_cmp.push_back(MaskComponent(MaskOperator::RESIDUE_LT, range, "<: (residues "
                                            "such that at all atoms lie within " +
                                            realToString(range) + " Angstrom)", op_init_pos, i));
        }
      }
      else {
        if (input_text[i + 1] == '@') {
          const double range = extractRangeValue(i + 2, end_position, &i);
          scope_cmp.push_back(MaskComponent(MaskOperator::ATOM_GT, range, ">@ (atoms beyond " +
                                            realToString(range) + "A)", op_init_pos, i));
        }
        else {
          const double range = extractRangeValue(i + 2, end_position, &i);
          scope_cmp.push_back(MaskComponent(MaskOperator::RESIDUE_GT, range, ">: (residues "
                                            "such that no atom lies within " +
                                            realToString(range) + " Angstrom)", op_init_pos, i));
        }
      }
      continue;
    }
    else if ((input_text[i] == '<' || input_text[i] == '>') && i < end_position - 2 &&
             input_text[i + 1] == '=' && (input_text[i + 2] == '@' || input_text[i + 2] == ':')) {
      
      // Assign a "<=" or "=>" ranged operator
      const int op_init_pos = i;
      if (input_text[i] == '<') {
        if (input_text[i + 2] == '@') {
          const double range = extractRangeValue(i + 3, end_position, &i);
          scope_cmp.push_back(MaskComponent(MaskOperator::ATOM_LE, range, "<@ (atoms at or "
                                            "within " + realToString(range) + " Angstrom)",
                                            op_init_pos, i));
        }
        else {
          const double range = extractRangeValue(i + 3, end_position, &i);
          scope_cmp.push_back(MaskComponent(MaskOperator::RESIDUE_LE, range, "<: (residues "
                                            "such that at all atoms lie at or within " +
                                            realToString(range) + " Angstrom)", op_init_pos, i));
        }
      }
      else {
        if (input_text[i + 2] == '@') {
          const double range = extractRangeValue(i + 3, end_position, &i);
          scope_cmp.push_back(MaskComponent(MaskOperator::ATOM_GE, range, ">@ (atoms at or "
                                            "beyond " + realToString(range) + "Angstrom)",
                                            op_init_pos, i));
        }
        else {
          const double range = extractRangeValue(i + 3, end_position, &i);
          scope_cmp.push_back(MaskComponent(MaskOperator::RESIDUE_GE, range, ">: (residues "
                                            "such that any atom lies at or beyond " +
                                            realToString(range) + " Angstrom)", op_init_pos, i));
        }
      }
      continue;
    }
    else if (input_text[i] == '&') {
      scope_cmp.push_back(MaskComponent(MaskOperator::AND, "& (AND)", i, i + 1));
      i++;
      continue;
    }
    else if (input_text[i] == '|') {
      scope_cmp.push_back(MaskComponent(MaskOperator::OR, "| (OR)", i, i + 1));
      i++;
      continue;
    }
    else if (input_text[i] == '!') {
      scope_cmp.push_back(MaskComponent(MaskOperator::NOT, "! (NOT)", i, i + 1));
      i++;
      continue;
    }

    // Identify groups of atoms and residues
    else if (input_text[i] == ':' || input_text[i] == '@') {
      
      // Associate all residues or atoms in the subsequent list
      int j = i;
      while (j < end_position &&
             input_text[j] != '!' && input_text[j] != '&' && input_text[j] != '|' &&
             input_text[j] != '(' && input_text[j] != ')' && input_text[j] != '{' &&
             input_text[j] != '}' && input_text[j] != '[' && input_text[j] != ']') {
        j++;
      }
      if (j == i + 1) {
        rtErr("Selection initiator " + std::to_string(input_text[i]) + " with no inclusion list "
              "at position " + std::to_string(i + 1) + " of ambmask " + input_text + ".",
              "AtomMask");
      }
      std::string incl_list = input_text.substr(i, j - i);
      
      // Evaluate this inclusion list and add it to the scope components
      scope_cmp.push_back(MaskComponent(evaluateInclusions(incl_list), ag_pointer, chemfe,
                                        incl_list, i, j));

      // Advance the overall counter to reflect the progress of reading this inclusion list
      i = j;
      continue;
    }

    // Recursive call to handle higher level scopes
    else if (input_text[i] == '(' || input_text[i] == '{' || input_text[i] == '[') {
      *position = i + 1;
      const std::vector<uint> subscope_mask = parseMask(scope_levels, position, cfr, chemfe);
      scope_cmp.push_back(MaskComponent(subscope_mask,
                                        input_text.substr(i, *position - i), i, *position));
      i = *position;
    }

    // Keep going
    else {
      i++;
    }
  }

  // Ensure that there is at least one mask in the list of tasks
  bool mask_found = false;
  for (size_t i = 0; i < scope_cmp.size(); i++) {
    mask_found = (mask_found || scope_cmp[i].getKind() == MaskComponentKind::MASK);
  }
  if (mask_found == false) {

    // Different errors for an entire mask string versus a specific scope of it
    if (start_position == 0 && end_position == static_cast<int>(input_text.size())) {
      rtErr("Atom mask string " + input_text + " contains nothing which can be parsed as a mask.",
            "parseMask");
    }
    else {
      rtErr("Sub-scope " + input_text.substr(start_position, end_position - start_position) +
            " beginning at position " + std::to_string(start_position) + " within atom mask "
            "string " + input_text + " contains nothing which can be parsed as a mask.",
            "parseMask");
    }
  }
  
  // Evaluate the accumulated list of masks and operations--there are no scopes anymore, just an
  // order of binding in the "&", "|", and "!" operations.  Begin with ! operations--apply them to
  // primitive masks directly after them in the list.
  const std::string mask_explained("either from a collection of atoms or the result of a prior "
                                   "evaluation enclosed by [], (), and {}.");
  for (size_t i = 0; i < scope_cmp.size(); i++) {
    if (scope_cmp[i].getKind() == MaskComponentKind::OPERATOR &&
        scope_cmp[i].getOperator() == MaskOperator::NOT) {
      if (i == scope_cmp.size() - 1 || scope_cmp[i + 1].getKind() != MaskComponentKind::MASK) {
        rtErr("Operator " + scope_cmp[i].getText() + " at position " +
              std::to_string(scope_cmp[i].getTextLimits().x) + " of mask string " + input_text +
              " must be followed by a mask, " + mask_explained, "parseMask");
      }

      // Reverse the bits in the mask immediately following the "!" operator
      scope_cmp[i + 1].applyNotOperator();
      scope_cmp.erase(scope_cmp.begin() + i);
    }
  }
  
  // Process AND operators.
  for (size_t i = 0; i < scope_cmp.size(); i++) {
    if (scope_cmp[i].getKind() == MaskComponentKind::OPERATOR &&
        scope_cmp[i].getOperator() == MaskOperator::AND) {
      if (i == 0 || i == scope_cmp.size() - 1 ||
          scope_cmp[i + 1].getKind() != MaskComponentKind::MASK ||
          scope_cmp[i - 1].getKind() != MaskComponentKind::MASK) {
        rtErr("Operator " + scope_cmp[i].getText() + " at position " +
              std::to_string(scope_cmp[i].getTextLimits().x) + " of mask string " + input_text +
              " must be preceded and followed by masks, " + mask_explained, "parseMask");
      }

      // Apply the AND operation between two masks
      scope_cmp[i + 1].applyAndOperator(scope_cmp[i - 1].getMask());
      scope_cmp.erase(scope_cmp.begin() + i - 1, scope_cmp.begin() + i + 1);
      i--;
    }
  }

  // Process OR operators.  They bind before AND.
  for (size_t i = 0; i < scope_cmp.size(); i++) {
    if (scope_cmp[i].getKind() == MaskComponentKind::OPERATOR &&
        scope_cmp[i].getOperator() == MaskOperator::OR) {
      if (i == 0 || i == scope_cmp.size() - 1 ||
          scope_cmp[i + 1].getKind() != MaskComponentKind::MASK ||
          scope_cmp[i - 1].getKind() != MaskComponentKind::MASK) {
        rtErr("Operator " + scope_cmp[i].getText() + " at position " +
              std::to_string(scope_cmp[i].getTextLimits().x) + " of mask string " + input_text +
              " must be preceded and followed by masks, " + mask_explained, "parseMask");
      }

      // Apply the OR operation between two masks
      scope_cmp[i + 1].applyOrOperator(scope_cmp[i - 1].getMask());
      scope_cmp.erase(scope_cmp.begin() + i - 1, scope_cmp.begin() + i + 1);
      i--;
    }
  }
  
  // Process ranged operators.  Each of them should be the last component of their respective
  // scopes, following a single primitive mask.
  for (size_t i = 0; i < scope_cmp.size(); i++) {
    if (scope_cmp[i].getKind() != MaskComponentKind::OPERATOR) {
      continue;
    }
    MaskOperator thisop = scope_cmp[i].getOperator();
    if (thisop == MaskOperator::ATOM_LT || thisop == MaskOperator::ATOM_GT ||
        thisop == MaskOperator::ATOM_LE || thisop == MaskOperator::ATOM_GE ||
        thisop == MaskOperator::RESIDUE_LT || thisop == MaskOperator::RESIDUE_GT ||
        thisop == MaskOperator::RESIDUE_LE || thisop == MaskOperator::RESIDUE_GE) {

      // The ranged operator should be the last component following a single mask by this point
      if (i != 1 || scope_cmp[0].getKind() != MaskComponentKind::MASK || scope_cmp.size() != 2) {
        rtErr("A ranged operation of mask string " + scope_cmp[i].getText() + " at position " +
              std::to_string(scope_cmp[i].getTextLimits().x) + " must ultimately be the final "
              "component of its scope, following a single evaluated mask expression for all other "
              "components.", "parseMask");
      }

      // Form a mask out of the range operator and its preceding MASK-kind MaskComponent, push
      // that to the back of the list, then pop the first two elements.
      scope_cmp.push_back(scope_cmp[i].applyRangeOperator(scope_cmp[0].getMask(), ag_pointer,
                                                          cfr));
      scope_cmp.erase(scope_cmp.begin(), scope_cmp.begin() + 2);
    }
  }
  
  // There should be one MaskComponent remaining at this point, and it should be MASK-kind
  if (scope_cmp.size() != 1 || scope_cmp[0].getKind() != MaskComponentKind::MASK) {
    std::string unparsed_list;
    const size_t ncmp = scope_cmp.size();
    bool first_mask_found = false;
    for (size_t i = 0; i < ncmp; i++) {
      switch (scope_cmp[i].getKind()) {
      case MaskComponentKind::OPERATOR:
        unparsed_list += "op ";
        break;
      case MaskComponentKind::MASK:
        unparsed_list += "mask ";
        break;
      }
      const int2 s_limits = scope_cmp[i].getTextLimits();
      unparsed_list += input_text.substr(s_limits.x, s_limits.y - s_limits.x);
      if (i < ncmp - 1) {
        unparsed_list += ", ";
      }
    }
    rtErr("Mask string " + input_text + " contained unparsed components " + unparsed_list + ".",
          "parseMask");
  }

  // Set the position counter to the end of the scope
  *position = end_position;
  
  // Return the last mask component's primitive mask
  return scope_cmp[0].getMask();
}

} // namespace chemistry
} // namespace stormm
