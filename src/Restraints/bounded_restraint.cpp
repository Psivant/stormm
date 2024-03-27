#include "copyright.h"
#include "Chemistry/atommask.h"
#include "Constants/symbol_values.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "Structure/local_arrangement.h"
#include "bounded_restraint.h"
#include "restraint_util.h"

namespace stormm {
namespace restraints {

using chemistry::AtomMask;
using chemistry::MaskInputMode;
using parse::realToString;
using parse::char4ToString;
using parse::NumberFormat;
using structure::imageValue;
using structure::ImagingMethod;
using symbols::pi;
using symbols::twopi;

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const AtomGraph *ag_in) :
    atom_i{-1}, atom_j{-1}, atom_k{-1}, atom_l{-1}, order{0},
    initial_step{0},
    final_step{0},
    initial_keq{0.0, 0.0},
    initial_r{0.0, 0.0, 0.0, 0.0},
    final_keq{0.0, 0.0},
    final_r{0.0, 0.0, 0.0, 0.0},
    initial_center{0.0, 0.0, 0.0},
    final_center{0.0, 0.0, 0.0},
    ag_pointer{ag_in}
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const std::string &mask_i_in, const std::string &mask_j_in,
                                   const std::string &mask_k_in, const std::string &mask_l_in,
                                   const AtomGraph *ag_in, const ChemicalFeatures &chemfe,
                                   const CoordinateFrameReader &cfr, const int init_step_in,
                                   const int final_step_in, const double init_k2_in,
                                   const double init_k3_in, const double init_r1_in,
                                   const double init_r2_in, const double init_r3_in,
                                   const double init_r4_in, const double final_k2_in,
                                   const double final_k3_in, const double final_r1_in,
                                   const double final_r2_in, const double final_r3_in,
                                   const double final_r4_in, const double3 init_ref_crd_in,
                                   const double3 final_ref_crd_in) :
    atom_i{-1}, atom_j{-1}, atom_k{-1}, atom_l{-1}, order{0},
    initial_step{init_step_in},
    final_step{final_step_in},
    initial_keq{init_k2_in, init_k3_in},
    initial_r{init_r1_in, init_r2_in, init_r3_in, init_r4_in},
    final_keq{final_k2_in, final_k3_in},
    final_r{final_r1_in, final_r2_in, final_r3_in, final_r4_in},
    initial_center{init_ref_crd_in},
    final_center{final_ref_crd_in},
    ag_pointer{ag_in}
{
  // At least one atom is required
  if (mask_i_in.size() == 0LLU) {
    rtErr("At least one atom must be be specified in order to create a bounded restraint.",
          "BoundedRestraint");
  }

  // Obtain the atoms to which this restraint applies based on the input atom masks.
  std::vector<AtomMask> all_masks = {
    AtomMask(mask_i_in, ag_pointer, chemfe, cfr, MaskInputMode::AMBMASK,
             "Bounded restraint atom I mask"),
    AtomMask(mask_j_in, ag_pointer, chemfe, cfr, MaskInputMode::AMBMASK,
             "Bounded restraint atom J mask"),
    AtomMask(mask_k_in, ag_pointer, chemfe, cfr, MaskInputMode::AMBMASK,
             "Bounded restraint atom K mask"),
    AtomMask(mask_l_in, ag_pointer, chemfe, cfr, MaskInputMode::AMBMASK,
             "Bounded restraint atom : mask") };
  std::vector<int> n_mask_atom(4);
  bool running = true;
  for (int i = 0; i < 4; i++) {
    n_mask_atom[i] = all_masks[i].getMaskedAtomCount();
    if (n_mask_atom[i] > 1 || (all_masks[i].getInputText().size() > 0 && n_mask_atom[i] == 0)) {
      rtErr("Atom mask \"" + all_masks[i].getInputText() + "\" specifies " +
            std::to_string(n_mask_atom[i]) + " atoms in topology " + ag_pointer->getFileName() +
            ".  The mask must specify one and only one atom.", "BoundedRestraint");
    }
    running = (running && (n_mask_atom[i] == 1));
    order += running;
  }

  // Report if the first atom mask is invalid
  if (n_mask_atom[0] == 1) {
    atom_i = all_masks[0].getMaskedAtomList()[0];
  }
  else {
    rtErr("Atom mask \"" + mask_i_in + "\" specifies " + std::to_string(n_mask_atom[0]) +
          " atoms in topology " + ag_pointer->getFileName() + ".  The mask must specify one and "
          "only one atom.", "BoundedRestraint");
  }

  // Set other atoms
  if (order > 1) {
    atom_j = n_mask_atom[1];
  }
  if (order > 2) {
    atom_k = n_mask_atom[2];
  }
  if (order > 3) {
    atom_l = n_mask_atom[3];
  }

  // Correct the displacment values if necessary
  checkDisplacementLimits(&initial_r);
  checkDisplacementLimits(&final_r);
}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const std::string &mask_i_in, const std::string &mask_j_in,
                                   const std::string &mask_k_in, const AtomGraph *ag_in,
                                   const ChemicalFeatures &chemfe,
                                   const CoordinateFrameReader &cfr, const int init_step_in,
                                   const int final_step_in, const double init_k2_in,
                                   const double init_k3_in, const double init_r1_in,
                                   const double init_r2_in, const double init_r3_in,
                                   const double init_r4_in, const double final_k2_in,
                                   const double final_k3_in, const double final_r1_in,
                                   const double final_r2_in, const double final_r3_in,
                                   const double final_r4_in) :
    BoundedRestraint(mask_i_in, mask_j_in, mask_k_in, std::string(""), ag_in, chemfe, cfr,
                     init_step_in, final_step_in, init_k2_in, init_k3_in, init_r1_in, init_r2_in,
                     init_r3_in, init_r4_in, final_k2_in, final_k3_in, final_r1_in, final_r2_in,
                     final_r3_in, final_r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const std::string &mask_i_in, const std::string &mask_j_in,
                                   const AtomGraph *ag_in, const ChemicalFeatures &chemfe,
                                   const CoordinateFrameReader &cfr, const int init_step_in,
                                   const int final_step_in, const double init_k2_in,
                                   const double init_k3_in, const double init_r1_in,
                                   const double init_r2_in, const double init_r3_in,
                                   const double init_r4_in, const double final_k2_in,
                                   const double final_k3_in, const double final_r1_in,
                                   const double final_r2_in, const double final_r3_in,
                                   const double final_r4_in) :
    BoundedRestraint(mask_i_in, mask_j_in, std::string(""), std::string(""), ag_in, chemfe, cfr,
                     init_step_in, final_step_in, init_k2_in, init_k3_in, init_r1_in, init_r2_in,
                     init_r3_in, init_r4_in, final_k2_in, final_k3_in, final_r1_in, final_r2_in,
                     final_r3_in, final_r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const std::string &mask_i_in, const AtomGraph *ag_in,
                                   const ChemicalFeatures &chemfe,
                                   const CoordinateFrameReader &cfr, const int init_step_in,
                                   const int final_step_in, const double init_k2_in,
                                   const double init_k3_in, const double init_r1_in,
                                   const double init_r2_in, const double init_r3_in,
                                   const double init_r4_in, const double final_k2_in,
                                   const double final_k3_in, const double final_r1_in,
                                   const double final_r2_in, const double final_r3_in,
                                   const double final_r4_in, const double3 init_ref_crd_in,
                                   const double3 final_ref_crd_in) :
    BoundedRestraint(mask_i_in, std::string(""), std::string(""), std::string(""), ag_in, chemfe,
                     cfr, init_step_in, final_step_in, init_k2_in, init_k3_in, init_r1_in,
                     init_r2_in, init_r3_in, init_r4_in, final_k2_in, final_k3_in, final_r1_in,
                     final_r2_in, final_r3_in, final_r4_in, init_ref_crd_in, final_ref_crd_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const std::string &mask_i_in, const std::string &mask_j_in,
                                   const std::string &mask_k_in, const std::string &mask_l_in,
                                   const AtomGraph *ag_in, const ChemicalFeatures &chemfe,
                                   const CoordinateFrameReader &cfr, const double k2_in,
                                   const double k3_in, const double r1_in, const double r2_in,
                                   const double r3_in, const double r4_in) :
    BoundedRestraint(mask_i_in, mask_j_in, mask_k_in, mask_l_in, ag_in, chemfe, cfr, 0, 0,
                     k2_in, k3_in, r1_in, r2_in, r3_in, r4_in, k2_in, k3_in, r1_in, r2_in, r3_in,
                     r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const std::string &mask_i_in, const std::string &mask_j_in,
                                   const std::string &mask_k_in, const AtomGraph *ag_in,
                                   const ChemicalFeatures &chemfe,
                                   const CoordinateFrameReader &cfr, const double k2_in,
                                   const double k3_in, const double r1_in, const double r2_in,
                                   const double r3_in, const double r4_in) :
    BoundedRestraint(mask_i_in, mask_j_in, mask_k_in, std::string(""), ag_in, chemfe, cfr, 0, 0,
                     k2_in, k3_in, r1_in, r2_in, r3_in, r4_in, k2_in, k3_in, r1_in, r2_in, r3_in,
                     r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const std::string &mask_i_in, const std::string &mask_j_in,
                                   const AtomGraph *ag_in, const ChemicalFeatures &chemfe,
                                   const CoordinateFrameReader &cfr, const double k2_in,
                                   const double k3_in, const double r1_in, const double r2_in,
                                   const double r3_in, const double r4_in) :
    BoundedRestraint(mask_i_in, mask_j_in, std::string(""), std::string(""), ag_in, chemfe, cfr, 0,
                     0, k2_in, k3_in, r1_in, r2_in, r3_in, r4_in, k2_in, k3_in, r1_in, r2_in,
                     r3_in, r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const std::string &mask_i_in, const AtomGraph *ag_in,
                                   const ChemicalFeatures &chemfe,
                                   const CoordinateFrameReader &cfr, const double k2_in,
                                   const double k3_in, const double r1_in, const double r2_in,
                                   const double r3_in, const double r4_in,
                                   const double3 ref_crd_in) :
    BoundedRestraint(mask_i_in, std::string(""), std::string(""), std::string(""), ag_in, chemfe,
                     cfr, 0, 0, k2_in, k3_in, r1_in, r2_in, r3_in, r4_in, k2_in, k3_in, r1_in,
                     r2_in, r3_in, r4_in, ref_crd_in, ref_crd_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const int atom_i_in, const int atom_j_in, const int atom_k_in,
                                   const int atom_l_in, const AtomGraph *ag_in,
                                   const int init_step_in, const int final_step_in,
                                   const double init_k2_in, const double init_k3_in,
                                   const double init_r1_in, const double init_r2_in,
                                   const double init_r3_in, const double init_r4_in,
                                   const double final_k2_in, const double final_k3_in,
                                   const double final_r1_in, const double final_r2_in,
                                   const double final_r3_in, const double final_r4_in,
                                   const double3 init_ref_crd_in, const double3 final_ref_crd_in) :
    atom_i{atom_i_in}, atom_j{atom_j_in}, atom_k{atom_k_in}, atom_l{atom_l_in}, order{0},
    initial_step{init_step_in},
    final_step{final_step_in},
    initial_keq{init_k2_in, init_k3_in},
    initial_r{init_r1_in, init_r2_in, init_r3_in, init_r4_in},
    final_keq{final_k2_in, final_k3_in},
    final_r{final_r1_in, final_r2_in, final_r3_in, final_r4_in},
    initial_center{init_ref_crd_in},
    final_center{final_ref_crd_in},
    ag_pointer{ag_in}
{
  // At least one atom is required
  if (atom_i < 0) {
    rtErr("At least one atom must be be specified in order to create a bounded restraint.",
          "BoundedRestraint");
  }
  
  // Get the order of the restraint
  order = (atom_i >= 0) + (atom_j >= 0) + (atom_k >= 0) + (atom_l >= 0);

  // Confirm that all atom selections are valid
  const int natom = ag_pointer->getAtomCount();
  if (atom_i >= natom || (order > 1 && atom_j >= natom) || (order > 2 && atom_k >= natom) ||
      (order > 3 && atom_l >= natom)) {
    std::string atom_number_list = std::to_string(atom_i);
    if (order > 1) {
      atom_number_list += ", " + std::to_string(atom_j);
    }
    if (order > 2) {
      atom_number_list += ", " + std::to_string(atom_k);
    }
    if (order > 3) {
      atom_number_list += ", " + std::to_string(atom_l);
    }
    rtErr("The atom selection " + atom_number_list + " does not fall within the atom indexing of "
          "topology " + ag_pointer->getFileName() + ".", "BoundedRestraint");
  }

  // Correct the displacment values if necessary
  checkDisplacementLimits(&initial_r);
  checkDisplacementLimits(&final_r);
}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const int atom_i_in, const int atom_j_in, const int atom_k_in,
                                   const AtomGraph *ag_in, const int init_step_in,
                                   const int final_step_in, const double init_k2_in,
                                   const double init_k3_in, const double init_r1_in,
                                   const double init_r2_in, const double init_r3_in,
                                   const double init_r4_in, const double final_k2_in,
                                   const double final_k3_in, const double final_r1_in,
                                   const double final_r2_in, const double final_r3_in,
                                   const double final_r4_in) :
    BoundedRestraint(atom_i_in, atom_j_in, atom_k_in, -1, ag_in, init_step_in, final_step_in,
                     init_k2_in, init_k3_in, init_r1_in, init_r2_in, init_r3_in, init_r4_in,
                     final_k2_in, final_k3_in, final_r1_in, final_r2_in, final_r3_in, final_r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const int atom_i_in, const int atom_j_in,
                                   const AtomGraph *ag_in, const int init_step_in,
                                   const int final_step_in, const double init_k2_in,
                                   const double init_k3_in, const double init_r1_in,
                                   const double init_r2_in, const double init_r3_in,
                                   const double init_r4_in, const double final_k2_in,
                                   const double final_k3_in, const double final_r1_in,
                                   const double final_r2_in, const double final_r3_in,
                                   const double final_r4_in) :
    BoundedRestraint(atom_i_in, atom_j_in, -1, -1, ag_in, init_step_in, final_step_in,
                     init_k2_in, init_k3_in, init_r1_in, init_r2_in, init_r3_in, init_r4_in,
                     final_k2_in, final_k3_in, final_r1_in, final_r2_in, final_r3_in, final_r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const int atom_i_in, const AtomGraph *ag_in,
                                   const int init_step_in, const int final_step_in,
                                   const double init_k2_in, const double init_k3_in,
                                   const double init_r1_in, const double init_r2_in,
                                   const double init_r3_in, const double init_r4_in,
                                   const double final_k2_in, const double final_k3_in,
                                   const double final_r1_in, const double final_r2_in,
                                   const double final_r3_in, const double final_r4_in,
                                   const double3 init_ref_crd_in, const double3 final_ref_crd_in) :
    BoundedRestraint(atom_i_in, -1, -1, -1, ag_in, init_step_in, final_step_in, init_k2_in,
                     init_k3_in, init_r1_in, init_r2_in, init_r3_in, init_r4_in, final_k2_in,
                     final_k3_in, final_r1_in, final_r2_in, final_r3_in, final_r4_in,
                     init_ref_crd_in, final_ref_crd_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const int atom_i_in, const int atom_j_in, const int atom_k_in,
                                   const int atom_l_in, const AtomGraph *ag_in, const double k2_in,
                                   const double k3_in, const double r1_in, const double r2_in,
                                   const double r3_in, const double r4_in) :
    BoundedRestraint(atom_i_in, atom_j_in, atom_k_in, atom_l_in, ag_in, 0, 0, k2_in, k3_in, r1_in,
                     r2_in, r3_in, r4_in, k2_in, k3_in, r1_in, r2_in, r3_in, r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const int atom_i_in, const int atom_j_in, const int atom_k_in,
                                   const AtomGraph *ag_in, const double k2_in, const double k3_in,
                                   const double r1_in, const double r2_in, const double r3_in,
                                   const double r4_in) :
    BoundedRestraint(atom_i_in, atom_j_in, atom_k_in, -1, ag_in, 0, 0, k2_in, k3_in, r1_in,
                     r2_in, r3_in, r4_in, k2_in, k3_in, r1_in, r2_in, r3_in, r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const int atom_i_in, const int atom_j_in,
                                   const AtomGraph *ag_in, const double k2_in, const double k3_in,
                                   const double r1_in, const double r2_in, const double r3_in,
                                   const double r4_in) :
    BoundedRestraint(atom_i_in, atom_j_in, -1, -1, ag_in, 0, 0, k2_in, k3_in, r1_in,
                     r2_in, r3_in, r4_in, k2_in, k3_in, r1_in, r2_in, r3_in, r4_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const int atom_i_in, const AtomGraph *ag_in, const double k2_in,
                                   const double k3_in, const double r1_in, const double r2_in,
                                   const double r3_in, const double r4_in, double3 ref_crd_in) :
    BoundedRestraint(atom_i_in, -1, -1, -1, ag_in, 0, 0, k2_in, k3_in, r1_in,
                     r2_in, r3_in, r4_in, k2_in, k3_in, r1_in, r2_in, r3_in, r4_in,
                     ref_crd_in, ref_crd_in)
{}

//-------------------------------------------------------------------------------------------------
BoundedRestraint::BoundedRestraint(const int atom_index, const AtomGraph *ag_in,
                                   const CoordinateFrameReader &cfr, const double k2_in,
                                   const double k3_in, const double r1_in, const double r2_in,
                                   const double r3_in, const double r4_in, const int refr_index) :
  BoundedRestraint(atom_index, -1, -1, -1, ag_in, 0, 0, k2_in, k3_in, r1_in, r2_in, r3_in, r4_in,
                   k2_in, k3_in, r1_in, r2_in, r3_in, r4_in, { 0.0, 0.0, 0.0 }, { 0.0, 0.0, 0.0 })
{
  // Make a bounds check on the input before setting the positional restraint center
  if (refr_index >= cfr.natom) {
    rtErr("The requested atom index " + std::to_string(refr_index) + " was invalid for a "
          "coordinate frame with " + std::to_string(cfr.natom) + " atoms.", "BoundedRestraint");
  }
  if (refr_index < 0) {

    // Default to restraining the atom to its known coordinates
    initial_center = { cfr.xcrd[atom_index], cfr.ycrd[atom_index], cfr.zcrd[atom_index] };
  }
  else {
    initial_center = { cfr.xcrd[refr_index], cfr.ycrd[refr_index], cfr.zcrd[refr_index] };
  }
  final_center = initial_center;
}

//-------------------------------------------------------------------------------------------------
int BoundedRestraint::getAtomIndex(const int restrained_atom_number) const {
  if (restrained_atom_number > order) {
    rtErr("A restraint of order " + std::to_string(order) + " cannot produce a valid atom index "
          "for participating atom " + std::to_string(restrained_atom_number) + ".",
          "BoundedRestraint", "getAtomIndex");
  }
  switch (restrained_atom_number) {
  case 1:
    return atom_i;
  case 2:
    return atom_j;
  case 3:
    return atom_k;
  case 4:
    return atom_l;
  default:
    rtErr("Valid atoms include 1, 2, 3, and 4, not " + std::to_string(restrained_atom_number) +
          ".", "BoundedRestraint", "getAtomIndex");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int BoundedRestraint::getOrder() const {
  return order;
}

//-------------------------------------------------------------------------------------------------
int BoundedRestraint::getInitialStep() const {
  return initial_step;
}

//-------------------------------------------------------------------------------------------------
int BoundedRestraint::getFinalStep() const {
  return final_step;
}

//-------------------------------------------------------------------------------------------------
double2 BoundedRestraint::getStiffness(const int step_number) const {
  const Vec2<double> mixwt = computeRestraintMixture<double>(step_number, initial_step,
                                                             final_step);
  double2 td_keq;
  td_keq.x = (mixwt.x * initial_keq.x) + (mixwt.y * final_keq.x);
  td_keq.y = (mixwt.x * initial_keq.y) + (mixwt.y * final_keq.y);
  return td_keq;
}

//-------------------------------------------------------------------------------------------------
double2 BoundedRestraint::getInitialStiffness() const {
  return initial_keq;
}

//-------------------------------------------------------------------------------------------------
double2 BoundedRestraint::getFinalStiffness() const {
  return final_keq;
}

//-------------------------------------------------------------------------------------------------
double4 BoundedRestraint::getDisplacements(const int step_number) const {
  const Vec2<double> mixwt = computeRestraintMixture<double>(step_number, initial_step,
                                                             final_step);
  double4 td_r;
  td_r.x = (mixwt.x * initial_r.x) + (mixwt.y * final_r.x);
  td_r.y = (mixwt.x * initial_r.y) + (mixwt.y * final_r.y);
  td_r.z = (mixwt.x * initial_r.z) + (mixwt.y * final_r.z);
  td_r.w = (mixwt.x * initial_r.w) + (mixwt.y * final_r.w);
  return td_r;
}

//-------------------------------------------------------------------------------------------------
double4 BoundedRestraint::getInitialDisplacements() const {
  return initial_r;
}

//-------------------------------------------------------------------------------------------------
double4 BoundedRestraint::getFinalDisplacements() const {
  return final_r;
}

//-------------------------------------------------------------------------------------------------
double3 BoundedRestraint::getTargetSite(const int step_number) const {

  // Check the order of the restraint and its time dependence
  if (order != 1) {
    rtErr("The target site is meaningless outside of positional restraints on individual atoms.",
          "BoundedRestraint", "getInitialTargetSite");
  }
  const Vec2<double> mixwt = computeRestraintMixture<double>(step_number, initial_step,
                                                             final_step);
  double3 td_center;
  td_center.x = (mixwt.x * initial_center.x) + (mixwt.y * final_center.x);
  td_center.y = (mixwt.x * initial_center.y) + (mixwt.y * final_center.y);
  td_center.z = (mixwt.x * initial_center.z) + (mixwt.y * final_center.z);
  return td_center;
}

//-------------------------------------------------------------------------------------------------
double3 BoundedRestraint::getInitialTargetSite() const {

  // Check the order of the restraint 
  if (order != 1) {
    rtErr("The target site is meaningless outside of positional restraints on individual atoms.",
          "BoundedRestraint", "getInitialTargetSite");
  }
  return initial_center;
}

//-------------------------------------------------------------------------------------------------
double3 BoundedRestraint::getFinalTargetSite() const {

  // Check the order of the restraint 
  if (order != 1) {
    rtErr("The target site is meaningless outside of positional restraints on individual atoms.",
          "BoundedRestraint", "getFinalTargetSite");
  }
  return final_center;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* BoundedRestraint::getTopologyPointer() const {
  return ag_pointer;
}

//-------------------------------------------------------------------------------------------------
void BoundedRestraint::setInitialStep(const int new_init_step) {
  initial_step = new_init_step;
}

//-------------------------------------------------------------------------------------------------
void BoundedRestraint::setFinalStep(const int new_final_step) {
  final_step = new_final_step;
}

//-------------------------------------------------------------------------------------------------
void BoundedRestraint::setStiffness(const double new_keq) {
  if (initial_step != final_step && initial_step != 0) {
    rtErr("The time-dependence of the restraint forbids use of static parameter setters.  Use "
          "set(Initial,Final)Stiffness() instead.", "BoundedRestraint", "setStiffness");
  }
  initial_keq = { new_keq, new_keq };
}

//-------------------------------------------------------------------------------------------------
void BoundedRestraint::setStiffnesses(const double new_k2, const double new_k3) {
  if (initial_step != final_step && initial_step != 0) {
    rtErr("The time-dependence of the restraint forbids use of static parameter setters.  Use "
          "set(Initial,Final)Stiffnesses() instead.", "BoundedRestraint", "setStiffnesses");
  }
  initial_keq = { new_k2, new_k3 };
}

//-------------------------------------------------------------------------------------------------
void BoundedRestraint::setInitialStiffness(const double new_init_keq) {
  if (initial_step == final_step && initial_step == 0) {
    rtErr("The time-independence of the restraint requires use of the dynamic parameter setter.  "
          "Use setStiffness() instead.", "BoundedRestraint", "setInitialStiffness");
  }
  initial_keq = { new_init_keq, new_init_keq };
}

//-------------------------------------------------------------------------------------------------
void BoundedRestraint::setInitialStiffnesses(const double new_init_k2, const double new_init_k3) {
  if (initial_step == final_step && initial_step == 0) {
    rtErr("The time-independence of the restraint requires use of the dynamic parameter setter.  "
          "Use setStiffnesses() instead.", "BoundedRestraint", "setInitialStiffnesses");
  }
  initial_keq = { new_init_k2, new_init_k3 };
}

//-------------------------------------------------------------------------------------------------
void BoundedRestraint::setFinalStiffness(const double new_final_keq) {
  if (initial_step == final_step && initial_step == 0) {
    rtErr("The time-independence of the restraint requires use of the dynamic parameter setter.  "
          "Use setStiffness() instead.", "BoundedRestraint", "setFinalStiffness");
  }
  final_keq = { new_final_keq, new_final_keq };
}

//-------------------------------------------------------------------------------------------------
void BoundedRestraint::setFinalStiffnesses(const double new_final_k2, const double new_final_k3) {
  if (initial_step == final_step && initial_step == 0) {
    rtErr("The time-independence of the restraint requires use of the dynamic parameter setter.  "
          "Use setStiffnesses() instead.", "BoundedRestraint", "setFinalStiffnesses");
  }
  final_keq = { new_final_k2, new_final_k3 };
}

//-------------------------------------------------------------------------------------------------
void BoundedRestraint::setDisplacements(const double new_r1, const double new_r2,
                                        const double new_r3, const double new_r4) {
  if (initial_step != final_step && initial_step != 0) {
    rtErr("The time-dependence of the restraint forbids use of static parameter setters.  Use "
          "set(Initial,Final)Displacements() instead.", "BoundedRestraint", "setDisplacments");
  }
  initial_r = { new_r1, new_r2, new_r3, new_r4 };
}

//-------------------------------------------------------------------------------------------------
void BoundedRestraint::setInitialDisplacements(const double new_r1, const double new_r2,
                                               const double new_r3, const double new_r4) {
  if (initial_step == final_step && initial_step == 0) {
    rtErr("The time-independence of the restraint requires use of static parameter setters.  Use "
          "setDisplacements() instead.", "BoundedRestraint", "setInitialDisplacments");
  }
  initial_r = { new_r1, new_r2, new_r3, new_r4 };
}

//-------------------------------------------------------------------------------------------------
void BoundedRestraint::setFinalDisplacements(const double new_r1, const double new_r2,
                                             const double new_r3, const double new_r4) {
  if (initial_step == final_step && initial_step == 0) {
    rtErr("The time-independence of the restraint requires use of static parameter setters.  Use "
          "setDisplacements() instead.", "BoundedRestraint", "setFinalDisplacments");
  }
  final_r = { new_r1, new_r2, new_r3, new_r4 };
}

//-------------------------------------------------------------------------------------------------
void BoundedRestraint::setTargetSite(const double new_ref_x, const double new_ref_y,
                                     const double new_ref_z) {
  if (order != 1) {
    rtErr("Target sites are only applicable to positional restraints.  This restraint order = " +
          std::to_string(order) + ".", "BoundedRestraint", "setTargetSite");
  }
  if (initial_step != final_step && initial_step != 0) {
    rtErr("The time-dependence of the restraint forbids use of static parameter setters.  Use "
          "set(Initial,Final)TargetSite() instead.", "BoundedRestraint", "setTargetSite");
  }
  initial_center = { new_ref_x, new_ref_y, new_ref_z };
}

//-------------------------------------------------------------------------------------------------
void BoundedRestraint::setTargetSite(const double3 new_ref_crd) {
  setTargetSite(new_ref_crd.x, new_ref_crd.y, new_ref_crd.z);
}

//-------------------------------------------------------------------------------------------------
void BoundedRestraint::setInitialTargetSite(const double new_ref_x, const double new_ref_y,
                                            const double new_ref_z) {
  if (order != 1) {
    rtErr("Target sites are only applicable to positional restraints.  This restraint order = " +
          std::to_string(order) + ".", "BoundedRestraint", "setInitialTargetSite");
  }
  if (initial_step == final_step && initial_step == 0) {
    rtErr("The time-independence of the restraint requires use of static parameter setters.  Use "
          "setTargetSite() instead.", "BoundedRestraint", "setInitialTargetSite");
  }
  initial_center = { new_ref_x, new_ref_y, new_ref_z };
}

//-------------------------------------------------------------------------------------------------
void BoundedRestraint::setInitialTargetSite(const double3 new_ref_crd) {
  setInitialTargetSite(new_ref_crd.x, new_ref_crd.y, new_ref_crd.z);
}

//-------------------------------------------------------------------------------------------------
void BoundedRestraint::setFinalTargetSite(const double new_ref_x, const double new_ref_y,
                                          const double new_ref_z) {
  if (order != 1) {
    rtErr("Target sites are only applicable to positional restraints.  This restraint order = " +
          std::to_string(order) + ".", "BoundedRestraint", "setFinalTargetSite");
  }
  if (initial_step == final_step && initial_step == 0) {
    rtErr("The time-independence of the restraint requires use of static parameter setters.  Use "
          "setTargetSite() instead.", "BoundedRestraint", "setFinalTargetSite");
  }
  final_center = { new_ref_x, new_ref_y, new_ref_z };
}

//-------------------------------------------------------------------------------------------------
void BoundedRestraint::setFinalTargetSite(const double3 new_ref_crd) {
  setFinalTargetSite(new_ref_crd.x, new_ref_crd.y, new_ref_crd.z);
}

//-------------------------------------------------------------------------------------------------
std::string BoundedRestraint::reportAtomList() {
  std::string result;
  for (int i = 0; i < order; i++) {
    result += std::to_string(ag_pointer->getStructuralAtomNumber(i)) + " ";
    result += char4ToString(ag_pointer->getAtomName(i)) + " ";
    result += char4ToString(ag_pointer->getResidueName(ag_pointer->getResidueIndex(i))) + " ";
    result += std::to_string(ag_pointer->getResidueNumber(i));
    if (i < order - 1) {
      result += ", ";
    }
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
void BoundedRestraint::checkDisplacementLimits(double4 *rval) {

  // The parameters r1, r2, r3, and r4 must be monotonically increasing.
  if (rval->x > rval->y || rval->y > rval->z || rval->z > rval->w) {
    rtErr("Restraints must operate with monotonically increasing displacement values.  Invalid "
          "displacements for a restraint between atoms " + reportAtomList() + ": " +
          realToString(rval->x, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(rval->y, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(rval->z, 4, NumberFormat::STANDARD_REAL) + ", " +
          realToString(rval->w, 4, NumberFormat::STANDARD_REAL) + ".", "BoundedRestraint",
          "checkDisplacementLimits");
  }
  
  // The range of r2 to r3 must be within [ 0, pi ) for an angle and [ -pi, pi ) for a dihedral.
  // The bound r1 and r4 can be set with whatever remaining room exists in the space.
  if (order == 3 && (rval->y < 0 || rval->z >= pi)) {
    rtErr("Three-point angle restraints require that r2 and r3 stay within the range [0, pi).  "
          "Invalid displacements for a restraint between atoms " + reportAtomList() + ": " +
          "r2 = " + realToString(rval->y, 4, NumberFormat::STANDARD_REAL) + ", " +
          "r3 = " + realToString(rval->z, 4, NumberFormat::STANDARD_REAL) + ".",
          "BoundedRestraint", "checkDisplacementLimits");
  }
  if (order == 4 && rval->z - rval->y >= twopi) {
    rtErr("Four-point dihedral angle restraints require that the difference between r2 and r3 "
          "stay within the range [-pi, pi).  Invalid displacements for a restraint between "
          "atoms " + reportAtomList() + ": r2 = " +
          realToString(rval->y, 4, NumberFormat::STANDARD_REAL) + ", r3 = " +
          realToString(rval->z, 4, NumberFormat::STANDARD_REAL) + ".", "BoundedRestraint",
          "checkDisplacementLimits");
  }
  if (order == 3) {
    if (rval->x < 0.0) {
      rval->x = 0.0;
    }
    else if (rval->w >= pi) {
      rval->w = pi;
    }
  }
  else if (order == 4 && rval->w - rval->x > twopi) {
    rval->x = rval->y - (0.5 * (twopi - (rval->z - rval->y)));
    rval->w = rval->z + (0.5 * (twopi - (rval->z - rval->y)));
  }

  // Ensure that the displacements for a dihedral restraint sit in the range [0, 2*pi) or higher.
  if (order == 4) {
    const double midpoint = 0.5 * (rval->y + rval->z);
    double midpt_copy = imageValue(midpoint, twopi, ImagingMethod::PRIMARY_UNIT_CELL);
    rval->x += midpt_copy - midpoint;
    rval->y += midpt_copy - midpoint;
    rval->z += midpt_copy - midpoint;
    rval->w += midpt_copy - midpoint;
  }
}

} // namespace restraints
} // namespace stormm
