// -*-c++-*-
#ifndef STORMM_NMR_RESTRAINT_H
#define STORMM_NMR_RESTRAINT_H

#include <string>
#include "copyright.h"
#include "Chemistry/chemical_features.h"
#include "DataTypes/stormm_vector_types.h"
#include "Trajectory/coordinateframe.h"
#include "Topology/atomgraph.h"

namespace stormm {
namespace restraints {

using chemistry::ChemicalFeatures;
using topology::AtomGraph;
using trajectory::CoordinateFrameReader;

class BoundedRestraint {
public:

  /// Constructors take either four atom masks (each of whcih must evaluate to exactly one atom)
  /// or four atom numbers (numbers are given for a series 1... atom_count in the topology pointer,
  /// decremented when constructing the mask to index into the actual memory)
  ///
  /// Overloaded:
  ///   - Empty constructor (creates a restraint of order zero but for a particular topology)
  ///   - Accept up to four atom masks and all other inputs
  ///   - Accept up to four atom numbers and all other inputs
  ///   - Accept up to four atom masks and basic inputs
  ///   - Accept up to four atom numbers and basic inputs
  ///   - If only a single atom is specified by mask or by number, reference coordinates should
  ///     also be supplied
  ///
  /// \param mask_i_in         Atom mask identifying the first atom
  /// \param mask_j_in         Atom mask identifying the second atom (empty string for no atom)
  /// \param mask_k_in         Atom mask identifying the third atom (empty string for no atom)
  /// \param mask_l_in         Atom mask identifying the fourth atom (empty string for no atom)
  /// \param atom_i_in         Topological index of the first atom (starts from 0, skip atom mask)
  /// \param atom_i_in         Topological index of the second atom (starts from 0, skip atom mask)
  /// \param atom_i_in         Topological index of the third atom (starts from 0, skip atom mask)
  /// \param atom_i_in         Topological index of the fourth atom (starts from 0, skip atom mask)
  /// \param ag_in             Topology of the system in question
  /// \param chemfe            Chemical perception output for the system in question 
  /// \param cfr               Coordinates of the system in question (for positional restraints)
  /// \param init_step_in      Step at which to begin to apply the restraints
  /// \param final_step_in     Step at which the retraint is to take on its mature values
  /// \param init_k2_in        Start value of the stiffness of the left-hand parabola
  /// \param init_k3_in        Start value of the stiffness of the right-hand parabola
  /// \param init_r1_in        Start value of the leftmost point of the left-hand half parabola
  /// \param init_r2_in        Start value of the rightmost point of the left-hand half parabola
  /// \param init_r3_in        Start value of the leftmost point of the right-hand half parabola
  /// \param init_r4_in        Start value of the rightmost point of the right-hand half parabola
  /// \param final_k2_in       Mature value of the stiffness of the left-hand parabola
  /// \param final_k3_in       Mature value of the stiffness of the right-hand parabola
  /// \param final_r1_in       Mature value of the leftmost point of the left-hand half parabola
  /// \param final_r2_in       Mature value of the rightmost point of the left-hand half parabola
  /// \param final_r3_in       Mature value of the leftmost point of the right-hand half parabola
  /// \param final_r4_in       Mature value of the rightmost point of the right-hand half parabola
  /// \param init_ref_crd_in   The initial target location for an atomic positional restraint
  /// \param final_ref_crd_in  The final target location for an atomic positional restraint
  /// \{
  BoundedRestraint(const AtomGraph *ag_in);
  
  BoundedRestraint(const std::string &mask_i_in, const std::string &mask_j_in, 
                   const std::string &mask_k_in, const std::string &mask_l_in,
                   const AtomGraph *ag_in, const ChemicalFeatures &chemfe,
                   const CoordinateFrameReader &cfr, int init_step_in, int final_step_in,
                   double init_k2_in, double init_k3_in, double init_r1_in, double init_r2_in,
                   double init_r3_in, double init_r4_in, double final_k2_in,
                   double final_k3_in, double final_r1_in, double final_r2_in,
                   double final_r3_in, double final_r4_in,
                   const double3 init_ref_crd_in = { 0.0, 0.0, 0.0 },
                   const double3 final_ref_crd_in = { 0.0, 0.0, 0.0 });

  BoundedRestraint(const std::string &mask_i_in, const std::string &mask_j_in, 
                   const std::string &mask_k_in, const AtomGraph *ag_in,
                   const ChemicalFeatures &chemfe, const CoordinateFrameReader &cfr,
                   int init_step_in, int final_step_in, double init_k2_in, double init_k3_in,
                   double init_r1_in, double init_r2_in, double init_r3_in, double init_r4_in,
                   double final_k2_in, double final_k3_in, double final_r1_in,
                   double final_r2_in, double final_r3_in, double final_r4_in);

  BoundedRestraint(const std::string &mask_i_in, const std::string &mask_j_in, 
                   const AtomGraph *ag_in, const ChemicalFeatures &chemfe,
                   const CoordinateFrameReader &cfr, int init_step_in, int final_step_in,
                   double init_k2_in, double init_k3_in, double init_r1_in, double init_r2_in,
                   double init_r3_in, double init_r4_in, double final_k2_in,
                   double final_k3_in, double final_r1_in, double final_r2_in,
                   double final_r3_in, double final_r4_in);

  BoundedRestraint(const std::string &mask_i_in, const AtomGraph *ag_in,
                   const ChemicalFeatures &chemfe, const CoordinateFrameReader &cfr,
                   int init_step_in, int final_step_in, double init_k2_in, double init_k3_in,
                   double init_r1_in, double init_r2_in, double init_r3_in, double init_r4_in,
                   double final_k2_in, double final_k3_in, double final_r1_in, double final_r2_in,
                   double final_r3_in, double final_r4_in, const double3 init_ref_crd_in,
                   const double3 final_ref_crd_in);
  
  BoundedRestraint(const std::string &mask_i_in, const std::string &mask_j_in, 
                   const std::string &mask_k_in, const std::string &mask_l_in,
                   const AtomGraph *ag_in, const ChemicalFeatures &chemfe,
                   const CoordinateFrameReader &cfr, double k2_in, double k3_in,
                   double r1_in, double r2_in, double r3_in, double r4_in);

  BoundedRestraint(const std::string &mask_i_in, const std::string &mask_j_in,
                   const std::string &mask_k_in, const AtomGraph *ag_in,
                   const ChemicalFeatures &chemfe, const CoordinateFrameReader &cfr, double k2_in,
                   double k3_in, double r1_in, double r2_in, double r3_in, double r4_in);

  BoundedRestraint(const std::string &mask_i_in, const std::string &mask_j_in,
                   const AtomGraph *ag_in, const ChemicalFeatures &chemfe,
                   const CoordinateFrameReader &cfr, double k2_in, double k3_in,
                   double r1_in, double r2_in, double r3_in, double r4_in);

  BoundedRestraint(const std::string &mask_i_in, const AtomGraph *ag_in,
                   const ChemicalFeatures &chemfe, const CoordinateFrameReader &cfr, double k2_in,
                   double k3_in, double r1_in, double r2_in, double r3_in, double r4_in,
                   const double3 ref_crd_in);

  BoundedRestraint(const std::string &mask_i_in, const AtomGraph *ag_in,
                   const ChemicalFeatures &chemfe, const CoordinateFrameReader &cfr, double k2_in,
                   double k3_in, double r1_in, double r2_in, double r3_in, double r4_in,
                   const std::string & mask_ref_in);

  BoundedRestraint(int atom_i_in, int atom_j_in, int atom_k_in, int atom_l_in,
                   const AtomGraph *ag_in, int init_step_in, int final_step_in, double init_k2_in,
                   double init_k3_in, double init_r1_in, double init_r2_in, double init_r3_in,
                   double init_r4_in, double final_k2_in, double final_k3_in,
                   double final_r1_in, double final_r2_in, double final_r3_in,
                   double final_r4_in, const double3 init_ref_crd_in = { 0.0, 0.0, 0.0},
                   const double3 final_ref_crd_in = { 0.0, 0.0, 0.0 });

  BoundedRestraint(int atom_i_in, int atom_j_in, int atom_k_in, const AtomGraph *ag_in,
                   int init_step_in, int final_step_in, double init_k2_in, double init_k3_in,
                   double init_r1_in, double init_r2_in, double init_r3_in, double init_r4_in,
                   double final_k2_in, double final_k3_in, double final_r1_in, double final_r2_in,
                   double final_r3_in, double final_r4_in);

  BoundedRestraint(int atom_i_in, int atom_j_in, const AtomGraph *ag_in, int init_step_in,
                   int final_step_in, double init_k2_in, double init_k3_in, double init_r1_in,
                   double init_r2_in, double init_r3_in, double init_r4_in, double final_k2_in,
                   double final_k3_in, double final_r1_in, double final_r2_in, double final_r3_in,
                   double final_r4_in);

  BoundedRestraint(int atom_i_in, const AtomGraph *ag_in, int init_step_in, int final_step_in,
                   double init_k2_in, double init_k3_in, double init_r1_in, double init_r2_in,
                   double init_r3_in, double init_r4_in, double final_k2_in, double final_k3_in,
                   double final_r1_in, double final_r2_in, double final_r3_in, double final_r4_in,
                   const double3 init_ref_crd_in, const double3 final_ref_crd_in);

  BoundedRestraint(int atom_i_in, int atom_j_in, int atom_k_in, int atom_l_in,
                   const AtomGraph *ag_in, double k2_in, double k3_in, double r1_in,
                   double r2_in, double r3_in, double r4_in);

  BoundedRestraint(int atom_i_in, int atom_j_in, int atom_k_in, const AtomGraph *ag_in,
                   double k2_in, double k3_in, double r1_in, double r2_in, double r3_in,
                   double r4_in);

  BoundedRestraint(int atom_i_in, int atom_j_in, const AtomGraph *ag_in, double k2_in,
                   double k3_in, double r1_in, double r2_in, double r3_in, double r4_in);
  
  BoundedRestraint(int atom_i_in, const AtomGraph *ag_in, double k2_in, double k3_in, double r1_in,
                   double r2_in, double r3_in, double r4_in, const double3 ref_crd_in);

  BoundedRestraint(int atom_index, const AtomGraph *ag_in, const CoordinateFrameReader &cfr,
                   double k2_in, double k3_in, double r1_in, double r2_in, double r3_in,
                   double r4_in, int refr_index = -1);
  /// \}

  /// \brief Take the default copy, copy assignment, move, and move assignment constructors.
  /// \{
  BoundedRestraint(const BoundedRestraint &original) = default;
  BoundedRestraint(BoundedRestraint &&original) = default;
  BoundedRestraint& operator=(const BoundedRestraint &original) = default;
  BoundedRestraint& operator=(BoundedRestraint &&original) = default;
  /// \}
  
  /// Obtain the topology index of an atom in this restraint
  ///
  /// \param restrained_atom_number  The 1st, 2nd, 3rd, or 4th atom (specify 1, 2, 3, or 4)
  int getAtomIndex(int restrained_atom_number) const;

  /// \brief Get the order of this restraint
  int getOrder() const;

  /// \brief Get the initial step at which to begin applying this restraint
  int getInitialStep() const;

  /// \brief Get the simulation step at which to finish applying this restraint
  int getFinalStep() const;

  /// \brief Get the stiffnesses of a restraint at a given step in the simulation
  ///
  /// \param step_number  The step at which to compute the target site
  double2 getStiffness(int step_number = 0) const;

  /// \brief Get the initial stiffnesses to use when applying this restraint
  double2 getInitialStiffness() const;

  /// \brief Get the final stiffnesses of the restraint in its complete form
  double2 getFinalStiffness() const;

  /// \brief Get the displacement parameters of a restraint at a given step in the simulation
  ///
  /// \param step_number  The step at which to compute the target site
  double4 getDisplacements(int step_number = 0) const;

  /// \brief Get the initial displacement parameters to use in applying this restraint
  double4 getInitialDisplacements() const;

  /// \brief Get the final displacement parameters of the restraint in its complete form
  double4 getFinalDisplacements() const;

  /// \brief Get the target site of a positional restraint at a given step in the simulation
  ///
  /// \param step_number  The step at which to compute the target site
  double3 getTargetSite(int step_number = 0) const;
  
  /// \brief Get the initial target of a positional restraint
  double3 getInitialTargetSite() const;

  /// \brief Get the final target of a positional restraint
  double3 getFinalTargetSite() const;

  /// \brief Get the topology pointer
  const AtomGraph* getTopologyPointer() const;

  /// \brief Modify the initial step at which the restraint is applied.
  ///
  /// \param new_init_step  The new value for the initial application step
  void setInitialStep(int new_init_step);

  /// \brief Modify the final step at which the restraint application reaches its final value.
  ///
  /// \param new_final_step  The new value for the final application step
  void setFinalStep(int new_final_step);

  /// \brief Modify the stiffness parameters, set with a single value.  This member function is
  ///        only applicable if there is no time dependence in the restraint and will throw an
  ///        error otherwise.
  ///
  /// \param new_keq  New value of both k2 and k3.  Both values of the initial stiffness tuple
  ///                 will be set to these values.
  void setStiffness(double new_keq);

  /// \brief Modify the stiffness parameters.  This member function is only applicable if there is
  ///        no time dependence in the restraint and will throw an error otherwise.
  ///
  /// \param new_k2  New value of k2
  /// \param new_k3  New value of k3
  void setStiffnesses(double new_k2, double new_k3);

  /// \brief Modify the initial stiffness parameters, setting both to a single value.  This member
  ///        function is only applicable if there is time dependence in the restraint.
  ///
  /// \param new_init_keq  New initial value of both k2 and k3
  void setInitialStiffness(double new_init_keq);

  /// \brief Modify the initial stiffness parameters.  This member function is only applicable if
  ///        there is time dependence in the restraint.
  ///
  /// \param new_init_k2  New initial value of k2
  /// \param new_init_k3  New initial value of k3
  void setInitialStiffnesses(double new_init_k2, double new_init_k3);

  /// \brief Modify the final stiffness parameters, setting both to a single value.  This member
  ///        function is only applicable if there is time dependence in the restraint.
  ///
  /// \param new_init_keq  New final value of both k2 and k3
  void setFinalStiffness(double new_init_keq);

  /// \brief Modify the final stiffness parameters.  This member function is only applicable if
  ///        there is time dependence in the restraint.
  ///
  /// \param new_init_k2  New final value of k2
  /// \param new_init_k3  New final value of k3
  void setFinalStiffnesses(double new_init_k2, double new_init_k3);

  /// \brief Modify the displacements over which each segment of the flat-bottom restraint applies.
  ///        This member function requires that the restraint be time-independent.
  ///
  /// \param new_r1  New value of r1, the leftmost point of a harmonic potential scaled by k2
  /// \param new_r2  New value of r2, the rightmost point of a harmonic potential scaled by k2
  /// \param new_r3  New value of r3, the leftmost point of a harmonic potential scaled by k3
  /// \param new_r4  New value of r4, the rightmost point of a harmonic potential scaled by k3
  void setDisplacements(double new_r1, double new_r2, double new_r3, double new_r4);

  /// \brief Modify the initial displacements over which each segment of the flat-bottom restraint
  ///        applies.  This member function requires that the restraint be time-dependent.
  ///
  /// \param new_r1  New initial r1, the leftmost point of a harmonic potential scaled by k2
  /// \param new_r2  New initial r2, the rightmost point of a harmonic potential scaled by k2
  /// \param new_r3  New initial r3, the leftmost point of a harmonic potential scaled by k3
  /// \param new_r4  New initial r4, the rightmost point of a harmonic potential scaled by k3
  void setInitialDisplacements(double new_r1, double new_r2, double new_r3, double new_r4);

  /// \brief Modify the final displacements over which each segment of the flat-bottom restraint
  ///        applies.  This member function requires that the restraint be time-dependent.
  ///
  /// \param new_r1  New final r1, the leftmost point of a harmonic potential scaled by k2
  /// \param new_r2  New final r2, the rightmost point of a harmonic potential scaled by k2
  /// \param new_r3  New final r3, the leftmost point of a harmonic potential scaled by k3
  /// \param new_r4  New final r4, the rightmost point of a harmonic potential scaled by k3
  void setFinalDisplacements(double new_r1, double new_r2, double new_r3, double new_r4);

  /// \brief Set the target site.  This function can only apply to a time-independent positional
  ///        restraint with only one active atom.
  ///
  /// Overloaded:
  ///   - Take separate x, y, and z coordinate arguments
  ///   - Take a 3-tuple of real numbers for x, y, and z coordinates
  ///
  /// \param new_ref_x    New reference Cartesian x coordinate
  /// \param new_ref_y    New reference Cartesian y coordinate
  /// \param new_ref_z    New reference Cartesian z coordinate
  /// \param new_ref_crd  New reference Cartesian coordinate tuple
  /// \{
  void setTargetSite(double new_ref_x, double new_ref_y, double new_ref_z);
  void setTargetSite(double3 new_ref_crd);
  /// \}

  /// \brief Set the initial target site.  This function can only apply to a time-dependent
  ///        positional restraint with only one active atom.
  ///
  /// Overloaded:
  ///   - Take separate x, y, and z coordinate arguments
  ///   - Take a 3-tuple of real numbers for x, y, and z coordinates
  ///
  /// \param new_ref_x    New initial reference Cartesian x coordinate
  /// \param new_ref_y    New initial reference Cartesian y coordinate
  /// \param new_ref_z    New initial reference Cartesian z coordinate
  /// \param new_ref_crd  New initial reference Cartesian coordinate tuple
  /// \{
  void setInitialTargetSite(double new_ref_x, double new_ref_y, double new_ref_z);
  void setInitialTargetSite(double3 new_ref_crd);
  /// \}

  /// \brief Set the final target site.  This function can only apply to a time-dependent
  ///        positional restraint with only one active atom.
  ///
  /// Overloaded:
  ///   - Take separate x, y, and z coordinate arguments
  ///   - Take a 3-tuple of real numbers for x, y, and z coordinates
  ///
  /// \param new_ref_x    New final reference Cartesian x coordinate
  /// \param new_ref_y    New final reference Cartesian y coordinate
  /// \param new_ref_z    New final reference Cartesian z coordinate
  /// \param new_ref_crd  New final reference Cartesian coordinate tuple
  /// \{
  void setFinalTargetSite(double new_ref_x, double new_ref_y, double new_ref_z);
  void setFinalTargetSite(double3 new_ref_crd);
  /// \}

private:
  int atom_i;             ///< Index of atom I in the corresponding topology
  int atom_j;             ///< Index of atom J in the corresponding topology
  int atom_k;             ///< Index of atom K in the corresponding topology
  int atom_l;             ///< Index of atom L in the corresponding topology
  int order;              ///< Order of the NMR restraint: 2 = distance restraint between atoms I
                          ///<   and J (implied if K and L are unspecified), 3 = angle restraint
                          ///<   between atoms I, J, and K (implied if L is unspecified),
                          ///<   4 = dihedral restraint between atoms I, J, K, and L
  int initial_step;       ///< Initial step of the simulation at which to begin applying the
                          ///<   restraint with initial keq and r parameters
  int final_step;         ///< Step of the simulation at which to finish applying the full
                          ///<   restraint, with final kew and r parameters
  double2 initial_keq;    ///< Initial stiffness constants for parabolic restraints between points
                          ///<   r1 and r2 (x member of the tuple) and points r3 and r4 (y member
                          ///<   of the tuple)
  double4 initial_r;      ///< Initial displacement parameters r1 (x), r2 (y), r3 (z), and r4 (w)
  double2 final_keq;      ///< Final stiffness constants
  double4 final_r;        ///< Final displacement parameters
  double3 initial_center; ///< Initial center of the restraint potential (positional restraints
                          ///<   only)
  double3 final_center;   ///< Final center of the restraint potential (positional restraints only)
  
  /// Pointer to the topology for which this restraint applies
  const AtomGraph *ag_pointer;  

  /// \brief Report the list of atoms in this restraint as a string, giving structural atom and
  ///        residue numbers as well as names from the topology in a format reminiscent of the PDB:
  ///        [Atom #] [Atom Name] [Residue Name] [Residue #]
  std::string reportAtomList();

  /// \brief Check the displacements of this restraint, in particular for angle restraints that
  ///        only apply within a range that the angle could possibly take and for dihedral
  ///        restraints for which the angle is periodic.
  void checkDisplacementLimits(double4 *rval);
};

} // namespace restraints
} // namespace stormm

#endif
