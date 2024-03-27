// -*-c++-*-
#ifndef STORMM_MESH_FORCEFIELD_H
#define STORMM_MESH_FORCEFIELD_H

#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
#    include <cuda.h>
#    include <cuda_runtime.h>
#  endif
#endif
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/symbol_values.h"
#include "Math/math_enumerators.h"
#include "Potential/energy_enumerators.h"
#include "Potential/soft_core_potentials.h"
#include "Reporting/error_format.h"
#include "Topology/atomgraph.h"

namespace stormm {
namespace structure {

using card::Hybrid;
using card::HybridKind;
using card::HybridTargetLevel;
using energy::cubicSoftCore;
using energy::quinticSoftCore;
using energy::VdwCombiningRule;
using stmath::Interpolant;
using symbols::amber_ancient_bioq;
using topology::AtomGraph;

/// \brief An object to store the probe's specific non-bonded parameters and softcore function
///        coefficients used in mesh construction or mechanics.
template <typename T> struct MeshFFKit {

  /// \brief The constructor takes arguments for all member variables.  One additional constructor
  ///        is provided to create a blank object, as this abstract may need to be mocked as a
  ///        placeholder to submit to some functions.
  /// \{
  MeshFFKit();
  MeshFFKit(VdwCombiningRule lj_rule_in, double coulomb_in, double clash_ratio_in,
            double clash_distance_in, const T* probe_lja_in, const T* probe_ljb_in,
            const T* probe_ljsig_in, const T* softcore_lja_in, const T* softcore_ljb_in,
            const T* softcore_ljc_in, const T* softcore_ljd_in, const T* softcore_lje_in,
            const T* softcore_ljf_in, const T* softcore_qq_in);
  /// \}

  /// \brief The default copy and move constructors will be valid for this object.  Const members
  ///        negate the use of default copy and move assignment operators.
  ///
  /// \param original  The object to copy or move
  /// \{
  MeshFFKit(const MeshFFKit<T> &original) = default;
  MeshFFKit(MeshFFKit<T> &&original) = default;
  /// \}

  // None of the member variables are mutable.
  const VdwCombiningRule ljrule;  ///< Combining rule used to determine Lennard-Jones parameters
  const double coulomb;           ///< The Coulomb constant for electrostatic interactions.
  const float coulomb_f;          ///< Single-precision variant of the Coulomb constant
  const double clash_ratio;       ///< Ratio of the inter-particle distance to the Lennard-Jones
                                  ///<   sigma parameter at which the softcore polynomial potential
                                  ///<   pairwise takes over
  const float clash_ratio_f;      ///< Single-precision variant of clash_ratio
  const double clash_distance;    ///< The inter-paricle distance beneath which an electrostatic
                                  ///<   interaction is declared to be in conflict and the softcore
                                  ///<   polynomial-based potential takes over
  const float clash_distance_f;   ///< Single-precision variant of the distance at which two
                                  ///<   particles' electrostatic interaction is deemed to conflict
  const T* probe_lja;             ///< Pairwise Lennard-Jones A coefficients between the mesh probe
                                  ///<   and each atom type in the underlying topology
  const T* probe_ljb;             ///< Pairwise Lennard-Jones B coefficients between the mesh probe
                                  ///<   and each atom type in the underlying topology
  const T* probe_ljsig;           ///< Pairwise Lennard-Jones sigma radii between the mesh probe
                                  ///<   and each atom type in the underlying topology
  const T* softcore_lja;          ///< Array of polynomial coefficients for the highest-order term
                                  ///<   in the softcore function for the mesh probe interacting
                                  ///<   with each atom type in the underlying topology
  const T* softcore_ljb;          ///< Array of polynomial coefficients for the next highest-order
                                  ///<   terms
  const T* softcore_ljc;          ///< Arrays of additional Lennard-Jones softcore polynomial
                                  ///<   coefficients
  const T* softcore_ljd;          ///< Arrays of additional Lennard-Jones softcore polynomial
                                  ///<   coefficients
  const T* softcore_lje;          ///< Arrays of additional Lennard-Jones softcore polynomial
                                  ///<   coefficients
  const T* softcore_ljf;          ///< Arrays of additional Lennard-Jones softcore polynomial
                                  ///<   coefficients
  const T* softcore_qq;           ///< A single array of coefficients for the electrostatic
                                  ///<   softcore polynomial
};

/// \brief A class to hold the rules by which the system underlying a mesh object interacts with
///        its surroundings.  The mesh object itself may have one probe or an array of them (for an
///        array of potential surfaces), but the constants, rules, and modifying potentials in this
///        object will determine how those probes interact with the underlying system.
template <typename T> class MeshForceField {
public:

  /// \brief The object is a convenient place to store constants and arrays of coefficients,
  ///        nothing more.  The constructor accepts these constants as well as sizing information
  ///        for its arrays and allocates the appropriate amount of memory.  Filling the arrays is
  ///        a matter for the particular mesh object "creating" (perhaps more precise to say
  ///        "populating") the force field.
  /// \}
  MeshForceField(VdwCombiningRule mixing_protocol_in = VdwCombiningRule::LORENTZ_BERTHELOT,
                 double coulomb_constant_in = amber_ancient_bioq,
                 double clash_ratio_in = 0.0, double clash_distance_in = 0.0,
                 int lj_type_count_in = 0);

  MeshForceField(VdwCombiningRule mixing_protocol_in, double clash_ratio_in,
                 double clash_distance_in, const AtomGraph *ag);

  MeshForceField(VdwCombiningRule mixing_protocol_in, double clash_ratio_in,
                 double clash_distance_in, const AtomGraph &ag);

  MeshForceField(VdwCombiningRule mixing_protocol_in, double coulomb_constant_in,
                 double clash_ratio_in, double clash_distance_in,
                 const std::vector<double> &probe_sigma,
                 const std::vector<double> &probe_epsilon);
  /// \}

  /// \brief The copy and move constructors as well as assignment operators must be written out to
  ///        repair various pointers.
  ///
  /// \param original  The original object to copy or move
  /// \param other     An existing object fulfilling the right hand side of an assignment statement
  /// \{
  MeshForceField(const MeshForceField<T> &original);
  MeshForceField(MeshForceField<T> &&original);
  MeshForceField& operator=(const MeshForceField<T> &other);
  MeshForceField& operator=(MeshForceField<T> &&other);
  /// \}

  /// \brief Get the Lennard-Jones mixing rule.
  VdwCombiningRule getCombiningRule() const;

  /// \brief Get the ratio of the van-der Waals (Lennard-Jones) sigma parameter for interacting
  ///        pairs of particles at which softcore van-der Waals interactions take over.
  double getClashRatio() const;

  /// \brief Get the absolute distance at which softcore electrostatic interactions take over.
  double getClashDistance() const;

  /// \brief Get the abstract in whatever templated type (float or double) the object is cast in.
  ///
  /// \param tier  Obtain the data at the level of the CPU host or GPU device
  const MeshFFKit<T> data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the refernce data in double precision, regardless of how the object is cast.
  ///
  /// \param tier  Obtain the data at the level of the CPU host or GPU device
  const MeshFFKit<double> referenceData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a template-free form of the object's abstract in its native type.
  ///
  /// \param tier  Obtain the data at the level of the CPU host or GPU device
  const MeshFFKit<void> templateFreeData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  /// \brief Get the abstract, visible on the device but mounted on the host, in the object's
  ///        native type.
  const MeshFFKit<T> deviceViewToHostData() const;

  /// \brief Get the reference data abstract, visible on the device but mounted on the host, in
  ///        the object's native type.
  const MeshFFKit<double> deviceViewToReferenceHostData() const;

  /// \brief Get the abstract, visible on the device but mounted on the host, in the object's
  ///        native type but with templating removed.
  const MeshFFKit<void> deviceViewToTemplateFreeHostData() const;
#  endif

  /// \brief Upload all data to the device layer.
  void upload();

  /// \brief Download all data from the device layer.
  void download();
#endif

  /// \brief Set the combining rule that will be used to make the probe interact with the receptor
  ///        on any mesh (with the exception of an electrostatic field).  Checks on the validity of
  ///        this setting must be made by the mesh that assembles the MeshForceField object.
  ///
  /// \param mixing_protocol_in  The method to use
  void setCombiningRule(VdwCombiningRule mixing_protocol_in);

  /// \brief Set the absolute distance at which electrostatic charges will begin to be damped
  ///        according to the softcore scheme (a harmonic potential with a relative maximum at
  ///        r = -1.0 takes over).
  ///
  /// \param clash_distance_in  The handoff range to set
  void setClashDistance(double clash_distance_in);

  /// \brief Set the ratio of the van-der Waals (Lennard-Jones) sigma parameters at which softcore
  ///        interactions take over (the interaction switches to a quartic potential with a
  ///        relative maximum at r = -1.0).
  ///
  /// \param clash_ratio_in  The handoff sigma ratio to set
  void setClashRatio(double clash_ratio_in);

  /// \brief Set a specific coefficient in the softcore electrostatic polynomial.
  ///
  /// \param coef  The coefficient value to apply
  /// \param pos   Order of the coefficient, e.g. in Ax^3 + Bx^2 + Cx + D, D has order zero and
  ///              C has order 1
  void setElecSoftcoreParameter(double coef, int pos);

  /// \brief Compute the polynomial for an electrostatic softcore potential suitable for a certain
  ///        interpolation strategy.
  ///
  /// \param stencil_kind  Prepare for computing energies that will feed into this type of
  ///                      interpolation stencil.  SMOOTHNESS interpolants require a fifth-order
  ///                      softcore function providing for three layers of continuous derivatives
  ///                      at the handoff, FUNCTION_VALUE interpolants, while more expensive to
  ///                      compute, require only third-order softcore polynomials with a single
  ///                      continuous derivative at the handoff.
  void setElecSoftcoreParameters(Interpolant stencil_kind);

  /// \brief Set the Lennard-Jones pairwise coefficients to describe the interaction of one atom
  ///        type in the mesh's underlying topology with its probe.
  ///
  /// \param index   Index of the atom type in the underlying topology
  /// \param lja_in  The Lennard-Jones A coefficient, as in U(LJ) = A/r^12 - B/r^6
  /// \param ljb_in  The Lennard-Jones B coefficient
  void setLJCoefficients(int index, double lja_in, double ljb_in);

  /// \brief Set the Lennard-Jones pairwise coefficients based on the mesh's underlying topology.
  ///
  /// Overloaded:
  ///   - Provide the topology as a const pointer or const reference
  ///   - Provide the a single value for the probe radius (sigma) and well depth (epsilon) to be
  ///     combined with the topology's parameters according to the MeshForceField object's rule
  ///
  /// \param ag                The mesh's underlying topology
  /// \param probe_radius      The spherical probe's radius (its sigma parameter if computing a
  ///                          Lennard-Jones non-bonded potential
  /// \param probe_well_depth  The probe epsilon parameter
  /// \{
  void setLJCoefficients(const AtomGraph *ag, const std::vector<double> &probe_radius,
                         const std::vector<double> &probe_well_depth);

  void setLJCoefficients(const AtomGraph &ag, const std::vector<double> &probe_radius,
                         const std::vector<double> &probe_well_depth);

  void setLJCoefficients(const AtomGraph *ag, double probe_radius, double probe_well_depth);

  void setLJCoefficients(const AtomGraph &ag, double probe_radius, double probe_well_depth);
  /// \}

  /// \brief Set a specific softcore polynomial coefficient.
  ///
  /// \param index  Lennard-Jones type index of a mesh atom type for which to apply the coefficient
  /// \param coef   The coefficient value to apply
  /// \param pos    Order of the coefficient, e.g. in Ax^3 + Bx^2 + Cx + D, D has order zero and
  ///               C has order 1
  void setLJSoftcoreParameter(int index, double coef, int pos);

  /// \brief Compute softcore polynomials for the loaded Lennard-Jones parameters, ratio of the
  ///        sigma parameter at which the handoff is to occur, and a given interpolation strategy.
  ///
  /// \param stencil_kind  Prepare for computing energies that will feed into this type of
  ///                      interpolation stencil.
  void setLJSoftcoreParameters(Interpolant stencil_kind);
  
  /// \brief Reallocate the object to hold a different number of Lennard-Jones types.  The
  ///        electrostatic softcore potential coefficients will be held in a temporary array during
  ///        this process so that they may be preserved.
  ///
  /// \param lj_type_count_in  The new number of Lennard-Jones types to prepare for
  void reallocate(const int lj_type_count_in);
  
private:

  /// The number of Lennard-Jones types for which this object is prepared to hold softcore
  /// coefficients.
  int lj_type_count;
  
  /// The combining rule used (or to use) in computing clashes and Lennard-Jones interactions.
  VdwCombiningRule mixing_protocol;

  /// The coulomb constant to use in electrostatics-based meshes, likely copied from the underlying
  /// topology.
  double coulomb_constant;

  /// Ratio of the sigma parameter at which to begin damping Lennard-Jones (van-der Waals)
  /// interactions between particles.  This follows the methods of softcore potentials used in
  /// energy minimizatons to make the potential mapped to various grid points tractable.
  double clash_ratio;

  /// Absolute range at which to begin damping electrostatic interactions.  This follows the same
  /// approach as dampening interactions in non-bonded softcore potentials for energy
  /// minimizations.
  double clash_distance;

  // The following arrays store atom type-specific Lennard-Jones parameters for meshes detailing
  // steric clashes or Lennard-Jones potentials.  If the combining rule is NBFIX, the pairwise
  // parameters for interactions of the probe with each type of atom in the underlying topology
  // must be provided.  Otherwise, these parameters will be inferred based on the topology and a
  // single set of parameters for the solvent probe using one of the other combining rules.
  Hybrid<T> probe_lja;           ///< Lennard-Jones A coefficients (4 epsilon * sigma^12) for the
                                 ///<   solvent probe interacting with each atom type in the
                                 ///<   underlying topology
  Hybrid<T> probe_ljb;           ///< Lennard-Jones B coefficients (4 epsilon * sigma^6) for the
                                 ///<   solvent probe interacting with each atom type in the
                                 ///<   underlying topology
  Hybrid<T> probe_lj_sigma;      ///< Pairwise Lennard-Jones sigma parameters for the solvent probe
                                 ///<   interacting with each atom type in the underlying topology
  Hybrid<T> probe_softcore_lja;  ///< An array of coefficients for softcore potentials applicable
                                 ///<   to each Lennard-Jones pair interaction.  This is the first
                                 ///<   array providing the A coefficients for a polynomial
                                 ///<   Ar^(n) + Br^(n-1) + ... + M.  The softcore polynomial takes
                                 ///<   over at inter-particle distances of less than clash_ratio
                                 ///<   times the pairwise sigma parameter found in probe_ljsig.
  Hybrid<T> probe_softcore_ljb;  ///< A second array of softcore Lennard-Jones polynomial
                                 ///<   coefficients.
  Hybrid<T> probe_softcore_ljc;  ///< Another array of softcore Lennard-Jones polynomial
                                 ///<   coefficients.
  Hybrid<T> probe_softcore_ljd;  ///< Another array of softcore Lennard-Jones polynomial
                                 ///<   coefficients.
  Hybrid<T> probe_softcore_lje;  ///< Another array of softcore Lennard-Jones polynomial
                                 ///<   coefficients.
  Hybrid<T> probe_softcore_ljf;  ///< A final array of softcore Lennard-Jones polynomial
                                 ///<   coefficients.  The softcore polynomial may take up to
                                 ///<   quintic order.
  Hybrid<T> elec_softcore;       ///< Coefficients for the electrostatic softcore polynomial,
                                 ///<   taking over from the Coulomb interaction at
                                 ///<   inter-particle distances of less than clash_distance

  /// ARRAY-kind Hybrid targeted by the POINTER-kind softcore Hybrid objects above
  Hybrid<T> softcore_data;
  
  /// In addition to the non-bonded parameters for mesh-based calculations, it is important to have
  /// a double-precision reference for each of the softcore coefficient arrays above, as they may
  /// be used in computing the mesh itself and the precision model for that could be single- or
  /// double-precision.  The following arrays will be loaded into an explicitly requested
  /// <double>-typed MeshFFKit.
  /// \{
  Hybrid<double> ref_probe_lja;
  Hybrid<double> ref_probe_ljb;
  Hybrid<double> ref_probe_lj_sigma;
  Hybrid<double> ref_probe_softcore_lja;
  Hybrid<double> ref_probe_softcore_ljb;
  Hybrid<double> ref_probe_softcore_ljc;
  Hybrid<double> ref_probe_softcore_ljd;
  Hybrid<double> ref_probe_softcore_lje;
  Hybrid<double> ref_probe_softcore_ljf;
  Hybrid<double> ref_elec_softcore;
  Hybrid<double> ref_softcore_data;
  /// \}

  /// Allocate memory for the object (this will also set the various pointers).  Memory for all
  /// coefficient arrays, including electrostatics and van-der Waals parameters with the requested
  /// number of types, will be allocated.
  ///
  /// \param lj_type_count_in  Allocate storage for this number of Lennard-Jones coefficients,
  ///                          sigma parameters, and softcore polynomials
  void allocate(int lj_type_count_in);
  
  /// Repair the object's pointers after a deep copy.
  void rebasePointers();

  /// Validate the requested Lennard-Jones type index.
  ///
  /// \param index  The Lennard-Jones type of interest
  void validateTypeIndex(int index);
};

/// \brief Restore the true data type of a void-casted MeshForceField abstract.  This will be
///        compiled by both C++ and HPC compilers that build executable objects dependent on the
///        MeshForceField class, such that executable objects built by either compiler will have
///        the means to remove and restore templated character to abstracts of MeshForceField.
///
/// Overloaded:
///   - Accept the original object by pointer or by reference (the original object will not be
///     altered)
///   - Overloads in other libraries produce outputs to match their own input types
///
/// \param rasa  Abstract of the original mesh object, free of templated pointers (inspired by the
///              phrase "tabula rasa", a blank slate)
/// \{
template <typename T> MeshFFKit<T> restoreType(const MeshFFKit<void> *rasa);
template <typename T> MeshFFKit<T> restoreType(const MeshFFKit<void> &rasa);
/// \}

} // namespace structure
} // namespace stormm

#include "mesh_forcefield.tpp"

#endif
