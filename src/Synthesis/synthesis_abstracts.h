// -*-c++-*-
#ifndef STORMM_SYNTHESIS_ABSTRACTS_H
#define STORMM_SYNTHESIS_ABSTRACTS_H

#include "copyright.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Topology/atomgraph_enumerators.h"
#include "synthesis_enumerators.h"

namespace stormm {
namespace synthesis {

using topology::UnitCellType;
using topology::ImplicitSolventModel;
  
/// \brief Collect the critical valence parameters and indexing information for work unit-based
///        evaluation of the systems in an AtomGraphSynthesis.
template <typename T> struct SyValenceKit {

  /// \brief The constructor takes a straight list of arguments for each member variable.
  explicit SyValenceKit(int nvwu_in, const T coulomb_in, const T* bond_keq_in,
                        const T* bond_leq_in, const T* angl_keq_in, const T* angl_theta_in,
                        const T* dihe_amp_in, const T* dihe_freq_in, const T* dihe_phi_in,
                        const T* attn14_elec_in, const T* attn14_vdw_in, const T* charges_in,
                        const T* lja_14_coeff_in, const T* ljb_14_coeff_in,
                        const T* ljc_14_coeff_in, const T* lj_14_sigma_in, const int* lj_idx_in,
                        const int* n_lj_types_in, const int* ljabc_offsets_in,
                        const T* ubrd_keq_in, const T* ubrd_leq_in, const T* cimp_keq_in,
                        const T* cimp_phi_in, const int* cmap_dim_in, const T* cmap_patches_in,
                        const int* cmap_patch_bounds_in, const int2* vwu_abstracts_in,
                        const int* vwu_imports_in, const uint2* cbnd_insr_in,
                        const uint2* angl_insr_in, const uint2* cdhe_insr_in,
                        const uint* cdhe_ovrt_insr_in, const uint2* cmap_insr_in,
                        const uint* infr14_insr_in, const uint* cbnd_acc_in,
                        const uint* angl_acc_in, const uint* cdhe_acc_in, const uint* cmap_acc_in,
                        const uint* infr14_acc_in);

  /// \brief The copy and move constructors are taken at their default values for this abstract
  ///        containing const elements.
  /// \{
  SyValenceKit(const SyValenceKit &original) = default;
  SyValenceKit(SyValenceKit &&original) = default;
  /// \}

  // Member variables (all public)
  const int nvwu;                ///< The number of valence work units
  const T coulomb;               ///< Coulomb's constant in kcal-A2/mol-e2

  // Condensed valence parameter tables
  const T* bond_keq;             ///< Consensus table of harmonic bond stiffness constants
  const T* bond_leq;             ///< Consensus table of harmonic bond equilibrium lengths
  const T* angl_keq;             ///< Consensus table of harmonic angle stiffness constants
  const T* angl_theta;           ///< Consensus table of harmonic angle equilibrium settings
  const T* dihe_amp;             ///< Consensus table of cosine-based dihedral amplitudes
  const T* dihe_freq;            ///< Consensus table of cosine-based dihedral frequencies 
  const T* dihe_phi;             ///< Consensus table of cosine-based dihedral phase angles
  const T* attn14_elec;          ///< Consensus table of 1:4 electrostatic scaling factors
  const T* attn14_vdw;           ///< Consensus table of 1:4 van-der Waals scaling factors
  const T* charges;              ///< Charges on all atoms
  const T* lja_14_coeff;         ///< Lennard-Jones 1:4 interaction A coefficients
  const T* ljb_14_coeff;         ///< Lennard-Jones 1:4 interaction B coefficients
  const T* ljc_14_coeff;         ///< Lennard-Jones 1:4 interaction C coefficients
  const T* lj_14_sigma;          ///< Lennard-Jones 1:4 sigma parameters for all atom type pairs
  const int* lj_idx;             ///< Lennard-Jones type indices for all atoms
  const int* n_lj_types;         ///< Lennard-Jones type counts for all systems
  const int* ljabc_offsets;      ///< Offsets for Lennard-Jones A, B, and C coefficient tables for
                                 ///<   all systems
  const T* ubrd_keq;             ///< Consensus table of Urey-Bradley stiffnesses
  const T* ubrd_leq;             ///< Consensus table of Urey-Bradley I:K separation distances
  const T* cimp_keq;             ///< Consensus table of CHARMM improper dihedral stiffness
                                 ///<   constants
  const T* cimp_phi;             ///< Consensus table of CHARMM improper dihedral phase angles
  const int* cmap_dim;           ///< Dimensions for each member of the consensus table of CMAPs
  const T* cmap_patches;         ///< Consensus table of CMAP pre-computed surface patch data
  const int* cmap_patch_bounds;  ///< Bounds array for cmap_patches
  
  // Details of each work unit and instruction sets
  const int2* vwu_abstracts;     ///< Abstracts for all valence work units
  const int* vwu_imports;        ///< Imported atom tables for all valence work units
  const uint2* cbnd_insr;        ///< Composite bond and Urey-Bradley instructions implemented by
                                 ///<   valence work units
  const uint2* angl_insr;        ///< Harmonic angle instructions implemented by valence work units
  const uint2* cdhe_insr;        ///< Composite cosine-based and CHARMM improper dihedral
                                 ///<   instructions implemented by valence work units
  const uint* cdhe_ovrt_insr;    ///< Instructions for composite cosine-based dihedral overtones
  const uint2* cmap_insr;        ///< CMAP instructions implemented by valence work units
  const uint* infr14_insr;       ///< Inferred 1:4 attenuated interaction instructions implemented
                                 ///<   by valence work units
  const uint* cbnd_acc;          ///< Composite bond energy accumulation masks
  const uint* angl_acc;          ///< Harmonic angle energy accumulation masks
  const uint* cdhe_acc;          ///< Composite dihedral energy accumulation masks
  const uint* cmap_acc;          ///< CMAP energy accumulation masks
  const uint* infr14_acc;        ///< Inferred 1:4 interaction energy accumulation masks

};

/// \brief Collect the critical restraint parameters and masking information for work unit-based
///        evaluation of the systems in an AtomGraphSynthesis.
template <typename T, typename T2, typename T4> struct SyRestraintKit {

  /// \brief The constructor takes a straight list of arguments for each member variable.
  explicit SyRestraintKit(const int2* rposn_step_bounds_in, const int2* rbond_step_bounds_in,
                          const int2* rangl_step_bounds_in, const int2* rdihe_step_bounds_in,
                          const T2* rposn_init_k_in, const T2* rposn_finl_k_in,
                          const T4* rposn_init_r_in, const T4* rposn_finl_r_in,
                          const T2* rposn_init_xy_in, const T* rposn_init_z_in,
                          const T2* rposn_finl_xy_in, const T* rposn_finl_z_in,
                          const T2* rbond_init_k_in, const T2* rbond_finl_k_in,
                          const T4* rbond_init_r_in, const T4* rbond_finl_r_in,
                          const T2* rangl_init_k_in, const T2* rangl_finl_k_in,
                          const T4* rangl_init_r_in, const T4* rangl_finl_r_in,
                          const T2* rdihe_init_k_in, const T2* rdihe_finl_k_in,
                          const T4* rdihe_init_r_in, const T4* rdihe_finl_r_in,
                          const uint2* rposn_insr_in, const uint2* rbond_insr_in,
                          const uint2* rangl_insr_in, const uint2* rdihe_insr_in,
                          const uint* rposn_acc_in, const uint* rbond_acc_in,
                          const uint* rangl_acc_in, const uint* rdihe_acc_in);

  /// \brief The copy and move constructors are taken at their default values for this abstract
  ///        containing const elements.
  /// \{
  SyRestraintKit(const SyRestraintKit &original) = default;
  SyRestraintKit(SyRestraintKit &&original) = default;
  /// \}

  // Member variables (all public)
  const int2* rposn_step_bounds;  ///< Steps for initiating and completely applying positional
                                  ///<   restraints
  const int2* rbond_step_bounds;  ///< Steps for initiating and completely applying distance
                                  ///<   restraints
  const int2* rangl_step_bounds;  ///< Steps for initiating and completely applying angle
                                  ///<   restraints
  const int2* rdihe_step_bounds;  ///< Steps for initiating and completely applying dihedral
                                  ///<   restraints
  const T2* rposn_init_k;         ///< Initial positional restraint stiffness parameters
  const T2* rposn_finl_k;         ///< Final positional restraint stiffness parameters
  const T4* rposn_init_r;         ///< Initial positional restraint displacement parameters
  const T4* rposn_finl_r;         ///< Final positional restraint displacement parameters
  const T2* rposn_init_xy;        ///< Initial positional restraint X and Y reference coordinates
  const T* rposn_init_z;          ///< Initial positional restraint Z reference coordinates
  const T2* rposn_finl_xy;        ///< Final positional restraint X and Y reference coordinates
  const T* rposn_finl_z;          ///< Final positional restraint Z reference coordinates
  const T2* rbond_init_k;         ///< Initial distance restraint stiffness parameters
  const T2* rbond_finl_k;         ///< Final distance restraint stiffness parameters
  const T4* rbond_init_r;         ///< Initial distance restraint displacement parameters
  const T4* rbond_finl_r;         ///< Final distance restraint displacement parameters
  const T2* rangl_init_k;         ///< Initial angle restraint stiffness parameters
  const T2* rangl_finl_k;         ///< Final angle restraint stiffness parameters
  const T4* rangl_init_r;         ///< Initial angle restraint displacement parameters
  const T4* rangl_finl_r;         ///< Final angle restraint displacement parameters
  const T2* rdihe_init_k;         ///< Initial dihedral angel restraint stiffness parameters
  const T2* rdihe_finl_k;         ///< Final dihedral angle restraint stiffness parameters
  const T4* rdihe_init_r;         ///< Initial dihedral angle restraint displacement parameters
  const T4* rdihe_finl_r;         ///< Final dihedral angle restraint displacement parameters
  const uint2* rposn_insr;        ///< Instructions for positional restraints
  const uint2* rbond_insr;        ///< Instructions for distance restraints
  const uint2* rangl_insr;        ///< Instructions for three-point angle restraints
  const uint2* rdihe_insr;        ///< Instructions for four-point dihedral restraints
  const uint* rposn_acc;          ///< Masks of whether to accumulate positional restraint effects
  const uint* rbond_acc;          ///< Masks of whether to accumulate distance restraint effects
  const uint* rangl_acc;          ///< Masks of whether to accumulate angle restraint effects
  const uint* rdihe_acc;          ///< Masks of whether to accumulate dihedral restraint effects
};

/// \brief Collect the critical non-bonded parameters and masking information for work unit-based
///        evaluation of the systems in an AtomGraphSynthesis.  Compared to the abstracts for a
///        single system's AtomGraph, this merges the ImplicitSolventKit and NonbondedKit.  For
///        comparison, the SyValenceKit object merges the information for valence interactions,
///        constraints, and virtual sites abstracts from the single topology's set of abstracts. 
template <typename T, typename T2> struct SyNonbondedKit {

  /// \brief The constructor takes a straight list of arguments for each member variable.
  explicit SyNonbondedKit(int nsys_in, UnitCellType unit_cell_in, int nnbwu_in,
                          const int* nbwu_abstracts_in, const uint2* nbwu_insr_in,
                          const int* atom_offsets, const int* atom_counts, T coulomb_in,
                          ImplicitSolventModel igb_in, int neck_table_size_in, T dielectric_in,
                          T kappa_in, T saltcon_in, T gb_offset_in, T gb_neckscale_in,
                          T gb_neckcut_in, const T* charge_in, const int* q_idx_in,
                          const T* q_params_in, const int* lj_idx_in, const int* n_lj_types_in,
                          const int* ljabc_offsets_in, const T* lja_coeff_in,
                          const T* ljb_coeff_in, const T* ljc_coeff_in, const T* lj_sigma_in,
                          const int* neck_gb_idx_in, const T* pb_radii_in, const T* gb_screen_in,
                          const T* gb_alpha_in, const T* gb_beta_in, const T* gb_gamma_in,
                          const T2* neck_limits_in);

  /// \brief The copy and move constructors are taken at their default values for this abstract
  ///        containing const elements.
  /// \{
  SyNonbondedKit(const SyNonbondedKit &original) = default;
  SyNonbondedKit(SyNonbondedKit &&original) = default;
  /// \}

  // Member variables (all public)
  const int nsys;                 ///< The total number of systems in the synthesis
  const UnitCellType unit_cell;   ///< The type of unit cell for all systems
  const int nnbwu;                ///< The number of non-bonded work units in the synthesis
  const int* nbwu_abstracts;      ///< Abstracts for all non-bonded work units
  const uint2* nbwu_insr;         ///< Instructions for all non-bonded work units
  const int* atom_offsets;        ///< Offsets for the atom indices of each system (needed to
                                  ///<   determine whether tiles overrun the end of any one system)
  const int* atom_counts;         ///< Atom counts for each system (needed to determine whether
                                  ///<   tiles overrun the end of any one system)
  const T coulomb;                ///< Coulomb's constant (charges are stored in atomic units and
                                  ///<   converted to factor this in at the time they are cached)
  const ImplicitSolventModel igb; ///< The flavor of Generalized Born to use, i.e. Hawkins /
                                  ///<   Cramer / Truhlar
  const int neck_table_size;      ///< Length of each side of the square "neck" GB limit tables
  const T dielectric;             ///< Solvent dielectric constant
  const T kappa;                  ///< Inverse of the Debye-Huckel length (based on the dielectic)
  const T saltcon;                ///< The salt concentration to use in GB calculations
  const T gb_offset;              ///< Offset for baseline Generalized Born radii (default 0.9,
                                  ///<   larger for neck GB model II)
  const T gb_neckscale;           ///< Neck function scaling parameter for neck GB
  const T gb_neckcut;             ///< Cutoff for the "neck" GB function
  const T* charge;                ///< Partial charges for all atoms in the synthesis (one
                                  ///<   concatenated array, with each system's atoms padded by the
                                  ///<   warp size)
  const int* q_idx;               ///< Indices of atomic charge parameters.  The charge for atom i
                                  ///<   can be obtained by accessing the ith element of the array
                                  ///<   charge, or by accessing the ith element of this array and
                                  ///<   then following that index into the q_params array.
  const T* q_params;              ///< List of unique atomic partial charge parameters, indexed by
                                  ///<   q_idx above.
  const int* lj_idx;              ///< Lennard-Jones indices for all atoms.  Each system may have
                                  ///<   its own specific parameter matricies, with offsets given
                                  ///<   in the ljabc_offsets array (see below).
  const int* n_lj_types;          ///< Lennard-Jones type counts for all systems
  const int* ljabc_offsets;       ///< Offsets for Lennard-Jones A, B, and C coefficient tables for
                                  ///<   all systems
  const T* lja_coeff;             ///< Lennard-Jones interaction A coefficients
  const T* ljb_coeff;             ///< Lennard-Jones interaction B coefficients
  const T* ljc_coeff;             ///< Lennard-Jones interaction C coefficients
  const T* lj_sigma;              ///< Lennard-Jones sigma parameters for all atom type pairs
  const int* neck_gb_idx;         ///< Neck GB indicies for all atoms in each system, applicable
                                  ///<   for Mongan's "neck" GB models
  const T* pb_radii;              ///< Atomic PB radii for all atoms in all systems
  const T* gb_screen;             ///< Generalized Born screening factors for all atoms and systems
  const T* gb_alpha;              ///< Generalized born alpha parameters (one parameter per atom)
  const T* gb_beta;               ///< Generalized born beta parameters (one parameter per atom)
  const T* gb_gamma;              ///< Generalized born gamma parameters (one parameter per atom)
  const T2* neck_limits;          ///< Interlaced neck GB maximum separations and values tables
};

/// \brief Collect the virtual site details and constraint parameters of the topology synthesis
///        into a single abstract.  This is designed to work in the context of the valence work
///        units from the same topology synthesis, whether in a standalone function or as an extra
///        step in a workflow that compute valence interactions and then moves particles.
template <typename T, typename T2, typename T4> struct SyAtomUpdateKit {

  /// \brief The constructor takes a list of parameter arrays.  Like other AtomGraphSynthesis
  ///        abstracts, the actual numbers of parameters in each array, which are irrelevant to
  ///        the implementation as long as the work unit instructions do not overrun the bounds,
  ///        are omitted to save space on the constants imported as arguments to each kernel.  The
  ///        numbers of each instruction are found in the valence work units.
  SyAtomUpdateKit(const T* masses_in, const T* inv_masses_in, int largest_group_in,
                  const T4* vs_params_in, const T4* settle_geom_in, const T2* settle_mass_in,
                  const T2* cnst_grp_params_in, const uint2* vste_insr_in,
                  const uint2* sett_insr_in, const uint2* cnst_insr_in, const uint2* vwu_manip_in);

  /// \brief The copy and move constructors are taken at their default values for this abstract
  ///        containing const elements.
  /// \{
  SyAtomUpdateKit(const SyAtomUpdateKit &original) = default;
  SyAtomUpdateKit(SyAtomUpdateKit &&original) = default;
  /// \}

  // General information relevant to moving particles
  const T* masses;            ///< Masses of all particles
  const T* inv_masses;        ///< Inverse masses of all particles

  // Parameters and instruction sets for virtual sites and constraints
  const int largest_group;    ///< Number of constrained bonds in the largest hub-and-spoke group
                              ///<   found in any topology
  const T4* vs_params;        ///< Unique virtual site frames, with the first, second, and third
                              ///<   dimension parameters in the 
  const T4* settle_geom;      ///< Geometric considerations for SETTLE-constrained groups.  The
                              ///<   SETTLE instructions will indicate which parameter set, and
                              ///<   thus which index of this array to take.
  const T2* settle_mass;      ///< Inverse mass considerations for SETTLE-constrained groups.  The
                              ///<   SETTLE instructions will indicate which parameter set, and
                              ///<   thus which index of this array to take.  
  const T2* cnst_grp_params;  ///< Squared bond length (x member) and the sum of inverse masses for
                              ///<   both atoms (y member) in one constraint (of perhaps several)
                              ///<   of a hub-and-spoke constraint group.
  const uint2* vste_insr;     ///< Virtual site placement instructions
  const uint2* sett_insr;     ///< SETLLE group constraints instructions
  const uint2* cnst_insr;     ///< Hub-and-spoke constraint instructions  
  const uint2* vwu_manip;     ///< Manipulation masks (movement in x member, update in y member)
                              ///<   for all valence work units
};
  
} // namespace synthesis
} // namespace stormm

#include "synthesis_abstracts.tpp"

#endif

