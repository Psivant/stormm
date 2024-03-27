// -*-c++-*-
#ifndef STORMM_ATOMGRAPH_ABSTRACTS_H
#define STORMM_ATOMGRAPH_ABSTRACTS_H

#include "copyright.h"
#include "atomgraph_enumerators.h"

namespace stormm {
namespace topology {

/// \brief Unguarded struct to store the ingredients of a single Urey-Bradley interaction
template <typename T> struct UreyBradleyTerm {
  int i_atom;     ///< Atom i in the term
  int k_atom;     ///< Atom k in the term (there is no atom j in this harmonic two-body potential)
  int param_idx;  ///< Parameter index (for tracing purposes--critical values are given next)
  T keq;          ///< Harmonic stretching stiffness
  T leq;          ///< Harmonic equilibrium length (Angstroms)
};

/// \brief Unguarded struct to store the ingredients of a single CHARMM improper dihedral
template <typename T> struct CharmmImprTerm {
  int i_atom;     ///< First atom in the term
  int j_atom;     ///< Second (central) atom in the term
  int k_atom;     ///< Third atom in the term
  int l_atom;     ///< Fourth atom in the term
  int param_idx;  ///< Parameter index (for tracing purposes--critical values are given next)
  T keq;          ///< Harmonic planar angle stiffness
  T phi_eq;       ///< Harmonic equilibrium angle (radians)
};

/// \brief Unguarded struct to store the ingredients of a single CMAP interaction
template <typename T> struct CmapTerm {
  int i_atom;     ///< First atom in the first dihedral term
  int j_atom;     ///< Second atom in the first dihedral term, or first atom in the second
  int k_atom;     ///< Third atom in the first dihedral term, or second atom in the second
  int l_atom;     ///< Fourth atom in the first dihedral term, or third atom in the second
  int m_atom;     ///< Fourth atom in the second dihedral term
  int surf_idx;   ///< Surface term index (for tracing purposes--dimensions of a CMAP and a
                  ///<   pointer to it are given next)
  int surf_dim;   ///< Dimension of the (square, periodic) CMAP potential grid
  T* surf;        ///< Pointer to the surface potential values (stored in column-major format)
};

/// \brief Unguarded struct to store the ingredients of a single bond stretching interaction
template <typename T> struct BondTerm {
  int i_atom;     ///< Atom i in the term
  int j_atom;     ///< Atom j in the term
  int param_idx;  ///< Parameter index (for tracing purposes--critical values are given next)
  T keq;          ///< Harmonic stretching stiffness
  T leq;          ///< Harmonic equilibrium length (Angstroms)
};

/// \brief Unguarded struct to store the ingredients of a single bond angle interaction
template <typename T> struct AngleTerm {
  int i_atom;     ///< First atom in the term
  int j_atom;     ///< Second (central) atom in the term
  int k_atom;     ///< Third atom in the term
  int param_idx;  ///< Parameter index (for tracing purposes--critical values are given next)
  T keq;          ///< Harmonic angle stiffness
  T theta_eq;     ///< Harmonic equilibrium angle (radians)
};

/// \brief Unguarded struct to store the ingredients of a single dihedral (proper or improper)
///        torsion interaction.  The trefoil improper torsion scheme used by OpenMM is implemented
///        by applying multiple improper dihedrals of the Amber format, where the third atom is
///        taken to be the central atom.
template <typename T> struct DihedralTerm {
  int i_atom;     ///< First atom in the dihedral
  int j_atom;     ///< Second atom in the dihedral
  int k_atom;     ///< Third atom in the dihedral (central atom if this is an improper)
  int l_atom;     ///< Fourth atom in the dihedral
  int param_idx;  ///< Parameter index (for tracing purposes--critical values are given next)
  T amplitude;    ///< Cosine series amplitude
  T phase;        ///< Phase angle for the dihedral
  T periodicity;  ///< Periodicity of the dihedral (integral value, but given as a real number)
  T elec_screen;  ///< Electrostatic screening (electrostatic 1:4 interactions will be computed at
                  ///<   this intensity relative to their normal values)
  T vdw_screen;   ///< van-der Waals screening (van-der Waals 1:4 interactions will be computed at
                  ///<   this intensity relative to their normal values)
};

/// \brief Information need for bonded calculations.  Templating is used to serve either of two
///        levels of precision: single (which, on almost any machine and HPC accelerator card,
///        will be fp32) or double (fp64).  Parameters for bonds, the highest frequency terms in
///        any system by a factor of roughly 5, will always be stored in double precision, as will
///        parameters for bond angles, the next highest frequency group of terms.  The precision of
///        other parameters varies according to the template parameter.
template <typename T> struct ValenceKit {

  /// \brief The constructor takes a Fortran77-worthy list of arguments.  It would not be any
  ///        cleaner to make this object a friend of or nested struct within AtomGraph, due to its
  ///        templated nature.
  ///
  ///        See the descriptions of eponymous member variables (minus _in) within the ValenceKit
  ///        object for descriptions of each parameter.
  explicit ValenceKit(int natom_in, int nbond_in, int nangl_in, int ndihe_in, int nbond_param_in,
                      int nangl_param_in, int ndihe_param_in, int ninfr14_in, int nattn14_param_in,
                      int nubrd_in, int ncimp_in, int ncmap_in, int nubrd_param_in,
                      int ncimp_param_in, int ncmap_surf_in, const T* bond_keq_in,
                      const T* bond_leq_in, const T* angl_keq_in,
                      const T* angl_theta_in, const T* dihe_amp_in, const T* dihe_freq_in,
                      const T* dihe_phi_in, const T* attn14_elec_in, const T* attn14_vdw_in,
                      const int* bond_i_atoms_in, const int* bond_j_atoms_in,
                      const int* bond_param_idx_in, const char4* bond_modifiers_in,
                      const int* angl_i_atoms_in, const int* angl_j_atoms_in,
                      const int* angl_k_atoms_in, const int* angl_param_idx_in,
                      const char4* angl_modifiers_in, const int* dihe_i_atoms_in,
                      const int* dihe_j_atoms_in, const int* dihe_k_atoms_in,
                      const int* dihe_l_atoms_in, const int* dihe_param_idx_in,
                      const int* dihe14_param_idx_in, const char4* dihe_modifiers_in,
                      const int* infr14_i_atoms_in, const int* infr14_l_atoms_in,
                      const int* infr14_param_idx_in, const int* ubrd_i_atoms_in,
                      const int* ubrd_k_atoms_in, const int* ubrd_param_idx_in,
                      const int* cimp_i_atoms_in, const int* cimp_j_atoms_in,
                      const int* cimp_k_atoms_in, const int* cimp_l_atoms_in,
                      const int* cimp_param_idx_in, const int* cmap_i_atoms_in,
                      const int* cmap_j_atoms_in, const int* cmap_k_atoms_in,
                      const int* cmap_l_atoms_in, const int* cmap_m_atoms_in,
                      const int* cmap_dim_in, const int* cmap_surf_bounds_in,
                      const int* cmap_patch_bounds_in, const int* cmap_surf_idx_in,
                      const T* ubrd_keq_in, const T* ubrd_leq_in, const T* cimp_keq_in,
                      const T* cimp_phi_in, const T* cmap_surf_in, const T* cmap_dphi_in,
                      const T* cmap_dpsi_in, const T* cmap_dphi_dpsi_in, const T* cmap_patches_in,
                      const int* bond_asgn_atoms_in, const int* bond_asgn_index_in,
                      const int* bond_asgn_terms_in, const int* bond_asgn_bounds_in,
                      const int* angl_asgn_atoms_in, const int* angl_asgn_index_in,
                      const int* angl_asgn_terms_in, const int* angl_asgn_bounds_in,
                      const int* dihe_asgn_atoms_in, const int* dihe_asgn_index_in,
                      const int* dihe_asgn_terms_in, const int* dihe_asgn_bounds_in,
                      const int* ubrd_asgn_atoms_in, const int* ubrd_asgn_index_in,
                      const int* ubrd_asgn_terms_in, const int* ubrd_asgn_bounds_in,
                      const int* cimp_asgn_atoms_in, const int* cimp_asgn_index_in,
                      const int* cimp_asgn_terms_in, const int* cimp_asgn_bounds_in,
                      const int* cmap_asgn_atoms_in, const int* cmap_asgn_index_in,
                      const int* cmap_asgn_terms_in, const int* cmap_asgn_bounds_in);

  /// \brief Take the default copy and move constructors.  The move assignment operator will get
  ///        implicitly deleted as this is just a collection of constants.
  /// \{
  ValenceKit(const ValenceKit &original) = default;
  ValenceKit(ValenceKit &&other) = default;
  /// \}
  
  // The purpose of this struct is to store a collection of pointers for HPC kernels.  As such, it
  // does not have any private member variables.  As a provider of parameters, it also does not
  // allow modification of the data that it points to.
  const int natom;              ///< The number of atoms in the system (needed for the overall
                                ///<   length of bond_asgn_abounds and similar arrays towards the
                                ///<   end, if not just for general convenience)
  const int nbond;              ///< Number of bonds in the system
  const int nangl;              ///< Number of bond angles in the system
  const int ndihe;              ///< Number of dihedrals or torsions (includes cosine-based
                                ///<   improper torsions) in the system
  const int nbond_param;        ///< Number of unique bond parameters
  const int nangl_param;        ///< Number of unqiue bond angle parameters
  const int ndihe_param;        ///< Number of unique cosine dihedral parameters
  const int ninfr14;            ///< The number of inferred 1:4 exclusions
  const int nattn14_param;      ///< Number of 1:4 attenuation factor pairs (vdW, electrostatic)
  const int nubrd;              ///< Number of Urey-Bradley 1:3 harmonic angle interactions
  const int ncimp;              ///< Number of CHARMM harmonic improper interactions
  const int ncmap;              ///< Number of CHARMM CMAP surface spline-based interactions
  const int nubrd_param;        ///< Number of unique Urey-Bradley parameters
  const int ncimp_param;        ///< Number of unique CHARMM harmonic improper parameters
  const int ncmap_surf;         ///< Number of unique CMAP surfaces
  const T* bond_keq;            ///< Equilibrium stiffness constants for all unique bonds
  const T* bond_leq;            ///< Equilibrium lengths for all unique bonds
  const T* angl_keq;            ///< Equilibrium stiffness constants for all unique bond angles
  const T* angl_theta;          ///< Equilibrium lengths for all unique bond angles
  const T* dihe_amp;            ///< Amplitudes for all unique cosine-based dihedrals
  const T* dihe_freq;           ///< Periodicities for cosine-based dihedral / torsion terms
  const T* dihe_phi;            ///< Phase angles for cosine-based dihedral / torsion terms
  const T* attn14_elec;         ///< Electrostatic 1:4 attenuation parameters, including a zero for
                                ///<   no 1:4 interaction.  The arrays dihe14_param_idx and
                                ///<   infr14_param_idx both index into this array.
  const T* attn14_vdw;          ///< van-der Waals 1:4 attenuation parameters, including a zero for
                                ///<   no 1:4 interaction.  The arrays dihe14_param_idx and
                                ///<   infr14_param_idx both index into this array.
  const int* bond_i_atoms;      ///< Array of first atoms in each bond
  const int* bond_j_atoms;      ///< Array of second atoms in each bond
  const int* bond_param_idx;    ///< Parameter indices for each bond in the system
  const char4* bond_modifiers;  ///< Modifying details of each bond stretching term
  const int* angl_i_atoms;      ///< Array of first atoms in each bond angle
  const int* angl_j_atoms;      ///< Array of second atoms in each bond angle
  const int* angl_k_atoms;      ///< Array of third atoms in each bond angle
  const int* angl_param_idx;    ///< Parameter indices for each bond angle in the system
  const char4* angl_modifiers;  ///< Modifying details of each angle bending term
  const int* dihe_i_atoms;      ///< Array of first atoms in each cosine-based dihedral
  const int* dihe_j_atoms;      ///< Array of second atoms in each cosine-based dihedral
  const int* dihe_k_atoms;      ///< Array of third atoms in each cosine-based dihedral (these are
                                ///<   the center atoms if the dihedral describes a cosine-based
                                ///<   improper)
  const int* dihe_l_atoms;      ///< Array of fourth atoms in each cosine-based dihedral
  const int* dihe_param_idx;    ///< Parameter indices for each cosine-based dihedral in the system
  const int* dihe14_param_idx;  ///< Parameter indices for 1:4 attenuated interactions handled by
                                ///<   each dihedral through their I and L atoms
  const char4* dihe_modifiers;  ///< Modifying details of each dihedral term, including whether it
                                ///<   is a proper or improper dihedral with or without a 1:4
                                ///<   attenuated interaction
  const int* infr14_i_atoms;    ///< Inferred 1:4 interactions' I atoms.  These cannot be handled
                                ///<   directly by any dihedrals and so must be treated as separate
                                ///<   terms.
  const int* infr14_l_atoms;    ///< Inferred 1:4 interactions' J atoms
  const int* infr14_param_idx;  ///< Indices into the arrays of 1:4 attenuation parameters for the
                                ///<   set of 1:4 interactions with inferred parameters
  const int* ubrd_i_atoms;      ///< First atoms in each Urey-Bradley harmonic angle interaction
  const int* ubrd_k_atoms;      ///< Second atoms in each Urey-Bradley harmonic angle interaction
  const int* ubrd_param_idx;    ///< Parameter indices for each Urey-Bradley interaction
  const int* cimp_i_atoms;      ///< First atoms in each CHARMM harmonic improper dihedral
  const int* cimp_j_atoms;      ///< Second atoms in each CHARMM harmonic improper dihedral
  const int* cimp_k_atoms;      ///< Third (center) atoms in each CHARMM harmonic improper dihedral
  const int* cimp_l_atoms;      ///< Fourth atoms in each CHARMM harmonic improper dihedral
  const int* cimp_param_idx;    ///< Parameter indices for each CHARMM harmonic improper dihedral
  const int* cmap_i_atoms;      ///< First atoms of the first dihedral in each CHARMM CMAP term
  const int* cmap_j_atoms;      ///< Second atoms of the first dihedral, first atoms of the second
                                ///<   dihedral, in each CHARMM CMAP term
  const int* cmap_k_atoms;      ///< Third atoms of the first dihedral, second atoms of the second
                                ///<   dihedral, in each CHARMM CMAP term
  const int* cmap_l_atoms;      ///< Fourth atoms of the first dihedral, third atoms of the second
                                ///<   dihedral, in each CHARMM CMAP term
  const int* cmap_m_atoms;      ///< Fourth atoms of the second dihedral in each CHARMM CMAP term
  const int* cmap_dim;          ///< Dimension of each CMAP surface (each surface is a square grid)
  const int* cmap_surf_bounds;  ///< Bounds array for individual CMAP surfaces
  const int* cmap_patch_bounds; ///< Bounds array for individual CMAP surfaces composed of patch
                                ///<   matrices (this array is basically 16 times cmap_surf_bounds)
  const int* cmap_surf_idx;     ///< Parameter (surface) indices for CHARMM CMAP interactions
  const T* ubrd_keq;            ///< Array of unique Urey-Bradley 1:3 harmonic stiffnesses
  const T* ubrd_leq;            ///< Array of unique Urey-Bradley 1:3 harmonic equilibrium lengths
  const T* cimp_keq;            ///< Array of unique CHARMM harmonic improper dihedral stiffnesses
  const T* cimp_phi;            ///< Array of unique CHARMM harmonic improper dihedral phase angles
  const T* cmap_surf;           ///< Array of CHARMM CMAP surface values
  const T* cmap_dphi;           ///< Array of CHARMM CMAP surface "phi" derivative values
  const T* cmap_dpsi;           ///< Array of CHARMM CMAP surface "psi" derivative values
  const T* cmap_dphi_dpsi;      ///< Array of CHARMM CMAP surface cross derivative values
  const T* cmap_patches;        ///< Array of interlaced CHARMM CMAP patch matrices

  // The following arrays provide a different perspective on the valence terms.  Each bond, angle,
  // dihedral, Urey-Bradley, CMAP, and other valence term is assigned to one of its constituent
  // atoms.  Atoms are then responsible for executing any of their assigned terms.  Atom a_idx
  // is responsible for terms in each array indexed by [(term)_asgn_bounds[a_idx] ...
  // (term)_asgn_bounds[a_idx] + 1), which provides the relevant elements of the associated
  // (term)_asgn_index array (for the parameter index) or the (term)_assigned_atoms array
  // (multiply by (number of atoms in the term - 1), but this lists all other atoms in the term).
  // This allows instant lookup of all relevant bonds, angles, CMAPs, etc. based on any given atom
  // number, although to get all bonds in a given molecule one still has to loop over all atoms in
  // the molecule.  To find all bonds associated with a given atom a_idx, look for all bonds that
  // a_idx itself is responsible for and then cycle through its 1:2 exclusions neighbor list to see
  // which of those atom are responsible for bonds back to a_idx.  That will enumerate all bonds,
  // with parameters, affecting atom a_idx.
  const int* bond_asgn_atoms;  ///< Other atoms for bond terms (one atom per term)
  const int* bond_asgn_index;  ///< Parameter indices for bond terms
  const int* bond_asgn_terms;  ///< Bond term indices assigned to each atom (i.e. the 5th atom is
                               ///<   assigned bond terms with indices 18, 20, and 21)
  const int* bond_asgn_bounds; ///< Bounds for each atom's assigned bond list
  const int* angl_asgn_atoms;  ///< Other atoms for bond angle terms (two atoms per term)
  const int* angl_asgn_index;  ///< Parameter indices for bond angle terms
  const int* angl_asgn_terms;  ///< Angle term indices assigned to each atom (i.e. the 9th atom is
                               ///<   assigned angle terms with indices 46, 48, and 50)
  const int* angl_asgn_bounds; ///< Bounds for each atom's assigned bond angle list
  const int* dihe_asgn_atoms;  ///< Other atoms for dihedral terms (three atoms per term)
  const int* dihe_asgn_index;  ///< Parameter indices for dihedral terms
  const int* dihe_asgn_terms;  ///< Dihedral term indices assigned to each atom (i.e. the 4th atom
                               ///<   is assigned dihedral terms with indices 27, 28, and 29)
  const int* dihe_asgn_bounds; ///< Bounds for each atom's assigned dihedral list
  const int* ubrd_asgn_atoms;  ///< Other atoms for Urey-Bradley terms (one atom per term)
  const int* ubrd_asgn_index;  ///< Parameter indices for Urey-Bradley terms
  const int* ubrd_asgn_terms;  ///< Urey-Bradley term indices assigned to each atom
  const int* ubrd_asgn_bounds; ///< Bounds for each atom's assigned Urey-Bradley list
  const int* cimp_asgn_atoms;  ///< Other atoms for CHARMM impropers (three atoms per term)
  const int* cimp_asgn_index;  ///< Parameter indices for CHARMM improper dihedral terms
  const int* cimp_asgn_terms;  ///< CHARMM improper term indices assigned to each atom
  const int* cimp_asgn_bounds; ///< Bounds for each atom's assigned CHARMM improper list
  const int* cmap_asgn_atoms;  ///< Other atoms for CMAP terms (four atoms per term)
  const int* cmap_asgn_index;  ///< Parameter indices for CMAP terms
  const int* cmap_asgn_terms;  ///< CMAP term indices assigned to each atom
  const int* cmap_asgn_bounds; ///< Bounds for each atom's assigned CMAP list
};

/// \brief Information needed for non-bonded real-space calculations.  Templating is used as above,
///        to provide different levels of precision in the real number representation.
template <typename T> struct NonbondedKit {

  /// \brief As with most other astracts, the constructor is the only member function of
  ///        significance.  It takes a long list of arguments one for each of its member variables.
  explicit NonbondedKit(int natom_in, int n_q_types_in, int n_lj_types_in,
                        const T coulomb_constant_in, const T* charge_in, const int* q_idx_in,
                        const int* lj_idx_in, const T* q_parameter_in, const T* lja_coeff_in,
                        const T* ljb_coeff_in, const T* ljc_coeff_in, const T* lja_14_coeff_in,
                        const T* ljb_14_coeff_in, const T* ljc_14_coeff_in, const T* lj_sigma,
                        const T* lj_14_sigma, const int* nb11x_in, const int* nb11_bounds_in,
                        const int* nb12x_in, const int* nb12_bounds_in, const int* nb13x_in,
                        const int* nb13_bounds_in, const int* nb14x_in, const int* nb14_bounds_in,
                        const T* lj_type_corr_in);

  /// \brief Take the default copy and move constructors.  The assignment operators will get
  ///        implicitly deleted as this is just a collection of constants.
  /// \{
  NonbondedKit(const NonbondedKit &original) = default;
  NonbondedKit(NonbondedKit &&other) = default;
  /// \}

  // Member variables again store a collection of atomic parameters
  const int natom;          ///< The number of atoms in the system
  const int n_q_types;      ///< The number of unique charge types in the system
  const int n_lj_types;     ///< The number of unique Lennard-Jones atom types in the system
  const T coulomb_constant; ///< Coulomb's constant in units of kcal-A/mol-e^2
  const T* charge;          ///< Partial atomic charges on all atoms
  const int* q_idx;         ///< Partial charge type indices for all atoms
  const int* lj_idx;        ///< Lennard-Jones type indices of all atoms
  const T* q_parameter;     ///< Partial atomic charges for each charge type (this will almost
                            ///<   certainly be smaller than the array of charges for every atom,
                            ///<   as there are only about 350 unique chargees in a protein force
                            ///<   field and two in most water models).  The representation here is
                            ///<   not more memory-efficient than accessing the single-precision
                            ///<   charge data, however, as the 32-bit index is a memory access of
                            ///<   its own.  However, in the AtomGraphSynthesis an array of
                            ///<   bit-packed unsigned integers that delivers both the charge and
                            ///<   Lennard-Jones parameter indices for each atom will offer the
                            ///<   most performant solution.
  const T* lja_coeff;       ///< Lennard-Jones A coefficients for all atoms, tabulated as a square
                            ///<   matrix of rank (number of Lennard-Jones types)
  const T* ljb_coeff;       ///< Lennard_jones B coefficients, again tabulated as a square matrix
  const T* ljc_coeff;       ///< Lennard_jones C coefficients, again tabulated as a square matrix
  const T* lja_14_coeff;    ///< Lennard-Jones A coefficients for 1:4 interactions, tabulated as a
                            ///<   square matrix in the same form as lja_coeff
  const T* ljb_14_coeff;    ///< Lennard_jones B coefficients for 1:4 interactions
  const T* ljc_14_coeff;    ///< Lennard_jones C coefficients for 1:4 interactions
  const T* lj_sigma;        ///< Lennard_jones sigma parameters
  const T* lj_14_sigma;     ///< Lennard_jones sigma parameters for 1:4 interactions
  const int* nb11x;         ///< Non-bonded 1:1 exclusions, applicable to virtual sites
  const int* nb11_bounds;   ///< Non-bonded 1:1 exclusion array bounds for each atom
  const int* nb12x;         ///< Non-bonded 1:2 exclusions, applicable to bonded atoms
  const int* nb12_bounds;   ///< Non-bonded 1:2 exclusion array bounds for each atom
  const int* nb13x;         ///< Non-bonded 1:3 exclusions, applicable to atoms in bond angles and
                            ///<   Urey-Bradley terms
  const int* nb13_bounds;   ///< Non-bonded 1:3 exclusion array bounds for each atom
  const int* nb14x;         ///< Non-bonded 1:4 exclusions, applicable to atoms in proper torsions
  const int* nb14_bounds;   ///< Non-bonded 1:4 exclusion array bounds for each atom
  const T* lj_type_corr;    ///< Lennard-Jones energy corrections for each atom type (used only
                            ///<   in energy calculations)
};

/// \brief Information needed for Generalized Born (and perhaps other) implicit solvent methods.
///        This information is collected into an object separate from the non-bonded kit because
///        only a subset of calculations will use GB for the solvent conditions.
template <typename T> struct ImplicitSolventKit {

  explicit ImplicitSolventKit(int natom_in, ImplicitSolventModel igb_in, T dielectric_in,
                              T saltcon_in, const int* neck_gb_idx_in, const T* pb_radii_in,
                              const T* gb_screen_in, const T* gb_alpha_in,
                              const T* gb_beta_in, const T* gb_gamma_in);

  /// \brief Take the default copy and move constructors as well as assignment operators
  /// \{
  ImplicitSolventKit(const ImplicitSolventKit &original) = default;
  ImplicitSolventKit(ImplicitSolventKit &&other) = default;
  /// \}

  // Member variables again store a collection of atomic parameters
  const int natom;                ///< The number of atoms in the system
  const ImplicitSolventModel igb; ///< The flavor of Generalized Born to use, i.e. Hawkins /
                                  ///<   Cramer / Truhlar
  const T dielectric;             ///< The dielectric constant to use (default 80.0)
  const T saltcon;                ///< The salt concentration to use in GB calculations
  const int* neck_gb_idx;         ///< Neck GB indicies for all atoms, applicable for Mongan's
                                  ///<   "neck" GB models
  const T* pb_radii;              ///< Atomic PB radii
  const T* gb_screen;             ///< Generalized Born screening factors
  const T* gb_alpha;              ///< Generalized born alpha parameters (one parameter per atom)
  const T* gb_beta;               ///< Generalized born beta parameters (one parameter per atom)
  const T* gb_gamma;              ///< Generalized born gamma parameters (one parameter per atom)
};

/// \brief Information on atoms and residues which may be useful for applying atom masks or
///        identifying specific parts of the sytem
struct ChemicalDetailsKit {

  /// \brief Simple constructor based on many detais.  This, like other abstracts, will most likely
  ///        be produced by some AtomGraph member function.
  ChemicalDetailsKit(int natom_in, int nres_in, int nmol_in, int free_dof_in, int cnst_dof_in,
                     const char4* atom_names_in, const char4* res_names_in,
                     const char4* atom_types_in, const int* z_numbers_in, const int* res_limits_in,
                     const int* atom_numbers_in, const int* res_numbers_in, const int* mol_home_in,
                     const int* mol_contents_in, const int* mol_limits_in, const double* masses_in,
                     const float* sp_masses_in, const double* inv_masses_in,
                     const float* sp_inv_masses_in);

  /// \brief Take the default copy and move constructors.  The assignment operators will get
  ///        implicitly deleted as this is just a collection of constants.
  /// \{
  ChemicalDetailsKit(const ChemicalDetailsKit &original) = default;
  ChemicalDetailsKit(ChemicalDetailsKit &&other) = default;
  /// \}

  // Member variables store the atom count, residue count, and other ways to quantify the system
  // in addition to many useful pointers
  const int natom;             ///< The number of atoms in the system
  const int nres;              ///< The number of residues in the system
  const int nmol;              ///< The number of molecules in the system
  const int free_dof;          ///< Number of degrees of freedom, without consideration to
                               ///<   geometric constraints
  const int cnst_dof;          ///< Number of degrees of freedom, when geometric constraints are
                               ///<   in effect
  const char4* atom_names;     ///< Names of all atoms in the system
  const char4* res_names;      ///< Names of all residues in the system
  const char4* atom_types;     ///< Atom type names for all atoms in the system
  const int* z_numbers;        ///< Atomic numbers for all atoms in the system
  const int* res_limits;       ///< Residue limits, a capped (exclusive) prefix sum
  const int* atom_numbers;     ///< Structural atom numbers for every atom
  const int* res_numbers;      ///< Structural residue numbers for every atom (atom, not residue,
                               ///<   even though the numbers refer to residues)
  const int* mol_home;         ///< Molecule index to which each atom belongs
  const int* mol_contents;     ///< Contents of every molecule in the system, as lists of atom
                               ///<   indices
  const int* mol_limits;       ///< Molecule limits, the bounds by which to read mol_contents
  const double* masses;        ///< Masses of atoms in the system
  const float* sp_masses;      ///< Masses of atoms in the system (single precision)
  const double* inv_masses;    ///< Masses of atoms in the system
  const float* sp_inv_masses;  ///< Masses of atoms in the system (single precision)
};

/// \brief Information needed for the placement of virtual sites and transmission of forces on
///        these sites to their frame atoms which have mass.
template <typename T> struct VirtualSiteKit {

  /// \brief The constructor follows other abstracts and is produced based on pointers from an
  ///        AtomGraph object.
  explicit VirtualSiteKit(int nsite_in, int nframe_set_in, const int* vs_atoms_in,
                          const int* frame1_idx_in, const int* frame2_idx_in,
                          const int* frame3_idx_in, const int* frame4_idx_in,
                          const int* vs_param_idx_in, const int* vs_types_in, const T* dim1_in,
                          const T* dim2_in, const T* dim3_in);

  /// \brief Take the default copy and move constructors.  The assignment operators will get
  ///        implicitly deleted as this is just a collection of constants.
  /// \{
  VirtualSiteKit(const VirtualSiteKit &original) = default;
  VirtualSiteKit(VirtualSiteKit &&other) = default;
  /// \}

  const int nsite;          ///< The number of virtual sites in the topology
  const int nframe_set;     ///< The number of unique frame parameter sets in the topology
  const int* vs_atoms;      ///< Topological indicies of frame atoms
  const int* frame1_idx;    ///< Topological indices of frame atom 1 (the parent atom)
  const int* frame2_idx;    ///< Topological indices of frame atom 2
  const int* frame3_idx;    ///< Topological indices of frame atom 3
  const int* frame4_idx;    ///< Topological indices of frame atom 4
  const int* vs_param_idx;  ///< Parameter indices for each virtual site, offering the correct
                            ///<   elements of vs_types, dim1, dim2, and dim3 to access in order
                            ///<   to understand the frame
  const int* vs_types;      ///< Virtual site frame types
  const T* dim1;            ///< Frame dimension 1
  const T* dim2;            ///< Frame dimension 2
  const T* dim3;            ///< Frame dimension 3
};

/// \brief Information needed to manage constraint groups.  This additional abstract is needed
///        due to the way that some collections of rigid bonds all connect to the same atom, and
///        because analytic SETTLE constraints are a thing.
template <typename T> struct ConstraintKit {

  /// \brief The constructor follows other abstracts and is produced based on pointers from an
  ///        AtomGraph object.
  explicit ConstraintKit(const int nsettle_in, const int nsett_param_in, const int ngroup_in,
                         const int ncnst_param_in, const int* settle_ox_atoms_in,
                         const int* settle_h1_atoms_in, const int* settle_h2_atoms_in,
                         const int* settle_param_idx_in, const int* group_list_in,
                         const int* group_bounds_in, const int* group_param_idx_in,
                         const int* group_param_bounds_in, const T* settle_mormt_in,
                         const T* settle_mhrmt_in, const T* settle_ra_in, const T* settle_rb_in,
                         const T* settle_rc_in, const T* settle_invra_in,
                         const T* group_sq_lengths_in, const T* group_inv_masses_in);

  /// \brief Take the default copy and move constructors.  The assignment operators will get
  ///        implicitly deleted as this is just a collection of constants.
  /// \{
  ConstraintKit(const ConstraintKit &original) = default;
  ConstraintKit(ConstraintKit &&other) = default;
  /// \}
  
  const int nsettle;              ///< Number of SETTLE (analytic) constrained groups of bonds.
  const int nsett_param;          ///< The number of distinct SETTLE parameter sets
  const int ngroup;               ///< Number of "hub and spoke" constrained groups of bonds.
  const int ncnst_param;          ///< The number of distinct constraint group parameter sets
  const int* settle_ox_atoms;     ///< Central 'oxygen' atoms for analytic SETTLE constraints
  const int* settle_h1_atoms;     ///< First 'hydrogen' atoms for analytic SETTLE constraints
  const int* settle_h2_atoms;     ///< Second 'hydrogen' atoms for analytic SETTLE constraints
  const int* settle_param_idx;    ///< Parameter set indices for each SETLLE group
  const int* group_list;          ///< List of all atoms involved in "hub and spoke" constraint
                                  ///<   groups. In each group, the central atom, to which all
                                  ///<   others bind, is listed first.  It is the first atom of
                                  ///<   any of the constrained bonds.  All other atoms are the
                                  ///<   distal termini of those constrained bonds.
  const int* group_bounds;        ///< Bounds array for group_list above.  Every segment of
                                  ///<   group_list will be at least two atoms, and generally two
                                  ///<   to five atoms.
  const int* group_param_idx;     ///< Parameters array for each constraint group.  To reference
                                  ///<   parameter set k is an instruction to access indices k and
                                  ///<   k + 1 from the group_param_bounds array, read the
                                  ///<   group_lengths and group_inv_masses arrays between those
                                  ///<   limits, and apply the lengths and inverse masses to each
                                  ///<   particle in the group.
  const int* group_param_bounds;  ///< Bounds array for the constraint group parameter arrays below
  const T* settle_mormt;          ///< Proportional mass of "oxygen" in SETTLE systems
  const T* settle_mhrmt;          ///< Proportional mass of "hydrogen" in SETTLE systems
  const T* settle_ra;             ///< Internal distance measurement of SETTLE groups
  const T* settle_rb;             ///< Internal distance measurement of SETTLE groups
  const T* settle_rc;             ///< Internal distance measurement of SETTLE groups
  const T* settle_invra;          ///< Internal distance measurement of SETTLE groups
  const T* group_sq_lengths;      ///< Bond length targets for the group atoms
  const T* group_inv_masses;      ///< Inverse masses for the particles in each group
};

} // namespace topology
} // namespace stormm

#include "atomgraph_abstracts.tpp"

#endif
