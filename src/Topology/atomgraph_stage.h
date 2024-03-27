// -*-c++-*-
#ifndef STORMM_ATOMGRAPH_STAGE_H
#define STORMM_ATOMGRAPH_STAGE_H

#include <vector>
#include <string>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "MoleculeFormat/mdlmol.h"
#include "atomgraph.h"
#include "atomgraph_abstracts.h"
#include "atomgraph_enumerators.h"

namespace stormm {
namespace topology {

using structure::MdlMol;

/// \brief Default settings associated with topology staging
/// \{
constexpr int default_minimum_solute_atoms = 15;
constexpr char4 wildcard_atom_type = { 'W', 'I', 'L', 'D' };
constexpr char4 unknown_atom_name = { 'U', 'N', 'K', ' ' };
constexpr char4 unknown_residue_name = { 'U', 'N', 'K', ' ' };
constexpr char4 unknown_atom_type = { 'U', 'N', 'K', ' ' };
/// \}

/// \brief A struct to encode the non-bonded atom types in a dihedral and associated electrostatic
///        and van-der Waals 1:4 interactions.
struct ScreeningFactorTag {

  /// \brief The simplest constructor labels the middle two atom types as "wildcard" and the end
  ///        types as "unknown," with 1.0 screening factors.
  /// \{
  ScreeningFactorTag();
  ScreeningFactorTag(int dihe_idx_in, double elec14_in, double vdw14_in);
  ScreeningFactorTag(const char4 i_at_in, const char4 l_at_in, double elec14_in, double vdw14_in);
  ScreeningFactorTag(const char4 i_at_in, const char4 j_at_in, const char4 k_at_in,
                     const char4 l_at_in, double elec14_in, double vdw14_in);
  /// \}

  int dihe_idx;   ///< The dihedral index to which this scaling factor pair applies.  If greater
                  ///<   than or equal to zero, this will take precedence over the atom types in
                  ///<   determining how to apply the scaling factors.  If less that zero, the atom
                  ///<   type codes will take precedence.
  char4 i_at;     ///< Atom type of the Ith atom in the proper dihedral (the dihedral will be
                  ///<   compared in both orders { I, J, K, L } and { L, K, J, I }) 
  char4 j_at;     ///< Atom type of the Jth atom in the proper dihedral
  char4 k_at;     ///< Atom type of the Kth atom in the proper dihedral
  char4 l_at;     ///< Atom type of atom L in the proper dihedral
  double elec14;  ///< Electrostatic 1:4 screening factor
  double vdw14;   ///< Van-der Waals 1:4 screening factor
};
  
/// \brief The topology stage is a CPU--bound object composed of Standard Template Library vectors
///        for its data arrays, with each member corresponding to a (probably eponymous) member of
///        the AtomGraph class.  This object allows a topology to be built atom by atom,
///        deconstructed, rebuilt, and even exported to a (much less mutable) AtomGraph.
class AtomGraphStage {
public:

  /// \brief The constructor can accept a number of atoms, residues, and molecules, an existing
  ///        AtomGraph, or a ChemicalFeatures object.
  /// \{
  AtomGraphStage(int atom_count_in = 0, const std::vector<int> &residue_limits_in = {},
                 ExceptionResponse policy_in = ExceptionResponse::DIE);
  AtomGraphStage(const AtomGraph *ag_in, ExceptionResponse policy_in = ExceptionResponse::DIE);
  AtomGraphStage(const AtomGraph &ag_in, ExceptionResponse policy_in = ExceptionResponse::DIE);
  AtomGraphStage(const MdlMol *cmpd_in, ExceptionResponse policy_in = ExceptionResponse::DIE);
  AtomGraphStage(const MdlMol &cmpd_in, ExceptionResponse policy_in = ExceptionResponse::DIE);
  /// \}

  /// \brief The utility of this object lies in its simple construction and independent arrays,
  ///        which can be resized as needed.  This implies no pointers to repair and no
  ///        POINTER-kind Hybrid objects which could require special copy or move constructors and
  ///        assignment operators.  All of the defaults apply.
  ///
  /// \param original  The object to copy or move
  /// \param other     An existing object placed on the left hand side of the assignment statement
  /// \{
  AtomGraphStage(const AtomGraphStage &original) = default;
  AtomGraphStage(AtomGraphStage &&original) = default;
  AtomGraphStage& operator=(const AtomGraphStage &other) = default;
  AtomGraphStage& operator=(AtomGraphStage &&other) = default;
  /// \}

  /// \brief Get the title of the topology as it is staged.
  const std::string& getTitle() const;

  /// \brief Get the name of the output file to be referenced by the resulting topology.  Printed
  ///        results will go to this file.
  const std::string& getOutputFile() const;
  
  /// \brief Get the total number of atoms.
  int getAtomCount() const;

  /// \brief Get the total number of residues.
  int getResidueCount() const;

  /// \brief Get the total number of molecules.
  int getMoleculeCount() const;

  /// \brief Get the general 1:4 screening factor on electrostatic interactions.
  double getElectrostatic14Screening() const;
  
  /// \brief Get the general 1:4 screening factor on van-der Waals interactions.
  double getVanDerWaals14Screening() const;
  
  /// \brief Export the topology based on the internal data arrays.  A final list of checks and
  ///        fixes will be applied to ensure consistent parameterization.
  ///
  /// \param minimum_solute_size             The minimum size of a solute molecule.  Anything
  ///                                        larger will be automatically included in the solute.
  /// \param solute_included_residues_names  Names of additional residues that should be included
  ///                                        in the solute
  /// \param hmass_repartition_factor        Hydrogen mass repartitioning factor.  The mass of
  ///                                        hydrogen atoms in the resulting topology will be
  ///                                        scaled by this amount.  The default of 1.0 implies no
  ///                                        repartitioning.  Simulations seeking to use a 4fs time
  ///                                        step typically repartition with a factor of 3.0.  The
  ///                                        extra mass added to hydrogen atoms will be borrowed
  ///                                        from their parent heavy atoms.
  AtomGraph exportTopology(int minimum_solute_size = default_minimum_solute_atoms,
                           const std::vector<char4> &solute_included_residue_names = {},
                           double hmass_repartition_factor = 1.0) const;

  /// \brief Set the title of the topology.
  ///
  /// \param title_in  The title to apply (will be checked for length, max 80 characters)
  void setTitle(const std::string &title_in);

  /// \brief Set the output file.
  ///
  /// \param output_file_in  Name of the output file to associate with the topology, the file to
  ///                        write
  void setOutputFile(const std::string &output_file_in);

  /// \brief Add new Urey-Bradley stretching parameters, accessible either by an index determined
  ///        by the order of addition or non-bonded atom types of the participating atoms.
  ///
  /// \param k_eq_in      Force constant (stiffness) for the harmonic bond interaction, in units of
  ///                     kcal/mol-A^2
  /// \param l_eq_in      Equilibrium length of the bond, in Angstroms
  /// \param atom_i_type  Non-bonded atom type of the first atom in the term.  The atom types are
  ///                     reflexive, a Urey-Bradley stretching term between types tI and tK being
  ///                     equivalent to a term between types tK and tI.
  /// \param atom_k_type  Non-bonded atom type of the second atom in the bond
  void addUbrdParameters(double k_eq_in, double l_eq_in,
                         const char4 atom_i_type = unknown_atom_type,
                         const char4 atom_k_type = unknown_atom_type);

  /// \brief Add new parameters for a CHARMM improper dihedral, accessible either by an index
  ///        determined by the order of addition or non-bonded atom types of the participating
  ///        atoms.
  ///
  /// \param ampl_in      Force constant (amplitude) for the harmonic improper torsion term, in
  ///                     units of kcal/mol-rad^2
  /// \param phi_in       Phase angle for the improper term, in unis of radians
  /// \param atom_i_type  Non-bonded atom type of the first atom in the improper torsion.  The atom
  ///                     types are again reflexive, a harmonic improper spanning types tI, tJ, tK,
  ///                     and tL being equivalent to an angle between types tL, tK, tJ, and tI.
  /// \param atom_j_type  Non-bonded atom type of the second atom in the improper torsion
  /// \param atom_k_type  Non-bonded atom type of the third atom in the improper torsion
  /// \param atom_l_type  Non-bonded atom type of the fourth atom in the improper torsion
  void addCImpParameters(double k_eq_in, double phi_in,
                         const char4 atom_i_type = unknown_atom_type,
                         const char4 atom_j_type = unknown_atom_type,
                         const char4 atom_k_type = unknown_atom_type,
                         const char4 atom_l_type = unknown_atom_type);

  /// \brief Add new parameters for a CHARMM correction map two-dimensional surface term,
  ///        accessible either by an index determined by the order of addition or non-bonded atom
  ///        types of the participating atoms.
  ///
  /// \param dimension_in  Dimension of the two-dimensional surface, assumed to be square and
  ///                      periodic on the range [-pi, pi)
  /// \param surf_in       Details of the two-dimensional energy surface, which will be solved by
  ///                      Catmull-Rom splines wherein the first and second cross-derivatives at
  ///                      each grid point are set to be equal between segements which exactly
  ///                      interpolate the potential values as provided.  Units of kcal/mol.
  /// \param atom_i_type   Non-bonded atom type of the first atom in the CMAP.  The atom types are
  ///                      again reflexive, a harmonic improper spanning types tI, tJ, tK, tL, and
  ///                      tM being equivalent to an angle between types tM, tL, tK, tJ, and tI.
  /// \param atom_j_type   Non-bonded atom type of the second atom in the CMAP
  /// \param atom_k_type   Non-bonded atom type of the third atom in the CMAP
  /// \param atom_l_type   Non-bonded atom type of the fourth atom in the CMAP
  /// \param atom_m_type   Non-bonded atom type of the fifth atom in the CMAP
  void addCmapParameters(const int dimension_in, const std::vector<double> &surf_in,
                         const char4 atom_i_type = unknown_atom_type,
                         const char4 atom_j_type = unknown_atom_type,
                         const char4 atom_k_type = unknown_atom_type,
                         const char4 atom_l_type = unknown_atom_type,
                         const char4 atom_m_type = unknown_atom_type);

  /// \brief Add new bond parameters, accessible either by an index determined by the order of
  ///        addition or non-bonded atom types of the participating atoms.
  ///
  /// \param k_eq_in      Force constant (stiffness) for the harmonic bond interaction, in units of
  ///                     kcal/mol-A^2
  /// \param l_eq_in      Equilibrium length of the bond, in Angstroms
  /// \param atom_i_type  Non-bonded atom type of the first atom in the bond.  The atom types are
  ///                     reflexive, a bond between types tI and tJ being equivalent to a bond
  ///                     between types tJ and tI.
  /// \param atom_j_type  Non-bonded atom type of the second atom in the bond
  void addBondParameters(double k_eq_in, double l_eq_in,
                         const char4 atom_i_type = unknown_atom_type,
                         const char4 atom_j_type = unknown_atom_type);

  /// \brief Add new angle parameters, accessible either by an index determined by the order of
  ///        addition or non-bonded atom types of the participating atoms.
  ///
  /// \param k_eq_in      Force constant (stiffness) for the harmonic bond angle interaction, in
  ///                     units of kcal/mol-radian^2
  /// \param th_eq_in      Equilibrium angle, in radians
  /// \param atom_i_type  Non-bonded atom type of the first atom in the angle.  The atom types are
  ///                     reflexive, an angle spanning types tI, tJ, and tK being equivalent to an
  ///                     angle between types tK, tJ, and tI.
  /// \param atom_j_type  Non-bonded atom type of the second atom in the angle
  /// \param atom_k_type  Non-bonded atom type of the third atom in the angle
  void addAnglParameters(double k_eq_in, double th_eq_in,
                         const char4 atom_i_type = unknown_atom_type,
                         const char4 atom_j_type = unknown_atom_type,
                         const char4 atom_k_type = unknown_atom_type);

  /// \brief Add new dihedral parameters, accessible either by an index determined by the order of
  ///        addition or non-bonded atom types of the participating atoms.
  ///
  /// \param ampl_in      Force constant (amplitude) for the cosine-based torsion term, in units of
  ///                     kcal/mol
  /// \param period_in    Periodicity of the cosine argument (unitless)
  /// \param phi_in       Phase angle for the cosine term, in unis of radians
  /// \param atom_i_type  Non-bonded atom type of the first atom in the torsion.  The atom types
  ///                     are again reflexive, a torsion spanning types tI, tJ, tK, and tL being
  ///                     equivalent to an angle between types tL, tK, tJ, and tI.
  /// \param atom_j_type  Non-bonded atom type of the second atom in the torsion
  /// \param atom_k_type  Non-bonded atom type of the third atom in the torsion
  /// \param atom_l_type  Non-bonded atom type of the fourth atom in the torsion
  void addDiheParameters(double ampl_in, double period_in, double phi_in,
                         const char4 atom_i_type = unknown_atom_type,
                         const char4 atom_j_type = unknown_atom_type,
                         const char4 atom_k_type = unknown_atom_type,
                         const char4 atom_l_type = unknown_atom_type);

  /// \brief Add atoms to the system.  This will create space for new atoms at the end of each
  ///        array or at the position specified, moving all existing atoms back and promoting their
  ///        topological indices as needed.  Residue and molecule limits will be adjust such that
  ///        the added atoms become part of the residue containing the placement point (including
  ///        placement points at the tail of the residue).  If the residue limits are [0, 5),
  ///        [5, 14), and [14, 18) and six atoms are added at placement 5, the first residue gains
  ///        the atoms and the rsidue limits become [0, 11), [11, 20), and [20, 24).  An array of
  ///        indices indicating the positions at which the atoms were added will be returned.
  ///
  /// \param z_numbers                    Atomic numbers of the atoms to add
  /// \param placement                    Topological index to give the first added atom.  Negative
  ///                                     values will add atoms to the back of the list.  Other
  ///                                     added atoms will be inserted in sequence.
  /// \param added_residue_limits         Optional definition of the boundaries of residues within
  ///                                     the added atom set.  Indices apply to the list of added
  ///                                     atoms, not the topology itself, although they will affect
  ///                                     how the residue bounds of the topology develop after the
  ///                                     additions.
  /// \param added_atom_names             Optional names of the added atoms
  /// \param added_atom_types             Optional charater codes for added atom non-bonded types
  /// \param added_residue_names          Optional names of the residues composed by the added
  ///                                     atoms.  Must meet the dimension of the added residue
  ///                                     limits array.
  /// \param added_charge_indices         Optional charge parameter indices for added atoms
  /// \param added_atomic_charges         Optional partial charges for added atoms
  /// \param added_lennard_jones_indices  Optional Lennard-Jones parameter indices to impart to the
  ///                                     added atoms
  std::vector<int> addAtoms(const std::vector<int> &z_numbers, int placement = -1,
                            const std::vector<int> &added_residue_limits = {},
                            const std::vector<char4> &added_atom_names = {},
                            const std::vector<char4> &added_atom_types = {},
                            const std::vector<char4> &added_residue_names = {},
                            const std::vector<int> &added_charge_indices = {},
                            const std::vector<int> &added_lennard_jones_indices = {},
                            const std::vector<double> &added_atomic_charges = {});

  /// \brief Add an atom to the topology.  This is a single-atom variant of the addAtoms() function
  ///        above, which it calls internally, and returns a single index where the atom was
  ///        placed.
  ///
  /// \param z_number             Atomic number of the atom to add, required to specify the element
  /// \param placement            Topological index at which to insert the atom into the topology
  ///                             as it exists when the function is called.  If negative or beyond
  ///                             the end of the list, the atom will be palced at the end of the
  ///                             list.
  /// \param atom_name            Name of the atom
  /// \param atom_type            Name of the non-bonded atom type
  /// \param res_name             Name of the residue that the atom is to be part of, if the atom
  ///                             is seeding a new residue
  /// \param charge_index         Partial charge parameter index of the atom
  /// \param lennard_jones_index  Lennard-Jones parameter index of the atom
  /// \param atomic_charge        Partial charge of the atom
  int addAtom(int z_number, int placement = -1, const char4 atom_name = { 'U', 'N', 'K', ' ' },
              const char4 atom_type = { 'U', 'N', 'K', ' ' },
              const char4 res_name = { 'U', 'N', 'K', ' ' }, int charge_index = -1,
              int lennard_jones_index = -1, double atomic_charge = 0.0);
  
  /// \brief Add a series of virtual sites to a topology.  A frame parameter set of each specified
  ///        index must be added in advance, and any frame atoms referenced by the sites must also
  ///        exist in advance.  If other atoms or other virtual sites are added later, the atom
  ///        indexing of existing virtual sites will update along with the growing topology.  The
  ///        indices of the added virtual sites in the topology are returned so that future
  ///        additions, whether of virtual sites or other atoms, can be adjusted accordingly.
  ///
  /// \param parameter_indices  Indices of the virtual site parameter set to apply (this must
  ///                           match the number of frame atoms provided)
  /// \param frame1_atoms       The parent atoms of each virtual site
  /// \param frame2_atoms       Second atoms in each virtual site frame (all site types have at
  ///                           least two atoms)
  /// \param frame3_atoms       Optional third atoms in the virtual site frames
  /// \param frame4_atoms       Optional fourth atoms in the virtual site frames
  std::vector<int> addVirtualSites(const std::vector<int> &parameter_indices,
                                   const std::vector<int> &frame1_atoms,
                                   const std::vector<int> &frame2_atoms,
                                   const std::vector<int> &frame3_atoms,
                                   const std::vector<int> &frame4_atoms);
  
  /// \brief Add a virtual site to the topology.  This is a single-particle variant of the
  ///        addVirtualSites() function above, which it calls internally, and returns a single
  ///        index where the virtual site was placed.
  ///
  /// \param parameter_index  Index of the virtual site parameter set to apply (this must match the
  ///                         number of frame atoms provided)
  /// \param frame1_atom      The parent atom of the virtual site
  /// \param frame2_atom      Second atom in the frame (all site types have at least two atoms)
  /// \param frame3_atom      Optional third atom in the frame
  /// \param frame4_atom      Optional fourth atom in the frame
  int addVirtualSite(int parameter_index, int frame1_atom, int frame2_atom, int frame3_atom = -1,
                     int frame4_atom = -1);

  /// \brief Set multiple bonds throughout the structure.  The bond term lists, bond parameter
  ///        lists, and number of bonds will be updated as appropriate.
  ///
  /// \param atom_i           The first atom in each of the bonds
  /// \param atom_j           The second atom in each of the bonds
  /// \param parameter_index  Indices of the parameter sets to use for each bond
  void setBonds(const std::vector<int> &atom_i, const std::vector<int> &atom_j,
                const std::vector<int> &parameter_index);

  /// \brief Set a bond between two atoms.  The bond term lists, bond parameter lists, and number
  ///        of bonds will be updated as appropriate.
  ///
  /// \param atom_i           The first atom in the bond
  /// \param atom_j           The second atom in the bond
  /// \param parameter_index  Index of the parameter set to use for the bond
  void setBond(int atom_i, int atom_j, int parameter_index);

  /// \brief Set the general screening factor on 1:4 electrostatic interactions.
  ///
  /// \param screening_factor_in  The screening factor to apply
  void setElectrostatic14Screening(double screening_factor_in);

  /// \brief Set the general screening factor on 1:4 van-der Waals interactions.
  ///
  /// \param screening_factor_in  The screening factor to apply
  void setVanDerWaals14Screening(double screening_factor_in);
  
private:

  // The title and source can be specified.
  ExceptionResponse policy;  ///< Indicate the way to handle bad input or manipulations
  std::string title;         ///< Title of the topology
  std::string output_file;   ///< Name of the file best associated with this topology.  If
                             ///<   constructed based on an existing AtomGraph, the source file of
                             ///<   that object will be taken.  This field can be modified
                             ///<   post-construction and set to hold the name of a file to write.

  /// A list of identifiers for any number of force fields that the topology file references
  std::vector<Citation> force_fields;
  
  // Sizing constants follow the AtomGraph class
  int atom_count;                 ///< Total number of atoms and virtual sites
  int residue_count;              ///< Total number of residues, including solvent molecules
  int molecule_count;             ///< Total number of molecules in the system

  // Arrays of limits define the structures and composition of the molecular system.
  std::vector<int> residue_limits;       ///< Bounds array for each residue in the system
  std::vector<int> molecule_limits;      ///< Bounds array for each molecule in the system,
                                         ///<   indexing molecule_contents (which will likely track
                                         ///<   topological indices, although this is not required)
  std::vector<int> molecule_membership;  ///< Indices of the molecules to which each atom belongs,
                                         ///<   the inverse of the information in molecule_contents
  std::vector<int> molecule_contents;    ///< Topological indices of each atom in any given
                                         ///<   molecule, bounded by molecule_limits

  // Arrays of atomic descriptors define the chemical properties of the system
  std::vector<int> atomic_numbers;     ///< Atomic numbers identifying the elements of each atom

  // Arrays of atom and residue names
  std::vector<char4> atom_names;     ///< Names of individual atoms, atom_count in length
  std::vector<char4> atom_types;     ///< Names of non-bonded atom property types (these atom types
                                     ///<   may also guide the application of non-bonded
                                     ///<   parameters)
  std::vector<char4> residue_names;  ///< Names of individual residues, residue_count in length

  // Sizing constants on CHARMM force field valence terms
  int urey_bradley_term_count;                ///< Total number of Urey-Bradley angle stretch terms
  int charmm_impr_term_count;                 ///< Total number of CHARMM impropers
  int cmap_term_count;                        ///< Total number of CMAP terms
  int urey_bradley_parameter_count;           ///< Number of unique Urey-Bradley parameter pairs
  int charmm_impr_parameter_count;            ///< Number of unique CHARMM improper parameter pairs
  int cmap_surface_count;                     ///< Number of unique CMAP surfaces

  // Atom indexing on CHARMM force field terms
  std::vector<int> urey_bradley_i_atoms;           ///< First atom in each Urey-Bradley term  
  std::vector<int> urey_bradley_k_atoms;           ///< Second atom in each Urey-Bradley term 
  std::vector<int> urey_bradley_parameter_indices; ///< Index of Urey-Bradley parameters for each
                                                   ///<   term   
  std::vector<int> charmm_impr_i_atoms;            ///< First atom in each CHARMM improper term
  std::vector<int> charmm_impr_j_atoms;            ///< Second atom in each CHARMM improper term
  std::vector<int> charmm_impr_k_atoms;            ///< Third atom in each CHARMM improper term
  std::vector<int> charmm_impr_l_atoms;            ///< Fourth atom in each CHARMM improper term
  std::vector<int> charmm_impr_parameter_indices;  ///< Index of each CHARMM improper's parameter
                                                   ///<   pair
  std::vector<int> cmap_i_atoms;                   ///< 1st:UNK atom of each CMAP 2-D surface term
  std::vector<int> cmap_j_atoms;                   ///< 2nd:1st atom of each CMAP 2-D surface term
  std::vector<int> cmap_k_atoms;                   ///< 3rd:2nd atom of each CMAP 2-D surface term
  std::vector<int> cmap_l_atoms;                   ///< 4th:3rd atom of each CMAP 2-D surface term
  std::vector<int> cmap_m_atoms;                   ///< UNK:4th atom of each CMAP 2-D surface term
  std::vector<int> cmap_surface_indices;           ///< Index of each CMAP term's parameter surface

  // Valence parameters for CHARMM force field terms
  std::vector<int> cmap_surface_dimensions;      ///< Sizes of each CMAP surface (this array has a
                                                 ///<   size based on the number of unique CMAP
                                                 ///<   terms)  
  std::vector<int> cmap_surface_bounds;          ///< Starting points for each CMAP surface (this
                                                 ///<   array has a size based on the number of
                                                 ///<   unique CMAP terms)
  std::vector<double> urey_bradley_stiffnesses;  ///< Stiffness constant of each Urey-Bradley
                                                 ///<   stretching term
  std::vector<double> urey_bradley_equilibria;   ///< Equilibrium length of each Urey-Bradley
                                                 ///<   stretching term
  std::vector<double> charmm_impr_stiffnesses;   ///< CHARMM impropers are harmonic, too!
  std::vector<double> charmm_impr_phase_angles;  ///< The "equilibria" for CHARMM impropers
  std::vector<double> cmap_surfaces;             ///< Concatenated, column-major format matrices
                                                 ///<   for every CMAP surface term
  std::vector<char4> urey_bradley_i_atom_types;  ///< Atom types for the first atom in each
                                                 ///<   Urey-Bradley parameter set.  These type
                                                 ///<   names may be left as "unknown" ('UNK ') if
                                                 ///<   the bonds are not assigned by non-bonded
                                                 ///<   atom types.  If bond terms are assigned by
                                                 ///<   non-bonded atom types, these types will be
                                                 ///<   searched in order to assign parameters and
                                                 ///<   parameter indices, irrespective of any
                                                 ///<   other parameter index inputs.
  std::vector<char4> urey_bradley_k_atom_types;  ///< Atom types for the second atom in each
                                                 ///<   Urey-Bradley parameter set
  std::vector<char4> charmm_impr_i_atom_types;   ///< Atom types for the first atom in each CHARMM
                                                 ///<   improper dihedral parameter set
  std::vector<char4> charmm_impr_j_atom_types;   ///< Atom types for the second atom in each CHARMM
                                                 ///<   improper dihedral parameter set
  std::vector<char4> charmm_impr_k_atom_types;   ///< Atom types for the third atom in each CHARMM
                                                 ///<   improper dihedral parameter set
  std::vector<char4> charmm_impr_l_atom_types;   ///< Atom types for the fourth atom in each CHARMM
                                                 ///<   improper dihedral parameter set
  std::vector<char4> cmap_i_atom_types;          ///< Atom types for the first atom in each CMAP
                                                 ///<   parameter set
  std::vector<char4> cmap_j_atom_types;          ///< Atom types for the second atom in each CMAP
                                                 ///<   parameter set
  std::vector<char4> cmap_k_atom_types;          ///< Atom types for the third atom in each CMAP
                                                 ///<   parameter set
  std::vector<char4> cmap_l_atom_types;          ///< Atom types for the fourth atom in each CMAP
                                                 ///<   parameter set
  std::vector<char4> cmap_m_atom_types;          ///< Atom types for the fifth atom in each CMAP
                                                 ///<   parameter set

  // Sizing constants on common valence terms (AMBER, GROMACS, CHARMM)
  int bond_term_count;           ///< Total number of bonded interactions
  int angl_term_count;           ///< Total number of bond angle interactions
  int dihe_term_count;           ///< Total number of dihedral cosine terms in the system
  int bond_parameter_count;      ///< The number of unique bond parameter sets
  int angl_parameter_count;      ///< The number of unique angle parameter sets
  int dihe_parameter_count;      ///< Number of unique dihedral Fourier basis functions
  int attenuated_14_type_count;  ///< The number of distinct 1:4 scaling factor pairs
  int inferred_14_attenuations;  ///< Most 1:4 exclusions are actually attenuations according to
                                 ///<   parameters found in dihe_elec_screenings and
                                 ///<   dihe_vdw_screenings.  This is the number of 1:4 exclusions
                                 ///<   whose attenuated strengths were not taken directly from
                                 ///<   dihedral interactions spanning the two atoms.  In these
                                 ///<   cases the 1:4 exclusion involves a virtual site which is
                                 ///<   not, itself, subject to dihedral interactions but inherits
                                 ///<   its exclusions from a parent atom which does paticipate in
                                 ///<   dihedral interactions.

  // Atom indexing on common valence terms
  std::vector<int> bond_i_atoms;              ///< The first atom of each bond stretching term
  std::vector<int> bond_j_atoms;              ///< The second atom of each bond stretching term
  std::vector<int> bond_parameter_indices;    ///< Parameter index of each bond stretching term
  std::vector<int> angl_i_atoms;              ///< The first atom of each angle bending term
  std::vector<int> angl_j_atoms;              ///< The second atom of each angle bending term
  std::vector<int> angl_k_atoms;              ///< The third atom of each angle bending term
  std::vector<int> angl_parameter_indices;    ///< Parameter index of each angle bending term
  std::vector<int> dihe_i_atoms;              ///< First atom of each dihedral / torsion term   
  std::vector<int> dihe_j_atoms;              ///< Second atom of each dihedral / torsion term  
  std::vector<int> dihe_k_atoms;              ///< Third atom of each dihedral / torsion term   
  std::vector<int> dihe_l_atoms;              ///< Fourth atom of each dihedral / torsion term  
  std::vector<int> dihe_parameter_indices;    ///< Parameter index of each dihedral twisting term
  std::vector<int> dihe14_parameter_indices;  ///< Parameter index for each dihedral's 1:4
                                              ///<   screening factor pair
  std::vector<int> infr14_i_atoms;            ///< First atoms in 1:4 attenuated pair interactions
                                              ///<   with inferred screening factors
  std::vector<int> infr14_l_atoms;            ///< Second atoms in 1:4 attenuated pair interactions
                                              ///<   with inferred screening factors
  std::vector<int> infr14_parameter_indices;  ///< Parameter indices of dihedral interactions
                                              ///<   guiding 1:4 exclusions with inferred
                                              ///<   attenuations (this index goes into a pair of
                                              ///<   arrays of up to 32 unique pairs of 1:4 scaling
                                              ///<   factors, the same as indexed by
                                              ///<   dihe14_parameter_indices).
  std::vector<char4> bond_modifiers;          ///< Enumerations for special aspects of bonds
  std::vector<char4> angl_modifiers;          ///< Enumerations for special aspects of angles    
  std::vector<char4> dihe_modifiers;          ///< Enumerations for special aspects of dihedrals

  // Valence parameters for common force field terms
  std::vector<double> bond_stiffnesses;    ///< Stiffness of each bond stretch, kcal/mol-A^2
  std::vector<double> bond_equilibria;     ///< Equilibrium lengths of all bonds, A
  std::vector<double> angl_stiffnesses;    ///< Stiffness of each angle bend, kcal/mol-rad^2
  std::vector<double> angl_equilibria;     ///< Equilibrium angle for all bending terms, radians
  std::vector<double> dihe_amplitudes;     ///< Amplitudes of each dihedral cosine term, kcal/mol
  std::vector<double> dihe_periodicities;  ///< Periodicity of each dihedral / torsion cosine term
  std::vector<double> dihe_phase_angles;   ///< Phase angle of each dihedral / torsion cosine term
  std::vector<char4> bond_i_atom_types;    ///< Atom types of the first atom in each bond parameter
                                           ///<   set.  These type names may be left as "unknown"
                                           ///<   ('UNK ') if the bonds are not assigned by
                                           ///<   non-bonded atom types.  If bond terms are
                                           ///<   assigned by non-bonded atom types, these types
                                           ///<   will be searched in order to assign parameters
                                           ///<   and parameter indices, irrespective of any other
                                           ///<   parameter index inputs.
  std::vector<char4> bond_j_atom_types;    ///< Atom types of the second atom in each bond
                                           ///<   parameter set
  std::vector<char4> angl_i_atom_types;    ///< Atom types of the first atom in each angle
                                           ///<   parameter set
  std::vector<char4> angl_j_atom_types;    ///< Atom types of the second atom in each angle
                                           ///<   parameter set
  std::vector<char4> angl_k_atom_types;    ///< Atom types of the third atom in each angle
                                           ///<   parameter set
  std::vector<char4> dihe_i_atom_types;    ///< Atom types of the first atom in each dihedral
                                           ///<   parameter set
  std::vector<char4> dihe_j_atom_types;    ///< Atom types of the second atom in each dihedral
                                           ///<   parameter set
  std::vector<char4> dihe_k_atom_types;    ///< Atom types of the third atom in each dihedral
                                           ///<   parameter set
  std::vector<char4> dihe_l_atom_types;    ///< Atom types of the fourth atom in each dihedral
                                           ///<   parameter set

  // Information relevant to virtual site placement
  int virtual_site_count;                           ///< Number of v-sites / extra points
  int virtual_site_parameter_set_count;             ///< Number of parameter sets describing unique
                                                    ///<   virtual site frames (the length of
                                                    ///<   virtual_site_frame_dim1 and related
                                                    ///<   arrays)
  std::vector<int> virtual_site_atoms;              ///< List of atoms which are massless virtual
                                                    ///<   sites
  std::vector<VirtualSiteKind> vs_frame_types;      ///< Frame types for each virtual site
                                                    ///<   parameter set
  std::vector<int> virtual_site_frame1_atoms;       ///< Parent atoms (frame I) for each virtual
                                                    ///<   site
  std::vector<int> virtual_site_frame2_atoms;       ///< Frame atom 2 for each virtual site
  std::vector<int> virtual_site_frame3_atoms;       ///< Frame atom 3 (optional) for each virtual
                                                    ///<   site   
  std::vector<int> virtual_site_frame4_atoms;       ///< Frame atom 4 (optional) for each virtual
                                                    ///<   site
  std::vector<int> virtual_site_parameter_indices;  ///< Parameter indices for virtual sites,
                                                    ///<   indicating the frame type and the values
                                                    ///<   of up to three dimensional measurements.
                                                    ///<   Virtual sites operate much like valence
                                                    ///<   parameters: each has a unique set of two
                                                    ///<   to four frame atoms and references a set
                                                    ///<   of parameters.           
  std::vector<double> virtual_site_frame_dim1;      ///< First frame dimension for each v-site
  std::vector<double> virtual_site_frame_dim2;      ///< Second (optional) frame dimension for each
                                                    ///<   v-site
  std::vector<double> virtual_site_frame_dim3;      ///< Third (optional) frame dimension for each
                                                    ///<   v-site

  // Sizing constants and settings for the non-bonded calculation
  int charge_type_count;            ///< Number of distinct atomic partial charges
  int lj_type_count;                ///< Number of distinct Lennard-Jones types
  int total_exclusions;             ///< Total number of non-bonded exclusions, including 1:4
  UnitCellType periodic_box_class;  ///< Type of periodic box: 0 = nonperiodic, 1 = rectilinear
                                    ///<   prism, 2 = triclinic prism
  double coulomb_constant;          ///< Coulomb's constant in units of kcal-A/mol-e^2 (Amber
                                    ///<   differs from other programs in terms of what this is,
                                    ///<   so it can be set here)                
  double elec14_screening_factor;   ///< General screening factor for electrostatic 1:4
                                    ///<   interactions, applied if no pair-specific factor is
                                    ///<   available
  double vdw14_screening_factor;    ///< General screening factor for 1:4 van-der Waals
                                    ///<   interactions
  
  // Arrays of non-bonded parameters
  std::vector<int> charge_indices;          ///< Atomic charge indices, 0 to charge_type_count - 1
  std::vector<int> lennard_jones_indices;   ///< Lennard-Jones indices, 0 to atom_type_count - 1
  std::vector<int> nb11_exclusion_bounds;   ///< Bounds of non-bonded 1:1 exclusions to each atom
  std::vector<int> nb11_exclusion_list;     ///< Non-bonded 1:1 exclusions (virtual site::parent
                                            ///<   atom)
  std::vector<int> nb12_exclusion_bounds;   ///< Bounds of non-bonded 1:2 exclusions to each atom
  std::vector<int> nb12_exclusion_list;     ///< Non-bonded 1:2 exclusions (separated by one bond)
  std::vector<int> nb13_exclusion_bounds;   ///< Bounds of non-bonded 1:3 exclusions to each atom
  std::vector<int> nb13_exclusion_list;     ///< Non-bonded 1:3 exclusions (separated by two bonds)
  std::vector<int> nb14_exclusion_bounds;   ///< Bounds of non-bonded 1:4 exclusions to each atom
  std::vector<int> nb14_exclusion_list;     ///< Non-bonded 1:4 exclusions (separated by three
                                            ///<   bonds)
  std::vector<double> atomic_charges;       ///< Partial charges assigned to each atom
  std::vector<double> charge_parameters;    ///< Atomic partial charges, condensed to a list of
                                            ///<   unique values 0 ... charge_type_count - 1 and
                                            ///<   indexed by the charge_indices array
  std::vector<double> lj_a_values;          ///< Lennard-Jones A coefficients, U = A/r^12 - B/r^6
  std::vector<double> lj_b_values;          ///< Lennard-Jones B coefficients, U = A/r^12 - B/r^6
  std::vector<double> lj_c_values;          ///< Lennard-Jones C coefficients, U = A/r^12 - B/r^6
  std::vector<double> lj_14_a_values;       ///< Lennard-Jones A coefficients for 1:4 interactions
  std::vector<double> lj_14_b_values;       ///< Lennard-Jones B coefficients for 1:4 interactions
  std::vector<double> lj_14_c_values;       ///< Lennard-Jones C coefficients for 1:4 interactions
  std::vector<double> lj_sigma_values;      ///< Lennard-Jones sigma parameters, supplied directly
                                            ///<   and determinants of the pair coefficients unless
                                            ///<   the Lennard-Jones rule is set to NBFIX
  std::vector<double> lj_epsilon_values;    ///< Lennard-Jones epsilon parameters, supplied
                                            ///<   directly and determinants of the pair
                                            ///<   coefficients unless the Lennard-Jones rule is
                                            ///<   set to NBFIX
  std::vector<double> lj_14_sigma_values;   ///< Lennard-Jones 1:4 sigma parameters, again supplied
                                            ///<   directly
  std::vector<double> lj_14_epsilon_values; ///< Lennard-Jones 1:4 sigma parameters, again supplied
                                            ///<   directly

  // Whereas the AtomGraph contains condensed arrays of the distinct 1:4 non-bonded screening
  // parameter pairs, the staging object must collect instructions and create the 1:4 screening
  // factors for every dihedral in the topology.  This is accomplished by accumulating a series of
  // atom type / screening factor pair records, then building the 1:4 screenings for all dihedrals
  // in a flash just before exporting a topology.  This array covers pairs of electrostatic and
  // van-der Waals screening factors for up to four distinct atom types. 
  std::vector<ScreeningFactorTag> attn14_tags;

  // Atom, atom type, and residue name overflow keys
  std::vector<char4> atom_overflow_names;    ///< Codified names of atoms which were too long
  std::vector<char4> atom_overflow_types;    ///< Codified names of atom types which were too long
  std::vector<char4> residue_overflow_names; ///< Codified names of residues which were too long

  /// \brief Validate the requeste atom index.
  ///
  /// \param atom_index  Index of the atom of interest
  void validateAtomIndex(int atom_index) const;

  /// \brief Validate the index of a virtual site frame's parameters.
  ///
  /// \param set_index  The index of the requested virtual site parameter set
  void validateFrameIndex(int set_index) const;

  /// \brief Insert a particle, atom or virtual site, into the topological order.
  ///
  /// \param z_numbers              The atomic number of the particle (0 indicates it is a virtual
  ///                               site, any positive value a real atom, other values are errors)
  /// \param residue_homes          The residue that is to absorb the new particle (must work in
  ///                               concert with z_numbers--if there are no entries of
  ///                               residue_homes to cover some or all entries of z_numbers, the
  ///                               remaining atoms will be added to the back of the atom list in
  ///                               a unique residue containing all unassigned atoms)
  /// \param added_vs_frame_types   Frame types for any added virtual sites, must be long enough to
  ///                               cover any virtual sites among the added atoms (an atomic number
  ///                               of 0 in the z_numbers array) with NONE as the kind otherwise
  /// \param added_vs_frame1_atoms  Atoms (in the ordering found when the function is first called)
  ///                               that will serve as the parent atoms for each added virtual site
  ///                               (the actual frame1 index recorded will reflect changes made due
  ///                               to the particle insertions).  This array and others must be
  ///                               large enough to cover any virtual sites in the input.
  /// \param added_vs_frame2_atoms  The second frame atom for each added virtual site
  /// \param added_vs_frame3_atoms  The optional third atom for each added virtual site
  /// \param added_vs_frame4_atoms  The optional fourth frame atom for each added virtual site
  std::vector<int> insertParticles(const std::vector<int> &z_numbers,
                                   const std::vector<int> &residue_homes, 
                                   const std::vector<int> &added_vs_parameter_indices = {},
                                   const std::vector<int> &added_vs_frame1_atoms = {},
                                   const std::vector<int> &added_vs_frame2_atoms = {},
                                   const std::vector<int> &added_vs_frame3_atoms = {},
                                   const std::vector<int> &added_vs_frame4_atoms = {});

  /// \brief Add a series of connections between pairs of atoms.  The exclusions list will be
  ///        updated as appropriate.  The number of bond terms, angle terms, and dihedral terms
  ///        will not be updated and must be modified in step using the setBond(), setAngle(), and
  ///        setDihedral() functions.
  std::vector<int> addConnections(const std::vector<int> &atom_i, const std::vector<int> &atom_j);
  
  /// \brief Add a connection between two atoms.  The exclusions lists will be updated as
  ///        appropriate.  The number of bond terms, angle terms, and dihedral terms will not be
  ///        updated and must be modified in step using the setBond(), setAngle(), and
  ///        setDihedral() functions.
  ///
  /// Overloaded:
  ///   - Add a connection between two atoms
  ///   - Add connections between two lists of atoms (list I:0 -> list J:0, list I:1 -> list J:1,
  ///     list I:2 -> list J:2 in an element-wise fashion, not all combinations from each list)
  ///
  /// \param atom_i  The first atom to connect
  /// \param atom_j  The second atom to connect
  int addConnection(int atom_i, int atom_j);

  /// \brief Check that an array of supplied properties matches the number of atoms to be added,
  ///        if the array has any length at all.  Return TRUE if the length of the array is nonzero
  ///        and raise an exception if there is any mismatch otherwise.
  ///
  /// \param array_size   The length of the atom property array
  /// \param added_count  The number of atoms to be added to the system
  /// \param desc         Description of the array provided
  bool applyProperty(int array_size, int added_count, const std::string &desc);

  /// \brief Mask the solute atoms in the topology, returning a boolean mask of qualifying atoms.
  ///
  /// \param minimum_solute_size             The minimum size of a solute molecule.  Anything
  ///                                        larger will be automatically included in the solute.
  /// \param solute_included_residues_names  Names of additional residues that should be included
  ///                                        in the solute
  std::vector<bool>
  maskSoluteAtoms(int minimum_solute_size = default_minimum_solute_atoms,
                  const std::vector<char4> &solute_included_residue_names = {}) const;

  /// \brief Create the array of Amber prmtop format descriptors, a series of about 38 integers
  ///        indicating the numbers of atoms, various energy terms, box shape, and
  ///        simulation-related partitions.  This list appears near the top of an Amber topology
  ///        file.
  ///
  /// \param largest_residue_atoms  The number of atoms in the largest residue (calculated prior
  ///                               to calling this routine and passed in for convenience)
  /// \param total_exclusions       The total number of exclusions that will be listed in the
  ///                               Amber prmtop format, including values of "zero" (no exclusion)
  ///                               for atoms that have no exclusions of higher index atoms. 
  std::vector<int> createPrmtopDescriptors(int largest_residue_atoms,
                                           int total_exclusions) const;

  /// \brief Compute the atomic masses based on the atomic numbers and any hydrogen mass
  ///        repartitioning requested.
  ///
  /// \param hmass_repartition_factor  Hydrogen mass repartitioning factor, 1.0 implies no change
  std::vector<double> computeAtomicMasses(double hmass_repartition_factor) const;

  /// \brief Prepare arrays of 1:4 non-bonded screening factors for the accumulated set of dihedral
  ///        parameters.  The manner in which valence parameters are applied (i.e. by chemical
  ///        perception or lookup based on assigned atom types) immaterial: dihedral parameters
  ///        will span various 1:4 interactions and link different combinations of non-bonded atom
  ///        types, which this function will sort into coherent arrays of scaling factors based on
  ///        the provided general electrostatic or van-der Waals 1:4 scaling factors and a list of
  ///        additional parameters tailored for particular atom types.  If there are no
  ///        type-specific parameters, no dihedral-specific factors will be developed so as to let
  ///        subsequent routines rely on the general screening factors.
  ///
  /// \param attn14_elec_factors  Electrostatic 1:4 screening factors for each dihedral parameter
  ///                             set, modified and returned
  /// \param attn14_vdw_factors   Van-der Waals 1:4 screening factors for each dihedral parameter
  ///                             set, modified and returned
  void buildNonbonded14Screen(std::vector<double> *attn14_elec_factors,
                              std::vector<double> *attn14_vdw_factors) const;
};

/// \brief Create a vector with the topological indices of an atom and any virtual sites to which
///        the atom is a parent.
std::vector<int> atomWithVirtualSiteChildren(const int atom_index,
                                             const std::vector<int> &nb11_exclusion_list,
                                             const std::vector<int> &nb11_exclusion_bounds);

} // namespace topology
} // namespace stormm

#endif
