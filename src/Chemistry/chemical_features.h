// -*-c++-*-
#ifndef STORMM_CHEMICAL_FEATURES_H
#define STORMM_CHEMICAL_FEATURES_H

#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/vector_ops.h"
#include "Numerics/split_fixed_precision.h"
#include "Topology/atomgraph.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"
#include "UnitTesting/stopwatch.h"
#include "chemistry_enumerators.h"
#include "match_bonding_pattern.h"

namespace stormm {
namespace chemistry {

using card::Hybrid;
using card::HybridTargetLevel;
using data_types::isFloatingPointScalarType;
using stmath::crossProduct;
using stmath::dot;
using stmath::project;
using stmath::roundUp;
using synthesis::Condensate;
using synthesis::CondensateReader;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisReader;
using testing::StopWatch;
using topology::AtomGraph;
using topology::ChemicalDetailsKit;
using topology::NonbondedKit;
using topology::ValenceKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesReader;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceReader;

/// \brief The maximum number of nodes permitted in the tree searches of fused ring systems
constexpr int max_fused_ring_tree_size = 33554432;
  
/// \brief A struct to serve as a tracker of progress through a molecule in the search for loops.
///        Proceeding forward in the search, every link will have come from one and only one
///        previous link, but could go in multiple directions thereafter.
class BondedNode {
public:
  
  /// \brief Basic constructor simply initializes the members to blank values
  BondedNode();

  /// \brief Take a pre-existing std::vector<int> in order to set the object's own next_links
  ///        pointer.  This is in the interest of not repeatedly allocating and freeing many tiny
  ///        std::vector<int> objects as BondedNodes are created and destroyed to grow a chain of
  ///        atoms.
   ///
  /// \param vi            Vector of integers in which this BondedNode will store some of its data
  /// \param pos           Indicator of the space in the storage vector allocate to this object
  /// \param max_branches  The maximum number of branches that this chain link shall support
  void setBranchPointer(std::vector<int> *vi, size_t pos, size_t max_branches);

  /// \brief Set this as the next link in the chain, or as the chain initiator.
  ///
  /// Overloaded:
  ///   - Optionally accept the index of the previous node for traceback purposes.
  ///
  /// \param previous_in      Index of the previous BondedNode leading to this one, or -1 if this
  ///                         is intended to be the chain initiator
  /// \param current_atom     Topological index of the current atom this BondedNode describes
  /// \param current_layer    The current layer in the growing tree
  /// \param nbk              Nonbonded interactions abstract for the overarching topology (needed
  ///                         for list of 1:2 exclusions, which indicate bonds)
  /// \param valid_atom_mask  Reference table of which atoms are valid to include in the tree
  void addToTree(int previous_in, int current_atom, int current_layer,
                 const NonbondedKit<double> &nbk, const std::vector<bool> &valid_atom_mask);
  void addToTree(int previous_in, int current_atom, int current_layer, int previous_node_in,
                 const NonbondedKit<double> &nbk, const std::vector<bool> &valid_atom_mask);

  /// \brief Add the bond order between this atom and its previous atom
  ///
  /// \param vk           Valence interactions abstract from the original topology
  /// \param bond_orders  Orders of bonds determined by Lewis structure drawing
  void addBondOrder(const ValenceKit<double> &vk, const Hybrid<double> &bond_orders);
  
  /// \brief Get the previous atom index
  int getPreviousAtom() const;

  /// \brief Get the previous node's list index
  int getPreviousNode() const;

  /// \brief Get the central atom for this BondedNode.
  int getAtom() const;

  /// \brief Get the layer of which this BondedNode is a part.
  int getLayer() const;

  /// \brief Get the number of branches associated with this BondedNode's central atom (this is the
  ///        total number of other atoms that it bonds to, less the previous atom in the chain)
  int getBranchCount() const;

  /// \brief Get the topological atom index of a specific branch that can come off of this
  ///        BondedNode.
  ///
  /// \param index  The branch index of interest (could be any topological atom)
  int getBranchAtom(int index) const;

  /// \brief Get the branch index of a specific atom in a BondedNode.
  ///
  /// \param search_index  The topological atom index to search for within the node's branches
  int findBranchAtom(int search_index) const;

  /// \brief Get the ring completion status for a particular node
  ///
  /// \param branch_index  The branch index of the partner atom with which this node forms a ring
  ///                      (this index is not a topological index, and must be ascertained in
  ///                      advance by findBranchAtom() above)
  uint getRingCompletion(int branch_index) const;

  /// \brief Get the bond order of the bond between this atom and the previous atom in the tree.
  double getRootBondOrder() const;

  /// \brief Set the ring completion for this object (this relies on no atom having more bonds than
  ///        the number of bits in an unsigned int)
  ///
  /// \param partner_index  Topological index of the atom to which a ring has been formed
  void setRingCompletion(int partner_index);

  /// \brief Wipe clean the ring completion for this object
  void wipeRingCompletion();

private:
  int previous_atom_index;   ///< Topological index of the previous atom
  int previous_node_index;   ///< Tree list index of the previous node
  int atom_index;            ///< Topological index of this atom
  int layer_index;           ///< Index of the layer to which this BondedNode belongs.  The first
                             ///<   atom in the tree is layer zero.
  double root_bond_order;    ///< Order by which this atom bonds to the previous atom (initalized
                             ///<   to zero by default, filled with addBondOrder())
  int branch_count;          ///< The number of branching particles coming off of this atom,
                             ///<   excluding the previous atom
  int* branch_atoms;         ///< Pointer to data array containing indices of all branch atoms
  uint rings_completed;      ///< Bitmask indicating, in the nth bit, whether a ring has been
                             ///<   completed by following the nth branch
};

/// \brief A shorthand form of a conformational degree of freedom, which could be a rotable bond,
///        cis-trans invertible bond, or chiral center.
class IsomerPlan {
public:

  /// \brief The constructor takes all applicable arguments.  Different overloads check the value
  ///        of the ConformationEdit parameter to determine whether they have all of the necessary
  ///        information.
  ///
  /// Overloaded:
  ///   - Accept a list of moving atoms
  ///   - Leave the field blank to be filled later
  ///   - Specify the chiral inversion protocol (used only in the case of motion_in being
  ///     CHIRAL_INVERSION)
  /// \{
  IsomerPlan(ConformationEdit motion_in, ChiralInversionProtocol chiral_plan_in, int root_atom_in,
             int pivot_atom_in, const std::vector<int> &moving_atoms_in,
             const AtomGraph* ag_pointer_in);

  IsomerPlan(ConformationEdit motion_in, int root_atom_in, int pivot_atom_in,
             const std::vector<int> &moving_atoms_in, const AtomGraph* ag_pointer_in);

  IsomerPlan(ConformationEdit motion_in, ChiralInversionProtocol chiral_plan_in, int root_atom_in,
             int pivot_atom_in, const AtomGraph* ag_pointer_in);

  IsomerPlan(ConformationEdit motion_in, int root_atom_in, int pivot_atom_in,
             const AtomGraph* ag_pointer_in);
  /// \}

  /// \brief The default copy and move constructors will suffice, and with no const members the
  ///        copy and move assignment operators can also take on their default forms.
  /// \{
  IsomerPlan(const IsomerPlan &original) = default;
  IsomerPlan(IsomerPlan &&original) = default;
  IsomerPlan& operator=(const IsomerPlan &other) = default;
  IsomerPlan& operator=(IsomerPlan &&other) = default;
  /// \}

  /// \brief Get the isomerizing motion.
  ConformationEdit getMotion() const;

  /// \brief Get the plan for inverting a chiral center.
  ChiralInversionProtocol getChiralPlan() const;

  /// \brief Get the root atom for the motion.
  int getRootAtom() const;

  /// \brief Get the pivot atom for the motion.
  int getPivotAtom() const;

  /// \brief Get the root atom handle (for a rotatable bond J-K, this might be I of the sequence
  ///        I-J-K-L).
  int getRootHandle() const;

  /// \brief Get the pivot atom handle (for a rotatable bond J-K, this might be L of the sequence
  ///        I-J-K-L).
  int getPivotHandle() const;

  /// \brief Get the number of moving atoms.
  int getMovingAtomCount() const;

  /// \brief Get the index of an individual moving atom from within this plan.
  ///
  /// \param index  Index of the atom of interest within the plan (the topological index of the
  ///               atom will be returned)
  int getMovingAtom(const size_t index) const;
  
  /// \brief Get a const reference to the stored list of atoms that move as a consequence of the
  ///        isomerization.
  const std::vector<int>& getMovingAtoms() const;

  /// \brief Add atoms to the moving atoms list.
  ///
  /// \param new_atom_idx  The list of new atoms
  void addMovingAtoms(const std::vector<int> &new_atom_idx);

  /// \brief Erase the list of moving atoms to start over.  This will resize the vector to zero.
  void eraseMovingAtoms();
  
private:
  ConformationEdit motion;              ///< The motion to make that will change the conformation
                                        ///<   (this includes bond rotations as well as true
                                        ///<   isomerizations)
  ChiralInversionProtocol chiral_plan;  ///< Plan for inverting a chiral center, if that is the
                                        ///<   motion involved in this isomer plan
  int root_atom;                        ///< Root of the rotatable bond, chiral center atom
  int pivot_atom;                       ///< The end atom of the rotatable bond (unused if motion
                                        ///<   is set to CHIRAL_INVERSION)
  int root_handle;                      ///< The most significant atom connected to the root atom,
                                        ///<   decided by the highest Z number or, in the case of a
                                        ///<   tie, the atom with the lowest topological index.
                                        ///<   Only relevant in the case of a rotatable or
                                        ///<   cis-trans bond, set to -1 for chiral centers.
  int pivot_handle;                     ///< The most significan atom connected to the pivot atom
  std::vector<int> moving_atoms;        ///< List of all atoms that turn as a consequence of
                                        ///<   twisting about the rotatable bond axis
  AtomGraph *ag_pointer;                ///< Pointer to the AtomGraph object to which this
                                        ///<   isomerization applies
};

/// \brief A pairing of some process (a ConformationEdit, being a rotation about a bond, flipping
///        groups across an isomeric cis-trans bond, or inversion of a chiral center) with an
///        integer index to indicate its association with some list or group.  This struct serves
///        the SynthesisPermutor (see synthesis_permutor.h and the synthesis namespace) as well as
///        the IsomerPlan object (see above).  In Hybrid objects, it is transmutable with int2.
struct CoupledEdit {
public:

  /// An explicit statement of the default constructor, plus additional constructors based on the
  /// related int2 type.
  /// \{
  CoupledEdit();
  CoupledEdit(ConformationEdit edit_in, int index_in);
  CoupledEdit(int2 data_in);
  /// \}
  
  /// \brief The copy and move constructors, as well as assignment operators, can take their
  ///        default forms.
  /// \{
  CoupledEdit(const CoupledEdit &original) = default;
  CoupledEdit(CoupledEdit &&original) = default;
  CoupledEdit& operator=(const CoupledEdit &other) = default;
  CoupledEdit& operator=(CoupledEdit &&other) = default;
  /// \}

  /// \brief Overload the assignment operator to permit communication back and forth with int2.
  CoupledEdit operator=(const int2 &other);
  
  ConformationEdit edit;  ///< The process (this will put the index in context: does it apply to
      	      	      	  ///<   a connected rotatable bond, cis-trans bond, or a chiral center?
  int index;              ///< Index of the connected feature from whatever appropriate list
};

/// \brief Abstract of the chemical features object: this is about as complex as an object with a
///        single abstract should get.
struct ChemicalFeaturesReader {

  /// \brief The constructor takes a straight list of critical constants and pointers.
  ChemicalFeaturesReader(int atom_count_in, int planar_atom_count_in, int ring_count_in,
                         int fused_ring_count_in, int twistable_ring_count_in,
                         int conjugated_group_count_in, int aromatic_group_count_in,
                         int polar_hydrogen_count_in, int hbond_donor_count_in,
                         int hbond_acceptor_count_in, int chiral_center_count_in,
                         int rotatable_bond_count_in, int cis_trans_bond_count_in,
                         int double_bond_count_in, int triple_bond_count_in,
                         int max_ring_size_in, double temperature_in,
                         bool rotating_groups_mapped_in, bool chiralities_computed_in,
                         const int* planar_centers_in, const int* ring_atoms_in,
                         const int* ring_atom_bounds_in, const int* aromatic_pi_electrons_in,
                         const int* aromatic_groups_in, const int* aromatic_group_bounds_in,
                         const int* polar_hydrogens_in, const int* hydrogen_bond_donors_in,
                         const int* hydrogen_bond_acceptors_in, const int* chiral_centers_in,
                         const int* chiral_inversion_methods_in, const int* rotatable_groups_in,
                         const int* rotatable_group_bounds_in, const int* cis_trans_groups_in,
                         const int* cis_trans_group_bounds_in, const int* invertible_groups_in,
                         const int* invertible_group_bounds_in, const int4* chiral_arm_atoms_in,
                         const double* formal_charges_in, const double* bond_orders_in,
                         const double* free_electrons_in, const ullint* ring_inclusion_in,
                         const AtomGraph *ag_pointer_in);

  /// \brief The default copy and move constructors are acceptable, but assignment operators are
  ///        implicitly deleted due to the presence of const elements.
  ///
  /// \param original  The original object to copy or move
  /// \{
  ChemicalFeaturesReader(const ChemicalFeaturesReader &original) = default;
  ChemicalFeaturesReader(ChemicalFeaturesReader &&original) = default;
  /// \}

  // Counts of atoms and features (see the parent object for more detailed information)
  const int atom_count;               ///< Total number of atoms in the system
  const int planar_atom_count;        ///< Number of atoms with enforced planarity
  const int ring_count;               ///< Number of rings
  const int fused_ring_count;         ///< Number of fused rings
  const int twistable_ring_count;     ///< Number of twistable (i.e. boat-chair) rings
  const int conjugated_group_count;   ///< Number of separate conjugated pi-systems
  const int aromatic_group_count;     ///< Number of distinct aromatic groups
  const int polar_hydrogen_count;     ///< Number of polar hydrogens
  const int hbond_donor_count;        ///< Number of hydrogen bond donors
  const int hbond_acceptor_count;     ///< Number of hydrogen bond acceptors
  const int chiral_center_count;      ///< Number of chiral centers
  const int rotatable_bond_count;     ///< Number of rotatable single bonds
  const int cis_trans_bond_count;     ///< Number of cis-trans flippable bonds
  const int double_bond_count;        ///< Number of double bonds
  const int triple_bond_count;        ///< Number of triple bonds
  const int max_ring_size;            ///< The maximum ring size found
  const double temperature;           ///< Temperature at which resonance structures were computed
  const bool rotating_groups_mapped;  ///< Indicate whether the rotatable groups of the system
                                      ///<   (as well as chiral inversion groups) have been mapped
  const bool chiralities_computed;    ///< Indicate whether the chiral centers of the system have
                                      ///<   been determined

  // Lists of atom indices and atom properties (see the parent object for detailed descriptions)
  const int* planar_centers;            ///< Topological indices of atoms with enforced planarity
  const int* ring_atoms;                ///< Lists of topological indices for atoms involved in
                                        ///<   rings, arranged back-to-back
  const int* ring_atom_bounds;          ///< Bounds array for ring_atoms
  const int* aromatic_pi_electrons;     ///< Numbers of pi electrons in each aromatic group (not
                                        ///<   a per-atom quantity, and not associated with limits
                                        ///<   in the aromatic_group_bounds array)
  const int* aromatic_groups;           ///< Lists of atoms involved in distinct aromatic groups,
                                        ///<   arranged back-to-back
  const int* aromatic_group_bounds;     ///< Bounds array for aromatic_groups
  const int* polar_hydrogens;           ///< List of topological indices for polar hydrogen atoms
  const int* hydrogen_bond_donors;      ///< List of topological indices for H-bond donors
                                        ///<   (atoms bearin gpolar hydrogens)
  const int* hydrogen_bond_acceptors;   ///< List of topological indices for H-bond acceptors
  const int* chiral_centers;            ///< List of topological indices for chiral centers
  const int* chiral_inversion_methods;  ///< Integer translations of the ChiralInversionProtocol
                                        ///<   data types instructing 
  const int* rotatable_groups;          ///< Lists of atoms that move when rotating about single
                                        ///<   bonds, arranged back-to-back
  const int* rotatable_group_bounds;    ///< Bounds array for rotatable groups
  const int* cis_trans_groups;          ///< Lists of atoms that move when creating cis-trans
                                        ///<  isomers, arranged back-to-back
  const int* cis_trans_group_bounds;    ///< Bounds array for cis_trans_groups
  const int* invertible_groups;         ///< Lists of atoms that move based on inversion of a
                                        ///<   particular chiral center (various groups are
                                        ///<   arranged back-to-back in this array)
  const int* invertible_group_bounds;   ///< Bounds array for invertible_groups
  const int4* chiral_arm_atoms;         ///< Topological indices of the base atoms for all four
                                        ///<   arms around each chiral center (the y and z members
                                        ///<   are identical to entries in the anchor_a_branches
                                        ///<   and anchor_b_branches arrays)
  const double* formal_charges;         ///< Formal charges of each atom (real numbers taken from
                                        ///<   the resonance structure, often between 0.0 and 1.0)
  const double* bond_orders;            ///< The order of each bond (real numbers which frequently
                                        ///<   fall between 1.0 and 2.0)
  const double* free_electrons;         ///< Free electron content of each atom
  const ullint* ring_inclusion;         ///< Bitmask indicating the various sizes of rings that
                                        ///<  each atom is found in
  const AtomGraph *ag_pointer;          ///< Pointer to the associated topology (only useable in
                                        ///<   the HOST version of the abstract)
};
  
/// \brief An object to store information about chemical motifs: participation in rings, planarity,
///        chirality, aromaticity, conjugation, planarity, and bonds with different rotatability.
class ChemicalFeatures {
public:

  /// \brief The constructor requires a topology and some coordinate set.
  ///
  /// Overloaded:
  ///   - Create a blank object
  ///   - Create a with a topology only (provided by const pointer or by const reference)
  ///   - Create with a topology and one of the basic (single-system) coordinate objects
  ///
  /// \param ag_in           Pointer to the system topology.  This topology will not be modified by
  ///                        submitting it to this constructor, but it is needed as a constant
  ///                        pointer so that the object itself can store a valid pointer to the
  ///                        original topology (passing by const reference would not create a valid
  ///                        pointer).
  /// \param cfr             Coordinates of the system
  /// \param ps              Coordinates of the system (a CoordinateFrameReader will be extracted)
  /// \param map_groups_in   Indicator of whether to map rotatable groups (this is an O(N^2)
  ///                        algorithm in memory as well as computation, as larger structures will
  ///                        have more rotatable bonds and the entirety of the structure must be
  ///                        traced in relation to each of them)
  /// \param temperature_in  Temperature at which to take Boltzmann weights of different resonance
  ///                        states
  /// \{
  ChemicalFeatures(const AtomGraph *ag_in = nullptr,
                   MapRotatableGroups map_groups_in = MapRotatableGroups::NO,
                   double temperature_in = 300.0, StopWatch *timer_in = nullptr);

  ChemicalFeatures(const AtomGraph &ag_in,
                   MapRotatableGroups map_groups_in = MapRotatableGroups::NO,
                   double temperature_in = 300.0, StopWatch *timer_in = nullptr);
  
  ChemicalFeatures(const AtomGraph *ag_in, const CoordinateFrameReader &cfr,
                   MapRotatableGroups map_groups_in = MapRotatableGroups::NO,
                   double temperature_in = 300.0, StopWatch *timer_in = nullptr);

  ChemicalFeatures(const AtomGraph *ag_in, const CoordinateFrame &cf,
                   MapRotatableGroups map_groups_in = MapRotatableGroups::NO,
                   double temperature_in = 300.0, StopWatch *timer_in = nullptr);

  ChemicalFeatures(const AtomGraph *ag_in, const PhaseSpace &ps,
                   MapRotatableGroups map_groups_in = MapRotatableGroups::NO,
                   double temperature_in = 300.0, StopWatch *timer_in = nullptr);

  ChemicalFeatures(const AtomGraph &ag_in, const CoordinateFrame &cf,
                   MapRotatableGroups map_groups_in = MapRotatableGroups::NO,
                   double temperature_in = 300.0, StopWatch *timer_in = nullptr);

  ChemicalFeatures(const AtomGraph &ag_in, const PhaseSpace &ps,
                   MapRotatableGroups map_groups_in = MapRotatableGroups::NO,
                   double temperature_in = 300.0, StopWatch *timer_in = nullptr);
  /// \}

  /// \brief Copy and move constructors
  ///
  /// \param original  The ChemicalFeatures object to copy
  /// \{
  ChemicalFeatures(const ChemicalFeatures &original);
  ChemicalFeatures(ChemicalFeatures &&original);
  /// \}
  
  /// \brief Copy assignment and move assignment operators
  ///
  /// \param other     The ChemicalFeatures object to copy (a different name for a better semantic
  ///                  fit in the context of the = sign)
  /// \{
  ChemicalFeatures& operator=(const ChemicalFeatures &other);
  ChemicalFeatures& operator=(ChemicalFeatures &&other);
  /// \}
  
#ifdef STORMM_USE_HPC
  /// \brief Upload data to the GPU.
  void upload();

  /// \brief Download data from the GPU.
  void download();
#endif

  /// \brief Get the total number of atoms.
  int getAtomCount() const;
  
  /// \brief Get the number of planar atoms.
  int getPlanarAtomCount() const;

  /// \brief Get the number of rings in the system.
  int getRingCount() const;

  /// \brief Get the number of fused rings in the system.
  int getFusedRingCount() const;

  /// \brief Get the number of malleable, non-aromatic twistable rings in the system.
  int getMutableRingCount() const;

  /// \brief Get the number of aromatic groups in the system.
  int getAromaticGroupCount() const;

  /// \brief Get the number of polar hydrogens in the system.
  int getPolarHydrogenCount() const;

  /// \brief Get the number of hydrogen bond donors in the system.
  int getHydrogenBondDonorCount() const;

  /// \brief Get the number of hydrogen bond acceptors in the system.
  int getHydrogenBondAcceptorCount() const;

  /// \brief Get the number of chiral centers in the system.
  int getChiralCenterCount() const;

  /// \brief Get the number of rotatable bonds in the system.
  int getRotatableBondCount() const;

  /// \brief Get the number of bonds involved in cis-trans isomerization.
  int getCisTransBondCount() const;

  /// \brief Indicate whether the rotatable groups (including rotatable bonds as well as cis-trans
  ///        bonds) have been computed.
  bool rotatableGroupsMapped() const;
  
  /// \brief Indicate whether chiralities have been determined for this set of chemical features.
  bool chiralitiesComputed() const;
  
  /// \brief Return a mask of rings within a given size range for this system.
  ///
  /// \param min_ring_size  The minimum number of atoms in the rings that will be reported
  /// \param max_ring_size  The maximum number of atoms in the rings that will be reported
  std::vector<uint> getRingMask(int min_ring_size, int max_ring_size) const;

  /// \brief Return a mask of atomatic atoms in the system.
  ///
  /// \param min_pi_electrons  Minimum number of electrons in the aromatic ring system to report
  /// \param max_pi_electrons  Maximum number of electrons in the aromatic ring system to report
  std::vector<uint> getAromaticMask(int min_pi_electrons, int max_pi_electrons) const;

  /// \brief Get a list of all polar hydrogen atoms.
  std::vector<int> getPolarHydrogenList() const;

  /// \brief Get a list of all hydrogen bond donor atoms.
  std::vector<int> getHydrogenBondDonorList() const;

  /// \brief Get a list of all hydrogen bond acceptor atoms.
  std::vector<int> getHydrogenBondAcceptorList() const;

  /// \brief Get a bit mask of all polar hydrogen atoms in the system, acceptable for inputs to
  ///        creating new atom masks.
  std::vector<uint> getPolarHydrogenMask() const;

  /// \brief Get a bit mask of all hydrogen bond donors in the system, acceptable for inputs to
  ///        creating new atom masks.
  std::vector<uint> getHydrogenBondDonorMask() const;

  /// \brief Get a bit mask of all hydrogen bond acceptors in the system, acceptable for inputs to
  ///        creating new atom masks.
  std::vector<uint> getHydrogenBondAcceptorMask() const;

  /// \brief List the chiral centers in a system, using topological indices.
  ///
  /// \param direction  Preferred chiral orientation of the centers to return (D-, L-, or both)
  std::vector<int> getChiralCenters(ChiralOrientation direction = ChiralOrientation::NONE) const;

  /// \brief Get the bases of each arm for all chiral centers.  The result returns the lowest
  ///        priority arm in the "x" member, the highest in the "y" member, and the second- and
  ///        third-highest priority arms in the "z" and "w" members of each tuple, respectively.
  std::vector<int4> getChiralArmBaseAtoms() const;
  
  /// \brief Return the chiral orientations for one or more atoms in the system.
  ///
  /// Overloaded:
  ///   - Get the chiral orientation for a specific atom
  ///   - Get the chiral orientations for a range of atoms
  ///   - Get the chiral orientations for all atoms
  ///
  /// \param atom_index  The one atom of interest, based on the topological ordering
  /// \param low_index   The start of a series of atoms for which to get ring inclusions
  /// \param high_index  The upper limit of a series of atoms for which to get ring inclusions
  /// \{
  ChiralOrientation getAtomChirality(int atom_index) const;
  std::vector<ChiralOrientation> getAtomChirality(int low_index, int high_index) const;
  std::vector<ChiralOrientation> getAtomChirality() const;  
  /// \}

  /// \brief Return a mask of chiral centers in the system.
  ///
  /// \param direction  Allows one to select R- (D-), S- (L-), or both chiralities for the mask
  std::vector<uint> getChiralityMask(ChiralOrientation direction) const;

  /// \brief Return a vector containing the formal charges on all particles in the system (this
  ///        includes virtual sites, which will have formal charges of zero since they are not real
  ///        atoms).
  std::vector<double> getFormalCharges() const;

  /// \brief Return a vector containing the orders of all bonds in the system (this includes bonds
  ///        to virtual sites only if they are defined in the topology's connectivity, and such
  ///        bonds will have order zero).
  std::vector<double> getBondOrders() const;

  /// \brief Return the (resonance-averaged) free electron content for one or more atoms in the
  ///        system.
  ///
  /// Overloaded:
  ///   - Get the free electron content for a specific atom
  ///   - Get the free electron content for a range of atoms
  ///   - Get the free electron content for all atoms
  ///
  /// \param atom_index  The one atom of interest
  /// \param low_index   The start of a series of atoms for which to get ring inclusions
  /// \param high_index  The upper limit of a series of atoms for which to get ring inclusions
  /// \{
  double getFreeElectrons(int atom_index) const;
  std::vector<double> getFreeElectrons(int low_index, int high_index) const;
  std::vector<double> getFreeElectrons() const;  
  /// \}

  /// \brief Return the formal charges for a representative Lewis structure of the molecules in
  ///        this topology.
  const std::vector<int>& getZeroKelvinFormalCharges() const;

  /// \brief Return the bond orders for a representative Lewis structure of the molecules in this
  ///        topology.
  const std::vector<int>& getZeroKelvinBondOrders() const;
  
  /// \brief Return the free electron content of atoms in a representative Lewis structure of the
  ///        molecules in this topology.
  const std::vector<int>& getZeroKelvinFreeElectrons() const;

  /// \brief Return the ring inclusion specifications for one or more atoms in the system.
  ///
  /// Overloaded:
  ///   - Get the ring inclusion for a specific atom
  ///   - Get the ring inclusion for a range of atoms
  ///   - Get the ring inclusion for all atoms
  ///
  /// \param atom_index  The one atom of interest
  /// \param low_index   The start of a series of atoms for which to get ring inclusions
  /// \param high_index  The upper limit of a series of atoms for which to get ring inclusions
  /// \{
  ullint getRingInclusion(int atom_index) const;
  std::vector<ullint> getRingInclusion(int low_index, int high_index) const;
  std::vector<ullint> getRingInclusion() const;
  /// \}

  /// \brief Return whether a bond in the system is in a ring group or not.
  ///
  /// Overloaded:
  ///   - Provide the two atoms at either endpoint of the bond (checked for existence of an actual
  ///     bond)
  ///   - Provide the topological index of the bond (faster method if a means of interpreting the
  ///     context of each bond is prepared)
  ///
  /// \param atom_i      Topological index of the first atom in the bond
  /// \param atom_j      Topological index of the second atom in the bond
  /// \param bond_index  Topological index of the bond of interest
  /// \{
  bool bondIsInRing(int atom_i, int atom_j) const;
  bool bondIsInRing(int bond_index) const;
  /// \}
  
  /// \brief Get the atom endpoints of a rotatable bond.  The bond root atom is returned in the x
  ///        member of the tuple, the pivot atom (the second atom, closest to atoms that will turn)
  ///        is listed in the y member.
  ///
  /// Overloaded:
  ///   - Get all rotatable bonds without re-ordering the list.
  ///   - Get a list of all rotatable bonds upon which a minimum number of atoms turn.
  ///   - Get a list of rotatable bonds for which the pivot is within a specific cutoff of the
  ///     current conformation's center of mass.
  ///   - Get a list of rotatable bonds for which the pivot is within a specific cutoff of the
  ///     center of mass of some atom mask (the mask must be supplied as a raw bitmask of the
  ///     entire system, a std::vector of unsigned ints, to prevent making a circular dependency
  ///     whereby the AtomMask object depends on ChemicalFeatures and vice-versa).
  ///   - Get a list of the N largest rotatable groups of atoms.
  ///
  /// \param cutoff     The threshold at which to accept rotatable bond groups.  The meaning
  ///                   depends on the value of the choice enumeration (see below).  If the choice
  ///                   is COM_PROXIMITY, then cutoff is a distance with units of Angstroms.  If
  ///                   the choice is GROUP_SIZE, then cutoff is a minimium number of rotating
  ///                   atoms.
  /// \param mol_index  The molecule of interest (the system may have multiple molecules).
  /// \{
  std::vector<IsomerPlan> getRotatableBondGroups() const;
  std::vector<IsomerPlan> getRotatableBondGroups(int cutoff, int mol_index = 0) const;
  /// \}

  /// \brief Get the moving atom groups and bond endpoints involved in cis-trans isomerization.
  ///        This function clones the simple form of getRotatableBondGroups but the reliance on
  ///        so many different arrays and counters makes it fruitless to abstract.
  std::vector<IsomerPlan> getCisTransIsomerizationGroups() const;
  
  /// \brief Get the group of atoms that can invert a chiral center by a C2 symmetry rotation.
  ///        This re-uses the IsomerPlan struct, this time casting the root and pivot atoms as
  ///        the origins of the heaviest and second-heaviest branches, respectively, which will
  ///        again not move even as the rest of the atoms rotate.
  std::vector<IsomerPlan> getChiralInversionGroups() const;
  
  /// \brief Get the means for inverting one or more chiral centers.
  ///
  /// Overloaded:
  ///   - Get all chiral inversion instructions
  ///   - Get the instruction for a specific center, numbered according to the list of all chiral
  ///     centers in the molecule.
  ///
  /// \param index  The index of the chiral center of interest
  /// \{
  std::vector<ChiralInversionProtocol> getChiralInversionMethods() const;
  ChiralInversionProtocol getChiralInversionMethods(int index) const;
  /// \}

  /// \brief Get a pointer to the AtomGraph which built this object.
  const AtomGraph* getTopologyPointer() const;

  /// \brief Get a const pointer to the object itself in host memory.
  const ChemicalFeatures* getSelfPointer() const;
  
  /// \brief Get the abstract.
  ///
  /// \param tier  Extract pointers to data on either the CPU host or GPU device
  ChemicalFeaturesReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  
  /// \brief Find the chiral orientations of the system's chiral centers.  This function will be
  ///        called once by the constructor but can be called additional times if the coordinates
  ///        of the system change.  It will update the object's internal array of chiral centers,
  ///        which combines orientational information with the centers' topological indices.
  ///
  /// \param cfr  Present coordinates of the system
  void findChiralOrientations(const CoordinateFrameReader &cfr);

  /// \brief Find the rotatable (and cis-trans invertible, and chiral invertible) atom groups in
  ///        the topology based on rotatable bonds or chiral centers found by the rest of the
  ///        chemical perception calculations.
  ///
  /// \param timer  Timekeeping object for performance profiling
  void findRotatableBondGroups(StopWatch *timer = nullptr);
  
private:
  int atom_count;              ///< Number of atoms in the system
  int planar_atom_count;       ///< Number of atoms at the centers of improper dihedrals
  int ring_count;              ///< Number of unique ring systems found in the topology
  int fused_ring_count;        ///< Number of fused rings in the system
  int twistable_ring_count;    ///< Number of rings which can undergo conformational changes, i.e.
                               ///<   boat / chair.  The conformations are not called boat or chair
                               ///<   per se, but provided in tabulated lists based on the size of
                               ///<   the ring.  This indicates the number of such rings, which
                               ///<   will be a subset of all rings (excluding those which are
                               ///<   aromatic).
  int conjugated_group_count;  ///< Number of conjugated systems found in the topology
  int aromatic_group_count;    ///< Number of aromatic groups found in the topology
  int polar_hydrogen_count;    ///< Number of polar hydrogens in the system, across all molecules
  int hbond_donor_count;       ///< Number of hydrogen bond donors, across all molecules
  int hbond_acceptor_count;    ///< Number of hydrogen bond acceptors, across all molecules
  int chiral_center_count;     ///< Number of chiral centers in the topology
  int rotatable_bond_count;    ///< Number of fully rotatable single bonds
  int cis_trans_bond_count;    ///< Number of cis-trans invertible double bonds
  int double_bond_count;       ///< Number of detected double bonds
  int triple_bond_count;       ///< Number of detected triple bonds
  int max_ring_size;           ///< The maximum size of a ring, based on the size of a long long
                               ///<   unsigned integer (stored for reference)
  double temperature;          ///< The temperature at which these chemical features were
                               ///<   determined (this can influence the Lewis structure)
  bool rotating_groups_mapped; ///< Flag to indicate that rotating groups have been mapped
  bool chiralities_computed;   ///< Flag to indicate that chiral centers have been identified and
                               ///<   designated (R- or S-)
  
  /// List of atoms which constitute planar centers (needs no bounds array, simply a list of
  /// unique atoms at the centers of improper dihedrals)
  Hybrid<int> planar_centers;
  
  /// Indicates, for each atom i in the nth bit, whether the atom is part of an n-membered ring.
  /// Rings of up to 60 atoms can be tracked in this manner.  The highest four bits are reserved
  /// to indicate how many unique rings an atom takes part in (up to 15, which should be chemically
  /// impossible except in the case of a bonded crystal).  For atoms in chemically bonded crystal
  /// lattices, the ring numbers could technically explode (every atom being considered part of a
  /// multiple higher-order rings.  For this reason, the algorithm for detecting rings will be
  /// recursive but always refer back to the array itself for the current status of the ring
  /// occupancy.  It will cut off if all 60 slots have been filled, and declare the number of
  /// simultaneously occupied rings to be 15 (a sort of inf representation).
  Hybrid<ullint> ring_inclusion;

  /// Bounds for the rings that have been found
  Hybrid<int> ring_atom_bounds;
  
  /// A list of all rings in the system.  The starting point of the kth ring is given by
  /// ring_atom_bounds[k] and the size of the ring can be computed by referencing
  /// ring_atom_bonds[k + 1].
  Hybrid<int> ring_atoms;
  
  /// Bounds array for the aromatic atoms list below
  Hybrid<int> aromatic_group_bounds;

  /// Counts of pi electrons in each detected aromatic group
  Hybrid<int> aromatic_pi_electrons;
  
  /// List of atoms indicating groups of aromaticity in the system
  Hybrid<int> aromatic_groups;

  /// List of polar hydrogen atoms
  Hybrid<int> polar_hydrogens;
  
  /// List of hydrogen bond donor atoms (while the order of donors and actual polar hydrogens
  /// may be similar, the fact that any one donor could have more than one polar hydrogen attached
  /// to it prevents there being a 1:1 mapping between the arrays)
  Hybrid<int> hydrogen_bond_donors;

  /// List of hydrogen bond acceptor atoms
  Hybrid<int> hydrogen_bond_acceptors;

  /// List of topological indices for each of four atoms at the roots of chiral arms for each
  /// chiral center.  The x member contains the lowest priority arm (the one which orients the
  /// wheel), the y, z, and w members contain the highest, second highest, and third highest
  /// priority arms.
  Hybrid<int4> chiral_arm_atoms;
  
  /// List of chiral centers (needs no bounds array).  Positive values in this list indicate that
  /// an atom is L-chiral, while negative values indicate that an atom with the index equal to the
  /// absolute value is D-chiral.
  Hybrid<int> chiral_centers;

  /// List of chiral inversion capabilities.  The integers are tied to the chemistry enumeration
  /// ChiralInversionProtocol, spanning "do not invert" to "invert by full reflection and reverse
  /// all other centers to compensate."
  Hybrid<int> chiral_inversion_methods;
  
  /// List of rotating atoms, and the endpoints of bonds defining the axes of rotation.  The first
  /// atom in each group is the "root," distal to the atoms that move from the "pivot," the second
  /// atom listed in each group.  Neither the root atom or the pivot actually moves, but subsequent
  /// atoms in each group do.  Use the rotatable_group_bounds array to determine the extent of each
  /// group.
  Hybrid<int> rotatable_groups;

  /// Bounds array for rotatable groups.  Each group will have at least three members: the root
  /// atom, the pivot atom, and one or more atoms beyond the pivot that rotate as a consequence of
  /// twisting about the rotatable bond.
  Hybrid<int> rotatable_group_bounds;

  /// List of atoms involved in cis-trans isomerization about relevant double bonds.  The format
  /// follows rotatable_groups, above.
  Hybrid<int> cis_trans_groups;

  /// Bounds array for cis_trans_groups.  The format and size considerations follow form
  /// rotatable_group_bounds, above.
  Hybrid<int> cis_trans_group_bounds;
  
  /// List of atoms that move (by a C2 symmetry operation) to invert a chiral center.  These groups
  /// are chosen to minimize the number of moving atoms involved in flipping a chiral center.  For
  /// some molecules like proteins, the list can get very long, so a flag is provided for disabling
  /// the search if it is not needed.
  Hybrid<int> invertible_groups;
  
  /// Bounds array for chiral inversion groups.  Each group will have at least two members,
  /// covering two branches of the chiral center that will undergo a C2 symmetry operation in order
  /// to invert the center.
  Hybrid<int> invertible_group_bounds;

  /// First atoms of the heaviest chain branching from each chiral center
  Hybrid<int> anchor_a_branches;
  
  /// First atoms of the second-heaviest chain branching from each chiral center
  Hybrid<int> anchor_b_branches;
  
  /// Formal charges determined for all atoms in the system, based on a Lewis structure drawing
  /// with the Indigo method
  Hybrid<double> formal_charges;

  /// Bond orders determined for all bonds in the system, based on a Lewis structure drawing with
  /// the Indigo method
  Hybrid<double> bond_orders;

  /// Free electron content of each atom in the molecule, averaged over resonance states
  Hybrid<double> free_electrons;
  
  /// Integer data (the Hybrid<int> arrays above are POINTER-kind and will targer this storage)
  Hybrid<int> int_data;

  /// Mutable group data.  This is also integer data, targeted by arrays that pertain to moving
  /// groups of atoms.  The presence of this array, separate from other integer data, makes it
  /// feasible to remap the rotatable bond groups and inversion centers.
  Hybrid<int> mutable_group_data;

  /// Double-precision data, targeted by the POINTER-kind Hybrid<double> objects above
  Hybrid<double> double_data;

  /// Formal charges of a representative Lewis structure preserved from the octet-rule-based
  /// methods used to construct the resonance-averaged formal charges, bond orders, and free
  /// electron counts on each atom.  This is useful in printing SD files and other formats that
  /// require formal charges and bond orders, but inaccessible to the GPU as it is much less
  /// useful in subsequent calculations.
  std::vector<int> zerok_formal_charges;

  /// Bond orders of the representative Lewis structure
  std::vector<int> zerok_bond_orders;

  /// Free electron content, per atom, of the representative Lewis structure
  std::vector<int> zerok_free_electrons;

  /// Array to indicate whether each bond in the topology is part of any ring system.
  std::vector<bool> bond_in_ring;
  
  /// Pointer to the topology that this object describes (needed, if not just for reference, to
  /// obtain the actual atom indices by the various bonds enumerated in some of the lists above)
  AtomGraph *ag_pointer;

  /// Pointer to a timings object that, if not nullptr, will record the operations of this object
  /// and help profile its contributions to the total wall time.
  StopWatch *timer;

  /// \brief Impart the details of a topology onto the object.  This includes detecting which atoms
  ///        will have chirality and how to invert those centers if necessary, but does not compute
  ///        the chirality of those atoms (use the findChiralOrientations() member function).
  ///
  /// \param map_groups_in  Indicate whether to map rotable bond (and invertible center) groups
  void analyzeTopology(MapRotatableGroups map_groups_in);
  
  /// \brief Find all planar atoms in the system based on improper dihedral terms.  This prepares
  ///        data arrays that will later be loaded into the ChemicalFeatures object, but does not
  ///        directly modify the object's members.
  ///
  /// \param vk  Abstract of valence terms from the topology
  std::vector<int> findPlanarAtoms(const ValenceKit<double> &vk) const;

  /// \brief Trace all rings in the system based on a tree structure linked list.
  ///
  /// \param nbk                   Nonbonded abstract from the original topology
  /// \param cdk                   Atom, residue, and molecule details taken from the original
  ///                              topology
  /// \param tmp_ring_inclusions   Developing list of ring inlusions.  The nth bit of the kth
  ///                              array element indicates whether the kth atom participates in a
  ///                              ring of size n.
  /// \param tmp_ring_atoms        Growing list of atoms in all rings found thus far
  /// \param tmp_ring_atom_bounds  Bounds array for tmp_ring_atoms.  Must be initialized with a
  ///                              leading zero so that the push_back method can add the boundaries
  ///                              of each successive group and arrive at a capped, exclusive
  ///                              prefix sum of all ring sizes.
  void traceTopologicalRings(const NonbondedKit<double> &nbk, const ChemicalDetailsKit &cdk,
                             std::vector<ullint> *tmp_ring_inclusion,
                             std::vector<int> *tmp_ring_atoms,
                             std::vector<int> *tmp_ring_atom_bounds);
  
  /// \brief Draw a ring containing two atoms, based on the histories of each atom's branch.
  ///
  /// \param j_atom                The first atom that is definitely in the ring, links to k_atom
  /// \param k_atom                The second atom that is definitely in the ring, links to j_atom
  /// \param tree_positions        Positions of various atoms in the linked list of BondedNode
  ///                              objects
  /// \param node_count            The number of nodes in the linked list object
  /// \param links                 The linked list object storing information on what atoms bond
  /// \param tmp_ring_inclusion    Array recording whether each atom takes part in a ring of size
  ///                              n, as shown by the nth bit of each ullint 
  /// \param tmp_ring_atoms        Growing list of known ring-forming atoms
  /// \param tmp_ring_atom_bounds  Bounds array for tmp_ring_atoms
  void markRingAtoms(int j_atom, int k_atom, const std::vector<int> &tree_positions,
                     int node_count, std::vector<BondedNode> *links,
                     std::vector<ullint> *tmp_ring_inclusion, std::vector<int> *tmp_ring_atoms,
                     std::vector<int> *tmp_ring_atom_bounds,
                     const ChemicalDetailsKit &cdk);

  /// \brief Draw Lewis structures over the entire topology using the Indigo method.  Lewis
  ///        structures will be drawn for each unique molecule and copied otherwise, but that
  ///        requires a list of all unique molecules.
  ///
  /// \param vk   Information on valence interactions, taken from the original topology
  /// \param nbk  Information on non-bonded interactions, taken from the original topology
  /// \param cdk  Atom, residue, and molecule details taken from the original topology
  void drawLewisStructures(const ValenceKit<double> &vk, const NonbondedKit<double> &nbk,
                           const ChemicalDetailsKit &cdk);

  /// \brief Find groups of aromatic atoms.  This will also enumerate the number of fused rings in
  ///        the system as such systems must be considered for aromaticity even if their component
  ///        rings fall short.
  ///
  /// \param cdk                        Chemical details abstract from the original topology, with
  ///                                   atomic numbers of all atoms
  /// \param vk                         Valence term abstract from the original topology
  /// \param tmp_ring_atoms             List of atoms involved in rings.  The associated bounds
  ///                                   array is not needed in this case, as the function just
  ///                                   steps over all atoms that are in some ring.
  /// \param tmp_ring_atom_bounds       Bounds array for each ring in the topology
  /// \param tmp_aromatic_groups        Growing list of aromatic groups detected in the topology
  /// \param tmp_aromatic_group_bounds  Developing bounds array for aromatic groups.  Must be
  ///                                   pre-intialized with a leading zero to allow use of the
  ///                                   push_back method to grow the array.
  void findAromaticGroups(const ChemicalDetailsKit &cdk, const ValenceKit<double> &vk,
                          const std::vector<int> &tmp_ring_atoms,
                          const std::vector<int> &tmp_ring_atom_bounds,
                          std::vector<int> *tmp_aromatic_group_bounds,
                          std::vector<int> *tmp_aromatic_pi_electrons,
                          std::vector<int> *tmp_aromatic_groups);

  /// \brief Find chiral centers in the system.  This will use methods similar to the ring
  ///        detection system.  Apply IUPAC rules based on substituent atomic numbers and bond
  ///        orders.  Isotopes will not be applied.  This returns a list of detected chiral
  ///        centers and assumes them all to be L-.  Atoms at the bases of their four distinct
  ///        arms are recorded and posted to one of the internal object's ARRAY-kind Hybrids, for
  ///        later reference.  The true orientations of each center are detected in a subsequent
  ///        call to the public member function findChiralOrientations().
  ///
  /// \param nbk                   Nonbonded system details abstract from the original topology
  /// \param vk                    Valence term abstract from the original topology
  /// \param cdk                   Chemical details abstract from the original topology
  std::vector<int> detailChiralCenters(const NonbondedKit<double> &nbk,
                                       const ValenceKit<double> &vk,
                                       const ChemicalDetailsKit &cdk);

  /// \brief Determine whether each chiral center can be inverted, and whether that is best done by
  ///        rotating two chiral branches by C2 symmetry or reflecting the entire molecule across a
  ///        plane.  Reflection across a plane inverts all chiral centers, so other centers which
  ///        can be flipped with C2 rotations must be reset to their original configurations
  ///        following a reflection.  If there are multiple chiral atoms requiring complete
  ///        reflection, their centers can only be flipped simultaneously, so only the first such
  ///        center in a molecule will be considered invertible at all.
  ///
  /// \param tmp_chiral_centers    List of previously detected chiral centers
  /// \param tmp_ring_atoms        List of all atoms in rings
  /// \param tmp_ring_atom_bounds  Bounds array for tmp_ring_atoms
  std::vector<int> findChiralInversionMethods(const std::vector<int> &tmp_chiral_centers,
                                              const std::vector<int> &tmp_ring_atoms,
                                              const std::vector<int> &tmp_ring_atom_bounds);

  /// \brief Find rotatable bonds in the system, those with bond order of 1.0 and nontrivial groups
  ///        sprouting from either end, and return a vector of the atom indices at either end.
  ///        This routine also traces cis- and trans-groups that can invert if the bond "rotates"
  ///        180 degrees.
  ///
  /// \param vk                          Valence term abstract from the original topology
  /// \param cdk                         Chemical details of the system (for atomic numbers)
  /// \param nbk                         Nonbonded system details, abstract taken from the
  ///                                    original topology
  /// \param ring_atoms                  List of atoms in rings, indexing the original topology
  /// \param ring_atom_bounds            Bounds array for the ring_atoms array
  /// \param tmp_rotatable_groups        List of all atoms involved in rotation about a rotatable
  ///                                    bond, including the endpoints of the bond itself in the
  ///                                    first two slots.  This vector is assembled and returned.
  /// \param tmp_rotatable_group_bounds  Bounds array for tmp_rotatable_groups, assembled and
  ///                                    returned
  /// \param tmp_cis_trans_groups        List of all atoms involved in flipping the cis- or trans-
  ///                                    nature of a double bond.  This vector is assembled and
  ///                                    returned.
  /// \param tmp_cis_trans_group_bounds  Bounds array for tmp_cis_trans_groups, assembled and
  ///                                    returned
  void findRotatableBonds(const ValenceKit<double> &vk, const ChemicalDetailsKit &cdk,
                          const NonbondedKit<double> &nbk, const std::vector<int> &ring_atoms,
                          const std::vector<int> &ring_atom_bounds,
                          std::vector<int> *tmp_rotatable_groups,
                          std::vector<int> *tmp_rotatable_group_bounds,
                          std::vector<int> *tmp_cis_trans_groups,
                          std::vector<int> *tmp_cis_trans_group_bounds);

  /// \brief Find invertible groups in the system, those comprising two branches of a chiral center
  ///        that have the fewest possible atoms.  The number of invertible groups is the number of
  ///        chiral centers.
  ///
  /// \param tmp_chiral_centers           Array of pre-determined chiral centers
  /// \param tmp_inversion_methods        Array of pre-determined inversion protocols
  /// \param tmp_anchor_a_branches        Roots of the largest branch on each chiral center (in
  ///                                     terms of the number of atoms), not to be inverted
  /// \param tmp_anchor_b_branches        Roots of the second largest branch on each chiral center,
  ///                                     not to be inverted
  /// \param tmp_invertible_groups        Array of atoms in invertible groups (assembled and
  ///                                     returned)
  /// \param tmp_invertible_group_bounds  Bounds array for tmp_invertible_groups (assembled and
  ///                                     returned)
  void findInvertibleGroups(const std::vector<int> &tmp_chiral_centers,
                            const std::vector<int> &tmp_inversion_methods,
                            std::vector<int> *tmp_anchor_a_branches,
                            std::vector<int> *tmp_anchor_b_branches,
                            std::vector<int> *tmp_invertible_groups,
                            std::vector<int> *tmp_invertible_group_bounds);
  
  /// \brief Find polar heavy atoms that can act as hydrogen bond donors, and label polar hydrogens
  ///        in the process.  Unlike findRotatableBonds above, this is a fast evaluation and will
  ///        be performed automatically with the ChemicalFeatures object construction.
  ///
  /// \param nbk          Non-bonded interaction abstract from the original topology
  /// \param cdk          Chemical details of the system (for atomic numbers)
  /// \param tmp_polar_h  Array for listing polar hydrogen atoms (modified and returned)
  /// \param tmp_hb_don   Array for listing hydrogen bond donor atoms (modified and returned)
  /// \param tmp_hb_acc   Array for listing hydrogen bond acceptor atoms (modified and returned)
  void findHydrogenBondElements(const NonbondedKit<double> &nbk, const ChemicalDetailsKit &cdk,
                                std::vector<int> *tmp_polar_h, std::vector<int> *tmp_hb_don,
                                std::vector<int> *tmp_hb_acc);

  /// \brief Allocate and assign POINTER-kind Hybrid objects to the mutable groups data array.  The
  ///        function accepts various concatenated lists of atom groups (plus their respective
  ///        bounds arrays), allocates the appropriate amount of Hybrid int data, then stores the
  ///        arrays in eponymous POINTER-kind Hybrid members of the ChemicalFeatures object
  ///        targeting the newly allocated storage array.
  ///
  /// \param tmp_rotatable_groups         List of atoms involved in rotating groups
  /// \param tmp_rotatable_group_bounds   Bounds array for tmp_rotatable_groups
  /// \param tmp_cis_trans_groups         List of atoms involved in cis-trans flips
  /// \param tmp_cis_trans_group_bounds   Bounds array for tmp_cis_trans_groups
  /// \param tmp_invertible_groups        List of atoms involve din chiral inversions
  /// \param tmp_invertible_group_bounds  Bounds array for tmp_invertible_groups
  /// \param tmp_anchor_a_branches        List of chiral centers' highest-priority branches
  /// \param tmp_anchor_b_branches        List of chiral centers' second-highest branches
  void allocateMutableData(const std::vector<int> &tmp_rotatable_groups,
                           const std::vector<int> &tmp_rotatable_group_bounds,
                           const std::vector<int> &tmp_cis_trans_groups,
                           const std::vector<int> &tmp_cis_trans_group_bounds,
                           const std::vector<int> &tmp_cis_invertible_groups,
                           const std::vector<int> &tmp_cis_invertible_group_bounds,
                           const std::vector<int> &tmp_anchor_a_branches,
                           const std::vector<int> &tmp_anchor_b_branches);

  /// \brief Reset POINTER-kind Hybrid objects to target the appropriate ARRAY-kind object in
  ///        the copy constructor and copy assignment operator.
  void repairPointers();
};

/// \brief Concatenate the lists of rotating groups that define rotatable bonds as well as
///        cis-trans isomerization.
///
/// \param bond_endpoints    List of endpoints for the relevant bonds, root atom (x member) and
///                          pivot atom (y member)
/// \param moving_lists      The list of lists for minimal rotating groups pertaining to all
///                          relevant bonds (rotatable or cis-trans invertible)
/// \param tmp_groups        Concatenated list of atoms involved in each rotatable group,
///                          assembled and returned
/// \param tmp_group_bounds  Bounds of rotatable or cis-trans invertible groups, assembled and
///                          returned
void unpackRotatableGroups(const std::vector<int2> &bond_endpoints,
                           const std::vector<std::vector<int>> &moving_lists,
                           std::vector<int> *tmp_groups, std::vector<int> *tmp_group_bounds);

/// \brief Score the four branches of a chiral molecule.  This is called by the findChiralCenters()
///        member function of the ChemicalFeatures object, but written as a free function as it
///        would not benefit from any of the object's member variables.  Returns true if scoring
///        indicates that more branch exploration is needed and might successfully discriminate
///        among the four branches, false otherwise.
///
/// \param links             Quartet of trees describing branches out of the putative chiral
///                          center.  Also contains information on bond orders.
/// \param layer_llim        Lower bounds for each branch's atoms added since the previous pass
/// \param layer_hlim        Upper bounds for each branch's atoms added since the previous pass
/// \param cdk               Chemical details of the system, an abstract taken from the original
///                          topology.  Contains atomic numbers for all atoms.
/// \param chiral_dominance  Dominance matrix describing whether one branch beats another.
///                          Modified and returned.
/// \param parallel_growth   Matrix indicating whether each branch grows in parallel with another
bool scoreChiralBranches(const std::vector<std::vector<BondedNode>> &links,
                         const std::vector<int> &layer_llim, const std::vector<int> &layer_hlim,
                         const ChemicalDetailsKit &cdk, std::vector<int> *chiral_dominance,
                         std::vector<int> *parallel_growth);

/// \brief Determine the orientation of a chiral center given the priorities of its various
///        branches.  Only the first atom of each branch is critical at this point.  Return +1 for
///        L-chiral centers (S-) or -1 for D-chiral (R-) centers.  This is then used as a
///        multiplier for the chiral center atom index in subsequent masks or atom retrievals.
///
/// Overloaded:
///   - Supply any type of coordinate object, by const pointer or by const reference
///   - Compute in the native precision level of the object (the demarcation between L- and D-
///     chirality should be sufficiently large in nearly all cases that there is no chance of
///     ambiguity based on the precision model)
///
/// \param cfr            Coordinates of the entire system
/// \param cf             Coordinates of the entire system
/// \param psr            Coordinates of the entire system
/// \param ps             Coordinates of the entire system
/// \param center_atom    Index of the center atom
/// \param root_atom      Index of the "root" atom--the very lowest priority branch
/// \param branch_a_atom  Index of the highest priority branch
/// \param branch_b_atom  Index of the next highest priority branch
/// \param branch_c_atom  Index of the lowest priority branch, aside from the root branch
/// \{
template <typename T>
int getChiralOrientation(const T* xcrd, const T* ycrd, const T* zcrd, int center_atom,
                         int root_atom, int branch_a_atom, int branch_b_atom, int branch_c_atom);
                         
int getChiralOrientation(const CoordinateFrameReader &cfr, int center_atom, int root_atom,
                         int branch_a_atom, int branch_b_atom, int branch_c_atom);

int getChiralOrientation(const CoordinateFrame *cf, int center_atom, int root_atom,
                         int branch_a_atom, int branch_b_atom, int branch_c_atom);

int getChiralOrientation(const CoordinateFrame &cf, int center_atom, int root_atom,
                         int branch_a_atom, int branch_b_atom, int branch_c_atom);

int getChiralOrientation(const PhaseSpaceReader &psr, int center_atom, int root_atom,
                         int branch_a_atom, int branch_b_atom, int branch_c_atom);

int getChiralOrientation(const PhaseSpace *ps, int center_atom, int root_atom, int branch_a_atom,
                         int branch_b_atom, int branch_c_atom);

int getChiralOrientation(const PhaseSpace &ps, int center_atom, int root_atom, int branch_a_atom,
                         int branch_b_atom, int branch_c_atom);

template <typename T>
int getChiralOrientation(const CoordinateSeriesReader<T> &csr, size_t frame_index, int center_atom,
                         int root_atom, int branch_a_atom, int branch_b_atom, int branch_c_atom);

template <typename T>
int getChiralOrientation(const CoordinateSeries<T> *cs, size_t frame_index, int center_atom,
                         int root_atom, int branch_a_atom, int branch_b_atom, int branch_c_atom);

template <typename T>
int getChiralOrientation(const CoordinateSeries<T> &cs, size_t frame_index, int center_atom,
                         int root_atom, int branch_a_atom, int branch_b_atom, int branch_c_atom);

int getChiralOrientation(const CondensateReader &cdnsr, int system_index, int center_atom,
                         int root_atom, int branch_a_atom, int branch_b_atom, int branch_c_atom);

int getChiralOrientation(const Condensate *cdns, int system_index, int center_atom,
                         int root_atom, int branch_a_atom, int branch_b_atom, int branch_c_atom);

int getChiralOrientation(const Condensate &cdns, int system_index, int center_atom,
                         int root_atom, int branch_a_atom, int branch_b_atom, int branch_c_atom);

int getChiralOrientation(const PsSynthesisReader &poly_psr, int system_index, int center_atom,
                         int root_atom, int branch_a_atom, int branch_b_atom, int branch_c_atom);

int getChiralOrientation(const PhaseSpaceSynthesis *poly_ps, int system_index, int center_atom,
                         int root_atom, int branch_a_atom, int branch_b_atom, int branch_c_atom);

int getChiralOrientation(const PhaseSpaceSynthesis &poly_ps, int system_index, int center_atom,
                         int root_atom, int branch_a_atom, int branch_b_atom, int branch_c_atom);
/// \}

/// \brief Beginning with two distinct atoms in a molecule within a topology, proceed throughout
///        the molecular bonding pattern verifying that the atoms one encounters when stepping
///        outward from each atom have identical properties.
///
/// Overloaded:
///   - Accept a ChemicalFeatures object from which to extract critical information.  This calls
///     previous definitions of the function found in match_bonding_pattern.h.
///
/// \param ag              System topology
/// \param chemfe          Chemical details of the system, containing formal charges, bond orders,
///                        free electron counts, and ring inclusions
/// \param atom_a          The first atom to compare   
/// \param atom_b          The second atom to compare
/// \{
bool matchBondingPattern(const AtomGraph &ag, const ChemicalFeatures &chemfe, int atom_a,
                         int atom_b);
/// \}

/// \brief Determine whether two ways of isomerizing a molecule are linked, whether by sharing an
///        atom in the chiral center or rotatable bond, or by having any of those atoms be, in
///        turn, bonded to one another.
///
/// \param isomerizers  List of atom groups involved in each isomer creation for the molecule
/// \param permi        The first origin of isomerization and thus permutations
/// \param permj        The second origin of isomerization and thus permutations
/// \param nbk          Contains the non-bonded exclusions, including bonded exclusions
bool permutationsAreLinked(const std::vector<IsomerPlan> &isomerizers, int permi, int permj,
                           const NonbondedKit<double> &nbk);

} // namespace chemistry
} // namespace stormm

#include "chemical_features.tpp"

#endif
