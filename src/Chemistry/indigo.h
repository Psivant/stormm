// -*-c++-*-
#ifndef STORMM_INDIGO_H
#define STORMM_INDIGO_H

#include <map>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "DataTypes/common_types.h"
#include "DataTypes/mixed_types.h"
#include "Topology/atomgraph.h"

namespace stormm {
namespace chemistry {

using card::Hybrid;
using card::HybridTargetLevel;
using topology::AtomGraph;

/// \brief The range of formal charges that STORMM will handle.  This is well beyond the range of 
///        formal charges that can exist on bio-organic molecules.
/// \{
constexpr int minimum_formal_charge = -8;
constexpr int maximum_formal_charge =  7;
constexpr int formal_charge_range   = maximum_formal_charge - minimum_formal_charge + 1;
constexpr int formal_charge_bits    = 4;
/// \}

/// \brief The maximum energy spacing for Indigo to consider is set at 10 kcal/mol by default.
constexpr int default_indigo_energy_gap = 3187;

/// \brief The maximum number of negatively charged carbon atoms that a fragment can have (some
///        corner case might lead to this limit being increased under special circumstances)
constexpr int maximum_negative_carbons = 1;
  
/// \brief The maximum bond order is a triple bond
constexpr int maximum_bond_order = 3;
constexpr int bond_order_range   = maximum_bond_order + 1;
  
/// \brief Indigo is only geared to handle elements through Bromine at present
constexpr int maximum_indigo_atomic_number = 35;
constexpr int indigo_atomic_number_range   = maximum_indigo_atomic_number + 1;
  
/// \brief Create an unsigned int key for an atomic formal charge.  These numbers are taken from
///        Table S2 of the Supplement for the original Indigo paper:
///
/// Welsh ID and Allison JR. "Automated simultaneous assignment of bond orders and formal charges."
/// (2019) J. Chem. Informatics 11:18.  D.O.I. 10.1186/s13321-019-0340-0
///
/// https://jcheminf.biomedcentral.com/articles/10.1186/s13321-019-0340-0#Sec2
///
/// \param atomic_number  The atomic number under consideration
/// \param formal_charge  The charge of the element
int IndigoFormalChargeKey(int atomic_number, int formal_charge);

/// \brief Create an unsigned int key for a bond of a given order between two elements.  These
///        numbers are taken from Table S6 of the Supplement for the original Indigo paper:
///
/// Welsh ID and Allison JR. "Automated simultaneous assignment of bond orders and formal charges."
/// (2019) J. Chem. Informatics 11:18.  D.O.I. 10.1186/s13321-019-0340-0
///
/// https://jcheminf.biomedcentral.com/articles/10.1186/s13321-019-0340-0
///
/// \param atomic_number_i  The atomic number of the first element in the bond
/// \param atomic_number_j  The atomic number of the second element in the bond
/// \param bond_order_in    The bond order
int IndigoBondOrderKey(int atomic_number_i, int atomic_number_j, int bond_order);

/// \brief Create a table for Indigo's atomic formal charge assignments
std::vector<int> indigoFormalChargeScores();
  
/// \brief Create a table for Indigo's bond order assignments
std::vector<int> indigoBondOrderScores();

/// \brief Enumerate the possible combinations of formal charge state and the orders of bonds
///        relevant to a particular atom.  Also list the other atoms linked by those bonds.
class IndigoAtomCenter {
public:

  /// \brief The constructor takes the parent IndigoTable's lists of atoms and bonds
  IndigoAtomCenter(int table_index_in, int z_number_in, int bond_count_in,
                   int valence_electrons_in, const std::vector<int> &relevant_bonds_in,
                   const std::vector<int> &partner_atoms_in, const std::vector<int> &atom_scores,
                   const std::vector<int> &bond_scores);

  /// \brief Let it be known which atom this center pertains to
  int getIndigoTableIndex() const;
  
  /// \brief Return the number of states this atomc center can take on
  int getStateCount() const;

  /// \brief Return the bitmask describing a particular state of this atom center.
  ///
  /// \param state_index  Index of the state in question
  uint getState(int state_index) const;

  /// \brief Get the charge of an atom center in a particular state
  ///
  /// \param state_index  Index of the state in question
  int getCharge(int state_index) const;

  /// \brief Get the atomic number of the atom
  int getAtomicNumber() const;
  
  /// \brief Return the score associated with a particular state of this atom center.
  ///
  /// \param state_index  Index of the state in question
  int getScore(int state_index) const;

  /// \brief Sort the states to put the lowest energy (lowest scoring) state first
  void sortStates();
  
  /// \brief Get the number of partners in this atom center
  int getPartnerCount() const;

  /// \brief Get the index of the relevant bond connecting this atom to the partner of the same
  ///        index.  The index is produced into the original topology.
  ///
  /// \param index  Index of the atom center's partner, to which the bond of interest connects
  int getRelevantBond(int index) const;

  /// \brief Get one of the partners for a specific atom center.  Return the atom center index
  ///        of the partner in the overarching IndigoTable.
  ///
  /// \param index  Index of the specific partner within the center's local list
  int getPartner(int index) const;

  /// \brief Find a partner atom within the local list of this atom center.  Return the index of
  ///        the partner atom in the local list (also corresponds to the index of the bond that
  ///        this atom will use to connect to that partner atom).
  int findPartnerIndex(int partner_atom_number) const;

  /// \brief Get the order of a bond made by this atom center to a specified partner when the atom
  ///        center lies in a particular state.
  ///
  /// \param bond_index   Index of the bond in question, in the local list kept by this atom center
  ///                     This also indicates the local index of the partner of interest.
  /// \param state_index  Index of the state of the atom center to query.  Bond orders can change
  ///                     depending on different states of the atom center.
  int getBondOrderOfState(int bond_index, int state_index) const;  
  
  /// \brief Determine whether an atom center makes a bond of a certain order to a specific
  ///        partner.
  ///
  /// \param partner_index  The index of the partner atom center in the overarching IndigoTable
  /// \param partner_bo     The bond order that must be found in some state
  bool queryBondToPartner(int partner_atom_number, int partner_bo) const;

  /// \brief Cull states from a particular atom center if its partners cannot reciprocate with the
  ///        correct bond order.  Return the number of states culled.
  ///
  /// \param acen  Array of all atom centers from the overarching IndigoTable (this array contains
  ///              the atom center that is currently being operated on, but is passed as read-only)
  int cullByPartners(const std::vector<IndigoAtomCenter> &acen);
  
private:
  int table_index;                  ///< Index of this atom center's atom in the parent IndigoTable
  int z_number;                     ///< The atomic number of the atom center (this, in addition to
                                    ///<   the number of bonds and required number of valence
                                    ///<   electrons, will determine the number of states that the
                                    ///<   center can take on)
  int bond_count;                   ///< The number of bonds to and from the atom center,
                                    ///<   regardless of their orders
  int valence_electrons;            ///< Number of valence electrons that this center must maintain
                                    ///<   (2 for hydrogen, 10 for phosphorous with >=3 bonds,
                                    ///<   12 for sulfur with >= 3 bonds)
  int possible_states;              ///< The number of states that satisfy the valence requirements
                                    ///<   of the atom center (regardless of their energy)
  std::vector<int> relevant_bonds;  ///< Indices of bonds in the overarching IndigoTable whose
                                    ///<   orders are controlled by this atom center
  std::vector<int> partner_atoms;   ///< Atoms at the other end of each bond (same order as the
                                    ///<   relevant_bonds array)
  std::vector<uint> states;         ///< Bit-packed states of the atom center: the first four bits
                                    ///<   indicate the formal charge state of the atom, while
                                    ///<   subsequent segments of two bits each indicate the bond
                                    ///<   order state of the corresponding bonds from
                                    ///<   relevant_bonds
  std::vector<int> scores;          ///< Scores for the atom center, indicating twice the atom
                                    ///<   formal charge score and the sum of all bond order scores
                                    ///<   (double-counting the bond order scores in this way, as
                                    ///<   each will be shared by two centers, allows the
                                    ///<   configuration's total score to be divided by two to
                                    ///<   recover the estimated energy in Hartrees)

  /// \brief Compute the state bitmask for a particular formal charge and bonding configuration
  ///
  /// \param fc_value     Formal charge of the atom (it will be adjusted to a scale of [0, 16)
  ///                     before storing it in a 4-bit format)
  /// \param bond_orders  Settings for each bond going in and out of the atom center
  uint computeState(const int fc_value, const std::vector<uint> &bond_orders);

  /// \brief Compute the score for a particular formal charge and bonding configuration.  The
  ///        atom score contributes twice and each bond score contributes its full amount, so that
  ///        the final result of all atom centers' scores in particular (mutually compatible)
  ///        states can be divided by 2 to give the total molecular energy of that Lewis Structure.
  ///
  /// \param state        State of the system, encoded in a bit mask
  /// \param atom_scores  Scores for each atom's possible formal charge states
  /// \param bond_scores  Scores for each bond's possible orders
  int computeScore(const uint state, const std::vector<int> &atom_scores,
                   const std::vector<int> &bond_scores);
};

/// \brief A fragment of a structure detailed with the Indigo scoring function
class IndigoFragment {
public:

  /// \brief Constructor takes a list of atoms, then determines the bonds between them and
  ///        whatever connections to other fragments.  The identities of other fragments are not
  ///        determined at the time of construction, but space to store them is allocated.
  ///
  /// \param idg_tab
  /// \param atom_list    Subset of atoms that will make up this fragment
  /// \param score_delta  The maximum score above some ground state, for a particular charge
  ///                     increment, to continue catalogging new states.  Other above this
  ///                     threshold will be ignored.
  IndigoFragment(const std::vector<int> &centers_list_in,
                 const std::vector<IndigoAtomCenter> &all_centers,
                 int score_delta = default_indigo_energy_gap);

  /// \brief Get the number of centers in this fragment
  int getCenterCount() const;
  
  /// \brief Get the number of states in this fragment
  int getStateCount() const;
  
  /// \brief Return the vector of {atom center number, atom center state} tuples describing a
  ///        particular state of this fragment.  This list will have to be interpreted in light of
  ///        the complete list of atom centers to know what formal charge states and bond orders
  ///        the fragment state actually implies.
  ///
  /// \param state_index  Index of the fragment state in question
  std::vector<int2> getState(int state_index) const;

  /// \brief Get the charge of this fragment in a particular state.  There total charges were
  ///        compute when the object was constructed, hence no further references to a master list
  ///        of atom centers and their individual states is necessary.
  ///
  /// \param state_index  Index of the state in question
  int getCharge(int state_index) const;

  /// \brief Return the score associated with a particular state of this atom center.
  ///
  /// \param state_index  Index of the state in question
  int getScore(int state_index) const;

  /// \brief Create a vector of the unique charge states in a fragment
  std::vector<int> getChargeStates() const;

  /// \brief Get the range of charges and the minimum energy states that can produce each charge
  ///        offered by this fragment.  The vector of results is a list of tuples providing the
  ///        fragment charge (in the x member), the minimum energy state (in the y member), and
  ///        the index of the fragment state for that minimum energy and charge combination (in
  ///        the z member).
  ///
  /// \param state_index  Index of the state in question
  std::vector<int3> getChargesAndBestEnergies() const;

  /// \brief Get all states in a fragment bearing a particular charge.  The energy of each state
  ///        satisfying the net charge condition is returned in the x member of each tuple.  The
  ///        index of the state is returned in the y member.
  ///
  /// \param charge_value  The charge criterion for selecting states
  std::vector<int2> getStatesBearingCharge(int charge_value) const;

  /// \brief Assess whether two fragments are equivalent, just with different atom indices.
  ///        Return the result as a boolean expression.
  ///
  /// \param other          Another fragment
  /// \param real_atom_map  Map of fragment atoms into the original topology
  /// \param ag_pointer     Pointer to the original topology
  /// \param atom_centers   Array of all atom centers in the IndigoTable containing these fragments
  bool testEquivalence(const IndigoFragment &other, const std::vector<int> &real_atom_map,
                       const AtomGraph *ag_pointer,
                       const std::vector<IndigoAtomCenter> &atom_centers) const;

  /// \brief Cull states of a fragment that bear a specific charge.
  ///
  /// \param charge_value  Cull states bearing this net charge
  int cullStatesBearingCharge(int charge_value);  

  /// \brief Cull non-optimal states of a fragment that bear a specific charge.
  ///
  /// \param charge_value  Cull states bearing this net charge with sub-optimal energies
  /// \param score_delta   Threshold at which to cull non-optimal states (default 3187 implies an
  ///                      energy difference of 10 kcal/mol)
  int cullHigherEnergyStatesByCharge(int charge_value,
                                     int score_delta = default_indigo_energy_gap);
  
private:
  int center_count;               ///< The number of atoms in this fragment
  int possible_states;            ///< The number of states that this fragment, as a whole, can
                                  ///<   adopt.  This number will likely be considerably less than
                                  ///<   the product of all its component centers' state counts, as
                                  ///<   not all centers may be able to adopt any of their states
                                  ///<   independently.
  std::vector<int> centers_list;  ///< List of the atom centers making up this fragment
  std::vector<int*> states;       ///< List of the states that this fragment can adopt.  Each
                                  ///<   state is given as a series of unsigned integers detailing
                                  ///<   the states of its component atom centers.  This is a list
                                  ///<   of pointers to a larger vector of the actual unsigned int
                                  ///<   data, which makes indexing and sorting easier.  This array
                                  ///<   is possible_states in length.
  std::vector<int> states_data;   ///< The actual data for each atom center in each fragment state.
                                  ///<   This array is center_count x possible states in size.
  std::vector<int> net_charges;   ///< List of net charges for each fragment state (possible_states
                                  ///<   in length)
  std::vector<int> scores;        ///< List of scores for each fragment state (possible_states in
                                  ///<   length)
};
  
/// \brief Table of options and optimization targets for an Indigo formal charge and bond order
///        assignment.
class IndigoTable {
public:

  /// \brief The constructor takes a topology and options for the min and max formal charge (these
  ///        are merely to accelerate the optimization and search for possible values from the
  ///        respective maps).
  ///
  /// \param ag                     System topology
  /// \param molecule_index         Molecule within the topology to develop a Lewis structure for
  ///                               (if there are solvent molecules, this is a good way to ignore
  ///                               them)
  /// \param temperature            System temperature in Kelvin
  IndigoTable(const AtomGraph *ag_in, int molecule_index = 0, double temperature_in = 300.0);

  /// \brief Get the number of atoms  
  int getAtomCount() const;

  /// \brief Get the number of bonds
  int getBondCount() const;

  /// \brief Get the net charge on the system
  int getNetCharge() const;

  /// \brief Return the system's ground state formal charges in a form amenable to the original
  ///        topology.  The output is a list of tuples containg the topology atom index in the x
  ///        member (round and convert to nearest integer) and the formal charge in the y member.
  std::vector<CombineIDp> getGroundStateFormalCharges() const;

  /// \brief Return the system's ground state bond orders in a form amenable to the original
  ///        topology.  The output is a list of tuples containg the topology bond index in the x
  ///        member (round and convert to nearest integer) and the order in the y member.
  std::vector<CombineIDp> getGroundStateBondOrders() const;

  /// \brief Return the number of free electrons on each atom in the system's ground state.
  std::vector<CombineIDp> getGroundStateFreeElectrons() const;

  /// \brief Return the system's ground state formal charges, taking only a single representation
  ///        of each bond and formal charge state.  This equates to a single Lewis structure (out
  ///        of many), with the lowest or tied-for-lowest energy state at 0K.  This will produce
  ///        integral formal charges, suitable for an SD file.  The topological index of each atom
  ///        is given in the "x" member while the formal charge is given in the "y" member of each
  ///        tuple in the output array.
  std::vector<int2> getZeroKelvinFormalCharges() const;

  /// \brief Return the system's ground state bond orders, taking only a single representation of
  ///        each bond and formal charge state.  As with getZeroKelvinFormalCharges() above, this
  ///        will produce a single Lewis structure with integral bond orders.  The topological
  ///        index of each atom is given in the "x" member while the bond order is given in the "y"
  ///        member of each tuple in the output array.
  std::vector<int2> getZeroKelvinBondOrders() const;

  /// \brief Return the system's ground state electron content, taking only a single representation
  ///        of each bond and formal charge state.  As with getZeroKelvinFormalCharges() above,
  ///        this will produce a single Lewis structure with integral bond orders.  The topological
  ///        index of each atom is given in the "x" member while the bond order is given in the "y"
  ///        member of each tuple in the output array.
  std::vector<int2> getZeroKelvinFreeElectrons() const;
  
  /// \brief Return the system's ground state energy
  double getGroundStateEnergy() const;  
  
private:
  int atom_count;                       ///< Number of real H/C/N/O/S/P/F/Cl/Br atoms in the system
  int bond_count;                       ///< Number of bonds in the system
  int net_charge;                       ///< Net charge of the molecule (indicates the total number
                                        ///<   of excess electrons to place on the system)
  int fragment_count;                   ///< The number of fragments into which the molecule is
                                        ///<   broken
  double temperature;                   ///< Temperature at which to evaluate distributions of
                                        ///<   alternative states
  double ground_state_energy;           ///< Estimate of the ground state energy
  std::vector<int> real_atom_map;       ///< Map of atoms in the table to real atom indices in the
                                        ///<   parent topology (this allows the formal charge and
                                        ///<   bond order estimate to see through virtual sites)
  std::vector<int> real_bond_map;       ///< Map of bonds in the table to real bond indices in the
                                        ///<   parent topology
  std::vector<int> bond_i_atoms;        ///< First atoms within this object to which each bond
                                        ///<   pertains
  std::vector<int> bond_j_atoms;        ///< Second atoms within this object to which each bond
                                        ///<   pertains
  std::vector<int> atom_connect_bounds; ///< Bounds for each atom's entries in the array of bonded
                                        ///<   neighbors
  std::vector<int> atom_relevant_bonds; ///< Indices (in the system defined for this IndigoTable,
                                        ///<   not necessarily the same ordering as the original
                                        ///<   topology) for all bonds that affect a given atom
                                        ///<   center.  The atom_nb12_bounds array indexes the
                                        ///<   start and stop of any given atom center's list of
                                        ///<   bonds in this and the atom_nb12_partners array.
  std::vector<int> atom_bond_partners;  ///< Indices (in the system defined for this IndigoTable,
                                        ///<   not necessarily the same ordering as the original
                                        ///<   topology) for all of any given atom's bond partners
  std::vector<int> atom_scores;         ///< Scores for each atom's formal charge options
  std::vector<int> bond_scores;         ///< Scores for each bond's bond order options
  std::vector<int> valence_electrons;   ///< Targets for valence electron populations in each atom
                                        ///<   (must be met for the configuration to be viable)

  /// Given the above, atom centers can be constructed, one per atom, that combine the formal
  /// charge state and the bond orders of all bonds connecting the atom to any others.  The
  /// available states for these collective variables can thereby be enumerated, all within the
  /// object construction.
  std::vector<IndigoAtomCenter> atom_centers;

  /// Atom centers that have only one possible state can be ignored when seraching for the minimum
  /// energy structure.  Atom centers with more than one possible state but exist in isolation
  /// (surrounded by atom centers with only one possible state) will contribute to the minimum
  /// energy state by adopting their minimum energy configuration sbuject to overall charge
  /// constraints.  Clusters of contiguous atom centers with more than one possible state present
  /// the hardest challenge: how many distinct states might they adopt, and what net charges would
  /// those states contribute to the final result?  For this reason, one or more contiguous groups
  /// of atoms which all have multiple states will be grouped together into fragments and all
  /// (possible) states of each fragment will be evaluated for their combined energy
  std::vector<IndigoFragment> mutable_fragments;

  /// Formal charges on all atoms, in atomic units (+1 = the charge of a proton), without any
  /// adjustment for the formal charge scale used in scoring.
  std::vector<double> ground_state_formal_charges;

  /// Bond orders on all atoms in the ground states
  std::vector<double> ground_state_bond_orders;

  /// Bond orders on all atoms in the ground states
  std::vector<double> ground_state_free_electrons;

  /// Formal charges on all atoms in the single lowest-energy state of each molecule
  std::vector<double> zerok_formal_charges;

  /// Orders of all bonds in the single lowest-energy state of each molecule
  std::vector<double> zerok_bond_orders;

  /// Free electron content all atoms in the single lowest-energy state of each molecule
  std::vector<double> zerok_free_electrons;

  /// Pointer to the original topology
  const AtomGraph *ag_pointer;

  /// \brief Place an atom center in a particular state as part of the ground state
  ///
  /// \param ac             The atom center to work with
  /// \param state_index    State of the atom center that will contribute
  /// \param accumulation   Number of the iteration over which the average state is accumulated.
  ///                       If zero, the values of the atomic formal charge and associated bonds
  ///                       will be set to the respective state values times the Boltzmann weight.
  ///                       Otherwise the values times the probability will be added to the
  ///                       existing totals.
  /// \param probability    Boltzmann weight of this state among others that might contribute
  double addToGroundState(const IndigoAtomCenter &ac, int state_index,
                          std::vector<double> *acc_formal_charges,
                          std::vector<double> *acc_bond_orders,
                          std::vector<double> *acc_free_electrons, int accumulation = 0,
                          double probability = 1.0);
};

} // namespace chemistry
} // namespace stormm

#endif
