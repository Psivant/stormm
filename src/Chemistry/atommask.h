// -*-c++-*-
#ifndef STORMM_AMBMASK_H
#define STORMM_AMBMASK_H

#include <string>
#include <vector>
#include <climits>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"
#include "chemical_features.h"

namespace stormm {
namespace chemistry {

using constants::ExceptionResponse;
using parse::WildCardKind;
using topology::AtomGraph;
using trajectory::CoordinateFrame;
using trajectory::PhaseSpace;

/// \brief Enumerate the operators that can compose an atom mask.  This list is a union of
///        operators used by all types of masks.
enum class MaskOperator {
  NONE,            ///< No operation (this MaskComponent does not entail an operation)
  ATOM_LT,         ///< All atoms within some range of a selected group
  ATOM_LE,         ///< All atoms lying at or within some range of a selected group
  ATOM_GT,         ///< All atoms outside of some range of a selected group
  ATOM_GE,         ///< All atoms lying at or greater than a given range from a selected group
  RESIDUE_LT,      ///< All residues such that all atoms lie within range of a selected group
  RESIDUE_GT,      ///< All residues with no atoms in some range of a selected group
  RESIDUE_LE,      ///< Residues with at least one atom lying in range of a selected group
  RESIDUE_GE,      ///< Residues with at least one atom lying beyond some range of a selected group
  AND,             ///< One selection AND another
  OR,              ///< One selection OR another (or both)
  NOT              ///< Not atoms in some selection
};

/// \brief Different mask components will target different pieces of topology: atoms (by name or
///        index), atomic elements, atom types, or residues (by name or index).
enum class SelectionItemKind {
  ATOM_NAME,      ///< Any atoms in the topology, filtered by name
  ATOM_NUMBER,    ///< Any atoms in the topology, filtered by number in the sequence, starting at 1
  ELEMENT_NAME,   ///< Elements, filtered by name
  ELEMENT_NUMBER, ///< Elements, filtered by atomic number (0 denotes virtual sites)
  ATOM_TYPE,      ///< Atom types, the abstract and artificial molecular mechanics classifications
                  ///<   that make force fields what they are.  Filtered by name.
  RESIDUE_NAME,   ///< Residues, the repeating units of biological heteropolymers or also
                  ///<   individual molecules like water or ions, filtered by name in the topology.
  RESIDUE_NUMBER, ///< Residues, filtered by number in the topology, starting at 1.
#if 0
  RING_SIZE,      ///< Filter by the sizes of rings in which the atom participates.
#endif
  NONE            ///< No filtering criteria (default)
};

/// \brief The quantum of an atom mask: a specific name, or a numerical range, of atoms or
///        residues.  A vector of these composes one part of an atom mask.  Multiple masks combine
///        according to additional operators (not '!', and '&', as well as or '|') in order to
///        compile a single ambmask.  This struct is unguarded.
class SelectionItem {
public:

  /// \brief Initialization for this struct makes it simpler to create an array of such objects
  ///        and be assured that they will not direct mask formulation until filled.
  ///
  /// Overloaded:
  ///   - Takes a kind and a name
  ///   - Takes a kind and up to two numbers (the second defaults to the C-standard minimum int
  ///     value, to ensure that it is less than or equal to the first number which will indicate
  ///     that it should then be set to the last atom in the system)
  /// \{
  SelectionItem(const SelectionItemKind kind_in, const char4 name_in,
                const std::vector<WildCardKind> &wildcards_in);
  SelectionItem(const SelectionItemKind kind_in, const int begin_in, const int end_in = INT_MIN);
  /// \}
  
  SelectionItemKind kind;               ///< Indicates the type of selection, i.e. atom name
  char4 name;                           ///< The (perhaps codified) name to look for, e.g. of an
                                        ///<   atom or residue
  std::vector<WildCardKind> wildcards;  ///< Indicators of wildcards in the name
  int scan_begin;                       ///< Point at which to begin the scan
  int scan_end;                         ///< Point at which to end the scan
};

/// \brief Enumerate the types of components to a mask
enum class MaskComponentKind {
  OPERATOR,  ///< This mask component is an operator that joins or modifies MASK components on
             ///<   either side of it
  MASK       ///< This mask component is a bitmask spanning all atoms in the system
};

/// \brief A mask is a series of operators and more primitive masks.  This struct comprises both
///        things, with an indicator as to what it contains.  A complete mask for a given scope is
///        composed by stepping through a list of these objects.
class MaskComponent {
public:

  /// \brief Constructor takes either an operator designation or a series of atoms or residues
  ///
  /// Overloaded:
  ///   - Build a basic "&", "|", or "!" operator component
  ///   - Build a ranged operator component
  ///   - Build a primitive mask component from a list of named or numbered atoms or residues
  ///   - Take a primitive mask component output by a recursive call to parseMask and store it
  ///
  /// \param op_in              Operator designation
  /// \param parts_in           A series of atom or residue names and indices to find
  /// \param ag                 System topology containing atom and residue names and indices
  /// \param primitive_mask_in  Bit mask for all atoms, as output by parseMask (see below)
  /// \{
  MaskComponent(MaskOperator op_in, const std::string &basis_in, int start_idx, int end_idx);
  MaskComponent(MaskOperator op_in, double range, const std::string &basis_in, int start_idx,
                int end_idx);
  MaskComponent(const std::vector<SelectionItem> &parts_in, const AtomGraph *ag,
                const ChemicalFeatures &chemfe, const std::string &basis_in, int start_idx,
                int end_idx);
  MaskComponent(const std::vector<uint> &primitive_mask_in, const std::string &basis_in,
                int start_idx, int end_idx);
  /// \}

  /// \brief Get the kind of component
  MaskComponentKind getKind() const;

  /// \brief Get the operator, if relevant
  MaskOperator getOperator() const;

  /// \brief Get the bitmask, if relevant
  std::vector<uint> getMask() const;

  /// \brief Get the textual basis of this mask component
  std::string getText() const;

  /// \brief Get the indices from within the original mask string where this component's textual
  ///         basis is to be found.
  int2 getTextLimits() const;

  /// \brief Apply the NOT operator to a primitive mask
  void applyNotOperator();

  /// \brief Apply the AND operator to a primitive mask, given another primitive mask
  ///
  /// \param other  Pre-determined primitive mask of atoms
  void applyAndOperator(const std::vector<uint> &other);

  /// \brief Apply the OR operator to a primitive mask, given another primitive mask
  ///
  /// \param other  Pre-determined primitive mask of atoms
  void applyOrOperator(const std::vector<uint> &other);

  /// \brief Apply one of the ranged operators.  This will produce a new MASK-kind MaskComponent.
  ///        Ranged operators are expected to occur last in their scopes.  Once the scope is
  ///        broken into MaskComponents of operators and other masks, the list is processed by
  ///        applying operators to eliminate MaskComponents until there is a single MASK-kind
  ///        MaskComponent left.  With range operators, the list whittles down to a MASK-kind
  ///        followed by an OPERATOR MaskComponent.  Both of these will be eliminated after the
  ///        output of this function is tacked onto the end of the list, where it will then become
  ///        the only element.
  ///
  /// \param other  Pre-determined primitive mask of atoms
  /// \param ag     Pointer to the original topology
  /// \param cfr    Coordinates of the system
  MaskComponent applyRangeOperator(const std::vector<uint> &other, const AtomGraph *ag,
                                   const CoordinateFrameReader &cfr);

private:
  MaskComponentKind kind;
  MaskOperator op;
  double range;
  std::vector<uint> primitive_mask;
  std::string text_basis;
  int2 text_limits;
};

/// \brief An atom mask may comprise nearly the entire system or so many atoms dispersed throughout
///        the list that the best approach is test every atom to see whether it is part of the
///        mask.  This is not common, however.  In most cases, the intensive search can be narrowed
///        to one part of the list, a handful of contiguous stretches of masked atoms, or even a
///        single stretch of masked atoms.
enum class MaskTraversalMode {
  COMPLETE,    ///< Scan all atoms in the entire list, testing each one
  PARTIAL,     ///< Scan all atoms in part of the list, testing each one
  SEGMENTED    ///< Scan between a series of lower and upper limits.  Every atom between these
               ///<   limits is guaranteed to be part of the mask.
};

/// \brief The input mode for an AtomMask, indicating how it was constructed.
enum class MaskInputMode {
  AMBMASK,  ///< Interpreting Amber's ambmask format 
  VMD       ///< Interpreting VMD's atom selection format
};
  
/// \brief An atom selection within a system.  Internally, this stores a bitmask, one bit for every
///        particle, to denote whether it is in the mask.  The object also stores the original
///        string upon which the mask was constructed and a description of the mask, as it is a
///        significant computation to build.
class AtomMask {
public:

  /// \brief Constructor takes a mask string, topology, and coordinates
  ///
  /// Overloaded:
  ///   - Blank constructor (build later from a provided topology and optional coordinates)
  ///   - Require no coordinates (only valid if the mask has no range operators)
  ///   - Take coordinates as a PhaseSpace object
  ///   - Take coordinates as a CoordinateFrame object
  ///   - Take coordinates as three std::vectors for Cartesian X, Y, and Z positions, with or
  ///     without unit cell transforms (if the transform vectors do not exist, the system is
  ///     assumed to have isolated boundary conditions)
  ///   - Take coordinates as three C-style arrays, with unit cell transforms (if the transform
  ///     arrays are nullptr, the system is assumed to exist in isolated boundary conditions)
  ///
  /// \param input_text_in  The mask string
  /// \param ag_in          System topology
  /// \param chemfe         Chemical features computed for the topology
  /// \param crd            Coordinates obtained from a CoordinateFrame object
  /// \param xcrd           Cartesian X coordinates of all atoms
  /// \param ycrd           Cartesian Y coordinates of all atoms
  /// \param zcrd           Cartesian Z coordinates of all atoms
  /// \param description    [Optional] description of the mask and its purpose
  /// \{
  AtomMask(const AtomGraph *ag_in = nullptr);

  AtomMask(const std::string &input_text_in, const AtomGraph *ag_in,
           const ChemicalFeatures &chemfe, const CoordinateFrameReader &cfr,
           MaskInputMode mode = MaskInputMode::AMBMASK,
           const std::string &description_in = std::string("No description provided"));

  AtomMask(const std::string &input_text_in, const AtomGraph *ag_in,
           const ChemicalFeatures &chemfe, const CoordinateFrame &cf,
           MaskInputMode mode = MaskInputMode::AMBMASK,
           const std::string &description_in = std::string("No description provided"));

  AtomMask(const std::string &input_text_in, const AtomGraph *ag_in,
           const ChemicalFeatures &chemfe, const PhaseSpace &ps,
           MaskInputMode mode = MaskInputMode::AMBMASK,
           const std::string &description_in = std::string("No description provided"));
  /// \}

  /// The standard copy and move constructors will be effective for this object, which has no
  /// pointers to repair.  Because there is a const member, the copy and move assignment operators
  /// will be implicitly deleted.
  /// \{
  AtomMask(const AtomMask &original) = default;
  AtomMask(AtomMask &&original) = default;
  /// \}
  
  /// \brief Get the recommendation for scanning this mask
  MaskTraversalMode getRecommendation() const;

  /// \brief Get the raw mask as a vector of unsigned long integers (bit strings storing whether
  ///        each atom is masked at one atom per bit)
  const std::vector<uint>& getRawMask() const;

  /// \brief Get a count of the number of masked atoms
  int getMaskedAtomCount() const;

  /// \brief Get the mask as a vector of boolean values (not just a bit mask as the object
  ///        stores privately)
  std::vector<bool> getMask() const;

  // \brief Get a list of the masked atoms
  std::vector<int> getMaskedAtomList() const;

  /// \brief Get the segments of atoms over which the mask applies (only filled out if the
  ///        recommended scanning method is by segment)
  std::vector<int2> getSegments() const;

  /// \brief Get the input text used to generate this atom mask
  std::string getInputText() const;

  /// \brief Get the input text mode, e.g. AMBMASK to indicate that an Amber-style atom mask
  ///        string was used to make this mask.
  MaskInputMode getInputKind() const;
  
  /// \brief Get the description of the mask, if provided
  std::string getDescription() const;

  /// \brief Get a pointer to the topology that this mask describes
  const AtomGraph* getTopologyPointer() const;

  /// \brief Answer whether an atom is in the mask.
  ///
  /// Overloaded:
  ///   - Provide the atom name only (any instance of the atom name found in the mask will suffice
  ///     to say that an atom with the specified name is present in the mask)
  ///   - Provide the atom residue name and atom name (any instance of an atom conforming to this
  ///     pair of names will result in the function returning TRUE)
  ///   - Provide the atom index number (indexing starting at zero)
  ///
  /// \param atom_name   Name of the atom of interest
  /// \param res_name    Name of the residue of interest
  /// \param atom_index  Topological index number of the atom of interest
  /// \{
  bool isAtomInMask(const char4 atom_name) const;
  bool isAtomInMask(const char4 res_name, const char4 atom_name) const;
  bool isAtomInMask(int atom_index) const;
  /// \}
  
  /// \brief Add atoms to the atom mask.
  ///
  /// Overloaded:
  ///   - Add a series of atom indices, with indexing beginning at 0.
  ///   - Add a series of atom names.  All examples of any supplied name will be matched and added
  ///     to the mask (an O(NM) operation for N atom names and a topology of M atoms)
  ///   - Add the atoms of another mask (the topologies will be checked for a reasonable degree of
  ///     congruity)
  ///   - Add atoms implied by another mask string.
  ///   - Supply the system's pre-computed chemical features, or have it constructed temporarily.
  ///
  /// \param new_indices  The atom indices to add to the current mask
  /// \param new_names    The names of atoms to add to the current mask
  /// \param new_mask     Another atom mask (or mask string) to add to the current mask
  /// \param policy       The policy to take in the event that an error is encountered
  /// \{
  void addAtoms(const std::vector<int> &new_indices,
                ExceptionResponse policy = ExceptionResponse::DIE);

  void addAtoms(const std::vector<char4> &new_names);

  void addAtoms(const AtomMask &new_mask, const CoordinateFrame &cf,
                const ChemicalFeatures &chemfe);

  void addAtoms(const std::string &new_mask, const CoordinateFrame &cf,
		const ChemicalFeatures &chemfe);

  void addAtoms(const AtomMask &new_mask, const CoordinateFrame &cf);
  
  void addAtoms(const std::string &new_mask, const CoordinateFrame &cf);
  /// \}
  
private:
  MaskTraversalMode recommended_scan;      ///< Indicates the most efficient way to read this mask
  MaskInputMode style;                     ///< Mask syntax (default AMBMASK, see above)
  int masked_atom_count;                   ///< Number of atoms included in the mask
  std::vector<uint> raw_mask;              ///< Compact bitmask for the topology
  std::vector<int2> segments;              ///< Segments over which the mask has nonzero elements
  std::string input_text;                  ///< Input text for the mask-generating string, i.e.
                                           ///<   the ambmask string ':5-8 & @N,CA,C,O'
  std::string description;                 ///< Description of this AtomMask object, indicating its
                                           ///<   place in a program
  const AtomGraph *ag_pointer;             ///< Pointer to the topology this mask is based on

  /// \brief Function to extract a ranged value when prompted by range operators in an atom mask
  ///
  /// \param start_char  The first character after a detected range operator
  /// \param final_char  The last character in the scope
  /// \param position    Used to return the position of the last character used in finding the
  ///                    range value of interest
  double extractRangeValue(int start_char, int final_char, int *position);

  /// \brief Function to evaluate a list of atom or residue inclusions within a more complex
  ///        atom mask.  Returns a primitive bitmask for all atoms, which can become the raw_mask
  ///        of the AtomMask object itself, but it subject to any other operations "&", "|", and
  ///        "!" in the scope and then any other operations acting on masks evaluated in other
  ///        scopes.  
  ///
  /// \param inclusion_list  The list of inclusions, i.e. ":1-5,TYR,10-14,ARG,SER" or "@CA,CB"
  std::vector<SelectionItem> evaluateInclusions(const std::string &inclusion_list);
  
  /// \brief Recursive function for parsing an Amber ambmask at multiple levels
  ///
  /// \param scope_levels  The scope level of each character in the input mask string
  /// \param position      Used to return the position of the last character used in finding the
  ///                      range value of interest
  std::vector<uint> parseMask(const std::vector<int> &scope_levels, int *position,
                              const CoordinateFrameReader &cfr, const ChemicalFeatures &chemfe);
};

} // namespace chemistry
} // namespace stormm

#endif
