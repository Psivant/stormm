// -*-c++-*-
#ifndef STORMM_MDL_FILE_H
#define STORMM_MDL_FILE_H

#include <fstream>
#include <string>
#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Chemistry/znumber.h"
#include "Chemistry/chemical_features.h"
#include "Constants/behavior.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "FileManagement/file_util.h"
#include "Parsing/ascii_numbers.h"
#include "Parsing/parse.h"
#include "Potential/energy_enumerators.h"
#include "Restraints/restraint_apparatus.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"
#include "Trajectory/trajectory_enumerators.h"
#include "molecule_file_io.h"
#include "molecule_format_enumerators.h"
#include "mdlmol_atomlist.h"
#include "mdlmol_bond.h"
#include "mdlmol_dataitem.h"
#include "mdlmol_property.h"

namespace stormm {
namespace structure {

using card::HybridTargetLevel;
using chemistry::symbolToZNumber;
using chemistry::ChemicalFeatures;
using constants::CartesianDimension;
using constants::CaseSensitivity;
using constants::ExceptionResponse;
using data_types::isFloatingPointScalarType;
using diskutil::PrintSituation;
using energy::StateVariable;
using parse::TextFile;
using restraints::RestraintApparatus;
using synthesis::Condensate;
using synthesis::CondensateReader;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisReader;
using topology::AtomGraph;
using topology::ChemicalDetailsKit;
using topology::NonbondedKit;
using trajectory::CoordinateCycle;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateFrameWriter;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesReader;
using trajectory::CoordinateSeriesWriter;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceReader;
using trajectory::PhaseSpaceWriter;

/// \brief Default settings for the MDL MOL object atom initializations
/// \{
constexpr char4 default_mdl_atomic_symbol = { ' ', ' ', ' ', ' ' };
constexpr int default_mdl_atomic_number = -1;
constexpr int default_mdl_formal_charge = 0;
constexpr RadicalState default_mdl_radical_state = RadicalState::NONE;
constexpr int default_mdl_isotopic_shift = 0;
constexpr MolObjAtomStereo default_mdl_stereo_parity = MolObjAtomStereo::NOT_STEREO;
constexpr HydrogenAssignment default_hydrogenation = HydrogenAssignment::DO_NOT_HYDROGENATE;
constexpr int default_mdl_implicit_hydrogen = 0;
constexpr int default_mdl_valence_connections = 0;
constexpr int default_mdl_map_count = 0;
constexpr bool default_mdl_stereo_considerations = false;
constexpr bool default_mdl_exact_change = false;
constexpr StereoRetention default_mdl_stereo_retention = StereoRetention::NOT_APPLIED;
/// \}

/// \brief An SText group from an MDL MOL (.mol) or SDF file.  This information is used only by
///        old ISIS / Desktop programs and is otherwise deprecated.
class MolObjSTextGroup {
public:
private:
};
  
/// \brief A molecular three-dimensional feature.  This special class of MOL object properies has
///        its own data lines.
  
/// \brief A molecule read from an MDL .mol file, or one of many read from a concatenated SDF file.
///        Many of the enumerators above are translated according to member functions of this
///        object based on the documentation in the ctfileformats.pdf file in this library's
///        directory, also available at:
///
///        http://help.accelrysonline.com/ (...)
///          ulm/onelab/1.0/content/ulm_pdfs/direct/reference/ctfileformats2016.pdf
class MdlMol {
public:

  /// \brief Constructors for the MDL molecule format object (known as MolObj in RDKit)
  ///
  /// Overloaded:
  ///   - Basic constructor for creating a blank MdlMol, referenced by all other constructors'
  ///     initializer lists
  ///   - Constructors based on a file name (string or const char* array)
  ///   - Constructor based on a TextFile object from an SDF container previously committed to RAM
  ///
  /// \param filename        Name of the file to read, containing one MDL file
  /// \param tf              Text file data previously read into memory (this is the way to handle
  ///                        an SDF file)
  /// \param line_start      The irst line at which to begin reading the TextFile object
  /// \param line_end        Last relevant line of the TextFile object at which reading will stop
  /// \param capitalization  Indicate whether atomic symbol capitalization can be ignored when
  ///                        assigning elements to each atom
  /// \param policy_in       Course of action to take if errors are encountered when inferring
  ///                        atomic elements or other critical features of the structure
  /// \param dimod_policy    Preferred course of action if the names of some data items are found
  ///                        to break the strict Biovia standard: should they be repaired
  ///                        (modified) or left as they are?
  /// \param dimod_notify    Indicate whether to alert the user if modifications are made to data
  ///                        items.
  /// \{
  MdlMol(ExceptionResponse policy_in = ExceptionResponse::WARN);

  MdlMol(const std::string &filename, ExceptionResponse policy_in = ExceptionResponse::WARN,
         ModificationPolicy dimod_policy = ModificationPolicy::DO_NOT_MODIFY,
         ExceptionResponse dimod_notify = ExceptionResponse::WARN);

  MdlMol(const char* filename, ExceptionResponse policy_in = ExceptionResponse::WARN,
         ModificationPolicy dimod_policy = ModificationPolicy::DO_NOT_MODIFY,
         ExceptionResponse dimod_notify = ExceptionResponse::WARN);

  MdlMol(const TextFile &tf, int line_start = 0, int line_end = -1,
         CaseSensitivity capitalization = CaseSensitivity::YES,
         ExceptionResponse policy_in = ExceptionResponse::WARN,
         ModificationPolicy dimod_policy = ModificationPolicy::DO_NOT_MODIFY,
         ExceptionResponse dimod_notify = ExceptionResponse::WARN);

  template <typename T>
  MdlMol(const ChemicalFeatures *chemfe, const T* xcrd, const T* ycrd, const T* zcrd,
         double inv_scale, int molecule_index = 0);

  template <typename T>
  MdlMol(const ChemicalFeatures *chemfe, const T* xcrd, const T* ycrd, const T* zcrd,
         int molecule_index = 0);

  MdlMol(const ChemicalFeatures *chemfe, const CoordinateFrameReader cfr, int molecule_index = 0);
  
  MdlMol(const ChemicalFeatures *chemfe, const CoordinateFrame *cf, int molecule_index = 0);

  MdlMol(const ChemicalFeatures &chemfe, const CoordinateFrame &cf, int molecule_index = 0);

  MdlMol(const ChemicalFeatures *chemfe, const PhaseSpaceReader psr, int molecule_index = 0);

  MdlMol(const ChemicalFeatures *chemfe, const PhaseSpace *ps, int molecule_index = 0);

  MdlMol(const ChemicalFeatures &chemfe, const PhaseSpace &ps, int molecule_index = 0);

  template <typename T>
  MdlMol(const ChemicalFeatures *chemfe, const CoordinateSeriesReader<T> csr, int frame_index,
         int molecule_index = 0);

  template <typename T>
  MdlMol(const ChemicalFeatures *chemfe, const CoordinateSeries<T> *cs, int frame_index,
         int molecule_index = 0);

  template <typename T>
  MdlMol(const ChemicalFeatures &chemfe, const CoordinateSeries<T> &cs, int frame_index,
         int molecule_index = 0);

  MdlMol(const ChemicalFeatures *chemfe, const PsSynthesisReader poly_psr, int system_index,
         int molecule_index = 0);

  MdlMol(const ChemicalFeatures *chemfe, const PhaseSpaceSynthesis *poly_ps, int system_index,
         int molecule_index = 0);

  MdlMol(const ChemicalFeatures &chemfe, const PhaseSpaceSynthesis &poly_ps, int system_index,
         int molecule_index = 0);

  MdlMol(const ChemicalFeatures *chemfe, const CondensateReader *cdnsr, int system_index,
         int molecule_index = 0);

  MdlMol(const ChemicalFeatures *chemfe, const Condensate *cdns, int system_index,
         int molecule_index = 0);

  MdlMol(const ChemicalFeatures &chemfe, const Condensate &cdns, int system_index,
         int molecule_index = 0);
  /// \}

  /// \brief Default copy and move constructors, as well as assignment operators, are appropriate
  ///        for this object which consists entirely of scalar data types, Standard Template
  ///        Library objects, and no const members.
  /// \{
  MdlMol(const MdlMol &original) = default;
  MdlMol(MdlMol &&original) = default;
  MdlMol& operator=(const MdlMol &original) = default;
  MdlMol& operator=(MdlMol &&original) = default;
  /// \}

  /// \brief Get the title of the MDL MOL entry (for its first three lines)
  const std::string& getTitle() const;
  
  /// \brief Get the system's atom count.
  int getAtomCount() const;
  
  /// \brief Get the number of bonds in the system.
  int getBondCount() const;

  /// \brief Get the number of properties found in the MDL MOL entry.
  int getPropertiesCount() const;

  /// \brief Get the number of data items (XML-like note found after the MDL section of the entry
  ///        in a BIOVIA Structure Data (SD)-format file).
  int getDataItemCount() const;

  /// \brief Get the { X, Y, Z } coordinate tuple for a particular atom, or for all atoms.
  ///
  /// Overloaded:
  ///   - Get a unique, modifiable tuple for a single atom.
  ///   - Get a const reference to the array of tuples for coordinates of all atoms.
  ///   - Get a modifiable vector of all Cartesian X, Y, or Z coordinates.
  ///
  /// \param index  Index of the atom of interest
  /// \param dim    The Cartesian dimension of interest
  /// \{
  double3 getCoordinates(int index) const;
  const std::vector<double3>& getCoordinates() const;
  std::vector<double> getCoordinates(CartesianDimension dim) const;
  /// \}

  /// \brief Get a const pointer to the object, useful if a pointer is needed and the object is
  ///        available by const reference.
  const MdlMol* getSelfPointer() const;
  
  /// \brief Export the coordinates as a PhaseSpace object suitable for basic molecular mechanics
  ///        force computations.
  PhaseSpace exportPhaseSpace() const ;

  /// \brief Export the coordinates as a stripped-down CoordinateFrame object suitable for
  ///        molecular mechanics energy computations.
  CoordinateFrame exportCoordinateFrame() const;

  /// \brief Get the atomic symbol for a particular atom.
  ///
  /// \param index  Index of the atom of interest
  char4 getAtomSymbol(int index) const;

  /// \brief Get a const reference to the atomic symbols for all atoms.
  const std::vector<char4>& getAtomSymbols() const;
  
  /// \brief Get the atomic number of a particular atom.
  ///
  /// \param index  Index number of the atom of interest
  int getAtomicNumber(int index) const;

  /// \brief Get a const reference to the vector of all atomic numbers in the system.
  const std::vector<int>& getAtomicNumbers() const;

  /// \brief Get the formal charge on a particular atom.
  ///
  /// \param index  Index number of the atom of interest
  int getFormalCharge(int index) const;

  /// \brief Get a const reference to the vector of all formal charges.
  const std::vector<int>& getFormalCharges() const;
  
  /// \brief Impart a new set of coordinates to the atoms based on one of STORMM's major coordinate
  ///        classes.  This will permanently modify the underlying coordinates of the object and
  ///        also any properties dependent of the structure.  If providing a standalone coordinate
  ///        object, the atom count of the incoming coordinate set will be checked.  In all cases,
  ///        the order of atoms in the incoming coordinate set is expect to match this object's
  ///        own arrangement.  If presenting an abstract of a PhaseSpace object, the coordinates of
  ///        the WHITE stage in the time cycle will be taken.  When presenting the original
  ///        objects, the developer may select the stage in the time cycle.
  ///
  /// Overloaded:
  ///   - Present three C-style arrays for X, Y, and Z coordinates with trusted lengths
  ///   - Present a PhaseSpace object, or abstract thereof
  ///   - Present a CoordinateFrame object, or abstract thereof
  ///   - Present a CoordinateSeries object, or abstract thereof, with a frame index number
  ///
  /// \param xcrd            Cartesian X coordinates of all particles, trusted to describe a number
  ///                        of atoms equal to that found in the MdlMol and in the same order
  ///                        unless supplemented with a ChemicalDetailsKit and molecule index,
  ///                        which can indicate the proper subset and indices of the atomic
  ///                        coordinates to use
  /// \param ycrd            Cartesian Y coordinates of all particles
  /// \param zcrd            Cartesian Z coordinates of all particles
  /// \param xcrd_ovrf       Overflow bits for Cartesian X coordinates of all particles
  /// \param ycrd_ovrf       Overflow bits for Cartesian Y coordinates of all particles
  /// \param zcrd_ovrf       Overflow bits for Cartesian Z coordinates of all particles
  /// \param scale_factor    The scaling factor by which to multiply coordinates in order to take
  ///                        them into units of Angstroms (in general, this is the inverse scaling
  ///                        factor found in the abstracts of fixed-precision coordinate objects)
  /// \param cdk             Contains indices of the atoms of interest within the (likely larger)
  ///                        arrays provided
  /// \param molecule_index  Index of the molecule of interest within some larger topology
  /// \param ps              Coordinates to transfer
  /// \param cf              Coordinates to transfer
  /// \param cs              A series of frames, one containing the coordinates to transfer
  /// \param poly_ps         A complex collection of systems containing coordinates to transfer
  /// \param frame_index     Frame within the coordinate series to transfer
  /// \param system_index    System within the phase space synthesis to transfer
  /// \param tier            Indicate whether to obtain coordinates from either the CPU host or
  ///                        GPU device memory
  /// \{
  template <typename T> void impartCoordinates(const T* xcrd, const T* ycrd, const T* zcrd,
                                               double scale_factor);

  template <typename T> void impartCoordinates(const T* xcrd, const T* ycrd, const T* zcrd,
                                               double scale_factor, const ChemicalDetailsKit &cdk,
                                               int molecule_index);
  
  void impartCoordinates(const llint* xcrd, const llint* ycrd, const llint* zcrd,
                         const int* xcrd_ovrf, const int* ycrd_ovrf, const int* zcrd_ovrf,
                         double scale_factor);

  void impartCoordinates(const llint* xcrd, const llint* ycrd, const llint* zcrd,
                         const int* xcrd_ovrf, const int* ycrd_ovrf, const int* zcrd_ovrf,
                         double scale_factor, const ChemicalDetailsKit &cdk, int molecule_index);
  
  void impartCoordinates(const PhaseSpaceReader &psr);

  void impartCoordinates(const PhaseSpaceWriter &psw);

  void impartCoordinates(const PhaseSpace *ps, CoordinateCycle orientation,
                         HybridTargetLevel tier = HybridTargetLevel::HOST);

  void impartCoordinates(const PhaseSpace *ps,
                         HybridTargetLevel tier = HybridTargetLevel::HOST);

  void impartCoordinates(const PhaseSpace &ps, CoordinateCycle orientation,
                         HybridTargetLevel tier = HybridTargetLevel::HOST);

  void impartCoordinates(const PhaseSpace &ps,
                         HybridTargetLevel tier = HybridTargetLevel::HOST);

  void impartCoordinates(const CoordinateFrameReader &cfr);

  void impartCoordinates(const CoordinateFrameWriter &cfw);

  void impartCoordinates(const CoordinateFrame *cf,
                         HybridTargetLevel tier = HybridTargetLevel::HOST);

  void impartCoordinates(const CoordinateFrame &cf,
                         HybridTargetLevel tier = HybridTargetLevel::HOST);

  template <typename T>
  void impartCoordinates(const CoordinateSeriesReader<T> &cs, int frame_index,
                         HybridTargetLevel tier = HybridTargetLevel::HOST);

  template <typename T>
  void impartCoordinates(const CoordinateSeriesWriter<T> &cs, int frame_index,
                         HybridTargetLevel tier = HybridTargetLevel::HOST);

  template <typename T>
  void impartCoordinates(const CoordinateSeries<T> *cs, int frame_index,
                         HybridTargetLevel tier = HybridTargetLevel::HOST);

  template <typename T>
  void impartCoordinates(const CoordinateSeries<T> &cs, int frame_index,
                         HybridTargetLevel tier = HybridTargetLevel::HOST);

  void impartCoordinates(const PhaseSpaceSynthesis *poly_ps, int system_index,
                         HybridTargetLevel tier = HybridTargetLevel::HOST);

  void impartCoordinates(const PhaseSpaceSynthesis &poly_ps, int system_index,
                         HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Add a property to the object.  Properties must follow a list of pre-approved BIOVIA
  ///        codes, but provide a means of extending the V2000 format to include information found
  ///        in the V3000 format.
  ///
  /// Overloaded:
  ///   - Provide the property by pointer
  ///   - Provide the property by reference
  ///
  /// \param mmprop  The property to add
  /// /{
  void addProperty(const MdlMolProperty *mmprop);
  void addProperty(const MdlMolProperty &mmprop);
  /// \}

  /// \brief Add a data item to the object.  These data items follow the classifications set forth
  ///        in the DataRequestKind enumerator (see molecule_format_enumerators.h).  Each overload
  ///        is designed to serve a particular case.
  ///
  /// Overloaded:
  ///   - Add a molecular mechanics energy term, with a value supplied by evaluating the structure
  ///     coordinates with the topology and restraint apparatus.
  ///   - Add a summary of the valence interaction influences (including restraints) on each atom
  ///     in a mask, with forces and energies supplied by evaluating the topology and restraints.
  ///   - Track down a specific restraint or topological energy term based on a series of atom
  ///     types and display its settings.
  ///   - Add a user-defined string as a message in the SD output.
  ///
  /// \param ask  The request for an item to impart ot the SD file
  /// \param ag   The topology guiding the motion of the system
  /// \param ra   The system of restraints guiding the structure alongside the topology
  /// \{
  void addDataItem(const MdlMolDataRequest &ask, const AtomGraph &ag,
                   const RestraintApparatus &ra);
  /// \}

  /// \brief Get the classification of a data item from within the MDL MOL entry.
  ///
  /// \param item_index  Position of the data item of interest in the data_items array
  MdlMolDataItemKind getDataItemKind(int item_index) const;

  /// \brief Get the state variable that one of the data items is tracking.
  StateVariable getTrackedState(int item_index) const;

  /// \brief Add a line of text to a specific data item within the MDL MOL entry.  This text may
  ///        contain a pre-formatted string displaying some energy quantity, or more complex
  ///        syntax.
  ///
  /// \param text        The text to append
  /// \param item_index  Index of the data item to update
  void addLineToDataItem(const std::string &text, int item_index);

  /// \brief Get the name of a data item based on its number
  ///
  /// \param item_index  Index of the data item of interest
  const std::string& getDataItemName(int item_index) const;

  /// \brief Get the name of a data item based on its number, fit for output in the proper Biovia
  ///        SD file format standard.
  ///
  /// \param item_index  Index of the data item of interest
  const std::string& getDataItemOutputName(int item_index) const;

  /// \brief Get the index of a data item based on a name.
  ///
  /// \param item_name   Name of the data item of interest
  int getDataItemIndex(const std::string &item_name) const;
  
  /// \brief Get the lines from a data item as an array of strings.
  ///
  /// Overloaded:
  ///   - Accept the index of the data item and return a const reference to its content (this is
  ///     possible because a data item with a valid index will have content to reference)
  ///   - Accept the name of the data item and return a copy of the content, if any data item is
  ///     found (if the search fails, a const reference to a temporary would be dangerous)
  ///
  /// \param item_index  Index of the data item of interest
  /// \param item_name   Name of the data item of interest
  /// \{
  const std::vector<std::string>& getDataItemContent(int item_index) const;
  std::vector<std::string> getDataItemContent(const std::string &item_name) const;
  /// \}

  /// \brief Write a set of molecular coordinates, bonds, and their annotations in MDL MOL format.
  ///        Apply all properties already stored in the object, such that the result is not an
  ///        exact reversal of the operations for reading such a file but correctly conveys the
  ///        relevant information to any other SDF / MDL MOL reader.  The data items are not
  ///        part of this output--for that section of the SD file format, use writeDataItems()
  ///        below.  The intention is that different coordinates can then be edited for a new
  ///        structure while the (smaller) data section is re-written separately for a more
  ///        efficient file I/O process.
  ///
  /// Overloaded:
  ///   - Write to an output file
  ///   - Write to a string
  ///
  /// \param foutp        Output file stream
  /// \param fname        Name of the output file
  /// \param vformat      Format version of the MDL MOL standard to use.  If provided as either
  ///                     V2000 or V3000, this will override information stored in the object.
  /// \param expectation  Anticipated (or required) condition of the output file that is to be
  ///                     opened, if only its name is provided
  /// \{
  void writeMdl(std::ofstream *foutp, MdlMolVersion vformat) const;

  void writeMdl(const std::string &fname, MdlMolVersion vformat, PrintSituation expectation) const;

  std::string writeMdl(MdlMolVersion vformat = MdlMolVersion::V3000) const;
  /// \}

  /// \brief Write the non-MDL components of an SD file, including all applicable data items, to a
  ///        string or file.  This is to be used in concert with the writeMdl() member function
  ///        above.
  ///
  /// Overloaded:
  ///   - Write to an output file
  ///   - Write to a string
  ///
  /// \param foutp        Output file stream
  /// \param fname        Name of the output file
  /// \param expectation  Anticipated (or required) condition of the output file that is to be
  ///                     opened, if only its name is provided
  /// \param mol_index    Index of the molecule within the SD file (a detail which may be added to
  ///                     each data item's header line)
  /// \{
  void writeDataItems(std::ofstream *foutp, int mol_index = 0) const;

  void writeDataItems(const std::string &fname, PrintSituation expectation,
                      int mol_index = 0) const;

  std::string writeDataItems(int mol_index = 0) const;
  /// \}
  
private:

  // Items describing quantities of information (most of them from the counts line)
  ExceptionResponse policy;   ///< Action to take if errors in an input file are encountered
  MdlMolVersion version_no;   ///< The format in which this entry was read (does not necessarily
                              ///<   dictate the version in which it will be written)
  int atom_count;             ///< The number of atoms in the molecule
  int bond_count;             ///< Number of bonds of all types between atoms in the system
  int list_count;             ///< The number of atom (element) lists
  int stext_entry_count;      ///< The number of S-text entries
  int properties_count;       ///< The number of additional properties (the default, and expected,
                              ///<   value is 999)
  int sgroup_count;           ///< The number of S-groups
  int constraint_count;       ///< The number of three-dimensional constraints
  MolObjChirality chirality;  ///< The molecule's chirality (assumes only one significant center)
  int registry_number;        ///< The molecule registry number
  int data_item_count;        ///< The number of data items

  // Atomic properties
  bool property_formal_charges;             ///< Flag to indicate that the formal charge properties
                                            ///<   were set by properties.  This will cause any
                                            ///<   printing in V2000 format to list a zero for the
                                            ///<   formal charge property of every atom in the
                                            ///<   atom block.  The actual value will be given in
                                            ///<   one of the MDL MOL properties.
  bool property_radicals;                   ///< Flag to indicate that the radical character of
                                            ///<   one or more atoms was set by a property,
                                            ///<   removing such information from the atom block
                                            ///<   of a V2000 printout.
  bool property_isotopes;                   ///< Flag to indicate that the isotopic shift of one or
                                            ///<   more atoms was set by a property, removing such
                                            ///<   information from the atom block of a V2000
                                            ///<   printout.
  bool property_element_lists;              ///< Flag to indicate that properties control the "atom
                                            ///<   lists," which are really lists of elements with
                                            ///<   an attachment point.  If this is the case, no
                                            ///<   atom list block information from the deprecated
                                            ///<   format will be written to a V2000 MOL file.
  std::vector<double3> coordinates;         ///< Cartesian coordinates of all atoms
  std::vector<char4> atomic_symbols;        ///< Symbols for all atoms
  std::vector<int> atomic_numbers;          ///< Atomic numbers for all atoms
  std::vector<int> formal_charges;          ///< Formal charges for all atoms
  std::vector<RadicalState> radicals;       ///< Indications that any atom contains a radical in
                                            ///<   singlet, doublet, or triplet excitation state
  std::vector<int> isotopic_shifts;         ///< Isotope numbers shifting each atoms nuclear mass
                                            ///<   (~1 Dalton from the most common isotopic mass)
  std::vector<MolObjAtomStereo> parities;   ///< Stereochemical parity of each atom
  std::vector<int> implicit_hydrogens;      ///< The number of implicit hydrogens that may be
                                            ///<   considered to reside bonded to a particular atom
                                            ///<   (in addition to those that are explicitly drawn)
  std::vector<bool> stereo_considerations;  ///< Indications of whether each atom has steroisomeric
                                            ///<   considerations.  Having stereochemistry at both
                                            ///<   ends of a double bond indicates that the
                                            ///<   double-bond affects the structural properties
                                            ///<   of the molecule.
  std::vector<int> valence_connections;     ///< The number of bonds made by each atom to others
                                            ///<   in the same molecule
  std::vector<int> atom_atom_mapping_count; ///< The number of atoms mapped to this one in a
                                            ///<   chemical reaction
  std::vector<bool> exact_change_enforced;  ///< Flags to indicate that the changes on an atom
                                            ///<   must be exactly as described in the reaction

  /// The assigned method of adding implicit hydrogens: add at least the number defined in the
  /// field (up to the amount needed to satisfy the anticipated valence shell electron content),
  /// add the number needed to satisfy the valence shell, or do not add hydrogens. 
  std::vector<HydrogenAssignment> hydrogenation_protocol;

  /// Indication of whether stereochemistry is inverted or retained in a chemical reaction
  std::vector<StereoRetention> orientation_stability;

  /// Bonds between atoms
  std::vector<MdlMolBond> bonds;

  /// Lists of atomic elements to be used in arbitrary operations
  std::vector<MdlMolAtomList> element_lists;

  /// Stext entries (these are deprecated and used by ISIS / Desktop applications only)
  std::vector<MolObjSTextGroup> stext_entries;
  
  /// Properties
  std::vector<MdlMolProperty> properties;

  /// Data items: this is what makes the SD file (.sdf) format so extensible.  While technically
  /// not part of the MDL MOL format, they are stored in the MdlMol for association.  Data items
  /// will be written back to new .sdf files but not .mol files based on the data in this object.
  std::vector<MdlMolDataItem> data_items;
  
  /// Title (first line from the file header)
  std::string title;

  /// Details of the software that generated this MDL .mol file (second line of the header)
  std::string software_details;

  /// General comments for the file (third line of the file header)
  std::string general_comment;

  /// An external registry number for the molecule.  This can be imparted explicitly, or imparted
  /// by adding a data item indicating an external registry number.
  std::string external_regno;

  /// \brief Validate the index of an atom query.
  ///
  /// \param index   The atom of interest, numbering starts from zero
  /// \param caller  Name of the calling function
  void validateAtomIndex(int index, const char* caller) const;
  
  /// \brief Allocate space for information in the object.
  ///
  /// Overloaded:
  ///   - Allocate as described on the MOL entry's counts line
  ///   - Allocate as described in topology abstracts
  ///
  /// \param cdk      Chemical details of the system to model
  /// \param vk       Valence parameters and connections of the system to model
  /// \param mol_idx  Index of the molecule within the topology
  /// \{
  void allocate();
  void allocate(const ChemicalDetailsKit &cdk, const NonbondedKit<double> &nbk, int mol_idx = 0);
  /// \}
  
  /// \brief Produce the correct code for an atom's isotopic shift.  While the number in the atoms
  ///        block of the V2000 format most often corresponds to the value in the actual array,
  ///        there are extreme isotopes (> +4 or < -3) that place a zero in the atom block and
  ///        must be handled by property lines.  Furthermore, if any atoms' shifts are handled in
  ///        this way, it invalidates the atom block information for all of them, so report zero
  ///        if there are any isotope-related properties.
  ///
  /// \param atom_index  Index of the atom in question
  int getIsotopicShiftCode(int atom_index) const;
  
  /// \brief Interpret the formal charge of an atom based on an integer code.  A doublet radical
  ///        can also emerge from this analysis.
  ///
  /// \param charge_in   The charge code to interpret
  /// \param atom_index  Index of the atom to which the charge applies (the radical state may also
  ///                    be set by this function)
  void interpretFormalCharge(int charge_in, int atom_index);

  /// \brief Reverse the work of interpretFormalCharge() above to obtain a code suitable for the
  ///        atom block of the V2000 format.
  ///
  /// \param atom_index  Index of the atom in question
  int getFormalChargeCode(int atom_index) const;
  
  /// \brief Interpret the stereochemical parity of an atom based on an integral numeric code.
  ///
  /// \param setting_in  The code to parse
  MolObjAtomStereo interpretStereoParity(int setting_in);

  /// \brief Interpret the implicit hydrogen count for an atom.  At least this many additional
  ///        hydrogens are implied to be present around an atom center, in addition to any hydrogen
  ///        atoms explicitly placed in the structure.
  ///
  /// \param nh_in       The number of hydrogens that are to be inferred around the atom's
  ///                    structure.  This number is decremented by one, meaning that 1 implies no
  ///                    implicit hydrogen content.
  /// \param atom_index  Atom to which the hydrogen content applies.  This is provided in order to
  ///                    set both the actual number and the action fields for the atom.
  void interpretImplicitHydrogenContent(int nh_in, int atom_index);

  /// \brief Reverse the operation of interpretImplicitHydrogenContent() above, to produce the
  ///        appropriate code for implicit hydrogens in a V2000 MDL MOL format entry atom block.
  ///
  /// \param atom_index  Atom to which the hydrogen content applies.
  int getImplicitHydrogenCode(const int atom_index) const;
  
  /// \brief Interpret a boolean value from an integer.
  ///
  /// \param code_in  The code to translate.  This one is simple: 1 = TRUE, 0 = FALSE
  /// \param desc     Description of the activity, for error tracing purposes
  bool interpretBooleanValue(int value_in, const std::string &desc);

  /// \brief Translate the count of valence interactions for a particular atom.
  ///
  /// \param count_in  The count to translate.  The input is typically unchanged, but for the fact
  ///                  that 0 and 15 both translate to no valence interactions.
  int interpretValenceNumber(int count_in);

  /// \brief Interpret the stability of stereochemical arrangements listed for each atom.
  ///
  /// \param code_in  A numeric code to be translated into inversion or retention options
  StereoRetention interpretStereoStability(int code_in);

  /// \brief Scan properties (of a V2000-format molecule) and update information that may have been
  ///        read from the atom block.
  void updateV2kAtomAttributes();
  
  /// \brief Add hydrogens to fill out valence shells based on the stated bond orders, formal
  ///        charges, and elements of the molecule.
  void hydrogenate();

  /// \brief Compare the external registry number (if any) coming from an added data item with
  ///        that stored in the molecule object itself.  This provides a check that all data items
  ///        will display a consistent external registry number.
  ///
  /// \param regno_in  The registry number to compare to anything that exists.  If this string is
  ///                  blank, no comparison will be made.  If this string is substantial and the
  ///                  object's own number is blank, take this as the object's registry number.
  ///                  If both the object's own number and this string are substantial and they
  ///                  disagree, produce a warning or error.
  void compareExternalRegistryNumbers(const std::string &regno_in);

  /// \brief Check the atom count of some external object against the internal number.
  ///
  /// \param ext_atom_count  The number of atoms coming from the external struct
  void checkAtomCount(int ext_atom_count) const;

  /// \brief Check the data item count of the object against a requested index.
  ///
  /// \param item_index  The item of interest
  /// \param caller      Name of the calling function (for backtracing purposes)
  void checkDataItemIndex(int item_index, const char* caller) const;

  /// \brief Create a property for one of the common atom / integer combination lists (i.e. M  CHG)
  ///        and add it to the MDL MOL entry.
  ///
  /// \param notable_data  List of atom indices and property pairs.  The atom indices need not be
  ///                      inflated by one, as they are going "straight into the bloodstream"
  ///                      rather than passing through filters in the MdlMolProperty constructor.
  /// \param prcode        The property code, i.e. { 'C', 'H', 'G', 'M' }
  void addAtomIntProperty(const std::vector<int2> &notable_data, const char4 prcode);
  
  /// \brief Transfer information from topology and chemical features objects to complete details
  ///        of the MDL MOL entry.  This function assumes that the object's internal arrays for
  ///        atomic properties have been allocated (resized) to their full lengths, and that the
  ///        objec'ts array of bonds have been reserved to accommodate its final length.
  ///
  /// \param chemfe          The chemical features of the system, containing a topology pointer
  /// \param molecule_index  Index of the molecule of interest within the topology
  void transferTopologicalDetails(const ChemicalFeatures *chemfe, int molecule_index = 0);
};

/// \brief Read a structure data file (.sdf extension) containing one or more MDL MOL entries.
///        Return the results as a Standard Template Library vector of MDL MOL objects.
///
/// Overloaded:
///   - Take a file name
///   - Take a pre-converted TextFile object
///
/// \param tf       
/// \param file_name
/// \param capitalization  Indicate whether atomic symbol capitalization can be ignored when
///                        assigning elements to each atom
/// \param policy          Course of action to take if errors are encountered when inferring
///                        atomic elements
/// \{
std::vector<MdlMol> readStructureDataFile(const TextFile &tf, int low_frame_limit,
                                          int high_frame_limit,
                                          CaseSensitivity capitalization = CaseSensitivity::YES,
                                          ExceptionResponse policy = ExceptionResponse::WARN);
std::vector<MdlMol> readStructureDataFile(const std::string &file_name,
                                          CaseSensitivity capitalization = CaseSensitivity::YES,
                                          ExceptionResponse policy = ExceptionResponse::WARN);

std::vector<MdlMol> readStructureDataFile(const TextFile &tf,
                                          CaseSensitivity capitalization = CaseSensitivity::YES,
                                          ExceptionResponse policy = ExceptionResponse::WARN);
/// \}
  
} // namespace structure
} // namespace stormm

#include "mdlmol.tpp"

#endif
