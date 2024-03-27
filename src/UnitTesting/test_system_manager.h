// -*-c++-*-
#ifndef STORMM_TEST_SYSTEM_MANAGER_H
#define STORMM_TEST_SYSTEM_MANAGER_H

#include "copyright.h"
#include <string>
#include <vector>
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "Constants/scaling.h"
#include "FileManagement/file_listing.h"
#include "Math/rounding.h"
#include "Numerics/split_fixed_precision.h"
#include "Random/random.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/systemcache.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"
#include "test_environment.h"

namespace stormm {
namespace testing {

using constants::ExceptionResponse;
using stmath::roundUp;
using numerics::hostDoubleToInt95;
using numerics::hostSplitFPSum;
using numerics::default_velocity_scale_bits;
using numerics::default_force_scale_bits;
using random::Xoshiro256ppGenerator;
using synthesis::Condensate;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisWriter;
using synthesis::SystemCache;
using topology::AtomGraph;
using topology::UnitCellType;
using trajectory::CoordinateFileKind;
using trajectory::CoordinateFrame;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesWriter;
using trajectory::PhaseSpace;
  
/// \brief Provide the means to read a series of topology and coordinate files, then organize them
///        into corresponding lists of AtomGraph and PhaseSpace objects.  The object can also
///        export coordinates in other formats.  Error checking is provided to abort the (test)
///        program or issue a warning, and the object can signal whether dependent tests should
///        be run or aborted.
class TestSystemManager {
public:

  /// \brief The constructor takes lists of topologies, then lists of coordinates, with a
  ///        one-to-one correspondence of files in each list.
  ///
  /// Overloaded:
  ///   - Create a blank object (so that this object can be easily containerized)
  ///   - Accept straight lists of file names for both topologies and the coordinates
  ///   - Accept a common base name and extension for the topologies and the coordinates
  ///   - Accept a common base name for the topologies and the coordinates
  /// \{
  TestSystemManager();

  TestSystemManager(const std::string &topology_base_in, const std::string &topology_extn_in,
                    const std::vector<std::string> &topology_names_in,
                    const std::string &coordinate_base_in, const std::string &coordinate_extn_in,
                    const std::vector<std::string> &coordinate_names_in,
                    ExceptionResponse policy = ExceptionResponse::WARN,
                    TestPriority fault_response_in = TestPriority::ABORT,
                    TestPriority all_go_response_in = TestPriority::CRITICAL);

  TestSystemManager(const std::string &topology_base_in,
                    const std::vector<std::string> &topology_names_in,
                    const std::string &coordinate_base_in,
                    const std::vector<std::string> &coordinate_names_in,
                    ExceptionResponse policy = ExceptionResponse::WARN,
                    TestPriority fault_response_in = TestPriority::ABORT,
                    TestPriority all_go_response_in = TestPriority::CRITICAL);

  TestSystemManager(const std::vector<std::string> &topology_names_in,
                    const std::vector<std::string> &coordinate_names_in,
                    ExceptionResponse policy = ExceptionResponse::WARN,
                    TestPriority fault_response_in = TestPriority::ABORT,
                    TestPriority all_go_response_in = TestPriority::CRITICAL);
  /// \}

  /// \brief Get the number of systems in the object.
  ///
  /// Overloaded:
  ///   - Obtain the total number of systems
  ///   - Obtain the number of systems adhering to selected boundary conditions
  ///
  /// \param query_bc  The boundary conditions of interest
  /// \{
  int getSystemCount() const;
  int getSystemCount(UnitCellType query_bc) const;
  int getSystemCount(const std::vector<UnitCellType> &query_bc) const;
  /// \}

  /// \brief Get a list of all topologies meeting a specific unit cell type.
  ///
  /// Overloaded:
  ///   - Return a list of topologies matching one selected unit cell type.
  ///   - Return a list of topologies matching any of a list of unit cell types.
  ///
  /// \param uc_choice  The unit cell type of interest
  /// \{
  std::vector<int> getQualifyingSystems(UnitCellType uc_choice) const;
  std::vector<int> getQualifyingSystems(const std::vector<UnitCellType> &uc_choice) const;
  /// \}
  
  /// \brief Return the base path for topologies
  std::string getTopologyBasePath() const;
  
  /// \brief Return the base path for coordinates
  std::string getCoordinateBasePath() const;

  /// \brief Return the extension for topologies
  std::string getTopologyExtension() const;

  /// \brief Return the extension for coordinates
  std::string getCoordinateExtension() const;
  
  /// \brief Return the full name of a topology file according to some index in the list.
  std::string getTopologyFile(int index) const;

  /// \brief Return the full name of a coordinate file according to some index in the list.
  std::string getCoordinateFile(int index) const;

  /// \brief Get the planned course of action for subsequent tests in the event that a file is
  ///        non-existent, or for any reason unreadable.
  ///
  /// Overloaded:
  ///   - Return a general status for any and all systems covered by the manager
  ///   - Return a specific status for an individual system by index
  ///   - Return the overall status for a collections of systems, specified by index
  ///
  /// \param index    Index of the system of interest
  /// \param indices  Indices of the systems of interest
  /// \{
  TestPriority getTestingStatus() const;
  TestPriority getTestingStatus(int index) const;
  TestPriority getTestingStatus(const std::vector<int> &indices) const;
  /// \}

  /// \brief Get a copy of the coordinates for one system as a CoordinateFrame object.
  ///
  /// \param index  The system of interest
  CoordinateFrame exportCoordinateFrame(int index) const;

  /// \brief Get a copy of the coordinates for one system as a PhaseSpace object.
  ///
  /// \param index  The system of interest
  PhaseSpace exportPhaseSpace(int index) const;

  /// \brief Export the topology for one system.
  ///
  /// \param index  The system of interest.
  AtomGraph exportAtomGraph(int index) const;
  
  /// \brief Get a const reference to the coordinates for one or more systems.
  ///
  /// Overloaded:
  ///   - Get a reference to a single coordinate set
  ///   - Get a const reference to the vector of all coordinate sets
  ///
  /// \param index  The system of interest
  /// \{
  PhaseSpace& viewCoordinates(int index);
  const std::vector<PhaseSpace>& viewCoordinates();
  /// \}

  /// \brief Get a const reference to the topology for one or more systems.
  ///
  /// Overloaded:
  ///   - Get a reference to a single topology
  ///   - Get a reference to the vector of all topologies in order
  ///
  /// \param index  The system of interest
  /// \{
  const AtomGraph& getTopologyReference(int index) const;
  const std::vector<AtomGraph>& getTopologyReference() const;
  /// \}
  
  /// \brief Get a pointer to one or more topologies
  ///
  /// Overloaded:
  ///   - Get a pointer to a single topology
  ///   - Get a vector of pointers to all topologies in order
  ///   - Return const results for a const object or non-const results for a non-const object
  ///
  /// \param index  Identifier of a specific topology of interest
  /// \{
  AtomGraph* getTopologyPointer(int index);
  const AtomGraph* getTopologyPointer(int index) const;
  std::vector<AtomGraph*> getTopologyPointer();
  const std::vector<AtomGraph*> getTopologyPointer() const;
  /// \}

  /// \brief Get a coordinate series out of the object, based on one system from its repository.
  ///
  /// \param base_system         
  /// \param frame_count         The number of frames in the series, each a copy of the original
  ///                            coordinate set
  /// \param perturbation_sigma  The Gaussian width by which to perturb coordinates
  /// \param xrs_seed            Random number generator seed
  /// \param scale_bits          The number of bits after the decimal to include, if the data type
  ///                            of the coordinate series is integral (fixed precision)
  template <typename T>
  CoordinateSeries<T> exportCoordinateSeries(const int base_system, int frame_count,
                                             double perturbation_sigma = 0.0,
                                             int xrs_seed = 915083, int scale_bits = 0) const;

  /// \brief Get a synthesis of coordinates with fixed-precision storage ready for positions,
  ///        velocities, and forces.
  ///
  /// Overloaded:
  ///   - Select a list of systems based on their indices within the manager
  ///   - Select all systems conforming to one or more specific types of unit cell (a critical
  ///     factor in building a synthesis is that all systems use similar boundary conditions)
  ///
  /// Parameter descriptions follow from getCoordinateSeries() above, with the addition of: 
  ///
  /// \param index_key  A series of system indices from within the TestSystemManager with which
  ///                   to compose the synthesis.  Topologies and coordinates will be drawn upon.
  /// \param uc_choice  The choice or choices of unit cell to select (one instance of each matching
  ///                   system will be included in the result)
  /// \param vel_bits   The number of bits in the fraction for the velocity representation
  /// \param frc_bits   The number of bits in the fraction for the force representation
  /// \{
  PhaseSpaceSynthesis exportPhaseSpaceSynthesis(const std::vector<int> &index_key,
                                                double perturbation_sigma = 0.0,
                                                int xrs_seed = 7108262, int scale_bits = 40,
                                                int vel_bits = default_velocity_scale_bits,
                                                int frc_bits = default_force_scale_bits) const;
  PhaseSpaceSynthesis exportPhaseSpaceSynthesis(UnitCellType uc_choice) const;
  PhaseSpaceSynthesis exportPhaseSpaceSynthesis(const std::vector<UnitCellType> &uc_choice) const;
  /// \}
  
  /// \brief Get a synthesis of topologies based on a series of system indices from the manager.
  ///
  /// Overloaded:
  ///   - Select a list of systems based on their indices within the manager
  ///   - Select all systems conforming to one or more specific types of unit cell (a critical
  ///     factor in building a synthesis is that all systems use similar boundary conditions)
  ///
  /// Descriptions of parameters follow from getPhaseSpaceSynthesis() above, with the addition of:
  ///
  /// \param policy  Specify whether to report abnormal input.  In most cases, the
  ///                AtomGraphSynthesis will report being given a longer list of topologies than
  ///                was required (unused topologies).  For most testing purposes, this is benign,
  ///                and should not be raised to the end user's attention.
  /// \{
  AtomGraphSynthesis
  exportAtomGraphSynthesis(const std::vector<int> &index_key,
                           ExceptionResponse policy = ExceptionResponse::SILENT) const;

  AtomGraphSynthesis
  exportAtomGraphSynthesis(UnitCellType uc_choice,
                           ExceptionResponse policy = ExceptionResponse::SILENT) const;

  AtomGraphSynthesis
  exportAtomGraphSynthesis(const std::vector<UnitCellType> &uc_choice,
                           ExceptionResponse policy = ExceptionResponse::SILENT) const;
  /// \}

  /// \brief Produce a SystemCache based on a compilation of the systems read into the manager.
  ///
  /// Overloaded:
  ///   - Accept a list of system indices, and optionally labels to apply
  ///   - Accept a string input for a &files namelist
  ///
  /// \param index_key    A series of system indices from within the TestSystemManager with which
  ///                     to compose the cache.  Topologies and coordinates will be drawn upon.
  /// \param x_name_list  A list of trajectory file names for each cache item
  /// \param x_type_list  A list of file types for each trajectory file to write
  /// \param r_name_list  A list of restart file names for each cache item
  /// \param r_type_list  A list of file types for each restart file to write
  /// \param label_list   A list of labels to apply to each cache item.  If provided, this array
  ///                     have the same length as index_key.
  /// \param fnml         A string containing the contents of a &files namelist
  /// \{
  SystemCache exportSystemCache(const std::vector<int> &index_key, const TestEnvironment &oe,
                                const std::vector<std::string> &x_name_list = {},
                                const std::vector<CoordinateFileKind> &x_type_list = {},
                                const std::vector<std::string> &r_name_list = {},
                                const std::vector<CoordinateFileKind> &r_type_list = {},
                                const std::vector<std::string> &label_list = {});

  SystemCache exportSystemCache(const std::string &fnml);
  /// \}
  
private:
  int system_count;              ///< The number of systems managed by this object
  std::string topology_base;     ///< Common base name for all topology files (a directory
                                 ///<   character will be placed between this and any subsequent
                                 ///<   topology file names)
  std::string topology_extn;     ///< Common extension for all topology files (a dot '.' will come
                                 ///<   between this and any preceding coordinate file names)
  std::string coordinate_base;   ///< Common base name for all coordinate files
  std::string coordinate_extn;   ///< Common extension for all coordinate files
  TestPriority fault_response;   ///< Action to take if files from either list are missing
  TestPriority all_go_response;  ///< Action to take if all files are present and parseable
  bool fault_found;              ///< Set to TRUE if any files are unreadable for any reason (if
                                 ///<   any element of the subsequent <bool>-type vectors are TRUE)

  // Vectors of names and flags describing the file parsing
  std::vector<std::string> topology_names;    ///< Names for all topology files
  std::vector<std::string> coordinate_names;  ///< Names for all coordinate files
  std::vector<bool> topology_exist;           ///< Indications that topology files were found
  std::vector<bool> topology_success;         ///< Indications of successful topology parsing
  std::vector<bool> coordinate_exist;         ///< Indications that coordinate files were found
  std::vector<bool> coordinate_success;       ///< Indications of successful coordinate parsing
  std::vector<bool> compatibility;            ///< Indications that the topologies and coordinates
                                              ///<   for corresponding systems are compatible

  /// All topologies (this array contains empty topologies if the corresponding files were
  /// unreadable)
  std::vector<AtomGraph> all_topologies;

  /// All coordinates (this array contains empty coordinates objects if the corresponding files
  /// were unreadable)
  std::vector<PhaseSpace> all_coordinates;

  /// \brief Test that an index does not exceed the valid range for this object.
  ///
  /// \param index   The index of a topology or coordinate file to access
  /// \param caller  Name of the calling function, for error tracing purposes
  void checkIndexing(int index, const char* caller) const;
};

} // namespace testing
} // namespace stormm

#include "test_system_manager.tpp"

#endif
