// -*-c++-*-
#ifndef STORMM_NETCDF_LINK_H
#define STORMM_NETCDF_LINK_H

#ifdef STORMM_INCLUDE_NETCDF
#  include <netcdf.h>
#endif
#include <string>
#include <vector>
#include "copyright.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Parsing/parsing_enumerators.h"
#include "Parsing/polynumeric.h"

#if STORMM_INCLUDE_NETCDF
namespace stormm {
namespace trajectory {

using parse::NumberFormat;
using parse::PolyNumeric;

/// \brief Store the critical details of a single number in a type-agnostic format using a union,
///        a type enumeration, and the character string to tag the variable.
class EcumenicalScalar {
public:

  /// \brief The constructor requires a value, type, and tag.
  template <typename T> EcumenicalScalar(T value_in, size_t type_code_in,
                                         const std::string &key_string_in);

  /// \brief The class is written in C and C++, with no const elements.  The default copy and move
  ///        constructors, as well as assignment operators, will apply.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object found on the right hand side of the assignment statement
  /// \{
  EcumenicalScalar(const EcumenicalScalar &original) = default;
  EcumenicalScalar(EcumenicalScalar &&original) = default;
  EcumenicalScalar& operator=(const EcumenicalScalar &original) = default;
  EcumenicalScalar& operator=(EcumenicalScalar &&original) = default;

  /// \brief Get the name identifier of the array.  This will be the identifier used to tag the
  ///        array in the NetCDF file.
  const std::string& getKeyString() const;

  /// \brief Get the value of the constant.  The return type will be checked for validity.
  template <typename T> T getValue() const;

  /// \brief Get the type code, for inspection.
  size_t getTypeCode() const;
  
  /// \brief Get the type name, for inspection.
  size_t getTypeName() const;

  /// \brief Assign the NetCDF variable code for the scalar value, with codes for the size of the
  ///        type.  This requires an active NetCDF file.
  ///
  /// \param netcdf_file_id  The identifier for the active NetCDF file into which this array will
  ///                        be written
  /// \param filename        Name of the file being written (for error tracing purposes)
  void assignNCVar(int netcdf_file_id,
                   const std::string &filename = std::string("not conveyed by the program's "
                                                             "developer"));

  /// \brief Write the scalar variable to a NetCDF file.
  void write(int netcdf_file_id, const std::string &filename = std::string("not conveyed by the "
                                                                           "program's developer"));
  
private:

  std::vector<uint8_t> value;  ///< The value of the constant, encoded as a byte string
  int value_length;            ///< An explicit note of the size of the data type contained in the
                               ///<   "scalar", enabling distinctions between various tuple values
                               ///<   and other atomic types.  This is obtained from sizeof(T).
  size_t type_code;            ///< The type index
  std::string key_string;      ///< The tag to associate with this constant in the NetCDF file
  std::string type_name;       ///< Human-readable name of the original data type, for error
                               ///<   reporting purposes
  int nc_dim_code[1];          ///< NetCDF code for storing or retrieving the length of data in
                               ///<   the scalar value definition
  int nc_var_code;             ///< NetCDF code for storing or retrieving the value from a file
};

/// \brief Store the critical details of an array in a type-agnostic format, using C's
///        pointer-array duality, a type identifier code, and sizing constants.  While the
///        underlying array is expected to be contiguous in memory, the array can be cast as having
///        up to four dimensions.
class EcumenicalArray {
public:

  /// \brief The constructor assumes a one-dimensional array, but more dimensions may be added.
  ///        The format will accept up to four dimensions, with the fastest varying expected first
  ///        and slower varying dimensions included later.
  EcumenicalArray(const void* tip_in, size_t type_code_in, const std::string &key_string_in,
                  size_t first_dimension_in, size_t second_dimension_in = 0,
                  size_t third_dimension_in = 0, size_t fourth_dimension_in = 0,
                  size_t io_batch_size_in = 4194304);

  /// \brief The class is written in C and C++, with no const elements.  The default copy and move
  ///        constructors, as well as assignment operators, will apply.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object found on the right hand side of the assignment statement
  /// \{
  EcumenicalArray(const EcumenicalArray &original) = default;
  EcumenicalArray(EcumenicalArray &&original) = default;
  EcumenicalArray& operator=(const EcumenicalArray &original) = default;
  EcumenicalArray& operator=(EcumenicalArray &&original) = default;
  /// \}
  
  /// \brief Get the name identifier of the array.  This will be the identifier used to tag the
  ///        array in the NetCDF file.
  const std::string& getKeyString() const;

  /// \brief Get the number of dimensions in the underlying array
  int getDimensionality() const;
  
  /// \brief Get the length of the fastest-incrementing dimension
  size_t getFirstIncrementingSize() const;

  /// \brief Get the length of the second fastest-incrementing dimension
  size_t getSecondIncrementingSize() const;

  /// \brief Get the length of the third fastest-incrementing dimension
  size_t getThirdIncrementingSize() const;

  /// \brief Get the length of the fourth fastest-incrementing dimension
  size_t getFourthIncrementingSize() const;

  /// \brief Get the length of one of the array dimensions.
  ///
  /// \param order  The order of the index to retrieve, beginning at 0 (fastest incrementing) and
  ///               reaching up to 3 (slowest incrementing).  A check will be applied to ensure
  ///               the array has as many dimensions as are being requested.
  size_t size(int order = 0) const;

  /// \brief Get the IO batch size used by the array.
  size_t getIOBatchSize() const;

  /// \brief Get the maximum dimensions of the IO area / volume to be read or written in any one
  ///        disk transaction involving NetCDF.  This is for inspection purposes.
  std::vector<size_t> getIOLengths() const;
  
  /// \brief Get a pointer to the array's data at a particular point along any of its higher
  ///        indices.
  ///
  /// Overloaded:
  ///   - Obtain a const-qualified pointer from a const-qualified object
  ///   - Obtain a mutable pointer from a non-cosnt object
  ///
  /// \param second_index  Optional index along the array's second-fastest incrementing dimension
  /// \param third_index   Optional index along the array's third-fastest incrementing dimension
  /// \param fourth_index  Optional index along the slowest incrementing dimension
  /// \{
  template <typename T> T* getPointer(int second_index = 0, int third_index = 0,
                                      int fourth_index = 0);

  template <typename T>
  const T* getPointer(int second_index = 0, int third_index = 0, int fourth_index = 0) const;
  /// \}
  
  /// \brief Read a value from the underlying array.
  ///
  /// \param pos_first   The index to read from the fastest-incrementing dimension of the array
  /// \param pos_second  Optional index to read from the second fastest-incrementing dimension
  /// \param pos_third   Optional index to read from the second fastest-incrementing dimension
  /// \param pos_fourth  Optional index to read from the second fastest-incrementing dimension
  template <typename T> T readValue(size_t pos_first, size_t pos_second = 0, size_t pos_third = 0,
                                    size_t pos_fourth = 0) const;
  
  /// \brief Set the value of an element of the underlying array
  template <typename T>
  void putValue(T value, size_t pos_first, size_t pos_second = 0, size_t pos_third = 0,
                size_t pos_fourth = 0);

  /// \brief Assign the NetCDF variable code for the array, and codes for all dimensions (the size
  ///        of the array along each dimension gets its own identifier in the NetCDF meta data).
  ///        This requires an active NetCDF file.
  ///
  /// \param netcdf_file_id  The identifier for the active NetCDF file into which this array will
  ///                        be written
  /// \param filename        Name of the file being written (for error tracing purposes)
  void assignNCVar(const int netcdf_file_id,
                   const std::string &filename = std::string("not conveyed by the program's "
                                                             "developer"));

  /// \brief Write the array variable, or a subset of it, to a NetCDF file.
  ///
  /// \param netcdf_file_id       The NetCDF file identifier
  /// \param filename             Name of the NetCDF file being written.  This is provided for
  ///                             error tracing purposes.
  /// \param cornerstone_indices  The starting indices of a subsection of the array to write.  One
  ///                             index for each dimension must be provided.  The entire array will
  ///                             be written to disk if this holds no data.
  /// \param edge_lengths         Lengths of the subset to write.  The entire array will be written
  ///                             to disk if this holds no data.
  void write(int netcdf_file_id, const std::string &filename = std::string("not conveyed by the "
                                                                           "program's developer"),
             const std::vector<size_t> &cornerstone_indices = {},
             const std::vector<size_t> &edge_lengths = {});

  /// \brief Read the array variable, or a subset of it, from a NetCDF file.  Descriptions of input
  ///        parameters follow from write(), above.
  void read(int netcdf_file_id, const std::string &filename = std::string("not conveyed by the "
                                                                          "program's developer"),
            const std::vector<size_t> &cornerstone_indices = {},
            const std::vector<size_t> &edge_lengths = {});
  
private:

  size_t type_code;        ///< The codified type index for the data type in the array, according
                           ///<   to STORMM's internal definitions created at runtime.  See the
                           ///<   files src/DataTypes/common_types.h and
                           ///<   src/DataTypes/stormm_vector_types.h.   Like Hybrid objects, this
                           ///<   is limited to certain basic data types, such as express integer
                           ///<   or floating point data, possibly tuples thereof.
  int nc_type_code;        ///< The codified type index for the data type in the array.  This is
                           ///<   assigned according to the NetCDF package's native definitions.
  int nc_var_code;         ///< Code assigned to the array variable in a NetCDF file
  int nc_dim_codes[4];     ///< Codes assigned to each dimension of the array in a NetCDF file
  std::string key_string;  ///< The human-readable string name for this array, conveyed in the
                           ///<   NetCDF file meta data
  size_t length[4];        ///< The length(s) of the underlying array.  The object stores a pointer
                           ///<   to the head of the array and understands its length.  It has the
                           ///<   ability to change elements of the underlying array but not its
                           ///<   size or dimensionality.  These are the lengths of the related
                           ///<   arrays in STORMM.
  size_t nc_length[4];     ///< The length(s) of the array as will be written to a NetCDF file.
                           ///<   Tuple data types (e.g. double2) will count 2, 3, or 4 times as
                           ///<   many array elements to pack the data into an array of the
                           ///<   base atomic type.
  void* tip;               ///< A pointer to the start of the underlying array.  This must be
                           ///<   interpreted in the context of the data type code (type_code).

  /// Multi-dimensional arrays may need to be written in pieces, as NetCDF expects C-ordered
  /// arrays whereas STORMM keeps arrays in Fortran order.  This batch size will determine how to
  /// size each patch.  This is a number of units of the original data type, not a total byte
  /// count.
  size_t io_batch_size;

  /// The maximum portion of the array which will be read or written at any given time.  This is
  /// relevant only for multi-dimensional arrays.  STORMM stores such data in Fortran order,
  /// whereas NetCDF expects the data to be in C order.  A meticulous dance is needed to let NetCDF
  /// read and write data into Fortran-format arrays.  This arrays tracks with the length array,
  /// not nc_length: it shows the number of elements of the underlying data type that will be read
  /// read or written in each batch.
  std::vector<size_t> io_length;

  /// \brief Check the stated dimensionality.  While size_t values are unsigned, an unreasonably
  ///        large number can still be caught.
  ///
  /// \param length_in  The length to check
  void validateLength(size_t length_in);

  /// \brief Validate arrays specifying that a subset of the array be written to or read from disk.
  ///
  /// \param cornerstone_indices  Origin indices of the subset to be written, in each dimension
  /// \param edge_lengths         The length of the subset to be written, in each dimension
  bool validateSubsetDirective(const std::vector<size_t> &cornerstone_indices,
                               const std::vector<size_t> &edge_lengths);
  
  /// \brief A templated function will create pointers to a buffer of data as well as the
  ///        underlying array itself, to stash the Fortran-ordered data and protect it from
  ///        unintentional overwriting.
  ///
  ///
  /// \param buffer  A void-casted pointer to memory allocated for temporary storage of the array's
  ///                original values
  template <typename T> void stashSubset(const std::vector<size_t> &cornerstone,
                                         const std::vector<size_t> &sbs_lengths, void* buffer);

  /// \brief Unroll the templated function for preparing a multi-dimensional array to write part
  ///        of its contents to a NetCDF file.
  ///
  void unrollWritingPrep(const std::vector<size_t> &cornerstone,
                         const std::vector<size_t> &sbs_lengths,
                         void* buffer);
};

/// \brief This class stores a series of ecumenical pointers (void*, with type IDs and length
///        constants to define the true content behind them) to translate data from STORMM objects
///        to Network Commond Data Form (NetCDF).  The class also serves to transfer data from
///        NetCDF to STORMM class objects.  In a manner bearing vague similarities to the
///        NamelistEmulator, the workflow is to configure an object of this class and then fill it
///        with the contents (or, more precisely, pointers and other critical constants based on
///        the abstract) of a class object in STORMM.  One free function will take another STORMM
///        class object and produce an object of this class, while another free function will take
///        an object of this class and fill another STORMM class object.
class NetCDFLink {
public:

  /// \brief To avoid dependencies on every class that the link might serve, the class constructor
  ///        will take arrays of ecumenical pointers, type IDs, and length constants.  Member
  ///        functions will serve to add new array pointers, individual constants, and other
  ///        metadata.
  /// \{
  NetCDFLink();

  NetCDFLink(const std::vector<EcumenicalScalar> &scalar_contents_in,
             const std::vector<EcumenicalArray> &tensor_contents_in,
             const std::string &filename_in);
  /// \}

  /// \brief The class is written in C and C++, with no const elements.  The default copy and move
  ///        constructors, as well as assignment operators, will apply.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object found on the right hand side of the assignment statement
  /// \{
  NetCDFLink(const NetCDFLink &original) = default;
  NetCDFLink(NetCDFLink &&original) = default;
  NetCDFLink& operator=(const NetCDFLink &original) = default;
  NetCDFLink& operator=(NetCDFLink &&original) = default;

  /// Get a pointer to one of the arrays comprised by the interface.
  ///
  /// Overloaded:
  ///   - Get the pointer to mutable data for a non-const object.
  ///   - Get a const pointer to data held within a const-qualified object.
  ///
  /// \param 
  /// \{
  template <typename T> T* getContentPointer();
  template <typename T> const T* getContentPointer() const;
  /// \}

  /// \brief Set the file name to read or write
  ///
  /// \param filename_in  The name of the file in question
  void setFileName(const std::string &filename_in);

  /// \brief Add a scalar variable to the list.
  ///
  /// Overloaded:
  ///   - Add the pre-formed object
  ///   - Provide the necessary parameters to construct a new object at the back of the current
  ///     list
  ///
  /// \param value  The value to append to the growing list
  /// \param kind   The kind of number represented by the input parameter value
  /// \{
  void addScalar(EcumenicalScalar);
  void addScalar(PolyNumeric value, NumberFormat kind, const std::string &var_name);
  /// \}
  
  /// \brief Add an array variable to the list.
  ///
  /// \param value  The array to append to the growing list
  
private:
  std::vector<EcumenicalScalar> scalar_contents;  ///< A list of scalars that will be encoded in
                                                  ///<   the resulting NetCDF file
  std::vector<EcumenicalArray> array_contents;    ///< A list of arrays that will be encoded in the
                                                  ///<   resulting NetCDF file
  std::string filename;                           ///< The name of the NetCDF file to read or write
};

/// \brief Convert one of STORMM's data type index codes (recorded at runtime) into a NetCDF data
///        type identifier.
///
/// \param stormm_type_index  The runtime-derived type index to convert
int assignNCTypeIndex(size_t stormm_type_code);

/// \brief Check the return code from a NetCDF definition call.  This will print detailed error
///        messages and help trace the problem to a specific file if multiple files are in
///        production.
///
/// \param chk_code         The return code of the NetCDF defining function
/// \param varname          Name of the variable being defined
/// \param filename         The file be written
/// \param array_lengths    Lengths of the array along each of its dimensions
/// \param dimension_count  The number of dimensions in the underlying array
/// \param class_caller     Name of a class calling the routine.  This is expected to be either
///                         EcumenicalScalar or EcumenicalArray, but if the developer is calling
///                         the routine from a free function, the name of the free function should
///                         be supplied as the "class" caller.
/// \param method_caller    Name of the member function within any class calling the validator.
///                         Like class_caller, this serves error tracing purposes.
void validateNCDefinitionCode(const int chk_code, const std::string &varname,
                              const std::string &filename = std::string("not conveyed by the "
                                                                        "program's developer"),
                              const int* array_lengths = nullptr, const int dimension_count = 1,
                              const char* class_caller = nullptr,
                              const char* method_caller = nullptr);

/// \brief Check the return code from a NetCDF write instruction.  This will print detailed error
///        messages and help trace the problem to a specific file if multiple files are in
///        production.  Descriptions of input parameters follow from validateNCDefinitionCode(),
///        above.
void validateNCWriteCode(const int chk_code, const std::string &varname,
                         const std::string &filename = std::string("not conveyed by the program's "
                                                                   "developer"),
                         const int* array_lengths = nullptr, const int dimension_count = 1,
                         const char* class_caller = nullptr,
                         const char* method_caller = nullptr);

} // namespace trajectory
} // namespace stormm

#include "netcdf_link.tpp"
#endif // STORMM_INCLUDE_NETCDF

#endif
