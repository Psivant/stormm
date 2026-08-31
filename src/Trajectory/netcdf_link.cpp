#include "copyright.h"
#include "Reporting/error_format.h"
#include "netcdf_link.h"

#if STORMM_USE_NETCDF
namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
const std::string& EcumencialScalar::getKeyString() const {
  return key_string;
}

//-------------------------------------------------------------------------------------------------
size_t EcumenicalScalar::getTypeCode() const {
  return type_code;
}

//-------------------------------------------------------------------------------------------------
const std::string& getTypeName() const {
  return type_name;
}

//-------------------------------------------------------------------------------------------------
void EcumenicalScalar::assignNCVar(const int netcdf_file_id, const std::string &filename) {
  int chk_code;
  const int ndims = this->getDimensionality();
  const std::string dim_str = key_string + "_bytes";
  chk_code = nc_def_dim(netcdf_file_id, dim_str, value_length, nc_dim_code);
  const int* vl_ptr = &value_length;
  validateDefinitionCode(chk_code, filename, vl_ptr, 1, "EcumenicalScalar", "assignNCVar");
  chk_code = nc_def_var(netcdf_file_id, key_string.c_str(), NC_UBYTE, 1, nc_dim_code,
                        &nc_var_code);
  validateNCDefinitionCode(chk_code, filename, vl_ptr, 1, "EcumenicalScalar", "assignNCVar");
}

//-------------------------------------------------------------------------------------------------
void EcumenicalScalar::write(int netcdf_file_id, const std::string &filename) {
  const int chk_code = nc_put_var_ubyte(netcdf_file_id, nc_var_code, value.data());
  validateNCWriteCode(chk_code, key_string, filename, "EcumenicalScalar", "write");
}

//-------------------------------------------------------------------------------------------------
void EcumenicalScalar::read(int netcdf_file_id, const std::string &filename) {
  const size_t start_place = 0;
  const size_t total_count = value_length;
  const int chk_code = nc_get_vara_ubyte(netcdf_file_id, nc_var_code, &start_place, &total_count,
                                         value.data());
  validateNCReadCode(chk_code, key_string, filename, "EcumenicalScalar", "read", &start_place,
                     &total_count, 1);
}
  
//-------------------------------------------------------------------------------------------------
EcumenicalArray::EcumenicalArray(const void* tip_in, const size_t type_code_in,
                                 const std::string &key_string_in, const size_t first_dimension_in,
                                 const size_t second_dimension_in, const size_t third_dimension_in,
                                 const size_t fourth_dimension_in, const size_t io_batch_size_in) :
    type_code{type_code_in},
    nc_type_code{assignNCTypeIndex(type_code_in)},
    nc_var_code{0}, nc_dim_codes{}, key_string{key_string_in}, length{}, nc_length{},
    io_batch_size{io_batch_size_in}, io_length{}, dimension_count{1},
    tip{const_cast<void*>(tip_in)}
{
  validateLength(first_dimension_in);
  validateLength(second_dimension_in);
  validateLength(third_dimension_in);
  validateLength(fourth_dimension_in);
  length[0] = first_dimension_in;
  length[1] = second_dimension_in;
  length[2] = third_dimension_in;
  length[3] = fourth_dimension_in;
  for (int i = 0; i < 4; i++) {
    nc_dim_codes[0] = 0;
    nc_length[i] = length[i];
  }
  dimension_count = this->getDimensionality();
  
  // Prepare to pack tuple data types in the NetCDF file
  if (type_code ==   int2_type_index || type_code ==   double2_type_index ||
      type_code == float2_type_index || type_code ==     char2_type_index ||
      type_code == short2_type_index || type_code ==   ushort2_type_index ||
      type_code == uchar2_type_index || type_code == longlong2_type_index ||
      type_code == ulonglong2_type_index) {
    nc_length[dimension_count - 1] *= 2;
  }
  else if (type_code ==   int3_type_index || type_code ==   double3_type_index ||
           type_code == float3_type_index || type_code ==     char3_type_index ||
           type_code == short3_type_index || type_code ==   ushort3_type_index ||
           type_code == uchar3_type_index || type_code == longlong3_type_index ||
           type_code == ulonglong3_type_index) {
    nc_length[dimension_count - 1] *= 3;
  }
  else if (type_code ==   int4_type_index || type_code ==   double4_type_index ||
           type_code == float4_type_index || type_code ==     char4_type_index ||
           type_code == short4_type_index || type_code ==   ushort4_type_index ||
           type_code == uchar4_type_index || type_code == longlong4_type_index ||
           type_code == ulonglong4_type_index) {
    nc_length[dimension_count - 1] *= 4;
  }

  // Prepare to read and write a multi-dimensional array
  io_length.resize(dimension_count);
  if (dimension_count > 1) {
    size_t tmp_volume = 1;
    int dim_depth = 0;
    while (dim_depth < ndims) {
      if (length[dim_depth] * tmp_volume < io_batch_size) {
        io_length[dim_depth] = length[dim_depth];
        tmp_volume *= length[dim_depth];
        dim_depth++;
      }
      else {
        io_length[dim_depth] = io_batch_size / tmp_volume;
        tmp_volume *= length[dim_depth];
        for (int i = dim_depth + 1; i < ndims; i++) {
          io_length[i] = 1;
        }
        dim_depth = ndims;
      }
    }
  }
  else {
    io_length[0] = length[0];
  }
}

//-------------------------------------------------------------------------------------------------
const std::string& EcumenicalArray::getKeyString() const {
  return key_string;
}

//-------------------------------------------------------------------------------------------------
size_t EcumenicalArray::getTypeCode() const {
  return type_code;
}

//-------------------------------------------------------------------------------------------------
const std::string& EcumenicalArray::getTypeName() const {
  return type_name;
}

//-------------------------------------------------------------------------------------------------
int EcumenicalArray::getDimensionality() const {
  if (length[3] > 0) {
    if (length[2] == 0 || length[1] == 0 || length[0] == 0) {
      rtErr("The fourth dimension has nonzero length, but the lower dimensions have lengths " +
            std::to_string(length[0]) + ", " + std::to_string(length[1]) + ", and " +
            std::to_string(length[2]) + ".", "EcumenicalArray", "getDimensionality");
    }
    return 4;
  }
  else if (length[2] > 0) {
    if (length[1] == 0 || length[0] == 0) {
      rtErr("The third dimension has nonzero length, but the lower dimensions have lengths " +
            std::to_string(length[0]) + " and " + std::to_string(length[1]) + ".",
            "EcumenicalArray", "getDimensionality");
    }
    return 3;
  }
  else if (length[1] > 0) {
    if (length[0] == 0) {
      rtErr("The second dimension has nonzero length, but the first dimension has length zero.",
            "EcumenicalArray", "getDimensionality");
    }
    return 2;
  }
  else {
    return 1;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
size_t EcumenicalArray::getFirstIncrementingSize() const {
  return length[0];
}

//-------------------------------------------------------------------------------------------------
size_t EcumenicalArray::getSecondIncrementingSize() const {
  if (length[1] == 0) {
    rtErr("The underlying array does not have a second dimension.", "EcumenicalArray",
          "getSecondIncrementingSize");
  }
  return length[1];
}

//-------------------------------------------------------------------------------------------------
size_t EcumenicalArray::getThirdIncrementingSize() const {
  if (length[2] == 0) {
    rtErr("The underlying array does not have a third dimension.", "EcumenicalArray",
          "getThirdIncrementingSize");
  }
  return length[2];
}

//-------------------------------------------------------------------------------------------------
size_t EcumenicalArray::getFourthIncrementingSize() const {
  if (length[3] == 0) {
    rtErr("The underlying array does not have a fourth dimension.", "EcumenicalArray",
          "getFourthIncrementingSize");
  }
  return length[3];
}

//-------------------------------------------------------------------------------------------------
size_t EcumenicalArray::size(const int order) const {
  if (order < 4) {
    return length[order];
  }
  else {
    rtErr("The underlying array cannot have " + std::to_string(order + 1) + " dimensions.",
          "EcumenicalArray", "size");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
size_t EcumenicalArray::getIOBatchSize() const {
  return io_batch_size;
}

//-------------------------------------------------------------------------------------------------
std::vector<size_t> EcumenicalArray::getIOLengths() const {
  std::vector<size_t> result(dimension_count);
  for (int i = 0; i < dimension_count; i++) {
    result[i] = io_length[i];
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
void EcumenicalArray::assignNCVar(const int netcdf_file_id, const std::string &filename) {
  int chk_code;
  const int ndims = this->getDimensionality();
  for (int i = 0; i < ndims; i++) {
    const std::string dim_str = key_string + "_" + std::to_string(i);
    chk_code = nc_def_dim(netcdf_file_id, dim_str, nc_length[i], &nc_dim_codes[i]);
    validateNCDefinitionCode(chk_code, filename, nc_length, ndims, "EcumenicalArray",
                             "assignNCVar");
  }
  chk_code = nc_def_var(netcdf_file_id, key_string.c_str(), nc_type_code, ndims, nc_dim_codes,
                        nc_var_code);
  validateNCDefinitionCode(chk_code, filename, nc_length, ndims, "EcumenicalArray", "assignNCVar");
}
  
//-------------------------------------------------------------------------------------------------
void EcumenicalArray::write(const int netcdf_file_id, const std::string &filename,
                            const std::vector<size_t> &cornerstone_indices,
                            const std::vector<size_t> &edge_lengths) {

  // Check for patch printing commands
  const bool subset = validateSubsetDirective();
  if (dimension_count > 1) {
    if (type_code == char_type_index) {
      writeBySubsets<char, char>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == short_type_index) {
      writeBySubsets<short, short>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == int_type_index) {
      writeBySubsets<int, int>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == llint_type_index) {
      writeBySubsets<llint, llint>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == float_type_index) {
      writeBySubsets<float, float>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == double_type_index) {
      writeBySubsets<double, double>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == uchar_type_index) {
      writeBySubsets<uchar, uchar>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == ushort_type_index) {
      writeBySubsets<ushort, ushort>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == uint_type_index) {
      writeBySubsets<uint, uint>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == ullint_type_index) {
      writeBySubsets<ullint, ullint>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == char2_type_index) {
      writeBySubsets<char2, char>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == short2_type_index) {
      writeBySubsets<short2, short>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == int2_type_index) {
      writeBySubsets<int2, int>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == llint2_type_index) {
      writeBySubsets<llint2, llint>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == float2_type_index) {
      writeBySubsets<float2, float>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == double2_type_index) {
      writeBySubsets<double2, double>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == ushort2_type_index) {
      writeBySubsets<ushort2, ushort>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == uint2_type_index) {
      writeBySubsets<uint2, uint>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == ullint2_type_index) {
      writeBySubsets<ullint2, ullint>(netcdf_file_id, filename, cornerstone_indices,
                                       edge_lengths);
    }
    else if (type_code == char3_type_index) {
      writeBySubsets<char3, char>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == short3_type_index) {
      writeBySubsets<short3, short>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == int3_type_index) {
      writeBySubsets<int3, int>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == llint3_type_index) {
      writeBySubsets<llint3, llint>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == float3_type_index) {
      writeBySubsets<float3, float>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == double3_type_index) {
      writeBySubsets<double3, double>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == ushort3_type_index) {
      writeBySubsets<ushort3, ushort>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == uint3_type_index) {
      writeBySubsets<uint3, uint>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == ullint3_type_index) {
      writeBySubsets<ullint3, ullint>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == char4_type_index) {
      writeBySubsets<char4, char>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == short4_type_index) {
      writeBySubsets<short4, short>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == int4_type_index) {
      writeBySubsets<int4, int>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == llint4_type_index) {
      writeBySubsets<llint4, llint>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == float4_type_index) {
      writeBySubsets<float4, float>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == double4_type_index) {
      writeBySubsets<double4_16a, double>(netcdf_file_id, filename, cornerstone_indices,
                                          edge_lengths);
    }
    else if (type_code == ushort4_type_index) {
      writeBySubsets<ushort4, ushort>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == uint4_type_index) {
      writeBySubsets<uint4, uint>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
    else if (type_code == ullint4_type_index) {
      writeBySubsets<ullint4, ullint>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
  }
  else {
    const size_t tuple_multiplier = nc_length[dimension_count - 1] / length[dimension_count - 1];
    std::vector<size_t> actual_cornerstone_indices = cornerstone_indices;
    std::vector<size_t> actual_edge_lengths = edge_lengths;
    actual_cornerstone_indices[dimension_count - 1] *= tuple_multiplier;
    actual_edge_lengths[dimension_count - 1] *= tuple_multiplier;
    int chk_code;
    if (subset) {

      // Write a subset of the entire variable to the file.  Because the array has only one
      // dimension, the difference between C-ordered arrays and Fortran-ordered arrays does not
      // apply.
      switch (nc_type_code) {
      case NC_CHAR:
        chk_code = nc_put_vara_schar(netcdf_file_id, nc_var_code,
                                     actual_cornerstone_indices.data(), actual_edge_lengths.data(),
                                     reinterpret_cast<char*>(tip));
        break;
      case NC_SHORT:
        chk_code = nc_put_vara_short(netcdf_file_id, nc_var_code,
                                     actual_cornerstone_indices.data(), actual_edge_lengths.data(),
                                     reinterpret_cast<short*>(tip));
        break;
      case NC_INT:
        chk_code = nc_put_vara_int(netcdf_file_id, nc_var_code, actual_cornerstone_indices.data(),
                                   actual_edge_lengths.data(), reinterpret_cast<int*>(tip));
        break;
      case NC_INT64:
        chk_code = nc_put_vara_int64(netcdf_file_id, nc_var_code,
                                     actual_cornerstone_indices.data(), actual_edge_lengths.data(),
                                     reinterpret_cast<llint*>(tip));
        break;
      case NC_FLOAT:
        chk_code = nc_put_vara_float(netcdf_file_id, nc_var_code,
                                     actual_cornerstone_indices.data(), actual_edge_lengths.data(),
                                     reinterpret_cast<float*>(tip));
        break;
      case NC_DOUBLE:
        chk_code = nc_put_vara_double(netcdf_file_id, nc_var_code,
                                      actual_cornerstone_indices.data(),
                                      actual_edge_lengths.data(), reinterpret_cast<double*>(tip));
        break;
      case NC_UBYTE:
        chk_code = nc_put_vara_ubyte(netcdf_file_id, nc_var_code,
                                     actual_cornerstone_indices.data(), actual_edge_lengths.data(),
                                     reinterpret_cast<uint8_t*>(tip));
        break;
      case NC_USHORT:
        chk_code = nc_put_vara_ushort(netcdf_file_id, nc_var_code,
                                      actual_cornerstone_indices.data(),
                                      actual_edge_lengths.data(), reinterpret_cast<ushort*>(tip));
        break;
      case NC_UINT:
        chk_code = nc_put_vara_uint(netcdf_file_id, nc_var_code, actual_cornerstone_indices.data(),
                                    actual_edge_lengths.data(), reinterpret_cast<uint*>(tip));
        break;
      case NC_UINT64:
        chk_code = nc_put_vara_uint64(netcdf_file_id, nc_var_code,
                                      actual_cornerstone_indices.data(),
                                      actual_edge_lengths.data(), reinterpret_cast<ullint*>(tip));
        break;
      default:
        rtErr("Attempting to write an usupported data type for variable " + key_string +
              " in file " + filename + ".", "EcumenicalArray", "write");
      }
    }
    else {

      // Write the entire variable to the file.
      switch (nc_type_code) {
      case NC_CHAR:
        chk_code = nc_put_var_schar(netcdf_file_id, nc_var_code, reinterpret_cast<char*>(tip));
        break;
      case NC_SHORT:
        chk_code = nc_put_var_short(netcdf_file_id, nc_var_code, reinterpret_cast<short*>(tip));
        break;
      case NC_INT:
        chk_code = nc_put_var_int(netcdf_file_id, nc_var_code, reinterpret_cast<int*>(tip));
        break;
      case NC_INT64:
        chk_code = nc_put_var_int64(netcdf_file_id, nc_var_code, reinterpret_cast<llint*>(tip));
        break;
      case NC_FLOAT:
        chk_code = nc_put_var_float(netcdf_file_id, nc_var_code, reinterpret_cast<float*>(tip));
        break;
      case NC_DOUBLE:
        chk_code = nc_put_var_double(netcdf_file_id, nc_var_code, reinterpret_cast<double*>(tip));
        break;
      case NC_UBYTE:
        chk_code = nc_put_var_ubyte(netcdf_file_id, nc_var_code, reinterpret_cast<uint8_t*>(tip));
        break;
      case NC_USHORT:
        chk_code = nc_put_var_ushort(netcdf_file_id, nc_var_code, reinterpret_cast<ushort*>(tip));
        break;
      case NC_UINT:
        chk_code = nc_put_var_uint(netcdf_file_id, nc_var_code, reinterpret_cast<uint*>(tip));
        break;
      case NC_UINT64:
        chk_code = nc_put_var_uint64(netcdf_file_id, nc_var_code, reinterpret_cast<ullint*>(tip));
        break;
      default:
        rtErr("Attempting to write an usupported data type for variable " + key_string +
              " in file " + filename + ".", "EcumenicalArray", "write");
      }
    }
    validateNCWriteCode(chk_code, key_string, filename, "EcumenicalScalar", "write");
  }
}

//-------------------------------------------------------------------------------------------------
void EcumenicalArray::read(const int netcdf_file_id, const std::string &filename,
                           const std::vector<size_t> &cornerstone_indices,
                           const std::vector<size_t> &edge_lengths) {

  // Check for patch reading commands
  const bool subset = validateSubsetDirective();
  if (dimension_count > 1) {
    if (type_code == char_type_index) {
      writeBySubsets<char, char>(netcdf_file_id, filename, cornerstone_indices, edge_lengths);
    }
  }
  else {
    const size_t tuple_multiplier = nc_length[dimension_count - 1] / length[dimension_count - 1];
    std::vector<size_t> actual_cornerstone_indices = cornerstone_indices;
    std::vector<size_t> actual_edge_lengths = edge_lengths;
    actual_cornerstone_indices[dimension_count - 1] *= tuple_multiplier;
    actual_edge_lengths[dimension_count - 1] *= tuple_multiplier;
    int chk_code;
    if (subset) {

      // Read a subset of the entire variable from the file.
      switch (nc_type_code) {
      case NC_CHAR:
        chk_code = nc_put_vara_schar(netcdf_file_id, nc_var_code,
                                     actual_cornerstone_indices.data(), actual_edge_lengths.data(),
                                     reinterpret_cast<char*>(tip));
        break;
      }
    }
    else {
    }
  }
}

//-------------------------------------------------------------------------------------------------
void EcumenicalArray::validateLength(const size_t length_in) {
  if (length_in > 281474976710656LLU) {
    rtErr("A length of " + std::to_string(length_in) + " is too large.", "EcumenicalArray",
          "validateLength");
  }
}

//-------------------------------------------------------------------------------------------------
void EcumenicalArray::validateSubsetDirective(const std::vector<size_t> &cornerstone_indices,
                                              const std::vector<size_t> &edge_lengths) {
  bool result = false;
  if (cornerstone_indices.size() > 0) {
    if (cornerstone_indices.size() != dimension_count) {
      rtErr("In order to write a subset of the data, the starting index of the subset along each "
            "dimension must be specified.  " + std::to_string(cornerstone_indices.size()) +
            " starting indices were provided for " + std::to_string(dimension_count) + ".",
            "EcumenicalArray", "write");
    }
    if (edge_lengths.size() != dimension_count) {
      rtErr("In order to write a subset of the data, the edge length of the subset along each "
            "dimension must be specified.  " + std::to_string(edge_lengths.size()) +
            " edge lengths were provided for " + std::to_string(dimension_count) + ".",
            "EcumenicalArray", "write");
    }
    result = true;
  }
  if ((cornerstone_indices.size() > 0 && edge_lengths.size() == 0) ||
      (cornerstone_indices.size() == 0 && edge_lengths.size() > 0)) {
    rtErr("In order to write a subset of the data, a set of cornerstone starting indices and a "
          "set of edge lengths must be provided.", "EcumenicalArray", "write");
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
NetCDFLink::NetCDFLink() :
    scalar_contents{}, array_contents{}, filename{}
{}

//-------------------------------------------------------------------------------------------------
NetCDFLink::NetCDFLink(const std::vector<EcumenicalScalar> &scalar_contents_in,
                       const std::vector<EcumenicalArray> &tensor_contents_in,
                       const std::string &filename_in) :
    scalar_contents{scalar_contents_in}, array_contents{array_contents_in}, filename{filename_in}
{}

//-------------------------------------------------------------------------------------------------
int assignNCTypeIndex(size_t stormm_type_code) {

  // As is the case elsewhere, because the index codes of various data types are not defined
  // until runtime, a branched if statement must be used rather than a switch.
  if (stormm_type_code ==  int_type_index || stormm_type_code == int2_type_index ||
      stormm_type_code == int3_type_index || stormm_type_code == int4_type_index) {
    return NC_INT;
  }
  if (stormm_type_code ==  uint_type_index || stormm_type_code == uint2_type_index ||
      stormm_type_code == uint3_type_index || stormm_type_code == uint4_type_index) {
    return NC_UINT;
  }
  else if (stormm_type_code ==     llint_type_index || stormm_type_code == longlong2_type_index ||
           stormm_type_code == longlong3_type_index || stormm_type_code == longlong4_type_index) {
    return NC_INT64;
  }
  else if (stormm_type_code ==     ullint_type_index ||
           stormm_type_code == ulonglong2_type_index ||
           stormm_type_code == ulonglong3_type_index ||
           stormm_type_code == ulonglong4_type_index) {
    return NC_UINT64;
  }
  else if (stormm_type_code ==  short_type_index || stormm_type_code == short2_type_index ||
           stormm_type_code == short3_type_index || stormm_type_code == short4_type_index) {
    return NC_SHORT;
  }
  else if (stormm_type_code ==  ushort_type_index || stormm_type_code == ushort2_type_index ||
           stormm_type_code == ushort3_type_index || stormm_type_code == ushort4_type_index) {
    return NC_USHORT;
  }
  else if (stormm_type_code ==  char_type_index || stormm_type_code == char2_type_index ||
           stormm_type_code == char3_type_index || stormm_type_code == char4_type_index) {
    return NC_CHAR;
  }
  else if (stormm_type_code ==  uchar_type_index || stormm_type_code == uchar2_type_index ||
           stormm_type_code == uchar3_type_index || stormm_type_code == uchar4_type_index) {
    return NC_UBYTE;
  }
  else if (stormm_type_code ==  float_type_index || stormm_type_code == float2_type_index ||
           stormm_type_code == float3_type_index || stormm_type_code == float4_type_index) {
    return NC_FLOAT;
  }
  else if (stormm_type_code ==  double_type_index || stormm_type_code == double2_type_index ||
           stormm_type_code == double3_type_index || stormm_type_code == double4_type_index) {
    return NC_DOUBLE;
  }
  else {
    rtErr("Unrecognized type code " + std::to_string(stormm_type_code) + ".", "assignNCTypeIndex");
  }
}

//-------------------------------------------------------------------------------------------------
void validateNCDefinitionCode(const int chk_code, const std::string &varname,
                              const std::string &filename, const int* array_lengths,
                              const int dimension_count, const char* class_caller,
                              const char* method_caller) {
  switch (chk_code) {
  case NC_NOERR:
    break;
  case NC_EBADID:
    rtErr("The NetCDF file identifier for " + filename + " is invalid. (Attempting to write "
          "variable " + varname + ".)", class_caller, method_caller);
  case NC_EMAXNAME:
    rtErr("The NetCDF dimension name based on variable \"" + varname + "\" is too long.",
          class_caller, method_caller);
  case NC_EBADNAME:
    rtErr("The NetCDF naming conventions do not allow \"" + varname + "\" as a name.  Some "
          "characters may be problematic.", class_caller, method_caller);
  case NC_EINVAL:
    rtErr("Invalid input parameters to file " + filename + ".", class_caller, method_caller);
  case NC_ENOTINDEFINE:
    rtErr("The NetCDF output file " + filename + " is not in define mode.  Meta data can no "
          "longer be written for variable " + varname + ".", class_caller, method_caller);
  case NC_EDIMSIZE:
    {
      std::string all_dims = std::to_string(array_lengths[0]);
      for (int j = 1; j < ndims; j++) {
        all_dims += std::string(", ") + std::to_string(array_lengths[j]);
      }
      rtErr("An invalid dimension size (among " + all_dims + ") was specified while writing "
            "meta data for " + varname + " in file " + filename + ".", class_caller,
            method_caller);
    }
    break;
  case NC_EVARSIZE:
    {
      size_t stormm_est_size = array_lengths[0];
      for (int i = 0; i < dimension_count; i++) {
        stormm_est_size *= array_lengths[i];
      }
      rtErr("The size of variable " + varname + " (STORMM estimate " +
            std::to_string(stormm_est_size) + " bytes) exceeds the constraints of the NetCDF "
            "file format.", class_caller, method_caller);
    }
    break;
  case NC_EFILEMETA:
    rtErr("An error occurred while writing " + varname + " to file " + filename + ".",
          class_caller, method_caller);
  case NC_EUNLIMIT:
    rtErr("The one permitted \"unlimited\" dimension is already in use for the NetCDF output "
          "file, although variable " + varname + " calls for another.", class_caller,
          method_caller);
  case NC_EMAXDIMS:
    rtErr("The maximum number of array dimensions has been exceeded in variable " + varname +
          " when writing " + filename + ".  This should never happen with modern NetCDF.",
          class_caller, method_caller);
  case NC_ENAMEINUSE:
    rtErr("The variable name " + varname + " is already in use while trrying to write file " +
          filename + ".", class_caller, method_caller);
  case NC_ENOMEM:
    rtErr("A memory allocation failure occured in the NetCDF routines while writing variable " +
          varname + " in file " + filename + ".", class_caller, method_caller);
  case NC_EPERM:
    rtErr("The program is attempting to write data from variable " + varname + " into a read-only "
          "file " + filename + ".", class_caller, method_caller);
  default:
    rtErr("Unknown NetCDF error while writing variable " + varname + " to file " + filename + ".",
          class_caller, method_caller);
  }
}

//-------------------------------------------------------------------------------------------------
void validateNCWriteCode(const int chk_code, const std::string &varname,
                         const std::string &filename, const char* class_caller,
                         const char* method_caller) {
  switch (chk_code) {
  case NC_NOERR:
    break;
  case NC_EHDFERR:
    rtErr("An error was reported by the HDF-5 layer when attempting to write variable " + varname +
          " to file " + filename + ".", class_caller, method_caller);
  case NC_EINDEFINE:
    rtErr("The NetCDF file " + filename + " is in define mode while attempting to write " +
          varname + ".", class_caller, method_caller);
  case NC_EBADID:
    rtErr("The specified NetCDF identifier does not refer to an active file.", class_caller,
          method_caller);
  case NC_ENOTVAR:
    rtErr("The variable identifier provided does not correspond to a known variable in the NetCDF "
          "layer.", class_caller, method_caller);
  default:
    rtErr("Unknown NetCDF error while writing file " + filename + ".", class_caller,
          method_caller);
  }
}

//-------------------------------------------------------------------------------------------------
void validateNCReadCode(const int chk_code, const std::string &varname,
                        const std::string &filename, const char* class_caller,
                        const char* method_caller, const int* start_indices, const int* ranges,
                        const int* var_lengths, const int dimension_count) {
  if (chk_code == NC_NOERR) {
    return;
  }

  // Create a message about the variable size in the event that one of the errors will need it.
  std::string lengths_msg;
  if (var_lengths != nullptr) {
    if (dimension_count == 1) {
      lengths_msg = std::string(" (size ") + std::to_string(var_lengths[0]) + std::string(")");
    }
    else {
      lengths_msg = " (size ";
      for (int i = 0; i < dimension_count; i++) {
        lengths_msg += std::to_string(var_lengths[i]);
        if (i < dimension_count - 1) {
          lengths_msg += ", ";
        }
      }
      lengths_msg += ")";
    }
  }

  // Branch over all error codes
  switch (chk_code) {
  case NC_NOERR:
    break;
  case NC_EINVALIDCOORDS:
    if (dimension_count == 1) {
      std::string base_index_msg;
      if (start_indices != nullptr) {
        base_index_msg = std::string("(") + std::to_string(start_indices[0]) + std::string(") ");
      }
      rtErr("The requested base array index " + base_index_msg + "is invalid for variable " +
            varname + lengths_msg + " while reading file " + filename + ".", class_caller,
            method_caller);
    }
    else {
      std::string crnr_index_msg;
      if (start_indices != nullptr) {
        crnr_index_msg = "(";
        for (int i = 0; i < dimension_count; i++) {
          crnr_index_msg += std::to_string(start_indices[i]);
          if (i < dimension_count - 1) {
            crnr_index_msg += ", ";
          }
        }
        crnr_index_msg += ")";
      }
      rtErr("The requested cornerstone array indices " + crnr_index_msg + "are invalid for "
            "variable " + varname + lengths_msg + " while reading file " + filename + ".",
            class_caller, method_caller);
    }
    break;
  case NC_EEDGE:
    {
      if (dimension_count == 1) {
        std::string range_msg, lengths_msg;
        if (start_indices != nullptr && ranges != nullptr) {
          range_msg = std::string("(") + std::to_string(ranges[0]) + " starting at index " +
                      std::to_string(start_indices[0]) + std::string(") ");
        }
      }
      else {
        std::string range_msg, lengths_msg;
        if (start_indices != nullptr && ranges != nullptr) {
          range_msg = "(";
          for (int i = 0; i < dimension_count; i++) {
            range_msg += std::to_string(ranges[i]);
            if (i < dimension_count - 1) {
              range_msg += ", ";
            }
          }
          range_msg += " starting at index ";
          for (int i = 0; i < dimension_count; i++) {
            range_msg += std::to_string(start_indices[i]);
            if (i < dimension_count - 1) {
              range_msg += ", ";
            }
          }
          range_msg += ") ";
        }
      }
      rtErr("The requested range " + range_msg + "would violate the bounds of variable " +
            varname + lengths_msg + " while reading file " + filename + ".", class_caller,
            method_caller);
    }
    break;
  case NC_EINDEFINE:
    rtErr("The NetCDF file " + filename + " is in define mode while attempting to read " +
          varname + ".", class_caller, method_caller);
  case NC_EBADID:
    rtErr("The specified NetCDF identifier does not refer to an active file.", class_caller,
          method_caller);
  case NC_ENOTVAR:
    rtErr("The variable identifier provided does not correspond to a known variable in the NetCDF "
          "layer.", class_caller, method_caller);
  default:
    rtErr("Unknown NetCDF error while reading file " + filename + ".", class_caller,
          method_caller);
  }
}

} // namespace trajectory
} // namespace stormm
#endif // STORMM_USE_NETCDF
