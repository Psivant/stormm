// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
template <typename T>
EcumenicalScalar::EcumenicalScalar(const T value_in, const size_t type_code_in,
                                   const std::string &key_string_in) :
    value{}, value_length{static_cast<int>(sizeof(T))}, type_code{type_code_in},
    key_string{key_string_in}
{
  // Transform the value into a byte string
  const T* ptr = &value_in;
  const uint8_t* ui_ptr = reinterpret_cast<const uint8_t*>(ptr);
  value.resize(value_length);
  for (int i = 0; i < value_length; i++) {
    value[i] = ui_ptr[i];
  }

  // Store the type name of the original value, for error reporting purposes
  if (isHpcVectorType<T>()) {
    type_name = getHpcVectorTypeName<T>();
  }
  else if (isScalarType<T>()) {
    type_name = getStormmScalarTypeName<T>();
  }
  else {
    rtErr("Only known scalar and HPC vector types are permitted.", "EcumenicalScalar");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> T EcumenicalScalar::getValue() const {

  // Test that the requested return type matches the value stored.
  if (std::type_index(typeid(T)) != type_code) {
    std::string return_type_name;
    if (isHpcVectorType<T>()) {
      return_type_name = getHpcVectorTypeName<T>();
    }
    else if (isScalarType<T>()) {
      return_type_name = getStormmScalarTypeName<T>();
    }
    else {
      return_type_name = "unknown";
    }
    rtErr("The value stored is of type " + type_name + ", but the requested return value is of "
          "type " + return_type_name, "EcumenicalScalar", "getValue");
  }
  T result;
  T* res_ptr;
  uint8_t* ui_ptr = reinterpret_cast<uint8_t*>(res_ptr);
  for (int i = 0; i < value_length; i++) {
    ui_ptr[i] = value[i];
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
template <typename T>
T* EcumenicalArray::getPointer(const int second_index, const int third_index,
                               const int fourth_index) {
  T* rc_tip = reinterpret_cast<T*>(tip);
  return &rc_tip[((((fourth_index * length[3]) + third_index) * length[2]) + second_index) *
                 length[1]];
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const T* EcumenicalArray::getPointer(const int second_index, const int third_index,
                                     const int fourth_index) const {
  const T* rc_tip = reinterpret_cast<T*>(tip);
  return &rc_tip[((((fourth_index * length[3]) + third_index) * length[2]) + second_index) *
                 length[1]];
}

//-------------------------------------------------------------------------------------------------
template <typename T>
T EcumenicalArray::readValue(const size_t pos_first, const size_t pos_second,
                             const size_t pos_third, const size_t pos_fourth) const {
  const T* rc_tip = reinterpret_cast<T*>(tip);
  return rc_tip[(((((pos_fourth * length[3]) + pos_third) * length[2]) + pos_second) *
                 length[1]) + pos_first];
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void EcumenicalArray::putValue(const T value, const size_t pos_first, const size_t pos_second,
                               const size_t pos_third, const size_t pos_fourth) {
  T* rc_tip = reinterpret_cast<T*>(tip);
  rc_tip[(((((pos_fourth * length[3]) + pos_third) * length[2]) + pos_second) *
          length[1]) + pos_first] = value;
}

//-------------------------------------------------------------------------------------------------
template <typename Tstormm, typename Tncdf>
void EcumenicalArray::writeSubset(const int netcdf_file_id, const std::string &filename,
                                  const std::vector<size_t> &cornerstone,
                                  const std::vector<size_t> &sbs_lengths, void* f_order_buffer,
                                  void* c_order_buffer) {
  Tstormm* rc_buffer = reinterpret_cast<Tstormm*>(c_order_buffer);
  Tstormm* rf_buffer = reinterpret_cast<Tstormm*>(f_order_buffer);
  const Tstormm* r_tip = reinterpret_cast<Tstormm*>(tip);

  // Stash the appropriate subset of the array for safe keeping
  for (size_t i = cornerstone[0]; i < i_hlim; i++) {
    const int i_del = i - cornerstone[0];
    for (size_t j = cornerstone[1]; j < j_hlim; j++) {
      const int j_del = j - cornerstone[1];
      if (dimension_count > 2) {
        const size_t k_hlim = cornerstone[2] + sbs_lengths[2];
        for (size_t k = cornerstone[2]; k < k_hlim; k++) {
          const int k_del = k - cornerstone[2];
          if (dimension_count == 4) {
            const size_t m_hlim = cornerstone[3] + sbs_lengths[3];
            for (size_t m = cornerstone[3]; m < m_hlim; m++) {
              const int m_del = m - cornerstone[3];

              // Both buffers remain in Fortran order, like other STORMM arrays.
              const size_t buffer_idx = (((((m_del * sbs_lengths[2]) + k_del) *
                                           sbs_lengths[1]) + j_del) * sbs_lengths[0]) + i_del;

              // Compute the Fortran order indices in the array.  These are indices that will send
              // values to other parts of the array.
              const size_t f_gbl_idx = (((((m * lengths[2]) + k) * lengths[1]) + j) *
                                        lengths[0]) + i; 
              rf_buffer[buffer_idx] = r_tip[f_gbl_idx];
              
	      // Compute the C order index in the array.  These are indices that might take on
              // new values. 
              const size_t c_gbl_idx = (((((i * lengths[1]) + j) * lengths[2]) + k) *
                                        lengths[3]) + m;
              rc_buffer[buffer_idx] = r_tip[c_gbl_idx];
            }
          }
          else {

            // For three dimensions, the array contains slabs, rows, and columns, with the index
            // of array arr for either C or Fortran reading arr[slab_id][row_id][column_id].  In
            // C, the column count increments fastest, slab count increments slowest.  In Fortran,
            // the slab count increments fastest, while the column count increments slowest.
            const size_t buffer_idx = (((k_del * sbs_lengths[1]) + j_del) * sbs_lengths[0]) +
                                      i_del;
            const size_t f_gbl_idx = (((k * lengths[1]) + j) * lengths[0]) + i;
            rf_buffer[buffer_idx] = r_tip[f_gbl_idx];
            const size_t c_gbl_idx = (((i * lengths[1]) + j) * lengths[2]) + k;
            rc_buffer[buffer_idx] = r_tip[c_gbl_idx];
          }
        }
      }
      else {
        const size_t buffer_idx = (j_del * sbs_lengths[0]) + i_del;
        const size_t f_gbl_idx = (j * sbs_lengths[0]) + i;
        const size_t c_gbl_idx = (i * sbs_lengths[1]) + j;
        rc_buffer[buffer_idx] = rc_tip[global_idx];
      }
    }
  }

  // Move the Fortran-ordered data into the C-ordered indices.  In order to ensure that the value
  // of (i, j, k, m) in the Fortran-ordered array has not been corrupted by some prior transfer,
  // take the copy stashed in the Fortran buffer.
  for (size_t i = cornerstone[0]; i < i_hlim; i++) {
    const int i_del = i - cornerstone[0];
    for (size_t j = cornerstone[1]; j < j_hlim; j++) {
      const int j_del = j - cornerstone[1];
      if (dimension_count > 2) {
        const size_t k_hlim = cornerstone[2] + sbs_lengths[2];
        for (size_t k = cornerstone[2]; k < k_hlim; k++) {
          const int k_del = k - cornerstone[2];
          if (dimension_count == 4) {
            const size_t m_hlim = cornerstone[3] + sbs_lengths[3];
            for (size_t m = cornerstone[3]; m < m_hlim; m++) {
              const int m_del = m - cornerstone[3];

              // The index i, j, k, m of the Fortran-ordered array (STORMM's native format, and
              // the format in which most of the array will probably remain) will go into index
              // i, j, k, m of the C-ordered array (NetCDF's native format).  The indices of the
              // strung-out, one-dimensional array are what differ.  In order to ensure that the
              // value of (i, j, k, m) in the Fortran-ordered array has not been corrupted by some
              // prior transfer, take the copy stashed in the Fortran buffer.
              const size_t buffer_idx = (((((m_del * sbs_lengths[2]) + k_del) *
                                           sbs_lengths[1]) + j_del) * sbs_lengths[0]) + i_del;
              const size_t c_gbl_idx = (((((i * lengths[1]) + j) * lengths[2]) + k) *
                                        lengths[3]) + m;
              r_tip[c_gbl_idx] = rf_buffer[buffer_idx];
            }
          }
          else {
            const size_t buffer_idx = (((k_del * sbs_lengths[1]) + j_del) * sbs_lengths[0]) +
                                      i_del;
            const size_t c_gbl_idx = (((i * lengths[1]) + j) * lengths[2]) + k;
            r_tip[c_gbl_idx] = rf_buffer[buffer_idx];
          }
        }
      }
      else {
        const size_t buffer_idx = (j_del * sbs_lengths[0]) + i_del;
        const size_t c_gbl_idx = (i * sbs_lengths[1]) + j;
        r_tip[c_gbl_idx] = rf_buffer[buffer_idx];
      }
    }
  }

  // Write the NetCDF file.  This is the point at which to translate between two-, three-, and
  // four-tuples and the corresponding atomic type stored by the NetCDF file.
  const size_t tuple_multiplier = sizeof(Tstormm) / sizeof(Tncdf);
  std::vector<size_t> nc_sbs_lengths(dimension_count), nc_cornerstone_indices(dimension_count);
  for (int i = 0; i < dimension_count; i++) {
    nc_sbs_lengths[i] = sbs_lengths[i] * tuple_multiplier;
    nc_cornerstone_indices[i] = cornerstone_indices[i];
  }
  nc_cornerstone_indices[dimension_count - 1] *= tuple_multiplier;
  int chk_code;
  switch (nc_type_code) {
  case NC_CHAR:
    chk_code = nc_put_vara_schar(netcdf_file_id, nc_var_code, nc_cornerstone_indices.data(),
                                 nc_sbs_lengths.data(), reinterpret_cast<char*>(tip));
    break;
  case NC_SHORT:
    chk_code = nc_put_vara_short(netcdf_file_id, nc_var_code, nc_cornerstone_indices.data(),
                                 nc_sbs_lengths.data(), reinterpret_cast<short*>(tip));
    break;
  case NC_INT:
    chk_code = nc_put_vara_int(netcdf_file_id, nc_var_code, nc_cornerstone_indices.data(),
                               nc_sbs_lengths.data(), reinterpret_cast<int*>(tip));
    break;
  case NC_INT64:
    chk_code = nc_put_vara_int64(netcdf_file_id, nc_var_code, nc_cornerstone_indices.data(),
                                 nc_sbs_lengths.data(), reinterpret_cast<llint*>(tip));
    break;
  case NC_FLOAT:
    chk_code = nc_put_vara_float(netcdf_file_id, nc_var_code, nc_cornerstone_indices.data(),
                                 nc_sbs_lengths.data(), reinterpret_cast<float*>(tip));
    break;
  case NC_DOUBLE:
    chk_code = nc_put_vara_double(netcdf_file_id, nc_var_code, nc_cornerstone_indices.data(),
                                  nc_sbs_lengths.data(), reinterpret_cast<double*>(tip));
    break;
  case NC_UBYTE:
    chk_code = nc_put_vara_ubyte(netcdf_file_id, nc_var_code, nc_cornerstone_indices.data(),
                                 nc_sbs_lengths.data(), reinterpret_cast<uint8_t*>(tip));
    break;
  case NC_USHORT:
    chk_code = nc_put_vara_ushort(netcdf_file_id, nc_var_code, nc_cornerstone_indices.data(),
                                  nc_sbs_lengths.data(), reinterpret_cast<ushort*>(tip));
    break;
  case NC_UINT:
    chk_code = nc_put_vara_uint(netcdf_file_id, nc_var_code, nc_cornerstone_indices.data(),
                                nc_sbs_lengths.data(), reinterpret_cast<uint*>(tip));
    break;
  case NC_UINT64:
    chk_code = nc_put_vara_uint64(netcdf_file_id, nc_var_code, nc_cornerstone_indices.data(),
                                  nc_sbs_lengths.data(), reinterpret_cast<ullint*>(tip));
    break;
  default:
    rtErr("Attempting to write an unsupported data type for variable " + key_string +
          " in file " + filename + ".", "EcumenicalArray", "writeSubset");
  }
  validateNCWriteCode(chk_code, key_string, filename, "EcumenicalScalar", "writeSubset");
  
  // Replace the data in the C-ordered indices so that the array, once again, holds its original
  // Fortran-ordered content.
  for (size_t i = cornerstone[0]; i < i_hlim; i++) {
    const int i_del = i - cornerstone[0];
    for (size_t j = cornerstone[1]; j < j_hlim; j++) {
      const int j_del = j - cornerstone[1];
      if (dimension_count > 2) {
        const size_t k_hlim = cornerstone[2] + sbs_lengths[2];
        for (size_t k = cornerstone[2]; k < k_hlim; k++) {
          const int k_del = k - cornerstone[2];
          if (dimension_count == 4) {
            const size_t m_hlim = cornerstone[3] + sbs_lengths[3];
            for (size_t m = cornerstone[3]; m < m_hlim; m++) {
              const int m_del = m - cornerstone[3];
              const size_t buffer_idx = (((((m_del * sbs_lengths[2]) + k_del) *
                                           sbs_lengths[1]) + j_del) * sbs_lengths[0]) + i_del;
              r_tip[c_gbl_idx] = rc_buffer[buffer_idx];
            }
          }
          else {
            const size_t buffer_idx = (((k_del * sbs_lengths[1]) + j_del) * sbs_lengths[0]) +
                                      i_del;
            const size_t c_gbl_idx = (((i * lengths[1]) + j) * lengths[2]) + k;
            r_tip[c_gbl_idx] = rc_buffer[buffer_idx];
          }
        }
      }
      else {
        const size_t buffer_idx = (j_del * sbs_lengths[0]) + i_del;
        const size_t c_gbl_idx = (i * sbs_lengths[1]) + j;
        r_tip[c_gbl_idx] = rc_buffer[buffer_idx];
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tstormm, typename Tncdf>
void EcumenicalArray::writeBySubsets(int netcdf_file_id, const std::string &filename,
                                     const std::vector<size_t> &cornerstone_indices,
                                     const std::vector<size_t> &edge_lengths) {

  // Allocate a pair of buffers for re-arranging this array while writing chunks of it to the disk.
  // The first array will store the values of indices in the Fortran order which will send
  // their contents to other indices in order to make the array look as if it is in C order.  The
  // second array will store the values which would be overwritten by this process.  The two arrays
  // may share some indices in common, but both are necessary.
  std::vector<Tstormm> fortran_buffer(io_batch_size), c_buffer(io_batch_size);
  
  // Write pieces of the variable to the file, one at a time, until the entirety of the variable
  // or requested subset is written.
  std::vector<size_t> tmp_crnr(4, 0), tmp_apex(4, 0), local_crnr(4, 0), local_lengths(4, 0);
  for (int i = 0; i < ndims; i++) {
    tmp_crnr[i] = (subset) ? cornerstone_indices[i] : 0;
    tmp_apex[i]  = (subset) ? cornerstone_indices[i] + edge_lengths[i] : length[i];
  }
  for (size_t i_llim = tmp_crnr[0]; i_llim < tmp_apex[0]; i_llim += io_length[0]) {
    local_lengths[0] = (i_llim + io_length[0] < tmp_apex[0]) ? io_length[0] : tmp_apex[0] - i_llim;
    for (size_t j_llim = tmp_crnr[1]; j_llim < tmp_apex[1]; j_llim += io_length[1]) {
      local_lengths[1] = (j_llim + io_length[1] < tmp_apex[1]) ? io_length[1] :
                                                                 tmp_apex[1] - j_llim;
      if (ndims > 2) {
        for (size_t k_llim = tmp_crnr[2]; k_llim < tmp_apex[2]; k_llim += io_length[2]) {
         local_lengths[2] = (k_llim + io_length[2] < tmp_apex[2]) ? io_length[2] :
                                                                    tmp_apex[2] - k_llim;
          if (ndims > 3) {
            for (size_t m_llim = tmp_crnr[3]; m_llim < tmp_apex[3]; m_llim += io_length[3]) {
              local_lengths[3] = (m_llim + io_length[3] < tmp_apex[3]) ? io_length[3] :
                                                                         tmp_apex[3] - m_llim;
              writeSubset<Tstormm, Tncdf>(netcdf_file_id, filename, local_crnr, local_lengths,
                                          fortran_buffer.data(), c_buffer.data());
            }
          }
          else {
            writeSubset<Tstormm, Tncdf>(netcdf_file_id, filename, local_crnr, local_lengths,
                                        fortran_buffer.data(), c_buffer.data());
          }
        }
        else {
          writeSubset<Tstormm, Tncdf>(netcdf_file_id, filename, local_crnr, local_lengths,
                                      fortran_buffer.data(), c_buffer.data());

        }
      }
      else {
        writeSubset<Tstormm, Tncdf>(netcdf_file_id, filename, local_crnr, local_lengths,
                                    fortran_buffer.data(), c_buffer.data());
      }
    }	
  }
}

} // namespace trajectory
} // namespace stormm
