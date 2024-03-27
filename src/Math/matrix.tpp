// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
template <typename T>
HpcMatrix<T>::HpcMatrix(const char* tag_in, const HybridFormat format_in) :
    n_rows{n_rows_pr},
    n_cols{n_cols_pr},
    n_elem{n_elem_pr},
    n_rows_pr{0},
    n_cols_pr{0},
    n_elem_pr{0},
    contents{0, tag_in},
#ifdef STORMM_USE_HPC
    transfer_buffer{0, "matrix_buffer", HybridFormat::EXPEDITED},
    devc_data{contents.data(HybridTargetLevel::DEVICE)},
#endif
    host_data{contents.data(HybridTargetLevel::HOST)}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
HpcMatrix<T>::HpcMatrix(const size_t rows_in , const size_t cols_in, const char* tag_in,
                        const HybridFormat format_in, const MatrixFillMode fill_method,
                        const T fill_value) :
    n_rows{n_rows_pr},
    n_cols{n_cols_pr},
    n_elem{n_elem_pr},
    n_rows_pr{rows_in},
    n_cols_pr{cols_in},
    n_elem_pr{rows_in * cols_in},
    contents{n_elem_pr, tag_in},
#ifdef STORMM_USE_HPC
    transfer_buffer{rows_in, "matrix_buffer", HybridFormat::EXPEDITED},
    devc_data{contents.data(HybridTargetLevel::DEVICE)},
#endif
    host_data{contents.data(HybridTargetLevel::HOST)}
{
  switch (fill_method) {
  case MatrixFillMode::ZEROS:
    break;
  case MatrixFillMode::ONES:
    {
      const T tval = static_cast<T>(1);
      for (size_t i = 0; i < n_elem_pr; i++) {
        host_data[i] = tval;
      }
    }
    break;
  case MatrixFillMode::EYE:
    {
      const T tval = static_cast<T>(1);
      for (size_t i = 0; i < n_rows_pr; i++) {
        if (i < n_cols_pr) {
          host_data[(i * n_rows_pr) + i] = tval;
        }
      }
    }
    break;
  case MatrixFillMode::RANDU:
  case MatrixFillMode::RANDN:
    rtErr("Matrix intialization to random values requires a pointer to a random number "
          "generator.  Use an alternative constructor.", "HpcMatrix");
  case MatrixFillMode::VALUE:
    for (size_t i = 0; i < n_elem_pr; i++) {
      host_data[i] = fill_value;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
HpcMatrix<T>::HpcMatrix(const size_t rows_in , const size_t cols_in,
                        const MatrixFillMode fill_method, const T fill_value, const char* tag_in,
                        const HybridFormat format_in) :
  HpcMatrix(rows_in ,cols_in, tag_in, format_in, fill_method, fill_value)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
HpcMatrix<T>::HpcMatrix(const size_t rows_in , const size_t cols_in,
                        const MatrixFillMode fill_method, Ran2Generator *matrng,
                        const char* tag_in, const HybridFormat format_in) :
    n_rows{n_rows_pr},
    n_cols{n_cols_pr},
    n_elem{n_elem_pr},
    n_rows_pr{rows_in},
    n_cols_pr{cols_in},
    n_elem_pr{rows_in * cols_in},
    contents{n_elem_pr, tag_in},
#ifdef STORMM_USE_HPC
    transfer_buffer{rows_in, "matrix_buffer", HybridFormat::EXPEDITED},
    devc_data{contents.data(HybridTargetLevel::DEVICE)},
#endif
    host_data{contents.data(HybridTargetLevel::HOST)}
{
  switch (fill_method) {
  case MatrixFillMode::ZEROS:
  case MatrixFillMode::ONES:
  case MatrixFillMode::EYE:
  case MatrixFillMode::VALUE:
    rtErr("Matrix intialization with a random number generator must proceed in MatrixFillMode "
          "RANDU or RANDN.", "HpcMatrix");
  case MatrixFillMode::RANDU:
    for (size_t i = 0; i < n_elem_pr; i++) {
      host_data[i] = matrng->uniformRandomNumber();
    }
    break;
  case MatrixFillMode::RANDN:
    for (size_t i = 0; i < n_elem_pr; i++) {
      host_data[i] = matrng->gaussianRandomNumber();
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
HybridLabel HpcMatrix<T>::getLabel() const {
  return contents.getLabel();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void HpcMatrix<T>::resize(const size_t new_rows, const size_t new_cols) {
  const size_t new_elem = new_rows * new_cols;
  if (new_rows < n_rows_pr) {
    const size_t relevant_elem = std::min(new_cols, n_cols_pr) * std::min(new_rows, n_rows_pr);
    for (size_t counter = 0LLU; counter < relevant_elem; counter++) {
      const size_t col_idx = counter / new_rows;
      const size_t row_idx = counter - (col_idx * new_rows);
      contents[counter] = contents[(col_idx * n_rows_pr) + row_idx];
    }
  }
  else if (new_rows > n_rows_pr) {

    // Lengthen each column to accommodate the new number of rows
    contents.resize(new_rows * n_cols_pr);
    for (size_t counter = n_elem; counter <= n_elem; counter--) {
      const size_t col_idx = counter / n_rows_pr;
      const size_t row_idx = counter - (col_idx * n_rows_pr);
      contents[(col_idx * new_rows) + row_idx] = contents[counter];
    }

    // Fill in the new space with zeros
    for (size_t col_idx = 0LLU; col_idx < n_cols_pr; col_idx++) {
      for (size_t row_idx = n_rows_pr; row_idx < new_rows; row_idx++) {
        contents[(col_idx * new_rows) + row_idx] = 0.0;
      }
    }
  }

  // Reallocate the contents--if more space is allocated, it will be filled with zeros.
  contents.reallocate(new_elem);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void HpcMatrix<T>::reset() {
  n_rows_pr = 0LLU;
  n_cols_pr = 0LLU;
  n_elem_pr = 0LLU;
  contents.reallocate(0);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const T* HpcMatrix<T>::memptr(const HybridTargetLevel tier) const {
  return contents.data(tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
T* HpcMatrix<T>::memptr(const HybridTargetLevel tier) {
  return contents.data(tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
T HpcMatrix<T>::operator()(const size_t row_idx, const size_t col_idx,
                           const HybridTargetLevel tier) {
  if (row_idx >= n_rows_pr) {
    rtErr("Request for row " + std::to_string(row_idx) + " is out of bounds for a matrix of "
          "rank " + std::to_string(n_rows_pr) + ".", "HpcMatrix");
  }
  if (col_idx >= n_cols_pr) {
    rtErr("Request for column " + std::to_string(col_idx) + " is out of bounds for a matrix of " +
          std::to_string(n_cols_pr) + " columns.", "HpcMatrix");
  }
#ifdef STORMM_USE_HPC
  switch (tier) {
  case HybridTargetLevel::HOST:
    return contents.readHost((col_idx * n_rows_pr) + row_idx);
  case HybridTargetLevel::DEVICE:
    return contents.readDevice((col_idx * n_rows_pr) + row_idx);
  }
#else
  return contents.readHost((col_idx * n_rows_pr) + row_idx);
#endif
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
T HpcMatrix<T>::atHost(const size_t row_idx, const size_t col_idx) {
  return host_data[(col_idx * n_rows_pr) + row_idx];
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
template <typename T>
T HpcMatrix<T>::atDevice(const size_t row_idx, const size_t col_idx) {
  return contents.readDevice((col_idx * n_rows_pr) + row_idx);
}
#endif

//-------------------------------------------------------------------------------------------------
template <typename T>
Hybrid<T> HpcMatrix<T>::colptr(const size_t col_idx, const char* tag_in) {
  if (col_idx >= n_cols_pr) {
    rtErr("Matrix has only " + std::to_string(n_cols_pr) + ", cannot access column " +
          std::to_string(col_idx) + ".", "HpcMatrix", "colptr");
  }
  Hybrid<T> hptr(HybridKind::POINTER, tag_in, contents.getFormat());
  hptr.setPointer(&contents, n_rows_pr * col_idx, n_rows_pr);
  return hptr;
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
template <typename T>
void HpcMatrix<T>::upload() {
  contents.upload();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void HpcMatrix<T>::download() {
  contents.download();
}
#endif

} // namespace stmath
} // namespace stormm
