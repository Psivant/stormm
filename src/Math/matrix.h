// -*-c++-*-
#ifndef STORMM_MATRIX_H
#define STORMM_MATRIX_H

#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Random/random.h"

namespace stormm {
namespace stmath {

using card::default_hpc_format;
using card::Hybrid;
using card::HybridFormat;
using card::HybridKind;
using card::HybridLabel;
using card::HybridTargetLevel;
using random::Ran2Generator;

enum class MatrixFillMode {
  ZEROS,    ///< Fill with all values set to zero (equivalent to not providing a fill method)
  ONES,     ///< Fill with all values set to 1.0
  EYE,      ///< Fill with the identity matrix
  RANDU,    ///< Fill with random numbers from a uniform distribution in the range [0, 1)
  RANDN,    ///< Fill with random numbers from a normal distribution with sigma 1.0 and mean 0.0
  VALUE,    ///< Fill with a specified value (given in a subsequent argument)
};

/// \brief A matrix object that can traverse CPU and GPU memory
template <typename T> class HpcMatrix {
public:

  // Public-facing, read-only forms of the private variables (allows syntactically simpler access
  // than getter functions)
  const size_t& n_rows;
  const size_t& n_cols;
  const size_t& n_elem;

  /// \brief Constructors must set public-facing reference variables in addition to private member
  ///        variables.
  ///
  /// Overloaded:
  ///   - Initialize an empty matrix with zero rows and zero columns
  ///   - Initialize a matrix with a given numer of rows and columns, all elements set to zero
  ///
  /// \param rows_in      The number of rows to allocate
  /// \param cols_in      The number of columns to allocate
  /// \param fill_method  Manner in which to fill the indices
  /// \param fill_value   A scalar value for initializing all elements if fill_method is VALUE, or
  ///                     the random seed for initializing a random number generator if fill_method
  ///                     is RANDU or RANDN
  /// \{
  HpcMatrix(const char* tag_in = nullptr, HybridFormat format_in = default_hpc_format);
  HpcMatrix(size_t rows_in, size_t cols_in, const char* tag_in = nullptr,
            HybridFormat format_in = default_hpc_format,
            MatrixFillMode fill_method = MatrixFillMode::VALUE, T fill_value = static_cast<T>(0));
  HpcMatrix(size_t rows_in, size_t cols_in, MatrixFillMode fill_method,
            T fill_value = static_cast<T>(0), const char* tag_in = nullptr,
            HybridFormat format_in = default_hpc_format);
  HpcMatrix(size_t rows_in, size_t cols_in, MatrixFillMode fill_method,
            Ran2Generator *matrng, const char* tag_in = nullptr,
            HybridFormat format_in = default_hpc_format);
  /// \}

  /// \brief Produce the label for this matrix (obtained from the label of its primary contents)
  HybridLabel getLabel() const;

  /// \brief Change the number of rows and columns.  Elements outside the bounds of the new matrix
  ///        will be discarded.  If the new matrix is larger than the old one, all elements will be
  ///        preserved and new elements will be initialized to zero.
  ///
  /// \param new_rows  New number of rows to allocate
  /// \param new_cols  New number of columns to allocate
  void resize(size_t new_rows, size_t new_cols);

  /// \brief Empty the matrix and set its size to 0 by 0.
  void reset();

  /// \brief Pointer to the memory itself, on the host or the device.  This takes its name from
  ///        the Armadillo vec, mat, and cube classes.
  ///
  /// Overloaded:
  ///   - Get a const pointer if the object is const
  ///   - Get a pointer that can modify the underlying data if the object is non-const
  ///
  /// \param tier  Indicator of whether to offer a pointer to host or GPU device memory
  /// \{
  const T* memptr(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  T* memptr(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}

  /// \brief Element-wise access via the function call () operator.  This entails two bounds
  ///        checks.
  ///
  /// \param row_idx   Index of the row to access, starting at 0
  /// \param col_idx   Index of the column to access, starting at 0
  /// \param tier      Indicator of whether to access data on the host or on the GPU
  T operator()(size_t row_idx, size_t col_idx, HybridTargetLevel tier = HybridTargetLevel::HOST);

  /// \brief Element-wise access to host data without bounds checks, following on the Armadillo
  ///        naming conventions.
  ///
  /// \param row_idx   Index of the row to access, starting at 0
  /// \param col_idx   Index of the column to access, starting at 0
  T atHost(size_t row_idx, size_t col_idx);

  /// \brief Element-wise access to data on the GPU without bounds checks, following on the
  ///        Armadillo naming conventions.
  ///
  /// \param row_idx  Index of the row to access, starting at 0
  /// \param col_idx  Index of the column to access, starting at 0
  T atDevice(size_t row_idx, size_t col_idx);

  /// \brief Get a Hybrid POINTER-kind object that is set directly to the column of interest.
  ///
  /// \param col_idx  The column of interest
  /// \param tag_in   A note to help identify the column pointer
  Hybrid<T> colptr(size_t col_idx, const char* tag_in = nullptr);

#ifdef STORMM_USE_HPC
  /// \brief Upload the matrix to the device
  void upload();

  /// \brief Download the matrix from the device
  void download();
#endif

private:
  size_t n_rows_pr;           ///< Number of rows
  size_t n_cols_pr;           ///< Number of columns
  size_t n_elem_pr;           ///< Total number of elements
  Hybrid<T> contents;         ///< The main, contiguous memory array for column-major data
#ifdef STORMM_USE_HPC
  Hybrid<T> transfer_buffer;  ///< Buffer with space allocated equal to the size of one column
                              ///<   for sparse or scattered data exchange between the host and
                              ///<   the device.
  T* devc_data;               ///< Additional pointer maintained by this object to track its
                              ///<   underlying Hybrid array contents
#endif
  T* host_data;               ///< Additional pointer maintained by this object to track its
                              ///<   underlying Hybrid array contents
};

} // namespace stmath
} // namespace stormm

#include "matrix.tpp"

#endif
