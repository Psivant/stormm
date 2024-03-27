// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
template <typename T>
void matrixMultiply(const T* a, const size_t row_a, const size_t col_a, const T* b,
                    const size_t row_b, const size_t col_b, T* c, const T scale_a,
                    const T scale_b, const T scale_c, const TransposeState x_a,
                    const TransposeState x_b) {

  // Check dimensions
  size_t relevant_row_a, relevant_col_a, relevant_row_b, relevant_col_b;
  switch (x_a) {
  case TransposeState::AS_IS:
    relevant_row_a = row_a;
    relevant_col_a = col_a;
    break;
  case TransposeState::TRANSPOSE:
    relevant_row_a = col_a;
    relevant_col_a = row_a;
    break;
  }
  switch (x_b) {
  case TransposeState::AS_IS:
    relevant_row_b = row_b;
    relevant_col_b = col_b;
    break;
  case TransposeState::TRANSPOSE:
    relevant_row_b = col_b;
    relevant_col_b = row_b;
    break;
  }
  if (relevant_col_a != relevant_row_b) {
    rtErr("The effective number of columns in matrix A (" +
          std::to_string(relevant_col_a) + ") must equal the effective number of rows in matrix "
          "B (" + std::to_string(relevant_row_b) + ").", "matrixMultiply");
  }
  const size_t tile_size = ((row_a >= 128 && col_a >= 128) ||
                            (row_b >= 128 && col_b >= 128)) ? 128 : 
                           ((row_a >= 32 && col_a >= 32) ||
                            (row_b >= 32 && col_b >= 32)) ? 32 : 8;
  std::vector<T> buffer_a(tile_size * tile_size, 0.0);
  std::vector<T> buffer_b(tile_size * tile_size, 0.0);
  const size_t nt_row_a = (relevant_row_a + tile_size - 1) / tile_size;
  const size_t nt_col_a = (relevant_col_a + tile_size - 1) / tile_size;
  const size_t nt_row_b = (relevant_row_b + tile_size - 1) / tile_size;
  const size_t nt_col_b = (relevant_col_b + tile_size - 1) / tile_size;
  
  // Perform the multiplication.  Step through tiles of each matrix, first with each tile row
  // of matrix A and pairing it with the correct tile from the corresponding tile column of
  // matrix B.
  const double scale_ab = scale_a * scale_b;
  for (size_t i = 0; i < relevant_row_a * relevant_col_b; i++) {
    c[i] *= scale_c;
  }
  for (size_t i = 0; i < nt_row_a; i++) {
    const size_t row_llim_a = i * tile_size;
    const size_t row_hlim_a = std::min((i + 1) * tile_size, relevant_row_a);
    const size_t ii_hlim = row_hlim_a - row_llim_a;
    for (size_t j = 0; j < nt_col_b; j++) {
      const size_t col_llim_b = j * tile_size;
      const size_t col_hlim_b = std::min((j + 1) * tile_size, relevant_col_b);

      // nt_col_a == nt_row_b by an implication of the check above.  This is analogous to what
      // will happen in the innermost loops.  We are moving up the columns of A while moving
      // down the rows of B.  Likewise, row_llim_b == col_llim_a and row_hlim_b == col_hlim_a.
      for (size_t k = 0; k < nt_col_a; k++) {
        const size_t col_llim_a = k * tile_size;
        const size_t col_hlim_a = std::min((k + 1) * tile_size, relevant_col_a);
        const size_t kk_hlim = col_hlim_a - col_llim_a;
        switch (x_a) {
        case TransposeState::AS_IS:

          // Buffer matrix A, transposed
          for (size_t ii = row_llim_a; ii < row_hlim_a; ii++) {
            const size_t iits = (ii - row_llim_a) * tile_size;
            for (size_t jj = col_llim_a; jj < col_hlim_a; jj++) {
              buffer_a[iits + jj - col_llim_a] = a[(jj * row_a) + ii];
            }
          }
          switch (x_b) {
          case TransposeState::AS_IS:

            // Buffer matrix B, untransposed
            for (size_t ii = col_llim_b; ii < col_hlim_b; ii++) {
              const size_t iits = (ii - col_llim_b) * tile_size;
              for (size_t jj = col_llim_a; jj < col_hlim_a; jj++) {
                buffer_b[iits + jj - col_llim_a] = b[(ii * row_b) + jj];
              }
            }
            for (size_t ii = 0; ii < ii_hlim; ii++) {
              const size_t iits = ii * tile_size;
              for (size_t jj = col_llim_b; jj < col_hlim_b; jj++) {
                double dsum = 0.0;
                const size_t jjts = (jj - col_llim_b) * tile_size;
                for (size_t kk = 0; kk < kk_hlim; kk++) {
                  dsum += buffer_a[iits + kk] * buffer_b[jjts + kk];
                }
                c[(jj * row_a) + ii + row_llim_a] += dsum * scale_ab;
              }
            }
            break;
          case TransposeState::TRANSPOSE:

            // Buffer matrix B, transposed
            for (size_t ii = col_llim_b; ii < col_hlim_b; ii++) {
              const size_t iits = (ii - col_llim_b) * tile_size;
              for (size_t jj = col_llim_a; jj < col_hlim_a; jj++) {
                buffer_b[iits + jj - col_llim_a] = b[(jj * row_b) + ii];
              }
            }
            for (size_t ii = 0; ii < ii_hlim; ii++) {
              const size_t iits = ii * tile_size;
              for (size_t jj = col_llim_b; jj < col_hlim_b; jj++) {
                const size_t jjts = (jj - col_llim_b) * tile_size;
                double dsum = 0.0;
                for (size_t kk = 0; kk < kk_hlim; kk++) {
                  dsum += buffer_a[iits + kk] * buffer_b[jjts + kk];
                }
                c[(jj * row_a) + ii + row_llim_a] += dsum * scale_ab;
              }
            }
            break;
          }
          break;
        case TransposeState::TRANSPOSE:
          switch (x_b) {
          case TransposeState::AS_IS:

            // Buffer matrix B, untransposed
            for (size_t ii = col_llim_b; ii < col_hlim_b; ii++) {
              const size_t iits = (ii - col_llim_b) * tile_size;
              for (size_t jj = col_llim_a; jj < col_hlim_a; jj++) {
                buffer_b[iits + jj - col_llim_a] = b[(ii * row_b) + jj];
              }
            }
            for (size_t ii = row_llim_a; ii < row_hlim_a; ii++) {
              const size_t ii_rowa = ii * row_a;
              for (size_t jj = col_llim_b; jj < col_hlim_b; jj++) {
                double dsum = 0.0;
                const size_t jjts = (jj - col_llim_b) * tile_size;
                for (size_t kk = col_llim_a; kk < col_hlim_a; kk++) {
                  dsum += a[ii_rowa + kk] * buffer_b[jjts + kk - col_llim_a];
                }
                c[(jj * col_a) + ii] += dsum * scale_ab;
              }
            }
            break;
          case TransposeState::TRANSPOSE:

            // Buffer matrix B, transposed
            for (size_t ii = col_llim_b; ii < col_hlim_b; ii++) {
              const size_t iits = (ii - col_llim_b) * tile_size;
              for (size_t jj = col_llim_a; jj < col_hlim_a; jj++) {
                buffer_b[iits + jj - col_llim_a] = b[(jj * row_b) + ii];
              }
            }
            for (size_t ii = row_llim_a; ii < row_hlim_a; ii++) {
              const size_t ii_rowa = ii * row_a;
              for (size_t jj = col_llim_b; jj < col_hlim_b; jj++) {
                const size_t jjts = (jj - col_llim_b) * tile_size;
                double dsum = 0.0;
                for (size_t kk = col_llim_a; kk < col_hlim_a; kk++) {
                  dsum += a[ii_rowa + kk] * buffer_b[jjts + kk - col_llim_a];
                }
                c[(jj * col_a) + ii] += dsum * scale_ab;
              }
            }
            break;
          }
          break;
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void matrixMultiply(const std::vector<T> &a, const size_t row_a, const size_t col_a,
                    const std::vector<T> &b, const size_t row_b, const size_t col_b,
                    std::vector<T> *c, const T scale_a, const T scale_b, const T scale_c,
                    const TransposeState x_a, const TransposeState x_b) {
  if (a.size() != row_a * col_a) {
    rtErr("The size of the input matrix A (" + std::to_string(a.size()) + ") must comport with "
          "its stated row and column counts (" + std::to_string(row_a) + ", " +
          std::to_string(col_a) + ").", "matrixMultiply");
  }
  if (b.size() != row_b * col_b) {
    rtErr("The size of the input matrix B (" + std::to_string(b.size()) + ") must comport with "
          "its stated row and column counts (" + std::to_string(row_b) + ", " +
          std::to_string(col_b) + ").", "matrixMultiply");
  }
  const size_t relevant_row_c = (x_a == TransposeState::AS_IS) ? row_a : col_a;
  const size_t relevant_col_c = (x_b == TransposeState::AS_IS) ? col_b : row_b;
  if (c->size() != relevant_row_c * relevant_col_c) {
    rtErr("The size of the output matrix C (" + std::to_string(c->size()) + ") must comport with "
          "its stated row and column counts (" + std::to_string(relevant_row_c) + ", " +
          std::to_string(relevant_col_c) + ").", "matrixMultiply");
  }
  matrixMultiply(a.data(), row_a, col_a, b.data(), row_b, col_b, c->data(), scale_a, scale_b,
                 scale_c, x_a, x_b);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void matrixMultiply(const Hybrid<T> &a, const size_t row_a, const size_t col_a, const Hybrid<T> &b,
                    const size_t row_b, const size_t col_b, Hybrid<T> *c, const T scale_a,
                    const T scale_b, const T scale_c, const TransposeState x_a,
                    const TransposeState x_b) {
  if (a.size() != row_a * col_a) {
    rtErr("The size of the input matrix A (" + std::to_string(a.size()) + ", Hybrid " +
          std::string(a.getLabel().name) + ") must comport with its stated row and column counts "
          "(" + std::to_string(row_a) + ", " + std::to_string(col_a) + ").", "matrixMultiply");
  }
  if (b.size() != row_b * col_b) {
    rtErr("The size of the input matrix B (" + std::to_string(b.size()) + ", Hybrid " +
          std::string(b.getLabel().name) + ") must comport with its stated row and column counts "
          "(" + std::to_string(row_b) + ", " + std::to_string(col_b) + ").", "matrixMultiply");
  }
  const size_t relevant_row_c = (x_a == TransposeState::AS_IS) ? row_a : col_a;
  const size_t relevant_col_c = (x_b == TransposeState::AS_IS) ? col_b : row_b;
  if (c->size() != relevant_row_c * relevant_col_c) {
    rtErr("The size of the output matrix C (" + std::to_string(c->size()) + ", Hybrid " +
          std::string(c->getLabel().name) + ") must comport with its stated row and column counts "
          "(" + std::to_string(relevant_row_c) + ", " + std::to_string(relevant_col_c) + ").",
          "matrixMultiply");
  }
  matrixMultiply(a.data(), row_a, col_a, b.data(), row_b, col_b, c->data(), scale_a, scale_b,
                 scale_c, x_a, x_b);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void matrixMultiply(const HpcMatrix<T> &a, const HpcMatrix<T> &b, HpcMatrix<T> *c, const T scale_a,
                    const T scale_b, const T scale_c, const TransposeState x_a,
                    const TransposeState x_b) {
  switch (x_a) {
  case TransposeState::AS_IS:
    if (a.n_rows != c->n_rows) {
      rtErr("The number of rows in the output matrix must equal the number of rows in the "
            "left-hand input matrix.  Currently " + std::to_string(a.n_rows) + " (input) vs " +
            std::to_string(c->n_rows) + " (output).", "matrixMultiply");
    }
    break;
  case TransposeState::TRANSPOSE:
    if (a.n_cols != c->n_rows) {
      rtErr("The number of rows in the output matrix must equal the number of columns in the "
            "left-hand input matrix (which is to be transposed).  Currently " +
            std::to_string(a.n_cols) + " (input) vs " + std::to_string(c->n_rows) + " (output).",
            "matrixMultiply");
    }
    break;
  }
  switch (x_b) {
  case TransposeState::AS_IS:
    if (b.n_cols != c->n_cols) {
      rtErr("The number of columns in the output matrix must equal the number of columns in the "
            "right-hand input matrix.  Currently " + std::to_string(b.n_cols) + " (input) vs " +
            std::to_string(c->n_cols) + " (output).", "matrixMultiply");
    }
    break;
  case TransposeState::TRANSPOSE:
    if (b.n_rows != c->n_cols) {
      rtErr("The number of columns in the output matrix must equal the number of rows in the "
            "right-hand input matrix (which is to be transposed).  Currently " +
            std::to_string(b.n_rows) + " (input) vs " + std::to_string(c->n_cols) + " (output).",
            "matrixMultiply");
    }
    break;
  }
  matrixMultiply(a.memptr(), a.n_rows, a.n_cols, b.memptr(), b.n_rows, b.n_cols, c->memptr(),
                 scale_a, scale_b, scale_c, x_a, x_b);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void matrixVectorMultiply(const T* a, const T* x, T* b, const size_t row_a, const size_t col_a,
                          const T scale_a, const T scale_x, const T scale_b,
                          const TransposeState x_a) {

  // Scale the vector b as needed
  for (size_t i = 0; i < row_a; i++) {
    b[i] *= scale_b;
  }
  const T scale_ax = scale_a * scale_x;

  // Handle transpose cases
  switch (x_a) {
  case TransposeState::AS_IS:

    // Step along the matrix's rows to have each one act upon the vector x.
    for (size_t i = 0; i < col_a; i++) {
      const size_t i_rowa = i * row_a;
      T dsum = 0.0;
      const T xval = scale_ax * x[i];
      for (size_t j = 0; j < row_a; j++) {
        b[j] += xval * a[i_rowa + j];
      }
    }
    break;
  case TransposeState::TRANSPOSE:

    // Step along the matrix's columns, treating them as rows to act upon the vector x.
    for (size_t i = 0; i < col_a; i++) {
      const size_t i_rowa = i * row_a;
      T dsum = 0.0;
      for (size_t j = 0; j < row_a; j++) {
        dsum += scale_ax * a[i_rowa + j] * x[j];
      }
      b[i] = dsum;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void matrixVectorMultiply(const std::vector<T> &a, const std::vector<T> &x, std::vector<T> *b,
                          const size_t row_a, const size_t col_a, const T scale_a, const T scale_x,
                          const T scale_b, const TransposeState x_a) {
  if (row_a * col_a != a.size()) {
    rtErr("The number of elements in the input matrix (" + std::to_string(a.size()) +
          ") must equal the product of its stated row and column counts (" +
          std::to_string(row_a) + " and " + std::to_string(col_a) + ", respectively).",
          "matrixVectorMultiply");
  }
  switch (x_a) {
  case TransposeState::AS_IS:
    if (row_a != b->size()) {
      rtErr("The number of elements in the output vector must equal the number of rows in the "
            "left-hand input matrix.  Currently " + std::to_string(row_a) + " (input) vs " +
            std::to_string(b->size()) + " (output).", "matrixVectorMultiply");
    }
    if (col_a != x.size()) {
      rtErr("The number of elements in the right-hand input vector must equal the number of "
            "columns in the left-hand input matrix.  Currently " + std::to_string(col_a) +
            " (matrix) vs " + std::to_string(x.size()) + " (vector).", "matrixVectorMultiply");
    }
    break;
  case TransposeState::TRANSPOSE:
    if (col_a != b->size()) {
      rtErr("The number of elements in the output vector must equal the number of columns in the "
            "left-hand input matrix (which is to be transposed).  Currently " +
            std::to_string(col_a) + " (input) vs " + std::to_string(b->size()) + " (output).",
            "matrixVectorMultiply");
    }
    break;
  }
  matrixVectorMultiply(a.data(), x.data(), b->data(), row_a, col_a, scale_a, scale_x, scale_b,
                       x_a);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void matrixVectorMultiply(const Hybrid<T> &a, const Hybrid<T> &x, Hybrid<T> *b,
                          const size_t row_a, const size_t col_a, const T scale_a, const T scale_x,
                          const T scale_b, const TransposeState x_a) {
  if (row_a * col_a != a.size()) {
    rtErr("The number of elements in the input matrix (" + std::to_string(a.size()) +
          ", Hybrid " + std::string(a.getLabel().name) + ") must equal the product of its stated "
          "row and column counts (" + std::to_string(row_a) + " and " + std::to_string(col_a) +
          ", respectively).", "matrixVectorMultiply");
  }
  switch (x_a) {
  case TransposeState::AS_IS:
    if (row_a != b->size()) {
      rtErr("The number of elements in the output vector must equal the number of rows in the "
            "left-hand input matrix.  Currently " + std::to_string(row_a) + " (input matrix, "
            "Hybrid " + std::string(a.getLabel().name) + ") vs " + std::to_string(b->size()) +
            " (output vector, Hybrid " + std::string(b->getLabel().name) + ").",
            "matrixVectorMultiply");
    }
    if (col_a != x.size()) {
      rtErr("The number of elements in the right-hand input vector must equal the number of "
            "columns in the left-hand input matrix.  Currently " + std::to_string(col_a) +
            " (matrix, Hybrid " + std::string(a.getLabel().name) + ") vs " +
            std::to_string(x.size()) + " (vector, Hybrid " + std::string(x.getLabel().name) + ").",
            "matrixVectorMultiply");
    }
    break;
  case TransposeState::TRANSPOSE:
    if (col_a != b->size()) {
      rtErr("The number of elements in the output vector must equal the number of columns in the "
            "left-hand input matrix (which is to be transposed).  Currently " +
            std::to_string(col_a) + " (input) vs " + std::to_string(b->size()) + " (output).",
            "matrixVectorMultiply");
    }
    break;
  }
  matrixVectorMultiply(a.data(), x.data(), b->data(), row_a, col_a, scale_a, scale_x, scale_b,
                       x_a);  
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void matrixVectorMultiply(const HpcMatrix<T> &a, const Hybrid<T> &x, Hybrid<T> *b, const T scale_a,
                          const T scale_x, const T scale_b, const TransposeState x_a) {
  switch (x_a) {
  case TransposeState::AS_IS:
    if (a.n_rows != b->size()) {
      rtErr("The number of elements in the output vector must equal the number of rows in the "
            "left-hand input matrix.  Currently " + std::to_string(a.n_rows) + " (input) vs " +
            std::to_string(b->size()) + " (output).", "matrixVectorMultiply");
    }
    if (a.n_cols != x.size()) {
      rtErr("The number of elements in the right-hand input vector must equal the number of "
            "columns in the left-hand input matrix.  Currently " + std::to_string(a.n_cols) +
            " (matrix) vs " + std::to_string(x.size()) + " (vector).", "matrixVectorMultiply");
    }
    break;
  case TransposeState::TRANSPOSE:
    if (a.n_cols != b->size()) {
      rtErr("The number of elements in the output vector must equal the number of columns in the "
            "left-hand input matrix (which is to be transposed).  Currently " +
            std::to_string(a.n_cols) + " (input) vs " + std::to_string(b->size()) + " (output).",
            "matrixVectorMultiply");
    }
    break;
  }
  matrixVectorMultiply(a.memptr(), x.data(), b->data(), a.n_rows, a.n_cols, scale_a, scale_x,
                       scale_b, x_a);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void transpose(T* a, size_t rows, size_t columns) {
  const size_t space = (rows * columns) - 1LLU;
  std::vector<bool> visited(rows * columns, false);
  visited[0] = true;
  visited[space] = true;
  size_t i = 1;
  while (i < space) {

    // Proceed to swap one cycle
    size_t focus_index = i;
    size_t cycle_begin = focus_index;
    T t = a[focus_index];
    do {
      size_t next = (focus_index * columns) % space;
      std::swap(a[next], t);
      visited[focus_index] = true;
      focus_index = next;
    } while (focus_index != cycle_begin);

    // Advance to the next cycle based on elements that have not been moved
    while (i < space && visited[i]) {
      i++;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void transpose(const T* a, T* a_transpose, size_t rows, size_t columns) {
  for (size_t i = 0; i < rows; i++) {
    for (size_t j = 0; j < columns; j++) {
      a_transpose[(i * columns) + j] = a[(j * rows) + i];
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void transpose(std::vector<T> *a, size_t rows, size_t columns) {
  if (a->size() != rows * columns) {
    rtErr("The stated dimensions of the matrix (" + std::to_string(rows) + ", " +
          std::to_string(columns) + ") do not match the array size (" +
          std::to_string(a->size()) + ").", "transpose");
  }
  transpose(a->data(), rows, columns);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void transpose(const std::vector<T> &a, std::vector<T> *a_transpose, size_t rows, size_t columns) {
  if (a.size() != rows * columns) {
    rtErr("The stated dimensions of the matrix (" + std::to_string(rows) + ", " +
          std::to_string(columns) + ") do not match the array size (" +
          std::to_string(a.size()) + ").", "transpose");
  }
  if (a.size() != a_transpose->size()) {
    rtErr("The transposed matrix array (" + std::to_string(a_transpose->size()) + " elements) is "
          "not of the same size as the original (" + std::to_string(a.size()) + " elements).",
          "transpose");
  }
  transpose(a.data(), a_transpose->data(), rows, columns);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void transpose(Hybrid<T> *a, size_t rows, size_t columns) {
  if (a->size() != rows * columns) {
    rtErr("The stated dimensions of the matrix (" + std::to_string(rows) + ", " +
          std::to_string(columns) + ") do not match the array size (" +
          std::to_string(a->size()) + ").", "transpose");
  }
  transpose(a->data(), rows, columns);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void transpose(const Hybrid<T> &a, Hybrid<T> *a_transpose, size_t rows, size_t columns) {
  if (a.size() != rows * columns) {
    rtErr("The stated dimensions of the matrix (" + std::to_string(rows) + ", " +
          std::to_string(columns) + ") do not match the array size (" +
          std::to_string(a.size()) + ").", "transpose");
  }
  if (a.size() != a_transpose->size()) {
    rtErr("The transposed matrix array (" + std::to_string(a_transpose->size()) + " elements) is "
          "not of the same size as the original (" + std::to_string(a.size()) + " elements).",
          "transpose");
  }
  transpose(a.data(), a_transpose->data(), rows, columns);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void invertSquareMatrix(T* matrix, T* inverse, const size_t rank) {

  // Return immediately if the rank is zero
  if (rank == 0LLU) {
    return;
  }

  // Assign a set of pointers to the input and output matrices.  Form the identity matrix in the
  // output matrix (inverse).
  std::vector<TwinPointer<T>> s;
  s.resize(rank);
  for (size_t i = 0; i < rank; i++) {
    s[i].ptr_a = &matrix[i * rank];
    s[i].ptr_b = &inverse[i * rank];
    for (size_t j = 0; j < rank; j++) {
      s[i].ptr_b[j] = 0.0;
    }
    s[i].ptr_b[i] = 1.0;
  }

  // Proceed with partial-pivoting Gauss-Jordan elimination, stepping along the ith row to obtain
  // a lower-triangular matrix
  for (size_t i = 0; i < rank; i++) {
 
    // Find the largest value in the ith row
    T pivot_val = fabs(s[i].ptr_a[i]);
    size_t pivot_col = i;
    for (size_t j = i; j < rank; j++) {
      const T pivot_test = fabs(s[j].ptr_a[i]);
      if (pivot_test > pivot_val) {
        pivot_val = pivot_test;
        pivot_col = j;
      }
    }

    // Swap the pivot row into the prime position
    std::swap(s[i], s[pivot_col]);

    // Destroy everything in the ith row with the pivot column, modifying other columns as needed
    const T point_val = s[i].ptr_a[i];
    for (size_t j = i + 1; j < rank; j++) {
      if (fabs(s[j].ptr_a[i]) < constants::verytiny) {
        continue;
      }
      const T factor = s[j].ptr_a[i] / point_val;
      for (size_t k = 0; k < rank; k++) {
        s[j].ptr_a[k] -= factor * s[i].ptr_a[k];
        s[j].ptr_b[k] -= factor * s[i].ptr_b[k];
      }
    }
  }

  // Starting with the final column, crush backwards to obtain the identity matrix in place of the
  // original matrix.  The output matrix will then contain the inverse of the input matrix.
  for (size_t i = rank; i > 0; i--) {
    const size_t ieff = i - 1;
    const T point_val = 1.0 / s[ieff].ptr_a[ieff]; 
    for (size_t j = 0; j < rank; j++) {
      s[ieff].ptr_b[j] *= point_val;
    }
    s[ieff].ptr_a[ieff] = 1.0;
    for (size_t j = 0; j < ieff; j++) {
      const T factor = s[j].ptr_a[ieff];
      s[j].ptr_a[ieff] = 0.0;
      for (size_t k = 0; k < rank; k++) {
        s[j].ptr_b[k] -= s[ieff].ptr_b[k] * factor;
      }
    }
  }

  // While the result is now set, the arrays of pointers used to perform row swaps in the
  // Gauss-Jordan elimination make the matrix incomprehensible.  The original matrix has been
  // mangled anyway-use it as scratch space to re-arrange the columns of the inverse matrix for
  // output.
  for (size_t i = 0; i < rank; i++) {
    for (size_t j = 0; j < rank; j++) {
      matrix[(rank * i) + j] = s[i].ptr_b[j];
    }
  }
  const size_t rank2 = rank * rank;
  for (size_t i = 0; i < rank2; i++) {
    inverse[i] = matrix[i];
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> void invertSquareMatrix(const T* matrix, T* inverse, const size_t rank) {

  // Create a copy of the matrix, thus accomplishing out-of-place matrix inversion
  const size_t rank2 = rank * rank;
  std::vector<T> matrix_copy(rank2, 0.0); 
  for (size_t i = 0; i < rank2; i++) {
    matrix_copy[i] = matrix[i];    
  }
  invertSquareMatrix(matrix_copy.data(), inverse, rank);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void invertSquareMatrix(const std::vector<T> &matrix, std::vector<T> *inverse) {
  std::vector<T> mat_copy = matrix;
  invertSquareMatrix(&mat_copy, inverse);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void invertSquareMatrix(std::vector<T> *matrix, std::vector<T> *inverse) {
  const size_t rank = round(sqrt(static_cast<double>(matrix->size())));
  if (matrix->size() != rank * rank) {
    rtErr("An array of " + std::to_string(matrix->size()) + " cannot represent a square matrix.",
          "invertSquareMatrix");
  }
  if (inverse->size() != matrix->size()) {
    rtErr("The inverse matrix must be of the same rank as the original.", "invertSquareMatrix");
  }
  invertSquareMatrix(matrix->data(), inverse->data(), rank);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void invertSquareMatrix(const Hybrid<T> &matrix, Hybrid<T> *inverse) {
  Hybrid<T> mat_copy = matrix;
  invertSquareMatrix(&mat_copy, inverse);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void invertSquareMatrix(Hybrid<T> *matrix, Hybrid<T> *inverse) {
  const size_t rank = round(sqrt(static_cast<double>(matrix->size())));
  if (matrix->size() != rank * rank) {
    rtErr("An array of " + std::to_string(matrix->size()) + " cannot represent a square matrix.",
          "invertSquareMatrix");
  }
  if (inverse->size() != matrix->size()) {
    rtErr("The inverse matrix must be of the same rank as the original.", "invertSquareMatrix");
  }
  invertSquareMatrix(matrix->data(), inverse->data(), rank);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void jacobiEigensolver(T* a, T* v, T* d, const size_t rank, const ExceptionResponse policy) {
  double tresh, theta, tau, t, sm, s, h, g, c;

  // Allocate local arrays for workspace
  std::vector<T> b(rank, 0.0);
  std::vector<T> z(rank, 0.0);

  // Ensure that the eigenvector matrix is set to the identity matrix
  for (size_t ip = 0LLU; ip < rank * rank; ip++) {
    v[ip] = 0.0;
  }
  for (size_t ip = 0LLU; ip < rank; ip++) {
    v[(ip * rank) + ip] = 1.0;
    d[ip] = a[(ip * rank) + ip];
    b[ip] = d[ip];
    z[ip] = 0.0;
  }
  
  // Ensure that the matrix is symmetric
  for (size_t ip = 1LLU; ip < rank; ip++) {
    for (size_t iq = 0; iq < ip; iq++) {
      if (fabs(a[(ip * rank) + iq] - a[(iq * rank) + ip]) > constants::tiny) {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("Matrix of rank " + std::to_string(rank) + " is not symmetric: " +
                realToString(a[(ip * rank) + iq]) + " != " + realToString(a[(iq * rank) + ip]) +
                ".", "jacobiEigensolver");
        case ExceptionResponse::WARN:
          rtWarn("Matrix of rank " + std::to_string(rank) + " is not symmetric: " +
                 realToString(a[(ip * rank) + iq]) + " != " + realToString(a[(iq * rank) + ip]) +
                 ".", "jacobiEigensolver");
          return;
        case ExceptionResponse::SILENT:
          return;
        }
      }
    }
  }
  
  // Iteratively solve the eigenvalue problem
  for (int i = 1; i <= 50; i++) {
    sm = 0.0;
    for (size_t ip = 0; ip < rank - 1LLU; ip++) {
      for (size_t iq = ip + 1LLU; iq < rank; iq++) {
        sm += fabs(a[iq * rank + ip]);
      }
    }
    if (sm == 0.0) {

      // The eigenvalues were computed successfully
      return;
    }
    if (i < 4) {
      tresh = 0.2 * sm / (rank * rank);
    }
    else {
      tresh = 0.0;
    }
    for (size_t ip = 0; ip < rank - 1LLU; ip++) {
      for (size_t iq = ip + 1LLU; iq < rank; iq++) {
        g = 100.0 * fabs(a[(iq * rank) + ip]);
        if (i > 4 && fabs(d[ip]) + g == fabs(d[ip]) &&
            fabs(d[iq]) + g == fabs(d[iq])) {
          a[(iq * rank) + ip] = 0.0;
        }
        else {
          if (fabs(a[(iq * rank) + ip]) > tresh) {
            h = d[iq] - d[ip];
            if (fabs(h) + g == fabs(h)) {
              t = (a[(iq * rank) + ip]) / h;
            }
            else {
              theta = 0.5 * h / (a[(iq * rank) + ip]);
              t = 1.0 / (fabs(theta) + sqrt(1.0 + (theta * theta)));
              if (theta < 0.0) {
                t = -t;
              }
            }
            c = 1.0 / sqrt(1.0 + t * t);
            s = t * c;
            tau = s / (1.0 + c);
            h = t * a[(iq * rank) + ip];
            z[ip] -= h;
            z[iq] += h;
            d[ip] -= h;
            d[iq] += h;
            a[(iq * rank) + ip] = 0.0;
            if (ip > 0LLU) {
              for (size_t j = 0LLU; j <= ip - 1LLU; j++) {
                g = a[(ip * rank) + j];
                h = a[(iq * rank) + j];
                a[(ip * rank) + j] = g - s * (h + g * tau);
                a[(iq * rank) + j] = h + s * (g - h * tau);
              }
            }
            if (iq > 0LLU) {
              for (size_t j = ip + 1LLU; j <= iq - 1LLU; j++) {
                g = a[(j * rank) + ip];
                h = a[(iq * rank) + j];
                a[(j * rank) + ip] = g - s * (h + g * tau);
                a[(iq * rank) + j] = h + s * (g - h * tau);
              }
            }
            for (size_t j = iq + 1LLU; j < rank; j++) {
              g = a[(j * rank) + ip];
              h = a[(j * rank) + iq];
              a[(j * rank) + ip] = g - s * (h + g * tau);
              a[(j * rank) + iq] = h + s * (g - h * tau);
            }
            for (size_t j = 0LLU; j < rank; j++) {
              g = v[(ip * rank) + j];
              h = v[(iq * rank) + j];
              v[(ip * rank) + j] = g - s * (h + g * tau);
              v[(iq * rank) + j] = h + s * (g - h * tau);
            }
          }
        }
      }
    }
    for (size_t ip = 0; ip < rank; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }
  
  switch (policy) {
  case ExceptionResponse::DIE:
    rtErr("Jacobi iterations on a " + std::to_string(rank) + "-rank matrix did not converge.",
          "jacobiEigensolver");
  case ExceptionResponse::WARN:
    rtWarn("Jacobi iterations on a " + std::to_string(rank) + "-rank matrix did not converge.",
           "jacobiEigensolver");
    break;
  case ExceptionResponse::SILENT:
    break;
  }
}
  
//-------------------------------------------------------------------------------------------------
template <typename T>
void jacobiEigensolver(std::vector<T> *a, std::vector<T> *v, std::vector<T> *d, const size_t rank,
                       const ExceptionResponse policy) {
  const size_t actual_rank = (rank == 0LLU) ? d->size() : rank;
  checkArraySizeToMatrixRank(actual_rank, rank, a->size(), v->size(), policy);
  jacobiEigensolver(a->data(), v->data(), d->data(), actual_rank, policy);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void jacobiEigensolver(Hybrid<T> *a, Hybrid<T> *v, Hybrid<T> *d, const size_t rank,
                       const ExceptionResponse policy) {
  const size_t actual_rank = (rank == 0LLU) ? d->size() : rank;
  checkArraySizeToMatrixRank(actual_rank, rank, a->size(), v->size(), policy);
  jacobiEigensolver(a->data(), v->data(), d->data(), actual_rank, policy);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void jacobiEigensolver(HpcMatrix<T> *a, HpcMatrix<T> *v, Hybrid<T> *d,
                       const ExceptionResponse policy) {
  if (a->n_rows != a->n_cols) {
    rtErr("Matrix must be square (current dimensions " + std::to_string(a->n_rows) + " by " +
          std::to_string(a->n_cols) + ").", "jacobiEigensolver");
  }
  else if (v->n_rows != v->n_cols || v->n_rows != a->n_rows) {
    rtErr("Eigenvector matrix must be square and match the rank of the input matrix (current "
          "dimensions " + std::to_string(v->n_rows) + " by " +
          std::to_string(v->n_cols) + ") against a rank " + std::to_string(a->n_rows) +
          "input matrix.", "jacobiEigensolver");
  }
  else if (d->size() != a->n_rows) {
    rtErr("Eigenvalue results require a vector matching the rank of the input matrix (" +
          std::to_string(a->n_rows) + "), currently " + std::to_string(d->size()) + ".",
          "jacobiEigensolver");
  }
  jacobiEigensolver(a->memptr(), v->memptr(), d->data(), a->n_rows, policy);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void realSymmEigensolver(T* amat, const int rank, T* eigv, T* sdiag, const EigenWork process) {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const bool t_is_double = (ct == double_type_index);
  const T zero = 0.0;
  const T one = 1.0;
  const T two = 2.0;
  for (int i = rank - 1; i >= 1; i--) {
    int l = i - 1;
    T h = zero;
    T scale = zero;
    if (l > 0) {
      for (int k = 0; k <= l; k++) {
        scale += std::abs(amat[(k * rank) + i]);
      }
      if (scale == zero) {
        sdiag[i] = amat[(l * rank) + i];
      }
      else {
        for (int k = 0; k <= l; k++) {
          const T aik = amat[(k * rank) + i] / scale;
          h += aik * aik;
          amat[(k * rank) + i] = aik;
        }
        T f = amat[(l * rank) + i];
        T g;
        if (t_is_double) {
          g = (f >= zero) ? -sqrt(h) : sqrt(h);
        }
        else {
          g = (f >= zero) ? -sqrtf(h) : sqrtf(h);
        }
        sdiag[i] = scale * g;
	h -= f * g;
        amat[(l * rank) + i] = f - g;
        f = 0.0;
        for (int j = 0; j <= l; j++) {
          switch (process) {
          case EigenWork::EIGENVALUES:
            break;
          case EigenWork::EIGENVECTORS:
            amat[(i * rank) + j] = amat[(j * rank) + i] / h;
            break;
          }
          g = 0.0;
          for (int k = 0; k <= j; k++) {
            g += amat[(k * rank) + j] * amat[(k * rank) + i];
          }
          for (int k = j+1; k <=l; k++) {
            g += amat[(j * rank) + k] * amat[(k * rank) + i];
          }
          sdiag[j] = g / h;
          f += sdiag[j] * amat[(j * rank) + i];
        }
        const T hh = f / (h + h);
	for (int j = 0; j <= l; j++) {
          f = amat[(j * rank) + i];
          g = sdiag[j] - (hh * f);
          sdiag[j] = g;
          for (int k = 0; k <= j; k++) {
            amat[(k * rank) + j] -= (f * sdiag[k]) + (g * amat[(k * rank) + i]);
          }
        }
      }
    }
    else {
      sdiag[i] = amat[(l * rank) + i];
    }
    eigv[i] = h;
  }

  // Accumulate the eigenvalues.
  sdiag[0] = zero;
  switch (process) {
  case EigenWork::EIGENVALUES:
    for (int i = 0; i < rank; i++) {
      eigv[i] = amat[(i * rank) + i];
    }
    break;
  case EigenWork::EIGENVECTORS:
    eigv[0] = zero;
    for (int i = 0; i < rank; i++) {
      int l = i - 1;
      if (eigv[i]) {
        for (int j = 0; j <= l; j++) {
          T g = zero;
          for (int k = 0; k <= l; k++) {
            g += amat[(k * rank) + i] * amat[(j * rank) + k];
          }
          for (int k = 0; k <= l; k++) {
            amat[(j * rank) + k] -= g * amat[(i * rank) + k];
          }
        }
      }
      eigv[i] = amat[(i * rank) + i];
      amat[(i * rank) + i] = one;
      for (int j = 0; j <= l; j++) {
        amat[(j * rank) + i] = zero;
        amat[(i * rank) + j] = zero;
      }
    }
    break;
  }
  for (int i = 1; i < rank; i++) {
    sdiag[i - 1] = sdiag[i];
  }
  sdiag[rank - 1] = zero;
  for (int l = 0; l < rank; l++) {
    int iter = 0;
    int m = l - 1;
    while (m != l) {
      for (m = l; m < rank - 1; m++) {
        T dd = std::abs(eigv[m]) + std::abs(eigv[m + 1]);
        if (std::abs(sdiag[m] + dd) == dd) {
          break;
        }
      }
      if (m != l) {
        if (iter++ == maximum_ql_iterations) {
          rtErr("Too many iterations in the QL solver.", "realSymmEigensolver");
        }
        T g = (eigv[l + 1] - eigv[l]) / (two * sdiag[l]);
        T r = hypotenuse<T>(g, one);
        T sign_result;
        sign_result = (g >= zero) ? std::abs(r) : -std::abs(r);
        g = eigv[m] - eigv[l] + (sdiag[l] / (g + sign_result));
        T c = one;
        T s = one;
        T p = zero;
        bool early_finish = false;
        for (int i = m - 1; i >= l; i--) {
          T f = s * sdiag[i];
          T b = c * sdiag[i];
          r = hypotenuse<T>(f, g);
          sdiag[i + 1] = r;
          if (r == zero) {
            eigv[i + 1] -= p;
            sdiag[m] = zero;
            early_finish = true;
            break;
          }
          s = f / r;
          c = g / r;
          g = eigv[i + 1] - p;
          r = ((eigv[i] - g) * s) + (two * c * b);
          eigv[i + 1] = g + (p = s * r);
          g = (c * r) - b;
          switch (process) {
          case EigenWork::EIGENVALUES:
            break;
          case EigenWork::EIGENVECTORS:
            for (int k = 0; k < rank; k++) {
              f = amat[((i + 1) * rank) + k];
              amat[((i + 1) * rank) + k] = (s * amat[(i * rank) + k]) + (c * f);
              amat[(i * rank) + k] = (c * amat[(i * rank) + k]) - (s * f);
            }
            break;
          }
        }
        if (r == zero && early_finish) {
          continue;
        }
        eigv[l] -= p;
        sdiag[l] = g;
        sdiag[m] = zero;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void realSymmEigensolver(std::vector<T> *amat, std::vector<T> *eigv, std::vector<T> *sdiag,
                         const EigenWork process) {
  const size_t rank = eigv->size();
  if (amat->size() != rank * rank) {
    rtErr("Eigenvalue storage is ready for a matrix of rank " + std::to_string(rank) +
          ", but the matrix rank does not appear to match.", "realSymmEigensolver");
  }
  realSymmEigensolver(amat->data(), rank, eigv->data(), sdiag->data());
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void realSymmEigensolver(const T* amat, T* vmat, const int rank, T* eigv, T* sdiag,
                         const EigenWork process) {
  const size_t rank2 = rank * rank;
  for (size_t i = 0LLU; i < rank2; i++) {
    vmat[i] = amat[i];
  }
  realSymmEigensolver(vmat, rank, eigv, sdiag);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void realSymmEigensolver(const std::vector<T> &amat, std::vector<T> *vmat, std::vector<T> *eigv,
                         std::vector<T> *sdiag, const EigenWork process) {
  const size_t rank = eigv->size();
  const size_t rank2 = rank * rank;
  if (amat.size() != rank2) {
    rtErr("Eigenvalue storage is ready for a matrix of rank " + std::to_string(rank) +
          ", but the matrix rank does not appear to match.", "realSymmEigensolver");
  }
  T* vmat_ptr = vmat->data();
  for (size_t i = 0LLU; i < rank2; i++) {
    vmat_ptr[i] = amat[i];
  }
  realSymmEigensolver(vmat->data(), rank, eigv->data(), sdiag->data());
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void pseudoInverse(T* amat, T* workspace, T* inv_workspace, const size_t rows,
                   const size_t columns) {

  // If the matrix is already square, just invert it
  if (rows == columns) {
    invertSquareMatrix(amat, workspace, columns);
    for (size_t i = 0; i < columns * columns; i++) {
      amat[i] = workspace[i];
    }
    return;
  }
  if (rows > columns) {

    // Compute A' * A for a "left" pseudo-inverse if a "tall" matrix
    matrixMultiply(amat, rows, columns, amat, rows, columns, workspace, 1.0, 1.0, 0.0,
                   TransposeState::TRANSPOSE, TransposeState::AS_IS);
    invertSquareMatrix(workspace, inv_workspace, columns);
  }
  else {

    // Compute A * A' for a "right" pseudo-inverse of a "wide" matrix
    matrixMultiply(amat, rows, columns, amat, rows, columns, workspace, 1.0, 1.0, 0.0,
                   TransposeState::AS_IS, TransposeState::TRANSPOSE);
    invertSquareMatrix(workspace, inv_workspace, rows);
  }

  // Perform an in-place transpose of the matrix A, then use the workspace array to perform an
  // in-place replacement of the matrix with inv(A' * A) * A' or A' * inv(A * A'), depending on
  // whether the original matrix was tall or wide.
  if (rows > columns) {
    transpose(amat, rows, columns);
    for (size_t i = 0; i < rows; i += columns) {
      T* achunk = &amat[i * columns];
      const size_t chunk_columns = std::min(columns, rows - i);
      matrixMultiply(inv_workspace, columns, columns, achunk, columns, chunk_columns, workspace,
                     1.0, 1.0, 0.0);
      for (size_t j = 0; j < chunk_columns * columns; j++) {
        achunk[j] = workspace[j];
      }
    }
  }
  else {

    // To compute A' * inv(A * A'), use the identity (A * B)' = B' * A'.  Thus, A' * inv(A * A') =
    // (inv(A * A')' * A)'.  By first transposing the inverse matrix, the in-place multiplication
    // can be accomplished in a manner similar to that in the case of the tall matrix.  The result
    // can then be transposed into the "right" pseudo-inverse.
    transpose(inv_workspace, rows, rows);
    for (size_t i = 0; i < columns; i += rows) {
      T* achunk = &amat[i * rows];
      const size_t chunk_columns = std::min(rows, columns - i);
      matrixMultiply(inv_workspace, rows, rows, achunk, rows, chunk_columns, workspace,
                     1.0, 1.0, 0.0);
      for (size_t j = 0; j < chunk_columns * rows; j++) {
        achunk[j] = workspace[j];
      }
    }
    transpose(amat, rows, columns);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void pseudoInverse(T* amat, const size_t rows, const size_t columns) {
  const size_t lesser = std::min(rows, columns);
  std::vector<T> workspace(lesser * lesser, 0.0), inv_workspace(lesser * lesser, 0.0);
  pseudoInverse(amat, workspace.data(), inv_workspace.data(), rows, columns);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void pseudoInverse(const T* amat, T* result, T* workspace, T* inv_workspace, const size_t rows,
                   const size_t columns) { 
  const size_t rc = rows * columns;
  for (size_t i = 0; i < rc; i++) {
    result[i] = amat[i];
  }
  pseudoInverse(result, workspace, inv_workspace, rows, columns);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void pseudoInverse(std::vector<T> *amat, std::vector<T> *workspace, std::vector<T> *inv_workspace,
                   size_t rows, size_t columns) {
  if (amat->size() != rows * columns) {
    rtErr("The size of the input matrix (" + std::to_string(amat->size()) + ") must equal the "
          "product of the stated numbers of rows (" + std::to_string(rows) + ") and columns (" +
          std::to_string(columns) + ").", "pseudoInverse");
  }
  const size_t lesser = std::min(rows, columns);
  const std::string pinv_type = (rows > columns) ? "left" : "right";
  if (workspace->size() != lesser * lesser) {
    rtErr("The workspace (size " + std::to_string(workspace->size()) + ") must accommodate a "
          "square matrix of rank " + std::to_string(lesser) + " to compute the " + pinv_type +
          " pseudo-inverse of a " + std::to_string(rows) + " x " + std::to_string(columns) +
          " matrix.", "pseudoInverse");
  }
  if (inv_workspace->size() != lesser * lesser) {
    rtErr("The inverse workspace (size " + std::to_string(workspace->size()) + ") must "
          "accommodate a square matrix of rank " + std::to_string(lesser) + " to compute the " +
          pinv_type + " pseudo-inverse of a " + std::to_string(rows) + " x " +
          std::to_string(columns) + " matrix.", "pseudoInverse");
  }
  pseudoInverse(amat->data(), workspace->data(), inv_workspace->data(), rows, columns);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void pseudoInverse(std::vector<T> *amat, size_t rows, size_t columns) {
  const size_t lesser = std::min(rows, columns);
  std::vector<T> workspace(lesser * lesser, 0.0), inv_workspace(lesser * lesser, 0.0);
  pseudoInverse(amat, &workspace, &inv_workspace, rows, columns);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void pseudoInverse(const std::vector<T> &amat, std::vector<T> *result, std::vector<T> *workspace,
                   std::vector<T> *inv_workspace, size_t rows, size_t columns) {
  if (amat->size() != rows * columns) {
    rtErr("The size of the input matrix (" + std::to_string(amat->size()) + ") must equal the "
          "product of the stated numbers of rows (" + std::to_string(rows) + ") and columns (" +
          std::to_string(columns) + ").", "pseudoInverse");
  }
  if (result->size() != amat->size()) {
    result->resize(amat->size());
  }
  const size_t rc = rows * columns;
  T* res_ptr = result->data();
  for (size_t i = 0; i < rc; i++) {
    res_ptr[i] = amat[i];
  }
  pseudoInverse(result, workspace, inv_workspace, rows, columns);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> pseudoInverse(const std::vector<T> &amat, size_t rows, size_t columns) {
  std::vector<T> result = amat;
  const size_t lesser = std::min(rows, columns);
  std::vector<T> workspace(lesser * lesser, 0.0), inv_workspace(lesser * lesser, 0.0);
  pseudoInverse(&result, &workspace, &inv_workspace, rows, columns);
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<T> pseudoInverse(const std::vector<T> &amat, std::vector<T> *workspace,
                             std::vector<T> *inv_workspace, size_t rows, size_t columns) {
  std::vector<T> result = amat;
  pseudoInverse(&result, workspace, inv_workspace, rows, columns);
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void pseudoInverse(Hybrid<T> *amat, Hybrid<T> *workspace, Hybrid<T> *inv_workspace, size_t rows,
                   size_t columns) {
  if (amat->size() != rows * columns) {
    rtErr("The size of the input matrix (" + std::to_string(amat->size()) + ", Hybrid " +
          std::string(amat->getLabel().name) + ") must equal the product of the stated numbers of "
          "rows (" + std::to_string(rows) + ") and columns (" + std::to_string(columns) + ").",
          "pseudoInverse");
  }
  const size_t lesser = std::min(rows, columns);
  const std::string pinv_type = (rows > columns) ? "left" : "right";
  if (workspace->size() != lesser * lesser) {
    rtErr("The workspace (size " + std::to_string(workspace->size()) + ", Hybrid " +
          std::string(inv_workspace->getLabel().name) + ") must accommodate a square matrix of "
          "rank " + std::to_string(lesser) + " to compute the " + pinv_type +
          " pseudo-inverse of a " + std::to_string(rows) + " x " + std::to_string(columns) +
          " matrix.", "pseudoInverse");
  }
  if (inv_workspace->size() != lesser * lesser) {
    rtErr("The inverse workspace (size " + std::to_string(workspace->size()) + ", Hybrid " +
          std::string(workspace->getLabel().name) + ") must accommodate a square matrix of rank " +
          std::to_string(lesser) + " to compute the " + pinv_type + " pseudo-inverse of a " +
          std::to_string(rows) + " x " + std::to_string(columns) + " matrix.", "pseudoInverse");
  }
  pseudoInverse(amat->data(), workspace->data(), inv_workspace->data(), rows, columns);  
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void pseudoInverse(Hybrid<T> *amat, size_t rows, size_t columns) {

  // Perform the size check on the amat object--the mixture of std::vector<T> and Hybrid<T> must
  // go directly to the baseline C-style array overload.
  if (amat->size() != rows * columns) {
    rtErr("The size of the input matrix (" + std::to_string(amat->size()) + ", Hybrid " +
          std::string(amat->getLabel().name) + ") must equal the product of the stated numbers of "
          "rows (" + std::to_string(rows) + ") and columns (" + std::to_string(columns) + ").",
          "pseudoInverse");
  }
  std::vector<T> workspace(columns * columns);
  std::vector<T> inv_workspace(columns * columns);
  pseudoInverse(amat->data(), workspace->data(), inv_workspace->data(), rows, columns);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void pseudoInverse(const Hybrid<T> &amat, Hybrid<T> *result, Hybrid<T> *workspace,
                   Hybrid<T> *inv_workspace, size_t rows, size_t columns) {
  if (amat->size() != rows * columns) {
    rtErr("The size of the input matrix (" + std::to_string(amat->size()) + ", Hybrid " +
          std::string(amat->getLabel().name) + ") must equal the product of the stated numbers of "
          "rows (" + std::to_string(rows) + ") and columns (" + std::to_string(columns) + ").",
          "pseudoInverse");
  }
  if (result->size() != amat->size()) {
    result->resize(amat->size());
  }
  const size_t rc = rows * columns;
  T* res_ptr = result->data();
  for (size_t i = 0; i < rc; i++) {
    res_ptr[i] = amat[i];
  }
  pseudoInverse(result, workspace, inv_workspace, rows, columns);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void qrSolver(T* amat, T* xvec, T* bvec, const size_t a_rows, const size_t a_cols) {

  // QR decomposition of rank-deficient matrices will not be attempted
  if (a_rows < a_cols) {
    rtErr("The matrix is rank deficient (" + std::to_string(a_rows) + " rows, " +
          std::to_string(a_cols) + " columns).", "qrSolver");
  }

  // Constants in the relevant data type
  const T value_negative_one = -1.0;
  const T value_zero         =  0.0;
  const T value_one          =  1.0;
  const T value_two          =  2.0;
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const bool tcalc_is_double = (ct == double_type_index);
  
  // Allocate workspace
  std::vector<T> v(a_rows), v_prime(a_rows);
  for (size_t k = 0; k < a_cols; k++) {

    // Compute the kth column of Q*
    T tnm_v2 = 0.0;
    const size_t arowk = (k * a_rows) + k;
    for (size_t i = 0; i < a_rows - k; i++) {
      v[i] = amat[arowk + i];
      tnm_v2 += v[i] * v[i];
    }
    const T sign_v = value_negative_one + (value_two * (v[0] >= value_zero));
    T tnm_v = (tcalc_is_double) ? sqrt(tnm_v2) : sqrtf(tnm_v2);
    tnm_v2 -= v[0] * v[0];
    v[0] += sign_v * tnm_v;
    tnm_v = (tcalc_is_double) ? value_one / sqrt(tnm_v2 + (v[0] * v[0])) :
                                value_one / sqrtf(tnm_v2 + (v[0] * v[0]));
    for (size_t i = 0; i < a_rows - k; i++) {
      v[i] = v[i] * tnm_v;
    }

    // Update A as R evolves
    for (size_t i = 0; i < a_cols - k; i++) {
      v_prime[i] = value_zero;
    }

    const size_t acolsmk = a_cols - k;
    const size_t arowsmk = a_rows - k;
    for (size_t j = 0; j < acolsmk; j++) {
      T vp_acc = 0.0;
      const size_t arowkjk = ((k + j) * a_rows) + k;
      for (size_t i = 0; i < arowsmk; i++) {
        vp_acc += v[i] * amat[arowkjk + i];
      }
      v_prime[j] += vp_acc;
    }
    for (size_t j = 0; j < acolsmk; j++) {
      const T tempval = value_two * v_prime[j];
      const size_t arowkjk = ((k + j) * a_rows) + k;
      for (size_t i = 0; i < arowsmk; i++) {
        amat[arowkjk + i] -= tempval * v[i];        
      }
    }

    // Update b as Q* evolves
    const T tempval = value_two * dot(v.data(), &bvec[k], a_rows - k);
    for (size_t i = 0; i < a_rows - k; i++) {
      bvec[k + i] -= tempval * v[i];
    }
  }

  // Use back-substitution to solve the system of equations R * x = b.
  for (size_t i = a_cols - 1; i > 0; i--) {
    const T pivot = value_one / amat[(i * a_rows) + i];
    const T temp_b = bvec[i];

    // Use the criterion that j stays less than the number of rows to anticipate the fact that
    // the always positive size_t will wrap after passing below zero.
    const size_t arowi = i * a_rows;
    for (size_t j = i - 1; j < a_rows; j--) {
      const T multval = amat[arowi + j] * pivot;
      bvec[j] -= multval * temp_b;
    }
    bvec[i] *= pivot;
  }
  bvec[0] /= amat[0];

  // Read the entries of x out of b
  for (size_t i = 0; i < a_cols; i++) {
    xvec[i] = bvec[i];
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void qrSolver(std::vector<T> *amat, std::vector<T> *xvec, std::vector<T> *bvec) {
  if (amat->size() != xvec->size() * bvec->size()) {
    rtErr("The stiffness matrix must hold data commensurate with the sizes of the solution "
          "vector (" + std::to_string(bvec->size()) + ") and coefficients vector (" +
          std::to_string(xvec->size()) + ") counts.  The stiffness matrix holds " +
          std::to_string(amat->size()) + " elements.", "qrSolver");
  }
  qrSolver(amat->data(), xvec->data(), bvec->data(), bvec->size(), xvec->size());
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void qrSolver(const std::vector<T> &amat, std::vector<T> *xvec, const std::vector<T> &bvec) {
  std::vector<T> amat_copy = amat;
  std::vector<T> bvec_copy = bvec;
  qrSolver(amat_copy.data(), xvec->data(), bvec_copy.data(), bvec.size(), xvec->size());
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void qrSolver(Hybrid<T> *amat, Hybrid<T> *xvec, Hybrid<T> *bvec) {
  if (amat->size() != xvec->size() * bvec->size()) {
    rtErr("The stiffness matrix must hold data commensurate with the sizes of the solution "
          "vector (" + std::to_string(bvec->size()) + ") and coefficients vector (" +
          std::to_string(xvec->size()) + ") counts.  The stiffness matrix holds " +
          std::to_string(amat->size()) + " elements.", "qrSolver");
  }
  qrSolver(amat->data(), xvec->data(), bvec->data(), bvec->size(), xvec->size());
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void qrSolver(const Hybrid<T> &amat, Hybrid<T> *xvec, const Hybrid<T> &bvec) {
  std::vector<T> amat_copy = amat;
  std::vector<T> bvec_copy = bvec;
  qrSolver(amat_copy.data(), xvec->data(), bvec_copy.data(), bvec.size(), xvec->size());
}
  
//-------------------------------------------------------------------------------------------------
template <typename T>
double leibnizDeterminant(const T* amat, const size_t rank) {
  if (rank == 1) {
    return amat[0];
  }
  else if (rank == 2) {
    return (amat[0] * amat[3]) - (amat[1] * amat[2]);
  }
  else if (rank == 3) {
    return (amat[0] * ((amat[4] * amat[8]) - (amat[5] * amat[7]))) -
           (amat[3] * ((amat[1] * amat[8]) - (amat[2] * amat[7]))) +
           (amat[6] * ((amat[1] * amat[5]) - (amat[2] * amat[4])));
  }
  else if (rank <= 8) {
    double detsum = 0.0;
    double sgnfac = 1.0;
    const int irank = rank;
    const int subrank = rank - 1;
    std::vector<T> submat(subrank * subrank);
    for (int i = 0; i < irank; i++) {
      int jcon = 0;
      for (int j = 0; j < irank; j++) {
        if (j == i) {
          continue;
        }
        int kcon = 0;
        for (int k = 1; k < rank; k++) {
          submat[(jcon * subrank) + kcon] = amat[(j * rank) + k];
          kcon++;
        }
        jcon++;
      }
      detsum += sgnfac * amat[i * rank] * leibnizDeterminant(submat.data(), subrank);
      sgnfac = -sgnfac;
    }
    return detsum;
  }
  else {
    rtErr("The Leibniz formula is numerically ill-conditioned and requires O(N! * N) operations.  "
          "Choose a QR or LU decomposition to obtain the determinant of larger matrices.",
          "leibNizDeterminant");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
double leibnizDeterminant(const std::vector<T> &amat) {
  const size_t rank = sqrt(static_cast<double>(amat.size()));
  if (rank * rank != amat.size()) {
    rtErr("An array of length " + std::to_string(amat.size()) + " cannot represent a square "
          "matrix.", "leibnizDeterminant");
  }
  return leibnizDeterminant(amat.data(), rank);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
double leibnizDeterminant(const Hybrid<T> &amat) {
  const size_t rank = sqrt(static_cast<double>(amat.size()));
  if (rank * rank != amat.size()) {
    rtErr("An array of length " + std::to_string(amat.size()) + " cannot represent a square "
          "matrix.", "leibnizDeterminant");
  }
  return leibnizDeterminant(amat.data(), rank);
}
  
//-------------------------------------------------------------------------------------------------
template <typename T>
void computeBoxTransform(const T lx, const T ly, const T lz, const T alpha, const T beta,
                         const T gamma, T* umat, T* invu) {
  const T dx = ((cos(beta) * cos(gamma)) - cos(alpha)) / (sin(beta) * sin(gamma));
  const T dy = sqrt(1.0 - (dx * dx));
  invu[0] =  lx;
  invu[1] = 0.0;
  invu[2] = 0.0;
  invu[3] =  ly * cos(gamma);
  invu[5] = 0.0;
  invu[6] =  lz * cos(beta);
  invu[4] =  ly * sin(gamma);
  invu[7] = -lz * sin(beta) * dx;
  invu[8] =  lz * sin(beta) * dy;

  std::vector<T> invu_copy(9, 0.0);
  for (int i = 0; i < 9; i++) {
    invu_copy[i] = invu[i];
  }
  invertSquareMatrix(invu_copy.data(), umat, 3);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void computeBoxTransform(const std::vector<T> &dims, std::vector<T> *umat,
                                               std::vector<T> *invu) {
  computeBoxTransform(dims[0], dims[1], dims[2], dims[3], dims[4], dims[5], umat->data(),
                      invu->data());
}

//-------------------------------------------------------------------------------------------------
template <typename T> void computeBoxTransform(const T* dims, T* umat, T* invu) {
  computeBoxTransform(dims[0], dims[1], dims[2], dims[3], dims[4], dims[5], umat, invu);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void extractBoxDimensions(T *lx, T *ly, T *lz, T *alpha, T *beta, T *gamma, const T* invu) {

  // The first box length is still in the matrix
  *lx = invu[0];

  // Determine the second box length based on (sin^2)(x) + (cos^2)(x) = 1.
  const T ly_cos_gamma = invu[3];
  const T ly_sin_gamma = invu[4];
  *ly = sqrt((ly_cos_gamma * ly_cos_gamma) + (ly_sin_gamma * ly_sin_gamma));

  // Now determine the third box angle
  *gamma = acos(ly_cos_gamma / (*ly));

  // Go after beta and lz next.
  // (lz * sin(beta) * dy)^2 = lz^2 * (sin^2)(beta) * (1.0 - dx * dx)
  //                         = lz^2 * (sin^2)(beta) - lz^2 * (sin^2)(beta) * dx * dx
  // (lz + cos(beta))^2      = lz^2 * (cos^2)(beta)
  // invu[8] * invu[8] + invu[6] * invu[6] = lz^2 - lz^2 * (sin^2)(beta) * dx * dx
  // invu[8] * invu[8] + invu[6] * invu[6] = lz^2 - invu[7] * invu[7];
  *lz = sqrt((invu[6] * invu[6]) + (invu[7] * invu[7]) + (invu[8] * invu[8]));
  *beta = acos(invu[6] / (*lz));

  // Finally, go after alpha
  const T dx = -invu[7] / (sin(*beta) * (*lz));
  *alpha = acos((cos(*beta) * cos(*gamma)) - (dx * sin(*beta) * sin(*gamma)));
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> hessianNormalWidths(const T* invu) {
  T x, y, z;
  hessianNormalWidths(invu, &x, &y, &z);
  return { x, y, z };
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<T> hessianNormalWidths(const std::vector<T> &invu) {
  T x, y, z;
  hessianNormalWidths(invu.data(), &x, &y, &z);
  return { x, y, z };
}

//-------------------------------------------------------------------------------------------------
template <typename T> void hessianNormalWidths(const std::vector<T> &invu, T *x, T *y, T *z) {
  hessianNormalWidths(invu.data(), x, y, z);
}
  
//-------------------------------------------------------------------------------------------------
template <typename T> void hessianNormalWidths(const T* invu, T *x, T *y, T *z) {
  double xyz_n[3], thx[3], thy[3], thz[3], column_a[3], column_b[3], column_c[3];
  xyz_n[0] = invu[0] + invu[3] + invu[6];
  xyz_n[1] = invu[4] + invu[7];
  xyz_n[2] = invu[8];
  for (int i = 0; i < 3; i++) {
    column_a[i] = invu[i    ];
    column_b[i] = invu[i + 3];
    column_c[i] = invu[i + 6];
  }
  crossProduct(column_b, column_c, thx);
  crossProduct(column_a, column_c, thy);
  crossProduct(column_a, column_b, thz);
  normalize(thx, 3);
  normalize(thy, 3);
  normalize(thz, 3);
  *x = fabs(dot(thx, xyz_n, 3));
  *y = fabs(dot(thy, xyz_n, 3));
  *z = fabs(dot(thz, xyz_n, 3));
}
  
//-------------------------------------------------------------------------------------------------
template <typename T>
void printMatrix(const T* matrix, const int rows, const int cols, const std::string &varname,
                 const std::string &file_name, const PrintSituation expectation) {
  if (file_name.size() == 0LLU) {

    // Print to screen
    printf("%s = [\n", varname.c_str());
    for (int i = 0; i < rows; i++) {
      for (int j = 0; j < cols; j++) {
        printf("%16.8e", matrix[(j * rows) + i]);
      }
      printf("\n");
    }
    printf("];\n");
  }
  else {

    // Print to a file
    std::vector<PolyNumeric> smat;
    smat.resize(rows * cols);
    for (int i = 0; i < rows * cols; i++) {
      smat[i].d = matrix[i];
    }
    writeSnapshot(file_name, smat, varname, 1.0e-8, NumberFormat::SCIENTIFIC, expectation);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void printMatrix(const std::vector<T> &matrix, const int rows, const int cols,
                 const std::string &varname, const std::string &file_name,
                 const PrintSituation expectation) {
  printMatrix(matrix.data(), rows, cols, varname, file_name, expectation);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void printMatrix(const Hybrid<T> &matrix, const int rows, const int cols,
                 const std::string &varname, const std::string &file_name,
                 const PrintSituation expectation) {
  printMatrix(matrix.data(), rows, cols, varname, file_name, expectation);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void printMatrix(const HpcMatrix<T> &matrix, const std::string &varname,
                 const std::string &file_name, const PrintSituation expectation) {
  printMatrix(matrix.memptr(), matrix.n_rows, matrix.n_cols, varname, file_name, expectation);
}

} // namespace stormm
} // namespace stormm
