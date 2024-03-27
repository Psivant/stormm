// -*-c++-*-
#ifndef STORMM_MATRIX_OPS_H
#define STORMM_MATRIX_OPS_H

#include <algorithm>
#include <cmath>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/scaling.h"
#include "Accelerator/hybrid.h"
#include "FileManagement/file_enumerators.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "UnitTesting/file_snapshot.h"
#include "matrix.h"
#include "vector_ops.h"

namespace stormm {
namespace stmath {

using constants::ExceptionResponse;
using card::Hybrid;
using card::HybridTargetLevel;
using diskutil::PrintSituation;
using parse::NumberFormat;
using parse::PolyNumeric;
using parse::realToString;
using testing::writeSnapshot;

/// \brief The maximum number of QL iterations that will be allowed before the algorithm is
///        declared non-convergent.
constexpr int maximum_ql_iterations = 30;
  
/// \brief A twin pointer, useful for performing partial-pivoting Gaussian elimination on
///        row-major matrices.
template <typename T> struct TwinPointer {
  T* ptr_a;  ///< Pointer to the first stretch of data
  T* ptr_b;  ///< Pointer to the second stretch of data
};

/// \brief The enumerator makes the production of eigenvectors optional.
enum class EigenWork {
  EIGENVALUES,   ///< Compute eigenvalues only
  EIGENVECTORS,  ///< Compute both eigenvalues and eigenvectors
};
  
/// \brief Enumerate the transpose states of a matrix, for some basic matrix operations
enum class TransposeState {
  AS_IS,     ///< The matrix shall be taken as it is found
  TRANSPOSE  ///< The matrix shall be handled by first taking its transpose (the original matrix
             ///<   will not be disturbed)
};

/// \brief Multiply two matrices, with or without transposition of either, after multiplying each
///        by scalar prefactors and any pre-existing result matrix by its own prefactor.  This
///        routine is designed to emulate BLAS functionality, with a slightly different API.
///
/// Overloaded:
///   - Operate on raw pointers to C-style arrays, with trusted row and column bounds
///   - Accept Standard Template Library objects
///   - Accept Hybrid objects for the arrays
///   - Take HpcMatrices as inputs, with desired transposes
///
/// \param a        Hereafter "matrix A", the first of the operand matrices
/// \param row_a    Number of rows in matrix A
/// \param col_a    Number of columns in matrix A
/// \param b        Hereafter "matrix B", the second of the operand matrices
/// \param row_b    Number of rows in matrix B
/// \param col_b    Number of columns in matrix B
/// \param c        Hereafter "matrix C", the third of the operand matrices
/// \param scale_a  Scalar prefactor for matrix A
/// \param scale_b  Scalar prefactor for matrix B
/// \param scale_c  Scalar prefactor for matrix C
/// \param x_a      Transpose instruction for matrix A
/// \param x_b      Transpose instruction for matrix B
/// \{
template <typename T>
void matrixMultiply(const T* a, size_t row_a, size_t col_a, const T* b, size_t row_b, size_t col_b,
                    T* c, T scale_a = 1.0, T scale_b = 1.0, T scale_c = 1.0,
                    TransposeState x_a = TransposeState::AS_IS,
                    TransposeState x_b = TransposeState::AS_IS);

template <typename T>
void matrixMultiply(const std::vector<T> &a, size_t row_a, size_t col_a, const std::vector<T> &b,
                    size_t row_b, size_t col_b, std::vector<T> *c, T scale_a = 1.0,
                    T scale_b = 1.0, T scale_c = 1.0, TransposeState x_a = TransposeState::AS_IS,
                    TransposeState x_b = TransposeState::AS_IS);

template <typename T>
void matrixMultiply(const Hybrid<T> &a, size_t row_a, size_t col_a, const Hybrid<T> &b,
                    size_t row_b, size_t col_b, Hybrid<T> *c, T scale_a = 1.0, T scale_b = 1.0,
                    T scale_c = 1.0, TransposeState x_a = TransposeState::AS_IS,
                    TransposeState x_b = TransposeState::AS_IS);

template <typename T>
void matrixMultiply(const HpcMatrix<T> &a, const HpcMatrix<T> &b, HpcMatrix<T> *c, T scale_a = 1.0,
                    T scale_b = 1.0, T scale_c = 1.0, TransposeState x_a = TransposeState::AS_IS,
                    TransposeState x_b = TransposeState::AS_IS);
/// \}

/// \brief Multiply a vector by a matrix, to arrive at another vector.  The constants for scaling
///        the matrix, input vector, or output vector (prior to adding them matrix-vector product)
///        are carried over from the more general matrix-matrix multiply routine.  Transposing the
///        input matrix prior to the matrix-vector multiplication is also feasible.
///
/// Overloaded:
///   - Operate on raw pointers to C-style arrays, with trusted row and column bounds
///   - Accept Standard Template Library objects
///   - Accept Hybrid objects for the arrays
///   - Take an HpcMatrix, with either desired transpose, plus Hybrid object "vectors" as inputs 
///
/// \param a        Hereafter "matrix A", the matrix in b = scale_b * b + scale_a * scale_x * (A x)
/// \param x        The left-hand side vector to use in multiplication, hereafter "vector X"
/// \param b        The right-hand side vector, hereafter, "vector B"
/// \param row_a    The number of rows in matrix A (taken to be the number of elements in vector B,
///                 unless A is transposed in the calculation, in which case it is taken as the
///                 number of elements in vector X))
/// \param col_a    The number of columns in matrix A (taken to be the number of elements in
///                 vector X unless A is transposed, when it would be the number of elements in
///                 vector B)
/// \param scale_a  Scaling factor for matrix A
/// \param scale_x  Scaling factor for vector X
/// \param scale_b  Scaling factor for vector B
/// \param x_a      Indicator of whether to do a virtual transpose of A prior to the computation
/// \{
template <typename T>
void matrixVectorMultiply(const T* a, const T* x, T* b, size_t row_a, size_t col_a,
                          T scale_a = 1.0, T scale_x = 1.0, T scale_b = 1.0,
                          TransposeState x_a = TransposeState::AS_IS);

template <typename T>
void matrixVectorMultiply(const std::vector<T> &a, const std::vector<T> &x, std::vector<T> *b,
                          size_t row_a, size_t col_a, T scale_a = 1.0, T scale_x = 1.0,
                          T scale_b = 1.0, TransposeState x_a = TransposeState::AS_IS);

template <typename T>
void matrixVectorMultiply(const Hybrid<T> &a, const Hybrid<T> &x, Hybrid<T> *b, size_t row_a,
                          size_t col_a, T scale_a = 1.0, T scale_x = 1.0, T scale_b = 1.0,
                          TransposeState x_a = TransposeState::AS_IS);

template <typename T>
void matrixVectorMultiply(const HpcMatrix<T> &a, const Hybrid<T> &x, Hybrid<T> *b, T scale_a = 1.0,
                          T scale_x = 1.0, T scale_b = 1.0,
                          TransposeState x_a = TransposeState::AS_IS);
/// \}

/// \brief Compute the transpose of a matrix.
///
/// Overloaded:
///   - Compute the in-place transpose using a pre-allocated workspace of known length (this will
///     be checked for minimum space requirements)
///   - Compute the out-of-place transpose
///   - Operate on C-style arrays of trusted lengths
///   - Operate on Standard Template Library Vectors
///   - Operate on Hybrid objects
/// \{
template <typename T>
void transpose(T* a, size_t rows, size_t columns, size_t n_work);

template <typename T>
void transpose(const T* a, T* a_transpose, size_t rows, size_t columns);

template <typename T>
void transpose(std::vector<T> *a, size_t rows, size_t columns, size_t n_work);

template <typename T>
void transpose(const std::vector<T> &a, std::vector<T> *a_transpose, size_t rows, size_t columns);

template <typename T>
void transpose(Hybrid<T> *a, size_t rows, size_t columns, size_t n_work);

template <typename T>
void transpose(const Hybrid<T> &a, Hybrid<T> *a_transpose, size_t rows, size_t columns);
/// \}
  
/// \brief Invert a square matrix.  Simple function to avoid the heavy lift and compilation of
///        BLAS calls.  Templated for single- and double-precision functionality.
///
/// Overloaded:
///   - Provide C-style arrays of trusted length, Standard Template Library vectors, or Hybrid
///     objects
///   - Preserve the original matrix
///   - Destroy the original matrix 
///
/// \param matrix   The matrix to invert, in column-mjaor (contiguous columns) format
/// \param inverse  The inverse matrix result, pre-allocated and returned in column-major format
/// \param rank     Matrix rank, implying a trusted length if C-style arrays are supplied
/// \{
template <typename T> void invertSquareMatrix(const T* matrix, T* inverse, size_t rank);

template <typename T> void invertSquareMatrix(T* matrix, T* inverse, size_t rank);

template <typename T> void invertSquareMatrix(const std::vector<T> &matrix,
                                              std::vector<T> *inverse);

template <typename T> void invertSquareMatrix(std::vector<T> *matrix, std::vector<T> *inverse);

template <typename T> void invertSquareMatrix(const Hybrid<T> &matrix, Hybrid<T> *inverse);

template <typename T> void invertSquareMatrix(Hybrid<T> *matrix, Hybrid<T> *inverse);
/// \}

/// \brief Check the dimensions of matrices against the determined rank.  This function serves the
///        jacobiEigensolver routine.
///
/// \param actual_rank  The presumed rank of the matrix
/// \param s_a          Size of the array allocated for matrix A
/// \param s_v          Size of the array allocated for eigenvector matrix V (or some other
///                     matrix)
/// \param policy       Procedure in the event that the matrices are not of the right size
void checkArraySizeToMatrixRank(size_t actual_rank, size_t rank, size_t s_a, size_t s_v,
                                ExceptionResponse policy);
  
/// \brief Compute the eigenvalues and eigenvectors of a real, symmetric matrix.
///
///  Overloaded:
///    - Take templated type pointers for the matrices and vectors, with an explicitly stated rank
///    - Take std::vectors for the matrices and vector, with the rank to be inferred if it is not
///      given explicitly
///    - Take Hybrid objects for the matrices and vector, with the rank to be inferred if it is not
///      explicitly stated
///
/// \param a       The dense matrix of interest
/// \param v       Pre-allocated matrix which will eventually hold the eiegnvectors
/// \param d       Pre-allocated vector which will eventually hold the eigenvalues
/// \param rank    Rank of the matrix
/// \param policy  Procedure in the event that the iterative solver does not converge
/// \{
template <typename T>
void jacobiEigensolver(T* a, T* v, T* d, size_t rank,
                       ExceptionResponse policy = ExceptionResponse::WARN);

template <typename T>
void jacobiEigensolver(std::vector<T> *a, std::vector<T> *v, std::vector<T> *d, size_t rank = 0,
                       ExceptionResponse policy = ExceptionResponse::WARN);

template <typename T>
void jacobiEigensolver(Hybrid<T> *a, Hybrid<T> *v, Hybrid<T> *d, size_t rank = 0,
                       ExceptionResponse policy = ExceptionResponse::WARN);

template <typename T>
void jacobiEigensolver(HpcMatrix<T> *a, HpcMatrix<T> *v, Hybrid<T> *d,
                       ExceptionResponse policy = ExceptionResponse::WARN);
/// \}

/// \brief Solve the eigenvalues (and, optionally) eigenvectors of a symmetric, real matrix.
///
/// Overloaded:
///   - Provide C-style arrays for the matrix, diagonal, and proto-eigenvalue storage
///   - Provide Standard Template Library vectors for all of the above
///   - Provide Hybrid objects for all of the above
///
/// \param amat     The matrix from which to obtain eigenvalues and eigenvectors.  If eigenvectors
///                 are requested and amat is all that is supplied, they will replace the contents
///                 of amat on output.
/// \param vmat     Optional matrix for out-of-place eigenvector computation.  The original matrix
///                 amat will be preserved if this is provided.
/// \param rank     Rank of the matrix
/// \param diag     Storage space for the matrix diagonal
/// \param eigv     Storage space for computed eigenvalues
/// \param process  Toggle the computation of eigenvectors (active by default)
/// \{
template <typename T>
void realSymmEigensolver(T* amat, const int rank, T* diag, T* eigv,
                         EigenWork process = EigenWork::EIGENVECTORS);

template <typename T>
void realSymmEigensolver(std::vector<T> *amat, std::vector<T> *diag, std::vector<T> *eigv,
                         EigenWork process = EigenWork::EIGENVECTORS);

template <typename T>
void realSymmEigensolver(const T* amat, T* vmat, const int rank, T* diag, T* eigv,
                         EigenWork process = EigenWork::EIGENVECTORS);

template <typename T>
void realSymmEigensolver(const std::vector<T> &amat, std::vector<T> *vmat, std::vector<T> *diag,
                         std::vector<T> *eigv, EigenWork process = EigenWork::EIGENVECTORS);
/// \}

/// \brief Solve the pseudo-inverse of a real matrix A by computing inv(At * A) * At.  The matrix
///        A must have more rows than columns.  If the matrix is well-behaved, this produces the
///        Moore-Penrose Inverse.
///
/// Overloaded:
///   - Overwrite the original matrix A with its pseudo-inverse
///   - Allocate new memory for the pseudo-inverse
///   - Provide pre-allocated C-style arrays with trusted lengths
///   - Provide Standard Template Library vectors
///   - Provide Hybrid objects
///
/// \{
template <typename T>
void pseudoInverse(T* amat, T* workspace, T* inv_workspace, size_t rows, size_t columns);

template <typename T>
void pseudoInverse(T* amat, size_t rows, size_t columns);

template <typename T>
void pseudoInverse(const T* amat, T* result, T* workspace, T* inv_workspace, size_t rows,
                   size_t columns);

template <typename T>
void pseudoInverse(std::vector<T> *amat, std::vector<T> *workspace, std::vector<T> *inv_workspace,
                   size_t rows, size_t columns);

template <typename T>
void pseudoInverse(std::vector<T> *amat, size_t rows, size_t columns);

template <typename T>
void pseudoInverse(const std::vector<T> &amat, std::vector<T> *result, std::vector<T> *workspace,
                   std::vector<T> *inv_workspace, size_t rows, size_t columns);

template <typename T>
std::vector<T> pseudoInverse(const std::vector<T> &amat, size_t rows, size_t columns);

template <typename T>
std::vector<T> pseudoInverse(const std::vector<T> &amat, std::vector<T> *workspace,
                             std::vector<T> *inv_workspace, size_t rows, size_t columns);

template <typename T>
void pseudoInverse(Hybrid<T> *amat, Hybrid<T> *workspace, Hybrid<T> *inv_workspace, size_t rows,
                   size_t columns);

template <typename T>
void pseudoInverse(Hybrid<T> *amat, size_t rows, size_t columns);

template <typename T>
void pseudoInverse(const Hybrid<T> &amat, Hybrid<T> *result, Hybrid<T> *workspace,
                   Hybrid<T> *inv_workspace, size_t rows, size_t columns);
/// \}

/// \brief Solve a system of equations using the QR decomposition.  The matrix shall be presented
///        in column-major format.
///
/// Overloaded:
///   - Provide C-style arrays, all modifiable and with trusted row and column counts, for the
///     stiffness matrix, the solution vector, and the coefficients vector.
///   - Provide Standard Template Library vectors with compatible lengths
///   - Provide Hybrid objects (work will be done by the host CPU with host-side data)
///
/// \param amat    The stiffness matrix
/// \param xvec    The vector of coefficients to be solved
/// \param bvec    The solution vector for the system of equations
/// \param a_rows  Trusted number of rows in amat, and length of bvec
/// \param a_cols  Trusted number of columns in amat, and length of xvec
/// \{
template <typename T>
void qrSolver(T* amat, T* xvec, T* bvec, const size_t a_rows, const size_t a_cols);

template <typename T>
void qrSolver(std::vector<T> *amat, std::vector<T> *xvec, std::vector<T> *bvec);  

template <typename T>
void qrSolver(const std::vector<T> &amat, std::vector<T> *xvec, const std::vector<T> &bvec);

template <typename T>
void qrSolver(Hybrid<T> *amat, Hybrid<T> *xvec, Hybrid<T> *bvec);  

template <typename T>
void qrSolver(const Hybrid<T> &amat, Hybrid<T> *xvec, const Hybrid<T> &bvec);
/// \}
  
/// \brief Compute the determinant of a rank N matrix using the Leibniz formula.
///
/// Overloaded:
///   - Provide a C-style array with N^2 elements (N > 5 will be trapped as inefficient)
///   - Provide a Standard template library vector
///   - Provide a Hybrid object
///
/// \param amat  The matrix of interest, presented in column-major (Fortran) format
/// \param rank  Trusted square root of the length of amat
/// \{
template <typename T>
double leibnizDeterminant(const T* amat, const size_t rank);

template <typename T>
double leibnizDeterminant(const std::vector<T> &amat);

template <typename T>
double leibnizDeterminant(const Hybrid<T> &amat);
/// \}
  
/// \brief Compute the box space transformation matrix given a sequence of six real numbers.
///        Templated for single- and double-precision real forms.
///
/// Overloaded:
///   - Take each dimension as a separate value
///   - Take all dimensions as a const std::vector reference
///   - Take all dimensions as a const C-style array pointer
///
/// \param lx      Length along the first principal axis
/// \param ly      Length along the second principal axis
/// \param lz      Length along the third principal axis
/// \param alpha   First box angle, in radians
/// \param beta    Second box angle, in radians
/// \param gamma   Third box angle, in radians
/// \param matrix  Pre-allocated space for the matrix 
/// \param dims    Vector of all dimensions in the order {lx, ly, lz, alpha, beta, gamma}
/// \{
template <typename T>
void computeBoxTransform(T lx, T ly, T lz, T alpha, T beta, T gamma, T* umat, T* invu);

template <typename T> void computeBoxTransform(const T* dims, T* umat, T* invu);

template <typename T> void computeBoxTransform(const std::vector<T> &dims, std::vector<T> *umat,
                                               std::vector<T> *invu);
/// \}

/// \brief Pull the box dimensions back out of the inverse transformation matrix.
///
/// \param lx     The first box dimension (returned)
/// \param ly     The second box dimension (returned)
/// \param lz     The third box dimension (returned)
/// \param alpha  The first box angle (returned)
/// \param beta   The second box angle (returned)
/// \param gamma  The third box angle (returned)
/// \param invu   Inverse transformation matrix (takes fractional coordinates into real space)
template <typename T>
void extractBoxDimensions(T *lx, T *ly, T *lz, T *alpha, T *beta, T *gamma, const T* invu);

/// \brief Compute the Hessian normal form to obtain the thickness of a unit cell (or mesh element)
///        based on its "a", "b", and "c" box vectors.  The thicknesses between planes defined by
///        the "b" and "c" vectors (x in an orthorhombic unit cell), "a" and "c" vectors (y in an
///        orthorhombic unit cell) and "a" and "b" vectors (z in an orthorhombic unit cell) are
///        returned.
///
/// Overloaded:
///   - Accept inputs as a C-style array or as a Standard Template Library object
///   - Return the result as a vector in the data type used for computations
///   - Return the result in three mutable pointers for each thickness
///
/// \param invu  The 3x3 inverse transformation matrix taking coordinates in unit cell (box) space
///              back to real space.  Given in Fortran order (first three elements are the first
///              column, second three are the second column...).
/// \{
template <typename T> std::vector<T> hessianNormalWidths(const T* invu);

template <typename T> std::vector<T> hessianNormalWidths(const std::vector<T> &invu);

template <typename T> void hessianNormalWidths(const T* invu, T *xw, T *yw, T *zw);

template <typename T> void hessianNormalWidths(const std::vector<T> &invu, T *xw, T *yw, T *zw);
/// \}
  
/// \brief Print a matrix, either to the screen (if no file name is supplied) or to a file.  This
///        makes use of the writeSnapshot() function in 
///
/// Overloaded:
///   - Print a matrix based on a raw pointer to scalar data
///   - Print a matrix based on a std::vector with scalar data
///   - Print a matrix based on a Hybrid object with scalar data
///
/// \param matrix       The matrix to print
/// \param rows         Number of rows in the matrix
/// \param cols         Number of columns in the matrix
/// \param varname      Variable name to give the matrix when printed
/// \param file_name    File name to print the matrix into (stdout if this argument is empty)
/// \param expectation  Indicator of wheter to expect a file of the given name, and what to do
/// \{
template <typename T>
void printMatrix(const T* matrix, const int rows, const int cols, const std::string &varname,
                 const std::string &file_name = std::string(""),
                 const PrintSituation expectation = PrintSituation::APPEND);

template <typename T>
void printMatrix(const std::vector<T> &matrix, const int rows, const int cols,
                 const std::string &varname, const std::string &file_name = std::string(""),
                 const PrintSituation expectation = PrintSituation::APPEND);

template <typename T>
void printMatrix(const Hybrid<T> &matrix, const int rows, const int cols,
                 const std::string &varname, const std::string &file_name = std::string(""),
                 const PrintSituation expectation = PrintSituation::APPEND);

template <typename T>
void printMatrix(const HpcMatrix<T> &matrix, const std::string &varname,
                 const std::string &file_name = std::string(""),
                 const PrintSituation expectation = PrintSituation::APPEND);
/// \}

} // namespace stmath
} // namespace stormm

#include "matrix_ops.tpp"

#endif
