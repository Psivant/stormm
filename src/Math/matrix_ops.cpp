#include "copyright.h"
#include "matrix_ops.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
void checkArraySizeToMatrixRank(const size_t actual_rank, const size_t rank, const size_t s_a,
                                const size_t s_v, const ExceptionResponse policy) {
  if (actual_rank == 0LLU) {
    rtErr("A matrix of zero rank is invalid.", "jacobiEigensolver");
  }
  if (actual_rank * actual_rank != s_a || actual_rank * actual_rank != s_v) {
    if (s_a < actual_rank * actual_rank) {
      if (rank == 0LLU) {
        rtErr("Inferred rank " + std::to_string(actual_rank) + " suggests that the input matrix "
              "does not contain enough data (" + std::to_string(s_a) + " elements in the "
              "object).", "jacobiEigensolver");
      }
      else {
        rtErr("The input matrix does not contain enough data (" + std::to_string(s_a) +
              " elements) to be of rank " + std::to_string(rank) + ".", "jacobiEigensolver");
      }
    }
    else if (s_v < actual_rank * actual_rank) {
      rtErr("Eigenvector matrix contains too little space to store eigenvectors of a rank " +
            std::to_string(actual_rank) + " matrix.", "jacobiEigensolver");
    }
    else {
      if (rank == 0LLU) {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("Input matrix or eigenvector space are too large for a rank " +
                std::to_string(actual_rank) + " matrix.", "jacobiEigensolver");
        case ExceptionResponse::WARN:
          rtWarn("Input matrix or eigenvector space are too large for a rank " +
                 std::to_string(actual_rank) + " matrix.", "jacobiEigensolver");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
      }
      else {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("Input matrix or eigenvector space are too large for a matrix of inferred rank " +
                std::to_string(actual_rank) + ".", "jacobiEigensolver");
        case ExceptionResponse::WARN:
          rtWarn("Input matrix or eigenvector space are too large for a matrix of inferred rank " +
                 std::to_string(actual_rank) + ".", "jacobiEigensolver");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
      }
    }
  }
}
  
} // namespace stmath
} // namespace stormm
