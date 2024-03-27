// -*-c++-*-
#ifndef STORMM_FFT_H
#define STORMM_FFT_H

#include "copyright.h"
#include "Accelerator/hybrid.h"

namespace stormm {
namespace math {

/// \brief Object for managing an FFT, equivalent in some sense to FFTW-wisdom.
template <typename T> class FFTConfig {
public:

  /// \brief The constructor needs a dimension.
  FFTConfig(const int length);

private:
  bool do_inverse;     ///< Flag to have the inverse FFT computed
  int3 length;         ///< Length of the FFT along the X, Y, and Z dimensions
  ullint3 factors;     ///< Bit-packed unsigned long long integer storing up to sixteen FFT factors
                       ///<   along the X, Y, and Z dimensions
  int3 factor_counts;  ///< Numbers of factors along the X, Y, and Z dimensions
  T twiddles;   
};

/// \brief Perform a discrete transform of complex floating-point data with a dimension containing
///        a factor of 2.  The data type must be either float2 or double2.
///
/// \param data
template <typename T> void fftButterfly2(T* data, int fstride, const FFTConfig &state, int m);

/// \brief Function to loop over all factors and manage an FFT over one dimension of the data
  
/// \brief Perform a real-to-complex FFT on a data set of a dimension containing at least one
///        factor of 2, and other factors of 2, 3, 5, or 7.
///
/// \param data    The data to transform
/// \param buffer  Buffer for the out-of-place transpose
/// \param length  The length of the data
template <typename T> void r2c_dft(T* data, T* buffer, int length);

} // namespace math
} // namespace stormm

#include "fft.tpp"

#endif
