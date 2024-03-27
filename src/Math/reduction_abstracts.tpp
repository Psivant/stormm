// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
template <typename T>
GenericRdSubstrate<T>::GenericRdSubstrate(const T* x_read_in, double* x_buffer_in, T* x_write_in,
                                          const int scale_bits_in) :
    dim{1},
    scale_bits{scale_bits_in}, fp_scaling{pow(2.0, scale_bits)}, inv_fp_scaling{1.0 / fp_scaling},
    x_read{x_read_in}, y_read{nullptr}, z_read{nullptr},
    x_read_ovrf{nullptr}, y_read_ovrf{nullptr}, z_read_ovrf{nullptr},
    x_buffer{x_buffer_in}, y_buffer{nullptr}, z_buffer{nullptr},
    x_write{x_write_in}, y_write{nullptr}, z_write{nullptr},
    x_write_ovrf{nullptr}, y_write_ovrf{nullptr}, z_write_ovrf{nullptr}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
GenericRdSubstrate<T>::GenericRdSubstrate(const T* x_read_in, const T* y_read_in,
                                          const T* z_read_in, double* x_buffer_in,
                                          double* y_buffer_in, double* z_buffer_in, T* x_write_in,
                                          T* y_write_in, T* z_write_in, const int scale_bits_in) :
    dim{1 + (y_read_in != nullptr) + (z_read_in != nullptr)},
    scale_bits{scale_bits_in}, fp_scaling{pow(2.0, scale_bits)}, inv_fp_scaling{1.0 / fp_scaling},
    x_read{x_read_in}, y_read{y_read_in}, z_read{z_read_in},
    x_read_ovrf{nullptr}, y_read_ovrf{nullptr}, z_read_ovrf{nullptr},
    x_buffer{x_buffer_in}, y_buffer{y_buffer_in}, z_buffer{z_buffer_in},
    x_write{x_write_in}, y_write{y_write_in}, z_write{z_write_in},
    x_write_ovrf{nullptr}, y_write_ovrf{nullptr}, z_write_ovrf{nullptr}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
GenericRdSubstrate<T>::GenericRdSubstrate(const T* x_read_in, const int* x_read_ovrf_in,
                                          double* x_buffer_in, T* x_write_in, int* x_write_ovrf_in,
                                          const int scale_bits_in) :
    dim{1},
    scale_bits{scale_bits_in}, fp_scaling{pow(2.0, scale_bits)}, inv_fp_scaling{1.0 / fp_scaling},
    x_read{x_read_in}, y_read{nullptr}, z_read{nullptr},
    x_read_ovrf{x_read_ovrf_in}, y_read_ovrf{nullptr}, z_read_ovrf{nullptr},
    x_buffer{x_buffer_in}, y_buffer{nullptr}, z_buffer{nullptr},
    x_write{x_write_in}, y_write{nullptr}, z_write{nullptr},
    x_write_ovrf{x_write_ovrf_in}, y_write_ovrf{nullptr}, z_write_ovrf{nullptr}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
GenericRdSubstrate<T>::GenericRdSubstrate(const T* x_read_in, const int* x_read_ovrf_in,
                                          const T* y_read_in, const int* y_read_ovrf_in,
                                          const T* z_read_in, const int* z_read_ovrf_in,
                                          double* x_buffer_in, double* y_buffer_in,
                                          double* z_buffer_in, T* x_write_in, int* x_write_ovrf_in,
                                          T* y_write_in, int* y_write_ovrf_in, T* z_write_in,
                                          int* z_write_ovrf_in, const int scale_bits_in) :
    dim{1 + (y_read_in != nullptr) + (z_read_in != nullptr)},
    scale_bits{scale_bits_in}, fp_scaling{pow(2.0, scale_bits)}, inv_fp_scaling{1.0 / fp_scaling},
    x_read{x_read_in}, y_read{y_read_in}, z_read{z_read_in},
    x_read_ovrf{x_read_ovrf_in}, y_read_ovrf{y_read_ovrf_in}, z_read_ovrf{z_read_ovrf_in},
    x_buffer{x_buffer_in}, y_buffer{y_buffer_in}, z_buffer{z_buffer_in},
    x_write{x_write_in}, y_write{y_write_in}, z_write{z_write_in},
    x_write_ovrf{x_write_ovrf_in}, y_write_ovrf{y_write_ovrf_in}, z_write_ovrf{z_write_ovrf_in}
{}

} // namespace stmath
} // namespace stormm
