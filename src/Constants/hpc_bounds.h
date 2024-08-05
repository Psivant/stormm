// -*-c++-*-
#ifndef STORMM_HPC_BOUNDS_H
#define STORMM_HPC_BOUNDS_H

#include "copyright.h"
#include "DataTypes/common_types.h"

namespace stormm {
namespace constants {

/// \brief The warp size can be given in many number formats to minimize type conversions during
///        compute-intensive processes
/// \{
#if defined(STORMM_USE_HIP)
constexpr size_t warp_size_zu = 64;
constexpr int warp_size_int = 64;
constexpr unsigned long int warp_size_uint = 64LU;
constexpr long long int warp_size_lld = 64LL;
constexpr unsigned long long int warp_size_llu = 64LLU;
constexpr int warp_bits = 6;
#elif defined(STORMM_USE_INTEL)
constexpr size_t warp_size_zu = 16;
constexpr int warp_size_int = 16;
constexpr unsigned long int warp_size_uint = 16LU;
constexpr long long int warp_size_lld = 16LL;
constexpr unsigned long long int warp_size_llu = 16LLU;
constexpr int warp_bits = 4;
#else
constexpr size_t warp_size_zu = 32;
constexpr int warp_size_int = 32;
constexpr unsigned long int warp_size_uint = 32LU;
constexpr long long int warp_size_lld = 32LL;
constexpr unsigned long long int warp_size_llu = 32LLU;
constexpr int warp_bits = 5;
#endif
constexpr size_t twice_warp_size_zu = 2LLU * warp_size_zu;
constexpr int twice_warp_size_int = 2 * warp_size_int;
constexpr unsigned long int twice_warp_size_uint = 2LU * warp_size_uint;
constexpr long long int twice_warp_size_lld = 2LL * warp_size_lld;
constexpr unsigned long long int twice_warp_size_llu = 2LLU * warp_size_llu;
constexpr size_t warp_bits_mask_zu = warp_size_zu - 1LLU;
constexpr int warp_bits_mask_int = warp_size_int - 1;
constexpr unsigned long int warp_bits_mask_lu = warp_size_uint - 1LU;
constexpr long long int warp_bits_mask_lld = warp_size_lld - 1LL;
constexpr unsigned long long int warp_bits_mask_llu = warp_size_llu - 1LLU;
constexpr size_t twice_warp_bits_mask_zu = warp_bits_mask_zu + warp_size_zu;
constexpr int twice_warp_bits_mask_int = warp_bits_mask_int + warp_size_int;
constexpr unsigned long int twice_warp_bits_mask_lu = warp_bits_mask_lu + warp_size_uint;
constexpr long long int twice_warp_bits_mask_lld = warp_bits_mask_lld + warp_size_lld;
constexpr unsigned long long int twice_warp_bits_mask_llu = warp_bits_mask_llu + warp_size_llu;
constexpr int quarter_warp_size_int = (warp_size_int >> 2);
constexpr int half_warp_size_int = (warp_size_int >> 1);
constexpr int three_quarter_warp_size_int = quarter_warp_size_int + half_warp_size_int;
constexpr int quarter_warp_bits = 3;
constexpr int half_warp_bits = 4;
/// \}

/// \brief Enumerator for general block sizes
enum class BlockSize {
  TINY,    ///< Tiny block of 128 threads
  SMALL,   ///< Small block of 256 threads
  MEDIUM,  ///< Medium-sized block of 512 threads
  LARGE    ///< Largest block of 1024 threads
};

/// \brief Increment for block size tests used in determining the optimal size for kernel launches
constexpr int block_size_test_increment = 64;

/// \brief HPC block sizes, based on a maximum number of 1024 threads per block (common to both
///        CUDA and HIP).  The tiny block size represents a minimum size for any one thread block,
///        and the large block size represents a maximum.  Blocks requiring thread counts other
///        than those presented here are possible, within those limits.
/// \{
constexpr int tiny_block_size   = 128;
constexpr int small_block_size  = 256;
constexpr int medium_block_size = 512;
constexpr int large_block_size  = 1024;
/// \}

/// \brief The number of warps in any given block is a convienient number to have.
/// \{
constexpr int tiny_block_warps   = (tiny_block_size >> warp_bits);
constexpr int small_block_warps  = (small_block_size >> warp_bits);
constexpr int medium_block_warps = (medium_block_size >> warp_bits);
constexpr int large_block_warps  = (large_block_size >> warp_bits);
/// \}
  
} // namespace constants
} // namespace stormm

// Share the major block bounds
namespace stormm {
using constants::warp_size_zu;
using constants::warp_size_int;
using constants::warp_size_uint;
using constants::warp_size_lld;
using constants::warp_size_llu;
using constants::twice_warp_size_zu;
using constants::twice_warp_size_int;
using constants::twice_warp_size_uint;
using constants::twice_warp_size_lld;
using constants::twice_warp_size_llu;
using constants::twice_warp_bits_mask_zu;
using constants::twice_warp_bits_mask_int;
using constants::twice_warp_bits_mask_lu;
using constants::twice_warp_bits_mask_lld;
using constants::twice_warp_bits_mask_llu;
using constants::warp_bits;
using constants::warp_bits_mask_int;
using constants::warp_bits_mask_lu;
using constants::warp_bits_mask_lld;
using constants::warp_bits_mask_llu;
using constants::warp_bits_mask_zu;
using constants::quarter_warp_size_int;
using constants::half_warp_size_int;
using constants::three_quarter_warp_size_int;
using constants::quarter_warp_bits;
using constants::half_warp_bits;
using constants::BlockSize;
using constants::large_block_size;
using constants::medium_block_size;
using constants::small_block_size;
using constants::tiny_block_size;
using constants::large_block_warps;
using constants::medium_block_warps;
using constants::small_block_warps;
using constants::tiny_block_warps;
}

#endif
