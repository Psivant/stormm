// -*-c++-*-
#include "copyright.h"
#include "Accelerator/ptx_macros.h"
#include "Constants/hpc_bounds.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Numerics/split_fixed_precision.h"
#include "Potential/energy_abstracts.h"
#include "Random/random.h"
#include "Synthesis/implicit_solvent_workspace.h"
#include "Synthesis/nonbonded_workunit.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Topology/atomgraph_enumerators.h"
#include "hpc_nonbonded_potential.h"

namespace stormm {
namespace energy {

using constants::PrecisionModel;
using constants::twice_warp_bits_mask_int;
using constants::twice_warp_size_int;
using constants::warp_size_int;
using constants::warp_bits;
using constants::warp_bits_mask_int;
using numerics::chooseAccumulationMethod;
using synthesis::NbwuKind;
using synthesis::small_block_max_imports;
using synthesis::small_block_max_atoms;
using synthesis::tile_groups_wu_abstract_length;
using topology::ImplicitSolventModel;
using namespace random;
  
//-------------------------------------------------------------------------------------------------
// Get the number atoms in a particular tile stretch.
//
// Arguments:
//   nbwu_map:  Details of the non-bonded work unit, condensed into a simple array of integers
//   pos:       Thread position in the list of atoms
//-------------------------------------------------------------------------------------------------
__device__ __forceinline__ int getTileSideAtomCount(const int* nbwu_map, const int pos) {
  const int key_idx  = pos / 4;
  const int key_slot = pos - (key_idx * 4);
  return ((nbwu_map[small_block_max_imports + 1 + key_idx] >> (8 * key_slot)) & 0xff);
}

#include "Math/rounding.cui"
#include "Numerics/accumulation.cui"
#include "Random/xor_shift_rng.cui"

//-------------------------------------------------------------------------------------------------
// Load coordinates relating to atoms in a non-bonded tile suitable for isolated boundary
// conditions.
//
// Overloaded:
//   - Work with an array of long long integers appropriate for single-precision arithmetic
//   - Work with dual arrays of long long int and int types, appropriate for double-precision
//     arithmetic
//
// Arguments:
//   pos:              Position in the tile list (not the atom list)
//   import_count:     Number of groups of atoms imported to populate one side of one or more tiles
//   iter:             Number of passes made by this or related routines (incrementation of
//                     iter is essential to maintain the correct procession through all loads)
//   nbwu_map:         Non-bonded work unit details            
//   read_crd:         Array of coordinates to read from
//   write_crd:        Array of coordinates to write into
//   read_crd_ovrf:    Overflow buffers for coordinates to be read
//   write_crd_ovrf:   Overflow buffers for local copies of coordinates
//   sh_tile_cog:      Array holding mean values of the positions of each imported atom group (this
//                     later expedites computing the center of geometry for each complete tile)
//   gpos_scale:       Scaling factor for coordinates in the fixed-precision representation (this
//                     is needed only to place dummy atom coordinates for blank slots of a tile)
//-------------------------------------------------------------------------------------------------
__device__ int loadTileCoordinates(const int pos, const int iter, const int* nbwu_map,
                                   const llint* read_crd, llint* write_crd, float* sh_tile_cog,
                                   const float gpos_scale) {
  const int tile_sides_per_warp = (warp_size_int / tile_length);
  const int warps_per_block = blockDim.x >> warp_bits;
  const int tile_lane_idx = (threadIdx.x & tile_length_bits_mask);
  const int import_count = nbwu_map[0];
  const int padded_import_count = devcRoundUp(import_count, tile_sides_per_warp);
  int rel_pos = pos - (iter * padded_import_count);
  while (rel_pos < padded_import_count) {
    float fval;
    if (rel_pos < import_count) {
      const size_t read_idx = nbwu_map[rel_pos + 1] + tile_lane_idx;
      const size_t write_idx = (rel_pos * tile_length) + tile_lane_idx;
      if (tile_lane_idx < getTileSideAtomCount(nbwu_map, rel_pos)) {
        const llint ival = __ldcs(&read_crd[read_idx]);
        fval = (float)(ival);
        write_crd[write_idx] = ival;
      }
      else {
        fval = (float)(0.0);
        write_crd[write_idx] = (128 * (rel_pos + 8) * tile_lane_idx) * gpos_scale;
      }
    }
    else {
      fval = (float)(0.0);
    }
    for (int i = half_tile_length; i > 0; i >>= 1) {
      fval += SHFL_DOWN(fval, i);
    }
    if (tile_lane_idx == 0 && rel_pos < import_count) {
      sh_tile_cog[rel_pos] = fval;
    }
    rel_pos += tile_sides_per_warp * warps_per_block;
  }
  return rel_pos + (iter * padded_import_count);
}

__device__ int loadTileCoordinates(const int pos, const int iter, const int* nbwu_map,
                                   const llint* read_crd, llint* write_crd,
                                   const int* read_crd_ovrf, int* write_crd_ovrf,
                                   double* sh_tile_cog, const double gpos_scale) {
  const int tile_sides_per_warp = (warp_size_int / tile_length);
  const int warps_per_block = blockDim.x >> warp_bits;
  const int tile_lane_idx = (threadIdx.x & tile_length_bits_mask);
  const int import_count = nbwu_map[0];
  const int padded_import_count = devcRoundUp(import_count, tile_sides_per_warp);
  int rel_pos = pos - (iter * padded_import_count);
  while (rel_pos < padded_import_count) {
    double fval;
    if (rel_pos < import_count) {
      const size_t read_idx = nbwu_map[rel_pos + 1] + tile_lane_idx;
      const size_t write_idx = (rel_pos * tile_length) + tile_lane_idx;
      if (tile_lane_idx < getTileSideAtomCount(nbwu_map, rel_pos)) {
        const llint ival = __ldcs(&read_crd[read_idx]);
        fval = (double)(ival);
        write_crd[write_idx] = ival;
        const int ival_ovrf = __ldcs(&read_crd_ovrf[read_idx]);
        fval += (double)(ival_ovrf) * max_llint_accumulation;
        write_crd_ovrf[write_idx] = ival_ovrf;
      }
      else {
        fval = 0.0;
        const int95_t fake_val = doubleToInt95((128 * (rel_pos + 8) * tile_lane_idx) * gpos_scale);
        write_crd[write_idx] = fake_val.x;
        write_crd_ovrf[write_idx] = fake_val.y;
      }
    }
    else {
      fval = 0.0;
    }
    for (int i = half_tile_length; i > 0; i >>= 1) {
      fval += SHFL_DOWN(fval, i);
    }
    if (tile_lane_idx == 0 && rel_pos < import_count) {
      sh_tile_cog[rel_pos] = fval;
    }
    rel_pos += tile_sides_per_warp * warps_per_block;
  }
  return rel_pos + (iter * padded_import_count);
}

//-------------------------------------------------------------------------------------------------
// Load scalar values (integral or real) from global memory into local arrays, on a tile-by-tile
// basis for non-bonded work units involving isolated systems with all-to-all interaction matrices.
//
// Overloaded:
//   - Copy the values directly
//   - Fold in a scalar multiple
//   - Fold in a scalar addition
//
// Parameter descriptors follow from loadTileCoordinates() above, with alterations:
//   read_array:   Generic array of (global) information to read from
//   write_array:  Generic (local) array of information to write 
//-------------------------------------------------------------------------------------------------
template <typename T> __device__
int loadTileProperty(const int pos, const int iter, const int* nbwu_map, const T* read_array,
                     T* write_array) {
  const int tile_sides_per_warp = (warp_size_int / tile_length);
  const int warps_per_block = blockDim.x >> warp_bits;
  const int tile_lane_idx = (threadIdx.x & tile_length_bits_mask);
  const int import_count = nbwu_map[0];
  const int padded_import_count = devcRoundUp(import_count, tile_sides_per_warp);  
  int rel_pos = pos - (iter * padded_import_count);
  while (rel_pos < padded_import_count) {
    if (rel_pos < import_count) {
      const size_t read_idx = nbwu_map[rel_pos + 1] + tile_lane_idx;
      const size_t write_idx = (rel_pos * tile_length) + tile_lane_idx;
      if (tile_lane_idx < getTileSideAtomCount(nbwu_map, rel_pos)) {
        write_array[write_idx] = __ldcs(&read_array[read_idx]);
      }
      else {
        write_array[write_idx] = (T)(0);
      }
    }
    rel_pos += tile_sides_per_warp * warps_per_block;
  }
  return rel_pos + (iter * padded_import_count);
}

template <typename T> __device__
int loadTileProperty(const int pos, const int iter, const int* nbwu_map, const T* read_array,
                     T* write_array, T multiplier) {
  const int tile_sides_per_warp = (warp_size_int / tile_length);
  const int warps_per_block = blockDim.x >> warp_bits;
  const int tile_lane_idx = (threadIdx.x & tile_length_bits_mask);
  const int import_count = nbwu_map[0];
  const int padded_import_count = devcRoundUp(import_count, tile_sides_per_warp);  
  int rel_pos = pos - (iter * padded_import_count);
  while (rel_pos < padded_import_count) {
    if (rel_pos < import_count) {
      const size_t read_idx = nbwu_map[rel_pos + 1] + tile_lane_idx;
      const size_t write_idx = (rel_pos * tile_length) + tile_lane_idx;
      if (tile_lane_idx < getTileSideAtomCount(nbwu_map, rel_pos)) {
        write_array[write_idx] = __ldcs(&read_array[read_idx]) * multiplier;
      }
      else {
        write_array[write_idx] = (T)(0);
      }
    }
    rel_pos += tile_sides_per_warp * warps_per_block;
  }
  return rel_pos + (iter * padded_import_count);
}

template <typename T> __device__
int loadTileProperty(const int pos, const int iter, const int* nbwu_map, const T* read_array,
                     T increment, T* write_array) {
  const int tile_sides_per_warp = (warp_size_int / tile_length);
  const int warps_per_block = blockDim.x >> warp_bits;
  const int tile_lane_idx = (threadIdx.x & tile_length_bits_mask);
  const int import_count = nbwu_map[0];
  const int padded_import_count = devcRoundUp(import_count, tile_sides_per_warp);  
  int rel_pos = pos - (iter * padded_import_count);
  while (rel_pos < padded_import_count) {
    if (rel_pos < import_count) {
      const size_t read_idx = nbwu_map[rel_pos + 1] + tile_lane_idx;
      const size_t write_idx = (rel_pos * tile_length) + tile_lane_idx;
      if (tile_lane_idx < getTileSideAtomCount(nbwu_map, rel_pos)) {
        write_array[write_idx] = __ldcs(&read_array[read_idx]) + increment;
      }
      else {
        write_array[write_idx] = (T)(0);
      }
    }
    rel_pos += tile_sides_per_warp * warps_per_block;
  }
  return rel_pos + (iter * padded_import_count);
}

//-------------------------------------------------------------------------------------------------
// Write information about the atoms in tile groups back to global accumulators.  Relevant for
// systems with all-to-all interactions in isolated boundary conditions.
//
// Overloaded:
//   - Accept various combinations of single- or double-integer local accumulators to contribute
//     to the implied single- or double-integer global accumulators
//
// Arguments:
//   pos:                   Position in the tile list (not the atom list)
//   iter:                  Number of passes made by this or related routines (incrementation of
//                          iter is essential to maintain the correct procession through all loads)
//   nbwu_map:              Non-bonded work unit details
//   tile_prop:             Primary local accumulator for the tile-based computed property
//   tile_prop_ovrf:        Local overflow accumulator for the tile-based computed property
//   gbl_accumulator:       Primary (or, perhaps lone) global accumulator for the computed property
//   gbl_accumulator_ovrf:  Overflow global accumulator for the computed property
//-------------------------------------------------------------------------------------------------
__device__ int accumulateTileProperty(const int pos, const int iter, const int* nbwu_map,
                                      const int* tile_prop, const int* tile_prop_ovrf,
                                      llint* gbl_accumulator) {
  const int tile_sides_per_warp = (warp_size_int / tile_length);
  const int warps_per_block = blockDim.x >> warp_bits;
  const int tile_lane_idx = (threadIdx.x & tile_length_bits_mask);
  const int import_count = nbwu_map[0];
  const int padded_import_count = devcRoundUp(import_count, tile_sides_per_warp);
  int rel_pos = pos - (iter * padded_import_count);
  while (rel_pos < padded_import_count) {
    if (rel_pos < import_count) {
      const size_t write_idx = nbwu_map[rel_pos + 1] + tile_lane_idx;
      const size_t read_idx = (rel_pos * tile_length) + tile_lane_idx;
      if (tile_lane_idx < getTileSideAtomCount(nbwu_map, rel_pos)) {
        llint itp = tile_prop_ovrf[read_idx];
        itp *= max_int_accumulation_ll;
        itp += tile_prop[read_idx];
        atomicAdd((ullint*)&gbl_accumulator[write_idx], (ullint)(itp));
      }
    }
    rel_pos += tile_sides_per_warp * warps_per_block;
  }
  return rel_pos + (iter * padded_import_count);
}

__device__ int accumulateTileProperty(const int pos, const int iter, const int* nbwu_map,
                                      const llint* tile_prop, llint* gbl_accumulator) {
  const int tile_sides_per_warp = (warp_size_int / tile_length);
  const int warps_per_block = blockDim.x >> warp_bits;
  const int tile_lane_idx = (threadIdx.x & tile_length_bits_mask);
  const int import_count = nbwu_map[0];
  const int padded_import_count = devcRoundUp(import_count, tile_sides_per_warp);
  int rel_pos = pos - (iter * padded_import_count);
  while (rel_pos < padded_import_count) {
    if (rel_pos < import_count) {
      const size_t write_idx = nbwu_map[rel_pos + 1] + tile_lane_idx;
      const size_t read_idx = (rel_pos * tile_length) + tile_lane_idx;
      if (tile_lane_idx < getTileSideAtomCount(nbwu_map, rel_pos)) {
        atomicAdd((ullint*)&gbl_accumulator[write_idx], (ullint)(tile_prop[read_idx]));
      }
    }
    rel_pos += tile_sides_per_warp * warps_per_block;
  }
  return rel_pos + (iter * padded_import_count);
}

__device__ int accumulateTileProperty(const int pos, const int iter, const int* nbwu_map,
                                      const llint* tile_prop, const int* tile_prop_ovrf,
                                      llint* gbl_accumulator, int* gbl_accumulator_ovrf) {
  const int tile_sides_per_warp = (warp_size_int / tile_length);
  const int warps_per_block = blockDim.x >> warp_bits;
  const int tile_lane_idx = (threadIdx.x & tile_length_bits_mask);
  const int import_count = nbwu_map[0];
  const int padded_import_count = devcRoundUp(import_count, tile_sides_per_warp);
  int rel_pos = pos - (iter * padded_import_count);
  while (rel_pos < padded_import_count) {
    if (rel_pos < import_count) {
      const size_t write_idx = nbwu_map[rel_pos + 1] + tile_lane_idx;
      const size_t read_idx = (rel_pos * tile_length) + tile_lane_idx;
      if (tile_lane_idx < getTileSideAtomCount(nbwu_map, rel_pos)) {
        atomicSplit(tile_prop[read_idx], tile_prop_ovrf[read_idx], write_idx, gbl_accumulator,
                    gbl_accumulator_ovrf);
      }
    }
    rel_pos += tile_sides_per_warp * warps_per_block;
  }
  return rel_pos + (iter * padded_import_count);
}

// Single-precision non-bonded kernel floating point definitions
#define TCALC float
#  define TCALC2 float2
#  define TCALC_IS_SINGLE
#  ifdef STORMM_USE_CUDA
#    if (__CUDA_ARCH__ >= 750) && (__CUDA_ARCH__ < 800)
#      define NONBOND_KERNEL_BLOCKS_MULTIPLIER 4
#      define GBRADII_KERNEL_BLOCKS_MULTIPLIER 4
#      define GBDERIV_KERNEL_BLOCKS_MULTIPLIER 4
#    else
#      define NONBOND_KERNEL_BLOCKS_MULTIPLIER 5
#      define GBRADII_KERNEL_BLOCKS_MULTIPLIER 5
#      define GBDERIV_KERNEL_BLOCKS_MULTIPLIER 5
#    endif
#  else
#    define NONBOND_KERNEL_BLOCKS_MULTIPLIER 5
#    define GBRADII_KERNEL_BLOCKS_MULTIPLIER 5
#    define GBDERIV_KERNEL_BLOCKS_MULTIPLIER 5
#  endif
#  define LLCONV_FUNC __float2ll_rn
#  define SQRT_FUNC sqrtf
#  define CBRT_FUNC cbrtf
#  define LOG_FUNC  logf
#  define EXP_FUNC  expf
#  define TANH_FUNC tanhf
#  define FABS_FUNC fabsf
#  define SPLIT_FORCE_ACCUMULATION
#    define KERNEL_NAME ktgfsCalculateGBRadii
#      include "gbradii_tilegroups.cui"
#    undef KERNEL_NAME
#    define DO_NECK_CORRECTION
#      define KERNEL_NAME ktgfsCalculateGBNeckRadii
#        include "gbradii_tilegroups.cui"
#      undef KERNEL_NAME
#    undef DO_NECK_CORRECTION
#  undef SPLIT_FORCE_ACCUMULATION
#  define KERNEL_NAME ktgfCalculateGBRadii
#    include "gbradii_tilegroups.cui"
#  undef KERNEL_NAME
#  define DO_NECK_CORRECTION
#    define KERNEL_NAME ktgfCalculateGBNeckRadii
#      include "gbradii_tilegroups.cui"
#    undef KERNEL_NAME
#  undef DO_NECK_CORRECTION
#  define COMPUTE_FORCE
#    define SPLIT_FORCE_ACCUMULATION
#      define COMPUTE_ENERGY
#        define KERNEL_NAME ktgfsVacuumForceEnergy
#          include "nonbonded_potential_tilegroups.cui"
#        undef KERNEL_NAME
#        define DO_GENERALIZED_BORN
#          define KERNEL_NAME ktgfsGBForceEnergy
#            include "nonbonded_potential_tilegroups.cui"
#          undef KERNEL_NAME
#          define DO_NECK_CORRECTION
#            define KERNEL_NAME ktgfsGBNeckForceEnergy
#              include "nonbonded_potential_tilegroups.cui"
#            undef KERNEL_NAME
#          undef DO_NECK_CORRECTION
#        undef DO_GENERALIZED_BORN
#      undef COMPUTE_ENERGY
#      define KERNEL_NAME ktgfsVacuumForce
#        include "nonbonded_potential_tilegroups.cui"
#      undef KERNEL_NAME
#      define DO_GENERALIZED_BORN
#        define KERNEL_NAME ktgfsGBForce
#          include "nonbonded_potential_tilegroups.cui"
#        undef KERNEL_NAME
#        define KERNEL_NAME ktgfsCalculateGBDerivatives
#          include "gbderivative_tilegroups.cui"
#        undef KERNEL_NAME
#        define DO_NECK_CORRECTION
#          define KERNEL_NAME ktgfsGBNeckForce
#            include "nonbonded_potential_tilegroups.cui"
#          undef KERNEL_NAME
#          define KERNEL_NAME ktgfsCalculateGBNeckDerivatives
#            include "gbderivative_tilegroups.cui"
#          undef KERNEL_NAME
#        undef DO_NECK_CORRECTION
#      undef DO_GENERALIZED_BORN
#    undef SPLIT_FORCE_ACCUMULATION
#    define COMPUTE_ENERGY
#      define KERNEL_NAME ktgfVacuumForceEnergy
#        include "nonbonded_potential_tilegroups.cui"
#      undef KERNEL_NAME
#      define DO_GENERALIZED_BORN
#        define KERNEL_NAME ktgfGBForceEnergy
#          include "nonbonded_potential_tilegroups.cui"
#        undef KERNEL_NAME
#        define DO_NECK_CORRECTION
#          define KERNEL_NAME ktgfGBNeckForceEnergy
#            include "nonbonded_potential_tilegroups.cui"
#          undef KERNEL_NAME
#        undef DO_NECK_CORRECTION
#      undef DO_GENERALIZED_BORN
#    undef COMPUTE_ENERGY
#    define KERNEL_NAME ktgfVacuumForce
#      include "nonbonded_potential_tilegroups.cui"
#    undef KERNEL_NAME
#    define DO_GENERALIZED_BORN
#      define KERNEL_NAME ktgfGBForce
#        include "nonbonded_potential_tilegroups.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME ktgfCalculateGBDerivatives
#        include "gbderivative_tilegroups.cui"
#      undef KERNEL_NAME
#      define DO_NECK_CORRECTION
#        define KERNEL_NAME ktgfGBNeckForce
#          include "nonbonded_potential_tilegroups.cui"
#        undef KERNEL_NAME
#        define KERNEL_NAME ktgfCalculateGBNeckDerivatives
#          include "gbderivative_tilegroups.cui"
#        undef KERNEL_NAME
#      undef DO_NECK_CORRECTION
#    undef DO_GENERALIZED_BORN
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define KERNEL_NAME ktgfVacuumEnergy
#      include "nonbonded_potential_tilegroups.cui"
#    undef KERNEL_NAME
#    define DO_GENERALIZED_BORN
#      define KERNEL_NAME ktgfGBEnergy
#        include "nonbonded_potential_tilegroups.cui"
#      undef KERNEL_NAME
#      define DO_NECK_CORRECTION
#        define KERNEL_NAME ktgfGBNeckEnergy
#          include "nonbonded_potential_tilegroups.cui"
#        undef KERNEL_NAME
#      undef DO_NECK_CORRECTION
#    undef DO_GENERALIZED_BORN
#  undef COMPUTE_ENERGY
#  undef GBRADII_KERNEL_BLOCKS_MULTIPLIER
#  undef GBDERIV_KERNEL_BLOCKS_MULTIPLIER
#  undef NONBOND_KERNEL_BLOCKS_MULTIPLIER
#  define CLASH_FORGIVENESS
#    define COMPUTE_FORCE
#      define SPLIT_FORCE_ACCUMULATION
#        define COMPUTE_ENERGY
#          define NONBOND_KERNEL_BLOCKS_MULTIPLIER 4
#          define KERNEL_NAME ktgfsVacuumForceEnergyNonClash
#            include "nonbonded_potential_tilegroups.cui"
#          undef KERNEL_NAME
#          define DO_GENERALIZED_BORN
#            define KERNEL_NAME ktgfsGBForceEnergyNonClash
#              include "nonbonded_potential_tilegroups.cui"
#            undef KERNEL_NAME
#            define DO_NECK_CORRECTION
#              define KERNEL_NAME ktgfsGBNeckForceEnergyNonClash
#                include "nonbonded_potential_tilegroups.cui"
#              undef KERNEL_NAME
#            undef DO_NECK_CORRECTION
#          undef DO_GENERALIZED_BORN
#          undef NONBOND_KERNEL_BLOCKS_MULTIPLIER
#        undef COMPUTE_ENERGY
#        define NONBOND_KERNEL_BLOCKS_MULTIPLIER 5
#        define KERNEL_NAME ktgfsVacuumForceNonClash
#          include "nonbonded_potential_tilegroups.cui"
#        undef KERNEL_NAME
#        define DO_GENERALIZED_BORN
#          define KERNEL_NAME ktgfsGBForceNonClash
#            include "nonbonded_potential_tilegroups.cui"
#          undef KERNEL_NAME
#          define DO_NECK_CORRECTION
#            define KERNEL_NAME ktgfsGBNeckForceNonClash
#              include "nonbonded_potential_tilegroups.cui"
#            undef KERNEL_NAME
#          undef DO_NECK_CORRECTION
#        undef DO_GENERALIZED_BORN
#        undef NONBOND_KERNEL_BLOCKS_MULTIPLIER
#      undef SPLIT_FORCE_ACCUMULATION
#      define COMPUTE_ENERGY
#        define NONBOND_KERNEL_BLOCKS_MULTIPLIER 4
#        define KERNEL_NAME ktgfVacuumForceEnergyNonClash
#          include "nonbonded_potential_tilegroups.cui"
#        undef KERNEL_NAME
#        define DO_GENERALIZED_BORN
#          define KERNEL_NAME ktgfGBForceEnergyNonClash
#            include "nonbonded_potential_tilegroups.cui"
#          undef KERNEL_NAME
#          define DO_NECK_CORRECTION
#            define KERNEL_NAME ktgfGBNeckForceEnergyNonClash
#              include "nonbonded_potential_tilegroups.cui"
#            undef KERNEL_NAME
#          undef DO_NECK_CORRECTION
#        undef DO_GENERALIZED_BORN
#        undef NONBOND_KERNEL_BLOCKS_MULTIPLIER
#      undef COMPUTE_ENERGY
#      define NONBOND_KERNEL_BLOCKS_MULTIPLIER 5
#      define KERNEL_NAME ktgfVacuumForceNonClash
#        include "nonbonded_potential_tilegroups.cui"
#      undef KERNEL_NAME
#      define DO_GENERALIZED_BORN
#        define KERNEL_NAME ktgfGBForceNonClash
#          include "nonbonded_potential_tilegroups.cui"
#        undef KERNEL_NAME
#        define DO_NECK_CORRECTION
#          define KERNEL_NAME ktgfGBNeckForceNonClash
#            include "nonbonded_potential_tilegroups.cui"
#          undef KERNEL_NAME
#        undef DO_NECK_CORRECTION
#      undef DO_GENERALIZED_BORN
#      undef NONBOND_KERNEL_BLOCKS_MULTIPLIER
#    undef COMPUTE_FORCE
#    define COMPUTE_ENERGY
#      define NONBOND_KERNEL_BLOCKS_MULTIPLIER 5
#      define KERNEL_NAME ktgfVacuumEnergyNonClash
#        include "nonbonded_potential_tilegroups.cui"
#      undef KERNEL_NAME
#      define DO_GENERALIZED_BORN
#        define KERNEL_NAME ktgfGBEnergyNonClash
#          include "nonbonded_potential_tilegroups.cui"
#        undef KERNEL_NAME
#        define DO_NECK_CORRECTION
#          define KERNEL_NAME ktgfGBNeckEnergyNonClash
#            include "nonbonded_potential_tilegroups.cui"
#          undef KERNEL_NAME
#        undef DO_NECK_CORRECTION
#      undef DO_GENERALIZED_BORN
#      undef NONBOND_KERNEL_BLOCKS_MULTIPLIER
#    undef COMPUTE_ENERGY
#  undef CLASH_FORGIVENESS
#  undef LLCONV_FUNC
#  undef SQRT_FUNC
#  undef CBRT_FUNC
#  undef LOG_FUNC
#  undef EXP_FUNC
#  undef TANH_FUNC
#  undef FABS_FUNC
#  undef TCALC_IS_SINGLE
#  undef TCALC2
#undef TCALC

// Double-precision non-bonded kernel floating point definitions
#define TCALC double
#  define TCALC2 double2
#  define SPLIT_FORCE_ACCUMULATION
#  define NONBOND_KERNEL_BLOCKS_MULTIPLIER 3
#  define GBRADII_KERNEL_BLOCKS_MULTIPLIER 3
#  define GBDERIV_KERNEL_BLOCKS_MULTIPLIER 3
#  define LLCONV_FUNC __double2ll_rn
#  define SQRT_FUNC sqrt
#  define CBRT_FUNC cbrt
#  define LOG_FUNC  log
#  define EXP_FUNC  exp
#  define TANH_FUNC tanh
#  define FABS_FUNC fabs
#  define KERNEL_NAME ktgdsCalculateGBRadii
#    include "gbradii_tilegroups.cui"
#  undef KERNEL_NAME
#  define DO_NECK_CORRECTION
#    define KERNEL_NAME ktgdsCalculateGBNeckRadii
#      include "gbradii_tilegroups.cui"
#    undef KERNEL_NAME
#  undef DO_NECK_CORRECTION
#  define COMPUTE_FORCE
#    define COMPUTE_ENERGY
#      define KERNEL_NAME ktgdsVacuumForceEnergy
#        include "nonbonded_potential_tilegroups.cui"
#      undef KERNEL_NAME
#      define DO_GENERALIZED_BORN
#        define KERNEL_NAME ktgdsGBForceEnergy
#          include "nonbonded_potential_tilegroups.cui"
#        undef KERNEL_NAME
#        define KERNEL_NAME ktgdsCalculateGBDerivatives
#          include "gbderivative_tilegroups.cui"
#        undef KERNEL_NAME
#        define DO_NECK_CORRECTION
#          define KERNEL_NAME ktgdsGBNeckForceEnergy
#            include "nonbonded_potential_tilegroups.cui"
#          undef KERNEL_NAME
#          define KERNEL_NAME ktgdsCalculateGBNeckDerivatives
#            include "gbderivative_tilegroups.cui"
#          undef KERNEL_NAME
#        undef DO_NECK_CORRECTION
#      undef DO_GENERALIZED_BORN
#    undef COMPUTE_ENERGY
#    define KERNEL_NAME ktgdsVacuumForce
#      include "nonbonded_potential_tilegroups.cui"
#    undef KERNEL_NAME
#    define DO_GENERALIZED_BORN
#      define KERNEL_NAME ktgdsGBForce
#        include "nonbonded_potential_tilegroups.cui"
#      undef KERNEL_NAME
#      define DO_NECK_CORRECTION
#        define KERNEL_NAME ktgdsGBNeckForce
#          include "nonbonded_potential_tilegroups.cui"
#        undef KERNEL_NAME
#      undef DO_NECK_CORRECTION
#    undef DO_GENERALIZED_BORN
#  undef COMPUTE_FORCE
#  define COMPUTE_ENERGY
#    define KERNEL_NAME ktgdVacuumEnergy
#      include "nonbonded_potential_tilegroups.cui"
#    undef KERNEL_NAME
#    define DO_GENERALIZED_BORN
#      define KERNEL_NAME ktgdGBEnergy
#        include "nonbonded_potential_tilegroups.cui"
#      undef KERNEL_NAME
#      define DO_NECK_CORRECTION
#        define KERNEL_NAME ktgdGBNeckEnergy
#          include "nonbonded_potential_tilegroups.cui"
#        undef KERNEL_NAME
#      undef DO_NECK_CORRECTION
#    undef DO_GENERALIZED_BORN
#  undef COMPUTE_ENERGY
#  undef GBRADII_KERNEL_BLOCKS_MULTIPLIER
#  undef GBDERIV_KERNEL_BLOCKS_MULTIPLIER
#  undef NONBOND_KERNEL_BLOCKS_MULTIPLIER
#  define CLASH_FORGIVENESS
#    define COMPUTE_FORCE
#      define COMPUTE_ENERGY
#        define NONBOND_KERNEL_BLOCKS_MULTIPLIER 3
#        define KERNEL_NAME ktgdsVacuumForceEnergyNonClash
#          include "nonbonded_potential_tilegroups.cui"
#        undef KERNEL_NAME
#        undef NONBOND_KERNEL_BLOCKS_MULTIPLIER
#        define DO_GENERALIZED_BORN
#          define NONBOND_KERNEL_BLOCKS_MULTIPLIER 2
#          define KERNEL_NAME ktgdsGBForceEnergyNonClash
#            include "nonbonded_potential_tilegroups.cui"
#          undef KERNEL_NAME
#          define DO_NECK_CORRECTION
#            define KERNEL_NAME ktgdsGBNeckForceEnergyNonClash
#              include "nonbonded_potential_tilegroups.cui"
#            undef KERNEL_NAME
#          undef DO_NECK_CORRECTION
#          undef NONBOND_KERNEL_BLOCKS_MULTIPLIER
#        undef DO_GENERALIZED_BORN
#      undef COMPUTE_ENERGY
#      define NONBOND_KERNEL_BLOCKS_MULTIPLIER 3
#      define KERNEL_NAME ktgdsVacuumForceNonClash
#        include "nonbonded_potential_tilegroups.cui"
#      undef KERNEL_NAME
#      define DO_GENERALIZED_BORN
#        define KERNEL_NAME ktgdsGBForceNonClash
#          include "nonbonded_potential_tilegroups.cui"
#        undef KERNEL_NAME
#        define DO_NECK_CORRECTION
#          define KERNEL_NAME ktgdsGBNeckForceNonClash
#            include "nonbonded_potential_tilegroups.cui"
#          undef KERNEL_NAME
#        undef DO_NECK_CORRECTION
#      undef DO_GENERALIZED_BORN
#      undef NONBOND_KERNEL_BLOCKS_MULTIPLIER
#    undef COMPUTE_FORCE
#    define COMPUTE_ENERGY
#      define NONBOND_KERNEL_BLOCKS_MULTIPLIER 3
#      define KERNEL_NAME ktgdVacuumEnergyNonClash
#        include "nonbonded_potential_tilegroups.cui"
#      undef KERNEL_NAME
#      define DO_GENERALIZED_BORN
#        define KERNEL_NAME ktgdGBEnergyNonClash
#          include "nonbonded_potential_tilegroups.cui"
#        undef KERNEL_NAME
#        define DO_NECK_CORRECTION
#          define KERNEL_NAME ktgdGBNeckEnergyNonClash
#            include "nonbonded_potential_tilegroups.cui"
#          undef KERNEL_NAME
#        undef DO_NECK_CORRECTION
#      undef DO_GENERALIZED_BORN
#      undef NONBOND_KERNEL_BLOCKS_MULTIPLIER
#    undef COMPUTE_ENERGY
#  undef CLASH_FORGIVENESS
#  undef LLCONV_FUNC
#  undef SQRT_FUNC
#  undef CBRT_FUNC
#  undef LOG_FUNC
#  undef EXP_FUNC
#  undef TANH_FUNC
#  undef FABS_FUNC
#  undef SPLIT_FORCE_ACCUMULATION
#  undef TCALC2
#undef TCALC

//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes
queryNonbondedKernelRequirements(const PrecisionModel prec, const NbwuKind kind,
                                 const EvaluateForce eval_frc, const EvaluateEnergy eval_nrg,
                                 const AccumulationMethod acc_meth,
                                 const ImplicitSolventModel igb,
                                 const ClashResponse collision_handling) {
  
  // The kernel manager will have information about the GPU to use--look at the work units from
  // the perspective of overall occupancy on the GPU.
  cudaFuncAttributes attr;
  switch (collision_handling) {
  case ClashResponse::NONE:
    switch (igb) {
    case ImplicitSolventModel::NONE:
      switch (prec) {
      case PrecisionModel::DOUBLE:
        switch (eval_frc) {
        case EvaluateForce::YES:
          switch (eval_nrg) {
          case EvaluateEnergy::YES:
            if (cudaFuncGetAttributes(&attr, ktgdsVacuumForceEnergy) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgdsVacuumForceEnergy.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          case EvaluateEnergy::NO:
            if (cudaFuncGetAttributes(&attr, ktgdsVacuumForce) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgdsVacuumForce.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          }
          break;
        case EvaluateForce::NO:
          if (cudaFuncGetAttributes(&attr, ktgdVacuumEnergy) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel ktgdVacuumEnergy.",
                  "queryNonbondedKernelRequirements");
          }
          break;
        }
        break;
      case PrecisionModel::SINGLE:
        switch (kind) {
          case NbwuKind::TILE_GROUPS:
          switch (eval_frc) {
          case EvaluateForce::YES:
            switch (eval_nrg) {
            case EvaluateEnergy::YES:
              switch (acc_meth) {
              case AccumulationMethod::SPLIT:
                if (cudaFuncGetAttributes(&attr, ktgfsVacuumForceEnergy) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfsVacuumForceEnergy.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              case AccumulationMethod::WHOLE:
                if (cudaFuncGetAttributes(&attr, ktgfVacuumForceEnergy) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfVacuumForceEnergy.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              }
              break;
            case EvaluateEnergy::NO:
              switch (acc_meth) {
              case AccumulationMethod::SPLIT:
                if (cudaFuncGetAttributes(&attr, ktgfsVacuumForce) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfsVacuumForce.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              case AccumulationMethod::WHOLE:
                if (cudaFuncGetAttributes(&attr, ktgfVacuumForce) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfVacuumForce.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              }
              break;
            }
            break;
          case EvaluateForce::NO:
            if (cudaFuncGetAttributes(&attr, ktgfVacuumEnergy) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgfVacuumEnergy.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          }
          break;
        case NbwuKind::SUPERTILES:
          break;
        case NbwuKind::HONEYCOMB:
          break;
        }
        break;
      }
      break;
    case ImplicitSolventModel::HCT_GB:
    case ImplicitSolventModel::OBC_GB:
    case ImplicitSolventModel::OBC_GB_II:
      switch (prec) {
      case PrecisionModel::DOUBLE:
        switch (eval_frc) {
        case EvaluateForce::YES:
          switch (eval_nrg) {
          case EvaluateEnergy::YES:
            if (cudaFuncGetAttributes(&attr, ktgdsGBForceEnergy) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgdsGBForceEnergy.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          case EvaluateEnergy::NO:
            if (cudaFuncGetAttributes(&attr, ktgdsGBForce) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgdsGBForce.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          }
          break;
        case EvaluateForce::NO:
          if (cudaFuncGetAttributes(&attr, ktgdGBEnergy) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel ktgdGBEnergy.",
                  "queryNonbondedKernelRequirements");
          }
          break;
        }
        break;
      case PrecisionModel::SINGLE:
        switch (kind) {
          case NbwuKind::TILE_GROUPS:
          switch (eval_frc) {
          case EvaluateForce::YES:
            switch (eval_nrg) {
            case EvaluateEnergy::YES:
              switch (acc_meth) {
              case AccumulationMethod::SPLIT:
                if (cudaFuncGetAttributes(&attr, ktgfsGBForceEnergy) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfsGBForceEnergy.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              case AccumulationMethod::WHOLE:
                if (cudaFuncGetAttributes(&attr, ktgfGBForceEnergy) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfGBForceEnergy.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              }
              break;
            case EvaluateEnergy::NO:
              switch (acc_meth) {
              case AccumulationMethod::SPLIT:
                if (cudaFuncGetAttributes(&attr, ktgfsGBForce) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfsGBForce.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              case AccumulationMethod::WHOLE:
                if (cudaFuncGetAttributes(&attr, ktgfGBForce) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfGBForce.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              }
              break;
            }
            break;
          case EvaluateForce::NO:
            if (cudaFuncGetAttributes(&attr, ktgfGBEnergy) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgfGBEnergy.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          }
          break;
        case NbwuKind::SUPERTILES:
          break;
        case NbwuKind::HONEYCOMB:
          break;
        }
        break;
      }
      break;
    case ImplicitSolventModel::NECK_GB:
    case ImplicitSolventModel::NECK_GB_II:
      switch (prec) {
      case PrecisionModel::DOUBLE:
        switch (eval_frc) {
        case EvaluateForce::YES:
          switch (eval_nrg) {
          case EvaluateEnergy::YES:
            if (cudaFuncGetAttributes(&attr, ktgdsGBNeckForceEnergy) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgdsGBNeckForceEnergy.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          case EvaluateEnergy::NO:
            if (cudaFuncGetAttributes(&attr, ktgdsGBNeckForce) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgdsGBNeckForce.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          }
          break;
        case EvaluateForce::NO:
          if (cudaFuncGetAttributes(&attr, ktgdGBNeckEnergy) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel ktgdGBNeckEnergy.",
                  "queryNonbondedKernelRequirements");
          }
          break;
        }
        break;
      case PrecisionModel::SINGLE:
        switch (kind) {
          case NbwuKind::TILE_GROUPS:
          switch (eval_frc) {
          case EvaluateForce::YES:
            switch (eval_nrg) {
            case EvaluateEnergy::YES:
              switch (acc_meth) {
              case AccumulationMethod::SPLIT:
                if (cudaFuncGetAttributes(&attr, ktgfsGBNeckForceEnergy) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfsGBNeckForceEnergy.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              case AccumulationMethod::WHOLE:
                if (cudaFuncGetAttributes(&attr, ktgfGBNeckForceEnergy) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfGBNeckForceEnergy.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              }
              break;
            case EvaluateEnergy::NO:
              switch (acc_meth) {
              case AccumulationMethod::SPLIT:
                if (cudaFuncGetAttributes(&attr, ktgfsGBNeckForce) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfsGBNeckForce.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              case AccumulationMethod::WHOLE:
                if (cudaFuncGetAttributes(&attr, ktgfGBNeckForce) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfGBNeckForce.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              }
              break;
            }
            break;
          case EvaluateForce::NO:
            if (cudaFuncGetAttributes(&attr, ktgfGBNeckEnergy) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgfGBNeckEnergy.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          }
          break;
        case NbwuKind::SUPERTILES:
          break;
        case NbwuKind::HONEYCOMB:
          break;
        }
        break;
      }
      break;
    }
    break;
  case ClashResponse::FORGIVE:
    switch (igb) {
    case ImplicitSolventModel::NONE:
      switch (prec) {
      case PrecisionModel::DOUBLE:
        switch (eval_frc) {
        case EvaluateForce::YES:
          switch (eval_nrg) {
          case EvaluateEnergy::YES:
            if (cudaFuncGetAttributes(&attr, ktgdsVacuumForceEnergyNonClash) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgdsVacuumForceEnergyNonClash.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          case EvaluateEnergy::NO:
            if (cudaFuncGetAttributes(&attr, ktgdsVacuumForceNonClash) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgdsVacuumForceNonClash.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          }
          break;
        case EvaluateForce::NO:
          if (cudaFuncGetAttributes(&attr, ktgdVacuumEnergyNonClash) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel ktgdVacuumEnergyNonClash.",
                  "queryNonbondedKernelRequirements");
          }
          break;
        }
        break;
      case PrecisionModel::SINGLE:
        switch (kind) {
          case NbwuKind::TILE_GROUPS:
          switch (eval_frc) {
          case EvaluateForce::YES:
            switch (eval_nrg) {
            case EvaluateEnergy::YES:
              switch (acc_meth) {
              case AccumulationMethod::SPLIT:
                if (cudaFuncGetAttributes(&attr, ktgfsVacuumForceEnergyNonClash) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfsVacuumForceEnergyNonClash.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              case AccumulationMethod::WHOLE:
                if (cudaFuncGetAttributes(&attr, ktgfVacuumForceEnergyNonClash) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfVacuumForceEnergyNonClash.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              }
              break;
            case EvaluateEnergy::NO:
              switch (acc_meth) {
              case AccumulationMethod::SPLIT:
                if (cudaFuncGetAttributes(&attr, ktgfsVacuumForceNonClash) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfsVacuumForceNonClash.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              case AccumulationMethod::WHOLE:
                if (cudaFuncGetAttributes(&attr, ktgfVacuumForceNonClash) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfVacuumForceNonClash.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              }
              break;
            }
            break;
          case EvaluateForce::NO:
            if (cudaFuncGetAttributes(&attr, ktgfVacuumEnergyNonClash) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgfVacuumEnergyNonClash.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          }
          break;
        case NbwuKind::SUPERTILES:
          break;
        case NbwuKind::HONEYCOMB:
          break;
        }
        break;
      }
      break;
    case ImplicitSolventModel::HCT_GB:
    case ImplicitSolventModel::OBC_GB:
    case ImplicitSolventModel::OBC_GB_II:
      switch (prec) {
      case PrecisionModel::DOUBLE:
        switch (eval_frc) {
        case EvaluateForce::YES:
          switch (eval_nrg) {
          case EvaluateEnergy::YES:
            if (cudaFuncGetAttributes(&attr, ktgdsGBForceEnergyNonClash) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgdsGBForceEnergyNonClash.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          case EvaluateEnergy::NO:
            if (cudaFuncGetAttributes(&attr, ktgdsGBForceNonClash) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgdsGBForceNonClash.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          }
          break;
        case EvaluateForce::NO:
          if (cudaFuncGetAttributes(&attr, ktgdGBEnergyNonClash) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel ktgdGBEnergyNonClash.",
                  "queryNonbondedKernelRequirements");
          }
          break;
        }
        break;
      case PrecisionModel::SINGLE:
        switch (kind) {
          case NbwuKind::TILE_GROUPS:
          switch (eval_frc) {
          case EvaluateForce::YES:
            switch (eval_nrg) {
            case EvaluateEnergy::YES:
              switch (acc_meth) {
              case AccumulationMethod::SPLIT:
                if (cudaFuncGetAttributes(&attr, ktgfsGBForceEnergyNonClash) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfsGBForceEnergyNonClash.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              case AccumulationMethod::WHOLE:
                if (cudaFuncGetAttributes(&attr, ktgfGBForceEnergyNonClash) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfGBForceEnergyNonClash.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              }
              break;
            case EvaluateEnergy::NO:
              switch (acc_meth) {
              case AccumulationMethod::SPLIT:
                if (cudaFuncGetAttributes(&attr, ktgfsGBForceNonClash) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfsGBForceNonClash.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              case AccumulationMethod::WHOLE:
                if (cudaFuncGetAttributes(&attr, ktgfGBForceNonClash) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfGBForceNonClash.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              }
              break;
            }
            break;
          case EvaluateForce::NO:
            if (cudaFuncGetAttributes(&attr, ktgfGBEnergyNonClash) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgfGBEnergyNonClash.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          }
          break;
        case NbwuKind::SUPERTILES:
          break;
        case NbwuKind::HONEYCOMB:
          break;
        }
        break;
      }
      break;
    case ImplicitSolventModel::NECK_GB:
    case ImplicitSolventModel::NECK_GB_II:
      switch (prec) {
      case PrecisionModel::DOUBLE:
        switch (eval_frc) {
        case EvaluateForce::YES:
          switch (eval_nrg) {
          case EvaluateEnergy::YES:
            if (cudaFuncGetAttributes(&attr, ktgdsGBNeckForceEnergyNonClash) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgdsGBNeckForceEnergyNonClash.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          case EvaluateEnergy::NO:
            if (cudaFuncGetAttributes(&attr, ktgdsGBNeckForceNonClash) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgdsGBNeckForceNonClash.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          }
          break;
        case EvaluateForce::NO:
          if (cudaFuncGetAttributes(&attr, ktgdGBNeckEnergyNonClash) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel ktgdGBNeckEnergyNonClash.",
                  "queryNonbondedKernelRequirements");
          }
          break;
        }
        break;
      case PrecisionModel::SINGLE:
        switch (kind) {
          case NbwuKind::TILE_GROUPS:
          switch (eval_frc) {
          case EvaluateForce::YES:
            switch (eval_nrg) {
            case EvaluateEnergy::YES:
              switch (acc_meth) {
              case AccumulationMethod::SPLIT:
                if (cudaFuncGetAttributes(&attr, ktgfsGBNeckForceEnergyNonClash) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfsGBNeckForceEnergyNonClash.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              case AccumulationMethod::WHOLE:
                if (cudaFuncGetAttributes(&attr, ktgfGBNeckForceEnergyNonClash) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfGBNeckForceEnergyNonClash.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              }
              break;
            case EvaluateEnergy::NO:
              switch (acc_meth) {
              case AccumulationMethod::SPLIT:
                if (cudaFuncGetAttributes(&attr, ktgfsGBNeckForceNonClash) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfsGBNeckForceNonClash.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              case AccumulationMethod::WHOLE:
                if (cudaFuncGetAttributes(&attr, ktgfGBNeckForceNonClash) != cudaSuccess) {
                  rtErr("Error obtaining attributes for kernel ktgfGBNeckForceNonClash.",
                        "queryNonbondedKernelRequirements");
                }
                break;
              }
              break;
            }
            break;
          case EvaluateForce::NO:
            if (cudaFuncGetAttributes(&attr, ktgfGBNeckEnergyNonClash) != cudaSuccess) {
              rtErr("Error obtaining attributes for kernel ktgfGBNeckEnergyNonClash.",
                    "queryNonbondedKernelRequirements");
            }
            break;
          }
          break;
        case NbwuKind::SUPERTILES:
          break;
        case NbwuKind::HONEYCOMB:
          break;
        }
        break;
      }
      break;
    }
    break;
  }
  return attr;
}

//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes
queryBornRadiiKernelRequirements(const PrecisionModel prec, const NbwuKind kind,
                                 const AccumulationMethod acc_meth,
                                 const ImplicitSolventModel igb) {
  cudaFuncAttributes attr;
  switch (igb) {
  case ImplicitSolventModel::NONE:
    break;
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
    switch (prec) {
    case PrecisionModel::DOUBLE:
      switch (kind) {
      case NbwuKind::TILE_GROUPS:
        if (cudaFuncGetAttributes(&attr, ktgdsCalculateGBRadii) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel ktgdCalculateGBRadii.",
                "queryNonbondedKernelRequirements");
        }
        break;
      case NbwuKind::SUPERTILES:
      case NbwuKind::HONEYCOMB:
        break;
      }
      break;
    case PrecisionModel::SINGLE:
      switch (kind) {
      case NbwuKind::TILE_GROUPS:
        switch (acc_meth) {
        case AccumulationMethod::SPLIT:
          if (cudaFuncGetAttributes(&attr, ktgfsCalculateGBRadii) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel ktgfsCalculateGBRadii.",
                  "queryNonbondedKernelRequirements");
          }
          break;
        case AccumulationMethod::WHOLE:
          if (cudaFuncGetAttributes(&attr, ktgfCalculateGBRadii) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel ktgfCalculateGBRadii.",
                  "queryNonbondedKernelRequirements");
          }
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case NbwuKind::SUPERTILES:
      case NbwuKind::HONEYCOMB:
        break;
      }
      break;
    }
    break;
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    switch (prec) {
    case PrecisionModel::DOUBLE:
      switch (kind) {
      case NbwuKind::TILE_GROUPS:
        if (cudaFuncGetAttributes(&attr, ktgdsCalculateGBNeckRadii) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel ktgdCalculateGBNeckRadii.",
                "queryNonbondedKernelRequirements");
        }
        break;
      case NbwuKind::SUPERTILES:
      case NbwuKind::HONEYCOMB:
        break;
      }
      break;
    case PrecisionModel::SINGLE:
      switch (kind) {
      case NbwuKind::TILE_GROUPS:
        switch (acc_meth) {
        case AccumulationMethod::SPLIT:
          if (cudaFuncGetAttributes(&attr, ktgfsCalculateGBNeckRadii) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel ktgfsCalculateGBNeckRadii.",
                  "queryNonbondedKernelRequirements");
          }
          break;
        case AccumulationMethod::WHOLE:
          if (cudaFuncGetAttributes(&attr, ktgfCalculateGBNeckRadii) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel ktgfCalculateGBNeckRadii.",
                  "queryNonbondedKernelRequirements");
          }
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case NbwuKind::SUPERTILES:
      case NbwuKind::HONEYCOMB:
        break;
      }
      break;
    }
    break;
  }
  return attr;
}

//-------------------------------------------------------------------------------------------------
extern cudaFuncAttributes
queryBornDerivativeKernelRequirements(const PrecisionModel prec, const NbwuKind kind,
                                      const AccumulationMethod acc_meth,
                                      const ImplicitSolventModel igb) {
  cudaFuncAttributes attr;
  switch (igb) {
  case ImplicitSolventModel::NONE:
    break;
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
    switch (prec) {
    case PrecisionModel::DOUBLE:
      switch (kind) {
      case NbwuKind::TILE_GROUPS:
        if (cudaFuncGetAttributes(&attr, ktgdsCalculateGBDerivatives) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel ktgdCalculateGBDerivatives.",
                "queryNonbondedKernelRequirements");
        }
        break;
      case NbwuKind::SUPERTILES:
      case NbwuKind::HONEYCOMB:
        break;
      }
      break;
    case PrecisionModel::SINGLE:
      switch (kind) {
      case NbwuKind::TILE_GROUPS:
        switch (acc_meth) {
        case AccumulationMethod::SPLIT:
          if (cudaFuncGetAttributes(&attr, ktgfsCalculateGBDerivatives) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel ktgfsCalculateGBDerivatives.",
                  "queryNonbondedKernelRequirements");
          }
          break;
        case AccumulationMethod::WHOLE:
          if (cudaFuncGetAttributes(&attr, ktgfCalculateGBDerivatives) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel ktgfCalculateGBDerivatives.",
                  "queryNonbondedKernelRequirements");
          }
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case NbwuKind::SUPERTILES:
      case NbwuKind::HONEYCOMB:
        break;
      }
      break;
    }
    break;
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    switch (prec) {
    case PrecisionModel::DOUBLE:
      switch (kind) {
      case NbwuKind::TILE_GROUPS:
        if (cudaFuncGetAttributes(&attr, ktgdsCalculateGBNeckDerivatives) != cudaSuccess) {
          rtErr("Error obtaining attributes for kernel ktgdCalculateGBNeckDerivatives.",
                "queryNonbondedKernelRequirements");
        }
        break;
      case NbwuKind::SUPERTILES:
      case NbwuKind::HONEYCOMB:
        break;
      }
      break;
    case PrecisionModel::SINGLE:
      switch (kind) {
      case NbwuKind::TILE_GROUPS:
        switch (acc_meth) {
        case AccumulationMethod::SPLIT:
          if (cudaFuncGetAttributes(&attr, ktgfsCalculateGBNeckDerivatives) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel ktgfsCalculateGBNeckDerivatives.",
                  "queryNonbondedKernelRequirements");
          }
          break;
        case AccumulationMethod::WHOLE:
          if (cudaFuncGetAttributes(&attr, ktgfCalculateGBNeckDerivatives) != cudaSuccess) {
            rtErr("Error obtaining attributes for kernel ktgfCalculateGBNeckDerivatives.",
                  "queryNonbondedKernelRequirements");
          }
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case NbwuKind::SUPERTILES:
      case NbwuKind::HONEYCOMB:
        break;
      }
      break;
    }
    break;
  }    
  return attr;
}

//-------------------------------------------------------------------------------------------------
extern void launchNonbonded(const NbwuKind kind, const SyNonbondedKit<double, double2> &poly_nbk,
                            const SeMaskSynthesisReader &poly_ser, MMControlKit<double> *ctrl,
                            PsSynthesisWriter *poly_psw, ThermostatWriter<double> *tstw,
                            ScoreCardWriter *scw, CacheResourceKit<double> *gmem_r,
                            ISWorkspaceKit<double> *iswk, const EvaluateForce eval_force,
                            const EvaluateEnergy eval_energy, const int2 bt, const int2 gbr_bt,
                            const int2 gbd_bt, const double clash_minimum_distance,
                            const double clash_ratio) {

  // As with the valence kernels, detect the intention to intercept and dampen clashes by nonzero
  // values in the parameters, not an explicit input based on the ClashResponse enumerator.
  if (clash_minimum_distance >= 1.0e-6 || clash_ratio >= 1.0e-6) {
    switch (kind) {
    case NbwuKind::TILE_GROUPS:
      switch (poly_nbk.igb) {
      case ImplicitSolventModel::NONE:
        switch (eval_force) {
        case EvaluateForce::YES:
          switch (eval_energy) {
          case EvaluateEnergy::YES:
            ktgdsVacuumForceEnergyNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                           clash_minimum_distance, clash_ratio,
                                                           *scw, *tstw, *gmem_r);
            break;
          case EvaluateEnergy::NO:
            ktgdsVacuumForceNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                     clash_minimum_distance, clash_ratio, *tstw,
                                                     *gmem_r);
            break;
          }
          break;
        case EvaluateForce::NO:
          ktgdVacuumEnergyNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                   clash_minimum_distance, clash_ratio, *scw,
                                                   *tstw, *gmem_r);
          break;
        }
        break;
      case ImplicitSolventModel::HCT_GB:
      case ImplicitSolventModel::OBC_GB:
      case ImplicitSolventModel::OBC_GB_II:
        ktgdsCalculateGBRadii<<<gbr_bt.x, gbr_bt.y>>>(poly_nbk, *ctrl, *poly_psw, *tstw, *iswk,
                                                      *gmem_r);
        switch (eval_force) {
        case EvaluateForce::YES:
          switch (eval_energy) {
          case EvaluateEnergy::YES:
            ktgdsGBForceEnergyNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                       clash_minimum_distance, clash_ratio, *scw,
                                                       *iswk, *gmem_r);
            break;
          case EvaluateEnergy::NO:
            ktgdsGBForceNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                 clash_minimum_distance, clash_ratio, *iswk,
                                                 *gmem_r);
            break;
          }
          ktgdsCalculateGBDerivatives<<<gbd_bt.x, gbd_bt.y>>>(poly_nbk, *ctrl, *poly_psw, *iswk,
                                                              *gmem_r);
          break;
        case EvaluateForce::NO:
          ktgdGBEnergyNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                               clash_minimum_distance, clash_ratio, *scw, *iswk,
                                               *gmem_r);
          break;
        }
        break;
      case ImplicitSolventModel::NECK_GB:
      case ImplicitSolventModel::NECK_GB_II:
        ktgdsCalculateGBNeckRadii<<<gbr_bt.x, gbr_bt.y>>>(poly_nbk, *ctrl, *poly_psw, *tstw, *iswk,
                                                          *gmem_r);
        switch (eval_force) {
        case EvaluateForce::YES:
          switch (eval_energy) {
          case EvaluateEnergy::YES:
            ktgdsGBNeckForceEnergyNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                           clash_minimum_distance, clash_ratio,
                                                           *scw, *iswk, *gmem_r);
            break;
          case EvaluateEnergy::NO:
            ktgdsGBNeckForceNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                     clash_minimum_distance, clash_ratio, *iswk,
                                                     *gmem_r);
            break;
          }
          ktgdsCalculateGBNeckDerivatives<<<gbd_bt.x, gbd_bt.y>>>(poly_nbk, *ctrl, *poly_psw,
                                                                  *iswk, *gmem_r);
          break;
        case EvaluateForce::NO:
          ktgdGBNeckEnergyNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                   clash_minimum_distance, clash_ratio, *scw,
                                                   *iswk, *gmem_r);
          break;
        }
        break;
      }
      break;
    case NbwuKind::SUPERTILES:
    case NbwuKind::HONEYCOMB:
      break;
    }
  }
  else {
    switch (kind) {
    case NbwuKind::TILE_GROUPS:
      switch (poly_nbk.igb) {
      case ImplicitSolventModel::NONE:
        switch (eval_force) {
        case EvaluateForce::YES:
          switch (eval_energy) {
          case EvaluateEnergy::YES:
            ktgdsVacuumForceEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                                   *tstw, *gmem_r);
            break;
          case EvaluateEnergy::NO:
            ktgdsVacuumForce<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *tstw, *gmem_r);
            break;
          }
          break;
        case EvaluateForce::NO:
          ktgdVacuumEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw, *tstw,
                                           *gmem_r);
          break;
        }
        break;
      case ImplicitSolventModel::HCT_GB:
      case ImplicitSolventModel::OBC_GB:
      case ImplicitSolventModel::OBC_GB_II:
        ktgdsCalculateGBRadii<<<gbr_bt.x, gbr_bt.y>>>(poly_nbk, *ctrl, *poly_psw, *tstw, *iswk,
                                                      *gmem_r);
        switch (eval_force) {
        case EvaluateForce::YES:
          switch (eval_energy) {
          case EvaluateEnergy::YES:
            ktgdsGBForceEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                               *iswk, *gmem_r);
            break;
          case EvaluateEnergy::NO:
            ktgdsGBForce<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *iswk,
                                         *gmem_r);
            break;
          }
          ktgdsCalculateGBDerivatives<<<gbd_bt.x, gbd_bt.y>>>(poly_nbk, *ctrl, *poly_psw, *iswk,
                                                              *gmem_r);
          break;
        case EvaluateForce::NO:
          ktgdGBEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw, *iswk, *gmem_r);
          break;
        }
        break;
      case ImplicitSolventModel::NECK_GB:
      case ImplicitSolventModel::NECK_GB_II:
        ktgdsCalculateGBNeckRadii<<<gbr_bt.x, gbr_bt.y>>>(poly_nbk, *ctrl, *poly_psw, *tstw, *iswk,
                                                          *gmem_r);
        switch (eval_force) {
        case EvaluateForce::YES:
          switch (eval_energy) {
          case EvaluateEnergy::YES:
            ktgdsGBNeckForceEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                                   *iswk, *gmem_r);
            break;
          case EvaluateEnergy::NO:
            ktgdsGBNeckForce<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *iswk, *gmem_r);
            break;
          }
          ktgdsCalculateGBNeckDerivatives<<<gbd_bt.x, gbd_bt.y>>>(poly_nbk, *ctrl, *poly_psw,
                                                                  *iswk, *gmem_r);
          break;
        case EvaluateForce::NO:
          ktgdGBNeckEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw, *iswk,
                                           *gmem_r);
          break;
        }
        break;
      }
      break;
    case NbwuKind::SUPERTILES:
    case NbwuKind::HONEYCOMB:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchNonbonded(const NbwuKind kind, const SyNonbondedKit<float, float2> &poly_nbk,
                            const SeMaskSynthesisReader &poly_ser, MMControlKit<float> *ctrl,
                            PsSynthesisWriter *poly_psw, ThermostatWriter<float> *tstw,
                            ScoreCardWriter *scw, CacheResourceKit<float> *gmem_r,
                            ISWorkspaceKit<float> *iswk, const EvaluateForce eval_force,
                            const EvaluateEnergy eval_energy, const AccumulationMethod force_sum,
                            const int2 bt, const int2 gbr_bt, const int2 gbd_bt,
                            const float clash_minimum_distance, const float clash_ratio) {
  const AccumulationMethod actual_force_sum = (force_sum == AccumulationMethod::AUTOMATIC) ?
                                              chooseAccumulationMethod(poly_psw->frc_bits) :
                                              force_sum;
  // As with the valence kernels, detect the intention to intercept and dampen clashes by nonzero
  // values in the parameters, not an explicit input based on the ClashResponse enumerator.
  if (clash_minimum_distance > 0.000001f || clash_ratio > 0.000001f) {
    switch (kind) {
    case NbwuKind::TILE_GROUPS:
      switch (poly_nbk.igb) {
      case ImplicitSolventModel::NONE:
        switch (eval_force) {
        case EvaluateForce::YES:
          switch (actual_force_sum) {
          case AccumulationMethod::SPLIT:
            switch (eval_energy) {
            case EvaluateEnergy::YES:
              ktgfsVacuumForceEnergyNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                             clash_minimum_distance, clash_ratio,
                                                             *scw, *tstw, *gmem_r);
              break;
            case EvaluateEnergy::NO:
              ktgfsVacuumForceNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                       clash_minimum_distance, clash_ratio, *tstw,
                                                       *gmem_r);
              break;
            }
            break;
          case AccumulationMethod::WHOLE:
            switch (eval_energy) {
            case EvaluateEnergy::YES:
              ktgfVacuumForceEnergyNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                            clash_minimum_distance, clash_ratio,
                                                            *scw, *tstw, *gmem_r);
              break;
            case EvaluateEnergy::NO:
              ktgfVacuumForceNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                      clash_minimum_distance, clash_ratio, *tstw,
                                                      *gmem_r);
              break;
            }
            break;
          case AccumulationMethod::AUTOMATIC:

            // This case was converted into SPLIT or WHOLE by evaluating actual_force_sum
            break;
          }
          break;
        case EvaluateForce::NO:
          ktgfVacuumEnergyNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                   clash_minimum_distance, clash_ratio, *scw,
                                                   *tstw, *gmem_r);
          break;
        }
        break;
      case ImplicitSolventModel::HCT_GB:
      case ImplicitSolventModel::OBC_GB:
      case ImplicitSolventModel::OBC_GB_II:
        switch (actual_force_sum) {
        case AccumulationMethod::SPLIT:
          ktgfsCalculateGBRadii<<<gbr_bt.x, gbr_bt.y>>>(poly_nbk, *ctrl, *poly_psw, *tstw, *iswk,
                                                        *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          ktgfCalculateGBRadii<<<gbr_bt.x, gbr_bt.y>>>(poly_nbk, *ctrl, *poly_psw, *tstw, *iswk,
                                                       *gmem_r);
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        switch (eval_force) {
        case EvaluateForce::YES:
          switch (actual_force_sum) {
          case AccumulationMethod::SPLIT:
            switch (eval_energy) {
            case EvaluateEnergy::YES:
              ktgfsGBForceEnergyNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                         clash_minimum_distance, clash_ratio, *scw,
                                                         *iswk, *gmem_r);
              break;
            case EvaluateEnergy::NO:
              ktgfsGBForceNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                   clash_minimum_distance, clash_ratio, *iswk,
                                                   *gmem_r);
              break;
            }
            ktgfsCalculateGBDerivatives<<<gbd_bt.x, gbd_bt.y>>>(poly_nbk, *ctrl, *poly_psw, *iswk,
                                                                *gmem_r);
            break;
          case AccumulationMethod::WHOLE:
            switch (eval_energy) {
            case EvaluateEnergy::YES:
              ktgfGBForceEnergyNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                        clash_minimum_distance, clash_ratio, *scw,
                                                        *iswk, *gmem_r);
              break;
            case EvaluateEnergy::NO:
              ktgfGBForceNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                  clash_minimum_distance, clash_ratio, *iswk,
                                                  *gmem_r);
              break;
            }
            ktgfCalculateGBDerivatives<<<gbd_bt.x, gbd_bt.y>>>(poly_nbk, *ctrl, *poly_psw, *iswk,
                                                               *gmem_r);
            break;
          case AccumulationMethod::AUTOMATIC:
            break;
          }
          break;
        case EvaluateForce::NO:
          ktgfGBEnergyNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                               clash_minimum_distance, clash_ratio, *scw, *iswk,
                                               *gmem_r);
          break;
        }
        break;
      case ImplicitSolventModel::NECK_GB:
      case ImplicitSolventModel::NECK_GB_II:
        switch (actual_force_sum) {
        case AccumulationMethod::SPLIT:
          ktgfsCalculateGBNeckRadii<<<gbr_bt.x, gbr_bt.y>>>(poly_nbk, *ctrl, *poly_psw, *tstw,
                                                            *iswk, *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          ktgfCalculateGBNeckRadii<<<gbr_bt.x, gbr_bt.y>>>(poly_nbk, *ctrl, *poly_psw, *tstw,
                                                           *iswk, *gmem_r);
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        switch (eval_force) {
        case EvaluateForce::YES:
          switch (actual_force_sum) {
          case AccumulationMethod::SPLIT:
            switch (eval_energy) {
            case EvaluateEnergy::YES:
              ktgfsGBNeckForceEnergyNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                             clash_minimum_distance, clash_ratio,
                                                             *scw, *iswk, *gmem_r);
              break;
            case EvaluateEnergy::NO:
              ktgfsGBNeckForceNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                       clash_minimum_distance, clash_ratio, *iswk,
                                                       *gmem_r);
              break;
            }
            ktgfsCalculateGBNeckDerivatives<<<gbd_bt.x, gbd_bt.y>>>(poly_nbk, *ctrl, *poly_psw,
                                                                    *iswk, *gmem_r);
            break;
          case AccumulationMethod::WHOLE:
            switch (eval_energy) {
            case EvaluateEnergy::YES:
              ktgfGBNeckForceEnergyNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                            clash_minimum_distance, clash_ratio,
                                                            *scw, *iswk, *gmem_r);
              break;
            case EvaluateEnergy::NO:
              ktgfGBNeckForceNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                      clash_minimum_distance, clash_ratio, *iswk,
                                                      *gmem_r);
              break;
            }
            ktgfCalculateGBNeckDerivatives<<<gbd_bt.x, gbd_bt.y>>>(poly_nbk, *ctrl, *poly_psw,
                                                                   *iswk, *gmem_r);
            break;
          case AccumulationMethod::AUTOMATIC:
            break;
          }
          break;
        case EvaluateForce::NO:
          ktgfGBNeckEnergyNonClash<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw,
                                                   clash_minimum_distance, clash_ratio, *scw,
                                                   *iswk, *gmem_r);
          break;
        }
        break;
      }
      break;
    case NbwuKind::SUPERTILES:
    case NbwuKind::HONEYCOMB:
      break;
    }
  }
  else {
    switch (kind) {
    case NbwuKind::TILE_GROUPS:
      switch (poly_nbk.igb) {
      case ImplicitSolventModel::NONE:
        switch (eval_force) {
        case EvaluateForce::YES:
          switch (actual_force_sum) {
          case AccumulationMethod::SPLIT:
            switch (eval_energy) {
            case EvaluateEnergy::YES:
              ktgfsVacuumForceEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                                     *tstw, *gmem_r);
              break;
            case EvaluateEnergy::NO:
              ktgfsVacuumForce<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *tstw,
                                               *gmem_r);
              break;
            }
            break;
          case AccumulationMethod::WHOLE:
            switch (eval_energy) {
            case EvaluateEnergy::YES:
              ktgfVacuumForceEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                                    *tstw, *gmem_r);
              break;
            case EvaluateEnergy::NO:
              ktgfVacuumForce<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *tstw,
                                              *gmem_r);
              break;
            }
            break;
          case AccumulationMethod::AUTOMATIC:

            // This case was converted into SPLIT or WHOLE by evaluating actual_force_sum
            break;
          }
          break;
        case EvaluateForce::NO:
          ktgfVacuumEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw, *tstw,
                                           *gmem_r);
          break;
        }
        break;
      case ImplicitSolventModel::HCT_GB:
      case ImplicitSolventModel::OBC_GB:
      case ImplicitSolventModel::OBC_GB_II:
        switch (actual_force_sum) {
        case AccumulationMethod::SPLIT:
          ktgfsCalculateGBRadii<<<gbr_bt.x, gbr_bt.y>>>(poly_nbk, *ctrl, *poly_psw, *tstw, *iswk,
                                                        *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          ktgfCalculateGBRadii<<<gbr_bt.x, gbr_bt.y>>>(poly_nbk, *ctrl, *poly_psw, *tstw, *iswk,
                                                       *gmem_r);
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        switch (eval_force) {
        case EvaluateForce::YES:
          switch (actual_force_sum) {
          case AccumulationMethod::SPLIT:
            switch (eval_energy) {
            case EvaluateEnergy::YES:
              ktgfsGBForceEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw, *iswk,
                                                 *gmem_r);
              break;
            case EvaluateEnergy::NO:
              ktgfsGBForce<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *iswk, *gmem_r);
              break;
            }
            ktgfsCalculateGBDerivatives<<<gbd_bt.x, gbd_bt.y>>>(poly_nbk, *ctrl, *poly_psw, *iswk,
                                                                *gmem_r);
            break;
          case AccumulationMethod::WHOLE:
            switch (eval_energy) {
            case EvaluateEnergy::YES:
              ktgfGBForceEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw, *iswk,
                                                *gmem_r);
              break;
            case EvaluateEnergy::NO:
              ktgfGBForce<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *iswk, *gmem_r);
              break;
            }
            ktgfCalculateGBDerivatives<<<gbd_bt.x, gbd_bt.y>>>(poly_nbk, *ctrl, *poly_psw, *iswk,
                                                               *gmem_r);
            break;
          case AccumulationMethod::AUTOMATIC:
            break;
          }
          break;
        case EvaluateForce::NO:
          ktgfGBEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw, *iswk, *gmem_r);
          break;
        }
        break;
      case ImplicitSolventModel::NECK_GB:
      case ImplicitSolventModel::NECK_GB_II:
        switch (actual_force_sum) {
        case AccumulationMethod::SPLIT:
          ktgfsCalculateGBNeckRadii<<<gbr_bt.x, gbr_bt.y>>>(poly_nbk, *ctrl, *poly_psw, *tstw,
                                                            *iswk, *gmem_r);
          break;
        case AccumulationMethod::WHOLE:
          ktgfCalculateGBNeckRadii<<<gbr_bt.x, gbr_bt.y>>>(poly_nbk, *ctrl, *poly_psw, *tstw,
                                                           *iswk, *gmem_r);
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        switch (eval_force) {
        case EvaluateForce::YES:
          switch (actual_force_sum) {
          case AccumulationMethod::SPLIT:
            switch (eval_energy) {
            case EvaluateEnergy::YES:
              ktgfsGBNeckForceEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                                     *iswk, *gmem_r);
              break;
            case EvaluateEnergy::NO:
              ktgfsGBNeckForce<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *iswk,
                                               *gmem_r);
              break;
            }
            ktgfsCalculateGBNeckDerivatives<<<gbd_bt.x, gbd_bt.y>>>(poly_nbk, *ctrl, *poly_psw,
                                                                    *iswk, *gmem_r);
            break;
          case AccumulationMethod::WHOLE:
            switch (eval_energy) {
            case EvaluateEnergy::YES:
              ktgfGBNeckForceEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw,
                                                    *iswk, *gmem_r);
              break;
            case EvaluateEnergy::NO:
              ktgfGBNeckForce<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *iswk,
                                              *gmem_r);
              break;
            }
            ktgfCalculateGBNeckDerivatives<<<gbd_bt.x, gbd_bt.y>>>(poly_nbk, *ctrl, *poly_psw,
                                                                   *iswk, *gmem_r);
            break;
          case AccumulationMethod::AUTOMATIC:
            break;
          }
          break;
        case EvaluateForce::NO:
          ktgfGBNeckEnergy<<<bt.x, bt.y>>>(poly_nbk, poly_ser, *ctrl, *poly_psw, *scw, *iswk,
                                           *gmem_r);
          break;
        }
        break;
      }
      break;
    case NbwuKind::SUPERTILES:
    case NbwuKind::HONEYCOMB:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchNonbonded(const PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                            const StaticExclusionMaskSynthesis &poly_se,
                            MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                            Thermostat *heat_bath, ScoreCard *sc, CacheResource *tb_space,
                            ImplicitSolventWorkspace *ism_space, const EvaluateForce eval_force,
                            const EvaluateEnergy eval_energy, const AccumulationMethod force_sum,
                            const CoreKlManager &launcher, const double clash_minimum_distance,
                            const double clash_ratio) {
  const HybridTargetLevel tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(tier);
  ScoreCardWriter scw = sc->data(tier);
  const SeMaskSynthesisReader poly_ser = poly_se.data();
  const NbwuKind nb_work_type = poly_ag.getNonbondedWorkType();
  const ImplicitSolventModel ism_type = poly_ag.getImplicitSolventModel();
  const int2 bt = launcher.getNonbondedKernelDims(prec, nb_work_type, eval_force, eval_energy,
                                                  force_sum, ism_type, ClashResponse::NONE);
  const int2 gbr_bt = launcher.getBornRadiiKernelDims(prec, nb_work_type, force_sum, ism_type);
  const int2 gbd_bt = launcher.getBornDerivativeKernelDims(prec, nb_work_type, force_sum,
                                                           ism_type);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyNonbondedKit<double,
                           double2> poly_nbk = poly_ag.getDoublePrecisionNonbondedKit(tier);
      MMControlKit<double> ctrl = mmctrl->dpData(tier);
      CacheResourceKit<double> gmem_r = tb_space->dpData(tier);
      ThermostatWriter<double> tstw = heat_bath->dpData(tier);
      ISWorkspaceKit<double> iswk = ism_space->dpData(poly_ps->getCyclePosition(), tier);
      launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl, &poly_psw, &tstw, &scw, &gmem_r,
                      &iswk, eval_force, eval_energy, bt, gbr_bt, gbd_bt, clash_minimum_distance,
                      clash_ratio);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyNonbondedKit<float,
                           float2> poly_nbk = poly_ag.getSinglePrecisionNonbondedKit(tier);
      MMControlKit<float> ctrl = mmctrl->spData(tier);
      CacheResourceKit<float> gmem_r = tb_space->spData(tier);
      ThermostatWriter<float> tstw = heat_bath->spData(tier);
      ISWorkspaceKit<float> iswk = ism_space->spData(poly_ps->getCyclePosition(), tier);
      launchNonbonded(nb_work_type, poly_nbk, poly_ser, &ctrl, &poly_psw, &tstw, &scw, &gmem_r,
                      &iswk, eval_force, eval_energy, force_sum, bt, gbr_bt, gbd_bt,
                      clash_minimum_distance, clash_ratio);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
extern void launchNonbonded(const PrecisionModel prec, const AtomGraphSynthesis &poly_ag,
                            const StaticExclusionMaskSynthesis &poly_se,
                            MolecularMechanicsControls *mmctrl, PhaseSpaceSynthesis *poly_ps,
                            Thermostat *heat_bath, ScoreCard *sc, CacheResource *tb_space,
                            ImplicitSolventWorkspace *ism_space, const EvaluateForce eval_force,
                            const EvaluateEnergy eval_energy, const CoreKlManager &launcher,
                            const double clash_minimum_distance, const double clash_ratio) {
  if (prec == PrecisionModel::DOUBLE || poly_ps->getForceAccumulationBits() <= 24) {
    launchNonbonded(prec, poly_ag, poly_se, mmctrl, poly_ps, heat_bath, sc, tb_space, ism_space,
                    eval_force, eval_energy, AccumulationMethod::SPLIT, launcher,
                    clash_minimum_distance, clash_ratio);
  }
  else {
    launchNonbonded(prec, poly_ag, poly_se, mmctrl, poly_ps, heat_bath, sc, tb_space, ism_space,
                    eval_force, eval_energy, AccumulationMethod::WHOLE, launcher,
                    clash_minimum_distance, clash_ratio);
  }
}

} // namespace energy
} // namespace stormm
