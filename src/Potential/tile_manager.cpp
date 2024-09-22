#include "copyright.h"
#include "tile_manager.h"

#include "Math/formulas.h"
#include "Constants/hpc_bounds.h"

namespace stormm {
namespace energy {

using card::HybridKind;
using stmath::ipow;
  
//-------------------------------------------------------------------------------------------------
TilePlan::TilePlan(const int nchoice_in, const int* read_assign_in, const int* self_assign_in,
                   const int* reduce_prep_in, const int* self_prep_in, const float* scalings_in,
                   const int* nt_stencil_in, int* xfrc_ovrf_in, int* yfrc_ovrf_in,
                   int* zfrc_ovrf_in) :
    nchoice{nchoice_in}, read_assign{read_assign_in}, self_assign{self_assign_in},
    reduce_prep{reduce_prep_in}, self_prep{self_prep_in}, scalings{scalings_in},
    nt_stencil{nt_stencil_in}, xfrc_ovrf{xfrc_ovrf_in}, yfrc_ovrf{yfrc_ovrf_in},
    zfrc_ovrf{zfrc_ovrf_in}
{}

//-------------------------------------------------------------------------------------------------
TileManager::TileManager(const int2 launch_parameters, const int max_deg_in) :
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
    maximum_degeneracy{(max_deg_in <= 0) ? 2 : max_deg_in},
#  else
    maximum_degeneracy{(max_deg_in <= 0) ? 2 : max_deg_in},
#  endif
#else
    maximum_degeneracy{(max_deg_in <= 0) ? 2 : max_deg_in},
#endif
    block_count{launch_parameters.x},
    thread_count{launch_parameters.y},
    read_assignments{HybridKind::POINTER, "tile_recv_reads"},
    self_assignments{HybridKind::POINTER, "tile_self_reads"},
    reduce_preparations{HybridKind::POINTER, "tile_reduce_prep"},
    self_preparations{HybridKind::POINTER, "tile_slfrdc_prep"},
    thread_scalings{HybridKind::ARRAY, "tile_scalings"},
    tower_plate_stencil{HybridKind::POINTER, "tile_stencil"},
    x_force_overflow{HybridKind::POINTER, "tile_xfrc_ovrf"},
    y_force_overflow{HybridKind::POINTER, "tile_yfrc_ovrf"},
    z_force_overflow{HybridKind::POINTER, "tile_zfrc_ovrf"},
    int_data{HybridKind::ARRAY, "tile_int_storage"}
{
  // Check that the maximum degeneracy will not break the tile scheme.  The degeneracy cannot
  // grow larger than the point at which even one iteration of a tile would result in
  // double-counting.
  const int min_unique_atoms = warp_size_int / ipow(2, maximum_degeneracy);
  if (min_unique_atoms * min_unique_atoms < warp_size_int) {
    rtErr("The maximum degeneracy of a tile batch is too great to guarantee that each thread can "
          "have a unique interaction to compute.", "TileManager");
  }

  // Allocate array space
  allocate();

  // Loop over all possible subdivisions of the sending and receiving atom batches.
  int sending_division = 1;
  std::vector<int> batch_assignments(warp_size_int);
  int* red_prep_ptr = reduce_preparations.data();
  int* slf_prep_ptr = self_preparations.data();
  for (int i = 0; i <= maximum_degeneracy; i++) {
    const int snd_unique_atoms = warp_size_int / sending_division;
    int recving_division = 1;
    for (int j = 0; j <= maximum_degeneracy; j++) {
      const int rcv_unique_atoms = warp_size_int / recving_division;

      // Compute the layout by which the warp will read a batch of receiving atoms with the given
      // degeneracies of both sending and receiving atoms.
      for (int k = 0; k < recving_division; k++) {
        const std::vector<int> recv_order = computeStaggeredOrder(rcv_unique_atoms,
                                                                  sending_division, k);
        for (int m = 0; m < rcv_unique_atoms; m++) {
          batch_assignments[(k * rcv_unique_atoms) + m] = recv_order[m];
        }
      }
      const int ij_offset = (((maximum_degeneracy + 1) * i) + j) * warp_size_int;
      read_assignments.putHost(batch_assignments, ij_offset, warp_size_int);
      if (i != j) {
        self_assignments.putHost(batch_assignments, ij_offset, warp_size_int);
      }
      else {

        // A warp of 32 threads divided into quarters creates an 8 x 8 tile.  A warp of 64
        // threads divided into quarters creates a 16 x 16 tile. (A warp of 16 threads cannot be
        // divided into quarters for effective use.) A warp of any size that is not subdivided
        // creates a tile of warp_size_int x warp_size_int.  In all cases there are some
        // interactions which will be double-counted as the self-interactions tile is evaluated,
        // and it is desirable to have this double-counting occur in a predictable way.  The
        // common pattern for self-interaction tiles?  Shift the last subdivision of the warp by
        // half the subdivision width--this is a priority to ensure that all double-counted
        // interactions occur in the first iteration of the tile evaluation.  Beyond that, shift
        // the first subdivision of the warp +1 (if the warp is subdivided at all), then shift
        // subsequent subdivisions of the warp by an additional N positions each, where N is the
        // number of iterations that will be needed to fully evaluate the tile (shift by N so that
        // subsequent iterations do not run over previously computed interations).
        std::vector<int> self_batch_assignments(warp_size_int);
        const int niter = rcv_unique_atoms * rcv_unique_atoms / warp_size_int / 2;
        for (int k = 0; k < recving_division; k++) {
          int kmcon = (k == recving_division - 1) ? (rcv_unique_atoms / 2) : (niter * k) + 1;
          for (int m = 0; m < rcv_unique_atoms; m++) {
            self_batch_assignments[(k * rcv_unique_atoms) + m] = kmcon;
            kmcon++;
            if (kmcon == rcv_unique_atoms) {
              kmcon = 0;
            }
          }
        }
        self_assignments.putHost(self_batch_assignments, ij_offset, warp_size_int);
        
        // Compute the preparatory assignments for self-interacting tiles.  These assignments
        // depend on the special batch read assignments computed within this local scope, and a
        // reduced tile depth.
        const int self_tile_depth = (snd_unique_atoms * snd_unique_atoms) / twice_warp_size_int;
        for (int k = 0; k < self_tile_depth - 1; k++) {
          const int zero_idx = self_batch_assignments[0];
          for (int m = 0; m < warp_bits_mask_int; m++) {
            self_batch_assignments[m] = self_batch_assignments[m + 1];
          }
          self_batch_assignments[warp_bits_mask_int] = zero_idx;
        }
        
        // In tiles where the sending and receiving atoms are not the same, the receiving atoms
        // count from 0 to N-1 in the first N slots, where N is the number of unique receiving
        // atoms however the warp is subdivided.  This is a natural arrangement for thinking of
        // the way the writeback will go and makes it easy to define the other lanes that each
        // thread must look to in order to prepare for the reduction.  However, for the self
        // interactions, the receiving atoms in the first subdivision are shifted forward by one
        // lane.  Because it is convenient to mark down image indices on the GPU with based on
        // the initial assignments, those threads must collect back the relevant data from threads
        // in other subdivisions of the warp.
        for (int k = 0; k < snd_unique_atoms; k++) {
          const int rel_seek = self_assignments.readHost(ij_offset + k);
          int subdiv = 0;
          for (int m = 0; m < warp_size_int; m++) {
            if (self_batch_assignments[m] == rel_seek) {
              slf_prep_ptr[(subdiv * snd_unique_atoms) + k + ij_offset] = m;
              subdiv++;
            }
          }
        }
      }
      
      // Jog the arrangement forward, as if processing a tile for the required number of
      // iterations, in order to determine the final resting places of each atom replica and thus
      // the manner in which the results must be re-arranged in order to be ready for force
      // reduction.
      const int tile_depth = (snd_unique_atoms * rcv_unique_atoms) / warp_size_int;
      for (int k = 0; k < tile_depth - 1; k++) {

        // The warp bits mask is equal to the warp size minus one.
        const int zero_idx = batch_assignments[0];
        for (int m = 0; m < warp_bits_mask_int; m++) {
          batch_assignments[m] = batch_assignments[m + 1];
        }
        batch_assignments[warp_bits_mask_int] = zero_idx;
      }
      std::vector<int> copies_fulfilled(rcv_unique_atoms, 0);
      for (int k = 0; k < warp_size_int; k++) {
        const int rel_idx = batch_assignments[k];
        const int reduction_ready_idx = (copies_fulfilled[rel_idx] * rcv_unique_atoms) + rel_idx;
        copies_fulfilled[rel_idx] += 1;
        red_prep_ptr[reduction_ready_idx + ij_offset] = k;
        if (i != j) {
          slf_prep_ptr[reduction_ready_idx + ij_offset] = k;
        }
      }

      // Fill out the scaling factors
      std::vector<float> tmp_scale(warp_size_int);
      for (int k = 0; k < recving_division; k++) {
        const float tk = (k == recving_division - 1) ? 0.5 : 1.0;
        for (int m = 0; m < rcv_unique_atoms; m++) {
          tmp_scale[(k * rcv_unique_atoms) + m] = tk;
        }
      }
      thread_scalings.putHost(tmp_scale, ij_offset, warp_size_int);

      // Increment the number of subdivisions
      recving_division *= 2;
    }
    sending_division *= 2;
  }

  // Lay out the tower-plant neutral territory stencil
  planTowerPlateStencil();

  // Upload the plans immediately
#ifdef STORMM_USE_HPC
  this->upload();
#endif
}

//-------------------------------------------------------------------------------------------------
TileManager::TileManager(const TileManager &original) :
    maximum_degeneracy{original.maximum_degeneracy},
    block_count{original.block_count},
    thread_count{original.thread_count},
    read_assignments{original.read_assignments},
    self_assignments{original.self_assignments},
    reduce_preparations{original.reduce_preparations},
    self_preparations{original.self_preparations},
    thread_scalings{original.thread_scalings},
    tower_plate_stencil{original.tower_plate_stencil},
    x_force_overflow{original.x_force_overflow},
    y_force_overflow{original.y_force_overflow},
    z_force_overflow{original.z_force_overflow},
    int_data{original.int_data}
{
  // The allocation function will take care of pointer rebasing
  allocate();
}

//-------------------------------------------------------------------------------------------------
TileManager::TileManager(TileManager &&original) :
    maximum_degeneracy{original.maximum_degeneracy},
    block_count{original.block_count},
    thread_count{original.thread_count},
    read_assignments{std::move(original.read_assignments)},
    self_assignments{std::move(original.self_assignments)},
    reduce_preparations{std::move(original.reduce_preparations)},
    self_preparations{std::move(original.self_preparations)},
    thread_scalings{std::move(original.thread_scalings)},
    tower_plate_stencil{std::move(original.tower_plate_stencil)},
    x_force_overflow{std::move(original.x_force_overflow)},
    y_force_overflow{std::move(original.y_force_overflow)},
    z_force_overflow{std::move(original.z_force_overflow)},
    int_data{std::move(original.int_data)}
{}

//-------------------------------------------------------------------------------------------------
TileManager& TileManager::operator=(const TileManager &other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }
  maximum_degeneracy = other.maximum_degeneracy;
  block_count = other.block_count;
  thread_count = other.thread_count;
  read_assignments = other.read_assignments;
  self_assignments = other.self_assignments;
  reduce_preparations = other.reduce_preparations;
  self_preparations = other.self_preparations;
  thread_scalings = other.thread_scalings;
  tower_plate_stencil = other.tower_plate_stencil;
  x_force_overflow = other.x_force_overflow;
  y_force_overflow = other.y_force_overflow;
  z_force_overflow = other.z_force_overflow;
  int_data = other.int_data;

  // The allocation function will take care of pointer rebasing
  allocate();
  return *this;
}

//-------------------------------------------------------------------------------------------------
TileManager& TileManager::operator=(TileManager &&other) {

  // Guard against self assignment
  if (this == &other) {
    return *this;
  }
  maximum_degeneracy = other.maximum_degeneracy;
  block_count = other.block_count;
  thread_count = other.thread_count;
  read_assignments = std::move(other.read_assignments);
  self_assignments = std::move(other.self_assignments);
  reduce_preparations = std::move(other.reduce_preparations);
  self_preparations = std::move(other.self_preparations);
  thread_scalings = std::move(other.thread_scalings);
  tower_plate_stencil = std::move(other.tower_plate_stencil);
  x_force_overflow = std::move(other.x_force_overflow);
  y_force_overflow = std::move(other.y_force_overflow);
  z_force_overflow = std::move(other.z_force_overflow);
  int_data = std::move(other.int_data);

  // No pointer repair is needed in most move constructors or assignment operators
  return *this;
}

//-------------------------------------------------------------------------------------------------
int TileManager::getMaximumBatchDegeneracy() const {
  return maximum_degeneracy;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> TileManager::getSendingAtomLayout(const int sending_atom_degeneracy) const {
  if (warp_size_int % sending_atom_degeneracy != 0 || sending_atom_degeneracy > warp_size_int) {
    rtErr("A degeneracy of " + std::to_string(sending_atom_degeneracy) + " is invalid for a warp "
          "of size " + std::to_string(warp_size_int) + ".", "TileManager", "getSendingAtomLayout");
  }
  std::vector<int> result(warp_size_int);
  const int unique_natom = warp_size_int / sending_atom_degeneracy;
  const int nrep = warp_size_int / unique_natom;
  size_t k = 0;
  for (int i = 0; i < nrep; i++) {
    for (int j = 0; j < unique_natom; j++) {
      result[k] = j;
      k++;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> TileManager::getReadAssignments(const int sending_atom_degeneracy,
                                                 const int recving_atom_degeneracy) const {
  const int lsend_deg = round(log2(sending_atom_degeneracy));
  const int lrecv_deg = round(log2(recving_atom_degeneracy));
  return read_assignments.readHost(((lsend_deg * (maximum_degeneracy + 1)) + lrecv_deg) *
                                   warp_size_int, warp_size_int);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> TileManager::getSelfAssignments(const int sending_atom_degeneracy,
                                                 const int recving_atom_degeneracy) const {
  const int lsend_deg = round(log2(sending_atom_degeneracy));
  const int lrecv_deg = round(log2(recving_atom_degeneracy));
  return self_assignments.readHost(((lsend_deg * (maximum_degeneracy + 1)) + lrecv_deg) *
                                   warp_size_int, warp_size_int);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> TileManager::getReductionPreparations(const int sending_atom_degeneracy,
                                                       const int recving_atom_degeneracy) const {
  const int lsend_deg = round(log2(sending_atom_degeneracy));
  const int lrecv_deg = round(log2(recving_atom_degeneracy));
  return reduce_preparations.readHost(((lsend_deg * (maximum_degeneracy + 1)) + lrecv_deg) *
                                      warp_size_int, warp_size_int);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> TileManager::getSelfPreparations(const int sending_atom_degeneracy,
                                                  const int recving_atom_degeneracy) const {
  const int lsend_deg = round(log2(sending_atom_degeneracy));
  const int lrecv_deg = round(log2(recving_atom_degeneracy));
  return self_preparations.readHost(((lsend_deg * (maximum_degeneracy + 1)) + lrecv_deg) *
                                    warp_size_int, warp_size_int);
}

//-------------------------------------------------------------------------------------------------
std::vector<float> TileManager::getThreadScalings(const int atom_degeneracy) const {
  const int l_deg = round(log2(atom_degeneracy));
  return thread_scalings.readHost(l_deg * (maximum_degeneracy + 1) * warp_size_int, warp_size_int);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> TileManager::getNeutralTerritoryStencil() const {
  return tower_plate_stencil.readHost();
}

//-------------------------------------------------------------------------------------------------
TilePlan TileManager::data(const HybridTargetLevel tier) {
  return TilePlan(maximum_degeneracy + 1, read_assignments.data(tier), self_assignments.data(tier),
                  reduce_preparations.data(tier), self_preparations.data(tier),
                  thread_scalings.data(tier), tower_plate_stencil.data(tier),
                  x_force_overflow.data(tier), y_force_overflow.data(tier),
                  z_force_overflow.data(tier));
}
  
#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void TileManager::upload() {
  int_data.upload();
  thread_scalings.upload();
}

//-------------------------------------------------------------------------------------------------
void TileManager::download() {
  int_data.download();
  thread_scalings.download();
}
#endif

//-------------------------------------------------------------------------------------------------
void TileManager::allocate() {
  const int table_len = (maximum_degeneracy + 1) * (maximum_degeneracy + 1) * warp_size_int;
  const int grid_len = block_count * thread_count;
  int_data.resize((4 * table_len) + 32 + (3 * grid_len));
  read_assignments.setPointer(&int_data,                                        0, table_len);
  self_assignments.setPointer(&int_data,         table_len                       , table_len);
  reduce_preparations.setPointer(&int_data,  2 * table_len                       , table_len);
  self_preparations.setPointer(&int_data,    3 * table_len                       , table_len);
  tower_plate_stencil.setPointer(&int_data,  4 * table_len                       , 32);
  x_force_overflow.setPointer(&int_data,    (4 * table_len) +                  32, grid_len);
  y_force_overflow.setPointer(&int_data,    (4 * table_len) +      grid_len  + 32, grid_len);
  z_force_overflow.setPointer(&int_data,    (4 * table_len) + (2 * grid_len) + 32, grid_len);
  thread_scalings.resize(table_len);
}

//-------------------------------------------------------------------------------------------------
void TileManager::planTowerPlateStencil() {
  std::vector<int> pln(32, 0);

  // Each member of the 32-element plan (best for an NVIDIA 32-threaded warp, but also amenable
  // to other architectures) will be a bit-packed int encoding the moves relative to some central
  // cell.  The exact cell index for any given lane will be calculated in the context of the
  // central neighbor list cell.
  for (int i = 0; i < 5; i++) {
    pln[i] = (i << 16);
  }
  int pln_idx = 16; 
  for (int j = 0; j < 3; j++) {
    for (int i = 0; i < 5; i++) {
      if (j < 2 || i < 2) {
        pln[pln_idx] = ((2 << 16) | (j << 8) | i);
        pln_idx++;
      }
    }
  }
  tower_plate_stencil.putHost(pln);
}

//-------------------------------------------------------------------------------------------------
std::vector<int> TileManager::computeStaggeredOrder(const int natom, const int sending_reps,
                                                    const int iter) {
  const int start_index = iter * natom;
  const int sending_batch_size = warp_size_int / sending_reps;
  const int rep_idx = start_index / std::max(sending_batch_size, natom);
  const int offset = rep_idx * (natom / sending_reps);
  std::vector<int> result(natom);
  for (int i = 0; i < natom; i++) {
    result[i] = offset + i;
    result[i] -= (result[i] >= natom) * natom;
  }
  return result;
}

} // namespace energy
} // namespace stormm
