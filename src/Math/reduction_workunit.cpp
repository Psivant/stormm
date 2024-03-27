#include "copyright.h"
#include "Constants/hpc_bounds.h"
#include "Reporting/error_format.h"
#include "reduction_workunit.h"
#include "rounding.h"
#include "summation.h"
#include "vector_ops.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
ReductionWorkUnit::ReductionWorkUnit(const int atom_start_in, const int atom_end_in,
                                     const int result_index_in, const int dependency_start_in,
                                     const int dependency_end_in, const int system_index_in) :
  atom_start{atom_start_in}, atom_end{atom_end_in}, result_index{result_index_in},
  dependency_start{dependency_start_in}, dependency_end{dependency_end_in},
  system_index{system_index_in}
{}

//-------------------------------------------------------------------------------------------------
int ReductionWorkUnit::getAtomStart() const {
  return atom_start;
}

//-------------------------------------------------------------------------------------------------
int ReductionWorkUnit::getAtomEnd() const {
  return atom_end;
}

//-------------------------------------------------------------------------------------------------
int ReductionWorkUnit::getResultIndex() const {
  return result_index;
}

//-------------------------------------------------------------------------------------------------
int ReductionWorkUnit::getDependencyStart() const {
  return dependency_start;
}

//-------------------------------------------------------------------------------------------------
int ReductionWorkUnit::getDependencyEnd() const {
  return dependency_end;
}

//-------------------------------------------------------------------------------------------------
int ReductionWorkUnit::getSystemIndex() const {
  return system_index;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> ReductionWorkUnit::getAbstract() const {
  std::vector<int> result = { atom_start, atom_end, result_index, dependency_start,
                              dependency_end, system_index, 0, 0 };
  return result;
}

//-------------------------------------------------------------------------------------------------
int optReductionKernelSubdivision(const int* atom_counts, const int n_systems,
                                  const GpuDetails &gpu) {
  const int nsmp = gpu.getSMPCount();
  if (n_systems > nsmp * 16) {
    return 2;
  }
  const int total_atoms = sum<int>(atom_counts, n_systems);
  if (total_atoms > nsmp * 16 * 256) {
    return 2;
  }
  else {
    return 1;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
int optReductionKernelSubdivision(const std::vector<int> &atom_counts, const GpuDetails &gpu) {
  return optReductionKernelSubdivision(atom_counts.data(), atom_counts.size(), gpu);
}

//-------------------------------------------------------------------------------------------------
int optReductionKernelSubdivision(const Hybrid<int> &atom_counts, const GpuDetails &gpu) {
  return optReductionKernelSubdivision(atom_counts.data(), atom_counts.size(), gpu);
}

//-------------------------------------------------------------------------------------------------
std::vector<ReductionWorkUnit> buildReductionWorkUnits(const std::vector<int> &atom_starts,
                                                       const std::vector<int> &atom_counts,
                                                       const GpuDetails &gpu,
                                                       const int tasks_per_atom) {
  if (atom_starts.size() != atom_counts.size()) {
    rtErr("Starting indices were provided for the atoms of " + std::to_string(atom_starts.size()) +
          " systems, but atom counts for " + std::to_string(atom_counts.size()) + " systems were "
          "provided.", "buildReductionWorkUnits");
  }
  const int nsys = atom_starts.size();
  
  // Compare the number of systems to the number of streaming multiprocessors on the GPU, when
  // blocks of various sizes are used.  Compute the average number of tasks (atoms times
  // tasks_per_atom) for each system and assume that the system size stays close to the average
  // throughout the workload.
  double overall_tasks = 0.0;
  int overall_atoms = 0;
  for (int i = 0; i < nsys; i++) {
    overall_tasks += static_cast<double>(atom_counts[i]);
    overall_atoms += atom_counts[i];
  }
  overall_tasks *= static_cast<double>(tasks_per_atom);
  const int max_natom = maxValue(atom_counts);
  const double max_system_tasks = static_cast<double>(max_natom) *
                                  static_cast<double>(tasks_per_atom);
  const double average_tasks = overall_tasks / static_cast<double>(nsys);
  const int n_smp = gpu.getSMPCount();

  // It is assumed that each streaming multiprocessor can accept at least the large block size
  // worth of threads.  If this assumption is violated, the launch grid will not entirely fit on
  // the GPU at one time, but that is not, in principle, a problem.
  const int n_reduction_threads = n_smp * large_block_size;
  const int n_reduction_blocks = n_reduction_threads / small_block_size;
  
  // Reduction workloads, even if fused with other work as part of a more elaborate kernel, will
  // require extreme rates of global memory transactions.  Let the reduction thread blocks take a
  // size of 256 threads (small_block_size) and let up to maximum_gathering_results (1024) thread
  // blocks work on the reduction operations for any one system.  That will give the GPUs plenty
  // of room to grow, and will ensure that even a single very large system will be able to have a
  // card-saturating 262,144 threads handling its reductions.
  const double tasks_per_thread = overall_tasks / static_cast<double>(n_reduction_threads);
  
  // Say there are 2.25 tasks per thread.  Each system contains its atom count times the number
  // of tasks per atom (tasks_per_atom).  Therefore each block should be responsible for two such
  // tasks, and every fourth block should be responsible for three.  If, however, the number of
  // tasks per thread is less than twice the number of tasks per atom, it is more efficient to
  // forego a second kernel call and have one one block be responsible for the entirety of the
  // computations in any given system.  This block can accomplish both gathering and scattering
  // oprations in a single step.  Otherwise, assign a consistent number of atoms to each block and
  // a number of blocks to each system sufficient to cover all atoms in it.
  std::vector<ReductionWorkUnit> result;  
  if (max_system_tasks < 2 * small_block_size * tasks_per_atom) {
    result.reserve(nsys);
    for (int i = 0; i < nsys; i++) {
      result.emplace_back(atom_starts[i], atom_starts[i] + atom_counts[i], i, i, i + 1, i);
    }
  }
  else {

    // Each block must deal with what is probably a subset of the atoms in any given system, but
    // how many should it target?  If the warp size is 32, dealing with 23 atoms is the same as
    // dealing with 32, 39 might as well be 64.  Furthermore, each block will be able to assign
    // warps to different tasks related to each atom, but for a given number of reduction tasks
    // and atoms assigned to a block there may or may not be idle warps in each block.  Also,
    // each system can have a limited number of separate blocks working on it, to limit the
    // breadth of secondary reductions that must occur in a gather / scatter over two kernels.
    // The largest system informs the minimum number of atoms that each block must deal with.
    const int min_atoms_per_block = roundUp((max_natom + maximum_gathering_results - 1) /
                                            maximum_gathering_results, warp_size_int);
    int atoms_per_block = std::min(min_atoms_per_block,
                                   (overall_atoms + n_reduction_blocks - 1) / n_reduction_blocks);
    atoms_per_block = roundUp(atoms_per_block, warp_size_int);
    const int imin = std::max(atoms_per_block - warp_size_int, warp_size_int);
    const int imax = atoms_per_block + warp_size_int;
    int best_efficiency = 0.0;
    int best_atoms_per_block = atoms_per_block;
    const int warps_per_block = small_block_size / warp_size_int;
    for (int i = imin; i <= imax; i++) {

      // Compute the efficiency of a block assigned to reduce the standard number of atoms
      const int warps_occupied  = (i / warp_size_int) * tasks_per_atom;
      const int baseline_cycles = (warps_occupied + warps_per_block - 1) / warps_per_block;
      const double baseline_efficiency = static_cast<double>(warps_occupied) /
                                         static_cast<double>(baseline_cycles * warps_per_block);
      double active_efficiency = 0.0;
      int total_blocks_alloted = 0;
      for (int j = 0; j < nsys; j++) {
        const int blk_allotment  = (atom_counts[j] + i - 1) / i;
        const int last_blk_atoms = atom_counts[j] - ((blk_allotment - 1) * i);
        const int last_blk_warps_occupied = ((last_blk_atoms + warp_size_int - 1) /
                                             warp_size_int) * tasks_per_atom;
        const int last_blk_cycles = (last_blk_warps_occupied + warps_per_block - 1) /
                                    warps_per_block;
        const double last_blk_efficiency = static_cast<double>(last_blk_warps_occupied) /
                                           static_cast<double>(last_blk_cycles * warps_per_block);
        active_efficiency += (baseline_efficiency * static_cast<double>(blk_allotment - 1)) +
                             last_blk_efficiency;
        total_blocks_alloted += blk_allotment;
      }
      const int grid_cycles = (total_blocks_alloted + n_reduction_blocks - 1) / n_reduction_blocks;
      const double grid_efficiency = static_cast<double>(total_blocks_alloted) /
                                     static_cast<double>(grid_cycles * n_reduction_blocks);
      active_efficiency /= static_cast<double>(total_blocks_alloted);
      const double overall_efficiency = grid_efficiency * active_efficiency;
      if (overall_efficiency >= best_efficiency) {
        best_efficiency = overall_efficiency;
        best_atoms_per_block = i;
      }
    }

    // Use the (arduously computed) optimal block size to make work units for each system
    int nredwu = 0;
    for (int i = 0; i < nsys; i++) {
      nredwu += (atom_counts[i] + best_atoms_per_block - 1) / best_atoms_per_block;
    }
    result.reserve(nredwu);
    int result_con = 0;
    for (int i = 0; i < nsys; i++) {
      const int dep_start = result_con;
      const int dep_end   = result_con + ((atom_counts[i] + best_atoms_per_block - 1) /
                                         best_atoms_per_block);
      for (int j = 0; j < atom_counts[i]; j += best_atoms_per_block) {
        const int jlimit = atom_starts[i] + std::min(j + best_atoms_per_block, atom_counts[i]);
        result.emplace_back(atom_starts[i] + j, jlimit, result_con, dep_start, dep_end, i);
        result_con++;
      }
    }
  }
  
  return result;
}

} // namespace stmath
} // namespace stormm
