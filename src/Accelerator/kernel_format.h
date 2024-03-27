// -*-c++-*-
#ifndef STORMM_KERNEL_MANAGER_H
#define STORMM_KERNEL_MANAGER_H

#include <map>
#include <string>
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
#    include <cuda_runtime.h>
#  endif
#endif
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "gpu_details.h"

namespace stormm {
namespace card {

#ifndef STORMM_USE_HPC
using data_types::int2;
#endif

/// \brief The name of a generic kernel manager
const char generic_kernel_manager_name[] = "KernelManager";
  
/// \brief Encapsulate the operations to store and retrieve information about a kernel's format.
class KernelFormat {
public:

  /// \brief The constructor takes launch bounds and other information that can be plucked from a
  ///        cudaFuncAttributes object.
  ///
  /// Overloaded:
  ///   - Construct a blank object
  ///   - Provide explicit instructions on whether to consider breaking up the blocks into smaller
  ///     units
  ///   - Assume that the largest possible block size is always to be used
  ///
  /// \param lb_max_threads_per_block  Maximum threads per block, as stated in the launch bounds
  /// \param lb_min_blocks_per_smp     Minimum blocks per multiprocessor, from the launch bounds
  /// \param register_usage_in         Input register usage
  /// \param shared_usage_in           Input __shared__ memory usage
  /// \param block_subdivision         Preferred block multiplicity (this will compound the input
  ///                                  minimum number of blocks per multiprocessor)
  /// \param attr                      Result of a CUDA runtime query to get kernel specifications
  /// \param gpu                       Details of the available GPU (likely passed in from a
  ///                                  CoreKlManager struct containing many KernelFormat objects)
  /// \param kernel_name_in            Name of the kernel, for reporting purposes later (optional)
  /// \{
  KernelFormat();
  
  KernelFormat(int lb_max_threads_per_block, int lb_min_blocks_per_smp, int register_usage_in,
               int shared_usage_in, int block_subdivision, const GpuDetails &gpu,
               const std::string &kernel_name_in = std::string(""));

  KernelFormat(int lb_max_threads_per_block, int lb_min_blocks_per_smp, int register_usage_in,
               int shared_usage_in, const GpuDetails &gpu,
               const std::string &kernel_name_in = std::string(""));

#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  KernelFormat(const cudaFuncAttributes &attr, int lb_min_blocks_per_smp, int block_subdivision,
               const GpuDetails &gpu, const std::string &kernel_name_in = std::string(""));
#  endif
#endif
  /// \}

  /// \brief Take the default copy and move constructors as well as assignment operators.
  /// \{
  KernelFormat(const KernelFormat &original) = default;
  KernelFormat(KernelFormat &&original) = default;
  KernelFormat& operator=(const KernelFormat &other) = default;
  KernelFormat& operator=(KernelFormat &&other) = default;
  /// \}
  
  /// \brief Get the optimal block and grid sizes for kernel launches with the present GPU.
  int2 getLaunchParameters() const;

  /// \brief Get the register usage of the kernel.
  int getRegisterUsage() const;

  /// \brief Get the maximum thread count for a single block in the kernel launch.
  int getBlockSizeLimit() const;

  /// \brief Get the amount of __shared__ memory needed by any one block.
  int getSharedMemoryRequirement() const;

  /// \brief Get the name of this kernel
  const std::string& getKernelName() const;
  
private:
  int block_size_limit;        ///< The largest block size usable by the kernel launch (exceeding
                               ///<   this will cause the kernel launch to fail)
  int shared_usage;            ///< The maximum amount of __shared__ memory needed by each block
  int block_dimension;         ///< Computed optimal block dimension to use in kernel launches
  int grid_dimension;          ///< Computed optimal grid size to use in kernel launches
  int register_usage;          ///< The number of registers needed by each thread of the kernel
                               ///<   as it is compiled for the current executable
  std::string kernel_name;     ///< Name of the kernel, for reporting purposes
};

/// \brief Parent class for other kernel managers, incorporating the common dictionary of kernel
///        keys and GPU details.
class KernelManager {
public:

  /// \brief Get the GPU information for the active GPU.
  const GpuDetails& getGpu() const;

  /// \brief Print out the kernel launch parameters found for this workload.
  ///
  /// \param k_key  Identifier string of the kernel for which to print the parameters (if blank,
  ///               no kernels' parameters will be printed, and if "ALL" (case-insensitive), all
  ///               kernels' parameters will be printed)
  void printLaunchParameters(const std::string &k_key = std::string(""),
                             const char* true_class_name = generic_kernel_manager_name) const;

protected:

  /// \brief The constructor for this base class takes the GPU specifications
  KernelManager(const GpuDetails &gpu_in = null_gpu);

  /// \brief A virtual destructor ensures proper behavior in the destructors of derived classes
  ///        for managing particular groups of kernels.
  virtual ~KernelManager();

  /// The details of the GPU in use are simply copied into this object.
  GpuDetails gpu;

  /// Store the resource requirements and selected launch parameters for a variety of kernels.
  /// Keys are determined according to the free functions further on in this library.
  std::map<std::string, KernelFormat> k_dictionary;
};
  
} // namespace card
} // namespace stormm

#endif
