// -*-c++-*-
#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/hpc_bounds.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "kernel_format.h"

namespace stormm {
namespace card {

using constants::CaseSensitivity;
using constants::ExceptionResponse;
using parse::strcmpCased;
using stmath::findBin;

//-------------------------------------------------------------------------------------------------
KernelFormat::KernelFormat() :
  block_size_limit{1}, shared_usage{0}, block_dimension{1}, grid_dimension{1}, register_usage{0},
  kernel_name{std::string("")}
{}
  
//-------------------------------------------------------------------------------------------------
KernelFormat::KernelFormat(const int lb_max_threads_per_block, const int lb_min_blocks_per_smp,
                           const int register_usage_in, const int shared_usage_in,
                           const int block_subdivision, const GpuDetails &gpu,
                           const std::string &kernel_name_in) :
    block_size_limit{lb_max_threads_per_block},
    shared_usage{shared_usage_in},
    block_dimension{(block_size_limit / block_subdivision / warp_size_int) * warp_size_int},
    grid_dimension{block_subdivision * lb_min_blocks_per_smp * gpu.getSMPCount()},
    register_usage{register_usage_in},
    kernel_name{kernel_name_in}
{
  // Refine the register usage (this is unreliable with current cudart function calls, and should
  // not be trusted even after this step).
  const std::vector<int> register_break_points = {  0, 40, 48, 56, 64, 72, 80, 128, 256 };
  const std::vector<int> register_warp_counts  = { 48, 40, 36, 32, 28, 24, 16,   8 };
  const int register_bin = findBin(register_break_points, register_usage, ExceptionResponse::WARN);
  register_usage = register_warp_counts[register_bin];
}

//-------------------------------------------------------------------------------------------------
KernelFormat::KernelFormat(const int lb_max_threads_per_block, const int lb_min_blocks_per_smp,
                           const int register_usage_in, const int shared_usage_in,
                           const GpuDetails &gpu, const std::string &kernel_name_in) :
    KernelFormat(lb_max_threads_per_block, lb_min_blocks_per_smp, register_usage_in,
                 shared_usage_in, 1, gpu, kernel_name_in)
{}

#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
//-------------------------------------------------------------------------------------------------
KernelFormat::KernelFormat(const cudaFuncAttributes &attr, const int lb_min_blocks_per_smp,
                           const int block_subdivision, const GpuDetails &gpu,
                           const std::string &kernel_name_in) :
    KernelFormat(attr.maxThreadsPerBlock, lb_min_blocks_per_smp, attr.numRegs,
                 attr.sharedSizeBytes, block_subdivision, gpu, kernel_name_in)
{}
#  endif
#endif

//-------------------------------------------------------------------------------------------------
int2 KernelFormat::getLaunchParameters() const {
  return { grid_dimension, block_dimension };
}

//-------------------------------------------------------------------------------------------------
int KernelFormat::getRegisterUsage() const {
  return register_usage;
}

//-------------------------------------------------------------------------------------------------
int KernelFormat::getBlockSizeLimit() const {
  return block_size_limit;
}

//-------------------------------------------------------------------------------------------------
int KernelFormat::getSharedMemoryRequirement() const {
  return shared_usage;
}

//-------------------------------------------------------------------------------------------------
const std::string& KernelFormat::getKernelName() const {
  return kernel_name;
}

//-------------------------------------------------------------------------------------------------
KernelManager::KernelManager(const GpuDetails &gpu_in) :
    gpu{gpu_in}, k_dictionary{}
{}

//-------------------------------------------------------------------------------------------------
KernelManager::~KernelManager() {
}

//-------------------------------------------------------------------------------------------------
const GpuDetails& KernelManager::getGpu() const {
  return gpu;
}

//-------------------------------------------------------------------------------------------------
void KernelManager::printLaunchParameters(const std::string &k_key,
                                          const char* true_class_name) const {
  if (k_key.size() == 0) {
    return;
  }
  if (strcmpCased(k_key, "all", CaseSensitivity::NO)) {
    for (auto it = k_dictionary.begin(); it != k_dictionary.end(); it++) {
      printLaunchParameters(it->first);
    }
    return;
  }
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("No kernel with identifier " + k_key + " is known.", true_class_name,
          "printLaunchParameters");
  }
  const int2 lp = k_dictionary.at(k_key).getLaunchParameters();
  const int mtpb = k_dictionary.at(k_key).getBlockSizeLimit();
  const int nreg = k_dictionary.at(k_key).getRegisterUsage();
  printf("  %12.12s :: %4d blocks, %4d threads, %3d registers (%4d SMP, %4d max threads)\n",
         k_key.c_str(), lp.x, lp.y, nreg, gpu.getSMPCount(), mtpb);
}

} // namespace card
} // namespace stormm
