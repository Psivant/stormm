// -*-c++-*-
#include <cstdio>
#include <cstdlib>
#include <vector>
#ifdef STORMM_USE_CUDA
#include <cuda.h>
#include <curand_kernel.h>
#include <nvml.h>
#endif
#include "copyright.h"
#include "Reporting/error_format.h"
#include "hpc_config.h"

namespace stormm {
namespace card {

#ifdef STORMM_USE_CUDA
//-------------------------------------------------------------------------------------------------
GpuDetails::GpuDetails(const cudaDeviceProp &dev_prop, const int dev_index) :
    available{true},
    supported{true},
    arch_major{dev_prop.major},
    arch_minor{dev_prop.minor},
    smp_count{dev_prop.multiProcessorCount},
    card_ram{static_cast<int>(dev_prop.totalGlobalMem / static_cast<size_t>(constants::mega))},
    max_threads_per_block{dev_prop.maxThreadsPerBlock},
    max_threads_per_smp{dev_prop.maxThreadsPerMultiProcessor},
    max_blocks_per_smp{16},
    max_shared_per_block{static_cast<int>(dev_prop.sharedMemPerBlock)},
    max_shared_per_smp{static_cast<int>(dev_prop.sharedMemPerMultiprocessor)},
    global_cache_size{dev_prop.l2CacheSize},
    registers_per_block{dev_prop.regsPerBlock},
    registers_per_smp{dev_prop.regsPerMultiprocessor},
    card_name{std::string(dev_prop.name)}
{
  if (dev_prop.major < 3) {
    available = false;
    supported = false;
  }
  else if (dev_prop.major >= 3) {
    std::vector<nvmlProcessInfo_t> nvml_info(32);
    nvmlDevice_t nt_device;
    if (nvmlDeviceGetHandleByIndex_v2(dev_index, &nt_device) != NVML_SUCCESS) {
      rtWarn("Unable to get device handle for GPU + " + std::to_string(dev_index) +
             ".  This device (" + card_name + ") will not be accepted for program execution.",
             "GpuDetails");
      supported = false;
    }
    else {
      supported = true;
      unsigned int nvml_item_count = 0;
      nvmlReturn_t nv_status = nvmlDeviceGetComputeRunningProcesses(nt_device, &nvml_item_count,
                                                                    nvml_info.data());
      if (nv_status != NVML_SUCCESS && nv_status != NVML_ERROR_INSUFFICIENT_SIZE) {
        rtWarn("Unable to monitor activity on GPU " + std::to_string(dev_index) + " [error " +
               std::to_string(nv_status) + "]\n", "GpuDetails");
      }
      unsigned long long int mem_occ = 0;
      for (int i = 0; i < nvml_item_count; i++) {
        mem_occ += nvml_info[i].usedGpuMemory;
      }
      available = (mem_occ < significant_gpu_memory);
    }
  }
}

//-------------------------------------------------------------------------------------------------
extern const std::vector<GpuDetails> queryGpuStats(const ExceptionResponse policy) {

  // Test that there is a GPU in the system.  Initialize
  // a vector of specs to store all detected GPUs.
  int n_gpus;
  std::vector<GpuDetails> device_catalog;
  if (cudaGetDeviceCount(&n_gpus) != cudaSuccess || n_gpus == 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("No CUDA-capable devices were found.", "queryGpuStats");
    case ExceptionResponse::WARN:
      rtWarn("No CUDA-capable devices were found.  This will lead to errors if calls to an "
             "accelerator card are issued later in the program.", "queryGpuStats");
      return device_catalog;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  device_catalog.resize(n_gpus);
  
  // Activate zero-copy
  if (cudaSetDeviceFlags(cudaDeviceMapHost) != cudaSuccess) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Unable to establish cudaDeviceMapHost with cudaSetDeviceFlags().", "queryGpuStats");
    case ExceptionResponse::WARN:
      rtWarn("Unable to establish cudaDeviceMapHost with cudaSetDeviceFlags().  This may lead "
             "to errors later in the program.", "queryGpuStats");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }

  // Initialize the NVIDIA Management Library
  nvmlReturn_t nvml_init_result = nvmlInit_v2();
  std::string error_index;
  switch (nvml_init_result) {
  case NVML_ERROR_DRIVER_NOT_LOADED:
    error_index = "NVML_ERROR_DRIVER_NOT_LOADED";
  case NVML_ERROR_NO_PERMISSION:
    error_index = "NVML_ERROR_NO_PERMISSION";
  case NVML_ERROR_UNKNOWN:
    error_index = "NVML_ERROR_UNKNOWN";
  default:
    break;
  }
  if (error_index.size() > 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("nvmlInit_v2() failed with " + error_index + ".", "queryGpuStats");
    case ExceptionResponse::WARN:
      rtWarn("nvmlInit_v2() failed with " + error_index + ".  This may lead to errors later in "
             "the program.", "queryGpuStats");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }

  // Get device properties
  for (int i = 0; i < n_gpus; i++) {
    cudaDeviceProp device_properties;
    if (cudaGetDeviceProperties(&device_properties, i) != cudaSuccess) {
      rtWarn("Unable to query properties for GPU " + std::to_string(i) + ".  This device will not "
             "be accepted for program execution.", "queryGpuStats");
      continue;
    }

    // Transcribe information about this GPU
    device_catalog[i] = GpuDetails(device_properties, i);
  }

  // Shut down the NVIDIA Management Library
  if (nvmlShutdown() != NVML_SUCCESS) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Error executing nvmlShutdown().", "queryGpuStats");
    case ExceptionResponse::WARN:
      rtWarn("Error executing nvmlShutdown().", "queryGpuStats");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }

  return device_catalog;
}

//-------------------------------------------------------------------------------------------------
HpcConfig::HpcConfig(ExceptionResponse policy) :
    overall_gpu_count{0},
    available_gpu_count{0},
    supported_gpu_count{0},
    gpu_list{queryGpuStats(policy)}
{
  // Count valid and available GPUs
  overall_gpu_count = gpu_list.size();
  available_gpu_count = 0;
  supported_gpu_count = 0;
  for (int i = 0; i < overall_gpu_count; i++) {
    available_gpu_count += (gpu_list[i].getAvailability());
    supported_gpu_count += (gpu_list[i].getGpuSupported());
  }
  if (available_gpu_count == 0 && supported_gpu_count > 0) {
    const std::string user_msg = "No valid GPUs were available.  " +
                                 std::to_string(supported_gpu_count) +
                                 " GPUs were found to be occupied with other jobs.";
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr(user_msg, "HpcConfig");
    case ExceptionResponse::WARN:
      rtWarn(user_msg, "HpcConfig");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }

  // Crate handles for cuBLAS and cuSolver linear algebra packages
  cublasCreate(&cublas_handle);
  cusolverDnCreate(&cusolver_handle);
}

//-------------------------------------------------------------------------------------------------
HpcConfig::~HpcConfig() {

}
  
//-------------------------------------------------------------------------------------------------
int HpcConfig::getOverallGpuCount() const {
  return overall_gpu_count;
}

//-------------------------------------------------------------------------------------------------
int HpcConfig::getAvailableGpuCount() const {
  return available_gpu_count;
}

//-------------------------------------------------------------------------------------------------
int HpcConfig::getSupportedGpuCount() const {
  return supported_gpu_count;
}

//-------------------------------------------------------------------------------------------------
GpuDetails HpcConfig::getGpuInfo(const int gpu_index) const {
  return gpu_list[gpu_index];
}

//-------------------------------------------------------------------------------------------------
std::vector<int> HpcConfig::getGpuDevice(int requested_count) const {

  // Make a list of supported and available GPUs by their device indices
  std::vector<int> options;
  for (int i = 0; i < overall_gpu_count; i++) {
    if (gpu_list[i].getAvailability()) {
      options.push_back(i);
    }
  }
  
  // Loop over the list of available devices and mark them for CUDA calls
  if (cudaSetValidDevices(options.data(), available_gpu_count) != cudaSuccess) {
    cudaDeviceReset();
    rtErr("Error searching for compatible GPU.", "getGpuDevice");
  }

  // Establish the CUDA context
  if (cudaFree(0) != cudaSuccess) {
    cudaDeviceReset();
    rtErr("Error initializing the CUDA context with cudaFree(0).", "getGpuDevice");
  }

  // Get the device (this is a sanity check to ensure that the device can still be seen)
  int i = 0;
  std::vector<int> selections;
  while (i < requested_count) {
    int selected_device;
    if (cudaGetDevice(&selected_device) != cudaSuccess) {
      cudaDeviceReset();
      rtErr("Error selecting GPU.", "getGpuDevice");
    }

    // Set the device so that it will be used in all future calculations launched by this thread
    if (cudaSetDevice(selected_device) != cudaSuccess) {
      cudaDeviceReset();
      rtErr("Error setting GPU.", "getGpuDevice");
    }
    if (gpu_list[selected_device].getAvailability()) {
      selections.push_back(selected_device);
    }
    i++;
  }
  return selections;
}

//-------------------------------------------------------------------------------------------------
cublasHandle_t HpcConfig::getCuBlasHandle() const {
  return cublas_handle;
}

//-------------------------------------------------------------------------------------------------
cusolverDnHandle_t HpcConfig::getCuSolverHandle() const {
  return cusolver_handle;
}
#endif

} // namespace card
} // namespace stormm
