#include <vector>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include "copyright.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/UnitTesting/unit_test.h"

using stormm::constants::tiny;
using stormm::constants::ExceptionResponse;
using stormm::random::Ran2Generator;
using stormm::review::stormmSplash;
using stormm::review::stormmWatermark;
using namespace stormm::card;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Obtain environment variables or command-line input, if available
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Section 1
  section("Detecting GPU specifications");

  // Section 2
  section("cuBLAS bindings");

  // Section 3
  section("cuSolver bindings");

  // Get stats on the available GPUs
  section(1);
  HpcConfig gpu_config(ExceptionResponse::WARN);
  const int n_gpus = gpu_config.getOverallGpuCount();
  const int n_avail = gpu_config.getAvailableGpuCount();
  TestPriority gpu_calc_tests = (n_avail > 0) ? TestPriority::CRITICAL : TestPriority::ABORT;
  check(n_gpus, RelationalOperator::GE, n_avail, "The number of available GPU devices exceeds the "
        "total number of GPUs in the machine.");

  // Switch over the various architectures: check that various GPU stats are correct
  for (int i = 0; i < n_gpus; i++) {
    GpuDetails g_deets = gpu_config.getGpuInfo(i);
    if (g_deets.getGpuSupported() == false) {
      continue;
    }
    check(g_deets.getArchMajor(), RelationalOperator::GE, 3, "An unsupported, pre-Kepler device "
          "was approved for use.");
    if (g_deets.getArchMajor() == 7 && g_deets.getArchMinor() == 5) {
      check(g_deets.getMaxThreadsPerSMP() == 1024, "The maximum number of threads per SMP was "
            "detected as " + std::to_string(g_deets.getMaxThreadsPerSMP()) + " on a Turing device "
            "of compute capability 7.5: " + g_deets.getCardName());
    }
    else {
      check(g_deets.getMaxThreadsPerSMP(), RelationalOperator::GE, 1536, "The maximum number of "
            "threads per SMP was detected as " + std::to_string(g_deets.getMaxThreadsPerSMP()) +
            " on a device of compute capability " + std::to_string(g_deets.getArchMajor()) + "." +
            std::to_string(g_deets.getArchMinor()) + ": " + g_deets.getCardName(),
            TestPriority::NON_CRITICAL);
    }
    check(g_deets.getMaxThreadsPerBlock() == 1024, "The maximum number of threads per block was "
          "detected as " + std::to_string(g_deets.getMaxThreadsPerBlock()) + " on a device of "
          "compute capability " + std::to_string(g_deets.getArchMajor()) + "." +
          std::to_string(g_deets.getArchMinor()) + ": " + g_deets.getCardName());
    check(g_deets.getMaxBlocksPerSMP(), RelationalOperator::GE, 16, "The maximum number of blocks "
          "per SMP was detected as " + std::to_string(g_deets.getMaxBlocksPerSMP()) + " on a "
          "device of compute capability " + std::to_string(g_deets.getArchMajor()) + "." +
          std::to_string(g_deets.getArchMinor()) + ": " + g_deets.getCardName());
  }

  // Prepare a positive definite matrix based on random numbers
  section(2);
  Ran2Generator my_prng(oe.getRandomSeed());
  const int m_rank = 64;
  Hybrid<double> matrix_a(m_rank * m_rank, "Matrix A");
  Hybrid<double> matrix_b(m_rank * m_rank, "Matrix B");
  Hybrid<double> matrix_c(m_rank * m_rank, "Matrix C");
  double *ma_ptr = matrix_a.data(HybridTargetLevel::HOST);
  double *mb_ptr = matrix_b.data(HybridTargetLevel::HOST);
  for (int i = 0; i < m_rank * m_rank; i++) {
    ma_ptr[i] = 5.0 * (0.5 - my_prng.uniformRandomNumber());
    mb_ptr[i] = ma_ptr[i];
  }
  matrix_a.upload();
  matrix_b.upload();
  const double multiplier = 1.0;
  cublasDgemm(gpu_config.getCuBlasHandle(), CUBLAS_OP_T, CUBLAS_OP_N, m_rank, m_rank, m_rank,
              &multiplier, matrix_a.data(HybridTargetLevel::DEVICE), m_rank,
              matrix_b.data(HybridTargetLevel::DEVICE), m_rank, &multiplier,
              matrix_c.data(HybridTargetLevel::DEVICE), m_rank);
  cudaDeviceSynchronize();

  // A simple matrix-matrix multiply when one matrix's transpose is implicit
  double *mc_ptr = matrix_c.data(HybridTargetLevel::HOST);
  for (int i = 0; i < m_rank; i++) {
    for (int j = 0; j < m_rank; j++) {
      double mc_sum = 0.0;
      for (int k = 0; k < m_rank; k++) {
        mc_sum += ma_ptr[(i * m_rank) + k] * mb_ptr[(j * m_rank) + k];
      }
      mc_ptr[(i * m_rank) + j] = mc_sum;
    }
  }

  // Compare the results to confirm a successful DGEMM operation
  const std::vector<double> mc_host_image = matrix_c.readHost();
  const std::vector<double> mc_devc_image = matrix_c.readDevice();
  check(mc_host_image, RelationalOperator::EQUAL, Approx(mc_devc_image).margin(1.0e-6), "Matrix "
        "multiplication by cuBLAS and a hand-rolled algorithm yield different results.");

  // Repeat the process with an SGEMM
  const std::vector<double> ma_host_image = matrix_a.readHost();
  const std::vector<float> spma_host_image(ma_host_image.begin(), ma_host_image.end());
  check(spma_host_image, RelationalOperator::NE,
        Approx(ma_host_image, ComparisonType::MEAN_UNSIGNED_ERROR).margin(tiny),
        "Single- and double-precision real data represents numbers to unusually high precision.",
        TestPriority::NON_CRITICAL);
  Hybrid<float> sp_matrix_a(m_rank * m_rank, "fp32 Matrix A", HybridFormat::DEVICE_ONLY);
  Hybrid<float> sp_matrix_b(m_rank * m_rank, "fp32 Matrix B", HybridFormat::DEVICE_ONLY);
  Hybrid<float> sp_matrix_c(m_rank * m_rank, "fp32 Matrix C", HybridFormat::DEVICE_ONLY);
  sp_matrix_a.putDevice(spma_host_image);
  sp_matrix_b.putDevice(spma_host_image);
  const float sp_mult = 1.0;
  cublasSgemm(gpu_config.getCuBlasHandle(), CUBLAS_OP_T, CUBLAS_OP_N, m_rank, m_rank, m_rank,
              &sp_mult, sp_matrix_a.data(HybridTargetLevel::DEVICE), m_rank,
              sp_matrix_b.data(HybridTargetLevel::DEVICE), m_rank, &sp_mult,
              sp_matrix_c.data(HybridTargetLevel::DEVICE), m_rank);
  cudaDeviceSynchronize();
  const std::vector<float> spmc_devc_image = sp_matrix_c.readDevice();
  check(spmc_devc_image, RelationalOperator::EQUAL,
        Approx(mc_devc_image, ComparisonType::MEAN_UNSIGNED_ERROR).margin(1.0e-5),
        "Disagreement between cuBLAS dgemm and sgemm results is larger than expected.");

  // Print a summary of tests run
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
