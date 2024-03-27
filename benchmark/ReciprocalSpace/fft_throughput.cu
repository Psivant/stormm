// -*-c++-*-
#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <limits.h>
#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Accelerator/gpu_details.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/Constants/behavior.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/stopwatch.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::card;
using namespace stormm::constants;
using namespace stormm::data_types;
using namespace stormm::parse;
using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// Multiply a block of data by a constant.
//
// Arguments:
//   data:    Pointer to the array to scale
//   n:       Trusted length of data
//   mult:    Factor by which to multiply all elements of the array
//-------------------------------------------------------------------------------------------------
template <typename T>
__global__ void __launch_bounds__(large_block_size, 1)
kMult(T* data, const size_t n, const T mult) {
  size_t pos = threadIdx.x + (blockIdx.x * blockDim.x);
  while (pos < n) {
    data[pos] *= mult;
    pos += (blockDim.x * gridDim.x);
  }
}

//-------------------------------------------------------------------------------------------------
// Lay out an FFT of a particular size and perform it for a given number of repeats using in-place
// transformations.  Time the result.
//
// Arguments:
//   nx:      Dimension of the FFT along the unit cell A axis
//   ny:      Dimension of the FFT along the unit cell B axis
//   nz:      Dimension of the FFT along the unit cell C axis
//   nbatch:  The number of FFTs to combine into a batch
//   prec:    Precision in which to perform calculations
//   use_ip:  Set to TRUE to apply in-place transforms, FALSE for out-of-place transforms
//   chkdir:  Set to TRUE to have forward and backward FFTs independently timed, if the problem
//            size and number of iterations is reasonable
//   xrs:     Random number generator to use in preparing data
//   iter:    The number of times to repeat the FFT / inverse FFT cycle
//   timer:   Tracks the wall time
//   gpu:     Details of the available GPU
//-------------------------------------------------------------------------------------------------
template <typename T>
void runFFT(const int nx, const int ny, const int nz, const int nbatch, const bool use_ip,
            const bool chkdir, Xoshiro256ppGenerator *xrs, const int iter, StopWatch *timer,
            const GpuDetails &gpu) {
  const PrecisionModel prec = (std::type_index(typeid(T)).hash_code() == double_type_index) ?
                              PrecisionModel::DOUBLE : PrecisionModel::SINGLE;
  const std::string place_str = (use_ip) ? "IP" : "OOP";
  const int t_id = timer->addCategory(std::string("FFT(") + intToString(nx, 3) + ", " +
                                      intToString(ny, 3) + ", " + intToString(nz, 3) + "), " +
                                      getEnumerationName(prec) + ", " + place_str);
  const int m_id = timer->addCategory(std::string("Scl(") + intToString(nx, 3) + ", " +
                                      intToString(ny, 3) + ", " + intToString(nz, 3) + "), " +
                                      getEnumerationName(prec) + ", " + place_str);
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  const int nxyz     = nx * ny * nz;
  const int nz_pad   = 2 * ((nz / 2) + 1);
  const int nz_eff   = (use_ip) ? nz_pad : nz;
  const int nxyz_eff = nx * ny * nz_eff;
  std::vector<double> trial_load = gaussianRand(xrs, nxyz * nbatch, 1.0);
  cufftHandle frwd_plan, bkwd_plan;
  bool frwd_problem = false;
  bool bkwd_problem = false;

  // The following variables are computed for batch FFT setup.
  int dims[] = { nx, ny, nz };
  int real_embed[3];
  real_embed[0] = nx;
  real_embed[1] = ny;
  real_embed[2] = nz_eff;
  int cmpx_embed[3];
  cmpx_embed[0] = nx;
  cmpx_embed[1] = ny;
  cmpx_embed[2] = (use_ip) ? nz_pad / 2 : nz;
  const int real_length = real_embed[0] * real_embed[1] * real_embed[2];
  const int cmpx_length = cmpx_embed[0] * cmpx_embed[1] * cmpx_embed[2];
  const std::string batch_msg = (nbatch == 1) ? ")" : "), batch " + std::to_string(nbatch);
  const double nxyz_inv_scale = 1.0 / static_cast<double>(nxyz);

  // Create the FFT plans
  cufftType frwd_kind, bkwd_kind;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    frwd_kind = CUFFT_D2Z;
    bkwd_kind = CUFFT_Z2D;
    break;
  case PrecisionModel::SINGLE:
    frwd_kind = CUFFT_R2C;
    bkwd_kind = CUFFT_C2R;
    break;
  }
  if (nbatch == 1) {
    if (cufftPlan3d(&frwd_plan, nx, ny, nz, frwd_kind) != CUFFT_SUCCESS) {
      frwd_problem = true;
    }
    if (cufftPlan3d(&bkwd_plan, nx, ny, nz, bkwd_kind) != CUFFT_SUCCESS) {
      bkwd_problem = true;
    }
  }
  else {
    if (cufftPlanMany(&frwd_plan, 3, dims, real_embed, 1, real_length, cmpx_embed, 1,
                      cmpx_length, frwd_kind, nbatch) != CUFFT_SUCCESS) {
      frwd_problem = true;
    }
    if (cufftPlanMany(&bkwd_plan, 3, dims, cmpx_embed, 1, cmpx_length, real_embed, 1,
                      real_length, bkwd_kind, nbatch) != CUFFT_SUCCESS) {
      bkwd_problem = true;
    }
  }
  if (frwd_problem || bkwd_problem) {
    const std::string pprob = (frwd_problem) ? "forward" : "inverse";
    rtErr("Failed to create " + pprob + " plan for " + std::to_string(nx) + " x " +
          std::to_string(ny) + " x " + std::to_string(nz) + " points (" +
          getEnumerationName(prec) + batch_msg + ".", "runFFT");
  }
  const int nblocks  = gpu.getSMPCount();
  const int nthreads = large_block_size;

  // Allocate and fill a data array, upload, and perform the FFT cycle 100x
  Hybrid<T> trial(nbatch * nxyz_eff, "trial_fft");
  Hybrid<T> trial_t(HybridKind::ARRAY, "trial_fft_t");
  T* trial_ptr = trial.data();
  for (int n = 0; n < nbatch; n++) {
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        const int nij_t = (n * nxyz_eff) + (((i * ny) + j) * nz_eff);
        const int nij   = (n * nxyz)     + (((i * ny) + j) * nz);
        for (int k = 0; k < nz; k++) {
          trial_ptr[nij_t + k] = trial_load[nij + k];
        }
      }
    }
  }
  trial.upload();
  timer->assignTime(0);

  // Time the scaling kernel--this must be applied to ensure that the data does not become full of
  // NaN or Inf values after many iterations.
  for (int i = 0; i < iter; i++) {
    kMult<T><<<nblocks, nthreads>>>(trial.data(devc_tier), nbatch * real_length, 1.0);
  }
  if (cudaDeviceSynchronize() != cudaSuccess) {
    rtErr("CUDA device synchronize failed.", "runIpFFT");
  }
  timer->assignTime(m_id);

  // Set pointers to the data.  Pairs of pointers for both single- and double-precision data will
  // be created, but only one pair will be valid and used.
  cufftDoubleReal* real_data = reinterpret_cast<cufftDoubleReal*>(trial.data(devc_tier));
  cufftReal* sp_real_data    = reinterpret_cast<cufftReal*>(trial.data(devc_tier));
  cufftDoubleComplex* cmpx_data;
  cufftComplex* sp_cmpx_data;
  if (use_ip) {
    cmpx_data    = reinterpret_cast<cufftDoubleComplex*>(trial.data(devc_tier));        
    sp_cmpx_data = reinterpret_cast<cufftComplex*>(trial.data(devc_tier));        
  }
  else {
    trial_t.resize(2 * nbatch * nx * ny * nz);
    cmpx_data = reinterpret_cast<cufftDoubleComplex*>(trial_t.data(devc_tier));
    sp_cmpx_data = reinterpret_cast<cufftComplex*>(trial_t.data(devc_tier));        
  }

  // Run the FFT for the required number of iterations.  Rescale the result if needed.
  timer->assignTime(0);
  for (int i = 0; i < iter; i++) {
    switch (prec) {
    case PrecisionModel::DOUBLE:
      if (cufftExecD2Z(frwd_plan, real_data, cmpx_data) != CUFFT_SUCCESS) {
        frwd_problem = true;
        break;
      }
      if (cufftExecZ2D(bkwd_plan, cmpx_data, real_data) != CUFFT_SUCCESS) {
        bkwd_problem = true;
        break;
      }
      break;
    case PrecisionModel::SINGLE:
      if (cufftExecR2C(frwd_plan, sp_real_data, sp_cmpx_data) != CUFFT_SUCCESS) {
        frwd_problem = true;
        break;
      }
      if (cufftExecC2R(bkwd_plan, sp_cmpx_data, sp_real_data) != CUFFT_SUCCESS) {
        bkwd_problem = true;
        break;
      }
      break;
    }
    kMult<T><<<nblocks, nthreads>>>(trial.data(devc_tier), nbatch * real_length, nxyz_inv_scale);
  }
  
  // Synchronize the CPU to get timings on the FFT
  if (cudaDeviceSynchronize() != cudaSuccess) {
    rtErr("CUDA device synchronize failed.", "runIpFFT");
  }
  timer->assignTime(t_id);
  if (frwd_problem || bkwd_problem) {
    const std::string pprob = (frwd_problem) ? "forward" : "inverse";
    rtErr("Failed to execute " + pprob + " FFT for " + std::to_string(nx) + " x " +
          std::to_string(ny) + " x " + std::to_string(nz) + " points (" +
          getEnumerationName(prec) + batch_msg + ".", "runIpFFT");
  }

  // Check the results, if it is reasonable to do so.  Many iterations of the FFT may create some
  // degree of wander in the output, but working with just a few iterations should make it possible
  // to recover the original data to within some fairly tight bounds.
  trial.download();
  std::vector<double> trial_outcome(nbatch * nxyz);
  for (int n = 0; n < nbatch; n++) {
    for (int i = 0; i < nx; i++) {
      for (int j = 0; j < ny; j++) {
        const int nij_t = (n * nxyz_eff) + (((i * ny) + j) * nz_eff);
        const int nij   = (n * nxyz)     + (((i * ny) + j) * nz);
        for (int k = 0; k < nz; k++) {
          trial_outcome[nij + k] = trial_ptr[nij_t + k];
        }
      }
    }
  }
  if (iter < 5) {
    const double tol = static_cast<double>(iter) * 1.0e-5;
    check(trial_outcome, RelationalOperator::EQUAL, Approx(trial_load).margin(tol), "The result "
          "of forward and inverse FFTs on " + std::to_string(nx) + " x " + std::to_string(ny) +
          " x " + std::to_string(nz) + " points did not recover the original data.  Precision "
          "model: " + getEnumerationName(prec) + ".  Batch count: " + std::to_string(nbatch) +
          ".");
  }

  // Check the forward and backward FFT timings independently
  if (nx == ny && nx == nz && nbatch == 1 && iter < (INT_MAX / (nx * ny * nz)) - 1 &&
      iter < 100 && chkdir) {
    Hybrid<T> x_trial(iter * nx * ny * nz_eff, "replica_trial");
    Hybrid<T> x_trial_t(HybridKind::ARRAY, "replica_trial_t");
    cufftDoubleReal* x_real_data = reinterpret_cast<cufftDoubleReal*>(x_trial.data(devc_tier));
    cufftReal* sp_x_real_data    = reinterpret_cast<cufftReal*>(x_trial.data(devc_tier));
    cufftDoubleComplex* x_cmpx_data;
    cufftComplex* sp_x_cmpx_data;
    if (use_ip) {
      x_cmpx_data    = reinterpret_cast<cufftDoubleComplex*>(x_trial.data(devc_tier));        
      sp_x_cmpx_data = reinterpret_cast<cufftComplex*>(x_trial.data(devc_tier));        
    }
    else {
      x_trial_t.resize(2 * iter * nx * ny * nz);
      x_cmpx_data = reinterpret_cast<cufftDoubleComplex*>(x_trial_t.data(devc_tier));
      sp_x_cmpx_data = reinterpret_cast<cufftComplex*>(x_trial_t.data(devc_tier));        
    }
    T* x_trial_ptr = x_trial.data();
    for (int n = 0; n < iter; n++) {
      for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
          const int nij_t = (n * nxyz_eff) + (((i * ny) + j) * nz_eff);
          for (int k = 0; k < nz; k++) {
            x_trial_ptr[nij_t + k] = xrs->gaussianRandomNumber();
          }
        }
      }
    }
    x_trial.upload();
    const int frwd_id = timer->addCategory(std::string("Fwd(") + intToString(nx, 3) + ", " +
                                           intToString(ny, 3) + ", " + intToString(nz, 3) + "), " +
                                           getEnumerationName(prec) + ", " + place_str);
    const int bkwd_id = timer->addCategory(std::string("Inv(") + intToString(nx, 3) + ", " +
                                           intToString(ny, 3) + ", " + intToString(nz, 3) + "), " +
                                           getEnumerationName(prec) + ", " + place_str);
    timer->assignTime(0);
    for (int i = 0; i < iter; i++) {
      const size_t real_idx = i * nxyz_eff;
      const size_t cmpx_idx = (use_ip) ? i * nxyz_eff / 2 : i * nxyz_eff;
      switch (prec) {
      case PrecisionModel::DOUBLE:
        if (cufftExecD2Z(frwd_plan,
                         &x_real_data[real_idx], &x_cmpx_data[cmpx_idx]) != CUFFT_SUCCESS) {
          frwd_problem = true;
          break;
        }
        break;
      case PrecisionModel::SINGLE:
        if (cufftExecR2C(frwd_plan,
                         &sp_x_real_data[real_idx], &sp_x_cmpx_data[cmpx_idx]) != CUFFT_SUCCESS) {
          frwd_problem = true;
          break;
        }
        break;
      }
    }
    if (cudaDeviceSynchronize() != cudaSuccess) {
      rtErr("CUDA device synchronize failed.", "runIpFFT");
    }
    timer->assignTime(frwd_id);
    for (int i = 0; i < iter; i++) {
      const size_t real_idx = i * nxyz_eff;
      const size_t cmpx_idx = (use_ip) ? i * nxyz_eff / 2 : i * nxyz_eff;
      switch (prec) {
      case PrecisionModel::DOUBLE:
        if (cufftExecZ2D(bkwd_plan,
                         &x_cmpx_data[cmpx_idx], &x_real_data[real_idx]) != CUFFT_SUCCESS) {
          bkwd_problem = true;
          break;
        }
        break;
      case PrecisionModel::SINGLE:
        if (cufftExecC2R(bkwd_plan,
                         &sp_x_cmpx_data[cmpx_idx], &sp_x_real_data[real_idx]) != CUFFT_SUCCESS) {
          bkwd_problem = true;
          break;
        }
        break;
      }
    }
    if (cudaDeviceSynchronize() != cudaSuccess) {
      rtErr("CUDA device synchronize failed.", "runIpFFT");
    }
    timer->assignTime(bkwd_id);
  }
  
  // Destroy the plans to prevent having too many active handles
  if (cufftDestroy(frwd_plan) != CUFFT_SUCCESS) {
    frwd_problem = true;
  }
  if (cufftDestroy(bkwd_plan) != CUFFT_SUCCESS) {
    bkwd_problem = true;
  }
  if (frwd_problem || bkwd_problem) {
    const std::string pprob = (frwd_problem) ? "forward" : "inverse";
    rtErr("Failed to destroy the " + pprob + " FFT plan for " + std::to_string(nx) + " x " +
          std::to_string(ny) + " x " + std::to_string(nz) + " points (" +
          getEnumerationName(prec) + batch_msg + ".", "runIpFFT");
  }
}

//-------------------------------------------------------------------------------------------------
// Perform a simple, out-of-place, round-trip FFT for an 8 x 8 x 8 grid and check the result.
//
// Arguments:
//   xrs:     Random number generator to produce initial data
//-------------------------------------------------------------------------------------------------
void simpleOopFFT(Xoshiro256ppGenerator *xrs) {
  Hybrid<double> trial(512, "trial_fft"), trial_t(1024, "trial_fft_t");
  const std::vector<double> trial_fill = gaussianRand(xrs, trial.size(), 1.0);
  trial.putHost(trial_fill);
  trial.upload();
  cufftHandle trial_plan, trial_plan_inv;
  if (cufftPlan3d(&trial_plan, 8, 8, 8, CUFFT_D2Z) != CUFFT_SUCCESS) {
    rtErr("Unable to create forward plan for the trial.", "main");
  }
  if (cufftPlan3d(&trial_plan_inv, 8, 8, 8, CUFFT_Z2D) != CUFFT_SUCCESS) {
    rtErr("Unable to create inverse plan for the trial.", "main");
  }
  cufftExecD2Z(trial_plan,
               reinterpret_cast<cufftDoubleReal*>(trial.data(HybridTargetLevel::DEVICE)),
               reinterpret_cast<cufftDoubleComplex*>(trial_t.data(HybridTargetLevel::DEVICE)));
  cufftExecZ2D(trial_plan_inv,
               reinterpret_cast<cufftDoubleComplex*>(trial_t.data(HybridTargetLevel::DEVICE)),
               reinterpret_cast<cufftDoubleReal*>(trial.data(HybridTargetLevel::DEVICE)));
  cufftDestroy(trial_plan);
  cufftDestroy(trial_plan_inv);
  std::vector<double> trial_result = trial.readDevice();
  for (int i = 0; i < 512; i++) {
    trial_result[i] /= 512.0;
  }
  check(trial_result, RelationalOperator::EQUAL, trial_fill, "The result of forward and inverse "
        "out-of-place FFTs does not match the original input.");
}

//-------------------------------------------------------------------------------------------------
// Perform a simple, in-place, round-trip FFT for an 8 x 8 x 8 grid and check the result.
//
// Arguments:
//   xrs:     Random number generator to produce initial data
//-------------------------------------------------------------------------------------------------
void simpleIpFFT(Xoshiro256ppGenerator *xrs) {
  const int nx = 8;
  const int ny = 8;
  const int nz = 8;
  const int nz_pad = 2 * ((nz / 2) + 1);
  Hybrid<double> trial(nx * ny * nz_pad, "trial_fft");
  const std::vector<double> trial_fill = gaussianRand(xrs, nx * ny * nz, 1.0);
  double* trial_ptr = trial.data();
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        trial_ptr[(((i * ny) + j) * nz_pad) + k] = trial_fill[(((i * ny) + j) * nz) + k];
      }
    }
  }
  trial.upload();
  cufftHandle trial_plan, trial_plan_inv;
  if (cufftPlan3d(&trial_plan, 8, 8, 8, CUFFT_D2Z) != CUFFT_SUCCESS) {
    rtErr("Unable to create forward plan for the trial.", "main");
  }
  if (cufftPlan3d(&trial_plan_inv, 8, 8, 8, CUFFT_Z2D) != CUFFT_SUCCESS) {
    rtErr("Unable to create inverse plan for the trial.", "main");
  }
  cufftExecD2Z(trial_plan,
               reinterpret_cast<cufftDoubleReal*>(trial.data(HybridTargetLevel::DEVICE)),
               reinterpret_cast<cufftDoubleComplex*>(trial.data(HybridTargetLevel::DEVICE)));
  cufftExecZ2D(trial_plan_inv,
               reinterpret_cast<cufftDoubleComplex*>(trial.data(HybridTargetLevel::DEVICE)),
               reinterpret_cast<cufftDoubleReal*>(trial.data(HybridTargetLevel::DEVICE)));
  cufftDestroy(trial_plan);
  cufftDestroy(trial_plan_inv);
  const std::vector<double> trial_result_raw = trial.readDevice();
  std::vector<double> trial_result(nx * ny * nz);
  for (int i = 0; i < nx; i++) {
    for (int j = 0; j < ny; j++) {
      for (int k = 0; k < nz; k++) {
        trial_result[(((i * ny) + j) * nz) + k] = trial_result_raw[(((i * ny) + j) * nz_pad) + k];
      }
    }
  }
  for (int i = 0; i < 512; i++) {
    trial_result[i] /= 512.0;
  }
  check(trial_result, RelationalOperator::EQUAL, trial_fill, "The result of forward and inverse "
        "in-place FFTs does not match the original input.");
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Baseline variables
  TestEnvironment oe(argc, argv, ExceptionResponse::SILENT);

  // Take in additional command-line variables
  bool use_ip = false;
  int nbatch = 1;
  int iter = 100;
  int max_grid = 256;
  int max_pts = 1024 * 1024 * 1024;
  int min_grid = 8;
  int min_pts = 512;
  bool test_double = false;
  bool test_single = false;
  bool check_directional_ffts = false;
  std::vector<int> radices, non_radices;
  for (int i = 0; i < argc; i++) {
    if (strcmpCased(argv[i], "-ip", CaseSensitivity::NO)) {
      use_ip = true;
    }
    else if (strcmpCased(argv[i], "-double", CaseSensitivity::NO)) {
      test_double = true;
    }
    else if (strcmpCased(argv[i], "-single", CaseSensitivity::NO)) {
      test_single = true;
    }
    else if (strcmpCased(argv[i], "-directional", CaseSensitivity::NO)) {
      check_directional_ffts = true;
    }
    else if (i < argc - 1 && strcmpCased(argv[i], "-batch", CaseSensitivity::NO)) {
      if (verifyNumberFormat(argv[i + 1], NumberFormat::INTEGER)) {
        nbatch = atoi(argv[i + 1]);
      }
      else {
        rtErr("Invalid batch count \"" + std::string(argv[i + 1]) + "\".", "main");
      }
    }
    else if (i < argc - 1 && strcmpCased(argv[i], "-iter", CaseSensitivity::NO)) {
      if (verifyNumberFormat(argv[i + 1], NumberFormat::INTEGER)) {
        iter = atoi(argv[i + 1]);
      }
      else {
        rtErr("Invalid iteration count \"" + std::string(argv[i + 1]) + "\".", "main");
      }
    }
    else if (i < argc - 1 && strcmpCased(argv[i], "-max_grid", CaseSensitivity::NO)) {
      if (verifyNumberFormat(argv[i + 1], NumberFormat::INTEGER)) {
        max_grid = atoi(argv[i + 1]);
      }
      else {
        rtErr("Invalid maximum grid dimension \"" + std::string(argv[i + 1]) + "\".", "main");
      }
    }
    else if (i < argc - 1 && strcmpCased(argv[i], "-min_grid", CaseSensitivity::NO)) {
      if (verifyNumberFormat(argv[i + 1], NumberFormat::INTEGER)) {
        min_grid = atoi(argv[i + 1]);
      }
      else {
        rtErr("Invalid minimum grid dimension \"" + std::string(argv[i + 1]) + "\".", "main");
      }
    }
    else if (i < argc - 1 && strcmpCased(argv[i], "-max_pts", CaseSensitivity::NO)) {
      if (verifyNumberFormat(argv[i + 1], NumberFormat::INTEGER)) {
        max_pts = atoi(argv[i + 1]);
      }
      else {
        rtErr("Invalid maximum point count \"" + std::string(argv[i + 1]) + "\".", "main");
      }
    }
    else if (i < argc - 1 && strcmpCased(argv[i], "-min_pts", CaseSensitivity::NO)) {
      if (verifyNumberFormat(argv[i + 1], NumberFormat::INTEGER)) {
        min_pts = atoi(argv[i + 1]);
      }
      else {
        rtErr("Invalid minimum point count \"" + std::string(argv[i + 1]) + "\".", "main");
      }
    }
    else if (i < argc - 1 && strcmpCased(argv[i], "-radix", CaseSensitivity::NO)) {
      if (verifyNumberFormat(argv[i + 1], NumberFormat::INTEGER)) {
        radices.push_back(atoi(argv[i + 1]));
      }
      else {
        rtErr("Invalid radix \"" + std::string(argv[i + 1]) + "\".", "main");
      }
    }
    else if (i < argc - 1 && strcmpCased(argv[i], "-noradix", CaseSensitivity::NO)) {
      if (verifyNumberFormat(argv[i + 1], NumberFormat::INTEGER)) {
        non_radices.push_back(atoi(argv[i + 1]));
      }
      else {
        rtErr("Invalid excluded radix \"" + std::string(argv[i + 1]) + "\".", "main");
      }
    }
  }
  if (test_double == false && test_single == false) {
    test_double = true;
    test_single = true;
  }
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
  StopWatch timer;
#ifdef STORMM_USE_CUDA
  const HpcConfig gpu_config(ExceptionResponse::WARN);
  const std::vector<int> my_gpus = gpu_config.getGpuDevice(1);
  const GpuDetails gpu = gpu_config.getGpuInfo(my_gpus[0]);
  Hybrid<int> force_gpu_to_engage(1);
#endif

  // Initialize a random number generator and run some trials
  Xoshiro256ppGenerator xrs;
  simpleOopFFT(&xrs);
  simpleIpFFT(&xrs);
  const std::vector<int> ffts = { 24, 28, 30, 32, 33, 36, 40, 44, 45, 48, 54, 55, 56, 60, 64, 66,
                                  72, 75, 77, 80, 84, 88, 90, 96, 99, 100, 108, 110, 112, 120, 128,
                                  132, 135, 140, 144, 150, 154, 160, 165, 168, 176, 180, 192, 198,
                                  200, 210, 220, 224, 231, 240, 250, 256 };
  std::vector<PrecisionModel> all_prec;
  if (test_double) {
    all_prec.push_back(PrecisionModel::DOUBLE);
  }
  if (test_single) {
    all_prec.push_back(PrecisionModel::SINGLE);
  }
  int radix_product = 1;
  if (radices.size() > 0) {
    for (size_t i = 0; i < radices.size(); i++) {
      radix_product *= radices[i];
    }
  }
  for (size_t p = 0; p < all_prec.size(); p++) {
    for (size_t i = 0; i < ffts.size(); i++) {
      for (size_t j = i; j < ffts.size(); j++) {
        for (size_t k = j; k < ffts.size(); k++) {

          // Check that the aspect ratios are realistic for a typical MD simulation
          const int min_ij = std::min(ffts[i], ffts[j]);
          const int min_ik = std::min(ffts[i], ffts[k]);
          const int min_jk = std::min(ffts[j], ffts[k]);
          if (ffts[k] > max_grid || ffts[i] < min_grid ||
              (std::max(ffts[i], ffts[j]) + (min_ij - 1)) / min_ij > 2 ||
              (std::max(ffts[i], ffts[k]) + (min_ik - 1)) / min_ik > 2 ||
              (std::max(ffts[j], ffts[k]) + (min_jk - 1)) / min_jk > 2) {
            continue;
          }
          const int total_pts = ffts[i] * ffts[j] * ffts[k];

          // Skip if the total point count is outside the range of interest
          if (total_pts * nbatch < min_pts || total_pts * nbatch > max_pts) {
            continue;
          }

          // Skip is the requested radices are not present
          if (radix_product > 1 && (total_pts % radix_product) > 0) {
            continue;
          }

          // Skip if a forbidden radix is present, after accounting for requested radices
          const int rquot = total_pts / radix_product;
          bool has_bad_radix = false;
          for (size_t ir = 0; ir < non_radices.size(); ir++) {
            has_bad_radix = (has_bad_radix || (rquot % non_radices[ir]) == 0);
          }
          if (has_bad_radix) {
            continue;
          }
          
          // Perform the FFTs in the requested mode
          switch (all_prec[p]) {
          case PrecisionModel::DOUBLE:
            runFFT<double>(ffts[i], ffts[j], ffts[k], nbatch, use_ip, check_directional_ffts,
                           &xrs, iter, &timer, gpu);
            break;
          case PrecisionModel::SINGLE:
            runFFT<float>(ffts[i], ffts[j], ffts[k], nbatch, use_ip, check_directional_ffts,
                          &xrs, iter, &timer, gpu);
            break;
          }
        }
      }
    }
  }

  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  return 0;
}
