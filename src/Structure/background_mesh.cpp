#include "copyright.h"
#include "Numerics/host_popc.h"
#include "background_mesh.h"

namespace stormm {
namespace structure {

using numerics::hostPopcll;

//-------------------------------------------------------------------------------------------------
void accumulateOcclusionMesh(BackgroundMeshWriter<llint> *target,
                             const BackgroundMeshReader<llint> &contribution) {
  
  // Check mesh dimensions.  If the contribution and target match, determine the correspondence
  // of mesh elements and sum the bit counts.
  checkMeshCompatibility(*target, contribution);
  int period_factor;
  switch (contribution.dims.bounds) {
  case BoundaryCondition::ISOLATED:
    period_factor = 1;
    break;
  case BoundaryCondition::PERIODIC:
    period_factor = 0;
    break;
  }
  const int samp_a = 4 * contribution.dims.na / (target->dims.na + period_factor);
  const int samp_b = 4 * contribution.dims.nb / (target->dims.nb + period_factor);
  const int samp_c = 4 * contribution.dims.nc / (target->dims.nc + period_factor);
  const size_t cnt_na_zu = contribution.dims.na;
  const size_t cnt_nb_zu = contribution.dims.nb;
  const size_t cnt_nc_zu = contribution.dims.nc;
  const size_t trg_na_zu = target->dims.na;
  const size_t trg_nb_zu = target->dims.nb;
  const size_t trg_nc_zu = target->dims.nc;
  
  // There are 64 x (elements) in the contribution mesh, each of them a 64-bit integer.  Determine
  // which of them are relevant to each accumulator in the target, and accumulate the result.
  const std::vector<int> bit_counts = { 0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4 };
  for (int i = 0; i <= target->dims.na; i++) {
    for (int j = 0; j <= target->dims.nb; j++) {
      for (int k = 0; k <= target->dims.nc; k++) {

        // Loop over the limits of the relevant contribution
        const int min_ca = samp_a * i;
        const int min_cb = samp_b * j;
        const int min_cc = samp_c * k;
        const int max_ca = samp_a * (i + 1);
        const int max_cb = samp_b * (j + 1);
        const int max_cc = samp_c * (k + 1);
        int nbits = 0;
        for (int ci = min_ca; ci < max_ca; ci++) {
          size_t ci_cube = ci / 4;
          const size_t ci_cubelet = ci - (4 * ci_cube);
          if (ci_cube >= contribution.dims.na) {
            ci_cube -= contribution.dims.na;
          }
          for (int cj = min_cb; cj < max_cb; cj++) {
            size_t cj_cube = cj / 4;
            const size_t cj_cubelet = cj - (4 * cj_cube);
            if (cj_cube >= contribution.dims.nb) {
              cj_cube -= contribution.dims.nb;
            }
            for (int ck = min_cc; ck < max_cc; ck++) {

              // Extract the occlusion mesh value and test its bits
              size_t ck_cube = ck / 4;
              const size_t ck_cubelet = ck - (4 * ck_cube);
              if (ck_cube >= contribution.dims.nc) {
                ck_cube -= contribution.dims.nc;
              }
              const size_t ocidx = (64LLU * ((((ck_cube * cnt_nb_zu) + cj_cube) *
                                              cnt_na_zu) + ci_cube)) +
                                   static_cast<size_t>((((ck_cubelet * 4) + cj_cubelet) * 4) +
                                                       ci_cubelet);
              const ullint ocval = contribution.coeffs[ocidx];
              if (ocval > 0LLU) {
                nbits += hostPopcll(ocval);
              }
            }
          }
        }

        // Replicate the potential across up to eight mesh elements
        for (int ci = 0; ci <= 1; ci++) {
          int im = i - ci;
          if (im < 0) {
            switch (contribution.dims.bounds) {
            case BoundaryCondition::ISOLATED:
              continue;
            case BoundaryCondition::PERIODIC:
              im += target->dims.na;
              break;
            }
          }
          for (int cj = 0; cj <= 1; cj++) {
            int jm = j - cj;
            if (jm < 0) {
              switch (contribution.dims.bounds) {
              case BoundaryCondition::ISOLATED:
                continue;
              case BoundaryCondition::PERIODIC:
                jm += target->dims.nb;
                break;
              }
            }
            for (int ck = 0; ck <= 1; ck++) {
              int km = k - ck;
              if (km < 0) {
                switch (contribution.dims.bounds) {
                case BoundaryCondition::ISOLATED:
                  continue;
                case BoundaryCondition::PERIODIC:
                  km += target->dims.nc;
                  break;
                }
              }

              // The first eight slots of the coefficients array for this mesh element will hold
              // the potential at each of the element's corners.  When accumulating meshes on the
              // GPU, the d/dx, d/dy, and d/dz first derivatives will fall into slots [8, 16),
              // [16, 24), and [24, 32).  The d2/dxy, d2/dxz, and d2/dyz mixed particle derivatives
              // will be accumulated in slots [32, 40), [40, 48), and [48, 56).  The d3/dxyz
              // partial derivatives at each of the element's corners will be accumulated in slots
              // [56, 64).
              const size_t target_idx = (((km * target->dims.nb) + jm) * target->dims.na) + im;
              const size_t target_subidx = (((ck * 2) + cj) * 2) + ci;
              target->coeffs[(64LLU * target_idx) + target_subidx] = nbits;
            }
          }
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void accumulateOcclusionMesh(BackgroundMesh<llint> *target,
                             const BackgroundMesh<llint> &contribution, const GpuDetails &gpu,
                             const HybridTargetLevel availability) {
  BackgroundMeshWriter<llint> targetw = target->data(availability);
  const BackgroundMeshReader<llint> contribr = contribution.data(availability);
  if (gpu != null_gpu) {
    switch (availability) {
    case HybridTargetLevel::HOST:
      accumulateOcclusionMesh(&targetw, contribr);
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      launchAccOcclusionMesh(&targetw, contribr, gpu);
      break;
#endif
    }
    return;
  }
  else {
    switch (availability) {
    case HybridTargetLevel::HOST:
      accumulateOcclusionMesh(&targetw, contribr);
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      rtErr("Unable to perform mesh accumulation on the GPU device without providing valid GPU "
            "specifications.", "accumulateOcclusionMesh");
#endif
    }
  }
}

} // namespace structure
} // namespace stormm
