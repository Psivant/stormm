// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
template <typename T>
double gatherNormalization(const GenericRdSubstrate<T> rsbs, const int start_pos,
                           const int end_pos) {

  // Scale down the values of vector elements if fixed precision is in use (if it is not, the
  // value of the scaling factor will just one 1.0).  While the double-precision format could
  // likely handle the squared values of large forces scaled up by, say, 2^72, that is 2^144 extra
  // stress on the number format that is not needed.  If this were ever to go to single-precision,
  // it would be impossible.
  double tsum = 0.0;
  if (rsbs.y_read == nullptr && rsbs.z_read == nullptr) {
    for (int j = start_pos; j < end_pos; j++) {
      const double dx = static_cast<double>(rsbs.x_read[j]) * rsbs.inv_fp_scaling;
      tsum += (dx * dx);
    }
  }
  else if (rsbs.z_read == nullptr) {
    for (int j = start_pos; j < end_pos; j++) {
      const double dx = static_cast<double>(rsbs.x_read[j]) * rsbs.inv_fp_scaling;
      const double dy = static_cast<double>(rsbs.y_read[j]) * rsbs.inv_fp_scaling;
      tsum += (dx * dx) + (dy * dy);
    }
  }
  else {
    for (int j = start_pos; j < end_pos; j++) {
      const double dx = static_cast<double>(rsbs.x_read[j]) * rsbs.inv_fp_scaling;
      const double dy = static_cast<double>(rsbs.y_read[j]) * rsbs.inv_fp_scaling;
      const double dz = static_cast<double>(rsbs.z_read[j]) * rsbs.inv_fp_scaling;
      tsum += (dx * dx) + (dy * dy) + (dz * dz);
    }
  }
  return tsum;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
double3 gatherCenterOnZero(const GenericRdSubstrate<T> rsbs, const int start_pos,
                           const int end_pos) {
  double tsum_x = 0.0;
  double tsum_y = 0.0;
  double tsum_z = 0.0;
  const int nval = end_pos - start_pos;
  tsum_x = sum<double>(&rsbs.x_read[start_pos], nval);
  if (rsbs.y_read != nullptr) {
    tsum_y = sum<double>(&rsbs.y_read[start_pos], nval);
  }
  if (rsbs.z_read != nullptr) {
    tsum_z = sum<double>(&rsbs.z_read[start_pos], nval);
  }
  return { -tsum_x, -tsum_y, -tsum_z };
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void scatterNormalization(GenericRdSubstrate<T> rsbs, const double tsum, const int start_pos,
                          const int end_pos) {

  // When normalizing extended fixed-precision values, it is not necessary to scale them down with
  // whatever fixed-precision scaling factor only to scale them right back up.  Just divide by the
  // magnitude of the vector found earlier.
  const double nfactor = 1.0 / sqrt(tsum);
  if (rsbs.y_read == nullptr && rsbs.z_read == nullptr) {
    for (int j = start_pos; j < end_pos; j++) {
      const double dx = static_cast<double>(rsbs.x_write[j]);
      rsbs.x_write[j] = static_cast<T>(dx * nfactor);
    }
  }
  else if (rsbs.z_read == nullptr) {
    for (int j = start_pos; j < end_pos; j++) {
      const double dx = static_cast<double>(rsbs.x_write[j]);
      const double dy = static_cast<double>(rsbs.y_write[j]);
      rsbs.x_write[j] = static_cast<T>(dx * nfactor);
      rsbs.y_write[j] = static_cast<T>(dy * nfactor);
    }
  }
  else {
    for (int j = start_pos; j < end_pos; j++) {
      const double dx = static_cast<double>(rsbs.x_write[j]);
      const double dy = static_cast<double>(rsbs.y_write[j]);
      const double dz = static_cast<double>(rsbs.z_write[j]);
      rsbs.x_write[j] = static_cast<T>(dx * nfactor);
      rsbs.y_write[j] = static_cast<T>(dy * nfactor);
      rsbs.z_write[j] = static_cast<T>(dz * nfactor);
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void scatterCenterOnZero(GenericRdSubstrate<T> rsbs, const double tsum_x, const double tsum_y,
                         const double tsum_z, const int natom, const int start_pos,
                         const int end_pos) {
  
  // Values entering the calculation of the sum were never scaled down from their fixed-precision
  // values, so there is no need to rescale any fixed-precision representations here.
  const double inv_norm = 1.0 / static_cast<double>(natom);
  const T center_x = static_cast<T>(tsum_x * inv_norm);
  addScalarToVector(&rsbs.x_write[start_pos], end_pos - start_pos, center_x);
  if (rsbs.y_read != nullptr) {
    const T center_y = static_cast<T>(tsum_y * inv_norm);
    addScalarToVector(&rsbs.y_write[start_pos], end_pos - start_pos, center_y);
  }
  if (rsbs.z_read != nullptr) {
    const T center_z = static_cast<T>(tsum_z * inv_norm);
    addScalarToVector(&rsbs.z_write[start_pos], end_pos - start_pos, center_z);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void evalReduction(GenericRdSubstrate<T> *rsbs, const ReductionKit &redk,
                   const ReductionStage process, const ReductionGoal purpose) {

  // Make the presence of an overflow array for the first data dimension a bellwether for the
  // implementation of extended precision throughout.
  const int atom_start_id = static_cast<int>(RdwuAbstractMap::ATOM_START);
  const int atom_end_id = static_cast<int>(RdwuAbstractMap::ATOM_END);
  const int result_id = static_cast<int>(RdwuAbstractMap::RESULT_INDEX);
  const int dep_start_id = static_cast<int>(RdwuAbstractMap::DEPN_START);
  const int dep_end_id = static_cast<int>(RdwuAbstractMap::DEPN_END);
  const int system_id = static_cast<int>(RdwuAbstractMap::SYSTEM_ID);
  switch (process) {
  case ReductionStage::GATHER:
    for (int i = 0; i < redk.nrdwu; i++) {
      const int start_pos  = redk.rdwu_abstracts[(i * rdwu_abstract_length) + atom_start_id];
      const int end_pos    = redk.rdwu_abstracts[(i * rdwu_abstract_length) + atom_end_id];
      const int result_pos = redk.rdwu_abstracts[(i * rdwu_abstract_length) + result_id];
      const int nval       = end_pos - start_pos;
      switch (purpose) {
      case ReductionGoal::NORMALIZE:
        rsbs->x_buffer[result_pos] = gatherNormalization(*rsbs, start_pos, end_pos);
        break;
      case ReductionGoal::CENTER_ON_ZERO:
        {
          const double3 tsum3 = gatherCenterOnZero(*rsbs, start_pos, end_pos);
          rsbs->x_buffer[result_pos] = tsum3.x;
          if (rsbs->y_read != nullptr) {
            rsbs->y_buffer[result_pos] = tsum3.y;
          }
          if (rsbs->z_read != nullptr) {
            rsbs->z_buffer[result_pos] = tsum3.z;
          }
        }
        break;
      case ReductionGoal::CONJUGATE_GRADIENT:
        break;
      }
    }
    break;
  case ReductionStage::SCATTER:
    for (int i = 0; i < redk.nrdwu; i++) {
      const int start_pos      = redk.rdwu_abstracts[(i * rdwu_abstract_length) + atom_start_id];
      const int end_pos        = redk.rdwu_abstracts[(i * rdwu_abstract_length) + atom_end_id];
      const int result_pos     = redk.rdwu_abstracts[(i * rdwu_abstract_length) + result_id];
      const int depn_start_pos = redk.rdwu_abstracts[(i * rdwu_abstract_length) + dep_start_id];
      const int depn_end_pos   = redk.rdwu_abstracts[(i * rdwu_abstract_length) + dep_end_id];
      const int depn_nval      = depn_end_pos - depn_start_pos;
      const int system_pos     = redk.rdwu_abstracts[(i * rdwu_abstract_length) + system_id];
      
      // Branch for different reduction goals.
      switch (purpose) {
      case ReductionGoal::NORMALIZE:
        scatterNormalization(*rsbs, sum<double>(&rsbs->x_buffer[depn_start_pos], depn_nval),
                             start_pos, end_pos);
        break;
      case ReductionGoal::CENTER_ON_ZERO:
        {
          const double tsum_x = sum<double>(&rsbs->x_buffer[depn_start_pos], depn_nval);
          double tsum_y = 0.0;
          double tsum_z = 0.0;
          if (rsbs->y_read != nullptr) {
            tsum_y = sum<double>(&rsbs->y_buffer[depn_start_pos], depn_nval);
          }
          if (rsbs->z_read != nullptr) {
            tsum_z = sum<double>(&rsbs->z_buffer[depn_start_pos], depn_nval);
          }
          scatterCenterOnZero(*rsbs, tsum_x, tsum_y, tsum_z, redk.atom_counts[system_pos],
                              start_pos, end_pos);
        }
        break;
      case ReductionGoal::CONJUGATE_GRADIENT:
        break;
      }
    }
    break;
  case ReductionStage::RESCALE:
    for (int i = 0; i < redk.nrdwu; i++) {
      const int start_pos  = redk.rdwu_abstracts[(i * rdwu_abstract_length) + atom_start_id];
      const int end_pos    = redk.rdwu_abstracts[(i * rdwu_abstract_length) + atom_end_id];
      const int result_pos = redk.rdwu_abstracts[(i * rdwu_abstract_length) + system_id];
      const double xscale  = rsbs->x_buffer[result_pos];
      for (int j = start_pos; j < end_pos; j++) {
        const double dx = static_cast<double>(rsbs->x_read[j]);
        rsbs->x_write[j] = static_cast<T>(dx * xscale);
      }
      if (rsbs->y_read != nullptr) {
        const double yscale  = rsbs->y_buffer[result_pos];
        for (int j = start_pos; j < end_pos; j++) {
          const double dy = static_cast<double>(rsbs->y_read[j]);
          rsbs->y_write[j] = static_cast<T>(dy * yscale);
        }
      }
      if (rsbs->z_read != nullptr) {
        const double zscale  = rsbs->z_buffer[result_pos];
        for (int j = start_pos; j < end_pos; j++) {
          const double dz = static_cast<double>(rsbs->z_read[j]);
          rsbs->z_write[j] = static_cast<T>(dz * zscale);
        }
      }
    }
    break;
  case ReductionStage::ALL_REDUCE:
    switch (redk.rps) {
    case RdwuPerSystem::ONE:
      for (int i = 0; i < redk.nrdwu; i++) {
        const int start_pos      = redk.rdwu_abstracts[(i * rdwu_abstract_length) + atom_start_id];
        const int end_pos        = redk.rdwu_abstracts[(i * rdwu_abstract_length) + atom_end_id];
        const int system_pos     = redk.rdwu_abstracts[(i * rdwu_abstract_length) + system_id];
        switch (purpose) {
        case ReductionGoal::NORMALIZE:
          {
            const double tsum = gatherNormalization(*rsbs, start_pos, end_pos);
            scatterNormalization(*rsbs, tsum, start_pos, end_pos);
          }
          break;
        case ReductionGoal::CENTER_ON_ZERO:
          {
            const double3 tsum3 = gatherCenterOnZero(*rsbs, start_pos, end_pos);
            scatterCenterOnZero(*rsbs, tsum3.x, tsum3.y, tsum3.z, redk.atom_counts[system_pos],
                                start_pos, end_pos);
          }
          break;
        case ReductionGoal::CONJUGATE_GRADIENT:
          break;
        }
      }
      break;
    case RdwuPerSystem::MULTIPLE:
      evalReduction(rsbs, redk, ReductionStage::GATHER, purpose);
      evalReduction(rsbs, redk, ReductionStage::SCATTER, purpose);
      break;
    }
    break;
  }
}

} // namespace stmath
} // namespace stormm
