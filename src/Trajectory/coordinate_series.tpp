// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeriesWriter<T>::CoordinateSeriesWriter(const int natom_in, const int nframe_in,
                                                  const UnitCellType unit_cell_in,
                                                  const int gpos_bits_in,
                                                  const double gpos_scale_in,
                                                  const double inv_gpos_scale_in, T* xcrd_in,
                                                  T* ycrd_in, T* zcrd_in, double* umat_in,
                                                  double* invu_in, double* boxdim_in) :
    natom{natom_in}, nframe{nframe_in}, unit_cell{unit_cell_in}, gpos_bits{gpos_bits_in},
    gpos_scale{gpos_scale_in}, inv_gpos_scale{inv_gpos_scale_in}, xcrd{xcrd_in}, ycrd{ycrd_in},
    zcrd{zcrd_in}, umat{umat_in}, invu{invu_in}, boxdim{boxdim_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeriesReader<T>::CoordinateSeriesReader(const int natom_in, const int nframe_in,
                                                  const UnitCellType unit_cell_in,
                                                  const int gpos_bits_in,
                                                  const double gpos_scale_in,
                                                  const double inv_gpos_scale_in, const T* xcrd_in,
                                                  const T* ycrd_in, const T* zcrd_in,
                                                  const double* umat_in, const double* invu_in,
                                                  const double* boxdim_in) :
    natom{natom_in}, nframe{nframe_in}, unit_cell{unit_cell_in}, gpos_bits{gpos_bits_in},
    gpos_scale{gpos_scale_in}, inv_gpos_scale{inv_gpos_scale_in}, xcrd{xcrd_in}, ycrd{ycrd_in},
    zcrd{zcrd_in}, umat{umat_in}, invu{invu_in}, boxdim{boxdim_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeriesReader<T>::CoordinateSeriesReader(const CoordinateSeriesWriter<T> &csw) :
    natom{csw.natom}, nframe{csw.nframe}, unit_cell{csw.unit_cell}, gpos_bits{csw.gpos_bits},
    gpos_scale{csw.gpos_scale}, inv_gpos_scale{csw.inv_gpos_scale}, xcrd{csw.xcrd}, ycrd{csw.ycrd},
    zcrd{csw.zcrd}, umat{csw.umat}, invu{csw.invu}, boxdim{csw.boxdim}
{}
  
//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T>::CoordinateSeries(const int natom_in, const int nframe_in,
                                      const UnitCellType unit_cell_in,
                                      const int globalpos_scale_bits_in,
                                      const HybridFormat format_in) :
    format{format_in}, atom_count{natom_in}, frame_count{nframe_in}, frame_capacity{nframe_in},
    globalpos_scale_bits{globalpos_scale_bits_in}, unit_cell{unit_cell_in},
    globalpos_scale{pow(2.0, globalpos_scale_bits)},
    inverse_globalpos_scale{1.0 / globalpos_scale},
    x_coordinates{static_cast<size_t>(nframe_in) * roundUp<size_t>(natom_in, warp_size_zu),
                  "cser_x_coords", format},
    y_coordinates{static_cast<size_t>(nframe_in) * roundUp<size_t>(natom_in, warp_size_zu),
                  "cser_y_coords", format},
    z_coordinates{static_cast<size_t>(nframe_in) * roundUp<size_t>(natom_in, warp_size_zu),
                  "cser_z_coords", format},
    box_space_transforms{static_cast<size_t>(nframe_in) * roundUp<size_t>(9, warp_size_zu),
                         "cser_umat", format},
    inverse_transforms{static_cast<size_t>(nframe_in) * roundUp<size_t>(9, warp_size_zu),
                       "cser_invu", format},
    box_dimensions{static_cast<size_t>(nframe_in) * roundUp<size_t>(6, warp_size_zu),
                   "cser_boxdims", format}
{
  // Limits on the valid data types
  if (isFloatingPointScalarType<T>() == false && isSignedIntegralScalarType<T>() == false) {
    rtErr("A CoordinateSeries object is only valid with a signed integer or real number scalar "
          "data type.", "CoordinateSeries");
  }
  if (globalpos_scale_bits < 0) {
    rtErr("A fixed precision representation cannot have a negative bit count (" +
          std::to_string(globalpos_scale_bits) + ") after the decimal.", "CoordinateSeries");
  }
  if (isFloatingPointScalarType<T>() && globalpos_scale_bits != 0) {
    rtErr("A real-numbered representation of coordinates is incompatible with a fixed-precision "
          "representation of " + std::to_string(globalpos_scale_bits) + " bits after the decimal.",
          "CoordinateSeries");
  }
  if (isSignedIntegralScalarType<T>() && globalpos_scale_bits == 0) {
    globalpos_scale_bits = default_trajpos_scale_bits;
    globalpos_scale = pow(2.0, globalpos_scale_bits);
    inverse_globalpos_scale = 1.0 / globalpos_scale;
  }
  allocate(frame_count);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T>::CoordinateSeries(const std::string &file_name, const int atom_count_in,
                                      const CoordinateFileKind file_kind,
                                      const std::vector<int> &frame_numbers,
                                      const int replica_count, const UnitCellType unit_cell_in,
                                      const int globalpos_scale_bits_in,
                                      const HybridFormat format_in) :
    CoordinateSeries(atom_count_in, replica_count * static_cast<int>(frame_numbers.size()),
                     unit_cell_in, globalpos_scale_bits_in, format_in)
{
  importFromFile(file_name, file_kind, frame_numbers, replica_count, 0);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T>::CoordinateSeries(PhaseSpace *ps, const int nframe_in,
                                      const int globalpos_scale_bits_in,
                                      const HybridFormat format_in) :
    CoordinateSeries(ps->getAtomCount(), nframe_in, ps->getUnitCellType(), globalpos_scale_bits_in,
                     format_in)
{
  for (int i = 0; i < frame_count; i++) {
    import(ps, 0, ps->getAtomCount(), i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T>::CoordinateSeries(const PhaseSpace &ps, const int nframe_in,
                                      const int globalpos_scale_bits_in,
                                      const HybridFormat format_in) :
    CoordinateSeries(ps.getAtomCount(), nframe_in, ps.getUnitCellType(), globalpos_scale_bits_in,
                     format_in)
{
  for (int i = 0; i < frame_count; i++) {
    import(ps, 0, ps.getAtomCount(), i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T>::CoordinateSeries(CoordinateFrame *cf, const int nframe_in,
                                      const int globalpos_scale_bits_in,
                                      const HybridFormat format_in) :
    CoordinateSeries(cf->getAtomCount(), nframe_in, cf->getUnitCellType(), globalpos_scale_bits_in,
                     format_in)
{
  for (int i = 0; i < frame_count; i++) {
    import(cf, 0, cf->getAtomCount(), i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T>::CoordinateSeries(const CoordinateFrame &cf, const int nframe_in,
                                      const int globalpos_scale_bits_in,
                                      const HybridFormat format_in) :
    CoordinateSeries(cf.getAtomCount(), nframe_in, cf.getUnitCellType(), globalpos_scale_bits_in,
                     format_in)
{
  for (int i = 0; i < frame_count; i++) {
    import(cf, 0, cf.getAtomCount(), i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeries<T>::CoordinateSeries(const CoordinateSeries<T> &original,
                                      const HybridFormat format_in) :
    CoordinateSeries(original.atom_count, original.frame_count, original.unit_cell,
                     original.globalpos_scale_bits, format_in)
{
  // After delegating to the basic constructor, fill the new object with deep copies of the data
  // from the original object.
  deepCopy(&x_coordinates, original.x_coordinates);
  deepCopy(&y_coordinates, original.y_coordinates);
  deepCopy(&z_coordinates, original.z_coordinates);
  deepCopy(&box_space_transforms, original.box_space_transforms);
  deepCopy(&inverse_transforms, original.inverse_transforms);
  deepCopy(&box_dimensions, original.box_dimensions);
}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Toriginal>
CoordinateSeries<T>::CoordinateSeries(const CoordinateSeries<Toriginal> &original,
                                      const int globalpos_scale_bits_in,
                                      const HybridFormat format_in) :
    CoordinateSeries(original.getAtomCount(), original.getFrameCount(), original.getUnitCellType(),
                     globalpos_scale_bits_in, format_in)
{
  for (int i = 0; i < frame_count; i++) {
    const CoordinateFrame cf = original.exportFrame(i);
    import(cf, 0, cf.getAtomCount(), i);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T> HybridFormat CoordinateSeries<T>::getFormat() const {
  return format;
}

//-------------------------------------------------------------------------------------------------
template <typename T> int CoordinateSeries<T>::getAtomCount() const {
  return atom_count;
}
  
//-------------------------------------------------------------------------------------------------
template <typename T> int CoordinateSeries<T>::getFrameCount() const {
  return frame_count;
}

//-------------------------------------------------------------------------------------------------
template <typename T> int CoordinateSeries<T>::getFrameCapacity() const {
  return frame_capacity;
}

//-------------------------------------------------------------------------------------------------
template <typename T> UnitCellType CoordinateSeries<T>::getUnitCellType() const {
  return unit_cell;
}

//-------------------------------------------------------------------------------------------------
template <typename T> int CoordinateSeries<T>::getFixedPrecisionBits() const {
  if (isFloatingPointScalarType<T>()) {
    rtWarn("A CoordinateSeries object with a real-numbered data representation does not have a "
           "meaningful fixed-precision scaling factor.", "CoordinateSeries",
           "getFixedPrecisionBits");
  }
  return globalpos_scale_bits;
}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Treport> std::vector<Treport>
CoordinateSeries<T>::getInterlacedCoordinates(const int frame_index,
                                              const int globalpos_bits_out,
                                              const HybridTargetLevel tier) const {
  return getInterlacedCoordinates<Treport>(frame_index, 0, atom_count, globalpos_bits_out, tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Treport> std::vector<Treport>
CoordinateSeries<T>::getInterlacedCoordinates(const int frame_index, const int low_index,
                                              const int high_index, const int globalpos_bits_out,
                                              const HybridTargetLevel tier) const {

  // Range check as this will use pointers
  if (low_index < 0 || high_index > atom_count || high_index <= low_index) {
    rtErr("Invalid atom range " + std::to_string(low_index) + " to " + std::to_string(high_index) +
          " in an object with " + std::to_string(atom_count) + " atoms.", "CoordinateSeries",
          "getInterlacedCoordinates");
  }
  checkFormatCompatibility(tier, format, "CoordinateSeries", "getInterlacedCoordinates");
  const int actual_bits_out = (globalpos_bits_out < 0) ? globalpos_scale_bits : globalpos_bits_out;
  std::vector<Treport> result(3 * (high_index - low_index));
  const size_t fidx_zu  = static_cast<size_t>(frame_index);
  const size_t natom_zu = roundUp(static_cast<size_t>(atom_count), warp_size_zu);
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const T* xptr = x_coordinates.data();
      const T* yptr = y_coordinates.data();
      const T* zptr = z_coordinates.data();
      const size_t frame_offset = natom_zu * fidx_zu;
      if (isSignedIntegralScalarType<T>() && isSignedIntegralScalarType<Treport>()) {
        if (globalpos_scale_bits >= actual_bits_out) {
          llint divisor = 1LL;
          for (int i = actual_bits_out; i < globalpos_scale_bits; i++) {
            divisor *= 2LL;
          }
          for (int i = low_index; i < high_index; i++) {
            const int base_idx = 3 * (i - low_index);
            const size_t access_idx = frame_offset + static_cast<size_t>(i);
            result[base_idx    ] = static_cast<llint>(xptr[access_idx]) / divisor;
            result[base_idx + 1] = static_cast<llint>(yptr[access_idx]) / divisor;
            result[base_idx + 2] = static_cast<llint>(zptr[access_idx]) / divisor;
          }
        }
        else {
          const int shift_left = actual_bits_out - globalpos_scale_bits;
          llint multiplier = 1LL;
          for (int i = globalpos_scale_bits; i < actual_bits_out; i++) {
            multiplier *= 2LL;
          }
          for (int i = low_index; i < high_index; i++) {
            const int base_idx = 3 * (i - low_index);
            const size_t access_idx = frame_offset + static_cast<size_t>(i);
            result[base_idx    ] = static_cast<llint>(xptr[access_idx]) * multiplier;
            result[base_idx + 1] = static_cast<llint>(yptr[access_idx]) * multiplier;
            result[base_idx + 2] = static_cast<llint>(zptr[access_idx]) * multiplier;
          }
        }
      }
      else if (isSignedIntegralScalarType<T>() && isFloatingPointScalarType<Treport>()) {
        const Treport multiplier = 1.0 / pow(2.0, globalpos_scale_bits);
        for (int i = low_index; i < high_index; i++) {
          const int base_idx = 3 * (i - low_index);
          const size_t access_idx = frame_offset + static_cast<size_t>(i);
          result[base_idx    ]  = xptr[access_idx];
          result[base_idx + 1]  = yptr[access_idx];
          result[base_idx + 2]  = zptr[access_idx];
          result[base_idx    ] *= multiplier;
          result[base_idx + 1] *= multiplier;
          result[base_idx + 2] *= multiplier;
        }
      }
      else if (isFloatingPointScalarType<T>() && isSignedIntegralScalarType<Treport>()) {
        const Treport multiplier = pow(2.0, actual_bits_out);
        for (int i = low_index; i < high_index; i++) {
          const int base_idx = 3 * (i - low_index);
          const size_t access_idx = frame_offset + static_cast<size_t>(i);
          result[base_idx    ]  = llround(xptr[access_idx] * multiplier);
          result[base_idx + 1]  = llround(yptr[access_idx] * multiplier);
          result[base_idx + 2]  = llround(zptr[access_idx] * multiplier);
        }
      }
      else {
        for (int i = low_index; i < high_index; i++) {
          const int base_idx = 3 * (i - low_index);
          const size_t access_idx = frame_offset + static_cast<size_t>(i);
          result[base_idx    ] = xptr[access_idx];
          result[base_idx + 1] = yptr[access_idx];
          result[base_idx + 2] = zptr[access_idx];
        }
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const size_t llim = (fidx_zu * natom_zu) + static_cast<size_t>(low_index);
      const size_t hlim = (fidx_zu * natom_zu) + static_cast<size_t>(high_index);
      const std::vector<T> xval = x_coordinates.readDevice(llim, hlim);
      const std::vector<T> yval = y_coordinates.readDevice(llim, hlim);
      const std::vector<T> zval = z_coordinates.readDevice(llim, hlim);
      const int irange = high_index - low_index;
      if (isSignedIntegralScalarType<T>() && isSignedIntegralScalarType<Treport>()) {
        if (globalpos_scale_bits >= actual_bits_out) {
          const int shift_right = globalpos_scale_bits - actual_bits_out;
          for (int i = 0; i < irange; i++) {
            result[(3 * i)    ] = (static_cast<llint>(xval[i]) >> shift_right);
            result[(3 * i) + 1] = (static_cast<llint>(yval[i]) >> shift_right);
            result[(3 * i) + 2] = (static_cast<llint>(zval[i]) >> shift_right);
          }
        }
        else {
          const int shift_left = actual_bits_out - globalpos_scale_bits;
          for (int i = 0; i < high_index - low_index; i++) {
            result[(3 * i)    ] = (static_cast<llint>(xval[i]) << shift_left);
            result[(3 * i) + 1] = (static_cast<llint>(yval[i]) << shift_left);
            result[(3 * i) + 2] = (static_cast<llint>(zval[i]) << shift_left);
          }
        }
      }
      else if (isSignedIntegralScalarType<T>() && isFloatingPointScalarType<Treport>()) {
        const Treport multiplier = 1.0 / pow(2.0, globalpos_scale_bits);
        for (int i = 0; i < irange; i++) {
          result[(3 * i)    ]  = xval[i];
          result[(3 * i) + 1]  = yval[i];
          result[(3 * i) + 2]  = zval[i];
          result[(3 * i)    ] *= multiplier;
          result[(3 * i) + 1] *= multiplier;
          result[(3 * i) + 2] *= multiplier;
        }
      }
      else if (isFloatingPointScalarType<T>() && isSignedIntegralScalarType<Treport>()) {
        const Treport multiplier = pow(2.0, globalpos_scale_bits);
        for (int i = 0; i < irange; i++) {
          result[(3 * i)    ]  = llround(xval[i] * multiplier);
          result[(3 * i) + 1]  = llround(yval[i] * multiplier);
          result[(3 * i) + 2]  = llround(zval[i] * multiplier);
        }
      }
      else {
        for (int i = 0; i < irange; i++) {
          result[(3 * i)    ] = xval[i];
          result[(3 * i) + 1] = yval[i];
          result[(3 * i) + 2] = zval[i];
        }
      }
    }
    break;
#endif
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<double> CoordinateSeries<T>::getBoxSpaceTransform(const int frame_index,
                                                              const HybridTargetLevel tier) const {
  checkFormatCompatibility(tier, format, "CoordinateSeries", "getBoxSpaceTransform");
  const size_t read_start = static_cast<size_t>(frame_index) * roundUp<size_t>(9, warp_size_zu);
  switch (tier) {
  case HybridTargetLevel::HOST:
    return box_space_transforms.readHost(read_start, 9);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return box_space_transforms.readDevice(read_start, 9);
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<double> CoordinateSeries<T>::getInverseTransform(const int frame_index,
                                                             const HybridTargetLevel tier) const {
  checkFormatCompatibility(tier, format, "CoordinateSeries", "getInverseTransform");
  const size_t read_start = static_cast<size_t>(frame_index) * roundUp<size_t>(9, warp_size_zu);
  switch (tier) {
  case HybridTargetLevel::HOST:
    return inverse_transforms.readHost(read_start, 9);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return inverse_transforms.readDevice(read_start, 9);
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<double> CoordinateSeries<T>::getBoxDimensions(const int frame_index,
                                                          const HybridTargetLevel tier) const {
  checkFormatCompatibility(tier, format, "CoordinateSeries", "getBoxDimensions");
  const size_t read_start = static_cast<size_t>(frame_index) * roundUp<size_t>(6, warp_size_zu);
  switch (tier) {
  case HybridTargetLevel::HOST:
    return box_dimensions.readHost(read_start, 6);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    return box_dimensions.readDevice(read_start, 6);
    break;
#endif
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const Hybrid<double>& CoordinateSeries<T>::getBoxDimensions() const {
  return box_dimensions;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::extractFrame(CoordinateFrame *cf, const size_t frame_index,
                                       const GpuDetails &gpu) const {

  // Transfer coordinates
  Hybrid<double>* x_ptr = cf->getCoordinateHandle(CartesianDimension::X);
  Hybrid<double>* y_ptr = cf->getCoordinateHandle(CartesianDimension::Y);
  Hybrid<double>* z_ptr = cf->getCoordinateHandle(CartesianDimension::Z);
  const size_t fidx = frame_index * roundUp<size_t>(atom_count, warp_size_zu);
  const size_t natom = cf->getAtomCount();
  deepRecast<double, T>(x_ptr, &x_coordinates, cf->getAtomCount(), natom, 0, fidx, 0,
                        globalpos_scale_bits, gpu); 
  deepRecast<double, T>(y_ptr, &y_coordinates, cf->getAtomCount(), natom, 0, fidx, 0,
                        globalpos_scale_bits, gpu); 
  deepRecast<double, T>(z_ptr, &z_coordinates, cf->getAtomCount(), natom, 0, fidx, 0,
                        globalpos_scale_bits, gpu);

  // Transfer box dimensions and transformation matrices
  Hybrid<double>* umat_ptr = cf->getBoxTransformHandle();
  Hybrid<double>* invu_ptr = cf->getInverseTransformHandle();
  Hybrid<double>* bdim_ptr = cf->getBoxDimensionsHandle();
  const size_t matrix_stride = roundUp(9, warp_size_int);
  const size_t boxdim_stride = roundUp(6, warp_size_int);
  deepCopy<double>(umat_ptr, box_space_transforms, 9, 0, matrix_stride * frame_index, gpu);
  deepCopy<double>(invu_ptr, inverse_transforms, 9, 0, matrix_stride * frame_index, gpu);
  deepCopy<double>(bdim_ptr, box_dimensions, 6, 0, boxdim_stride * frame_index, gpu);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::extractFrame(CoordinateFrame *cf, const size_t frame_index,
                                       const HybridTargetLevel tier) const {
  checkFormatCompatibility(tier, format, "CoordinateSeries", "extractFrame");
  if (cf->getAtomCount() != atom_count) {
    rtErr("Atom count of the destination CoordinateFrame (" + std::to_string(cf->getAtomCount()) +
          ") does not match the object (" + std::to_string(atom_count) + ").", "CoordinateSeries",
          "extractFrame");
  }
  const size_t frame_offset = roundUp<size_t>(atom_count, warp_size_zu) * frame_index;
  const size_t bdim_offset  = roundUp<size_t>(6, warp_size_zu) * frame_index;
  switch (tier) {
  case HybridTargetLevel::HOST:
    cf->fill(&x_coordinates.data()[frame_offset], &y_coordinates.data()[frame_offset],
             &z_coordinates.data()[frame_offset], globalpos_scale_bits,
             &box_dimensions.data()[bdim_offset]);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const size_t natom_zu = atom_count;
      const std::vector<T> tmp_xcrd = x_coordinates.readDevice(frame_offset, natom_zu);
      const std::vector<T> tmp_ycrd = y_coordinates.readDevice(frame_offset, natom_zu);
      const std::vector<T> tmp_zcrd = z_coordinates.readDevice(frame_offset, natom_zu);
      const std::vector<double> tmp_bdim = box_dimensions.readDevice(bdim_offset, 6);
      cf->fill(tmp_xcrd, tmp_ycrd, tmp_zcrd, globalpos_scale_bits, tmp_bdim);
    }
    break;
#endif
  }

  // Set the frame number based on the position in the series
  cf->setFrameNumber(frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::extractFrame(PhaseSpace *ps, const size_t frame_index,
                                       const TrajectoryKind kind, const CoordinateCycle time_point,
                                       const HybridTargetLevel tier) const {
  checkFormatCompatibility(tier, format, "CoordinateSeries", "extractFrame");
  if (ps->getAtomCount() != atom_count) {
    rtErr("Atom count of the destination CoordinateFrame (" + std::to_string(ps->getAtomCount()) +
          ") does not match the object (" + std::to_string(atom_count) + ").", "CoordinateSeries",
          "extractFrame");
  }
  const size_t frame_offset = roundUp<size_t>(atom_count, warp_size_zu) * frame_index;
  const size_t bdim_offset  = roundUp<size_t>(6, warp_size_zu) * frame_index;
  switch (tier) {
  case HybridTargetLevel::HOST:
    ps->fill(&x_coordinates.data()[frame_offset], &y_coordinates.data()[frame_offset],
             &z_coordinates.data()[frame_offset], kind, time_point, globalpos_scale_bits,
             &box_dimensions.data()[bdim_offset]);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const size_t natom_zu = atom_count;
      const std::vector<T> tmp_xcrd = x_coordinates.readDevice(frame_offset, natom_zu);
      const std::vector<T> tmp_ycrd = y_coordinates.readDevice(frame_offset, natom_zu);
      const std::vector<T> tmp_zcrd = z_coordinates.readDevice(frame_offset, natom_zu);
      const std::vector<double> tmp_bdim = box_dimensions.readDevice(bdim_offset, 6);
      ps->fill(tmp_xcrd, tmp_ycrd, tmp_zcrd, kind, time_point, globalpos_scale_bits, tmp_bdim);
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::extractFrame(PhaseSpace *ps, const size_t frame_index,
                                       const TrajectoryKind kind,
                                       const HybridTargetLevel tier) const {
  extractFrame(ps, frame_index, kind, ps->getCyclePosition(), tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::extractFrame(PhaseSpace *ps, const size_t frame_index,
                                       const HybridTargetLevel tier) const {
  extractFrame(ps, frame_index, TrajectoryKind::POSITIONS, ps->getCyclePosition(), tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateFrame CoordinateSeries<T>::exportFrame(const size_t frame_index,
                                                 const HybridTargetLevel tier) const {
  checkFormatCompatibility(tier, format, "CoordinateSeries", "exportFrame");
  CoordinateFrame result(atom_count, unit_cell);
  const size_t frame_offset = roundUp<size_t>(atom_count, warp_size_zu) * frame_index;
  const size_t bdim_offset  = roundUp<size_t>(6, warp_size_zu) * frame_index;
  switch (tier) {
  case HybridTargetLevel::HOST:
    result.fill(&x_coordinates.data()[frame_offset], &y_coordinates.data()[frame_offset],
                &z_coordinates.data()[frame_offset], globalpos_scale_bits,
                &box_dimensions.data()[bdim_offset]);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const size_t natom_zu = atom_count;
      const std::vector<T> tmp_xcrd = x_coordinates.readDevice(frame_offset, natom_zu);
      const std::vector<T> tmp_ycrd = y_coordinates.readDevice(frame_offset, natom_zu);
      const std::vector<T> tmp_zcrd = z_coordinates.readDevice(frame_offset, natom_zu);
      const std::vector<double> tmp_bdim = box_dimensions.readDevice(bdim_offset, 6);
      result.fill(tmp_xcrd, tmp_ycrd, tmp_zcrd, globalpos_scale_bits, tmp_bdim);
    }
    break;
#endif
  }
  
  // Set the frame number based on the position in the series
  result.setFrameNumber(frame_index);
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
PhaseSpace CoordinateSeries<T>::exportPhaseSpace(const size_t frame_index,
                                                 const HybridTargetLevel tier) const {
  checkFormatCompatibility(tier, format, "CoordinateSeries", "exportPhaseSpace");
  PhaseSpace result(atom_count, unit_cell);
  const size_t frame_offset = roundUp<size_t>(atom_count, warp_size_zu) * frame_index;
  const size_t bdim_offset  = roundUp<size_t>(6, warp_size_zu) * frame_index;
  switch (tier) {
  case HybridTargetLevel::HOST:
    result.fill(&x_coordinates.data()[frame_offset], &y_coordinates.data()[frame_offset],
                &z_coordinates.data()[frame_offset], TrajectoryKind::POSITIONS,
                CoordinateCycle::WHITE, globalpos_scale_bits,
                &box_dimensions.data()[bdim_offset]);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const size_t natom_zu = atom_count;
      const std::vector<T> tmp_xcrd = x_coordinates.readDevice(frame_offset, natom_zu);
      const std::vector<T> tmp_ycrd = y_coordinates.readDevice(frame_offset, natom_zu);
      const std::vector<T> tmp_zcrd = z_coordinates.readDevice(frame_offset, natom_zu);
      const std::vector<double> tmp_bdim = box_dimensions.readDevice(bdim_offset, 6);
      result.fill(tmp_xcrd, tmp_ycrd, tmp_zcrd,	TrajectoryKind::POSITIONS,
                  CoordinateCycle::WHITE, globalpos_scale_bits, tmp_bdim);
    }
    break;
#endif
  }
  return result;  
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::exportToFile(const std::string &file_name, const CoordinateFileKind kind,
                                       const PrintSituation expectation, const int low_index,
                                       const int high_index, const HybridTargetLevel tier) const {
  checkFormatCompatibility(tier, format, "CoordinateSeries", "exportToFile");
  if (low_index < 0 || low_index >= frame_count || high_index < low_index ||
      high_index >= frame_count) {
    rtErr("The frame index range " + std::to_string(low_index) + " to " +
          std::to_string(high_index) + " is invalid for a series with " +
          std::to_string(frame_count) + " frames.", "CoordinateSeries", "exportToFile");
  }
  const int actual_high_index = (high_index > low_index) ? high_index : frame_count;
  const PrintSituation aexp = adjustTrajectoryOpeningProtocol(expectation, kind,
                                                              "CoordinateSeries", "exportToFile");
  const DataFormat style = getTrajectoryFormat(kind);
  const bool fi_exists = (getDrivePathType(file_name) == DrivePathType::FILE);
  std::ofstream foutp;
  foutp = openOutputFile(file_name, aexp, "Open an output file for writing CoordinateSeries "
                         "contents.", style);
  switch (kind) {
  case CoordinateFileKind::AMBER_CRD:
  case CoordinateFileKind::SDF:
  case CoordinateFileKind::AMBER_NETCDF:
  case CoordinateFileKind::AMBER_NETCDF_RST:
    if (fi_exists == false ||
        aexp == PrintSituation::OVERWRITE || aexp == PrintSituation::OPEN_NEW) {
      initializeTrajectory(&foutp, kind, atom_count);
    }
    break;
  case CoordinateFileKind::AMBER_INPCRD:
  case CoordinateFileKind::AMBER_ASCII_RST:
    initializeTrajectory(&foutp, kind, atom_count);
    break;
  case CoordinateFileKind::UNKNOWN:
    rtErr("A trajectory format must be specified.", "CoordinateSeries", "exportToFile");
  }
  switch (kind) {
  case CoordinateFileKind::AMBER_CRD:
  case CoordinateFileKind::AMBER_INPCRD:
    if (kind == CoordinateFileKind::AMBER_INPCRD && actual_high_index - low_index != 1) {
      rtErr("The " + getEnumerationName(kind) + " file format requires one and only "
            "one frame.  It cannot accept a series of " +
            std::to_string(actual_high_index - low_index) + " frames.", "CoordinateSeries",
            "exportToFile");
    }
    for (int i = low_index; i < actual_high_index; i++) {
      const CoordinateFrame cf = exportFrame(i, tier);
      const CoordinateFrameReader cfr = cf.data();
      writeFrame(&foutp, file_name, kind, atom_count, cfr.xcrd, cfr.ycrd, cfr.zcrd,
                 nullptr, nullptr, nullptr, cfr.unit_cell, cfr.boxdim);
    }
    break;
  case CoordinateFileKind::SDF:
    rtErr("The object does not have sufficient information to create an annotated SD file.  The "
          "program must use one of the writeFrame() overloads from the write_annotated_frame "
          "library instead.", "CoordinateSeries", "exportToFile");
    break;
  case CoordinateFileKind::AMBER_NETCDF:
    break;
  case CoordinateFileKind::AMBER_ASCII_RST:
  case CoordinateFileKind::AMBER_NETCDF_RST:
    rtErr("A restart file cannot be written based on a CoordinateSeries.  The object will not be "
          "able to store both coordinates and velocities needed for checkpointing.",
          "CoordinateFrame", "exportToFile");
    break;
  case CoordinateFileKind::UNKNOWN:
    break;
  }
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const Hybrid<T>& CoordinateSeries<T>::getCoordinateReference(CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return x_coordinates;
  case CartesianDimension::Y:
    return y_coordinates;
  case CartesianDimension::Z:
    return z_coordinates;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const Hybrid<T>* CoordinateSeries<T>::getCoordinatePointer(CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return &x_coordinates;
  case CartesianDimension::Y:
    return &y_coordinates;
  case CartesianDimension::Z:
    return &z_coordinates;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const Hybrid<double>& CoordinateSeries<T>::getBoxTransforms() const {
  return box_space_transforms;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const Hybrid<double>* CoordinateSeries<T>::getBoxTransformPointer() const {
  return &box_space_transforms;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const Hybrid<double>& CoordinateSeries<T>::getInverseTransforms() const {
  return inverse_transforms;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const Hybrid<double>* CoordinateSeries<T>::getInverseTransformPointer() const {
  return &inverse_transforms;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const Hybrid<double>* CoordinateSeries<T>::getBoxDimensionPointer() const {
  return &box_dimensions;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const Hybrid<T>* CoordinateSeries<T>::getFramesHandle(const CartesianDimension dim) const {
  switch (dim) {
  case CartesianDimension::X:
    return &x_coordinates;
  case CartesianDimension::Y:
    return &y_coordinates;
  case CartesianDimension::Z:
    return &z_coordinates;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
Hybrid<T>* CoordinateSeries<T>::getFramesHandle(const CartesianDimension dim) {
  switch (dim) {
  case CartesianDimension::X:
    return &x_coordinates;
  case CartesianDimension::Y:
    return &y_coordinates;
  case CartesianDimension::Z:
    return &z_coordinates;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const Hybrid<double>* CoordinateSeries<T>::getBoxTransformsHandle() const {
  return &box_space_transforms;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
Hybrid<double>* CoordinateSeries<T>::getBoxTransformsHandle() {
  return &box_space_transforms;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const Hybrid<double>* CoordinateSeries<T>::getInverseTransformsHandle() const {
  return &inverse_transforms;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
Hybrid<double>* CoordinateSeries<T>::getInverseTransformsHandle() {
  return &inverse_transforms;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const Hybrid<double>* CoordinateSeries<T>::getBoxDimensionsHandle() const {
  return &box_dimensions;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
Hybrid<double>* CoordinateSeries<T>::getBoxDimensionsHandle() {
  return &box_dimensions;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const CoordinateSeries<T>* CoordinateSeries<T>::getSelfPointer() const {
  return this;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const CoordinateSeriesReader<T> CoordinateSeries<T>::data(const HybridTargetLevel tier) const {
  checkFormatCompatibility(tier, format, "CoordinateSeries", "data");
  return CoordinateSeriesReader<T>(atom_count, frame_count, unit_cell, globalpos_scale_bits,
                                   globalpos_scale, inverse_globalpos_scale,
                                   x_coordinates.data(tier), y_coordinates.data(tier),
                                   z_coordinates.data(tier), box_space_transforms.data(tier),
                                   inverse_transforms.data(tier), box_dimensions.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeriesWriter<T> CoordinateSeries<T>::data(const HybridTargetLevel tier) {
  checkFormatCompatibility(tier, format, "CoordinateSeries", "data");
  return CoordinateSeriesWriter<T>(atom_count, frame_count, unit_cell, globalpos_scale_bits,
                                   globalpos_scale, inverse_globalpos_scale,
                                   x_coordinates.data(tier), y_coordinates.data(tier),
                                   z_coordinates.data(tier), box_space_transforms.data(tier),
                                   inverse_transforms.data(tier), box_dimensions.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const CoordinateSeriesReader<void>
CoordinateSeries<T>::templateFreeData(const HybridTargetLevel tier) const {
  checkFormatCompatibility(tier, format, "CoordinateSeries", "data");
  return CoordinateSeriesReader<void>(atom_count, frame_count, unit_cell, globalpos_scale_bits,
                                      globalpos_scale, inverse_globalpos_scale,
                                      reinterpret_cast<const void*>(x_coordinates.data(tier)),
                                      reinterpret_cast<const void*>(y_coordinates.data(tier)),
                                      reinterpret_cast<const void*>(z_coordinates.data(tier)),
                                      box_space_transforms.data(tier),
                                      inverse_transforms.data(tier), box_dimensions.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeriesWriter<void>
CoordinateSeries<T>::templateFreeData(const HybridTargetLevel tier) {
  checkFormatCompatibility(tier, format, "CoordinateSeries", "data");
  return CoordinateSeriesWriter<void>(atom_count, frame_count, unit_cell, globalpos_scale_bits,
                                      globalpos_scale, inverse_globalpos_scale,
                                      reinterpret_cast<void*>(x_coordinates.data(tier)),
                                      reinterpret_cast<void*>(y_coordinates.data(tier)),
                                      reinterpret_cast<void*>(z_coordinates.data(tier)),
                                      box_space_transforms.data(tier),
                                      inverse_transforms.data(tier), box_dimensions.data(tier));
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
template <typename T>
const CoordinateSeriesReader<T> CoordinateSeries<T>::deviceViewToHostData() const {
  confirmHostVisibleToGpu(format, "Host memory is not visible to the GPU in this object",
                          "CoordinateSeries", "deviceViewToHostData");
  const T* xcrd = x_coordinates.getDeviceValidHostPointer();
  const T* ycrd = y_coordinates.getDeviceValidHostPointer();
  const T* zcrd = z_coordinates.getDeviceValidHostPointer();
  const double* umat = box_space_transforms.getDeviceValidHostPointer();
  const double* invu = inverse_transforms.getDeviceValidHostPointer();
  const double* boxdim = box_dimensions.getDeviceValidHostPointer();
  return CoordinateSeriesReader<T>(atom_count, frame_count, unit_cell, globalpos_scale_bits,
                                   globalpos_scale, inverse_globalpos_scale, xcrd, ycrd, zcrd,
                                   umat, invu, boxdim);
}

//-------------------------------------------------------------------------------------------------
template <typename T> CoordinateSeriesWriter<T> CoordinateSeries<T>::deviceViewToHostData() {
  confirmHostVisibleToGpu(format, "Host memory is not visible to the GPU in this object",
                          "CoordinateSeries", "deviceViewToHostData");
  T* xcrd = x_coordinates.getDeviceValidHostPointer();
  T* ycrd = y_coordinates.getDeviceValidHostPointer();
  T* zcrd = z_coordinates.getDeviceValidHostPointer();
  double* umat = box_space_transforms.getDeviceValidHostPointer();
  double* invu = inverse_transforms.getDeviceValidHostPointer();
  double* boxdim = box_dimensions.getDeviceValidHostPointer();
  return CoordinateSeriesWriter<T>(atom_count, frame_count, unit_cell, globalpos_scale_bits,
                                   globalpos_scale, inverse_globalpos_scale, xcrd, ycrd, zcrd,
                                   umat, invu, boxdim);  
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const CoordinateSeriesReader<void> CoordinateSeries<T>::deviceViewToTemplateFreeHostData() const {
  confirmHostVisibleToGpu(format, "Host memory is not visible to the GPU in this object",
                          "CoordinateSeries", "deviceViewToTemplateFreeHostData");
  const T* xcrd = x_coordinates.getDeviceValidHostPointer();
  const T* ycrd = y_coordinates.getDeviceValidHostPointer();
  const T* zcrd = z_coordinates.getDeviceValidHostPointer();
  const void* v_xcrd = reinterpret_cast<const void*>(xcrd);
  const void* v_ycrd = reinterpret_cast<const void*>(ycrd);
  const void* v_zcrd = reinterpret_cast<const void*>(zcrd);
  const double* umat = box_space_transforms.getDeviceValidHostPointer();
  const double* invu = inverse_transforms.getDeviceValidHostPointer();
  const double* boxdim = box_dimensions.getDeviceValidHostPointer();
  return CoordinateSeriesReader<void>(atom_count, frame_count, unit_cell, globalpos_scale_bits,
                                      globalpos_scale, inverse_globalpos_scale, v_xcrd, v_ycrd,
                                      v_zcrd, umat, invu, boxdim);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeriesWriter<void> CoordinateSeries<T>::deviceViewToTemplateFreeHostData() {
  confirmHostVisibleToGpu(format, "Host memory is not visible to the GPU in this object",
                          "CoordinateSeries", "deviceViewToTemplateFreeHostData");
  T* xcrd = x_coordinates.getDeviceValidHostPointer();
  T* ycrd = y_coordinates.getDeviceValidHostPointer();
  T* zcrd = z_coordinates.getDeviceValidHostPointer();
  void* v_xcrd = reinterpret_cast<void*>(xcrd);
  void* v_ycrd = reinterpret_cast<void*>(ycrd);
  void* v_zcrd = reinterpret_cast<void*>(zcrd);
  double* umat = box_space_transforms.getDeviceValidHostPointer();
  double* invu = inverse_transforms.getDeviceValidHostPointer();
  double* boxdim = box_dimensions.getDeviceValidHostPointer();
  return CoordinateSeriesWriter<void>(atom_count, frame_count, unit_cell, globalpos_scale_bits,
                                      globalpos_scale, inverse_globalpos_scale, v_xcrd, v_ycrd,
                                      v_zcrd, umat, invu, boxdim);
}

//-------------------------------------------------------------------------------------------------
template <typename T> void CoordinateSeries<T>::upload() {
  x_coordinates.upload();
  y_coordinates.upload();
  z_coordinates.upload();
  box_space_transforms.upload();
  inverse_transforms.upload();
  box_dimensions.upload();
}

//-------------------------------------------------------------------------------------------------
template <typename T> void CoordinateSeries<T>::download() {
  x_coordinates.download();
  y_coordinates.download();
  z_coordinates.download();
  box_space_transforms.download();
  inverse_transforms.download();
  box_dimensions.download();
}
#endif

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::import(const CoordinateFrameReader &cfr, const int atom_start,
                                 const int atom_end, const int frame_index) {
  const int actual_atom_end = (atom_end > atom_start) ? atom_end : cfr.natom;
  if (actual_atom_end - atom_start != atom_count) {
    rtErr("A CoordinateSeries with frames of " + std::to_string(atom_count) + " atoms cannot "
          "accept " + std::to_string(actual_atom_end - atom_start) + " atoms from a system with " +
          std::to_string(cfr.natom) + " atoms.", "CoordinateSeries", "import");
  }
  
  // Compute the actual frame index that shall be written, and the upper limit of the atoms that
  // will be written.  If the frame index exceeds the current capacity, make more capacity.
  const int actual_frame_index = (frame_index == -1) ? frame_count : frame_index;
  const size_t atom_offset = static_cast<size_t>(actual_frame_index) *
                             roundUp(static_cast<size_t>(atom_count), warp_size_zu);
  const size_t frame_limit = atom_offset + static_cast<size_t>(actual_atom_end - atom_start);

  // Exapnd the capacity to at least 125% of its original size to avoid continuous resizing of the
  // arrays as more frames are added.
  if (actual_frame_index >= frame_capacity) {
    allocate(((actual_frame_index * 5) + 3) / 4);
  }
  T* xptr = x_coordinates.data();
  T* yptr = y_coordinates.data();
  T* zptr = z_coordinates.data();
  double* umat_ptr = box_space_transforms.data();
  double* invu_ptr = inverse_transforms.data();
  double* bdim_ptr = box_dimensions.data();
  size_t jcon = atom_start;
  if (globalpos_scale_bits == 0) {
    for (size_t j = atom_offset; j < frame_limit; j++) {
      xptr[j] = cfr.xcrd[jcon];
      yptr[j] = cfr.ycrd[jcon];
      zptr[j] = cfr.zcrd[jcon];
      jcon++;
    }
  }
  else {
    const double lgpos_scale = globalpos_scale;
    for (size_t j = atom_offset; j < frame_limit; j++) {
      xptr[j] = llround(cfr.xcrd[jcon] * lgpos_scale);
      yptr[j] = llround(cfr.ycrd[jcon] * lgpos_scale);
      zptr[j] = llround(cfr.zcrd[jcon] * lgpos_scale);
      jcon++;
    }
  }
  const size_t xfrm_offset = actual_frame_index * roundUp(static_cast<size_t>(9), warp_size_zu);
  const size_t bdim_offset = actual_frame_index * roundUp(static_cast<size_t>(6), warp_size_zu);
  const size_t xfrm_limit = xfrm_offset + 9LLU;
  const size_t bdim_limit = bdim_offset + 9LLU;
  jcon = 0LLU;
  for (size_t j = xfrm_offset; j < xfrm_limit; j++) {
    umat_ptr[j] = cfr.umat[jcon];
    invu_ptr[j] = cfr.invu[jcon];
    jcon++;
  }
  jcon = 0LLU;
  for (size_t j = bdim_offset; j < bdim_limit; j++) {
    bdim_ptr[j] = cfr.boxdim[jcon];
    jcon++;
  }

  // Update the number of actual frames, if appropriate
  if (actual_frame_index >= frame_count) {
    frame_count = actual_frame_index + 1;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::import(const CoordinateFrameReader &cfr, const int frame_index) {
  import(cfr, 0, cfr.natom, frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::import(const CoordinateFrameWriter &cfw, const int atom_start,
                                 const int atom_end, const int frame_index) {
  import(CoordinateFrameReader(cfw), atom_start, atom_end, frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::import(const CoordinateFrameWriter &cfw, const int frame_index) {
  import(CoordinateFrameReader(cfw), 0, cfw.natom, frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::import(const CoordinateFrame &cf, const int atom_start,
                                 const int atom_end, const int frame_index) {
  import(cf.data(), atom_start, atom_end, frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::import(const CoordinateFrame *cf, const int frame_index) {
  import(cf->data(), 0, cf->getAtomCount(), frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::import(const PhaseSpace &ps, const int frame_index,
                                 const TrajectoryKind kind, const CoordinateCycle orientation) {
  import(CoordinateFrameReader(ps, kind, orientation), 0, ps.getAtomCount(), frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::import(const PhaseSpace *ps, const int frame_index,
                                 const TrajectoryKind kind, const CoordinateCycle orientation) {
  import(CoordinateFrameReader(ps, kind, orientation), 0, ps->getAtomCount(), frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::import(const PhaseSpace &ps, const int atom_start, const int atom_end,
                                 const int frame_index, const TrajectoryKind kind,
                                 const CoordinateCycle orientation) {
  import(CoordinateFrameReader(ps, kind, orientation), atom_start, atom_end, frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::import(const PhaseSpace *ps, const int atom_start, const int atom_end,
                                 const int frame_index, const TrajectoryKind kind,
                                 const CoordinateCycle orientation) {
  import(CoordinateFrameReader(ps, kind, orientation), atom_start, atom_end, frame_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::importFromFile(const std::string &file_name,
                                         const CoordinateFileKind file_kind,
                                         const std::vector<int> &frame_numbers,
                                         const int replica_count, const int frame_index_start) {
  
  // Try to detect the file format if it is not already specified.  If it remains UNKNOWN, that
  // will ultimately lead to an error.
  CoordinateFileKind actual_kind = file_kind;
  if (file_kind == CoordinateFileKind::UNKNOWN) {
    actual_kind = detectCoordinateFileKind(file_name);
  }
  switch (actual_kind) {
  case CoordinateFileKind::AMBER_CRD:
    {
      // The number of atoms must be known a-priori in order to read from a .crd trajectory file.
      if (atom_count == 0) {
        rtErr("A number of atoms matching the trajectory must be known prior to reading a .crd "
              "file.", "CoordinateSeries", "buildFromFile");
      }
      TextFile tf(file_name);
      std::vector<int> actual_frame_numbers;
      if (frame_numbers.size() == 0LLU) {
        const size_t nframe_detected = countAmberCrdFormatFrames(tf, atom_count, unit_cell);
        actual_frame_numbers.resize(nframe_detected);
        for (int i = 0; i < nframe_detected; i++) {
          actual_frame_numbers[i] = i;
        }
      }
      else {
        actual_frame_numbers.insert(actual_frame_numbers.end(), frame_numbers.begin(),
                                    frame_numbers.end());
      }
      const int nimports = actual_frame_numbers.size();
      const int orig_frame_count = frame_count;
      if (frame_index_start < 0) {
        resize((nimports * replica_count) + frame_count);
      }
      else {
        resize((nimports * replica_count) + frame_index_start);
      }

      // The file must be read in as double precision, despite the multiple modes that a
      // CoordinateSeries can operate in.  Buffer the data in a CoordinateFrame object, feeding
      // its pointers to a low-level overload of the reading function.  Transfer the results to
      // the CoordinateSeries.
      CoordinateFrame tmp_cf(atom_count, unit_cell);
      CoordinateFrameWriter tmp_cfw = tmp_cf.data();
      const CoordinateFrameReader tmp_cfr(tmp_cfw);
      T* xptr = x_coordinates.data();
      T* yptr = y_coordinates.data();
      T* zptr = z_coordinates.data();
      double* umat_ptr = box_space_transforms.data();
      double* invu_ptr = inverse_transforms.data();
      double* bdim_ptr = box_dimensions.data();
      for (int i = 0; i < nimports; i++) {
        readAmberCrdFormat(tf, tmp_cfw.xcrd, tmp_cfw.ycrd, tmp_cfw.zcrd, tmp_cfw.natom, unit_cell,
                           tmp_cfw.umat, tmp_cfw.invu, tmp_cfw.boxdim, actual_frame_numbers[i]);
        for (int j = 0; j < replica_count; j++) {
          if (frame_index_start < 0) {
            import(tmp_cfr, orig_frame_count + i + (j * nimports));
          }
          else {
            import(tmp_cfr, frame_index_start + i + (j * nimports));            
          }
        }
      }
    }
    break;
  case CoordinateFileKind::AMBER_INPCRD:
  case CoordinateFileKind::AMBER_ASCII_RST:
    {
      // The number of atoms need not be known when first reading an Amber ASCII input coordinate
      // or restart file--it is present in those files, and any NetCDF coordinate file.  However,
      // if the number of atoms currently in the CoordinateSeries is not zero and there are one or
      // more frames, assume that it is a legitimate number and enforce that any subsequent files
      // match that number of atoms.
      TextFile tf(file_name);
      const int test_atom_count = getAmberRestartAtomCount(tf);
      if (atom_count > 0 && frame_count > 0 && test_atom_count != atom_count) {
        rtErr("The atom count detected in " + file_name + " (" + std::to_string(test_atom_count) +
              " does not agree with the atom count of " + std::to_string(frame_count) +
              " frames of the existing series (" + std::to_string(atom_count) + " atoms).",
              "CoordinateSeries", "importFromFile");
      }
      else {
        atom_count = test_atom_count;
      }
      if (frame_index_start < 0) {
        resize(replica_count + frame_count);
      }
      else {
        resize(replica_count + frame_index_start);
      }
      CoordinateFrame tmp_cf(atom_count, unit_cell);
      CoordinateFrameWriter tmp_cfw = tmp_cf.data();
      const CoordinateFrameReader tmp_cfr(tmp_cfw);
      T* xptr = x_coordinates.data();
      T* yptr = y_coordinates.data();
      T* zptr = z_coordinates.data();
      double* umat_ptr = box_space_transforms.data();
      double* invu_ptr = inverse_transforms.data();
      double* bdim_ptr = box_dimensions.data();
      getAmberInputCoordinates(tf, tmp_cfw.xcrd, tmp_cfw.ycrd, tmp_cfw.zcrd, tmp_cfw.natom,
                               tmp_cfw.umat, tmp_cfw.invu, tmp_cfw.boxdim);
      const int orig_frame_count = frame_count;
      for (int i = 0; i < replica_count; i++) {
        if (frame_index_start < 0) {
          import(tmp_cfr, orig_frame_count + i);
        }
        else {
          import(tmp_cfr, frame_index_start + i);
        }
      }
    }
    break;
  case CoordinateFileKind::SDF:
    break;
  case CoordinateFileKind::AMBER_NETCDF:
  case CoordinateFileKind::AMBER_NETCDF_RST:
    break;
  case CoordinateFileKind::UNKNOWN:
    rtErr("The file type of " + file_name + " could not be understood.", "CoordinateSeries",
          "importFromFile");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::reserve(const int new_frame_capacity) {
  allocate(new_frame_capacity);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::shrinkToFit() {
  const size_t natom_zu = roundUp<size_t>(atom_count, warp_size_zu);
  const size_t fc_zu    = static_cast<size_t>(frame_count);
  const size_t xfrm_zu  = roundUp<size_t>(9, warp_size_zu);
  const size_t bdim_zu  = roundUp<size_t>(6, warp_size_zu);
  x_coordinates.resize(natom_zu * fc_zu);
  y_coordinates.resize(natom_zu * fc_zu);
  z_coordinates.resize(natom_zu * fc_zu);
  box_space_transforms.resize(xfrm_zu * fc_zu);
  inverse_transforms.resize(xfrm_zu * fc_zu);
  box_dimensions.resize(xfrm_zu * fc_zu);
  x_coordinates.shrinkToFit();
  y_coordinates.shrinkToFit();
  z_coordinates.shrinkToFit();
  box_space_transforms.shrinkToFit();
  inverse_transforms.shrinkToFit();
  box_dimensions.shrinkToFit();
  frame_capacity = frame_count;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::resize(const int new_frame_count) {
  allocate(new_frame_count);
  frame_count = new_frame_count;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::resize(const int new_frame_count, const CoordinateFrameReader &cfr,
                                 const int atom_start, const int atom_end) {
  allocate(new_frame_count);
  const int orig_frame_count = frame_count;
  for (int i = orig_frame_count; i < new_frame_count; i++) {
    import(cfr, atom_start, atom_end, i);
  }
  frame_count = new_frame_count;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::resize(const int new_frame_count, const CoordinateFrameWriter &cfw,
                                 const int atom_start, const int atom_end) {
  resize(new_frame_count, CoordinateFrameReader(cfw), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::resize(const int new_frame_count, const CoordinateFrame &cf,
                                 const int atom_start, const int atom_end) {
  resize(new_frame_count, cf.data(), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::resize(const int new_frame_count, CoordinateFrame *cf,
                                 const int atom_start, const int atom_end) {
  resize(new_frame_count, CoordinateFrameReader(cf->data()), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::resize(const int new_frame_count, const PhaseSpace &ps,
                                 const int atom_start, const int atom_end) {
  resize(new_frame_count, CoordinateFrameReader(ps), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::resize(const int new_frame_count, PhaseSpace *ps,
                                 const int atom_start, const int atom_end) {
  resize(new_frame_count, CoordinateFrameReader(CoordinateFrameWriter(ps)), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::pushBack(const CoordinateFrameReader &cfr, const int atom_start,
                                   const int atom_end) {
  const int orig_frame_count = frame_count;
  if (frame_count >= frame_capacity) {
    allocate(((frame_count * 5) + 3) / 4);
  }
  import(cfr, atom_start, atom_end, orig_frame_count);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::pushBack(const CoordinateFrameWriter &cfw, const int atom_start,
                                   const int atom_end) {
  pushBack(CoordinateFrameReader(cfw), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::pushBack(const CoordinateFrame &cf, const int atom_start,
                                   const int atom_end) {
  pushBack(cf.data(), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::pushBack(CoordinateFrame *cf, const int atom_start, const int atom_end) {
  pushBack(CoordinateFrameReader(cf->data()), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::pushBack(const PhaseSpace &ps, const int atom_start,
                                   const int atom_end) {
  pushBack(CoordinateFrameReader(ps), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::pushBack(PhaseSpace *ps, const int atom_start, const int atom_end) {
  pushBack(CoordinateFrameReader(CoordinateFrameWriter(ps)), atom_start, atom_end);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void CoordinateSeries<T>::allocate(const int new_frame_capacity) {
  if (new_frame_capacity > frame_capacity) {
    frame_capacity = new_frame_capacity;
    const size_t fc_zu = static_cast<size_t>(frame_capacity);
    const size_t total_atoms = fc_zu * roundUp(static_cast<size_t>(atom_count), warp_size_zu);
    const size_t total_xfrm = fc_zu * roundUp<size_t>(9, warp_size_zu);
    const size_t total_bdim = fc_zu * roundUp<size_t>(6, warp_size_zu);
    
    // Allocate space for the new capacity.  This will reallocate each array, but one by one in
    // order to avoid keeping what could be nearly double the amount of data at any one time.
    x_coordinates.resize(total_atoms);
    y_coordinates.resize(total_atoms);
    z_coordinates.resize(total_atoms);
    box_space_transforms.resize(total_xfrm);
    inverse_transforms.resize(total_xfrm);
    box_dimensions.resize(total_bdim);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeriesWriter<T> restoreType(CoordinateSeriesWriter<void> *rasa) {
  return CoordinateSeriesWriter<T>(rasa->natom, rasa->nframe, rasa->unit_cell, rasa->gpos_bits,
                                   rasa->gpos_scale, rasa->inv_gpos_scale,
                                   reinterpret_cast<T*>(rasa->xcrd),
                                   reinterpret_cast<T*>(rasa->ycrd),
                                   reinterpret_cast<T*>(rasa->zcrd), rasa->umat, rasa->invu,
                                   rasa->boxdim);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
CoordinateSeriesWriter<T> restoreType(const CoordinateSeriesWriter<void> &rasa) {
  return CoordinateSeriesWriter<T>(rasa.natom, rasa.nframe, rasa.unit_cell, rasa.gpos_bits,
                                   rasa.gpos_scale, rasa.inv_gpos_scale,
                                   reinterpret_cast<T*>(const_cast<void*>(rasa.xcrd)),
                                   reinterpret_cast<T*>(const_cast<void*>(rasa.ycrd)),
                                   reinterpret_cast<T*>(const_cast<void*>(rasa.zcrd)), rasa.umat,
                                   rasa.invu, rasa.boxdim);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const CoordinateSeriesReader<T> restoreType(const CoordinateSeriesReader<void> *rasa) {
  return CoordinateSeriesReader<T>(rasa->natom, rasa->nframe, rasa->unit_cell, rasa->gpos_bits,
                                   rasa->gpos_scale, rasa->inv_gpos_scale,
                                   reinterpret_cast<const T*>(rasa->xcrd),
                                   reinterpret_cast<const T*>(rasa->ycrd),
                                   reinterpret_cast<const T*>(rasa->zcrd), rasa->umat, rasa->invu,
                                   rasa->boxdim);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
const CoordinateSeriesReader<T> restoreType(const CoordinateSeriesReader<void> &rasa) {
  return CoordinateSeriesReader<T>(rasa.natom, rasa.nframe, rasa.unit_cell, rasa.gpos_bits,
                                   rasa.gpos_scale, rasa.inv_gpos_scale,
                                   reinterpret_cast<const T*>(rasa.xcrd),
                                   reinterpret_cast<const T*>(rasa.ycrd),
                                   reinterpret_cast<const T*>(rasa.zcrd), rasa.umat, rasa.invu,
                                   rasa.boxdim);
}

//-------------------------------------------------------------------------------------------------
template <typename Torig, typename Tnew>
CoordinateSeries<Tnew> changeCoordinateSeriesType(const CoordinateSeries<Torig> &cs,
                                                  const int globalpos_scale_bits_in) {
  const CoordinateSeriesReader<Torig> csr = cs.data();
  CoordinateSeries<Tnew> result(csr.natom, csr.nframe, csr.unit_cell, globalpos_scale_bits_in);
  CoordinateSeriesWriter<Tnew> resultw = result.data();
  const size_t natom_zu        = csr.natom;
  const size_t padded_natom_zu = roundUp<size_t>(csr.natom, warp_size_zu);
  const size_t padded_xfrm_zu  = roundUp<size_t>(9, warp_size_zu);
  const size_t padded_bdim_zu  = roundUp<size_t>(6, warp_size_zu);
  const size_t fc_zu           = static_cast<size_t>(csr.nframe);
  if (isFloatingPointScalarType<Torig>() && isFloatingPointScalarType<Tnew>()) {
    for (size_t i = 0; i < fc_zu; i++) {
      for (size_t j = 0; j < natom_zu; j++) {
        const size_t pos = (i * padded_natom_zu) + j;
        resultw.xcrd[pos] = csr.xcrd[pos];
        resultw.ycrd[pos] = csr.ycrd[pos];
        resultw.zcrd[pos] = csr.zcrd[pos];
      }
    }
  }
  else if (isFloatingPointScalarType<Torig>() && isSignedIntegralScalarType<Tnew>()) {
    for (size_t i = 0; i < fc_zu; i++) {
      for (size_t j = 0; j < natom_zu; j++) {
        const size_t pos = (i * padded_natom_zu) + j;
        resultw.xcrd[pos] = llround(csr.xcrd[pos] * resultw.gpos_scale);
        resultw.ycrd[pos] = llround(csr.ycrd[pos] * resultw.gpos_scale);
        resultw.zcrd[pos] = llround(csr.zcrd[pos] * resultw.gpos_scale);
      }
    }
  }
  else if (isFloatingPointScalarType<Tnew>() && isSignedIntegralScalarType<Torig>()) {
    const double conv_factor = pow(2.0, -globalpos_scale_bits_in);
    for (size_t i = 0; i < fc_zu; i++) {
      for (size_t j = 0; j < natom_zu; j++) {
        const size_t pos = (i * padded_natom_zu) + j;
        resultw.xcrd[pos] = static_cast<double>(csr.xcrd[pos]) * csr.inv_gpos_scale;
        resultw.ycrd[pos] = static_cast<double>(csr.ycrd[pos]) * csr.inv_gpos_scale;
        resultw.zcrd[pos] = static_cast<double>(csr.zcrd[pos]) * csr.inv_gpos_scale;
      }
    }
  }
  else if (isSignedIntegralScalarType<Tnew>() && isSignedIntegralScalarType<Torig>()) {
    if (globalpos_scale_bits_in == csr.gpos_bits) {
      for (size_t i = 0; i < fc_zu; i++) {
        for (size_t j = 0; j < natom_zu; j++) {
          const size_t pos = (i * padded_natom_zu) + j;
          resultw.xcrd[pos] = static_cast<Tnew>(csr.xcrd[pos]);
          resultw.ycrd[pos] = static_cast<Tnew>(csr.ycrd[pos]);
          resultw.zcrd[pos] = static_cast<Tnew>(csr.zcrd[pos]);
        }
      }
    }
    else if (globalpos_scale_bits_in > csr.gpos_bits) {
      const int forward_shift = globalpos_scale_bits_in - csr.gpos_bits;
      for (size_t i = 0; i < fc_zu; i++) {
        for (size_t j = 0; j < natom_zu; j++) {
          const size_t pos = (i * padded_natom_zu) + j;
          const llint ixval = csr.xcrd[pos];
          const llint iyval = csr.ycrd[pos];
          const llint izval = csr.zcrd[pos];
          resultw.xcrd[pos] = static_cast<Tnew>(ixval << forward_shift);
          resultw.ycrd[pos] = static_cast<Tnew>(iyval << forward_shift);
          resultw.zcrd[pos] = static_cast<Tnew>(izval << forward_shift);
        }
      }
    }
    else {
      const int backward_shift = globalpos_scale_bits_in - csr.gpos_bits;
      for (size_t i = 0; i < fc_zu; i++) {
        for (size_t j = 0; j < natom_zu; j++) {
          const size_t pos = (i * padded_natom_zu) + j;
          const llint ixval = csr.xcrd[pos];
          const llint iyval = csr.ycrd[pos];
          const llint izval = csr.zcrd[pos];
          resultw.xcrd[pos] = static_cast<Tnew>(ixval >> backward_shift);
          resultw.ycrd[pos] = static_cast<Tnew>(iyval >> backward_shift);
          resultw.zcrd[pos] = static_cast<Tnew>(izval >> backward_shift);
        }
      }
    }
  }
  for (size_t i = 0; i < fc_zu; i++) {
    for (size_t j = 0; j < 9LLU; j++) {
      const size_t pos = (i * padded_xfrm_zu) + j;
      resultw.umat[pos] = csr.umat[pos];
      resultw.invu[pos] = csr.invu[pos];
    }
    for (size_t j = 0; j < 6LLU; j++) {
      const size_t pos = (i * padded_bdim_zu) + j;
      resultw.boxdim[pos] = csr.boxdim[pos];
    }
  }
  
  return result;
}

} // namespace trajectory
} // namespace stormm
