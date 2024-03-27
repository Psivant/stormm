// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
template <typename Tdest, typename Torig>
void copyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest, const Torig* xorig,
                       const Torig* yorig, const Torig* zorig, const int natom) {

  // Check the data types
  if (isFloatingPointScalarType<Tdest>() == false) {
    rtErr("Coordinates of type " + getStormmScalarTypeName<Tdest>() + " cannot be copied without "
          "supplying a conversion scaling factor.", "copyCoordinateXYZ");
  }
  if (isFloatingPointScalarType<Torig>() == false) {
    rtErr("Coordinates of type " + getStormmScalarTypeName<Torig>() + " cannot be filled without "
          "supplying a conversion scaling factor.", "copyCoordinateXYZ");
  }
  for (int i = 0; i < natom; i++) {
    xdest[i] = xorig[i];
    ydest[i] = yorig[i];
    zdest[i] = zorig[i];
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest, typename Torig>
void copyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest, const double dest_scale,
                       const Torig* xorig, const Torig* yorig, const Torig* zorig,
                       const int natom) {

  // The origin data must be of some floating point type
  if (isFloatingPointScalarType<Torig>() == false) {
    rtErr("Coordinates of type " + getStormmScalarTypeName<Torig>() + " cannot be filled without "
          "supplying a conversion scaling factor.", "copyCoordinateXYZ");
  }
  const size_t ct_dest = std::type_index(typeid(Tdest)).hash_code();
  if (ct_dest == llint_type_index) {
    for (int i = 0; i < natom; i++) {
      xdest[i] = llround(xorig[i] * dest_scale);
      ydest[i] = llround(yorig[i] * dest_scale);
      zdest[i] = llround(zorig[i] * dest_scale);
    }
  }
  else if (ct_dest == int_type_index || ct_dest == short_type_index) {
    for (int i = 0; i < natom; i++) {
      xdest[i] = round(xorig[i] * dest_scale);
      ydest[i] = round(yorig[i] * dest_scale);
      zdest[i] = round(zorig[i] * dest_scale);
    }
  }
  else {
    for (int i = 0; i < natom; i++) {
      xdest[i] = xorig[i] * dest_scale;
      ydest[i] = yorig[i] * dest_scale;
      zdest[i] = zorig[i] * dest_scale;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest, typename Torig>
void copyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest, const Torig* xorig,
                       const Torig* yorig, const Torig* zorig, const double orig_scale,
                       const int natom) {
  const double inv_orig_scale = 1.0 / orig_scale;
  for (int i = 0; i < natom; i++) {
    xdest[i] = static_cast<double>(xorig[i]) * inv_orig_scale;
    ydest[i] = static_cast<double>(yorig[i]) * inv_orig_scale;
    zdest[i] = static_cast<double>(zorig[i]) * inv_orig_scale;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest, typename Torig>
void copyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest, const double dest_scale,
                       const Torig* xorig, const Torig* yorig, const Torig* zorig,
                       const double orig_scale, const int natom) {
  const bool dest_is_int = isSignedIntegralScalarType<Tdest>();
  const bool orig_is_int = isSignedIntegralScalarType<Torig>();
  if (dest_is_int && orig_is_int) {
    const int dest_bits = round(log2(dest_scale));
    const int orig_bits = round(log2(orig_scale));
    if (orig_bits < dest_bits) {
      const int shft_bits = dest_bits - orig_bits;
      for (int i = 0; i < natom; i++) {
        const llint xtmp = xorig[i];
        const llint ytmp = yorig[i];
        const llint ztmp = zorig[i];
        xdest[i] = (xtmp << shft_bits);
        ydest[i] = (ytmp << shft_bits);
        zdest[i] = (ztmp << shft_bits);
      }
    }
    else if (orig_bits > dest_bits) {
      const int shft_bits = orig_bits - dest_bits;
      for (int i = 0; i < natom; i++) {
        const llint xtmp = xorig[i];
        const llint ytmp = yorig[i];
        const llint ztmp = zorig[i];
        xdest[i] = (xtmp >> shft_bits);
        ydest[i] = (ytmp >> shft_bits);
        zdest[i] = (ztmp >> shft_bits);
      }
    }
    else {
      for (int i = 0; i < natom; i++) {
        xdest[i] = xorig[i];
        ydest[i] = yorig[i];
        zdest[i] = zorig[i];
      }
    }
  }
  else {
    const double conv_factor = dest_scale / orig_scale;
    if (dest_is_int) {
      for (int i = 0; i < natom; i++) {
        xdest[i] = llround(static_cast<double>(xorig[i]) * conv_factor);
        ydest[i] = llround(static_cast<double>(yorig[i]) * conv_factor);
        zdest[i] = llround(static_cast<double>(zorig[i]) * conv_factor);
      }
    }
    else {
      for (int i = 0; i < natom; i++) {
        xdest[i] = static_cast<double>(xorig[i]) * conv_factor;
        ydest[i] = static_cast<double>(yorig[i]) * conv_factor;
        zdest[i] = static_cast<double>(zorig[i]) * conv_factor;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest>
void copyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest,
                       const llint* xorig, const int* xorig_ovrf, const llint* yorig,
                       const int* yorig_ovrf, const llint* zorig, const	int* zorig_ovrf,
                       const double orig_scale, const int natom) {
  if (isFloatingPointScalarType<Tdest>() == false) {
    rtErr("Coordinates of type " + getStormmScalarTypeName<Tdest>() + " cannot be copied without "
          "supplying a conversion scaling factor.", "copyCoordinateXYZ");
  }
  const double inv_orig_scale = 1.0 / orig_scale;
  for (int i = 0; i < natom; i++) {
    xdest[i] = hostInt95ToDouble(xorig[i], xorig_ovrf[i]) * inv_orig_scale;
    ydest[i] = hostInt95ToDouble(yorig[i], yorig_ovrf[i]) * inv_orig_scale;
    zdest[i] = hostInt95ToDouble(zorig[i], zorig_ovrf[i]) * inv_orig_scale;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest>
void copyCoordinateXYZ(Tdest* xdest, Tdest* ydest, Tdest* zdest, const double dest_scale,
                       const llint* xorig, const int* xorig_ovrf, const llint* yorig,
                       const int* yorig_ovrf, const llint* zorig, const	int* zorig_ovrf,
                       const double orig_scale, const int natom) {
  const size_t ct = std::type_index(typeid(Tdest)).hash_code();
  if (ct == llint_type_index) {

    // Short of the 95-bit integer format, a long long integer representation is the only format
    // that has any possibility of exceeding the information content of double-precision floating
    // point numbers.  Use a full int95_t bit conversion procedure and then take the 64-bit
    // component, which must hold all of the information or the llint format is busted anyway.
    const int dest_bits = round(log2(dest_scale));
    const int orig_bits = round(log2(orig_scale));
    for (int i = 0; i < natom; i++) {
      const int95_t xo = { xorig[i], xorig_ovrf[i] };
      const int95_t yo = { yorig[i], yorig_ovrf[i] };
      const int95_t zo = { zorig[i], zorig_ovrf[i] };
      const int95_t xn = hostChangeFPBits(xo, orig_bits, dest_bits);
      const int95_t yn = hostChangeFPBits(xo, orig_bits, dest_bits);
      const int95_t zn = hostChangeFPBits(xo, orig_bits, dest_bits);
      xdest[i] = xn.x;
      ydest[i] = yn.x;
      zdest[i] = zn.x;
    }
  }
  else {
    const double conv_scale = dest_scale / orig_scale;
    for (int i = 0; i < natom; i++) {
      xdest[i] = hostInt95ToDouble(xorig[i], xorig_ovrf[i]) * conv_scale;
      ydest[i] = hostInt95ToDouble(yorig[i], yorig_ovrf[i]) * conv_scale;
      zdest[i] = hostInt95ToDouble(zorig[i], zorig_ovrf[i]) * conv_scale;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Torig>
void copyCoordinateXYZ(llint* xdest, int* xdest_ovrf, llint* ydest, int* ydest_ovrf, llint* zdest,
                       int* zdest_ovrf, const double dest_scale, const Torig* xorig,
                       const Torig* yorig, const Torig* zorig, const int natom) {
  if (isFloatingPointScalarType<Torig>() == false) {
    rtErr("Coordinates of type " + getStormmScalarTypeName<Torig>() + " cannot be copied without "
          "supplying a conversion scaling factor.", "copyCoordinateXYZ");
  }
  for (int i = 0; i < natom; i++) {
    const int95_t xn = hostDoubleToInt95(xorig[i] * dest_scale);
    xdest[i] = xn.x;
    xdest_ovrf[i] = xn.y;
    const int95_t yn = hostDoubleToInt95(yorig[i] * dest_scale);
    ydest[i] = yn.x;
    ydest_ovrf[i] = yn.y;
    const int95_t zn = hostDoubleToInt95(zorig[i] * dest_scale);
    zdest[i] = zn.x;
    zdest_ovrf[i] = zn.y;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Torig>
void copyCoordinateXYZ(llint* xdest, int* xdest_ovrf, llint* ydest, int* ydest_ovrf, llint* zdest,
                       int* zdest_ovrf, const double dest_scale, const Torig* xorig,
                       const Torig* yorig, const Torig* zorig, const double orig_scale,
                       const int natom) {
  const size_t ct = std::type_index(typeid(Torig)).hash_code();
  if (ct == llint_type_index) {
    const int orig_bits = round(log2(orig_scale));
    const int dest_bits = round(log2(dest_scale));
    
    // Again, handle the long long integer format with a special case as it is the only type that
    // might lose bits in a conversion via double-precision numbers.
    for (int i = 0; i < natom; i++) {
      const int95_t xo = { static_cast<llint>(xorig[i]), 0 };
      const int95_t yo = { static_cast<llint>(yorig[i]), 0 };
      const int95_t zo = { static_cast<llint>(zorig[i]), 0 };
      const int95_t xn = hostChangeFPBits(xo, orig_bits, dest_bits);
      const int95_t yn = hostChangeFPBits(yo, orig_bits, dest_bits);
      const int95_t zn = hostChangeFPBits(zo, orig_bits, dest_bits);
      xdest[i] = xn.x;
      ydest[i] = yn.x;
      zdest[i] = zn.x;
      xdest_ovrf[i] = xn.y;
      ydest_ovrf[i] = yn.y;
      zdest_ovrf[i] = zn.y;
    }
  }
  else {
    const double conv_factor = dest_scale / orig_scale;
    for (int i = 0; i < natom; i++) {
      const int95_t xn = hostDoubleToInt95(xorig[i] * conv_factor);
      xdest[i] = xn.x;
      xdest_ovrf[i] = xn.y;
      const int95_t yn = hostDoubleToInt95(yorig[i] * conv_factor);
      ydest[i] = yn.x;
      ydest_ovrf[i] = yn.y;
      const int95_t zn = hostDoubleToInt95(zorig[i] * conv_factor);
      zdest[i] = zn.x;
      zdest_ovrf[i] = zn.y;
    }
  }
}
  
//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateFrame *destination, const CoordinateSeries<T> &origin,
               const int frame_orig, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopyValidateFrameIndex(frame_orig, origin.getFrameCount());
  coordCopyValidateAtomCounts(destination->getAtomCount(), origin.getAtomCount());
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      origin.extractFrame(destination, frame_orig);

      // Return after completing the host-to-host transfer
      return;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }

  // Any transfers involving device data fall through the switch above to execute the following
  // code.  A rare unrolling of the memory tier enumeration is used to suppress code bloat.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, destination_tier, origin_tier);
  CoordinateFrameWriter destr = (destination_tier == HybridTargetLevel::HOST) ?
                                destination->deviceViewToHostData() :
                                destination->data(destination_tier);
  const CoordinateSeriesReader<T> origr = (origin_tier == HybridTargetLevel::HOST) ?
                                          origin.deviceViewToHostData() : origin.data(origin_tier);
  const int orig_atom_start = roundUp(origin.getAtomCount(), warp_size_int) * frame_orig;
  const int orig_xfrm_start = roundUp(9, warp_size_int) * frame_orig;
  const int orig_bdim_start = roundUp(6, warp_size_int) * frame_orig;
  launchCopyCoordinateXYZBox(destr.xcrd, destr.ycrd, destr.zcrd, 1.0, double_type_index,
                             origr.xcrd, origr.ycrd, origr.zcrd, origr.gpos_scale,
                             std::type_index(typeid(T)).hash_code(), destination->getAtomCount(),
                             destr.umat, destr.invu, destr.boxdim, origr.umat, origr.invu,
                             origr.boxdim, 0, orig_atom_start, 0, orig_xfrm_start,
                             0, orig_bdim_start, gpu);
  launchResolution(sync, destination_tier, origin_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(PhaseSpace *destination, const TrajectoryKind kind,
               const CoordinateCycle orientation, const CoordinateSeries<T> &origin,
               const int frame_orig, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopyValidateFrameIndex(frame_orig, origin.getFrameCount());
  coordCopyValidateAtomCounts(destination->getAtomCount(), origin.getAtomCount());
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      origin.extractFrame(destination, frame_orig, kind, orientation);

      // Return after completing the host-to-host transfer
      return;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }

  // Any transfers involving device data fall through the switch above to execute the following
  // code.  A rare unrolling of the memory tier enumeration is used to suppress code bloat.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, destination_tier, origin_tier);
  PhaseSpaceWriter destr = (destination_tier == HybridTargetLevel::HOST) ?
                            destination->deviceViewToHostData(orientation) :
                            destination->data(orientation, destination_tier);
  const CoordinateSeriesReader<T> origr = (origin_tier == HybridTargetLevel::HOST) ?
                                          origin.deviceViewToHostData() : origin.data(origin_tier);
  const int orig_atom_start = roundUp(origin.getAtomCount(), warp_size_int) * frame_orig;
  const int orig_xfrm_start = roundUp(9, warp_size_int) * frame_orig;
  const int orig_bdim_start = roundUp(6, warp_size_int) * frame_orig;
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    launchCopyCoordinateXYZBox(destr.xcrd, destr.ycrd, destr.zcrd, 1.0, double_type_index,
                               origr.xcrd, origr.ycrd, origr.zcrd, origr.gpos_scale,
                               std::type_index(typeid(T)).hash_code(), destination->getAtomCount(),
                               destr.umat, destr.invu, destr.boxdim, origr.umat, origr.invu,
                               origr.boxdim, 0, orig_atom_start, 0, orig_xfrm_start,
                               0, orig_bdim_start, gpu);
    break;
  case TrajectoryKind::VELOCITIES:
    launchCopyCoordinateXYZ(destr.xvel, destr.yvel, destr.zvel, 1.0, double_type_index,
                            origr.xcrd, origr.ycrd, origr.zcrd, origr.gpos_scale,
                            std::type_index(typeid(T)).hash_code(), destination->getAtomCount(),
                            0, orig_atom_start, gpu);
    break;
  case TrajectoryKind::FORCES:
    launchCopyCoordinateXYZ(destr.xfrc, destr.yfrc, destr.zfrc, 1.0, double_type_index,
                            origr.xcrd, origr.ycrd, origr.zcrd, origr.gpos_scale,
                            std::type_index(typeid(T)).hash_code(), destination->getAtomCount(),
                            0, orig_atom_start, gpu);
    break;
  }
  launchResolution(sync, destination_tier, origin_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(PhaseSpace *destination, const TrajectoryKind kind,
               const CoordinateSeries<T> &origin, const int frame_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy<T>(destination, kind, destination->getCyclePosition(), origin, frame_orig,
               destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(PhaseSpace *destination, const CoordinateSeries<T> &origin,
               const int frame_orig, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopy<T>(destination, TrajectoryKind::POSITIONS, destination->getCyclePosition(), origin,
               frame_orig, destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeries<T> *destination, const int frame_dest,
               const CoordinateFrameReader &origin, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopyValidateFrameIndex(frame_dest, destination->getFrameCount());
  coordCopyValidateAtomCounts(destination->getAtomCount(), origin.natom);
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      destination->import(origin, frame_dest);
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }

  // The host-to-host import function involves checks on the size of the coordinate series as
  // well as the atom count.  Any procedures involving data on the device will fall through the
  // switch to execute the following code.  Place explicit checks here.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, destination_tier, origin_tier);
  const int dest_atom_start = roundUp(destination->getAtomCount(), warp_size_int) * frame_dest;
  const int dest_xfrm_start = roundUp(9, warp_size_int) * frame_dest;
  const int dest_bdim_start = roundUp(6, warp_size_int) * frame_dest;
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    {
      CoordinateSeriesWriter<T> csw = destination->deviceViewToHostData();
      launchCopyCoordinateXYZBox(csw.xcrd, csw.ycrd, csw.zcrd, csw.gpos_scale,
                                 std::type_index(typeid(T)).hash_code(), origin.xcrd, origin.ycrd,
                                 origin.zcrd, 1.0, double_type_index, origin.natom, csw.umat,
                                 csw.invu, csw.boxdim, origin.umat, origin.invu, origin.boxdim, 
                                 dest_atom_start, 0, dest_xfrm_start, 0, dest_bdim_start, 0, gpu);
    }
    break;
  case HybridTargetLevel::DEVICE:
    {
      CoordinateSeriesWriter<T> csw = destination->data(destination_tier);
      launchCopyCoordinateXYZBox(csw.xcrd, csw.ycrd, csw.zcrd, csw.gpos_scale,
                                 std::type_index(typeid(T)).hash_code(), origin.xcrd, origin.ycrd,
                                 origin.zcrd, 1.0, double_type_index, origin.natom, csw.umat,
                                 csw.invu, csw.boxdim, origin.umat, origin.invu, origin.boxdim, 
                                 dest_atom_start, 0, dest_xfrm_start, 0, dest_bdim_start, 0, gpu);
    }
    break;
  }
  launchResolution(sync, destination_tier, origin_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeries<T> *destination, const int frame_dest,
               const CoordinateFrame &origin, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {

  // As this routine merely delegates to another with checks on the frame count and number of
  // atoms in each system, no checks are needed here.
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    coordCopy(destination, frame_dest, origin.data(origin_tier), destination_tier, origin_tier,
              gpu);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      coordCopy(destination, frame_dest, origin.deviceViewToHostData(), destination_tier,
                origin_tier, gpu, sync);
      break;
    case HybridTargetLevel::DEVICE:      
      coordCopy(destination, frame_dest, origin.data(origin_tier), destination_tier, origin_tier,
                gpu);
      break;
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeries<T> *destination, const int frame_dest,
               const PhaseSpaceReader &origin, const TrajectoryKind kind,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopyValidateFrameIndex(frame_dest, destination->getFrameCount());
  coordCopyValidateAtomCounts(destination->getAtomCount(), origin.natom);
  const int dest_atom_start = roundUp(destination->getAtomCount(), warp_size_int) * frame_dest;
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        CoordinateSeriesWriter<T> csw = destination->data(destination_tier);
        switch (kind) {
        case TrajectoryKind::POSITIONS:
          copyBoxInformation(csw.boxdim, csw.umat, csw.invu, frame_dest, origin.boxdim,
                             origin.umat, origin.invu);
          copyCoordinateXYZ<T, double>(&csw.xcrd[dest_atom_start], &csw.ycrd[dest_atom_start],
                                       &csw.zcrd[dest_atom_start], csw.gpos_scale, origin.xcrd,
                                       origin.ycrd, origin.zcrd, 1.0, origin.natom);
          break;
        case TrajectoryKind::VELOCITIES:
          copyCoordinateXYZ<T, double>(&csw.xcrd[dest_atom_start], &csw.ycrd[dest_atom_start],
                                       &csw.zcrd[dest_atom_start], csw.gpos_scale, origin.xvel,
                                       origin.yvel, origin.zvel, 1.0, origin.natom);
          break;
        case TrajectoryKind::FORCES:          
          copyCoordinateXYZ<T, double>(&csw.xcrd[dest_atom_start], &csw.ycrd[dest_atom_start],
                                       &csw.zcrd[dest_atom_start], csw.gpos_scale, origin.xfrc,
                                       origin.yfrc, origin.zfrc, 1.0, origin.natom);
          break;
        }
      }

      // Return immediately after completing the host-to-host copy
      return;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }

  // Any copy procedures involving data on the GPU device will fall through the switch above to
  // execute the following code.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, destination_tier, origin_tier);
  CoordinateSeriesWriter<T> csw = (destination_tier == HybridTargetLevel::HOST) ?
                                  destination->deviceViewToHostData() :
                                  destination->data(destination_tier);
  const size_t ct_dest = std::type_index(typeid(T)).hash_code();
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    {
      const int dest_xfrm_start = roundUp(9, warp_size_int) * frame_dest;
      const int dest_bdim_start = roundUp(6, warp_size_int) * frame_dest;
      launchCopyCoordinateXYZBox(csw.xcrd, csw.ycrd, csw.zcrd, csw.gpos_scale, ct_dest,
                                 origin.xcrd, origin.ycrd, origin.zcrd, 1.0, double_type_index,
                                 origin.natom, csw.umat, csw.invu, csw.boxdim, origin.umat,
                                 origin.invu, origin.boxdim, dest_atom_start, 0, dest_xfrm_start,
                                 0, dest_bdim_start, 0, gpu);
    }
    break;
  case TrajectoryKind::VELOCITIES:
    launchCopyCoordinateXYZ(csw.xcrd, csw.ycrd, csw.zcrd, csw.gpos_scale, ct_dest, origin.xvel,
                            origin.yvel, origin.zvel, 1.0, double_type_index, origin.natom,
                            dest_atom_start, 0, gpu);
    break;
  case TrajectoryKind::FORCES:
    launchCopyCoordinateXYZ(csw.xcrd, csw.ycrd, csw.zcrd, csw.gpos_scale, ct_dest, origin.xfrc,
                            origin.yfrc, origin.zfrc, 1.0, double_type_index, origin.natom,
                            dest_atom_start, 0, gpu);
    break;
  }
  launchResolution(sync, destination_tier, origin_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeries<T> *destination, const int frame_dest,
               const PhaseSpaceReader &origin, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, frame_dest, origin, TrajectoryKind::POSITIONS, destination_tier,
            origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeries<T> *destination, const int frame_dest, const PhaseSpace &origin,
               const TrajectoryKind kind, const CoordinateCycle orientation,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    coordCopy(destination, frame_dest, origin.data(orientation, origin_tier), kind,
              destination_tier, origin_tier, gpu, sync);
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      coordCopy(destination, frame_dest, origin.deviceViewToHostData(orientation), kind,
                destination_tier, origin_tier, gpu, sync);
      break;
    case HybridTargetLevel::DEVICE:
      coordCopy(destination, frame_dest, origin.data(orientation, origin_tier), kind,
                destination_tier, origin_tier, gpu, sync);
      break;
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeries<T> *destination, const int frame_dest, const PhaseSpace &origin,
               const TrajectoryKind kind, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopy(destination, frame_dest, origin, kind, origin.getCyclePosition(), destination_tier,
            origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeries<T> *destination, const int frame_dest, const PhaseSpace &origin,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, frame_dest, origin, TrajectoryKind::POSITIONS, origin.getCyclePosition(),
            destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest, typename Torig>
void coordCopy(CoordinateSeriesWriter<Tdest> *destination, const size_t frame_dest,
               const CoordinateSeriesReader<Torig> &origin, const size_t frame_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {

  // The coordinate series abstract can still recover its own starting atom indices.
  coordCopyValidateFrameIndex(frame_dest, destination->nframe);
  coordCopyValidateFrameIndex(frame_orig, origin.nframe);
  coordCopyValidateAtomCounts(destination->natom, origin.natom);
  const size_t dest_atom_start = roundUp<size_t>(destination->natom, warp_size_zu) * frame_dest;
  const size_t orig_atom_start = roundUp<size_t>(origin.natom, warp_size_zu) * frame_orig;
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      copyBoxInformation(destination->boxdim, destination->umat, destination->invu, frame_dest,
                         origin.boxdim, origin.umat, origin.invu, frame_orig);
      copyCoordinateXYZ<Tdest, Torig>(&destination->xcrd[dest_atom_start],
                                      &destination->ycrd[dest_atom_start],
                                      &destination->zcrd[dest_atom_start], destination->gpos_scale,
                                      &origin.xcrd[orig_atom_start], &origin.ycrd[orig_atom_start],
                                      &origin.zcrd[orig_atom_start], origin.gpos_scale,
                                      origin.natom);
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }

  // Any copy operations involving data on the GPU device will fall through the above switch and
  // execute the following code.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, destination_tier, origin_tier);
  const int xfrm_w = roundUp(9, warp_size_int);
  const int bdim_w = roundUp(6, warp_size_int);
  const int dest_xfrm_start = xfrm_w * frame_dest;
  const int dest_bdim_start = bdim_w * frame_dest;
  const int orig_xfrm_start = xfrm_w * frame_orig;
  const int orig_bdim_start = bdim_w * frame_orig;
  launchCopyCoordinateXYZBox(destination->xcrd, destination->ycrd, destination->zcrd,
                             destination->gpos_scale, std::type_index(typeid(Tdest)).hash_code(),
                             origin.xcrd, origin.ycrd, origin.zcrd, origin.gpos_scale,
                             std::type_index(typeid(Torig)).hash_code(), origin.natom,
                             destination->umat, destination->invu, destination->boxdim,
                             origin.umat, origin.invu, origin.boxdim, dest_atom_start,
                             orig_atom_start, dest_xfrm_start, orig_xfrm_start, dest_bdim_start,
                             orig_bdim_start, gpu);
  launchResolution(sync, destination_tier, origin_tier);
#endif
}
  
//-------------------------------------------------------------------------------------------------
template <typename Tdest, typename Torig>
void coordCopy(CoordinateSeries<Tdest> *destination, const size_t frame_dest,
               const CoordinateSeries<Torig> &origin, const size_t frame_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        CoordinateSeriesWriter<Tdest> csw = destination->data(destination_tier);
        coordCopy(&csw, frame_dest, origin.data(origin_tier), frame_orig, destination_tier,
                  origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        CoordinateSeriesWriter<Tdest> csw = destination->deviceViewToHostData();
        coordCopy(&csw, frame_dest, origin.data(origin_tier), frame_orig, destination_tier,
                  origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CoordinateSeriesWriter<Tdest> csw = destination->data(destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&csw, frame_dest, origin.deviceViewToHostData(), frame_orig, destination_tier,
                  origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&csw, frame_dest, origin.data(origin_tier), frame_orig, destination_tier,
                  origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeriesWriter<T> *destination, const size_t frame_dest,
               const PsSynthesisReader &origin, const int orig_atom_start, const int index_orig,
               const TrajectoryKind kind, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopyValidateFrameIndex(frame_dest, destination->nframe);
  coordCopyValidateSystemIndex(index_orig, origin.system_count);
  const int natom = destination->natom;
  const size_t dest_atom_start = roundUp<size_t>(destination->natom, warp_size_zu) * frame_dest;
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      switch (kind) {
      case TrajectoryKind::POSITIONS:
        copyBoxInformation(destination->boxdim, destination->umat, destination->invu, frame_dest,
                           origin.boxdims, origin.umat, origin.invu, index_orig);
        copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->ycrd[dest_atom_start],
                          &destination->zcrd[dest_atom_start], destination->gpos_scale,
                          &origin.xcrd[orig_atom_start], &origin.xcrd_ovrf[orig_atom_start],
                          &origin.ycrd[orig_atom_start], &origin.ycrd_ovrf[orig_atom_start],
                          &origin.zcrd[orig_atom_start], &origin.zcrd_ovrf[orig_atom_start],
                          origin.gpos_scale, natom);
        break;
      case TrajectoryKind::VELOCITIES:
        copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->ycrd[dest_atom_start],
                          &destination->zcrd[dest_atom_start], destination->gpos_scale,
                          &origin.xvel[orig_atom_start], &origin.xvel_ovrf[orig_atom_start],
                          &origin.yvel[orig_atom_start], &origin.yvel_ovrf[orig_atom_start],
                          &origin.zvel[orig_atom_start], &origin.zvel_ovrf[orig_atom_start],
                          origin.vel_scale, natom);
        break;
      case TrajectoryKind::FORCES:
        copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->ycrd[dest_atom_start],
                          &destination->zcrd[dest_atom_start], destination->gpos_scale,
                          &origin.xfrc[orig_atom_start], &origin.xfrc_ovrf[orig_atom_start],
                          &origin.yfrc[orig_atom_start], &origin.yfrc_ovrf[orig_atom_start],
                          &origin.zfrc[orig_atom_start], &origin.zfrc_ovrf[orig_atom_start],
                          origin.frc_scale, natom);
        break;
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }

  // Coordinate copy procedures involving data on the GPU device will fall through the switch
  // above and execute the following code.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, destination_tier, origin_tier);
  const size_t ct_dest = std::type_index(typeid(T)).hash_code();
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    {
      const int xfrm_w = roundUp(9, warp_size_int);
      const int bdim_w = roundUp(6, warp_size_int);
      const int dest_xfrm_start = xfrm_w * frame_dest;
      const int dest_bdim_start = bdim_w * frame_dest;
      const int orig_xfrm_start = xfrm_w * index_orig;
      const int orig_bdim_start = bdim_w * index_orig;
      launchCopyCoordinateXYZBox(destination->xcrd, destination->ycrd, destination->zcrd,
                                 destination->gpos_scale, destination->gpos_bits, ct_dest,
                                 origin.xcrd, origin.ycrd, origin.zcrd, origin.xcrd_ovrf,
                                 origin.ycrd_ovrf, origin.zcrd_ovrf, origin.gpos_scale,
                                 origin.gpos_bits, natom, destination->umat, destination->invu,
                                 destination->boxdim, origin.umat, origin.invu, origin.boxdims,
                                 dest_atom_start, orig_atom_start, dest_xfrm_start,
                                 orig_xfrm_start, dest_bdim_start, orig_bdim_start, gpu);
    }
    break;
  case TrajectoryKind::VELOCITIES:
    launchCopyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd,
                            destination->gpos_scale, destination->gpos_bits, ct_dest,
                            origin.xvel, origin.yvel, origin.zvel, origin.xvel_ovrf,
                            origin.yvel_ovrf, origin.zvel_ovrf, origin.vel_scale, origin.vel_bits,
                            natom, dest_atom_start, orig_atom_start, gpu);
    break;
  case TrajectoryKind::FORCES:
    launchCopyCoordinateXYZ(destination->xcrd, destination->ycrd, destination->zcrd,
                            destination->gpos_scale, destination->gpos_bits, ct_dest,
                            origin.xfrc, origin.yfrc, origin.zfrc, origin.xfrc_ovrf,
                            origin.yfrc_ovrf, origin.zfrc_ovrf, origin.frc_scale, origin.frc_bits,
                            natom, dest_atom_start, orig_atom_start, gpu);
    break;
  }
  launchResolution(sync, destination_tier, origin_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeriesWriter<T> *destination, const size_t frame_dest,
               const PsSynthesisReader &origin, const int orig_atom_start, const int index_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, frame_dest, origin, orig_atom_start, index_orig,
            TrajectoryKind::POSITIONS, destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeries<T> *destination, const size_t frame_dest,
               const PhaseSpaceSynthesis &origin, const int index_orig, const TrajectoryKind kind,
               const CoordinateCycle orientation, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopyValidateFrameIndex(frame_dest, destination->getFrameCount());
  coordCopyValidateSystemIndex(index_orig, origin.getSystemCount());
  coordCopyValidateAtomCounts(destination->getAtomCount(), origin.getAtomCount(index_orig));
  const int orig_atom_start = origin.getAtomOffset(index_orig);
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        CoordinateSeriesWriter<T> csw = destination->data(destination_tier);
        coordCopy(&csw, frame_dest, origin.data(orientation, origin_tier), orig_atom_start,
                  index_orig, kind, destination_tier, origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        CoordinateSeriesWriter<T> csw = destination->deviceViewToHostData();
        coordCopy(&csw, frame_dest, origin.data(orientation, origin_tier), orig_atom_start,
                  index_orig, kind, destination_tier, origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CoordinateSeriesWriter<T> csw = destination->data(destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&csw, frame_dest, origin.deviceViewToHostData(orientation), orig_atom_start,
                  index_orig, kind, destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&csw, frame_dest, origin.data(orientation, origin_tier), orig_atom_start,
                  index_orig, kind, destination_tier, origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeries<T> *destination, const size_t frame_dest,
               const PhaseSpaceSynthesis &origin, const int index_orig, TrajectoryKind kind,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, frame_dest, origin, index_orig, kind, origin.getCyclePosition(),
            destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeries<T> *destination, const size_t frame_dest,
               const PhaseSpaceSynthesis &origin, const int index_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, frame_dest, origin, index_orig, TrajectoryKind::POSITIONS,
            origin.getCyclePosition(), destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeriesWriter<T> *destination, size_t frame_dest,
               const CondensateReader &origin, size_t orig_atom_start, size_t index_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_orig, origin.system_count);
  coordCopyValidateFrameIndex(frame_dest, destination->nframe);
  const int natom = destination->natom;
  const size_t dest_atom_start = frame_dest * roundUp<size_t>(natom, warp_size_zu);
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      copyBoxInformation(destination->boxdim, destination->umat, destination->invu, frame_dest,
                         origin.boxdims, origin.umat, origin.invu, index_orig);
      switch (origin.mode) {
      case PrecisionModel::DOUBLE:
        copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->ycrd[dest_atom_start],
                          &destination->zcrd[dest_atom_start], destination->gpos_scale,
                          &origin.xcrd[orig_atom_start], &origin.ycrd[orig_atom_start],
                          &origin.zcrd[orig_atom_start], natom);
        break;
      case PrecisionModel::SINGLE:
        copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->ycrd[dest_atom_start],
                          &destination->zcrd[dest_atom_start], destination->gpos_scale,
                          &origin.xcrd_sp[orig_atom_start], &origin.ycrd_sp[orig_atom_start],
                          &origin.zcrd_sp[orig_atom_start], natom);
        break;
      }

      // Return after completing the host-to-host copy procedure
      return;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }

  // Any copy procedures involving coordinates on the device will fall through the switch above to
  // execute the following code.  It is assumed that the Condensate contains positions data, such
  // that box information is relevant.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, destination_tier, origin_tier);
  const int xfrm_w = roundUp(9, warp_size_int);
  const int bdim_w = roundUp(6, warp_size_int);
  const int dest_xfrm_start = xfrm_w * frame_dest;
  const int dest_bdim_start = bdim_w * frame_dest;
  const int orig_xfrm_start = xfrm_w * index_orig;
  const int orig_bdim_start = bdim_w * index_orig;
  switch (origin.mode) {
  case PrecisionModel::DOUBLE:
    launchCopyCoordinateXYZBox(destination->xcrd, destination->ycrd, destination->zcrd,
                               destination->gpos_scale, std::type_index(typeid(T)).hash_code(),
                               origin.xcrd, origin.ycrd, origin.zcrd, 1.0, double_type_index,
                               natom, destination->umat, destination->invu, destination->boxdim,
                               origin.umat, origin.invu, origin.boxdims, dest_atom_start,
                               orig_atom_start, dest_xfrm_start, orig_xfrm_start, dest_bdim_start,
                               orig_bdim_start, gpu);
    break;
  case PrecisionModel::SINGLE:
    launchCopyCoordinateXYZBox(destination->xcrd, destination->ycrd, destination->zcrd,
                               destination->gpos_scale, std::type_index(typeid(T)).hash_code(),
                               origin.xcrd_sp, origin.ycrd_sp, origin.zcrd_sp, 1.0,
                               float_type_index, natom, destination->umat, destination->invu,
                               destination->boxdim, origin.umat, origin.invu, origin.boxdims,
                               dest_atom_start, orig_atom_start, dest_xfrm_start, orig_xfrm_start,
                               dest_bdim_start, orig_bdim_start, gpu);
    break;
  }
  launchResolution(sync, destination_tier, origin_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CoordinateSeries<T> *destination, const int frame_dest,
               const Condensate &origin, const int index_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  const int natom = destination->getAtomCount();
  coordCopyValidateAtomCounts(natom, origin.getAtomCount(index_orig));
  const size_t orig_atom_start = origin.getAtomOffset(index_orig);
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        CoordinateSeriesWriter<T> csw = destination->data(destination_tier);
        coordCopy(&csw, frame_dest, origin.data(origin_tier), orig_atom_start, index_orig,
                  destination_tier, origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        CoordinateSeriesWriter<T> csw = destination->data(destination_tier);
        coordCopy(&csw, frame_dest, origin.data(origin_tier), orig_atom_start, index_orig,
                  destination_tier, origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CoordinateSeriesWriter<T> csw = destination->data(destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&csw, frame_dest, origin.deviceViewToHostData(), orig_atom_start, index_orig,
                  destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&csw, frame_dest, origin.data(origin_tier), orig_atom_start, index_orig,
                  destination_tier, origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(PsSynthesisWriter *destination, const int dest_atom_start, const int index_dest,
               const TrajectoryKind kind, const CoordinateSeriesReader<T> &origin,
               const int frame_orig, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_dest, destination->system_count);
  coordCopyValidateFrameIndex(frame_orig, origin.nframe);
  const size_t orig_atom_start = roundUp<size_t>(origin.natom, warp_size_zu) *
                                 static_cast<size_t>(frame_orig);
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      switch (kind) {
      case TrajectoryKind::POSITIONS:
        {
          copyBoxInformation(destination->boxdims, destination->umat, destination->invu,
                             index_dest, origin.boxdim, origin.umat, origin.invu, frame_orig,
                             destination->boxvecs, destination->boxvec_ovrf);
          copyCoordinateXYZ(&destination->xcrd[dest_atom_start],
                            &destination->xcrd_ovrf[dest_atom_start],
                            &destination->ycrd[dest_atom_start],
                            &destination->ycrd_ovrf[dest_atom_start],
                            &destination->zcrd[dest_atom_start],
                            &destination->zcrd_ovrf[dest_atom_start], destination->gpos_scale,
                            &origin.xcrd[orig_atom_start], &origin.ycrd[orig_atom_start],
                            &origin.zcrd[orig_atom_start], origin.gpos_scale, origin.natom);
        }
        break;
      case TrajectoryKind::VELOCITIES:
        copyCoordinateXYZ(&destination->xvel[dest_atom_start],
                          &destination->xvel_ovrf[dest_atom_start],
                          &destination->yvel[dest_atom_start],
                          &destination->yvel_ovrf[dest_atom_start],
                          &destination->zvel[dest_atom_start],
                          &destination->zvel_ovrf[dest_atom_start],
                          destination->vel_scale, &origin.xcrd[orig_atom_start],
                          &origin.ycrd[orig_atom_start], &origin.zcrd[orig_atom_start], 
                          origin.gpos_scale, origin.natom);
        break;
      case TrajectoryKind::FORCES:
        copyCoordinateXYZ(&destination->xfrc[dest_atom_start],
                          &destination->xfrc_ovrf[dest_atom_start],
                          &destination->yfrc[dest_atom_start],
                          &destination->yfrc_ovrf[dest_atom_start],
                          &destination->zfrc[dest_atom_start],
                          &destination->zfrc_ovrf[dest_atom_start],
                          destination->frc_scale, &origin.xcrd[orig_atom_start],
                          &origin.ycrd[orig_atom_start], &origin.zcrd[orig_atom_start],
                          origin.gpos_scale, origin.natom);
        break;
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }

  // Copy procedures involving coordinate data on the GPU device will fall through the switch above
  // to execute the following code.  The switch over trajectory kinds needs to be repeated.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, destination_tier, origin_tier);
  const size_t ct_orig = std::type_index(typeid(T)).hash_code();
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    {
      const int xfrm_w = roundUp(9, warp_size_int);
      const int bdim_w = roundUp(9, warp_size_int);
      const int dest_xfrm_start = xfrm_w * index_dest;
      const int dest_bdim_start = bdim_w * index_dest;
      const int orig_xfrm_start = xfrm_w * frame_orig;
      const int orig_bdim_start = bdim_w * frame_orig;
      launchCopyCoordinateXYZBox(destination->xcrd, destination->ycrd, destination->zcrd,
                                 destination->xcrd_ovrf, destination->ycrd_ovrf,
                                 destination->zcrd_ovrf, destination->gpos_scale,
                                 destination->gpos_bits, origin.xcrd, origin.ycrd, origin.zcrd,
                                 origin.gpos_scale, origin.gpos_bits, ct_orig, origin.natom,
                                 destination->umat, destination->invu, destination->boxdims,
                                 origin.umat, origin.invu, origin.boxdim, destination->boxvecs,
                                 destination->boxvec_ovrf, dest_atom_start, orig_atom_start,
                                 dest_xfrm_start, orig_xfrm_start, dest_bdim_start,
                                 orig_bdim_start, gpu);
    }
    break;
  case TrajectoryKind::VELOCITIES:
    launchCopyCoordinateXYZ(destination->xvel, destination->yvel, destination->zvel,
                            destination->xvel_ovrf, destination->yvel_ovrf, destination->zvel_ovrf,
                            destination->vel_scale, destination->vel_bits, origin.xcrd,
                            origin.ycrd, origin.zcrd, origin.gpos_scale, origin.gpos_bits, ct_orig,
                            origin.natom, dest_atom_start, orig_atom_start, gpu);
    break;
  case TrajectoryKind::FORCES:
    launchCopyCoordinateXYZ(destination->xfrc, destination->yfrc, destination->zfrc,
                            destination->xfrc_ovrf, destination->yfrc_ovrf, destination->zfrc_ovrf,
                            destination->frc_scale, destination->frc_bits, origin.xcrd,
                            origin.ycrd, origin.zcrd, origin.gpos_scale, origin.gpos_bits, ct_orig,
                            origin.natom, dest_atom_start, orig_atom_start, gpu);
    break;
  }
  launchResolution(sync, destination_tier, origin_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(PsSynthesisWriter *destination, const int index_dest,
               const CoordinateSeriesReader<T> &origin, const int frame_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, index_dest, TrajectoryKind::POSITIONS, CoordinateCycle::WHITE, origin,
            frame_orig, destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest, const TrajectoryKind kind,
               const CoordinateCycle orientation, const CoordinateSeries<T> &origin,
               const int frame_orig, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_dest, destination->getSystemCount());
  coordCopyValidateFrameIndex(frame_orig, origin.getFrameCount());
  coordCopyValidateAtomCounts(destination->getAtomCount(index_dest), origin.getAtomCount());
  const int dest_atom_start = destination->getAtomOffset(index_dest);
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        PsSynthesisWriter poly_psw = destination->data(orientation, destination_tier);
        coordCopy(&poly_psw, dest_atom_start, index_dest, kind, origin.data(origin_tier),
                  frame_orig, destination_tier, origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        PsSynthesisWriter poly_psw = destination->deviceViewToHostData(orientation);
        coordCopy(&poly_psw, dest_atom_start, index_dest, kind, origin.data(origin_tier),
                  frame_orig, destination_tier, origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      PsSynthesisWriter poly_psw = destination->data(orientation, destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&poly_psw, dest_atom_start, index_dest, kind, origin.deviceViewToHostData(),
                  frame_orig, destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&poly_psw, dest_atom_start, index_dest, kind, origin.data(origin_tier),
                  frame_orig, destination_tier, origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest, const TrajectoryKind kind,
               const CoordinateSeries<T> &origin, const int frame_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, index_dest, kind, destination->getCyclePosition(), origin, frame_orig,
            destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(PhaseSpaceSynthesis *destination, const int index_dest,
               const CoordinateSeries<T> &origin, const int frame_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, index_dest, TrajectoryKind::POSITIONS, destination->getCyclePosition(),
            origin, frame_orig, destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(CondensateWriter *destination, const size_t dest_atom_start, const int index_dest,
               const CoordinateSeriesReader<T> &origin, const size_t frame_orig,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_dest, destination->system_count);
  coordCopyValidateFrameIndex(frame_orig, origin.nframe);
  const size_t orig_atom_start = roundUp<size_t>(origin.natom, warp_size_zu) * frame_orig;
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      copyBoxInformation(destination->boxdims, destination->umat, destination->invu, index_dest,
                         origin.boxdim, origin.umat, origin.invu);
      switch (destination->mode) {
      case PrecisionModel::DOUBLE:
        copyCoordinateXYZ(&destination->xcrd[dest_atom_start], &destination->ycrd[dest_atom_start],
                          &destination->zcrd[dest_atom_start], &origin.xcrd[orig_atom_start],
                          &origin.ycrd[orig_atom_start], &origin.zcrd[orig_atom_start],
                          origin.gpos_scale, origin.natom);
        break;
      case PrecisionModel::SINGLE:
        copyCoordinateXYZ(&destination->xcrd_sp[dest_atom_start],
                          &destination->ycrd_sp[dest_atom_start],
                          &destination->zcrd_sp[dest_atom_start], &origin.xcrd[orig_atom_start],
                          &origin.ycrd[orig_atom_start], &origin.zcrd[orig_atom_start],
                          origin.gpos_scale, origin.natom);
        break;
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }

  // Copy procedures involving coordinates on the GPU device will fall through the switch above to
  // execute the following code.
#ifdef STORMM_USE_HPC
  launchPreparation(sync, destination_tier, origin_tier);
  const int xfrm_w = roundUp(9, warp_size_int);
  const int bdim_w = roundUp(9, warp_size_int);
  const int dest_xfrm_start = xfrm_w * index_dest;
  const int dest_bdim_start = bdim_w * index_dest;
  const int orig_xfrm_start = xfrm_w * frame_orig;
  const int orig_bdim_start = bdim_w * frame_orig;
  switch (destination->mode) {
  case PrecisionModel::DOUBLE:
    launchCopyCoordinateXYZBox(destination->xcrd, destination->ycrd, destination->zcrd, 1.0,
                               double_type_index, origin.xcrd, origin.ycrd, origin.zcrd,
                               origin.gpos_scale, std::type_index(typeid(T)).hash_code(),
                               origin.natom, destination->umat, destination->invu,
                               destination->boxdims, origin.umat, origin.invu, origin.boxdim,
                               dest_atom_start, orig_atom_start, dest_xfrm_start, orig_xfrm_start,
                               dest_bdim_start, orig_bdim_start, gpu);
    break;
  case PrecisionModel::SINGLE:
    launchCopyCoordinateXYZBox(destination->xcrd_sp, destination->ycrd_sp, destination->zcrd_sp,
                               1.0, float_type_index, origin.xcrd, origin.ycrd, origin.zcrd,
                               origin.gpos_scale, std::type_index(typeid(T)).hash_code(),
                               origin.natom, destination->umat, destination->invu,
                               destination->boxdims, origin.umat, origin.invu, origin.boxdim,
                               dest_atom_start, orig_atom_start, dest_xfrm_start, orig_xfrm_start,
                               dest_bdim_start, orig_bdim_start, gpu);
    break;
  }
  launchResolution(sync, destination_tier, origin_tier);
#endif
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void coordCopy(Condensate *destination, const int index_dest, const CoordinateSeries<T> &origin,
               const size_t frame_orig, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopyValidateSystemIndex(index_dest, destination->getSystemCount());
  coordCopyValidateFrameIndex(frame_orig, origin.getFrameCount());
  coordCopyValidateAtomCounts(destination->getAtomCount(index_dest), origin.getAtomCount());
  const size_t dest_atom_start = destination->getAtomOffset(index_dest);
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        CondensateWriter cdnsw = destination->data(destination_tier);
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.data(origin_tier), frame_orig,
                  destination_tier, origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        CondensateWriter cdnsw = destination->deviceViewToHostData();
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.data(origin_tier), frame_orig,
                  destination_tier, origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      CondensateWriter cdnsw = destination->data(destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.deviceViewToHostData(), frame_orig,
                  destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&cdnsw, dest_atom_start, index_dest, origin.data(origin_tier), frame_orig,
                  destination_tier, origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest, typename Torig>
void coordCopy(CoordinateSeries<Tdest> *destination, const CoordinateSeries<Torig> &origin,
               const std::vector<int2> &system_pairs, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  const Hybrid<int2> hsys_tmp(system_pairs, "xfer_temporary", HybridFormat::HOST_ONLY);
  coordCopy(destination, origin, hsys_tmp, destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest, typename Torig>
void coordCopy(CoordinateSeries<Tdest> *destination, const CoordinateSeries<Torig> &origin,
               const Hybrid<int2> &system_pairs, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  const size_t ct_dest = std::type_index(typeid(Tdest)).hash_code();
  const size_t ct_orig = std::type_index(typeid(Torig)).hash_code();
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        CoordinateSeriesWriter<void> destw = destination->templateFreeData(destination_tier);
        coordCopy(&destw, ct_dest, origin.templateFreeData(origin_tier), ct_orig,
                  system_pairs.data(), system_pairs.size(), destination_tier, origin_tier, gpu,
                  sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        CoordinateSeriesWriter<void> destw = destination->deviceViewToTemplateFreeHostData();
        coordCopy(&destw, ct_dest, origin.templateFreeData(origin_tier), ct_orig,
                  system_pairs.getDeviceValidHostPointer(), system_pairs.size(), destination_tier,
                  origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const int2* sys_ptr = system_pairs.getDeviceValidHostPointer();
      const int npairs = system_pairs.size();
      CoordinateSeriesWriter<void> destw = destination->data(destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&destw, ct_dest, origin.deviceViewToTemplateFreeHostData(), ct_orig, sys_ptr,
                  npairs, destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&destw, ct_dest, origin.data(origin_tier), ct_orig, sys_ptr, npairs,
                  destination_tier, origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest>
void coordCopy(CoordinateSeries<Tdest> *destination, const PhaseSpaceSynthesis &origin,
               TrajectoryKind kind, CoordinateCycle orientation,
               const std::vector<int2> &system_pairs, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  const Hybrid<int2> hsys_tmp(system_pairs, "xfer_temporary", HybridFormat::HOST_ONLY);
  coordCopy(destination, origin, kind, orientation, hsys_tmp, destination_tier, origin_tier, gpu,
            sync);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest>
void coordCopy(CoordinateSeries<Tdest> *destination, const PhaseSpaceSynthesis &origin,
               TrajectoryKind kind, const std::vector<int2> &system_pairs,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, origin, kind, destination->getCyclePosition(), system_pairs,
            destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest>
void coordCopy(CoordinateSeries<Tdest> *destination, const PhaseSpaceSynthesis &origin,
               const std::vector<int2> &system_pairs, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopy(destination, origin, TrajectoryKind::POSITIONS, destination->getCyclePosition(),
            system_pairs, destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest>
void coordCopy(CoordinateSeries<Tdest> *destination, const PhaseSpaceSynthesis &origin,
               const TrajectoryKind kind, const CoordinateCycle orientation,
               const Hybrid<int2> &system_pairs, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  const size_t ct_dest = std::type_index(typeid(Tdest)).hash_code();
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        CoordinateSeriesWriter<void> destw = destination->templateFreeData(destination_tier);
        coordCopy(&destw, ct_dest, origin.data(orientation, origin_tier), kind,
                  system_pairs.data(), system_pairs.size(), destination_tier, origin_tier, gpu,
                  sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        CoordinateSeriesWriter<void> destw = destination->deviceViewToTemplateFreeHostData();
        coordCopy(&destw, ct_dest, origin.data(orientation, origin_tier), kind,
                  system_pairs.getDeviceValidHostPointer(), system_pairs.size(), destination_tier,
                  origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const int2* sys_ptr = system_pairs.getDeviceValidHostPointer();
      const int npairs = system_pairs.size();
      CoordinateSeriesWriter<void> destw = destination->data(destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&destw, ct_dest, origin.deviceViewToHostData(orientation), kind, sys_ptr,
                  npairs, destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&destw, ct_dest, origin.data(orientation, origin_tier), kind, sys_ptr, npairs,
                  destination_tier, origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest>
void coordCopy(CoordinateSeries<Tdest> *destination, const PhaseSpaceSynthesis &origin,
               const TrajectoryKind kind, const Hybrid<int2> &system_pairs,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, origin, kind, destination->getCyclePosition(), system_pairs,
            destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest>
void coordCopy(CoordinateSeries<Tdest> *destination, const PhaseSpaceSynthesis &origin,
               const Hybrid<int2> &system_pairs, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopy(destination, origin, TrajectoryKind::POSITIONS, destination->getCyclePosition(),
            system_pairs, destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest>
void coordCopy(CoordinateSeries<Tdest> *destination, const Condensate &origin,
               const std::vector<int2> &system_pairs, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  const Hybrid<int2> hsys_tmp(system_pairs, "xfer_temporary", HybridFormat::HOST_ONLY);
  coordCopy(destination, origin, hsys_tmp, destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest>
void coordCopy(CoordinateSeries<Tdest> *destination, const Condensate &origin,
               const Hybrid<int2> &system_pairs, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  const size_t ct_dest = std::type_index(typeid(Tdest)).hash_code();
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        CoordinateSeriesWriter<void> destw = destination->templateFreeData(destination_tier);
        coordCopy(&destw, ct_dest, origin.data(origin_tier), system_pairs.data(),
                  system_pairs.size(), destination_tier, origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        CoordinateSeriesWriter<void> destw = destination->deviceViewToTemplateFreeHostData();
        coordCopy(&destw, ct_dest, origin.data(origin_tier),
                  system_pairs.getDeviceValidHostPointer(), system_pairs.size(), destination_tier,
                  origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const int2* sys_ptr = system_pairs.getDeviceValidHostPointer();
      const int npairs = system_pairs.size();
      CoordinateSeriesWriter<void> destw = destination->data(destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&destw, ct_dest, origin.deviceViewToHostData(), sys_ptr, npairs,
                  destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&destw, ct_dest, origin.data(origin_tier), sys_ptr, npairs, destination_tier,
                  origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Torig>
void coordCopy(PhaseSpaceSynthesis *destination, const TrajectoryKind kind,
               const CoordinateCycle orientation, const CoordinateSeries<Torig> &origin,
               const std::vector<int2> &system_pairs, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  const Hybrid<int2> hsys_tmp(system_pairs, "xfer_temporary", HybridFormat::HOST_ONLY);
  coordCopy(destination, kind, orientation, origin, hsys_tmp, destination_tier, origin_tier, gpu,
            sync);
}

//-------------------------------------------------------------------------------------------------
template <typename Torig>
void coordCopy(PhaseSpaceSynthesis *destination, const TrajectoryKind kind,
               const CoordinateSeries<Torig> &origin, const std::vector<int2> &system_pairs,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, kind, destination->getCyclePosition(), origin, system_pairs,
            destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename Torig>
void coordCopy(PhaseSpaceSynthesis *destination, const CoordinateSeries<Torig> &origin,
               const std::vector<int2> &system_pairs, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopy(destination, TrajectoryKind::POSITIONS, destination->getCyclePosition(), origin,
            system_pairs, destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename Torig>
void coordCopy(PhaseSpaceSynthesis *destination, const TrajectoryKind kind,
               const CoordinateCycle orientation, const CoordinateSeries<Torig> &origin,
               const Hybrid<int2> &system_pairs, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  const size_t ct_orig = std::type_index(typeid(Torig)).hash_code();
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        PsSynthesisWriter destw = destination->data(orientation, destination_tier);
        coordCopy(&destw, kind, origin.templateFreeData(origin_tier), ct_orig, system_pairs.data(),
                  system_pairs.size(), destination_tier, origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        PsSynthesisWriter destw = destination->deviceViewToHostData(orientation);
        coordCopy(&destw, kind, origin.templateFreeData(origin_tier), ct_orig,
                  system_pairs.getDeviceValidHostPointer(), system_pairs.size(), destination_tier,
                  origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const int2* sys_ptr = system_pairs.getDeviceValidHostPointer();
      const int npairs = system_pairs.size();
      PsSynthesisWriter destw = destination->data(orientation, destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&destw, kind, origin.deviceViewToTemplateFreeHostData(), ct_orig, sys_ptr,
                  npairs, destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&destw, kind, origin.templateFreeData(origin_tier), ct_orig, sys_ptr, npairs,
                  destination_tier, origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Torig>
void coordCopy(PhaseSpaceSynthesis *destination, const TrajectoryKind kind,
               const CoordinateSeries<Torig> &origin, const Hybrid<int2> &system_pairs,
               const HybridTargetLevel destination_tier, const HybridTargetLevel origin_tier,
               const GpuDetails &gpu, const HpcKernelSync sync) {
  coordCopy(destination, kind, destination->getCyclePosition(), origin, system_pairs,
            destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename Torig>
void coordCopy(PhaseSpaceSynthesis *destination, const CoordinateSeries<Torig> &origin,
               const Hybrid<int2> &system_pairs, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  coordCopy(destination, TrajectoryKind::POSITIONS, destination->getCyclePosition(), origin,
            system_pairs, destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename Torig>
void coordCopy(Condensate *destination, const CoordinateSeries<Torig> &origin,
               const std::vector<int2> &system_pairs, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  const Hybrid<int2> hsys_tmp(system_pairs, "xfer_temporary", HybridFormat::HOST_ONLY);
  coordCopy(destination, origin, hsys_tmp, destination_tier, origin_tier, gpu, sync);
}

//-------------------------------------------------------------------------------------------------
template <typename Torig>
void coordCopy(Condensate *destination, const CoordinateSeries<Torig> &origin,
               const Hybrid<int2> &system_pairs, const HybridTargetLevel destination_tier,
               const HybridTargetLevel origin_tier, const GpuDetails &gpu,
               const HpcKernelSync sync) {
  const size_t ct_orig = std::type_index(typeid(Torig)).hash_code();
  switch (destination_tier) {
  case HybridTargetLevel::HOST:
    switch (origin_tier) {
    case HybridTargetLevel::HOST:
      {
        CondensateWriter destw = destination->data(destination_tier);
        coordCopy(&destw, origin.templateFreeData(origin_tier), ct_orig, system_pairs.data(),
                  system_pairs.size(), destination_tier, origin_tier, gpu, sync);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridTargetLevel::DEVICE:
      {
        CondensateWriter destw = destination->deviceViewToHostData();
        coordCopy(&destw, origin.templateFreeData(origin_tier), ct_orig,
                  system_pairs.getDeviceValidHostPointer(), system_pairs.size(), destination_tier,
                  origin_tier, gpu, sync);
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const int2* sys_ptr = system_pairs.getDeviceValidHostPointer();
      const int npairs = system_pairs.size();
      CondensateWriter destw = destination->data(destination_tier);
      switch (origin_tier) {
      case HybridTargetLevel::HOST:
        coordCopy(&destw, origin.deviceViewToTemplateFreeHostData(), ct_orig, sys_ptr,
                  npairs, destination_tier, origin_tier, gpu, sync);
        break;
      case HybridTargetLevel::DEVICE:
        coordCopy(&destw, origin.templateFreeData(origin_tier), ct_orig, sys_ptr, npairs,
                  destination_tier, origin_tier, gpu, sync);
        break;
      }
    }
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tdest>
void unrollCCXYZOrigin(CoordinateSeriesWriter<Tdest> *destination,
                       const CoordinateSeriesReader<void> &origin, const size_t ct_orig,
                       const int2* system_pairs, const int copy_count,
                       const HybridTargetLevel destination_tier,
                       const HybridTargetLevel origin_tier) {

  // Check that the memory tiers are both set to the CPU host
  if (destination_tier != HybridTargetLevel::HOST || origin_tier != HybridTargetLevel::HOST) {
    rtErr("This variant of the function can only be called for operations on the CPU host.  "
          "Coordinate copying between multiple systems is handled by a separate, overloaded "
          "variant, but it is the responsibility of the program developer to use the proper form "
          "of this function.", "unrollCCXYZOrigin");
  }
  if (ct_orig == double_type_index) {
    const CoordinateSeriesReader<double> tr_orig = restoreType<double>(origin);
    for (int i = 0; i < copy_count; i++) {
      coordCopy<Tdest, double>(destination, system_pairs[i].y, tr_orig, system_pairs[i].x);
    }
  }
  else if (ct_orig == float_type_index) {
    const CoordinateSeriesReader<float> tr_orig = restoreType<float>(origin);
    for (int i = 0; i < copy_count; i++) {
      coordCopy<Tdest, float>(destination, system_pairs[i].y, tr_orig, system_pairs[i].x);
    }
  }
  else if (ct_orig == short_type_index) {
    const CoordinateSeriesReader<short int> tr_orig = restoreType<short int>(origin);
    for (int i = 0; i < copy_count; i++) {
      coordCopy<Tdest, short int>(destination, system_pairs[i].y, tr_orig, system_pairs[i].x);
    }
  }
  else if (ct_orig == int_type_index) {
    const CoordinateSeriesReader<int> tr_orig = restoreType<int>(origin);
    for (int i = 0; i < copy_count; i++) {
      coordCopy<Tdest, int>(destination, system_pairs[i].y, tr_orig, system_pairs[i].x);
    }
  }
  else if (ct_orig == llint_type_index) {
    const CoordinateSeriesReader<llint> tr_orig = restoreType<llint>(origin);
    for (int i = 0; i < copy_count; i++) {
      coordCopy<Tdest, llint>(destination, system_pairs[i].y, tr_orig, system_pairs[i].x);
    }
  }
}

} // namespace trajectory
} // namespace stormm
