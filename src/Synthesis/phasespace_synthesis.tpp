// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace synthesis {

//-------------------------------------------------------------------------------------------------
template <typename T>
void PhaseSpaceSynthesis::importSystem(const T* x_import, const T* y_import, const T* z_import,
                                       const double* box_xform_in, const double* inverse_xform_in,
                                       const double* box_dimensions_in, const int system_index,
                                       const CoordinateCycle orientation,
                                       const double inverse_scaling_factor,
                                       const TrajectoryKind kind, const HybridTargetLevel tier) {
  checkFormatCompatibility(tier, format, "PhaseSpaceSynthesis", "import");
  llint *x_recv, *y_recv, *z_recv;
  int *x_recv_ovrf, *y_recv_ovrf, *z_recv_ovrf;
  double *box_xform_ptr, *inverse_xform_ptr, *box_dim_ptr;
  double conv_factor;
  bool needs_overflow;
  switch (kind) {
  case TrajectoryKind::POSITIONS:
    conv_factor = inverse_scaling_factor * globalpos_scale;
    needs_overflow = (globalpos_scale_bits > globalpos_scale_nonoverflow_bits);
    switch (orientation) {
    case CoordinateCycle::WHITE:
      x_recv            = x_coordinates.data(tier);
      y_recv            = y_coordinates.data(tier);
      z_recv            = z_coordinates.data(tier);
      x_recv_ovrf       = x_coordinate_overflow.data(tier);
      y_recv_ovrf       = y_coordinate_overflow.data(tier);
      z_recv_ovrf       = z_coordinate_overflow.data(tier);
      box_xform_ptr     = box_space_transforms.data(tier);
      inverse_xform_ptr = inverse_transforms.data(tier);
      box_dim_ptr       = box_dimensions.data(tier);
      break;
    case CoordinateCycle::BLACK:
      x_recv            = x_alt_coordinates.data(tier);
      y_recv            = y_alt_coordinates.data(tier);
      z_recv            = z_alt_coordinates.data(tier);
      x_recv_ovrf       = x_alt_coord_overflow.data(tier);
      y_recv_ovrf       = y_alt_coord_overflow.data(tier);
      z_recv_ovrf       = z_alt_coord_overflow.data(tier);
      box_xform_ptr     = alt_box_transforms.data(tier);
      inverse_xform_ptr = alt_inverse_transforms.data(tier);
      box_dim_ptr       = alt_box_dimensions.data(tier);      
      break;
    }
    break;
  case TrajectoryKind::VELOCITIES:
    conv_factor = inverse_scaling_factor * velocity_scale;
    needs_overflow = (velocity_scale_bits > velocity_scale_nonoverflow_bits);
    switch (orientation) {
    case CoordinateCycle::WHITE:
      x_recv      = x_velocities.data(tier);
      y_recv      = y_velocities.data(tier);
      z_recv      = z_velocities.data(tier);
      x_recv_ovrf = x_velocity_overflow.data(tier);
      y_recv_ovrf = y_velocity_overflow.data(tier);
      z_recv_ovrf = z_velocity_overflow.data(tier);
      break;
    case CoordinateCycle::BLACK:
      x_recv      = x_alt_velocities.data(tier);
      y_recv      = y_alt_velocities.data(tier);
      z_recv      = z_alt_velocities.data(tier);
      x_recv_ovrf = x_alt_velocity_overflow.data(tier);
      y_recv_ovrf = y_alt_velocity_overflow.data(tier);
      z_recv_ovrf = z_alt_velocity_overflow.data(tier);
      break;
    }
    break;
  case TrajectoryKind::FORCES:
    conv_factor = inverse_scaling_factor * force_scale;
    needs_overflow = (force_scale_bits > force_scale_nonoverflow_bits);
    switch (orientation) {
    case CoordinateCycle::WHITE:
      x_recv      = x_forces.data(tier);
      y_recv      = y_forces.data(tier);
      z_recv      = z_forces.data(tier);
      x_recv_ovrf = x_force_overflow.data(tier);
      y_recv_ovrf = y_force_overflow.data(tier);
      z_recv_ovrf = z_force_overflow.data(tier);
      break;
    case CoordinateCycle::BLACK:
      x_recv      = x_alt_forces.data(tier);
      y_recv      = y_alt_forces.data(tier);
      z_recv      = z_alt_forces.data(tier);
      x_recv_ovrf = x_alt_force_overflow.data(tier);
      y_recv_ovrf = y_alt_force_overflow.data(tier);
      z_recv_ovrf = z_alt_force_overflow.data(tier);
      break;
    }
    break;
  }
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const int pos_start   = atom_starts.readHost(system_index);
      const int pos_end     = pos_start + atom_counts.readHost(system_index);
      const int box_offset  = roundUp(9, warp_size_int) * system_index;
      const int dim_offset  = roundUp(6, warp_size_int) * system_index;
      switch (kind) {
      case TrajectoryKind::POSITIONS:
        for (size_t i = 0; i < 9LLU; i++) {
          box_xform_ptr[box_offset + i] = box_xform_in[i];
          inverse_xform_ptr[box_offset + i] = inverse_xform_in[i];
          const int95_t fpbv = hostDoubleToInt95(inverse_xform_in[i] * globalpos_scale);
          box_vectors.putHost(fpbv.x, box_offset + i);
          box_vector_overflow.putHost(fpbv.y, box_offset + i);
        }
        for (size_t i = 0; i < 6LLU; i++) {
          box_dim_ptr[dim_offset + i] = box_dimensions_in[i];
        }
        break;
      case TrajectoryKind::VELOCITIES:
      case TrajectoryKind::FORCES:
        break;
      }
      if (needs_overflow) {
        for (int i = pos_start; i < pos_end; i++) {
          const size_t ip = i - pos_start;
          const int95_t fpx = hostDoubleToInt95(static_cast<double>(x_import[ip]) * conv_factor);
          const int95_t fpy = hostDoubleToInt95(static_cast<double>(y_import[ip]) * conv_factor);
          const int95_t fpz = hostDoubleToInt95(static_cast<double>(z_import[ip]) * conv_factor);
          x_recv[i]      = fpx.x;
          x_recv_ovrf[i] = fpx.y;
          y_recv[i]      = fpy.x;
          y_recv_ovrf[i] = fpy.y;
          z_recv[i]      = fpz.x;
          z_recv_ovrf[i] = fpz.y;
        }
      }
      else {
        for (int i = pos_start; i < pos_end; i++) {
          const size_t ip = i - pos_start;
          const llint fpx = llround(static_cast<double>(x_import[ip]) * conv_factor);
          const llint fpy = llround(static_cast<double>(y_import[ip]) * conv_factor);
          const llint fpz = llround(static_cast<double>(z_import[ip]) * conv_factor);
          x_recv[i]      = fpx;
          y_recv[i]      = fpy;
          z_recv[i]      = fpz;
          x_recv_ovrf[i] = 0;
          y_recv_ovrf[i] = 0;
          z_recv_ovrf[i] = 0;
        }
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void PhaseSpaceSynthesis::importSystem(const CoordinateSeriesReader<T> &csr, const int frame_index,
                                       const int system_index, const CoordinateCycle orientation,
                                       const TrajectoryKind kind, const HybridTargetLevel tier) {
  const size_t xfrm_offset = frame_index * roundUp(9, warp_size_int);
  const size_t bdim_offset = frame_index * roundUp(6, warp_size_int);
  const size_t atom_offset = static_cast<size_t>(frame_index) *
                             static_cast<size_t>(roundUp(csr.natom, warp_size_int));
  importSystem(&csr.xcrd[atom_offset], &csr.ycrd[atom_offset], &csr.zcrd[atom_offset],
         &csr.umat[xfrm_offset], &csr.invu[xfrm_offset], &csr.boxdim[bdim_offset], system_index,
         orientation, csr.inv_gpos_scale, kind, tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void PhaseSpaceSynthesis::importSystem(const CoordinateSeriesReader<T> &csr, const int frame_index,
                                       const int system_index, const TrajectoryKind kind,
                                       const HybridTargetLevel tier) {
  importSystem(csr, frame_index, system_index, cycle_position, kind, tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void PhaseSpaceSynthesis::importSystem(const CoordinateSeriesWriter<T> &csw, const int frame_index,
                                       const int system_index, const CoordinateCycle orientation,
                                       const TrajectoryKind kind, const HybridTargetLevel tier) {
  importSystem(CoordinateSeriesReader<T>(csw), frame_index, system_index, orientation, kind, tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void PhaseSpaceSynthesis::importSystem(const CoordinateSeriesWriter<T> &csw, const int frame_index,
                                       const int system_index, const TrajectoryKind kind,
                                       const HybridTargetLevel tier) {
  importSystem(CoordinateSeriesReader<T>(csw), frame_index, system_index, cycle_position, kind,
               tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void PhaseSpaceSynthesis::importSystem(const CoordinateSeries<T> &cs, const int frame_index,
                                       const int system_index, const CoordinateCycle orientation,
                                       const TrajectoryKind kind, const HybridTargetLevel tier) {
  importSystem(cs.data(), frame_index, system_index, orientation, kind, tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void PhaseSpaceSynthesis::importSystem(const CoordinateSeries<T> &cs, const int frame_index,
                                       const int system_index, const TrajectoryKind kind,
                                       const HybridTargetLevel tier) {
  importSystem(cs.data(), frame_index, system_index, cycle_position, kind, tier);
}

//-------------------------------------------------------------------------------------------------
template <typename T, typename TReader>
void PhaseSpaceSynthesis::mountSystems(const std::vector<T> &crd_list,
                                       const std::vector<AtomGraph*> &ag_list,
                                       const GpuDetails &gpu) {

  // Check that all members in the list of coordinate objects are of the same format.
  const HybridFormat input_format = crd_list[0].getFormat();
  for (int i = 1; i < system_count; i++) {
    if (crd_list[i].getFormat() != input_format) {
      rtErr("The format of all input coordinate objects must be identical.  Coordinates at "
            "index " + std::to_string(i) + " have format " +
            getEnumerationName(crd_list[i].getFormat()) + ", but " +
            getEnumerationName(crd_list[0].getFormat()) + " is needed.", "PhaseSpaceSynthesis");
    }
  }

  // Allocate data and set internal pointers
  int atom_stride = 0;
  for (int i = 0; i < system_count; i++) {
    atom_stride += roundUp(crd_list[i].getAtomCount(), warp_size_int);
  }
  allocate(atom_stride);

  // For a device-only memory layout, it is best to lay out temporary arrays that will hold
  // system-wide descriptors and perform the uploads of each one at a time.
#ifdef STORMM_USE_HPC
  std::vector<int> tmp_atom_counts, tmp_atom_starts, tmp_stib_buffer, tmp_sti_buffer;
  std::vector<int> tmp_replica_buffer, tmp_utr_buffer;
  switch (format) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::HOST_MOUNTED:
    break;
  case HybridFormat::DEVICE_ONLY:
    tmp_atom_counts.resize(system_count);
    tmp_atom_starts.resize(system_count);
    tmp_sti_buffer.resize(system_count);
    tmp_stib_buffer.resize(unique_topology_count + 1, 0);
    tmp_replica_buffer.resize(system_count);
    tmp_utr_buffer.resize(system_count);
    break;
  }
#endif

  // Survey all systems and list all examples using each unique topology.
  int *sti_ptr, *stib_ptr, *replica_ptr, *utr_ptr;
  switch (format) {
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::HOST_MOUNTED:
    sti_ptr  = shared_topology_instances.data();
    stib_ptr = shared_topology_instance_bounds.data();
    replica_ptr = shared_topology_instance_index.data();
    utr_ptr = unique_topology_reference.data();
    break;
  case HybridFormat::DEVICE_ONLY:
    sti_ptr = tmp_sti_buffer.data();
    stib_ptr = tmp_stib_buffer.data();
    replica_ptr = tmp_replica_buffer.data();
    utr_ptr = tmp_utr_buffer.data();
    break;
#else
  case HybridFormat::HOST_ONLY:
    sti_ptr  = shared_topology_instances.data();
    stib_ptr = shared_topology_instance_bounds.data();
    replica_ptr = shared_topology_instance_index.data();
    utr_ptr = unique_topology_reference.data();
    break;
#endif
  }
  for (int i = 0; i < system_count; i++) {
    const AtomGraph* iag_ptr = topologies[i];
    for (int j = 0; j < unique_topology_count; j++) {
      if (iag_ptr == unique_topologies[j]) {
        stib_ptr[j] += 1;
      }
    }
  }
  prefixSumInPlace(stib_ptr, unique_topology_count + 1, PrefixSumType::EXCLUSIVE,
                   "PhaseSpaceSynthesis");
  std::vector<int> stib_counters;
  switch (format) {
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::HOST_MOUNTED:
    stib_counters = shared_topology_instance_bounds.readHost();
    break;
  case HybridFormat::DEVICE_ONLY:
    shared_topology_instance_bounds.putDevice(tmp_stib_buffer);
    stib_counters = tmp_stib_buffer;
    break;
#else
  case HybridFormat::HOST_ONLY:
    stib_counters = shared_topology_instance_bounds.readHost();
    break;
#endif
  }
  if (stib_counters.back() != system_count) {
    rtErr("Counts of systems linked to each unique topology are incorrect.",
          "PhaseSpaceSynthesis");
  }
  for (int i = 0; i < system_count; i++) {
    const AtomGraph* iag_ptr = topologies[i];
    for (int j = 0; j < unique_topology_count; j++) {
      if (iag_ptr == unique_topologies[j]) {
        sti_ptr[stib_counters[j]] = i;
        replica_ptr[i] = stib_counters[j] - stib_ptr[j];
        utr_ptr[i] = j;
        stib_counters[j] += 1;
      }
    }
  }

  // Check that coordinates match topologies.  Set atom starts and counts in the process.
  int acc_limit = 0;
  for (int i = 0; i < system_count; i++) {
    const int natom = crd_list[i].getAtomCount();
    if (natom != ag_list[i]->getAtomCount()) {
      rtErr("Input topology and coordinate sets disagree on atom counts (" +
            std::to_string(ag_list[i]->getAtomCount()) + " vs. " + std::to_string(natom) + ").",
            "PhaseSpaceSynthesis");
    }
    switch (format) {
#ifdef STORMM_USE_HPC
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
    case HybridFormat::UNIFIED:
    case HybridFormat::HOST_ONLY:
    case HybridFormat::HOST_MOUNTED:
      atom_counts.putHost(natom, i);
      atom_starts.putHost(acc_limit, i);
      break;
    case HybridFormat::DEVICE_ONLY:
      tmp_atom_counts[i] = natom;
      tmp_atom_starts[i] = acc_limit;
      break;
#else
    case HybridFormat::HOST_ONLY:
      atom_counts.putHost(natom, i);
      atom_starts.putHost(acc_limit, i);
      break;
#endif
    }
    acc_limit += roundUp(natom, warp_size_int);
  }
#ifdef STORMM_USE_HPC
  switch (format) {
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::HOST_MOUNTED:
    break;
  case HybridFormat::DEVICE_ONLY:
    atom_counts.putDevice(tmp_atom_counts);
    atom_starts.putDevice(tmp_atom_starts);
    shared_topology_instances.putDevice(tmp_sti_buffer);
    shared_topology_instance_bounds.putDevice(tmp_stib_buffer);
    shared_topology_instance_index.putDevice(tmp_replica_buffer);
    unique_topology_reference.putDevice(tmp_utr_buffer);
    break;
  }
#endif
  // Establish the unit cell type
  bool uc_none = false;
  bool uc_orth = false;
  bool uc_tric = false;
  for (int i = 0; i < system_count; i++) {
    switch (crd_list[i].getUnitCellType()) {
    case UnitCellType::NONE:
      uc_none = true;
      break;
    case UnitCellType::ORTHORHOMBIC:
      uc_orth = true;
      break;
    case UnitCellType::TRICLINIC:
      uc_tric = true;
      break;
    }
  }
  if (uc_none) {
    if (uc_orth || uc_tric) {
      rtErr("A coordinate synthesis cannot be formed with a combination of systems having "
            "periodic boundary conditions as well as systems having no boundary conditions.",
            "PhaseSpaceSynthesis");
    }
    unit_cell = UnitCellType::NONE;
  }
  if (uc_orth || uc_tric) {
    unit_cell = (uc_tric) ? UnitCellType::TRICLINIC : UnitCellType::ORTHORHOMBIC;
  }

  // Loop over all systems and import coordinates.  If the input format has memory on the host,
  // take this as authoritative and build the synthesis based on those structures.  Otherwise,
  // take information from the input's memory staged on the HPC device.
  switch (format) {
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
  case HybridFormat::DECOUPLED:
  case HybridFormat::UNIFIED:
  case HybridFormat::HOST_ONLY:
  case HybridFormat::HOST_MOUNTED:
#else
  case HybridFormat::HOST_ONLY:
#endif
    switch (input_format) {
#ifdef STORMM_USE_HPC
    case HybridFormat::EXPEDITED:
    case HybridFormat::DECOUPLED:
    case HybridFormat::UNIFIED:
    case HybridFormat::HOST_ONLY:
    case HybridFormat::HOST_MOUNTED:
#else
    case HybridFormat::HOST_ONLY:
#endif
      // Host-to-host object loading is mediated by C++ code.  All other types of loading, while
      // must less common, require device-to-host temporary copies or a kernel.
      for (int i = 0; i < system_count; i++) {
        loadHostCoordinates(crd_list[i], i);
      }
      break;
#ifdef STORMM_USE_HPC
    case HybridFormat::DEVICE_ONLY:

      // Only in the case of host memory not accessible to the device does a very cumbersome
      // download of the input objects' device memory memory need to occur.  Otherwise, a kernel
      // can handle the loading.
      if (format == HybridFormat::HOST_ONLY) {
        for (int i = 0; i < system_count; i++) {
          T tmp_crd(crd_list[i].getAtomCount(), crd_list[i].getUnitCellType(),
                    HybridFormat::HOST_ONLY);
          const Hybrid<double> *crd_storage = crd_list[i].getStorageHandle();
          deepCopy(tmp_crd.getStorageHandle(), *crd_storage);
          loadHostCoordinates(tmp_crd, i);
        }
      }
      else {

        // A kernel handles communication the communication of device-resident input coordinates
        // to device-accessible memory on the host.
        PsSynthesisWriter poly_psw = this->deviceViewToHostData();
        for (int i = 0; i < system_count; i++) {
          const TReader crdr = crd_list[i].data(HybridTargetLevel::DEVICE);
          loadXPciCoordinates(&poly_psw, i, crdr, gpu);
        }
      }
      break;
#endif
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridFormat::DEVICE_ONLY:
    {
      PsSynthesisWriter poly_psw = this->data(HybridTargetLevel::DEVICE);
      switch (input_format) {
      case HybridFormat::EXPEDITED:
      case HybridFormat::UNIFIED:
      case HybridFormat::HOST_MOUNTED:
        for (int i = 0; i < system_count; i++) {

          // A kernel handles communication between input coordinates held in device-accessible
          // host memory and device-resident object data.
          const TReader crdr = crd_list[i].deviceViewToHostData();
          loadXPciCoordinates(&poly_psw, i, crdr, gpu);
        }
        break;
      case HybridFormat::DECOUPLED:
      case HybridFormat::DEVICE_ONLY:
      case HybridFormat::HOST_ONLY:
        for (int i = 0; i < system_count; i++) {

          // The case of host-resident input memory inaccessible to the device and device-resident
          // memory in the object is the reverse of the tedious process of device-resident input
          // and host-exclusive object memory.  A copy of each PhaseSpace input will be made in a
          // format that the device can see, then uploaded to the device.  This will not place an
          // undue burden on page-locked memory resources to create one system at a time in this
          // manner.
          T tmp_crd(crd_list[i], HybridFormat::HOST_MOUNTED);
          const TReader crdr = tmp_crd.deviceViewToHostData();
          loadXPciCoordinates(&poly_psw, i, crdr, gpu);
        }
        break;
      }
    }
    break;
#endif
  }
}

} // namespace synthesis
} // namespace stormm
