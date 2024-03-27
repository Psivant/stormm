#include "copyright.h"
#include "Constants/fixed_precision.h"
#include "condensate.h"

namespace stormm {
namespace synthesis {

using card::HybridKind;
using numerics::max_llint_accumulation;
using numerics::globalpos_scale_nonoverflow_bits;

//-------------------------------------------------------------------------------------------------
CondensateBorders::CondensateBorders(const int system_count_in, const size_t* atom_starts_in,
                                     const int* atom_counts_in) :
    system_count{system_count_in}, atom_starts{atom_starts_in}, atom_counts{atom_counts_in}
{}

//-------------------------------------------------------------------------------------------------
CondensateWriter::CondensateWriter(const PrecisionModel mode_in, const StructureSource basis_in,
                                   const int system_count_in, const UnitCellType unit_cell_in,
                                   const size_t* atom_starts_in, const int* atom_counts_in,
                                   float* xcrd_sp_in, float* ycrd_sp_in, float* zcrd_sp_in,
                                   double* xcrd_in, double* ycrd_in, double* zcrd_in,
                                   double* umat_in, double* invu_in, double* boxdims_in) :
    mode{mode_in}, basis{basis_in}, system_count{system_count_in}, unit_cell{unit_cell_in},
    atom_starts{atom_starts_in}, atom_counts{atom_counts_in}, xcrd_sp{xcrd_sp_in},
    ycrd_sp{ycrd_sp_in}, zcrd_sp{zcrd_sp_in}, xcrd{xcrd_in}, ycrd{ycrd_in}, zcrd{zcrd_in},
    umat{umat_in}, invu{invu_in}, boxdims{boxdims_in}
{}

//-------------------------------------------------------------------------------------------------
CondensateReader::CondensateReader(const PrecisionModel mode_in, const StructureSource basis_in,
                                   const int system_count_in, const UnitCellType unit_cell_in,
                                   const size_t* atom_starts_in, const int* atom_counts_in,
                                   const float* xcrd_sp_in, const float* ycrd_sp_in,
                                   const float* zcrd_sp_in, const double* xcrd_in,
                                   const double* ycrd_in, const double* zcrd_in,
                                   const double* umat_in, const double* invu_in,
                                   const double* boxdims_in) :
    mode{mode_in}, basis{basis_in}, system_count{system_count_in}, unit_cell{unit_cell_in},
    atom_starts{atom_starts_in}, atom_counts{atom_counts_in}, xcrd_sp{xcrd_sp_in},
    ycrd_sp{ycrd_sp_in}, zcrd_sp{zcrd_sp_in}, xcrd{xcrd_in}, ycrd{ycrd_in}, zcrd{zcrd_in},
    umat{umat_in}, invu{invu_in}, boxdims{boxdims_in}
{}

//-------------------------------------------------------------------------------------------------
CondensateReader::CondensateReader(const CondensateWriter &cdw) :
    mode{cdw.mode}, basis{cdw.basis}, system_count{cdw.system_count}, unit_cell{cdw.unit_cell},
    atom_starts{cdw.atom_starts}, atom_counts{cdw.atom_counts}, xcrd_sp{cdw.xcrd_sp},
    ycrd_sp{cdw.ycrd_sp}, zcrd_sp{cdw.zcrd_sp}, xcrd{cdw.xcrd}, ycrd{cdw.ycrd}, zcrd{cdw.zcrd},
    umat{cdw.umat}, invu{cdw.invu}, boxdims{cdw.boxdims}
{}

//-------------------------------------------------------------------------------------------------
Condensate::Condensate() :
    mode{PrecisionModel::SINGLE}, basis{StructureSource::NONE},
    system_count{0}, unit_cell{UnitCellType::NONE},
    holds_own_data{false}, csptr_data_type{0},
    atom_starts{HybridKind::ARRAY, "cdns_atom_starts"},
    atom_counts{HybridKind::ARRAY, "cdns_atom_counts"},
    x_coordinates_sp{HybridKind::POINTER, "cdns_xcrd_sp"},
    y_coordinates_sp{HybridKind::POINTER, "cdns_ycrd_sp"},
    z_coordinates_sp{HybridKind::POINTER, "cdns_zcrd_sp"},
    x_coordinates{HybridKind::POINTER, "cdns_xcrd"},
    y_coordinates{HybridKind::POINTER, "cdns_ycrd"},
    z_coordinates{HybridKind::POINTER, "cdns_zcrd"},
    box_transforms{HybridKind::POINTER, "cdns_umat"},
    inv_transforms{HybridKind::POINTER, "cdns_invu"},
    box_dimensions{HybridKind::POINTER, "cdns_bdim"},
    pps_ptr{nullptr},
    cs_ptr{nullptr},
    float_data{HybridKind::ARRAY, "cdns_floats"},
    double_data{HybridKind::ARRAY, "cdns_doubles"}
{}

//-------------------------------------------------------------------------------------------------
Condensate::Condensate(const PhaseSpaceSynthesis *poly_ps_in, const PrecisionModel mode_in,
                       const GpuDetails &gpu) :
    Condensate()
{
  mode = mode_in;
  pps_ptr = const_cast<PhaseSpaceSynthesis*>(poly_ps_in);
  rebuild(pps_ptr, mode, gpu);
}

//-------------------------------------------------------------------------------------------------
Condensate::Condensate(const Condensate &original) :
    mode{original.mode},
    basis{original.basis},
    system_count{original.system_count},
    unit_cell{original.unit_cell},
    holds_own_data{original.holds_own_data},
    csptr_data_type{original.csptr_data_type},
    atom_starts{original.atom_starts},
    atom_counts{original.atom_counts},
    x_coordinates_sp{original.x_coordinates_sp},
    y_coordinates_sp{original.y_coordinates_sp},
    z_coordinates_sp{original.z_coordinates_sp},
    x_coordinates{original.x_coordinates},
    y_coordinates{original.y_coordinates},
    z_coordinates{original.z_coordinates},
    box_transforms{original.box_transforms},
    inv_transforms{original.inv_transforms},
    box_dimensions{original.box_dimensions},
    pps_ptr{original.pps_ptr},
    cs_ptr{original.cs_ptr},
    float_data{original.float_data},
    double_data{original.double_data}
{
  repairPointers();
}

//-------------------------------------------------------------------------------------------------
Condensate::Condensate(Condensate &&original) :
    mode{original.mode},
    basis{original.basis},
    system_count{original.system_count},
    unit_cell{original.unit_cell},
    holds_own_data{original.holds_own_data},
    csptr_data_type{original.csptr_data_type},
    atom_starts{std::move(original.atom_starts)},
    atom_counts{std::move(original.atom_counts)},
    x_coordinates_sp{std::move(original.x_coordinates_sp)},
    y_coordinates_sp{std::move(original.y_coordinates_sp)},
    z_coordinates_sp{std::move(original.z_coordinates_sp)},
    x_coordinates{std::move(original.x_coordinates)},
    y_coordinates{std::move(original.y_coordinates)},
    z_coordinates{std::move(original.z_coordinates)},
    box_transforms{std::move(original.box_transforms)},
    inv_transforms{std::move(original.inv_transforms)},
    box_dimensions{std::move(original.box_dimensions)},
    pps_ptr{original.pps_ptr},
    cs_ptr{original.cs_ptr},
    float_data{std::move(original.float_data)},
    double_data{std::move(original.double_data)}
{}

//-------------------------------------------------------------------------------------------------
Condensate& Condensate::operator=(const Condensate &other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }

  // Typical copying
  mode = other.mode;
  basis = other.basis;
  system_count = other.system_count;
  unit_cell = other.unit_cell;
  holds_own_data = other.holds_own_data;
  csptr_data_type = other.csptr_data_type;
  atom_starts = other.atom_starts;
  atom_counts = other.atom_counts;
  x_coordinates_sp = other.x_coordinates_sp;
  y_coordinates_sp = other.y_coordinates_sp;
  z_coordinates_sp = other.z_coordinates_sp;
  x_coordinates = other.x_coordinates;
  y_coordinates = other.y_coordinates;
  z_coordinates = other.z_coordinates;
  box_transforms = other.box_transforms;
  inv_transforms = other.inv_transforms;
  box_dimensions = other.box_dimensions;
  pps_ptr = other.pps_ptr;
  cs_ptr = other.cs_ptr;
  float_data = other.float_data;
  double_data = other.double_data;

  // Repair pointers and return the result
  repairPointers();
  return *this;
}

//-------------------------------------------------------------------------------------------------
Condensate& Condensate::operator=(Condensate &&other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }

  // Typical copying
  mode = other.mode;
  basis = other.basis;
  system_count = other.system_count;
  unit_cell = other.unit_cell;
  holds_own_data = other.holds_own_data;
  csptr_data_type = other.csptr_data_type;
  atom_starts = std::move(other.atom_starts);
  atom_counts = std::move(other.atom_counts);
  x_coordinates_sp = std::move(other.x_coordinates_sp);
  y_coordinates_sp = std::move(other.y_coordinates_sp);
  z_coordinates_sp = std::move(other.z_coordinates_sp);
  x_coordinates = std::move(other.x_coordinates);
  y_coordinates = std::move(other.y_coordinates);
  z_coordinates = std::move(other.z_coordinates);
  box_transforms = std::move(other.box_transforms);
  inv_transforms = std::move(other.inv_transforms);
  box_dimensions = std::move(other.box_dimensions);
  pps_ptr = other.pps_ptr;
  cs_ptr = other.cs_ptr;
  float_data = std::move(other.float_data);
  double_data = std::move(other.double_data);
  return *this;
}

//-------------------------------------------------------------------------------------------------
Condensate::Condensate(const PhaseSpaceSynthesis &poly_ps_in, const PrecisionModel mode_in,
                       const GpuDetails &gpu) :
    Condensate(poly_ps_in.getSelfPointer(), mode_in, gpu)
{}

//-------------------------------------------------------------------------------------------------
PrecisionModel Condensate::getMode() const {
  return mode;
}

//-------------------------------------------------------------------------------------------------
StructureSource Condensate::getBasis() const {
  return basis;
}

//-------------------------------------------------------------------------------------------------
int Condensate::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
bool Condensate::ownsCoordinates() const {
  return holds_own_data;
}

//-------------------------------------------------------------------------------------------------
size_t Condensate::getAtomOffset(const int system_index) const {
  return atom_starts.readHost(system_index);
}

//-------------------------------------------------------------------------------------------------
int Condensate::getAtomCount(const int system_index) const {
  return atom_counts.readHost(system_index);
}

//-------------------------------------------------------------------------------------------------
size_t Condensate::getCoordinateSeriesTypeID() const {
  return csptr_data_type;
}

//-------------------------------------------------------------------------------------------------
const PhaseSpaceSynthesis* Condensate::getSynthesisPointer() const {
  return pps_ptr;
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame Condensate::exportCoordinateFrame(const int system_index,
                                                  const HybridTargetLevel tier) const {
  validateSystemIndex(system_index);
  const size_t natom = atom_counts.readHost(system_index);
  const size_t atom_offset = atom_starts.readHost(system_index);
  const size_t xfrm_offset = static_cast<size_t>(system_index) * roundUp<size_t>(9, warp_size_zu);
  const size_t bdim_offset = static_cast<size_t>(system_index) * roundUp<size_t>(6, warp_size_zu);
  CoordinateFrame result(natom);
  CoordinateFrameWriter resw = result.data();
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      switch (mode) {
      case PrecisionModel::DOUBLE:
        {
          const double* x_ptr = x_coordinates.data();
          const double* y_ptr = y_coordinates.data();
          const double* z_ptr = z_coordinates.data();
          for (size_t i = 0; i < natom; i++) {
            resw.xcrd[i] = x_ptr[atom_offset + i];
            resw.ycrd[i] = y_ptr[atom_offset + i];
            resw.zcrd[i] = z_ptr[atom_offset + i];
          }
        }
        break;
      case PrecisionModel::SINGLE:
        {
          const float* x_ptr = x_coordinates_sp.data();
          const float* y_ptr = y_coordinates_sp.data();
          const float* z_ptr = z_coordinates_sp.data();
          for (size_t i = 0; i < natom; i++) {
            resw.xcrd[i] = x_ptr[atom_offset + i];
            resw.ycrd[i] = y_ptr[atom_offset + i];
            resw.zcrd[i] = z_ptr[atom_offset + i];
          }
        }
        break;
      }

      // Copy the transformation matrices.
      const double* fr_umat = box_transforms.data();
      const double* fr_invu = inv_transforms.data();
      const double* fr_bdim = box_dimensions.data();
      for (size_t i = 0; i < 9; i++) {
        resw.umat[i] = fr_umat[xfrm_offset + i];
        resw.invu[i] = fr_invu[xfrm_offset + i];
      }
      for (size_t i = 0; i < 6; i++) {
        resw.boxdim[i] = fr_bdim[bdim_offset + i];
      }
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      switch (mode) {
      case PrecisionModel::DOUBLE:
        {
          const std::vector<double> x_data = x_coordinates.readDevice(atom_offset, natom);
          const std::vector<double> y_data = y_coordinates.readDevice(atom_offset, natom);
          const std::vector<double> z_data = z_coordinates.readDevice(atom_offset, natom);
          for (size_t i = 0; i < natom; i++) {
            resw.xcrd[i] = x_data[i];
            resw.ycrd[i] = y_data[i];
            resw.zcrd[i] = z_data[i];
          }
        }
        break;
      case PrecisionModel::SINGLE:
        {
          const std::vector<float> x_data = x_coordinates_sp.readDevice(atom_offset, natom);
          const std::vector<float> y_data = y_coordinates_sp.readDevice(atom_offset, natom);
          const std::vector<float> z_data = z_coordinates_sp.readDevice(atom_offset, natom);
          for (size_t i = 0; i < natom; i++) {
            resw.xcrd[i] = x_data[i];
            resw.ycrd[i] = y_data[i];
            resw.zcrd[i] = z_data[i];
          }
        }
        break;
      }

      // Copy the transformation matrices.
      const std::vector<double> fr_umat = box_transforms.readDevice(xfrm_offset, 9);
      const std::vector<double> fr_invu = inv_transforms.readDevice(xfrm_offset, 9);
      const std::vector<double> fr_bdim = box_dimensions.readDevice(bdim_offset, 6);
      for (int i = 0; i < 9; i++) {
        resw.umat[i] = fr_umat[i];
        resw.invu[i] = fr_invu[i];
      }
      for (size_t i = 0; i < 6; i++) {
        resw.boxdim[i] = fr_bdim[i];
      }
    }
    break;
#endif
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> Condensate::getInterlacedCoordinates(const int system_index,
                                                         const HybridTargetLevel tier) const {
  validateSystemIndex(system_index);
  const size_t llim = atom_starts.readHost(system_index);
  const size_t hlim = llim + atom_counts.readHost(system_index);
  std::vector<double> result(3LLU * (hlim - llim));
  switch (tier) {
  case HybridTargetLevel::HOST:
    switch (mode) {
    case PrecisionModel::DOUBLE:
      {
        const double* xptr = x_coordinates.data();
        const double* yptr = y_coordinates.data();
        const double* zptr = z_coordinates.data();
        for (size_t i = llim; i < hlim; i++) {
          result[(3 * (i - llim))    ] = xptr[i];
          result[(3 * (i - llim)) + 1] = yptr[i];
          result[(3 * (i - llim)) + 2] = zptr[i];
        }
      }
      break;
    case PrecisionModel::SINGLE:
      {
        const float* xptr = x_coordinates_sp.data();
        const float* yptr = y_coordinates_sp.data();
        const float* zptr = z_coordinates_sp.data();
        for (size_t i = llim; i < hlim; i++) {
          result[(3 * (i - llim))    ] = xptr[i];
          result[(3 * (i - llim)) + 1] = yptr[i];
          result[(3 * (i - llim)) + 2] = zptr[i];
        }
      }
      break;
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    {
      const size_t natom = hlim - llim;
      switch (mode) {
      case PrecisionModel::DOUBLE:
        {
          const std::vector<double> x_copy = x_coordinates.readDevice(llim, natom);
          const std::vector<double> y_copy = y_coordinates.readDevice(llim, natom);
          const std::vector<double> z_copy = z_coordinates.readDevice(llim, natom);
          for (size_t i = 0; i < natom; i++) {
            result[(3 * i)    ] = x_copy[i];
            result[(3 * i) + 1] = y_copy[i];
            result[(3 * i) + 2] = z_copy[i];
          }
        }
        break;
      case PrecisionModel::SINGLE:
        {
          const std::vector<float> x_copy = x_coordinates_sp.readDevice(llim, natom);
          const std::vector<float> y_copy = y_coordinates_sp.readDevice(llim, natom);
          const std::vector<float> z_copy = z_coordinates_sp.readDevice(llim, natom);
          for (size_t i = 0; i < natom; i++) {
            result[(3 * i)    ] = x_copy[i];
            result[(3 * i) + 1] = y_copy[i];
            result[(3 * i) + 2] = z_copy[i];
          }
        }
        break;
      }
    }
    break;
#endif
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const Condensate* Condensate::getSelfPointer() const {
  return this;
}

//-------------------------------------------------------------------------------------------------
const CondensateReader Condensate::data(const HybridTargetLevel tier) const {
  return CondensateReader(mode, basis, system_count, unit_cell, atom_starts.data(tier),
                          atom_counts.data(tier), x_coordinates_sp.data(tier),
                          y_coordinates_sp.data(tier), z_coordinates_sp.data(tier),
                          x_coordinates.data(tier), y_coordinates.data(tier),
                          z_coordinates.data(tier), box_transforms.data(tier),
                          inv_transforms.data(tier), box_dimensions.data(tier));
}

//-------------------------------------------------------------------------------------------------
CondensateWriter Condensate::data(const HybridTargetLevel tier) {
  return CondensateWriter(mode, basis, system_count, unit_cell, atom_starts.data(tier),
                          atom_counts.data(tier), x_coordinates_sp.data(tier),
                          y_coordinates_sp.data(tier), z_coordinates_sp.data(tier),
                          x_coordinates.data(tier), y_coordinates.data(tier),
                          z_coordinates.data(tier), box_transforms.data(tier),
                          inv_transforms.data(tier), box_dimensions.data(tier));
}

//-------------------------------------------------------------------------------------------------
const CondensateBorders Condensate::borders(const HybridTargetLevel tier) const {
  return CondensateBorders(system_count, atom_starts.data(tier), atom_counts.data(tier));
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
const CondensateReader Condensate::deviceViewToHostData() const {
  const size_t* astarts = atom_starts.getDeviceValidHostPointer();
  const int* acounts = atom_counts.getDeviceValidHostPointer();
  const float* xcrd_sp = x_coordinates_sp.getDeviceValidHostPointer();
  const float* ycrd_sp = y_coordinates_sp.getDeviceValidHostPointer();
  const float* zcrd_sp = z_coordinates_sp.getDeviceValidHostPointer();
  const double* xcrd = x_coordinates.getDeviceValidHostPointer();
  const double* ycrd = y_coordinates.getDeviceValidHostPointer();
  const double* zcrd = z_coordinates.getDeviceValidHostPointer();
  const double* umat_ptr = box_transforms.getDeviceValidHostPointer();
  const double* invu_ptr = inv_transforms.getDeviceValidHostPointer();
  const double* bdim_ptr = box_dimensions.getDeviceValidHostPointer();
  return CondensateReader(mode, basis, system_count, unit_cell, astarts, acounts, xcrd_sp, ycrd_sp,
                          zcrd_sp, xcrd, ycrd, zcrd, umat_ptr, invu_ptr, bdim_ptr);
}

//-------------------------------------------------------------------------------------------------
CondensateWriter Condensate::deviceViewToHostData() {
  size_t* astarts = atom_starts.getDeviceValidHostPointer();
  int* acounts = atom_counts.getDeviceValidHostPointer();
  float* xcrd_sp = x_coordinates_sp.getDeviceValidHostPointer();
  float* ycrd_sp = y_coordinates_sp.getDeviceValidHostPointer();
  float* zcrd_sp = z_coordinates_sp.getDeviceValidHostPointer();
  double* xcrd = x_coordinates.getDeviceValidHostPointer();
  double* ycrd = y_coordinates.getDeviceValidHostPointer();
  double* zcrd = z_coordinates.getDeviceValidHostPointer();
  double* umat_ptr = box_transforms.getDeviceValidHostPointer();
  double* invu_ptr = inv_transforms.getDeviceValidHostPointer();
  double* bdim_ptr = box_dimensions.getDeviceValidHostPointer();
  return CondensateWriter(mode, basis, system_count, unit_cell, astarts, acounts, xcrd_sp, ycrd_sp,
                          zcrd_sp, xcrd, ycrd, zcrd, umat_ptr, invu_ptr, bdim_ptr);
}

//-------------------------------------------------------------------------------------------------
void Condensate::upload() {
  atom_starts.upload();
  atom_counts.upload();
  float_data.upload();
  double_data.upload();
}

//-------------------------------------------------------------------------------------------------
void Condensate::download() {
  atom_starts.download();
  atom_counts.download();
  float_data.download();
  double_data.download();
}
#endif

//-------------------------------------------------------------------------------------------------
void Condensate::rebuild(const PhaseSpaceSynthesis *poly_ps_in, const PrecisionModel mode_in,
                         const GpuDetails &gpu) {

  // The basis must be set to a synthesis, even if that synthesis is the nullptr.
  basis = StructureSource::SYNTHESIS;
  
  // Exit if the synthesis is the null pointer.  Only PhaseSpaceSynthesis objects will come in as
  // nullptr
  if (poly_ps_in == nullptr) {
    system_count = 0;
    unit_cell = UnitCellType::NONE;
    atom_starts.resize(0);
    atom_counts.resize(0);
    float_data.resize(0);
    float_data.shrinkToFit();
    double_data.resize(0);
    double_data.shrinkToFit();
    x_coordinates.setPointer(&double_data, 0, 0);
    y_coordinates.setPointer(&double_data, 0, 0);
    z_coordinates.setPointer(&double_data, 0, 0);
    x_coordinates_sp.setPointer(&float_data, 0, 0);
    y_coordinates_sp.setPointer(&float_data, 0, 0);
    z_coordinates_sp.setPointer(&float_data, 0, 0);
    box_transforms.setPointer(&double_data, 0, 0);
    inv_transforms.setPointer(&double_data, 0, 0);
    box_dimensions.setPointer(&double_data, 0, 0);
    return;
  }

  // Build based on the synthesis
  pps_ptr = const_cast<PhaseSpaceSynthesis*>(poly_ps_in);
  cs_ptr = nullptr;
  const PsSynthesisReader poly_psr = pps_ptr->data();
  system_count = poly_psr.system_count;
  unit_cell = poly_psr.unit_cell;
  atom_starts.resize(system_count);
  atom_counts.resize(system_count);
  size_t* atom_starts_ptr = atom_starts.data();
  int*    atom_counts_ptr = atom_counts.data();
  for (int i = 0; i < system_count; i++) {
    atom_starts_ptr[i] = poly_psr.atom_starts[i];
    atom_counts_ptr[i] = poly_psr.atom_counts[i];
  }
  const int last_sys = poly_psr.system_count - 1;
  const size_t padded_atoms = poly_psr.atom_starts[last_sys] +
                              roundUp(poly_psr.atom_counts[last_sys], warp_size_int);
  const size_t xfrm_spacing = system_count * roundUp<size_t>(9, warp_size_zu);
  mode = mode_in;
  holds_own_data = true;
  switch (mode) {
  case PrecisionModel::DOUBLE:
    float_data.resize(0);
    float_data.shrinkToFit();
    x_coordinates_sp.setPointer(&float_data, 0, 0);
    y_coordinates_sp.setPointer(&float_data, 0, 0);
    z_coordinates_sp.setPointer(&float_data, 0, 0);
    double_data.resize((3LLU * padded_atoms) + (3 * xfrm_spacing));
    double_data.shrinkToFit();
    x_coordinates.setPointer(&double_data,                   0, padded_atoms);
    y_coordinates.setPointer(&double_data,        padded_atoms, padded_atoms);
    z_coordinates.setPointer(&double_data, 2LLU * padded_atoms, padded_atoms);
    box_transforms.setPointer(&double_data,  3LLU * padded_atoms                , xfrm_spacing);
    inv_transforms.setPointer(&double_data, (3LLU * padded_atoms) + xfrm_spacing, xfrm_spacing);
    box_dimensions.setPointer(&double_data, (3LLU * padded_atoms) + (2LLU *xfrm_spacing),
                              xfrm_spacing);
    break;
  case PrecisionModel::SINGLE:
    float_data.resize(3LLU * padded_atoms);
    float_data.shrinkToFit();
    x_coordinates_sp.setPointer(&float_data,                   0, padded_atoms);
    y_coordinates_sp.setPointer(&float_data,        padded_atoms, padded_atoms);
    z_coordinates_sp.setPointer(&float_data, 2LLU * padded_atoms, padded_atoms);
    double_data.resize(3LLU * xfrm_spacing);
    double_data.shrinkToFit();
    x_coordinates.setPointer(&double_data, 0, 0);
    y_coordinates.setPointer(&double_data, 0, 0);
    z_coordinates.setPointer(&double_data, 0, 0);
    box_transforms.setPointer(&double_data, 0, xfrm_spacing);
    inv_transforms.setPointer(&double_data,         xfrm_spacing,  xfrm_spacing);
    box_dimensions.setPointer(&double_data, (2LLU * xfrm_spacing), xfrm_spacing);
    break;
  }

  // Map the work units.  Condensates based on CoordinateSeries objects have a similar process,
  // but with many topologies to consider the PhaseSpaceSynthesis is different enough that the
  // builder has its own code.
  update();
}

//-------------------------------------------------------------------------------------------------
void Condensate::rebuild(const PhaseSpaceSynthesis &poly_ps_in, const PrecisionModel mode_in,
                         const GpuDetails &gpu) {
  rebuild(poly_ps_in.getSelfPointer(), mode_in, gpu);
}

//-------------------------------------------------------------------------------------------------
void Condensate::update(const HybridTargetLevel tier, const GpuDetails &gpu) {

  // Updating is only relevant if the object holds its own, separate copy of the coordinates.
  if (holds_own_data) {
    if (pps_ptr != nullptr) {
      switch (tier) {
      case HybridTargetLevel::HOST:
        {
          const PsSynthesisReader poly_psr = pps_ptr->data(tier);
          CondensateWriter cdw = this->data(tier);
          for (int i = 0; i < poly_psr.system_count; i++) {
            const size_t llim = poly_psr.atom_starts[i];
            const size_t hlim = llim + poly_psr.atom_counts[i];
            for (size_t j = llim; j < hlim; j++) {
              double xij = poly_psr.xcrd[j];
              double yij = poly_psr.ycrd[j];
              double zij = poly_psr.zcrd[j];
              if (poly_psr.gpos_bits >= globalpos_scale_nonoverflow_bits) {
                xij += static_cast<double>(poly_psr.xcrd_ovrf[j]) * max_llint_accumulation;
                yij += static_cast<double>(poly_psr.ycrd_ovrf[j]) * max_llint_accumulation;
                zij += static_cast<double>(poly_psr.zcrd_ovrf[j]) * max_llint_accumulation;
              }
              xij *= poly_psr.inv_gpos_scale;
              yij *= poly_psr.inv_gpos_scale;
              zij *= poly_psr.inv_gpos_scale;
              switch (mode) {
              case PrecisionModel::DOUBLE:
                cdw.xcrd[j] = xij;
                cdw.ycrd[j] = yij;
                cdw.zcrd[j] = zij;
                break;
              case PrecisionModel::SINGLE:
                cdw.xcrd_sp[j] = xij;
                cdw.ycrd_sp[j] = yij;
                cdw.zcrd_sp[j] = zij;
                break;
              }
            }
            const size_t xllim = static_cast<size_t>(i) * roundUp<size_t>(9, warp_size_zu);
            const size_t xhlim = xllim + 9LLU;
            for (size_t j = xllim; j < xhlim; j++) {
              cdw.umat[j] = poly_psr.umat[j];
              cdw.invu[j] = poly_psr.invu[j];
            }
            const size_t bllim = static_cast<size_t>(i) * roundUp<size_t>(6, warp_size_zu);
            const size_t bhlim = xllim + 6LLU;
            for (size_t j = xllim; j < xhlim; j++) {
              cdw.boxdims[j] = poly_psr.boxdims[i];
            }
          }
        }
        break;
#ifdef STORMM_USE_HPC
      case HybridTargetLevel::DEVICE:
        launchCondensateUpdate(gpu);
        break;
#endif
      }
    }
    else if (cs_ptr != nullptr) {

      // If not based on a PhaseSpaceSynthesis, the Condensate must be based on a CoordinateSeries.
      // Assess the type of the CoordinateSeries, making a switch-like apparatus to fire off the
      // templated overload of this function with the correct CoordinateSeries<T> pointer.
      if (csptr_data_type == double_type_index) {
        update(reinterpret_cast<CoordinateSeries<double>*>(cs_ptr), tier, gpu);
      }
      else if (csptr_data_type == float_type_index) {
        update(reinterpret_cast<CoordinateSeries<float>*>(cs_ptr), tier, gpu);
      }
      else if (csptr_data_type == llint_type_index) {
        update(reinterpret_cast<CoordinateSeries<llint>*>(cs_ptr), tier, gpu);
      }
      else if (csptr_data_type == int_type_index) {
        update(reinterpret_cast<CoordinateSeries<int>*>(cs_ptr), tier, gpu);
      }
      else if (csptr_data_type == short_type_index) {
        update(reinterpret_cast<CoordinateSeries<int>*>(cs_ptr), tier, gpu);
      }
      else if (csptr_data_type == char_type_index) {
        update(reinterpret_cast<CoordinateSeries<char>*>(cs_ptr), tier, gpu);
      }
      else {
        rtErr("A CoordinateSeries must be typed as double, float, or some signed integer in "
              "order to submit for analysis.  Check the data type, or call the update() function "
              "by directly supplying a pointer to the original CoordinateSeries.", "Condensate",
              "update");
      }
    }
    else {
      rtErr("There is no current coordinate synthesis or series to base the object upon.",
            "Condensate", "update");
    }
  }
}

//-------------------------------------------------------------------------------------------------
void Condensate::repairPointers() {

  // Repairs only occur if the object holds its own data.  Otherwise, a Condensate that points at
  // the contents of a CoordinateSeries can be copied with no need for pointer manipulation.
  if (holds_own_data) {
    x_coordinates_sp.swapTarget(&float_data);
    y_coordinates_sp.swapTarget(&float_data);
    z_coordinates_sp.swapTarget(&float_data);
    x_coordinates.swapTarget(&double_data);
    y_coordinates.swapTarget(&double_data);
    z_coordinates.swapTarget(&double_data);
    box_transforms.swapTarget(&double_data);
    inv_transforms.swapTarget(&double_data);
    box_dimensions.swapTarget(&double_data);
  }
}

//-------------------------------------------------------------------------------------------------
void Condensate::validateSystemIndex(const int index) const {
  if (index < 0 || index >= system_count) {
    rtErr("System index " + std::to_string(index) + " is invalid for a synthesis of " +
          std::to_string(system_count) + " systems.", "Condensate", "validateSystemIndex");
  }
}

} // namespace synthesis
} // namespace stormm
