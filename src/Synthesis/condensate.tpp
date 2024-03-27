// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace synthesis {

//-------------------------------------------------------------------------------------------------
template <typename T>
Condensate::Condensate(const CoordinateSeries<T> *cs_in, const PrecisionModel mode_in,
                       const GpuDetails &gpu) :
    Condensate(nullptr, mode_in, gpu)
{
  rebuild(cs_in, mode, gpu);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
Condensate::Condensate(const CoordinateSeries<T> &cs_in, const PrecisionModel mode_in,
                       const GpuDetails &gpu) :
    Condensate(cs_in.getSelfPointer(), mode_in, gpu)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
const CoordinateSeries<T>* Condensate::getSeriesPointer() const {
  if (cs_ptr == nullptr) {
    rtErr("No CoordinateSeries is referenced.", "Condensate", "getSeriesPointer");
  }
  return reinterpret_cast<CoordinateSeries<T>*>(cs_ptr);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void Condensate::rebuild(const CoordinateSeries<T> *cs_in, const PrecisionModel mode_in,
                         const GpuDetails &gpu) {
  basis = StructureSource::SERIES;

  // Reinterpret the templated CoordinatesSeries<T> pointer as a CoordinateSeries of an arbitrary
  // type.  This prevents the Condensate class as a whole from taking on a template requirement.
  // Also correct the basis to indicate the source of the data.
  cs_ptr = reinterpret_cast<CoordinateSeries<int>*>(const_cast<CoordinateSeries<T>*>(cs_in));
  pps_ptr = nullptr;
  csptr_data_type = std::type_index(typeid(T)).hash_code();
  basis = StructureSource::SERIES;
  const CoordinateSeriesReader<T> csr = cs_in->data();
  system_count = csr.nframe;
  unit_cell = csr.unit_cell;
  atom_starts.resize(system_count);
  atom_counts.resize(system_count);
  size_t* atom_starts_ptr = atom_starts.data();
  int*    atom_counts_ptr = atom_counts.data();
  size_t curr_start = 0;
  const size_t padded_atom_count = roundUp(csr.natom, warp_size_int);
  for (int i = 0; i < system_count; i++) {
    atom_starts_ptr[i] = curr_start;
    atom_counts_ptr[i] = csr.natom;
    curr_start += padded_atom_count;
  }
  const size_t padded_atoms = static_cast<size_t>(csr.nframe) *
                              roundUp(static_cast<size_t>(csr.natom), warp_size_zu);
  const size_t xfrm_spacing = system_count * roundUp<size_t>(9, warp_size_zu);
  const CartesianDimension xdim = CartesianDimension::X;
  const CartesianDimension ydim = CartesianDimension::Y;
  const CartesianDimension zdim = CartesianDimension::Z;
  mode = mode_in;
  const size_t ct = std::type_index(typeid(T)).hash_code();

  // If the coordinate series is of a suitable type, set the pointers of this object to its own
  // data arrays.
  CoordinateSeries<float> *cs_in_f = (ct == float_type_index) ?
    reinterpret_cast<CoordinateSeries<float>*>(const_cast<CoordinateSeries<T>*>(cs_in)) : nullptr;
  CoordinateSeries<double> *cs_in_d = (ct == double_type_index) ?
    reinterpret_cast<CoordinateSeries<double>*>(const_cast<CoordinateSeries<T>*>(cs_in)) : nullptr;
  holds_own_data = (cs_in_f == nullptr && cs_in_d == nullptr);
  switch (mode) {
  case PrecisionModel::DOUBLE:
    float_data.resize(0);
    float_data.shrinkToFit();
    x_coordinates_sp.setPointer(&float_data, 0, 0);
    y_coordinates_sp.setPointer(&float_data, 0, 0);
    z_coordinates_sp.setPointer(&float_data, 0, 0);
    if (cs_in_d == nullptr) {
      double_data.resize((3LLU * padded_atoms) + (2LLU * xfrm_spacing));
      double_data.shrinkToFit();
      x_coordinates.setPointer(&double_data,                   0, padded_atoms);
      y_coordinates.setPointer(&double_data,        padded_atoms, padded_atoms);
      z_coordinates.setPointer(&double_data, 2LLU * padded_atoms, padded_atoms);
      box_transforms.setPointer(&double_data,  3LLU * padded_atoms                , xfrm_spacing);
      inv_transforms.setPointer(&double_data, (3LLU * padded_atoms) + xfrm_spacing, xfrm_spacing);
    }
    else {
      double_data.resize(0);
      double_data.shrinkToFit();
      x_coordinates.setPointer(const_cast<Hybrid<double>*>(cs_in_d->getCoordinatePointer(xdim)),
                               0, padded_atoms);
      y_coordinates.setPointer(const_cast<Hybrid<double>*>(cs_in_d->getCoordinatePointer(ydim)),
                               0, padded_atoms);
      z_coordinates.setPointer(const_cast<Hybrid<double>*>(cs_in_d->getCoordinatePointer(zdim)),
                               0, padded_atoms);
      box_transforms.setPointer(const_cast<Hybrid<double>*>(cs_in_d->getBoxTransformPointer()));
      const Hybrid<double>* cs_in_dptr = cs_in_d->getInverseTransformPointer();
      inv_transforms.setPointer(const_cast<Hybrid<double>*>(cs_in_dptr));
      box_dimensions.setPointer(const_cast<Hybrid<double>*>(cs_in_d->getBoxDimensionPointer()));
    }
    break;
  case PrecisionModel::SINGLE:
    double_data.resize(2LLU * xfrm_spacing);
    double_data.shrinkToFit();
    x_coordinates.setPointer(&double_data, 0, 0);
    y_coordinates.setPointer(&double_data, 0, 0);
    z_coordinates.setPointer(&double_data, 0, 0);
    if (cs_in_f == nullptr) {
      float_data.resize(3LLU * padded_atoms);
      float_data.shrinkToFit();
      x_coordinates_sp.setPointer(&float_data,                   0, padded_atoms);
      y_coordinates_sp.setPointer(&float_data,        padded_atoms, padded_atoms);
      z_coordinates_sp.setPointer(&float_data, 2LLU * padded_atoms, padded_atoms);
      box_transforms.setPointer(&double_data,            0, xfrm_spacing);
      inv_transforms.setPointer(&double_data, xfrm_spacing, xfrm_spacing);
    }
    else {
      float_data.resize(0);
      float_data.shrinkToFit();
      x_coordinates_sp.setPointer(const_cast<Hybrid<float>*>(cs_in_f->getCoordinatePointer(xdim)),
                                  0, padded_atoms);
      y_coordinates_sp.setPointer(const_cast<Hybrid<float>*>(cs_in_f->getCoordinatePointer(ydim)),
                                  0, padded_atoms);
      z_coordinates_sp.setPointer(const_cast<Hybrid<float>*>(cs_in_f->getCoordinatePointer(zdim)),
                                  0, padded_atoms);

      // Even the float-type CoordinateSeries will emit double-precision transform pointers.
      box_transforms.setPointer(const_cast<Hybrid<double>*>(cs_in_f->getBoxTransformPointer()));
      const Hybrid<double>* cs_in_fptr = cs_in_f->getInverseTransformPointer();
      inv_transforms.setPointer(const_cast<Hybrid<double>*>(cs_in_fptr));
      box_dimensions.setPointer(const_cast<Hybrid<double>*>(cs_in_f->getBoxDimensionPointer()));
    }
    break;
  }
  update(cs_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void Condensate::update(const CoordinateSeries<T> *cs_basis, const HybridTargetLevel tier,
                        const GpuDetails &gpu) {

  // Do not perform any updates if the Condensate already points to the CoordinateSeries
  if (holds_own_data == false) {
    return;
  }

  // Produce an error if the Condensate is not ready to accept the CoordinateSeries in question.
  if (cs_basis != reinterpret_cast<CoordinateSeries<T>*>(cs_ptr)) {
    rtErr("The pointer to the CoordinateSeries supplied does not match the pointer stored "
          "internally.", "Condensate", "Update");
  }

  // Copy the data
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const CoordinateSeriesReader<T> csr = cs_basis->data(tier);
      CondensateWriter cdw = this->data(tier);
      const size_t natom_zu = csr.natom;
      const size_t padded_natom = roundUp(natom_zu, warp_size_zu);
      const size_t nframe_zu = csr.nframe;
      for (size_t i = 0; i < nframe_zu; i++) {
        const size_t llim = i * padded_natom;
        const size_t hlim = llim + natom_zu;
        for (size_t j = llim; j < hlim; j++) {
          double xij = csr.xcrd[j] * csr.inv_gpos_scale;
          double yij = csr.ycrd[j] * csr.inv_gpos_scale;
          double zij = csr.zcrd[j] * csr.inv_gpos_scale;
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

//-------------------------------------------------------------------------------------------------
template <typename T>
void Condensate::update(const CoordinateSeries<T> &cs_basis, const HybridTargetLevel tier,
                        const GpuDetails &gpu) {
  update(cs_basis.getSelfPointer(), tier, gpu);
}

} // namespace synthesis
} // namespace stormm
