#include "copyright.h"
#include "Constants/fixed_precision.h"
#include "structure_ops.h"
#include "radius_gyration.h"

namespace stormm {
namespace structure {

using numerics::globalpos_scale_nonoverflow_bits;

//-------------------------------------------------------------------------------------------------
double radiusOfGyration(const CoordinateFrameReader &cfr, const ChemicalDetailsKit &cdk,
                        const HybridTargetLevel tier) {
  double result = 0.0;
  switch (tier) {
  case HybridTargetLevel::HOST:
    {
      const double3 c_of_m = centerOfMass(cfr.xcrd, cfr.ycrd, cfr.zcrd, cdk.masses,
                                          0, cdk.natom);
      double total_mass = 0.0;
      for (int i = 0; i < cdk.natom; i++) {
        const double dx = cfr.xcrd[i] - c_of_m.x;
        const double dy = cfr.ycrd[i] - c_of_m.y;
        const double dz = cfr.zcrd[i] - c_of_m.z;
        const double tm = cdk.masses[i];
        result += ((dx * dx) + (dy * dy) + (dz * dz)) * tm;
        total_mass += tm;        
      }
      result = sqrt(result / total_mass);
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double radiusOfGyration(const CoordinateFrame *cf, const AtomGraph &ag,
                        const HybridTargetLevel tier) {
  return radiusOfGyration(cf->data(tier), ag.getChemicalDetailsKit(tier), tier);
}

//-------------------------------------------------------------------------------------------------
double radiusOfGyration(const CoordinateFrame &cf, const AtomGraph &ag,
                        const HybridTargetLevel tier) {
  return radiusOfGyration(cf.data(tier), ag.getChemicalDetailsKit(tier), tier);
}

//-------------------------------------------------------------------------------------------------
double radiusOfGyration(const PhaseSpace *ps, const AtomGraph &ag, const HybridTargetLevel tier) {
  return radiusOfGyration(CoordinateFrameReader(ps), ag.getChemicalDetailsKit(tier), tier);
}

//-------------------------------------------------------------------------------------------------
double radiusOfGyration(const PhaseSpace &ps, const AtomGraph &ag, const HybridTargetLevel tier) {
  return radiusOfGyration(CoordinateFrameReader(ps), ag.getChemicalDetailsKit(tier), tier);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> radiusOfGyration(const PsSynthesisReader &poly_psr,
                                     const SyAtomUpdateKit<double, double2, double4> &poly_auk,
                                     const HybridTargetLevel tier) {
  std::vector<double> result(poly_psr.system_count);
  switch (tier) {
  case HybridTargetLevel::HOST:
    for (int i = 0; i < poly_psr.system_count; i++) {
      double3 c_of_m = { 0.0, 0.0, 0.0 };
      double total_mass = 0.0;
      const int jlim = poly_psr.atom_counts[i] + poly_psr.atom_starts[i];
      if (poly_psr.gpos_bits <= globalpos_scale_nonoverflow_bits) {
        for (int j = poly_psr.atom_starts[i]; j < jlim; j++) {
          c_of_m.x += static_cast<double>(poly_psr.xcrd[j]);
          c_of_m.y += static_cast<double>(poly_psr.ycrd[j]);
          c_of_m.z += static_cast<double>(poly_psr.zcrd[j]);
          total_mass += poly_auk.masses[j];
        }
      }
      else {
        for (int j = poly_psr.atom_starts[i]; j < jlim; j++) {
          c_of_m.x += hostInt95ToDouble(poly_psr.xcrd[j], poly_psr.xcrd_ovrf[j]);
          c_of_m.y += hostInt95ToDouble(poly_psr.ycrd[j], poly_psr.ycrd_ovrf[j]);
          c_of_m.z += hostInt95ToDouble(poly_psr.zcrd[j], poly_psr.zcrd_ovrf[j]);
          total_mass += poly_auk.masses[j];
        }
      }
      c_of_m.x /= total_mass;
      c_of_m.y /= total_mass;
      c_of_m.z /= total_mass;
      double rslti = 0.0;
      if (poly_psr.gpos_bits <= globalpos_scale_nonoverflow_bits) {
        for (int j = poly_psr.atom_starts[i]; j < jlim; j++) {
          const double dx = static_cast<double>(poly_psr.xcrd[i]) - c_of_m.x;
          const double dy = static_cast<double>(poly_psr.ycrd[i]) - c_of_m.y;
          const double dz = static_cast<double>(poly_psr.zcrd[i]) - c_of_m.z;
          rslti += ((dx * dx) + (dy * dy) + (dz * dz)) * poly_auk.masses[j];
        }
      }
      else {
        for (int j = poly_psr.atom_starts[i]; j < jlim; j++) {
          const double dx = hostInt95ToDouble(poly_psr.xcrd[i], poly_psr.xcrd_ovrf[i]) - c_of_m.x;
          const double dy = hostInt95ToDouble(poly_psr.ycrd[i], poly_psr.ycrd_ovrf[i]) - c_of_m.y;
          const double dz = hostInt95ToDouble(poly_psr.zcrd[i], poly_psr.zcrd_ovrf[i]) - c_of_m.z;
          rslti += ((dx * dx) + (dy * dy) + (dz * dz)) * poly_auk.masses[j];
        }
      }
      result[i] = sqrt(rslti / total_mass) * poly_psr.inv_gpos_scale;
    }
    break;
#ifdef STORMM_USE_HPC
  case HybridTargetLevel::DEVICE:
    break;
#endif
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> radiusOfGyration(const PhaseSpaceSynthesis *poly_ps,
                                     const AtomGraphSynthesis &poly_ag,
                                     const HybridTargetLevel tier) {
  return radiusOfGyration(poly_ps->data(tier), poly_ag.getDoublePrecisionAtomUpdateKit(), tier);
}

//-------------------------------------------------------------------------------------------------
std::vector<double> radiusOfGyration(const PhaseSpaceSynthesis &poly_ps,
                                     const AtomGraphSynthesis &poly_ag,
                                     const HybridTargetLevel tier) {
  return radiusOfGyration(poly_ps.data(tier), poly_ag.getDoublePrecisionAtomUpdateKit(), tier);
}

} // namespace structure
} // namespace stormm
