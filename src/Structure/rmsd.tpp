// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename T>
void checkCompatibility(const ComparisonGuide &cg, const RMSDPlan &rplan,
                        const PhaseSpaceSynthesis &poly_ps, const Hybrid<T> *result,
                        const RMSDTask process, const SystemGrouping organization) {
  const PsSynthesisReader poly_psr = poly_ps.data();
  const int ntop = poly_psr.unique_topology_count;
  if (ntop != rplan.getPlanCount()) {
    rtErr("The number of unique topologies in the synthesis (" + std::to_string(ntop) + ") must "
          "match the number of individual RMSD calculation plans (" +
          std::to_string(rplan.getPlanCount()) + ").", "checkCompatibility");
  }
  switch (process) {
  case RMSDTask::REFERENCE:
    if (result->size() != cg.getAllToReferenceOutputSize()) {
      rtErr("The anticipated total size of RMSD-to-reference calculations (" +
            std::to_string(cg.getAllToReferenceOutputSize()) + ") is not "
            "matched by the size of the results array (" + std::to_string(result->size()) + ").",
            "checkCompatibility");
    }
    break;
  case RMSDTask::MATRIX:
    if (result->size() != cg.getSymmetryEquivalentPairOutputSize(organization)) {
      rtErr("The anticipated total size of RMSD-to-reference calculations (" +
            std::to_string(cg.getSymmetryEquivalentPairOutputSize(organization)) + ") is not "
            "matched by the size of the results array (" + std::to_string(result->size()) + ").",
            "checkCompatibility");
    }
    break;
  }
  if (rplan.getSynthesisMapPointer() != nullptr) {
    const SynthesisCacheMap *scmap = rplan.getSynthesisMapPointer();
    if (scmap->getCoordinateSynthesisPointer() != poly_ps.getSelfPointer()) {
      bool problem = (scmap->getTopologySynthesisPointer() == nullptr);

      // Check the syntheses for a reasonable degree of similarity
      if (problem == false) {
        const SyNonbondedKit<double, double2> multi_nbk =
          scmap->getTopologySynthesisPointer()->getDoublePrecisionNonbondedKit();
        const PsSynthesisReader poly_psr = poly_ps.data();
        problem = (multi_nbk.nsys != poly_psr.system_count ||
                   scmap->getTopologySynthesisPointer()->getUniqueTopologyCount() !=
                   poly_ps.getUniqueTopologyCount());
        if (problem == false) {
          for (int i = 0; i < poly_psr.system_count; i++) {
            problem = (problem || poly_psr.atom_counts[i] != multi_nbk.atom_counts[i]);
          }
        }
      }
      if (problem) {
        rtErr("The condensed coordinate representation does not reference the original synthesis.",
              "checkCompatibility");
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void checkCompatibility(const ComparisonGuide &cg, const RMSDPlan &rplan,
                        const PhaseSpaceSynthesis &poly_ps, const Condensate &cdns,
                        const Hybrid<T> *result, const RMSDTask process,
                        const SystemGrouping organization) {
  checkCompatibility(cg, rplan, poly_ps, result, process, organization);
  if (cdns.getSynthesisPointer() != poly_ps.getSelfPointer()) {
    rtErr("The condensed coordinate representation does not reference the original synthesis.",
          "checkCompatibility");
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc rmsd(const Tcoord* xcrd_a, const Tcoord* ycrd_a, const Tcoord* zcrd_a, const Tcoord* xcrd_b,
           const Tcoord* ycrd_b, const Tcoord* zcrd_b, const Tcalc* masses,
           const RMSDMethod method, const int lower_limit, const int upper_limit,
           const Tcalc inv_gpos_scale_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const Tcalc value_one = 1.0;
  const Tcalc value_zero = 0.0;
  Tcalc result = value_zero;
  Tcalc inv_mass_divisor;
  switch (method) {
  case RMSDMethod::ALIGN_MASS:
  case RMSDMethod::NO_ALIGN_MASS:
    inv_mass_divisor = value_one / sum<Tcalc>(&masses[lower_limit], upper_limit - lower_limit);
    break;
  case RMSDMethod::ALIGN_GEOM:
  case RMSDMethod::NO_ALIGN_GEOM:
    inv_mass_divisor = value_one / static_cast<Tcalc>(upper_limit - lower_limit);
    break;
  }
  switch (method) {
  case RMSDMethod::ALIGN_MASS:
  case RMSDMethod::ALIGN_GEOM:
    {
      // In order to compute the aligned RMSD without moving the coordinates, the movement must
      // be computed in temporary variables only.
      Tcalc sab_xx = value_zero;
      Tcalc sab_xy = value_zero;
      Tcalc sab_xz = value_zero;
      Tcalc sab_yx = value_zero;
      Tcalc sab_yy = value_zero;
      Tcalc sab_yz = value_zero;
      Tcalc sab_zx = value_zero;
      Tcalc sab_zy = value_zero;
      Tcalc sab_zz = value_zero;
      Tcalc sa_x = value_zero;
      Tcalc sa_y = value_zero;
      Tcalc sa_z = value_zero;
      Tcalc sb_x = value_zero;
      Tcalc sb_y = value_zero;
      Tcalc sb_z = value_zero;
      const bool use_mass = (method == RMSDMethod::ALIGN_MASS);
      if (tcoord_is_sgnint) {
        for (int i = lower_limit; i < upper_limit; i++) {
          const Tcalc tmass = (use_mass) ? masses[i] : 1.0;
          const Tcalc xca = static_cast<Tcalc>(xcrd_a[i]) * inv_gpos_scale_factor;
          const Tcalc yca = static_cast<Tcalc>(ycrd_a[i]) * inv_gpos_scale_factor;
          const Tcalc zca = static_cast<Tcalc>(zcrd_a[i]) * inv_gpos_scale_factor;
          const Tcalc xcb = static_cast<Tcalc>(xcrd_b[i]) * inv_gpos_scale_factor;
          const Tcalc ycb = static_cast<Tcalc>(ycrd_b[i]) * inv_gpos_scale_factor;
          const Tcalc zcb = static_cast<Tcalc>(zcrd_b[i]) * inv_gpos_scale_factor;
          sab_xx += tmass * xca * xcb;
          sab_xy += tmass * xca * ycb;
          sab_xz += tmass * xca * zcb;
          sab_yx += tmass * yca * xcb;
          sab_yy += tmass * yca * ycb;
          sab_yz += tmass * yca * zcb;
          sab_zx += tmass * zca * xcb;
          sab_zy += tmass * zca * ycb;
          sab_zz += tmass * zca * zcb;
          sa_x += tmass * xca;
          sa_y += tmass * yca;
          sa_z += tmass * zca;
          sb_x += tmass * xcb;
          sb_y += tmass * ycb;
          sb_z += tmass * zcb;
        }
      }
      else {
        for (int i = lower_limit; i < upper_limit; i++) {
          const Tcalc tmass = (use_mass) ? masses[i] : 1.0;
          sab_xx += tmass * xcrd_a[i] * xcrd_b[i];
          sab_xy += tmass * xcrd_a[i] * ycrd_b[i];
          sab_xz += tmass * xcrd_a[i] * zcrd_b[i];
          sab_yx += tmass * ycrd_a[i] * xcrd_b[i];
          sab_yy += tmass * ycrd_a[i] * ycrd_b[i];
          sab_yz += tmass * ycrd_a[i] * zcrd_b[i];
          sab_zx += tmass * zcrd_a[i] * xcrd_b[i];
          sab_zy += tmass * zcrd_a[i] * ycrd_b[i];
          sab_zz += tmass * zcrd_a[i] * zcrd_b[i];
          sa_x += tmass * xcrd_a[i];
          sa_y += tmass * ycrd_a[i];
          sa_z += tmass * zcrd_a[i];
          sb_x += tmass * xcrd_b[i];
          sb_y += tmass * ycrd_b[i];
          sb_z += tmass * zcrd_b[i];
        }
      }
      const Tcalc coma_x = sa_x * inv_mass_divisor;
      const Tcalc coma_y = sa_y * inv_mass_divisor;
      const Tcalc coma_z = sa_z * inv_mass_divisor;
      const Tcalc comb_x = sb_x * inv_mass_divisor;
      const Tcalc comb_y = sb_y * inv_mass_divisor;
      const Tcalc comb_z = sb_z * inv_mass_divisor;

      // Assemble the Kabsch matrix and diagonalize it.  The weighted product of the centers of
      // mass is equal to the weighted sum sa_(x,y,z) or sb_(x,y,z) times the center of mass
      // comb_(x,y,z) or coma_(x,y,z), so rather than evaluate the full F.O.I.L. to complete
      // mass * ((x,y,z)crd_a - coma_(x,y,z)) * ((x,y,z)crd_b - comb_(x,y,z)), the final two
      // terms cancel.
      const Tcalc aa = sab_xx - (sb_x * coma_x);
      const Tcalc ab = sab_xy - (sb_y * coma_x);
      const Tcalc ac = sab_xz - (sb_z * coma_x);
      const Tcalc ba = sab_yx - (sb_x * coma_y);
      const Tcalc bb = sab_yy - (sb_y * coma_y);
      const Tcalc bc = sab_yz - (sb_z * coma_y);
      const Tcalc ca = sab_zx - (sb_x * coma_z);
      const Tcalc cb = sab_zy - (sb_y * coma_z);
      const Tcalc cc = sab_zz - (sb_z * coma_z);
      std::vector<Tcalc> rmat(16);
      rmat[ 0] = aa + bb + cc;
      rmat[ 1] = cb - bc;
      rmat[ 2] = ac - ca;
      rmat[ 3] = ba - ab;
      rmat[ 5] = aa - bb - cc;
      rmat[ 6] = ab + ba;
      rmat[ 7] = ca + ac;
      rmat[10] = bb - aa - cc;
      rmat[11] = bc + cb;
      rmat[15] = cc - aa - bb;
      rmat[ 4] = rmat[ 1];
      rmat[ 8] = rmat[ 2];
      rmat[12] = rmat[ 3];
      rmat[ 9] = rmat[ 6];
      rmat[13] = rmat[ 7];
      rmat[14] = rmat[11];
      for (int i = 0; i < 16; i++) {
        rmat[i] *= inv_mass_divisor;
      }
      std::vector<Tcalc> eigval(4, value_zero), sdiag(4, value_zero);
      realSymmEigensolver(&rmat, &eigval, &sdiag, EigenWork::EIGENVECTORS);
      int max_eig_loc = 0;
      for (int i = 1; i < 4; i++) {
        if (eigval[i] > eigval[max_eig_loc]) {
          max_eig_loc = i;
        }
      }
      
      // Form the rotation matrix.  This is prone to some instability when crd_a and crd_b are
      // close to identical, or indeed entirely identical.  As a result, when computed in single
      // precision, the RMSD of identical structures computed in 32 bit floating-point arithmetic
      // can be reported to be on the order of 0.001 Angstroms, an error much higher than that
      // computed for structures with perturbations that cause their true RMSD to sit at 0.2-0.5
      // Angstroms.
      const Tcalc a = rmat[(4 * max_eig_loc)    ];
      const Tcalc x = rmat[(4 * max_eig_loc) + 1];
      const Tcalc y = rmat[(4 * max_eig_loc) + 2];
      const Tcalc z = rmat[(4 * max_eig_loc) + 3];
      std::vector<Tcalc> umat(9);
      umat[0] = (a * a) + (x * x) - (y * y) - (z * z);
      umat[3] = 2.0 * ((x * y) + (a * z));
      umat[6] = 2.0 * ((x * z) - (a * y));
      umat[1] = 2.0 * ((x * y) - (a * z));
      umat[4] = (a * a) - (x * x) + (y * y) - (z * z);
      umat[7] = 2.0 * ((y * z) + (a * x));
      umat[2] = 2.0 * ((x * z) + (a * y));
      umat[5] = 2.0 * ((y * z) - (a * x));
      umat[8] = (a * a) - (x * x) - (y * y) + (z * z);      

      // Shift and rotate the coordinates of the first frame (in temporary variables only) and
      // compare them to the shifted, unrotated coordinates of the second frame.
      if (tcoord_is_sgnint) {
        for (int i = lower_limit; i < upper_limit; i++) {
          const Tcalc tmass = (use_mass) ? masses[i] : 1.0;
          const Tcalc shfta_x = (static_cast<Tcalc>(xcrd_a[i]) * inv_gpos_scale_factor) - coma_x;
          const Tcalc shfta_y = (static_cast<Tcalc>(ycrd_a[i]) * inv_gpos_scale_factor) - coma_y;
          const Tcalc shfta_z = (static_cast<Tcalc>(zcrd_a[i]) * inv_gpos_scale_factor) - coma_z;
          const Tcalc shftb_x = (static_cast<Tcalc>(xcrd_b[i]) * inv_gpos_scale_factor) - comb_x;
          const Tcalc shftb_y = (static_cast<Tcalc>(ycrd_b[i]) * inv_gpos_scale_factor) - comb_y;
          const Tcalc shftb_z = (static_cast<Tcalc>(zcrd_b[i]) * inv_gpos_scale_factor) - comb_z;
          const Tcalc rota_x = (umat[0] * shfta_x) + (umat[3] * shfta_y) + (umat[6] * shfta_z);
          const Tcalc rota_y = (umat[1] * shfta_x) + (umat[4] * shfta_y) + (umat[7] * shfta_z);
          const Tcalc rota_z = (umat[2] * shfta_x) + (umat[5] * shfta_y) + (umat[8] * shfta_z);
          const Tcalc dx = shftb_x - rota_x;
          const Tcalc dy = shftb_y - rota_y;
          const Tcalc dz = shftb_z - rota_z;
          result += tmass * ((dx * dx) + (dy * dy) + (dz * dz));
        }
      }
      else {
        for (int i = lower_limit; i < upper_limit; i++) {
          const Tcalc tmass = (use_mass) ? masses[i] : 1.0;
          const Tcalc shfta_x = xcrd_a[i] - coma_x;
          const Tcalc shfta_y = ycrd_a[i] - coma_y;
          const Tcalc shfta_z = zcrd_a[i] - coma_z;
          const Tcalc shftb_x = xcrd_b[i] - comb_x;
          const Tcalc shftb_y = ycrd_b[i] - comb_y;
          const Tcalc shftb_z = zcrd_b[i] - comb_z;
          const Tcalc rota_x = (umat[0] * shfta_x) + (umat[3] * shfta_y) + (umat[6] * shfta_z);
          const Tcalc rota_y = (umat[1] * shfta_x) + (umat[4] * shfta_y) + (umat[7] * shfta_z);
          const Tcalc rota_z = (umat[2] * shfta_x) + (umat[5] * shfta_y) + (umat[8] * shfta_z);
          const Tcalc dx = shftb_x - rota_x;
          const Tcalc dy = shftb_y - rota_y;
          const Tcalc dz = shftb_z - rota_z;
          result += tmass * ((dx * dx) + (dy * dy) + (dz * dz));
        }
      }
    }
    break;
  case RMSDMethod::NO_ALIGN_MASS:
    if (tcoord_is_sgnint) {
      for (int i = lower_limit; i < upper_limit; i++) {
        const Tcalc dx = static_cast<Tcalc>(xcrd_b[i] - xcrd_a[i]) * inv_gpos_scale_factor;
        const Tcalc dy = static_cast<Tcalc>(ycrd_b[i] - ycrd_a[i]) * inv_gpos_scale_factor;
        const Tcalc dz = static_cast<Tcalc>(zcrd_b[i] - zcrd_a[i]) * inv_gpos_scale_factor;
        result += masses[i] * ((dx * dx) + (dy * dy) + (dz * dz));
      }      
    }
    else {
      for (int i = lower_limit; i < upper_limit; i++) {
        const Tcalc dx = xcrd_b[i] - xcrd_a[i];
        const Tcalc dy = ycrd_b[i] - ycrd_a[i];
        const Tcalc dz = zcrd_b[i] - zcrd_a[i];
        result += masses[i] * ((dx * dx) + (dy * dy) + (dz * dz));
      }
    }
    break;
  case RMSDMethod::NO_ALIGN_GEOM:
    if (tcoord_is_sgnint) {
      for (int i = lower_limit; i < upper_limit; i++) {
        const Tcalc dx = static_cast<Tcalc>(xcrd_b[i] - xcrd_a[i]) * inv_gpos_scale_factor;
        const Tcalc dy = static_cast<Tcalc>(ycrd_b[i] - ycrd_a[i]) * inv_gpos_scale_factor;
        const Tcalc dz = static_cast<Tcalc>(zcrd_b[i] - zcrd_a[i]) * inv_gpos_scale_factor;
        result += (dx * dx) + (dy * dy) + (dz * dz);
      }      
    }
    else {
      for (int i = lower_limit; i < upper_limit; i++) {
        const Tcalc dx = xcrd_b[i] - xcrd_a[i];
        const Tcalc dy = ycrd_b[i] - ycrd_a[i];
        const Tcalc dz = zcrd_b[i] - zcrd_a[i];
        result += (dx * dx) + (dy * dy) + (dz * dz);
      }
    }
    break;
  }
  return (tcalc_is_double) ? sqrt(result * inv_mass_divisor) : sqrtf(result * inv_mass_divisor);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc rmsd(const CoordinateSeriesReader<Tcoord> &csr, const size_t frame_a, const size_t frame_b,
           const ChemicalDetailsKit &cdk, const RMSDMethod method, const int lower_limit,
           const int upper_limit) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const size_t padded_natom_zu = roundUp<size_t>(csr.natom, warp_size_zu);
  const size_t fa_offset = frame_a * padded_natom_zu;
  const size_t fb_offset = frame_b * padded_natom_zu;
  const Tcalc* mass_ptr = (tcalc_is_double) ? (Tcalc*)(cdk.masses) : (Tcalc*)(cdk.sp_masses);
  return rmsd<Tcoord, Tcalc>(&csr.xcrd[fa_offset], &csr.ycrd[fa_offset], &csr.zcrd[fa_offset],
                             &csr.xcrd[fb_offset], &csr.ycrd[fb_offset], &csr.zcrd[fb_offset],
                             mass_ptr, method, lower_limit, upper_limit, csr.inv_gpos_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc rmsd(const CoordinateSeriesWriter<Tcoord> &csw, const size_t frame_a, const size_t frame_b,
           const ChemicalDetailsKit &cdk, const RMSDMethod method, const int lower_limit,
           const int upper_limit) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const size_t padded_natom_zu = roundUp<size_t>(csw.natom, warp_size_zu);
  const size_t fa_offset = frame_a * padded_natom_zu;
  const size_t fb_offset = frame_b * padded_natom_zu;
  const Tcalc* mass_ptr = (tcalc_is_double) ? (Tcalc*)(cdk.masses) : (Tcalc*)(cdk.sp_masses);
  return rmsd<Tcoord, Tcalc>(&csw.xcrd[fa_offset], &csw.ycrd[fa_offset], &csw.zcrd[fa_offset],
                             &csw.xcrd[fb_offset], &csw.ycrd[fb_offset], &csw.zcrd[fb_offset],
                             mass_ptr, method, lower_limit, upper_limit, csw.inv_gpos_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
Tcalc rmsd(const CoordinateSeries<Tcoord> &cs, const int frame_a, const int frame_b,
           const AtomGraph &ag, const RMSDMethod method, const int lower_limit,
           const int upper_limit) {
  const int nframes = cs.getFrameCount();
  if (frame_a < 0 || frame_a >= nframes || frame_b < 0 || frame_b >= nframes) {
    rtErr("Frames " + std::to_string(frame_a) + " and " + std::to_string(frame_b) + " are not "
          "accessible in a series of " + std::to_string(nframes) + " frames.", "rmsd");
  }
  return rmsd(cs.data(), frame_a, frame_b, ag.getChemicalDetailsKit(), method, lower_limit,
              upper_limit);
}

//-------------------------------------------------------------------------------------------------
double rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan,
            const CoordinateFrameReader &reference, const CoordinateFrameReader &snapshot) {
  const ChemicalDetailsKit cdk = rplan.getTopologyPointer(0)->getChemicalDetailsKit();
  const int natom = cdk.natom;
  if (rplan.getPlanCount() == 1) {
    if (natom != reference.natom && natom != snapshot.natom) {
      rtErr("The system and topology do not have matching atom counts (" + std::to_string(natom) +
            " vs. " + std::to_string(reference.natom) + ").", "rmsd");
    }
  }
  else {
    rtErr("A plan must have only one topology in order to serve a specific coordinate frame.",
          "RMSDPlan", "rmsd");
  }
  switch (rplan.getAlignmentProtocol(0)) {
  case RMSDAlignmentProtocol::BUILD_CORE:
  case RMSDAlignmentProtocol::ALIGN_CORE:
  case RMSDAlignmentProtocol::ALIGN_ALL:
    return rmsd(reference, snapshot, cdk, rplan.getGeneralStrategy());
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
double rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const CoordinateFrame &reference,
            const CoordinateFrame &snapshot) {
  return rmsd(cg, rplan, reference.data(), snapshot.data());
}

//-------------------------------------------------------------------------------------------------
double rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const PhaseSpace &reference,
            const PhaseSpace &snapshot) {
  return rmsd(cg, rplan, CoordinateFrameReader(reference), CoordinateFrameReader(snapshot));
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<double> rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan,
                         const CoordinateFrame &reference, const CoordinateSeries<T> &snapshots) {
  const int frame_count = snapshots.getFrameCount();
  std::vector<double> result(frame_count);
  const ChemicalDetailsKit cdk = rplan.getTopologyPointer(0)->getChemicalDetailsKit();
  const CoordinateFrameReader refr = reference.data();
  for (int i = 0; i < frame_count; i++) {
    const CoordinateFrame cfi = snapshots.exportFrame(i);
    if (i == 0) {
      
      // Invoke the primary overload to engage its error checking
      result[i] = rmsd(cg, rplan, refr, cfi.data());
    }
    else {
      switch (rplan.getAlignmentProtocol(0)) {
      case RMSDAlignmentProtocol::BUILD_CORE:
      case RMSDAlignmentProtocol::ALIGN_CORE:
      case RMSDAlignmentProtocol::ALIGN_ALL:
        result[i] = rmsd(refr, cfi.data(), cdk, rplan.getGeneralStrategy());
        break;
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<double> rmsd(const RMSDPlan &rplan, const CoordinateSeries<T> &snapshots,
                         const int reference_frame) {
  const CoordinateFrame reference = snapshots.exportFrame(reference_frame);
  return rmsd(rplan, reference, snapshots);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const CoordinateFrame &reference,
          const CoordinateSeries<T> &snapshots, Hybrid<T> *result) {
  T* res_ptr = result->data();
  const int nval = result->size();
  if (nval != snapshots.getFrameCount()) {
    rtErr("A result vector of size " + std::to_string(nval) + " is not set up to accept "
          "results from a comparison of " + std::to_string(snapshots.getFrameCount()) +
          " snapshots to the reference structure.");
  }
  std::vector<double> output = rmsd(cg, rplan, reference, snapshots);
  for (int i = 0; i < nval; i++) {
    res_ptr[i] = output[i];
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const CoordinateSeries<T> &snapshots,
          const int reference_frame, Hybrid<T> *result) {
  const CoordinateFrame reference = snapshots.exportFrame(reference_frame);
  rmsd(cg, rplan, reference, snapshots, result);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const CoordinateSeries<T> &snapshots,
          Hybrid<T> *result) {
  T* res_ptr = result->data();
  const size_t nval = result->size() ;
  const size_t nfrm = snapshots.getFrameCount();
  if (nval != roundUp<size_t>((nfrm - 1LLU) * nfrm / 2LLU, warp_size_zu)) {
    rtErr("A result vector of size " + std::to_string(nval) + " is not set up to accept "
          "results from a comparison of " + std::to_string(snapshots.getFrameCount()) +
          " snapshots to the reference structure.");        
  }
  CoordinateSeriesReader<T> snapr = snapshots.data();
  const ChemicalDetailsKit cdk = rplan.getTopologyPointer(0)->getChemicalDetailsKit();
  const double* mass_ptr = cdk.masses;
  const size_t natom = snapr.natom;
  const size_t padded_natom = roundUp(snapr.natom, warp_size_int);
  switch (rplan.getAlignmentProtocol(0)) {
  case RMSDAlignmentProtocol::BUILD_CORE:
  case RMSDAlignmentProtocol::ALIGN_CORE:
  case RMSDAlignmentProtocol::ALIGN_ALL:
    for (size_t i = 1; i < nfrm; i++) {
      const T* ixcrd = &snapr.xcrd[i * padded_natom];
      const T* iycrd = &snapr.ycrd[i * padded_natom];
      const T* izcrd = &snapr.zcrd[i * padded_natom];
      for (size_t j = 0; j < i; j++) {
        const T* jxcrd = &snapr.xcrd[j * padded_natom];
        const T* jycrd = &snapr.ycrd[j * padded_natom];
        const T* jzcrd = &snapr.zcrd[j * padded_natom];
        res_ptr[((i - 1LLU) * i / 2LLU) + j] = rmsd(ixcrd, iycrd, izcrd, jxcrd, jycrd, jzcrd,
                                                    mass_ptr, rplan.getGeneralStrategy(), 0, natom,
                                                    snapr.inv_gpos_scale);
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const PhaseSpaceSynthesis &snapshots,
          const Hybrid<int> &reference_frames, Hybrid<T> *result,
          const SystemGrouping organization) {
  checkCompatibility(cg, rplan, snapshots, result, RMSDTask::REFERENCE, organization);
  const PsSynthesisReader snapr = snapshots.data();
  const RMSDPlanReader<double> rplanr = rplan.dpData();
  T* res_ptr = result->data();
  const int n_result_groups = cg.getPartitionCount(organization);
  for (int i = 0; i < n_result_groups; i++) {
    const std::vector<int> frames_in_group = cg.getGroupStructureIndices(i, organization);
    const int nfrm = frames_in_group.size();
    if (nfrm == 0) {
      continue;
    }
    const int plan_index = rplan.getPlanIndex(i, organization);
    const int plan_offset = cg.getAllToOneResultOffset(i, organization);
    const double *mass_ptr = &rplanr.masses[rplanr.atom_starts[plan_index]];
    const CoordinateFrame ref_i = snapshots.exportCoordinates(reference_frames.readHost(i));
    const CoordinateFrameReader ref_ir = ref_i.data();
    switch (rplan.getAlignmentProtocol(plan_index)) {
    case RMSDAlignmentProtocol::BUILD_CORE:
    case RMSDAlignmentProtocol::ALIGN_CORE:
    case RMSDAlignmentProtocol::ALIGN_ALL:
      for (int j = 0; j < nfrm; j++) {
        const CoordinateFrame cf_ij = snapshots.exportCoordinates(frames_in_group[j]);
        const CoordinateFrameReader cf_ijr = cf_ij.data();
        res_ptr[plan_offset + j] = rmsd<double, double>(ref_ir.xcrd, ref_ir.ycrd, ref_ir.zcrd,
                                                        cf_ijr.xcrd, cf_ijr.ycrd, cf_ijr.zcrd,
                                                        mass_ptr, rplanr.strategy, 0,
                                                        ref_ir.natom);
      }
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const PhaseSpaceSynthesis &snapshots,
          Hybrid<T> *result, const SystemGrouping organization) {
  checkCompatibility(cg, rplan, snapshots, result, RMSDTask::MATRIX, organization);
  const PsSynthesisReader snapr = snapshots.data();
  const RMSDPlanReader<double> rplanr = rplan.dpData();
  T* res_ptr = result->data();
  const int n_result_groups = cg.getPartitionCount(organization);
  for (int i = 0; i < n_result_groups; i++) {
    const std::vector<int> frames_in_group = cg.getGroupStructureIndices(i, organization);
    const int nfrm = frames_in_group.size();
    if (nfrm == 0) {
      continue;
    }
    const int plan_index = rplan.getPlanIndex(i, organization);
    const size_t plan_offset = cg.getSymmetryEquivalentResultOffset(i, organization);
    const double *mass_ptr = &rplanr.masses[rplanr.atom_starts[plan_index]];
    switch (rplan.getAlignmentProtocol(plan_index)) {
    case RMSDAlignmentProtocol::BUILD_CORE:
    case RMSDAlignmentProtocol::ALIGN_CORE:
    case RMSDAlignmentProtocol::ALIGN_ALL:
      for (size_t j = 1; j < nfrm; j++) {
        const CoordinateFrame cf_ij = snapshots.exportCoordinates(frames_in_group[j]);
        const CoordinateFrameReader cf_ijr = cf_ij.data();
        for (size_t k = 0; k < j; k++) {
          const CoordinateFrame cf_ik = snapshots.exportCoordinates(frames_in_group[k]);
          const CoordinateFrameReader cf_ikr = cf_ik.data();
          const size_t result_idx = plan_offset + ((j - 1) * j / 2) + k;
          res_ptr[result_idx] = rmsd<double, double>(cf_ijr.xcrd, cf_ijr.ycrd, cf_ijr.zcrd,
                                                     cf_ikr.xcrd, cf_ikr.ycrd, cf_ikr.zcrd,
                                                     mass_ptr, rplanr.strategy, 0, cf_ijr.natom);
        }
      }
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const PhaseSpaceSynthesis &snapshots,
          const Condensate &cdns, const Hybrid<int> &reference_frames, Hybrid<T> *result,
          const SystemGrouping organization) {
  checkCompatibility(cg, rplan, snapshots, cdns, result, RMSDTask::REFERENCE, organization);
  T* res_ptr = result->data();
  const PsSynthesisReader snapr = snapshots.data();
  const CondensateReader cdr = cdns.data();
  

  // Obtain both RMSD plan abstracts, anticipating one to match either mode for the Condensate
  const RMSDPlanReader<double> rplandr = rplan.dpData();
  const RMSDPlanReader<float> rplanfr = rplan.spData();
  const SynthesisMapReader scmapr = rplan.getSynthesisMapPointer()->data();
  const std::vector<int4> atr_insr = cg.getATRInstructionMembers(organization);
  const std::vector<int> atr_group = cg.getATRInstructionGroups(organization);
  const int instruction_count = cg.getATRInstructionCount(organization);
  for (int i = 0; i < instruction_count; i++) {

    // Rather than check the plan for the number of frames, check the Condensate for the work
    // unit.  On the CPU, this will indicate an amount of work that spans all replicas for each
    // system, but in parallel computating situations (HPC) the work unit may indicate a subset
    // of the systems to work on.
    const int4 tinsr = atr_insr[i];
    const int  tgrp  = atr_group[i];
    const size_t result_offset = cg.getAllToOneResultOffset(tgrp, organization);
    const int ref_idx = reference_frames.readHost(tinsr.y);
    const size_t ref_crd_idx = snapr.atom_starts[ref_idx];

    // Assign all pointers at this level to accept their const-ness.  Adding a switch would
    // require that the pointers be cast as non-const, or that they be inherently confined to
    // the scope of one case in the switch.
    const double *ref_xcrd = &cdr.xcrd[ref_crd_idx];
    const double *ref_ycrd = &cdr.ycrd[ref_crd_idx];
    const double *ref_zcrd = &cdr.zcrd[ref_crd_idx];
    const float *ref_xcrd_sp = &cdr.xcrd_sp[ref_crd_idx];
    const float *ref_ycrd_sp = &cdr.ycrd_sp[ref_crd_idx];
    const float *ref_zcrd_sp = &cdr.zcrd_sp[ref_crd_idx];
    const double *ref_masses = &rplandr.masses[rplandr.atom_starts[tinsr.z]];
    const float *ref_masses_sp = &rplanfr.masses[rplandr.atom_starts[tinsr.z]];
    const int ni_atom = snapr.atom_counts[ref_idx];
    for (int j = 0; j < tinsr.w; j++) {
      int system_idx;
      switch (organization) {
      case SystemGrouping::SOURCE:
        system_idx = scmapr.csystem_proj[scmapr.csystem_bounds[tgrp] + tinsr.x + j];
        break;
      case SystemGrouping::TOPOLOGY:
        system_idx = scmapr.ctopol_proj[scmapr.ctopol_bounds[tgrp] + tinsr.x + j];
        break;
      case SystemGrouping::LABEL:
        system_idx = scmapr.clabel_proj[scmapr.clabel_bounds[tgrp] + tinsr.x + j];
        break;
      }
      const size_t system_crd_idx = snapr.atom_starts[system_idx];
      const size_t result_idx = static_cast<size_t>(tinsr.x + j) + result_offset;
      switch (cdr.mode) {
      case PrecisionModel::DOUBLE:
        {
          const double *system_xcrd = &cdr.xcrd[system_crd_idx];
          const double *system_ycrd = &cdr.ycrd[system_crd_idx];
          const double *system_zcrd = &cdr.zcrd[system_crd_idx];
          res_ptr[result_idx] = rmsd<double, double>(ref_xcrd, ref_ycrd, ref_zcrd, system_xcrd,
                                                     system_ycrd, system_zcrd, ref_masses,
                                                     rplandr.strategy, 0, ni_atom);
        }
        break;
      case PrecisionModel::SINGLE:
        {
          const float *system_xcrd_sp = &cdr.xcrd_sp[system_crd_idx];
          const float *system_ycrd_sp = &cdr.ycrd_sp[system_crd_idx];
          const float *system_zcrd_sp = &cdr.zcrd_sp[system_crd_idx];
          res_ptr[result_idx] = rmsd<float, float>(ref_xcrd_sp, ref_ycrd_sp, ref_zcrd_sp,
                                                   system_xcrd_sp, system_ycrd_sp,
                                                   system_zcrd_sp, ref_masses_sp,
                                                   rplanfr.strategy, 0, ni_atom);
        }
        break;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const PhaseSpaceSynthesis &snapshots,
          const Condensate &cdns, Hybrid<T> *result, const SystemGrouping organization) {
  checkCompatibility(cg, rplan, snapshots, cdns, result, RMSDTask::MATRIX, organization);
  T* res_ptr = result->data();
  const PsSynthesisReader snapr = snapshots.data();
  const CondensateReader cdr = cdns.data();
  const int half_int_bits = sizeof(uint) * 4;
  int half_mask = 0;
  for (int i = 0; i < half_int_bits; i++) {
    half_mask |= (0x1 << i);
  }
  
  // Again, obtain both RMSD plan abstracts, and always use them rather than getter functions
  // querying the underlying object in the inner loops to mimic GPU kernel activity.
  const RMSDPlanReader<double> rplandr = rplan.dpData();
  const RMSDPlanReader<float> rplanfr = rplan.spData();
  const SynthesisMapReader scmapr = rplan.getSynthesisMapPointer()->data();
  const std::vector<int4> ata_insr = cg.getATAInstructionMembers(organization);
  const std::vector<int> ata_group = cg.getATAInstructionGroups(organization);
  const int instruction_count = cg.getATAInstructionCount(organization);
  for (int i = 0; i < instruction_count; i++) {
    const int4 tinsr = ata_insr[i];
    const int  tgrp  = ata_group[i];
    const size_t result_offset = cg.getSymmetryEquivalentResultOffset(tgrp, organization);
    const double *i_masses   = &rplandr.masses[rplandr.atom_starts[tinsr.z]];
    const float *i_masses_sp = &rplanfr.masses[rplanfr.atom_starts[tinsr.z]];
    const int jside = (tinsr.w & half_mask);
    const int kside = ((tinsr.w >> half_int_bits) & half_mask);
    size_t inst_offset;
    switch (organization) {
    case SystemGrouping::SOURCE:
      inst_offset = scmapr.csystem_bounds[tgrp];
      break;
    case SystemGrouping::TOPOLOGY:
      inst_offset = scmapr.ctopol_bounds[tgrp];
      break;
    case SystemGrouping::LABEL:
      inst_offset = scmapr.clabel_bounds[tgrp];
      break;
    }
    for (int j = 0; j < jside; j++) {
      const size_t jrepno = tinsr.x + j;
      int jrep_idx;
      switch (organization) {
      case SystemGrouping::SOURCE:
        jrep_idx = scmapr.csystem_proj[inst_offset + jrepno];
        break;
      case SystemGrouping::TOPOLOGY:
        jrep_idx = scmapr.ctopol_proj[inst_offset + jrepno];
        break;
      case SystemGrouping::LABEL:
        jrep_idx = scmapr.clabel_proj[inst_offset + jrepno];
        break;
      }

      // Obtain pointers to the coordinates of the "reference" frame for this row of the matrix
      const size_t jrep_crd_idx = snapr.atom_starts[jrep_idx];
      const int njk_atom = snapr.atom_counts[jrep_idx];

      // Again, assign all pointers at this level to accept their const-ness.
      const double *jxcrd = &cdr.xcrd[jrep_crd_idx];
      const double *jycrd = &cdr.ycrd[jrep_crd_idx];
      const double *jzcrd = &cdr.zcrd[jrep_crd_idx];
      const float *jxcrd_sp = &cdr.xcrd_sp[jrep_crd_idx];
      const float *jycrd_sp = &cdr.ycrd_sp[jrep_crd_idx];
      const float *jzcrd_sp = &cdr.zcrd_sp[jrep_crd_idx];
      for (int k = 0; k < kside; k++) {
        const size_t krepno = tinsr.y + k;
        int krep_idx;
        switch (organization) {
        case SystemGrouping::SOURCE:
          krep_idx = scmapr.csystem_proj[inst_offset + krepno];
          break;
        case SystemGrouping::TOPOLOGY:
          krep_idx = scmapr.ctopol_proj[inst_offset + krepno];
          break;
        case SystemGrouping::LABEL:
          krep_idx = scmapr.clabel_proj[inst_offset + krepno];
          break;
        }

        // Skip the calculation if outside the lower triangular region.
        if (krepno < jrepno) {
          const size_t krep_crd_idx = snapr.atom_starts[krep_idx];
          const size_t result_idx = (jrepno * (jrepno - 1) / 2) + krepno + result_offset;
          switch (cdr.mode) {
          case PrecisionModel::DOUBLE:
            {
              const double *kxcrd = &cdr.xcrd[krep_crd_idx];
              const double *kycrd = &cdr.ycrd[krep_crd_idx];
              const double *kzcrd = &cdr.zcrd[krep_crd_idx];
              res_ptr[result_idx] = rmsd<double, double>(jxcrd, jycrd, jzcrd, kxcrd, kycrd,
                                                         kzcrd, i_masses, rplandr.strategy, 0,
                                                         njk_atom);
            }
            break;
          case PrecisionModel::SINGLE:
            {
              const float *kxcrd_sp = &cdr.xcrd_sp[krep_crd_idx];
              const float *kycrd_sp = &cdr.ycrd_sp[krep_crd_idx];
              const float *kzcrd_sp = &cdr.zcrd_sp[krep_crd_idx];
              res_ptr[result_idx] = rmsd<float, float>(jxcrd_sp, jycrd_sp, jzcrd_sp, kxcrd_sp,
                                                       kycrd_sp, kzcrd_sp, i_masses_sp,
                                                       rplanfr.strategy, 0, njk_atom);
            }
            break;
          }
        }
      }
    }
  }
}

} // namespace structure
} // namespace stormm
