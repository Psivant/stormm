// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void rotateAboutBond(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, const int atom_i, const int atom_j,
                     const int* moving_atoms, const int moving_atom_count,
                     const Tcalc rotation_angle, const Tcalc globalpos_scale_factor) {

  // Center and image the coordinates
  const Tcoord center_x = xcrd[atom_j];
  const Tcoord center_y = ycrd[atom_j];
  const Tcoord center_z = zcrd[atom_j];
  xcrd[atom_i] -= center_x;
  ycrd[atom_i] -= center_y;
  zcrd[atom_i] -= center_z;
  xcrd[atom_j] = static_cast<Tcoord>(0.0);
  ycrd[atom_j] = static_cast<Tcoord>(0.0);
  zcrd[atom_j] = static_cast<Tcoord>(0.0);
  for (int i = 0; i < moving_atom_count; i++) {
    const int mk = moving_atoms[i];
    xcrd[mk] -= center_x;
    ycrd[mk] -= center_y;
    zcrd[mk] -= center_z;
  }
  
  // Define the vector of rotation, then the matrix.  The coordinate data type is assumed to be
  // either a real numbered type (float, double), or a signed integral type (int, long long int),
  // a fact which will be enforced by the CoordinateSeries object, the only templated coordinate
  // representation.  All other coordinate representations work in double or long long int.
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  Tcalc dx, dy, dz, invdr, cos_ra, sin_ra;
  const Tcalc value_one = 1.0;
  const Tcalc inv_globalpos_scale_factor = value_one / globalpos_scale_factor;  
  if (tcoord_is_sgnint) {
    dx = static_cast<Tcalc>(xcrd[atom_j] - xcrd[atom_i]) * inv_globalpos_scale_factor;
    dy = static_cast<Tcalc>(ycrd[atom_j] - ycrd[atom_i]) * inv_globalpos_scale_factor;
    dz = static_cast<Tcalc>(zcrd[atom_j] - zcrd[atom_i]) * inv_globalpos_scale_factor;
  }
  else {
    dx = xcrd[atom_j] - xcrd[atom_i];
    dy = ycrd[atom_j] - ycrd[atom_i];
    dz = zcrd[atom_j] - zcrd[atom_i];
  }
  if (tcalc_is_double) {
    invdr = value_one / sqrt((dx * dx) + (dy * dy) + (dz * dz));
    cos_ra = cos(rotation_angle);
    sin_ra = sin(rotation_angle);
  }
  else {
    invdr = value_one / sqrtf((dx * dx) + (dy * dy) + (dz * dz));
    cos_ra = cosf(rotation_angle);
    sin_ra = sinf(rotation_angle);
  }
  dx *= invdr;
  dy *= invdr;
  dz *= invdr;
  const std::vector<Tcalc> rmat = { cos_ra + (dx * dx * (value_one - cos_ra)),
                                    (dy * dx * (value_one - cos_ra)) + (dz * sin_ra),
                                    (dz * dx * (value_one - cos_ra)) - (dy * sin_ra),
                                    (dx * dy * (value_one - cos_ra)) - (dz * sin_ra),
                                    cos_ra + (dy * dy * (value_one - cos_ra)),
                                    (dz * dy * (value_one - cos_ra)) + (dx * sin_ra),
                                    (dx * dz * (value_one - cos_ra)) + (dy * sin_ra),
                                    (dy * dz * (value_one - cos_ra)) - (dx * sin_ra),
                                    cos_ra + (dz * dz * (value_one - cos_ra)) };
  
  // Restore the original imaging and location of the bond atoms
  xcrd[atom_j] = center_x;
  ycrd[atom_j] = center_y;
  zcrd[atom_j] = center_z;
  xcrd[atom_i] += center_x;
  ycrd[atom_i] += center_y;
  zcrd[atom_i] += center_z;

  // Loop over all moving particles, rotate about the vector, and re-apply their original
  // translations relative to the origin.
  for (int i = 0; i < moving_atom_count; i++) {
    const int mk = moving_atoms[i];
    const Tcalc nx = (rmat[0] * xcrd[mk]) + (rmat[3] * ycrd[mk]) + (rmat[6] * zcrd[mk]);
    const Tcalc ny = (rmat[1] * xcrd[mk]) + (rmat[4] * ycrd[mk]) + (rmat[7] * zcrd[mk]);
    const Tcalc nz = (rmat[2] * xcrd[mk]) + (rmat[5] * ycrd[mk]) + (rmat[8] * zcrd[mk]);
    if (tcoord_is_sgnint) {
      xcrd[mk] = llround(nx * globalpos_scale_factor);
      ycrd[mk] = llround(ny * globalpos_scale_factor);
      zcrd[mk] = llround(nz * globalpos_scale_factor);
    }
    else {
      xcrd[mk] = nx;
      ycrd[mk] = ny;
      zcrd[mk] = nz;
    }
    xcrd[mk] += center_x;
    ycrd[mk] += center_y;
    zcrd[mk] += center_z;    
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void rotateAboutBond(CoordinateSeries<Tcoord> *cs, const int frame_index, const int atom_i,
                     const int atom_j, const std::vector<int> &moving_atoms,
                     const Tcalc rotation_angle) {
  rotateAboutBond<Tcoord, Tcalc>(cs->data(), frame_index, atom_i, atom_j, moving_atoms.data(),
                                 moving_atoms.size(), rotation_angle);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void rotateAboutBond(CoordinateSeriesWriter<Tcoord> csw, const int frame_index, const int atom_i,
                     const int atom_j, const int* moving_atoms, const int moving_atom_count,
                     const Tcalc rotation_angle) {
  const size_t fidx_zu  = static_cast<size_t>(frame_index);
  const size_t natom_zu = roundUp<size_t>(csw.natom, warp_size_zu);
  rotateAboutBond<Tcoord, Tcalc>(&csw.xcrd[fidx_zu * natom_zu], &csw.ycrd[fidx_zu * natom_zu],
                                 &csw.zcrd[fidx_zu * natom_zu], atom_i, atom_j, moving_atoms,
                                 moving_atom_count, rotation_angle, csw.gpos_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void rotateAboutBond(PhaseSpaceSynthesis *psynth, const int system_index, const int atom_i,
                     const int atom_j, const std::vector<int> &moving_atoms,
                     const Tcalc rotation_angle) {
  rotateAboutBond<Tcalc>(psynth->data(), system_index, atom_i, atom_j, moving_atoms.data(),
                         moving_atoms.size(), rotation_angle);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void rotateAboutBond(PsSynthesisWriter psynthw, const int system_index, const int atom_i,
                     const int atom_j, const int* moving_atoms, const int moving_atom_count,
                     const Tcalc rotation_angle) {
  const int system_offset = psynthw.atom_starts[system_index];
  std::vector<int> local_atoms(moving_atom_count);
  for (int i = 0; i < moving_atom_count; i++) {
    local_atoms[i] = i + 2;
  }
  std::vector<Tcalc> xcrd(moving_atom_count + 2);
  std::vector<Tcalc> ycrd(moving_atom_count + 2);
  std::vector<Tcalc> zcrd(moving_atom_count + 2);
  const Tcalc inv_scl = psynthw.inv_gpos_scale;
  if (psynthw.gpos_bits <= globalpos_scale_nonoverflow_bits) {
    xcrd[0] = static_cast<Tcalc>(psynthw.xcrd[atom_i + system_offset]) * inv_scl;
    ycrd[0] = static_cast<Tcalc>(psynthw.ycrd[atom_i + system_offset]) * inv_scl;
    zcrd[0] = static_cast<Tcalc>(psynthw.zcrd[atom_i + system_offset]) * inv_scl;
    xcrd[1] = static_cast<Tcalc>(psynthw.xcrd[atom_j + system_offset]) * inv_scl;
    ycrd[1] = static_cast<Tcalc>(psynthw.ycrd[atom_j + system_offset]) * inv_scl;
    zcrd[1] = static_cast<Tcalc>(psynthw.zcrd[atom_j + system_offset]) * inv_scl;
    for (int i = 0; i < moving_atom_count; i++) {
      const int apos_n = moving_atoms[i] + system_offset;
      xcrd[i + 2] = static_cast<Tcalc>(psynthw.xcrd[apos_n]) * inv_scl;
      ycrd[i + 2] = static_cast<Tcalc>(psynthw.ycrd[apos_n]) * inv_scl;
      zcrd[i + 2] = static_cast<Tcalc>(psynthw.zcrd[apos_n]) * inv_scl;
    }
  }
  else {
    const size_t apos_i = atom_i + system_offset;
    const size_t apos_j = atom_j + system_offset;
    xcrd[0] = hostInt95ToDouble(psynthw.xcrd[apos_i], psynthw.xcrd_ovrf[apos_i]) * inv_scl;
    ycrd[0] = hostInt95ToDouble(psynthw.ycrd[apos_i], psynthw.ycrd_ovrf[apos_i]) * inv_scl;
    zcrd[0] = hostInt95ToDouble(psynthw.zcrd[apos_i], psynthw.zcrd_ovrf[apos_i]) * inv_scl;
    xcrd[1] = hostInt95ToDouble(psynthw.xcrd[apos_j], psynthw.xcrd_ovrf[apos_j]) * inv_scl;
    ycrd[1] = hostInt95ToDouble(psynthw.ycrd[apos_j], psynthw.ycrd_ovrf[apos_j]) * inv_scl;
    zcrd[1] = hostInt95ToDouble(psynthw.zcrd[apos_j], psynthw.zcrd_ovrf[apos_j]) * inv_scl;
    for (int i = 0; i < moving_atom_count; i++) {
      const int apos_n = moving_atoms[i] + system_offset;
      xcrd[i + 2] = hostInt95ToDouble(psynthw.xcrd[apos_n], psynthw.xcrd_ovrf[apos_n]) * inv_scl;
      ycrd[i + 2] = hostInt95ToDouble(psynthw.ycrd[apos_n], psynthw.ycrd_ovrf[apos_n]) * inv_scl;
      zcrd[i + 2] = hostInt95ToDouble(psynthw.zcrd[apos_n], psynthw.zcrd_ovrf[apos_n]) * inv_scl;
    }
  }
  rotateAboutBond<Tcalc>(xcrd.data(), ycrd.data(), zcrd.data(), 0, 1, local_atoms.data(),
                         local_atoms.size(), rotation_angle);
  const Tcalc pos_scl = psynthw.gpos_scale;
  if (psynthw.gpos_bits <= globalpos_scale_nonoverflow_bits) {
    for (int i = 0; i < moving_atom_count; i++) {
      const int apos_n = moving_atoms[i] + system_offset;
      psynthw.xcrd[apos_n] = llround(xcrd[i + 2] * pos_scl);
      psynthw.ycrd[apos_n] = llround(ycrd[i + 2] * pos_scl);
      psynthw.zcrd[apos_n] = llround(zcrd[i + 2] * pos_scl);
    }
  }
  else {
    for (int i = 0; i < moving_atom_count; i++) {
      const int apos_n = moving_atoms[i] + system_offset;
      const int95_t nx = hostDoubleToInt95(xcrd[i + 2] * pos_scl);
      const int95_t ny = hostDoubleToInt95(ycrd[i + 2] * pos_scl);
      const int95_t nz = hostDoubleToInt95(zcrd[i + 2] * pos_scl);
      psynthw.xcrd[apos_n] = nx.x;
      psynthw.ycrd[apos_n] = ny.x;
      psynthw.zcrd[apos_n] = nz.x;
      psynthw.xcrd_ovrf[apos_n] = nx.y;
      psynthw.ycrd_ovrf[apos_n] = ny.y;
      psynthw.zcrd_ovrf[apos_n] = nz.y;
    }
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void rotateAboutBond(Condensate *cdns, const int system_index, const int atom_i, const int atom_j,
                     const std::vector<int> &moving_atoms, const Tcalc rotation_angle) {
  rotateAboutBond<Tcalc>(cdns->data(), system_index, atom_i, atom_j, moving_atoms,
                         rotation_angle);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void rotateAboutBond(CondensateWriter cdnsw, const int system_index, const int atom_i,
                     const int atom_j, const std::vector<int> &moving_atoms,
                     const Tcalc rotation_angle) {
  const size_t atom_offset = cdnsw.atom_starts[system_index];
  switch (cdnsw.mode) {
  case PrecisionModel::DOUBLE:
    rotateAboutBond(&cdnsw.xcrd[atom_offset], &cdnsw.ycrd[atom_offset], &cdnsw.zcrd[atom_offset],
                    atom_i, atom_j, moving_atoms.data(), moving_atoms.size(), rotation_angle);
    break;
  case PrecisionModel::SINGLE:
    rotateAboutBond(&cdnsw.xcrd_sp[atom_offset], &cdnsw.ycrd_sp[atom_offset],
                    &cdnsw.zcrd_sp[atom_offset], atom_i, atom_j, moving_atoms.data(),
                    moving_atoms.size(), rotation_angle);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void flipChiralCenter(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, const int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups,
                      const Tcalc globalpos_scale_factor) {
  const size_t tcalc_ct = std::type_index(typeid(Tcalc)).hash_code();
  const bool tcalc_is_double = (tcalc_ct == double_type_index);
  const bool tcoord_is_sgnint = isSignedIntegralScalarType<Tcoord>();
  const Tcalc value_one = 1.0;
  const Tcalc inv_globalpos_scale_factor = value_one / globalpos_scale_factor;
  switch (chiral_protocols[center_idx]) {
  case ChiralInversionProtocol::ROTATE:
    {
      // Unpack the appropriate inversion group
      const int root_idx = inversion_groups[center_idx].getRootAtom();
      const int pivt_idx = inversion_groups[center_idx].getPivotAtom();
      const int ccen   = chiral_centers[center_idx];

      // Find the bisector of the root_idx : chiral_center : pivt_idx angle.  Shift the pivt_idx
      // atom to lie along the line of the bisector, rotate the moving atoms 180 degrees about
      // this "bond," then replace the pivt_idx atom.
      const Tcoord orig_bx = xcrd[pivt_idx];
      const Tcoord orig_by = ycrd[pivt_idx];
      const Tcoord orig_bz = zcrd[pivt_idx];
      Tcalc dax, day, daz, dbx, dby, dbz, ccenx, cceny, ccenz, invra, invrb;
      if (tcoord_is_sgnint) {
        dax = static_cast<Tcalc>(xcrd[root_idx] - xcrd[ccen]) * inv_globalpos_scale_factor;
        day = static_cast<Tcalc>(ycrd[root_idx] - ycrd[ccen]) * inv_globalpos_scale_factor;
        daz = static_cast<Tcalc>(zcrd[root_idx] - zcrd[ccen]) * inv_globalpos_scale_factor;
        dbx = static_cast<Tcalc>(orig_bx - xcrd[ccen]) * inv_globalpos_scale_factor;
        dby = static_cast<Tcalc>(orig_by - ycrd[ccen]) * inv_globalpos_scale_factor;
        dbz = static_cast<Tcalc>(orig_bz - zcrd[ccen]) * inv_globalpos_scale_factor;
        ccenx = static_cast<Tcalc>(xcrd[ccen]) * inv_globalpos_scale_factor;
        cceny = static_cast<Tcalc>(ycrd[ccen]) * inv_globalpos_scale_factor;
        ccenz = static_cast<Tcalc>(zcrd[ccen]) * inv_globalpos_scale_factor;
      }
      else {

        // Make no assumptions about the precision of the calculation relative to the precision
        // of the coordinates.  It is possible, and even common in MD programs, to store
        // coordinates in a higher degree of precision than the calculation.  It would be faster,
        // otherwise, to set ccen(x,y,z) at the front and re-use the values in computing da or db.
        dax = xcrd[root_idx] - xcrd[ccen];
        day = ycrd[root_idx] - ycrd[ccen];
        daz = zcrd[root_idx] - zcrd[ccen];
        dbx = orig_bx - xcrd[ccen];
        dby = orig_by - ycrd[ccen];
        dbz = orig_bz - zcrd[ccen];
        ccenx = xcrd[ccen];
        cceny = ycrd[ccen];
        ccenz = zcrd[ccen];
      }
      if (tcalc_is_double) {
        invrb = 1.0 / sqrt((dbx * dbx) + (dby * dby) + (dbz * dbz));
        invra = 1.0 / sqrt((dax * dax) + (day * day) + (daz * daz));
      }
      else {
        invrb = 1.0f / sqrtf((dbx * dbx) + (dby * dby) + (dbz * dbz));
        invra = 1.0f / sqrtf((dax * dax) + (day * day) + (daz * daz));
      }
      dbx = ccenx + (dbx * invrb);
      dby = cceny + (dby * invrb);
      dbz = ccenz + (dbz * invrb);
      dax = ccenx + (dax * invra);
      day = cceny + (day * invra);
      daz = ccenz + (daz * invra);
      const Tcalc midpoint_x = 0.5 * (dbx + dax);
      const Tcalc midpoint_y = 0.5 * (dby + day);
      const Tcalc midpoint_z = 0.5 * (dbz + daz);
      if (tcoord_is_sgnint) {
        xcrd[pivt_idx] = llround(midpoint_x * globalpos_scale_factor);
        ycrd[pivt_idx] = llround(midpoint_y * globalpos_scale_factor);
        zcrd[pivt_idx] = llround(midpoint_z * globalpos_scale_factor);
      }
      else {
        xcrd[pivt_idx] = midpoint_x;
        ycrd[pivt_idx] = midpoint_y;
        zcrd[pivt_idx] = midpoint_z;
      }
      rotateAboutBond(xcrd, ycrd, zcrd, pivt_idx, chiral_centers[center_idx],
                      inversion_groups[center_idx].getMovingAtoms().data(),
                      inversion_groups[center_idx].getMovingAtoms().size(),
                      static_cast<Tcalc>(symbols::pi), globalpos_scale_factor);
      xcrd[pivt_idx] = orig_bx;
      ycrd[pivt_idx] = orig_by;
      zcrd[pivt_idx] = orig_bz;
    }
    break;
  case ChiralInversionProtocol::REFLECT:
    {
      // Find the molecule home of the present center and flip those atoms only.
      const int natom = inversion_groups[center_idx].getMovingAtomCount();
      const int* moving_atoms = inversion_groups[center_idx].getMovingAtoms().data();
      for (int i = 0; i < natom; i++) {
        const int atom_idx = moving_atoms[i];
        xcrd[atom_idx] = -xcrd[atom_idx];
      }

      // All centers have been flipped by the reflection.  Loop over all other chiral centers and
      // perform chiral inversions in order to flip the others back, where possible.  Infinite
      // recursion is limited by the fact that at most one chiral center in a molecule can be given
      // the protocol "REFLECT," and reflection only triggers subsequent rotations.
      const int nchirals = chiral_protocols.size();
      for (int i = 0; i < nchirals; i++) {
        if (chiral_protocols[i] == ChiralInversionProtocol::ROTATE) {
          flipChiralCenter(xcrd, ycrd, zcrd, i, chiral_centers, chiral_protocols,
                           inversion_groups, globalpos_scale_factor);
        }
      }
    }
    break;
  case ChiralInversionProtocol::DO_NOT_INVERT:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void flipChiralCenter(CoordinateSeries<Tcoord> *cs, const int frame_index,
                      const int center_idx, const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups) {
  flipChiralCenter<Tcoord, Tcalc>(cs->data(), frame_index, center_idx, chiral_centers,
                                  chiral_protocols, inversion_groups);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tcalc>
void flipChiralCenter(CoordinateSeriesWriter<Tcoord> csw, const int frame_index,
                      const int center_idx, const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups) {
  const size_t fidx_zu  = static_cast<size_t>(frame_index);
  const size_t natom_zu = roundUp<size_t>(csw.natom, warp_size_zu);
  flipChiralCenter<Tcoord, Tcalc>(&csw.xcrd[fidx_zu * natom_zu], &csw.ycrd[fidx_zu * natom_zu],
                                  &csw.zcrd[fidx_zu * natom_zu], center_idx, chiral_centers,
                                  chiral_protocols, inversion_groups, csw.gpos_scale);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void flipChiralCenter(PsSynthesisWriter psynthw, const int system_index, const int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups) {

  // Prepare an array of local atom indices to indicate moving atoms.
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  switch (chiral_protocols[center_idx]) {
  case ChiralInversionProtocol::ROTATE:
    {
      const int system_offset = psynthw.atom_starts[system_index];
      const int nmove = inversion_groups[center_idx].getMovingAtomCount();
      std::vector<int> local_atoms(nmove);
      for (int i = 0; i < nmove; i++) {
        local_atoms[i] = i + 2;
      }
      std::vector<Tcalc> xcrd(nmove + 2), ycrd(nmove + 2), zcrd(nmove + 2);
      const Tcalc inv_scl = psynthw.inv_gpos_scale;
      const Tcalc value_half = 0.5;
      const Tcalc value_one  = 1.0;
      
      // Prepare the appropriate representations of the relevant coordinates, including the
      // midpoint of the two root atoms for the heaviest chiral branches and the chiral atom
      // itself.
      const int root_idx = inversion_groups[center_idx].getRootAtom() + system_offset;
      const int pivt_idx = inversion_groups[center_idx].getPivotAtom() + system_offset;
      const int chir_idx = chiral_centers[center_idx] + system_offset;
      Tcalc ccenx, cceny, ccenz, dax, day, daz, dbx, dby, dbz;
      if (psynthw.gpos_bits <= globalpos_scale_nonoverflow_bits) {
        ccenx = static_cast<Tcalc>(psynthw.xcrd[chir_idx]) * inv_scl;
        cceny = static_cast<Tcalc>(psynthw.ycrd[chir_idx]) * inv_scl;
        ccenz = static_cast<Tcalc>(psynthw.zcrd[chir_idx]) * inv_scl;
        dax = static_cast<Tcalc>(psynthw.xcrd[root_idx] - psynthw.xcrd[chir_idx]) * inv_scl;
        day = static_cast<Tcalc>(psynthw.ycrd[root_idx] - psynthw.ycrd[chir_idx]) * inv_scl;
        daz = static_cast<Tcalc>(psynthw.zcrd[root_idx] - psynthw.zcrd[chir_idx]) * inv_scl;
        dbx = static_cast<Tcalc>(psynthw.xcrd[pivt_idx] - psynthw.xcrd[chir_idx]) * inv_scl;
        dby = static_cast<Tcalc>(psynthw.ycrd[pivt_idx] - psynthw.ycrd[chir_idx]) * inv_scl;
        dbz = static_cast<Tcalc>(psynthw.zcrd[pivt_idx] - psynthw.zcrd[chir_idx]) * inv_scl;
      }
      else {
        const int95_t tmp_ccen_x  = { psynthw.xcrd[chir_idx], psynthw.xcrd_ovrf[chir_idx] };
        const int95_t tmp_ccen_y  = { psynthw.ycrd[chir_idx], psynthw.ycrd_ovrf[chir_idx] };
        const int95_t tmp_ccen_z  = { psynthw.zcrd[chir_idx], psynthw.zcrd_ovrf[chir_idx] };
        const int95_t tmp_roota_x = { psynthw.xcrd[root_idx], psynthw.xcrd_ovrf[root_idx] };
        const int95_t tmp_roota_y = { psynthw.ycrd[root_idx], psynthw.ycrd_ovrf[root_idx] };
        const int95_t tmp_roota_z = { psynthw.zcrd[root_idx], psynthw.zcrd_ovrf[root_idx] };
        const int95_t tmp_rootb_x = { psynthw.xcrd[pivt_idx], psynthw.xcrd_ovrf[pivt_idx] };
        const int95_t tmp_rootb_y = { psynthw.ycrd[pivt_idx], psynthw.ycrd_ovrf[pivt_idx] };
        const int95_t tmp_rootb_z = { psynthw.zcrd[pivt_idx], psynthw.zcrd_ovrf[pivt_idx] };
        ccenx = hostInt95ToDouble(tmp_ccen_x) * inv_scl;
        cceny = hostInt95ToDouble(tmp_ccen_y) * inv_scl;
        ccenz = hostInt95ToDouble(tmp_ccen_z) * inv_scl;
        dax = hostInt95ToDouble(hostSplitFPSubtract(tmp_roota_x, tmp_ccen_x.x, tmp_ccen_x.y)) *
              inv_scl;
        day = hostInt95ToDouble(hostSplitFPSubtract(tmp_roota_y, tmp_ccen_y.x, tmp_ccen_y.y)) *
              inv_scl;
        daz = hostInt95ToDouble(hostSplitFPSubtract(tmp_roota_z, tmp_ccen_z.x, tmp_ccen_z.y)) *
              inv_scl;
        dbx = hostInt95ToDouble(hostSplitFPSubtract(tmp_rootb_x, tmp_ccen_x.x, tmp_ccen_x.y)) *
              inv_scl;
        dby = hostInt95ToDouble(hostSplitFPSubtract(tmp_rootb_y, tmp_ccen_y.x, tmp_ccen_y.y)) *
              inv_scl;
        dbz = hostInt95ToDouble(hostSplitFPSubtract(tmp_rootb_z, tmp_ccen_z.x, tmp_ccen_z.y)) *
              inv_scl;
      }
      Tcalc invra, invrb;
      if (tcalc_is_double) {
        invra = value_one / sqrt((dax * dax) + (day * day) + (daz * daz));
        invrb = value_one / sqrt((dbx * dbx) + (dby * dby) + (dbz * dbz));
      }
      else {
        invra = value_one / sqrtf((dax * dax) + (day * day) + (daz * daz));
        invrb = value_one / sqrtf((dbx * dbx) + (dby * dby) + (dbz * dbz));
      }
      dbx = ccenx + (dbx * invrb);
      dby = cceny + (dby * invrb);
      dbz = ccenz + (dbz * invrb);
      dax = ccenx + (dax * invra);
      day = cceny + (day * invra);
      daz = ccenz + (daz * invra);
      const Tcalc midpoint_x = value_half * (dbx + dax);
      const Tcalc midpoint_y = value_half * (dby + day);
      const Tcalc midpoint_z = value_half * (dbz + daz);
      xcrd[0] = midpoint_x;
      ycrd[0] = midpoint_y;
      zcrd[0] = midpoint_z;
      xcrd[1] = ccenx;
      ycrd[1] = cceny;
      zcrd[1] = ccenz;
      const int* moving_atoms = inversion_groups[center_idx].getMovingAtoms().data();
      if (psynthw.gpos_bits <= globalpos_scale_nonoverflow_bits) {
        for (int i = 0; i < nmove; i++) {
          const int atom_idx = moving_atoms[i] + system_offset;
          xcrd[i + 2] = static_cast<Tcalc>(psynthw.xcrd[atom_idx]) * inv_scl;
          ycrd[i + 2] = static_cast<Tcalc>(psynthw.ycrd[atom_idx]) * inv_scl;
          zcrd[i + 2] = static_cast<Tcalc>(psynthw.zcrd[atom_idx]) * inv_scl;
        }
      }
      else {
        for (int i = 0; i < nmove; i++) {
          const int atom_idx = moving_atoms[i] + system_offset;
          xcrd[i + 2] = hostInt95ToDouble(psynthw.xcrd[atom_idx], psynthw.xcrd_ovrf[atom_idx]) *
                        inv_scl;
          ycrd[i + 2] = hostInt95ToDouble(psynthw.ycrd[atom_idx], psynthw.ycrd_ovrf[atom_idx]) *
                        inv_scl;
          zcrd[i + 2] = hostInt95ToDouble(psynthw.zcrd[atom_idx], psynthw.zcrd_ovrf[atom_idx]) *
                        inv_scl;
        }
      }
      rotateAboutBond<Tcalc>(xcrd.data(), ycrd.data(), zcrd.data(), 0, 1, local_atoms.data(),
                             local_atoms.size(), symbols::pi);
      const Tcalc pos_scl = psynthw.gpos_scale;
      if (psynthw.gpos_bits <= globalpos_scale_nonoverflow_bits) {
        for (int i = 0; i < nmove; i++) {
          const int ixyz_dest = moving_atoms[i] + system_offset;
          psynthw.xcrd[ixyz_dest] = llround(xcrd[i + 2] * pos_scl);
          psynthw.ycrd[ixyz_dest] = llround(ycrd[i + 2] * pos_scl);
          psynthw.zcrd[ixyz_dest] = llround(zcrd[i + 2] * pos_scl);
        }
      }
      else {
        for (int i = 0; i < nmove; i++) {
          const int ixyz_dest = moving_atoms[i] + system_offset;
          const int95_t ixcrd = hostDoubleToInt95(xcrd[i + 2] * pos_scl);
          const int95_t iycrd = hostDoubleToInt95(ycrd[i + 2] * pos_scl);
          const int95_t izcrd = hostDoubleToInt95(zcrd[i + 2] * pos_scl);
          psynthw.xcrd[ixyz_dest] = ixcrd.x;
          psynthw.ycrd[ixyz_dest] = iycrd.x;
          psynthw.zcrd[ixyz_dest] = izcrd.x;
          psynthw.xcrd_ovrf[ixyz_dest] = ixcrd.y;
          psynthw.ycrd_ovrf[ixyz_dest] = iycrd.y;
          psynthw.zcrd_ovrf[ixyz_dest] = izcrd.y;
        }
      }
    }
    break;
  case ChiralInversionProtocol::REFLECT:
    {
      const int system_offset = psynthw.atom_starts[system_index];
      const int natom = inversion_groups[center_idx].getMovingAtomCount();
      const int* moving_atoms = inversion_groups[center_idx].getMovingAtoms().data();
      if (psynthw.gpos_bits <= globalpos_scale_nonoverflow_bits) {
        for (int i = 0; i < natom; i++) {
          const int atom_idx = moving_atoms[i] + system_offset;
          psynthw.xcrd[atom_idx] = -psynthw.xcrd[atom_idx];
        }
      }
      else {
        for (int i = 0; i < natom; i++) {
          const int atom_idx = moving_atoms[i] + system_offset;
          psynthw.xcrd[atom_idx]      = -psynthw.xcrd[atom_idx];
          psynthw.xcrd_ovrf[atom_idx] = -psynthw.xcrd_ovrf[atom_idx];
        }
      }

      // All centers have been flipped by the reflection.  Loop over all other chiral centers and
      // perform chiral inversions in order to flip the others back, where possible.
      const int nchiral = chiral_protocols.size();
      for (int i = 0; i < nchiral; i++) {
        if (chiral_protocols[i] == ChiralInversionProtocol::ROTATE) {
          flipChiralCenter<Tcalc>(psynthw, system_index, i, chiral_centers, chiral_protocols,
                                  inversion_groups);
        }
      }
    }
    break;
  case ChiralInversionProtocol::DO_NOT_INVERT:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void flipChiralCenter(PhaseSpaceSynthesis *psynth, const int system_index, const int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups) {
  flipChiralCenter<Tcalc>(psynth->data(), system_index, center_idx, chiral_centers,
                          chiral_protocols, inversion_groups);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void flipChiralCenter(Condensate *cdns, const int system_index, const int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups) {
  flipChiralCenter<Tcalc>(cdns->data(), system_index, center_idx, chiral_centers, chiral_protocols,
                          inversion_groups);
}

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
void flipChiralCenter(CondensateWriter cdnsw, const int system_index, const int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups) {
  const size_t atom_offset = cdnsw.atom_starts[system_index];
  switch (cdnsw.mode) {
  case PrecisionModel::DOUBLE:
    flipChiralCenter<double, Tcalc>(&cdnsw.xcrd[atom_offset], &cdnsw.ycrd[atom_offset],
                                    &cdnsw.zcrd[atom_offset], center_idx, chiral_centers,
                                    chiral_protocols, inversion_groups);
    break;
  case PrecisionModel::SINGLE:
    flipChiralCenter<float, Tcalc>(&cdnsw.xcrd_sp[atom_offset], &cdnsw.ycrd_sp[atom_offset],
                                   &cdnsw.zcrd_sp[atom_offset], center_idx, chiral_centers,
                                   chiral_protocols, inversion_groups);
    break;
  }
}

} // namespace structure
} // namespace stormm
