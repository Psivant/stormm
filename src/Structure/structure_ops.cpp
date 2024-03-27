#include "copyright.h"
#include "Topology/topology_util.h"
#include "structure_ops.h"

namespace stormm {
namespace structure {

using topology::ChemicalDetailsKit;
using topology::moleculeValidityCheck;
using trajectory::PhaseSpaceReader;
using trajectory::CoordinateFrameReader;
  
//-------------------------------------------------------------------------------------------------
double3 centerOfMass(const AtomGraph &ag, const PhaseSpace &ps, const int mol_index) {
  const ChemicalDetailsKit cdk = ag.getChemicalDetailsKit();
  moleculeValidityCheck(mol_index, cdk.nmol, "centerOfMass");
  const PhaseSpaceReader psr = ps.data();
  return centerOfMass(psr.xcrd, psr.ycrd, psr.zcrd, cdk.masses, cdk.mol_limits[mol_index],
                      cdk.mol_limits[mol_index + 1], cdk.mol_contents);
}

//-------------------------------------------------------------------------------------------------
double3 centerOfMass(const AtomGraph *ag, const PhaseSpace *ps, const int mol_index) {
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  moleculeValidityCheck(mol_index, cdk.nmol, "centerOfMass");
  const PhaseSpaceReader psr = ps->data();
  return centerOfMass(psr.xcrd, psr.ycrd, psr.zcrd, cdk.masses, cdk.mol_limits[mol_index],
                      cdk.mol_limits[mol_index + 1], cdk.mol_contents);
}

//-------------------------------------------------------------------------------------------------
double3 centerOfMass(const AtomGraph &ag, const CoordinateFrame &cf, const int mol_index) {
  const ChemicalDetailsKit cdk = ag.getChemicalDetailsKit();
  moleculeValidityCheck(mol_index, cdk.nmol, "centerOfMass");
  const CoordinateFrameReader cfr = cf.data();
  return centerOfMass(cfr.xcrd, cfr.ycrd, cfr.zcrd, cdk.masses, cdk.mol_limits[mol_index],
                      cdk.mol_limits[mol_index + 1],  cdk.mol_contents);
}

//-------------------------------------------------------------------------------------------------
double3 centerOfMass(const AtomGraph *ag, const CoordinateFrame *cf, const int mol_index) {
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  moleculeValidityCheck(mol_index, cdk.nmol, "centerOfMass");
  const CoordinateFrameReader cfr = cf->data();
  return centerOfMass(cfr.xcrd, cfr.ycrd, cfr.zcrd, cdk.masses, cdk.mol_limits[mol_index],
                      cdk.mol_limits[mol_index + 1], cdk.mol_contents);
}

//-------------------------------------------------------------------------------------------------
double3 centerOfMass(const double* xcrd, const double* ycrd, const double* zcrd,
                     const double* masses, const int mol_start, const int mol_end,
                     const int* mol_contents) {
  double3 result = { 0.0, 0.0, 0.0 };
  double total_mass = 0.0;
  for (int i = mol_start; i < mol_end; i++) {
    const int atom_idx = (mol_contents != nullptr) ? mol_contents[i] : i;
    const double atom_mass = masses[atom_idx];
    result.x += atom_mass * xcrd[i];
    result.y += atom_mass * ycrd[i];
    result.z += atom_mass * zcrd[i];
    total_mass += atom_mass;
  }
  total_mass = 1.0 / total_mass;
  result.x *= total_mass;
  result.y *= total_mass;
  result.z *= total_mass;
  return result;
}

//-------------------------------------------------------------------------------------------------
double3 molecularTorque(const AtomGraph &ag, const PhaseSpace &ps, const int mol_index) {
  const ChemicalDetailsKit cdk = ag.getChemicalDetailsKit();
  moleculeValidityCheck(mol_index, cdk.nmol, "molecularTorque");
  const PhaseSpaceReader psr = ps.data();
  return molecularTorque(psr.xcrd, psr.ycrd, psr.zcrd, psr.xfrc, psr.yfrc, psr.zfrc, cdk.masses,
                         cdk.mol_contents, cdk.mol_limits[mol_index],
                         cdk.mol_limits[mol_index + 1]);
}

//-------------------------------------------------------------------------------------------------
double3 molecularTorque(const AtomGraph *ag, const PhaseSpace *ps, const int mol_index) {
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  moleculeValidityCheck(mol_index, cdk.nmol, "molecularTorque");
  const PhaseSpaceReader psr = ps->data();
  return molecularTorque(psr.xcrd, psr.ycrd, psr.zcrd, psr.xfrc, psr.yfrc, psr.zfrc, cdk.masses,
                         cdk.mol_contents, cdk.mol_limits[mol_index],
                         cdk.mol_limits[mol_index + 1]);
}

//-------------------------------------------------------------------------------------------------
double3 molecularTorque(const double* xcrd, const double* ycrd, const double* zcrd,
                        const double* xfrc, const double* yfrc, const double* zfrc,
                        const double* masses, const int* mol_contents, const int mol_start,
                        const int mol_end) {
  double3 result = { 0.0, 0.0, 0.0 };
  const double3 c_of_m = centerOfMass(xcrd, ycrd, zcrd, masses, mol_start, mol_end, mol_contents);
  for (int i = mol_start; i < mol_end; i++) {

    // Compute the displacement from the center of mass.  This relies on cross-products, but
    // do the accumulation directly for the most efficient code.
    const double dx = xcrd[i] - c_of_m.x;
    const double dy = ycrd[i] - c_of_m.y;
    const double dz = zcrd[i] - c_of_m.z;
    const double fx = xfrc[i];
    const double fy = yfrc[i];
    const double fz = zfrc[i];
    result.x += (dy * fz) - (dz * fy);
    result.y -= (dx * fz) - (dz * fx);
    result.z += (dx * fy) - (dy * fx);
  }
  return result;
}

} // namespace structure
} // namespace stormm
