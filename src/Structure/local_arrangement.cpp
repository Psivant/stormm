#include <cmath>
#include "copyright.h"
#include "Math/vector_ops.h"
#include "Reporting/error_format.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "local_arrangement.h"

namespace stormm {
namespace structure {

using stmath::crossProduct;
using trajectory::CoordinateFrameWriter;

//-------------------------------------------------------------------------------------------------
void imageCoordinates(CoordinateFrameWriter cfw, const ImagingMethod style) {
  imageCoordinates<double, double>(cfw.xcrd, cfw.ycrd, cfw.zcrd, cfw.natom, cfw.umat, cfw.invu,
                                   cfw.unit_cell, style);
}

//-------------------------------------------------------------------------------------------------
void imageCoordinates(CoordinateFrame *cf, const ImagingMethod style) {
  imageCoordinates(cf->data(), style);
}

//-------------------------------------------------------------------------------------------------
void imageCoordinates(PhaseSpaceWriter psw, const ImagingMethod style) {
  imageCoordinates<double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.natom, psw.umat, psw.invu,
                                   psw.unit_cell, style);
}

//-------------------------------------------------------------------------------------------------
void imageCoordinates(PhaseSpace *ps, const ImagingMethod style) {
  imageCoordinates(ps->data(), style);
}

//-------------------------------------------------------------------------------------------------
double distance(const int atom_i, const int atom_j, const CoordinateFrameReader &cfr) {
  return distance<double, double>(atom_i, atom_j, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat, cfr.invu,
                                  cfr.unit_cell);
}

//-------------------------------------------------------------------------------------------------
double distance(const int atom_i, const int atom_j, const CoordinateFrame *cf) {
  return distance(atom_i, atom_j, cf->data());
}

//-------------------------------------------------------------------------------------------------
double distance(const int atom_i, const int atom_j, const CoordinateFrame &cf) {
  return distance(atom_i, atom_j, cf.data());
}

//-------------------------------------------------------------------------------------------------
double distance(const int atom_i, const int atom_j, const PhaseSpaceReader &psr) {
  return distance<double, double>(atom_i, atom_j, psr.xcrd, psr.ycrd, psr.zcrd, psr.umat, psr.invu,
                                  psr.unit_cell);
}

//-------------------------------------------------------------------------------------------------
double distance(const int atom_i, const int atom_j, const PhaseSpace *ps) {

  // When possible, convert from PhaseSpace objects directly to the CoordinateFrameReader.  This
  // is more efficient than making a PhaseSpace abstract.
  return distance(atom_i, atom_j, CoordinateFrameReader(ps));
}

//-------------------------------------------------------------------------------------------------
double distance(const int atom_i, const int atom_j, const PhaseSpace &ps) {
  return distance(atom_i, atom_j, CoordinateFrameReader(ps));
}

//-------------------------------------------------------------------------------------------------
double angle(const int atom_i, const int atom_j, const int atom_k,
             const CoordinateFrameReader &cfr) {
  return angle<double, double>(atom_i, atom_j, atom_k, cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat,
                               cfr.invu, cfr.unit_cell);
}

//-------------------------------------------------------------------------------------------------
double angle(const int atom_i, const int atom_j, const int atom_k, const CoordinateFrame *cf) {
  return angle(atom_i, atom_j, atom_k, cf->data());
}

//-------------------------------------------------------------------------------------------------
double angle(const int atom_i, const int atom_j, const int atom_k, const CoordinateFrame &cf) {
  return angle(atom_i, atom_j, atom_k, cf.data());
}

//-------------------------------------------------------------------------------------------------
double angle(const int atom_i, const int atom_j, const int atom_k,
             const PhaseSpaceReader &psr) {
  return angle<double, double>(atom_i, atom_j, atom_k, psr.xcrd, psr.ycrd, psr.zcrd, psr.umat,
                               psr.invu, psr.unit_cell);
}

//-------------------------------------------------------------------------------------------------
double angle(const int atom_i, const int atom_j, const int atom_k, const PhaseSpace *ps) {
  return angle(atom_i, atom_j, atom_k, CoordinateFrameReader(ps));
}

//-------------------------------------------------------------------------------------------------
double angle(const int atom_i, const int atom_j, const int atom_k, const PhaseSpace &ps) {
  return angle(atom_i, atom_j, atom_k, CoordinateFrameReader(ps));
}

//-------------------------------------------------------------------------------------------------
double dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l,
                     const CoordinateFrameReader &cfr) {
  return dihedralAngle<double, double>(atom_i, atom_j, atom_k, atom_l, cfr.xcrd, cfr.ycrd,
                                       cfr.zcrd, cfr.umat, cfr.invu, cfr.unit_cell);
}

//-------------------------------------------------------------------------------------------------
double dihedralAngle(const int atom_i, const int atom_j, const int atom_k, const int atom_l,
                     const CoordinateFrame *cf) {
  return dihedralAngle(atom_i, atom_j, atom_k, atom_l, cf->data());
}

//-------------------------------------------------------------------------------------------------
double dihedralAngle(const int atom_i, const int atom_j, const int atom_k, const int atom_l,
                     const CoordinateFrame &cf) {
  return dihedralAngle(atom_i, atom_j, atom_k, atom_l, cf.data());
}

//-------------------------------------------------------------------------------------------------
double dihedralAngle(const int atom_i, const int atom_j, const int atom_k, const int atom_l,
                     const PhaseSpaceReader &psr) {
  return dihedralAngle<double, double>(atom_i, atom_j, atom_k, atom_l, psr.xcrd, psr.ycrd,
                                       psr.zcrd, psr.umat, psr.invu, psr.unit_cell);
}
  
//-------------------------------------------------------------------------------------------------
double dihedralAngle(const int atom_i, const int atom_j, const int atom_k, const int atom_l,
                     const PhaseSpace *ps) {
  return dihedralAngle(atom_i, atom_j, atom_k, atom_l, CoordinateFrameReader(ps));
}

//-------------------------------------------------------------------------------------------------
double dihedralAngle(const int atom_i, const int atom_j, const int atom_k, const int atom_l,
                     const PhaseSpace &ps) {
  return dihedralAngle(atom_i, atom_j, atom_k, atom_l, CoordinateFrameReader(ps));
}

} // namespace structure
} // namespace stormm
