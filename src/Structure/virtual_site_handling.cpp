#include <cmath>
#include "copyright.h"
#include "Math/vector_ops.h"
#include "local_arrangement.h"
#include "virtual_site_handling.h"

namespace stormm {
namespace structure {

using stmath::dot;
using stmath::project;
using stmath::crossProduct;
using topology::VirtualSiteKind;

//-------------------------------------------------------------------------------------------------
void placeVirtualSites(PhaseSpace *ps, const AtomGraph *ag, const CoordinateCycle affix,
                       const PrecisionModel prec) {
  PhaseSpaceWriter psw = ps->data(affix);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    placeVirtualSites<double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                      psw.unit_cell, ag->getDoublePrecisionVirtualSiteKit());
    break;
  case PrecisionModel::SINGLE:
    placeVirtualSites<double, float>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                     psw.unit_cell, ag->getSinglePrecisionVirtualSiteKit());
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void placeVirtualSites(PhaseSpace *ps, const AtomGraph &ag, const CoordinateCycle affix,
                       const PrecisionModel prec) {
  placeVirtualSites(ps, ag.getSelfPointer(), affix, prec);
}

//-------------------------------------------------------------------------------------------------
void placeVirtualSites(PhaseSpace *ps, const AtomGraph *ag, const PrecisionModel prec) {
  placeVirtualSites(ps, ag, ps->getCyclePosition(), prec);
}

//-------------------------------------------------------------------------------------------------
void placeVirtualSites(PhaseSpace *ps, const AtomGraph &ag, const PrecisionModel prec) {
  placeVirtualSites(ps, ag.getSelfPointer(), ps->getCyclePosition(), prec);
}

//-------------------------------------------------------------------------------------------------
void placeVirtualSites(CoordinateFrame *cf, const AtomGraph *ag, const PrecisionModel prec) {
  CoordinateFrameWriter cfw = cf->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    placeVirtualSites<double, double>(cfw.xcrd, cfw.ycrd, cfw.zcrd, cfw.umat, cfw.invu,
                                      cfw.unit_cell, ag->getDoublePrecisionVirtualSiteKit());
    break;
  case PrecisionModel::SINGLE:
    placeVirtualSites<double, float>(cfw.xcrd, cfw.ycrd, cfw.zcrd, cfw.umat, cfw.invu,
                                     cfw.unit_cell, ag->getSinglePrecisionVirtualSiteKit());
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void placeVirtualSites(CoordinateFrame *cf, const AtomGraph &ag, const PrecisionModel prec) {
  placeVirtualSites(cf, ag.getSelfPointer(), prec);
}

//-------------------------------------------------------------------------------------------------
void placeVirtualSites(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                       const CoordinateCycle affix, const PrecisionModel prec) {
  PsSynthesisWriter poly_psw = poly_ps->data(affix);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyValenceKit<double> poly_vk = poly_ag->getDoublePrecisionValenceKit();
      const SyAtomUpdateKit<double,
                            double2,
                            double4_16a> poly_auk = poly_ag->getDoublePrecisionAtomUpdateKit();
      placeVirtualSites<double, double2, double4_16a>(&poly_psw, poly_vk, poly_auk);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag->getSinglePrecisionValenceKit();
      const SyAtomUpdateKit<float,
                            float2, float4> poly_auk = poly_ag->getSinglePrecisionAtomUpdateKit();
      placeVirtualSites<float, float2, float4>(&poly_psw, poly_vk, poly_auk);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void placeVirtualSites(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                       const CoordinateCycle affix, const PrecisionModel prec) {
  placeVirtualSites(poly_ps, poly_ag.getSelfPointer(), affix, prec);
}

//-------------------------------------------------------------------------------------------------
void placeVirtualSites(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                       const PrecisionModel prec) {
  placeVirtualSites(poly_ps, poly_ag, poly_ps->getCyclePosition(), prec);
}

//-------------------------------------------------------------------------------------------------
void placeVirtualSites(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                       const PrecisionModel prec) {
  placeVirtualSites(poly_ps, poly_ag.getSelfPointer(), poly_ps->getCyclePosition(), prec);
}

//-------------------------------------------------------------------------------------------------
void transmitVirtualSiteForces(PhaseSpace *ps, const AtomGraph *ag, const PrecisionModel prec) {
  PhaseSpaceWriter psw = ps->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    transmitVirtualSiteForces<double,
                              double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.xfrc, psw.yfrc,
                                              psw.zfrc, psw.umat, psw.invu, psw.unit_cell,
                                              ag->getDoublePrecisionVirtualSiteKit());
    break;
  case PrecisionModel::SINGLE:
    transmitVirtualSiteForces<double,
                              double, float>(psw.xcrd, psw.ycrd, psw.zcrd, psw.xfrc, psw.yfrc,
                                             psw.zfrc, psw.umat, psw.invu, psw.unit_cell,
                                             ag->getSinglePrecisionVirtualSiteKit());
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void transmitVirtualSiteForces(PhaseSpace *ps, const AtomGraph &ag, const PrecisionModel prec) {
  transmitVirtualSiteForces(ps, ag.getSelfPointer(), prec);
}

//-------------------------------------------------------------------------------------------------
void transmitVirtualSiteForces(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                               const PrecisionModel prec) {
  PsSynthesisWriter poly_psw = poly_ps->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const SyValenceKit<double> poly_vk = poly_ag->getDoublePrecisionValenceKit();
      const SyAtomUpdateKit<double,
                            double2,
                            double4_16a> poly_auk = poly_ag->getDoublePrecisionAtomUpdateKit();
      transmitVirtualSiteForces(&poly_psw, poly_vk, poly_auk);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const SyValenceKit<float> poly_vk = poly_ag->getSinglePrecisionValenceKit();
      const SyAtomUpdateKit<float,
                            float2, float4> poly_auk = poly_ag->getSinglePrecisionAtomUpdateKit();
      transmitVirtualSiteForces(&poly_psw, poly_vk, poly_auk);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void transmitVirtualSiteForces(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &ag,
                               const PrecisionModel prec) {

}

} // namespace structure
} // namespace stormm
