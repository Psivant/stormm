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
void placeVirtualSites(PhaseSpace *ps, const AtomGraph &ag) {
  PhaseSpaceWriter psw = ps->data();
  placeVirtualSites<double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                    psw.unit_cell, ag.getDoublePrecisionVirtualSiteKit());
}

//-------------------------------------------------------------------------------------------------
void placeVirtualSites(PhaseSpace *ps, const AtomGraph *ag) {
  PhaseSpaceWriter psw = ps->data();
  placeVirtualSites<double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.umat, psw.invu,
                                    psw.unit_cell, ag->getDoublePrecisionVirtualSiteKit());
}

//-------------------------------------------------------------------------------------------------
void placeVirtualSites(CoordinateFrame *cf, const AtomGraph &ag) {
  CoordinateFrameWriter cfw = cf->data();
  placeVirtualSites<double, double>(cfw.xcrd, cfw.ycrd, cfw.zcrd, cfw.umat, cfw.invu,
                                    cfw.unit_cell, ag.getDoublePrecisionVirtualSiteKit());
}

//-------------------------------------------------------------------------------------------------
void placeVirtualSites(CoordinateFrame *cf, const AtomGraph *ag) {
  CoordinateFrameWriter cfw = cf->data();
  placeVirtualSites<double, double>(cfw.xcrd, cfw.ycrd, cfw.zcrd, cfw.umat, cfw.invu,
                                    cfw.unit_cell, ag->getDoublePrecisionVirtualSiteKit());
}

//-------------------------------------------------------------------------------------------------
void transmitVirtualSiteForces(PhaseSpace *ps, const AtomGraph &ag) {
  PhaseSpaceWriter psw = ps->data();
  transmitVirtualSiteForces<double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.xfrc, psw.yfrc,
                                            psw.zfrc, psw.umat, psw.invu, psw.unit_cell,
                                            ag.getDoublePrecisionVirtualSiteKit());
}

//-------------------------------------------------------------------------------------------------
void transmitVirtualSiteForces(PhaseSpace *ps, const AtomGraph *ag) {
  PhaseSpaceWriter psw = ps->data();
  transmitVirtualSiteForces<double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.xfrc, psw.yfrc,
                                            psw.zfrc, psw.umat, psw.invu, psw.unit_cell,
                                            ag->getDoublePrecisionVirtualSiteKit());
}

} // namespace structure
} // namespace stormm
