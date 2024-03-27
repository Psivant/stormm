#include "copyright.h"
#include "Numerics/split_fixed_precision.h"
#include "Synthesis/synthesis_enumerators.h"
#include "global_manipulation.h"
#include "structure_utils.h"

namespace stormm {
namespace structure {

using numerics::hostDoubleToInt95;
using numerics::hostInt95ToDouble;
using numerics::hostSplitFPSum;
using synthesis::StructureSource;
using topology::UnitCellType;
  
//-------------------------------------------------------------------------------------------------
void rotateCoordinates(CoordinateFrameWriter *cfw, const VirtualSiteKit<double> &vsk,
                       const double alpha, const double beta, const double gamma,
                       const int lower_limit, const int upper_limit) {
  coordinateBoundsCheck(lower_limit, upper_limit, cfw->natom, "rotateCoordinates");
  const int actual_upper_limit = (upper_limit > lower_limit) ? upper_limit : cfw->natom;
  rotateCoordinates<double, double>(cfw->xcrd, cfw->ycrd, cfw->zcrd, alpha, beta, gamma,
                                    lower_limit, actual_upper_limit);
  placeVirtualSites<double, double>(cfw->xcrd, cfw->ycrd, cfw->zcrd, nullptr, nullptr,
                                    UnitCellType::NONE, vsk);
}

//-------------------------------------------------------------------------------------------------
void rotateCoordinates(CoordinateFrame *cf, const AtomGraph &ag, const double alpha,
                       const double beta, const double gamma, const int lower_limit,
                       const int upper_limit) {
  CoordinateFrameWriter cfw = cf->data();
  coordinateBoundsCheck(lower_limit, upper_limit, cfw.natom, "rotateCoordinates");
  const int actual_upper_limit = (upper_limit > lower_limit) ? upper_limit : cfw.natom;
  const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
  rotateCoordinates<double, double>(cfw.xcrd, cfw.ycrd, cfw.zcrd, alpha, beta, gamma, lower_limit,
                                    actual_upper_limit, &vsk);
}

//-------------------------------------------------------------------------------------------------
void rotateCoordinates(PhaseSpaceWriter *psw, const VirtualSiteKit<double> &vsk,
                       const double alpha, const double beta, const double gamma,
                       const int lower_limit, const int upper_limit) {
  coordinateBoundsCheck(lower_limit, upper_limit, psw->natom, "rotateCoordinates");
  const int actual_upper_limit = (upper_limit > lower_limit) ? upper_limit : psw->natom;
  rotateCoordinates<double, double>(psw->xcrd, psw->ycrd, psw->zcrd, alpha, beta, gamma,
                                    lower_limit, actual_upper_limit);
  placeVirtualSites<double, double>(psw->xcrd, psw->ycrd, psw->zcrd, nullptr, nullptr,
                                    UnitCellType::NONE, vsk);
}

//-------------------------------------------------------------------------------------------------
void rotateCoordinates(PhaseSpace *ps, const AtomGraph &ag, const double alpha, const double beta,
                       const double gamma, const int lower_limit, const int upper_limit) {
  PhaseSpaceWriter psw = ps->data();
  coordinateBoundsCheck(lower_limit, upper_limit, psw.natom, "rotateCoordinates");
  const int actual_upper_limit = (upper_limit > lower_limit) ? upper_limit : psw.natom;
  const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
  rotateCoordinates<double, double>(psw.xcrd, psw.ycrd, psw.zcrd, alpha, beta, gamma, lower_limit,
                                    actual_upper_limit, &vsk);
}

//-------------------------------------------------------------------------------------------------
void rotateCoordinates(Condensate *cdns, const int system_index, const double alpha,
                       const double beta, const double gamma, const int lower_limit,
                       const int upper_limit) {
  if (cdns->getBasis() != StructureSource::SYNTHESIS) {
    rtErr("When calling the overload with a Condensate, provide the topology if the coordinate "
          "object is based on a series, otherwise let the program trace back to the original "
          "topology using the Condensate's foundational synthesis.", "rotateCoordinates");
  }
  
  // Use precision commensurate with the object's encoding.  When calling this function with a
  // Condensate and without a topology, the Condensate is assumed to be based on a synthesis which
  // will contain the appropriate topology pointers.
  CondensateWriter cdnsw = cdns->data();
  const AtomGraph *ag = cdns->getSynthesisPointer()->getSystemTopologyPointer(system_index);
  switch (cdns->getMode()) {
  case PrecisionModel::DOUBLE:
    {
      const VirtualSiteKit<double> vsk = ag->getDoublePrecisionVirtualSiteKit();
      const int actual_upper_limit = (upper_limit > lower_limit) ? upper_limit :
                                                                   ag->getAtomCount();
      rotateCoordinates<double>(&cdnsw, system_index, vsk, alpha, beta, gamma, lower_limit,
                                actual_upper_limit);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const VirtualSiteKit<float> vsk = ag->getSinglePrecisionVirtualSiteKit();
      const int actual_upper_limit = (upper_limit > lower_limit) ? upper_limit :
                                                                   ag->getAtomCount();
      rotateCoordinates<float>(&cdnsw, system_index, vsk, alpha, beta, gamma, lower_limit,
                               actual_upper_limit);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void rotateCoordinates(Condensate *cdns, const int system_index, const AtomGraph &ag,
                       const double alpha, const double beta, const double gamma,
                       const int lower_limit, const int upper_limit) {
  if (cdns->getBasis() != StructureSource::SERIES) {
    rtErr("When calling the overload with a Condensate, provide the topology if the coordinate "
          "object is based on a series, otherwise let the program trace back to the original "
          "topology using the Condensate's foundational synthesis.", "rotateCoordinates");
  }

  // Use precision commensurate with the object's encoding.
  CondensateWriter cdnsw = cdns->data();
  switch (cdns->getMode()) {
  case PrecisionModel::DOUBLE:
    {
      const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
      rotateCoordinates<double>(&cdnsw, system_index, vsk, alpha, beta, gamma, lower_limit,
                                upper_limit);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const VirtualSiteKit<float> vsk = ag.getSinglePrecisionVirtualSiteKit();
      rotateCoordinates<float>(&cdnsw, system_index, vsk, alpha, beta, gamma, lower_limit,
                               upper_limit);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void translateCoordinates(CoordinateFrameWriter *cfw, const VirtualSiteKit<double> &vsk,
                          const double xmove, const double ymove, const double zmove,
                          const int lower_limit, const int upper_limit) {
  coordinateBoundsCheck(lower_limit, upper_limit, cfw->natom, "translateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? cfw->natom : upper_limit;
  translateCoordinates<double, double>(cfw->xcrd, cfw->ycrd, cfw->zcrd, xmove, ymove, zmove,
                                       lower_limit, actual_upper_limit);
  placeVirtualSites<double, double>(cfw->xcrd, cfw->ycrd, cfw->zcrd, nullptr, nullptr,
                                    UnitCellType::NONE, vsk);
}

//-------------------------------------------------------------------------------------------------
void translateCoordinates(CoordinateFrame *cf, const AtomGraph &ag, const double xmove,
                          const double ymove, const double zmove, const int lower_limit,
                          const int upper_limit) {
  CoordinateFrameWriter cfw = cf->data();
  coordinateBoundsCheck(lower_limit, upper_limit, cfw.natom, "translateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? cfw.natom : upper_limit;
  const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
  translateCoordinates<double, double>(cfw.xcrd, cfw.ycrd, cfw.zcrd, xmove, ymove, zmove,
                                       lower_limit, actual_upper_limit, &vsk);
}

//-------------------------------------------------------------------------------------------------
void translateCoordinates(PhaseSpaceWriter *psw, const VirtualSiteKit<double> &vsk,
                          const double xmove, const double ymove, const double zmove,
                          const int lower_limit, const int upper_limit) {
  coordinateBoundsCheck(lower_limit, upper_limit, psw->natom, "translateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? psw->natom : upper_limit;
  translateCoordinates<double, double>(psw->xcrd, psw->ycrd, psw->zcrd, xmove, ymove, zmove,
                                       lower_limit, actual_upper_limit);
  placeVirtualSites<double, double>(psw->xcrd, psw->ycrd, psw->zcrd, nullptr, nullptr,
                                    UnitCellType::NONE, vsk);
}

//-------------------------------------------------------------------------------------------------
void translateCoordinates(PhaseSpace *ps, const AtomGraph &ag, const double xmove,
                          const double ymove, const double zmove, const int lower_limit,
                          const int upper_limit) {
  PhaseSpaceWriter psw = ps->data();
  coordinateBoundsCheck(lower_limit, upper_limit, psw.natom, "translateCoordinates");
  const int actual_upper_limit = (upper_limit <= lower_limit) ? psw.natom : upper_limit;
  const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
  translateCoordinates<double, double>(psw.xcrd, psw.ycrd, psw.zcrd, xmove, ymove, zmove,
                                       lower_limit, actual_upper_limit, &vsk);
}

//-------------------------------------------------------------------------------------------------
void translateCoordinates(Condensate *cdns, const int system_index, const double xmove,
                          const double ymove, const double zmove, const int lower_limit,
                          const int upper_limit) {
  if (cdns->getBasis() != StructureSource::SYNTHESIS) {
    rtErr("When calling the overload with a Condensate, provide the topology if the coordinate "
          "object is based on a series, otherwise let the program trace back to the original "
          "topology using the Condensate's foundational synthesis.", "translateCoordinates");
  }
  CondensateWriter cdnsw = cdns->data();
  const AtomGraph *ag = cdns->getSynthesisPointer()->getSystemTopologyPointer(system_index);
  switch (cdns->getMode()) {
  case PrecisionModel::DOUBLE:
    {
      const VirtualSiteKit<double> vsk = ag->getDoublePrecisionVirtualSiteKit();
      translateCoordinates<double>(&cdnsw, system_index, vsk, xmove, ymove, zmove, lower_limit,
                                   upper_limit);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const VirtualSiteKit<float> vsk = ag->getSinglePrecisionVirtualSiteKit();
      translateCoordinates<float>(&cdnsw, system_index, vsk, xmove, ymove, zmove, lower_limit,
                                  upper_limit);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void translateCoordinates(Condensate *cdns, const int system_index, const AtomGraph &ag,
                          const double xmove, const double ymove, const double zmove,
                          const int lower_limit, const int upper_limit) {
  if (cdns->getBasis() != StructureSource::SERIES) {
    rtErr("When calling the overload with a Condensate, provide the topology if the coordinate "
          "object is based on a series, otherwise let the program trace back to the original "
          "topology using the Condensate's foundational synthesis.", "translateCoordinates");
  }
  CondensateWriter cdnsw = cdns->data();
  switch (cdns->getMode()) {
  case PrecisionModel::DOUBLE:
    {
      const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
      translateCoordinates<double>(&cdnsw, system_index, vsk, xmove, ymove, zmove, lower_limit,
                                   upper_limit);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const VirtualSiteKit<float> vsk = ag.getSinglePrecisionVirtualSiteKit();
      translateCoordinates<float>(&cdnsw, system_index, vsk, xmove, ymove, zmove, lower_limit,
                                  upper_limit);
    }
    break;
  }  
}

} // namespace structure
} // namespace stormm
