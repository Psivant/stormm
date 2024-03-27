#include "copyright.h"
#include "minimization.h"

namespace stormm {
namespace mm {

using energy::EvaluateForce;
using stmath::invertSquareMatrix;
using stmath::matrixVectorMultiply;
using structure::placeVirtualSites;
using structure::transmitVirtualSiteForces;
  
//-------------------------------------------------------------------------------------------------
ScoreCard minimize(PhaseSpace *ps, const AtomGraph &ag, const RestraintApparatus &ra,
                   const StaticExclusionMask &se, const MinimizeControls &mincon,
                   const int nrg_scale_bits) {
  const NeckGeneralizedBornTable ngb_tab;
  return minimize(ps, ag, ngb_tab, ra, se, mincon, nrg_scale_bits);
}

//-------------------------------------------------------------------------------------------------
ScoreCard minimize(PhaseSpace *ps, const AtomGraph &ag, const NeckGeneralizedBornTable &ngb_tab,
                   const RestraintApparatus &ra, const StaticExclusionMask &se,
                   const MinimizeControls &mincon, const int nrg_scale_bits) {
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  const ImplicitSolventKit<double> isk = ag.getDoublePrecisionImplicitSolventKit();
  const NeckGeneralizedBornKit<double> ngbk = ngb_tab.dpData();
  const ValenceKit<double> vk = ag.getDoublePrecisionValenceKit();
  const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
  PhaseSpaceWriter psw = ps->data();
  return minimize<double, double,
                  double, double2, double4>(psw.xcrd, psw.ycrd, psw.zcrd, psw.xfrc, psw.yfrc,
                                            psw.zfrc, psw.xvel, psw.yvel, psw.zvel, psw.xalt,
                                            psw.yalt, psw.zalt, vk, nbk, isk, ngbk, ra.dpData(),
                                            vsk, se.data(), mincon, nrg_scale_bits, 1.0, 1.0);
}

//-------------------------------------------------------------------------------------------------
ScoreCard minimize(PhaseSpace *ps, const AtomGraph *ag, const StaticExclusionMask &se,
                   const MinimizeControls &mincon, const int nrg_scale_bits) {
  const NeckGeneralizedBornTable ngb_tab;
  return minimize(ps, ag, ngb_tab, se, mincon, nrg_scale_bits);
}

//-------------------------------------------------------------------------------------------------
ScoreCard minimize(PhaseSpace *ps, const AtomGraph *ag, const NeckGeneralizedBornTable &ngb_tab,
                   const StaticExclusionMask &se, const MinimizeControls &mincon,
                   const int nrg_scale_bits) {
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  const ImplicitSolventKit<double> isk = ag->getDoublePrecisionImplicitSolventKit();
  const NeckGeneralizedBornKit<double> ngbk = ngb_tab.dpData();
  const ValenceKit<double> vk = ag->getDoublePrecisionValenceKit();
  const VirtualSiteKit<double> vsk = ag->getDoublePrecisionVirtualSiteKit();
  const RestraintApparatus ra(ag);
  PhaseSpaceWriter psw = ps->data();
  return minimize<double, double,
                  double, double2, double4>(psw.xcrd, psw.ycrd, psw.zcrd, psw.xfrc, psw.yfrc,
                                            psw.zfrc, psw.xvel, psw.yvel, psw.zvel, psw.xalt,
                                            psw.yalt, psw.zalt, vk, nbk, isk, ngbk, ra.dpData(),
                                            vsk, se.data(), mincon, nrg_scale_bits, 1.0, 1.0);
}

//-------------------------------------------------------------------------------------------------
ScoreCard minimize(PhaseSpaceWriter psw, const ValenceKit<double> &vk,
                   const NonbondedKit<double> &nbk, const ImplicitSolventKit<double> &isk,
                   const NeckGeneralizedBornKit<double> &ngbk,
                   const RestraintKit<double, double2, double4> &rar,
                   const VirtualSiteKit<double> &vsk, const StaticExclusionMaskReader &ser,
                   const MinimizeControls &mincon, int nrg_scale_bits) {
  return minimize<double, double,
                  double, double2, double4>(psw.xcrd, psw.ycrd, psw.zcrd, psw.xfrc, psw.yfrc,
                                            psw.zfrc, psw.xvel, psw.yvel, psw.zvel, psw.xalt,
                                            psw.yalt, psw.zalt, vk, nbk, isk, ngbk, rar, vsk, ser,
                                            mincon, nrg_scale_bits, 1.0, 1.0);
}
  
} // namespace mm
} // namespace stormm
