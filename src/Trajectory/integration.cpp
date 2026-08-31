#include "copyright.h"
#include "Constants/symbol_values.h"
#include "integration.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
void velocityVerletVelocityUpdate(PhaseSpaceWriter *psw, const ChemicalDetailsKit &cdk,
                                  const ThermostatReader<double> &tstr) {
  velocityVerletVelocityUpdate<double, double>(psw->xvel, psw->yvel, psw->zvel, psw->xfrc,
                                               psw->yfrc, psw->zfrc, psw->natom, cdk.masses,
                                               psw->vxalt, psw->vyalt, psw->vzalt, tstr);
}

//-------------------------------------------------------------------------------------------------
void velocityVerletVelocityUpdate(PhaseSpace *ps, const AtomGraph *ag, const Thermostat *tst) {
  PhaseSpaceWriter psw = ps->data();
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  const ThermostatReader<double> tstr = tst->dpData();
  velocityVerletVelocityUpdate(&psw, cdk, tstr);
}
  
//-------------------------------------------------------------------------------------------------
void velocityVerletVelocityUpdate(PhaseSpace *ps, const AtomGraph &ag, const Thermostat &tst) {
  PhaseSpaceWriter psw = ps->data();
  const ChemicalDetailsKit cdk = ag.getChemicalDetailsKit();
  const ThermostatReader<double> tstr = tst.dpData();
  velocityVerletVelocityUpdate(&psw, cdk, tstr);
}

//-------------------------------------------------------------------------------------------------
void velocityVerletVelocityUpdate(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis *poly_ag,
                                  const Thermostat *tst, const PrecisionModel prec) {
  PsSynthesisWriter poly_psw = poly_ps->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const ThermostatReader tstr = tst->dpData();
      const SyAtomUpdateKit<double,
                            double2,
                            double4_16a> poly_auk = poly_ag->getDoublePrecisionAtomUpdateKit();
      velocityVerletVelocityUpdate(&poly_psw, poly_auk, tstr);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const ThermostatReader tstr = tst->spData();
      const SyAtomUpdateKit<float,
                            float2, float4> poly_auk = poly_ag->getSinglePrecisionAtomUpdateKit();
      velocityVerletVelocityUpdate(&poly_psw, poly_auk, tstr);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void velocityVerletVelocityUpdate(PhaseSpaceSynthesis *poly_ps, const AtomGraphSynthesis &poly_ag,
                                  const Thermostat &tst, const PrecisionModel prec) {
  velocityVerletVelocityUpdate(poly_ps, poly_ag.getSelfPointer(), tst.getSelfPointer(), prec);
}

//-------------------------------------------------------------------------------------------------
void velocityVerletCoordinateUpdate(PhaseSpaceWriter *psw, const ChemicalDetailsKit &cdk,
                                    const ThermostatReader<double> &tstw) {
  velocityVerletCoordinateUpdate<double, double>(psw->xcrd, psw->ycrd, psw->zcrd, psw->xfrc,
                                                 psw->yfrc, psw->zfrc, psw->natom, cdk.masses,
                                                 psw->xalt, psw->yalt, psw->zalt, psw->vxalt,
                                                 psw->vyalt, psw->vzalt, tstw);
}

//-------------------------------------------------------------------------------------------------
void velocityVerletCoordinateUpdate(PhaseSpace *ps, const AtomGraph *ag, const Thermostat *tst) {
  PhaseSpaceWriter psw = ps->data();
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  velocityVerletCoordinateUpdate<double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.xfrc, psw.yfrc,
                                                 psw.zfrc, psw.natom, cdk.masses, psw.xalt,
                                                 psw.yalt, psw.zalt, psw.vxalt, psw.vyalt,
                                                 psw.vzalt, tst->dpData());
}

//-------------------------------------------------------------------------------------------------
void velocityVerletCoordinateUpdate(PhaseSpace *ps, const AtomGraph &ag, const Thermostat &tst) {
  PhaseSpaceWriter psw = ps->data();
  const ChemicalDetailsKit cdk = ag.getChemicalDetailsKit();
  velocityVerletCoordinateUpdate<double, double>(psw.xcrd, psw.ycrd, psw.zcrd, psw.xfrc, psw.yfrc,
                                                 psw.zfrc, psw.natom, cdk.masses, psw.xalt,
                                                 psw.yalt, psw.zalt, psw.vxalt, psw.vyalt,
                                                 psw.vzalt, tst.dpData());
}

//-------------------------------------------------------------------------------------------------
void velocityVerletCoordinateUpdate(PhaseSpaceSynthesis *poly_ps,
                                    const AtomGraphSynthesis *poly_ag, const Thermostat *tst,
                                    const PrecisionModel prec) {
  PsSynthesisWriter poly_psw = poly_ps->data();
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      const ThermostatReader tstr = tst->dpData();
      const SyAtomUpdateKit<double,
                            double2,
                            double4_16a> poly_auk = poly_ag->getDoublePrecisionAtomUpdateKit();
      velocityVerletCoordinateUpdate(&poly_psw, poly_auk, tstr);
    }
    break;
  case PrecisionModel::SINGLE:
    {
      const ThermostatReader tstr = tst->spData();
      const SyAtomUpdateKit<float,
                            float2, float4> poly_auk = poly_ag->getSinglePrecisionAtomUpdateKit();
      velocityVerletCoordinateUpdate(&poly_psw, poly_auk, tstr);
    }
    break;
  }
}

} // namespace structure
} // namespace stormm
