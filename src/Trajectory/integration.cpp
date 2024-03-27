#include "copyright.h"
#include "Constants/symbol_values.h"
#include "integration.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
void velocityVerletVelocityUpdate(PhaseSpaceWriter *psw, const ChemicalDetailsKit &cdk,
                                  const ThermostatReader<double> &tstw) {
  velocityVerletVelocityUpdate<double, double>(psw->xvel, psw->yvel, psw->zvel, psw->xfrc,
                                               psw->yfrc, psw->zfrc, psw->natom, cdk.masses,
                                               psw->vxalt, psw->vyalt, psw->vzalt, tstw);
}

//-------------------------------------------------------------------------------------------------
void velocityVerletVelocityUpdate(PhaseSpace *ps, const AtomGraph *ag, const Thermostat *tst) {
  PhaseSpaceWriter psw = ps->data();
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  const ThermostatReader<double> tstw = tst->dpData();
  velocityVerletVelocityUpdate(&psw, cdk, tstw);
}
  
//-------------------------------------------------------------------------------------------------
void velocityVerletVelocityUpdate(PhaseSpace *ps, const AtomGraph &ag, const Thermostat &tst) {
  PhaseSpaceWriter psw = ps->data();
  const ChemicalDetailsKit cdk = ag.getChemicalDetailsKit();
  const ThermostatReader<double> tstw = tst.dpData();
  velocityVerletVelocityUpdate(&psw, cdk, tstw);
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

} // namespace structure
} // namespace stormm
