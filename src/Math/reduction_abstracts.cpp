#include "copyright.h"
#include "Constants/behavior.h"
#include "reduction_abstracts.h"

namespace stormm {
namespace stmath {

using constants::CartesianDimension;
  
//-------------------------------------------------------------------------------------------------
ReductionKit::ReductionKit(const int nrdwu_in, const RdwuPerSystem rps_in,
                           const int* rdwu_abstracts_in, const int* atom_counts_in) :
  nrdwu{nrdwu_in}, rps{rps_in}, rdwu_abstracts{rdwu_abstracts_in}, atom_counts{atom_counts_in}
{}

//-------------------------------------------------------------------------------------------------
ReductionKit::ReductionKit(const AtomGraphSynthesis &poly_ag, const HybridTargetLevel tier) :
  nrdwu{poly_ag.getReductionWorkUnitCount()},
  rps{poly_ag.getRdwuPerSystem()},
  rdwu_abstracts{poly_ag.getReductionWorkUnitAbstracts().data(tier)},
  atom_counts{poly_ag.getSystemAtomCounts().data(tier)}
{}

//-------------------------------------------------------------------------------------------------
ConjGradSubstrate::ConjGradSubstrate(PsSynthesisWriter poly_psw, ReductionBridge *rbg,
                                     const HybridTargetLevel tier) :
    inv_frc_scale{poly_psw.inv_frc_scale},
    xfrc{poly_psw.xfrc}, yfrc{poly_psw.yfrc}, zfrc{poly_psw.zfrc},
    xfrc_ovrf{poly_psw.xfrc_ovrf}, yfrc_ovrf{poly_psw.yfrc_ovrf}, zfrc_ovrf{poly_psw.zfrc_ovrf},
    xprv{poly_psw.xalt}, yprv{poly_psw.yalt}, zprv{poly_psw.zalt},
    xprv_ovrf{poly_psw.xalt_ovrf}, yprv_ovrf{poly_psw.yalt_ovrf}, zprv_ovrf{poly_psw.zalt_ovrf},
    x_cg_temp{poly_psw.xvel}, y_cg_temp{poly_psw.yvel}, z_cg_temp{poly_psw.zvel},
    x_cg_temp_ovrf{poly_psw.xvel_ovrf},
    y_cg_temp_ovrf{poly_psw.yvel_ovrf},
    z_cg_temp_ovrf{poly_psw.zvel_ovrf},
    gg_buffer{rbg->getPointer(CartesianDimension::X, tier)},
    dgg_buffer{rbg->getPointer(CartesianDimension::Y, tier)},
    msum_buffer{rbg->getPointer(CartesianDimension::Z, tier)}
{}

//-------------------------------------------------------------------------------------------------
ConjGradSubstrate::ConjGradSubstrate(PhaseSpaceSynthesis *poly_ps, ReductionBridge *rbg,
                                     const HybridTargetLevel tier) :
    ConjGradSubstrate(poly_ps->data(tier), rbg, tier)
{}

} // namespace stmath
} // namespace stormm
