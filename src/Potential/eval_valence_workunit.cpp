#include "copyright.h"
#include "Math/vector_ops.h"
#include "Topology/atomgraph_enumerators.h"
#include "eval_valence_workunit.h"
#include "valence_potential.h"

namespace stormm {
namespace energy {

using stmath::crossProduct;
using stmath::readBitFromMask;
using topology::TorsionKind;
using trajectory::PhaseSpaceReader;
using trajectory::PhaseSpaceWriter;

//-------------------------------------------------------------------------------------------------
void evalVwuInitEnergy(ScoreCard *ecard, const VwuTask activity, const int sysid) {
  switch (activity) {
  case VwuTask::BOND:
    ecard->initialize(StateVariable::BOND, sysid);
    break;
  case VwuTask::ANGL:
    ecard->initialize(StateVariable::ANGLE, sysid);
    break;
  case VwuTask::DIHE:
    ecard->initialize(StateVariable::PROPER_DIHEDRAL, sysid);
    ecard->initialize(StateVariable::IMPROPER_DIHEDRAL, sysid);
    break;
  case VwuTask::UBRD:
    ecard->initialize(StateVariable::UREY_BRADLEY, sysid);
    break;
  case VwuTask::CIMP:
    ecard->initialize(StateVariable::CHARMM_IMPROPER, sysid);
    break;
  case VwuTask::CMAP:
    ecard->initialize(StateVariable::CMAP, sysid);
    break;
  case VwuTask::INFR14:
    ecard->initialize(StateVariable::ELEC_ONE_FOUR, sysid);
    ecard->initialize(StateVariable::VDW_ONE_FOUR, sysid);
    break;
  case VwuTask::RPOSN:
  case VwuTask::RBOND:
  case VwuTask::RANGL:
  case VwuTask::RDIHE:
    ecard->initialize(StateVariable::RESTRAINT, sysid);
    break;
  case VwuTask::SETTLE:
  case VwuTask::CGROUP:
  case VwuTask::VSITE:
  case VwuTask::CDHE:
  case VwuTask::CBND:
    break;
  case VwuTask::ALL_TASKS:
    ecard->initialize({ StateVariable::BOND, StateVariable::ANGLE, StateVariable::PROPER_DIHEDRAL,
                        StateVariable::IMPROPER_DIHEDRAL, StateVariable::UREY_BRADLEY,
                        StateVariable::CHARMM_IMPROPER }, sysid);
    break;
  }  
}

//-------------------------------------------------------------------------------------------------
void commitVwuEnergies(const llint bond_acc, const llint angl_acc, const llint dihe_acc,
                       const llint impr_acc, const llint ubrd_acc, const llint cimp_acc,
                       const llint cmap_acc, const llint qq14_acc, const llint lj14_acc,
                       const llint rest_acc, const int sysid, const VwuTask activity,
                       ScoreCard *ecard) {
  switch (activity) {
  case VwuTask::BOND:
    ecard->add(StateVariable::BOND, bond_acc, sysid);
    break;
  case VwuTask::ANGL:
    ecard->add(StateVariable::ANGLE, angl_acc, sysid);
    break;
  case VwuTask::DIHE:
    ecard->add(StateVariable::PROPER_DIHEDRAL, dihe_acc, sysid);
    ecard->add(StateVariable::IMPROPER_DIHEDRAL, impr_acc, sysid);
    break;
  case VwuTask::UBRD:
    ecard->add(StateVariable::UREY_BRADLEY, ubrd_acc, sysid);
    break;
  case VwuTask::CIMP:
    ecard->add(StateVariable::CHARMM_IMPROPER, cimp_acc, sysid);
    break;
  case VwuTask::CMAP:
    ecard->add(StateVariable::CMAP, cmap_acc, sysid);
    break;
  case VwuTask::INFR14:
    ecard->add(StateVariable::ELEC_ONE_FOUR, qq14_acc, sysid);
    ecard->add(StateVariable::VDW_ONE_FOUR, lj14_acc, sysid);
    break;
  case VwuTask::RPOSN:
  case VwuTask::RBOND:
  case VwuTask::RANGL:
  case VwuTask::RDIHE:
    ecard->add(StateVariable::RESTRAINT, rest_acc, sysid);
    break;
  case VwuTask::SETTLE:
  case VwuTask::CGROUP:
  case VwuTask::VSITE:
  case VwuTask::CDHE:
  case VwuTask::CBND:
    break;
  case VwuTask::ALL_TASKS:
    ecard->add(StateVariable::BOND, bond_acc, sysid);
    ecard->add(StateVariable::ANGLE, angl_acc, sysid);
    ecard->add(StateVariable::PROPER_DIHEDRAL, dihe_acc, sysid);
    ecard->add(StateVariable::IMPROPER_DIHEDRAL, impr_acc, sysid);
    ecard->add(StateVariable::UREY_BRADLEY, ubrd_acc, sysid);
    ecard->add(StateVariable::CHARMM_IMPROPER, cimp_acc, sysid);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void localVwuEvaluation(const ValenceKit<double> vk, const VirtualSiteKit<double> vsk,
                        const NonbondedKit<double> nbk,
                        const RestraintKit<double, double2, double4> rar,
                        const double* sh_charges, const int* sh_lj_idx, double* sh_xcrd,
                        double* sh_ycrd, double* sh_zcrd, double* sh_xfrc, double* sh_yfrc,
                        double* sh_zfrc, ScoreCard *ecard, const int sysid,
                        const ValenceWorkUnit &my_vwu, const EvaluateForce eval_force,
                        const VwuTask activity, const VwuGoal purpose, const int step_number) {

  // Initialize energy accumulators
  const double nrg_scale_factor = ecard->getEnergyScalingFactor<double>();
  llint bond_acc = 0LL;
  llint angl_acc = 0LL;
  llint dihe_acc = 0LL;
  llint impr_acc = 0LL;
  llint ubrd_acc = 0LL;
  llint cimp_acc = 0LL;
  llint cmap_acc = 0LL;
  llint qq14_acc = 0LL;
  llint lj14_acc = 0LL;
  llint rest_acc = 0LL;

  // Get all tasks from the work unit
  const std::vector<int> task_counts = my_vwu.getTaskCounts();

  // Evaluate bonds and Urey-Bradley terms
  if (activity == VwuTask::BOND || activity == VwuTask::UBRD || activity == VwuTask::ALL_TASKS) {
    const int ncbnd = task_counts[static_cast<int>(VwuTask::CBND)];
    const std::vector<uint> cbnd_acc_mask = my_vwu.getAccumulationFlags(VwuTask::CBND);
    for (int pos = 0; pos < ncbnd; pos++) {
      bool log_term = true;
      switch (purpose) {
      case VwuGoal::ACCUMULATE:
        if (readBitFromMask(cbnd_acc_mask, pos) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(cbnd_acc_mask, pos);
        break;
      }
      const uint2 tinsr = my_vwu.getCompositeBondInstruction(pos);
      const bool is_urey_bradley = ((tinsr.x >> 20) & 0x1);
      
      // Skip Urey-Bradley interactions if only bonds are desired, and skip bonds if only
      // Urey-Bradley interactions are of interest.
      if ((activity == VwuTask::BOND && is_urey_bradley) ||
          (activity == VwuTask::UBRD && (! is_urey_bradley))) {
        continue;
      }
      const int i_atom = (tinsr.x & 0x3ff);
      const int j_atom = ((tinsr.x >> 10) & 0x3ff);
      const int param_idx = tinsr.y;
      const double keq = (is_urey_bradley) ? vk.ubrd_keq[param_idx] : vk.bond_keq[param_idx];
      const double leq = (is_urey_bradley) ? vk.ubrd_leq[param_idx] :
                                             std::abs(vk.bond_leq[param_idx]);
      const double du =
        evalHarmonicStretch<double, double, double>(i_atom, j_atom, keq, leq, sh_xcrd, sh_ycrd,
                                                    sh_zcrd, nullptr, nullptr, UnitCellType::NONE,
                                                    sh_xfrc, sh_yfrc, sh_zfrc, eval_force);
      if (log_term) {
        if (is_urey_bradley) {
          ubrd_acc += llround(du * nrg_scale_factor);
        }
        else {
          bond_acc += llround(du * nrg_scale_factor);
        }
      }
    }
  }
    
  // Evaluate harmonic bond angles
  if (activity == VwuTask::ANGL || activity == VwuTask::ALL_TASKS) {
    const int nangl = task_counts[static_cast<int>(VwuTask::ANGL)];
    const std::vector<uint> angl_acc_mask = my_vwu.getAccumulationFlags(VwuTask::ANGL);
    for (int pos = 0; pos < nangl; pos++) {
      bool log_term = true;
      switch (purpose) {
      case VwuGoal::ACCUMULATE:
        if (readBitFromMask(angl_acc_mask, pos) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(angl_acc_mask, pos);
        break;
      }
      const uint2 tinsr = my_vwu.getAngleInstruction(pos);
      const int i_atom = (tinsr.x & 0x3ff);
      const int j_atom = ((tinsr.x >> 10) & 0x3ff);
      const int k_atom = ((tinsr.x >> 20) & 0x3ff);
      const int param_idx = tinsr.y;
      const double keq = vk.angl_keq[param_idx];
      const double theta0 = vk.angl_theta[param_idx];
      const double du =
        evalHarmonicBend<double, double, double>(i_atom, j_atom, k_atom, keq, theta0, sh_xcrd,
                                                 sh_ycrd, sh_zcrd, nullptr, nullptr,
                                                 UnitCellType::NONE, sh_xfrc, sh_yfrc, sh_zfrc,
                                                 eval_force);
      if (log_term) {
        angl_acc += llround(du * nrg_scale_factor);
      }
    }
  }

  // Evaluate cosine-based dihedrals and CHARMM improper dihedrals
  if (activity == VwuTask::DIHE || activity == VwuTask::CIMP || activity == VwuTask::INFR14 ||
      activity == VwuTask::ALL_TASKS) {
    const int ncdhe = task_counts[static_cast<int>(VwuTask::CDHE)];
    const std::vector<uint> cdhe_acc_mask = my_vwu.getAccumulationFlags(VwuTask::CDHE);
    for (int pos = 0; pos < ncdhe; pos++) {
      bool log_term = true;
      switch (purpose) {
      case VwuGoal::ACCUMULATE:
        if (readBitFromMask(cdhe_acc_mask, pos) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(cdhe_acc_mask, pos);
        break;
      }
      const uint3 tinsr = my_vwu.getCompositeDihedralInstruction(pos);

      // Compute 1:4 interactions and go on if that is all that is required.  Here, the
      // 1:4 interactions are being *inferred* from the dihedral that controls them, as this
      // routine can compute all van-der Waals and electrostatic 1:4 interactions as specific
      // energy quantities.
      const int i_atom = (tinsr.x & 0x3ff);
      const int l_atom = (tinsr.y & 0x3ff);
      if (activity == VwuTask::INFR14 || activity == VwuTask::ALL_TASKS) {
        int attn_idx = ((tinsr.y >> 10) & 0x1f);
        if (attn_idx == 0 && ((tinsr.y >> 15) & 0x1)) {

          // Check whether the secondary dihedral contains the 1:4 parameter set
          attn_idx = ((tinsr.z >> 20) & 0xfff);
        }
        if (attn_idx > 0) {
          const Vec2<double> uc =
            evaluateAttenuated14Pair(i_atom, l_atom, attn_idx, nbk.coulomb_constant, sh_charges,
                                     sh_lj_idx, vk.attn14_elec, vk.attn14_vdw, nbk.lja_14_coeff,
                                     nbk.ljb_14_coeff, nbk.lj_14_sigma, nbk.n_lj_types, sh_xcrd,
                                     sh_ycrd, sh_zcrd, nullptr, nullptr, UnitCellType::NONE,
                                     sh_xfrc, sh_yfrc, sh_zfrc, eval_force, eval_force);
          if (log_term) {
            qq14_acc += llround(uc.x * nrg_scale_factor);
            lj14_acc += llround(uc.y * nrg_scale_factor);
          }
        }
        if (activity == VwuTask::INFR14) {
          continue;
        }
      }
      const bool is_charmm_improper = ((tinsr.x >> 30) & 0x1);

      // Skip CHARMM improper interactions if only standard dihedrals are of interest, and
      // skip standard dihedrals if only CHARMM impropers are of interest.
      if ((activity == VwuTask::DIHE && is_charmm_improper) ||
          (activity == VwuTask::CIMP && (! is_charmm_improper))) {
        continue;
      }
      const int j_atom = ((tinsr.x >> 10) & 0x3ff);
      const int k_atom = ((tinsr.x >> 20) & 0x3ff);
      const int param_idx = ((tinsr.y >> 16) & 0xffff);
      const bool second_term = ((tinsr.y >> 15) & 0x1);
      const TorsionKind kind = (tinsr.x >> 31) ? TorsionKind::IMPROPER : TorsionKind::PROPER;

      // Compute displacements
      double ab[3], bc[3], cd[3], crabbc[3], crbccd[3], scr[3];
      ab[0] = sh_xcrd[j_atom] - sh_xcrd[i_atom];
      ab[1] = sh_ycrd[j_atom] - sh_ycrd[i_atom];
      ab[2] = sh_zcrd[j_atom] - sh_zcrd[i_atom];
      bc[0] = sh_xcrd[k_atom] - sh_xcrd[j_atom];
      bc[1] = sh_ycrd[k_atom] - sh_ycrd[j_atom];
      bc[2] = sh_zcrd[k_atom] - sh_zcrd[j_atom];
      cd[0] = sh_xcrd[l_atom] - sh_xcrd[k_atom];
      cd[1] = sh_ycrd[l_atom] - sh_ycrd[k_atom];
      cd[2] = sh_zcrd[l_atom] - sh_zcrd[k_atom];

      // Compute cross products and then the angle between the planes
      crossProduct(ab, bc, crabbc);
      crossProduct(bc, cd, crbccd);
      crossProduct(crabbc, crbccd, scr);
      double costheta = crabbc[0]*crbccd[0] + crabbc[1]*crbccd[1] + crabbc[2]*crbccd[2];
      costheta /= sqrt((crabbc[0]*crabbc[0] + crabbc[1]*crabbc[1] + crabbc[2]*crabbc[2]) *
                       (crbccd[0]*crbccd[0] + crbccd[1]*crbccd[1] + crbccd[2]*crbccd[2]));
      costheta = (costheta < -1.0) ? -1.0 : (costheta > 1.0) ? 1.0 : costheta;
      const double theta = (scr[0]*bc[0] + scr[1]*bc[1] + scr[2]*bc[2] > 0.0) ?  acos(costheta) :
                                                                                -acos(costheta);
      double sangle, sangle_ii, stiffness, stiffness_ii;
      if (is_charmm_improper) {
        stiffness = vk.cimp_keq[param_idx];
        sangle = theta - vk.cimp_phi[param_idx];
        double contrib = stiffness * sangle * sangle;
        if (log_term) {
          cimp_acc += llround(contrib * nrg_scale_factor);
        }
      }
      else {
        const double ampl = (param_idx < 65535) ? vk.dihe_amp[param_idx] : 0.0;
        const double freq = vk.dihe_freq[param_idx];
        const double phi  = vk.dihe_phi[param_idx];
        stiffness = ampl * freq;
        sangle = (freq * theta) - phi;
        double contrib = ampl * (1.0 + cos(sangle));
        if (second_term) {
          const int param_ii_idx = (tinsr.z & 0xfffff);
          const double ampl_ii = vk.dihe_amp[param_ii_idx];
          const double freq_ii = vk.dihe_freq[param_ii_idx];
          const double phi_ii  = vk.dihe_phi[param_ii_idx];
          stiffness_ii = ampl_ii * freq_ii;
          sangle_ii = (freq_ii * theta) - phi_ii;
          contrib += ampl_ii * (1.0 + cos(sangle_ii));
        }
        if (log_term) {
          if (kind == TorsionKind::PROPER) {
            dihe_acc += llround(contrib * nrg_scale_factor);
          }
          else {
            impr_acc += llround(contrib * nrg_scale_factor);
          }
        }
      }

      // Compute forces, if requested
      if (eval_force == EvaluateForce::YES) {
        double fr;
        if (is_charmm_improper) {
          fr = -2.0 * stiffness * sangle;
        }
        else {
          fr = stiffness * sin(sangle);
          if (second_term) {
            fr += stiffness_ii * sin(sangle_ii);
          }
        }
        const double mgab = sqrt(ab[0]*ab[0] + ab[1]*ab[1] + ab[2]*ab[2]);
        const double invab = 1.0 / mgab;
        const double mgbc = sqrt(bc[0]*bc[0] + bc[1]*bc[1] + bc[2]*bc[2]);
        const double invbc = 1.0 / mgbc;
        const double mgcd = sqrt(cd[0]*cd[0] + cd[1]*cd[1] + cd[2]*cd[2]);
        const double invcd = 1.0 / mgcd;
        const double cosb = -(ab[0]*bc[0] + ab[1]*bc[1] + ab[2]*bc[2]) * invab * invbc;
        const double isinb2 = (cosb * cosb < asymptotic_to_one_lf) ?
                              fr / (1.0 - (cosb * cosb)) : fr * inverse_one_minus_asymptote_lf;
        const double cosc = -(bc[0]*cd[0] + bc[1]*cd[1] + bc[2]*cd[2]) * invbc * invcd;
        const double isinc2 = (cosc * cosc < asymptotic_to_one_lf) ?
                              fr / (1.0 - (cosc * cosc)) : fr * inverse_one_minus_asymptote_lf;
        const double invabc = invab * invbc;
        const double invbcd = invbc * invcd;
        for (int i = 0; i < 3; i++) {
          crabbc[i] *= invabc;
          crbccd[i] *= invbcd;
        }

        // Transform the rotational derivatives to cartesian coordinates
        const double fa = -invab * isinb2;
        const double fb1 = (mgbc - (mgab * cosb)) * invabc * isinb2;
        const double fb2 = cosc * invbc * isinc2;
        const double fc1 = (mgbc - (mgcd * cosc)) * invbcd * isinc2;
        const double fc2 = cosb * invbc * isinb2;
        const double fd = -invcd * isinc2;
        sh_xfrc[i_atom] += crabbc[0] * fa;
        sh_xfrc[j_atom] += (fb1 * crabbc[0]) - (fb2 * crbccd[0]);
        sh_xfrc[k_atom] += (fc2 * crabbc[0]) - (fc1 * crbccd[0]);
        sh_xfrc[l_atom] -= fd * crbccd[0];
        sh_yfrc[i_atom] += crabbc[1] * fa;
        sh_yfrc[j_atom] += (fb1 * crabbc[1]) - (fb2 * crbccd[1]);
        sh_yfrc[k_atom] += (fc2 * crabbc[1]) - (fc1 * crbccd[1]);
        sh_yfrc[l_atom] -= fd * crbccd[1];
        sh_zfrc[i_atom] += crabbc[2] * fa;
        sh_zfrc[j_atom] += (fb1 * crabbc[2]) - (fb2 * crbccd[2]);
        sh_zfrc[k_atom] += (fc2 * crabbc[2]) - (fc1 * crbccd[2]);
        sh_zfrc[l_atom] -= fd * crbccd[2];
      }
    }
  }

  // Evaluate CMAP-based terms 
  if (activity == VwuTask::CMAP || activity == VwuTask::ALL_TASKS) {
    const int ncmap = task_counts[static_cast<int>(VwuTask::CMAP)];
    const std::vector<uint> cmap_acc_mask = my_vwu.getAccumulationFlags(VwuTask::CMAP);
    for (int pos = 0; pos < ncmap; pos++) {
      bool log_term = true;
      switch (purpose) {
      case VwuGoal::ACCUMULATE:
        if (readBitFromMask(cmap_acc_mask, pos) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(cmap_acc_mask, pos);
        break;
      }
      const uint2 tinsr = my_vwu.getCmapInstruction(pos);
      const int i_atom = (tinsr.x & 0x3ff);
      const int j_atom = ((tinsr.x >> 10) & 0x3ff);
      const int k_atom = ((tinsr.x >> 20) & 0x3ff);
      const int l_atom = (tinsr.y & 0x3ff);
      const int m_atom = ((tinsr.y >> 10) & 0x3ff);
      const int surf_idx = (tinsr.y >> 20);
      
      // Assume that, by this point, the imported atoms in the ValenceWorkUnit have been
      // arranged into an image that can be trusted for all interactions.  Avoid further
      // re-imaging of displacements inside the CMAP calculation.
      const double contrib = evalCmap(vk.cmap_patches, vk.cmap_patch_bounds, surf_idx,
                                      vk.cmap_dim[surf_idx], i_atom, j_atom, k_atom, l_atom,
                                      m_atom, sh_xcrd, sh_ycrd, sh_zcrd,
                                      nullptr, nullptr, UnitCellType::NONE, sh_xfrc,
                                      sh_yfrc, sh_zfrc, eval_force);
      if (log_term) {
        cmap_acc += llround(contrib * nrg_scale_factor);
      }
    }
  }
   
  // Evaluate remaining 1:4 interactions.  Calling this routine by the "INFR14" enumeration will
  // compute all 1:4 interactions, a redefinition of "inferred 1:4" relative to other situations.
  // In calling this function, "inferred 1:4" means *inferring* 1:4 interactions from dihedrals
  // that control them as well as extra 1:4 interactions that cannot be attached to a dihedral
  // interaction.  Elsewhere, including the following conditional statement, "inferred 1:4"
  // refers to only the extra interactions (which had to be *inferred* from virtual site
  // exclusion rules).  This minor bit of confusion was accepted for convenience, to reuse
  // VwuTask and avoid making a totally separate enumerator.  The StateVariable enumerator
  // doesn't serve quite the right definitions, either.
  if (activity == VwuTask::INFR14 || activity == VwuTask::ALL_TASKS) {
    const int ninfr14 = task_counts[static_cast<int>(VwuTask::INFR14)];
    const std::vector<uint> infr_acc_mask = my_vwu.getAccumulationFlags(VwuTask::INFR14);
    for (int pos = 0; pos < ninfr14; pos++) {
      bool log_term = true;
      switch (purpose) {
      case VwuGoal::ACCUMULATE:
        if (readBitFromMask(infr_acc_mask, pos) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(infr_acc_mask, pos);
        break;
      }
      const uint tinsr = my_vwu.getInferred14Instruction(pos);
      const int i_atom = (tinsr & 0x3ff);
      const int l_atom = ((tinsr >> 10) & 0x3ff);
      const int attn_idx = (tinsr >> 20);
      if (attn_idx == 0) {
        continue;
      }
      const Vec2<double> uc =
        evaluateAttenuated14Pair(i_atom, l_atom, attn_idx, nbk.coulomb_constant, sh_charges,
                                 sh_lj_idx, vk.attn14_elec, vk.attn14_vdw, nbk.lja_14_coeff,
                                 nbk.ljb_14_coeff, nbk.lj_14_sigma, nbk.n_lj_types, sh_xcrd,
                                 sh_ycrd, sh_zcrd, nullptr, nullptr, UnitCellType::NONE, sh_xfrc,
                                 sh_yfrc, sh_zfrc, eval_force, eval_force);
      if (log_term) {
        qq14_acc += llround(uc.x * nrg_scale_factor);
        lj14_acc += llround(uc.y * nrg_scale_factor);
      }
    }
  }

  // Evaluate positional restraints (all restraint energies contribute to the same accumulator)
  const double2 steady_restraint = { 1.0, 0.0 };
  if (activity == VwuTask::RPOSN || activity == VwuTask::ALL_TASKS) {
    const int nrposn = task_counts[static_cast<int>(VwuTask::RPOSN)];
    const std::vector<uint> rposn_acc_mask = my_vwu.getAccumulationFlags(VwuTask::RPOSN);
    for (int pos = 0; pos < nrposn; pos++) {
      bool log_term = true;
      switch (purpose) {
      case VwuGoal::ACCUMULATE:
        if (readBitFromMask(rposn_acc_mask, pos) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(rposn_acc_mask, pos);
        break;
      }
      const uint2 tinsr = my_vwu.getPositionalRestraintInstruction(pos);
      const int p_atom = (tinsr.x & 0x3ff);
      const int kr_param_idx = ((tinsr.x >> 10) & 0x1fffff);
      const int xyz_param_idx = tinsr.y;
      const double contrib =
        evalPosnRestraint(p_atom, step_number, rar.rposn_init_step[kr_param_idx],
                          rar.rposn_finl_step[kr_param_idx], rar.rposn_init_xy[xyz_param_idx],
                          rar.rposn_finl_xy[xyz_param_idx], rar.rposn_init_z[xyz_param_idx],
                          rar.rposn_finl_z[xyz_param_idx], rar.rposn_init_keq[kr_param_idx],
                          rar.rposn_finl_keq[kr_param_idx], rar.rposn_init_r[kr_param_idx],
                          rar.rposn_finl_r[kr_param_idx], sh_xcrd, sh_ycrd, sh_zcrd, nullptr,
                          nullptr, UnitCellType::NONE, sh_xfrc, sh_yfrc, sh_zfrc, eval_force);
      if (log_term) {
        rest_acc += llround(contrib * nrg_scale_factor);
      }
    }
  }

  // Evaluate distance restraints
  if (activity == VwuTask::RBOND || activity == VwuTask::ALL_TASKS) {
    const int nrbond = task_counts[static_cast<int>(VwuTask::RBOND)];
    const std::vector<uint> rbond_acc_mask = my_vwu.getAccumulationFlags(VwuTask::RBOND);
    for (int pos = 0; pos < nrbond; pos++) {
      bool log_term = true;
      switch (purpose) {
      case VwuGoal::ACCUMULATE:
        if (readBitFromMask(rbond_acc_mask, pos) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(rbond_acc_mask, pos);
        break;
      }
      const uint2 tinsr = my_vwu.getDistanceRestraintInstruction(pos);
      const int i_atom = (tinsr.x & 0x3ff);
      const int j_atom = ((tinsr.x >> 10) & 0x3ff);
      const int param_idx = tinsr.y;
      const double contrib =
        evalBondRestraint<double,
                          double,
                          double,
                          double2,
                          double4>(i_atom, j_atom, step_number, rar.rbond_init_step[param_idx],
                                   rar.rbond_finl_step[param_idx], rar.rbond_init_keq[param_idx],
                                   rar.rbond_finl_keq[param_idx], rar.rbond_init_r[param_idx],
                                   rar.rbond_finl_r[param_idx], sh_xcrd, sh_ycrd, sh_zcrd, nullptr,
                                   nullptr, UnitCellType::NONE, sh_xfrc, sh_yfrc, sh_zfrc,
                                   eval_force);
      if (log_term) {
        rest_acc += llround(contrib * nrg_scale_factor);
      }
    }
  }

  // Evaluate three-point angle restraints
  if (activity == VwuTask::RANGL || activity == VwuTask::ALL_TASKS) {
    const int nrangl = task_counts[static_cast<int>(VwuTask::RANGL)];
    const std::vector<uint> rangl_acc_mask = my_vwu.getAccumulationFlags(VwuTask::RANGL);
    for (int pos = 0; pos < nrangl; pos++) {
      bool log_term = true;
      switch (purpose) {
      case VwuGoal::ACCUMULATE:
        if (readBitFromMask(rangl_acc_mask, pos) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(rangl_acc_mask, pos);
        break;
      }
      const uint2 tinsr = my_vwu.getAngleRestraintInstruction(pos);
      const int i_atom = (tinsr.x & 0x3ff);
      const int j_atom = ((tinsr.x >> 10) & 0x3ff);
      const int k_atom = ((tinsr.x >> 20) & 0x3ff);
      const int param_idx = tinsr.y;
      const double contrib =
        evalAnglRestraint<double,
                          double,
                          double,
                          double2,
                          double4>(i_atom, j_atom, k_atom, step_number,
                                   rar.rangl_init_step[param_idx], rar.rangl_finl_step[param_idx],
                                   rar.rangl_init_keq[param_idx], rar.rangl_finl_keq[param_idx],
                                   rar.rangl_init_r[param_idx], rar.rangl_finl_r[param_idx],
                                   sh_xcrd, sh_ycrd, sh_zcrd, nullptr, nullptr, UnitCellType::NONE,
                                   sh_xfrc, sh_yfrc, sh_zfrc, eval_force);
      if (log_term) {
        rest_acc += llround(contrib * nrg_scale_factor);
      }
    }
  }

  // Evaluate four-point dihedral restraints
  if (activity == VwuTask::RDIHE || activity == VwuTask::ALL_TASKS) {
    const int nrdihe = task_counts[static_cast<int>(VwuTask::RDIHE)];
    const std::vector<uint> rdihe_acc_mask = my_vwu.getAccumulationFlags(VwuTask::RDIHE);
    for (int pos = 0; pos < nrdihe; pos++) {
      bool log_term = true;
      switch (purpose) {
      case VwuGoal::ACCUMULATE:
        if (readBitFromMask(rdihe_acc_mask, pos) == 0) {
          continue;
        }
        break;
      case VwuGoal::MOVE_PARTICLES:
        log_term = readBitFromMask(rdihe_acc_mask, pos);
        break;
      }
      const uint2 tinsr = my_vwu.getDihedralRestraintInstruction(pos);
      const int i_atom = (tinsr.x & 0x3ff);
      const int j_atom = ((tinsr.x >> 10) & 0x3ff);
      const int k_atom = ((tinsr.x >> 20) & 0x3ff);
      const int l_atom = (tinsr.y & 0x3ff);
      const int param_idx = ((tinsr.y >> 10) & 0x3fffff);
      const double contrib =
        evalDiheRestraint<double,
                          double,
                          double,
                          double2,
                          double4>(i_atom, j_atom, k_atom, l_atom, step_number,
                                   rar.rdihe_init_step[param_idx], rar.rdihe_finl_step[param_idx],
                                   rar.rdihe_init_keq[param_idx], rar.rdihe_finl_keq[param_idx],
                                   rar.rdihe_init_r[param_idx], rar.rdihe_finl_r[param_idx],
                                   sh_xcrd, sh_ycrd, sh_zcrd, nullptr, nullptr, UnitCellType::NONE,
                                   sh_xfrc, sh_yfrc, sh_zfrc, eval_force);
      if (log_term) {
        rest_acc += llround(contrib * nrg_scale_factor);
      }
    }
  }  

  // Contribute results as the instantaneous states
  commitVwuEnergies(bond_acc, angl_acc, dihe_acc, impr_acc, ubrd_acc, cimp_acc, cmap_acc,
                    qq14_acc, lj14_acc, rest_acc, sysid, activity, ecard);
}

//-------------------------------------------------------------------------------------------------
void evalValenceWorkUnits(const ValenceKit<double> vk, const VirtualSiteKit<double> vsk,
                          const NonbondedKit<double> nbk,
                          const RestraintKit<double, double2, double4> rar, double* xcrd,
                          double* ycrd, double* zcrd, const double* umat, const double* invu,
                          const UnitCellType unit_cell, double* xfrc, double* yfrc, double* zfrc,
                          ScoreCard *ecard, const int sysid,
                          const std::vector<ValenceWorkUnit> &vwu_list,
                          const EvaluateForce eval_force, const VwuTask activity,
                          const VwuGoal purpose, const int step_number) {  

  // Initialize accumulators
  evalVwuInitEnergy(ecard, activity, sysid);

  // Declare arrays outside the loop to avoid unecessary reallocations
  std::vector<double> sh_xcrd, sh_ycrd, sh_zcrd, sh_xfrc, sh_yfrc, sh_zfrc, sh_charges;
  std::vector<int> sh_lj_idx;
  
  // Loop over each work unit
  for (size_t vidx = 0LLU; vidx < vwu_list.size(); vidx++) {

    // Make local arrays of the coordinates and forces.  Import coordinates and extant forces
    // from the global arrays.
    const int max_atoms = vwu_list[vidx].getMaxAtoms();
    sh_xcrd.resize(max_atoms);
    sh_ycrd.resize(max_atoms);
    sh_zcrd.resize(max_atoms);
    sh_charges.resize(max_atoms);
    sh_lj_idx.resize(max_atoms);
    if (eval_force == EvaluateForce::YES) {
      sh_xfrc.resize(max_atoms);
      sh_yfrc.resize(max_atoms);
      sh_zfrc.resize(max_atoms);
    }
    const int n_imp_atoms = vwu_list[vidx].getImportedAtomCount();
    for (int i = 0; i < n_imp_atoms; i++) {
      const size_t atom_idx = vwu_list[vidx].getImportedAtomIndex(i);
      sh_xcrd[i] = xcrd[atom_idx];
      sh_ycrd[i] = ycrd[atom_idx];
      sh_zcrd[i] = zcrd[atom_idx];
      sh_charges[i] = nbk.charge[atom_idx];
      sh_lj_idx[i] = nbk.lj_idx[atom_idx];
    }
    if (eval_force == EvaluateForce::YES) {
      for (int i = 0; i < n_imp_atoms; i++) {    
        sh_xfrc[i] = 0.0;
        sh_yfrc[i] = 0.0;
        sh_zfrc[i] = 0.0;
      }
    }
    
    // Evaluate the work unit using the locally cached data.  If force evaluation is not requested,
    // it still works to feed the unallocated vectors' data pointers as they will not be used.
    localVwuEvaluation(vk, vsk, nbk, rar, sh_charges.data(), sh_lj_idx.data(), sh_xcrd.data(),
                       sh_ycrd.data(), sh_zcrd.data(), sh_xfrc.data(), sh_yfrc.data(),
                       sh_zfrc.data(), ecard, sysid, vwu_list[vidx], eval_force, activity, purpose,
                       step_number);
    switch (purpose) {
    case VwuGoal::ACCUMULATE:

      // Add accumulated forces back to the global arrays (this is not done by all GPU kernels, as
      // in some cases the ValenceWorkUnits also move atoms and then leave the global force arrays
      // zero'ed for atoms they are responsible for moving).  It is only possible to write these
      // forces back because one and only one work unit is assigned to log the force due to any
      // particular interaction.  If any two work units import the same atoms for any reason (i.e.
      // they are needed to support certain forces that both work units must know in order to move
      // some of their imported atoms) it is not possible to move atoms in a valid way and write
      // back the accumulated forces.  GPU kernels will generally move atoms and not write back
      // forces.  This function writes back forces but does not (cannot, safely) move atoms.
      for (int i = 0; i < n_imp_atoms; i++) {
        const size_t atom_idx = vwu_list[vidx].getImportedAtomIndex(i);
        xfrc[atom_idx] += sh_xfrc[i];
        yfrc[atom_idx] += sh_yfrc[i];
        zfrc[atom_idx] += sh_zfrc[i];
      }
      break;
    case VwuGoal::MOVE_PARTICLES:

      // Copy back the modified, cached coordinates if this work unit is responsible for moving
      // the particle.
      const std::vector<uint> updates = vwu_list[vidx].getAtomUpdateFlags();
      for (int i = 0; i < n_imp_atoms; i++) {
        if (readBitFromMask(updates, i) == 1) {
          const size_t atom_idx = vwu_list[vidx].getImportedAtomIndex(i);
          xcrd[atom_idx] = sh_xcrd[i];
          ycrd[atom_idx] = sh_ycrd[i];
          zcrd[atom_idx] = sh_zcrd[i];
        }
      }
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void evalValenceWorkUnits(const AtomGraph *ag, PhaseSpace *ps, const RestraintApparatus *ra,
                          ScoreCard *ecard, const int sysid,
                          const std::vector<ValenceWorkUnit> &vwu_list,
                          const EvaluateForce eval_force, const VwuTask activity,
                          const VwuGoal goal, const int step_number) {
  PhaseSpaceWriter psw = ps->data();
  evalValenceWorkUnits(ag->getDoublePrecisionValenceKit(), ag->getDoublePrecisionVirtualSiteKit(),
                       ag->getDoublePrecisionNonbondedKit(), ra->dpData(), psw.xcrd, psw.ycrd,
                       psw.zcrd, psw.umat, psw.invu, psw.unit_cell, psw.xfrc, psw.yfrc, psw.zfrc,
                       ecard, sysid, vwu_list, eval_force, activity, step_number);
}

//-------------------------------------------------------------------------------------------------
void evalValenceWorkUnits(const ValenceKit<double> vk, const VirtualSiteKit<double> vsk,
                          const NonbondedKit<double> nbk,
                          const RestraintKit<double, double2, double4> rar, const double* xcrd,
                          const double* ycrd, const double* zcrd, const double* umat,
                          const double* invu, const UnitCellType unit_cell, double* xfrc,
                          double* yfrc, double* zfrc, ScoreCard *ecard, const int sysid,
                          const std::vector<ValenceWorkUnit> &vwu_list,
                          const EvaluateForce eval_force, const VwuTask activity,
                          const int step_number) {  

  // Initialize accumulators
  evalVwuInitEnergy(ecard, activity, sysid);

  // Declare arrays outside the loop to avoid unecessary reallocations
  std::vector<double> sh_xcrd, sh_ycrd, sh_zcrd, sh_xfrc, sh_yfrc, sh_zfrc, sh_charges;
  std::vector<int> sh_lj_idx;

  // Loop over each work unit
  for (size_t vidx = 0LLU; vidx < vwu_list.size(); vidx++) {

    // Make local arrays of the coordinates and forces.  Import coordinates and extant forces
    // from the global arrays.
    const int max_atoms = vwu_list[vidx].getMaxAtoms();
    sh_xcrd.resize(max_atoms);
    sh_ycrd.resize(max_atoms);
    sh_zcrd.resize(max_atoms);
    sh_xfrc.resize(max_atoms);
    sh_yfrc.resize(max_atoms);
    sh_zfrc.resize(max_atoms);
    sh_charges.resize(max_atoms);
    sh_lj_idx.resize(max_atoms);
    const int n_imp_atoms = vwu_list[vidx].getImportedAtomCount();
    for (int i = 0; i < n_imp_atoms; i++) {
      const size_t atom_idx = vwu_list[vidx].getImportedAtomIndex(i);
      sh_xcrd[i] = xcrd[atom_idx];
      sh_ycrd[i] = ycrd[atom_idx];
      sh_zcrd[i] = zcrd[atom_idx];
      sh_xfrc[i] = 0.0;
      sh_yfrc[i] = 0.0;
      sh_zfrc[i] = 0.0;
      sh_charges[i] = nbk.charge[atom_idx];
      sh_lj_idx[i] = nbk.lj_idx[atom_idx];
    }

    // Evaluate the work unit using the locally cached data
    localVwuEvaluation(vk, vsk, nbk, rar, sh_charges.data(), sh_lj_idx.data(), sh_xcrd.data(),
                       sh_ycrd.data(), sh_zcrd.data(), sh_xfrc.data(), sh_yfrc.data(),
                       sh_zfrc.data(), ecard, sysid, vwu_list[vidx], eval_force, activity,
                       VwuGoal::ACCUMULATE, step_number);
    
    // Add accumulated forces back to the global arrays (this is not done by all GPU kernels, as
    // in some cases the ValenceWorkUnits also move atoms and then leave the global force arrays
    // zero'ed for atoms they are responsible for moving).  It is only possible to write these
    // forces back because one and only one work unit is assigned to log the force due to any
    // particular interaction.  If any two work units import the same atoms for any reason (i.e.
    // they are needed to support certain forces that both work units must know in order to move
    // some of their imported atoms) it is not possible to move atoms in a valid way and write
    // back the accumulated forces.  GPU kernels will generally move atoms and not write back
    // forces.  This function writes back forces but does not (cannot, safely) move atoms.
    for (int i = 0; i < n_imp_atoms; i++) {
      const size_t atom_idx = vwu_list[vidx].getImportedAtomIndex(i);
      xfrc[atom_idx] += sh_xfrc[i];
      yfrc[atom_idx] += sh_yfrc[i];
      zfrc[atom_idx] += sh_zfrc[i];
    }
  }
}

//-------------------------------------------------------------------------------------------------
void evalValenceWorkUnits(const AtomGraph &ag, const PhaseSpace &ps, const RestraintApparatus &ra,
                          ScoreCard *ecard, const int sysid,
                          const std::vector<ValenceWorkUnit> &vwu_list,
                          const VwuTask activity, const int step_number) {
  PhaseSpaceReader psr = ps.data();
  evalValenceWorkUnits(ag.getDoublePrecisionValenceKit(), ag.getDoublePrecisionVirtualSiteKit(),
                       ag.getDoublePrecisionNonbondedKit(), ra.dpData(), psr.xcrd, psr.ycrd,
                       psr.zcrd, psr.umat, psr.invu, psr.unit_cell, nullptr, nullptr, nullptr,
                       ecard, sysid, vwu_list, EvaluateForce::NO, activity, step_number);
}

} // namespace energy
} // namespace stormm
