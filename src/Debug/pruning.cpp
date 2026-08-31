#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Parsing/parse.h"
#include "Potential/valence_potential.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/coordinateframe.h"
#include "pruning.h"

namespace stormm {
namespace debug {

using card::HybridFormat;
using energy::evalHarmonicBend;
using energy::evalHarmonicStretch;
using energy::EvaluateForce;
using parse::char4ToString;
using topology::AtomGraph;
using topology::ChemicalDetailsKit;
using topology::NonbondedKit;
using topology::ValenceKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;

//-------------------------------------------------------------------------------------------------
std::vector<int> checkBondAndAngleSanity(const PhaseSpaceSynthesis *poly_ps,
                                         const DebugControls &dbgcon) {
  const double max_bond_strain = dbgcon.getBondStrainTolerance();
  const double max_angle_strain = dbgcon.getAngleStrainTolerance();
  std::vector<int> result;
  result.reserve(poly_ps->getSystemCount());
  for (int sysid = 0; sysid < poly_ps->getSystemCount(); sysid++) {
    const AtomGraph* ag = poly_ps->getSystemTopologyPointer(sysid);
    const CoordinateFrame cf = poly_ps->exportCoordinates(sysid, HybridFormat::HOST_ONLY);
    const ValenceKit vk = ag->getDoublePrecisionValenceKit();
    const NonbondedKit nbk = ag->getDoublePrecisionNonbondedKit();
    const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
    const CoordinateFrameReader cfr = cf.data();

    // Check bonds
    std::vector<double> bond_violations;
    std::vector<int2> violating_pairs;
    for (int pos = 0; pos < vk.nbond; pos++) {
      const int atom_i = vk.bond_i_atoms[pos];
      const int atom_j = vk.bond_j_atoms[pos];
      const int parm_idx = vk.bond_param_idx[pos];
      const double k_eq = vk.bond_keq[parm_idx];
      const double l_eq = vk.bond_leq[parm_idx];
      const double u_t = evalHarmonicStretch<double,
                                             double, double>(atom_i, atom_j, k_eq, l_eq, cfr.xcrd,
                                                             cfr.ycrd, cfr.zcrd, cfr.umat,
                                                             cfr.invu, cfr.unit_cell, nullptr,
                                                             nullptr, nullptr, EvaluateForce::NO);
      if (u_t > max_bond_strain) {
        bond_violations.push_back(u_t);
        violating_pairs.push_back({ atom_i, atom_j });
      }
    }
    const int nbond_violations = bond_violations.size();

    // Check angles in a similar fashion
    std::vector<double> angle_violations;
    std::vector<int3> violating_triples;
    for (int pos = 0; pos < vk.nangl; pos++) {
      const int atom_i = vk.angl_i_atoms[pos];
      const int atom_j = vk.angl_j_atoms[pos];
      const int atom_k = vk.angl_k_atoms[pos];
      const int parm_idx = vk.angl_param_idx[pos];
      const double k_eq = vk.angl_keq[parm_idx];
      const double t_eq = vk.angl_theta[parm_idx];
      const double u_t = evalHarmonicBend<double,
                                          double, double>(atom_i, atom_j, atom_k, k_eq, t_eq,
                                                          cfr.xcrd, cfr.ycrd, cfr.zcrd, cfr.umat,
                                                          cfr.invu, cfr.unit_cell, nullptr,
                                                          nullptr, nullptr, EvaluateForce::NO);
      if (u_t > max_angle_strain) {
        angle_violations.push_back(u_t);
        violating_triples.push_back({ atom_i, atom_j, atom_k });
      }
    }
    for (int pos = 0; pos < vk.nubrd; pos++) {
      const int atom_i = vk.ubrd_i_atoms[pos];
      const int atom_k = vk.ubrd_k_atoms[pos];
      const int parm_idx = vk.ubrd_param_idx[pos];
      const double k_eq = vk.ubrd_keq[parm_idx];
      const double l_eq = vk.ubrd_leq[parm_idx];
      const double u_t = evalHarmonicStretch<double,
                                             double, double>(atom_i, atom_k, k_eq, l_eq, cfr.xcrd,
                                                             cfr.ycrd, cfr.zcrd, cfr.umat,
                                                             cfr.invu, cfr.unit_cell, nullptr,
                                                             nullptr, nullptr, EvaluateForce::NO);
      if (u_t > max_angle_strain) {
        angle_violations.push_back(u_t);

        // Find the atom between i and k
        int atom_j = 0;
        for (int j = nbk.nb12_bounds[atom_i]; j < nbk.nb12_bounds[atom_i + 1]; j++) {
          const int guess_j = nbk.nb12x[j];
          bool link_k = false;
          for (int k = nbk.nb12_bounds[guess_j]; k < nbk.nb12_bounds[guess_j + 1]; k++) {
            link_k = (nbk.nb12x[j] == atom_k || link_k);
          }
          if (link_k) {
            atom_j = guess_j;
          }
        }

        // Even if the above search for the connecting atom fails, the tuple can still be created.
        violating_triples.push_back({ atom_i, atom_j, atom_k });
      }
    }
    const int nangl_violations = angle_violations.size();

    // Analyze the violations, if there are any.  Are the violations co-localized?  Do the atoms
    // share anything in common, such as being members of the same ring group?
    if (nbond_violations > 0) {
      for (int i = 0; i < nbond_violations; i++) {
        const int atom_i = violating_pairs[i].x;
        const int atom_j = violating_pairs[i].y;
        const int res_i = ag->getResidueIndex(atom_i);
        const int res_j = ag->getResidueIndex(atom_j);
        const std::string atom_i_name = char4ToString(cdk.atom_names[atom_i]);
        const std::string atom_j_name = char4ToString(cdk.atom_names[atom_j]);
        const std::string res_i_name = char4ToString(cdk.res_names[res_i]);
        const std::string res_j_name = char4ToString(cdk.res_names[res_j]);
        printf("Highly strained bond in system %2d : %4d %4.4s %4.4s -- %4d %4.4s %4.4s  "
               "%12.4lf\n", sysid, atom_i, atom_i_name.c_str(), res_i_name.c_str(), atom_j,
               atom_j_name.c_str(), res_j_name.c_str(), bond_violations[i]);
      }
    }
    if (nangl_violations > 0) {
      for (int i = 0; i < nangl_violations; i++) {
        const int atom_i = violating_triples[i].x;
        const int atom_j = violating_triples[i].y;
        const int atom_k = violating_triples[i].z;
        const int res_i = ag->getResidueIndex(atom_i);
        const int res_j = ag->getResidueIndex(atom_j);
        const int res_k = ag->getResidueIndex(atom_k);
        const std::string atom_i_name = char4ToString(cdk.atom_names[atom_i]);
        const std::string atom_j_name = char4ToString(cdk.atom_names[atom_j]);
        const std::string atom_k_name = char4ToString(cdk.atom_names[atom_k]);
        const std::string res_i_name = char4ToString(cdk.res_names[res_i]);
        const std::string res_j_name = char4ToString(cdk.res_names[res_j]);
        const std::string res_k_name = char4ToString(cdk.res_names[res_k]);
        printf("Highly strained angle in system %2d : %4d %4.4s %4.4s -- %4d %4.4s %4.4s -- "
               "%4d %4.4s %4.4s  %12.4lf\n", sysid, atom_i, atom_i_name.c_str(),
               res_i_name.c_str(), atom_j, atom_j_name.c_str(), res_j_name.c_str(),
               atom_k, atom_k_name.c_str(), res_k_name.c_str(), angle_violations[i]);
      }
    }
    
    // If there are no violations, include this system in the final list.
    if (nbond_violations == 0 && nangl_violations == 0) {
      result.push_back(sysid);
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> checkBondAndAngleSanity(const PhaseSpaceSynthesis &poly_ps,
                                         const DebugControls &dbgcon) {
  return checkBondAndAngleSanity(poly_ps.getSelfPointer(), dbgcon);
}

} // namespace debug
} // namespace stormm
