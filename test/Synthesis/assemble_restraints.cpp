#include "copyright.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/Random/random.h"
#include "../../src/Topology/atomgraph_abstracts.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "assemble_restraints.h"

#ifndef STORMM_USE_HPC
using stormm::data_types::double3;
#endif
using stormm::random::Xoroshiro128pGenerator;
using stormm::restraints::BoundedRestraint;
using stormm::topology::AtomGraph;
using stormm::topology::ValenceKit;
using stormm::trajectory::CoordinateFrameReader;
using stormm::trajectory::PhaseSpace;

//-------------------------------------------------------------------------------------------------
// Create a collection of positional, distance, angle, and dihedral restraints acting on a system
// which will exert mild to moderate forces.  All of the restraints, together, will be checked for
// total energy and forces.  This feeds into tests that are not designed to weed out individual
// problems but instead take a "pooled sample" and verify that there are no problems in it.
//
// Arguments:
//   ag:      System topology
//   ps:      System coordinates
//-------------------------------------------------------------------------------------------------
RestraintApparatus assembleRestraints(const AtomGraph *ag, const PhaseSpace &ps) {
  const CoordinateFrameReader cfr(ps);
  std::vector<BoundedRestraint> rlist;
  rlist.reserve((cfr.natom / 8) + 12);
  Xoroshiro128pGenerator xrs(87293);
  int nrst = 0;
  for (int i = 0; i < cfr.natom; i += cfr.natom / 8) {
    rlist.emplace_back(i, ag, cfr, 1.1, 1.4, 0.2, 0.6, 0.8, 1.4);
    double3 ts = rlist[nrst].getTargetSite();
    ts.x += 0.25 * (0.5 - xrs.uniformRandomNumber());
    ts.y += 0.25 * (0.5 - xrs.uniformRandomNumber());
    ts.z += 0.25 * (0.5 - xrs.uniformRandomNumber());
    rlist[nrst].setTargetSite(ts);
    nrst++;
  }
  const ValenceKit<double> vk = ag->getDoublePrecisionValenceKit();
  for (int pos = 0; pos < vk.nbond; pos += vk.nbond / 4) {
    const int i_atom = vk.bond_i_atoms[pos];
    const int j_atom = vk.bond_j_atoms[pos];
    rlist.emplace_back(i_atom, j_atom, ag, 2.4, 2.7, 0.5, 2.0, 2.4, 3.9);
  }
  for (int pos = 0; pos < vk.nangl; pos += vk.nangl / 4) {
    const int i_atom = vk.angl_i_atoms[pos];
    const int j_atom = vk.angl_j_atoms[pos];
    const int k_atom = vk.angl_k_atoms[pos];
    rlist.emplace_back(i_atom, j_atom, k_atom, ag, 3.1, 1.7, 1.5, 2.1, 2.1, 2.7);
  }
  for (int pos = 0; pos < vk.ndihe; pos += vk.ndihe / 4) {
    const int i_atom = vk.dihe_i_atoms[pos];
    const int j_atom = vk.dihe_j_atoms[pos];
    const int k_atom = vk.dihe_k_atoms[pos];
    const int l_atom = vk.dihe_l_atoms[pos];
    rlist.emplace_back(i_atom, j_atom, k_atom, l_atom, ag, 3.1, 1.7, 1.5, 2.1, 2.1, 2.7);
  }
  return RestraintApparatus(rlist, ag);
}
