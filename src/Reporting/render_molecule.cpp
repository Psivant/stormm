#include "copyright.h"
#include "render_molecule.h"

namespace stormm {
namespace review {

//-------------------------------------------------------------------------------------------------
std::vector<int> traceBondLines(const NonbondedKit<double> &nbk, const ChemicalDetailsKit &cdk,
                                const int mol_index, std::vector<uint> *atom_coverage,
                                std::vector<uint> *bond_coverage) {
  std::vector<int> result, history;
  const int mol_size = cdk.mol_limits[mol_index + 1] - cdk.mol_limits[mol_index];
  history.reserve(mol_size);
  result.reserve(2 * mol_size);
  int present_atom = cdk.mol_contents[cdk.mol_limits[mol_index]];
  uint* atom_coverage_ptr = atom_coverage->data();
  uint* bond_coverage_ptr = bond_coverage->data();
  accumulateBitmask(atom_coverage_ptr, present_atom);
  result.push_back(present_atom);
  bool all_bonds_tried = false;
  int atoms_covered = 1;
  while (atoms_covered < mol_size || history.size() > 0) {

    // The first priority is to find an unused bond to an atom that the path has not yet touched.
    // In lieu of that, it is preferable to travel next to an atom that has been included in the
    // path already, but is not the previous atom.  As a last resort, return to the previous atom.
    int new_atom_path = -1;
    int known_atom_path = -1;
    for (int i = nbk.nb12_bounds[present_atom]; i < nbk.nb12_bounds[present_atom + 1]; i++) {
      if (readBitFromMask(bond_coverage_ptr, i) == 0) {
        if (readBitFromMask(atom_coverage_ptr, nbk.nb12x[i]) == 0) {
          new_atom_path = i;
        }
        else if (history.size() == 0 || nbk.nb12x[i] != history.back()) {
          known_atom_path = i;
        }
      }
    }
    if (new_atom_path >= 0) {

      // Extend the path to include a new atom.
      accumulateBitmask(bond_coverage_ptr, new_atom_path);
      
      // The bond has been logged, but the history can only move forward if the atom to which
      // the bond extends is not yet covered.  Otherwise, begin to backtrack on the next cycle.
      history.push_back(present_atom);
      atoms_covered++;
      present_atom = nbk.nb12x[new_atom_path];
      accumulateBitmask(atom_coverage_ptr, present_atom);
      result.push_back(present_atom);
    }
    else if (known_atom_path >= 0) {

      // Extend the path to the atom, which is already included on the path in some other capacity,
      // then immediately backtrack.
      accumulateBitmask(bond_coverage_ptr, known_atom_path);
      const int known_atom = nbk.nb12x[known_atom_path];
      result.push_back(known_atom);
      for (int i = nbk.nb12_bounds[known_atom]; i < nbk.nb12_bounds[known_atom + 1]; i++) {
        if (nbk.nb12x[i] == present_atom) {
          accumulateBitmask(bond_coverage_ptr, i);
        }
      }
      result.push_back(present_atom);
    }
    else {

      // Backtrack.
      const int last_atom = history.back();
      for (int i = nbk.nb12_bounds[present_atom]; i < nbk.nb12_bounds[present_atom + 1]; i++) {
        if (nbk.nb12x[i] == last_atom) {
          accumulateBitmask(bond_coverage_ptr, i);
        }
      }
      present_atom = last_atom;
      result.push_back(present_atom);
      history.pop_back();
    }
  }
  return result;
}
  
//-------------------------------------------------------------------------------------------------
void renderMolecule(std::ofstream *foutp, const CoordinateFrame *cf, const AtomGraph *ag,
                    const RenderOptions *ropt, const GridFileSyntax syntax) {
  const CoordinateFrameReader cfr = cf->data();
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  renderMolecule<double>(foutp, cfr.xcrd, cfr.ycrd, cfr.zcrd, nbk, cdk, *ropt, syntax, 1);
}

//-------------------------------------------------------------------------------------------------
void renderMolecule(std::ofstream *foutp, const CoordinateFrame &cf, const AtomGraph &ag,
                    const RenderOptions &ropt, const GridFileSyntax syntax) {
  renderMolecule(foutp, cf.getSelfPointer(), ag.getSelfPointer(), ropt.getSelfPointer(), syntax);
}

} // namespace review
} // namespace stormm
