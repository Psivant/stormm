#include <cmath>
#include "copyright.h"
#include "Constants/symbol_values.h"
#include "DataTypes/stormm_vector_types.h"
#include "Numerics/split_fixed_precision.h"
#include "isomerization.h"
#include "local_arrangement.h"
#include "structure_enumerators.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(CoordinateFrame *cf, const int atom_i, const int atom_j,
                     const std::vector<int> &moving_atoms, const double rotation_angle) {
  rotateAboutBond(cf->data(), atom_i, atom_j, moving_atoms, rotation_angle);
}

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(CoordinateFrameWriter cfw, const int atom_i, const int atom_j,
                     const std::vector<int> &moving_atoms, const double rotation_angle) {
  rotateAboutBond(cfw.xcrd, cfw.ycrd, cfw.zcrd, atom_i, atom_j, moving_atoms.data(),
                  moving_atoms.size(), rotation_angle);
}

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(PhaseSpace *ps, const int atom_i, const int atom_j,
                     const std::vector<int> &moving_atoms, const double rotation_angle) {
  rotateAboutBond(ps->data(), atom_i, atom_j, moving_atoms, rotation_angle);
}

//-------------------------------------------------------------------------------------------------
void rotateAboutBond(PhaseSpaceWriter psw, const int atom_i, const int atom_j,
                     const std::vector<int> &moving_atoms, const double rotation_angle) {
  rotateAboutBond(psw.xcrd, psw.ycrd, psw.zcrd, atom_i, atom_j, moving_atoms.data(),
                  moving_atoms.size(), rotation_angle);
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(CoordinateFrame *cf, const int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups) {
  flipChiralCenter(cf->data(), center_idx, chiral_centers, chiral_protocols, inversion_groups);
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(CoordinateFrameWriter cfw, const int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups) {
  flipChiralCenter<double, double>(cfw.xcrd, cfw.ycrd, cfw.zcrd, center_idx, chiral_centers,
                                   chiral_protocols, inversion_groups);
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(PhaseSpace *ps, const int center_idx, const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups) {
  flipChiralCenter(ps->data(), center_idx, chiral_centers, chiral_protocols, inversion_groups);
}

//-------------------------------------------------------------------------------------------------
void flipChiralCenter(PhaseSpaceWriter psw, const int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups) {
  flipChiralCenter<double, double>(psw.xcrd, psw.ycrd, psw.zcrd, center_idx, chiral_centers,
                                   chiral_protocols, inversion_groups);
}

} // namespace structure
} // namespace stormm
