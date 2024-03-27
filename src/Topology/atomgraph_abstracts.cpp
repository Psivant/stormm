#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"
#include "atomgraph_abstracts.h"

namespace stormm {
namespace topology {

//-------------------------------------------------------------------------------------------------
ChemicalDetailsKit::ChemicalDetailsKit(const int natom_in, const int nres_in, const int nmol_in,
                                       const int free_dof_in, const int cnst_dof_in,
                                       const char4* atom_names_in, const char4* res_names_in,
                                       const char4* atom_types_in, const int* z_numbers_in,
                                       const int* res_limits_in, const int* atom_numbers_in,
                                       const int* res_numbers_in, const int* mol_home_in,
                                       const int* mol_contents_in, const int* mol_limits_in,
                                       const double* masses_in, const float* sp_masses_in,
                                       const double* inv_masses_in,
                                       const float* sp_inv_masses_in) :
    natom{natom_in}, nres{nres_in}, nmol{nmol_in}, free_dof{free_dof_in}, cnst_dof{cnst_dof_in},
    atom_names{atom_names_in}, res_names{res_names_in}, atom_types{atom_types_in},
    z_numbers{z_numbers_in}, res_limits{res_limits_in}, atom_numbers{atom_numbers_in},
    res_numbers{res_numbers_in}, mol_home{mol_home_in}, mol_contents{mol_contents_in},
    mol_limits{mol_limits_in}, masses{masses_in}, sp_masses{sp_masses_in},
    inv_masses{inv_masses_in}, sp_inv_masses{sp_inv_masses_in}
{}

} // namespace topology
} // namespace stormm
