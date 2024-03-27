// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace structure {

//-------------------------------------------------------------------------------------------------
template <typename T>
RMSDPlanReader<T>::RMSDPlanReader(const int plan_count_in, const RMSDMethod strategy_in,
                                  const double mass_fraction_in, const T* masses_in,
                                  const int* atom_counts_in, const int* atom_starts_in,
                                  const int* alignment_steps_in, const int* core_atoms_in,
                                  const int* core_counts_in, const int* core_starts_in,
                                  const int* symm_atoms_in, const int4* symm_bounds_in,
                                  const int2* symm_ranges_in) :
    plan_count{plan_count_in}, strategy{strategy_in}, mass_fraction{mass_fraction_in},
    masses{masses_in}, atom_counts{atom_counts_in}, atom_starts{atom_starts_in},
    alignment_steps{alignment_steps_in}, core_atoms{core_atoms_in}, core_counts{core_counts_in},
    core_starts{core_starts_in}, symm_atoms{symm_atoms_in}, symm_bounds{symm_bounds_in},
    symm_ranges{symm_ranges_in}
{}

} // namespace structure
} // namespace stormm
