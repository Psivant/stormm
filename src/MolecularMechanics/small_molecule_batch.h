// -*-c++-*-
#ifndef STORMM_SMALL_MOLECULE_BATCH_H
#define STORMM_SMALL_MOLECULE_BATCH_H

#include "copyright.h"

namespace stormm {
namespace mm {

/// \brief Atom limits for small molecules in various thread block sizes
/// \{
constexpr int max_smmol_atom_count_tiny_block   = 96;
constexpr int max_smmol_atom_count_small_block  = 192;
constexpr int max_smmol_atom_count_medium_block = 384;
constexpr int max_smmol_atom_count_large_block  = 720;
/// \}

/// \brief Minimize a batch of small molecules.  The CPU templated function mimics the GPU
///        templated kernels to minimize the small molecules, subject to restraints, one at a
///        time whereas the GPU will assign one molecule to a thread block at a time and continue
///        until all thread blocks have finished with the work.


} // namespace mm
} // namespace stormm

#include "small_molecule_batch.tpp"

#endif
