// -*-c++-*-
#include "copyright.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"

namespace stormm {
namespace mm {

//-------------------------------------------------------------------------------------------------
template <typename T> void minimizeSmallMoleculeBatch<T>(PhaseSpaceSynthesis *psynth,
                                                         const AtomGraphSynthesis &agsynth,
                                                         const 
                                                         const BlockSize block_simulator) {

  // Local data for each minimization.  Each thread block on the card will allocate similar
  // arrays.  
  std::vector<int> descriptors(ValenceWorkUnitSpecs::DESCRIPTOR_COUNT);
  int max_atom_count;
  switch (block_simulator) {
  case BlockSize::TINY:
    max_atom_count = max_smmol_atom_count_tiny_block;
    break;
  case BlockSize::SMALL:
    max_atom_count = max_smmol_atom_count_small_block;
    break;
  case BlockSize::MEDIUM:
    max_atom_count = max_smmol_atom_count_medium_block;
    break;
  case BlockSize::LARGE:
    max_atom_count = max_smmol_atom_count_large_block;
    break;
  }

  // The coordinates will be stored in high precision, 28 bits after the decimal (268 million
  // parts per Angstrom).
  std::vector<llint> xcrd(max_atom_count);
  std::vector<llint> ycrd(max_atom_count);
  std::vector<llint> zcrd(max_atom_count);

  // Forces will be accumulated in two tiers: most of the time, the first 32 bits will be
  // sufficient, particularly if the forces are accumulated with the default 23 bits of precision
  // after the decimal.  If the force is large, or if the precision model is tighter, the overflow
  // bits will be brought to bear.  The overflow will be stored off-chip, in GMEM, but the first
  // 32 bits will be stored in L1 __shared__ memory.  With the forces known, the particle moves
  // can be computed and reversed exactly for sampling along the gradient.
  std::vector<int> xfrc(max_atom_count);
  std::vector<int> yfrc(max_atom_count);
  std::vector<int> zfrc(max_atom_count);
  std::vector<int> xfrc_overflow(max_atom_count);
  std::vector<int> yfrc_overflow(max_atom_count);
  std::vector<int> zfrc_overflow(max_atom_count);

  // Arrays for stroing aspects of the Generalized Born calculation
  std::vector<T> psi(max_atom_count);
  std::vector<T> subdeijda(max_atom_count);
  std::vector<T> reff(max_atom_count);

  // Prepare abstracts for reading the topology and writing the coordinate synthesis objects
  PsSynthesisWriter psynthw = psynth->data();
  AgSynthesisReader agsynthr = agsynth->data();
  
  // Loop over all systems and accomplish each minimization one by one.  On the accelerator card,
  // one thread block will handle each system until its minimization is completed, in order to
  // conserve L1 cache content.
  for (int syscon = 0; syscon < psynth.system_count; syscon++) {
  
    // Load system stats.  A __syncthreads() must happen after this step in HPC execution.
    for (int i = 0; i < ValenceWorkUnitSpecs::DESCRIPTOR_COUNT; i++) {
      descriptors[i] = psynth.system_specs[(syscon * ValenceWorkUnitSpecs::DESCRIPTOR_COUNT) + i];
    }

    // Take in atom coordinates and initialize forces.  Another __syncthreads() must happen after
    // this step.
    const int atom_offset = descriptors[ValenceWorkUnitSpecs::ATOM_DIRECT_READ_START];
    const int atom_count = descriptors[ValenceWorkUnitSpecs::ATOM_COUNT];
    for (int i = 0; i < atom_count; i++) {
      xcrd[i] = psynthw.xcrd[atom_offset + i];
      ycrd[i] = psynthw.ycrd[atom_offset + i];
      zcrd[i] = psynthw.zcrd[atom_offset + i];
    }

    // These flags are held in __shared__ memory to conserve registers in an HPC scenario.
    bool overflow_used = true;
    bool minimization_converged = false;

    // Minimization cycles begin here.
    while (minimization_converged = false) {
      
      // Place virtual sites for the initial particle orientation.  Zero all forces.
      // A __syncthreads() must occur after this step in an HPC scenario.
      xfrc[j] = 0;
      yfrc[j] = 0;
      zfrc[j] = 0;
      if (overflow_used) {
        xfrc_overflow[j] = 0;
        yfrc_overflow[j] = 0;
        zfrc_overflow[j] = 0;
      }
      overflow_used = false;
      for (int i = 0; i < descriptors[ValenceWorkUnitSpecs::VSITE_COUNT]; i++) {
        
      }
      
      // Do all valence force field terms, rogue 1:4 interactions, plus restraints.
      const int bond_offset = descriptors[ValenceWorkUnitSpecs::BOND_INSR_START];
      for (int i = 0; i < descriptors[ValenceWorkUnitSpecs::BOND_COUNT]; i++) {

      }
      const int angl_offset = descriptors[ValenceWorkUnitSpecs::ANGL_INSR_START];
      for (int i = 0; i < descriptors[ValenceWorkUnitSpecs::ANGL_COUNT]; i++) {

      }
      const int dihe_offset = descriptors[ValenceWorkUnitSpecs::DIHE_INSR_START];
      for (int i = 0; i < descriptors[ValenceWorkUnitSpecs::DIHE_COUNT]; i++) {

      }
      const int ubrd_offset = descriptors[ValenceWorkUnitSpecs::UBRD_INSR_START];
      for (int i = 0; i < descriptors[ValenceWorkUnitSpecs::UBRD_COUNT]; i++) {

      }
      const int cimp_offset = descriptors[ValenceWorkUnitSpecs::CIMP_INSR_START];
      for (int i = 0; i < descriptors[ValenceWorkUnitSpecs::CIMP_COUNT]; i++) {

      }
      const int cmap_offset = descriptors[ValenceWorkUnitSpecs::CMAP_INSR_START];
      for (int i = 0; i < descriptors[ValenceWorkUnitSpecs::CMAP_COUNT]; i++) {

      }
      const int rg14_offset = descriptors[ValenceWorkUnitSpecs::ROGUE_14_INSR_START];
      for (int i = 0; i < descriptors[ValenceWorkUnitSpecs::ROGUE_14_COUNT]; i++) {

      }
      const int rest_offset = descriptors[ValenceWorkUnitSpecs::RESTRAINT_INSR_START];
      for (int i = 0; i < descriptors[ValenceWorkUnitSpecs::RESTRAINT_COUNT]; i++) {

      }

      // Do all non-bonded tiles.  This will complete the force accumulation and a
      // __syncthreads() will then be required in an HPC scenario.
      for (int i = 0; i < descriptors[ValenceWorkUnitSpecs::NONBONDED_TILE_COUNT]; i++) {
        
      }

      // Compute gradients, using the Generalized Born accumulator arrays for the X, Y, and Z
      // components, respectively.  This will also occur on the GPU to conserve L1 __shared__
      // cache space.
      for (int i = 0; i < descriptors[ValenceWorkUnitSpecs::ATOM_COUNT]; i++) {

      }

      // Move particles along the gradients
    }
  }
}
  
} // namespace mm
} // namespace stormm
