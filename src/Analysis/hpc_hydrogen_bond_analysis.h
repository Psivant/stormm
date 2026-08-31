// -*-c++-*-
#ifndef STORMM_HPC_HYDROGEN_BOND_ANALYSIS_H
#define STORMM_HPC_HYDROGEN_BOND_ANALYSIS_H

#include "copyright.h"
#include "Synthesis/phasespace_synthesis.h"
#include "hydrogen_bond_analysis.h"

namespace stormm {
namespace analysis {

using synthesis::PsSynthesisReader;
  
/// \brief Launch the kernel to evaluate hydrogen bonds in the current simulation snapshot
///        on the GPU.
///
/// \param step      The current simulation step number, needed to assign data to a particular
///                  statistical block if block averaging is in effect
/// \param hbw       Abstract of the hydrogen bonding analysis, with pointers to memory on the GPU
///                  device
/// \param poly_psr  Abstract of the coordinate synthesis, with pointers to memory on the GPU
void launchEvalHydrogenBondAnalysis(int step, HBondWriter *hbw, const PsSynthesisReader &poly_psr);

} // namespace analysis
} // namespace stormm

#endif 
