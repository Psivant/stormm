#ifdef STORMM_USE_HPC
#  include <cuda_runtime.h>
#endif
#include "copyright.h"
#include "atomgraph.h"

namespace stormm {
namespace topology {

//-------------------------------------------------------------------------------------------------
void AtomGraph::validateAtomIndex(const int index, const char* caller) {
  if (index < 0 || index >= atom_count) {
    rtErr("Atom index " + std::to_string(index) + " is invalid for a topology of " +
          std::to_string(atom_count) + " atoms.", "AtomGraph", caller);
  }
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void AtomGraph::upload() {
  int_data.upload();
  double_data.upload();
  float_data.upload();
  char4_data.upload();
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::download() {
  int_data.download();
  double_data.download();
  float_data.download();
  char4_data.download();
}
#endif

} // namespace topology
} // namespace stormm
