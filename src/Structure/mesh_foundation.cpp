#include "copyright.h"
#include "Math/vector_ops.h"
#include "mesh_foundation.h"

namespace stormm {
namespace structure {

using card::HybridKind;
using stmath::accumulateBitmask;

//-------------------------------------------------------------------------------------------------
MeshBasicsKit::MeshBasicsKit(const int* ngbr_in, const size_t* ngbr_bounds_in,
                             const uint* frozen_atoms_in) :
    ngbr{ngbr_in}, ngbr_bounds{ngbr_bounds_in}, frozen_atoms{frozen_atoms_in}
{}

//-------------------------------------------------------------------------------------------------
MeshFoundation::MeshFoundation(const CoordinateFrame *cf_in, const AtomGraph *ag_in,
                               const std::vector<std::string> &comments_in) :
    frames_ready{0},
    cf_pointer{nullptr},
    cf_ensemble{nullptr},
    ag_pointer{nullptr},
    neighbor_list{HybridKind::ARRAY, "mesh_neighbor_list"},
    neighbor_list_bounds{HybridKind::ARRAY, "mesh_nl_bounds"},
    frozen_atoms{HybridKind::ARRAY, "mesh_frozen_atoms"},
    comments{comments_in}
{
  setTopology(ag_in);
  setCoordinates(cf_in);
  testReadiness();
}

//-------------------------------------------------------------------------------------------------
MeshFoundation::MeshFoundation(const CoordinateFrame &cf_in, const AtomGraph &ag_in,
                               const std::vector<std::string> &comments) :
    MeshFoundation(cf_in.getSelfPointer(), ag_in.getSelfPointer(), comments)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshFoundation::MeshFoundation(const CoordinateSeries<T> *cs_in, const AtomGraph *ag_in,
                               const std::vector<std::string> &comments_in) :
    frames_ready{0},
    cf_pointer{nullptr},
    cf_ensemble{nullptr},
    ag_pointer{nullptr},
    neighbor_list{HybridKind::ARRAY, "mesh_neighbor_list"},
    neighbor_list_bounds{HybridKind::ARRAY, "mesh_nl_bounds"},
    frozen_atoms{HybridKind::ARRAY, "mesh_frozen_atoms"},
    comments{comments_in}
{
  setTopology(ag_in);
  setEnsemble(cs_in);
  testReadiness();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
MeshFoundation::MeshFoundation(const CoordinateSeries<T> &cs, const AtomGraph &ag_in,
                               const std::vector<std::string> &comments) :
    MeshFoundation(cs.getSelfPointer(), ag_in.getSelfPointer(), comments)
{}

//-------------------------------------------------------------------------------------------------
int MeshFoundation::getReadyFrameCount() const {
  return frames_ready;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* MeshFoundation::getTopologyPointer() const {
  return ag_pointer;
}

//-------------------------------------------------------------------------------------------------
const CoordinateFrame* MeshFoundation::getCoordinatePointer() const {
  return cf_pointer;
}

//-------------------------------------------------------------------------------------------------
size_t MeshFoundation::getEnsembleTypeCode() const {
  return ensemble_data_type;
}

//-------------------------------------------------------------------------------------------------
std::vector<uint> MeshFoundation::getFrozenAtomMask() const {
  return frozen_atoms.readHost();
}

//-------------------------------------------------------------------------------------------------
int MeshFoundation::getCommentCount() const {
  return static_cast<int>(comments.size());
}

//-------------------------------------------------------------------------------------------------
const std::string& MeshFoundation::getComment(const int comm_index) const {
  if (comm_index < static_cast<int>(comments.size())) {
    return comments[comm_index];
  }
  else {
    rtErr("The mesh has only " + std::to_string(comments.size()) + ".  Comment " +
          std::to_string(comm_index) + " was requested.", "MeshFoundation", "getComment");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
MeshBasicsKit MeshFoundation::data(const HybridTargetLevel tier) const {
  return MeshBasicsKit(neighbor_list.data(tier), neighbor_list_bounds.data(tier),
                       frozen_atoms.data(tier));
}

#ifdef STORMM_USE_HPC
//-------------------------------------------------------------------------------------------------
void MeshFoundation::upload() {
  if (neighbor_list.size() > 0) {
    neighbor_list.upload();
    neighbor_list_bounds.upload();
  }
  frozen_atoms.upload();
}

//-------------------------------------------------------------------------------------------------
void MeshFoundation::download() {
  if (neighbor_list.size() > 0) {
    neighbor_list.download();
    neighbor_list_bounds.download();
  }
  frozen_atoms.download();
}

#  ifdef STORMM_USE_CUDA
//-------------------------------------------------------------------------------------------------
MeshBasicsKit MeshFoundation::deviceViewToHostData(const HybridTargetLevel tier) const {
  return MeshBasicsKit(neighbor_list.getDeviceValidHostPointer(),
                       neighbor_list_bounds.getDeviceValidHostPointer(),
                       frozen_atoms.getDeviceValidHostPointer());
}
#  endif
#endif

//-------------------------------------------------------------------------------------------------
void MeshFoundation::setTopology(const AtomGraph *ag_in) {
  ag_pointer = const_cast<AtomGraph*>(ag_in);
  testReadiness();  

  // Keep a record of the frozen atoms as conveyed by the topology.
  if (ag_pointer != nullptr) {
    const int nbits = static_cast<int>(sizeof(uint)) * 8;
    frozen_atoms.resize((ag_pointer->getAtomCount() + nbits - 1) / nbits);
    std::vector<bool> moving_atoms = ag_pointer->getAtomMobility();
    const int natom = ag_pointer->getAtomCount();
    std::vector<uint> tmp_frozen(frozen_atoms.size(), 0U);
    for (int i = 0; i < natom; i++) {
      if (moving_atoms[i] == false) {
        accumulateBitmask(&tmp_frozen, i);
      }
    }  
    frozen_atoms.putHost(tmp_frozen);
  }
}

//-------------------------------------------------------------------------------------------------
void MeshFoundation::setCoordinates(const CoordinateFrame *cf_in) {
  cf_pointer = const_cast<CoordinateFrame*>(cf_in);
  testReadiness();
}

//-------------------------------------------------------------------------------------------------
void MeshFoundation::setCoordinates(const CoordinateFrame &cf_in) {
  setCoordinates(cf_in.getSelfPointer());
  testReadiness();
}

//-------------------------------------------------------------------------------------------------
void MeshFoundation::computeNeighborLists(const MeshParameters &mps, const MeshRulers &rlrs,
                                          const MeshKlManager &launcher, const PrecisionModel prec,
                                          const HybridTargetLevel availability) {
}

//-------------------------------------------------------------------------------------------------
void MeshFoundation::addComment(const std::string &verbiage) {
  comments.push_back(verbiage);
}

//-------------------------------------------------------------------------------------------------
void MeshFoundation::clearComments() {
  comments.resize(0);
  comments.shrink_to_fit();
}

//-------------------------------------------------------------------------------------------------
int MeshFoundation::validateAtomCounts() const {
  int natom, nframe;
  if (cf_ensemble != nullptr) {
    if (ensemble_data_type == double_type_index) {
      CoordinateSeries<double> *tmp = reinterpret_cast<CoordinateSeries<double>*>(cf_ensemble);
      natom = tmp->getAtomCount();
      nframe = tmp->getFrameCount();
    }
    else if (ensemble_data_type == float_type_index) {
      CoordinateSeries<float> *tmp = reinterpret_cast<CoordinateSeries<float>*>(cf_ensemble);
      natom = tmp->getAtomCount();
      nframe = tmp->getFrameCount();
    }
    else if (ensemble_data_type == llint_type_index) {
      CoordinateSeries<llint> *tmp = reinterpret_cast<CoordinateSeries<llint>*>(cf_ensemble);
      natom = tmp->getAtomCount();
      nframe = tmp->getFrameCount();
    }
    else if (ensemble_data_type == int_type_index) {
      CoordinateSeries<int> *tmp = reinterpret_cast<CoordinateSeries<int>*>(cf_ensemble);
      natom = tmp->getAtomCount();
      nframe = tmp->getFrameCount();
    }
    else if (ensemble_data_type == short_type_index) {
      CoordinateSeries<short> *tmp = reinterpret_cast<CoordinateSeries<short>*>(cf_ensemble);
      natom = tmp->getAtomCount();
      nframe = tmp->getFrameCount();
    }
    else {
      rtErr("A CoordinateSeries object must store its data in one of a selected group of "
            "formats.  It should not be possible to have formed a CoordinateSeries of any other "
            "type in the first place.", "BackgroundMesh", "validateAtomCounts");
    }
  }
  else if (cf_pointer != nullptr) {
    natom = cf_pointer->getAtomCount();
    nframe = 1;
  }
  else {
    natom = 0;
    nframe = 0;
  }
  if (ag_pointer != nullptr) {
    if ((cf_ensemble != nullptr || cf_pointer != nullptr) && natom != ag_pointer->getAtomCount()) {
      rtErr("The number of atoms in the structure coordinates (" + std::to_string(natom) +
            ") must equal that in the topology describing the system (" +
            std::to_string(ag_pointer->getAtomCount()) + ").", "MeshFoundation",
            "validateAtomCounts");
    }
  }
  else {
    nframe = 0;
  }
  return nframe;
}

//-------------------------------------------------------------------------------------------------
void MeshFoundation::testReadiness() {
  if (ag_pointer != nullptr) {
    frames_ready = validateAtomCounts();
  }
  else {
    frames_ready = 0;
  }
}
  
} // namespace structure
} // namespace stormm
