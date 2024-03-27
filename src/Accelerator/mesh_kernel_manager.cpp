#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/parse.h"
#include "Structure/hpc_mesh_support.h"
#include "mesh_kernel_manager.h"

namespace stormm {
namespace card {

using constants::CaseSensitivity;
using parse::strcmpCased;
  
//-------------------------------------------------------------------------------------------------
MeshKlManager::MeshKlManager(const GpuDetails &gpu_in) :
    KernelManager(gpu_in),
    occ_mesh_block_multiplier_dp{meshBlockMultiplier(PrecisionModel::DOUBLE,
                                                     GridDetail::OCCLUSION)},
    occ_mesh_block_multiplier_sp{meshBlockMultiplier(PrecisionModel::SINGLE,
                                                     GridDetail::OCCLUSION)},
    nbf_mesh_block_multiplier_dp{meshBlockMultiplier(PrecisionModel::DOUBLE,
                                                     GridDetail::NONBONDED_FIELD)},
    nbf_mesh_block_multiplier_sp{meshBlockMultiplier(PrecisionModel::SINGLE,
                                                     GridDetail::NONBONDED_FIELD)}
{
  catalogMeshKernel(PrecisionModel::DOUBLE, MappingActivity::PARTICLE_TO_MESH,
                    GridDetail::OCCLUSION, BoundaryCondition::ISOLATED, "kdColorOcclusionMesh");
  catalogMeshKernel(PrecisionModel::SINGLE, MappingActivity::PARTICLE_TO_MESH,
                    GridDetail::OCCLUSION, BoundaryCondition::ISOLATED, "kfColorOcclusionMesh");
  const std::vector<PrecisionModel> build_mdl = { PrecisionModel::DOUBLE, PrecisionModel::SINGLE };
  const std::vector<BoundaryCondition> pbc_mdl = { BoundaryCondition::ISOLATED,
                                                   BoundaryCondition::PERIODIC };
  const std::vector<UnitCellType> uc_mdl = { UnitCellType::ORTHORHOMBIC, UnitCellType::TRICLINIC };
  const std::vector<Interpolant> intp_mdl = { Interpolant::SMOOTHNESS,
                                              Interpolant::FUNCTION_VALUE };
  const std::vector<std::string> build_mdl_str = { "d", "f" };
  const std::vector<std::string> pbc_mdl_str = { "", "Pbc" };
  const std::vector<std::string> uc_mdl_str = { "", "Tric" };
  const std::vector<std::string> intp_mdl_str = { "", "Value" };
  for (size_t i = 0; i < build_mdl.size(); i++) {
    for (size_t j = 0; j < pbc_mdl.size(); j++) {
      for (size_t k = 0; k < uc_mdl.size(); k++) {
        for (size_t m = 0; m < intp_mdl.size(); m++) {
          const std::string knl_name = "k" + build_mdl_str[i] + "ColorNBFieldMesh" +
                                       pbc_mdl_str[j] + uc_mdl_str[k] + intp_mdl_str[m];
          catalogMeshKernel(build_mdl[i], MappingActivity::PARTICLE_TO_MESH,
                            GridDetail::NONBONDED_FIELD, pbc_mdl[j], uc_mdl[k], intp_mdl[m],
                            knl_name);
        }
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void MeshKlManager::catalogMeshKernel(const PrecisionModel prec, const MappingActivity process,
                                      const GridDetail picture, const BoundaryCondition bounds,
                                      const std::string &kernel_name) {
  const std::string k_key = meshKernelKey(prec, picture, bounds, process);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Exclusion grid coloring kernel identifier " + k_key + " already exists in the kernel "
          "map.", "MeshKlManager", "catalogColorExclusionKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryGridKernelRequirements(prec, picture, bounds, process);
  int block_multiplier;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    block_multiplier = occ_mesh_block_multiplier_dp;
    break;
  case PrecisionModel::SINGLE:
    block_multiplier = occ_mesh_block_multiplier_sp;
    break;
  }
  k_dictionary[k_key] = KernelFormat(attr, block_multiplier, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
void MeshKlManager::catalogMeshKernel(const PrecisionModel prec, const MappingActivity process,
                                      const GridDetail picture, const BoundaryCondition bounds,
                                      const UnitCellType unit_cell, const Interpolant stencil_kind,
                                      const std::string &kernel_name) {
  const std::string k_key = meshKernelKey(prec, picture, bounds, unit_cell, stencil_kind, process);
  std::map<std::string, KernelFormat>::iterator it = k_dictionary.find(k_key);
  if (it != k_dictionary.end()) {
    rtErr("Exclusion grid coloring kernel identifier " + k_key + " already exists in the kernel "
          "map.", "MeshKlManager", "catalogColorExclusionKernel");
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  const cudaFuncAttributes attr = queryGridKernelRequirements(prec, picture, bounds, unit_cell,
                                                              stencil_kind, process);
  int block_multiplier;
  switch (prec) {
  case PrecisionModel::DOUBLE:
    block_multiplier = nbf_mesh_block_multiplier_dp;
    break;
  case PrecisionModel::SINGLE:
    block_multiplier = nbf_mesh_block_multiplier_sp;
    break;
  }
  k_dictionary[k_key] = KernelFormat(attr, block_multiplier, 1, gpu, kernel_name);
#  endif
#else
  k_dictionary[k_key] = KernelFormat();
#endif
}

//-------------------------------------------------------------------------------------------------
int2 MeshKlManager::getMeshKernelDims(const PrecisionModel prec, const GridDetail picture,
                                      const BoundaryCondition bounds,
                                      const MappingActivity process) const {
  const std::string k_key = meshKernelKey(prec, picture, bounds, process);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Mesh kernel identifier " + k_key + " was not found in the kernel map.", "MeshKlManager",
          "getMeshKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int2 MeshKlManager::getMeshKernelDims(const PrecisionModel prec, const GridDetail picture,
                                      const BoundaryCondition bounds, const UnitCellType unit_cell,
                                      const Interpolant stencil_kind,
                                      const MappingActivity process) const {
  const std::string k_key = meshKernelKey(prec, picture, bounds, unit_cell, stencil_kind, process);
  if (k_dictionary.find(k_key) == k_dictionary.end()) {
    rtErr("Mesh kernel identifier " + k_key + " was not found in the kernel map.", "MeshKlManager",
          "getMeshKernelDims");
  }
  return k_dictionary.at(k_key).getLaunchParameters();
}

//-------------------------------------------------------------------------------------------------
int meshBlockMultiplier(const PrecisionModel prec, const GridDetail picture) {
#ifdef STORMM_USE_HPC
  switch (prec) {
  case PrecisionModel::DOUBLE:
    switch (picture) {
    case GridDetail::OCCLUSION:
      return 2;
    case GridDetail::OCCLUSION_FIELD:
      break;
    case GridDetail::NONBONDED_FIELD:
      return 1;
    case GridDetail::NONBONDED_ATOMIC:
      break;
    }
    return 1;
  case PrecisionModel::SINGLE:
    switch (picture) {
    case GridDetail::OCCLUSION:
      return 4;
    case GridDetail::OCCLUSION_FIELD:
      break;
    case GridDetail::NONBONDED_FIELD:
      return 1;
    case GridDetail::NONBONDED_ATOMIC:
      break;
    }
    return 1;
  }
  __builtin_unreachable();
#else
  return 1;
#endif
}

//-------------------------------------------------------------------------------------------------
std::string meshKernelKey(const PrecisionModel prec, const GridDetail picture,
                          const BoundaryCondition bounds, const MappingActivity process) {
  std::string k_key("mesh_");
  switch (prec) {
  case PrecisionModel::DOUBLE:
    k_key += "d";
    break;
  case PrecisionModel::SINGLE:
    k_key += "f";
    break;
  }
  switch (picture) {
  case GridDetail::OCCLUSION:
    k_key += "_oc";
    break;
  case GridDetail::OCCLUSION_FIELD:
    k_key += "_of";
    break;
  case GridDetail::NONBONDED_FIELD:
    k_key += "_nf";
    break;
  case GridDetail::NONBONDED_ATOMIC:
    k_key += "_na";
    break;
  }
  switch (bounds) {
  case BoundaryCondition::ISOLATED:
    k_key += "_iso";
    break;
  case BoundaryCondition::PERIODIC:
    k_key += "_pbc";
    break;
  }
  switch (process) {
  case MappingActivity::PARTICLE_TO_MESH:
    k_key += "_pm";
    break;
  case MappingActivity::MESH_TO_PARTICLE:
    k_key += "_mp";
    break;
  }
  return k_key;
}

//-------------------------------------------------------------------------------------------------
std::string meshKernelKey(const PrecisionModel prec, const GridDetail picture,
                          const BoundaryCondition bounds, const UnitCellType unit_cell,
                          const Interpolant stencil_kind, const MappingActivity process) {
  std::string k_key = meshKernelKey(prec, picture, bounds, process);
  k_key = k_key.substr(0, k_key.size() - 3);
  switch (unit_cell) {
  case UnitCellType::ORTHORHOMBIC:
    k_key += "_rc";
    break;
  case UnitCellType::TRICLINIC:
    k_key += "_tr";
    break;
  case UnitCellType::NONE:
    break;
  }
  switch (stencil_kind) {
  case Interpolant::SMOOTHNESS:
    k_key += "_sm";
    break;
  case Interpolant::FUNCTION_VALUE:
    k_key += "_va";
    break;
  }
  switch (process) {
  case MappingActivity::PARTICLE_TO_MESH:
    k_key += "_pm";
    break;
  case MappingActivity::MESH_TO_PARTICLE:
    k_key += "_mp";
    break;
  }
  return k_key;
}

} // namespace card
} // namespace stormm
