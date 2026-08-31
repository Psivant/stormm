// -*-c++-*-
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Accelerator/ptx_macros.h"
#include "Constants/symbol_values.h"
#include "Potential/energy_enumerators.h"
#include "Potential/scorecard.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Synthesis/valence_workunit.h"
#include "Topology/atomgraph_enumerators.h"
#include "hpc_integration.h"

namespace stormm {
namespace trajectory {

using card::HybridTargetLevel;
using energy::ScoreCardWriter;
using energy::StateVariable;
using energy::ValenceKernelSize;
using numerics::AccumulationMethod;
using numerics::getEnumerationName;
using symbols::kcal_to_gafs_f;
using symbols::gafs_to_kcal_f;
using symbols::boltzmann_constant_f;
using synthesis::AtomGraphSynthesis;
using synthesis::eighth_valence_work_unit_atoms;
using synthesis::half_valence_work_unit_atoms;
using synthesis::maximum_valence_work_unit_atoms;
using synthesis::quarter_valence_work_unit_atoms;
using synthesis::vwu_abstract_length;
using synthesis::VwuAbstractMap;
using topology::VirtualSiteKind;
  
#include "Accelerator/syncwarp.cui"
#include "Math/rounding.cui"
#include "Math/vector_formulas.cui"
#include "Numerics/accumulation.cui"

#define VERLET_STANDALONE

// Double-precision floating point arithmetic
#define TCALC double
#  define SPLIT_FORCE_ACCUMULATION
#  define TCALC2 double2
#  define TCALC3 double3
#  define TCALC4 double4_16a
#  define SQRT_FUNC sqrt
#  define COS_FUNC cos
#  define SIN_FUNC sin
#  define ABS_FUNC fabs
#  define LLCONV_FUNC __double2ll_rn
#    define INTEG_KERNEL_THREAD_COUNT 320
#    define INTEG_BLOCK_MULTIPLICITY 2
#      define KERNEL_NAME kdsIntegVelAdv
#        include "verlet_i.cui"
#      undef KERNEL_NAME
#      define PME_COMPATIBLE
#        define TCOORD double
#        define TACC llint
#        define TCOORD4 double4_16a
#        define TCOORD_IS_LONG
#          define KERNEL_NAME kdsdPmeIntegVelAdv
#            include "verlet_i.cui"
#          undef KERNEL_NAME
#          define DUAL_GRIDS
#            define KERNEL_NAME kdsdPmeDualIntegVelAdv
#              include "verlet_i.cui"
#            undef KERNEL_NAME
#          undef DUAL_GRIDS
#        undef TCOORD
#        undef TACC
#        undef TCOORD4
#        undef TCOORD_IS_LONG
#        define TCOORD float
#        define TACC int
#        define TCOORD4 float4
#          define KERNEL_NAME kdsfPmeIntegVelAdv
#            include "verlet_i.cui"
#          undef KERNEL_NAME
#          define DUAL_GRIDS
#            define KERNEL_NAME kdsfPmeDualIntegVelAdv
#              include "verlet_i.cui"
#            undef KERNEL_NAME
#          undef DUAL_GRIDS
#        undef TCOORD
#        undef TACC
#        undef TCOORD4
#        define KERNEL_NAME kdsPmeIntegVelCnst
#          include "velocity_constraints.cui"
#        undef KERNEL_NAME
#        define KERNEL_NAME kdsPmeIntegGeomCnst
#          include "geometry_constraints.cui"
#        undef KERNEL_NAME
#      undef PME_COMPATIBLE
#      define KERNEL_NAME kdsIntegVelCnst
#        include "velocity_constraints.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kdsKineticEnergy
#        include "MolecularMechanics/kinetic.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kdsIntegPosAdv
#        include "verlet_ii.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kdsIntegGeomCnst
#        include "geometry_constraints.cui"
#      undef KERNEL_NAME
#    undef INTEG_BLOCK_MULTIPLICITY
#    undef INTEG_KERNEL_THREAD_COUNT
#  undef SPLIT_FORCE_ACCUMULATION
#  undef TCALC2
#  undef TCALC3
#  undef TCALC4
#  undef SQRT_FUNC
#  undef COS_FUNC
#  undef SIN_FUNC
#  undef ABS_FUNC
#  undef LLCONV_FUNC
#undef TCALC

// Single-precision floating point arithmetic, split accumulation followed by fused accumulators
#define TCALC float
#  define TCALC_IS_SINGLE
#  define TCALC2 float2
#  define TCALC3 float3
#  define TCALC4 float4
#  define SQRT_FUNC sqrtf
#  define COS_FUNC cosf
#  define SIN_FUNC sinf
#  define ABS_FUNC fabsf
#  define LLCONV_FUNC __float2ll_rn
#  define SPLIT_FORCE_ACCUMULATION
#    define INTEG_KERNEL_THREAD_COUNT 512
#    define INTEG_BLOCK_MULTIPLICITY 2
#      define KERNEL_NAME kfsIntegVelAdvXL
#        include "verlet_i.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfsIntegVelCnstXL
#        include "velocity_constraints.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfsKineticEnergyXL
#        include "MolecularMechanics/kinetic.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfsIntegPosAdvXL
#        include "verlet_ii.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfsIntegGeomCnstXL
#        include "geometry_constraints.cui"
#      undef KERNEL_NAME
#      define PME_COMPATIBLE
#        define TCOORD double
#        define TACC llint
#        define TCOORD4 double4_16a
#        define TCOORD_IS_LONG
#          define KERNEL_NAME kfsdPmeIntegVelAdv
#            include "verlet_i.cui"
#          undef KERNEL_NAME
#          define DUAL_GRIDS
#            define KERNEL_NAME kfsdPmeDualIntegVelAdv
#              include "verlet_i.cui"
#            undef KERNEL_NAME
#          undef DUAL_GRIDS
#        undef TCOORD
#        undef TACC
#        undef TCOORD4
#        undef TCOORD_IS_LONG
#        define TCOORD float
#        define TACC int
#        define TCOORD4 float4
#          define KERNEL_NAME kfsfPmeIntegVelAdv
#            include "verlet_i.cui"
#          undef KERNEL_NAME
#          define DUAL_GRIDS
#            define KERNEL_NAME kfsfPmeDualIntegVelAdv
#              include "verlet_i.cui"
#            undef KERNEL_NAME
#          undef DUAL_GRIDS
#        undef TCOORD
#        undef TACC
#        undef TCOORD4
#        define KERNEL_NAME kfsPmeIntegVelCnst
#          include "velocity_constraints.cui"
#        undef KERNEL_NAME
#        define KERNEL_NAME kfsPmeIntegGeomCnst
#          include "geometry_constraints.cui"
#        undef KERNEL_NAME
#      undef PME_COMPATIBLE
#    undef INTEG_BLOCK_MULTIPLICITY
#    undef INTEG_KERNEL_THREAD_COUNT
#    define INTEG_KERNEL_THREAD_COUNT 256
#    define INTEG_BLOCK_MULTIPLICITY 4
#      define KERNEL_NAME kfsIntegVelAdvLG
#        include "verlet_i.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfsIntegVelCnstLG
#        include "velocity_constraints.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfsKineticEnergyLG
#        include "MolecularMechanics/kinetic.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfsIntegPosAdvLG
#        include "verlet_ii.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfsIntegGeomCnstLG
#        include "geometry_constraints.cui"
#      undef KERNEL_NAME
#    undef INTEG_BLOCK_MULTIPLICITY
#    undef INTEG_KERNEL_THREAD_COUNT
#    define INTEG_KERNEL_THREAD_COUNT 128
#    define INTEG_BLOCK_MULTIPLICITY 8
#      define KERNEL_NAME kfsIntegVelAdvMD
#        include "verlet_i.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfsIntegVelCnstMD
#        include "velocity_constraints.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfsKineticEnergyMD
#        include "MolecularMechanics/kinetic.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfsIntegPosAdvMD
#        include "verlet_ii.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfsIntegGeomCnstMD
#        include "geometry_constraints.cui"
#      undef KERNEL_NAME
#    undef INTEG_BLOCK_MULTIPLICITY
#    undef INTEG_KERNEL_THREAD_COUNT
#    define INTEG_KERNEL_THREAD_COUNT 64
#    define INTEG_BLOCK_MULTIPLICITY 16
#      define KERNEL_NAME kfsIntegVelAdvSM
#        include "verlet_i.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfsIntegVelCnstSM
#        include "velocity_constraints.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfsKineticEnergySM
#        include "MolecularMechanics/kinetic.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfsIntegPosAdvSM
#        include "verlet_ii.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfsIntegGeomCnstSM
#        include "geometry_constraints.cui"
#      undef KERNEL_NAME
#    undef INTEG_BLOCK_MULTIPLICITY
#    undef INTEG_KERNEL_THREAD_COUNT
#  undef SPLIT_FORCE_ACCUMULATION
#  define INTEG_KERNEL_THREAD_COUNT 512
#  define INTEG_BLOCK_MULTIPLICITY 2
#    define KERNEL_NAME kfIntegVelAdvXL
#      include "verlet_i.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfIntegVelCnstXL
#      include "velocity_constraints.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfKineticEnergyXL
#      include "MolecularMechanics/kinetic.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfIntegPosAdvXL
#      include "verlet_ii.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfIntegGeomCnstXL
#      include "geometry_constraints.cui"
#    undef KERNEL_NAME
#    define PME_COMPATIBLE
#      define TCOORD double
#      define TACC llint
#      define TCOORD4 double4_16a
#      define TCOORD_IS_LONG
#        define KERNEL_NAME kfdPmeIntegVelAdv
#          include "verlet_i.cui"
#        undef KERNEL_NAME
#        define DUAL_GRIDS
#          define KERNEL_NAME kfdPmeDualIntegVelAdv
#            include "verlet_i.cui"
#          undef KERNEL_NAME
#        undef DUAL_GRIDS
#      undef TCOORD
#      undef TACC
#      undef TCOORD4
#      undef TCOORD_IS_LONG
#      define TCOORD float
#      define TACC int
#      define TCOORD4 float4
#        define KERNEL_NAME kffPmeIntegVelAdv
#          include "verlet_i.cui"
#        undef KERNEL_NAME
#        define DUAL_GRIDS
#          define KERNEL_NAME kffPmeDualIntegVelAdv
#            include "verlet_i.cui"
#          undef KERNEL_NAME
#        undef DUAL_GRIDS
#      undef TCOORD
#      undef TACC
#      undef TCOORD4
#      define KERNEL_NAME kfPmeIntegVelCnst
#        include "velocity_constraints.cui"
#      undef KERNEL_NAME
#      define KERNEL_NAME kfPmeIntegGeomCnst
#        include "geometry_constraints.cui"
#      undef KERNEL_NAME
#    undef PME_COMPATIBLE
#  undef INTEG_BLOCK_MULTIPLICITY
#  undef INTEG_KERNEL_THREAD_COUNT
#  define INTEG_KERNEL_THREAD_COUNT 256
#  define INTEG_BLOCK_MULTIPLICITY 4
#    define KERNEL_NAME kfIntegVelAdvLG
#      include "verlet_i.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfIntegVelCnstLG
#      include "velocity_constraints.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfKineticEnergyLG
#      include "MolecularMechanics/kinetic.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfIntegPosAdvLG
#      include "verlet_ii.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfIntegGeomCnstLG
#      include "geometry_constraints.cui"
#    undef KERNEL_NAME
#  undef INTEG_BLOCK_MULTIPLICITY
#  undef INTEG_KERNEL_THREAD_COUNT
#  define INTEG_KERNEL_THREAD_COUNT 128
#  define INTEG_BLOCK_MULTIPLICITY 8
#    define KERNEL_NAME kfIntegVelAdvMD
#      include "verlet_i.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfIntegVelCnstMD
#      include "velocity_constraints.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfKineticEnergyMD
#      include "MolecularMechanics/kinetic.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfIntegPosAdvMD
#      include "verlet_ii.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfIntegGeomCnstMD
#      include "geometry_constraints.cui"
#    undef KERNEL_NAME
#  undef INTEG_BLOCK_MULTIPLICITY
#  undef INTEG_KERNEL_THREAD_COUNT
#  define INTEG_KERNEL_THREAD_COUNT 64
#  define INTEG_BLOCK_MULTIPLICITY 16
#    define KERNEL_NAME kfIntegVelAdvSM
#      include "verlet_i.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfIntegVelCnstSM
#      include "velocity_constraints.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfKineticEnergySM
#      include "MolecularMechanics/kinetic.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfIntegPosAdvSM
#      include "verlet_ii.cui"
#    undef KERNEL_NAME
#    define KERNEL_NAME kfIntegGeomCnstSM
#      include "geometry_constraints.cui"
#    undef KERNEL_NAME
#  undef INTEG_BLOCK_MULTIPLICITY
#  undef INTEG_KERNEL_THREAD_COUNT
#  undef TCALC_IS_SINGLE
#  undef TCALC2
#  undef TCALC3
#  undef TCALC4
#  undef SQRT_FUNC
#  undef COS_FUNC
#  undef SIN_FUNC
#  undef ABS_FUNC
#  undef LLCONV_FUNC
#undef TCALC

#undef VERLET_STANDALONE
  
//-------------------------------------------------------------------------------------------------
cudaFuncAttributes queryIntegrationKernelRequirements(const PrecisionModel calc_prec,
                                                      const PrecisionModel neighbor_prec,
                                                      const UnitCellType unit_cell,
                                                      const AccumulationMethod acc_meth,
                                                      const ValenceKernelSize kwidth,
                                                      const IntegrationStage process) {
  cudaFuncAttributes result;
  cudaError_t cfa;

  switch (unit_cell) {
  case UnitCellType::NONE:
    switch (calc_prec) {
    case PrecisionModel::DOUBLE:
      switch (process) {
      case IntegrationStage::CALC_FORCES:
        break;
      case IntegrationStage::VELOCITY_ADVANCE:
        cfa = cudaFuncGetAttributes(&result, kdsIntegVelAdv);
        break;
      case IntegrationStage::VELOCITY_CONSTRAINT:
        cfa = cudaFuncGetAttributes(&result, kdsIntegVelCnst);
        break;
      case IntegrationStage::CALC_KINETIC:
        cfa = cudaFuncGetAttributes(&result, kdsKineticEnergy);
        break;
      case IntegrationStage::POSITION_ADVANCE:
        cfa = cudaFuncGetAttributes(&result, kdsIntegPosAdv);
        break;
      case IntegrationStage::GEOMETRY_CONSTRAINT:
        cfa = cudaFuncGetAttributes(&result, kdsIntegGeomCnst);
        break;
      }
      break;
    case PrecisionModel::SINGLE:
      switch (acc_meth) {
      case AccumulationMethod::SPLIT:
        switch (kwidth) {
        case ValenceKernelSize::XL:
          switch (process) {
          case IntegrationStage::CALC_FORCES:
            break;
          case IntegrationStage::VELOCITY_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfsIntegVelAdvXL);
            break;
          case IntegrationStage::VELOCITY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfsIntegVelCnstXL);
            break;
          case IntegrationStage::CALC_KINETIC:
            cfa = cudaFuncGetAttributes(&result, kfsKineticEnergyXL);
            break;
          case IntegrationStage::POSITION_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfsIntegPosAdvXL);
            break;
          case IntegrationStage::GEOMETRY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfsIntegGeomCnstXL);
            break;
          }
          break;
        case ValenceKernelSize::LG:
          switch (process) {
          case IntegrationStage::CALC_FORCES:
            break;
          case IntegrationStage::VELOCITY_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfsIntegVelAdvLG);
            break;
          case IntegrationStage::VELOCITY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfsIntegVelCnstLG);
            break;
          case IntegrationStage::CALC_KINETIC:
            cfa = cudaFuncGetAttributes(&result, kfsKineticEnergyLG);
            break;
          case IntegrationStage::POSITION_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfsIntegPosAdvLG);
            break;
          case IntegrationStage::GEOMETRY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfsIntegGeomCnstLG);
            break;
          }
          break;
        case ValenceKernelSize::MD:
          switch (process) {
          case IntegrationStage::CALC_FORCES:
            break;
          case IntegrationStage::VELOCITY_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfsIntegVelAdvMD);
            break;
          case IntegrationStage::VELOCITY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfsIntegVelCnstMD);
            break;
          case IntegrationStage::CALC_KINETIC:
            cfa = cudaFuncGetAttributes(&result, kfsKineticEnergyMD);
            break;
          case IntegrationStage::POSITION_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfsIntegPosAdvMD);
            break;
          case IntegrationStage::GEOMETRY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfsIntegGeomCnstMD);
            break;
          }
          break;
        case ValenceKernelSize::SM:
          switch (process) {
          case IntegrationStage::CALC_FORCES:
            break;
          case IntegrationStage::VELOCITY_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfsIntegVelAdvSM);
            break;
          case IntegrationStage::VELOCITY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfsIntegVelCnstSM);
            break;
          case IntegrationStage::CALC_KINETIC:
            cfa = cudaFuncGetAttributes(&result, kfsKineticEnergySM);
            break;
          case IntegrationStage::POSITION_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfsIntegPosAdvSM);
            break;
          case IntegrationStage::GEOMETRY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfsIntegGeomCnstSM);
            break;
          }
          break;
        }
        break;
      case AccumulationMethod::WHOLE:
        switch (kwidth) {
        case ValenceKernelSize::XL:
          switch (process) {
          case IntegrationStage::CALC_FORCES:
            break;
          case IntegrationStage::VELOCITY_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfIntegVelAdvXL);
            break;
          case IntegrationStage::VELOCITY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfIntegVelCnstXL);
            break;
          case IntegrationStage::CALC_KINETIC:
            cfa = cudaFuncGetAttributes(&result, kfKineticEnergyXL);
            break;
          case IntegrationStage::POSITION_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfIntegPosAdvXL);
            break;
          case IntegrationStage::GEOMETRY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfIntegGeomCnstXL);
            break;
          }
          break;
        case ValenceKernelSize::LG:
          switch (process) {
          case IntegrationStage::CALC_FORCES:
            break;
          case IntegrationStage::VELOCITY_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfIntegVelAdvLG);
            break;
          case IntegrationStage::VELOCITY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfIntegVelCnstLG);
            break;
          case IntegrationStage::CALC_KINETIC:
            cfa = cudaFuncGetAttributes(&result, kfKineticEnergyLG);
            break;
          case IntegrationStage::POSITION_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfIntegPosAdvLG);
            break;
          case IntegrationStage::GEOMETRY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfIntegGeomCnstLG);
            break;
          }
          break;
        case ValenceKernelSize::MD:
          switch (process) {
          case IntegrationStage::CALC_FORCES:
            break;
          case IntegrationStage::VELOCITY_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfIntegVelAdvMD);
            break;
          case IntegrationStage::VELOCITY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfIntegVelCnstMD);
            break;
          case IntegrationStage::CALC_KINETIC:
            cfa = cudaFuncGetAttributes(&result, kfKineticEnergyMD);
            break;
          case IntegrationStage::POSITION_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfIntegPosAdvMD);
            break;
          case IntegrationStage::GEOMETRY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfIntegGeomCnstMD);
            break;
          }
          break;
        case ValenceKernelSize::SM:
          switch (process) {
          case IntegrationStage::CALC_FORCES:
            break;
          case IntegrationStage::VELOCITY_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfIntegVelAdvSM);
            break;
          case IntegrationStage::VELOCITY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfIntegVelCnstSM);
            break;
          case IntegrationStage::CALC_KINETIC:
            cfa = cudaFuncGetAttributes(&result, kfKineticEnergySM);
            break;
          case IntegrationStage::POSITION_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfIntegPosAdvSM);
            break;
          case IntegrationStage::GEOMETRY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfIntegGeomCnstSM);
            break;
          }
          break;
        }
        break;
      case AccumulationMethod::AUTOMATIC:
        break;
      }
      break;
    }
    break;
  case UnitCellType::ORTHORHOMBIC:
  case UnitCellType::TRICLINIC:
    switch (calc_prec) {
    case PrecisionModel::DOUBLE:
      switch (neighbor_prec) {
      case PrecisionModel::DOUBLE:
        switch (process) {
        case IntegrationStage::CALC_FORCES:
          break;
        case IntegrationStage::VELOCITY_ADVANCE:
          cfa = cudaFuncGetAttributes(&result, kdsdPmeIntegVelAdv);
          break;
        case IntegrationStage::VELOCITY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kdsPmeIntegVelCnst);
          break;
        case IntegrationStage::CALC_KINETIC:
          cfa = cudaFuncGetAttributes(&result, kdsKineticEnergy);
          break;
        case IntegrationStage::POSITION_ADVANCE:
          cfa = cudaFuncGetAttributes(&result, kdsIntegPosAdv);
          break;
        case IntegrationStage::GEOMETRY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kdsPmeIntegGeomCnst);
          break;
        }
        break;
      case PrecisionModel::SINGLE:
        switch (process) {
        case IntegrationStage::CALC_FORCES:
          break;
        case IntegrationStage::VELOCITY_ADVANCE:
          cfa = cudaFuncGetAttributes(&result, kdsfPmeIntegVelAdv);
          break;
        case IntegrationStage::VELOCITY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kdsPmeIntegVelCnst);
          break;
        case IntegrationStage::CALC_KINETIC:
          cfa = cudaFuncGetAttributes(&result, kdsKineticEnergy);
          break;
        case IntegrationStage::POSITION_ADVANCE:
          cfa = cudaFuncGetAttributes(&result, kdsIntegPosAdv);
          break;
        case IntegrationStage::GEOMETRY_CONSTRAINT:
          cfa = cudaFuncGetAttributes(&result, kdsPmeIntegGeomCnst);
          break;
        }
        break;
      }
      break;
    case PrecisionModel::SINGLE:
      switch (neighbor_prec) {
      case PrecisionModel::DOUBLE:
        switch (acc_meth) {
        case AccumulationMethod::SPLIT:
          switch (process) {
          case IntegrationStage::CALC_FORCES:
            break;
          case IntegrationStage::VELOCITY_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfsdPmeIntegVelAdv);
            break;
          case IntegrationStage::VELOCITY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfsPmeIntegVelCnst);
            break;
          case IntegrationStage::CALC_KINETIC:
            cfa = cudaFuncGetAttributes(&result, kfsKineticEnergyXL);
            break;
          case IntegrationStage::POSITION_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfsIntegPosAdvXL);
            break;
          case IntegrationStage::GEOMETRY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfsPmeIntegGeomCnst);
            break;
          }
          break;
        case AccumulationMethod::WHOLE:
          switch (process) {
          case IntegrationStage::CALC_FORCES:
            break;
          case IntegrationStage::VELOCITY_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfdPmeIntegVelAdv);
            break;
          case IntegrationStage::VELOCITY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfPmeIntegVelCnst);
            break;
          case IntegrationStage::CALC_KINETIC:
            cfa = cudaFuncGetAttributes(&result, kfsKineticEnergyXL);
            break;
          case IntegrationStage::POSITION_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfIntegPosAdvXL);
            break;
          case IntegrationStage::GEOMETRY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfPmeIntegGeomCnst);
            break;
          }
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      case PrecisionModel::SINGLE:
        switch (acc_meth) {
        case AccumulationMethod::SPLIT:
          switch (process) {
          case IntegrationStage::CALC_FORCES:
            break;
          case IntegrationStage::VELOCITY_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfsfPmeIntegVelAdv);
            break;
          case IntegrationStage::VELOCITY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfsPmeIntegVelCnst);
            break;
          case IntegrationStage::CALC_KINETIC:
            cfa = cudaFuncGetAttributes(&result, kfsKineticEnergyXL);
            break;
          case IntegrationStage::POSITION_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfsIntegPosAdvXL);
            break;
          case IntegrationStage::GEOMETRY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfsPmeIntegGeomCnst);
            break;
          }
          break;
        case AccumulationMethod::WHOLE:
          switch (process) {
          case IntegrationStage::CALC_FORCES:
            break;
          case IntegrationStage::VELOCITY_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kffPmeIntegVelAdv);
            break;
          case IntegrationStage::VELOCITY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfPmeIntegVelCnst);
            break;
          case IntegrationStage::CALC_KINETIC:
            cfa = cudaFuncGetAttributes(&result, kfsKineticEnergyXL);
            break;
          case IntegrationStage::POSITION_ADVANCE:
            cfa = cudaFuncGetAttributes(&result, kfIntegPosAdvXL);
            break;
          case IntegrationStage::GEOMETRY_CONSTRAINT:
            cfa = cudaFuncGetAttributes(&result, kfPmeIntegGeomCnst);
            break;
          }
          break;
        case AccumulationMethod::AUTOMATIC:
          break;
        }
        break;
      }
      break;
    }
    break;
  }

  // Check for errors
  if (cfa != cudaSuccess) {

    // Construct the appropriate error message
    std::string error_message("Error obtaining attributes for kernel k");
    switch (calc_prec) {
    case PrecisionModel::DOUBLE:
      error_message += "ds";
      switch (acc_meth) {
      case AccumulationMethod::SPLIT:
        break;
      case AccumulationMethod::WHOLE:
      case AccumulationMethod::AUTOMATIC:
        rtErr("No kernel is available for " + getEnumerationName(calc_prec) + "-precision " +
              getEnumerationName(acc_meth) + " force accumulation.",
              "queryIntegrationKernelRequirementts");
        break;
      }
      break;
    case PrecisionModel::SINGLE:
      error_message += "f";
      switch (acc_meth) {
      case AccumulationMethod::SPLIT:
        error_message += "s";
        break;
      case AccumulationMethod::WHOLE:
      case AccumulationMethod::AUTOMATIC:
        break;
      }
      break;
    }
    switch (unit_cell) {
    case UnitCellType::NONE:
      break;
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      switch (process) {
      case IntegrationStage::CALC_FORCES:
        break;
      case IntegrationStage::VELOCITY_ADVANCE:
      case IntegrationStage::VELOCITY_CONSTRAINT:
      case IntegrationStage::GEOMETRY_CONSTRAINT:
        switch (neighbor_prec) {
        case PrecisionModel::DOUBLE:
          error_message += "d";
          break;
        case PrecisionModel::SINGLE:
          error_message += "f";
          break;
        }
        error_message += "Pme";
        break;
      case IntegrationStage::CALC_KINETIC:
      case IntegrationStage::POSITION_ADVANCE:
        break;
      }
      break;
    }
    error_message += "Integ";
    switch (process) {
    case IntegrationStage::CALC_FORCES:
      break;
    case IntegrationStage::VELOCITY_ADVANCE:
      error_message += "VelAdv";
      break;
    case IntegrationStage::VELOCITY_CONSTRAINT:
      error_message += "VelCnst";
      break;
    case IntegrationStage::CALC_KINETIC:
      error_message += "KineticEnergy";
      break;
    case IntegrationStage::POSITION_ADVANCE:
      error_message += "PosAdv";
      break;
    case IntegrationStage::GEOMETRY_CONSTRAINT:
      error_message += "GeomCnst";
      break;
    }
    switch (unit_cell) {
    case UnitCellType::NONE:
      switch (kwidth) {
      case ValenceKernelSize::XL:
        error_message += "XL";
        break;
      case ValenceKernelSize::LG:
        error_message += "LG";
        break;
      case ValenceKernelSize::MD:
        error_message += "MD";
        break;
      case ValenceKernelSize::SM:
        error_message += "SM";
        break;
      }
      break;
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      switch (process) {
      case IntegrationStage::CALC_FORCES:
        break;
      case IntegrationStage::VELOCITY_ADVANCE:
      case IntegrationStage::VELOCITY_CONSTRAINT:
      case IntegrationStage::CALC_KINETIC:
      case IntegrationStage::GEOMETRY_CONSTRAINT:
        break;
      case IntegrationStage::POSITION_ADVANCE:
        error_message += "XL";
        break;
      }
      break;
    }
    error_message += ".";

    // Report the error
    rtErr(error_message, "queryIntegrationKernelRequirements");
  }  
  return result;
}

//-------------------------------------------------------------------------------------------------
cudaFuncAttributes queryIntegrationKernelRequirements(const PrecisionModel calc_prec,
                                                      const AccumulationMethod acc_meth,
                                                      const ValenceKernelSize kwidth,
                                                      const IntegrationStage process) {
  return queryIntegrationKernelRequirements(calc_prec, PrecisionModel::SINGLE, UnitCellType::NONE,
                                            acc_meth, kwidth, process);
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<double> *tb_resw,
                              MMControlKit<double> *ctrl, ScoreCardWriter *scw,
                              const SyValenceKit<double> &poly_vk,
                              const SyAtomUpdateKit<double, double2, double4_16a> &poly_auk,
                              const ThermostatWriter<double> &tstw, const int2 lp,
                              const IntegrationStage process) {
  switch (process) {
  case IntegrationStage::CALC_FORCES:

    // Force calculation is handled by multiple kernels.  This enumeration exists for completeness,
    // so that such a stage of the time step may be identified when intervening with customized
    // function pointers.
    break;
  case IntegrationStage::VELOCITY_ADVANCE:

    // If the unit cell type indicates periodic boundary conditions, a neighbor list object is
    // required ro run the proper velocity update.
    switch (poly_psw->unit_cell) {
    case UnitCellType::NONE:
      kdsIntegVelAdv<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      rtErr("A neighbor list abstract (CellGridReader) is required to launch the proper Verlet "
            "velocity update kernel.", "launchIntegrationProcess");
    }
    break;
  case IntegrationStage::VELOCITY_CONSTRAINT:

    // Systems with periodic boundary conditions will call versions of each constraint kernel that
    // include SETTLE routines for rigid water molecules.
    switch (poly_psw->unit_cell) {
    case UnitCellType::NONE:
      kdsIntegVelCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      kdsPmeIntegVelCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    }
    break;
  case IntegrationStage::CALC_KINETIC:
    if (scw == nullptr) {
      rtErr("This stage of the integration process must be called with a ScoreCard abstract for "
            "energy tracking.  Use a different overloaded function variant.",
            "launchIntegrationProcess.");
    }
    kdsKineticEnergy<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
    break;
  case IntegrationStage::POSITION_ADVANCE:

    // Systems with either periodic or isolated boundary conditions call the same Verlet position
    // advancement kernel.
    kdsIntegPosAdv<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
    break;
  case IntegrationStage::GEOMETRY_CONSTRAINT:
    switch (poly_psw->unit_cell) {
    case UnitCellType::NONE:
      kdsIntegGeomCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      kdsPmeIntegGeomCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<double> *tb_resw,
                              MMControlKit<double> *ctrl, const SyValenceKit<double> &poly_vk,
                              const SyAtomUpdateKit<double, double2, double4_16a> &poly_auk,
                              const ThermostatWriter<double> &tstw, const int2 lp,
                              const IntegrationStage process) {
  launchIntegrationProcess(poly_psw, tb_resw, ctrl, nullptr, poly_vk, poly_auk, tstw, lp, process);
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<float> *tb_resw,
                              MMControlKit<float> *ctrl, ScoreCardWriter *scw,
                              const SyValenceKit<float> &poly_vk,
                              const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                              const ThermostatWriter<float> &tstw, const int2 lp,
                              const AccumulationMethod acc_meth, const ValenceKernelSize kwidth,
                              const IntegrationStage process) {
  switch (process) {
  case IntegrationStage::CALC_FORCES:
    break;
  case IntegrationStage::VELOCITY_ADVANCE:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      switch (poly_psw->unit_cell) {
      case UnitCellType::NONE:
        switch (kwidth) {
        case ValenceKernelSize::XL:
          kfsIntegVelAdvXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::LG:
          kfsIntegVelAdvLG<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::MD:
          kfsIntegVelAdvMD<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::SM:
          kfsIntegVelAdvSM<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        }
        break;
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        rtErr("A neighbor list abstract (CellGridReader) is required to launch the proper Verlet "
              "velocity update kernel.", "launchIntegrationProcess");
      }
      break;
    case AccumulationMethod::WHOLE:
      switch (poly_psw->unit_cell) {
      case UnitCellType::NONE:
        switch (kwidth) {
        case ValenceKernelSize::XL:
          kfIntegVelAdvXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::LG:
          kfIntegVelAdvLG<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::MD:
          kfIntegVelAdvMD<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::SM:
          kfIntegVelAdvSM<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        }
        break;
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        rtErr("A neighbor list abstract (CellGridReader) is required to launch the proper Verlet "
              "velocity update kernel.", "launchIntegrationProcess");
      }
      break;
    }
    break;
  case IntegrationStage::VELOCITY_CONSTRAINT:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      switch (poly_psw->unit_cell) {
      case UnitCellType::NONE:
        switch (kwidth) {
        case ValenceKernelSize::XL:
          kfsIntegVelCnstXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::LG:
          kfsIntegVelCnstLG<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::MD:
          kfsIntegVelCnstMD<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::SM:
          kfsIntegVelCnstSM<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        }
        break;
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        kfsPmeIntegVelCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
        break;
      }
      break;
    case AccumulationMethod::WHOLE:
      switch (poly_psw->unit_cell) {
      case UnitCellType::NONE:
        switch (kwidth) {
        case ValenceKernelSize::XL:
          kfIntegVelCnstXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::LG:
          kfIntegVelCnstLG<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::MD:
          kfIntegVelCnstMD<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::SM:
          kfIntegVelCnstSM<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        }
        break;
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        kfPmeIntegVelCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
        break;
      }
      break;
    }
    break;
  case IntegrationStage::CALC_KINETIC:
    if (scw == nullptr) {
      rtErr("This stage of the integration process must be called with a ScoreCard abstract for "
            "energy tracking.  Use a different overloaded function variant.",
            "launchIntegrationProcess.");
    }
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:

      // A system's unit cell type does not have direct bearing on the choice of kinetic energy
      // computation kernel.  However, if simulations are taking place in isolated boundary
      // conditions, the kernel launch grid may be optimized in one way or another.  Therefore,
      // check the unit cell type and launch the correct kernel variant.
      switch (poly_psw->unit_cell) {
      case UnitCellType::NONE:
        switch (kwidth) {
        case ValenceKernelSize::XL:
          kfsKineticEnergyXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
          break;
        case ValenceKernelSize::LG:
          kfsKineticEnergyLG<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
          break;
        case ValenceKernelSize::MD:
          kfsKineticEnergyMD<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
          break;
        case ValenceKernelSize::SM:
          kfsKineticEnergySM<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
          break;
        }
        break;
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        kfsKineticEnergyXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
        break;
      }
      break;
    case AccumulationMethod::WHOLE:
      switch (poly_psw->unit_cell) {
      case UnitCellType::NONE:
        switch (kwidth) {
        case ValenceKernelSize::XL:
          kfKineticEnergyXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
          break;
        case ValenceKernelSize::LG:
          kfKineticEnergyLG<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
          break;
        case ValenceKernelSize::MD:
          kfKineticEnergyMD<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
          break;
        case ValenceKernelSize::SM:
          kfKineticEnergySM<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
          break;
        }
        break;
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        kfKineticEnergyXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
        break;
      }
    }
    break;
  case IntegrationStage::POSITION_ADVANCE:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      switch (poly_psw->unit_cell) {
      case UnitCellType::NONE:
        switch (kwidth) {
        case ValenceKernelSize::XL:
          kfsIntegPosAdvXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::LG:
          kfsIntegPosAdvLG<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::MD:
          kfsIntegPosAdvMD<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::SM:
          kfsIntegPosAdvSM<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        }
        break;
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        kfsIntegPosAdvXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
        break;
      }
      break;
    case AccumulationMethod::WHOLE:
      switch (poly_psw->unit_cell) {
      case UnitCellType::NONE:
        switch (kwidth) {
        case ValenceKernelSize::XL:
          kfIntegPosAdvXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::LG:
          kfIntegPosAdvLG<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::MD:
          kfIntegPosAdvMD<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::SM:
          kfIntegPosAdvSM<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        }
        break;
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        kfIntegPosAdvXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
        break;
      }
      break;
    }
    break;
  case IntegrationStage::GEOMETRY_CONSTRAINT:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      switch (poly_psw->unit_cell) {
      case UnitCellType::NONE:
        switch (kwidth) {
        case ValenceKernelSize::XL:
          kfsIntegGeomCnstXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::LG:
          kfsIntegGeomCnstLG<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::MD:
          kfsIntegGeomCnstMD<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::SM:
          kfsIntegGeomCnstSM<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        }
        break;
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        kfsPmeIntegGeomCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
        break;
      }
      break;
    case AccumulationMethod::WHOLE:
      switch (poly_psw->unit_cell) {
      case UnitCellType::NONE:
        switch (kwidth) {
        case ValenceKernelSize::XL:
          kfIntegGeomCnstXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::LG:
          kfIntegGeomCnstLG<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::MD:
          kfIntegGeomCnstMD<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        case ValenceKernelSize::SM:
          kfIntegGeomCnstSM<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
          break;
        }
        break;
      case UnitCellType::ORTHORHOMBIC:
      case UnitCellType::TRICLINIC:
        kfPmeIntegGeomCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
        break;
      }
      break;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<float> *tb_resw,
                              MMControlKit<float> *ctrl, const SyValenceKit<float> &poly_vk,
                              const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                              const ThermostatWriter<float> &tstw, const int2 lp,
                              const AccumulationMethod acc_meth, const ValenceKernelSize kwidth,
                              const IntegrationStage process) {
  launchIntegrationProcess(poly_psw, tb_resw, ctrl, nullptr, poly_vk, poly_auk, tstw, lp, acc_meth,
                           kwidth, process);
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<double> *tb_resw,
                              MMControlKit<double> *ctrl, ScoreCardWriter *scw,
                              const CellGridReader<double, llint, double, double4_16a> &cgr,
                              const SyValenceKit<double> &poly_vk,
                              const SyAtomUpdateKit<double, double2, double4_16a> &poly_auk,
                              const ThermostatWriter<double> &tstw, const int2 lp,
                              const IntegrationStage process) {
  switch (process) {
  case IntegrationStage::CALC_FORCES:
    break;
  case IntegrationStage::VELOCITY_ADVANCE:
    kdsdPmeIntegVelAdv<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, cgr, poly_auk, tstw, *tb_resw);
    break;
  case IntegrationStage::VELOCITY_CONSTRAINT:
    kdsPmeIntegVelCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
    break;
  case IntegrationStage::CALC_KINETIC:
    if (scw == nullptr) {
      rtErr("This stage of the integration process must be called with a ScoreCard abstract for "
            "energy tracking.  Use a different overloaded function variant.",
            "launchIntegrationProcess.");
    }
    kdsKineticEnergy<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
    break;
  case IntegrationStage::POSITION_ADVANCE:
    kdsIntegPosAdv<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
    break;
  case IntegrationStage::GEOMETRY_CONSTRAINT:
    kdsPmeIntegGeomCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<double> *tb_resw,
                              MMControlKit<double> *ctrl,
                              const CellGridReader<double, llint, double, double4_16a> &cgr,
                              const SyValenceKit<double> &poly_vk,
                              const SyAtomUpdateKit<double, double2, double4_16a> &poly_auk,
                              const ThermostatWriter<double> &tstw, const int2 lp,
                              const IntegrationStage process) {
  launchIntegrationProcess(poly_psw, tb_resw, ctrl, nullptr, cgr, poly_vk, poly_auk, tstw, lp,
                           process);
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<double> *tb_resw,
                              MMControlKit<double> *ctrl, ScoreCardWriter *scw,
                              const CellGridReader<float, int, float, float4> &cgr,
                              const SyValenceKit<double> &poly_vk,
                              const SyAtomUpdateKit<double, double2, double4_16a> &poly_auk,
                              const ThermostatWriter<double> &tstw, const int2 lp,
                              const IntegrationStage process) {
  switch (process) {
  case IntegrationStage::CALC_FORCES:
    break;
  case IntegrationStage::VELOCITY_ADVANCE:
    kdsfPmeIntegVelAdv<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, cgr, poly_auk, tstw, *tb_resw);
    break;
  case IntegrationStage::VELOCITY_CONSTRAINT:
    kdsPmeIntegVelCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
    break;
  case IntegrationStage::CALC_KINETIC:
    if (scw == nullptr) {
      rtErr("This stage of the integration process must be called with a ScoreCard abstract for "
            "energy tracking.  Use a different overloaded function variant.",
            "launchIntegrationProcess.");
    }
    kdsKineticEnergy<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
    break;
  case IntegrationStage::POSITION_ADVANCE:
    kdsIntegPosAdv<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
    break;
  case IntegrationStage::GEOMETRY_CONSTRAINT:
    kdsPmeIntegGeomCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<double> *tb_resw,
                              MMControlKit<double> *ctrl,
                              const CellGridReader<float, int, float, float4> &cgr,
                              const SyValenceKit<double> &poly_vk,
                              const SyAtomUpdateKit<double, double2, double4_16a> &poly_auk,
                              const ThermostatWriter<double> &tstw, const int2 lp,
                              const IntegrationStage process) {
  launchIntegrationProcess(poly_psw, tb_resw, ctrl, nullptr, cgr, poly_vk, poly_auk, tstw, lp,
                           process);
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<double> *tb_resw,
                              MMControlKit<double> *ctrl, ScoreCardWriter *scw,
                              const CellGridReader<double, llint, double, double4_16a> &cgr_qq,
                              const CellGridReader<double, llint, double, double4_16a> &cgr_lj,
                              const SyValenceKit<double> &poly_vk,
                              const SyAtomUpdateKit<double, double2, double4_16a> &poly_auk,
                              const ThermostatWriter<double> &tstw, const int2 lp,
                              const IntegrationStage process) {
  switch (process) {
  case IntegrationStage::CALC_FORCES:
    break;
  case IntegrationStage::VELOCITY_ADVANCE:
    kdsdPmeDualIntegVelAdv<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, cgr_qq, cgr_lj, poly_auk,
                                           tstw, *tb_resw);
    break;
  case IntegrationStage::VELOCITY_CONSTRAINT:
    kdsPmeIntegVelCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
    break;
  case IntegrationStage::CALC_KINETIC:
    if (scw == nullptr) {
      rtErr("This stage of the integration process must be called with a ScoreCard abstract for "
            "energy tracking.  Use a different overloaded function variant.",
            "launchIntegrationProcess.");
    }
    kdsKineticEnergy<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
    break;
  case IntegrationStage::POSITION_ADVANCE:
    kdsIntegPosAdv<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
    break;
  case IntegrationStage::GEOMETRY_CONSTRAINT:
    kdsPmeIntegGeomCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<double> *tb_resw,
                              MMControlKit<double> *ctrl,
                              const CellGridReader<double, llint, double, double4_16a> &cgr_qq,
                              const CellGridReader<double, llint, double, double4_16a> &cgr_lj,
                              const SyValenceKit<double> &poly_vk,
                              const SyAtomUpdateKit<double, double2, double4_16a> &poly_auk,
                              const ThermostatWriter<double> &tstw, const int2 lp,
                              const IntegrationStage process) {
  launchIntegrationProcess(poly_psw, tb_resw, ctrl, nullptr, cgr_qq, cgr_lj, poly_vk, poly_auk,
                           tstw, lp, process);
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<double> *tb_resw,
                              MMControlKit<double> *ctrl, ScoreCardWriter *scw,
                              const CellGridReader<float, int, float, float4> &cgr_qq,
                              const CellGridReader<float, int, float, float4> &cgr_lj,
                              const SyValenceKit<double> &poly_vk,
                              const SyAtomUpdateKit<double, double2, double4_16a> &poly_auk,
                              const ThermostatWriter<double> &tstw, const int2 lp,
                              const IntegrationStage process) {
  switch (process) {
  case IntegrationStage::CALC_FORCES:
    break;
  case IntegrationStage::VELOCITY_ADVANCE:
    kdsfPmeDualIntegVelAdv<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, cgr_qq, cgr_lj, poly_auk,
                                           tstw, *tb_resw);
    break;
  case IntegrationStage::VELOCITY_CONSTRAINT:
    kdsPmeIntegVelCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
    break;
  case IntegrationStage::CALC_KINETIC:
    if (scw == nullptr) {
      rtErr("This stage of the integration process must be called with a ScoreCard abstract for "
            "energy tracking.  Use a different overloaded function variant.",
            "launchIntegrationProcess.");
    }
    kdsKineticEnergy<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
    break;
  case IntegrationStage::POSITION_ADVANCE:
    kdsIntegPosAdv<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
    break;
  case IntegrationStage::GEOMETRY_CONSTRAINT:
    kdsPmeIntegGeomCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<double> *tb_resw,
                              MMControlKit<double> *ctrl,
                              const CellGridReader<float, int, float, float4> &cgr_qq,
                              const CellGridReader<float, int, float, float4> &cgr_lj,
                              const SyValenceKit<double> &poly_vk,
                              const SyAtomUpdateKit<double, double2, double4_16a> &poly_auk,
                              const ThermostatWriter<double> &tstw, const int2 lp,
                              const IntegrationStage process) {
  launchIntegrationProcess(poly_psw, tb_resw, ctrl, nullptr, cgr_qq, cgr_lj, poly_vk, poly_auk,
                           tstw, lp, process);
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<float> *tb_resw,
                              MMControlKit<float> *ctrl, ScoreCardWriter *scw,
                              const CellGridReader<double, llint, double, double4_16a> &cgr,
                              const SyValenceKit<float> &poly_vk,
                              const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                              const ThermostatWriter<float> &tstw, const int2 lp,
                              const AccumulationMethod acc_meth, const IntegrationStage process) {
  switch (process) {
  case IntegrationStage::CALC_FORCES:
    break;
  case IntegrationStage::VELOCITY_ADVANCE:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsdPmeIntegVelAdv<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, cgr, poly_auk, tstw, *tb_resw);
      break;
    case AccumulationMethod::WHOLE:
      kfdPmeIntegVelAdv<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, cgr, poly_auk, tstw, *tb_resw);
      break;
    }
    break;
  case IntegrationStage::VELOCITY_CONSTRAINT:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsPmeIntegVelCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    case AccumulationMethod::WHOLE:
      kfPmeIntegVelCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    }
    break;
  case IntegrationStage::CALC_KINETIC:
    if (scw == nullptr) {
      rtErr("This stage of the integration process must be called with a ScoreCard abstract for "
            "energy tracking.  Use a different overloaded function variant.",
            "launchIntegrationProcess.");
    }
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsKineticEnergyXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
      break;
    case AccumulationMethod::WHOLE:
      kfKineticEnergyXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
      break;
    }
    break;
  case IntegrationStage::POSITION_ADVANCE:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsIntegPosAdvXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    case AccumulationMethod::WHOLE:
      kfIntegPosAdvXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    }
    break;
  case IntegrationStage::GEOMETRY_CONSTRAINT:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsPmeIntegGeomCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    case AccumulationMethod::WHOLE:
      kfPmeIntegGeomCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<float> *tb_resw,
                              MMControlKit<float> *ctrl,
                              const CellGridReader<double, llint, double, double4_16a> &cgr,
                              const SyValenceKit<float> &poly_vk,
                              const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                              const ThermostatWriter<float> &tstw, const int2 lp,
                              const AccumulationMethod acc_meth, const IntegrationStage process) {
  launchIntegrationProcess(poly_psw, tb_resw, ctrl, nullptr, cgr, poly_vk, poly_auk, tstw, lp,
                           acc_meth, process);
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<float> *tb_resw,
                              MMControlKit<float> *ctrl, ScoreCardWriter *scw,
                              const CellGridReader<float, int, float, float4> &cgr,
                              const SyValenceKit<float> &poly_vk,
                              const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                              const ThermostatWriter<float> &tstw, const int2 lp,
                              const AccumulationMethod acc_meth, const IntegrationStage process) {
  switch (process) {
  case IntegrationStage::CALC_FORCES:
    break;
  case IntegrationStage::VELOCITY_ADVANCE:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsfPmeIntegVelAdv<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, cgr, poly_auk, tstw, *tb_resw);
      break;
    case AccumulationMethod::WHOLE:
      kffPmeIntegVelAdv<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, cgr, poly_auk, tstw, *tb_resw);
      break;
    }
    break;
  case IntegrationStage::VELOCITY_CONSTRAINT:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsPmeIntegVelCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    case AccumulationMethod::WHOLE:
      kfPmeIntegVelCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    }
    break;
  case IntegrationStage::CALC_KINETIC:
    if (scw == nullptr) {
      rtErr("This stage of the integration process must be called with a ScoreCard abstract for "
            "energy tracking.  Use a different overloaded function variant.",
            "launchIntegrationProcess.");
    }
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsKineticEnergyXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
      break;
    case AccumulationMethod::WHOLE:
      kfKineticEnergyXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
      break;
    }
    break;
  case IntegrationStage::POSITION_ADVANCE:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsIntegPosAdvXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    case AccumulationMethod::WHOLE:
      kfIntegPosAdvXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    }
    break;
  case IntegrationStage::GEOMETRY_CONSTRAINT:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsPmeIntegGeomCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    case AccumulationMethod::WHOLE:
      kfPmeIntegGeomCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<float> *tb_resw,
                              MMControlKit<float> *ctrl,
                              const CellGridReader<float, int, float, float4> &cgr,
                              const SyValenceKit<float> &poly_vk,
                              const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                              const ThermostatWriter<float> &tstw, const int2 lp,
                              const AccumulationMethod acc_meth, const IntegrationStage process) {
  launchIntegrationProcess(poly_psw, tb_resw, ctrl, nullptr, cgr, poly_vk, poly_auk, tstw, lp,
                           acc_meth, process);
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<float> *tb_resw,
                              MMControlKit<float> *ctrl, ScoreCardWriter *scw,
                              const CellGridReader<double, llint, double, double4_16a> &cgr_qq,
                              const CellGridReader<double, llint, double, double4_16a> &cgr_lj,
                              const SyValenceKit<float> &poly_vk,
                              const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                              const ThermostatWriter<float> &tstw, const int2 lp,
                              const AccumulationMethod acc_meth, const IntegrationStage process) {
  switch (process) {
  case IntegrationStage::CALC_FORCES:
    break;
  case IntegrationStage::VELOCITY_ADVANCE:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsdPmeDualIntegVelAdv<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, cgr_qq, cgr_lj, poly_auk,
                                             tstw, *tb_resw);
      break;
    case AccumulationMethod::WHOLE:
      kfdPmeDualIntegVelAdv<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, cgr_qq, cgr_lj, poly_auk,
                                            tstw, *tb_resw);
      break;
    }
    break;
  case IntegrationStage::VELOCITY_CONSTRAINT:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsPmeIntegVelCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    case AccumulationMethod::WHOLE:
      kfPmeIntegVelCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    }
    break;
  case IntegrationStage::CALC_KINETIC:
    if (scw == nullptr) {
      rtErr("This stage of the integration process must be called with a ScoreCard abstract for "
            "energy tracking.  Use a different overloaded function variant.",
            "launchIntegrationProcess.");
    }
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsKineticEnergyXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
      break;
    case AccumulationMethod::WHOLE:
      kfKineticEnergyXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
      break;
    }
    break;
  case IntegrationStage::POSITION_ADVANCE:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsIntegPosAdvXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    case AccumulationMethod::WHOLE:
      kfIntegPosAdvXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    }
    break;
  case IntegrationStage::GEOMETRY_CONSTRAINT:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsPmeIntegGeomCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    case AccumulationMethod::WHOLE:
      kfPmeIntegGeomCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<float> *tb_resw,
                              MMControlKit<float> *ctrl,
                              const CellGridReader<double, llint, double, double4_16a> &cgr_qq,
                              const CellGridReader<double, llint, double, double4_16a> &cgr_lj,
                              const SyValenceKit<float> &poly_vk,
                              const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                              const ThermostatWriter<float> &tstw, const int2 lp,
                              const AccumulationMethod acc_meth, const IntegrationStage process) {
  launchIntegrationProcess(poly_psw, tb_resw, ctrl, nullptr, cgr_qq, cgr_lj, poly_vk, poly_auk,
                           tstw, lp, acc_meth, process);
}
  
//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<float> *tb_resw,
                              MMControlKit<float> *ctrl, ScoreCardWriter *scw,
                              const CellGridReader<float, int, float, float4> &cgr_qq,
                              const CellGridReader<float, int, float, float4> &cgr_lj,
                              const SyValenceKit<float> &poly_vk,
                              const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                              const ThermostatWriter<float> &tstw, const int2 lp,
                              const AccumulationMethod acc_meth, const IntegrationStage process) {
  switch (process) {
  case IntegrationStage::CALC_FORCES:
    break;
  case IntegrationStage::VELOCITY_ADVANCE:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsfPmeDualIntegVelAdv<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, cgr_qq, cgr_lj, poly_auk,
                                             tstw, *tb_resw);
      break;
    case AccumulationMethod::WHOLE:
      kffPmeDualIntegVelAdv<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, cgr_qq, cgr_lj, poly_auk,
                                            tstw, *tb_resw);
      break;
    }
    break;
  case IntegrationStage::VELOCITY_CONSTRAINT:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsPmeIntegVelCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    case AccumulationMethod::WHOLE:
      kfPmeIntegVelCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    }
    break;
  case IntegrationStage::CALC_KINETIC:
    if (scw == nullptr) {
      rtErr("This stage of the integration process must be called with a ScoreCard abstract for "
            "energy tracking.  Use a different overloaded function variant.",
            "launchIntegrationProcess.");
    }
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsKineticEnergyXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
      break;
    case AccumulationMethod::WHOLE:
      kfKineticEnergyXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, *scw);
      break;
    }
    break;
  case IntegrationStage::POSITION_ADVANCE:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsIntegPosAdvXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    case AccumulationMethod::WHOLE:
      kfIntegPosAdvXL<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    }
    break;
  case IntegrationStage::GEOMETRY_CONSTRAINT:
    switch (acc_meth) {
    case AccumulationMethod::SPLIT:
    case AccumulationMethod::AUTOMATIC:
      kfsPmeIntegGeomCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    case AccumulationMethod::WHOLE:
      kfPmeIntegGeomCnst<<<lp.x, lp.y>>>(poly_vk, *ctrl, *poly_psw, poly_auk, tstw, *tb_resw);
      break;
    }
    break;
  }
}
  
//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PsSynthesisWriter *poly_psw, CacheResourceKit<float> *tb_resw,
                              MMControlKit<float> *ctrl,
                              const CellGridReader<float, int, float, float4> &cgr_qq,
                              const CellGridReader<float, int, float, float4> &cgr_lj,
                              const SyValenceKit<float> &poly_vk,
                              const SyAtomUpdateKit<float, float2, float4> &poly_auk,
                              const ThermostatWriter<float> &tstw, const int2 lp,
                              const AccumulationMethod acc_meth, const IntegrationStage process) {
  launchIntegrationProcess(poly_psw, tb_resw, ctrl, nullptr, cgr_qq, cgr_lj, poly_vk, poly_auk,
                           tstw, lp, acc_meth, process);
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PhaseSpaceSynthesis *poly_ps, CacheResource *tb_res, Thermostat *tst,
                              MolecularMechanicsControls *mmctrl, ScoreCard *sc,
                              const AtomGraphSynthesis &poly_ag, const CoreKlManager &launcher,
                              const PrecisionModel prec, const AccumulationMethod acc_meth,
                              const IntegrationStage process) {
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(devc_tier);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      CacheResourceKit<double> tb_resw = tb_res->dpData(devc_tier);
      MMControlKit<double> ctrl = mmctrl->dpData(devc_tier);
      SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(devc_tier);
      SyAtomUpdateKit<double,
                      double2,
                      double4_16a> poly_auk = poly_ag.getDoublePrecisionAtomUpdateKit(devc_tier);
      ThermostatWriter<double> tstw = tst->dpData(devc_tier);
      const int2 lp = launcher.getIntegrationKernelDims(prec, acc_meth, process);
      if (sc == nullptr) {
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, nullptr, poly_vk, poly_auk, tstw, lp,
                                 process);
      }
      else {
        ScoreCardWriter scw = sc->data(devc_tier);
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, &scw, poly_vk, poly_auk, tstw, lp,
                                 process);
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      CacheResourceKit<float> tb_resw = tb_res->spData(devc_tier);
      MMControlKit<float> ctrl = mmctrl->spData(devc_tier);
      SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(devc_tier);
      SyAtomUpdateKit<float,
                      float2,
                      float4> poly_auk = poly_ag.getSinglePrecisionAtomUpdateKit(devc_tier);
      ThermostatWriter<float> tstw = tst->spData(devc_tier);
      const int2 lp = launcher.getIntegrationKernelDims(prec, acc_meth, process);
      const ValenceKernelSize kwidth = poly_ag.getValenceThreadBlockSize();
      if (sc == nullptr) {
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, nullptr, poly_vk, poly_auk, tstw, lp,
                                 acc_meth, kwidth, process);
      }
      else {
        ScoreCardWriter scw = sc->data(devc_tier);
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, &scw, poly_vk, poly_auk, tstw, lp,
                                 acc_meth, kwidth, process);
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PhaseSpaceSynthesis *poly_ps, CacheResource *tb_res, Thermostat *tst,
                              MolecularMechanicsControls *mmctrl,
                              const AtomGraphSynthesis &poly_ag, const CoreKlManager &launcher,
                              const PrecisionModel prec, const AccumulationMethod acc_meth,
                              const IntegrationStage process) {
  launchIntegrationProcess(poly_ps, tb_res, tst, mmctrl, nullptr, poly_ag, launcher, prec,
                           acc_meth, process);
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PhaseSpaceSynthesis *poly_ps, CacheResource *tb_res, Thermostat *tst,
                              MolecularMechanicsControls *mmctrl, ScoreCard *sc,
                              const CellGrid<double, llint, double, double4_16a> &cg,
                              const AtomGraphSynthesis &poly_ag, const CoreKlManager &launcher,
                              const PrecisionModel prec, const AccumulationMethod acc_meth,
                              const IntegrationStage process) {
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(devc_tier);
  const CellGridReader<double, llint, double, double4_16a> cgr = cg.data(devc_tier);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      CacheResourceKit<double> tb_resw = tb_res->dpData(devc_tier);
      MMControlKit<double> ctrl = mmctrl->dpData(devc_tier);
      SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(devc_tier);
      SyAtomUpdateKit<double,
                      double2,
                      double4_16a> poly_auk = poly_ag.getDoublePrecisionAtomUpdateKit(devc_tier);
      ThermostatWriter<double> tstw = tst->dpData(devc_tier);
      const int2 lp = launcher.getIntegrationKernelDims(prec, acc_meth, process);
      if (sc == nullptr) {
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, nullptr, cgr, poly_vk, poly_auk, tstw,
                                 lp, process);
      }
      else {
        ScoreCardWriter scw = sc->data(devc_tier);
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, &scw, cgr, poly_vk, poly_auk, tstw,
                                 lp, process);
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      CacheResourceKit<float> tb_resw = tb_res->spData(devc_tier);
      MMControlKit<float> ctrl = mmctrl->spData(devc_tier);
      SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(devc_tier);
      SyAtomUpdateKit<float,
                      float2,
                      float4> poly_auk = poly_ag.getSinglePrecisionAtomUpdateKit(devc_tier);
      ThermostatWriter<float> tstw = tst->spData(devc_tier);
      const int2 lp = launcher.getIntegrationKernelDims(prec, acc_meth, process);
      if (sc == nullptr) {
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, nullptr, cgr, poly_vk, poly_auk, tstw,
                                 lp, acc_meth, process);
      }
      else {
        ScoreCardWriter scw = sc->data(devc_tier);
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, &scw, cgr, poly_vk, poly_auk, tstw,
                                 lp, acc_meth, process);
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PhaseSpaceSynthesis *poly_ps, CacheResource *tb_res, Thermostat *tst,
                              MolecularMechanicsControls *mmctrl,
                              const CellGrid<double, llint, double, double4_16a> &cg,
                              const AtomGraphSynthesis &poly_ag, const CoreKlManager &launcher,
                              const PrecisionModel prec, const AccumulationMethod acc_meth,
                              const IntegrationStage process) {
  launchIntegrationProcess(poly_ps, tb_res, tst, mmctrl, nullptr, cg, poly_ag, launcher, prec,
                           acc_meth, process);
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PhaseSpaceSynthesis *poly_ps, CacheResource *tb_res, Thermostat *tst,
                              MolecularMechanicsControls *mmctrl, ScoreCard *sc,
                              const CellGrid<float, int, float, float4> &cg,
                              const AtomGraphSynthesis &poly_ag, const CoreKlManager &launcher,
                              const PrecisionModel prec, const AccumulationMethod acc_meth,
                              const IntegrationStage process) {
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(devc_tier);
  const CellGridReader<float, int, float, float4> cgr = cg.data(devc_tier);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      CacheResourceKit<double> tb_resw = tb_res->dpData(devc_tier);
      MMControlKit<double> ctrl = mmctrl->dpData(devc_tier);
      SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(devc_tier);
      SyAtomUpdateKit<double,
                      double2,
                      double4_16a> poly_auk = poly_ag.getDoublePrecisionAtomUpdateKit(devc_tier);
      ThermostatWriter<double> tstw = tst->dpData(devc_tier);
      const int2 lp = launcher.getIntegrationKernelDims(prec, acc_meth, process);
      if (sc == nullptr) {
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, nullptr, cgr, poly_vk, poly_auk,
                                 tstw, lp, process);
      }
      else {
        ScoreCardWriter scw = sc->data(devc_tier);
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, &scw, cgr, poly_vk, poly_auk, tstw,
                                 lp, process);
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      CacheResourceKit<float> tb_resw = tb_res->spData(devc_tier);
      MMControlKit<float> ctrl = mmctrl->spData(devc_tier);
      SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(devc_tier);
      SyAtomUpdateKit<float,
                      float2,
                      float4> poly_auk = poly_ag.getSinglePrecisionAtomUpdateKit(devc_tier);
      ThermostatWriter<float> tstw = tst->spData(devc_tier);
      const int2 lp = launcher.getIntegrationKernelDims(prec, acc_meth, process);
      if (sc == nullptr) {
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, nullptr, cgr, poly_vk, poly_auk, tstw,
                                 lp, acc_meth, process);
      }
      else {
        ScoreCardWriter scw = sc->data(devc_tier);
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, &scw, cgr, poly_vk, poly_auk, tstw,
                                 lp, acc_meth, process);
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PhaseSpaceSynthesis *poly_ps, CacheResource *tb_res, Thermostat *tst,
                              MolecularMechanicsControls *mmctrl,
                              const CellGrid<float, int, float, float4> &cg,
                              const AtomGraphSynthesis &poly_ag, const CoreKlManager &launcher,
                              const PrecisionModel prec, const AccumulationMethod acc_meth,
                              const IntegrationStage process) {
  launchIntegrationProcess(poly_ps, tb_res, tst, mmctrl, nullptr, cg, poly_ag, launcher, prec,
                           acc_meth, process);
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PhaseSpaceSynthesis *poly_ps, CacheResource *tb_res, Thermostat *tst,
                              MolecularMechanicsControls *mmctrl, ScoreCard *sc,
                              const CellGrid<double, llint, double, double4_16a> &cg_qq,
                              const CellGrid<double, llint, double, double4_16a> &cg_lj,
                              const AtomGraphSynthesis &poly_ag, const CoreKlManager &launcher,
                              const PrecisionModel prec, const AccumulationMethod acc_meth,
                              const IntegrationStage process) {
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(devc_tier);
  const CellGridReader<double, llint, double, double4_16a> cgr_qq = cg_qq.data(devc_tier);
  const CellGridReader<double, llint, double, double4_16a> cgr_lj = cg_lj.data(devc_tier);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      CacheResourceKit<double> tb_resw = tb_res->dpData(devc_tier);
      MMControlKit<double> ctrl = mmctrl->dpData(devc_tier);
      SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(devc_tier);
      SyAtomUpdateKit<double,
                      double2,
                      double4_16a> poly_auk = poly_ag.getDoublePrecisionAtomUpdateKit(devc_tier);
      ThermostatWriter<double> tstw = tst->dpData(devc_tier);
      const int2 lp = launcher.getIntegrationKernelDims(prec, acc_meth, process);
      if (sc == nullptr) {
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, nullptr, cgr_qq, cgr_lj, poly_vk,
                                 poly_auk, tstw, lp, process);
      }
      else {
        ScoreCardWriter scw = sc->data(devc_tier);
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, &scw, cgr_qq, cgr_lj, poly_vk,
                                 poly_auk, tstw, lp, process);
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      CacheResourceKit<float> tb_resw = tb_res->spData(devc_tier);
      MMControlKit<float> ctrl = mmctrl->spData(devc_tier);
      SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(devc_tier);
      SyAtomUpdateKit<float,
                      float2,
                      float4> poly_auk = poly_ag.getSinglePrecisionAtomUpdateKit(devc_tier);
      ThermostatWriter<float> tstw = tst->spData(devc_tier);
      const int2 lp = launcher.getIntegrationKernelDims(prec, acc_meth, process);
      if (sc == nullptr) {
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, nullptr, cgr_qq, cgr_lj, poly_vk,
                                 poly_auk, tstw, lp, acc_meth, process);
      }
      else {
        ScoreCardWriter scw = sc->data(devc_tier);
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, &scw, cgr_qq, cgr_lj, poly_vk,
                                 poly_auk, tstw, lp, acc_meth, process);
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PhaseSpaceSynthesis *poly_ps, CacheResource *tb_res, Thermostat *tst,
                              MolecularMechanicsControls *mmctrl,
                              const CellGrid<double, llint, double, double4_16a> &cg_qq,
                              const CellGrid<double, llint, double, double4_16a> &cg_lj,
                              const AtomGraphSynthesis &poly_ag, const CoreKlManager &launcher,
                              const PrecisionModel prec, const AccumulationMethod acc_meth,
                              const IntegrationStage process) {
  launchIntegrationProcess(poly_ps, tb_res, tst, mmctrl, nullptr, cg_qq, cg_lj, poly_ag, launcher,
                           prec, acc_meth, process);
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PhaseSpaceSynthesis *poly_ps, CacheResource *tb_res, Thermostat *tst,
                              MolecularMechanicsControls *mmctrl, ScoreCard *sc,
                              const CellGrid<float, int, float, float4> &cg_qq,
                              const CellGrid<float, int, float, float4> &cg_lj,
                              const AtomGraphSynthesis &poly_ag, const CoreKlManager &launcher,
                              const PrecisionModel prec, const AccumulationMethod acc_meth,
                              const IntegrationStage process) {
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  PsSynthesisWriter poly_psw = poly_ps->data(devc_tier);
  const CellGridReader<float, int, float, float4> cgr_qq = cg_qq.data(devc_tier);
  const CellGridReader<float, int, float, float4> cgr_lj = cg_lj.data(devc_tier);
  switch (prec) {
  case PrecisionModel::DOUBLE:
    {
      CacheResourceKit<double> tb_resw = tb_res->dpData(devc_tier);
      MMControlKit<double> ctrl = mmctrl->dpData(devc_tier);
      SyValenceKit<double> poly_vk = poly_ag.getDoublePrecisionValenceKit(devc_tier);
      SyAtomUpdateKit<double,
                      double2,
                      double4_16a> poly_auk = poly_ag.getDoublePrecisionAtomUpdateKit(devc_tier);
      ThermostatWriter<double> tstw = tst->dpData(devc_tier);
      const int2 lp = launcher.getIntegrationKernelDims(prec, acc_meth, process);
      if (sc == nullptr) {
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, nullptr, cgr_qq, cgr_lj, poly_vk,
                                 poly_auk, tstw, lp, process);
      }
      else {
        ScoreCardWriter scw = sc->data(devc_tier);
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, &scw, cgr_qq, cgr_lj, poly_vk,
                                 poly_auk, tstw, lp, process);
      }
    }
    break;
  case PrecisionModel::SINGLE:
    {
      CacheResourceKit<float> tb_resw = tb_res->spData(devc_tier);
      MMControlKit<float> ctrl = mmctrl->spData(devc_tier);
      SyValenceKit<float> poly_vk = poly_ag.getSinglePrecisionValenceKit(devc_tier);
      SyAtomUpdateKit<float,
                      float2,
                      float4> poly_auk = poly_ag.getSinglePrecisionAtomUpdateKit(devc_tier);
      ThermostatWriter<float> tstw = tst->spData(devc_tier);
      const int2 lp = launcher.getIntegrationKernelDims(prec, acc_meth, process);
      const ValenceKernelSize kwidth = poly_ag.getValenceThreadBlockSize();
      if (sc == nullptr) {
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, nullptr, cgr_qq, cgr_lj, poly_vk,
                                 poly_auk, tstw, lp, acc_meth, process);
      }
      else {
        ScoreCardWriter scw = sc->data(devc_tier);
        launchIntegrationProcess(&poly_psw, &tb_resw, &ctrl, &scw, cgr_qq, cgr_lj, poly_vk,
                                 poly_auk, tstw, lp, acc_meth, process);
      }
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void launchIntegrationProcess(PhaseSpaceSynthesis *poly_ps, CacheResource *tb_res, Thermostat *tst,
                              MolecularMechanicsControls *mmctrl,
                              const CellGrid<float, int, float, float4> &cg_qq,
                              const CellGrid<float, int, float, float4> &cg_lj,
                              const AtomGraphSynthesis &poly_ag, const CoreKlManager &launcher,
                              const PrecisionModel prec, const AccumulationMethod acc_meth,
                              const IntegrationStage process) {
  launchIntegrationProcess(poly_ps, tb_res, tst, mmctrl, nullptr, cg_qq, cg_lj, poly_ag, launcher,
                           prec, acc_meth, process);
}

} // namespace trajectory
} // namespace stormm
