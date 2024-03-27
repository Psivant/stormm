// -*-c++-*-
#ifndef STORMM_CORE_KERNEL_MANAGER_H
#define STORMM_CORE_KERNEL_MANAGER_H

#include <map>
#include <string>
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
#    include <cuda_runtime.h>
#  endif
#endif
#include "copyright.h"
#include "Constants/behavior.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/reduction.h"
#include "Numerics/split_fixed_precision.h"
#include "Potential/energy_enumerators.h"
#include "Structure/structure_enumerators.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/trajectory_enumerators.h"
#include "gpu_details.h"
#include "kernel_format.h"

namespace stormm {
namespace card {

using constants::PrecisionModel;
using energy::ClashResponse;
using energy::EvaluateForce;
using energy::EvaluateEnergy;
using energy::getEnumerationName;
using energy::QMapMethod;
using energy::ValenceKernelSize;
using numerics::AccumulationMethod;
using structure::GridDetail;
using structure::RMSDTask;
using structure::VirtualSiteActivity;
using synthesis::AtomGraphSynthesis;
using synthesis::NbwuKind;
using synthesis::VwuGoal;
using stmath::ReductionGoal;
using stmath::ReductionStage;
using synthesis::VwuGoal;
using topology::ImplicitSolventModel;
using topology::UnitCellType;
using trajectory::IntegrationStage;
  
/// \brief A class to guide the implementation of GPU kernels, with selected thread counts per
///        block and block counts per launch grid for a specific GPU based on the workload.  This
///        object is constructed in a stepwise manner, with each kind of work unit contributing
///        new launch specifications.
class CoreKlManager : public KernelManager {
public:

  /// \brief The constructor will fill in values as if this were a single-threaded CPU "launch."
  ///
  /// \param gpu_in   Details of the GPU in use (this is relevant, as it will be used to interpret
  ///                 the layout of any kernels)
  /// \param poly_ag  Topologies for all systems, offering details of the workload
  CoreKlManager(const GpuDetails &gpu_in, const AtomGraphSynthesis &poly_ag);

  /// \brief Get the architecture-specific block multiplier.  This will run a minimum number of
  ///        blocks per streaming multiprocessor on some cards, specifically NVIDIA's GTX 1080-Ti,
  ///        when large blocks cannot use more than 32k registers in all.
  int getArchBlockMultiplier() const;
  
  /// \brief Get the block and thread counts for a valence kernel.
  ///
  /// \param prec                The type of floating point numbers in which the kernel shall work
  /// \param eval_force          Indication of whether the kernel will evaluate forces on atoms
  /// \param eval_nrg            Indication of whether to evaluate the energy of the system as a
  ///                            whole
  /// \param acc_meth            The force accumulation method (SPLIT or WHOLE, AUTOMATIC will
  ///                            produce an error in this context)
  /// \param collision_handling  Indication of whether clashes are to be dampened
  int2 getValenceKernelDims(PrecisionModel prec, EvaluateForce eval_force, EvaluateEnergy eval_nrg,
                            AccumulationMethod acc_meth, VwuGoal purpose,
                            ClashResponse collision_handling) const;

  /// \brief Get the block and thread counts for a stand-alone velocity constraint kernel.
  ///
  /// \param prec      The precision model for arithmetic to be used in applying constraints
  /// \param acc_meth  The method in which the velocities will be stored in __shared__ memory
  ///                  arrays.  The distinction comes about due to the multi-functional nature of
  ///                  the core code in this kernel: stand-alone or included in a valence
  ///                  interaction kernel.
  /// \param process   The part of the integration workflow to perform
  int2 getIntegrationKernelDims(PrecisionModel prec, AccumulationMethod acc_meth,
                                IntegrationStage process) const;

  /// \brief Get the block and thread counts for a non-bonded kernel.  Parameter descriptions
  ///        for this function follow from getValenceKernelDims() above, with the addition of:
  ///
  /// \param kind  The type of non-bonded work unit: tile groups, supertiles, or honeycomb
  /// \param igb   Type of implicit solvent model to use ("NONE" in periodic boundary conditions)
  int2 getNonbondedKernelDims(PrecisionModel prec, NbwuKind kind, EvaluateForce eval_force,
                              EvaluateEnergy eval_nrg, AccumulationMethod acc_meth,
                              ImplicitSolventModel igb, ClashResponse collision_handling) const;

  /// \brief Get the block and thread counts for a Born radii computation kernel.  Parameter
  ///        descriptions for this function follow from getValenceKernelDims() and
  ///        getNonbondedKernelDims(), above.
  int2 getBornRadiiKernelDims(PrecisionModel prec, NbwuKind kind, AccumulationMethod acc_meth,
                              ImplicitSolventModel igb) const;

  /// \brief Get the block and thread counts for a Born derivative computation kernel.  Parameter
  ///        descriptions for this function follow from getValenceKernelDims() and
  ///        getNonbondedKernelDims(), above.
  int2 getBornDerivativeKernelDims(PrecisionModel prec, NbwuKind kind, AccumulationMethod acc_meth,
                                   ImplicitSolventModel igb) const;

  /// \brief Get the launch parameters for a general-purpose particle-mesh density mapping kernel.
  ///        These kernels loop over atoms in a cell grid (see cellgrid.h) and issues atomic
  ///        instructions to accumulate density on a set of particle-mesh interaction grids.
  ///
  /// Overloaded:
  ///   - Provide only the calculation precision (valid for GENERAL_PURPOSE and ACC_REGISTER
  ///     kernels)
  ///   - Also provide the accumulation precision (needed for ACC_SHARED kernels)
  ///
  /// \param approach   The mapping approach guiding which kernel to choose
  /// \param prec       Precision in which all calculations will be performed.  The precision of
  ///                   the coordinate representation is a separate parameter, but does not tend to
  ///                   affect the number of registers used by the cards and is thus not a factor
  ///                   when determining the launch grid size.
  /// \param calc_prec  Precision in which all calculations will be performed (supply this for the
  ///                   ACC_SHARED approach)
  /// \param acc_prec   Precision in which accumulations will be expressed.  95-bit accumulators
  ///                   are used for real-valued double-precision grids, 63-bit accumulators are
  ///                   used for real-valued single-precision grids)
  /// \param overflow   Indicate whether the overflow bits might be engaged by the fixed-precision
  ///                   model
  /// \param cg_tmat    Type identifier code for the coordinate representation in the associated
  ///                   cell grid (a neighbor list spatial decomposition of the coordinate
  ///                   synthesis)
  /// \param order      Order of the B-spline based interpolation
  /// \{
  int2 getDensityMappingKernelDims(QMapMethod approach, PrecisionModel prec, size_t cg_tmat,
                                   int order) const;

  int2 getDensityMappingKernelDims(QMapMethod approach, PrecisionModel calc_prec,
                                   PrecisionModel acc_prec, bool overflow, size_t cg_tmat,
                                   int order) const;
  /// \}

  /// \brief Get the block and thread counts for a reduction kernel.
  ///
  /// \param prec     The type of floating point numbers represented by the kernel's substrate
  ///                 (all reductions happenin double-precision real numbers)
  /// \param purpose  The procedure requested, for which the appropriate kernel shall be found 
  /// \param process  The reduction step to perform
  int2 getReductionKernelDims(PrecisionModel prec, ReductionGoal purpose,
                              ReductionStage process) const;

  /// \brief Get the block and thread counts for a virtual site placement or force transmission
  ///        kernel.
  ///
  /// \param prec     Implies a level of detail in the fixed-precision coordinate representation,
  ///                 and specifies the precision in which to perform the placement or force
  ///                 transmission operations
  /// \param purpose  The process to perform with standalone virtual site kernels: place virtual
  ///                 sites or transmit forces from them to frame atoms with mass
  int2 getVirtualSiteKernelDims(PrecisionModel prec, VirtualSiteActivity purpose) const;

  /// \brief Get the launch parameters for an RMSD calculation kernel.
  ///
  /// \param prec   The precision model in which to perform calculations and present results
  /// \param order  The order of the calculation (all to reference, or all to all)
  int2 getRMSDKernelDims(PrecisionModel prec, RMSDTask order) const;
  
private:

  /// This parameter controls the valence kernel block multiplier.  There will always be at least
  /// two blocks of the valence work unit kernel resident on each SM, and it can have at most 1024
  /// threads which gives it no architecture dependence.  However, for different valence
  /// interaction kernels there can be different thread counts per block and different block
  /// multipliers associated with each kernel's size.
  ValenceKernelSize valence_kernel_width;

  /// The architecture-specific block multiplier for non-bonded kernels, essentially a provision
  /// for NVIDIA Turing cards.  The non-bonded kernels are not subject to the same constraints and
  /// block multiplicity as the valence kernels, as the non-bonded kernels will always run at
  /// least two blocks per streaming multiprocessor.  This value does depend on the unit cell type
  /// of the topology synthesis at hand.
  /// \{
  int nonbond_block_multiplier_dp;
  int nonbond_block_multiplier_sp;
  /// \}
  
  /// Architecture-specific block multipliers for Generalized Born radii computation kernels, again
  /// a provision for NVIDIA Turing cards.
  /// \{
  int gbradii_block_multiplier_dp;
  int gbradii_block_multiplier_sp;
  /// \}

  /// Architecture-specific block multipliers for Generalized Born radii derivative computation
  /// kernels, again a provision for NVIDIA Turing cards.
  /// \{
  int gbderiv_block_multiplier_dp;
  int gbderiv_block_multiplier_sp;
  /// \}

  /// Architecture-specific block multipliers for mapping kernels.  The specific multiplier for
  /// interpolation order k at a particular level of precision is given by the kth index of either
  /// array.  The "gen" prefix indicates a general-purpose mapping kernel, the "rac" prefix
  /// indicates register-centric accumulation, and "sac" indicates __shared__ memory accumulation.
  /// \{
  std::vector<int> gen_qmap_block_multiplier_dp;
  std::vector<int> gen_qmap_block_multiplier_sp;
  std::vector<int> sac_qmap_block_multiplier_dp;
  std::vector<int> sac_qmap_block_multiplier_sp;
  /// \}
  
  /// The workload-specific block multipliers for reduction kernels.  Like the valence kernels, the
  /// thread count per streaming multiprocessor will not go above 1024 (this time out of bandwidth
  /// limitations), but the block multiplicity (which starts at 4) could be increased.
  int reduction_block_multiplier;
  
  /// The architecture-specific block multipliers for virtual site handling kernels are implied by
  /// the sizes of valence work units and their recommended launch parameters, which contain the
  /// instructions for placing virtual sites.
  ValenceKernelSize virtual_site_kernel_width;

  /// The architecture-specific block multipliers for RMSD computation kernels.  The thread count
  /// for each kernel is set to 128 (tiny).  Each kernel operates strictly as a collection of
  /// warps.
  /// \{
  int rmsd_block_multiplier_dp;
  int rmsd_block_multiplier_sp;
  /// \}

  /// \brief Set the register, maximum block size, and thread counts for one of the valence
  ///        kernels.  This function complements getValenceKernelDims(), although the other
  ///        function reports numbers based on this functions input information and some further
  ///        analysis.
  ///
  /// \param prec                Type of floating point numbers in which the kernel shall work
  /// \param eval_force          Indication of whether the kernel will evaluate forces on atoms
  /// \param eval_nrg            Indication of whether to evaluate the energy of the system as a
  ///                            whole
  /// \param acc_meth            The force accumulation method (SPLIT or WHOLE, AUTOMATIC will
  ///                            produce an error in this context)
  /// \param collision_handling  Indication of whether collisions are to be forgiven
  /// \param kwidth              The size selection of the valence kernel, an explicit control over
  ///                            the kernel subdivision
  /// \param kernel_name         [Optional] Name of the kernel in the actual code
  void catalogValenceKernel(PrecisionModel prec, EvaluateForce eval_force, EvaluateEnergy eval_nrg,
                            AccumulationMethod acc_meth, VwuGoal purpose,
                            ClashResponse collision_handling, ValenceKernelSize kwidth,
                            const std::string &kernel_name = std::string(""));

  /// \brief Set the maximum block size and thread counts for one of the stand-alone trajectory
  ///        integration kernels.  Descriptions of input parameters follow from
  ///        catalogValenceKernel(), above, in addition to:
  ///
  /// \param process  The integration process to query, e.g. position advancement
  void catalogIntegrationKernel(PrecisionModel prec, AccumulationMethod acc_meth,
                                ValenceKernelSize kwidth, IntegrationStage process,
                                const std::string &kernel_name = std::string(""));
  
  /// \brief Set the register, maximum block size, and thread counts for one of the non-bonded
  ///        kernels.  Parameter descriptions for this function follow from
  ///        catalogValenceKernel() above, with the addition of:
  ///
  /// \param kind  Type of non-bonded work unit: tile groups, supertiles, or honeycomb
  /// \param igb   Type of implicit solvent model to use ("NONE" in periodic boundary conditions)
  void catalogNonbondedKernel(PrecisionModel prec, NbwuKind kind, EvaluateForce eval_force,
                              EvaluateEnergy eval_nrg, AccumulationMethod acc_meth,
                              ImplicitSolventModel igb, ClashResponse collision_handling,
                              const std::string &kernel_name = std::string(""));

  /// \brief Set the register, maximum block size, and thread counts for one of the Generalized
  ///        Born radius computation kernels.  Parameter descriptions for this function follow
  ///        from catalogValenceKernel() and catalogNonbondedKernel() above.
  void catalogBornRadiiKernel(PrecisionModel prec, NbwuKind kind, AccumulationMethod acc_meth,
                              ImplicitSolventModel igb,
                              const std::string &kernel_name = std::string(""));

  /// \brief Set the register, maximum block size, and thread counts for one of the Generalized
  ///        Born derivative computation kernels.  Parameter descriptions for this function follow
  ///        from catalogValenceKernel() and catalogNonbondedKernel() above.
  void catalogBornDerivativeKernel(PrecisionModel prec, NbwuKind kind, AccumulationMethod acc_meth,
                                   ImplicitSolventModel igb,
                                   const std::string &kernel_name = std::string(""));

  /// \brief Set the block multiplier for a general-purpose density mapping kernel.  The term
  ///        "general-purpose" stands in opposition to other methods for accomplishing the same
  ///        task involving special-purpose kernels which may not be applicable in all cases.
  ///
  /// \param prec         Indicate whether the kernel will perform arithmetic in single- or
  ///                     double-precision floating point numbers
  /// \param cg_tmat      Type index code of the coordinate (and spatial decomposition cell
  ///                     dimensions matrix) representation in the CellGrid object that will be
  ///                     submitted to the kernel
  /// \param order        The order of B-spline interpolation, defining the extent to which density
  ///                     is spread onto each system's particle-mesh interaction grid
  /// \param kernel_name  [Optional] Name of the kernel in the actual code
  void catalogGeneralQMapKernel(PrecisionModel prec, size_t cg_tmat, int order,
                                const std::string &kernel_name);

  /// \brief Set the block multiplier for a register-centered density mapping kernel.  These
  ///        kernels can take the place of the general mapping kernel in select, useful cases.
  ///        Descriptions of input parameters follow from catalogGeneralQMapKernel(), above, with
  ///        the distinctions:
  ///
  /// \param calc_prec  Indicate whether the kernel will perform arithmetic in single- or
  ///                   double-precision floating point numbers
  /// \param acc_prec   Indicate whether to accumulate in 95-bit or 63-bit split fixed-precision
  ///                   for eventual conversion to double- or single-precision real numbers,
  ///                   respectively
  /// \param overflow   Indicate whether the fixed-precision accumulation model might need the
  ///                   overflow accumulator bits
  void catalogShrAccQMapKernel(PrecisionModel calc_prec, PrecisionModel acc_prec,
                               bool overflow, size_t cg_tmat, int order,
                               const std::string &kernel_name);

  /// \brief Set the register, maximum block size, and thread counts for one of the reduction
  ///        kernels.  Parameter descriptions for this function follow from
  ///        setValenceKernelAttributes() above, with the addition of:
  ///
  /// \param purpose      Reason for doing the reduction, i.e. conjugate gradient transformation
  /// \param process      How far to take the reduction operation
  void catalogReductionKernel(PrecisionModel prec, ReductionGoal purpose, ReductionStage process,
                              int subdivision, const std::string &kernel_name = std::string(""));

  /// \brief Set the register, maximum block size, and thread counts for one of the virtual site
  ///        placement kernels.  The first two parameters follow from getVirtualSiteKernelDims(),
  ///        above.
  ///
  /// \param prec         The type of floating point numbers in which the kernel shall work
  /// \param purpose      The process to perform with standalone virtual site kernels
  /// \param kwidth       Width of the kernel required by the valence work units it will process.
  ///                     This has some bearing on the efficiency of virtual site manipulation in
  ///                     the standalone kernels.
  /// \param kernel_name  [Optional] Name of the kernel in the actual code
  void catalogVirtualSiteKernel(PrecisionModel prec, VirtualSiteActivity purpose,
                                ValenceKernelSize kwidth,
                                const std::string &kernel_name = std::string(""));

  /// \brief Set the maximum block size and thread counts for one of the RMSD calculation kernels.
  ///
  /// \param prec         The type of floating point numbers in which the kernel shall work
  /// \param order        The order of the calculation (all to reference, or all to all)
  /// \param kernel_name  [Optional] Name of the kernel in the actual code
  void catalogRMSDKernel(PrecisionModel prec, RMSDTask order,
                         const std::string &kernel_name = std::string(""));
};

/// \brief Obtain the architecture-specific block multiplier for non-bonded interaction kernels.
///
/// \param gpu        Details of the GPU that will perform the calculations
/// \param unit_cell  The unit cell type of the systems to evaluate
/// \param prec       The type of floating point numbers in which the kernel shall work
/// \param igb        The implicit solvent model that the kernel shall implement
int nonbondedBlockMultiplier(const GpuDetails &gpu, UnitCellType unit_cell, PrecisionModel prec,
                             ImplicitSolventModel igb);

/// \brief Obtain the architecture-specific block multiplier for Generalized Born radii kernels.
///
/// \param gpu   Details of the GPU that will perform the calculations
/// \param prec  The type of floating point numbers in which the kernel shall work
int gbRadiiBlockMultiplier(const GpuDetails &gpu, PrecisionModel prec);

/// \brief Obtain the architecture-specific block multiplier for Generalized Born derivative
///        computation kernels.
///
/// \param gpu   Details of the GPU that will perform the calculations
/// \param prec  The type of floating point numbers in which the kernel shall work
int gbDerivativeBlockMultiplier(const GpuDetails &gpu, PrecisionModel prec);

/// \brief Obtain the architecture-specific block multipler for general-purpose particle density
///        (charge or dispersion) mapping kernels.
///
/// \param gpu       Details of the GPU that will perform the calculations
/// \param prec      The type of floating point numbers in which the kernel shall work
/// \param approach  The method of density mapping, each of which defines a unique set of kernels
std::vector<int> densityMappingBlockMultiplier(const GpuDetails &gpu, PrecisionModel prec,
                                               QMapMethod approach);
  
/// \brief Obtain the workload-specific block multiplier for reduction kernels.
int reductionBlockMultiplier();

/// \brief Obtain the workload-specific block multiplier for virtual site handling kernels.  The
///        typical kernel will run on 256 threads per block.
///
/// \param prec  The type of floating point numbers in which the kernel shall work
int virtualSiteBlockMultiplier(PrecisionModel prec);

/// \brief Obtain the block multiplier for RMSD calculation kernels.  Each block will run on 128
///        threads.
///
/// \param prec  The type of floating point numbers in which the kernel shall work
int rmsdBlockMultiplier(PrecisionModel prec);

/// \brief Codify the kernel extension for various valence-related kernels.
///
/// \param prec    The precision in which the kernel operates.  Double-precision kernels are
///                generally not subject to modifications of the size of the thread blocks.
/// \param kwidth  The recommended width of kernels that will process valence work units
std::string valenceKernelWidthExtension(PrecisionModel prec, ValenceKernelSize kwidth);
  
/// \brief Obtain a unique string identifier for one of the valence kernels.  Each identifier
///        begins with "vale_" and is then appended with letter codes for different aspects
///        according to the following system:
///        - { d, f }            Perform calculations in double (d) or float (f) arithmetic
///        - { e, f, fe }        Compute energies (e), forces (f), or both (ef)
///        - { s, w }            Accumulate forces in split integers (s) or whole integers (w)
///        - { m, a }            Move atoms (m) or accumulate forces or energies (a
///        - { cl, nc }          Take no action (cl) in the event of a collision, or implement an
///                              increase in the perceived distance between particles (nc) 
///        - { xl, lg, md, sm }  Use an extra large (xl), large (lg), medium (md), or small (sm)
///                              kernel thread block
///
/// \param prec                Type of floating point numbers in which the kernel shall work
/// \param eval_force          Indication of whether the kernel will evaluate forces on atoms
/// \param eval_nrg            Indication of whether to evaluate the energy of the system as a
///                            whole
/// \param acc_meth            The force accumulation method (SPLIT or WHOLE, AUTOMATIC will
///                            produce an error in this context)
/// \param purpose             The intended action to take with computed forces
/// \param collision_handling  Indication of whether to dampen collision effects between particles
/// \param kwidth              The desired kernel width, from "extra large" to "small"
std::string valenceKernelKey(PrecisionModel prec, EvaluateForce eval_force,
                             EvaluateEnergy eval_nrg, AccumulationMethod acc_meth,
                             VwuGoal purpose, ClashResponse collision_handling,
                             ValenceKernelSize kwidth);

/// \brief Obtain a unique string identifier for one of the integration kernels.  Each identifier
///        begins with "intg_" and is then appended with letter codes for different attributes
///        according to the following system:
///        - { d, f }            Perform calculations in double (d) or float (f) arithmetic
///        - { s, w }            Accumulate forces in split integers (s) or whole integers (w)
///        - { xl, lg, md, sm }  Use an extra large (xl), large (lg), medium (md), or small (sm)
///                              kernel thread block
///        - { va, vc, pa, gc }  Four possible stages of the integration time step, between force
///                              calculations: velocity advancement, velocity constraints, position
///                              advancement, and geometric constraints
///                               
/// Descriptions of input parameters follow from valenceKernelKey(), above.
std::string integrationKernelKey(PrecisionModel prec, AccumulationMethod acc_meth,
                                 ValenceKernelSize kwidth, IntegrationStage process);
  
/// \brief Obtain a unique string identifier for one of the non-bonded kernels.  Each identifier
///        begins with "nonb_" and is then appended with letter codes for different aspects
///        according to the following system:
///        - { d, f }           Perform calculations in double (d) or float (f) arithmetic
///        - { tg, st, nl }     Use a "tile groups" or "supertiles" strategy for breaking down
///                             systems with isolated boundary conditions, or a "neighbor list"
///                             strategy for breaking down systems in periodic boundary conditions
///        - { vac, gbs, gbn }  Perform calculations in vacuum, standard GB implicit solvent, or
///                             "neck" GB implicit solvent (the latter two being relevant only for
///                             systems in isolated boundary conditions)
///        - { e, f, fe }       Compute energies (e), forces (f), or both (fe)
///        - { s, w }           Accumulate forces in split integers (s) or whole integers (w)
///        - { cl, nc }         Take no action (cl) in the event of a collision, or implement an
///                             increase in the perceived distance between particles (nc) 
///
/// \param prec                The type of floating point numbers in which the kernel shall work
/// \param kind                The type of non-bonded work unit to evaluate
/// \param eval_force          Indication of whether the kernel will evaluate forces on atoms
/// \param eval_nrg            Indication of whether to evaluate the energy of the system as a
///                            whole
/// \param acc_meth            The force accumulation method (SPLIT or WHOLE, AUTOMATIC will
///                            produce an error in this context)
/// \param igb                 Type of implicit solvent model ("NONE" in periodic boundary
///                            conditions)
/// \param collision_handling  Indication of whether to dampen collision effects between particles
std::string nonbondedKernelKey(PrecisionModel prec, NbwuKind kind, EvaluateForce eval_force,
                               EvaluateEnergy eval_nrg, AccumulationMethod acc_meth,
                               ImplicitSolventModel igb, ClashResponse collision_handling);

/// \brief Encapsulate the work of encoding a Generalized Born computation kernel key, shared
///        across radii and derivative computations.  Parameter descriptions, and the features of
///        the key, follow from nonbondedKernelKey() above.
std::string appendBornKernelKey(PrecisionModel prec, NbwuKind kind, AccumulationMethod acc_meth,
                                ImplicitSolventModel igb);
  
/// \brief Obtain a unique string identifier for one of the Born radii computation kernels.  Each
///        identifier begins with "gbrd_" and is then appended with letter codes for different
///        aspects of the kernel, following the codes set forth in nonbondedKernelKey() above.
///        Parameter descriptions also follow from nonbondedKernelKey().
std::string bornRadiiKernelKey(PrecisionModel prec, NbwuKind kind, AccumulationMethod acc_meth,
                               ImplicitSolventModel igb);

/// \brief Obtain a unique string identifier for one of the Born radii computation kernels.  Each
///        identifier begins with "gbrd_" and is then appended with letter codes for different
///        aspects of the kernel, following the codes set forth in nonbondedKernelKey() above.
///        Parameter descriptions also follow from nonbondedKernelKey().
std::string bornDerivativeKernelKey(PrecisionModel prec, NbwuKind kind,
                                    AccumulationMethod acc_meth, ImplicitSolventModel igb);

/// \brief Produce a descriptive, coded string which can be appended to a mapping kernel name
///        prefix.  The result contains letter codes for different precision models and the order
///        of interpolation:
///        - { s, l }  Short or long representations of the particle coordinates
///        - { i, r }  Integer (fixed precision) or real coordinate representations
///        - { d, f }  Double- or single-precision calculations
///
/// \param prec     The type of floating-point numbers in which the kernel performs calculations
/// \param cg_tmat  Type identifier code for the coordinate representation in the underlying
///                 CellGrid neighbor list 
/// \param order    Order of B-spline interpolation taking particle density onto the particle-mesh
///                 interaction grids
std::string appendQMapKernelKey(PrecisionModel prec, size_t cg_tmat, int order);

/// \brief Obtain a unique string identifier for one of the general-purpose density mapping
///        kernels.  Each identifier begins with "qmap_" and is then appended with a coded
///        string produced by appendQMapKernelKey(), above.  Descriptions of input parameters also
///        follow from that function.
std::string generalQMapKernelKey(PrecisionModel prec, size_t cg_tmat, int order);

/// \brief Obtain a unique string identifier for one of the density mapping kernels based on
///        accumulation in _shared__ memory.  Each identifier begins with "qmap_sacc_" and is then
///        appended with letter codes generated by appendQMapKernelKey(), above.  Descriptions of
///        input parameters also follow from that function, with the distinctions:
///
/// \param calc_prec  Indicate whether the kernel will perform arithmetic in single- or
///                   double-precision floating point numbers
/// \param acc_prec   Indicate whether to accumulate in 95-bit or 63-bit split fixed-precision for
///                   eventual conversion to double- or single-precision real numbers, respectively
/// \param overflow   Indicate whether the fixed-precision model is such that overflow bits might
///                   be needed
std::string shrAccQMapKernelKey(PrecisionModel calc_prec, PrecisionModel acc_prec,
                                bool overflow, size_t cg_tmat,  int order);

/// \brief Obtain a unique string identifier for one of the reduction kernels.  Each identifier
///        begins with "redc_" and is then appended with letter codes for different aspects
///        according to the following system:
///        - { d, f }        Perform calculations in double (d) or float (f) arithmetic
///        - { cg }          Calculate the conjugate gradient vector
///        - { gt, sc, rd }  Perform a gathering, scattering, or combined all-reduce operation
///
/// \param prec     The type of floating point numbers in which the kernel shall work
/// \param purpose  The reason for doing the reduction--the purpose of each kernel is unique
/// \param process  The reduction stage to perform
std::string reductionKernelKey(PrecisionModel prec, ReductionGoal purpose, ReductionStage process);

/// \brief Obtain a unique string identifier for one of the virtual site handling kernels.  Each
///        identifier begins with "vste_" and is then appended with letter codes for different
///        activities according to the following system:
///        - { d, f }            Perform calculations in double (d) or float (f) arithmetic
///        - { pl, xm }          Place particles or transmit forces to atoms with mass
///        - { xl, lg, md, sm }  Use an extra large (xl), large (lg), medium (md), or small (sm)
///                              kernel thread block
///
/// \param prec     The type of floating point numbers in which the kernel shall work
/// \param purpose  The process to perform with standalone virtual site kernels
/// \param kwidth   The recommended thread block size for valence-related kernels
std::string virtualSiteKernelKey(PrecisionModel prec, VirtualSiteActivity process,
                                 ValenceKernelSize kwidth);

/// \brief Obtain a unique string identifier for one of the RMSD calculation kernels.  Each
///        identifier begins with "rmsd_" and is then appended with letter codes for different
///        activities according to the following system:
///        - { d, f }  Perform calculations in double (d) or float (f) arithmetic
///        - { r, m }  Compute RMSD with respect to a reference structure or create matrices of
///                    all-to-all RMSD values
///
/// \param prec   The type of floating point numbers in which the kernel shall work
/// \param order  The order of the calculation (all to reference, or all to all)
std::string rmsdKernelKey(PrecisionModel prec, RMSDTask order);

} // namespace card
} // namespace stormm

#endif
