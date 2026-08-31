// -*-c++-*-
#ifndef STORMM_CONVOLUTION_MANAGER_H
#define STORMM_CONVOLUTION_MANAGER_H

#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
#    include <cuda_runtime.h>
#    include <cufft.h>
#  endif
#endif
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "DataTypes/stormm_type_name.h"
#include "Math/fft_stage.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "pme_util.h"
#include "pmigrid.h"
#include "scorecard.h"

namespace stormm {
namespace energy {

using constants::CartesianDimension;
using constants::PrecisionModel;
using constants::UnitCellAxis;
using card::Hybrid;
using card::HybridTargetLevel;
using data_types::getStormmTypeName;
using stmath::FFTStage;
using synthesis::AtomGraphSynthesis;

/// \brief Any abstract of the convolution manager is read-only.
template <typename T, typename T2> struct ConvolutionWriter {

  /// \brief As with other abstracts, the convolution kit is constructed with a series of values
  ///        for each associated pointer and critical constant.  See Essmann's 1995 paper (below)
  ///        for  the inspiration to some of the nomenclature.
  ConvolutionWriter(int system_count_in, int wu_count_in, T ew_coeff_in, const int* prf_offsets_in,
                    const T* self_ecorr_in, const T* bmesh_a_in, const T* bmesh_b_in,
                    const T* bmesh_c_in, const T* mval_a_in, const T* mval_b_in,
                    const T* mval_c_in, const T* msval_a_in, const T* msval_b_in,
                    const T* msval_c_in, const T* cmesh_a_in, const T* cmesh_b_in,
                    const T* cmesh_c_in, T2* freq_data_in, const uint* frq_offsets_in,
                    const uint4* wu_list_in, std::vector<FFTStage> *fft_ops_in);

  /// \brief Like other abstracts, the presence of const members implicitly forbids the copy and
  ///        move assignment operators.  Copy and move constructors may be taken in their deafult
  ///        forms.
  ///
  /// \param original  The original object to copy or move
  /// \{
  ConvolutionWriter(const ConvolutionWriter &original) = default;
  ConvolutionWriter(ConvolutionWriter &&original) = default;
  /// \}

  const int system_count;   ///< The number of indepedent systems
  const int wu_count;       ///< The number of HPC work units
  const T ew_coeff;         ///< The Ewald coefficien in use by the splitting function
  const int* prf_offsets;   ///< Offsets for grid elements and prefactors in each system.  The
                            ///<   maximum of all grid dimensions for a particular system, rounded
                            ///<   to the nearest multiple of the warp size, determines each
                            ///<   successive system offset in the arrays for prefactors and M
                            ///<   values along any of the A, B, or C unit cell dimensions.
  const T* self_ecorr;      ///< Self-correlation energy of the charges (electrostatic or
                            ///<    otherwise) in each system.
  const T* bmesh_a;         ///< "B" mesh prefactors for each system along the unit cell A axis
  const T* bmesh_b;         ///< "B" mesh prefactors for each system along the unit cell B axis
  const T* bmesh_c;         ///< "B" mesh prefactors for each system along the unit cell C axis
  const T* mval_a;          ///< "M" values for each system along the unit cell A axis
  const T* mval_b;          ///< "M" values for each system along the unit cell B axis
  const T* mval_c;          ///< "M" values for each system along the unit cell C axis
  const T* msval_a;         ///< Shifted "M" values for each system along the unit cell A axis
  const T* msval_b;         ///< Shifted "M" values for each system along the unit cell B axis
  const T* msval_c;         ///< Shifted "M" values for each system along the unit cell C axis
  const T* cmesh_a;         ///< "C" mesh prefactors for each system along the unit cell A axis.
                            ///<   These are only computed if all systems' unit cells are
                            ///<   orthorhombic.
  const T* cmesh_b;         ///< "C" mesh prefactors for each system along the unit cell B axis.
  const T* cmesh_c;         ///< "C" mesh prefactors for each system along the unit cell C axis.
  T2* freq_data;            ///< Pointer to a contiguous array of frequency transforms for all
                            ///<   systems.  This data may reside within the underlying PMIGrid,
                            ///<   if in-place FFTs are being performed, or within the source
                            ///<   ConvolutionManager if out-of-place FFTs are being performed.
  const uint* frq_offsets;  ///< Offsets for the frequency grids of each system.  These offsets are
                            ///<   set based on the real-space offsets of the underlying PMIGrid
                            ///<   object, if in-place FFTs are in called for, or based on a padded
                            ///<   layout of all problems' frequency-space grids in private data
                            ///<   allocated by the source ConvolutionManager.
  const uint4* wu_list;     ///< The list of work units guiding HPC applications of the Green's
                            ///<   function in a load-balanced and secure manner
  
  /// The array of FFT operations is accessible from the convolution manager's abstract.  While
  /// this information cannot be accessed in a meaningful way by a GPU kernel, the functions that
  /// accept the abstract and ultimately launch kernels to apply the Green's function will thus be
  /// able to launch the FFT operations before and afterwards.
  std::vector<FFTStage> *fft_ops;
};
  
/// \brief Any abstract of the convolution manager is read-only.
template <typename T, typename T2> struct ConvolutionReader {

  /// \brief As with other abstracts, the convolution kit is constructed with a series of values
  ///        for each associated pointer and critical constant.  See Essmann's 1995 paper (below)
  ///        for  the inspiration to some of the nomenclature.
  ConvolutionReader(int system_count_in, T ew_coeff_in, const int* prf_offsets_in,
                    const T* self_ecorr_in, const T* bmesh_a_in, const T* bmesh_b_in,
                    const T* bmesh_c_in, const T* mval_a_in, const T* mval_b_in,
                    const T* mval_c_in, const T* msval_a_in, const T* msval_b_in,
                    const T* msval_c_in, const T* cmesh_a_in, const T* cmesh_b_in,
                    const T* cmesh_c_in, const T2* freq_data_in, const uint* frq_offsets_in,
                    const std::vector<FFTStage> *fft_ops_in);

  /// \brief As with other abstracts, the reader can be constructed based on an equivalent writer.
  ///
  /// Overloaded:
  ///   - Provide the writer by reference
  ///   - Provide a pointer to the writer
  ///
  /// \param w  The original writer to transmute into a read-only object
  /// \{
  ConvolutionReader(const ConvolutionWriter<T, T2> *w);
  ConvolutionReader(const ConvolutionWriter<T, T2> &w);
  /// \}
  
  /// \brief Like other abstracts, the presence of const members implicitly forbids the copy and
  ///        move assignment operators.  Copy and move constructors may be taken in their deafult
  ///        forms.
  ///
  /// \param original  The original object to copy or move
  /// \{
  ConvolutionReader(const ConvolutionReader &original) = default;
  ConvolutionReader(ConvolutionReader &&original) = default;
  /// \}

  const int system_count;   ///< The number of indepedent systems
  const T ew_coeff;         ///< The Ewald coefficien in use by the splitting function
  const int* prf_offsets;   ///< Offsets for grid elements and prefactors in each system.  The
                            ///<   maximum of all grid dimensions for a particular system, rounded
                            ///<   to the nearest multiple of the warp size, determines each
                            ///<   successive system offset in the arrays for prefactors and M
                            ///<   values along any of the A, B, or C unit cell dimensions.
  const T* self_ecorr;      ///< Self-correlation energy of the charges (electrostatic or
                            ///<    otherwise) in each system.
  const T* bmesh_a;         ///< "B" mesh prefactors for each system along the unit cell A axis
  const T* bmesh_b;         ///< "B" mesh prefactors for each system along the unit cell B axis
  const T* bmesh_c;         ///< "B" mesh prefactors for each system along the unit cell C axis
  const T* mval_a;          ///< "M" values for each system along the unit cell A axis
  const T* mval_b;          ///< "M" values for each system along the unit cell B axis
  const T* mval_c;          ///< "M" values for each system along the unit cell C axis
  const T* msval_a;         ///< Shifted "M" values for each system along the unit cell A axis
  const T* msval_b;         ///< Shifted "M" values for each system along the unit cell B axis
  const T* msval_c;         ///< Shifted "M" values for each system along the unit cell C axis
  const T* cmesh_a;         ///< "C" mesh prefactors for each system along the unit cell A axis.
                            ///<   These are only computed if all systems' unit cells are
                            ///<   orthorhombic.
  const T* cmesh_b;         ///< "C" mesh prefactors for each system along the unit cell B axis.
  const T* cmesh_c;         ///< "C" mesh prefactors for each system along the unit cell C axis.
  const T2* freq_data;      ///< Pointer to a contiguous array of frequency transforms for all
                            ///<   systems.  This data may reside within the underlying PMIGrid,
                            ///<   if in-place FFTs are being performed, or within the source
                            ///<   ConvolutionManager if out-of-place FFTs are being performed.
  const uint* frq_offsets;  ///< Offsets for the frequency grids of each system.  These offsets are
                            ///<   set based on the real-space offsets of the underlying PMIGrid
                            ///<   object, if in-place FFTs are in called for, or based on a padded
                            ///<   layout of all problems' frequency-space grids in private data
                            ///<   allocated by the source ConvolutionManager.
  
  /// The array of FFT operations is accessible from the convolution manager's abstract.  While
  /// the array's individual elements cannot be used to launch kernels or perform CPU-based FFT
  /// operations due to the const qualification, the dimensions of each FFT may be of interest to
  /// a function using the read-only abstract.
  const std::vector<FFTStage> *fft_ops;
};
  
/// \brief Collect elements for performing the reciprocal space convolution in many systems.  This
///        class works most directly in conjunction with a PMIGrid object, and is expected to
///        reference the same PhaseSpaceSynthesis and CellGrid objects.  This object will also
///        point to the underlying AtomGraphSynthesis (topology synthesis).  Terminology in this
///        object follows from the 1995 Smooth Particle Mesh Ewald publication:
///
/// Ulrich Essmann, Lalith Perera, Max L. Berkowitz, Tom Darden, Hsing Lee, and Lee G. Pedersen.
/// (1995) "A Smooth Particle Mesh Ewald Method." Journal of Chemical Physics, 103:8577-8593.
class ConvolutionManager {
public:

  /// \brief The constructor depends on a PMIGrid and will refer back to the PMIGrid's associated
  ///        CellGrid object to retrieve the pointer to its topology synthesis.
  /// \{
  ConvolutionManager(const PMIGrid *pmig_in, const double ewald_coefficient_in,
                     const GpuDetails &gpu = null_gpu);

  ConvolutionManager(const PMIGrid &pmig_in, const double ewald_coefficient_in,
                     const GpuDetails &gpu = null_gpu);
  /// \}

  /// \brief The copy and move constructors, as well as copy and move assignemnt operators, must
  ///        all be given explicit definitions due to the presence of POINTER-kind Hybrid objects.
  ///
  /// \param original  The original object to copy or move
  /// \param other     Another object placed on the right hand side of the assignment operation
  /// \{
  ConvolutionManager(const ConvolutionManager &original);
  ConvolutionManager(ConvolutionManager &&original);
  ConvolutionManager& operator=(const ConvolutionManager &other);
  ConvolutionManager& operator=(ConvolutionManager &&other);
  /// \}
  
  /// \brief Get the number of systems in the associated synthesis.
  int getSystemCount() const;

  /// \brief Get the number of HPC work units applicable to the collection of all systems.
  int getHPCWorkUnitCount() const;
  
  /// \brief Get the Ewald coefficient used by the convolution.
  double getEwaldCoefficient() const;

  /// \brief Get the definition of Coulomb's constant used by the convolution.  This will return
  ///        Coulomb's constant as defined in the associated topology synthesis.
  double getCoulombConstant() const;

  /// \brief Get the Particle-Mesh Ewald grid dimensions for any one system.
  const PMIGrid* getPMIGridPointer() const;

  /// \brief Get a pointer to the topology synthesis.
  const AtomGraphSynthesis* getTopologySynthesisPointer() const;

  /// \brief Report the self energies of charges for each system, for inspection.
  ///
  /// \param prec  Indicate whether to draw from the SINGLE- or DOUBLE-precision array
  std::vector<double> getSelfEcorr(PrecisionModel prec) const;

  /// \brief Get the number of FFT batch groups that the object will perform in order to process
  ///        FFTs on all systems.
  int getFFTGroupCount() const;

  /// \brief Get the list of systems comprised by a specific FFT group.
  ///
  /// \param group_index  Index of the group of interest.  This will be checked for validity.
  std::vector<int> getFFTGroupSystemList(int group_index) const;
  
  /// \brief Get the level at which FFTs and the convolution will be performed.
  HybridTargetLevel getOperatingTier() const;

  /// \brief Get the double-precision abstract at the object's selected operating tier.  The tier
  ///        selection is automatic, unlike many other STORMM classes.
  ///
  /// Overloaded:
  ///   - Get an immutable abstract for a const-qualified object
  ///   - Get a writeable abstract for a non-const object
  /// \{
  const ConvolutionReader<double, double2> dpData() const;
  ConvolutionWriter<double, double2> dpData();
  /// \}

  /// \brief Get the double-precision abstract at the object's selected operating tier.  The tier
  ///        selection is automatic, unlike many other STORMM classes.
  ///
  /// Overloaded:
  ///   - Get an immutable abstract for a const-qualified object
  ///   - Get a writeable abstract for a non-const object
  /// \{
  const ConvolutionReader<float, float2> spData() const;
  ConvolutionWriter<float, float2> spData();
  /// \}

  /// \brief Get the object's abstract with explicit templating.  This function can be called
  ///        based on the templated type of another templated function.  Overloading follows from
  ///        dpData() and spData(), above, and like those functions the pointers in the abstract
  ///        will be oriented towards data on the CPU host or GPU device as appropriate.
  /// \{
  template <typename T, typename T2> const ConvolutionReader<T, T2> data() const;
  template <typename T, typename T2> ConvolutionWriter<T, T2> data();
  /// \}

  /// \brief Get a pointer to the ConvolutionManager object itself.
  const ConvolutionManager* getSelfPointer() const;
  
  /// \brief Run the forward FFT operations for all systems.  This will call the eponymous member
  ///        function for all FFTStage objects underneath it.
  void forwardFFT();

  /// \brief Apply the B and C meshes to the transformed density in order to get the frequency
  ///        spread of the non-bonded potential in each system.  This is provided as a "standalone"
  ///        method for applying the Green's function, although other methods are available which
  ///        make use of pre-established abstracts to the relevant objects.
  ///
  /// \param gpu  Details of the GPU that will carry out the calculation
  /// \param sc   Optional energy tracking object
  void applyGreensFunction(const GpuDetails &gpu = null_gpu, ScoreCard *sc = nullptr);

  /// \brief Run the backward FFT operations for all systems.  This will call the eponymous member
  ///        function for all FFTStage objects underneath it.
  void backwardFFT();
  
#ifdef STORMM_USE_HPC
  /// \brief Upload data to the GPU.  This will upload all supporting arrays, such as prefactors
  ///        and frequency grid data offsets, but not the frequency data itself.  The FFT outputs
  ///        are stored at only one memory tier, the level at which the FFTs are performed, and
  ///        may not even be located in data that this object has allocated for itself (the
  ///        underlying PMIGrid will contain the outputs in the case of in-place FFTs).
  void upload();
  
  /// \brief Download data from the GPU.  See above for details of what is downloaded.
  void download();
#endif
  
private:
  int system_count;                 ///< The number of systems in the associated synthesis
  int hpc_work_unit_count;          ///< The number of work units used in HPC implementations for
                                    ///<   applying the Green's function to all systems
  double ewald_coefficient;         ///< The Ewald coefficient, or half the inverse Gaussian
                                    ///<   width used to spread charges on the mesh.
  Hybrid<int> prefactor_offsets;    ///< Offsets for each system's arrays of prefactors in
                                    ///<   convolution computation (B, m, and m')
  Hybrid<double> self_ecorr;        ///< Total self energies of charges in each system in 64-bit
                                    ///<   precision
  Hybrid<double> b_prefactor_a;     ///< "B mesh" prefactors computed for each system along their
                                    ///<   respective unit cell A axes
  Hybrid<double> b_prefactor_b;     ///< "B mesh" prefactors computed for each system along their
                                    ///<   respective unit cell B axes
  Hybrid<double> b_prefactor_c;     ///< "B mesh" prefactors computed for each system along their
                                    ///<   respective unit cell C axes
  Hybrid<double> m_values_a;        ///< Unnormalized m coefficients for each system along their
                                    ///<   respective unit cell A axes
  Hybrid<double> m_values_b;        ///< Unnormalized m coefficients for each system along their
                                    ///<   respective unit cell B axes
  Hybrid<double> m_values_c;        ///< Unnormalized m coefficients for each system along their
                                    ///<   respective unit cell C axes
  Hybrid<double> mshift_values_a;   ///< Unnormalized m coefficients for each system along their
                                    ///<   respective unit cell A axes
  Hybrid<double> mshift_values_b;   ///< Unnormalized m coefficients for each system along their
                                    ///<   respective unit cell B axes
  Hybrid<double> mshift_values_c;   ///< Unnormalized m coefficients for each system along their
                                    ///<   respective unit cell C axes
  Hybrid<double> c_prefactor_a;     ///< "C mesh" prefactors computed for each system along their
                                    ///<   respective unit cell A axes.  These and other "C mesh"
                                    ///<   prefactors are computed only for orthorhombic unit cell
                                    ///<   cases.
  Hybrid<double> c_prefactor_b;     ///< "C mesh" prefactors computed for each system along their
                                    ///<   respective unit cell B axes
  Hybrid<double> c_prefactor_c;     ///< "C mesh" prefactors computed for each system along their
                                    ///<   respective unit cell C axes
  Hybrid<double> double_data;       ///< ARRAY-kind Hybrid targeted by each of the POINTER-kind
                                    ///<   Hybrid<float> objects in the class object

  // Single-precision (32-bit) variants of each of the above arrays
  Hybrid<float> sp_self_ecorr;       ///< Total self energies of charges in each system
  Hybrid<float> sp_b_prefactor_a;    ///< System "B mesh" prefactors computed along the A axis
  Hybrid<float> sp_b_prefactor_b;    ///< System "B mesh" prefactors computed along the B axis
  Hybrid<float> sp_b_prefactor_c;    ///< System "B mesh" prefactors computed along the C axis
  Hybrid<float> sp_m_values_a;       ///< Unnormalized system m coefficients along the A axis
  Hybrid<float> sp_m_values_b;       ///< Unnormalized system m coefficients along the B axis
  Hybrid<float> sp_m_values_c;       ///< Unnormalized system m coefficients along the C axis
  Hybrid<float> sp_mshift_values_a;  ///< Unnormalized, shifted system m coefficients on the A axis
  Hybrid<float> sp_mshift_values_b;  ///< Unnormalized, shifted system m coefficients on the B axis
  Hybrid<float> sp_mshift_values_c;  ///< Unnormalized, shifted system m coefficients on the C axis
  Hybrid<float> sp_c_prefactor_a;    ///< "C mesh" prefactors computed for each system along their
                                     ///<   respective unit cell A axes.  These and other "C mesh"
                                     ///<   prefactors are computed only for orthorhombic unit cell
                                     ///<   cases.
  Hybrid<float> sp_c_prefactor_b;    ///< "C mesh" prefactors computed for each system along their
                                     ///<   respective unit cell B axes
  Hybrid<float> sp_c_prefactor_c;    ///< "C mesh" prefactors computed for each system along their
                                     ///<   respective unit cell C axes
  Hybrid<float> float_data;          ///< ARRAY-kind Hybrid targeted by each of the POINTER-kind
                                     ///<   Hybrid<float> objects in the class object

  // Create plans for FFT operations.  In addition to the cuFFT plans needed for each system's
  // grid, there is an overarching plan created by STORMM itself to optimize the grouping of all
  // FFT operations on the card.
  HybridTargetLevel operating_tier;   ///< The level at which FFTs and the convolution shall be
                                      ///<   performed: the CPU host or the GPU device
  int fft_group_count;                ///< Number of FFT groups collected for C++ computations
                                      ///<   with PocketFFT or HPC calculations on the GPU
  std::vector<int> fft_groups;        ///< A list of collections of systems for related FFT
                                      ///<   operations carried out by the CPU using the
                                      ///<   C++-compatible PocketFFT or by the GPU using the
                                      ///<   appropriate HPC FFT package.  This information is
                                      ///<   held only on the CPU, for purposes of querying the
                                      ///<   details of an FFT call such as in backtracing.
  std::vector<int> fft_group_bounds;  ///< Bounds array for fft_groups.  This information is held
                                      ///<   on the CPU.
  Hybrid<double2> dp_frequency_data;  ///< Double-precision data allocated by the object in
                                      ///<   support of out-of-place FFTs set by the underlying
                                      ///<   particle-mesh interaction grid (PMIGrid) object
  Hybrid<float2> sp_frequency_data;   ///< Single-precision data allocated by the object in
                                      ///<   support of out-of-place FFTs set by the underlying
                                      ///<   particle-mesh interaction grid (PMIGrid) object
  Hybrid<uint> frequency_offsets;     ///< Starting points for the frequency-space grids of each
                                      ///<   system.  These offsets are calculated based on
                                      ///<   information from the underlying PMIGrid object, if
                                      ///<   in-place FFTs are called for, or to provide spacings
                                      ///<   in a padded arrangement of all such grids in the
                                      ///<   object's private allocations (dp_frequency_data or
                                      ///<   sp_frequency_data).

  /// The FFT operations themselves are handled by an array of nested class objects (an example of
  /// composition, rather than inheritance, in C++).  Each object can be queried for details about
  /// its specific procedures.
  std::vector<FFTStage> fft_operations;

  /// Work units are once again used to load-balance applications of the Green's function across
  /// multiple systems, which may be of disparate sizes.  In these simple work unit instructions,
  /// the "x" member of the tuple indicates the system upon which the work unit will focus.  The
  /// "y" member indicates the position in the system's grid where the work will start (this is a
  /// relative index, and must be offset by the starting point of the grid, in frequency space,
  /// on the contiguous array containing all frequency space representations).  The "z" member of
  /// the tuple contains the upper limit of grid elements that the work unit will span.  The "w"
  /// member is unused.
  Hybrid<uint4> work_units;
  
  // Pointers to associated objects
  PMIGrid *pmig_ptr;                ///< Pointer to the associated Particle-Mesh Interaction Grid
  AtomGraphSynthesis *poly_ag_ptr;  ///< Pointer to the associated topology synthesis

  /// \brief Allocate basic memory for the object.  This will not allocate fft_groups,
  ///        fft_group_bounds, or the frequency data and associated offsets.  What it does allocate
  ///        are the ARRAY-kind Hybrid objects that will be targeted by various POINTER-kind Hybrid
  ///        objects
  void allocate();

  /// \brief Compute the self energies of each system given their particle densities
  ///        (electrostatic or dispersion charge) and the chosen splitting constant.
  void computeSystemSelfEnergies();

  /// \brief Detect groups of systems that have the same grid dimensions.  These systems can all
  ///        use the same FFT staging object for C++ computations or, if they are continguous in
  ///        the list of systems, be batched for a single call to the HPC API.
  void makeFFTGroups(const GpuDetails &gpu = null_gpu);

  /// \brief Create a series of work units to load-balance work on the available GPU.
  void makeWorkUnits(const GpuDetails &gpu = null_gpu);
};

/// \brief Apply the Particle-Mesh Ewald Green's Function in the form of the "B" and "C" meshes to
///        a transformed density in order to get the frequency spread of the non-bonded potential
///        in each of a synthesis of systems.
///
/// Overloaded:
///   - Templated functions available for CPU-bound work
///   - GPU-enabled variants are enumerated along the non-bonded calculation types
///   - Provide abstracts or original objects for the convolution manager, particle-mesh
///     interaction grids, and system coordinate sets (unit cell transformation matrices)
///
/// \param cvolw  Abstract of the object itself, taken with the appropriate precision
/// \param pmigr  Abstract of the particle-mesh interactions for each system
/// \param pssb   Abbreviated abstract of the coordinate synthesis containing the transformation
///               matrices for each system
/// \param scw    Abstract to an energy tracking object.  Unlike earlier functions that compute
///               potential energy, this abstract (not the underlying energy tracking object) is
///               passed to functions that launch GPU kernels as well as function that run on the
///               CPU, despite some tedium in entering the results into the storage arrays.
/// \param gpu    Details of the GPU that will carry out the calculation
/// \{
template <typename T, typename T2>
std::vector<double> pmeGreensFunction(ConvolutionWriter<T, T2> *cvolw,
                                      const PsSynthesisBorders &pssb, const PMIGridReader &pmigr,
                                      ScoreCardWriter *scw = nullptr);

#ifdef STORMM_USE_HPC
void pmeGreensFunction(ConvolutionWriter<double, double2> *cvolw, const PsSynthesisBorders &pssb,
                       const PMIGridReader &pmigr, const GpuDetails &gpu,
                       ScoreCardWriter *scw = nullptr);

void pmeGreensFunction(ConvolutionWriter<float, float2> *cvolw, const PsSynthesisBorders &pssb,
                       const PMIGridReader &pmigr, const GpuDetails &gpu,
                       ScoreCardWriter *scw = nullptr);
#endif
/// \}

/// \brief Complete the entire series of convolutions for a synthesis of systems.
///
/// Overloaded:
///   - Templated function operating on the CPU
///   - Provide abstracts containing pre-calculated constants for the convolutions, unit cell
///     measurements, and the particle-mesh interaction grids
///   - Provide original class objects for all of the above
///
/// \param colvw     Abstract of the convolution management object
/// \param colv      Convolution management with pre-calculated constants for applying the
///                  Particle-Mesh Ewald Green's function, FFT plans, the Ewald coefficient
///                  describing the spread of Gaussian charges, and indicators about whether to
///                  carry out calculations on the CPU host or GPU device
/// \param pssb      Abstract of the coordinate synthesis containing measurements of unit cells,
///                  including the transformation matrices needed to normalize various aspects of
///                  the convolution and compute unit cell volumes
/// \param poly_ps   Synthesis of multiple systems' coordinates (atomic positions, velocities, and
///                  force accumulators)
/// \param pmigr     Read-only abstract of the particle-mesh interaction grid.  The convolution
///                  manager contains other pointers to the particle-mesh interaction grid's data,
///                  in mutable format.  This is needed for information such as the grid dimensions
///                  and FFT mode.
/// \param pmig      Read-only object providing the particle-mesh interaction grids.  The
///                  convolution manager contains other pointers to the particle-mesh interaction
///                  grid's data, in mutable format.  This is needed for information such as the
///                  grid dimensions and FFT mode.
/// \param scw       Writeable abstract of the ScoreCard energy tacking object
/// \param sc        Energy tracking object
/// \param gpu       Details of the GPU that may carry out the calculations
/// \{
template <typename T, typename T2>
std::vector<double> applyConvolution(ConvolutionWriter<float, float2> *cvolw,
                                     const PsSynthesisBorders &pssb,
                                     const PMIGridReader &pmigr, ScoreCardWriter *scw = nullptr);

std::vector<double> applyConvolution(ConvolutionManager *cvol, const PhaseSpaceSynthesis *poly_ps,
                                     const PMIGrid *pmig, ScoreCard *sc = nullptr);

std::vector<double> applyConvolution(ConvolutionManager *cvol, const PhaseSpaceSynthesis &poly_ps,
                                     const PMIGrid &pmig, ScoreCard *sc = nullptr);
#ifdef STORMM_USE_HPC
void applyConvolution(ConvolutionWriter<double, double2> *cvolw,
                      const PsSynthesisBorders &pssb, const PMIGridReader &pmigr,
                      const GpuDetails &gpu, ScoreCardWriter *scw = nullptr);

void applyConvolution(ConvolutionWriter<float, float2> *cvolw,
                      const PsSynthesisBorders &pssb, const PMIGridReader &pmigr,
                      const GpuDetails &gpu, ScoreCardWriter *scw = nullptr);

void applyConvolution(ConvolutionManager *cvol, const PhaseSpaceSynthesis *poly_ps,
                      const PMIGrid *pmig, const GpuDetails &gpu, ScoreCard *sc = nullptr);

void applyConvolution(ConvolutionManager *cvol, const PhaseSpaceSynthesis &poly_ps,
                      const PMIGrid &pmig, const GpuDetails &gpu, ScoreCard *sc = nullptr);
#endif
/// \}
  
} // namespace energy
} // namespace stormm

#include "convolution_manager.tpp"

#endif
