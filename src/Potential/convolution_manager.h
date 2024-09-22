// -*-c++-*-
#ifndef STORMM_CONVOLUTION_MANAGER_H
#define STORMM_CONVOLUTION_MANAGER_H

#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "pmigrid.h"

namespace stormm {
namespace energy {

using constants::PrecisionModel;
using card::Hybrid;
using card::HybridTargetLevel;
using synthesis::AtomGraphSynthesis;

/// \brief Any abstract of the convolution manager is read-only.
template <typename T> struct ConvolutionKit {

  /// \brief As with other abstracts, the convolution kit is constructed with a series of values
  ///        for each associated pointer and critical constant.  See Essmann's 1995 paper (below)
  ///        for  the inspiration to some of the nomenclature.
  ConvolutionKit(int system_count_in, T ew_coeff_in, const int* sys_offsets_in,
                 const T* self_ecorr_in, const T* bmesh_a_in, const T* bmesh_b_in,
                 const T* bmesh_c_in, const T* mval_a_in, const T* mval_b_in, const T* mval_c_in,
                 const T* msval_a_in, const T* msval_b_in, const T* msval_c_in,
                 const T* cmesh_a_in, const T* cmesh_b_in, const T* cmesh_c_in);

  /// \brief Like other abstracts, the presence of const members implicitly forbids the copy and
  ///        move assignment operators.  Copy and move constructors may be taken in their deafult
  ///        forms.
  ///
  /// \param original  The original object to copy or move
  /// \{
  ConvolutionKit(const ConvolutionKit &original) = default;
  ConvolutionKit(ConvolutionKit &&original) = default;
  /// \}

  const int system_count;  ///< The number of indepedent systems
  const T ew_coeff;        ///< The Ewald coefficien in use by the splitting function
  const int* sys_offsets;  ///< Offsets for grid elements and prefactors in each system.  The
                           ///<   maximum of all grid dimensions for a particular system, rounded
                           ///<   to the nearest multiple of the warp size, determines each
                           ///<   successive system offset in the arrays for prefactors and M
                           ///<   values along any of the A, B, or C unit cell dimensions.
  const T* self_ecorr;     ///< Self-correlation energy of the charges (electrostatic or otherwise)
                           ///<   in each system.
  const T* bmesh_a;        ///< "B" mesh prefactors for each system along the unit cell A axis
  const T* bmesh_b;        ///< "B" mesh prefactors for each system along the unit cell B axis
  const T* bmesh_c;        ///< "B" mesh prefactors for each system along the unit cell C axis
  const T* mval_a;         ///< "M" values for each system along the unit cell A axis
  const T* mval_b;         ///< "M" values for each system along the unit cell B axis
  const T* mval_c;         ///< "M" values for each system along the unit cell C axis
  const T* msval_a;        ///< Shifted "M" values for each system along the unit cell A axis
  const T* msval_b;        ///< Shifted "M" values for each system along the unit cell B axis
  const T* msval_c;        ///< Shifted "M" values for each system along the unit cell C axis
  const T* cmesh_a;        ///< "C" mesh prefactors for each system along the unit cell A axis.
                           ///<   These are only computed if all systems' unit cells are
                           ///<   orthorhombic.
  const T* cmesh_b;        ///< "C" mesh prefactors for each system along the unit cell B axis.
  const T* cmesh_c;        ///< "C" mesh prefactors for each system along the unit cell C axis.
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
  ConvolutionManager(const PMIGrid *pmig_in, const double ewald_coefficient_in);
  ConvolutionManager(const PMIGrid &pmig_in, const double ewald_coefficient_in);
  /// \}

  /// \brief The copy and move constructors, as well as copy andmove assignemnt operators, must all
  ///        be given explicit definitions due to the presence of POINTER-kind Hybrid objects.
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

#ifdef STORMM_USE_HPC
  /// \brief Upload data to the GPU
  void upload();
  
  /// \brief Download data from the GPU
  void download();
#endif
  
private:
  int system_count;                 ///< The number of systems in the associated synthesis
  double ewald_coefficient;         ///< The Ewald coefficient, or half the inverse Gaussian
                                    ///<   width used to spread charges on the mesh.
  Hybrid<int> system_offsets;       ///< Offsets for each system's arrays of prefators in
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

  // Pointers to associated objects
  PMIGrid *pmig_ptr;                ///< Pointer to the associated Particle-Mesh Interaction Grid
  AtomGraphSynthesis *poly_ag_ptr;  ///< Pointer to the associated topology synthesis

  /// \brief Allocate memory for the object.
  void allocate();

  /// \brief Compute the self energies of each system given their particle densities
  ///        (electrostatic or dispersion charge) and the chosen splitting constant.
  void computeSystemSelfEnergies();
};

} // namespace energy
} // namespace stormm

#include "convolution_manager.tpp"

#endif
