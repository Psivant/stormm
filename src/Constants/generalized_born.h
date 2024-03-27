// -*-c++-*-
#ifndef STORMM_GENERALIZED_BORN_H
#define STORMM_GENERALIZED_BORN_H

#include "copyright.h"
#include "Accelerator/hybrid.h"

namespace stormm {
namespace generalized_born_defaults {

using card::Hybrid;
using card::HybridTargetLevel;

/// \brief Scaling factors to be used in computing fixed-precision Generalized Born effective radii
/// \{
constexpr double gb_radii_scale_lf = 16277216.0;
constexpr float  gb_radii_scale_f  = (float)gb_radii_scale_lf;
constexpr double inverse_gb_radii_scale_lf = 1.0 / gb_radii_scale_lf;
constexpr float  inverse_gb_radii_scale_f  = (float)1.0 / gb_radii_scale_f;
/// \}

/// \brief Salt concentration dependence for Generalized Born kappa scaling parameter
constexpr double default_salt_kappa_dependence = 0.10806;
  
/// \brief Double-precision Generalized Born Taylor series expansion coefficients
/// \{
constexpr double gb_taylor_a_lf  = 0.33333333333333333333;
constexpr double gb_taylor_b_lf  = 0.4;
constexpr double gb_taylor_c_lf  = 0.42857142857142857143;
constexpr double gb_taylor_d_lf  = 0.44444444444444444444;
constexpr double gb_taylor_dd_lf = 0.45454545454545454545;
constexpr double gb_taylor_e_lf  = 1.33333333333333333333;
constexpr double gb_taylor_f_lf  = 2.4;
constexpr double gb_taylor_g_lf  = 3.42857142857142857143;
constexpr double gb_taylor_h_lf  = 4.44444444444444444444;
constexpr double gb_taylor_hh_lf = 5.45454545454545454545;
/// \}

/// \brief Single-precision Generalized Born Taylor series expansion coefficients
/// \{
constexpr float gb_taylor_a_f  = gb_taylor_a_lf;
constexpr float gb_taylor_b_f  = gb_taylor_b_lf;
constexpr float gb_taylor_c_f  = gb_taylor_c_lf;
constexpr float gb_taylor_d_f  = gb_taylor_d_lf;
constexpr float gb_taylor_dd_f = gb_taylor_dd_lf;
constexpr float gb_taylor_e_f  = gb_taylor_e_lf;
constexpr float gb_taylor_f_f  = gb_taylor_f_lf;
constexpr float gb_taylor_g_f  = gb_taylor_g_lf;
constexpr float gb_taylor_h_f  = gb_taylor_h_lf;
constexpr float gb_taylor_hh_f = gb_taylor_hh_lf;
/// \}

/// \brief Onufriev / Bashford / Case Generalized Born model I alpha, beta, and gamma values
///        (these are universal, applied across all atoms regardless of their element)
/// \{
constexpr double gb_obc_i_alpha = 0.8;
constexpr double gb_obc_i_beta  = 0.0;
constexpr double gb_obc_i_gamma = 2.909125;
/// \}

/// \brief Onufriev / Bashford / Case Generalized Born model II alpha, beta, and gamma values
///        (these are universal, applied across all atoms regardless of their element)
/// \{
constexpr double gb_obc_ii_alpha = 1.0;
constexpr double gb_obc_ii_beta  = 0.8;
constexpr double gb_obc_ii_gamma = 4.85;
/// \}

/// \brief "Neck" Generalized Born model I alpha, beta, and gamma parameters, plus some elements'
///        screening factors
///        
/// \{
constexpr double gb_neck_i_alpha    = 1.09511284;
constexpr double gb_neck_i_beta     = 1.90792938;
constexpr double gb_neck_i_gamma    = 2.50798245;
constexpr double gb_neck_i_screen_h = 1.09085413633;
constexpr double gb_neck_i_screen_c = 0.484353823306;
constexpr double gb_neck_i_screen_n = 0.700147318409;
constexpr double gb_neck_i_screen_o = 1.06557401132;
constexpr double gb_neck_i_screen_s = 0.602256336067;
constexpr double gb_neck_i_screen_default = 0.5;
/// \}

/// \brief "Neck" Generalized Born model II screening factors and element-specific alpha, beta,
///        and gamma parameters
/// \{
constexpr double gb_neck_ii_screen_h =  1.425952;
constexpr double gb_neck_ii_screen_c =  1.058554;
constexpr double gb_neck_ii_screen_n =  0.733599;
constexpr double gb_neck_ii_screen_o =  1.061039;
constexpr double gb_neck_ii_screen_s = -0.703469;
constexpr double gb_neck_ii_screen_p =  0.500000;
constexpr double gb_neck_ii_alpha_h  =  0.788440;
constexpr double gb_neck_ii_beta_h   =  0.798699;
constexpr double gb_neck_ii_gamma_h  =  0.437334;
constexpr double gb_neck_ii_alpha_c  =  0.733756;
constexpr double gb_neck_ii_beta_c   =  0.506378;
constexpr double gb_neck_ii_gamma_c  =  0.205844;
constexpr double gb_neck_ii_alpha_n  =  0.503364;
constexpr double gb_neck_ii_beta_n   =  0.316828;
constexpr double gb_neck_ii_gamma_n  =  0.192915;
constexpr double gb_neck_ii_alpha_os =  0.867814;
constexpr double gb_neck_ii_beta_os  =  0.876635;
constexpr double gb_neck_ii_gamma_os =  0.387882;
constexpr double gb_neck_ii_alpha_p  =  1.000000;
constexpr double gb_neck_ii_beta_p   =  0.800000;
constexpr double gb_neck_ii_gamma_p  =  4.850000;
/// \}

/// \brief The offset used to translate Poisson-Boltzmann atomic radii into baseline GB radii.
/// \{
constexpr double default_gb_radii_offset = 0.09;
constexpr double default_neck_ii_gb_radii_offset = 0.195141;
/// \}

/// \brief Other "neck" GB parameters.  The neck cutoff has a physical interpretation which would
///        put it near 2.8 Angstroms, the diameter of a water molecule, but it is set higher (6.8)
///        to reduce the discontinuity at the cutoff.  gb_kscale is a "Kappa" scaling parameter.
constexpr double default_gb_neck_cut      =  6.800000;
constexpr double default_gb_kscale        =  0.730000;
constexpr double default_gb_neck_scale    = 0.361825;
constexpr double default_gb_neck_ii_scale = 0.826836;
constexpr float default_gb_neck_cut_f      = default_gb_neck_cut;
constexpr float default_gb_kscale_f        = default_gb_kscale;
constexpr float default_gb_neck_scale_f    = default_gb_neck_scale;
constexpr float default_gb_neck_ii_scale_f = default_gb_neck_ii_scale;
/// \}

/// \brief Abstract for the NeckGeneralizedBornTable object, in single- or double-precision.
template <typename T> struct NeckGeneralizedBornKit {

  /// \brief Constructor takes a simple list of the relevant constants and pointers
  explicit NeckGeneralizedBornKit(int table_size, T neck_cut_in, T kscale_in,
                                  const T* max_separation_in, const T* max_value_in);

  const int table_size;    ///< Size of each table (they are both assumed to be square)
  const T neck_cut;        ///< The Generalized Born neck cutoff (see above for the default)
  const T kscale;          ///< Generalized Born Kappa scaling parameter (see the default above).
                           ///<   This parameter is carried through to this object to be loaded
                           ///<   from the CPU cache or HPC constants cache, rather than taken at
                           ///<   its default setting, to allow flexibility in the GB
                           ///<   implementation.  In contrast, the Generalzied Born Taylor series
                           ///<   coefficients are not set by the user and left as constexpr.
  const T* max_separation; ///< Maximum separations of spheres resulting in neck corrections
  const T* max_value;      ///< Maximum values of the neck correction in each case
};

/// \brief Object to hold a complex array of constants referenced by various GB calculations using
///        the "neck" formalism for the union of spheres.  This is an object of its own to
///        encapsulate the underlying Hybrid data structures.
class NeckGeneralizedBornTable {
public:

  /// The basic constructor hard-codes the large table of constants into the relevant Hybrid
  /// objects, with no arguments needed.  An overload of the constructor take custom arrays
  ///
  /// \param table_size_in  Length of the tables (they are assumed to be square)
  /// \param neck_cut_in    The GB neck cutoff distance (see above)
  /// \param kscale_in      The GB Kappa scaling parameter (see above)
  /// \param max_sep_in     Data for maximum separations (must be table_size_in * table_size_in)
  /// \param max_val_in     Data for maximum values (must be table_size_in * table_size_in)
  /// \{
  NeckGeneralizedBornTable();
  NeckGeneralizedBornTable(int table_size_in, double neck_cut_in, double kscale_in,
                           const std::vector<double> &max_sep_in,
                           const std::vector<double> &max_val_in);
  /// \}

  /// \brief Get the table size
  int size() const;

  /// \brief Get the maximum separation of neck placement for a pair of Generalized Born atom
  ///        types.  This reads from the host only, as there is no way to modify the device-side
  ///        data without first modifying it on the host.
  ///
  /// \param i_type  Generalized Born type of the first atom
  /// \param j_type  Generalized Born type of the second atom
  double getMaxSeparation(int i_type, int j_type) const;

  /// \brief Get the maximum value of the neck factor for a pair of Generalized Born atom types.
  ///        This reads from the host only, as there is no way to modify the device-side data
  ///        without first modifying it on the host.
  ///
  /// \param i_type  Generalized Born type of the first atom
  /// \param j_type  Generalized Born type of the second atom
  double getMaxValue(int i_type, int j_type) const;

  /// \brief Get a collection of double-precision pointers for this object
  ///
  /// \param tier  Level at which to get the pointers (CPU or GPU, HOST or DEVICE)
  const NeckGeneralizedBornKit<double>
  dpData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a collection of single-precision pointers for this object
  ///
  /// \param tier  Level at which to get the pointers (CPU or GPU, HOST or DEVICE)
  const NeckGeneralizedBornKit<float>
  spData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get a const pointer to the object itself, useful if a pointer is needed after passing
  ///        he object by const reference.
  const NeckGeneralizedBornTable* getSelfPointer() const;
  
private:
  int table_size;                        ///< Size of the neck GB lookup tables (if applicable)
  double neck_cut;                       ///< Neck GB cutoff (if applicable)
  double kscale;                         ///< GB Kappa scaling parameter
  Hybrid<double> neck_max_separation;    ///< Maximum separations between particles of given sizes
                                         ///<   for which the neck function has a value
  Hybrid<double> neck_max_value;         ///< Maximum values of the neck function for particles of
                                         ///<    two given, indexed sizes
  Hybrid<float> sp_neck_max_separation;  ///< Single-precision version of the eponymous array
  Hybrid<float> sp_neck_max_value;       ///< Single-precision version of the eponymous array
};

} // namespace generalized_born_defaults
} // namespace stormm

#include "generalized_born.tpp"

#endif
