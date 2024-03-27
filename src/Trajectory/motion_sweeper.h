// -*-c++-*-
#ifndef STORMM_MOTION_SWEEPER_H
#define STORMM_MOTION_SWEEPER_H

#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/fixed_precision.h"
#include "DataTypes/stormm_vector_types.h"
#include "Synthesis/phasespace_synthesis.h"

namespace stormm {
namespace trajectory {

using card::HybridTargetLevel;
using numerics::default_momentum_scale_bits;
using numerics::default_com_scale_bits;
using numerics::default_inertia_scale_bits;
using synthesis::PhaseSpaceSynthesis;

/// \brief Abstract of the MotionSweeper class, below, with non-const pointers to modifiable data.
struct MotionSweepWriter {
public:

  /// \brief As with other abstracts, the constructor takes critical constants and pointers to
  ///        initialize each of its member variables.
  MotionSweepWriter(int nwu_in, const int4* work_units_in, const double* total_mass_in,
                    double com_scale_in, double mv_scale_in, double inrt_scale_in, llint* xcom_in,
                    llint* ycom_in, llint* zcom_in, llint* xcom_nxt_in, llint* ycom_nxt_in,
                    llint* zcom_nxt_in, int* xcom_ovrf_in, int* ycom_ovrf_in, int* zcom_ovrf_in,
                    int* xcom_nxt_ovrf_in, int* ycom_nxt_ovrf_in, int* zcom_nxt_ovrf_in,
                    llint* xmv_in, llint* ymv_in, llint* zmv_in, llint* xmv_nxt_in,
                    llint* ymv_nxt_in, llint* zmv_nxt_in, int* xmv_ovrf_in, int* ymv_ovrf_in,
                    int* zmv_ovrf_in, int* xmv_nxt_ovrf_in, int* ymv_nxt_ovrf_in,
                    int* zmv_nxt_ovrf_in, llint* rxmv_in, llint* rymv_in, llint* rzmv_in,
                    llint* rxmv_nxt_in, llint* rymv_nxt_in, llint* rzmv_nxt_in, int* rxmv_ovrf_in,
                    int* rymv_ovrf_in, int* rzmv_ovrf_in, int* rxmv_nxt_ovrf_in,
                    int* rymv_nxt_ovrf_in, int* rzmv_nxt_ovrf_in, llint* inrt_in,
                    llint* inrt_nxt_in, int* inrt_ovrf_in, int* inrt_nxt_ovrf_in);

  /// \brief As with other abstracts, the presence of any const members invalidates the copy and
  ///        move assignment operators (unless using special features of more recent C++
  ///        standards).  The default copy and move constructors may be taken.
  ///
  /// \param original  The original abstract to copy or move
  /// \{
  MotionSweepWriter(const MotionSweepWriter &original) = default;
  MotionSweepWriter(MotionSweepWriter &&original) = default;
  /// \}
  
  const int nwu;             ///< The number of work units
  const int4* work_units;    ///< Work units spanning all systems.  Members of each tuple contain
                             ///<   the index of the starting atom for the work unit to accumulate,
                             ///<   upper limit of atom indices, and the system index to which the
                             ///<   accumulation pertains in the "x", "y", and "z" members.  The
                             ///<   "w" member of each tuple indicates other tasks assigned to any
                             ///<   given work unit.
  const double* total_mass;  ///< Total masses of each system in the synthesis
  const double com_scale;    ///< Scaling factor taking quantities in units of g/mol - Angstroms
                             ///<   into the fixed-precision format of the momentum accumulators
  const double mv_scale;     ///< Scaling factor taking quantities in units of g/mol - Angstroms
                             ///<   per fs into the fixed-precision format of the momentum
                             ///<   accumulators
  const double inrt_scale;   ///< Scaling factor taking quantities in units of g/mol - squared
                             ///<   Angstroms into the fixed-precision format of the inertial
                             ///<   tensor accumulators

  // Each of the following arrays has a dimension equal to the number of systems in the associated
  // synthesis.  The number of systems is not stored in this abstract as the correct index of each
  // array can be found in the work unit itself.
  llint* xcom;         ///< Accumulators for center of mass along the X dimension
  llint* ycom;         ///< Accumulators for center of mass along the Y dimension
  llint* zcom;         ///< Accumulators for center of mass along the Z dimension
  llint* xcom_nxt;     ///< Accumulators for COM along X in the next iteration 
  llint* ycom_nxt;     ///< Accumulators for COM along Y in the next iteration 
  llint* zcom_nxt;     ///< Accumulators for COM along Z in the next iteration 
  int* xcom_ovrf;      ///< Overflow bits for center of mass in the X dimension
  int* ycom_ovrf;      ///< Overflow bits for center of mass in the Y dimension
  int* zcom_ovrf;      ///< Overflow bits for center of mass in the Z dimension
  int* xcom_nxt_ovrf;  ///< Overflow bits for COM along X in the next iteration 
  int* ycom_nxt_ovrf;  ///< Overflow bits for COM along Y in the next iteration 
  int* zcom_nxt_ovrf;  ///< Overflow bits for COM along Z in the next iteration 
  llint* xmv;          ///< Accumulators for net momentum along the X dimension
  llint* ymv;          ///< Accumulators for net momentum along the Y dimension
  llint* zmv;          ///< Accumulators for net momentum along the Z dimension
  llint* xmv_nxt;      ///< Accumulators for net momentum in X in the next iteration
  llint* ymv_nxt;      ///< Accumulators for net momentum in Y in the next iteration
  llint* zmv_nxt;      ///< Accumulators for net momentum in Z in the next iteration
  int* xmv_ovrf;       ///< Overflow bits for net momentum in the X dimension
  int* ymv_ovrf;       ///< Overflow bits for net momentum in the Y dimension
  int* zmv_ovrf;       ///< Overflow bits for net momentum in the Z dimension
  int* xmv_nxt_ovrf;   ///< Overflow bits for net momentum in X in the next iteration
  int* ymv_nxt_ovrf;   ///< Overflow bits for net momentum in Y in the next iteration
  int* zmv_nxt_ovrf;   ///< Overflow bits for net momentum in Z in the next iteration
  llint* rxmv;         ///< Accumulators for angular momentum about the X axis
  llint* rymv;         ///< Accumulators for angular momentum about the Y axis
  llint* rzmv;         ///< Accumulators for angular momentum about the Z axis
  llint* rxmv_nxt;     ///< Angular momentum about the X axis in the next iteration
  llint* rymv_nxt;     ///< Angular momentum about the Y axis in the next iteration
  llint* rzmv_nxt;     ///< Angular momentum about the Z axis in the next iteration
  int* rxmv_ovrf;      ///< Overflow bits for angular momentum about the X axis
  int* rymv_ovrf;      ///< Overflow bits for angular momentum about the Y axis
  int* rzmv_ovrf;      ///< Overflow bits for angular momentum about the Z axis
  int* rxmv_nxt_ovrf;  ///< Overflow bits for X axis angular momentum in the next iteration
  int* rymv_nxt_ovrf;  ///< Overflow bits for Y axis angular momentum in the next iteration
  int* rzmv_nxt_ovrf;  ///< Overflow bits for Z axis angular momentum in the next iteration
  llint* inrt;         ///< Accumulators for the inertial tensor.  Slots for each system are stored
                       ///<   back-to-back, i.e. the slots for system i are (6 * i, (6 * i) + 1,
                       ///<   ..., (6 * i) + 5, and slots for system i + 1 are (6 * i) + 6,
                       ///<   (6 * i) + 7, ..., (6 * i) + 11.
  llint* inrt_nxt;     ///< Accumulators for the inertial tensor in the next iteration
  int* inrt_ovrf;      ///< Overflow bits for the inertial tensor
  int* inrt_nxt_ovrf;  ///< Overflow bits for the inertial tensor in the next iteration
};

/// \brief Read-only abstract of the MotionSweeper class, below, with const pointers to
///        non-modifiable data.
struct MotionSweepReader {
public:

  /// \brief As with other abstracts, the constructor takes critical constants and pointers to
  ///        initialize each of its member variables.  Constructors to covert the writeable
  ///        abstract to the read-only abstract are also provided.
  /// \{
  MotionSweepReader(int nwu_in, const int4* work_units_in, const double* total_mass_in,
                    double com_scale_in, double mv_scale_in, double inrt_scale,
                    const llint* xcom_in, const llint* ycom_in, const llint* zcom_in,
                    const llint* xcom_nxt_in, const llint* ycom_nxt_in, const llint* zcom_nxt_in,
                    const int* xcom_ovrf_in, const int* ycom_ovrf_in, const int* zcom_ovrf_in,
                    const int* xcom_nxt_ovrf_in, const int* ycom_nxt_ovrf_in,
                    const int* zcom_nxt_ovrf_in, const llint* xmv_in, const llint* ymv_in,
                    const llint* zmv_in, const llint* xmv_nxt_in, const llint* ymv_nxt_in,
                    const llint* zmv_nxt_in, const int* xmv_ovrf_in, const int* ymv_ovrf_in,
                    const int* zmv_ovrf_in, const int* xmv_nxt_ovrf_in, const int* ymv_nxt_ovrf_in,
                    const int* zmv_nxt_ovrf_in, const llint* rxmv_in, const llint* rymv_in,
                    const llint* rzmv_in, const llint* rxmv_nxt_in, const llint* rymv_nxt_in,
                    const llint* rzmv_nxt_in, const int* rxmv_ovrf_in, const int* rymv_ovrf_in,
                    const int* rzmv_ovrf_in, const int* rxmv_nxt_ovrf_in,
                    const int* rymv_nxt_ovrf_in, const int* rzmv_nxt_ovrf_in, const llint* inrt_in,
                    const llint* inrt_nxt_in, const int* inrt_ovrf_in,
                    const int* inrt_nxt_ovrf_in);

  MotionSweepReader(const MotionSweepWriter *w);
  
  MotionSweepReader(const MotionSweepWriter &w);
  /// \}

  /// \brief As with other abstracts, the presence of any const members invalidates the copy and
  ///        move assignment operators (unless using special features of more recent C++
  ///        standards).  The default copy and move constructors may be taken.
  ///
  /// \param original  The original abstract to copy or move
  /// \{
  MotionSweepReader(const MotionSweepReader &original) = default;
  MotionSweepReader(MotionSweepReader &&original) = default;
  /// \}
  
  const int nwu;             ///< The number of work units
  const int4* work_units;    ///< Work units spanning all systems.  Members of each tuple contain
                             ///<   the index of the starting atom for the work unit to accumulate,
                             ///<   upper limit of atom indices, and the system index to which the
                             ///<   accumulation pertains in the "x", "y", and "z" members.  The
                             ///<   "w" member of each tuple indicates other tasks assigned to any
                             ///<   given work unit.
  const double* total_mass;  ///< Total masses of each system in the synthesis
  const double com_scale;    ///< Scaling factor taking quantities in units of g/mol - Angstroms
                             ///<   into the fixed-precision format of the momentum accumulators
  const double mv_scale;     ///< Scaling factor taking quantities in units of g/mol - Angstroms
                             ///<   per fs into the fixed-precision format of the momentum
                             ///<   accumulators
  const double inrt_scale;   ///< Scaling factor taking quantities in units of g/mol - squared
                             ///<   Angstroms into the fixed-precision format of the inertial
                             ///<   tensor accumulators

  // Each of the following arrays has a dimension equal to the number of systems in the associated
  // synthesis.  The number of systems is not stored in this abstract as the correct index of each
  // array can be found in the work unit itself.
  const llint* xcom;         ///< Accumulators for center of mass along the X dimension
  const llint* ycom;         ///< Accumulators for center of mass along the Y dimension
  const llint* zcom;         ///< Accumulators for center of mass along the Z dimension
  const llint* xcom_nxt;     ///< Accumulators for COM along X in the next iteration 
  const llint* ycom_nxt;     ///< Accumulators for COM along Y in the next iteration 
  const llint* zcom_nxt;     ///< Accumulators for COM along Z in the next iteration 
  const int* xcom_ovrf;      ///< Overflow bits for center of mass in the X dimension
  const int* ycom_ovrf;      ///< Overflow bits for center of mass in the Y dimension
  const int* zcom_ovrf;      ///< Overflow bits for center of mass in the Z dimension
  const int* xcom_nxt_ovrf;  ///< Overflow bits for COM along X in the next iteration 
  const int* ycom_nxt_ovrf;  ///< Overflow bits for COM along Y in the next iteration 
  const int* zcom_nxt_ovrf;  ///< Overflow bits for COM along Z in the next iteration 
  const llint* xmv;          ///< Accumulators for net momentum along the X dimension
  const llint* ymv;          ///< Accumulators for net momentum along the Y dimension
  const llint* zmv;          ///< Accumulators for net momentum along the Z dimension
  const llint* xmv_nxt;      ///< Accumulators for net momentum in X in the next iteration
  const llint* ymv_nxt;      ///< Accumulators for net momentum in Y in the next iteration
  const llint* zmv_nxt;      ///< Accumulators for net momentum in Z in the next iteration
  const int* xmv_ovrf;       ///< Overflow bits for net momentum in the X dimension
  const int* ymv_ovrf;       ///< Overflow bits for net momentum in the Y dimension
  const int* zmv_ovrf;       ///< Overflow bits for net momentum in the Z dimension
  const int* xmv_nxt_ovrf;   ///< Overflow bits for net momentum in X in the next iteration
  const int* ymv_nxt_ovrf;   ///< Overflow bits for net momentum in Y in the next iteration
  const int* zmv_nxt_ovrf;   ///< Overflow bits for net momentum in Z in the next iteration
  const llint* rxmv;         ///< Accumulators for angular momentum about the X axis
  const llint* rymv;         ///< Accumulators for angular momentum about the Y axis
  const llint* rzmv;         ///< Accumulators for angular momentum about the Z axis
  const llint* rxmv_nxt;     ///< Angular momentum about the X axis in the next iteration
  const llint* rymv_nxt;     ///< Angular momentum about the Y axis in the next iteration
  const llint* rzmv_nxt;     ///< Angular momentum about the Z axis in the next iteration
  const int* rxmv_ovrf;      ///< Overflow bits for angular momentum about the X axis
  const int* rymv_ovrf;      ///< Overflow bits for angular momentum about the Y axis
  const int* rzmv_ovrf;      ///< Overflow bits for angular momentum about the Z axis
  const int* rxmv_nxt_ovrf;  ///< Overflow bits for X axis angular momentum in the next iteration
  const int* rymv_nxt_ovrf;  ///< Overflow bits for Y axis angular momentum in the next iteration
  const int* rzmv_nxt_ovrf;  ///< Overflow bits for Z axis angular momentum in the next iteration
  const llint* inrt;         ///< Accumulators for the intertial tensor.  Slots for each system are
                             ///<   stored back-to-back, i.e. the slots for system i are (6 * i,
                             ///<   (6 * i) + 1, ..., (6 * i) + 5, and slots for system i + 1 are
                             ///<   (6 * i) + 6, (6 * i) + 7, ..., (6 * i) + 11.
  const llint* inrt_nxt;     ///< Accumulators for the intertial tensor in the next iteration
  const int* inrt_ovrf;      ///< Overflow bits for the intertial tensor
  const int* inrt_nxt_ovrf;  ///< Overflow bits for the intertial tensor in the next iteration
};

/// \brief Object to manage and stage recentering and momentum removal for a coordinate synthesis.
///        The object does not support asynchronous work unit scheduling but does contain two
///        alternating sets of accumulators so that one set can be initialized when the other is
///        in use.
class MotionSweeper {
public:

  /// \brief Class constructors take as their argument the PhaseSpaceSynthesis which the object
  ///        will serve.  A fixed-precision model can also be specified.
  /// \{
  MotionSweeper(const PhaseSpaceSynthesis *poly_ps,
                int momentum_bit_count_in = default_momentum_scale_bits,
                int center_of_mass_bit_count_in = default_com_scale_bits,
                int inertia_bit_count_in = default_inertia_scale_bits);

  MotionSweeper(const PhaseSpaceSynthesis &poly_ps,
                int momentum_bit_count_in = default_momentum_scale_bits,
                int center_of_mass_bit_count_in = default_com_scale_bits,
                int inertia_bit_count_in = default_inertia_scale_bits);
  /// \}

  /// \brief The copy and move constructors as wel as assignment operators require pointer repair
  ///        and must therefore be coded manually.
  ///
  /// \param original  The original object to copy or move
  /// \param other     The object on the right hand side of the assignment operator
  /// \{
  MotionSweeper(const MotionSweeper &original);
  MotionSweeper(MotionSweeper &&original);
  MotionSweeper& operator=(const MotionSweeper &other);  
  MotionSweeper& operator=(MotionSweeper &&other);  
  /// \}
  
  /// \brief Get the number of systems that the object is set up to serve.
  int getSystemCount() const;

  /// \brief Get the total number of work units spanning all systems.
  int getWorkUnitCount() const;

  /// \brief Get the object's current place in the time cycle.
  CoordinateCycle getCyclePosition() const;

  /// \brief Get the number of bits used for fixed-precision accumulation of the momentum.
  int getMomentumBitCount() const;

  /// \brief Get the number of bits used for fixed-precision accumulation in the center of mass
  ///        computation.
  int getCenterOfMassBitCount() const;

  /// \brief Get the number of bits used for fixed-precision accumulation of the rotational moment
  ///        of inertia tensor.
  int getInertialTensorBitCount() const;

  /// \brief Get the net velocity of a specific system, relative to the coordinate origin.
  ///
  /// \param idx   The index of the system of interest
  /// \param tier  The memory level from which to take the tensor components
  double3 getNetVelocity(int idx, HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the center of mass as found for a specific system.  Like the other accessors,
  ///        above, this function can be used to inspect the accumulations and determine the
  ///        results before all such quantities are purged.  Descriptions of input parameters
  ///        follow from getNetMomentum(), above.
  double3 getCenterOfMass(int idx, HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the inertial tensor accumulated for a specific system.  Descriptions of input
  ///        parameters follow from getNetVelocity(), above.
  std::vector<double> getInertialTensor(int idx,
                                        HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the angular velocity of the system about its center of mass.  Descriptions of
  ///        input parameters follow from getNetVelocity(), above.
  double3 getAngularVelocity(int idx, HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Modify the object's stage in the coordinate cycle.
  ///
  /// Overloaded:
  ///   - Increment the cycle position, as would be done from one momentum purge to the next over
  ///     the course of molecular dynamics calculations
  ///   - Set the time cycle to a specific point
  ///
  /// \param time_point  The point in the time cycle which the object is to take as "current"
  /// \{
  void updateCyclePosition();
  void updateCyclePosition(CoordinateCycle time_point);
  /// \}

  /// \brief Get the abstract.
  ///
  /// Overloaded:
  ///   - Produce a writeable abstract for a non-const object.
  ///   - Produce a read-only abstract for a const object.
  ///   - Take the abstract at a particular point in the time cycle, or accept the object's
  ///     current state
  ///
  /// \param orientation  Specify a point in the time cycle at which to take the abstract
  /// \param tier         Indicate whether to set pointers to memory on the CPU host or GPU device
  /// \{
  MotionSweepWriter data(CoordinateCycle orientation,
                         HybridTargetLevel tier = HybridTargetLevel::HOST);

  MotionSweepWriter data(HybridTargetLevel tier = HybridTargetLevel::HOST);

  const MotionSweepReader data(CoordinateCycle orientation,
                               HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  const MotionSweepReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  /// \}
  
#ifdef STORMM_USE_HPC
  /// \brief Upload the work unit instructions to the GPU
  void uploadWorkUnits();

  /// \brief Download the work unit instructions from the GPU
  void downloadWorkUnits();
  
  /// \brief Upload the object's data (including work unit instructions) to the GPU
  void uploadAll();
  
  /// \brief Download the object's data (including work unit instructions) from the GPU
  void downloadAll();
#endif
  
private:
  
  CoordinateCycle cycle_position;   ///< Indicates whether the object is performing the WHITE or
                                    ///<   BLACK move, and thus which accumulators to work with or
                                    ///<   initialize
  int center_of_mass_bit_count;     ///  Number of fixed-precision bits with which to accumulate
                                    ///<   the center of mass of each system
  int momentum_bit_count;           ///< Number of fixed-precision bits with which to accumulate
                                    ///<   the momentum of each system
  int inertia_bit_count;            ///< Number of fixed-precision bits with which to accumulate
                                    ///<   the inertial tensor components
  Hybrid<llint> xcom_white;         ///< Center of mass accumulators for the Cartesian X
                                    ///<   dimension, when the object is set on WHITE's move
  Hybrid<llint> ycom_white;         ///< Center of mass accumulators along Y, WHITE
  Hybrid<llint> zcom_white;         ///< Center of mass accumulators along Z, WHITE
  Hybrid<llint> xcom_black;         ///< Center of mass accumulators for the Cartesian X
                                    ///<   dimension, when the object is set on BLACK's move
  Hybrid<llint> ycom_black;         ///< Center of mass accumulators along Y, BLACK
  Hybrid<llint> zcom_black;         ///< Center of mass accumulators along Z, BLACK
  Hybrid<int> xcom_overflow_white;  ///< Center of mass overflow bits, Cartesian X dimension, for
                                    ///<   WHITE's move
  Hybrid<int> ycom_overflow_white;  ///< Center of mass overflow bits along Y for WHITE's move
  Hybrid<int> zcom_overflow_white;  ///< Center of mass overflow bits along Z for WHITE's move
  Hybrid<int> xcom_overflow_black;  ///< Center of mass overflow bits along X for BLACK's move
  Hybrid<int> ycom_overflow_black;  ///< Center of mass overflow bits along Y for BLACK's move
  Hybrid<int> zcom_overflow_black;  ///< Center of mass overflow bits along Z for BLACK's move
  Hybrid<double> total_mass;        ///< Total masses of all systems, in a.m.u. (g / mol) 
  Hybrid<llint> xmv_white;          ///< Net momentum accumulators for the Cartesian X dimension,
                                    ///<   when the object is set on WHITE's move
  Hybrid<llint> ymv_white;          ///< Net momentum accumulators for Cartesian Y, WHITE's move
  Hybrid<llint> zmv_white;          ///< Net momentum accumulators for Cartesian Z, WHITE's move
  Hybrid<llint> xmv_black;          ///< Net momentum accumulators for the Cartesian X dimension,
                                    ///<   when the object is set on BLACK's move
  Hybrid<llint> ymv_black;          ///< Net momentum accumulators for Cartesian Y, BLACK's move
  Hybrid<llint> zmv_black;          ///< Net momentum accumulators for Cartesian Z, BLACK's move
  Hybrid<int> xmv_overflow_white;   ///< Net momentum overflow accumulator bits for the Cartesian
                                    ///<   X dimension, when the object is performing WHITE's move
  Hybrid<int> ymv_overflow_white;   ///< Net momentum overflow bits along Y, WHITE's move
  Hybrid<int> zmv_overflow_white;   ///< Net momentum overflow bits along Z, WHITE's move
  Hybrid<int> xmv_overflow_black;   ///< Net momentum overflow bits along X, BLACK's move
  Hybrid<int> ymv_overflow_black;   ///< Net momentum overflow bits along Y, BLACK's move
  Hybrid<int> zmv_overflow_black;   ///< Net momentum overflow bits along Z, BLACK's move
  Hybrid<llint> rxmv_white;         ///< Angular momentum accumulators about the Cartesian X
                                    ///<   axis, when the object is set on WHITE's move
  Hybrid<llint> rymv_white;         ///< Angular momentum accumulators for the Y axis, WHITE's move
  Hybrid<llint> rzmv_white;         ///< Angular momentum accumulators for the Z axis, WHITE's move
  Hybrid<llint> rxmv_black;         ///< Angular momentum accumulators about the Cartesian X
                                    ///<   axis, when the object is set on BLACK's move
  Hybrid<llint> rymv_black;         ///< Angular momentum accumulators for the Y axis, BLACK's move
  Hybrid<llint> rzmv_black;         ///< Angular momentum accumulators for the Z axis, BLACK's move
  Hybrid<int> rxmv_overflow_white;  ///< Angular momentum overflow accumulator bits for the X axis,
                                    ///<   when the object is performing WHITE's move
  Hybrid<int> rymv_overflow_white;  ///< Angular momentum overflow bits about Y, WHITE's move
  Hybrid<int> rzmv_overflow_white;  ///< Angular momentum overflow bits about Z, WHITE's move
  Hybrid<int> rxmv_overflow_black;  ///< Angular momentum overflow bits about X, BLACK's move
  Hybrid<int> rymv_overflow_black;  ///< Angular momentum overflow bits about Y, BLACK's move
  Hybrid<int> rzmv_overflow_black;  ///< Angular momentum overflow bits about Z, BLACK's move

  /// Inertial moment accumulators used when the object is performing either WHITE's move or
  /// BLACK's.  The accumulators for each system are stored back-to-back, every six elements being
  /// devoted to one system in the order {x, x}, {x, y}, {x, z}, {y, y}, {y, z}, and {z, z}.
  /// \{
  Hybrid<llint> inertial_tensor_white;
  Hybrid<llint> inertial_tensor_black;
  /// \}

  /// Overflow accumulators for the intertial moment accumulators, used when the object is
  /// performing either WHITE's or BLACK's move.
  /// \{
  Hybrid<int> inertial_tensor_overflow_white;
  Hybrid<int> inertial_tensor_overflow_black;
  /// \}

  /// ARRAY-kind Hybrid object targeted by the POINTER-kind long long int Hybrid arrays above
  Hybrid<llint> llint_data;

  /// ARRAY-kind Hybrid object targeted by the POINTER-kind int Hybrid arrays above
  Hybrid<int> int_data;
  
  /// The total number of work unit instructions spanning all systems in the synthesis.
  int work_unit_count;
  
  /// Work units for each thread block.  All blocks will take 256 threads, with four blocks per
  /// streaming multiprocessor.  Each thread will read one atom's X, Y, and Z coordinates
  /// (positions and velocities), then compute the center of mass, total momentum, and (if needed)
  /// components of the intertial moment tensor.  The members of each tuple contain the index of
  /// the starting atom for the work unit to accumulate, upper limit of atom indices, and the
  /// system index to which the accumulation pertains in the "x", "y", and "z" members.
  Hybrid<int4> work_units;

  /// Pointer to the coordinate synthesis that each object of the class serves
  PhaseSpaceSynthesis *poly_ps_ptr;

  // Allocate memory in ARRAY-kind Hybrid objects targeted by POINTER-kind Hybrid accumulators.
  // Memory for work units and total system masses is allocated separately in the constructor.
  void allocate();
};
  
} // namespace trajectory
} // namespace stormm
#endif
