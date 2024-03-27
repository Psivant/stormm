// -*-c++-*-
#ifndef STORMM_CLASH_DETECTION_H
#define STORMM_CLASH_DETECTION_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/rounding.h"
#include "Math/series_ops.h"
#include "Namelists/nml_minimize.h"
#include "Potential/static_exclusionmask.h"
#include "Reporting/error_format.h"
#include "Synthesis/condensate.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"
#include "structure_enumerators.h"

namespace stormm {
namespace structure {

using constants::PrecisionModel;
using data_types::getStormmScalarTypeName;
using data_types::isFloatingPointScalarType;
using energy::StaticExclusionMask;
using energy::StaticExclusionMaskReader;
using energy::supertile_length;
using energy::tile_length;
using energy::tile_lengths_per_supertile;
using stmath::roundUp;
using stmath::indexingArray;
using namelist::default_minimize_clash_ratio;
using namelist::default_minimize_clash_r0;
using numerics::globalpos_scale_nonoverflow_bits;
using synthesis::Condensate;
using synthesis::CondensateReader;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisReader;
using topology::AtomGraph;
using topology::NonbondedKit;
using topology::ValenceKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesReader;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceReader;

/// \brief The minimum system size needed in order to engage a neighbor list based search.
constexpr int clash_direct_calculation_size_limit = 32;

/// \brief An object for listing the clashes between atoms of a system.
class ClashReport {
public:

  /// \brief The constructor accepts optional parameters for the clash parameters, or just a
  ///        pointer to the desired topology.
  /// \{
  ClashReport(double clash_distance_in = default_minimize_clash_r0,
              double clash_ratio_in = default_minimize_clash_ratio,
              const AtomGraph *ag_pointer_in = nullptr);

  ClashReport(const AtomGraph *ag_pointer_in);
  /// \}

  /// \brief Copy and move constructors, as well as copy and move assignment operators, can all
  ///        take their default forms for this object with no POINTER-kind Hybrids or other
  ///        pointers needing repair.
  ///
  /// \param original  Some other object to base the construction upon
  /// \param other     Some other object to use in assignment
  /// \{
  ClashReport(const ClashReport &original) = default;
  ClashReport(ClashReport &&original) = default;
  ClashReport& operator=(const ClashReport &original) = default;
  ClashReport& operator=(ClashReport &&original) = default;
  /// \}
  
  /// \brief Get the number of detected clashes.
  int getClashCount() const;
  
  /// \brief Get the minimum clash distance.
  double getMinimumDistance() const;

  /// \brief Get the minimum ratio of inter-particle distance to the pairwise sigma parameter.
  double getMinimumSigmaRatio() const;

  /// \brief Get a const pointer to the topology describing the system.
  const AtomGraph* getTopologyPointer() const;
  
  /// \brief Get one of the clashing pairs of atoms.
  ///
  /// \param index  The index of the clashing pair
  int2 getClashingPair(int index) const;

  /// \brief Get the distance at which two particles clash
  ///
  /// \param index  The index of the clashing pair
  double getClashDistance(int index) const;
  
  /// \brief Get the type of clash made by two particles
  ///
  /// \param index  The index of the clashing pair
  ClashKind getClashKind(int index) const;

  /// \brief Produce a string describing the atoms, indices, and relevant names of one of the
  ///        clashes catalogged in this object.
  ///
  /// Overloaded:
  ///   - Get a string with just the atoms and the inter-particle distance
  ///   - Get a string with the atoms, inter-particle distance, and current coordinates of each
  ///     particle
  ///
  /// \param clash_index   The index of the clashing pair
  /// \param crd_object    Object containing the exact coordinates within which the clash occurs
  /// \param frame         Index of the frame within a series
  /// \param system_index  Index of the system within a synthesis
  /// \{
  std::string getClashDescription(int clash_index) const;

  std::string getClashDescription(int clash_index, const CoordinateFrameReader &crd_object) const;

  std::string getClashDescription(int clash_index, const CoordinateFrame &crd_object) const;

  std::string getClashDescription(int clash_index, const CoordinateFrame *crd_object) const;

  std::string getClashDescription(int clash_index, const PhaseSpaceReader &crd_object) const;

  std::string getClashDescription(int clash_index, const PhaseSpace &crd_object) const;

  std::string getClashDescription(int clash_index, const PhaseSpace *crd_object) const;

  template <typename T>
  std::string getClashDescription(int clash_index, const CoordinateSeriesReader<T> &crd_object,
                                  int frame) const;

  template <typename T>
  std::string getClashDescription(int clash_index, const CoordinateSeries<T> &crd_object,
                                  int frame) const;

  template <typename T>
  std::string getClashDescription(int clash_index, const CoordinateSeries<T> *crd_object,
                                  int frame) const;

  std::string getClashDescription(int clash_index, const PsSynthesisReader &crd_object,
                                  int system_index) const;

  std::string getClashDescription(int clash_index, const PhaseSpaceSynthesis &crd_object,
                                  int system_index) const;

  std::string getClashDescription(int clash_index, const PhaseSpaceSynthesis *crd_object,
                                  int system_index) const;

  std::string getClashDescription(int clash_index, const CondensateReader &crd_object,
                                  int system_index) const;

  std::string getClashDescription(int clash_index, const Condensate &crd_object,
                                  int system_index) const;

  std::string getClashDescription(int clash_index, const Condensate *crd_object,
                                  int system_index) const;
  /// \}

  /// \brief Add notes about a clash between two particles.
  ///
  /// \param atom_i       The first atom of the clashing pair
  /// \param atom_j       The second atom of the clashing pair
  /// \param distance_in  Distance between the two particles in the clashing pair
  void addClash(int atom_i, int atom_j, double distance_in);

  /// \brief Clear the clashes to recompile the report.
  void clear();

  /// \brief Set the topology that the clash report shall use.
  ///
  /// Overloaded:
  ///   - Supply a pointer to the topology of interest
  ///   - Pass the topology of interest by reference
  ///
  /// \param ag_pointer_in  Pointer to the desired topology
  /// \param ag_in          Reference to the desired topology
  /// \{
  void setTopologyPointer(const AtomGraph *ag_pointer_in);
  void setTopologyPointer(const AtomGraph &ag_in);
  /// \}

  /// \brief Set the minimum absolute distance for two particles' separation before a clash is
  ///        declared.
  ///
  /// \param clash_distance_in
  void setMinimumDistance(double clash_distance_in);
  
  /// \brief Set the minimum ratio of interparticle distance to the pairwise sigma value, below
  ///        which a clash is declared.
  ///
  /// \param clash_ratio_in
  void setMinimumSigmaRatio(double clash_ratio_in);
  
private:
  int clash_count;
  double clash_distance;
  double clash_ratio;
  std::vector<int2> pairs;
  std::vector<double> distances;
  std::vector<ClashKind> kinds;
  AtomGraph* ag_pointer;

  /// \brief Validate the index of some catalogged element, to ensure that it is within the scope
  ///        of this report.
  ///
  /// \param index  The index of interest
  void validateClashIndex(int index) const;

  /// \brief Validate the system index provided against the coordinate object available.
  void validateSystemIndex(int system_index, int system_count) const;
  
  /// \brief Update the description of a clash with the locations of two atoms.
  ///
  /// \param x_i  Cartesian X coordinate of the first atom
  /// \param y_i  Cartesian Y coordinate of the first atom
  /// \param z_i  Cartesian Z coordinate of the first atom
  /// \param x_j  Cartesian X coordinate of the second atom
  /// \param y_j  Cartesian Y coordinate of the second atom
  /// \param z_j  Cartesian Z coordinate of the second atom
  std::string atomPairCoordinates(double x_i, double y_i, double z_i, double x_j, double y_j,
                                  double z_j) const;
};

/// \brief Compute the maximum distance at which two particles can be considered to participate in
///        a van-der Waals clash.
///
/// \param nbk         Non-bonded abstract for the system topology
/// \param elec_limit  The limiting distance at which two particles (with presumed electrostatic
///                    properties) will be considered to clash
/// \param vdw_ratio   The minimum required ratio of the distance between any pair of particles and
///                    mutual (non-bonded) sigma parameters
template <typename Tcalc>
Tcalc maxClashingDistance(const NonbondedKit<Tcalc> &nbk, const Tcalc elec_limit,
                          const Tcalc vdw_ratio);

/// \brief Implement a trivial test to see whether the (non-imaged) Cartesian ranges of two
///        cached sets of atom coordinates might overlap.  Return TRUE if there is a possibility,
///        FALSE if not.
///
/// \param cachi_xcrd  Cartesian X coordinates of the first set of particles
/// \param cachj_xcrd  Cartesian X coordinates of the second set of particles
/// \param cachi_ycrd  Cartesian Y coordinates of the first set of particles
/// \param cachj_ycrd  Cartesian Y coordinates of the second set of particles
/// \param cachi_zcrd  Cartesian Z coordinates of the first set of particles
/// \param cachj_zcrd  Cartesian Z coordinates of the second set of particles
/// \param ni_atoms    Number of atoms in the first set
/// \param nj_atoms    Number of atoms in the second set
/// \param max_clash   Maximum distance at which two particles might be deemed to clash, by either
///                    a raw distance or a van-der Waals sigma ratio criterion
template <typename Tcoord, typename Tcalc>
bool trivialClashCheck(const std::vector<Tcalc> &cachi_xcrd, const std::vector<Tcalc> &cachj_xcrd,
                       const std::vector<Tcalc> &cachi_ycrd, const std::vector<Tcalc> &cachj_ycrd,
                       const std::vector<Tcalc> &cachi_zcrd, const std::vector<Tcalc> &cachj_zcrd,
                       int ni_atoms, int nj_atoms, Tcalc max_clash);
  
/// \brief Perform a direct comparison of particle positions in a molecule to determine whether
///        there are any clashes.
///
/// \param xcrd        Cartesian X coordinates of all particles
/// \param ycrd        Cartesian Y coordinates of all particles
/// \param zcrd        Cartesian Z coordinates of all particles
/// \param nbk         Non-bonded abstract for the system topology
/// \param maskr       Read-only abstract for the exclusion mask
/// \param vdw_ratio   The minimum required ratio of the distance between any pair of particles and
///                    mutual (non-bonded) sigma parameters
/// \param elec_limit  The limiting distance at which two particles (with presumed electrostatic
///                    properties) will be considered to clash
/// \param inv_scale   Inverse scaling factor for converting coordinates into Angstrom units
template <typename Tcoord, typename Tcalc>
bool directClashTesting(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                        const NonbondedKit<Tcalc> &nbk, const StaticExclusionMaskReader &maskr,
                        Tcalc elec_limit = default_minimize_clash_r0,
                        Tcalc vdw_ratio = default_minimize_clash_ratio, Tcalc inv_scale = 1.0,
                        ClashReport *summary = nullptr);

/// \brief Compute the number of grid cells (and, incidentally, the origin and box lengths) of an
///        orthorhombic grid to be used in neighbbor list decomposition for clash detection within
///        non-periodic systems.
///
/// \param xcrd          Cartesian X coordinates of all particles
/// \param ycrd          Cartesian Y coordinates of all particles
/// \param zcrd          Cartesian Z coordinates of all particles
/// \param natom         The number of atoms in the system, trusted length of xcrd, ycrd, and zcrd
/// \param grid_origin   The Cartesian coordinates of the grid origin
/// \param grid_lengths  Lengths of each orthonormal grid edge
template <typename Tcoord>
int3 clashGridDecomposition(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                            const int natom, double3 *grid_origin, double3 *grid_lengths);
  
/// \brief Detect a van-der Waals clash between particles based on a minimum required ratio against
///        any given pair's non-bonded sigma ratio.  Return TRUE if a clash is found, FALSE if not.
///
/// Overloaded:
///   - Provide various coordinate objects, by const pointer or by const reference, or their
///     read-only abstracts
///   - Provide the original topology and exclusion mask, or their abstracts
///
/// Parameter descriptions in this routine follow from directClashTesting() above, with the
/// additions of:
///
/// \param cf          The system coordinates
/// \param ps          The system coordinates (the current point in the time cycle will be queried)
/// \param cs          The system coordinates (a frame number is required)
/// \param frame       Frame index to examine, if the coordinate object provided is a
///                    CoordinateSeries or a Condensate
/// \param cfr         Read-only abstract of the coordinates
/// \param psr         Read-only abstract of the coordinates
/// \param ag          System topology, containing van-der Waals parameters and each atom's type
///                    index
/// \param mask        Static exclusion mask for the system, indicating exclusions in an all-to-all
///                    interaction context
/// \{
template <typename Tcoord, typename Tcalc>
bool detectClash(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                 const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                 const StaticExclusionMask *mask, Tcalc elec_limit = default_minimize_clash_r0,
                 Tcalc vdw_ratio = default_minimize_clash_ratio, Tcalc inv_scale = 1.0,
                 ClashReport *summary = nullptr);

bool detectClash(const CoordinateFrameReader &cfr, const ValenceKit<double> &vk,
                 const NonbondedKit<double> &nbk, const StaticExclusionMask *mask,
                 double elec_limit = default_minimize_clash_r0,
                 double vdw_ratio = default_minimize_clash_ratio, ClashReport *summary = nullptr);

bool detectClash(const CoordinateFrame &cf, const AtomGraph &ag, const StaticExclusionMask &mask,
                 double elec_limit = default_minimize_clash_r0,
                 double vdw_ratio = default_minimize_clash_ratio, ClashReport *summary = nullptr);

bool detectClash(const CoordinateFrame *cf, const AtomGraph *ag, const StaticExclusionMask *mask,
                 double elec_limit = default_minimize_clash_r0,
                 double vdw_ratio = default_minimize_clash_ratio, ClashReport *summary = nullptr);

bool detectClash(const CoordinateFrame &cf, const AtomGraph &ag, const StaticExclusionMask &mask,
                 ClashReport *summary = nullptr);

bool detectClash(const CoordinateFrame *cf, const AtomGraph *ag, const StaticExclusionMask *mask,
                 ClashReport *summary = nullptr);

bool detectClash(const PhaseSpaceReader &psr, const ValenceKit<double> &vk,
                 const NonbondedKit<double> &nbk, const StaticExclusionMask *mask,
                 double elec_limit = default_minimize_clash_r0,
                 double vdw_ratio = default_minimize_clash_ratio, ClashReport *summary = nullptr);

bool detectClash(const PhaseSpace &ps, const AtomGraph &ag, const StaticExclusionMask &mask,
                 double elec_limit = default_minimize_clash_r0,
                 double vdw_ratio = default_minimize_clash_ratio, ClashReport *summary = nullptr);

bool detectClash(const PhaseSpace *ps, const AtomGraph *ag, const StaticExclusionMask *mask,
                 double elec_limit = default_minimize_clash_r0,
                 double vdw_ratio = default_minimize_clash_ratio, ClashReport *summary = nullptr);

bool detectClash(const PhaseSpace &ps, const AtomGraph &ag, const StaticExclusionMask &mask,
                 ClashReport *summary = nullptr);

bool detectClash(const PhaseSpace *ps, const AtomGraph *ag, const StaticExclusionMask *mask,
                 ClashReport *summary = nullptr);

template <typename Tcoord, typename Tcalc>
bool detectClash(const CoordinateSeriesReader<Tcoord> &csr, size_t frame,
                 const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                 const StaticExclusionMask *mask, Tcalc elec_limit = default_minimize_clash_r0,
                 Tcalc vdw_ratio = default_minimize_clash_ratio, ClashReport *summary = nullptr);

template <typename Tcoord, typename Tcalc>
bool detectClash(const CoordinateSeries<Tcoord> *cs, int frame, const AtomGraph *ag,
                 const StaticExclusionMask *mask, Tcalc elec_limit = default_minimize_clash_r0,
                 Tcalc vdw_ratio = default_minimize_clash_ratio, ClashReport *summary = nullptr);
  
template <typename Tcoord, typename Tcalc>
bool detectClash(const CoordinateSeries<Tcoord> &cs, int frame, const AtomGraph &ag,
                 const StaticExclusionMask &mask, Tcalc elec_limit = default_minimize_clash_r0,
                 Tcalc vdw_ratio = default_minimize_clash_ratio, ClashReport *summary = nullptr);

template <typename Tcoord, typename Tcalc>
bool detectClash(const CoordinateSeries<Tcoord> *cs, const int frame, const AtomGraph *ag,
                 const StaticExclusionMask *mask, ClashReport *summary = nullptr);

template <typename Tcoord, typename Tcalc>
bool detectClash(const CoordinateSeries<Tcoord> &cs, const int frame, const AtomGraph &ag,
                 const StaticExclusionMask &mask, ClashReport *summary = nullptr);

template <typename Tcalc>
bool detectClash(const PsSynthesisReader &poly_psr, const int system_index,
                 const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                 const StaticExclusionMask *mask, const Tcalc elec_limit, const Tcalc vdw_ratio,
                 ClashReport *summary = nullptr);

template <typename Tcalc>
bool detectClash(const PhaseSpaceSynthesis *poly_ps, const int system_index,
                 const StaticExclusionMask *mask, const Tcalc elec_limit, const Tcalc vdw_ratio,
                 ClashReport *summary = nullptr);

template <typename Tcalc>
bool detectClash(const PhaseSpaceSynthesis &poly_ps, const int system_index,
                 const StaticExclusionMask *mask, const Tcalc elec_limit, const Tcalc vdw_ratio,
                 ClashReport *summary = nullptr);

template <typename Tcalc>
bool detectClash(const PhaseSpaceSynthesis *poly_ps, const int system_index,
                 const StaticExclusionMask *mask, ClashReport *summary = nullptr);

template <typename Tcalc>
bool detectClash(const PhaseSpaceSynthesis &poly_ps, const int system_index,
                 const StaticExclusionMask *mask, ClashReport *summary = nullptr);

template <typename Tcalc>
bool detectClash(const CondensateReader &cdnsr, const int system_index,
                 const ValenceKit<Tcalc> &vk, const NonbondedKit<Tcalc> &nbk,
                 const StaticExclusionMask *mask, Tcalc elec_limit = default_minimize_clash_r0,
                 Tcalc vdw_ratio = default_minimize_clash_ratio, ClashReport *summary = nullptr);

template <typename Tcalc>
bool detectClash(const Condensate *cdns, const int system_index, const AtomGraph *ag,
                 const StaticExclusionMask *mask, double elec_limit = default_minimize_clash_r0,
                 double vdw_ratio = default_minimize_clash_ratio, ClashReport *summary = nullptr);

template <typename Tcalc>
bool detectClash(const Condensate &cdns, const int system_index, const AtomGraph &ag,
                 const StaticExclusionMask &mask, double elec_limit = default_minimize_clash_r0,
                 double vdw_ratio = default_minimize_clash_ratio, ClashReport *summary = nullptr);

template <typename Tcalc>
bool detectClash(const Condensate *cdns, const int system_index, const AtomGraph *ag,
                 const StaticExclusionMask *mask, ClashReport *summary = nullptr);

template <typename Tcalc>
bool detectClash(const Condensate &cdns, const int system_index, const AtomGraph &ag,
                 const StaticExclusionMask &mask, ClashReport *summary = nullptr);
/// \}

} // namespace structure
} // namespace stormm

#include "clash_detection.tpp"

#endif
