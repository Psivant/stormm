// -*-c++-*-
#ifndef STORMM_LOCAL_ARRANGEMENT_H
#define STORMM_LOCAL_ARRANGEMENT_H

#include <cmath>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "Constants/symbol_values.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/rounding.h"
#include "Math/vector_ops.h"
#include "Numerics/split_fixed_precision.h"
#include "Reporting/error_format.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"
#include "structure_enumerators.h"

namespace stormm {
namespace structure {

using card::Hybrid;
using constants::PrecisionModel;
using data_types::isSignedIntegralScalarType;
using numerics::globalpos_scale_nonoverflow_bits;
using numerics::hostInt95ToDouble;
using numerics::hostInt95Sum;
using stmath::angleVerification;
using stmath::crossProduct;
using stmath::roundUp;
using synthesis::Condensate;
using synthesis::CondensateReader;
using synthesis::CondensateWriter;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisReader;
using synthesis::PsSynthesisWriter;
using topology::UnitCellType;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateFrameWriter;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesReader;
using trajectory::CoordinateSeriesWriter;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceReader;
using trajectory::PhaseSpaceWriter;

/// \brief Image a single value to a range.
///
/// \param x                  The value to re-image
/// \param range              Defines the interval over which to image the value
/// \param style              Determines whether to image the value to the interval [ 0, range )
///                           or [ -0.5 * range, 0.5 * range ).
/// \param gpos_scale_factor  Conversion factor for fixed precision coordinates
template <typename Tcoord, typename Tcalc>
Tcoord imageValue(Tcoord x, Tcalc range, ImagingMethod style, Tcalc gpos_scale_factor = 1.0);

/// \brief Image a coordinate (or distance between coordinates) to fit within the interval
///        [-0.5, 0.5], taking into account box specifications.
///
/// Overloaded:
///   - Operate on a single (x, y, z) tuple
///   - Operate on three C-style arrays of real numbers x, y, and z
///   - Operate on three std::vectors of real numbers x, y, and z
///   - Operate on three Hybrid arrays of real numbers x, y, and z
///
/// \param x                  Cartesian x coordinate of one or more particles
/// \param y                  Cartesian y coordinate of one or more particles
/// \param z                  Cartesian z coordinate of one or more particles
/// \param length             Number of particles or coordinate tuples to re-image
/// \param umat               Transformation matrix to go into fractional coordinates
/// \param invu               Transformation matrix from fractional coordinates back to real space
/// \param unit_cell          Shape of the unit cell
/// \param gpos_scale_factor  Conversion factor for fixed precision coordinates
/// \param x_ovrf             Overflow bits for Cartesian x coordinates
/// \param y_ovrf             Overflow bits for Cartesian y coordinates
/// \param z_ovrf             Overflow bits for Cartesian z coordinates
/// \{
template <typename Tcoord, typename Tcalc>
void imageCoordinates(Tcoord *x, Tcoord *y, Tcoord *z, const double* umat, const double* invu,
                      UnitCellType unit_cell, ImagingMethod style, Tcalc gpos_scale_factor = 1.0);

template <typename Tcoord, typename Tcalc>
void imageCoordinates(Tcoord* x, Tcoord* y, Tcoord* z, const int length, const double* umat,
                      const double* invu, UnitCellType unit_cell, ImagingMethod style,
                      Tcalc gpos_scale_factor = 1.0, int* xcrd_ovrf = nullptr,
                      int* ycrd_ovrf = nullptr, int* zcrd_ovrf = nullptr);

template <typename Tcoord, typename Tcalc>
void imageCoordinates(std::vector<Tcoord> *x, std::vector<Tcoord> *y, std::vector<Tcoord> *z,
                      const double* umat, const double* invu, UnitCellType unit_cell,
                      ImagingMethod style, Tcalc gpos_scale_factor = 1.0);

template <typename Tcoord, typename Tcalc>
void imageCoordinates(Hybrid<Tcoord> *x, Hybrid<Tcoord> *y, Hybrid<Tcoord> *z,
                      const double* umat, const double* invu, UnitCellType unit_cell,
                      ImagingMethod style, Tcalc gpos_scale_factor = 1.0);

void imageCoordinates(CoordinateFrameWriter cfw, ImagingMethod style);

void imageCoordinates(CoordinateFrame *cf, ImagingMethod style);

void imageCoordinates(PhaseSpaceWriter psw, ImagingMethod style);

void imageCoordinates(PhaseSpace *ps, ImagingMethod style);

template <typename Tcoord, typename Tcalc>
void imageCoordinates(CoordinateSeriesWriter<Tcoord> csw, size_t frame_index, ImagingMethod style);

template <typename Tcoord, typename Tcalc>
void imageCoordinates(CoordinateSeries<Tcoord> *cs, size_t frame_index, ImagingMethod style);

template <typename Tcalc>
void imageCoordinates(CondensateWriter cdnsw, int system_index, ImagingMethod style);

template <typename Tcalc>
void imageCoordinates(Condensate *cdns, int system_index, ImagingMethod style);

template <typename Tcalc>
void imageCoordinates(PsSynthesisWriter poly_psw, int system_index, ImagingMethod style);

template <typename Tcalc>
void imageCoordinates(PhaseSpaceSynthesis *poly_ps, int system_index, ImagingMethod style);
/// \}

/// \brief Compute the distance between two points in a coordinate set.  All forms of this
///        function work with double-precision coordinates on the host.
///
/// Overloaded:
///   - Accept two 
///   - Accept two atom indices and raw pointers to templated coordinate arrays (float, double, or
///     fixed precision signed integers) and templated box specifications (float or double, which
///     will determine the numerical precision of the internal calculations)
///   - Accept two atom indices and any of the coordinate storage objects
///
/// \param pti_x              Split fixed-precision, high accuracy representation of the ith
///                           particle's Cartesian X coordinate
/// \param pti_y              The ith particle's Cartesian Y coordinate
/// \param pti_z              The ith particle's Cartesian Z coordinate
/// \param ptj_x              The jth particle's Cartesian X coordinate
/// \param ptj_y              The jth particle's Cartesian Y coordinate
/// \param ptj_z              The jth particle's Cartesian Z coordinate
/// \param atom_i             Topological index of the first atom in the system
/// \param atom_j             Topological index of the second atom in the system
/// \param xcrd               Cartesian x coordinate of all particles in the system
/// \param ycrd               Cartesian y coordinate of all particles in the system
/// \param zcrd               Cartesian z coordinate of all particles in the system
/// \param umat               Transformation matrix to go into fractional coordinates
/// \param invu               Transformation matrix from fractional coordinates back to real space
/// \param unit_cell          Shape of the unit cell
/// \param gpos_scale_factor  Conversion factor for fixed precision coordinates
/// \param x_ovrf             Overflow bits for Cartesian x coordinates
/// \param y_ovrf             Overflow bits for Cartesian y coordinates
/// \param z_ovrf             Overflow bits for Cartesian z coordinates
/// \{
template <typename Tcalc4, typename Tcalc>
Tcalc4 distance(const int95_t pti_x, const int95_t pti_y, const int95_t pti_z, const int95_t ptj_x,
               	const int95_t ptj_y, const int95_t ptj_z, const	double*	umat, const double* invu,
                UnitCellType unit_cell, Tcalc gpos_scale_factor);

template <typename Tcoord, typename Tcalc>
Tcalc distance(int atom_i, int atom_j, const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
               const double* umat, const double* invu, UnitCellType unit_cell,
               Tcalc gpos_scale_factor = 1.0, const int* xcrd_ovrf = nullptr,
               const int* ycrd_ovrf = nullptr, const int* zcrd_ovrf = nullptr);

double distance(int atom_i, int atom_j, const CoordinateFrameReader &cfr);

double distance(int atom_i, int atom_j, const CoordinateFrame *cf);

double distance(int atom_i, int atom_j, const CoordinateFrame &cf);

double distance(int atom_i, int atom_j, const PhaseSpaceReader &psr);

double distance(int atom_i, int atom_j, const PhaseSpace *ps);

double distance(int atom_i, int atom_j, const PhaseSpace &ps);

template <typename Tcoord, typename Tcalc>
Tcalc distance(int atom_i, int atom_j, const CoordinateSeriesReader<Tcoord> &csr,
               size_t frame_index);

template <typename Tcoord, typename Tcalc>
Tcalc distance(int atom_i, int atom_j, const CoordinateSeries<Tcoord> *cs, size_t frame_index);

template <typename Tcoord, typename Tcalc>
Tcalc distance(int atom_i, int atom_j, const CoordinateSeries<Tcoord> &cs, size_t frame_index);

template <typename Tcalc>
Tcalc distance(int atom_i, int atom_j, const CondensateReader &cdnsr, int system_index);

template <typename Tcalc>
Tcalc distance(int atom_i, int atom_j, const Condensate *cdns, int system_index);

template <typename Tcalc>
Tcalc distance(int atom_i, int atom_j, const Condensate &cdns, int system_index);

template <typename Tcalc>
Tcalc distance(int atom_i, int atom_j, const PsSynthesisReader &poly_psr, int system_index);

template <typename Tcalc>
Tcalc distance(int atom_i, int atom_j, const PhaseSpaceSynthesis *poly_ps, int system_index);

template <typename Tcalc>
Tcalc distance(int atom_i, int atom_j, const PhaseSpaceSynthesis &poly_ps, int system_index);
/// \}

/// \brief Compute the angle between three points in a coordinate set.  All forms of this
///        function work with double-precision coordinates on the host.  The result is returned in
///        radians.
///
/// Overloaded:
///   - Accept three atom indices and raw pointers to templated coordinate arrays (float, double,
///     or fixed precision signed integers) and templated box specifications (float or double,
///     which will determine the numerical precision of the internal calculations)
///   - Accept three atom indices and any of the coordinate storage objects
///
/// Parameters for this function follow from descriptions in distance, above
/// \{
template <typename Tcoord, typename Tcalc>
Tcalc angle(int atom_i, int atom_j, int atom_k, const Tcoord* xcrd, const Tcoord* ycrd,
            const Tcoord* zcrd, const double* umat, const double* invu, UnitCellType unit_cell,
            Tcalc gpos_scale_factor = 1.0, const int* xcrd_ovrf = nullptr,
            const int* ycrd_ovrf = nullptr, const int* zcrd_ovrf = nullptr);

double angle(int atom_i, int atom_j, int atom_k, const CoordinateFrameReader &cfr);

double angle(int atom_i, int atom_j, int atom_k, const CoordinateFrame *cf);

double angle(int atom_i, int atom_j, int atom_k, const CoordinateFrame &cf);

double angle(int atom_i, int atom_j, int atom_k, const PhaseSpaceReader &psr);

double angle(int atom_i, int atom_j, int atom_k, const PhaseSpace *ps);

double angle(int atom_i, int atom_j, int atom_k, const PhaseSpace &ps);

template <typename Tcoord, typename Tcalc>
Tcalc angle(int atom_i, int atom_j, int atom_k, const CoordinateSeriesReader<Tcoord> &csr,
            size_t frame_index);

template <typename Tcoord, typename Tcalc>
Tcalc angle(int atom_i, int atom_j, int atom_k, const CoordinateSeries<Tcoord> *cs,
            size_t frame_index);

template <typename Tcoord, typename Tcalc>
Tcalc angle(int atom_i, int atom_j, int atom_k, const CoordinateSeries<Tcoord> &cs,
            size_t frame_index);

template <typename Tcalc>
Tcalc angle(int atom_i, int atom_j, int atom_k, const CondensateReader &cdnsr, int system_index);

template <typename Tcalc>
Tcalc angle(int atom_i, int atom_j, int atom_k, const Condensate *cdns, int system_index);

template <typename Tcalc>
Tcalc angle(int atom_i, int atom_j, int atom_k, const Condensate &cdns, int system_index);

template <typename Tcalc>
Tcalc angle(int atom_i, int atom_j, int atom_k, const PsSynthesisReader &poly_psr,
            int system_index);

template <typename Tcalc>
Tcalc angle(int atom_i, int atom_j, int atom_k, const PhaseSpaceSynthesis *poly_ps,
            int system_index);

template <typename Tcalc>
Tcalc angle(int atom_i, int atom_j, int atom_k, const PhaseSpaceSynthesis &poly_ps,
            int system_index);
/// \}

/// \brief Compute the dihedral angle between three points in a coordinate set.  All forms of this
///        function work with double-precision coordinates on the host.  The result is returned in
///        radians.
///
/// Overloaded:
///   - Accept four atom indices and raw pointers to templated coordinate arrays (float, double,
///     or fixed precision signed integers) and templated box specifications (float or double,
///     which will determine the numerical precision of the internal calculations)
///   - Accept four atom indices and any of the coordinate storage objects
///
/// Parameters for this function follow from descriptions in distance, above
/// \{
template <typename Tcoord, typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l, const Tcoord* xcrd,
                    const Tcoord* ycrd, const Tcoord* zcrd, const double* umat,
                    const double* invu, UnitCellType unit_cell, Tcalc gpos_scale_factor = 1.0,
                    const int* xcrd_ovrf = nullptr, const int* ycrd_ovrf = nullptr,
                    const int* zcrd_ovrf = nullptr);

double dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l,
                     const CoordinateFrameReader &cfr);

double dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l, const CoordinateFrame *cf);

double dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l, const CoordinateFrame &cf);

double dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l, const PhaseSpaceReader &cfr);

double dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l, const PhaseSpace *ps);

double dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l, const PhaseSpace &ps);

template <typename Tcoord, typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l,
                    const CoordinateSeriesReader<Tcoord> &csr, size_t frame_index);

template <typename Tcoord, typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l,
                    const CoordinateSeries<Tcoord> *cs, size_t frame_index);

template <typename Tcoord, typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l,
                    const CoordinateSeries<Tcoord> &cs, size_t frame_index);

template <typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l, const CondensateReader &cdnsr,
                    int system_index);

template <typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l, const Condensate *cdns,
                    int system_index);

template <typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l, const Condensate &cdns,
                    int system_index);

template <typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l,
                    const PsSynthesisReader &poly_psr, int system_index);

template <typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l,
                    const PhaseSpaceSynthesis *poly_ps, int system_index);

template <typename Tcalc>
Tcalc dihedralAngle(int atom_i, int atom_j, int atom_k, int atom_l,
                    const PhaseSpaceSynthesis &poly_ps, int system_index);
/// \}

} // namespace structure
} // namespace stormm

#include "local_arrangement.tpp"

#endif
