// -*-c++-*-
#ifndef STORMM_STRUCTURE_RMSD_H
#define STORMM_STRUCTURE_RMSD_H

#include "copyright.h"
#include "Analysis/comparison_guide.h"
#include "Constants/behavior.h"
#include "DataTypes/common_types.h"
#include "Math/matrix_ops.h"
#include "Math/rounding.h"
#include "Math/summation.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Synthesis/synthesis_cache_map.h"
#include "Synthesis/synthesis_enumerators.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"
#include "rmsd_plan.h"
#include "structure_enumerators.h"

namespace stormm {
namespace structure {

using analysis::ComparisonGuide;
using constants::PrecisionModel;
using data_types::isSignedIntegralScalarType;
using stmath::EigenWork;
using stmath::jacobiEigensolver;
using stmath::realSymmEigensolver;
using stmath::roundUp;
using stmath::sum;
using topology::AtomGraph;
using topology::ChemicalDetailsKit;
using synthesis::Condensate;
using synthesis::CondensateReader;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisReader;
using synthesis::SyNonbondedKit;
using synthesis::SynthesisCacheMap;
using synthesis::SynthesisMapReader;
using synthesis::SystemGrouping;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateFrameWriter;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesReader;
using trajectory::CoordinateSeriesWriter;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceReader;
using trajectory::PhaseSpaceWriter;

/// \brief Verify the compatibility of a plan for computing RMSD values and the synthesis from
///        which coordinates will be drawn.
///
/// Overloaded:
///   - Include only the plan and the synthesis itself
///   - Also include a condensed copy of the synthesis coordinates
///   - Also include a synthesis cache map relating systems back to original user specifications
///
/// \param cg            Guide for comparing structures within the synthesis
/// \param rplan         Plan for computing RMSD values for any topology included in the synthesis
/// \param snapshots     The synthesis of coordinates
/// \param cdns          Condensed form of the coordinate synthesis
/// \param result        Vector that shall hold the results (not modified herein, simply passed in
///                      to check its size)
/// \param process       The kind of RMSD calculations to perform
/// \param organization  The manner in which systems within the underlying synthesis are grouped
///                      for RMSD comparisons
/// \{
template <typename T>
void checkCompatibility(const ComparisonGuide &cg, const RMSDPlan &rplan,
                        const PhaseSpaceSynthesis &poly_ps, const Hybrid<T> *result,
                        RMSDTask process, SystemGrouping organization);

template <typename T>
void checkCompatibility(const ComparisonGuide &cg, const RMSDPlan &rplan,
                        const PhaseSpaceSynthesis &poly_ps, const Condensate &cdns,
                        const Hybrid<T> *result, RMSDTask process, SystemGrouping organization);
/// \}

/// \brief Compute the positional RMSD between two sets of coordinates.
///
/// Overloaded:
///   - Use const pointers to coordinates and mass
///   - Accept PhaseSpace objects and the topology as const pointers or references, or take
///     the relevant abstracts
///   - Accept CoordinateFrame objects and the topology as const pointers or references, or take
///     the relevant abstracts
///   - Accept CoordinateSeries objects and the topology as const pointers or references, or
///     take the relevant abstracts, and return a Standard Template Library vector of results
///     for all frames of one series to all frames of the other
///   - Use an RMSDPlan with various coordinate objects
///
/// \param xcrd_a            Cartesian X coordinates of the first frame
/// \param xcrd_b            Cartesian X coordinates of the second frame
/// \param ycrd_a            Cartesian Y coordinates of the first frame
/// \param ycrd_b            Cartesian Y coordinates of the second frame
/// \param zcrd_a            Cartesian Z coordinates of the first frame
/// \param zcrd_b            Cartesian Z coordinates of the second frame
/// \param masses            Masses of all particles
/// \param method            Alignment and mass-weighting scheme for computing RMSD
/// \param psr_a             Coordinate abstract for the first frame
/// \param psr_b             Coordinate abstract for the second frame
/// \param psw_a             Coordinate abstract for the first frame
/// \param psw_b             Coordinate abstract for the second frame
/// \param cdk               Chemical details from the topology describing both frames (for masses)
/// \param cfr_a             Coordinate abstract for the first frame
/// \param cfr_b             Coordinate abstract for the second frame
/// \param cfw_a             Coordinate abstract for the first frame
/// \param cfw_b             Coordinate abstract for the second frame
/// \param ps_a              Coordinates for the first frame
/// \param ps_b              Coordinates for the second frame
/// \param cf_a              Coordinates for the first frame
/// \param cf_b              Coordinates for the second frame
/// \param csw               Coordinate series containing both frames
/// \param cs                Coordinate series containing both frames
/// \param frame_a           Index of the first frame, if supplying a CoordinateSeries object
/// \param frame_b           Index of the second frame, if supplying a CoordinateSeries object
/// \param reference         The reference frame
/// \param reference_frame   Index of the reference frame within an array of coordinate sets
///                          pertaining to the same system
/// \param reference_frames  Indices of the reference frames within a synthesis of structures.
///                          Each index pertains to the absolute position of the structure in
///                          the overall synthesis, not its relative position among the list of
///                          structures that are instances of a particular topology.
/// \param snapshot          The frame of interest, to compare to the reference, or the
/// \param snapshots         Array of frames to compare to the reference(s)
/// \param ag                Topology for both frames (contains masses)
/// \param lower_limit       Lower limit of the atoms to cover in RMSD calculation
/// \param upper_limit       Upper limit of the atoms to cover in RMSD calculation
/// \{
template <typename Tcoord, typename Tcalc>
Tcalc rmsd(const Tcoord* xcrd_a, const Tcoord* ycrd_a, const Tcoord* zcrd_a, const Tcoord* xcrd_b,
           const Tcoord* ycrd_b, const Tcoord* zcrd_b, const Tcalc* masses, RMSDMethod method,
           int lower_limit, int upper_limit, Tcalc inv_gpos_scale_factor = 1.0);
  
double rmsd(const double* xcrd_a, const double* ycrd_a, const double* zcrd_a, const double* xcrd_b,
            const double* ycrd_b, const double* zcrd_b, const double* masses, RMSDMethod method,
            int lower_limit, int upper_limit);

double rmsd(const PhaseSpaceReader &psr_a, const PhaseSpaceReader &psr_b,
            const ChemicalDetailsKit &cdk, RMSDMethod method, int lower_limit = 0,
            int upper_limit = 0);

double rmsd(const PhaseSpaceWriter &psw_a, const PhaseSpaceWriter &psw_b,
            const ChemicalDetailsKit &cdk, RMSDMethod method, int lower_limit = 0,
            int upper_limit = 0);

double rmsd(const PhaseSpace &ps_a, const PhaseSpace &ps_b, const AtomGraph &ag, RMSDMethod method,
            int lower_limit = 0, int upper_limit = 0);

double rmsd(const CoordinateFrameReader &cfr_a, const CoordinateFrameReader &cfr_b,
            const ChemicalDetailsKit &cdk, RMSDMethod method, int lower_limit = 0,
            int upper_limit = 0);

double rmsd(const CoordinateFrameWriter &cfw_a, const CoordinateFrameWriter &cfw_b,
            const ChemicalDetailsKit &cdk, RMSDMethod method, int lower_limit = 0,
            int upper_limit = 0);

double rmsd(const CoordinateFrame &cf_a, const CoordinateFrame &cf_b, const AtomGraph &ag,
            RMSDMethod method, int lower_limit = 0, int upper_limit = 0);

template <typename Tcoord, typename Tcalc>
Tcalc rmsd(const CoordinateSeriesReader<Tcoord> &csr, const size_t frame_a, const size_t frame_b,
           const ChemicalDetailsKit &cdk, const RMSDMethod method, const int lower_limit,
           const int upper_limit);

template <typename Tcoord, typename Tcalc>
Tcalc rmsd(const CoordinateSeriesWriter<Tcoord> &csw, const size_t frame_a, const size_t frame_b,
           const ChemicalDetailsKit &cdk, const RMSDMethod method, const int lower_limit,
           const int upper_limit);

template <typename Tcoord, typename Tcalc>
Tcalc rmsd(const CoordinateSeries<Tcoord> &cs, const size_t frame_a, const size_t frame_b,
           const AtomGraph &ag, const RMSDMethod method, const int lower_limit,
           const int upper_limit);

double rmsd(const RMSDPlan &rplan, const PhaseSpaceReader &reference,
            const PhaseSpaceReader &snapshot);

double rmsd(const RMSDPlan &rplan, const PhaseSpace &reference, const PhaseSpace &snapshot);

double rmsd(const RMSDPlan &rplan,
            const CoordinateFrameReader &reference, const CoordinateFrameReader &snapshot);

double rmsd(const RMSDPlan &rplan, const CoordinateFrame &reference,
            const CoordinateFrame &snapshot);

template <typename T>
std::vector<double> rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan,
                         const CoordinateFrame &reference, const CoordinateSeries<T> &snapshots);

template <typename T>
std::vector<double> rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan,
                         const CoordinateSeries<T> &snapshots, int reference_frame);

template <typename T>
void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const CoordinateFrame &reference,
          const CoordinateSeries<T> &snapshots, Hybrid<T> *result);

template <typename T>
void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const CoordinateSeries<T> &snapshots,
          int reference_frame, Hybrid<T> *result);

template <typename T>
void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const CoordinateSeries<T> &snapshots,
          Hybrid<T> *result);

template <typename T>
void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const PhaseSpaceSynthesis &snapshots,
          const Hybrid<int> &reference_frames, Hybrid<T> *result,
          SystemGrouping organization = SystemGrouping::TOPOLOGY);

template <typename T>
void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const PhaseSpaceSynthesis &snapshots,
          Hybrid<T> *result, SystemGrouping organization = SystemGrouping::TOPOLOGY);

template <typename T>
void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const PhaseSpaceSynthesis &snapshots,
          const Condensate &cdns, const Hybrid<int> &reference_frames, Hybrid<T> *result,
          SystemGrouping organization = SystemGrouping::TOPOLOGY);

template <typename T>
void rmsd(const ComparisonGuide &cg, const RMSDPlan &rplan, const PhaseSpaceSynthesis &snapshots,
          const Condensate &cdns, Hybrid<T> *result,
          SystemGrouping organization = SystemGrouping::TOPOLOGY);
/// \}

} // namespace structure
} // namespace stormm

#include "rmsd.tpp"

#endif

