// -*-c++-*-
#ifndef STORMM_WRITE_SDF_H
#define STORMM_WRITE_SDF_H

#include <vector>
#include "copyright.h"
#include "FileManagement/file_util.h"
#include "Synthesis/phasespace_synthesis.h"
#include "coordinateframe.h"
#include "coordinate_series.h"
#include "coordinate_swap_plan.h"
#include "phasespace.h"

namespace stormm {
namespace trajectory {

using diskutil::PrintSituation;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisReader;
  
/// \brief Write one or more frames to a trajectory file in a format-agnostic manner.  This
///        function accepts text strings, with all of the information for the various lines of the
///        file concatenated together with the requisite carriage returns.
///
/// Overloaded:
///   - Accept a single TextFile object, assumed to be one frame
///   - Accept a vector of TextFile objects, each assumed to be one frame
///   - Accept one or more TextFile objects with a corresponding number of CoordinateFrame or
///     PhaseSpace objects, an appropriately sized PhaseSpaceSynthesis object, a CoordinateSeries
///     object, or abstracts thereof, and vectors of line and column indices at which the
///     coordinates exist.  In this case, the coordinates of the TextFile objects will be replaced
///     with those in the coordinate objects before writing and the TextFile objects themselves
///     will be unaltered.
///
/// \param foutp
/// \param filename
/// \param tf
/// \param tf_list
/// \param xcrd
/// \param ycrd
/// \param zcrd
/// \param excision
/// \param ps             Object holding coordinates to swap into the xcrd, ycrd, and zcrd arrays
///                       according to the supplied plan.
/// \param time_point     The point in the coordinate cycle at which to set the abstract, and thus
///                       the coordinates which are to be read.  If not provided, the object's
///                       internal cycle counter will be referenced.
/// \param psr            Abstract of the PhaseSpace coordinate object holding coordinates to
///                       swap into the supplied plan.  The WHITE coordinates found in this
///                       abstract are referenced.
/// \param cf             Object holding coordinates to swap into the xcrd, ycrd, and zcrd arrays
///                       according to the supplied plan.
/// \param cfr            Abstract of the CoordinateFrame object.  No time point applies.
/// \param cs             A templated series of coordinates, all assumed to relate to the same
///                       topology and thus applicable to the same annotation.
/// \param poly_ps        
/// \param excision_list  A list of
/// \param plan_indices
/// \{
void writeFrame(std::ofstream *foutp, const TextFile &tf);

void writeFrame(const std::string &filename, PrintSituation expectation, const TextFile &tf);

void writeFrame(std::ofstream *foutp, const TextFile &tf, const double* xcrd, const double* ycrd,
                const double* zcrd, const CoordinateSwapPlan &excision);

void writeFrame(const std::string &filename, PrintSituation expectation, const TextFile &tf,
                const double* xcrd, const double* ycrd, const double* zcrd,
                const CoordinateSwapPlan &excision);

void writeFrame(const std::string &filename, PrintSituation expectation, const TextFile &tf,
                const PhaseSpaceReader &psr, const CoordinateSwapPlan &excision);

void writeFrame(const std::string &filename, PrintSituation expectation, const TextFile &tf,
                const PhaseSpace &ps, const CoordinateSwapPlan &excision);

void writeFrame(const std::string &filename, PrintSituation expectation, const TextFile &tf,
                const PhaseSpace &ps, const CoordinateSwapPlan &excision,
                CoordinateCycle time_point = CoordinateCycle::WHITE);

void writeFrame(const std::string &filename, PrintSituation expectation, const TextFile &tf,
                const CoordinateFrameReader &cfr, const CoordinateSwapPlan &excision);

void writeFrame(const std::string &filename, PrintSituation expectation, const TextFile &tf,
                const CoordinateFrame &cf, const CoordinateSwapPlan &excision);

template <typename T>
void writeFrame(std::ofstream *foutp, PrintSituation expectation, const TextFile &tf,
                const CoordinateSeriesReader<T> &csr, const CoordinateSwapPlan &excision);

template <typename T>
void writeFrame(std::ofstream *foutp, const TextFile &tf, const CoordinateSeries<T> &cs,
                const CoordinateSwapPlan &excision);
  
template <typename T>
void writeFrame(const std::string &filename, PrintSituation expectation, const TextFile &tf,
                const CoordinateSeries<T> &cs, const CoordinateSwapPlan &excision);

void writeFrame(std::ofstream *foutp, const std::vector<TextFile> &tf_list,
                const PsSynthesisReader &poly_psr,
                const std::vector<CoordinateSwapPlan> &excision_list,
                const std::vector<int> &plan_indices);

void writeFrame(const std::string &filename, PrintSituation expectation,
                const std::vector<TextFile> &tf_list, const PhaseSpaceSynthesis &poly_ps,
                const std::vector<CoordinateSwapPlan> &excision_list,
                const std::vector<int> &plan_indices);

void writeFrame(const std::string &filename, PrintSituation expectation,
                const std::vector<TextFile> &tf_list, const PhaseSpaceSynthesis &poly_ps,
                const std::vector<CoordinateSwapPlan> &excision_list,
                const std::vector<int> &plan_indices,
                CoordinateCycle time_point = CoordinateCycle::WHITE);
/// \}

} // namespace trajectory
} // namespace stormm

#include "write_annotated_frame.tpp"

#endif
