// -*-c++-*-
#ifndef STORMM_RADIUS_GYRATION_H
#define STORMM_RADIUS_GYRATION_H

#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"

namespace stormm {
namespace structure {

using card::Hybrid;
using card::HybridTargetLevel;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisReader;
using synthesis::SyAtomUpdateKit;
using topology::AtomGraph;
using topology::ChemicalDetailsKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesReader;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceReader;

/// \brief Compute the radius of gyration for a system or collection of systems.
///
/// Overloaded:
///   - Compute the radius of gyration for one or more systems provided in various coordinate
///     objects
///   - Return the radius of gyration for multiple systems or multiple parts of one system as a
///     Standard Template Library vector
///   - Fill an existing array (Standard Template Library vector or Hybrid object) with the
///     requested radius of gyration calculation results
///
/// \param cf       The structure of interest
/// \param ps       The structure of interest
/// \param cs       The series of structures of interest
/// \param poly_ps  The structure of interest
/// \param cdns     The structure of interest
/// \param ag       Topology for the system of interest, containing particle masses
/// \param poly_ag  Topology synthesis for all systems, containing particle masses
/// \{
double radiusOfGyration(const CoordinateFrameReader &cfr, const ChemicalDetailsKit &cdk,
                        HybridTargetLevel tier = HybridTargetLevel::HOST);

double radiusOfGyration(const CoordinateFrame *cf, const AtomGraph &ag,
                        HybridTargetLevel tier = HybridTargetLevel::HOST);

double radiusOfGyration(const CoordinateFrame &cf, const AtomGraph &ag,
                        HybridTargetLevel tier = HybridTargetLevel::HOST);

double radiusOfGyration(const PhaseSpace *ps, const AtomGraph &ag,
                        HybridTargetLevel tier = HybridTargetLevel::HOST);

double radiusOfGyration(const PhaseSpace &ps, const AtomGraph &ag,
                        HybridTargetLevel tier = HybridTargetLevel::HOST);

std::vector<double> radiusOfGyration(const PsSynthesisReader &poly_psr,
                                     const SyAtomUpdateKit<double, double2, double4> &poly_auk,
                                     HybridTargetLevel tier = HybridTargetLevel::HOST);

std::vector<double> radiusOfGyration(const PhaseSpaceSynthesis *poly_ps,
                                     const AtomGraphSynthesis &poly_ag,
                                     HybridTargetLevel tier = HybridTargetLevel::HOST);

std::vector<double> radiusOfGyration(const PhaseSpaceSynthesis &poly_ps,
                                     const AtomGraphSynthesis &poly_ag,
                                     HybridTargetLevel tier = HybridTargetLevel::HOST);
/// \}

} // namespace structure
} // namespace stormm

#endif
