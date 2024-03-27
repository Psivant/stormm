// -*-c++-*-
#ifndef STORMM_RENDER_MOLECULE_H
#define STORMM_RENDER_MOLECULE_H

#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "copyright.h"
#include "Chemistry/periodic_table.h"
#include "Math/rounding.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "render_options.h"

namespace stormm {
namespace review {

using chemistry::element_maximum_count;
using chemistry::elemental_symbols;
using parse::char2ToString;
using parse::char4ToString;
using parse::intToString;
using parse::NumberFormat;
using parse::realToString;
using parse::removeTailingWhiteSpace;
using stmath::accumulateBitmask;
using stmath::readBitFromMask;
using stmath::roundUp;
using topology::AtomGraph;
using topology::ChemicalDetailsKit;
using topology::NonbondedKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameReader;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesReader;

/// \brief Trace a line through all bonds of a contiguous molecule, backtracking and retracing
///        steps as necessary.
///
/// \param nbk            Non-bonded parameters, including bond lists, for all atoms
/// \param cdk            Chemical details, including molecule limits, for the entire system
/// \param mol_index      Index of the molecule of interest
/// \param atom_coverage  Array of atoms covered thus far in the analysis, modified and returned
/// \param bond_coverage  Array of bonds covered thus far in the analysis, modified and returned
std::vector<int> traceBondLines(const NonbondedKit<double> &nbk, const ChemicalDetailsKit &cdk,
                                const int mol_index, std::vector<uint> *atom_coverage,
                                std::vector<uint> *bond_coverage);
  
/// \brief Plot the structure of a molecule by running a line from atom to atom, backtracking as
///        necessary to cover every atom and bond.
///
/// Overloaded:
///   - Provide the raw coordinate arrays and a frame count
///   - Provide a CoordinateFrame object
///   - Provide a CoordinateSeries object
/// 
/// \param foutp   Pointer to the output file
/// \param atom_x  Cartesian X coordinates of all frames of the molecular system
/// \param atom_y  Cartesian Y coordinates of all frames of the molecular system
/// \param atom_z  Cartesian Z coordinates of all frames of the molecular system
/// \param nbk     Non-bonded details of the molecule to plot, including the bonding structure
/// \param cdk     Chemical details of the molecule to plot, including atomic numbers of each atom
/// \param nframe  The number of coordinate frames.  All frames will be assumed to be located in
///                the coordinate arrays at intervals of the warp-padded number of atoms found in
///                cdk or nbk.
/// \{
template <typename T>
void renderMolecule(std::ofstream *foutp, const T* atom_x, const T* atom_y, const T* atom_z,
                    const NonbondedKit<double> &nbk, const ChemicalDetailsKit &cdk,
                    const RenderOptions &ropt, GridFileSyntax syntax, size_t nframe = 1);

void renderMolecule(std::ofstream *foutp, const CoordinateFrame *cf, const AtomGraph *ag,
                    const RenderOptions *ropt, GridFileSyntax syntax);

void renderMolecule(std::ofstream *foutp, const CoordinateFrame &cf, const AtomGraph &ag,
                    const RenderOptions &ropt, GridFileSyntax syntax);

template <typename T>
void renderMolecule(std::ofstream *foutp, const CoordinateSeries<T> *cs, const AtomGraph *ag,
                    const RenderOptions *ropt, GridFileSyntax syntax);

template <typename T>
void renderMolecule(std::ofstream *foutp, const CoordinateSeries<T> &cs, const AtomGraph &ag,
                    const RenderOptions &ropt, GridFileSyntax syntax);
/// \}

} // namespace review
} // namespace stormm

#include "render_molecule.tpp"

#endif
