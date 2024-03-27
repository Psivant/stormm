// -*-c++-*-
#ifndef STORMM_ISOMERIZATION_H
#define STORMM_ISOMERIZATION_H

#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/fixed_precision.h"
#include "Chemistry/chemical_features.h"
#include "Chemistry/chemistry_enumerators.h"
#include "DataTypes/common_types.h"
#include "Math/rounding.h"
#include "Topology/atomgraph.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"

namespace stormm {
namespace structure {

using chemistry::ChemicalFeatures;
using chemistry::ChiralInversionProtocol;
using chemistry::IsomerPlan;
using constants::PrecisionModel;
using data_types::isSignedIntegralScalarType;
using numerics::globalpos_scale_nonoverflow_bits;
using stmath::roundUp;
using synthesis::Condensate;
using synthesis::CondensateWriter;
using synthesis::PhaseSpaceSynthesis;
using synthesis::PsSynthesisWriter;
using topology::AtomGraph;
using topology::UnitCellType;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameWriter;
using trajectory::CoordinateSeries;
using trajectory::CoordinateSeriesWriter;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceWriter;

/// \brief Rotate a molecule such that selected internal coordinates achieve particular values.
///        The internal coordinates must define a rotatable bond (bonds in rings are prohibited
///        from rotation under this procedure).  This operation takes place only on the CPU.  The
///        results must be uploaded to the GPU but can serve to seed coordinates for more detailed
///        manipulations.
///
/// Overloaded:
///   - Take the X, Y, and Z particle coordinates directly, in a templated function that will
///     handle the fixed-precision conversions internally if needed
///   - Take a CoordinateFrame, PhaseSpace, PhaseSpaceSynthesis, or CoordinateSeries object, or
///     any of their writeable abstracts, for the coordinates
///
/// \param xcrd                    Cartesian X coordinates of particles
/// \param ycrd                    Cartesian Y coordinates of particles
/// \param zcrd                    Cartesian Z coordinates of particles
/// \param cf                      Coordinates to manipulate (will be modified upon return)
/// \param cfw                     CoordianteFrame writer abstract
/// \param ps                      PhaseSpace object with coordinates, velocities, and forces (the
///                                coordinates will be modified upon return)
/// \param psw                     PhaseSpace writer abstract
/// \param psynth                  Collection of systems, with coordinates, velocities, and forces
///                                represented in fixed precision
/// \param psynthw                 PhaseSpaceSynthesis writeable abstract
/// \param cs                      Collection of many replicas of one system, coordinates only
/// \param csw                     CoordinateSeries writeable abstract
/// \param system_index            System of interest, if coordinates come in a PhaseSpaceSynthesis
/// \param frame_index             Frame of interest, if coordinates come in a CoordinateSeries
/// \param atom_i                  Root of the rotatable bond
/// \param atom_j                  Second atom of the rotatable bond.
/// \param moving_atoms            Atoms branching from atom_j and distal to atom_i.  These will
///                                rotate.
/// \param rotation_angle          The angle about which to rotate atoms indexed in moving_atoms.
///                                The angle is oriented by the right hand rule with one's hand on
///                                atom_i and thumb pointed towards atom_j.
/// \param globalpos_scale_factor  Position scaling factor for coordinates, if the coordinates
///                                originated in a PhaseSpaceSynthesis or CoordinateSeries using
///                                an integer representation
/// \{
template <typename Tcoord, typename Tcalc>
void rotateAboutBond(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, int atom_i, int atom_j,
                     const int* moving_atoms, int moving_atom_count, Tcalc rotation_angle,
                     Tcalc globalpos_scale_factor = 1.0);

void rotateAboutBond(CoordinateFrame *cf, int atom_i, int atom_j,
                     const std::vector<int> &moving_atoms, double rotation_angle);

void rotateAboutBond(CoordinateFrameWriter cfw, int atom_i, int atom_j,
                     const std::vector<int> &moving_atoms, double rotation_angle);

void rotateAboutBond(PhaseSpace *ps, int atom_i, int atom_j, const std::vector<int> &moving_atoms,
                     double rotation_angle);

void rotateAboutBond(PhaseSpaceWriter psw, int atom_i, int atom_j,
                     const std::vector<int> &moving_atoms, double rotation_angle);

template <typename Tcoord, typename Tcalc>
void rotateAboutBond(CoordinateSeries<Tcoord> *cs, int frame_index, int atom_i, int atom_j,
                     const std::vector<int> &moving_atoms, Tcalc rotation_angle);

template <typename Tcoord, typename Tcalc>
void rotateAboutBond(CoordinateSeriesWriter<Tcoord> csw, int frame_index, int atom_i, int atom_j,
                     const int* moving_atoms, int moving_atom_count, Tcalc rotation_angle);

template <typename Tcalc>
void rotateAboutBond(PhaseSpaceSynthesis *psynth, int system_index, int atom_i, int atom_j,
                     const std::vector<int> &moving_atoms, Tcalc rotation_angle);

template <typename Tcalc>
void rotateAboutBond(PsSynthesisWriter psynthw, int system_index, int atom_i, int atom_j,
                     const int* moving_atoms, int moving_atom_count, Tcalc rotation_angle);

template <typename Tcalc>
void rotateAboutBond(Condensate *cdns, int system_index, int atom_i, int atom_j,
                     const std::vector<int> &moving_atoms, Tcalc rotation_angle);
  
template <typename Tcalc>
void rotateAboutBond(CondensateWriter cdnsw, int system_index, int atom_i, int atom_j,
                     const std::vector<int> &moving_atoms, Tcalc rotation_angle);
/// \}

/// \brief Rotate two branches of a chiral center 180 degrees, so as to invert the center.  
///
/// Overloaded:
///   - Take the X, Y, and Z particle coordinates directly, in a templated function that will
///     handle the fixed-precision conversions internally if ndeed
///   - Take a CoordinateFrame or writeable abstract for the coordinates
///   - Take a PhaseSpace object or writeable abstract for the coordinates
///
/// \param xcrd                    Cartesian X coordinates of particles
/// \param ycrd                    Cartesian Y coordinates of particles
/// \param zcrd                    Cartesian Z coordinates of particles
/// \param cf                      Coordinates to manipulate (will be modified upon return)
/// \param cfw                     CoordianteFrame writer abstract
/// \param ps                      PhaseSpace object with coordinates, velocities, and forces (the
///                                coordinates will be modified upon return)
/// \param psw                     PhaseSpace writer abstract
/// \param psynth                  Collection of systems, with coordinates, velocities, and forces
///                                represented in fixed precision
/// \param psynthw                 PhaseSpaceSynthesis writer abstract
/// \param system_index            System of interest, if coordinates come in a PhaseSpaceSynthesis
/// \param frame_index             Frame of interest, if coordinates come in a CoordinateSeries
/// \param center_idx              Index of the chiral center from within the following list
/// \param chiral_centers          List of chiral centers in the molecule
/// \param chiral_protocols        Instructions for how to invert each chiral center
/// \param inversion_groups        List of anchor atom and moving atoms which, combined with the
///                                information in chiral_protocols, set the atoms that move during
///                                a chiral center inversion
/// \param globalpos_scale_factor  Position scaling factor for coordinates, if the coordinates
///                                originated in a PhaseSpaceSynthesis or CoordinateSeries using
///                                an integer representation
/// \{
template <typename Tcoord, typename Tcalc>
void flipChiralCenter(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups,
                      Tcalc globalpos_scale_factor = 1.0);

void flipChiralCenter(CoordinateFrame *cf, int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups);

void flipChiralCenter(CoordinateFrameWriter cfw, int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups);

void flipChiralCenter(PhaseSpace *ps, int center_idx, const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups);

void flipChiralCenter(PhaseSpaceWriter psw, int center_idx, const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups);

template <typename Tcoord, typename Tcalc>
void flipChiralCenter(CoordinateSeries<Tcoord> *cs, int frame_index, int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups);

template <typename Tcoord, typename Tcalc>
void flipChiralCenter(CoordinateSeriesWriter<Tcoord> csw, int frame_index, int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups);

template <typename Tcalc>
void flipChiralCenter(PsSynthesisWriter psynthw, int system_index, int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups);

template <typename Tcalc>
void flipChiralCenter(PhaseSpaceSynthesis *psynth, int system_index, int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups);

template <typename Tcalc>
void flipChiralCenter(Condensate *cdns, int system_index, int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups);

template <typename Tcalc>
void flipChiralCenter(CondensateWriter cdnsw, int system_index, int center_idx,
                      const std::vector<int> &chiral_centers,
                      const std::vector<ChiralInversionProtocol> &chiral_protocols,
                      const std::vector<IsomerPlan> &inversion_groups);
/// \}
  
} // namespace structure
} // namespace stormm

#include "isomerization.tpp"

#endif
