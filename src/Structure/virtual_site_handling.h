// -*-c++-*-
#ifndef STORMM_VIRTUAL_SITE_PLACEMENT_H
#define STORMM_VIRTUAL_SITE_PLACEMENT_H

#include <cmath>
#include "copyright.h"
#include "Math/vector_ops.h"
#include "local_arrangement.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/phasespace.h"

namespace stormm {
namespace structure {

using stmath::dot;
using stmath::project;
using stmath::crossProduct;
using topology::AtomGraph;
using topology::UnitCellType;
using topology::VirtualSiteKind;
using topology::VirtualSiteKit;
using trajectory::CoordinateFrame;
using trajectory::CoordinateFrameWriter;
using trajectory::PhaseSpace;
using trajectory::PhaseSpaceWriter;

/// \brief Reference function for placing virtual sites, using double-precision math throughout.
///
/// Overloaded:
///   - Take templated pointers to raw coordinate and box transformation matrices, plus a virtual
///     sites abstract from the topology matching the precision of the transformation matrices
///   - Take modifiable forms of any of the coordinate objects or their abstracts (most will imply
///     a specific coordinate representation)
///
/// \param xcrd               Cartesian X coordinates of all particles
/// \param ycrd               Cartesian Y coordinates of all particles
/// \param zcrd               Cartesian Z coordinates of all particles
/// \param umat               Transformation matrix taking coordinates into fractional space
/// \param invu               Transformation matrix taking coordinates back to real space
/// \param unit_cell          The system's unit cell type
/// \param ps                 Coordinates of the system as a mutable PhaseSpace object
/// \param cf                 Coordinates of the system as a mutable CoordinateFrame object
/// \param vsk                Virtual sites details abstracted from the original topology
/// \param ag                 System topology containing virtual site specifications        
/// \param gpos_scale_factor  Scaling factor to convert real coordinates into a fixed-precision
///                           representation (the inverse is computed internally)
/// \{
template <typename Tcoord, typename Tcalc>
void placeVirtualSites(Tcoord* xcrd, Tcoord* ycrd, Tcoord* zcrd, const double* umat,
                       const double* invu, const UnitCellType unit_cell,
                       const VirtualSiteKit<Tcalc> &vsk, Tcalc gpos_scale_factor = 1.0);

void placeVirtualSites(PhaseSpace *ps, const AtomGraph &ag);

void placeVirtualSites(PhaseSpace *ps, const AtomGraph *ag);

void placeVirtualSites(CoordinateFrame *cf, const AtomGraph &ag);

void placeVirtualSites(CoordinateFrame *cf, const AtomGraph *ag);

template <typename Tcalc>
void placeVirtualSites(PhaseSpaceWriter psw, const VirtualSiteKit<Tcalc> vsk);

template <typename Tcalc>
void placeVirtualSites(CoordinateFrameWriter cfw, const VirtualSiteKit<Tcalc> vsk);
/// \}

/// \brief Transmit forces on virtual sites to their frame atoms, using double-precision math
///        throughout.
///
/// Overloaded:
///   - Take templated pointers to raw coordinate and box transformation matrices, plus a virtual
///     sites abstract from the topology matching the precision of the transformation matrices
///   - Take modifiable forms of appropriate coordinate objects or their abstracts (each will imply
///     a specific coordinate representation)
///
/// \param xcrd                Cartesian X coordinates of all particles
/// \param ycrd                Cartesian Y coordinates of all particles
/// \param zcrd                Cartesian Z coordinates of all particles
/// \param xfrc                Cartesian X forces acting on all particles
/// \param yfrc                Cartesian Y forces acting on all particles
/// \param zfrc                Cartesian Z forces acting on all particles
/// \param umat                Transformation matrix taking coordinates into fractional space
/// \param invu                Transformation matrix taking coordinates back to real space
/// \param unit_cell           The system's unit cell type
/// \param ps                  Coordinates of the system as a mutable PhaseSpace object
/// \param cf                  Coordinates of the system as a mutable CoordinateFrame object
/// \param vsk                 Virtual sites details abstracted from the original topology
/// \param ag                  System topology containing virtual site specifications
/// \param gpos_scale_factor   Scaling factor to convert real coordinates into a fixed-precision
///                            representation (the inverse is computed internally)
/// \param force_scale_factor  Scaling factor to convert real-valued force components into a
///                            fixed-precision representation (the inverse is computed internally)
/// \{
template <typename Tcoord, typename Tforce, typename Tcalc>
void transmitVirtualSiteForces(const Tcoord* xcrd, const Tcoord* ycrd, const Tcoord* zcrd,
                               Tforce* xfrc, Tforce* yfrc, Tforce* zfrc, const double* umat,
                               const double* invu, const UnitCellType unit_cell,
                               const VirtualSiteKit<Tcalc> &vsk, Tcalc gpos_scale_factor = 1.0,
                               Tcalc force_scale_factor = 1.0);

void transmitVirtualSiteForces(PhaseSpace *ps, const AtomGraph &ag);

void transmitVirtualSiteForces(PhaseSpace *ps, const AtomGraph *ag);

template <typename Tcalc>
void transmitVirtualSiteForces(PhaseSpace psw, VirtualSiteKit<Tcalc> vsk);
/// \}
  
} // namespace structure
} // namespace stormm

#include "virtual_site_handling.tpp"

#endif
