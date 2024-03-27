// -*-c++-*-
#ifndef STORMM_HONEYCOMB_H
#define STORMM_HONEYCOMB_H

#include "copyright.h"

namespace stormm {
namespace trajectory {

constexpr float default_scaling_resolution = (float)33554432.0;

/// \brief A compartmentalization of a molecule for all types of periodic boxes.  The name is
///        inspired by the shape of a honeycomb, with hexagonal cells tiling throughout space.  The
///        exact layout depends on the slant on the unit cell's alpha angle (alpha rotates in a
///        plane perpendicular to the x axis, along which all of the Honeycomb's pencil-shaped
///        hexagonal prisms are oriented).
class Honeycomb {
public:

private:
  int position_scaling_res;     ///< Scaling resolution (per Angstrom) for positional integers.
                                ///<   The default value of 33554432 guarantees precision of each
                                ///<   atom's coordinates to within +/-1.6 x 10^-8 Angstroms but
                                ///<   only permits evaluations of local systems that will always
                                ///<   fit within a 128 Angstrom cube.
  int force_scaling_res;        ///< Scaling resolution for minor forces.  The default value of
                                ///<   33554432 guarantees precision of each atom's force to within
                                ///<   +/-1.6 x 10^-8 kcal/mol-Angstrom for minor force components,
                                ///<   but overflows into bits of the major force components at
                                ///<   +/-32 kcal/mol-Angstrom (one bit of the minor force
                                ///<   component is sacrificed, but with ten bits to work with in
                                ///<   the major force component, forces in the range
                                ///<   [-16384, 16384) are represented).
  Hybrid<int4> cell_xyz_qljid;  ///< Cartesian coordinates of all atoms, re-imaged to their home
                                ///<   cells (this array is padded for every cell by the object's
                                ///<   preset maximum atom limit for any cell).  The x, y, and z
                                ///<   coordinates are found in the x, y, and z members of the
                                ///<   tuple, scaled by cell_xyz_resolution for the internal units
                                ///<   of Angstroms.  This information is followed immediately by a
                                ///<   bit-packed int indicating the atom's partial charge index
                                ///<   (up to 4194304 unique partial charges may be defined) in the
                                ///<   low 22 bits and the atom's Lennard-Jones type (up to 1024
                                ///<   unique Lennard-Jones atom types may be defined) in the high
                                ///<   10 bits.
  Hybrid<int4> cell_forces;     ///< Forces on all particles.  'Minor' x, y, and z force components
                                ///<   are stored in the x, y, and z members of the tuple, with
                                ///<   'major' forces for each component stored in the w member of
                                ///<   the tuple, for a precision of 41 bits overall.  The order of
                                ///<   forces in this array tracks the coordinates in xyz_qljid,
  Hybrid<int> cell_bin_index;   ///< Indices of the cells in which each atom resides (this refers
                                ///<   to a master cell list, in which multiple independent systems
                                ///<   may all have groups of non-interacting cells)
  Hybrid<int> cell_array_index; ///< Indices of the xyz_qljid and forces arrays at which to find
                                ///<   the positions and forces of each atom in any given system.
                                ///<   Atoms in this list are ordered according to the topologies
                                ///<   of their respective systems, following the system bounds
                                ///<   given in system_bounds.
  Hybrid<int> system_bounds;    ///< The bounds for each system's independent collection of atoms,
                                ///<   relevant to array_index and subsequently to the xyz_qljid
                                ///<   and forces arrays that 
  Hybrid<int2> atom_x_position; ///< Each atom's Cartesian x position, non-reimaged
  Hybrid<int2> atom_y_position; ///< Each atom's Cartesian y position, non-reimaged
  Hybrid<int2> atom_z_position; ///< Each atom's Cartesian z position, non-reimaged 
  Hybrid<int4> atom_forces;     ///< Each atom's Cartesian x, y, and z forces, with spillover from
                                ///<   the 'minor' quantities in the x, y, and z member variables
                                ///<   to the 'major' quantities packed into w

  // The following member variables are relevant for 
  int q_cell_atom_limit;        ///< The maximum number of partial-charge bearing particles that
                                ///<   any one cell can hold
  int lj_cell_atom_limit;       ///< The maximum number of Lennard-Jones bearing particles that any
                                ///<   one cell can hold
  Hybrid<int> cell_sizes;       ///< Cell sizes for each system, the distances between cell centers
                                ///<   given various translations (this indicates the amounts by
                                ///<   which to translate each particle when it exits one cell and
                                ///<   enters another--the unit cell dimensions are necessarily a
                                ///<   product of the individual cells' sizes).

  // The following information gives detailed instructions as to how to accomplish a set of work
  // units for this Honeycomb object.
  int work_unit_count;          ///< Number of work units in all systems comprised by this object

  /// The work unit intake tells the pair list building and nonbonded execution kernels which
  /// cells' atoms to take in.  Each work unit is outlined by a series of 128 numbers.  Up to 56
  /// cells may be part of the intake.  The numbers describe the work unit as follows:
  /// eight cells long with the first eight integers of every 64 indicating, in order:
  ///   (0): overall size of the work unit's cell contents (up to 56) in the low 6 bits and
  ///     reference level of the stack in the next lowest 4 bits.  Atoms from the cells will be
  ///     read into a local buffer, then subjected to the translations detailed in the subsequent
  ///     seven int values (see the next bullet points).  Finally, the x coordinates of atoms in
  ///     all cells will be translated to the reference frame of the center layer in the stack,
  ///     putting all coordinates into a single frame of reference for the purposes of computing
  ///     non-bonded interactions.  The high 22 bits offer a baseline index for the cell grid's
  ///     exclusion pair list buffers (multiply by 1024 to get the starting index).
  ///   (1): whether to translate pencil (A).  This is a bit-packed integer with the low six bits
  ///     indicating whether to translate the first cell read for the pencil by one box length in
  ///     -a, whether to translate the last cell read for the pencil one box length in +a, whether
  ///     to translate the entire pencil one box length in -b, then in +b, then in -c, then in +c
  ///     (the third and fourth bits are mutually exclusive, as are the fifth and sixth).
  ///   (2-7): similar to the above for pencil (B), then (C), (D), (E), (F), and (G).
  ///   (8-63): integer addresses of the cells to read
  ///   (64-?): lists of atom indices to test for particle-particle interactions, numbered
  ///     according to the local collection of atoms.  Up to sixteen sets of atom series, senders
  ///     [ atom_i0 ... atom_if ] vs receivers [ atom_j0 ... atom_jf, atom_k0 ... atom_kf,
  ///       atom_l0 ... atom_lf, atom_m0 ... atom_mf ] may be specified.  If there is only one
  ///     set of receiver atoms and they comprise the sending atoms this is a special case
  ///     and pairs will not be double-counted.
  Hybrid<int> work_unit_intake;

  /// The pairlist for a local frame is a complex thing, but requires a complete list of all
  /// interacting cells so that a series of unsigned integers can be laid out with bit masks
  /// indicating which cell-to-cell interactions involve exclusions.  The exclusion lists apply
  /// only to van-der Waals interactions: electrostatics can always be calculated and then,
  /// retroactively, removed.  In a periodic unit cell, there must be 47 cell-to-cell pair list
  /// buffers for every Honeycomb cell (each pencil is subdivided into a series of cells thick
  /// enough that atoms in one layer of cells all across the b/c plane cannot interact with others
  /// more than one plane of cells away).
  Hybrid<uint> pairlist_buffer;
};

} // namespace trajectory
} // namespace stormm

#endif
