// -*-c++-*-
#ifndef STORMM_ENERGY_ABSTRACTS_H
#define STORMM_ENERGY_ABSTRACTS_H

#include "copyright.h"
#include "Constants/generalized_born.h"
#include "Topology/atomgraph_abstracts.h"

namespace stormm {
namespace energy {

using topology::ImplicitSolventKit;
using topology::ImplicitSolventModel;
using namespace generalized_born_defaults;

template <typename T> struct ImplicitSolventRecipe {

  /// \brief The constructor takes two abstracts / kits, from the AtomGraph and from an
  ///        instantiation of the "neck" GB coefficient tables struct.  The results are combined
  ///        into a single struct, which is essentially a combined abstract but called a "recipe"
  ///        to denote its status as a fusion of multiple abstracts.
  ///
  /// \param isk   Implicit solvent details for the atoms of a specific system, obtained from an
  ///              AtomGraph object (see Topology/atomgraph.h)
  /// \param ngbk  Tables of constants for "neck" Generalized Born implementations--not needed for
  ///              all GB but essential for a complete solution.  Obtained from a
  ///              NeckGeneralizedBornTable object (see Constants/generalized_born.h).
  explicit ImplicitSolventRecipe(const ImplicitSolventKit<T> &isk,
                                 const NeckGeneralizedBornKit<T> &ngbk);

  // All member variables are public, just like any other abstract
  const int natom;                ///< Number of atoms in the system
  const ImplicitSolventModel igb; ///< Flavor of Generalized Born (see
                                  ///<   Topology/atomgraph_enumerators.h)
  const int table_size;           ///< Size of the neck Generalized Born tables (if applicable)
  const T dielectric;             ///< The dielectric constant to use in solvent screening
  const T kappa;                  ///< Inverse of the Debye-Huckel length (based on the dielectic)
  const T gb_offset;              ///< Offset for baseline Generalized Born radii (default 0.9,
                                  ///<   larger for neck GB model II)
  const T gb_neckscale;           ///< Neck function scaling parameter for neck GB
  const T gb_neckcut;             ///< Cutoff for the "neck" GB function
  const int* neck_gb_idx;         ///< Neck GB indices for each atom (if applicable--not for all GB
                                  ///<   models)
  const T* pb_radii;              ///< Poisson-Boltzmann (or Generalized Born) baseline radii
  const T* gb_screen;             ///< Generalized Born screening factors
  const T* gb_alpha;              ///< Generalized Born alpha parameters (obtained from the
                                  ///<   designated model and atomic numbers / elements of each
                                  ///<   atom)
  const T* gb_beta;               ///< Generalized Born beta parameters
  const T* gb_gamma;              ///< Generalized Born gamma parameters
  const T* neck_max_sep;          ///< Neck GB maximum separations table
  const T* neck_max_val;          ///< Neck GB maximum values table
};

} // namespace energy
} // namespace stormm

#include "energy_abstracts.tpp"

#endif
