// -*-c++-*-
#ifndef STORMM_TOPOLOGY_BOUNDS_CHECKS_H
#define STORMM_TOPOLOGY_BOUNDS_CHECKS_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/scaling.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/sorting_enumerators.h"
#include "Math/vector_ops.h"
#include "Potential/energy_enumerators.h"
#include "Reporting/error_format.h"
#include "atomgraph.h"
#include "atomgraph_abstracts.h"
#include "topology_limits.h"

namespace stormm {
namespace topology {

using card::Hybrid;
using card::HybridTargetLevel;
using data_types::getStormmScalarTypeName;
using data_types::getStormmHpcVectorTypeName;
using data_types::isScalarType;
using data_types::isHpcVectorType;
using energy::NonbondedTheme;
using energy::VdwCombiningRule;
using stmath::UniqueValueHandling;
using stmath::reduceUniqueValues;

/// \brief Extract a series of numbers from a bounded list, based on an index.  Return the result
///        as a std::vector of integers.  If data from multiple bins is requested, offer an option
///        to sort and prune the data for unique values.
///
/// Overloaded:
///   - Accept a single bin index
///   - Accept multiple bin indices
///   - Accept a pre-allocated result vector with one or more bin indices
///
/// \param result     The vector to pack the outcome into (should be pre-allocated, but will be
///                   resized as necessary)
/// \param va         The original list of entries with demarcations given in va_bounds
/// \param va_bounds  Bounds array for va
/// \param index      Index of the bin of interest (checked against the size of va_bounds)
/// \param indices    Vector of indices for the bins of interest (each will be checked against the
///                   size of va_bounds)
/// \param filter     Indication that the result should be reduced to unique values only (applies
///                   to cases with multiple collection bin indices)
/// \{
template <typename T> std::vector<T> extractBoundedListEntries(const std::vector<T> &va,
                                                               const std::vector<int> &va_bounds,
                                                               int index);

template <typename T>
void extractBoundedListEntries(std::vector<T> *result, const std::vector<T> &va,
                               const std::vector<int> &va_bounds, int index);

template <typename T> std::vector<T>
extractBoundedListEntries(const std::vector<T> &va, const std::vector<int> &va_bounds,
                          const std::vector<int> &indices,
                          UniqueValueHandling filter = UniqueValueHandling::UNIQUE_VALUES_ONLY);

template <typename T> void
extractBoundedListEntries(std::vector<T> *result, const std::vector<T> &va,
                          const std::vector<int> &va_bounds, const std::vector<int> &indices,
                          UniqueValueHandling filter = UniqueValueHandling::UNIQUE_VALUES_ONLY);
/// \}

/// \brief Obtain parameters in either single- or double-precision real data formats.
///
/// \param item        Double-precision representation of the data
/// \param sp_item     Single-precision representation of the data
/// \param tier        Indicator of whether to pull data from the host (CPU) or device (GPU)
/// \param low_index   Lower limit of the data to obtain
/// \param high_index  Upper limit of the data to obtain
/// \param caller      Optional name of the calling object (for error reporting purposes)
/// \param method      Optional member function within the calling object (for error reporting)
template <typename T>
std::vector<T> getRealParameters(const Hybrid<double> &item, const Hybrid<float> &sp_item,
                                 HybridTargetLevel tier, int low_index, int high_index,
                                 const char* caller = nullptr, const char* method = nullptr);
  
/// \brief Make a human-readable list of atoms based on a topology and a vector of indices,
///        providing atom numbers and name as well as residue numbers and names.  This is useful,
///        in particular, when printing error messages to the user.
///
/// \param atom_list  Topological indices for the atoms of interest
/// \param cdk        Chemical details abstract for the topology of interest
std::string writeAtomList(const std::vector<int> &atom_list, const ChemicalDetailsKit &cdk);

/// \brief Determine whether two atoms are bonded, based on non-bonded 1:2 exclusions.
///
/// Overloaded:
///   - Pass the topology's non-bonded abstract
///   - Pass the topology itself, as a const pointer or as a const reference
///
/// \param ag      Topology for the molecule of interest
/// \param nbk     Non-bonded abstract from ag
/// \param atom_i  Test atom I for bonding to atom J
/// \param atom_j  Test atom J for bonding to atom I
/// \{
bool isBonded(const NonbondedKit<double> &nbk, int atom_i, int atom_j);
bool isBonded(const AtomGraph &ag, int atom_i, int atom_j);
bool isBonded(const AtomGraph *ag, int atom_i, int atom_j);
/// \}

/// \brief Locate a match for a topology pointer among a list of topology pointers.  If the pointer
///        cannot be found, look for a topology originating from the same file.  Return the index
///        the array containing the matching topology.  If no match is found, return the length of
///        the array of topologies.
///
/// Overloaded:
///   - Match a pointer to a topology against an array of topology pointers
///   - Match a const reference to a topology against an array of topology pointers
///   - Match a pointer to a topology against an array of topologies
///   - Match a const reference to a topology against an array of topologies
///
/// \param query_ag   The topology pointer to match against the array
/// \param repo       The array of topology pointers
/// \{
int matchTopology(const AtomGraph *query_ag, const std::vector<AtomGraph*> &repo);
int matchTopology(const AtomGraph &query_ag, const std::vector<AtomGraph*> &repo);
int matchTopology(const AtomGraph *query_ag, const std::vector<AtomGraph> &repo);
int matchTopology(const AtomGraph &query_ag, const std::vector<AtomGraph> &repo);
/// \}

/// \brief Determine whether an atom in a topology has significant van-der Waals properties.  This
///        will check all interactions of the atom's type with other types, if necessary.
///
/// Overloaded:
///   - Provide a pointer or reference to the original abstract
///   - Provide the non-bonded parameter kit
///   - Provide the detected Lennard-Jones combining rule for optimization
///
/// \param ag          The topology to analyze
/// \param nbk         The non-bonded parameters obtained from the topology
/// \param lj_rule     The detected (and trusted) Lennard-Jones combining rule
/// \param atom_index  The index of the atom of interest
/// \{
template <typename T>
bool hasVdwProperties(const NonbondedKit<T> &nbk, int atom_index,
                      VdwCombiningRule lj_rule = VdwCombiningRule::NBFIX);

bool hasVdwProperties(const AtomGraph *ag, int atom_index,
                      VdwCombiningRule lj_rule = VdwCombiningRule::NBFIX);

bool hasVdwProperties(const AtomGraph &ag, int atom_index,
                      VdwCombiningRule lj_rule = VdwCombiningRule::NBFIX);
/// \}

/// \brief Determine whether an atom is relevant to a particular kind of non-bonded interactions.
///
/// \param nbk         Non-bonded parameters obtained from the topology
/// \param atom_index  Index of the atom of interest
/// \param theme       Indicate the kind of non-bonded properties to seek out in the atom (or
///                    virtual site)s
/// \param lj_rule     The perceived Lennard-Jones combining rule (important as NBFix combining
///                    rules require that all possible combinations of the atom with other types be
///                    checked, whereas Lorentz-Berthelot and geometric combining rules require
///                    only that the atom have significant self interations)
template <typename T>
bool hasRelevantProperties(const NonbondedKit<T> &nbk, int atom_index, NonbondedTheme theme,
                           VdwCombiningRule lj_rule = VdwCombiningRule::NBFIX);

} // namespace topology
} // namespace stormm

#include "topology_util.tpp"

#endif
