// -*-c++-*-
#ifndef STORMM_ATOMGRAPH_ANALYSIS_H
#define STORMM_ATOMGRAPH_ANALYSIS_H

#include <cmath>
#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/behavior.h"
#include "Constants/scaling.h"
#include "Potential/energy_enumerators.h"
#include "Reporting/error_format.h"
#include "atomgraph.h"
#include "atomgraph_abstracts.h"

namespace stormm {
namespace topology {

using card::Hybrid;
using constants::ExceptionResponse;
using energy::VdwCombiningRule;

/// \brief Apply a series of checks to identify the water model in use within a topology.  If there
///        is no water, the model will be set to NONE.
///
/// \param ag  The topology to analyze
WaterModel identifyWaterModel(const AtomGraph &ag);

/// \brief Identify the virtual site types present in a topology.
///
/// Overloaded:
///   - Take a const pointer to the topology
///   - Take a const reference to the topology
///   - Operate on the vector of virtual site types with a trusted indication of its length
/// 
/// \param ag        The topology of interest
/// \param vs_types  Array containing the frame types of all virtual sites in the system
/// \param nsite     The number of virtual sites in the system
/// \{
std::string listVirtualSiteFrameTypes(const AtomGraph &ag);

std::string listVirtualSiteFrameTypes(const AtomGraph *ag);

std::string listVirtualSiteFrameTypes(const int* vs_types, int nsite);
/// \}

/// \brief Color an atom graph based on connectivity.  Start at atom I, first color atom J, and do
///        not retrace any steps.  This will initialize and mark a pre-allocated array such that
///        connected atoms are marked TRUE.
///
/// \param nbk          Non-bonded details (needed for exclusion lists) of the original topology
/// \param atom_i       The root atom of the rotatable bond
/// \param atom_j       The second atom of the rotatable bond.  This can be taken as the center of
///                     coordinates for rotating the local system, if convenient.
/// \param marked       Array indicating whether the atoms have been colored, sized for the number
///                     of atoms in the system
/// \param ring_report  Flag to indicate that a ring was completed as a result of the coloring
///                     (returned)
int colorConnectivity(const NonbondedKit<double> &nbk, const ChemicalDetailsKit &cdk,
                      int atom_i, int atom_j, std::vector<uint> *marked, bool *ring_report);
  
/// \brief The essential work of the selectRotatingAtoms() function below, abstracted to allow
///        more flexibility in the way the function is called.
///
/// \param nbk       Non-bonded details (needed for exclusion lists) of the original topology
/// \param cdk       Chemical details (here, atom names) of the original topology
/// \param atom_i    The root atom of the rotatable bond
/// \param atom_j    The second atom of the rotatable bond.  This can be taken as the center of
///                  coordinates for rotating the local system, if convenient.
/// \param filename  The name of the original topology file (for error tracing purposes)
std::vector<int> mapRotatingGroup(const NonbondedKit<double> &nbk, const ChemicalDetailsKit &cdk,
                                  const int atom_i, const int atom_j, const std::string &filename);
  
/// \brief Select atoms for rotation, given a rotatable bond found in a molecule's chemical
///        features.
///
/// Overloaded:
///   - Take a const reference to an AtomGraph object
///   - Take a const pointer to an AtomGraph object
///
/// \param ag      System topology
/// \param atom_i  The root atom of the rotatable bond
/// \param atom_j  The second atom of the rotatable bond.  This can be taken as the center of
///                coordinates for rotating the local system, if convenient.
/// \{
std::vector<int> selectRotatingAtoms(const AtomGraph &ag, int atom_i, int atom_j);
std::vector<int> selectRotatingAtoms(const AtomGraph *ag, int atom_i, int atom_j);
/// \}

/// \brief Compute the number of Lennard-Jones atom types given two vectors of a particular length,
///        assuming them to be square matrices of the A and B coefficients.
///
/// \param length_a  The number of values in the matrix of A coefficients
/// \param length_b  The number of values in the matrix of B coefficients
/// \param caller    Name of the calling function
int inferLennardJonesTypeCount(const int length_a, const int length_b,
                               const char* caller = nullptr);

/// \brief Determine the Lennard-Jones combining rule in effect.  The function will accept either
///        single- or double-precision data, but internally it uses double-precision calculations.
///
/// Overloaded:
///   - Provide the A and B coefficient arrays by C-style arrays with trusted length
///   - Provide the A and B coefficient arrays as Standard Template Library vectors
///   - Provide the A and B coefficient arrays as Hybrid objects
///   - Provide a topology and extract the Lennard-Jones parameters from it
///
/// \param lja             Array of Lennard-Jones A coefficients
/// \param ljb             Array of Lennard-Jones B coefficients
/// \param lj_type_count   The number of Lennard-Jones types, indicating the trusted lengths of
///                        lja and ljb by the number's square
/// \param policy          Protocol in the event that there are only one (fewer) atom types
/// \param seek_prevalent  A rare boolean input that, if set to TRUE, will avoid declaring a
///                        parameter set "NBFIX" and instead declare it to be "GEOMETRIC" or
///                        "LORENTZ_BERTHELOT" depending on whether more cases of either combining
///                        rule can be found.  A parameter set might still be decalred "NBFIX" if
///                        no instances of the other combining rules can describe any off-diagonal
///                        interactions in the matrix.
/// \{
template <typename T>
VdwCombiningRule inferCombiningRule(const T* lja, const T* ljb, int lj_type_count,
                                    ExceptionResponse policy = ExceptionResponse::WARN,
                                    bool seek_prevalent = false);
  
template <typename T>
VdwCombiningRule inferCombiningRule(const std::vector<T> &lja, const std::vector<T> &ljb,
                                    ExceptionResponse policy = ExceptionResponse::WARN,
                                    bool seek_prevalent = false);

template <typename T>
VdwCombiningRule inferCombiningRule(const Hybrid<T> &lja, const Hybrid<T> &ljb,
                                    ExceptionResponse policy = ExceptionResponse::WARN,
                                    bool seek_prevalent = false);

VdwCombiningRule inferCombiningRule(const AtomGraph *ag,
                                    ExceptionResponse policy = ExceptionResponse::WARN,
                                    bool seek_prevalent = false);

VdwCombiningRule inferCombiningRule(const AtomGraph &ag,
                                    ExceptionResponse policy = ExceptionResponse::WARN,
                                    bool seek_prevalent = false);
/// \}

/// \brief Compute the number of degrees of freedom in some partition of a topology, delineated
///        by atom indices.  The subdivision could be protein and water, ligand, or specific parts
///        of a single biomolecule (although a warning will be raised if the constraints do not fit
///        precisely within the partitions as defined).
///
/// Overloaded:
///   - Provide the topology by const pointer
///   - Provide the topology by const reference
///
/// \param ag               Topology of the system, containing the atom count and constraint
///                         information
/// \param low_atom_index   The lower limit of atoms in the subset of interest
/// \param high_atom_index  The upper limit of atoms in the subset of interest
/// \{
int getConstrainedDegreesOfFreedom(const AtomGraph *ag, int low_atom_index = 0,
                                   int high_atom_index = 0);

int getConstrainedDegreesOfFreedom(const AtomGraph &ag, int low_atom_index = 0,
                                   int high_atom_index = 0);
/// \}

} // namespace topology
} // namespace stormm

#include "atomgraph_analysis.tpp"

#endif
