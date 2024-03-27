// -*-c++-*-
#ifndef CONFORMER_SETUP_H
#define CONFORMER_SETUP_H

#include "../../../src/Chemistry/atommask.h"
#include "../../../src/Chemistry/chemistry_enumerators.h"
#include "../../../src/Chemistry/chemical_features.h"
#include "../../../src/Constants/behavior.h"
#include "../../../src/Math/tickcounter.h"
#include "../../../src/MoleculeFormat/mdlmol.h"
#include "../../../src/Namelists/nml_conformer.h"
#include "../../../src/Namelists/nml_minimize.h"
#include "../../../src/Namelists/user_settings.h"
#include "../../../src/Potential/static_exclusionmask.h"
#include "../../../src/Random/random.h"
#include "../../../src/Restraints/restraint_apparatus.h"
#include "../../../src/Synthesis/phasespace_synthesis.h"
#include "../../../src/Synthesis/synthesis_cache_map.h"
#include "../../../src/Synthesis/synthesis_permutor.h"
#include "../../../src/Synthesis/systemcache.h"
#include "../../../src/Topology/atomgraph.h"
#include "../../../src/Topology/atomgraph_abstracts.h"
#include "../../../src/Trajectory/coordinateframe.h"
#include "../../../src/Trajectory/phasespace.h"
#include "../../../src/UnitTesting/stopwatch.h"

namespace conf_app {
namespace setup {

using stormm::chemistry::AtomMask;
using stormm::chemistry::ChemicalFeatures;
using stormm::chemistry::ChiralInversionProtocol;
using stormm::chemistry::IsomerPlan;
using stormm::constants::ExceptionResponse;
using stormm::energy::StaticExclusionMask;
using stormm::namelist::ConformerControls;
using stormm::namelist::MinimizeControls;
using stormm::namelist::UserSettings;
using stormm::restraints::RestraintApparatus;
using stormm::random::Xoshiro256ppGenerator;
using stormm::stmath::TickCounter;
using stormm::structure::MdlMol;
using stormm::synthesis::PhaseSpaceSynthesis;
using stormm::synthesis::SynthesisCacheMap;
using stormm::synthesis::SynthesisPermutor;
using stormm::synthesis::SystemCache;
using stormm::testing::StopWatch;
using stormm::topology::AtomGraph;
using stormm::topology::NonbondedKit;
using stormm::trajectory::PhaseSpace;
using stormm::trajectory::CoordinateFrame;

/// \brief Timings strings to track various setup procedures
/// \{
constexpr char tm_input_parsing[] = "Input parsing";
constexpr char tm_work_unit_building[] = "Work unit building";
constexpr char tm_coordinate_expansion[] = "Coordinate expansion";
constexpr char tm_energy_minimization[] = "Energy minimization";
constexpr char tm_conformer_selection[] = "Conformer selection";
/// \}
  
/// \brief Enumerate coarse-grained conformer sampling approaches: full if the combinatorial
///        permutations are few enough, randomized without replacement if the combinatorial
///        permutations are up to 32 times larger, and randomized with replacement if the
///        combinatorial permutations are huge.
enum class SamplingStrategy {
  FULL,     ///< Sample all permutations of conformers
  LIMITED,  ///< Randomly sample conformers without replacement
  SPARSE    ///< Randomly sample conformers with no regard to whether replacement could occur
};
  
/// \brief Get the core mask for one of the topologies in the cache of systems.  This will compute
///        the mask for one example of the systems using the topology, and as such should not
///        be used with core masks defined by coordinate-related details.
///
/// \param conf_input   The &conformer namelist user input
/// \param sdf_example  An MDL MOL structure containing coordinates and possibly a data item
///                     defining topology-specific core mask information
/// \param ps           An example of the coordinates, used if the MDL MOL entry is not available
/// \param ag           The topology of interest
/// \param chemfe       Chemical features computed for the topology of interest
/// \param policy       Indicate what to do if errors are encountered when making the mask (passed
///                     down from command-line input)
AtomMask getCoreMask(const ConformerControls &conf_input, const MdlMol &sdf_example,
                     const PhaseSpace &ps, const AtomGraph *ag, const ChemicalFeatures &chemfe,
                     const ExceptionResponse policy);

/// \brief Add implicit solvent models to topologies.  Set the position-restrained cores of each
///        proto-conformer in the system cache and return the list of such masks.
///
/// \param ui  User input settings, obtained from the input deck
/// \param sc  Cache of topologies and initial structures.  The topologies are modified temporarily
///            to immobilize certain atoms when determining rotatable bond groups, then returned to
///            their original mobility states.  Implicit solvent conditions are added to all
///            topologies.
/// \param tm  Timer to record the wall time spent on various setup procedures
void setGenerativeConditions(const UserSettings &ui, SystemCache *sc, StopWatch *tm);

/// \brief Return an estimate of the number of all permutations covering isomerizations involving
///        the same atoms or atoms that are bonded to one another.
///
/// \param limits       The number of possible settings for each permutation
/// \param isomerizers  Groups of atoms affecting different conformers and isomers, including bond
///                     rotations
/// \param ag           Topology for the molecule of interest
double computeLocalPermutations(const std::vector<int> &limits,
                                const std::vector<IsomerPlan> &isomerizers, const AtomGraph *ag);

/// \brief Create the sandbox for working with conformers.  All coordinates for each system in this
///        workspace will be initialized from the systems cache, to then be modified.
///
/// \param sc       The cache of systems, along with ChemicalFeatures for each of them
/// \param confcon  Conformer namelist controls from user input
/// \param mincon   Minimization controls, used to convey clash criteria
/// \param xrs      Random number generator for making permutations with LIGHT or HEAVY sampling of
///                 the user-specified systems
/// \param scmap    Map between the systems cache built from user input and the synthesis
///                 constructed inside this function (filled and returned)
/// \param syper    Object for guiding the permutations of each unique system in the user-specified
///                 cache of systems
/// \param tm       Timer to record the wall time spent on various setup procedures
PhaseSpaceSynthesis buildSamplingWorkspace(const SystemCache &sc, const ConformerControls &confcon,
                                           const MinimizeControls &mincon,
                                           Xoshiro256ppGenerator *xrs, SynthesisCacheMap *scmap,
                                           SynthesisPermutor *syper, StopWatch *tm);

/// \brief Build a list of restraint apparatus pointers for the synthesis of all system replicas.
///
/// \param scmap    Map between the systems cache built from user input and the synthesis
///                 constructed inside this function (provides a pointer to the original systems
///                 cache built from user input)
const std::vector<RestraintApparatus*> buildReplicaRestraints(const SynthesisCacheMap &scmap);
  
/// \brief Expand the initial list of systems into a complete list of initial states for the
///        population of conformers which conformer.omni will minimize in search of the
///        lowest-energy states.  This is done by parsing each topology into its chemical
///        details and then enumerating the rotatable bonds, cis-trans isomers, and chiral centers
///        which it could sample.
///
/// \param poly_ps          Pre-allocated synthesis of coordinates (modified and returned)
/// \param replica_counts   Numbers of replicas anticipated for each system
/// \param ui               User input settings, obtained from the input deck
/// \param sc               Cache of topologies and initial structures.  The list of topologies will
///                         not be expanded by this procedure, but the list of structures will
///                         undergo a radical expansion.
/// \param xrs              Random number generator (the master generator will guide the CPU-based,
///                         coarse-grained conformer selection)
/// \param tm               Timer to record the wall time spent on various setup procedures
void expandConformers(PhaseSpaceSynthesis *poly_ps, const std::vector<int> &replica_counts,
                      const UserSettings &ui, const SystemCache &sc, Xoshiro256ppGenerator *xrs,
                      StopWatch *tm);

} // namespace setup
} // namespace conf_app

#endif
