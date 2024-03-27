// -*-c++-*-
#ifndef STORMM_SYNTHESIS_PERMUTOR_H
#define STORMM_SYNTHESIS_PERMUTOR_H

#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Chemistry/chemical_features.h"
#include "Chemistry/chemistry_enumerators.h"
#include "Constants/fixed_precision.h"
#include "Math/tickcounter.h"
#include "Namelists/nml_conformer.h"
#include "Random/random.h"
#include "Structure/clash_detection.h"
#include "Structure/isomerization.h"
#include "Structure/local_arrangement.h"
#include "Structure/structure_enumerators.h"
#include "Topology/atomgraph.h"
#include "atomgraph_synthesis.h"
#include "phasespace_synthesis.h"
#include "synthesis_enumerators.h"

namespace stormm {
namespace synthesis {

using card::Hybrid;
using card::HybridTargetLevel;
using chemistry::ChemicalFeatures;
using chemistry::ChiralInversionProtocol;
using chemistry::ChiralOrientation;
using chemistry::getChiralOrientation;
using chemistry::ConformationEdit;
using chemistry::CoupledEdit;
using chemistry::IsomerPlan;
using namelist::ConformerControls;
using namelist::default_conf_clash_pairs;
using namelist::default_conf_max_seeding_attempts;
using numerics::default_globalpos_scale_bits;
using numerics::default_localpos_scale_bits;
using numerics::default_velocity_scale_bits;
using numerics::default_force_scale_bits;
using random::Xoshiro256ppGenerator;
using stmath::TickCounter;
using structure::ClashReport;
using structure::detectClash;
using structure::dihedralAngle;
using structure::flipChiralCenter;
using structure::rotateAboutBond;
using structure::SamplingIntensity;
using topology::AtomGraph;

/// \brief Read-only abstract for the SynthesisPermutor.  Access to all of the atom groups,
///        critical measurements, and settings for mutable components is contained herein.
template <typename T> struct SyPermutorKit {
public:

  /// \brief As with other abstracts, the constructor takes input arguments for every member
  ///        variable.
  explicit SyPermutorKit(int system_count_in, int perm_map_count_in,
                         const int* perm_map_idx_in, const int2* perm_elements_in,
                         const int* perm_element_bounds_in, const int* system_settings_in,
                         const int* system_settings_limits_in, const int* rot_grp_atoms_in,
                         const int* rot_grp_bounds_in, const int* prm_rot_grp_bounds_in,
                         const int* ctx_grp_atoms_in, const int* ctx_grp_bounds_in,
                         const int* prm_ctx_grp_bounds_in, const int* inv_grp_atoms_in,
                         const int* inv_grp_bounds_in, const int* prm_inv_grp_bounds_in,
                         const int* chiral_atoms_in, const int* chiral_protocols_in,
                         const int4* rot_bond_markers_in, const int4* ctx_bond_markers_in,
                         const int4* chiral_markers_in, const T* rot_bond_settings_in,
                         const T* ctx_bond_settings_in, const int* rot_bond_settings_bounds_in,
                         const int* ctx_bond_settings_bounds_in,
                         const int* chiral_settings_in, const int* chiral_settings_bounds_in);

  // The member variables store size constants for arrays (the numbers of systems in the synthesis
  // and the number of unique prmologies (translating to the number of unique permutor maps)
  // involved in manipulating each structure of the synthesis.
  int system_count;                     ///< Number of systems in the synthesis
  int perm_map_count;                   ///< Number of permutation maps (unique topologies in the
                                        ///<   synthesis)
  const int* perm_map_idx;              ///< Permutation map indices for every system in the
                                        ///<   synthesis
  const int2* perm_elements;            ///< Array identities and array indices of each element in
                                        ///<   each permutor map.  This information directs threads
                                        ///<   on where to search in rot_bond_settings,
                                        ///<   ctx_bond_settings, or chiral_settings after reading
                                        ///<   limits set forth in those arrays' bounds.
  const int* perm_element_bounds;       ///< Bounds array for perm_elements
  const int* system_settings;           ///< Array of settings indices for all elements of all
                                        ///<   systems in the synthesis.
  const int* system_settings_limits;    ///< Array of limits for system_settings--the ith element
                                        ///<   of system_settings must always fall within the range
                                        ///<   [ 0, system_settings_limits[i] ).
  const int* rot_grp_atoms;             ///< Rotatable bond group associated atoms
  const int* rot_grp_bounds;            ///< Bounds array for rot_grp_atoms
  const int* prm_rot_grp_bounds;        ///< Bounds array for rot_grp_bounds--giving the limits of
                                        ///<   rot_grp_atoms limits within each permutor map
  const int* ctx_grp_atoms;             ///< Cis-trans isomeric bond group atoms
  const int* ctx_grp_bounds;            ///< Bounds array for ctx_grp_atoms
  const int* prm_ctx_grp_bounds;        ///< Bounds array for ctx_grp_bounds--giving the limits of
                                        ///<   ctx_grp_atoms limits within each permutor map
  const int* inv_grp_atoms;             ///< Chiral center inversion group atoms
  const int* inv_grp_bounds;            ///< Bounds array for inv_grp_atoms
  const int* prm_inv_grp_bounds;        ///< Bounds array for inv_grp_bounds--giving the limits of
                                        ///<   inv_grp_atoms limits within each permutor map
  const int* chiral_atoms;              ///< Concatenated lists of chiral atoms in all systems,
                                        ///<   also bounded by prm_inv_grp_bounds
  const int* chiral_protocols;          ///< Protocols for inverting each chiral center in each
                                        ///<   system, also bounded by prm_inv_grp_bounds
  const int4* rot_bond_markers;         ///< Marker atoms for each rotatable bond, used to define
                                        ///<   the angle about that bond
  const int4* ctx_bond_markers;         ///< Marker atoms for each cis-trans isomeric bond, used to
                                        ///<   define the angle about that bond
  const int4* chiral_markers;           ///< Marker atoms for each chiral center, giving chiral arm
                                        ///<   priorities in the order .x = 0, .y = 3, .z = 2, and
                                        ///<   .w = 1.  See the ChemicalFeatures object
                                        ///<   documentation for more details.
  const T* rot_bond_settings;           ///< Specific values that each rotatable bond can take on
  const T* ctx_bond_settings;           ///< Specific values that each cis-trans isomeric bond can
                                        ///<   take on
  const int* chiral_settings;           ///< Specific values that each chiral center can take on
  const int* rot_bond_settings_bounds;  ///< Bounds array for rot_bond_settings
  const int* ctx_bond_settings_bounds;  ///< Bounds array for ctx_bond_settings
  const int* chiral_settings_bounds;    ///< Bounds array for chiral_settings
};
  
/// \brief An object for tracking the states and The object itself stores a series of permutor
///        maps, detailing the atoms that move as a consequence of rotating about some bond or
///        inverting a chiral center.  One such map is kept for each unique topology, and the
///        atom indices of a particular system are then obtained by adding the atom offset from
///        the PhaseSpaceSynthesis, Condensate, or AtomGraphSynthesis.
class SynthesisPermutor {
public:

  /// \brief The constructor accepts an array of pointers to existing chemical features objects,
  ///        in which case a pointer to the original object for each system will be retained.  The
  ///        constructor also accepts an array of topologies, in which case the chemical features
  ///        object will be temporarily created for the purpose of extracting critical arrays of
  ///        atoms used in manipulations while the pointer to the chemical features object itself
  ///        will go on to read null.
  ///
  /// \param poly_ag          Topology synthesis (sufficient for generating the object, but will
  ///                         calculate and discard chemical features for all individual
  ///                         topologies)
  /// \param poly_ps          Coordinate synthesis (also sufficient for generating the object,
  ///                         but will again imply creating and discarding chemical features)
  /// \param chemfe_in        A series of pointers to chemical features objects for each system
  /// \param retain_pointers  Indicate whether a pointer to the provided chemical features object
  ///                         should be retained.  This has a default value of true to make it
  ///                         transparent to future development, with any SynthesisPermutor created
  ///                         based on a topology passing a value of false.
  /// \param timer            Wall time tracking object, passed down to ChemicalFeatures
  ///                         computations
  /// \{
  SynthesisPermutor();

  SynthesisPermutor(const AtomGraphSynthesis *poly_ag, StopWatch *timer = nullptr);

  SynthesisPermutor(const AtomGraphSynthesis &poly_ag, StopWatch *timer = nullptr);

  SynthesisPermutor(const AtomGraphSynthesis *poly_ag,
                    const ConformerControls &confcon, StopWatch *timer = nullptr);

  SynthesisPermutor(const AtomGraphSynthesis &poly_ag,
                    const ConformerControls &confcon, StopWatch *timer = nullptr);
  
  SynthesisPermutor(const PhaseSpaceSynthesis *poly_ps, StopWatch *timer = nullptr);

  SynthesisPermutor(const PhaseSpaceSynthesis &poly_ps, StopWatch *timer = nullptr);

  SynthesisPermutor(const PhaseSpaceSynthesis *poly_ps,
                    const ConformerControls &confcon, StopWatch *timer = nullptr);

  SynthesisPermutor(const PhaseSpaceSynthesis &poly_ps,
                    const ConformerControls &confcon, StopWatch *timer = nullptr);
  
  SynthesisPermutor(const std::vector<ChemicalFeatures*> &chemfe_in, bool retain_pointers = true,
                    StopWatch *timer = nullptr);

  SynthesisPermutor(const std::vector<ChemicalFeatures> &chemfe_in, bool retain_pointers = true,
                    StopWatch *timer = nullptr);

  SynthesisPermutor(const std::vector<ChemicalFeatures*> &chemfe_in,
                    const PhaseSpaceSynthesis *poly_ps, const ConformerControls &confcon,
                    const bool retain_pointers = true, StopWatch *timer = nullptr);

  SynthesisPermutor(const std::vector<ChemicalFeatures*> &chemfe_in,
                    const PhaseSpaceSynthesis &poly_ps, const ConformerControls &confcon,
                    const bool retain_pointers = true, StopWatch *timer = nullptr);

  SynthesisPermutor(const std::vector<ChemicalFeatures> &chemfe_in,
                    const PhaseSpaceSynthesis &poly_ps, const ConformerControls &confcon,
                    const bool retain_pointers = true, StopWatch *timer = nullptr);
  /// \}

  /// \brief The presence of POINTER-kind Hybrid objects and multiple ARRAY-kind members means
  ///        that the copy and move constructors as well as assignment operators must be written
  ///        out explicitly.
  ///
  /// \param original  The object to copy or move
  /// \param other     Another object sitting on the right hand side of the assignment operator
  /// \{
  SynthesisPermutor(const SynthesisPermutor &original);
  SynthesisPermutor(SynthesisPermutor &&original);
  SynthesisPermutor& operator=(const SynthesisPermutor &other);
  SynthesisPermutor& operator=(SynthesisPermutor &&other);
  /// \}
  
  /// \brief Get the total number of systems tracked in the synthesis.
  int getSystemCount() const;

  /// \brief Get the number of unique permutor sets (the number of unique topologies served by the
  ///        object).
  int getPermutorSetCount() const;

  /// \brief Match a topology to the permutor map index, or a system in a synthesis to its
  ///        permutor map index.
  ///
  /// Overloaded:
  ///   - Provide the topology by pointer or by const reference
  ///   - Provide a coordinate synthesis by pointer or const reference with a system index
  ///     indicating a topology to search for
  ///
  /// \param query_ag      The topology to search for
  /// \param poly_ps       Coordinate synthesis containing many systems with topology pointers
  /// \param system_index  Index of system within the coordinate synthesis
  /// \{
  int getPermutorMapIndex(const AtomGraph *query_ag) const;
  int getPermutorMapIndex(const AtomGraph &query_ag) const;
  int getPermutorMapIndex(const PhaseSpaceSynthesis *poly_ps, int system_index) const;
  int getPermutorMapIndex(const PhaseSpaceSynthesis &poly_ps, int system_index) const;
  /// \}
  
  /// \brief Get the default number of rotatable bond angle samples.
  int getRotatableBondSampleCount() const;
  
  /// \brief Get the default number of cis-trans isomeric bond angle samples.
  int getCisTransBondSampleCount() const;

  /// \brief Get the number of relevant rotatable bonds.
  ///
  /// Overloaded:
  ///   - Get the number of rotatable bonds for a specific map
  ///   - Get the number of rotatable bonds across all maps
  ///
  /// \param permutor_map_index  Index of the map / topology of interest
  /// \{
  int getRotatableBondCount(const int permutor_map_index) const;
  int getRotatableBondCount() const;
  /// \}

  /// \brief Get the number of relevant cis-trans isomeric bonds.
  ///
  /// Overloaded:
  ///   - Get the number of cis-trans isomeric bonds for a specific map
  ///   - Get the number of cis-trans isomeric bonds across all maps
  ///
  /// \param permutor_map_index  Index of the map / topology of interest
  /// \{
  int getCisTransBondCount(const int permutor_map_index) const;
  int getCisTransBondCount() const;
  /// \}
  
  /// \brief Get the number of relevant chiral centers.
  ///
  /// Overloaded:
  ///   - Get the number of chiral centers for a specific map
  ///   - Get the number of chiral centers across all maps
  ///
  /// \param permutor_map_index  Index of the map / topology of interest
  /// \{
  int getChiralCenterCount(const int permutor_map_index) const;
  int getChiralCenterCount() const;
  /// \}
  
  /// \brief Get the numbers of options for each mutable element in a specific system.
  ///
  /// Overloaded:
  ///   - Search by system index
  ///   - Search by topology
  ///
  /// \param system_index  Index of the system of interest
  /// \param query_ag      Topology of interest (must match one of the topology pointers guiding a
  ///                      map, or the function will return an error)
  /// \{
  const std::vector<int>& getElementSampleCounts(int system_index) const;
  const std::vector<int>& getElementSampleCounts(const AtomGraph *query_ag) const;
  const std::vector<int>& getElementSampleCounts(const AtomGraph &query_ag) const;
  /// \}

  /// \brief Get the total number of variable elements in one of the systems: rotatable bonds,
  ///        cis-trans bonds, and chiral centers (even centers that cannot invert are counted).
  ///
  /// Overloaded:
  ///   - Search by system index
  ///   - Search by topology
  ///
  /// \param system_index  Index of the system of interest
  int getVariableElementCount(int system_index) const;

  /// \brief Get the number of settings to be sampled in one element, of a particular kind, from
  ///        the concatenated list spanning all systems.
  ///
  /// Overloaded:
  ///   - Provide the type and index of the element separately
  ///   - Provide a tuple of the two pieces of information
  ///
  /// \brief element_kind  The type of element, e.g. a rotatable bond
  /// \brief index         Index of the variable element from within the concatenated list
  /// \{
  int getElementSampleCount(ConformationEdit element_kind, int index) const;
  int getElementSampleCount(CoupledEdit ce) const;
  /// \}

  /// \brief Get the number of rotatable groups in one of the systems of the attached synthesis.
  ///
  /// \param system_index  Index of the system of interest
  int getSystemRotatableBondCount(int system_index) const;

  /// \brief Get the number of cis-trans isomeric groups in one of the systems of the attached
  ///        synthesis.
  ///
  /// \param system_index  Index of the system of interest
  int getSystemCisTransBondCount(int system_index) const;

  /// \brief Get the number of chiral centers in one of the systems of the attached synthesis.
  ///
  /// \param system_index  Index of the system of interest
  int getSystemChiralCenterCount(int system_index) const;

  /// \brief Get one of the rotatable groups.
  ///
  /// \param system_index  Index of the system of interest
  /// \param group_index   The rotatable group of interest
  const IsomerPlan& getRotatableGroup(int system_index, int group_index) const;

  /// \brief Get one of the cis-trans isomeric groups.
  ///
  /// \param system_index  Index of the system of interest within the attached synthesis
  /// \param group_index   The cis-trans isomeric group of interest
  const IsomerPlan& getCisTransGroup(int system_index, int group_index) const;

  /// \brief Get one of the invertible chiral centers.
  ///
  /// \param system_index  Index of the system of interest within the attached synthesis
  /// \param group_index   The chiral center of interest
  const IsomerPlan& getInvertibleGroup(int system_index, int group_index) const;

  /// \brief Get a const reference to the state tracker for the molecule.  Also useable to get a
  ///        copy of the state tracker which is then mutable.
  ///
  /// Overloaded:
  ///   - Get the state tracker based on an index in the permutor
  ///   - Get the state tracker based on a system's topology (this will raise an exception if the
  ///     search fails)
  ///
  /// \param system_index  Index of the system of interest
  /// \param ag            System topology
  /// \{
  const TickCounter<double>& getStateTracker(int system_index) const;
  const TickCounter<double>& getStateTracker(const AtomGraph *ag) const;
  const TickCounter<double>& getStateTracker(const AtomGraph &ag) const;
  /// \}

  /// \brief Get the number of replicas needed to fulfill a given level of sampling in one of the
  ///        permutor maps.
  ///
  /// Overloaded:
  ///   - Identify the map of interest by system index from the attached synthesis
  ///   - Identify the map of interest by topology pointer
  ///   - Obtain a list of replicas needed for all permutor maps
  ///
  /// \param query_index   Index of the system of interest within the attached synthesis
  /// \param query_ag      Pointer to the topology of interest
  /// \param effort        The requested level of sampling effort
  /// \{
  llint getReplicaCount(int query_index, SamplingIntensity effort) const;
  llint getReplicaCount(const AtomGraph *query_ag, SamplingIntensity effort) const;
  llint getReplicaCount(const AtomGraph &query_ag, SamplingIntensity effort) const;
  std::vector<llint> getReplicaCount(SamplingIntensity effort) const;
  /// \}

  /// \brief Get the pointer to the chemical features object for one of the permutor maps, which
  ///        can in turn reveal the original topology governing the map.  This function will apply
  ///        a check to ensure that the ChemicalFeatures pointers are valid.
  ///
  /// Overloaded:
  ///   - Get the pointer to one map's underlying features
  ///   - Get a const reference to the array of all such pointers
  ///
  /// \param map_index  Index of the map of interest (this will also be checked for validity) 
  /// \{
  const ChemicalFeatures* getChemicalFeaturesPointer(int map_index) const;
  const std::vector<ChemicalFeatures*> getChemicalFeaturesPointer() const;
  /// \}
  
  /// \brief Get the a pointer to the coordinate synthesis for the collection of permutor maps.
  const PhaseSpaceSynthesis* getSynthesisPointer() const;
  
  /// \brief Get the double-precision read-only abstract for the object.
  ///
  /// \param tier  Get data pointers relevant to the CPU, or to the GPU.
  const SyPermutorKit<double> dpData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

  /// \brief Get the single-precision read-only abstract for the object.
  ///
  /// \param tier  Get data pointers relevant to the CPU, or to the GPU.
  const SyPermutorKit<float> spData(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

#ifdef STORMM_USE_HPC
  /// \brief Upload data to the GPU device.
  void upload();

  /// \brief Download data from the GPU device.
  void download();
#endif

  /// \brief Load specific settings into the available states of rotatable bonds for a particular
  ///        system.  If the vectors differ from the originally stated rotatable bond sample
  ///        counts, the specific system will be updated.
  ///
  /// \param system_index  Index of the system of interest
  /// \param settings      A vector of vectors containing the acceptable settings as real numbers,
  ///                      in units of radians
  void defineRotatableBondSettings(int map_index,
                                   const std::vector<std::vector<double>> &settings);

  /// \brief Load specific settings into the available states of cis-trans bonds for a particular
  ///        system.  If the vectors are not of length 2 (the assumption for the number of choices
  ///        in a cis-trans rotatable bond), the number of states will be updated to accommodate
  ///        the provided choices.
  ///
  /// \param system_index  Index of the system of interest
  /// \param settings      A vector of vectors containing the acceptable settings as real numbers,
  ///                      in units of radians
  void defineCisTransBondSettings(int system_index,
                                  const std::vector<std::vector<double>> &settings);

  /// \brief Load specific settings into the available states of an invertible chiral center for a
  ///        particular system.
  ///
  /// \param system_index  Index of the system of interest
  /// \param settings      A vector of vectors containing the acceptable settings as 0 or 1
  void defineChiralCenterSettings(int system_index,
                                  const std::vector<std::vector<int>> &settings);

  /// \brief Apply a synthesis to the permutor maps, which will be used as the starting point for
  ///        defining new conformations of each molecule.  If cis-trans or chirality sampling are
  ///        shut off, the synthesis will be used to define the settings available for sampling.
  ///
  /// Overloaded:
  ///   - Provide the synthesis by const pointer
  ///   - Provide the synthesis by const reference
  ///
  /// \param poly_ps_in  The input synthesis of coordinates
  /// \param adj         The manner in which to adjust settings for the permutation of molecular
  ///                    geometric variables (e.g. rotatable bonds, chiral centers) in response to
  ///                    the initial state of the synthesis
  /// \{
  void applySynthesis(const PhaseSpaceSynthesis *poly_ps_in,
                      VariableTorsionAdjustment adj = VariableTorsionAdjustment::DO_NOT_CHANGE);

  void applySynthesis(const PhaseSpaceSynthesis &poly_ps_in,
                      VariableTorsionAdjustment adj = VariableTorsionAdjustment::DO_NOT_CHANGE);
  /// \}
  
  /// \brief Set the ranges of variables (rotatable bonds, cis-trans bonds, chiral centers) for
  ///        each permutor map.
  ///
  /// \param general_rot_settings_in  General-purpose settings for any rotatable bond
  /// \param general_ctx_settings_in  General-purpose settings for any cis-trans bond
  /// \param sample_cis_trans_states  Directive to sample cis- and trans-isomeric forms of double
  ///                                 bonds (likely passed down from a &conformer namelist, which
  ///                                 carries the same default value of FALSE)
  /// \param sample_chiral_states     Directive to sample D- and L-chiral states of invertible
  ///                                 atom centers (likely passed down from a &conformer namelist,
  ///                                 which carries the same default value of FALSE)
  void setVariableRanges(const std::vector<double> &general_rot_settings_in,
                         const std::vector<double> &general_ctx_settings_in,
                         bool sample_cis_trans_states = false, bool sample_chiral_states = false);

  /// \brief Transfer control data from a &conformer namelist to the object.  This will include
  ///        vectors of torsion angle settings as well as information on how to adjust the settings
  ///        based on a set of supplied structures.
  ///
  /// \param confcon  User input derived from a &conformer namelist
  void impartControlData(const ConformerControls &confcon);

  /// \brief Manipulate one of the systems inside a synthesis, based on the internal permutor maps.
  ///        The type of manipulation will be interpreted so that the appropriate group of atoms
  ///        can be set to the internal coordinates selected from the array of settings.  See
  ///        randomizeSystem() below for a discussion on the limited advantages of batching these
  ///        operations to work on a single, real-valued coordinate set and then converting all of
  ///        the coordinates back to the split fixed-precision format.
  ///
  /// Overloaded:
  ///   - Manipulate a structure in an arbitrary synthesis, such as one created by the object by
  ///     replicating structures in the synthesis it references internally.
  ///   - Manipulate a structure in the synthesis that this object references internally.
  ///
  /// \param psynth         The collection of systems to operate upon
  /// \param system_index   Index of the system to manipulate within synth
  /// \param map_index      Index of the permutor map (this can be inferred if the attached
  ///                       synthesis is the one being altered)
  /// \param ce             A tuple of the manipulation to perform and its index from within the
  ///                       concatenated lists spanning all permutor maps
  /// \param setting_index  Index of the setting to use, numbered 0, 1, ... N for a choice of N
  ///                       possible options
  /// \{
  void permuteSystem(PhaseSpaceSynthesis *psynth, int system_index, int map_index,
                     const CoupledEdit ce, int setting_index,
                     PrecisionModel prec = PrecisionModel::SINGLE) const;

  void permuteSystem(int system_index, const CoupledEdit ce, int setting_index,
                     PrecisionModel prec = PrecisionModel::SINGLE);
  /// \}

  /// \brief Randomize the system's degrees of freedom.  Up to three degrees of freedom can be
  ///        masked against any changes.  This function inherently converts the (possibly split)
  ///        fixed-precision coordinate arrays of the synthesis to real numbers for manipulation,
  ///        then converts the real numbers back.  However, due to the limited extent to which this
  ///        happens, it is not of much performance advantage to perform the full conversion of all
  ///        coordinates to real values, batch the manipulations, and then convert the results
  ///        back (primarily because of the reverse conversion).  Such an optimization would also
  ///        require additional code, and not be replicable on the GPU.
  ///
  /// Overloaded:
  ///   - Manipulate a structure in an arbitrary synthesis, such as one created by the object by
  ///     replicating structures in the synthesis it references internally.
  ///   - Manipulate a structure in the synthesis that this object references internally.
  ///
  /// \param psynth        The collection of systems to operate upon
  /// \param system_index  Index of the system to manipulate within synth
  /// \param map_index     Index of the permutor map to which the system corresponds (if the
  ///                      synthesis to manipulate differs from the one referenced in the object)
  /// \param xrs           Psuedo-random number generator for selecting random settings of each
  ///                      mutable element
  /// \param prec          Precision model under which to perform the structural operations
  /// \param mask_a        First of up to three features that will be left as they are found
  /// \param mask_c        Second of up to three features that will be left as they are found
  /// \param mask_b        Third of up to three features that will be left as they are found
  /// \{
  void randomizeSystem(PhaseSpaceSynthesis *psynth, int system_index, int map_index,
                       Xoshiro256ppGenerator *xrs, PrecisionModel prec = PrecisionModel::SINGLE,
                       const CoupledEdit mask_a = { ConformationEdit::BOND_ROTATION, -1 },
                       const CoupledEdit mask_b = { ConformationEdit::BOND_ROTATION, -1 },
                       const CoupledEdit mask_c = { ConformationEdit::BOND_ROTATION, -1 }) const;
  
  void randomizeSystem(int system_index, Xoshiro256ppGenerator *xrs,
                       PrecisionModel prec = PrecisionModel::SINGLE,
                       const CoupledEdit mask_a = { ConformationEdit::BOND_ROTATION, -1 },
                       const CoupledEdit mask_b = { ConformationEdit::BOND_ROTATION, -1 },
                       const CoupledEdit mask_c = { ConformationEdit::BOND_ROTATION, -1 });
  /// \}

  /// \brief Produce a new synthesis containing permutations of each structure in the current
  ///        synthesis linked to the object.
  ///
  /// Overloaded:
  ///   - Build the synthesis with or without specifying precision metrics for positions,
  ///     velocities, and force accumulators
  ///   - Provide an option to return the vector indicating which systems in the result originated
  ///     in each system of the synthesis attached to this SynthesisPermutor object.
  ///
  /// \param effort                The degree of sampling to be carried out on the structures of
  ///                              the current synthesis.
  /// \param xrs                   Random number generator to use in randomizing degrees of freedom
  ///                              not explicitly sampled
  /// \param system_limit          The maximum number of systems that the resulting synthesis can
  ///                              contain.  If this limit would be exceeded, structures in the
  ///                              current synthesis will have their respective permutation efforts
  ///                              dialed down until the result falls within the system limits.
  /// \param globalpos_scale_bits  Global position coordinate bits after the decimal
  /// \param localpos_scale_bits   Local position coordinate bits after the decimal
  /// \param velocity_scale_bits   Velocity coordinate bits after the decimal
  /// \param force_scale_bits      Force component bits after the decimal
  /// \param prec                  Precision model in whcih to do calculations for structural
  ///                              manipulations
  /// \param clrep                 Object to collect reports on clashes in each structure as it is
  ///                              produced.  A nullptr value will skip clash detection.  If
  ///                              clashes are detected during EXHAUSTIVE permutation evaluations,
  ///                              the clashing structures will not be included in the output.  In
  ///                              other contexts, if clash-free structures cannot be found within
  ///                              the allowed iterations, the structures will be removed from the
  ///                              output.
  /// \param iteration_limit       The maximum number of iterations for randomizing variable groups
  ///                              before the sturcture is declared a loss
  /// \param max_clashes           The maximum number of clashes between atoms before a structure
  ///                              is declared to clash
  /// \{
  PhaseSpaceSynthesis buildSynthesis(SamplingIntensity effort, Xoshiro256ppGenerator *xrs,
                                     std::vector<int> *correspondence, int system_limit,
                                     int globalpos_scale_bits = default_globalpos_scale_bits,
                                     int localpos_scale_bits = default_localpos_scale_bits,
                                     int velocity_scale_bits = default_velocity_scale_bits,
                                     int force_scale_bits = default_force_scale_bits,
                                     PrecisionModel prec = PrecisionModel::SINGLE,
                                     ClashReport *clrep = nullptr,
                                     int iteration_limit = default_conf_max_seeding_attempts,
                                     int max_clashes = default_conf_clash_pairs) const;

  PhaseSpaceSynthesis buildSynthesis(SamplingIntensity effort, Xoshiro256ppGenerator *xrs,
                                     std::vector<int> *correspondence, PrecisionModel prec,
                                     int system_limit, ClashReport *clrep = nullptr,
                                     int iteration_limit = default_conf_max_seeding_attempts,
                                     int max_clashes = default_conf_clash_pairs) const;

  PhaseSpaceSynthesis buildSynthesis(SamplingIntensity effort, Xoshiro256ppGenerator *xrs,
                                     int system_limit,
                                     int globalpos_scale_bits = default_globalpos_scale_bits,
                                     int localpos_scale_bits = default_localpos_scale_bits,
                                     int velocity_scale_bits = default_velocity_scale_bits,
                                     int force_scale_bits = default_force_scale_bits,
                                     PrecisionModel prec = PrecisionModel::SINGLE,
                                     ClashReport *clrep = nullptr,
                                     int iteration_limit = default_conf_max_seeding_attempts,
                                     int max_clashes = default_conf_clash_pairs) const;

  PhaseSpaceSynthesis buildSynthesis(SamplingIntensity effort, Xoshiro256ppGenerator *xrs,
                                     PrecisionModel prec, int system_limit,
                                     ClashReport *clrep = nullptr,
                                     int iteration_limit = default_conf_max_seeding_attempts,
                                     int max_clashes = default_conf_clash_pairs) const;
  /// \}

private:
  int system_count;                              ///< Total number of systems, and the length of
                                                 ///<   the permutor_indices array
  int permutor_map_count;                        ///< The number of permutor maps for different
                                                 ///<   systems in use by this object
  int rotatable_bond_samples;                    ///< General number of times to sample each
                                                 ///<   rotatable bond (default 3)
  int cis_trans_bond_samples;                    ///< General number of times to sample each
                                                 ///<   cis-trans isomeric bond (default 2)
  bool features_objects_accessible;              ///< Indicates whether the ChemicalFeatures
                                                 ///<   objects upon which the permutor maps are
                                                 ///<   based can be accessed by other functions
                                                 ///<   once the object is created.  Some of the
                                                 ///<   constructors create the objects temporarily
                                                 ///<   and the SynthesisPermutors thus created
                                                 ///<   cannot return valid pointers to the
                                                 ///<   features afterwards.
  double rotation_setting_snap;                  ///< Threshold at which rotatable bond settings
                                                 ///<   will be snapped to clustered values taken
                                                 ///<   from the synthesis at hand.  Units of
                                                 ///<   radians.
  double cis_trans_setting_snap;                 ///< Threshold at which cis-trans isomeric bond
                                                 ///<   settings will be snapped to clustered
                                                 ///<   values taken from the synthesis at hand.
                                                 ///<   Units of radians.
  int total_rotatable_groups;                    ///< The total number of rotatable bond groups
                                                 ///<   tracked by the object, spannig all systems
  int total_cis_trans_groups;                    ///< The total number of cis-trans isomeric bond
                                                 ///<   groups tracked by the object, spanning all
                                                 ///<   systems.
  int total_invertible_groups;                   ///< The total number of chiral centers (this
                                                 ///<   technically includes centers whose chiral
                                                 ///<   inversion protocol is DO_NO_INVERT) tracked
                                                 ///    by the object, spanning all centers
  Hybrid<int> permutor_map_indices;              ///< Index of the permutation map for each system
                                                 ///<   in the synthesis
  Hybrid<int2> permutor_elements;                ///< CoupledEdit objects, translated to int2,
                                                 ///<   indicating the appropriate arrays, and
                                                 ///<   indices therein, where the settings for
                                                 ///<   each mutable element of a given permtuor
                                                 ///<   map can be found.
  Hybrid<int> permutor_element_bounds;           ///< Bounds array for permutor_elements.  The
                                                 ///<   total number of variable elements in the
                                                 ///<   ith permutor map is given by
                                                 ///<   permutor_element_bounds[i + 1] -
                                                 ///<   permutor_element_bounds[i]
  Hybrid<int> system_settings;                   ///< Settings for each system in the synthesis
                                                 ///<   referencing the appropriate TickCounter
                                                 ///<   states.  One TickCounter, a CPU-only
                                                 ///<   object, is stored for each unique permutor
                                                 ///<   map, each unique topology in the synthesis.
                                                 ///<   Feed the settings in this array into the
                                                 ///<   TickCounter to see the actual angle values,
                                                 ///<   or reference the permutor_elements array,
                                                 ///<   then find the appropriate value series in
                                                 ///<   rotatable_bond_settings,
                                                 ///<   cis_trans_bond_settings, or
                                                 ///<   chiral_settings to get the values in memory
                                                 ///<   available to the GPU.
  Hybrid<int> system_settings_limits;            ///< Limits array critical for interpreting and
                                                 ///<   incrementing the values in system_settings,
                                                 ///<   replicating the operations of the
                                                 ///<   TickCounter object.
  Hybrid<int> synthesis_data;                    ///< An ARRAY-kind Hybrid targeted by the above
                                                 ///<   Hybrid<int> objects supporting the
                                                 ///<   synthesis.
  Hybrid<int> rotatable_group_atoms;             ///< Atoms present in each rotatable group of
                                                 ///<   each molecular system, starting with the
                                                 ///<   root and pivot atoms defining the rotatable
                                                 ///<   bond / axis of rotation and continuing to
                                                 ///<   list all that move as a consequence of
                                                 ///<   rotating about the bond.  All such groups
                                                 ///<   are listed back-to-back, with no padding
                                                 ///<   between them.
  Hybrid<int> rotatable_group_bounds;            ///< Bounds array for rotatable_group_atoms.  This
                                                 ///<   array knows where each list of topological
                                                 ///<   atom indices in a rotatable group starts
                                                 ///<   and stops, but not which topology that list
                                                 ///<   pertains to.  For that information, one
                                                 ///<   must follow the array below,
                                                 ///<   system_rotatable_group_bounds, stepping
                                                 ///<   over the range of indices for a particular
                                                 ///<   system.
  Hybrid<int> permutor_rotatable_group_bounds;   ///< Bounds array for rotatable_group_atoms.  This
                                                 ///<   array defines the range of rotatable groups
                                                 ///<   for the ith system in its ith and (i + 1)th
                                                 ///<   elements.
  Hybrid<int> cis_trans_group_atoms;             ///< Atoms present in each isomeric group of each
                                                 ///<   molecular system, beginning with the root
                                                 ///<   and pivot atoms of the bond itself.  All
                                                 ///<   such systems are listed back-to-back, with
                                                 ///<   no padding between them.
  Hybrid<int> cis_trans_group_bounds;            ///< Bounds array for cis_trans_group_atoms, laid
                                                 ///<   out in an analogous manner to
                                                 ///<   rotatable_group_bounds
  Hybrid<int> permutor_cis_trans_group_bounds;   ///< Bounds array for cis_trans_group_bounds, laid
                                                 ///<   out in an analogous manner to
                                                 ///<   system_rotatable_group_bounds
  Hybrid<int> invertible_group_atoms;            ///< Atoms present in each chiral inversion group.
                                                 ///<   The list of moving atoms is prefaced by
                                                 ///<   the indices of the anchor atoms for the two
                                                 ///<   smallest chains.
  Hybrid<int> invertible_group_bounds;           ///< Bounds array for invertible_group_atoms
  Hybrid<int> permutor_invertible_group_bounds;  ///< Bounds array for invertible_group_bounds
  Hybrid<int> invertible_atom_centers;           ///< Concatenated lists for all chiral centers in
                                                 ///<   each permutor map.  The bounds array
                                                 ///<   permutor_invertible_group_bounds also
                                                 ///<   serves to define system limits in this
                                                 ///<   array.
  Hybrid<int> invertible_group_protocols;        ///< Array with int representations of the
                                                 ///<   ChiralInversionProtocol values for each
                                                 ///<   chiral center in the map, also indexed by
                                                 ///<   permutor_invertible_group_bounds.
  Hybrid<int4> rotatable_bond_markers;           ///< Atoms that define the angle made by a
                                                 ///<   particular rotatable bond
  Hybrid<int4> cis_trans_bond_markers;           ///< Atoms that define the angle made by a
                                                 ///<   particular cis-trans isomeric bond
  Hybrid<int4> chiral_markers;                   ///< Atoms that define chirality of a particular
                                                 ///<   center, with the priority order going
                                                 ///<   y > z > w > x (see the chiral_arm_atoms
                                                 ///<   member variable in the ChemicalFeatures
                                                 ///<   object)
  Hybrid<double> rotatable_bond_settings;        ///< Angle settings for each rotatable bond, as
                                                 ///<   defined by atoms in rotatable_bond_markers
  Hybrid<double> cis_trans_bond_settings;        ///< Angle settings for each cis-trans bond, as
                                                 ///<   defined by atoms in cis_trans_bond_markers
  Hybrid<float> sp_rotatable_bond_settings;      ///< Angle settings for each rotatable bond, as
                                                 ///<   defined by atoms in rotatable_bond_markers
                                                 ///<   (single-precision)
  Hybrid<float> sp_cis_trans_bond_settings;      ///< Angle settings for each cis-trans bond, as
                                                 ///<   defined by atoms in cis_trans_bond_markers
                                                 ///<   (single-precision)
  Hybrid<int> rotatable_bond_settings_bounds;    ///< Bounds array for rotatable_bond_settings
  Hybrid<int> cis_trans_bond_settings_bounds;    ///< Bounds array for cis_trans_bond_settings
  Hybrid<int> chiral_settings_bounds;            ///< Bounds array for chiral_settings.  This is an
                                                 ///<   ARRAY-kind Hybrid object that does not
                                                 ///<   target group_data.
  Hybrid<int> group_data;                        ///< ARRAY-kind Hybrid object targeted by the
                                                 ///<   preceding POINTER-kind Hybrid<int> objects
  Hybrid<int4> marker_data;                      ///< ARRAY-kind Hybrid object targeted by the
                                                 ///<   preceding POINTER-kind Hybrid<int4> objects
  Hybrid<int> chiral_settings;                   ///< ChiralOrientation enumerations translated to
                                                 ///<   int values for each possible setting of
                                                 ///<   each chiral center.  If chiral sampling is
                                                 ///<   active for a specific topology, each center
                                                 ///<   (save for those with a transformation
                                                 ///<   protocol of DO_NOT_INVERT) will have two
                                                 ///<   options.  If chiral sampling is not active
                                                 ///<   for a specific topology, each center will
                                                 ///<   have one option for every orientation (D-
                                                 ///<   or L-) found in structures of the original
                                                 ///<   synthesis.  This is an ARRAY-kind Hybrid
                                                 ///<   object, like the rotatble and cis-trans
                                                 ///<   isomeric bond settings, that does not
                                                 ///<   target group_data.

  // The following arrays store couple variables needed to understand "light" and "heavy" sampling.
  // The case of "minimal" sampling is trivial--sample all states of each variable individually
  // while other variables get set to random values.  In "light" sampling the goal is to find all
  // combinations of two coupled variables, but some variables will exist only in isolation.  In
  // "heavy" sampling combinations of three variables are applied, but some variables may have no
  // couplings or only one other variable to connect with.  For these cases where pairs or triplets
  // cannot be fulfilled, the tuples representing them will carry a value of -1.  Exhaustive
  // sampling will place all variables in one massive tuple and do all permutations.
  Hybrid<int2> light_sampling_foci;        ///< Pairs of coupled variables, every two entries
                                           ///<   including one pair of coupled variables.  The "y"
                                           ///<   members of tuples in the array contain the
                                           ///<   indices of the relevant groups out of the
                                           ///<   concatenated arrays spanning all systems, while
                                           ///<   the "x" members of each tuple indicate which type
                                           ///<   of variable the "y" member is indexing.
  Hybrid<int2> heavy_sampling_foci;        ///< Triplets of linked variables, every three entries
                                           ///<   including one triplet.  The "x" and "y" members
                                           ///<   of the tuples work as in light_sampling_foci:
                                           ///<   reflecting the CoupledEdit object.
  Hybrid<int> light_sampling_foci_bounds;  ///< Bounds array for light_sampling_foci
  Hybrid<int> heavy_sampling_foci_bounds;  ///< Bounds array for heavy_sampling_foci
  Hybrid<int2> variables_data;             ///< ARRAY-kind Hybrid object targeted by three other
                                           ///<   member variables: permutor_elements,
                                           ///<   light_sampling_foci, and heavy_sampling_foci
  
  // The following CPU-based data arrays help to stage the Hybrid data.
  std::vector<double> general_rotation_settings;  ///< Settings for rotatable bonds, if they are
                                                  ///<   not given for a specific system
  std::vector<double> general_cis_trans_settings; ///< Settings for cis-trans rotatable bonds, if
                                                  ///<   they are not given for a specific system
  std::vector<IsomerPlan> rotatable_groups;       ///< A list of rotatable groups, copied from
                                                  ///<   each system's chemical features.  The
                                                  ///<   permutor_rotatable_group_bounds array can
                                                  ///<   also be thought of as bounding this.
  std::vector<IsomerPlan> cis_trans_groups;       ///< A list of cis-trans isomeric groups, copied
                                                  ///<   from each system's chemical features.  The
                                                  ///<   permutor_cis_trans_group_bounds array can
                                                  ///<   also be thought of as bounding this.
  std::vector<IsomerPlan> invertible_groups;      ///< A list of invertible groups for chiral
                                                  ///<   centers in each system, copied from their
                                                  ///<   respective chemical features objects.  The
                                                  ///<   permutor_invertible_group_bounds array can
                                                  ///<   also be thought of as bounding this.

  // These CPU-based arrays put dimensions on the permutations required to explore each system.
  std::vector<TickCounter<double>> state_trackers;  ///< Exemplary counter wheels for looping over
                                                    ///<   each unique system's mutable variables
  std::vector<int> minimal_sampling_replicas;       ///< Numer of permutations in each system
                                                    ///<   required for MINIMAL sampling
  std::vector<int> light_sampling_replicas;         ///< Number of permutations in each system
                                                    ///<   required for LIGHT (coupled pairs of
                                                    ///<   mutable variables) sampling
  std::vector<int> heavy_sampling_replicas;         ///< Number of permutations in each system
                                                    ///<   required for HEAVY (grouped triads of
                                                    ///<   mutable variables) sampling
  std::vector<llint> exhaustive_sampling_replicas;  ///< Number of permutations in each system
                                                    ///<   required for EXHAUSTIVE (all possible
                                                    ///<   combinations of mutable variables)
                                                    ///<   sampling

  /// A list of pointers to the unique topologies guiding each permutor map
  std::vector<AtomGraph*> topologies;
  
  /// Chemical features for each unique topology
  std::vector<ChemicalFeatures*> features;       

  /// The coordinate synthesis referenced by this permutor.  A coordinate synthesis must be present
  /// in order to properly account for sampling of cis-trans and chiral variables in each system:
  /// the permutor maps will accept settings based on what they see in the input structures.
  PhaseSpaceSynthesis *poly_ps_ptr;
  
  /// \brief Repair (rebase) POINTER-kind HYBRID objects to the object's present ARRAY-kind Hybrid
  ///        targets.
  void rebasePointers();
  
  /// \brief Create a temporary vector of ChemicalFeatures objects based on the unique topologies
  ///        at hand, for use in constructing the object.
  ///
  /// \param timer  Wall time tracker
  std::vector<ChemicalFeatures> temporaryFeatures(StopWatch *timer = nullptr) const;

  /// \brief Create pointers to each member of a list of ChemicalFeatures objects.
  ///
  /// \param feat_in  The list of ChemicalFeatures objects
  const std::vector<ChemicalFeatures*>
  temporaryFeaturesPointers(const std::vector<ChemicalFeatures> &feat_in) const;

  /// \brief Fill the permutor maps based on a defined collection of chemical features.
  ///
  /// \param chemfe_in  The chemical features for each unique system
  /// \param timer      Wall time tracker
  void fillPermutorMaps(const std::vector<ChemicalFeatures*> &chemfe_in, StopWatch *timer);
  
  /// \brief Validate a request for a specific system index.
  ///
  /// \param system_index  The index of interest
  void validateSystemIndex(int system_index) const;

  /// \brief Validate a request for a specific permutor map index.
  ///
  /// \param map_index  The index of interest
  void validateMapIndex(int map_index) const;

  /// \brief Alter the settings of a particular class of variables within one of the permutor maps.
  ///
  /// \param map_index                 Index of the permutor map of interest
  /// \param variable_settings         The settings for variables of the chosen type within the
  ///                                  map of interest.  Altered and returned.
  /// \param variable_settings_bounds  Bounds array for variable_settings
  /// \param permutor_map_bounds       Bounds array for variable_settings_bounds
  /// \param settings                  The new values to apply to variable_settings in the relevant
  ///                                  indices for permutor map map_index
  template <typename T> void alterMapSettings(int map_index, Hybrid<T> *variable_settings,
                                              const Hybrid<int> &variable_settings_bounds,
                                              const Hybrid<int> &permutor_map_bounds,
                                              const std::vector<std::vector<T>> &settings);

  /// \brief Private overload of the permuteSystem() function, called by the public member
  ///        functions if they are invoked.  Parameter descriptions follow from above, with the
  ///        addition of:
  ///
  /// \param cs  Coordinates of a single frame to alter.  The data type of the series will guide
  ///            the numerical precision used in calculations--it is a series, not a
  ///            CoordinateFrame object, to make use of the type formatting options.
  template <typename T>
  void permuteSystem(CoordinateSeries<T> *cs, int map_index, CoupledEdit ce,
                     int setting_index) const;

  /// \brief Private overload of the permuteSystem() function, called by the public member
  ///        functions if they are invoked.  Parameter descriptions follow from above, with the
  ///        addition of:
  ///
  /// \param cs  Coordinates of a single frame to alter.  The data type of the series will guide
  ///            the numerical precision used in calculations--it is a series, not a
  ///            CoordinateFrame object, to make use of the type formatting options.
  template <typename T>
  void randomizeSystem(CoordinateSeries<T> *cs, int map_index, Xoshiro256ppGenerator *xrs,
                       const CoupledEdit mask_a = { ConformationEdit::BOND_ROTATION, -1 },
                       const CoupledEdit mask_b = { ConformationEdit::BOND_ROTATION, -1 },
                       const CoupledEdit mask_c = { ConformationEdit::BOND_ROTATION, -1 }) const;

  /// \brief Determine whether atoms in the modified conformation contain clashes and continue to
  ///        randomize whatever elements are available until a suitable solution can be found.
  ///        The function returns TRUE if the clash check passes or was skipped, FALSE otherwise.
  ///
  /// \param cs               Coordinates of the conformation (one frame only)
  /// \param map_index        Index of the permutor map guiding the system randomizations
  /// \param xrs              Random number generator supplying random moves
  /// \param excl             Exclusion mask for all-to-all interactions
  /// \param summary          Contains clash criteria and tracks the number of clashing atom pairs
  /// \param iteration_limit  The maximum number of randomization iterations to attempt before
  ///                         declaring the conformation invalid
  /// \param max_clashes      The maximum number of atom pair clashes that will be tolerated before
  ///                         the system is declared invalid and must be randomized again
  /// \param mask_a           First feature to protect from alterations
  /// \param mask_b           Second feature to protect from alterations
  /// \param mask_c           Third feature to protect from alterations
  template <typename T>
  bool resolveClashes(CoordinateSeries<T> *cs, int map_index, Xoshiro256ppGenerator *xrs,
                      const StaticExclusionMask &excl, ClashReport *summary, int iteration_limit,
                      int max_clashes,
                      const CoupledEdit mask_a = { ConformationEdit::BOND_ROTATION, -1 },
                      const CoupledEdit mask_b = { ConformationEdit::BOND_ROTATION, -1 },
                      const CoupledEdit mask_c = { ConformationEdit::BOND_ROTATION, -1 }) const;
};

} // namespace synthesis
} // namespace stormm

#include "synthesis_permutor.tpp"

#endif

