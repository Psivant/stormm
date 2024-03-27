// -*-c++-*-
#ifndef STORMM_LOCAL_EXCLUSIONMASK_H
#define STORMM_LOCAL_EXCLUSIONMASK_H

#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "DataTypes/common_types.h"
#include "DataTypes/stormm_vector_types.h"
#include "Math/formulas.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "energy_enumerators.h"

namespace stormm {
namespace energy {

using card::Hybrid;
using card::HybridTargetLevel;
using stmath::ipow;
using synthesis::AtomGraphSynthesis;
using topology::AtomGraph;
using topology::ChemicalDetailsKit;
using topology::NonbondedKit;

/// \brief Modes for the profile, expressed a unsigned long long integers to eliminate the need
///        for an extra bit shift operation.
/// \{
static const ullint lmask_mode_a       = 0x0LLU;
static const ullint lmask_mode_b       = 0x2000000000000000LLU;
static const ullint lmask_mode_c       = 0x4000000000000000LLU;
static const ullint lmask_mode_d       = 0x6000000000000000LLU;
static const ullint lmask_mode_e       = 0x8000000000000000LLU;
static const ullint lmask_mode_f       = 0xa000000000000000LLU;
/// \}

/// Various masks useful for dissecting the exclusion profile
/// \{
static const ullint lmask_mode_bitmask = 0xe000000000000000LLU;
static const ullint lmask_a_excl       = 0x1fffffffffffffffLLU;
static const ullint lmask_b_excl       = 0x000000007fffffffLLU;
static const ullint lmask_b_lower_mask = 0x000001ff80000000LLU;
static const ullint lmask_b_upper_mask = 0x0007fe0000000000LLU;
static const ullint lmask_b_lower_shft = 0x00f8000000000000LLU;
static const ullint lmask_b_upper_shft = 0x1f00000000000000LLU;
static const ullint lmask_c_excl       = 0x000000007fffffffLLU;
static const ullint lmask_c_alt_mask   = 0x00001fff80000000LLU;
static const ullint lmask_c_shft       = 0x1fffe00000000000LLU;
static const ullint lmask_d_excl       = 0x000000007fffffffLLU;
static const ullint lmask_d_array_idx  = 0x003fffff80000000LLU;
static const ullint lmask_d_array_cnt  = 0x1fc0000000000000LLU;
static const ullint lmask_e_array_idx  = 0x000000ffffffffffLLU;
static const ullint lmask_e_array_cnt  = 0x1fffff0000000000LLU;
static const ullint lmask_f_idx        = 0x000000000000001fLLU;
/// \}

/// \brief Critical measurements for various exclusion profiles
/// \{
static const int lmask_long_local_span   = 30;
static const int lmask_short_local_span  = 15;
static const int lmask_long_pivot_point  = lmask_long_local_span + 1;
static const int lmask_short_pivot_point = lmask_short_local_span + 1;
static const int lmask_long_extra_span   = 14;
static const int lmask_short_extra_span  = 10;
/// \}
  
/// \brief Bit lengths of various components of different profiles.
/// \{
static const int lmask_mode_bits         =  3;
static const int lmask_b_lower_mask_bits = lmask_short_extra_span;
static const int lmask_b_upper_mask_bits = lmask_short_extra_span;
static const int lmask_b_shft_bits       =  5;
static const int lmask_b_max_reach       = lmask_short_local_span + ipow(2, lmask_b_shft_bits) +
                                           lmask_short_extra_span - 1;
static const int lmask_c_alt_mask_bits   = lmask_long_extra_span;
static const int lmask_c_shft_bits       = 16;
static const int lmask_d_array_idx_bits  = 23;
static const int lmask_d_array_cnt_bits  =  7;
static const int lmask_e_array_idx_bits  = 40;
static const int lmask_e_array_cnt_bits  = 21;
/// \}

/// \brief Bit counts by which to shift profiles to the right in order to put the important bits
///        of certain numbers in alignment for direct reading as integers.
/// \{
static const int lmask_mode_pos          = (2 * lmask_long_local_span) + 1;
static const int lmask_b_lower_mask_pos  = (2 * lmask_short_local_span) + 1;
static const int lmask_b_upper_mask_pos  = lmask_b_lower_mask_pos + lmask_short_extra_span;
static const int lmask_b_lower_shft_pos  = lmask_b_upper_mask_pos + lmask_short_extra_span;
static const int lmask_b_upper_shft_pos  = lmask_b_lower_shft_pos + lmask_b_shft_bits;
static const int lmask_c_alt_mask_pos    = (2 * lmask_short_local_span) + 1;
static const int lmask_c_shft_pos        = lmask_c_alt_mask_pos + lmask_c_alt_mask_bits;
static const int lmask_d_array_idx_pos   = (2 * lmask_short_local_span) + 1;
static const int lmask_d_array_cnt_pos   = lmask_d_array_idx_pos + lmask_d_array_idx_bits;
static const int lmask_e_array_cnt_pos   = 40;
/// \}

/// \brief A lean, read-only abstract for the LocalExclusionMask class.
struct LocalExclusionMaskReader {

  /// \brief The constructor for this and other abstracts takes input arguments for each member
  ///        variable.  It is lean, and best used in concert with the topology synthesis or cell
  ///        grid with which it is associated, so that the system atom limits are most accessible.
  LocalExclusionMaskReader(const int* prof_idx_in, const ullint* profiles_in,
                           const uint2* aux_masks_in);

  const int* prof_idx;     ///< Profile indices of each atom
  const ullint* profiles;  ///< Atom profiles, containing the mode and various other bit packed
                           ///<   elements to support multiple exclusion mask patterns
  const uint2* aux_masks;  ///< Auxiliary exclusion masks, with the offset from the index (central)
                           ///<   atom in the "x" member and the mask of 32 bits in the "y" member
};
  
/// \brief The local exclusion mask is an generalized form of the ForwardExclusionMask, listing
///        all exclusions in various ranges of an atom followed by additional masks.  Like the
///        ForwardExclusionMask, all atoms are assigned an index, but it points to a primary
///        array of profiles, not quite like the primary masks of the ForwardExclusionList.  There
///        are a finite number of profiles, each containing a three-bit code to indicate how it is
///        to be used.
///
/// When accessed, the profile is a 64-bit unsigned long long int, with several modes of operation
/// indicated by its high three bits.  For the purposes of the following discussion, we will take
/// atom (i) of the topology (or the atom with absolute index (i) out of a synthesis of topologies)
/// to access profile (j) and possibly secondary masks (k1) through (kN) when determining the
/// existence of an excluded interaction with some other atom (m).
///
/// The first mode is the simplest, and covers nearly all situations that will be found in a
/// typical MD simulation.  The low 61 bits describe a mask such that the 31st bit, always 0, is
/// the exclusion of the atom itself.  Bits 1-30 indicate exclusions for atoms 30, 29, ..., 2, 1
/// topological position prior to the atom of interest (i).  Bits 32-61 indicate exclusions for
/// atoms 1, 2, ..., 29, 30 positions after atom (i).  After retrieving the profile, the existence
/// of an exclusion between (i) and some other atom (m) is computed by:
///
///   int del_im = (m) - (i);
///   const bool excluded = (abs(del_im) < 31 && (((profile of i) >> (30 + del_im)) & 0x1));
///
/// Mode A:       000 ------------------------------ 1 ------------------------------
///
/// The mode of an atom with no exclusions, whether a monatomic ion or a water oxygen in the
/// context of a van-der Waals pairlist, is also mode A: every atom excludes itself.
///
/// In the second mode, when the 62nd bit is set to one, the 31 bits are used in much the same way
/// as the low 61 bits of the profile in Mode A.  However, the next higher 12 bits indicate a patch
/// of exclusions offset by -((6 * NNN) + 27) from atom (i) and the next higher 12 bits a similar
/// patch of exclusions offset +((6 * PPP) + 15) atom (i).  The code for evaluating exclusions is
/// somewhat more involved but these types of masks will be rare.
///
/// Mode B:       001 PPPPP NNNNN ---------- ---------- --------------- 1 ---------------
///
/// In the diagram above, NNNNN and PPPPP are unsigned integers extracted from bits 55-57 and
/// 58-60, respectively, and one bit is left unused.  This more complex footprint, which can
/// accommodate atoms with exclusions up to 56 topological indices away from the home atom, can
/// be easily evaluated in most modes, if the difference between (i) and (m) is greater than 56.
/// 
/// In the third mode of operation, the 63rd bit of the profile is set to 1.  The low 31 bits of
/// the profile continue to function as before, but the next 14 represent a mask that is shifted
/// relative to the original atom by -32768 to +32767 indices, as indicated by the high group of
/// 16 bits just below the three mode bits (read the S bits below as an unsigned integer in the
/// range [0, 65536) and subtract 32768).  This is effective for the largest known proteins with
/// disulfide bridges connecting one end of the chain to another.
///
/// Mode C:       010 SSSSSSSSSSSSSSSS -------------- --------------- 1 ---------------
///
/// In the fourth mode of operation, the 64th bit of the profile set to 1.  In this mode, the low
/// 31 bits continue to function as they did in modes A and B, but the next 23 bits indicate an
/// lower bound index into a secondary array, much as in the ForwardExclusionMask mode, where
/// additional groups of exclusions are to be found.  The highest 7 bits just behind the mode bits
/// are taken as an unsigned integer indicating an upper bound, the number of consecutive indices
/// to access from within the secondary array.  The secondary array is a collection of unsigned
/// int pairs, the "x" member being another local exclusion mask and the "y" member being a shift
/// relative to the original atom.  This mode would only be activated when there are extreme
/// numbers of far-flung and non-contiguous exclusions for a single atom.
///
/// Mode D:       011 KKKKKKK IIIIIIIIIIIIIIIIIIIIIIII --------------- 1 ---------------
///
/// In the sixth mode of operation, the 62nd and 64th bits of the profile are set to 1.  The
/// profile itself no longer provides local exclusion masks.  Instead, the low 40 bits operate as
/// the index into the array of secondary masks while the higher 21 bits (those just behind the
/// mode bits) provide the upper bound on the number of masks from the secondary array to access.
/// Any system operating in this mode would have to contain exclusions of extreme complexity and
/// the amount of available computer memory would become an issue before the format is broken.  It
/// could also be inefficient if taken to the extreme, where every atom has a large and far-flung
/// but unique exclusions profile.  Still, many thousands of atoms could be accommodated in this
/// manner.
///
/// Mode E:       100 KKKKKKKKKKKKKKKKKKKKKK IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
///
/// The final mode is applied when only a few atoms' exclusions (up to 32) are so astoundingly
/// complex that it is easier to simply enumerate them for all atoms of the system (or group of
/// systems) and perform a non-caching read of the pre-tabulated lists.  This is not expected to
/// occur in any ordinary MD simulation.  The five lowest bits of the profile indicate the
/// index of the pre-tabulated list to take in this mode.
///
/// Mode F:       101 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX IIIII
///
/// One final mode may be useful, but the list is quite thorough.
class LocalExclusionMask {
public:

  /// \brief Constructor requires a topology, and creates a blank object if nullptr is supplied.
  /// \{
  LocalExclusionMask(const AtomGraph *ag_in = nullptr,
                     NonbondedTheme theme_in = NonbondedTheme::ALL);

  LocalExclusionMask(const AtomGraph &ag_in, NonbondedTheme theme_in = NonbondedTheme::ALL);

  LocalExclusionMask(const AtomGraphSynthesis *poly_ag_in,
                     NonbondedTheme theme_in = NonbondedTheme::ALL);

  LocalExclusionMask(const AtomGraphSynthesis &poly_ag_in,
                     NonbondedTheme theme_in = NonbondedTheme::ALL);
  /// \}

  /// \brief The default copy and move constructors as well as the copy and move assignment
  ///        operator will suffice for this object, which has no POINTER-kind Hybrid objects among
  ///        its members.
  /// \{
  LocalExclusionMask(const LocalExclusionMask &original) = default;
  LocalExclusionMask(LocalExclusionMask &&original) = default;
  LocalExclusionMask& operator=(const LocalExclusionMask &other) = default;
  LocalExclusionMask& operator=(LocalExclusionMask &&other) = default;
  /// \}

  /// \brief Get the number of atoms in the system.
  int getAtomCount() const;

  /// \brief Get the number of unique profiles.
  ///
  /// Overloaded:
  /// - Get the total number of profiles 
  /// - Get the number of profiles using a particular mode
  ///
  /// \param mode_index  The index number of the mode of operation for the profile
  /// \{
  int getProfileCount() const;
  int getProfileCount(int mode_index) const;
  /// \}

  /// \brief Get the profile mode used by a particular atom as a positive integer.
  ///
  /// \param atom_index    Index of the atom of interest
  /// \param system_index  Index of the system of interest (applicable for topology syntheses)
  int getMode(int atom_index, int system_index = 0) const;

  /// \brief Get the number of secondary masks.
  size_t getSecondaryMaskCount() const;

  /// \brief Get the exclusion profile for an atom.
  ///
  /// \param atom_index    Topological index of the atom of interest
  /// \param system_index  Index of the system of interest (applicable for topology syntheses)
  ullint getProfile(int atom_index, int system_index = 0) const;

  /// \brief Get a POINTER-kind Hybrid object for viewing all secondary masks associated with an
  ///        atom.  Each presents the exclusions for a cluster of atoms following some reference
  ///        atom.  Descriptions of parameters follow from getProfile() above.
  const Hybrid<uint2> getSecondaryMaskView(int atom_index, int system_index = 0);

  /// \brief Get a copy of all secondary masks associated with an atom.  Descriptions of parameters
  ///        follow from getProfile() above.
  /// \{
  std::vector<uint2> getSecondaryMasks(int atom_index, int system_index = 0) const;
  /// \}

  /// \brief Get a pointer to the topology associated with this object.
  const AtomGraph* getTopologyPointer() const;

  /// \brief Get a pointer to the topology synthesis associated with this object.
  const AtomGraphSynthesis* getTopologySynthesisPointer() const;

  /// \brief Return an indication of whether two atoms in the same molecular system share an
  ///        excluded interaction.
  ///
  /// Overloaded:
  ///   - Provide two atom indices to test the exclusion of two atoms in a single topology
  ///   - Provide the system index in order to test the exclusion of two atoms in one system out
  ///     of a synthesis
  ///
  /// \param atom_i        Index of the first atom within the topology of interest
  /// \param atom_j        Index of the second atom within the topology of interest
  /// \param system_index  The system of interest (valid only for objects built on a topology
  ///                      synthesis)
  /// \{
  bool testExclusion(int atom_i, int atom_j) const;
  bool testExclusion(int atom_i, int atom_j, int system_index) const;
  /// \}

  /// \brief Get the abstract for the object.
  ///
  /// \param tier  Indicate whether to target pointers to memory on the CPU host or GPU device
  LocalExclusionMaskReader data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;

#ifdef STORMM_USE_HPC
  /// \brief Upload all data to the GPU device.
  void upload();
  
  /// \brief Download all data from the GPU device.
  void download();
#endif
  
  /// \brief Set (or reset) the topology for the object.  All exclusion masks and profiles will
  ///        be rebuilt if the topology changes.
  ///
  /// Overloaded:
  ///   - Provide a single topology by pointer or const reference
  ///   - Provide a synthesis of topologies by pointer or const reference
  ///
  /// \param ag_in  The new topology to set
  /// \{
  void setTopology(const AtomGraph *ag_in);
  void setTopology(const AtomGraph &ag_in);
  void setTopology(const AtomGraphSynthesis *poly_ag_in);
  void setTopology(const AtomGraphSynthesis &poly_ag_in);
  /// \}
  
private:

  /// The number of atoms with codes managed by this object.  This is the number of atoms in the
  /// topology, if this object is created with a single AtomGraph.  If it is created with a
  /// synthesis of topologies, it is the total padded number of atoms in the synthesis.
  int total_atom_count;

  /// The non-bonded theme to which the object is tailored.  The LocalExclusionMask can work with
  /// a CellGrid object, for example, which groups atoms by kinds of non-bonded properties.
  NonbondedTheme theme;
  
  /// The numbers of profiles operating in each of the modes described above
  /// \{
  int mode_a_profile_count;
  int mode_b_profile_count;
  int mode_c_profile_count;
  int mode_d_profile_count;
  int mode_e_profile_count;
  int mode_f_profile_count;
  /// \}
  
  /// Profile indices for all atoms
  Hybrid<int> atom_profile_codes;
  
  /// Primary exclusion profiles for all atoms.  There are four modes under which these profiles
  /// can be applied, as indicated by the high two bits of each 64-bit unsigned integer.  See the
  /// class documentation above for a thorough description.
  Hybrid<ullint> atom_profiles;

  /// Secondary exclusion masks, accessed by profiles set to use one of the latter two modes.
  /// Situations in which these would be necessary are not common, but this feature provides
  /// coverage for any conceivable case of molecular systems with extreme exclusion maps.
  Hybrid<uint2> secondary_masks;

  /// Pointer to the original topology
  AtomGraph *ag_pointer;

  /// Pointer to the original topology synthesis (this and ag_pointer are mutually exclusive--if
  /// one has a valid pointer the other will be set to nullptr)
  AtomGraphSynthesis *poly_ag_pointer;

  /// \brief Set the atom count and various profile counts all to zero, and erase any atom profiles
  ///        and secondary masks.  This will be called in response to both topology and topology
  ///        synthesis pointers being set to nullptr.
  void nullify();

  /// \brief Validate an atom index, whether against the entire atom count or the count of atoms in
  ///        a particular topology.
  ///
  /// Overloaded:
  ///   - Provide a system index (valid for topology syntheses)
  ///   - Provide just the atom index (valid for single topologies)
  ///
  /// \param atom_index    Index of the atom of interest
  /// \param system_index  Index of the system of interest (applicable for topology syntheses)
  /// \param caller        Name of the calling function (for backtracing purposes)
  /// \{
  void validateAtomIndex(int atom_index, int system_index = 0, const char* caller = nullptr) const;
  void validateAtomIndex(int atom_index, const char* caller) const;
  /// \}

  /// \brief Count the number of profiles of each type in the condensed list, based on the mode
  ///        bits.
  void countProfileKinds();
  
  /// \brief Extend the atom masks to cover a new topology.  This can be the one topology which
  ///        the object is designed to serve, or one of the series of topoologies within the
  ///        associated synthesis.  This function extends the lists of atom profiles and secondary
  ///        masks as needed, then returns a list of profile indices into the newly expanded
  ///        arrays.  The list of atom profile indices can be commuted directly into the
  ///        atom_profile_codes array if a single topology is in use, or applied iteratively to
  ///        all instances of the same topology in a synthesis.
  ///
  /// Overloaded:
  ///   - Provide a pointer to the system's topology
  ///   - Provide the nonbonded parameters and chemical details
  ///
  /// \param ag           The topology containing exclusions to track
  /// \param nbk          The system's non-bonded parameter abstract, containing exclusion lists
  /// \param cdk          The system's chemical details abstract, containing atom and residue names
  /// \{
  std::vector<int> extendMasks(const NonbondedKit<double> &nbk, const ChemicalDetailsKit &cdk);
  std::vector<int> extendMasks(const AtomGraph *ag);
  /// \}
};

/// \brief Create a list of exclusions for use in creating objects of the LocalExclusionMask class.
///        This function is free because there is no need to place it within the class itself, and
///        to facilitate unit testing with mock topologies, which are hard to create with arbitrary
///        connectivities.
///
/// \param nbk         Non-bonded parameters of the system, including 1:1, 1:2, 1:3, and 1:4
///                    exclusions
/// \param atom_index  Topological index of the central atom
/// \param theme       Indicate the kind of non-bonded interactions relevant to the associated
///                    neighbor list and cell grid
/// \param lj_rule     Combining rule determined for the topology 
std::vector<int> compileLocalExclusionList(const NonbondedKit<double> &nbk, int atom_index,
                                           NonbondedTheme theme, VdwCombiningRule lj_rule);

/// \brief Transform one of the mode codes from an atom profile in the LocalExclusionMask class to
///        a string.
std::string lMaskModeToString(const ullint mode);
  
/// \brief Determine whether the exclusion pattern fits mode A, the simplest mode after the case
///        of no exclusions (D).
///
/// \param atom_index      The index (central) atom for which the mask is being made
/// \param excluded_atoms  List of the topologicla indices of all atoms excluded by the index atom
bool lMaskFitsModeA(int atom_index, const std::vector<int> &excluded_atoms);

/// \brief Determine whether the exclusion pattern fits mode B, the most rapidly calculable profile
///        mode after mode A (in that no additional memory access is required).  Descriptions of
///        parameters follow from lMaskFitsModeA(), above.
bool lMaskFitsModeB(int atom_index, const std::vector<int> &excluded_atoms);

/// \brief Determine whether the exclusion pattern fits mode C, the next most rapidly calculable
///        profile requiring a single extra table lookup.  Mode C is typically the way to express
///        exclusions in disulfide bridges or other very large loops in macrocycles.  Descriptions
///        of parameters follow from lMaskFitsModeA(), above.
bool lMaskFitsModeC(int atom_index, const std::vector<int> &excluded_atoms);

/// \brief Determine whether the exclusion pattern fits mode D, an intensely laborious process
///        that will only be undertaken if there is no alternative, and is not known to be needed
///        by the molecules laid out in any conventional MD force field.  Many extra table lookups
///        are required in a compact loop of up to 128 auxiliary masks.  Descriptions of
///        parameters follow from lMaskFitsModeA(), above, in addition to:
///
/// \param tmp_secondary_masks  Temporary array of the secondary masks, used to check size limits
///                             in mode D profiles
bool lMaskFitsModeD(int atom_index, const std::vector<int> &excluded_atoms,
                    const std::vector<uint2> &tmp_secondary_masks);

/// \brief Create a profile in mode A for an object of the LocalExclusionMask class.  This will
///        return a 64-bit code for the atom of interest suitable for including in the array of
///        atom profiles.
///
/// Overloaded:
///   - Return the profile only for modes A, B, or C
///   - Return the array of secondary masks for modes D and E, with the profile returned by
///     modifying an input
///
/// \param atom_index       Topological index of the central atom
/// \param mode             The detected (and requested) mode for the profile
/// \param excluded_atoms   List of exclusions for the central atom, sorted
/// \param t_profile        The profile created by the function, modified and returned if mode E or
///                         F is requested.  It is returned in an incomplete state: the index
///                         portion of the array must be 
/// \{
ullint lMaskCreateProfile(int atom_index, ullint mode, const std::vector<int> &excluded_atoms);
  
std::vector<uint2> lMaskCreateProfile(int atom_index, ullint mode,
                                      const std::vector<int> &excluded_atoms, ullint *t_profile);
/// \}
  
/// \brief Set the profile index of a given atom, after checking for uniqueness of the profile
///        and contributing the profile to the list if necessary.
///
/// Overloaded:
///   - Include the 64-bit profile by itself
///   - Include also the arrays of secondary masks for the profile just found and those compiled
///     for all previous atoms
///
/// \param pos                    Topological index of the atom in question, which has exclusions
/// \param tprof                  The profile created for the atom at topological position pos
/// \param mode                   The mask mode is use (only requested when an array of secondary
///                               masks is important--this differentiates the way in which the
///                               profile parses and indexes into that array)
/// \param tprof_secondary_masks  Secondary masks array associated with the current profile
/// \param profile_indices        The developing array of profile indices for all atoms
/// \param tmp_profiles           The developing list of exclusion profiles for the system
/// \param tmp_secondary_masks    Growing array of secondary masks accumulated for all profiles
/// \{
void setProfile(int pos, ullint tprof, std::vector<int> *profile_indices,
                std::vector<ullint> *tmp_profiles);

void setProfile(int pos, ullint tprof, ullint mode,
                const std::vector<uint2> &tprof_secondary_masks, std::vector<int> *profile_indices,
                std::vector<ullint> *tmp_profiles, std::vector<uint2> *tmp_secondary_masks);
/// \}

} // namespace energy
} // namespace stormm

#endif
