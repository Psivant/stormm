// -*-c++-*-
#ifndef STORMM_NML_RANDOM_H
#define STORMM_NML_RANDOM_H

#include "copyright.h"
#include "Parsing/textfile.h"
#include "input.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using parse::WrapTextSearch;

/// \brief Default values for molecular dynamics
/// \{
constexpr int default_random_seed    = 30965871;
#ifdef STORMM_USE_HPC
constexpr int default_random_streams = 1048576;
#else
constexpr int default_random_streams = 1;
#endif
constexpr int default_random_stride  = 64;
constexpr int default_random_warmup  = 96;
/// \}

/// \brief Object to encapsulate random number generation controls.  Like other namelist
///        encapsualtors, this object can take input file data as part of its construction, or
///        by a series of setters.  Validation of each piece of data is handled as it appears
///        either in the contructor or via setters.  Getter functions dispense the internal
///        information to any application using STORMM libraries.
class RandomControls {
public:

  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &solvent namelist
  /// \param found_nml   Indication that the namelist was found (passed back to calling function)
  /// \param policy_in   Requested error handling behavior
  /// \param wrap        Indicate that the search for a &random namelist should carry on from
  ///                    the beginning of an input file if no such namelist is found starting
  ///                    from the original starting point
  /// \{
  RandomControls(ExceptionResponse policy_in = ExceptionResponse::DIE,
                 WrapTextSearch wrap = WrapTextSearch::NO);
  RandomControls(const TextFile &tf, int *start_line, bool *found_nml,
                 ExceptionResponse policy_in = ExceptionResponse::DIE,
                 WrapTextSearch wrap = WrapTextSearch::NO);
  /// \}

  /// \brief As with other control objects, copy and move constructors, plus copy and move
  ///        assignment operators, can all take their default forms.
  /// \{
  RandomControls(const RandomControls &original) = default;
  RandomControls(RandomControls &&original) = default;
  RandomControls& operator=(const RandomControls &original) = default;
  RandomControls& operator=(RandomControls &&original) = default;
  /// \}

  /// \brief Get the random seed supplied (and possibly, as corrected) from user input.
  int getRandomSeed() const;

  /// \brief Get the number of random streams.
  int getStreamCount() const;

  /// \brief Get the quantity of random numbers produced, per stream, each time a cache is filled.
  int getProductionStride() const;

  /// \brief Get the number of warmup cycles for the prime stream generator.
  int getWarmupCycleCount() const;

  /// \brief Get the original namelist emulator object as a transcript of the user input.
  const NamelistEmulator& getTranscript() const;
  
  /// \brief Set the random seed.
  ///
  /// \param igseed_in  The new choice of random seed (will be checked for validity)
  void setRandomSeed(const int igseed_in);

  /// \brief Set the quantity of random number streams.
  ///
  /// \param streams_in  The requested number of streams
  void setStreamCount(const int streams_in);

  /// \brief Set the random number production batch size, the quantity of random numbers that each
  ///        stream will produce each time a cache is filled.
  ///
  /// \param stride_in  The requested stride size
  void setProductionStride(const int stride_in);

  /// \brief Set the number of warmup cycles through which to process the prime stream's random
  ///        generator state prior to producing numbers in that stream or spawning new streams.
  ///
  /// \param cycles_in  Requested number of cycles
  void setWarmupCycles(const int cycles_in);
  
private:
  ExceptionResponse policy;  ///< Set the behavior when bad inputs are encountered.  DIE = abort
                             ///<   program, WARN = warn the user, and likely reset to the default
                             ///<   value if one is available, SILENT = do not warn the user, but
                             ///<   also likely reset to the default value if one is available.
  int igseed;                ///< Random number generator prime stream seed (all other streams
                             ///<   follow from this one seed)
  int stream_count;          ///< Quantity of random number state vectors to keep on hand.  This,
                             ///<   times the production stride, yields the size of a random
                             ///<   number batch.
  int production_stride;     ///< Quantity of random numbers to produce, per stream, when the
                             ///<   random number cache must be refreshed.
  int warmup_cycles;         ///< Number of warmup cycles to put the generator through prior to
                             ///<   producing numbers in the prime stream, or spawning any other
                             ///<   streams.

  /// Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;
  
  /// \brief Check the random number seed to ensure that it contains enough bits set to 1.
  ///        Adjust the number of warmup cycles as necessary.
  void validateRandomSeed();

  /// \brief Check the number of random streams for validity.
  void validateStreamCount();

  /// \brief Check the random number batch size (per stream) for validity.
  void validateProductionStride();
};
  
/// \brief Produce a namelist for specifying random number generation protocols.  This subsumes the
///        igseed keyword found in the &cntrl namelist of the Amber sander program, and much more.
///        STORMM random number streams created by any GPU can be reproduced exactly on the CPU.
///
/// \param tf          Input text file to scan immediately after the namelist has been created
/// \param start_line  Line at which to begin scanning the input file for the namelist (this
///                    function will wrap back to the beginning of the TextFile object, if
///                    necessary, to find a &random namelist)
/// \param found       Indication that the namelist was found (passed back to calling function)
/// \param policy      Reaction to exceptions encountered during namelist reading
/// \param wrap        Indicate that the search for a &random namelist should carry on from the
///                    beginning of an input file if no such namelist is found starting from the
///                    original starting point
NamelistEmulator randomInput(const TextFile &tf, int *start_line, bool *found,
                             ExceptionResponse policy = ExceptionResponse::DIE,
                             WrapTextSearch wrap = WrapTextSearch::NO);

/// \brief 
///
/// \param igseed         The provided random seed
/// \param warmup_cycles  Number of warmup cycles slated to process when creating the first
///                       random number stream generator (modified if necessary and returned)
void validateRandomSeed(int igseed, int *warmup_cycles);

} // namespace namelist
} // namespace stormm

#endif
