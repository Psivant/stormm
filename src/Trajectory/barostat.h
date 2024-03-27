// -*-c++-*-
#ifndef STORMM_BAROSTAT_H
#define STORMM_BAROSTAT_H

#include <string>
#include "copyright.h"

namespace stormm {
namespace trajectory {

/// \brief Enumerate the available barostat types
enum class BarostatKind {
  NONE, MONTE_CARLO, BERENDSEN
};

/// \brief Store the parameters for a simulation barostat.  Includes Monte-Carlo and Berendesen
///        barostats with atomic virials.
class Barostat {
public:

  /// \brief The constructor can be blank (implying a barostat of kind NONE), take a specific
  ///        barostat kind (implying default settings for that type of piston), or take specific
  ///        settings for all aspects of the barostat.
  /// \{
  Barostat();
  Barostat(BarostatKind kind);
  /// \}

  /// \brief Get the kind of barostat
  BarostatKind getKind() const;
  
private:
  BarostatKind kind;      ///< The type of barostat
};

/// \brief Return the name of the barostat choice (another enumerator string conversion function)
///
/// \param kind  The type of barostat
std::string getBarostatName(BarostatKind kind);

} // namespace trajectory
} // namespace stormm

#endif
