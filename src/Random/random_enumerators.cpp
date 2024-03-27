#include "copyright.h"
#include "random_enumerators.h"

namespace stormm {
namespace random {

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const RandomNumberKind input) {
  switch (input) {
  case RandomNumberKind::UNIFORM:
    return std::string("UNIFORM");
  case RandomNumberKind::GAUSSIAN:
    return std::string("GAUSSIAN");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const RandomAlgorithm input) {
  switch (input) {
  case RandomAlgorithm::XOROSHIRO_128P:
    return std::string("XOROSHIRO_128P");
  case RandomAlgorithm::XOSHIRO_256PP:
    return std::string("XOSHIRO_256PP");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const RngFillMode input) {
  switch (input) {
  case RngFillMode::COLUMNS:
    return std::string("COLUMNS");
  case RngFillMode::ROWS:
    return std::string("ROWS");
  }
  __builtin_unreachable();
}

} // namespace random
} // namespace stormm
