// -*-c++-*-
#ifndef STORMM_SCALING_H
#define STORMM_SCALING_H

#include "copyright.h"
#include "DataTypes/common_types.h"

namespace stormm {
namespace constants {

/// \brief Teh convenient, round number is 1024, not 1000.
/// \{
constexpr llint kilo = 1024LL;
constexpr llint mega = kilo * kilo;
constexpr llint giga = mega * kilo;
constexpr size_t kilo_zu = 1024LLU;
constexpr size_t mega_zu = kilo_zu * kilo_zu;
constexpr size_t giga_zu = mega_zu * kilo_zu;
/// \}

/// \brief a value bigger than tiny but still small enough to often be negligible
constexpr double small = 1.0e-8;
constexpr float small_f = small;
  
/// \brief The default tiny value, below which most quantities are unimportant
constexpr double tiny = 1.0e-10;
constexpr float tiny_f = tiny;

/// \brief An even tinier value, for the most stringent comparisons
constexpr double verytiny = 1.0e-12;
constexpr float verytiny_f = verytiny;

/// \brief Sizes of important data types
/// \{
constexpr int int_byte_count_int = sizeof(int);
constexpr int int_bit_count_int = int_byte_count_int * 8;
constexpr int uint_byte_count_int = sizeof(unsigned int);
constexpr int uint_bit_count_int = uint_byte_count_int * 8;
constexpr int llint_byte_count_int = sizeof(long long int);
constexpr int llint_bit_count_int = llint_byte_count_int * 8;
constexpr int ullint_byte_count_int = sizeof(unsigned long long int);
constexpr int ullint_bit_count_int = ullint_byte_count_int * 8;
constexpr int ushort_byte_count_int = sizeof(unsigned short int);
constexpr int ushort_bit_count_int = ushort_byte_count_int * 8;
constexpr size_t int_byte_count_zu = sizeof(int);
constexpr size_t int_bit_count_zu = int_byte_count_int * 8;
constexpr size_t uint_byte_count_zu = sizeof(unsigned int);
constexpr size_t uint_bit_count_zu = uint_byte_count_zu * 8LLU;
constexpr size_t llint_byte_count_zu = sizeof(long long int);
constexpr size_t llint_bit_count_zu = llint_byte_count_int * 8LLU;
constexpr size_t ullint_byte_count_zu = sizeof(unsigned long long int);
constexpr size_t ullint_bit_count_zu = ullint_byte_count_zu * 8LLU;
constexpr size_t ushort_byte_count_zu = sizeof(unsigned short int);
constexpr size_t ushort_bit_count_zu = ushort_byte_count_zu * 8LLU;
/// \}

/// \brief Conversion factor for obtaining the velocity change after a time interval (in units of
///        femtoseconds).
/// \{
constexpr double time_step_factor = 0.0204548282808729546544679323005766491405665874481201171875;
constexpr float time_step_factor_f = time_step_factor;
/// \}

} // namespace constants
} // namespace stormm

namespace stormm {
using constants::int_byte_count_int;
using constants::int_bit_count_int;
using constants::uint_byte_count_int;
using constants::uint_bit_count_int;
using constants::ullint_byte_count_int;
using constants::ullint_bit_count_int;
using constants::ushort_byte_count_int;
using constants::ushort_bit_count_int;
using constants::int_byte_count_zu;
using constants::int_bit_count_zu;
using constants::uint_byte_count_zu;
using constants::uint_bit_count_zu;
using constants::ullint_byte_count_zu;
using constants::ullint_bit_count_zu;
using constants::ushort_byte_count_zu;
using constants::ushort_bit_count_zu;
} // namespace stormm

#endif
