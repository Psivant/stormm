// -*-c++-*-
#ifndef STORMM_ZNUMBER_H
#define STORMM_ZNUMBER_H

#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "DataTypes/stormm_vector_types.h"

namespace stormm {
namespace chemistry {

using constants::ExceptionResponse;
using constants::CaseSensitivity;

/// \brief Obtain atomic Z-numbers for a list of atoms, assuming that the masses follow the
///        weighted natural abundances 
///
/// \param masses  Masses of the atoms in question
std::vector<int> massToZNumber(const std::vector<double> &masses);

/// \brief Obain the natural abundance-averaged mass of an element based on its atomic number.
///
/// \param z  The atomic number of interest
double zNumberToNaturalMass(int z);

/// \brief Convert a series of atomic numbers to element symbols from the periodic table
///
/// Overloaded:
///   - Obtain the symbol for a single Z-number
///   - Obtain symbols for a vector of Z-numbers
///
/// \param atomic_numbers  Z-numbers of the atoms in question (can in clude 0, for a virtual site)
/// \{
char2 zNumberToSymbol(const int atomic_number);
std::vector<char2> zNumberToSymbol(const std::vector<int> &atomic_numbers);
/// \}

/// \brief Convert a series of atomic symbols to elemental Z-numbers
///
/// Overloaded:
///   - Accept an array of char2 symbols
///   - Accept an array of char4 symbols and take just the first two characters of each tuple
///
/// \param atomic_symbols  Periodic table symbols of the elements in question
/// \param policy          Course of action if the symbol cannot be placed on the periodic table
/// \{
std::vector<int> symbolToZNumber(const std::vector<char2> &atomic_symbols,
                                 CaseSensitivity capitalization = CaseSensitivity::YES,
                                 ExceptionResponse policy = ExceptionResponse::WARN);

std::vector<int> symbolToZNumber(const std::vector<char4> &atomic_symbols,
                                 CaseSensitivity capitalization = CaseSensitivity::YES,
                                 ExceptionResponse policy = ExceptionResponse::WARN);
/// \}

} // namespace chemistry
} // namespace stormm

#endif

