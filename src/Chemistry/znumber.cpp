#include "copyright.h"
#include "Parsing/parse.h"
#include "UnitTesting/unit_test.h"
#include "znumber.h"
#include "periodic_table.h"

namespace stormm {
namespace chemistry {

using parse::uppercase;
using testing::Approx;
using testing::ComparisonType;  
  
//-------------------------------------------------------------------------------------------------
std::vector<int> massToZNumber(const std::vector<double> &masses) {

  // Create a vector of approximate masses for comparisons
  std::vector<Approx> natural_masses;
  for (int i = 0; i < element_maximum_count; i++) {
    natural_masses.push_back(Approx(elemental_masses[i], 0.01));
  }

  // Loop over all atoms in the system
  const int natom = masses.size();
  std::vector<int> result(natom, -1);
  for (int i = 0; i < natom; i++) {
    const double mass_i = masses[i];
    for (int j = 0; j < element_maximum_count; j++) {
      if (natural_masses[j].test(mass_i)) {
	result[i] = j;
	break;
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
double zNumberToNaturalMass(const int z) {
  if (z < 0 || z >= element_maximum_count) {
    rtErr("No mass is known for element " + std::to_string(z) + ".", "zNumberToNaturalMass");
  }
  return elemental_masses[z];
}

//-------------------------------------------------------------------------------------------------
char2 zNumberToSymbol(const int atomic_number) {
  if (atomic_number < 0) {
    rtWarn("Atomic number " + std::to_string(atomic_number) + " is invalid and will be "
           "indicated by symbol XX.", "zNumberToSymbol");
  }
  else if (atomic_number > element_maximum_count) {
    rtWarn("Atomic number " + std::to_string(atomic_number) + " is beyond the scope of "
           "stable elements covered by STORMM and will be indicated by symbol XX.",
           "zNumberToSymbol");      
  }
  return elemental_symbols[atomic_number];
}

//-------------------------------------------------------------------------------------------------
std::vector<char2> zNumberToSymbol(const std::vector<int> &atomic_numbers) {
  const int natom = atomic_numbers.size();
  std::vector<char2> symbs(natom);
  for (int i = 0; i < natom; i++) {
    symbs[i] = zNumberToSymbol(atomic_numbers[i]);
  }
  return symbs;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> symbolToZNumber(const std::vector<char2> &atomic_symbols,
                                 const CaseSensitivity capitalization,
                                 const ExceptionResponse policy) {
  const int natom = atomic_symbols.size();
  std::vector<int> result(natom, -1);
  for (int i = 0; i < natom; i++) {
    const char2 tmps = atomic_symbols[i];
    for (int j = 0; j < element_maximum_count; j++) {
      switch (capitalization) {
      case CaseSensitivity::YES:
      case CaseSensitivity::AUTOMATIC:
        if (tmps.x == elemental_symbols[j].x && tmps.y == elemental_symbols[j].y) {
          result[i] = j;
        }
        break;
      case CaseSensitivity::NO:
        if (uppercase(tmps.x) == uppercase(elemental_symbols[j].x) &&
            uppercase(tmps.y) == uppercase(elemental_symbols[j].y)) {
          result[i] = j;
        }
        break;
      }
    }
    if (result[i] < 0) {
      std::string tsymbol;
      tsymbol += tmps.x;
      tsymbol += tmps.y;
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("No atomic symbol " + tsymbol + " is known.  Atom " + std::to_string(i) + " cannot "
              "be represented as an element with mass.");
      case ExceptionResponse::WARN:
        rtWarn("No atomic symbol " + tsymbol + " is known.  Atom " + std::to_string(i) +
               "will be represented as a virtual site (symbol VS) with atomic number 0.",
               "symbolToZNumber");
        result[i] = 0;
        break;
      case ExceptionResponse::SILENT:
        result[i] = 0;
        break;
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> symbolToZNumber(const std::vector<char4> &atomic_symbols,
                                 const CaseSensitivity capitalization,
                                 const ExceptionResponse policy) {
  const size_t nsymb = atomic_symbols.size();
  std::vector<char2> abridged_symbols(nsymb);
  for (size_t i = 0LLU; i < nsymb; i++) {

    // Remove leading white space
    char4 asymb_tmp = atomic_symbols[i];
    int nws = 0;
    while (asymb_tmp.x == ' ' && nws < 3) {
      asymb_tmp.x = asymb_tmp.y;
      asymb_tmp.y = asymb_tmp.z;
      asymb_tmp.z = asymb_tmp.w;
      asymb_tmp.w = ' ';
      nws++;
    }
    abridged_symbols[i].x = atomic_symbols[i].x;
    abridged_symbols[i].y = atomic_symbols[i].y;
  }
  return symbolToZNumber(abridged_symbols, capitalization, policy);
}

} // namespace chemistry
} // namespace stormm
