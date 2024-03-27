#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/parse.h"
#include "math_enumerators.h"

namespace stormm {
namespace stmath {

using constants::CaseSensitivity;
using parse::strcmpCased;

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const FunctionLevel input) {
  switch (input) {
  case FunctionLevel::VALUE:
    return std::string("VALUE");
  case FunctionLevel::DX:
    return std::string("DX");
  case FunctionLevel::DY:
    return std::string("DY");
  case FunctionLevel::DZ:
    return std::string("DZ");
  case FunctionLevel::DXX:
    return std::string("DXX");
  case FunctionLevel::DXY:
    return std::string("DXY");
  case FunctionLevel::DXZ:
    return std::string("DXZ");
  case FunctionLevel::DYY:
    return std::string("DYY");
  case FunctionLevel::DYZ:
    return std::string("DYZ");
  case FunctionLevel::DZZ:
    return std::string("DZZ");
  case FunctionLevel::DXXX:
    return std::string("DXXX");
  case FunctionLevel::DXXY:
    return std::string("DXXY");
  case FunctionLevel::DXXZ:
    return std::string("DXXZ");
  case FunctionLevel::DXYY:
    return std::string("DXYY");
  case FunctionLevel::DXYZ:
    return std::string("DXYZ");
  case FunctionLevel::DXZZ:
    return std::string("DXZZ");
  case FunctionLevel::DYYY:
    return std::string("DYYY");
  case FunctionLevel::DYYZ:
    return std::string("DYYZ");
  case FunctionLevel::DYZZ:
    return std::string("DYZZ");
  case FunctionLevel::DZZZ:
    return std::string("DZZZ");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const Interpolant input) {
  switch (input) {
  case Interpolant::SMOOTHNESS:
    return std::string("SMOOTHNESS");
  case Interpolant::FUNCTION_VALUE:
    return std::string("FUNCTION_VALUE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const LimitApproach input) {
  switch (input) {
  case LimitApproach::BELOW:
    return std::string("BELOW");
  case LimitApproach::ABOVE:
    return std::string("ABOVE");    
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const LogSplineForm input) {
  switch (input) {
  case LogSplineForm::ELEC_PME_DIRECT:
    return std::string("ELEC_PME_DIRECT");
  case LogSplineForm::ELEC_PME_DIRECT_EXCL:
    return std::string("ELEC_PME_DIRECT_EXCL");
  case LogSplineForm::DELEC_PME_DIRECT:
    return std::string("DELEC_PME_DIRECT");
  case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    return std::string("DELEC_PME_DIRECT_EXCL");
  case LogSplineForm::CUSTOM:
    return std::string("CUSTOM");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const TableIndexing input) {
  switch (input) {
  case TableIndexing::ARG:
    return std::string("ARG");
  case TableIndexing::SQUARED_ARG:
    return std::string("SQUARED_ARG");
  case TableIndexing::ARG_OFFSET:
    return std::string("ARG_OFFSET");
  case TableIndexing::SQ_ARG_OFFSET:
    return std::string("SQ_ARG_OFFSET");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const BasisFunctions input) {
  switch (input) {
  case BasisFunctions::MIXED_FRACTIONS:
    return std::string("MIXED_FRACTIONS");
  case BasisFunctions::POLYNOMIAL:
    return std::string("POLYNOMIAL");
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const BSplineUnity input) {
  switch (input) {
  case BSplineUnity::CENTER_FILL:
    return std::string("CENTER_FILL");
  case BSplineUnity::NONE:
    return std::string("NONE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const FFTMode input) {
  switch (input) {
  case FFTMode::IN_PLACE:
    return std::string("IN_PLACE");
  case FFTMode::OUT_OF_PLACE:
    return std::string("OUT_OF_PLACE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
FFTMode translateFFTMode(const std::string &input) {
  if (strcmpCased(input, std::string("in_place"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("in-place"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("inplace"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("ip"), CaseSensitivity::NO)) {
    return FFTMode::IN_PLACE;
  }
  else if (strcmpCased(input, std::string("out_of_place"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("out-of-place"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("outofplace"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("oop"), CaseSensitivity::NO)) {
    return FFTMode::OUT_OF_PLACE;
  }
  else {
    rtErr("Unrecognized token \"" + input + "\".", "translateFFTMode");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
LogSplineForm translateLogSplineForm(const std::string &input) {
  if (strcmpCased(input, std::string("elec_pme"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("elecpme"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("elec_pme_direct"), CaseSensitivity::NO)) {
    return LogSplineForm::ELEC_PME_DIRECT;
  }
  else if (strcmpCased(input, std::string("delec_pme"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("delecpme"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("delec_pme_direct"), CaseSensitivity::NO)) {
    return LogSplineForm::DELEC_PME_DIRECT;
  }
  else if (strcmpCased(input, std::string("elec_pme_excl"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("elecpme_excl"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("elec_pme_direct_excl"), CaseSensitivity::NO)) {
    return LogSplineForm::ELEC_PME_DIRECT_EXCL;
  }
  else if (strcmpCased(input, std::string("delec_pme_excl"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("delecpme_excl"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("delec_pme_direct_excl"), CaseSensitivity::NO)) {
    return LogSplineForm::DELEC_PME_DIRECT_EXCL;
  }
  else if (strcmpCased(input, std::string("custom"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("customized"), CaseSensitivity::NO)) {
    return LogSplineForm::CUSTOM;
  }
  else {
    rtErr("Unrecognized token \"" + input + "\".", "translateLogSplineForm");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
TableIndexing translateTableIndexing(const std::string &input) {
  if (strcmpCased(input, std::string("arg"), CaseSensitivity::NO) ||
      strcmpCased(input, std::string("argument"), CaseSensitivity::NO)) {
    return TableIndexing::ARG;
  }
  else if (strcmpCased(input, std::string("squared_arg"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("sqarg"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("arg_squared"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("argsq"), CaseSensitivity::NO)) {
    return TableIndexing::SQUARED_ARG;
  }
  else if (strcmpCased(input, std::string("arg_offset"), CaseSensitivity::NO)) {
    return TableIndexing::ARG_OFFSET;
  }
  else if (strcmpCased(input, std::string("squared_arg_offset"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("sqarg_offset"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("sq_arg_offset"), CaseSensitivity::NO) ||
           strcmpCased(input, std::string("arg_squared_offset"), CaseSensitivity::NO)) {
    return TableIndexing::SQ_ARG_OFFSET;
  }
  else {
    rtErr("Unrecognized token \"" + input + "\".", "translateTableIndexing");
  }
  __builtin_unreachable();
}

} // namespace stormm
} // namespace stmath
