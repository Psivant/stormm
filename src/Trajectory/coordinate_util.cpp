#include "copyright.h"
#include "Reporting/error_format.h"
#include "coordinate_util.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
std::string nameCoordinateType(const size_t ct_coords) {
  if (ct_coords == cf_type_index) {
    return std::string("CoordinateFrame");
  }
  else if (ct_coords == ps_type_index) {
    return std::string("PhaseSpace");
  }
  else if (ct_coords == cs_dtype_index) {
    return std::string("CoordinateSeries<float64_t>");
  }  
  else if (ct_coords == cs_ftype_index) {
    return std::string("CoordinateSeries<float32_t>");
  }  
  else if (ct_coords == cs_stype_index) {
    return std::string("CoordinateSeries<int16_t>");
  }  
  else if (ct_coords == cs_itype_index) {
    return std::string("CoordinateSeries<int32_t>");
  }  
  else if (ct_coords == cs_ltype_index) {
    return std::string("CoordinateSeries<int64_t>");
  }  
  else if (ct_coords == poly_ps_type_index) {
    return std::string("PhaseSpaceSynthesis");
  }  
  else if (ct_coords == cdns_type_index) {
    return std::string("Condensate");
  }
  else {
    return std::string("Unrecognized coordinate object type.", "nameCoordinateType");
  }
  __builtin_unreachable();
}
  
} // namespace trajectory
} // namespace stormm
