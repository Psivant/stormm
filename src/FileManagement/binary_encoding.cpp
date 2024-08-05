#include "copyright.h"
#include "binary_encoding.h"

namespace stormm {
namespace diskutil {

//-------------------------------------------------------------------------------------------------
BinaryFileContent::BinaryFileContent(const size_t ct_vptr_in, const size_t length_in,
                                     const void* vptr_in, const bool repeating_in,
                                     const std::string &description_in,
                                     const HybridTargetLevel tier_in) :
    ct_vptr{ct_vptr_in}, type_code{codifyTypeIndex(ct_vptr, true, repeating_in)},
    length{length_in}, description{description_in}, vptr{const_cast<void*>(vptr_in)},
    tier{tier_in}, lower_bound{0}, upper_bound{0}
{}

//-------------------------------------------------------------------------------------------------
void BinaryFileContent::setLowerBound(const llint lower_bound_in) {
  lower_bound = lower_bound_in;
}

//-------------------------------------------------------------------------------------------------
void BinaryFileContent::setUpperBound(const llint upper_bound_in) {
  upper_bound = upper_bound_in;
}

//-------------------------------------------------------------------------------------------------
BinaryFileKey::BinaryFileKey(const std::string &description_in) :
    description{description_in}, tiers{}, singleton_items{}, array_items{}, toc_bytes{0}
{}

//-------------------------------------------------------------------------------------------------
BinaryFileKey::BinaryFileKey(const CoordinateFrame *cf,
                             const std::string &description_in, const HybridTargetLevel tier_in) :
    BinaryFileKey(description_in)
{
  const CoordinateFrameReader cfr = cf->data(tier_in);
  addArray(cfr.xcrd, cfr.natom, false, "Particle X coordinates", tier_in);
  addArray(cfr.ycrd, cfr.natom, false, "Particle Y coordinates", tier_in);
  addArray(cfr.zcrd, cfr.natom, false, "Particle Z coordinates", tier_in);
  addArray(cfr.umat, 9, false, "Transform into unit cell fractional space", tier_in);
  addArray(cfr.invu, 9, false, "Transform into Cartesian space", tier_in);
  addArray(cfr.boxdim, 6, false, "Box dimensions X/Y/Z + alpha/beta/gamma", tier_in);
}

//-------------------------------------------------------------------------------------------------
BinaryFileKey::BinaryFileKey(const CoordinateFrame &cf,
                             const std::string &description_in, const HybridTargetLevel tier_in) :
    BinaryFileKey(cf.getSelfPointer(), description_in, tier_in)
{}

//-------------------------------------------------------------------------------------------------
uint16_t codifyTypeIndex(const size_t ct, const bool is_array, const bool is_repeatable) {
  int number_kind, is_signed, tuple_count, element_size;
  if (ct == float_type_index) {
    number_kind = 1;
    is_signed = 1;
    tuple_count = 1;
    element_size = sizeof(float);
  }
  else if (ct == double_type_index) {
    number_kind = 1;
    is_signed = 1;
    tuple_count = 1;
    element_size = sizeof(double);
  }
  else if (ct == longdouble_type_index) {
    number_kind = 1;
    is_signed = 1;
    tuple_count = 1;
    element_size = sizeof(long double);
  }
  else if (ct == short_type_index) {
    number_kind = 0;
    is_signed = 1;
    tuple_count = 1;
    element_size = sizeof(short int);
  }
  else if (ct == int_type_index) {
    number_kind = 0;
    is_signed = 1;
    tuple_count = 1;
    element_size = sizeof(int);
  }
  else if (ct == llint_type_index) {
    number_kind = 0;
    is_signed = 1;
    tuple_count = 1;
    element_size = sizeof(llint);
  }
  else if (ct == ushort_type_index) {
    number_kind = 0;
    is_signed = 0;
    tuple_count = 1;
    element_size = sizeof(unsigned short int);
  }
  else if (ct == uint_type_index) {
    number_kind = 0;
    is_signed = 0;
    tuple_count = 1;
    element_size = sizeof(unsigned int);
  }
  else if (ct == ullint_type_index) {
    number_kind = 0;
    is_signed = 0;
    tuple_count = 1;
    element_size = sizeof(unsigned long long int);
  }
  else if (ct == char_type_index) {
    number_kind = 2;
    is_signed = 1;
    tuple_count = 1;
    element_size = sizeof(char);
  }
  else if (ct == uchar_type_index) {
    number_kind = 2;
    is_signed = 0;
    tuple_count = 1;
    element_size = sizeof(unsigned char);
  }
  else if (ct == float2_type_index) {
    number_kind = 1;
    is_signed = 1;
    tuple_count = 2;
    element_size = sizeof(float);
  }
  else if (ct == float3_type_index) {
    number_kind = 1;
    is_signed = 1;
    tuple_count = 3;
    element_size = sizeof(float);
  }
  else if (ct == float4_type_index) {
    number_kind = 1;
    is_signed = 1;
    tuple_count = 4;
    element_size = sizeof(float);
  }
  else if (ct == double2_type_index) {
    number_kind = 1;
    is_signed = 1;
    tuple_count = 2;
    element_size = sizeof(double);
  }
  else if (ct == double3_type_index) {
    number_kind = 1;
    is_signed = 1;
    tuple_count = 3;
    element_size = sizeof(double);
  }
  else if (ct == double4_type_index) {
    number_kind = 1;
    is_signed = 1;
    tuple_count = 4;
    element_size = sizeof(double);
  }
  else if (ct == short2_type_index) {
    number_kind = 0;
    is_signed = 1;
    tuple_count = 2;
    element_size = sizeof(short int);
  }
  else if (ct == short3_type_index) {
    number_kind = 0;
    is_signed = 1;
    tuple_count = 3;
    element_size = sizeof(short int);
  }
  else if (ct == short4_type_index) {
    number_kind = 0;
    is_signed = 1;
    tuple_count = 4;
    element_size = sizeof(short int);
  }
  else if (ct == int2_type_index) {
    number_kind = 0;
    is_signed = 1;
    tuple_count = 2;
    element_size = sizeof(int);
  }
  else if (ct == int3_type_index) {
    number_kind = 0;
    is_signed = 1;
    tuple_count = 3;
    element_size = sizeof(int);
  }
  else if (ct == int4_type_index) {
    number_kind = 0;
    is_signed = 1;
    tuple_count = 4;
    element_size = sizeof(int);
  }
  else if (ct == longlong2_type_index) {
    number_kind = 0;
    is_signed = 1;
    tuple_count = 2;
    element_size = sizeof(llint);
  }
  else if (ct == longlong3_type_index) {
    number_kind = 0;
    is_signed = 1;
    tuple_count = 3;
    element_size = sizeof(llint);
  }
  else if (ct == longlong4_type_index) {
    number_kind = 0;
    is_signed = 1;
    tuple_count = 4;
    element_size = sizeof(llint);
  }
  else if (ct == ushort2_type_index) {
    number_kind = 0;
    is_signed = 0;
    tuple_count = 2;
    element_size = sizeof(ushort);
  }
  else if (ct == ushort3_type_index) {
    number_kind = 0;
    is_signed = 0;
    tuple_count = 3;
    element_size = sizeof(ushort);
  }
  else if (ct == ushort4_type_index) {
    number_kind = 0;
    is_signed = 0;
    tuple_count = 4;
    element_size = sizeof(ushort);
  }
  else if (ct == uint2_type_index) {
    number_kind = 0;
    is_signed = 0;
    tuple_count = 2;
    element_size = sizeof(uint);
  }
  else if (ct == uint3_type_index) {
    number_kind = 0;
    is_signed = 0;
    tuple_count = 3;
    element_size = sizeof(uint);
  }
  else if (ct == uint4_type_index) {
    number_kind = 0;
    is_signed = 0;
    tuple_count = 4;
    element_size = sizeof(uint);
  }
  else if (ct == ulonglong2_type_index) {
    number_kind = 0;
    is_signed = 0;
    tuple_count = 2;
    element_size = sizeof(ullint);
  }
  else if (ct == ulonglong3_type_index) {
    number_kind = 0;
    is_signed = 0;
    tuple_count = 3;
    element_size = sizeof(ullint);
  }
  else if (ct == ulonglong4_type_index) {
    number_kind = 0;
    is_signed = 0;
    tuple_count = 4;
    element_size = sizeof(ullint);
  }
  else if (ct == char2_type_index) {
    number_kind = 2;
    is_signed = 1;
    tuple_count = 2;
    element_size = sizeof(char);
  }
  else if (ct == char3_type_index) {
    number_kind = 2;
    is_signed = 1;
    tuple_count = 3;
    element_size = sizeof(char);
  }
  else if (ct == char4_type_index) {
    number_kind = 2;
    is_signed = 1;
    tuple_count = 4;
    element_size = sizeof(char);
  }
  else if (ct == uchar2_type_index) {
    number_kind = 2;
    is_signed = 0;
    tuple_count = 2;
    element_size = sizeof(unsigned char);
  }
  else if (ct == uchar3_type_index) {
    number_kind = 2;
    is_signed = 0;
    tuple_count = 3;
    element_size = sizeof(unsigned char);
  }
  else if (ct == uchar4_type_index) {
    number_kind = 2;
    is_signed = 0;
    tuple_count = 4;
    element_size = sizeof(unsigned char);
  }
  const uint16_t result = number_kind + (is_signed << 2) + (tuple_count << 3) +
                          (static_cast<int>(is_array) << 7) +
                          (static_cast<int>(is_repeatable) << 8) + (element_size << 9);
  return result;
}

//-------------------------------------------------------------------------------------------------
void BinaryFileKey::addArray(const void* vptr, const size_t length, const size_t ct_vptr,
                             const bool repeating_in, const std::string &description_in,
                             const HybridTargetLevel tier_in) {
  array_items.emplace_back(ct_vptr, length, vptr, repeating_in, description_in, tier_in);
}
  
} // namespace disktuil
} // namespace stormm
