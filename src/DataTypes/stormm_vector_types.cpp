#include "copyright.h"
#include "stormm_vector_types.h"

namespace stormm {
namespace data_types {

//-------------------------------------------------------------------------------------------------
std::string getHpcVectorTypeName(const size_t index_code) {
  if (index_code == int2_type_index) {
    return getHpcVectorTypeName<int2>();
  }
  else if (index_code == int3_type_index) {
    return getHpcVectorTypeName<int3>();
  }
  else if (index_code == int4_type_index) {
    return getHpcVectorTypeName<int4>();
  }
  else if (index_code == double2_type_index) {
    return getHpcVectorTypeName<double2>();
  }
  else if (index_code == double3_type_index) {
    return getHpcVectorTypeName<double3>();
  }
  else if (index_code == double4_type_index) {
    return getHpcVectorTypeName<double4_16a>();
  }
  else if (index_code == float2_type_index) {
    return getHpcVectorTypeName<float2>();
  }
  else if (index_code == float3_type_index) {
    return getHpcVectorTypeName<float3>();
  }
  else if (index_code == float4_type_index) {
    return getHpcVectorTypeName<float4>();
  }
  else if (index_code == char2_type_index) {
    return getHpcVectorTypeName<char2>();
  }
  else if (index_code == char3_type_index) {
    return getHpcVectorTypeName<char3>();
  }
  else if (index_code == char4_type_index) {
    return getHpcVectorTypeName<char4>();
  }
  else if (index_code == uchar2_type_index) {
    return getHpcVectorTypeName<uchar2>();
  }
  else if (index_code == uchar3_type_index) {
    return getHpcVectorTypeName<uchar3>();
  }
  else if (index_code == uchar4_type_index) {
    return getHpcVectorTypeName<uchar4>();
  }
  else if (index_code == uint2_type_index) {
    return getHpcVectorTypeName<uint2>();
  }
  else if (index_code == uint3_type_index) {
    return getHpcVectorTypeName<uint3>();
  }
  else if (index_code == uint4_type_index) {
    return getHpcVectorTypeName<uint4>();
  }
  else if (index_code == longlong2_type_index) {
    return getHpcVectorTypeName<longlong2>();
  }
  else if (index_code == longlong3_type_index) {
    return getHpcVectorTypeName<longlong3>();
  }
  else if (index_code == longlong4_type_index) {
    return getHpcVectorTypeName<longlong4_16a>();
  }
  else if (index_code == ulonglong2_type_index) {
    return getHpcVectorTypeName<longlong3>();
  }
  else if (index_code == ulonglong3_type_index) {
    return getHpcVectorTypeName<ulonglong3>();
  }
  else if (index_code == ulonglong4_type_index) {
    return getHpcVectorTypeName<ulonglong4_16a>();
  }
  else if (index_code == short2_type_index) {
    return getHpcVectorTypeName<short2>();
  }
  else if (index_code == short3_type_index) {
    return getHpcVectorTypeName<short3>();
  }
  else if (index_code == short4_type_index) {
    return getHpcVectorTypeName<short4>();
  }
  else if (index_code == ushort2_type_index) {
    return getHpcVectorTypeName<ushort2>();
  }
  else if (index_code == ushort3_type_index) {
    return getHpcVectorTypeName<ushort3>();
  }
  else if (index_code == ushort4_type_index) {
    return getHpcVectorTypeName<ushort4>();
  }
  else if (index_code == int95t_type_index) {
    return getHpcVectorTypeName<int95_t>();
  }
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
  else if (index_code == cufftz_type_index) {
    return getHpcVectorTypeName<cufftDoubleComplex>();
  }
  else if (index_code == cufftc_type_index) {
    return getHpcVectorTypeName<cufftComplex>();
  }
#  endif
#endif
  else if (index_code == stdz_type_index) {
    return getHpcVectorTypeName<std::complex<double>>();
  }
  else if (index_code == stdc_type_index) {
    return getHpcVectorTypeName<std::complex<float>>();
  }
  __builtin_unreachable();
}
  
//-------------------------------------------------------------------------------------------------
uint short2ToUint(const short2 ituple) {
  return (static_cast<uint>(ituple.x) | (static_cast<uint>(ituple.y) << 16));
}

//-------------------------------------------------------------------------------------------------
uint ushort2ToUint(const ushort2 ituple) {
  return (static_cast<uint>(ituple.x) | (static_cast<uint>(ituple.y) << 16));
}

//-------------------------------------------------------------------------------------------------
uint char4ToUint(const char4 ctuple) {
  return (static_cast<uint>(ctuple.x) | (static_cast<uint>(ctuple.y) << 8) |
          (static_cast<uint>(ctuple.z) << 16) | (static_cast<uint>(ctuple.w) << 24));
}

//-------------------------------------------------------------------------------------------------
uint uchar4ToUint(const uchar4 ctuple) {
  return (static_cast<uint>(ctuple.x) | (static_cast<uint>(ctuple.y) << 8) |
          (static_cast<uint>(ctuple.z) << 16) | (static_cast<uint>(ctuple.w) << 24));
}

//-------------------------------------------------------------------------------------------------
short2 uintToShort2(const uint val) {
  return { static_cast<short>(val & 0xffff), static_cast<short>((val >> 16) & 0xffff) };
}

//-------------------------------------------------------------------------------------------------
ushort2 uintToUshort2(const uint val) {
  return { static_cast<ushort>(val & 0xffff), static_cast<ushort>((val >> 16) & 0xffff) };
}

//-------------------------------------------------------------------------------------------------
char4 uintToChar4(const uint val) {
  return { static_cast<char>(val & 0xff), static_cast<char>((val >> 8) & 0xff),
           static_cast<char>((val >> 16) & 0xff), static_cast<char>((val >> 24) & 0xff) };
}

//-------------------------------------------------------------------------------------------------
uchar4 uintToUchar4(const uint val) {
  return { static_cast<uchar>(val & 0xff), static_cast<uchar>((val >> 8) & 0xff),
           static_cast<uchar>((val >> 16) & 0xff), static_cast<uchar>((val >> 24) & 0xff) };
}

//-------------------------------------------------------------------------------------------------
bool operator==(const short2 lhs, const short2 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const short2 lhs, const short2 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y);
}

//-------------------------------------------------------------------------------------------------
bool operator==(const short3 lhs, const short3 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const short3 lhs, const short3 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y || lhs.z != rhs.z);
}

//-------------------------------------------------------------------------------------------------
bool operator==(const short4 lhs, const short4 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const short4 lhs, const short4 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y || lhs.z != rhs.z || lhs.w != rhs.w);
}

//-------------------------------------------------------------------------------------------------
bool operator==(const ushort2 lhs, const ushort2 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const ushort2 lhs, const ushort2 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y);
}

//-------------------------------------------------------------------------------------------------
bool operator==(const ushort3 lhs, const ushort3 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const ushort3 lhs, const ushort3 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y || lhs.z != rhs.z);
}

//-------------------------------------------------------------------------------------------------
bool operator==(const ushort4 lhs, const ushort4 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const ushort4 lhs, const ushort4 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y || lhs.z != rhs.z || lhs.w != rhs.w);
}

//-------------------------------------------------------------------------------------------------
bool operator==(const int2 lhs, const int2 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const int2 lhs, const int2 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y);
}

//-------------------------------------------------------------------------------------------------
bool operator==(const int3 lhs, const int3 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const int3 lhs, const int3 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y || lhs.z != rhs.z);
}

//-------------------------------------------------------------------------------------------------
bool operator==(const int4 lhs, const int4 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const int4 lhs, const int4 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y || lhs.z != rhs.z || lhs.w != rhs.w);
}

//-------------------------------------------------------------------------------------------------
bool operator==(const uint2 lhs, const uint2 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const uint2 lhs, const uint2 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y);
}

//-------------------------------------------------------------------------------------------------
bool operator==(const uint3 lhs, const uint3 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const uint3 lhs, const uint3 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y || lhs.z != rhs.z);
}

//-------------------------------------------------------------------------------------------------
bool operator==(const uint4 lhs, const uint4 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const uint4 lhs, const uint4 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y || lhs.z != rhs.z || lhs.w != rhs.w);
}

//-------------------------------------------------------------------------------------------------
bool operator==(const llint2 lhs, const llint2 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const llint2 lhs, const llint2 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y);
}

//-------------------------------------------------------------------------------------------------
bool operator==(const llint3 lhs, const llint3 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const llint3 lhs, const llint3 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y || lhs.z != rhs.z);
}

//-------------------------------------------------------------------------------------------------
bool operator==(const llint4 lhs, const llint4 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const llint4 lhs, const llint4 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y || lhs.z != rhs.z || lhs.w != rhs.w);
}

//-------------------------------------------------------------------------------------------------
bool operator==(const ullint2 lhs, const ullint2 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const ullint2 lhs, const ullint2 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y);
}

//-------------------------------------------------------------------------------------------------
bool operator==(const ullint3 lhs, const ullint3 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const ullint3 lhs, const ullint3 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y || lhs.z != rhs.z);
}

//-------------------------------------------------------------------------------------------------
bool operator==(const ullint4 lhs, const ullint4 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const ullint4 lhs, const ullint4 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y || lhs.z != rhs.z || lhs.w != rhs.w);
}

//-------------------------------------------------------------------------------------------------
bool operator<(const char2 lhs, const char2 rhs) {
  const char4 lhs4 = { lhs.x, lhs.y, ' ', ' ' };
  const char4 rhs4 = { rhs.x, rhs.y, ' ', ' ' };
  return (char4ToUint(lhs4) < char4ToUint(rhs4));
}

//-------------------------------------------------------------------------------------------------
bool operator<=(const char2 lhs, const char2 rhs) {
  const char4 lhs4 = { lhs.x, lhs.y, ' ', ' ' };
  const char4 rhs4 = { rhs.x, rhs.y, ' ', ' ' };
  return (char4ToUint(lhs4) <= char4ToUint(rhs4));
}

//-------------------------------------------------------------------------------------------------
bool operator>(const char2 lhs, const char2 rhs) {
  const char4 lhs4 = { lhs.x, lhs.y, ' ', ' ' };
  const char4 rhs4 = { rhs.x, rhs.y, ' ', ' ' };
  return (char4ToUint(lhs4) > char4ToUint(rhs4));
}

//-------------------------------------------------------------------------------------------------
bool operator>=(const char2 lhs, const char2 rhs) {
  const char4 lhs4 = { lhs.x, lhs.y, ' ', ' ' };
  const char4 rhs4 = { rhs.x, rhs.y, ' ', ' ' };
  return (char4ToUint(lhs4) >= char4ToUint(rhs4));
}

//-------------------------------------------------------------------------------------------------
bool operator==(const char2 lhs, const char2 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const char2 lhs, const char2 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y);
}

//-------------------------------------------------------------------------------------------------
bool operator<(const char3 lhs, const char3 rhs) {
  const char4 lhs4 = { lhs.x, lhs.y, lhs.z, ' ' };
  const char4 rhs4 = { rhs.x, rhs.y, rhs.z, ' ' };
  return (char4ToUint(lhs4) < char4ToUint(rhs4));
}

//-------------------------------------------------------------------------------------------------
bool operator<=(const char3 lhs, const char3 rhs) {
  const char4 lhs4 = { lhs.x, lhs.y, lhs.z, ' ' };
  const char4 rhs4 = { rhs.x, rhs.y, rhs.z, ' ' };
  return (char4ToUint(lhs4) <= char4ToUint(rhs4));
}

//-------------------------------------------------------------------------------------------------
bool operator>(const char3 lhs, const char3 rhs) {
  const char4 lhs4 = { lhs.x, lhs.y, lhs.z, ' ' };
  const char4 rhs4 = { rhs.x, rhs.y, rhs.z, ' ' };
  return (char4ToUint(lhs4) > char4ToUint(rhs4));
}

//-------------------------------------------------------------------------------------------------
bool operator>=(const char3 lhs, const char3 rhs) {
  const char4 lhs4 = { lhs.x, lhs.y, lhs.z, ' ' };
  const char4 rhs4 = { rhs.x, rhs.y, rhs.z, ' ' };
  return (char4ToUint(lhs4) >= char4ToUint(rhs4));
}

//-------------------------------------------------------------------------------------------------
bool operator==(const char3 lhs, const char3 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const char3 lhs, const char3 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y || lhs.z != rhs.z);
}

//-------------------------------------------------------------------------------------------------
bool operator<(const char4 lhs, const char4 rhs) {
  return (char4ToUint(lhs) < char4ToUint(rhs));
}

//-------------------------------------------------------------------------------------------------
bool operator<=(const char4 lhs, const char4 rhs) {
  return (char4ToUint(lhs) <= char4ToUint(rhs));
}

//-------------------------------------------------------------------------------------------------
bool operator>(const char4 lhs, const char4 rhs) {
  return (char4ToUint(lhs) > char4ToUint(rhs));
}

//-------------------------------------------------------------------------------------------------
bool operator>=(const char4 lhs, const char4 rhs) {
  return (char4ToUint(lhs) >= char4ToUint(rhs));
}

//-------------------------------------------------------------------------------------------------
bool operator==(const char4 lhs, const char4 rhs) {
  return (lhs.x == rhs.x && lhs.y == rhs.y && lhs.z == rhs.z && lhs.w == rhs.w);
}

//-------------------------------------------------------------------------------------------------
bool operator!=(const char4 lhs, const char4 rhs) {
  return (lhs.x != rhs.x || lhs.y != rhs.y || lhs.z != rhs.z || lhs.w != rhs.w);
}

} // namespace data_types
} // namespace stormm
