// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace data_types {

//-------------------------------------------------------------------------------------------------
template <typename T> bool isHpcVectorType() {
  return (std::is_same_v<T, int2>          || std::is_same_v<T, int3>           ||
          std::is_same_v<T, int4>          || std::is_same_v<T, double2>        ||
          std::is_same_v<T, double3>       || std::is_same_v<T, double4_16a>    ||
          std::is_same_v<T, float2>        || std::is_same_v<T, float3>         ||
          std::is_same_v<T, float4>        || std::is_same_v<T, char2>          ||
          std::is_same_v<T, char3>         || std::is_same_v<T, char4>          ||
          std::is_same_v<T, uchar2>        || std::is_same_v<T, uchar3>         ||
          std::is_same_v<T, uchar4>        || std::is_same_v<T, uint2>          ||
          std::is_same_v<T, uint3>         || std::is_same_v<T, uint4>          ||
          std::is_same_v<T, longlong2>     || std::is_same_v<T, longlong3>      ||
          std::is_same_v<T, longlong4_16a> || std::is_same_v<T, ulonglong2>     ||
          std::is_same_v<T, ulonglong3>    || std::is_same_v<T, ulonglong4_16a> ||
          std::is_same_v<T, short2>        || std::is_same_v<T, short3>         ||
          std::is_same_v<T, short4>        || std::is_same_v<T, ushort2>        ||
          std::is_same_v<T, ushort3>       || std::is_same_v<T, ushort4>);
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool isSignedIntegralHpcVectorType() {

  // The int95_t type is not counted among the signed integer vector types because it is not
  // efficient to make Hybrid objects out of such things, and it is not among the standard CUDA
  // or HIP HPC vector tuple types.
  return (std::is_same_v<T, int2>  || std::is_same_v<T, int3>  || std::is_same_v<T, int4> ||
          std::is_same_v<T, char2> || std::is_same_v<T, char3> || std::is_same_v<T, char4> ||
          std::is_same_v<T, longlong2>     || std::is_same_v<T, longlong3> ||
          std::is_same_v<T, longlong4_16a> || std::is_same_v<T, short2>    ||
          std::is_same_v<T, short3>        || std::is_same_v<T, short4>);
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool isUnsignedIntegralHpcVectorType() {
  return (std::is_same_v<T, uint2>  || std::is_same_v<T, uint3>  || std::is_same_v<T, uint4> ||
          std::is_same_v<T, uchar2> || std::is_same_v<T, uchar3> || std::is_same_v<T, uchar4> ||
          std::is_same_v<T, ulonglong2>     || std::is_same_v<T, ulonglong3> ||
          std::is_same_v<T, ulonglong4_16a> || std::is_same_v<T, ushort2>    ||
          std::is_same_v<T, ushort3>        || std::is_same_v<T, ushort4>);
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool isFloatingPointHpcVectorType() {
  return (std::is_same_v<T, double2>     || std::is_same_v<T, double3> ||
          std::is_same_v<T, double4_16a> || std::is_same_v<T, float2>  ||
          std::is_same_v<T, float3>      || std::is_same_v<T, float4>);
}

//-------------------------------------------------------------------------------------------------
template <typename T> int getHpcVectorTypeSize() {
  if (std::is_same_v<T, int2>      || std::is_same_v<T, uint2>      ||
      std::is_same_v<T, double2>   || std::is_same_v<T, float2>     ||
      std::is_same_v<T, char2>     || std::is_same_v<T, uchar2>     ||
      std::is_same_v<T, longlong2> || std::is_same_v<T, ulonglong2> ||
      std::is_same_v<T, short2>    || std::is_same_v<T, ushort2>) {
    return 2;
  }
  else if (std::is_same_v<T, int3> || std::is_same_v<T, uint3>      ||
      std::is_same_v<T, double3>   || std::is_same_v<T, float3>     ||
      std::is_same_v<T, char3>     || std::is_same_v<T, uchar3>     ||
      std::is_same_v<T, longlong3> || std::is_same_v<T, ulonglong3> ||
      std::is_same_v<T, short3>    || std::is_same_v<T, ushort3>) {
    return 3;
  }
  else if (std::is_same_v<T, int4>          || std::is_same_v<T, uint4>          ||
           std::is_same_v<T, double4_16a>   || std::is_same_v<T, float4>         ||
           std::is_same_v<T, char4>         || std::is_same_v<T, uchar4>         ||
           std::is_same_v<T, longlong4_16a> || std::is_same_v<T, ulonglong4_16a> ||
           std::is_same_v<T, short4>        || std::is_same_v<T, ushort4>) {
    return 4;
  }
  else {
    rtErr("Unknown data type " + std::string(std::type_index(typeid(T)).name()) + " encountered.",
          "getHpcVectorTypeSize");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::string getHpcVectorTypeName() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (std::is_same_v<T, int2>) return "int2";
  else if (std::is_same_v<T, int3>) return "int3";
  else if (std::is_same_v<T, int4>) return "int4";
  else if (std::is_same_v<T, uint2>) return "unsigned_int2";
  else if (std::is_same_v<T, uint3>) return "unsigned_int3";
  else if (std::is_same_v<T, uint4>) return "unsigned_int4";
  else if (std::is_same_v<T, double2>) return "double2";
  else if (std::is_same_v<T, double3>) return "double3";
  else if (std::is_same_v<T, double4_16a>) return "double4_16a";
  else if (std::is_same_v<T, float2>) return "float2";
  else if (std::is_same_v<T, float3>) return "float3";
  else if (std::is_same_v<T, float4>) return "float4";
  else if (std::is_same_v<T, char2>) return "char2";
  else if (std::is_same_v<T, char3>) return "char3";
  else if (std::is_same_v<T, char4>) return "char4";
  else if (std::is_same_v<T, uchar2>) return "unsigned_char2";
  else if (std::is_same_v<T, uchar3>) return "unsigned_char3";
  else if (std::is_same_v<T, uchar4>) return "unsigned_char4";
  else if (std::is_same_v<T, longlong2>) return "long_long_int2";
  else if (std::is_same_v<T, longlong3>) return "long_long_int3";
  else if (std::is_same_v<T, longlong4_16a>) return "long_long_int4";
  else if (std::is_same_v<T, ulonglong2>) return "unsigned_long_long_int2";
  else if (std::is_same_v<T, ulonglong3>) return "unsigned_long_long_int3";
  else if (std::is_same_v<T, ulonglong4_16a>) return "unsigned_long_long_int4";
  else if (std::is_same_v<T, short2>) return "short_int2";
  else if (std::is_same_v<T, short3>) return "short_int3";
  else if (std::is_same_v<T, short4>) return "short_int4";
  else if (std::is_same_v<T, ushort2>) return "unsigned_short_int2";
  else if (std::is_same_v<T, ushort3>) return "unsigned_short_int3";
  else if (std::is_same_v<T, ushort4>) return "unsigned_short_int4";
  else if (std::is_same_v<T, int95_t>) return "int95_t";
  else {
    rtErr("Data type " + std::string(std::type_index(typeid(T)).name()) + " is not a recognized "
          "HPC vector type.", "getHpcVectorTypeName");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> Vec2<T>::Vec2(const T x_in, const T y_in) :
    x{x_in}, y{y_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Tinput> Vec2<T>::Vec2(const Tinput v_in) :
    x{v_in.x}, y{v_in.y}
{}

//-------------------------------------------------------------------------------------------------
template <typename T> Vec3<T>::Vec3(const T x_in, const T y_in, const T z_in) :
    x{x_in}, y{y_in}, z{z_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Tinput> Vec3<T>::Vec3(const Tinput v_in) :
    x{v_in.x}, y{v_in.y}, z{v_in.z}
{}

//-------------------------------------------------------------------------------------------------
template <typename T> Vec4<T>::Vec4(const T x_in, const T y_in, const T z_in, const T w_in) :
    x{x_in}, y{y_in}, z{z_in}, w{w_in}
{}
  
//-------------------------------------------------------------------------------------------------
template <typename T> template <typename Tinput> Vec4<T>::Vec4(const Tinput v_in) :
    x{v_in.x}, y{v_in.y}, z{v_in.z}, w{v_in.w}
{}

//-------------------------------------------------------------------------------------------------
template <typename T> double2 vtConv2(const T rhs) {
  return { rhs.x, rhs.y };
}

//-------------------------------------------------------------------------------------------------
template <typename T> double3 vtConv3(const T rhs) {
  return { rhs.x, rhs.y, rhs.z };
}

//-------------------------------------------------------------------------------------------------
template <typename T> double4_16a vtConv4(const T rhs) {
  return { rhs.x, rhs.y, rhs.z, rhs.w };
}

//-------------------------------------------------------------------------------------------------
template <typename T> float2 vtConv2f(const T rhs) {
  return { static_cast<float>(rhs.x), static_cast<float>(rhs.y) };
}

//-------------------------------------------------------------------------------------------------
template <typename T> float3 vtConv3f(const T rhs) {
  return { static_cast<float>(rhs.x), static_cast<float>(rhs.y), static_cast<float>(rhs.z) };
}

//-------------------------------------------------------------------------------------------------
template <typename T> float4 vtConv4f(const T rhs) {
  return { static_cast<float>(rhs.x), static_cast<float>(rhs.y), static_cast<float>(rhs.z),
           static_cast<float>(rhs.w) };
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<double2> vtConv2(const std::vector<T> &rhs) {
  const size_t nelem = rhs.size();
  std::vector<double2> result(nelem);
  for (size_t i = 0; i < nelem; i++) {
    result[i] = { rhs[i].x, rhs[i].y };
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<double3> vtConv3(const std::vector<T> &rhs) {
  const size_t nelem = rhs.size();
  std::vector<double3> result(nelem);
  for (size_t i = 0; i < nelem; i++) {
    result[i] = { rhs[i].x, rhs[i].y, rhs[i].z };
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<double4_16a> vtConv4(const std::vector<T> &rhs) {
  const size_t nelem = rhs.size();
  std::vector<double4_16a> result(nelem);
  for (size_t i = 0; i < nelem; i++) {
    result[i] = { rhs[i].x, rhs[i].y, rhs[i].z, rhs[i].w };
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<float2> vtConv2f(const std::vector<T> &rhs) {
  const size_t nelem = rhs.size();
  std::vector<float2> result(nelem);
  for (size_t i = 0; i < nelem; i++) {
    result[i] = { static_cast<float>(rhs[i].x), static_cast<float>(rhs[i].y) };
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<float3> vtConv3f(const std::vector<T> &rhs) {
  const size_t nelem = rhs.size();
  std::vector<float3> result(nelem);
  for (size_t i = 0; i < nelem; i++) {
    result[i] = { static_cast<float>(rhs[i].x), static_cast<float>(rhs[i].y),
                  static_cast<float>(rhs[i].z) };
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::vector<float4> vtConv4f(const std::vector<T> &rhs) {
  const size_t nelem = rhs.size();
  std::vector<float4> result(nelem);
  for (size_t i = 0; i < nelem; i++) {
    result[i] = { static_cast<float>(rhs[i].x), static_cast<float>(rhs[i].y),
                  static_cast<float>(rhs[i].z), static_cast<float>(rhs[i].w) };
  }
  return result;
}

} // namespace data_types
} // namespace stormm
