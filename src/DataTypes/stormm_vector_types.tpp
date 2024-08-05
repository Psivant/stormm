// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace data_types {

//-------------------------------------------------------------------------------------------------
template <typename T> bool isHpcVectorType() {
  const size_t ct = std::type_index(typeid(T)).hash_code();

  // The int95_t type is not counted among the standard CUDA or HIP HPC vector tuple types.
  return (ct == int2_type_index || ct == int3_type_index || ct == int4_type_index ||
          ct == double2_type_index || ct == double3_type_index || ct == double4_type_index ||
          ct == float2_type_index || ct == float3_type_index || ct == float4_type_index ||
          ct == char2_type_index || ct == char3_type_index || ct == char4_type_index ||
          ct == uchar2_type_index || ct == uchar3_type_index || ct == uchar4_type_index ||
          ct == uint2_type_index || ct == uint3_type_index || ct == uint4_type_index ||
          ct == longlong2_type_index || ct == longlong3_type_index || ct == longlong4_type_index ||
          ct == ulonglong2_type_index || ct == ulonglong3_type_index ||
          ct == ulonglong4_type_index || ct == short2_type_index || ct == short3_type_index ||
          ct == short4_type_index || ct == ushort2_type_index || ct == ushort3_type_index ||
          ct == ushort4_type_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool isSignedIntegralHpcVectorType() {
  const size_t ct = std::type_index(typeid(T)).hash_code();

  // The int95_t type is not counted among the signed integer vector types because it is not
  // efficient to make Hybrid objects out of such things, and it is not among the standard CUDA
  // or HIP HPC vector tuple types.
  return (ct == int2_type_index || ct == int3_type_index || ct == int4_type_index ||
          ct == char2_type_index || ct == char3_type_index || ct == char4_type_index ||
          ct == longlong2_type_index || ct == longlong3_type_index || ct == longlong4_type_index ||
          ct == short2_type_index || ct == short3_type_index || ct == short4_type_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool isUnsignedIntegralHpcVectorType() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  return (ct == uint2_type_index || ct == uint3_type_index || ct == uint4_type_index ||
          ct == uchar2_type_index || ct == uchar3_type_index || ct == uchar4_type_index ||
          ct == ulonglong2_type_index || ct == ulonglong3_type_index ||
          ct == ulonglong4_type_index || ct == ushort2_type_index || ct == ushort3_type_index ||
          ct == ushort4_type_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T> bool isFloatingPointHpcVectorType() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  return (ct == double2_type_index || ct == double3_type_index || ct == double4_type_index ||
          ct == float2_type_index || ct == float3_type_index || ct == float4_type_index);
}

//-------------------------------------------------------------------------------------------------
template <typename T> int getHpcVectorTypeSize() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == int2_type_index || ct == uint2_type_index || ct == double2_type_index ||
      ct == float2_type_index || ct == char2_type_index || ct == uchar2_type_index ||
      ct == longlong2_type_index || ct == ulonglong2_type_index || ct == short2_type_index ||
      ct == ushort2_type_index) {
    return 2;
  }
  else if (ct == int3_type_index || ct == uint3_type_index || ct == double3_type_index ||
           ct == float3_type_index || ct == char3_type_index || ct == uchar3_type_index ||
           ct == longlong3_type_index || ct == ulonglong3_type_index || ct == short3_type_index ||
           ct == ushort3_type_index) {
    return 3;
  }
  else if (ct == int4_type_index || ct == uint4_type_index || ct == double4_type_index ||
           ct == float4_type_index || ct == char4_type_index || ct == uchar4_type_index ||
           ct == longlong4_type_index || ct == ulonglong4_type_index || ct == short4_type_index ||
           ct == ushort4_type_index) {
    return 4;
  }
  else {
    rtErr("Unknown data type " + std::string(std::type_index(typeid(T)).name()) + " encountered.",
          "getHpcVectorTypeSize");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
template <typename T> std::string getStormmHpcVectorTypeName() {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct == int2_type_index) return "int2";
  else if (ct == int3_type_index) return "int3";
  else if (ct == int4_type_index) return "int4";
  else if (ct == uint2_type_index) return "unsigned_int2";
  else if (ct == uint3_type_index) return "unsigned_int3";
  else if (ct == uint4_type_index) return "unsigned_int4";
  else if (ct == double2_type_index) return "double2";
  else if (ct == double3_type_index) return "double3";
  else if (ct == double4_type_index) return "double4";
  else if (ct == float2_type_index) return "float2";
  else if (ct == float3_type_index) return "float3";
  else if (ct == float4_type_index) return "float4";
  else if (ct == char2_type_index) return "char2";
  else if (ct == char3_type_index) return "char3";
  else if (ct == char4_type_index) return "char4";
  else if (ct == uchar2_type_index) return "unsigned_char2";
  else if (ct == uchar3_type_index) return "unsigned_char3";
  else if (ct == uchar4_type_index) return "unsigned_char4";
  else if (ct == longlong2_type_index) return "long_long_int2";
  else if (ct == longlong3_type_index) return "long_long_int3";
  else if (ct == longlong4_type_index) return "long_long_int4";
  else if (ct == ulonglong2_type_index) return "unsigned_long_long_int2";
  else if (ct == ulonglong3_type_index) return "unsigned_long_long_int3";
  else if (ct == ulonglong4_type_index) return "unsigned_long_long_int4";
  else if (ct == short2_type_index) return "short_int2";
  else if (ct == short3_type_index) return "short_int3";
  else if (ct == short4_type_index) return "short_int4";
  else if (ct == ushort2_type_index) return "unsigned_short_int2";
  else if (ct == ushort3_type_index) return "unsigned_short_int3";
  else if (ct == ushort4_type_index) return "unsigned_short_int4";
  else if (ct == int95t_type_index) return "int95_t";
  else {
    rtErr("Data type " + std::string(std::type_index(typeid(T)).name()) + " is not a recognized "
          "HPC vector type.", "getStormmHpcVectorTypeName");
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
template <typename T> double4 vtConv4(const T rhs) {
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
template <typename T> std::vector<double4> vtConv4(const std::vector<T> &rhs) {
  const size_t nelem = rhs.size();
  std::vector<double4> result(nelem);
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
