#include "copyright.h"
#include "stormm_vector_types.h"

namespace stormm {
namespace data_types {

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
