// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace diskutil {

//-------------------------------------------------------------------------------------------------
template <typename T>
BinaryFileContent::BinaryFileContent(const T content_in, const bool repeating_in,
                                     const std::string description_in,
                                     const HybridTargetLevel tier_in) :
    ct_vptr{std::type_index(typeid(T)).hash_code()},
    type_code{codifyTypeIndex(ct_vptr, false, repeating_in)}, 
    length{1}, description{description_in}, vptr{nullptr}, tier{tier_in}, lower_bound{0},
    upper_bound{0}
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BinaryFileKey::BinaryFileKey(const T* content, size_t length, const std::string &description_in,
                             const HybridTargetLevel tier_in) :
    description{std::string(default_binary_file_descriptor)}, tiers{}, singleton_items{},
    array_items{}, toc_bytes{0}
{
  const size_t tc = std::type_index(typeid(T)).hash_code();
  addArray(reinterpret_cast<const void*>(content), length, tc, description_in, tier_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
BinaryFileKey::BinaryFileKey(const std::vector<T> *content, const std::string &description_in) :
    BinaryFileKey(content->data(), content->size(), description_in, HybridTargetLevel::HOST)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BinaryFileKey::BinaryFileKey(const std::vector<T> &content, const std::string &description_in) :
    BinaryFileKey(content.data(), content.size(), description_in, HybridTargetLevel::HOST)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BinaryFileKey::BinaryFileKey(const Hybrid<T> *content, const std::string &description_in,
                             const HybridTargetLevel tier_in) :
    BinaryFileKey(content->data(tier_in), content->size(), description_in, tier_in)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
BinaryFileKey::BinaryFileKey(const Hybrid<T> &content, const std::string &description_in,
                             const HybridTargetLevel tier_in) :
    BinaryFileKey(content.data(tier_in), content.size(), description_in, tier_in)
{}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BinaryFileKey::addArray(const T* content, const size_t length, const bool repeating_in,
                             const std::string &description_in, const HybridTargetLevel tier_in) {
  addArray(reinterpret_cast<const void*>(content), length, std::type_index(typeid(T)).hash_code(),
           description_in, tier_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BinaryFileKey::addArray(const std::vector<T> *content, const bool repeating_in,
                             const std::string &description_in) {
  const size_t tc = std::type_index(typeid(T)).hash_code();
  addArray(content->data(), content->size(), tc, repeating_in, HybridTargetLevel::HOST,
           description_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BinaryFileKey::addArray(const std::vector<T> &content, const bool repeating_in,
                             const std::string &description_in) {
  const size_t tc = std::type_index(typeid(T)).hash_code();
  addArray(content.data(), content.size(), tc, repeating_in, HybridTargetLevel::HOST,
           description_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BinaryFileKey::addArray(const Hybrid<T> *content, const bool repeating_in,
                             const std::string &description_in, const HybridTargetLevel tier_in) {
  const size_t tc = std::type_index(typeid(T)).hash_code();
  addArray(content->data(tier_in), content->size(), tc, repeating_in, tier_in, description_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BinaryFileKey::addArray(const Hybrid<T> &content, const bool repeating_in,
                             const std::string &description_in, const HybridTargetLevel tier_in) {
  const size_t tc = std::type_index(typeid(T)).hash_code();
  addArray(content.data(tier_in), content.size(), tc, repeating_in, tier_in, description_in);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void BinaryFileKey::addScalar(const T parm, const std::string &description_in) {
  
}

} // namespace disktuil
} // namespace stormm
