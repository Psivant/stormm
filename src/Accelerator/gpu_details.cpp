// -*-c++-*-
#include <cstdio>
#include <cstdlib>
#include <vector>
#include "copyright.h"
#include "Reporting/error_format.h"
#include "gpu_details.h"

stormm::card::GpuDetails null_gpu;

namespace stormm {
namespace card {

//-------------------------------------------------------------------------------------------------
GpuDetails::GpuDetails() :
    available{false}, supported{false}, arch_major{0}, arch_minor{0}, smp_count{1}, card_ram{0},
    max_threads_per_block{1}, max_threads_per_smp{1}, max_blocks_per_smp{1},
    max_shared_per_block{0}, max_shared_per_smp{0}, global_cache_size{0}, registers_per_block{0},
    registers_per_smp{0},
    card_name{std::string("blank_gpu")}
{}

//-------------------------------------------------------------------------------------------------
bool GpuDetails::getAvailability() const {
  return available;
}

//-------------------------------------------------------------------------------------------------
bool GpuDetails::getGpuSupported() const {
  return supported;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getArchMajor() const {
  return arch_major;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getArchMinor() const {
  return arch_minor;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getSMPCount() const {
  return smp_count;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getCardRam() const {
  return card_ram;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getMaxThreadsPerBlock() const {
  return max_threads_per_block;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getMaxThreadsPerSMP() const {
  return max_threads_per_smp;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getMaxBlocksPerSMP() const {
  return max_blocks_per_smp;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getMaxSharedPerBlock() const {
  return max_shared_per_block;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getMaxSharedPerSMP() const {
  return max_shared_per_smp;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getGlobalCacheSize() const {
  return global_cache_size;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getRegistersPerBlock() const {
  return registers_per_block;
}

//-------------------------------------------------------------------------------------------------
int GpuDetails::getRegistersPerSMP() const {
  return registers_per_smp;
}

//-------------------------------------------------------------------------------------------------
std::string GpuDetails::getCardName() const {
  return card_name;
}

//-------------------------------------------------------------------------------------------------
void GpuDetails::setSMPCount(const int smp_count_in) {
  smp_count = smp_count_in;
}

//-------------------------------------------------------------------------------------------------
bool GpuDetails::operator==(const GpuDetails &right) const {
  return (available             == right.available &&
          supported             == right.supported &&
          arch_major            == right.arch_major &&
          arch_minor            == right.arch_minor &&
          smp_count             == right.smp_count &&
          card_ram              == right.card_ram &&
          card_name             == right.card_name &&
          max_threads_per_block == right.max_threads_per_block &&
          max_threads_per_smp   == right.max_threads_per_smp &&
          max_blocks_per_smp    == right.max_blocks_per_smp &&
          max_shared_per_block  == right.max_shared_per_block &&
          max_shared_per_smp    == right.max_shared_per_smp &&
          global_cache_size     == right.global_cache_size &&
          registers_per_block   == right.registers_per_block &&
          registers_per_smp     == right.registers_per_smp);
}

//-------------------------------------------------------------------------------------------------
bool GpuDetails::operator!=(const GpuDetails &right) const {
  return (! operator==(right));
}
  
} // namespace card
} // namespace stormm
