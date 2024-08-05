#include <vector>
#include <algorithm>
#include <string.h>
#include "copyright.h"
#include "Constants/scaling.h"
#include "Constants/symbol_values.h"
#include "hybrid.h"
#include "DataTypes/stormm_vector_types.h"
#include "DataTypes/common_types.h"
#include "Reporting/error_format.h"
#include "hybrid.h"

stormm::card::Ledger gbl_mem_balance_sheet;

namespace stormm {
namespace card {

using errors::rtErr;
using errors::terminalFormat;

//-------------------------------------------------------------------------------------------------
Ledger::Ledger() :
  total_expedited{0},
  total_decoupled{0},
  total_unified{0},
  total_host_mounted{0},
  total_devc_only{0},
  entries{},
  free_slots{}
{}

//-------------------------------------------------------------------------------------------------
int Ledger::getActiveEntryCount() const {
  return entries.size() - free_slots.size();
}

//-------------------------------------------------------------------------------------------------
int Ledger::getSerialNumber() {
  int serno;
  const int n_free_slots = free_slots.size();
  if (n_free_slots > 0) {
    serno = free_slots[0];
    free_slots[0] = free_slots[n_free_slots - 1];
    free_slots.resize(n_free_slots - 1);
  }
  else {
    serno = entries.size();
  }
  return serno;
}

//-------------------------------------------------------------------------------------------------
llint Ledger::getTotalExpedited() const {
  return total_expedited;
}

//-------------------------------------------------------------------------------------------------
llint Ledger::getTotalDecoupled() const {
  return total_decoupled;
}

//-------------------------------------------------------------------------------------------------
llint Ledger::getTotalUnified() const {
  return total_unified;
}

//-------------------------------------------------------------------------------------------------
llint Ledger::getTotalHostOnly() const {
  return total_host_only;
}

//-------------------------------------------------------------------------------------------------
llint Ledger::getTotalDevcOnly() const {
  return total_devc_only;
}

//-------------------------------------------------------------------------------------------------
llint Ledger::getTotalHostMounted() const {
  return total_host_mounted;
}

//-------------------------------------------------------------------------------------------------
void Ledger::setEntry(const int serno, const HybridKind kind, const HybridLabel hlbl,
                      const size_t length, const size_t element_size, const int allocations) {  
  LedgerEntry lgt;
  lgt.label = hlbl;
  lgt.length = length;
  lgt.element_size = element_size;
  lgt.kind = kind;
  lgt.allocations = allocations;
  lgt.active = true;
  if (static_cast<size_t>(serno) == entries.size()) {
    entries.push_back(lgt);
  }
  else if (static_cast<size_t>(serno) < entries.size()) {
    entries[serno] = lgt;
  }
  else {
    rtErr("A Hybrid object with serial number " + std::to_string(serno) + " cannot be part of a "
          "list with " + std::to_string(entries.size()) + " items.", "Ledger", "setEntry");
  }
}

//-------------------------------------------------------------------------------------------------
void Ledger::unsetEntry(const int serno) {
  if (static_cast<size_t>(serno) < entries.size()) {
    entries[serno].active = false;
    free_slots.push_back(serno);
  }
  else {
    rtErr("A Hybrid object with serial number " + std::to_string(serno) + " cannot be part of a "
          "list with " + std::to_string(entries.size()) + " items.", "Ledger", "unsetEntry");
  }
}

//-------------------------------------------------------------------------------------------------
LedgerEntry Ledger::getEntry(const int index) const {
  return entries[index];
}

//-------------------------------------------------------------------------------------------------
std::vector<LedgerEntry> Ledger::getEntry(const std::string &name) const {
  std::vector<LedgerEntry> named_entries;
  const size_t n_entries = entries.size();
  for (size_t i = 0; i < n_entries; i++) {
    if (strcmp(entries[i].label.name, name.c_str()) == 0) {
      named_entries.push_back(entries[i]);
    }
  }
  return named_entries;
}

//-------------------------------------------------------------------------------------------------
void Ledger::logMemory(const size_t capacity, const size_t element_size, const HybridFormat fmt,
                       const HybridLabel hlbl, const llint multiplier) {

  // Calculate the contribution, before knowing where to put it
  const llint contribution = static_cast<llint>(capacity) * static_cast<llint>(element_size) *
                             multiplier;

  // Using the serial number stored by the Hybrid object, adjust the recorded memory storage
  const int serno = hlbl.serial_number;
  if (strcmp(entries[serno].label.name, hlbl.name) != 0) {
    rtErr("Mismatch logging " + std::string(hlbl.name) + " into an entry for " +
          std::string(entries[serno].label.name) + ".", "Ledger", "logMemory");
  }
  if (entries[serno].element_size != element_size) {
    rtErr("Mismatch in element sizes for " + std::string(hlbl.name) + ": " +
          std::to_string(entries[serno].element_size) + " expected versus " +
          std::to_string(element_size) + " found.", "Ledger", "logMemory");
  }
  entries[serno].length = capacity;
  if (multiplier > 0LL) {
    entries[serno].allocations += 1;
  }

  // Log contribution in the pooled categories
  switch (fmt) {
#ifdef STORMM_USE_HPC
  case HybridFormat::EXPEDITED:
    total_expedited += contribution;
    break;
  case HybridFormat::DECOUPLED:
    total_decoupled += contribution;
    break;
  case HybridFormat::UNIFIED:
    total_unified += contribution;
    break;
  case HybridFormat::HOST_ONLY:
    total_host_only += contribution;
    break;
  case HybridFormat::DEVICE_ONLY:
    total_devc_only += contribution;
    break;
  case HybridFormat::HOST_MOUNTED:
    total_host_mounted += contribution;
    break;
#else
  case HybridFormat::HOST_ONLY:
    total_host_only += contribution;
    break;
#endif
  }
}

//-------------------------------------------------------------------------------------------------
bool Ledger::testActive(const int serno) {
  if (static_cast<size_t>(serno) >= entries.size()) {
    return false;
  }
  return entries[serno].active;
}

//-------------------------------------------------------------------------------------------------
int Ledger::getAllocations(int serno) {
  if (static_cast<size_t>(serno) >= entries.size()) {
    rtErr("a Hybrid object with serial number " + std::to_string(serno) + " cannot be part of a "
          "list with " + std::to_string(entries.size()) + " items.", "Ledger", "getAllocations");
  }
  if (entries[serno].active == false) {
    rtErr("no Hybrid object with serial number " + std::to_string(serno) + " is active.  Last "
          "known object name: " + std::string(entries[serno].label.name) + ".", "Ledger",
          "getAllocations");
  }
  return entries[serno].allocations;
}

//-------------------------------------------------------------------------------------------------
void Ledger::printMemoryProfile(const int n_display, const llint display_threshold) {

  // Print the summary across all formats
  printf("--------   (%c) Expedited arrays:       %12llu bytes\n", expedited_code,
         total_expedited);
  printf(" Hybrid    (%c) Decoupled arrays:       %12llu bytes\n", decoupled_code,
         total_decoupled);
  printf(" Ledger    (%c) Unified virtual memory: %12llu bytes\n", unified_code, total_unified);
  printf(" Memory    (%c) Host exclusive:         %12llu bytes\n", host_only_code,
         total_host_only);
  printf(" Memory    (%c) Host page-locked:       %12llu bytes\n", host_mounted_code,
         total_host_mounted);
  printf("--------   (%c) Device exclusive:       %12llu bytes\n\n", devc_only_code,
         total_devc_only);
  if (total_expedited == 0ll && total_decoupled == 0ll && total_unified == 0ll &&
      total_host_mounted == 0ll && total_devc_only == 0ll) {
    return;
  }
  printf("Significant allocations:\n");
  printf("Descriptor Tag       Format Kind  Serial Num.  Size (bytes)  Allocations\n");
  printf("---------------------- ---- ----  -----------  ------------  -----------\n");

  // Sort all arrays in decreasing order of size, looking at only the active entries
  const int n_entries = entries.size();
  const int n_active = n_entries - free_slots.size();
  std::vector<LedgerEntry> active_entries;
  active_entries.resize(n_active);
  std::vector<longlong2> memlist;
  memlist.resize(n_active);
  int j = 0;
  for (int i = 0; i < n_entries; i++) {
    if (entries[i].active) {
      active_entries[j] = entries[i];
      memlist[j].x = static_cast<llint>(entries[i].length) *
                     static_cast<llint>(entries[i].element_size);
      memlist[j].y = static_cast<llint>(j);
      j++;
    }
  }
  std::sort(memlist.begin(), memlist.end(), [](longlong2 a, longlong2 b) { return a.x > b.x; }); 

  // Sweep through the list, printing notable entries
  int n_printed = 0;
  bool ptr_printed = false;
  for (int i = 0; i < n_active; i++) {
    if (n_display < 0 || n_printed < n_display || memlist[i].x >= display_threshold) {
      const int idx = memlist[i].y;
      std::string kind;
      char astr_char;
      switch (active_entries[idx].kind) {
      case HybridKind::POINTER:
        kind = "Ptr.";
        astr_char = '*';
        ptr_printed = true;
        break;
      case HybridKind::ARRAY:
        kind = "Ary.";
        astr_char = ' ';
        break;
      }
      printf("%-22.22s    %c %4.4s  %11d  %12lld%c %11d%c\n", active_entries[idx].label.name,
             active_entries[idx].label.format, kind.c_str(),
             active_entries[idx].label.serial_number, memlist[i].x, astr_char,
             active_entries[idx].allocations, astr_char);
      active_entries[idx].active = false;
      n_printed++;
    }
  }

  // Print a summary of all remaining hybrid memory allocations
  if (n_printed < n_active) {
    llint remainder = 0;
    int agg_alloc = 0;
    std::string agg_types = "";
    for (int i = 0; i < n_active; i++) {
      if (active_entries[i].active == false || active_entries[i].kind == HybridKind::POINTER) {
        continue;
      }
      remainder += static_cast<llint>(active_entries[i].element_size) *
                   static_cast<llint>(active_entries[i].length);
      agg_alloc += active_entries[i].allocations;
      if (agg_types.find(active_entries[i].label.format) > agg_types.size()) {
        agg_types += active_entries[i].label.format;
      }
    }
    std::string new_arrays = std::to_string(n_active - n_printed) + " more arrays";
    printf("%-22.22s %4.4s                    %12lld  %11d\n", new_arrays.c_str(),
           agg_types.c_str(), remainder, agg_alloc);
  }

  // Print the explanation about POINTER-kind Hybrid objects
  if (ptr_printed) {
    std::string ptr_explanation = terminalFormat("[*] Sizes and allocations marked with an "
                                                 "asterisk indicate that POINTER-kind Hybrid "
                                                 "objects have implicit sizes based either on "
                                                 "specified bounds or on the extent of their "
                                                 "ARRAY-kind targets, and report allocations "
                                                 "based on their targets.", "", "", 0, 0, 4);
    printf("%s\n", ptr_explanation.c_str());
  }
}

} // namespace card
} // namespace stormm
