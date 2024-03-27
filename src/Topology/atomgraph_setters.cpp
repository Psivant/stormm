#include <cmath>
#include <cstdio>
#include <climits>
#include "copyright.h"
#include "Constants/generalized_born.h"
#include "FileManagement/file_listing.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "UnitTesting/approx.h"
#include "UnitTesting/unit_test_enumerators.h"
#include "atomgraph.h"

namespace stormm {
namespace topology {

using diskutil::getBaseName;
using parse::CaseSensitivity;
using parse::char4ToString;
using parse::realToString;
using parse::stringToChar4;
using parse::strncmpCased;
using data_types::operator==;
using testing::Approx;
using testing::ComparisonType;
using namespace generalized_born_defaults;

//-------------------------------------------------------------------------------------------------
void AtomGraph::modifyAtomMobility(const MobilitySetting movement) {

  // Just feed into the next, most general case.  The cost of the unneccessary bounds check is
  // trivial compared to the simplicity this affords in the code.
  modifyAtomMobility(0, atom_count, movement);
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::modifyAtomMobility(const int low_index, const int high_index,
                                   const MobilitySetting movement) {

  // Range check as this will use the pointer
  if (low_index < 0 || high_index > atom_count || high_index <= low_index) {
    rtErr("A topology with " + std::to_string(atom_count) + " atoms cannot change atom mobility "
          "for indices " + std::to_string(low_index) + " to " + std::to_string(high_index) + ".",
          "AtomGraph", "modifyAtomMobility");
  }
  const int int_bits = sizeof(int) * 8;
  int* m_ptr = mobile_atoms.data();
  for (int i = low_index; i < high_index; i++) {
    const int access_index = i / int_bits;
    const int bshift = i - (access_index * int_bits);
    const uint orig_mask = static_cast<uint>(m_ptr[access_index]);
    switch (movement) {
    case MobilitySetting::OFF:
      m_ptr[access_index] = (orig_mask & (~(0x1 << bshift)));
      break;
    case MobilitySetting::ON:
      m_ptr[access_index] = (orig_mask | (0x1 << bshift));
      break;
    case MobilitySetting::TOGGLE:
      m_ptr[access_index] = (orig_mask ^ (0x1 << bshift));
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::modifyAtomMobility(const int index, const MobilitySetting movement) {
  const int int_bits = sizeof(int) * 8;
  const int access_index = index / int_bits;
  const uint m_val = mobile_atoms.readHost(access_index);
  const int bshift = index - (access_index * int_bits);
  switch (movement) {
  case MobilitySetting::OFF:
    mobile_atoms.putHost(m_val & (~(0x1 << bshift)), access_index);
    break;
  case MobilitySetting::ON:
    mobile_atoms.putHost(m_val | (0x1 << bshift), access_index);
    break;
  case MobilitySetting::TOGGLE:
    mobile_atoms.putHost(m_val ^ (0x1 << bshift), access_index);
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::modifyAtomMobility(const std::vector<int> &mask, const MobilitySetting movement) {
  const size_t n_atoms = mask.size();
  for (size_t i = 0LLU; i < n_atoms; i++) {
    modifyAtomMobility(mask[i], movement);
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::setSource(const std::string &new_source) {
  source = new_source;
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::setBondParameters(const double new_keq, const double new_leq, const bool set_keq,
                                  const bool set_leq, const int parm_idx,
                                  const ExceptionResponse policy) {
  if (set_keq == false && set_leq == false) {
    return;
  }
  if (parm_idx < 0 || parm_idx >= bond_parameter_count) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A bond with parameter index " + std::to_string(parm_idx) + " does not exist in a "
            "topology with " + std::to_string(bond_parameter_count) + " unique bond types.",
            "AtomGraph", "setBondParameters");
    case ExceptionResponse::WARN:
      {
        std::string action_str("");
        if (set_keq && set_leq) {
          action_str += "stiffness to " + realToString(new_keq) + " and the equilibrium to " +
                        realToString(new_leq);
        }
        else if (set_keq) {
          action_str += "stiffness to " + realToString(new_keq);
        }
        else if (set_leq) {
          action_str += "equilibrium to " + realToString(new_leq);
        }
        rtWarn("A bond with parameter index " + std::to_string(parm_idx) + " does not exist in a "
               "topology with " + std::to_string(bond_parameter_count) + " unique bond types.  No "
               "action will be taken to change the " + action_str + ".", "AtomGraph",
               "setBondParameters");
      }
      break;
    case ExceptionResponse::SILENT:
      break;      
    }
  }
  else {
    if (set_keq) {
      bond_stiffnesses.putHost(new_keq, parm_idx);
      sp_bond_stiffnesses.putHost(new_keq, parm_idx);
    }
    if (set_leq) {
      bond_equilibria.putHost(new_leq, parm_idx);
      sp_bond_equilibria.putHost(new_leq, parm_idx);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::setBondParameters(const double new_keq, const double new_leq, const bool set_keq,
                                  const bool set_leq, const char4 type_a, const char4 type_b,
                                  const ExceptionResponse policy) {
  const char4* atom_type_ptr = atom_types.data();
  const int* bond_i_atom_ptr = bond_i_atoms.data();
  const int* bond_j_atom_ptr = bond_j_atoms.data();
  const int* atom_type_idx_ptr = lennard_jones_indices.data();
  int bond_parm_idx = -1;
  for (int pos = 0; pos < bond_term_count; pos++) {
    const int atom_i = bond_i_atom_ptr[pos];
    const int atom_j = bond_j_atom_ptr[pos];
    const int atyp_i = atom_type_idx_ptr[atom_i];
    const int atyp_j = atom_type_idx_ptr[atom_j];
    if ((atom_type_ptr[atyp_i] == type_a && atom_type_ptr[atyp_j] == type_b) ||
        (atom_type_ptr[atyp_j] == type_a && atom_type_ptr[atyp_i] == type_b)) {
      bond_parm_idx = bond_parameter_indices.readHost(pos);
      break;
    }
  }
  if (bond_parm_idx < 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("No bond between atom types " + char4ToString(type_a) + " and " +
            char4ToString(type_b) + " was found in topology " + source + ".", "AtomGraph",
            "setBondParameters");
    case ExceptionResponse::WARN:
      rtWarn("No bond between atom types " + char4ToString(type_a) + " and " +
             char4ToString(type_b) + " was found in topology " + source + ".  No action will be "
             "taken.", "AtomGraph", "setBondParameters");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  
  // Do not report the spurious parameter index if a warning about the types has already been
  // issued.
  setBondParameters(new_keq, new_leq, set_keq, set_leq, bond_parm_idx,
                    (policy == ExceptionResponse::WARN) ? ExceptionResponse::SILENT : policy); 
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::setAngleParameters(const double new_keq, const double new_teq, const bool set_keq,
                                   const bool set_teq, const int parm_idx,
                                   const ExceptionResponse policy) {
  if (set_keq == false && set_teq == false) {
    return;
  }
  if (parm_idx < 0 || parm_idx >= angl_parameter_count) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A harmonic angle term with parameter index " + std::to_string(parm_idx) + " does not "
            "exist in a topology with " + std::to_string(angl_parameter_count) + " unique angle "
            "types.", "AtomGraph", "setAngleParameters");
    case ExceptionResponse::WARN:
      {
        std::string action_str("");
        if (set_keq && set_teq) {
          action_str += "stiffness to " + realToString(new_keq) + " and the equilibrium to " +
                        realToString(new_teq);
        }
        else if (set_keq) {
          action_str += "stiffness to " + realToString(new_keq);
        }
        else if (set_teq) {
          action_str += "equilibrium to " + realToString(new_teq);
        }
        rtWarn("A harmonic angle with parameter index " + std::to_string(parm_idx) + " does not "
               "exist in a topology with " + std::to_string(angl_parameter_count) + " unique "
               "angle types.  No action will be taken to change the " + action_str + ".",
               "AtomGraph", "setAngleParameters");
      }
      break;
    case ExceptionResponse::SILENT:
      break;      
    }
  }
  else {
    if (set_keq) {
      angl_stiffnesses.putHost(new_keq, parm_idx);
      sp_angl_stiffnesses.putHost(new_keq, parm_idx);
    }
    if (set_teq) {
      angl_equilibria.putHost(new_teq, parm_idx);
      sp_angl_equilibria.putHost(new_teq, parm_idx);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::setAngleParameters(const double new_keq, const double new_teq, const bool set_keq,
                                   const bool set_teq, const char4 type_a, const char4 type_b,
                                   const char4 type_c, const ExceptionResponse policy) {
  const char4* atom_type_ptr = atom_types.data();
  const int* angl_i_atom_ptr = angl_i_atoms.data();
  const int* angl_j_atom_ptr = angl_j_atoms.data();
  const int* angl_k_atom_ptr = angl_k_atoms.data();
  const int* atom_type_idx_ptr = lennard_jones_indices.data();
  int angl_parm_idx = -1;
  for (int pos = 0; pos < angl_term_count; pos++) {
    const int atom_i = angl_i_atom_ptr[pos];
    const int atom_j = angl_j_atom_ptr[pos];
    const int atom_k = angl_k_atom_ptr[pos];
    const int atyp_i = atom_type_idx_ptr[atom_i];
    const int atyp_j = atom_type_idx_ptr[atom_j];
    const int atyp_k = atom_type_idx_ptr[atom_k];
    if (atom_type_ptr[atyp_j] == type_b &&
        ((atom_type_ptr[atyp_i] == type_a && atom_type_ptr[atyp_k] == type_c) ||
         (atom_type_ptr[atyp_k] == type_a && atom_type_ptr[atyp_i] == type_c))) {
      angl_parm_idx = angl_parameter_indices.readHost(pos);
      break;
    }
  }
  if (angl_parm_idx < 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("No angle between atom types " + char4ToString(type_a) + ", " + char4ToString(type_b) +
            " and " + char4ToString(type_c) + " was found in topology " + source + ".",
            "AtomGraph", "setAngleParameters");
    case ExceptionResponse::WARN:
      rtWarn("No angle between atom types " + char4ToString(type_a) + ", " +
             char4ToString(type_b) + ", and " + char4ToString(type_c) + " was found in topology " +
             source + ".  No action will be taken.", "AtomGraph", "setAngleParameters");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }

  // Do not report the spurious parameter index if a warning about the types has already been
  // issued.
  setAngleParameters(new_keq, new_teq, set_keq, set_teq, angl_parm_idx,
                     (policy == ExceptionResponse::WARN) ? ExceptionResponse::SILENT : policy); 
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::setDihedralParameters(const double new_amplitude, const double new_phase_angle,
                                      const bool set_amplitude, const bool set_phase_angle,
                                      const int parm_idx, const ExceptionResponse policy) {
  if (parm_idx < 0 || parm_idx >= dihe_parameter_count) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A cosine-based dihedral term with parameter index " + std::to_string(parm_idx) +
            " does not exist in a topology with " + std::to_string(dihe_parameter_count) +
            " unique dihedral types.", "AtomGraph", "setDihedralParameters");
    case ExceptionResponse::WARN:
      {
        std::string action_str("");
        if (set_amplitude && set_phase_angle) {
          action_str += "amplitude to " + realToString(new_amplitude) + " and the phase angle "
                        "to " + realToString(new_phase_angle);
        }
        else if (set_amplitude) {
          action_str += "amplitude to " + realToString(new_amplitude);
        }
        else if (set_phase_angle) {
          action_str += "phase angle to " + realToString(new_phase_angle);
        }
        rtWarn("A cosine-based dihedral term with parameter index " + std::to_string(parm_idx) +
               " does not exist in a topology with " + std::to_string(dihe_parameter_count) +
               " unique dihedral types.  No action will be taken to change the " + action_str +
               ".", "AtomGraph", "setDihedralParameters");
      }
      break;
    case ExceptionResponse::SILENT:
      break;      
    }
  }
  else {
    if (set_amplitude) {
      dihe_amplitudes.putHost(new_amplitude, parm_idx);
      sp_dihe_amplitudes.putHost(new_amplitude, parm_idx);
    }
    if (set_phase_angle) {
      dihe_phase_angles.putHost(new_phase_angle, parm_idx);
      sp_dihe_phase_angles.putHost(new_phase_angle, parm_idx);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::setDihedralParameters(const double new_amplitude, const double new_phase_angle,
                                      const bool set_amplitude, const bool set_phase_angle,
                                      const char4 type_a, const char4 type_b, const char4 type_c,
                                      const char4 type_d, const double periodicity,
                                      const ExceptionResponse policy) {
  const char4* atom_type_ptr = atom_types.data();
  const int* dihe_i_atom_ptr = dihe_i_atoms.data();
  const int* dihe_j_atom_ptr = dihe_j_atoms.data();
  const int* dihe_k_atom_ptr = dihe_k_atoms.data();
  const int* dihe_l_atom_ptr = dihe_l_atoms.data();
  const double* dihe_per_ptr = dihe_periodicities.data();
  const int* atom_type_idx_ptr = lennard_jones_indices.data();
  int dihe_parm_idx = -1;
  const Approx pdval(periodicity, ComparisonType::ABSOLUTE, 1.0e-4);
  for (int pos = 0; pos < dihe_term_count; pos++) {
    const int atom_i = dihe_i_atom_ptr[pos];
    const int atom_j = dihe_j_atom_ptr[pos];
    const int atom_k = dihe_k_atom_ptr[pos];
    const int atom_l = dihe_l_atom_ptr[pos];
    const int atyp_i = atom_type_idx_ptr[atom_i];
    const int atyp_j = atom_type_idx_ptr[atom_j];
    const int atyp_k = atom_type_idx_ptr[atom_k];
    const int atyp_l = atom_type_idx_ptr[atom_l];
    if (pdval.test(dihe_per_ptr[pos]) && 
        ((atom_type_ptr[atyp_i] == type_a && atom_type_ptr[atyp_j] == type_b &&
          atom_type_ptr[atyp_k] == type_c && atom_type_ptr[atyp_l] == type_d) ||
         (atom_type_ptr[atyp_l] == type_a && atom_type_ptr[atyp_k] == type_b &&
          atom_type_ptr[atyp_j] == type_c && atom_type_ptr[atyp_i] == type_d))) {
      dihe_parm_idx = dihe_parameter_indices.readHost(pos);
      break;
    }
  }
  if (dihe_parm_idx < 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("No dihedral angle between atom types " + char4ToString(type_a) + ", " +
            char4ToString(type_b) + ", " + char4ToString(type_c) + ", and " +
            char4ToString(type_d) + " was found in topology " + source + ".", "AtomGraph",
            "setDihedralParameters");
    case ExceptionResponse::WARN:
      rtWarn("No dihedral angle between atom types " + char4ToString(type_a) + ", " +
             char4ToString(type_b) + ", " + char4ToString(type_c) + ", and " +
             char4ToString(type_d) + " was found in topology " + source + ".  No action will be "
             "taken.", "AtomGraph", "setDihedralParameters");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }

  // Do not report the spurious parameter index if a warning about the types has already been
  // issued.
  setDihedralParameters(new_amplitude, new_phase_angle, set_amplitude, set_phase_angle,
                        dihe_parm_idx,
                        (policy == ExceptionResponse::WARN) ? ExceptionResponse::SILENT : policy); 
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::setUreyBradleyParameters(const double new_keq, const double new_leq,
                                         const bool set_keq, const bool set_leq,
                                         const int parm_idx, const ExceptionResponse policy) {
  if (set_keq == false && set_leq == false) {
    return;
  }
  if (parm_idx < 0 || parm_idx >= urey_bradley_parameter_count) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A Urey-Bradley term with parameter index " + std::to_string(parm_idx) + " does not "
            "exist in a topology with " + std::to_string(urey_bradley_parameter_count) +
            " unique Urey-Bradley types.", "AtomGraph", "setUreyBradleyParameters");
    case ExceptionResponse::WARN:
      {
        std::string action_str("");
        if (set_keq && set_leq) {
          action_str += "stiffness to " + realToString(new_keq) + " and the equilibrium to " +
                        realToString(new_leq);
        }
        else if (set_keq) {
          action_str += "stiffness to " + realToString(new_keq);
        }
        else if (set_leq) {
          action_str += "equilibrium to " + realToString(new_leq);
        }
        rtWarn("A Urey-Bradley term with parameter index " + std::to_string(parm_idx) + " does "
               "not exist in a topology with " + std::to_string(urey_bradley_parameter_count) +
               " unique Urey-Bradley types.  No action will be taken to change the " +
               action_str + ".", "AtomGraph", "setUreyBradleyParameters");
      }
      break;
    case ExceptionResponse::SILENT:
      break;      
    }
  }
  else {
    if (set_keq) {
      urey_bradley_stiffnesses.putHost(new_keq, parm_idx);
      sp_urey_bradley_stiffnesses.putHost(new_keq, parm_idx);
    }
    if (set_leq) {
      urey_bradley_equilibria.putHost(new_leq, parm_idx);
      sp_urey_bradley_equilibria.putHost(new_leq, parm_idx);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::setUreyBradleyParameters(const double new_keq, const double new_leq,
                                         const bool set_keq, const bool set_leq,
                                         const char4 type_a, const char4 type_b,
                                         const char4 type_c, const ExceptionResponse policy) {
  const char4* atom_type_ptr = atom_types.data();
  const int* ubrd_i_atom_ptr = urey_bradley_i_atoms.data();
  const int* ubrd_k_atom_ptr = urey_bradley_k_atoms.data();
  const int* atom_type_idx_ptr = lennard_jones_indices.data();
  int ubrd_parm_idx = -1;
  for (int pos = 0; pos < urey_bradley_term_count; pos++) {
    const int atom_i = ubrd_i_atom_ptr[pos];
    const int atom_k = ubrd_k_atom_ptr[pos];
    const int atyp_i = atom_type_idx_ptr[atom_i];
    const int atyp_k = atom_type_idx_ptr[atom_k];
    if ((atom_type_ptr[atyp_i] == type_a && atom_type_ptr[atyp_k] == type_c) ||
        (atom_type_ptr[atyp_k] == type_a && atom_type_ptr[atyp_i] == type_c)) {

      // Check that the central atom between atoms I and K has the correct type
      const int atom_i_llim = nb12_exclusion_bounds.readHost(atom_i);
      const int atom_i_hlim = nb12_exclusion_bounds.readHost(atom_i + 1);
      bool atyp_j_match = false;
      for (int i = atom_i_llim; i < atom_i_hlim; i++) {
        const int atom_j = nb12_exclusion_list.readHost(i);
        const int atom_j_llim = nb12_exclusion_bounds.readHost(atom_j);
        const int atom_j_hlim = nb12_exclusion_bounds.readHost(atom_j + 1);
        for (int j = atom_j_llim; j < atom_j_hlim; j++) {
          atyp_j_match = (atyp_j_match || (nb12_exclusion_list.readHost(j) == atom_k &&
                                           atom_type_ptr[atom_type_idx_ptr[atom_j]] == type_b));
        }
      }
      if (atyp_j_match) {
        ubrd_parm_idx = urey_bradley_parameter_indices.readHost(pos);
      }
      break;
    }
  }
  if (ubrd_parm_idx < 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("No Urey-Bradley angle between atom types " + char4ToString(type_a) + " and " +
            char4ToString(type_b) + " was found in topology " + source + ".", "AtomGraph",
            "setUreyBradleyParameters");
    case ExceptionResponse::WARN:
      rtWarn("No Urey-Bradley angle between atom types " + char4ToString(type_a) + " and " +
             char4ToString(type_b) + " was found in topology " + source + ".  No action will be "
             "taken.", "AtomGraph", "setUreyBradleyParameters");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }

  // Do not report the spurious parameter index if a warning about the types has already been
  // issued.
  setUreyBradleyParameters(new_keq, new_leq, set_keq, set_leq, ubrd_parm_idx,
                           (policy == ExceptionResponse::WARN) ? ExceptionResponse::SILENT :
                                                                 policy); 
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::setCharmmImprParameters(const double new_stiffness, const double new_phase_angle,
                                        const bool set_stiffness, const bool set_phase_angle,
                                        const int parm_idx, const ExceptionResponse policy) {
  if (parm_idx < 0 || parm_idx >= charmm_impr_parameter_count) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("A CHARMM improper dihedral term with parameter index " + std::to_string(parm_idx) +
            " does not exist in a topology with " + std::to_string(charmm_impr_parameter_count) +
            " unique parameter sets.", "AtomGraph", "setCharmmImproperParameters");
    case ExceptionResponse::WARN:
      {
        std::string action_str("");
        if (set_stiffness && set_phase_angle) {
          action_str += "stiffness to " + realToString(new_stiffness) + " and the phase angle "
                        "to " + realToString(new_phase_angle);
        }
        else if (set_stiffness) {
          action_str += "stiffness to " + realToString(new_stiffness);
        }
        else if (set_phase_angle) {
          action_str += "phase angle to " + realToString(new_phase_angle);
        }
        rtWarn("A CHARMM improper dihedral term with parameter index " + std::to_string(parm_idx) +
               " does not exist in a topology with " + std::to_string(dihe_parameter_count) +
               " unique dihedral types.  No action will be taken to change the " + action_str +
               ".", "AtomGraph", "setCharmmImproperParameters");
      }
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
  else {
    if (set_stiffness) {
      charmm_impr_stiffnesses.putHost(new_stiffness, parm_idx);
      sp_charmm_impr_stiffnesses.putHost(new_stiffness, parm_idx);
    }
    if (set_phase_angle) {
      charmm_impr_phase_angles.putHost(new_phase_angle, parm_idx);
      sp_charmm_impr_phase_angles.putHost(new_phase_angle, parm_idx);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::setCharmmImprParameters(const double new_stiffness, const double new_phase_angle,
                                        const bool set_stiffness, const bool set_phase_angle,
                                        const char4 type_a, const char4 type_b, const char4 type_c,
                                        const char4 type_d, const ExceptionResponse policy) {
  const char4* atom_type_ptr = atom_types.data();
  const int* cimp_i_atom_ptr = charmm_impr_i_atoms.data();
  const int* cimp_j_atom_ptr = charmm_impr_j_atoms.data();
  const int* cimp_k_atom_ptr = charmm_impr_k_atoms.data();
  const int* cimp_l_atom_ptr = charmm_impr_l_atoms.data();
  const int* atom_type_idx_ptr = lennard_jones_indices.data();
  int cimp_parm_idx = -1;
  for (int pos = 0; pos < charmm_impr_term_count; pos++) {
    const int atom_i = cimp_i_atom_ptr[pos];
    const int atom_j = cimp_j_atom_ptr[pos];
    const int atom_k = cimp_k_atom_ptr[pos];
    const int atom_l = cimp_l_atom_ptr[pos];
    const int atyp_i = atom_type_idx_ptr[atom_i];
    const int atyp_j = atom_type_idx_ptr[atom_j];
    const int atyp_k = atom_type_idx_ptr[atom_k];
    const int atyp_l = atom_type_idx_ptr[atom_l];
    if (((atom_type_ptr[atyp_i] == type_a && atom_type_ptr[atyp_j] == type_b &&
          atom_type_ptr[atyp_k] == type_c && atom_type_ptr[atyp_l] == type_d) ||
         (atom_type_ptr[atyp_l] == type_a && atom_type_ptr[atyp_k] == type_b &&
          atom_type_ptr[atyp_j] == type_c && atom_type_ptr[atyp_i] == type_d))) {
      cimp_parm_idx = charmm_impr_parameter_indices.readHost(pos);
      break;
    }
  }
  if (cimp_parm_idx < 0) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("No CHARMM improper dihedral between atom types " + char4ToString(type_a) + ", " +
            char4ToString(type_b) + ", " + char4ToString(type_c) + ", and " +
            char4ToString(type_d) + " was found in topology " + source + ".", "AtomGraph",
            "setCharmmImprParameters");
    case ExceptionResponse::WARN:
      rtWarn("No CHARMM improper dihedral between atom types " + char4ToString(type_a) + ", " +
             char4ToString(type_b) + ", " + char4ToString(type_c) + ", and " +
             char4ToString(type_d) + " was found in topology " + source + ".  No action will be "
             "taken.", "AtomGraph", "setCharmmImprParameters");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }

  // Do not report the spurious parameter index if a warning about the types has already been
  // issued.
  setCharmmImprParameters(new_stiffness, new_phase_angle, set_stiffness, set_phase_angle,
                          cimp_parm_idx,
                          (policy == ExceptionResponse::WARN) ? ExceptionResponse::SILENT :
                                                                policy); 
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::setImplicitSolventModel(const ImplicitSolventModel igb_in,
                                        const double dielectric_in, const double saltcon_in,
                                        const AtomicRadiusSet radii_set,
                                        const ExceptionResponse policy) {
  gb_style = igb_in;

  // Trap GB use cases with Parse radii (these are much older, for Poisson-Boltzmann calculations)
  switch (radii_set) {
  case AtomicRadiusSet::BONDI:
  case AtomicRadiusSet::AMBER6:
  case AtomicRadiusSet::MBONDI:
  case AtomicRadiusSet::MBONDI2:
  case AtomicRadiusSet::MBONDI3:
    break;
  case AtomicRadiusSet::PARSE:

    // The only use of Parse radii is PB calculations, which are currently not supported in STORMM.
    // Nonetheless, the radius set is included for completeness.  Just make sure that no one is
    // using it for Generalized Born calculations.
    switch (gb_style) {
    case ImplicitSolventModel::NONE:
      break;
    case ImplicitSolventModel::HCT_GB:
    case ImplicitSolventModel::OBC_GB:
    case ImplicitSolventModel::OBC_GB_II:
    case ImplicitSolventModel::NECK_GB:
    case ImplicitSolventModel::NECK_GB_II:
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("Parse radii are not intended for Generalized Born calculations and cannot "
              "therefore be applied to topology " + source + ".", "AtomGraph",
              "setImplicitSolventModel");
      case ExceptionResponse::WARN:
        rtWarn("Parse radii are not intended for Generalized Born calculations and should not "
               "therefore be applied to topology " + source + ".", "AtomGraph",
               "setImplicitSolventModel");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
    break;
  case AtomicRadiusSet::NONE:
    switch (gb_style) {
    case ImplicitSolventModel::NONE:
      break;
    case ImplicitSolventModel::HCT_GB:
    case ImplicitSolventModel::OBC_GB:
    case ImplicitSolventModel::OBC_GB_II:
    case ImplicitSolventModel::NECK_GB:
    case ImplicitSolventModel::NECK_GB_II:
      {
        // Radii are not being applied, so check that nonzero radii are at least present
        bool nonzero_radii_found = false;
        const double* radii_ptr = atomic_pb_radii.data();
        for (int i = 0; i < atom_count; i++) {
          nonzero_radii_found = (nonzero_radii_found || radii_ptr[i] > constants::tiny);
        }
        if (nonzero_radii_found == false) {

          // This is probably always an error, but use the policy switch in case someone is
          // trying something clever.
          switch (policy) {
          case ExceptionResponse::DIE:
            rtErr("No nonzero radii were found in topology file " + source + ", or applied "
                  "when setting a Generalized Born implicit solvent model.", "AtomGraph",
                  "setImplicitSolventModel");
          case ExceptionResponse::WARN:
            rtWarn("No nonzero radii were found in topology file " + source + ", or applied "
                   "when setting a Generalized Born implicit solvent model.  This is likely to "
                   "cause problems later.", "AtomGraph", "setImplicitSolventModel");
            break;
          case ExceptionResponse::SILENT:
            break;
          }
        }
      }
    }
    break;
  }

  // Check for offending radii, specifically those that might not fit within the bounds of the
  // "neck" GB models.  If such radii are found, recursively call this function and set the
  // appropriate (Bondi or mBondi3) radii and screening parameters.
  if ((gb_style == ImplicitSolventModel::NECK_GB    && radii_set != AtomicRadiusSet::BONDI) ||
      (gb_style == ImplicitSolventModel::NECK_GB_II && radii_set != AtomicRadiusSet::MBONDI3)) {
    const double* radii_ptr = atomic_pb_radii.data();
    bool bad_radii_found = false;
    for (int i = 0; i < atom_count; i++) {
      bad_radii_found = (bad_radii_found || radii_ptr [i] < 1.0 || radii_ptr[i] > 2.0);
    }
    if (bad_radii_found) {
      std::string correction;
      if (gb_style == ImplicitSolventModel::NECK_GB) {
        correction = "Radii will be replaced by BONDI parameters.";
      }
      else if (gb_style == ImplicitSolventModel::NECK_GB_II) {
        correction = "Radii will be replaced by Modified BONDI-3 parameters.";
      }
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("The \"neck\" GB models are incompatible with atomic radii smaller than 1.0 "
              "Anstroms or larger than 2.0 Angstroms, which were found in topology " +
              getBaseName(source) + ".  To automatically correct this situation (but use "
              "different parameters from what is written in the file) re-run with '-warn' as a "
              "command line argument.", "AtomGraph", "setImplicitSolventModel");
      case ExceptionResponse::WARN:
        rtWarn("The \"neck\" GB models are incompatible with atomic radii smaller than 1.0 "
               "Anstroms or larger than 2.0 Angstroms, which were found in topology " +
               getBaseName(source) + ".  " + correction, "AtomGraph", "setImplicitSolventModel");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
      if (gb_style == ImplicitSolventModel::NECK_GB) {
        setImplicitSolventModel(igb_in, dielectric_in, saltcon_in, AtomicRadiusSet::BONDI, policy);
      }
      else if (gb_style == ImplicitSolventModel::NECK_GB_II) {
        setImplicitSolventModel(igb_in, dielectric_in, saltcon_in, AtomicRadiusSet::MBONDI3,
                                policy);        
      }
      return;
    }
  }
  double* alpha_ptr   = gb_alpha_parameters.data();
  double* beta_ptr    = gb_beta_parameters.data();
  double* gamma_ptr   = gb_gamma_parameters.data();
  double* screen_ptr  = gb_screening_factors.data();
  float* sp_alpha_ptr  = sp_gb_alpha_parameters.data();
  float* sp_beta_ptr   = sp_gb_beta_parameters.data();
  float* sp_gamma_ptr  = sp_gb_gamma_parameters.data();
  float* sp_screen_ptr = sp_gb_screening_factors.data();
  const int* znum_ptr = atomic_numbers.data();

  // Set the dielectric constant and salt concentration, if appropriate
  switch (gb_style) {
  case ImplicitSolventModel::NONE:
    break;
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    dielectric_constant = dielectric_in;
    salt_concentration = saltcon_in;
    break;
  }

  // Set the radius set and screening parameters (a radius set of NONE will leave these
  // parameters unchanged in the topology)
  const int* nb12_bounds_ptr = nb12_exclusion_bounds.data();
  const int* nb12_excl_ptr = nb12_exclusion_list.data();
  const int* nb13_bounds_ptr = nb13_exclusion_bounds.data();
  const int* nb13_excl_ptr = nb13_exclusion_list.data();
  const char4* atom_type_ptr = atom_types.data();
  const char4* atom_name_ptr = atom_names.data();
  const char4* res_name_ptr = residue_names.data();
  const double* mass_ptr = atomic_masses.data();
  const std::vector<int> atom_residue_idx = getResidueIndex();
  for (int i = 0; i < atom_count; i++) {

    // Impart the radius set 
    double atom_rad = atomic_pb_radii.readHost(i);
    switch (radii_set) {
    case AtomicRadiusSet::BONDI:
    case AtomicRadiusSet::AMBER6:
    case AtomicRadiusSet::MBONDI:
    case AtomicRadiusSet::MBONDI2:
    case AtomicRadiusSet::MBONDI3:
      switch (znum_ptr[i]) {
      case 1:
        atom_rad = 1.2;
        if (nb12_bounds_ptr[i + 1] - nb12_bounds_ptr[i] > 0) {
          const int bonded_atom_idx = nb12_excl_ptr[nb12_bounds_ptr[i]];
          if (radii_set == AtomicRadiusSet::AMBER6 || radii_set == AtomicRadiusSet::MBONDI) {
            switch (znum_ptr[bonded_atom_idx]) {
            case 1:
              if (strncmpCased(char4ToString(atom_type_ptr[bonded_atom_idx]).c_str(), "HW", 2,
                               CaseSensitivity::NO)) {
                atom_rad = 0.8;
              }
              break;
            case 6:
              atom_rad = 1.3;
              break;
            case 7:
              if (radii_set == AtomicRadiusSet::MBONDI) {
                atom_rad = 1.3;
              }
              break;
            case 8:
            case 16:
              atom_rad = 0.8;
              break;
            default:
              break;
            }
          }
          else if (radii_set == AtomicRadiusSet::MBONDI2 ||
                   radii_set == AtomicRadiusSet::MBONDI3) {
            if (znum_ptr[bonded_atom_idx] == 7) {
              atom_rad = 1.3;
              const char4 atmc4 = atom_name_ptr[i];
              if (radii_set == AtomicRadiusSet::MBONDI3 &&
                  char4ToString(res_name_ptr[atom_residue_idx[bonded_atom_idx]]) == "ARG " &&
                  (atmc4.x == 'H' && (atmc4.y == 'H' || atmc4.y == 'E'))) {
                atom_rad = 1.17;
              }
            }
          }
        }
        else {
          switch (policy) {
          case ExceptionResponse::DIE:
            break;
          case ExceptionResponse::WARN:
            rtErr("Unbonded hydrogen atom " + char4ToString(atom_name_ptr[i]) + " detected in "
                  "residue " + char4ToString(res_name_ptr[atom_residue_idx[i]]) + ".  This "
                  "hydrogen's radius will keep the default Bondi value of 1.2 Angstroms.\n",
                  "AtomGraph", "setImplicitSolventModel");
            break;
          case ExceptionResponse::SILENT:
            break;
          }
        }
        break;
      case 6:
        {
          // Detect repartitioning in hydrogen masses when calculating the atom's true mass
          double atomi_mass = mass_ptr[i];
          for (int j = nb12_bounds_ptr[i]; j < nb12_bounds_ptr[i + 1]; j++) {
            const int neighbor_atom_idx = nb12_excl_ptr[j];
            if (znum_ptr[neighbor_atom_idx] == 1 && mass_ptr[neighbor_atom_idx] > 1.008) {
              atomi_mass += mass_ptr[neighbor_atom_idx] - 1.008;
            }
          }

          // Assuming all masses are based on the natural abundances, identify unified carbon
          // atoms as having (mass of C + mass of one or more H).
          const char4 attc4 = atom_type_ptr[i];
          if ((attc4.x == 'C' && attc4.y == '1' && atomi_mass > 13.018) ||
              (attc4.x == 'C' && attc4.y == '2' && atomi_mass > 14.026) ||
              (attc4.x == 'C' && attc4.y == '3' && atomi_mass > 15.034)) {

            // United atom carbon radius
            atom_rad = 2.2;
          }
          else {

            // Standard carbon radius
            atom_rad = 1.7;
          }
        }
        break;
      case 7:
        atom_rad = 1.55;
        break;
      case 8:
        atom_rad = 1.5;
        if (radii_set == AtomicRadiusSet::MBONDI3) {

          // Adjust carboxylic oxygens on proteins (side-chains as well as C-termini)
          const std::string resc4 = char4ToString(res_name_ptr[atom_residue_idx[i]]);
          const std::string atmc4 = char4ToString(atom_name_ptr[i]);
          if (((resc4 == "ASP " || resc4 == "AS4 ") && (atmc4 == "OD1 " || atmc4 == "OD2 ")) ||
              ((resc4 == "GLU " || resc4 == "GL4 ") && (atmc4 == "OE1 " || atmc4 == "OE2 ")) ||
              atmc4 == "OXT ") {
            atom_rad = 1.4;
          }
          if (znum_ptr[i] == 8 && nb12_bounds_ptr[i + 1] - nb12_bounds_ptr[i] == 1 &&
              znum_ptr[nb12_excl_ptr[nb12_bounds_ptr[i]]] == 6) {
            for (int j = nb13_bounds_ptr[i]; j < nb13_bounds_ptr[i + 1]; j++) {
              if (znum_ptr[nb13_excl_ptr[j]] == 8) {
                atom_rad = 1.4;
              }
            }
          }
        }
        break;
      case 9:
        atom_rad = 1.5;
        break;
      case 14:
        atom_rad = 2.1;
        break;
      case 15:
        atom_rad = 1.85;
        break;
      case 16:
        atom_rad = 1.8;
        break;
      case 17:
        atom_rad = 1.7;
        break;
      default:
        atom_rad = 1.5;
        break;
      }
      break;
    case AtomicRadiusSet::PARSE:
      switch (znum_ptr[i]) {
      case 1:
        atom_rad = 1.00;
        break;
      case 6:
        atom_rad = 1.7;
        break;
      case 7:
        atom_rad = 1.5;
        break;
      case 8:
        atom_rad = 1.4;
        break;
      case 16:
        atom_rad = 1.85;
        break;
      default:
        atom_rad = 1.5;
        break;
      }
      break;
    case AtomicRadiusSet::NONE:
      break;
    }
    atomic_pb_radii.putHost(atom_rad, i);

    // Impart the GB screening parameters
    double atom_screen = gb_screening_factors.readHost(i);
    switch (radii_set) {
    case AtomicRadiusSet::BONDI:
    case AtomicRadiusSet::AMBER6:
    case AtomicRadiusSet::MBONDI:
    case AtomicRadiusSet::MBONDI2:
    case AtomicRadiusSet::MBONDI3:
      switch (znum_ptr[i]) {
      case 1:
        atom_screen = 0.85;
        break;
      case 6:
        atom_screen = 0.72;
        break;
      case 7:
        atom_screen = 0.79;
        break;
      case 8:
        atom_screen = 0.85;
        break;
      case 9:
        atom_screen = 0.88;
        break;
      case 15:
        atom_screen = 0.86;
        break;
      case 16:
        atom_screen = 0.96;
        break;
      default:
        atom_screen = 0.8;
        break;
      }
    case AtomicRadiusSet::PARSE:
    case AtomicRadiusSet::NONE:
      break;
    }
    gb_screening_factors.putHost(atom_screen, i);
  }

  // Set alpha, beta, gamma, and screening factors if appropriate
  switch (gb_style) {
  case ImplicitSolventModel::NONE:
    break;
  case ImplicitSolventModel::HCT_GB:
    for (int i = 0; i < atom_count; i++) {
      alpha_ptr[i] = 0.0;
      beta_ptr[i]  = 0.0;
      gamma_ptr[i] = 0.0;
    }
    break;
  case ImplicitSolventModel::OBC_GB:
    for (int i = 0; i < atom_count; i++) {
      alpha_ptr[i] = gb_obc_i_alpha;
      beta_ptr[i]  = gb_obc_i_beta;
      gamma_ptr[i] = gb_obc_i_gamma;
    }
    break;
  case ImplicitSolventModel::OBC_GB_II:
    for (int i = 0; i < atom_count; i++) {
      alpha_ptr[i] = gb_obc_ii_alpha;
      beta_ptr[i]  = gb_obc_ii_beta;
      gamma_ptr[i] = gb_obc_ii_gamma;
    }
    break;
  case ImplicitSolventModel::NECK_GB:
    for (int i = 0; i < atom_count; i++) {
      alpha_ptr[i] = gb_neck_i_alpha;
      beta_ptr[i]  = gb_neck_i_beta;
      gamma_ptr[i] = gb_neck_i_gamma;
      switch(znum_ptr[i]) {
      case 1:
        screen_ptr[i] = gb_neck_i_screen_h;
        break;
      case 6:
        screen_ptr[i] = gb_neck_i_screen_c;
        break;
      case 7:
        screen_ptr[i] = gb_neck_i_screen_n;
        break;
      case 8:
        screen_ptr[i] = gb_neck_i_screen_o;
        break;
      case 16:
        screen_ptr[i] = gb_neck_i_screen_s;
        break;
      default:
        screen_ptr[i] = gb_neck_i_screen_default;
        break;
      }
    }
    break;
  case ImplicitSolventModel::NECK_GB_II:
    for (int i = 0; i < atom_count; i++) {
      switch(znum_ptr[i]) {
      case 1:
        alpha_ptr[i]  = gb_neck_ii_alpha_h;
        beta_ptr[i]   = gb_neck_ii_beta_h;
        gamma_ptr[i]  = gb_neck_ii_gamma_h;
        screen_ptr[i] = gb_neck_ii_screen_h;
        break;
      case 6:
        alpha_ptr[i]  = gb_neck_ii_alpha_c;
        beta_ptr[i]   = gb_neck_ii_beta_c;
        gamma_ptr[i]  = gb_neck_ii_gamma_c;
        screen_ptr[i] = gb_neck_ii_screen_c;
        break;
      case 7:
        alpha_ptr[i]  = gb_neck_ii_alpha_n;
        beta_ptr[i]   = gb_neck_ii_beta_n;
        gamma_ptr[i]  = gb_neck_ii_gamma_n;
        screen_ptr[i] = gb_neck_ii_screen_n;
        break;
      case 8:
        alpha_ptr[i]  = gb_neck_ii_alpha_os;
        beta_ptr[i]   = gb_neck_ii_beta_os;
        gamma_ptr[i]  = gb_neck_ii_gamma_os;
        screen_ptr[i] = gb_neck_ii_screen_o;
        break;
      case 15:
        alpha_ptr[i]  = gb_neck_ii_alpha_p;
        beta_ptr[i]   = gb_neck_ii_beta_p;
        gamma_ptr[i]  = gb_neck_ii_gamma_p;
        screen_ptr[i] = gb_neck_ii_screen_p;
        break;
      case 16:
        alpha_ptr[i]  = gb_neck_ii_alpha_os;
        beta_ptr[i]   = gb_neck_ii_beta_os;
        gamma_ptr[i]  = gb_neck_ii_gamma_os;
        screen_ptr[i] = gb_neck_ii_screen_s;
        break;
      default:
        alpha_ptr[i]  = 1.0;
        beta_ptr[i]   = 0.8;
        gamma_ptr[i]  = 4.85;
        screen_ptr[i] = 0.5;
        break;
      }
    }
    break;
  }

  // Copy double-precision values to single-precision arrays
  double* rad_ptr = atomic_pb_radii.data();
  float* sp_rad_ptr = sp_atomic_pb_radii.data();
  for (int i = 0; i < atom_count; i++) {
    sp_alpha_ptr[i]  = alpha_ptr[i];
    sp_beta_ptr[i]   = beta_ptr[i];
    sp_gamma_ptr[i]  = gamma_ptr[i];
    sp_screen_ptr[i] = screen_ptr[i];
    sp_rad_ptr[i]    = rad_ptr[i];
  }

  // Note the radius set in the topology, if it has indeed changed
  if (pb_radii_set.size() == 0 || radii_set != AtomicRadiusSet::NONE) {
    pb_radii_set = getEnumerationName(radii_set);
  }

  // Compute the neck GB indices based on the baseline atomic PB radii.  These values must later
  // be checked against the available table size.
  int* neck_idx_ptr = neck_gb_indices.data();
  double* radii_ptr = atomic_pb_radii.data();
  switch (gb_style) {
  case ImplicitSolventModel::NONE:
  case ImplicitSolventModel::HCT_GB:
  case ImplicitSolventModel::OBC_GB:
  case ImplicitSolventModel::OBC_GB_II:
    break;
  case ImplicitSolventModel::NECK_GB:
  case ImplicitSolventModel::NECK_GB_II:
    for (int i = 0; i < atom_count; i++) {
      neck_idx_ptr[i] = static_cast<int>(((radii_ptr[i] - 1.0) * 20.0) + 0.5);
    }
    break;
  }
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::setWaterResidueName(const char4 new_name) {
  water_residue_name = new_name;
}

//-------------------------------------------------------------------------------------------------
void AtomGraph::setWaterResidueName(const std::string &new_name) {
  if (new_name.size() > 4) {
    rtErr("The proposed water model name (" + new_name + ") cannot contain more than four "
          "characters.", "AtomGraph", "setWaterResidueName");
  }
  water_residue_name = stringToChar4(new_name);
}

} // namespace topology
} // namespace stormm
