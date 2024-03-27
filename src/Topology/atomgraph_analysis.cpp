#include <cmath>
#include "copyright.h"
#include "Constants/scaling.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "UnitTesting/approx.h"
#include "atomgraph_abstracts.h"
#include "atomgraph_analysis.h"
#include "atomgraph_enumerators.h"

namespace stormm {
namespace topology {

using stmath::findBin;
using stmath::accumulateBitmask;
using stmath::unsetBitInMask;
using stmath::readBitFromMask;
using parse::char4ToString;
using testing::Approx;

//-------------------------------------------------------------------------------------------------
WaterModel identifyWaterModel(const AtomGraph &ag) {

  // Obtain the necessary abstracts
  const ValenceKit<double> vnkit = ag.getDoublePrecisionValenceKit();
  const ChemicalDetailsKit cdkit = ag.getChemicalDetailsKit();
  const NonbondedKit<double> nbkit = ag.getDoublePrecisionNonbondedKit();

  // Search for molecules with the chemical formula H2O.
  std::vector<int> waters;
  for (int i = 0; i < cdkit.nmol; i++) {

    // Determine whether the molecule is water
    int n_hydrogen = 0;
    int n_oxygen = 0;
    int n_other_atoms = 0;
    for (int j = cdkit.mol_limits[i]; j < cdkit.mol_limits[i + 1]; j++) {
      const int element_number = cdkit.z_numbers[cdkit.mol_contents[j]];
      n_hydrogen += (element_number == 1);
      n_oxygen += (element_number == 8);
      n_other_atoms += (element_number > 0 && element_number != 1 && element_number != 8);
    }
    if (n_hydrogen == 2 && n_oxygen == 1 && n_other_atoms == 0) {
      waters.push_back(i);
    }
  }

  // Possible problems that may be encountered
  bool non_neutral_water = false;

  // Search the water molecules found for consistency and try to identify a particular model.
  // Return immediately if there is no water in the system.
  const int n_waters = waters.size();
  if (n_waters == 0) {
    return WaterModel::NONE;
  }
  std::vector<WaterModel> water_types(n_waters, WaterModel::UNKNOWN);
  for (int i = 0; i < n_waters; i++) {
    const int mol_contents_llim = cdkit.mol_limits[waters[i]];
    const int mol_contents_hlim = cdkit.mol_limits[waters[i] + 1];

    // Only three, four, and five-point water models are covered.
    if (mol_contents_hlim - mol_contents_llim < 3 || mol_contents_hlim - mol_contents_llim > 5) {
      continue;
    }

    // Find the non-bonded parameters, O-H equilibrium length, H-H equilibrium length (if
    // applicable), and H-O-H equilibrium angle (if applicable)
    double h_charge = 0.0;
    double h_ljsig  = 0.0;
    double h_ljeps  = 0.0;
    double o_charge = 0.0;
    double o_ljsig  = 0.0;
    double o_ljeps  = 0.0;
    double ep_charge = 0.0;
    double oh_leq = 0.0;
    double hh_leq = 0.0;
    double hoh_leq = 0.0;
    double oh_keq = 0.0;
    double hh_keq = 0.0;
    double hoh_keq = 0.0;
    double nep = 0.0;
    for (int j = mol_contents_llim; j < mol_contents_hlim; j++) {
      const int j_particle = cdkit.mol_contents[j];
      const int znum_j = cdkit.z_numbers[j_particle];

      // Seek hydrogen, oxygen, and massless site non-bonded parameters
      if (znum_j == 8) {
        o_charge = nbkit.charge[j_particle];
        const int o_ljidx = nbkit.lj_idx[j_particle] * (nbkit.n_lj_types + 1);
        const double o_lja = nbkit.lja_coeff[o_ljidx];
        const double o_ljb = nbkit.ljb_coeff[o_ljidx];
        if (o_ljb > constants::tiny) {
          o_ljsig = pow(o_lja / o_ljb, 1.0 / 6.0);
          o_ljeps = 0.25 * o_lja / pow(o_ljsig, 12.0);
        }
        else {
          o_ljsig = 0.0;
          o_ljeps = 0.0;
        }
      }
      else if (znum_j == 1) {
        h_charge = nbkit.charge[j_particle];
        const int h_ljidx = nbkit.lj_idx[j_particle] * (nbkit.n_lj_types + 1);
        const double h_lja = nbkit.lja_coeff[h_ljidx];
        const double h_ljb = nbkit.ljb_coeff[h_ljidx];
        if (h_ljb > constants::tiny) {
          h_ljsig = pow(h_lja / h_ljb, 1.0 / 6.0);
          h_ljeps = 0.25 * h_lja / pow(h_ljsig, 12.0);
        }
        else {
          h_ljsig = 0.0;
          h_ljeps = 0.0;
        }
      }
      else if (znum_j == 0) {
        ep_charge = nbkit.charge[j_particle];
        nep += 1.0;
      }

      // Find bond equilibrium lengths and stiffnesses
      for (int k = vnkit.bond_asgn_bounds[j_particle];
           k < vnkit.bond_asgn_bounds[j_particle + 1]; k++) {
        const int k_particle = vnkit.bond_asgn_atoms[k];
        const int znum_k = cdkit.z_numbers[k_particle];
        if ((znum_j == 8 && znum_k == 1) || (znum_j == 1 && znum_k == 8)) {
          const int oh_param_idx = vnkit.bond_asgn_index[k];
          oh_keq = vnkit.bond_keq[oh_param_idx];
          oh_leq = vnkit.bond_leq[oh_param_idx];
        }
        else if (znum_j == 1 && znum_k == 1) {
          const int hh_param_idx = vnkit.bond_asgn_index[k];
          hh_keq = vnkit.bond_keq[hh_param_idx];
          hh_leq = vnkit.bond_leq[hh_param_idx];
        }
      }
    }

    // Check that the charges make sense -- this will fail in cases like "Six-Site" water or
    // TIP7P, which do have more than one flavor of virtual site, but these are extremely rare
    // edge cases and the water model will just remain "UNKNOWN."
    if (o_charge + (2.0 * h_charge) + (nep * ep_charge) != Approx(0.0).margin(1.0e-5)) {
      non_neutral_water = true;
      continue;
    }

    // Identify this water molecule
    WaterModel charge_details = WaterModel::UNKNOWN;
    if (ep_charge == Approx(0.0).margin(1.0e-3)) {
      if (o_charge == Approx(-0.8952).margin(1.0e-3) &&
          h_charge == Approx(0.4476).margin(1.0e-3)) {
        charge_details = WaterModel::OPC3;
      }
      else if (o_charge == Approx(-0.82).margin(1.0e-3) &&
               h_charge == Approx(0.41).margin(1.0e-3)) {

        // Shared by SPC_FW
        charge_details = WaterModel::SPC;
      }
      else if (o_charge == Approx(-0.8476).margin(1.0e-3) &&
               h_charge == Approx(0.4238).margin(1.0e-3)) {

        // Shared by SPC_EB
        charge_details = WaterModel::SPC_E;
      }
      else if (o_charge == Approx(-0.870).margin(1.0e-3) &&
               h_charge == Approx(0.435).margin(1.0e-3)) {
        charge_details = WaterModel::SPC_HW;
      }
      else if (o_charge == Approx(-0.834).margin(1.0e-3) &&
               h_charge == Approx(0.417).margin(1.0e-3)) {

        // Shared by TIP3P_EW and TIP3P_FW
        charge_details = WaterModel::TIP3P;
      }
    }
    else {
      if (o_charge == Approx(0.0).margin(1.0e-3)) {
        if (ep_charge == Approx(-1.3582).margin(1.0e-3) &&
            h_charge == Approx(0.6791).margin(1.0e-3)) {
          charge_details = WaterModel::OPC;
        }
        else if (ep_charge == Approx(-1.04).margin(1.0e-3) &&
                 h_charge == Approx(0.52).margin(1.0e-3)) {
          charge_details = WaterModel::TIP4P;
        }
        else if (ep_charge == Approx(-1.04844).margin(1.0e-3) &&
                 h_charge == Approx(0.52422).margin(1.0e-3)) {
          charge_details = WaterModel::TIP4P_EW;
        }
        else if (ep_charge == Approx(-1.1128).margin(1.0e-3) &&
                 h_charge == Approx(0.5564).margin(1.0e-3)) {

          // Shared by TIP4P_2005F
          charge_details = WaterModel::TIP4P_2005;
        }
        else if (ep_charge == Approx(-1.1794).margin(1.0e-3) &&
                 h_charge == Approx(0.5897).margin(1.0e-3)) {
          charge_details = WaterModel::TIP4P_ICE;
        }
        else if (ep_charge == Approx(-1.054).margin(1.0e-3) &&
                 h_charge == Approx(0.527).margin(1.0e-3)) {
          charge_details = WaterModel::TIP4P_EPS;
        }
        else if (ep_charge == Approx(-1.16).margin(1.0e-3) &&
                 h_charge == Approx(0.58).margin(1.0e-3)) {
          charge_details = WaterModel::TIP4P_D;
        }
        else if (ep_charge == Approx(-0.241).margin(1.0e-3) &&
                 h_charge == Approx(0.241).margin(1.0e-3)) {

          // Shared by TIP5P_EW
          charge_details = WaterModel::TIP5P;
        }
      }
      else {

        // TIP5P-2018 is the lone model that puts charge on both the oxygen and the virtual sites
        if (o_charge == Approx(-0.641114).margin(1.0e-3) &&
            ep_charge == Approx(-0.07358).margin(1.0e-3) &&
            h_charge == Approx(0.394137).margin(1.0e-3)) {
          charge_details = WaterModel::TIP5P_2018;
        }
        else {
          charge_details = WaterModel::UNKNOWN;
        }
      }
    }

    // Continue to identify Lennard-Jones properties
    WaterModel lennard_jones_details = WaterModel::UNKNOWN;
    if (h_ljsig < constants::tiny) {
      if (o_ljsig == Approx(3.17427).margin(1.0e-3) &&
          o_ljeps == Approx(0.163406).margin(1.0e-3)) {
        lennard_jones_details = WaterModel::OPC3;
      }
      else if (o_ljsig == Approx(3.16572).margin(1.0e-3) &&
               o_ljeps == Approx(0.155354).margin(1.0e-3)) {

        // Also SPC_E, SPC_EB, SPC_HW, and SPC_FW
        lennard_jones_details = WaterModel::SPC;
      }
      else if (o_ljsig == Approx(3.1507).margin(1.0e-3) &&
               o_ljeps == Approx(0.1521).margin(1.0e-3)) {

        // Also TIP3P_FW
        lennard_jones_details = WaterModel::TIP3P;
      }
      else if (o_ljsig == Approx(3.19405).margin(1.0e-3) &&
               o_ljeps == Approx(0.0980).margin(1.0e-3)) {
        lennard_jones_details = WaterModel::TIP3P_EW;
      }
      else if (o_ljsig == Approx(3.16655).margin(1.0e-3) &&
               o_ljeps == Approx(0.21280).margin(1.0e-3)) {
        lennard_jones_details = WaterModel::OPC;
      }
      else if (o_ljsig == Approx(3.15365).margin(1.0e-3) &&
               o_ljeps == Approx(0.1550).margin(1.0e-3)) {
        lennard_jones_details = WaterModel::TIP4P;
      }
      else if (o_ljsig == Approx(3.16435).margin(1.0e-3) &&
               o_ljeps == Approx(0.16275).margin(1.0e-3)) {
        lennard_jones_details = WaterModel::TIP4P_EW;
      }
      else if (o_ljsig == Approx(3.1589).margin(1.0e-3) &&
               o_ljeps == Approx(0.18521).margin(1.0e-3)) {
        lennard_jones_details = WaterModel::TIP4P_2005;
      }
      else if (o_ljsig == Approx(3.1668).margin(1.0e-3) &&
               o_ljeps == Approx(0.21085).margin(1.0e-3)) {
        lennard_jones_details = WaterModel::TIP4P_ICE;
      }
      else if (o_ljsig == Approx(3.1644).margin(1.0e-3) &&
               o_ljeps == Approx(0.18521).margin(1.0e-3)) {
        lennard_jones_details = WaterModel::TIP4P_2005F;
      }
      else if (o_ljsig == Approx(3.165).margin(1.0e-3) &&
               o_ljeps == Approx(0.1848).margin(1.0e-3)) {
        lennard_jones_details = WaterModel::TIP4P_EPS;
      }
      else if (o_ljsig == Approx(3.165).margin(1.0e-3) &&
               o_ljeps == Approx(0.22384).margin(1.0e-3)) {
        lennard_jones_details = WaterModel::TIP4P_D;
      }
      else if (o_ljsig == Approx(3.12).margin(1.0e-3) &&
               o_ljeps == Approx(0.16).margin(1.0e-3)) {
        lennard_jones_details = WaterModel::TIP5P;
      }
      else if (o_ljsig == Approx(3.097).margin(1.0e-3) &&
               o_ljeps == Approx(0.17801).margin(1.0e-3)) {
        lennard_jones_details = WaterModel::TIP5P_EW;
      }
      else if (o_ljsig == Approx(3.145).margin(1.0e-3) &&
               o_ljeps == Approx(0.1888).margin(1.0e-3)) {
        lennard_jones_details = WaterModel::TIP5P_2018;
      }
    }
    else {
      if (o_ljsig == Approx(3.1507).margin(1.0e-3) && o_ljeps == Approx(0.1521).margin(1.0e-3) &&
          h_ljsig == Approx(0.4000).margin(1.0e-3) && h_ljeps == Approx(0.0460).margin(1.0e-3)) {
        lennard_jones_details = WaterModel::TIP3P_CHARMM;
      }
    }

    // Finally, identify the geometric properties
    WaterModel geometry_details = WaterModel::UNKNOWN;
    if (oh_leq == Approx(0.8724).margin(1.0e-3)) {
      geometry_details = WaterModel::OPC;
    }
    else if (oh_leq == Approx(1.0).margin(1.0e-3)) {

      // Also SPC_E and SPC_HW
      geometry_details = WaterModel::SPC;
    }
    else if (oh_leq == Approx(1.012).margin(1.0e-3)) {
      geometry_details = WaterModel::SPC_FW;
    }
    else if (oh_leq == Approx(1.01).margin(1.0e-3)) {
      geometry_details = WaterModel::SPC_EB;
    }
    else if (oh_leq == Approx(0.9572).margin(1.0e-3)) {

      // Also TIP3P_CHARMM, TIP3P_EW, TIP4P, TIP4P_EW, TIP4P_2005, TIP4P_ICE, TIP4P_EPS,
      // TIP4P_D, TIP4P_F, TIP4P_2005F, TIP5P, TIP5P_EW, TIP5P_2018
      geometry_details = WaterModel::TIP3P;
    }
    else if (oh_leq == Approx(0.96).margin(1.0e-3)) {
      geometry_details = WaterModel::TIP3P_FW;
    }

    // Verify that all of the details are consistent with one water model
    WaterModel this_water = WaterModel::UNKNOWN;
    switch (charge_details) {
    case WaterModel::NONE:
    case WaterModel::UNKNOWN:
    case WaterModel::CHIMERA:
    case WaterModel::MULTIPLE:
    case WaterModel::MULTI_CHIMERA:
      break;
    case WaterModel::OPC3:
      if (lennard_jones_details == WaterModel::OPC3 && geometry_details == WaterModel::OPC3) {
        this_water = WaterModel::OPC3;
      }
      else if (lennard_jones_details != WaterModel::UNKNOWN ||
               geometry_details != WaterModel::UNKNOWN) {
        this_water = WaterModel::CHIMERA;
      }
      break;
    case WaterModel::SPC:
    case WaterModel::SPC_FW:
      if (lennard_jones_details == WaterModel::SPC && geometry_details == WaterModel::SPC) {
        this_water = WaterModel::SPC;
      }
      else if (lennard_jones_details == WaterModel::SPC &&
               geometry_details == WaterModel::SPC_FW) {
        this_water = WaterModel::SPC_FW;
      }
      else if (lennard_jones_details != WaterModel::UNKNOWN ||
               geometry_details != WaterModel::UNKNOWN) {
        this_water = WaterModel::CHIMERA;
      }
      break;
    case WaterModel::SPC_E:
    case WaterModel::SPC_EB:
      if (lennard_jones_details == WaterModel::SPC && geometry_details == WaterModel::SPC) {
        this_water = WaterModel::SPC_E;
      }
      else if (lennard_jones_details == WaterModel::SPC &&
               geometry_details == WaterModel::SPC_EB) {
        this_water = WaterModel::SPC_EB;
      }
      else if (lennard_jones_details != WaterModel::UNKNOWN ||
               geometry_details != WaterModel::UNKNOWN) {
        this_water = WaterModel::CHIMERA;
      }
      break;
    case WaterModel::SPC_HW:
      if (lennard_jones_details == WaterModel::SPC && geometry_details == WaterModel::SPC) {
        this_water = WaterModel::SPC_HW;
      }
      else if (lennard_jones_details != WaterModel::UNKNOWN ||
               geometry_details != WaterModel::UNKNOWN) {
        this_water = WaterModel::CHIMERA;
      }
      break;
    case WaterModel::TIP3P:
    case WaterModel::TIP3P_FW:
    case WaterModel::TIP3P_EW:
    case WaterModel::TIP3P_CHARMM:
      if (lennard_jones_details == WaterModel::TIP3P && geometry_details == WaterModel::TIP3P) {
        this_water = WaterModel::TIP3P;
      }
      else if (lennard_jones_details == WaterModel::TIP3P &&
               geometry_details == WaterModel::TIP3P_FW) {
        this_water = WaterModel::TIP3P_FW;
      }
      else if (lennard_jones_details == WaterModel::TIP3P_EW &&
               geometry_details == WaterModel::TIP3P) {
        this_water = WaterModel::TIP3P_EW;
      }
      else if (lennard_jones_details == WaterModel::TIP3P_CHARMM &&
               geometry_details == WaterModel::TIP3P) {
        this_water = WaterModel::TIP3P_CHARMM;
      }
      else if (lennard_jones_details != WaterModel::UNKNOWN ||
               geometry_details != WaterModel::UNKNOWN) {
        this_water = WaterModel::CHIMERA;
      }
      break;
    case WaterModel::OPC:
      if (lennard_jones_details == WaterModel::OPC && geometry_details == WaterModel::OPC) {
        this_water = WaterModel::OPC;
      }
      else if (lennard_jones_details != WaterModel::UNKNOWN ||
               geometry_details != WaterModel::UNKNOWN) {
        this_water = WaterModel::CHIMERA;
      }
      break;
    case WaterModel::TIP4P:
      if (lennard_jones_details == WaterModel::TIP4P && geometry_details == WaterModel::TIP3P) {
        this_water = WaterModel::TIP4P;
      }
      else if (lennard_jones_details != WaterModel::UNKNOWN ||
               geometry_details != WaterModel::UNKNOWN) {
        this_water = WaterModel::CHIMERA;
      }
      break;
    case WaterModel::TIP4P_EW:
      if (lennard_jones_details == WaterModel::TIP4P_EW && geometry_details == WaterModel::TIP3P) {
        this_water = WaterModel::TIP4P_EW;
      }
      else if (lennard_jones_details != WaterModel::UNKNOWN ||
               geometry_details != WaterModel::UNKNOWN) {
        this_water = WaterModel::CHIMERA;
      }
      break;
    case WaterModel::TIP4P_2005:
      if (lennard_jones_details == WaterModel::TIP4P_2005 &&
          geometry_details == WaterModel::TIP3P) {
        this_water = WaterModel::TIP4P_2005;
      }
      else if (lennard_jones_details != WaterModel::UNKNOWN ||
               geometry_details != WaterModel::UNKNOWN) {
        this_water = WaterModel::CHIMERA;
      }
      break;
    case WaterModel::TIP4P_ICE:
      if (lennard_jones_details == WaterModel::TIP4P_ICE &&
          geometry_details == WaterModel::TIP3P) {
        this_water = WaterModel::TIP4P_ICE;
      }
      else if (lennard_jones_details != WaterModel::UNKNOWN ||
               geometry_details != WaterModel::UNKNOWN) {
        this_water = WaterModel::CHIMERA;
      }
      break;
    case WaterModel::TIP4P_EPS:
      if (lennard_jones_details == WaterModel::TIP4P_EPS &&
          geometry_details == WaterModel::TIP3P) {
        this_water = WaterModel::TIP4P_EPS;
      }
      else if (lennard_jones_details != WaterModel::UNKNOWN ||
               geometry_details != WaterModel::UNKNOWN) {
        this_water = WaterModel::CHIMERA;
      }
      break;
    case WaterModel::TIP4P_D:
      if (lennard_jones_details == WaterModel::TIP4P_D && geometry_details == WaterModel::TIP3P) {
        this_water = WaterModel::TIP4P_D;
      }
      else if (lennard_jones_details != WaterModel::UNKNOWN ||
               geometry_details != WaterModel::UNKNOWN) {
        this_water = WaterModel::CHIMERA;
      }
      break;
    case WaterModel::TIP4P_2005F:
      if (lennard_jones_details == WaterModel::TIP4P_2005F &&
          geometry_details == WaterModel::TIP3P) {
        this_water = WaterModel::TIP4P_2005F;
      }
      else if (lennard_jones_details != WaterModel::UNKNOWN ||
               geometry_details != WaterModel::UNKNOWN) {
        this_water = WaterModel::CHIMERA;
      }
      break;
    case WaterModel::TIP5P:
      if (lennard_jones_details == WaterModel::TIP5P && geometry_details == WaterModel::TIP3P) {
        this_water = WaterModel::TIP5P;
      }
      else if (lennard_jones_details != WaterModel::UNKNOWN ||
               geometry_details != WaterModel::UNKNOWN) {
        this_water = WaterModel::CHIMERA;
      }
      break;
    case WaterModel::TIP5P_EW:
      if (lennard_jones_details == WaterModel::TIP5P_EW && geometry_details == WaterModel::TIP3P) {
        this_water = WaterModel::TIP5P_EW;
      }
      else if (lennard_jones_details != WaterModel::UNKNOWN ||
               geometry_details != WaterModel::UNKNOWN) {
        this_water = WaterModel::CHIMERA;
      }
      break;
    case WaterModel::TIP5P_2018:
      if (lennard_jones_details == WaterModel::TIP5P_2018 &&
          geometry_details == WaterModel::TIP3P) {
        this_water = WaterModel::TIP5P_2018;
      }
      else if (lennard_jones_details != WaterModel::UNKNOWN ||
               geometry_details != WaterModel::UNKNOWN) {
        this_water = WaterModel::CHIMERA;
      }
      break;
    }

    // Catalog the result
    water_types[i] = this_water;
  }

  // Scan all waters and seek consensus
  bool water_homogeneous = true;
  for (int i = 1; i < n_waters; i++) {
    water_homogeneous = (water_homogeneous && water_types[i] == water_types[0]);
  }
  if (water_homogeneous) {
    return water_types[0];
  }
  else {

    // There are multiple water models.  See if there any of them are chimeric.
    bool has_chimera = false;
    for (int i = 0; i < n_waters; i++) {
      has_chimera = (has_chimera || water_types[i] == WaterModel::CHIMERA);
    }
    if (has_chimera) {
      return WaterModel::MULTI_CHIMERA;
    }
    else {
      return WaterModel::MULTIPLE;
    }
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string listVirtualSiteFrameTypes(const AtomGraph &ag) {
  const VirtualSiteKit<double> vsk = ag.getDoublePrecisionVirtualSiteKit();
  return listVirtualSiteFrameTypes(vsk.vs_types, vsk.nsite);
}

//-------------------------------------------------------------------------------------------------
std::string listVirtualSiteFrameTypes(const AtomGraph *ag) {
  const VirtualSiteKit<double> vsk = ag->getDoublePrecisionVirtualSiteKit();
  return listVirtualSiteFrameTypes(vsk.vs_types, vsk.nsite);
}

//-------------------------------------------------------------------------------------------------
std::string listVirtualSiteFrameTypes(const int* vs_types, const int nsite) {

  // Compile a list of all frame types found in this system.
  int n_unique = 0;
  std::vector<int> unique_frm;
  std::vector<bool> coverage(nsite);
  for (int i = 0; i < nsite; i++) {
    if (coverage[i]) {
      continue;
    }
    for (int j = i; j < nsite; j++) {
      coverage[j] = (coverage[j] || vs_types[j] == vs_types[i]);
    }
    unique_frm.push_back(vs_types[i]);
    n_unique++;
  }
  std::string frame_type_list;
  for (int i = 0; i < n_unique; i++) {
    frame_type_list += getEnumerationName(static_cast<VirtualSiteKind>(unique_frm[i]));
    if (i < n_unique - 1) {
      frame_type_list += ", ";
    }
  }
  return frame_type_list;
}

//-------------------------------------------------------------------------------------------------
int colorConnectivity(const NonbondedKit<double> &nbk, const ChemicalDetailsKit &cdk,
                      const int atom_i, const int atom_j, std::vector<uint> *marked,
                      bool *ring_report) {
  bool ring_completion = false;
  uint* marked_ptr = marked->data();
  const int nm_elem = marked->size();
  for (int i = 0; i < nm_elem; i++) {
    marked_ptr[i] = 0U;
  }
  accumulateBitmask(marked_ptr, atom_i);
  accumulateBitmask(marked_ptr, atom_j);
  std::vector<int> prev_atoms(16);
  std::vector<int> new_atoms(16);
  prev_atoms.resize(1);
  new_atoms.resize(0);
  prev_atoms[0] = atom_j;
  bool more_added = true;
  int nrot = 0;
  while (more_added) {
    more_added = false;
    new_atoms.resize(0);
    const int nprev = prev_atoms.size();
    for (int i = 0; i < nprev; i++) {
      const int jlim = nbk.nb12_bounds[prev_atoms[i] + 1];
      for (int j = nbk.nb12_bounds[prev_atoms[i]]; j < jlim; j++) {
	const int candidate_atom = nbk.nb12x[j];
        if (cdk.z_numbers[candidate_atom] == 0) {
          continue;
        }
        else if (readBitFromMask(marked_ptr, candidate_atom)) {
          ring_completion = (ring_completion ||
                             (candidate_atom == atom_i && prev_atoms[i] != atom_j));
          continue;
        }
        else {
          new_atoms.push_back(candidate_atom);
          accumulateBitmask(marked_ptr, candidate_atom);
        }
      }
    }
    const int nnew = new_atoms.size();
    nrot += nnew;
    if (nnew > 0) {
      more_added = true;
      prev_atoms.resize(nnew);
      for (int i = 0; i < nnew; i++) {
        prev_atoms[i] = new_atoms[i];
      }
    }
  }

  // Report whether a ring was completed by this coloring
  *ring_report = ring_completion;
  return nrot;
}
  
//-------------------------------------------------------------------------------------------------
std::vector<int> mapRotatingGroup(const NonbondedKit<double> &nbk, const ChemicalDetailsKit &cdk,
                                  const int atom_i, const int atom_j,
                                  const std::string &filename) {
  std::vector<uint> marked((nbk.natom + uint_bit_count_int - 1) / uint_bit_count_int, false);
  bool ring_completed;
  int nrot = colorConnectivity(nbk, cdk, atom_i, atom_j, &marked, &ring_completed);
  if (ring_completed) {
    rtWarn("The rotatable bond between atoms " + std::to_string(atom_i + 1) + " and " +
           std::to_string(atom_j + 1) + ", " +
           char4ToString(cdk.res_names[findBin(cdk.res_limits, atom_i, cdk.nres)]) + " " +
           std::to_string(cdk.res_numbers[atom_i]) + " :: " +
           char4ToString(cdk.atom_names[atom_i]) + " and " +
           char4ToString(cdk.res_names[findBin(cdk.res_limits, atom_j, cdk.nres)]) + " " +
           std::to_string(cdk.res_numbers[atom_j]) + " :: " +
           char4ToString(cdk.atom_names[atom_j]) + ", appears to be part of a ring system.  This "
           "should not have happened if the bond was selected from a ChemicalFeatures object, but "
           "may have occurred if there is a very complex fused ring system that could not be "
           "fully mapped.  No atoms will rotate about this bond.  Original topology: " +
           filename + ".", "selectRotatingAtoms");
    std::vector<int> tmp_result;
    return tmp_result;
  }
  if (nrot == 0 || nrot == nbk.natom - 2) {
    rtWarn("The rotatable bond between atoms " + std::to_string(atom_i + 1) + " and " +
           std::to_string(atom_j + 1) + ", " +
           char4ToString(cdk.res_names[findBin(cdk.res_limits, atom_i, cdk.nres)]) + " " +
           std::to_string(cdk.res_numbers[atom_i]) + " :: " +
           char4ToString(cdk.atom_names[atom_i]) + " and " +
           char4ToString(cdk.res_names[findBin(cdk.res_limits, atom_j, cdk.nres)]) + " " +
           std::to_string(cdk.res_numbers[atom_j]) + " :: " +
           char4ToString(cdk.atom_names[atom_j]) + ", does not appear to be worth rotating.  This "
           "should not have happened if the bond was selected from a ChemicalFeatures object.  "
           "Original topology: " + filename + ".", "selectRotatingAtoms");
  }
  std::vector<int> result(nrot);
  int counter = 0;
  unsetBitInMask(&marked, atom_i);
  unsetBitInMask(&marked, atom_j);
  for (int i = 0; i < nbk.natom; i++) {
    if (readBitFromMask(marked,i)) {
      result[counter] = i;
      counter++;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> selectRotatingAtoms(const AtomGraph &ag, const int atom_i, const int atom_j) {
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  const ChemicalDetailsKit cdk = ag.getChemicalDetailsKit();
  return mapRotatingGroup(nbk, cdk, atom_i, atom_j, ag.getFileName());
}

//-------------------------------------------------------------------------------------------------
std::vector<int> selectRotatingAtoms(const AtomGraph *ag, const int atom_i, const int atom_j) {
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  return mapRotatingGroup(nbk, cdk, atom_i, atom_j, ag->getFileName());
}

//-------------------------------------------------------------------------------------------------
int inferLennardJonesTypeCount(const int length_a, const int length_b, const char* caller) {
  if (length_a != length_b) {
    rtErr("Parameter matrices must have identical sizes (" + std::to_string(length_a) + " and " +
          std::to_string(length_b) + " provided).", caller);
  }
  const int n_lj_types = round(sqrt(length_a));
  if (n_lj_types * n_lj_types != static_cast<int>(length_a)) {
    rtErr("A number of atom types can only be inferred for square parameter matrices.", caller);
  }
  return n_lj_types;
}

//-------------------------------------------------------------------------------------------------
VdwCombiningRule inferCombiningRule(const AtomGraph *ag, const ExceptionResponse policy,
                                    const bool seek_prevalent) {
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  return inferCombiningRule<double>(nbk.lja_coeff, nbk.ljb_coeff, nbk.n_lj_types, policy,
                                    seek_prevalent);
}

//-------------------------------------------------------------------------------------------------
VdwCombiningRule inferCombiningRule(const AtomGraph &ag, const ExceptionResponse policy,
                                    const bool seek_prevalent) {
  const NonbondedKit<double> nbk = ag.getDoublePrecisionNonbondedKit();
  return inferCombiningRule<double>(nbk.lja_coeff, nbk.ljb_coeff, nbk.n_lj_types, policy,
                                    seek_prevalent);
}

//-------------------------------------------------------------------------------------------------
int getConstrainedDegreesOfFreedom(const AtomGraph *ag, const int low_atom_index,
                                   const int high_atom_index) {
  const int actual_high_atom_index = (high_atom_index <= 0) ? ag->getAtomCount() : high_atom_index;
  if (low_atom_index >= actual_high_atom_index) {
    rtErr("The atom range " + std::to_string(low_atom_index) + " to " +
          std::to_string(actual_high_atom_index) + " is invalid.",
          "getConstrainedDegreesOfFreedom");
  }
  if (low_atom_index < 0 || actual_high_atom_index > ag->getAtomCount()) {
    rtErr("The atom range " + std::to_string(low_atom_index) + " to " +
          std::to_string(actual_high_atom_index) + " is invalid for a topology of " +
          std::to_string(ag->getAtomCount()) + " atoms.", "getConstrainedDegreesOfFreedom");
  }
  int result = (3 * (actual_high_atom_index - low_atom_index) - 6);
  const ConstraintKit<double> cnk = ag->getDoublePrecisionConstraintKit();

  // Subtract one degree of freedom for every constrained bond between two atoms that are both part
  // of the subset in question.
  for (int i = 0; i < cnk.ngroup; i++) {
    const size_t llim = cnk.group_bounds[i];
    const size_t hlim = cnk.group_bounds[i + 1];
    if (cnk.group_list[llim] >= low_atom_index && cnk.group_list[llim] < high_atom_index) {
      for (size_t j = llim + 1; j < hlim; j++) {
        if (cnk.group_list[j] >= low_atom_index && cnk.group_list[j] < high_atom_index) {
          result--;
        }
      }
    }
  }

  // Subtract three degrees of freedom for every three-point rigid water molecule completely
  // comprised by the subset.
  for (int i = 0; i < cnk.nsettle; i++) {
    if (cnk.settle_ox_atoms[i] >= low_atom_index && cnk.settle_ox_atoms[i] < high_atom_index &&
        cnk.settle_h1_atoms[i] >= low_atom_index && cnk.settle_h1_atoms[i] < high_atom_index &&
        cnk.settle_h2_atoms[i] >= low_atom_index && cnk.settle_h2_atoms[i] < high_atom_index) {
      result -= 3;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int getConstrainedDegreesOfFreedom(const AtomGraph &ag, const int low_atom_index,
                                   const int high_atom_index) {
  return getConstrainedDegreesOfFreedom(ag.getSelfPointer(), low_atom_index, high_atom_index);
}

} // namespace topology
} // namespace stormm
