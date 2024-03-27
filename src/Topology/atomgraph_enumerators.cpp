#include <string>
#include "copyright.h"
#include "Parsing/parse.h"
#include "atomgraph_enumerators.h"

namespace stormm {
namespace topology {

using parse::strcmpCased;

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const TopologyKind input) {
  switch (input) {
  case TopologyKind::AMBER:
    return std::string("AMBER");
  case TopologyKind::CHARMM:
    return std::string("CHARMM");
  case TopologyKind::GROMACS:
    return std::string("GROMACS");
  case TopologyKind::OPENMM:
    return std::string("OPENMM");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const TopologyDescriptor input) {
  switch (input) {
  case TopologyDescriptor::ATOM_COUNT:
    return std::string("ATOM_COUNT");
  case TopologyDescriptor::ATOM_TYPE_COUNT:
    return std::string("ATOM_TYPE_COUNT");
  case TopologyDescriptor::BONDS_WITH_HYDROGEN:
    return std::string("BONDS_WITH_HYDROGEN");
  case TopologyDescriptor::BONDS_WITHOUT_HYDROGEN:
    return std::string("BONDS_WITHOUT_HYDROGEN");
  case TopologyDescriptor::ANGLES_WITH_HYDROGEN:
    return std::string("ANGLES_WITH_HYDROGEN");
  case TopologyDescriptor::ANGLES_WITHOUT_HYDROGEN:
    return std::string("ANGLES_WITHOUT_HYDROGEN");
  case TopologyDescriptor::DIHEDRALS_WITH_HYDROGEN:
    return std::string("DIHEDRALS_WITH_HYDROGEN");
  case TopologyDescriptor::DIHEDRALS_WITHOUT_HYDROGEN:
    return std::string("DIHEDRALS_WITHOUT_HYDROGEN");
  case TopologyDescriptor::NHPARM_UNUSED:
    return std::string("NHPARM_UNUSED");
  case TopologyDescriptor::ADDLES_CREATED:
    return std::string("ADDLES_CREATED");
  case TopologyDescriptor::TOTAL_EXCLUDED_ATOMS:
    return std::string("TOTAL_EXCLUDED_ATOMS");
  case TopologyDescriptor::RESIDUE_COUNT:
    return std::string("RESIUDE_COUNT");
  case TopologyDescriptor::NBONA_UNUSED:
    return std::string("NBONA_UNUSED");
  case TopologyDescriptor::NTHETA_UNUSED:
    return std::string("NTHETA_UNUSED");
  case TopologyDescriptor::NPHIA_UNUSED:
    return std::string("NPHIA_UNUSED");
  case TopologyDescriptor::BOND_TYPE_COUNT:
    return std::string("BOND_TYPE_COUNT");
  case TopologyDescriptor::ANGLE_TYPE_COUNT:
    return std::string("ANGLE_TYPE_COUNT");
  case TopologyDescriptor::DIHEDRAL_TYPE_COUNT:
    return std::string("DIHEDRAL_TYPE_COUNT");
  case TopologyDescriptor::NATYP_UNUSED:
    return std::string("NATYPE_UNUSED");
  case TopologyDescriptor::NPHB_UNUSED:
    return std::string("NPHB_UNUSED");
  case TopologyDescriptor::PERTURBATION:
    return std::string("PERTURBATION");
  case TopologyDescriptor::BOND_PERTURBATIONS:
    return std::string("BOND_PERTURBATIONS");
  case TopologyDescriptor::ANGLE_PERTURBATIONS:
    return std::string("ANGLE_PERTURBATIONS");
  case TopologyDescriptor::DIHEDRAL_PERTURBATIONS:
    return std::string("DIHEDRAL_PERTURBATIONS");
  case TopologyDescriptor::BONDS_IN_PERTURBED_GROUP:
    return std::string("BONDS_IN_PERTURBED_GROUP");
  case TopologyDescriptor::ANGLES_IN_PERTURBED_GROUP:
    return std::string("ANGLES_IN_PERTURBED_GROUP");
  case TopologyDescriptor::DIHEDRALS_IN_PERTURBED_GROUP:
    return std::string("DIHEDRALS_IN_PERTURBED_GROUP");
  case TopologyDescriptor::BOX_TYPE_INDEX:
    return std::string("BOX_TYPE_INDEX");
  case TopologyDescriptor::ATOM_COUNT_LARGEST_RESIDUE:
    return std::string("ATOM_COUNT_LARGEST_RESIDUE");
  case TopologyDescriptor::CAP:
    return std::string("CAP");
  case TopologyDescriptor::EXTRA_POINT_COUNT:
    return std::string("EXTRA_POINT_COUNT");
  case TopologyDescriptor::PIMD_SLICE_COUNT:
    return std::string("PIMD_SLICE_COUNT");
  case TopologyDescriptor::N_VALUES:
    return std::string("N_VALUES");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const SanderDescriptor input) {
  switch (input) {
  case SanderDescriptor::NATOM:
    return std::string("NATOM");
  case SanderDescriptor::NTYPES:
    return std::string("NTYPES");
  case SanderDescriptor::NBONH:
    return std::string("NBONH");
  case SanderDescriptor::MBONA:
    return std::string("MBONA");
  case SanderDescriptor::NTHETH:
    return std::string("NTHETH");
  case SanderDescriptor::MTHETA:
    return std::string("MTHETA");
  case SanderDescriptor::NPHIH:
    return std::string("NPHIH");
  case SanderDescriptor::MPHIA:
    return std::string("MPHIA");
  case SanderDescriptor::NHPARM:
    return std::string("NHPARM");
  case SanderDescriptor::NPARM:
    return std::string("NPARM");
  case SanderDescriptor::NNB:
    return std::string("NNB");
  case SanderDescriptor::NRES:
    return std::string("NRES");
  case SanderDescriptor::NBONA:
    return std::string("NBONA");
  case SanderDescriptor::NTHETA:
    return std::string("NTHETA");
  case SanderDescriptor::NPHIA:
    return std::string("NPHIA");
  case SanderDescriptor::NUMBND:
    return std::string("NUMBND");
  case SanderDescriptor::NUMANG:
    return std::string("NUMANG");
  case SanderDescriptor::NPTRA:
    return std::string("NPTRA");
  case SanderDescriptor::NATYP:
    return std::string("NATYP");
  case SanderDescriptor::NPHB:
    return std::string("NPHB");
  case SanderDescriptor::IFPERT:
    return std::string("IFPERT");
  case SanderDescriptor::NBPER:
    return std::string("NBPER");
  case SanderDescriptor::NGPER:
    return std::string("NGPER");
  case SanderDescriptor::NDPER:
    return std::string("NDPER");
  case SanderDescriptor::MBPER:
    return std::string("NBPER");
  case SanderDescriptor::MGPER:
    return std::string("MGPER");
  case SanderDescriptor::MDPER:
    return std::string("MDPER");
  case SanderDescriptor::IFBOX:
    return std::string("IFBOX");
  case SanderDescriptor::NMXRS:
    return std::string("NMXRS");
  case SanderDescriptor::IFCAP:
    return std::string("IFCAP");
  case SanderDescriptor::NUMEXTRA:
    return std::string("NMEXTRA");
  case SanderDescriptor::NCOPY:
    return std::string("NCOPY");
  case SanderDescriptor::N_VALUES:
    return std::string("N_VALUES");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const UnitCellType input) {
  switch (input) {
  case UnitCellType::NONE:
    return std::string("Isolated system (free in space)");
  case UnitCellType::ORTHORHOMBIC:
    return std::string("Rectilinear");
  case UnitCellType::TRICLINIC:
    return std::string("Non-rectilinear prism (at least one unit cell angle != 90 deg.)");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const MobilitySetting input) {
  switch (input) {
  case MobilitySetting::OFF:
    return std::string("OFF");
  case MobilitySetting::ON:
    return std::string("ON");
  case MobilitySetting::TOGGLE:
    return std::string("TOGGLE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ShakeSetting input) {
  switch (input) {
  case ShakeSetting::OFF:
    return std::string("OFF");
  case ShakeSetting::ON:
    return std::string("ON");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const SettleSetting input) {
  switch (input) {
  case SettleSetting::OFF:
    return std::string("OFF");
  case SettleSetting::ON:
    return std::string("ON");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const PerturbationSetting input) {
  switch (input) {
  case PerturbationSetting::OFF:
    return std::string("OFF");
  case PerturbationSetting::ON:
    return std::string("ON");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const SolventCapSetting input) {
  switch (input) {
  case SolventCapSetting::OFF:
    return std::string("OFF");
  case SolventCapSetting::ON:
    return std::string("ON");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const PolarizationSetting input) {
  switch (input) {
  case PolarizationSetting::OFF:
    return std::string("OFF");
  case PolarizationSetting::ON:
    return std::string("ON");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const TopologyRequirement input) {
  switch (input) {
  case TopologyRequirement::OPTIONAL:
    return std::string("OPTIONAL");
  case TopologyRequirement::ESSENTIAL:
    return std::string("ESSENTIAL");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const MassForm input) {
  switch (input) {
  case MassForm::ORDINARY:
    return std::string("ORDINARY");
  case MassForm::INVERSE:
    return std::string("INVERSE");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ConstraintStatus input) {
  switch (input) {
  case ConstraintStatus::FREE:
    return std::string("FREE");
  case ConstraintStatus::CONSTRAINED:
    return std::string("CONSTRAINED");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const HydrogenContent input) {
  switch (input) {
  case HydrogenContent::NO_HYDROGEN:
    return std::string("NO_HYDROGEN");
  case HydrogenContent::HAS_HYDROGEN:
    return std::string("HAS_HYDROGEN");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ForceFieldFamily input) {
  switch (input) {
  case ForceFieldFamily::BASIC:
    return std::string("BASIC");
  case ForceFieldFamily::AMBER:
    return std::string("AMBER");
  case ForceFieldFamily::CHARMM:
    return std::string("CHARMM");
  case ForceFieldFamily::OPENMM:
    return std::string("OPENMM");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const TorsionKind input) {
  switch (input) {
  case TorsionKind::PROPER:
    return std::string("PROPER");
  case TorsionKind::PROPER_NO_14:
    return std::string("PROPER_NO_14");
  case TorsionKind::IMPROPER:
    return std::string("IMPROPER");
  case TorsionKind::IMPROPER_NO_14:
    return std::string("IMPROPER_NO_14");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const ImplicitSolventModel input) {
  switch (input) {
  case ImplicitSolventModel::NONE:
    return std::string("no implicit solvent");
  case ImplicitSolventModel::HCT_GB:
    return std::string("Hawkins / Cramer / Truhlar Generalized Born");
  case ImplicitSolventModel::OBC_GB:
    return std::string("Onufriev / Bashford / Case Generalized Born (model I)");
  case ImplicitSolventModel::OBC_GB_II:
    return std::string("Onufriev / Bashford / Case Generalized Born (model II)");
  case ImplicitSolventModel::NECK_GB:
    return std::string("Mogan \"Neck\" Generalized Born (model I)");
  case ImplicitSolventModel::NECK_GB_II:
    return std::string("Mogan \"Neck\" Generalized Born (model II)");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const AtomicRadiusSet input) {
  switch (input) {
  case AtomicRadiusSet::NONE:
    return std::string("no atomic radii");
  case AtomicRadiusSet::BONDI:
    return std::string("Bondi radii");
  case AtomicRadiusSet::AMBER6:
    return std::string("Amber6 modified Bondi radii");
  case AtomicRadiusSet::MBONDI:
    return std::string("modified Bondi radii");
  case AtomicRadiusSet::MBONDI2:
    return std::string("amide H-modified Bondi (mBondi2) radii");
  case AtomicRadiusSet::MBONDI3:
    return std::string("Arg hydrogen and Asp / Glu oxygen modified mBondi2 (mBondi3) radii");
  case AtomicRadiusSet::PARSE:
    return std::string("Parse Poisson-Boltzmann radii");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const WaterModel input) {
  switch (input) {
  case WaterModel::NONE:
    return std::string("No water model");
  case WaterModel::UNKNOWN:
    return std::string("Unknown water model");
  case WaterModel::CHIMERA:
    return std::string("CHIMERIC water model (warning)");
  case WaterModel::MULTIPLE:
    return std::string("MULTIPLE water models (warning)");
  case WaterModel::MULTI_CHIMERA:
    return std::string("MULTIPLE water models with at least one CHIMERIC (warning)");
  case WaterModel::OPC3:
    return std::string("Optimized Point Charge, 3-Site");
  case WaterModel::SPC:
    return std::string("Simple Point Charge");
  case WaterModel::SPC_E:
    return std::string("Simple Point Charge, Extended");
  case WaterModel::SPC_EB:
    return std::string("Simple Point Charge, Extended-B");
  case WaterModel::SPC_HW:
    return std::string("Simple Point Charge, Heavy Water");
  case WaterModel::TIP3P:
    return std::string("Transferable Interaction Potential, 3-Site");
  case WaterModel::TIP3P_EW:
    return std::string("Transferable Interaction Potential, 3-Site, PME Adaptation F");
  case WaterModel::TIP3P_CHARMM:
    return std::string("Transferable Interaction Potential, 3-Site, CHARMM Adaptation");
  case WaterModel::SPC_FW:
    return std::string("Simple Point Charge, Flexible");
  case WaterModel::TIP3P_FW:
    return std::string("Transferable Interaction Potential, 3-Site, Flexible");
  case WaterModel::OPC:
    return std::string("Optimized Point Charge, 4-site");
  case WaterModel::TIP4P:
    return std::string("Transferable Interaction Potential, 4-Site");
  case WaterModel::TIP4P_EW:
    return std::string("Transferable Interaction Potential, 4-Site, PME Adaptation by IBM");
  case WaterModel::TIP4P_2005:
    return std::string("Transferable Interaction Potential, 4-Site, 2005 PME Adaptation");
  case WaterModel::TIP4P_ICE:
    return std::string("Transferable Interaction Potential, 4-Site, Ice Phase Diagram");
  case WaterModel::TIP4P_EPS:
    return std::string("Transferable Interaction Potential, 4-Site, 2016 Epsilon Adaptation");
  case WaterModel::TIP4P_D:
    return std::string("Transferable Interaction Potential, 4-Site, Dispersion Corrected");
  case WaterModel::TIP4P_2005F:
    return std::string("Transferable Interaction Potential, 4-Site, PME Adaptation, Flexible");
  case WaterModel::TIP5P:
    return std::string("Transferable Interaction Potential, 5-Site");
  case WaterModel::TIP5P_EW:
    return std::string("Transferable Interaction Potential, 5-Site, PME Adaptation");
  case WaterModel::TIP5P_2018:
    return std::string("Transferable Interaction Potential, 5-Site, 2018 Adaptation");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::string getEnumerationName(const VirtualSiteKind input) {
  switch (input) {
  case VirtualSiteKind::FLEX_2:
    return std::string("Flex-2 / FlexDis2");
  case VirtualSiteKind::FIXED_2:
    return std::string("FD-2 / FixedDis2");
  case VirtualSiteKind::FLEX_3:
    return std::string("Flex-3 / FlexDis3");
  case VirtualSiteKind::FIXED_3:
    return std::string("FD-3 / FixedDis3");
  case VirtualSiteKind::FAD_3:
    return std::string("FAD-3 / FixAnglDis");
  case VirtualSiteKind::OUT_3:
    return std::string("Out-3 / OutOfPlane");
  case VirtualSiteKind::FIXED_4:
    return std::string("FD-4 / FourPoint");
  case VirtualSiteKind::NONE:
    break;
  }
  return std::string("No detected frame type.");
}

//-------------------------------------------------------------------------------------------------
ImplicitSolventModel translateImplicitSolventModel(const int igb_val,
                                                   const ExceptionResponse policy) {
  switch (igb_val) {
  case 0:
  case 6:
    return ImplicitSolventModel::NONE;
  case 1:
    return ImplicitSolventModel::HCT_GB;
  case 2:
    return ImplicitSolventModel::OBC_GB;
  case 5:
    return ImplicitSolventModel::OBC_GB_II;
  case 7:
    return ImplicitSolventModel::NECK_GB;
  case 8:
    return ImplicitSolventModel::NECK_GB_II;
  default:
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Unrecognized implicit solvent model, igb = " + std::to_string(igb_val) + ".",
            "translateImplicitSolventModel");
    case ExceptionResponse::WARN:
      rtWarn("Unrecognized implicit solvent model, igb = " + std::to_string(igb_val) + ".  The "
             "model will be set to NONE instead.", "translateImplicitSolventModel");
      return ImplicitSolventModel::NONE;
    case ExceptionResponse::SILENT:
      return ImplicitSolventModel::NONE;
    }
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
AtomicRadiusSet translateAtomicRadiusSet(const std::string &pb_radii_in,
                                         const ExceptionResponse policy) {
  if (strcmpCased(pb_radii_in, std::string("bondi"), CaseSensitivity::NO)) {
    return AtomicRadiusSet::BONDI;
  }
  else if (strcmpCased(pb_radii_in, std::string("amber6"), CaseSensitivity::NO)) {
    return AtomicRadiusSet::AMBER6;
  }
  else if (strcmpCased(pb_radii_in, std::string("mbondi"), CaseSensitivity::NO)) {
    return AtomicRadiusSet::MBONDI;
  }
  else if (strcmpCased(pb_radii_in, std::string("mbondi2"), CaseSensitivity::NO)) {
    return AtomicRadiusSet::MBONDI2;
  }
  else if (strcmpCased(pb_radii_in, std::string("mbondi3"), CaseSensitivity::NO)) {
    return AtomicRadiusSet::MBONDI3;
  }
  else if (strcmpCased(pb_radii_in, std::string("parse"), CaseSensitivity::NO)) {
    return AtomicRadiusSet::PARSE;
  }
  else if (strcmpCased(pb_radii_in, std::string("none"), CaseSensitivity::NO)) {
    return AtomicRadiusSet::NONE;
  }
  else {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Unrecognized atomic radius set " + pb_radii_in + ".", "translatePBRadiiSet");
    case ExceptionResponse::WARN:
      rtWarn("Unrecognized atomic radius set " + pb_radii_in + ".  The radius set will be NONE.",
             "translatePBRadiiSet");
      return AtomicRadiusSet::NONE;
    case ExceptionResponse::SILENT:
      return AtomicRadiusSet::NONE;
    }
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
Citation getWaterModelCitation(const WaterModel wm) {
  Citation cite("Primary reference for " + getEnumerationName(wm));
  switch (wm) {
  case WaterModel::NONE:
  case WaterModel::UNKNOWN:
  case WaterModel::CHIMERA:
  case WaterModel::MULTIPLE:
  case WaterModel::MULTI_CHIMERA:
    cite.addJournal("");
    return cite;
  case WaterModel::OPC3:
    cite.addAuthors({"Saeed Izadi", "Alexey V. Onufriev"});
    cite.addTitle("Accuracy limit of rigid 3-point water models");
    cite.addJournal("J. Chem. Phys.");
    cite.addYear(2016);
    cite.addVolume(145);
    cite.addIssue(7);
    cite.addStartPage(074501);
    cite.addDoI("10.1063", "1.4960175");
    return cite;
  case WaterModel::SPC:
    cite.addAuthors({"H.J.C. Berendsen", "J.P.M. Postma", "W.F. van Gunsteren", "J. Hermans"});
    cite.addEditor("B. Pullman");
    cite.addTitle("Intermolecular Forces");
    cite.addPublisher("Springer, Dordrecht");
    cite.addYear(1981);
    cite.addStartPage(331);
    cite.addEndPage(342);
    cite.addDoI("10.1007", "978-94-015-7658-1_21");
    return cite;
  case WaterModel::SPC_E:
    cite.addAuthor("H.J.C. Berendsen");
    cite.addTitle("The missing term in effective pair potentials");
    cite.addJournal("J. Phys. Chem.");
    cite.addYear(1987);
    cite.addVolume(91);
    cite.addStartPage(6269);
    cite.addEndPage(6271);
    cite.addDoI("10.1021", "j100308a038");
    return cite;
  case WaterModel::SPC_EB:
    cite.addAuthors({"Kazuhiro Takemura", "Akio Kitao"});
    cite.addTitle("Water Model Tuning for Improved Reproduction of Rotational Diffusion and NMR "
                  "Spectral Density");
    cite.addJournal("J. Phys. Chem. B");
    cite.addYear(2012);
    cite.addVolume(116);
    cite.addIssue(22);
    cite.addStartPage(6279);
    cite.addEndPage(6287);
    cite.addDoI("10.1021", "jp301100g");
    return cite;
  case WaterModel::SPC_HW:
    cite.addAuthor("J. Raul Grigera");
    cite.addTitle("An effective pair potential for heavy water");
    cite.addJournal("J. Chem. Phys.");
    cite.addYear(2001);
    cite.addVolume(114);
    cite.addStartPage(8064);
    cite.addDoI("10.1063", "1.1359183");
    return cite;
  case WaterModel::TIP3P:
    cite.addAuthors({"William L. Jorgensen", "Jayaraman Chandrasekhar", "Jeffry D. Madura"});
    cite.addTitle("Comparison of simple potential functions for simulating liquid water");
    cite.addJournal("J. Chem. Phys.");
    cite.addYear(1983);
    cite.addVolume(79);
    cite.addStartPage(926);
    cite.addDoI("10.1063", "1.445869");
    return cite;
  case WaterModel::TIP3P_EW:
    cite.addAuthors({"Daniel J. Price", "Charles L. Brooks III"});
    cite.addTitle("A modified TIP3P water potential for simulation with Ewald summation");
    cite.addJournal("J. Chem. Phys.");
    cite.addYear(2004);
    cite.addVolume(121);
    cite.addStartPage(10096);
    cite.addDoI("10.1063", "1.1808117");
    return cite;
  case WaterModel::TIP3P_CHARMM:
    cite.addAuthors({"A.D. MacKerell", "D. Bashford", "M. Bellott", "R.L. Dunbrack",
                     "J.D. Evanseck", "M.J. Field", "S. Fischer", "J. Gao", "H. Guo", "S. Ha",
                     "D. Joseph-McCarthy", "L. Kuchnir", "K. Kuczera", "F.T. Lau", "C. Mattos",
                     "S. Michnick", "T. Ngo", "D.T. Nguyen", "B. Prodhom", "W.E. Reiher",
                     "B. Roux", "M. Schlenkrich", "J.C. Smith", "R. Stote", "J. Straub",
                     "M. Watanabe", "J. Wiorkiewicz-Kuczera", "D. Yin", "M. Karplus"});
    cite.addTitle("All-atom empirical potential for molecular modeling and dynamics studies of "
                  "proteins");
    cite.addJournal("J. Phys. Chem. B");
    cite.addYear(1998);
    cite.addVolume(102);
    cite.addIssue(18);
    cite.addStartPage(3586);
    cite.addEndPage(3616);
    cite.addDoI("10.1021", "jp973084f");
    return cite;
  case WaterModel::SPC_FW:
    cite.addAuthor("Yujie Wu");
    cite.addTitle("Flexible simple point-charge water model with improved liquid-state "
                  "properties");
    cite.addJournal("J. Chem. Phys.");
    cite.addYear(2006);
    cite.addVolume(124);
    cite.addStartPage(024503);
    cite.addDoI("10.1063", "1.2136877");
    return cite;
  case WaterModel::TIP3P_FW:
    cite.addAuthors({"Udo W. Schmitt", "Gregory A. Voth"});
    cite.addTitle("The computer simulation of proton transport in water");
    cite.addJournal("J. Chem. Phys.");
    cite.addYear(1999);
    cite.addVolume(111);
    cite.addStartPage(9361);
    cite.addDoI("10.1063", "1.480032");
    return cite;
  case WaterModel::OPC:
    cite.addAuthors({"Saeed Izadi", "Ramu Anandakrishnan", "Alexey V. Onufriev"});
    cite.addTitle("Building Water Models: A Different Approach");
    cite.addJournal("J. Phys. Chem. Lett.");
    cite.addYear(2014);
    cite.addVolume(5);
    cite.addIssue(21);
    cite.addStartPage(3863);
    cite.addEndPage(3871);
    cite.addDoI("10.1021", "jz501780a");
    return cite;
  case WaterModel::TIP4P:
    cite.addAuthors({"William L. Jorgensen", "Jeffry D. Madura"});
    cite.addTitle("Temperature and size dependence for Monte Carlo simulations of TIP4P water");
    cite.addJournal("Mol. Phys.");
    cite.addYear(1985);
    cite.addVolume(56);
    cite.addIssue(6);
    cite.addStartPage(1381);
    cite.addEndPage(1392);
    cite.addDoI("10.1080", "00268978500103111");
    return cite;
  case WaterModel::TIP4P_EW:
    cite.addAuthors({"Hans W. Horn", "William C. Swope", "Jed W. Pitera", "Jeffry D. Madura",
                     "Thomas J. Dick", "Greg L. Hura", "Teresa Head-Gordon"});
    cite.addTitle("Development of an improved four-site water model for biomolecular simulations: "
                  "TIP4P-Ew");
    cite.addJournal("J. Chem. Phys.");
    cite.addYear(2004);
    cite.addVolume(120);
    cite.addIssue(20);
    cite.addStartPage(9665);
    cite.addEndPage(9678);
    cite.addDoI("10.1063", "1.1683075");
    return cite;
  case WaterModel::TIP4P_2005:
    cite.addAuthors({"J. L. F. Abascal", "C. Vega"});
    cite.addTitle("A general purpose model for the condensed phases of water: TIP4P/2005");
    cite.addJournal("J. Chem. Phys.");
    cite.addYear(2005);
    cite.addVolume(123);
    cite.addStartPage(234505);
    cite.addDoI("10.1063", "1.2121687");
    return cite;
  case WaterModel::TIP4P_ICE:
    cite.addAuthors({"J. L. F. Abascal", "E. Sanz", "R. Garcia Fernandez", "C. Vega"});
    cite.addTitle("A potential model for the study of ices and amorphous water: TIP4P/Ice");
    cite.addJournal("J. Chem. Phys.");
    cite.addYear(2005);
    cite.addVolume(122);
    cite.addStartPage(234511);
    cite.addDoI("10.1063", "1.1931662");
    return cite;
  case WaterModel::TIP4P_EPS:
    cite.addAuthors({"Raul Fuentes-Azcatl", "Marcia C. Barbosa"});
    cite.addTitle("Thermodynamic and dynamic anomalous behavior in the TIP4P/eps water model");
    cite.addJournal("Physica A.");
    cite.addYear(2016);
    cite.addVolume(444);
    cite.addStartPage(86);
    cite.addEndPage(94);
    cite.addDoI("10.1016", "j.physa.2015.10.027");
    return cite;
  case WaterModel::TIP4P_D:
    cite.addAuthors({"Stefano Piana", "Alexander G. Donchev", "Paul Robustelli", "David E. Shaw"});
    cite.addTitle("Water Dispersion Interactions Strongly Influence Simulated Structural "
                  "Properties of Disordered Protein States");
    cite.addJournal("J. Phys. Chem. B");
    cite.addYear(2015);
    cite.addVolume(119);
    cite.addIssue(16);
    cite.addStartPage(5113);
    cite.addEndPage(5123);
    cite.addDoI("10.1021", "jp508971m");
    return cite;
  case WaterModel::TIP4P_2005F:
    cite.addAuthors({"Miguel A. Gonzalez", "Jose L. F. Abascal"});
    cite.addTitle("A flexible model for water based on TIP4P/2005");
    cite.addJournal("J. Chem. Phys.");
    cite.addYear(2011);
    cite.addVolume(135);
    cite.addStartPage(224516);
    cite.addDoI("10.1063", "1.3663219");
    return cite;
  case WaterModel::TIP5P:
    cite.addAuthor("Michael W. Mahoney");
    cite.addTitle("A five-site model for liquid water and the reproduction of the density anomaly "
                  "by rigid, nonpolarizable potential functions");
    cite.addJournal("J. Chem. Phys.");
    cite.addYear(2000);
    cite.addVolume(112);
    cite.addStartPage(8910);
    cite.addDoI("10.1063", "1.481505");
    return cite;
  case WaterModel::TIP5P_EW:
    cite.addAuthor("Steven E. Rick");
    cite.addTitle("A reoptimization of the five-site water potential (TIP5P) for use with Ewald "
                  "sums");
    cite.addJournal("J. Chem. Phys.");
    cite.addYear(2004);
    cite.addVolume(120);
    cite.addStartPage(6085);
    cite.addDoI("10.1063", "1.1652434");
    return cite;
  case WaterModel::TIP5P_2018:
    cite.addAuthors({"Yuriy Khalak", "Bjorn Baumeier", "Mikko Karttunen"});
    cite.addTitle("Improved general-purpose five-point model for water: TIP5P/2018");
    cite.addJournal("J. Chem. Phys.");
    cite.addYear(2018);
    cite.addVolume(149);
    cite.addStartPage(224507);
    cite.addDoI("10.1063", "1.5070137");
    return cite;
  }
  __builtin_unreachable();
}

} // namespace topology
} // namespace stormm
