#include <cmath>
#include "copyright.h"
#include "Chemistry/chemistry_enumerators.h"
#include "Chemistry/znumber.h"
#include "FileManagement/file_listing.h"
#include "Topology/atomgraph_enumerators.h"
#include "Trajectory/write_annotated_frame.h"
#include "mdlmol.h"

namespace stormm {
namespace structure {

using chemistry::ChemicalFeaturesReader;
using chemistry::ChiralOrientation;
using chemistry::zNumberToNaturalMass;
using chemistry::zNumberToSymbol;
using diskutil::getBaseName;
using parse::char4ToString;
using parse::NumberFormat;
using parse::readIntegerValue;
using parse::readRealValue;
using parse::separateText;
using parse::stringToChar4;
using parse::strncmpCased;
using parse::TextFileReader;
using parse::TextOrigin;
using parse::verifyContents;
using data_types::operator==;
using topology::TorsionKind;
using topology::ValenceKit;
using trajectory::writeFrame;

//-------------------------------------------------------------------------------------------------
MdlMol::MdlMol(const ExceptionResponse policy_in):
    policy{policy_in}, version_no{MdlMolVersion::V2000}, atom_count{0}, bond_count{0},
    list_count{0}, stext_entry_count{0}, properties_count{0}, sgroup_count{0}, constraint_count{0},
    chirality{MolObjChirality::ACHIRAL}, registry_number{-1}, data_item_count{0},
    property_formal_charges{false}, property_radicals{false}, property_isotopes{false},
    property_element_lists{false}, coordinates{}, atomic_symbols{}, atomic_numbers{},
    formal_charges{}, isotopic_shifts{}, parities{}, implicit_hydrogens{}, stereo_considerations{},
    valence_connections{}, atom_atom_mapping_count{}, orientation_stability{}, bonds{},
    element_lists{}, stext_entries{}, properties{}, title{""}, software_details{""},
    general_comment{""}
{}

//-------------------------------------------------------------------------------------------------
MdlMol::MdlMol(const TextFile &tf, const int line_start, const int line_end_in,
               const CaseSensitivity capitalization, const ExceptionResponse policy_in,
               const ModificationPolicy dimod_policy, const ExceptionResponse dimod_notify):
  MdlMol(policy_in)
{
  const TextFileReader tfr = tf.data();
  
  // Default line end of -1 indicates reading to the end of the file.  Otherwise, identify the
  // end of the formatting ("M  END").
  const int mdl_section_end = getMdlMolSectionEnd(tfr, line_start, line_end_in);
  const int sd_compound_end = getCompoundSectionEnd(tfr, line_start, line_end_in);
  
  // The range of data now extends from line_start to mdl_section_end.  Sift through that
  // information for a V2000 or V3000 specification.  This should be found on the fourth line.
  version_no = findMolObjVersion(tf, line_start + 3);
  
  // Begin by reading the molecule name (title), generating software details, and any general
  // comment (always three and only three distinct lines, even if left blank).
  if (line_start + 2 < tfr.line_count) {
    title            = tf.extractString(line_start);
    software_details = tf.extractString(line_start + 1);
    general_comment  = tf.extractString(line_start + 2);
  }
  
  // Read the counts line
  switch (version_no) {
  case MdlMolVersion::V2000:
    {
      const int counts_line_idx = line_start + 3;
      const char* counts_line_ptr = &tfr.text[tfr.line_limits[line_start + 3]];
      atom_count   = readIntegerValue(counts_line_ptr, 0, 3);
      bond_count   = readIntegerValue(counts_line_ptr, 3, 3);
      if (verifyContents(tf, counts_line_idx, 6, 3, NumberFormat::INTEGER)) {
        list_count   = readIntegerValue(counts_line_ptr, 6, 3);
      }
      if (verifyContents(tf, counts_line_idx, 12, 3, NumberFormat::INTEGER)) {
        const int chir_num = readIntegerValue(counts_line_ptr, 12, 3);
        if (chir_num == 0) {
          chirality = MolObjChirality::ACHIRAL;
        }
        else if (chir_num == 1) {
          chirality = MolObjChirality::CHIRAL;
        }
        else {
          rtErr("Invalid chirality setting detected at line " + std::to_string(counts_line_idx) +
                " in .sdf or MDL MOL file " + getBaseName(tf.getFileName()) + ".", "MdlMol");
        }
      }
      if (verifyContents(tf, counts_line_idx, 15, 3, NumberFormat::INTEGER)) {
        stext_entry_count = readIntegerValue(counts_line_ptr, 15, 3);
      }
      if (verifyContents(tf, counts_line_idx, 30, 3, NumberFormat::INTEGER)) {
        properties_count  = readIntegerValue(counts_line_ptr, 30, 3);
      }
      
      // Validation of the atom counts line
      if (atom_count <= 0 || bond_count < 0) {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("A V2000 MOL format entry in file " + getBaseName(tf.getFileName()) + " at line " +
                std::to_string(line_start) + " contains invalid numbers of atoms (" +
                std::to_string(atom_count) + ") and/or bonds (" + std::to_string(bond_count) +
                ").", "MdlMol");
        case ExceptionResponse::WARN:
          rtWarn("A V2000 MOL format entry in file " + getBaseName(tf.getFileName()) + " at " +
                 "line " + std::to_string(line_start) + " contains invalid numbers of atoms (" +
                 std::to_string(atom_count) + ") and/or bonds (" + std::to_string(bond_count) +
                 ").  This will become a blank entry", "MdlMol");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
        atom_count = 0;
        bond_count = 0;
        list_count = 0;
        chirality = MolObjChirality::ACHIRAL;
        stext_entry_count = 0;
        properties_count = 0;
        return;
      }
    }
    break;
  case MdlMolVersion::V3000:
    break;
  case MdlMolVersion::UNKNOWN:

    // This case would be unprocessable.
    rtErr("No valid MDL MOL version was detected.  Parsing in " + getBaseName(tf.getFileName()) +
          " cannot proceed.");
  }

  // Allocate space for information to be read
  allocate();
    
  // Read the atom block
  int iatm = 0;
  switch (version_no) {
  case MdlMolVersion::V2000:
    for (int i = line_start + 4; i < line_start + 4 + atom_count; i++) {
      const char* atom_line_ptr = &tfr.text[tfr.line_limits[i]];
      coordinates[iatm] = { readRealValue(atom_line_ptr,  0, 10),
                            readRealValue(atom_line_ptr, 10, 10),
                            readRealValue(atom_line_ptr, 20, 10) };
      if (verifyContents(tf, i, 31, 3, NumberFormat::CHAR4)) {
        atomic_symbols[iatm] = tf.extractChar4(i, 31, 3);
      }
      if (verifyContents(tf, i, 34, 2, NumberFormat::INTEGER)) {
        isotopic_shifts[iatm] = readIntegerValue(atom_line_ptr, 34, 2);
        if (isotopic_shifts[iatm] > 4 || isotopic_shifts[iatm] < -3) {
          rtErr("A V2000 MOL format entry should not describe an isotopic shift outside the range "
                "[-3, 4].  Shift found: " + std::to_string(isotopic_shifts[iatm]) +
                ".  Title of entry: \"" + title + "\".", "MdlMol");
        }
      }
      
      // Standard Template Library vector<bool> works differently from other vectors.  Set its
      // contents in a different manner.
      if (verifyContents(tf, i, 36, 3, NumberFormat::INTEGER)) {
        interpretFormalCharge(readIntegerValue(atom_line_ptr, 36, 3), iatm);
      }
      if (verifyContents(tf, i, 39, 3, NumberFormat::INTEGER)) {
        parities[iatm] = interpretStereoParity(readIntegerValue(atom_line_ptr, 39, 3));
      }
      if (verifyContents(tf, i, 42, 3, NumberFormat::INTEGER)) {
        interpretImplicitHydrogenContent(readIntegerValue(atom_line_ptr, 42, 3), iatm);
      }
      if (verifyContents(tf, i, 45, 3, NumberFormat::INTEGER)) {
        stereo_considerations[iatm] =
          interpretBooleanValue(readIntegerValue(atom_line_ptr, 45, 3),
                                "interpreting stereochemical considerations");
      }
      if (verifyContents(tf, i, 48, 3, NumberFormat::INTEGER)) {
        valence_connections[iatm] = interpretValenceNumber(readIntegerValue(atom_line_ptr, 48, 3));
      }
      if (verifyContents(tf, i, 51, 3, NumberFormat::INTEGER)) {
        if (readIntegerValue(atom_line_ptr, 51, 3) == 1 && implicit_hydrogens[iatm] > 0) {
          rtErr("The H0 designation, indicating that implicit hydrogens are not allowed on atom " +
                std::to_string(iatm + 1) + " of MDL MOL entry \"" +  title + "\", is present but "
                "the number of implicit hydrogens has also been indicated as " +
                std::to_string(implicit_hydrogens[iatm]) + ".", "MdlMol");
        }
      }
      if (verifyContents(tf, i, 60, 3, NumberFormat::INTEGER)) {
        atom_atom_mapping_count[iatm] = readIntegerValue(atom_line_ptr, 60, 3);
      }
      if (verifyContents(tf, i, 63, 3, NumberFormat::INTEGER)) {
        orientation_stability[iatm] =
          interpretStereoStability(readIntegerValue(atom_line_ptr, 63, 3));
      }
      if (verifyContents(tf, i, 66, 3, NumberFormat::INTEGER)) {
        exact_change_enforced[iatm] = interpretBooleanValue(readIntegerValue(atom_line_ptr, 66, 3),
                                                            "interpreting exact change flag");
      }
      iatm++;
    }
    break;
  case MdlMolVersion::V3000:
    break;
  case MdlMolVersion::UNKNOWN:
    break;
  }

  // Read the bonds block
  switch (version_no) {
  case MdlMolVersion::V2000:
    {
      const int bond_line_start = line_start + 4 + atom_count;
      for (int pos = 0; pos < bond_count; pos++) {
        bonds.emplace_back(tf, bond_line_start + pos, title);
      }
    }
    break;
  case MdlMolVersion::V3000:
    break;
  case MdlMolVersion::UNKNOWN:
    break;
  }

  // Read the atom lists (this information is superceded by the presence of "M  ALS" properties)
  switch (version_no) {
  case MdlMolVersion::V2000:
    {
      // Scan for atom lists imbedded in properties
      const std::string als_tag("M  ALS");
      for (int pos = line_start; pos < mdl_section_end; pos++) {
        if (tfr.line_lengths[pos] >= 6 && strncmpCased(tf.getLinePointer(pos), als_tag)) {
          element_lists.emplace_back(tf, pos, title);
        }
      }
      if (element_lists.size() == 0LLU) {
        const int alst_line_start = line_start + 4 + atom_count + bond_count;
        for (int pos = 0; pos < list_count; pos++) {
          element_lists.emplace_back(tf, alst_line_start + pos, title);          
        }
      }
    }
    break;
  case MdlMolVersion::V3000:
    break;
  case MdlMolVersion::UNKNOWN:
    break;
  }

  // Read the stext entries
  switch (version_no) {
  case MdlMolVersion::V2000:
    {
      const int stxt_line_start = line_start + 4 + atom_count + bond_count + list_count;
      for (int pos = 0; pos < stext_entry_count; pos += 2) {
      }
    }
    break;
  case MdlMolVersion::V3000:
    break;
  case MdlMolVersion::UNKNOWN:
    break;
  }

  // Read various properties
  switch (version_no) {
  case MdlMolVersion::V2000:
    {
      const int prop_line_start = line_start + 4 + atom_count + bond_count + list_count +
                                  (2 * stext_entry_count);
      for (int pos = prop_line_start; pos < mdl_section_end; pos++) {
        int adv_pos;
        properties.emplace_back(tf, pos, &adv_pos, title);
      }
      
      // Update the properties count
      properties_count = properties.size();
    }
    break;
  case MdlMolVersion::V3000:
    break;
  case MdlMolVersion::UNKNOWN:
    break;
  }
    
  // Read the data items
  for (int pos = mdl_section_end; pos < sd_compound_end; pos++) {
    if (tf.getLineLength(pos) >= 2 && tf.getChar(tf.getLineLimits(pos)) == '>') {
      int adv_pos;
      data_items.emplace_back(tf, pos, &adv_pos, sd_compound_end, title, dimod_policy,
                              dimod_notify);
      compareExternalRegistryNumbers(data_items.back().getExternalRegistryNumber());
      pos = adv_pos;
    }
  }
  data_item_count = data_items.size();
  if (data_item_count == 0 && sd_compound_end - mdl_section_end > 1) {
    rtErr("If there are no data items, the compound section must terminate immediately after the "
          "MDL MOL format section.  File " + getBaseName(tf.getFileName()) + " violates SD file "
          "conventions at lines " + std::to_string(mdl_section_end) + " - " +
          std::to_string(sd_compound_end));
  }
  
  // Update the atom attributes based on properties data.  This provides V3000 functionality and
  // backwards compatibility for the V2000 format.
  updateV2kAtomAttributes();
  
  // Make some basic inferences.  The version number is irrelevant by now: all information has been
  // converted into version-agnostic internal representations and the version number is kept only
  // as a footnote for reference when re-printing the file later.
  const std::vector<int> tmp_znum = symbolToZNumber(atomic_symbols, capitalization, policy);
  int nvs = 0;
  for (int i = 0; i < atom_count; i++) {
    atomic_numbers[i] = tmp_znum[i];
    nvs += (atomic_numbers[i] == 0);
  }
  if (nvs > 0) {
    
  }

  // Add implicit hydrogens.  This may re-allocate the data arrays and extend the bonding patterns.
  hydrogenate();
}

//-------------------------------------------------------------------------------------------------
MdlMol::MdlMol(const std::string &filename, const ExceptionResponse policy_in,
               const ModificationPolicy dimod_policy, const ExceptionResponse dimod_notify):
  MdlMol(TextFile(filename), 0, -1, CaseSensitivity::YES, policy_in, dimod_policy, dimod_notify)
{}

//-------------------------------------------------------------------------------------------------
MdlMol::MdlMol(const char* filename, const ExceptionResponse policy_in,
               const ModificationPolicy dimod_policy, const ExceptionResponse dimod_notify):
    MdlMol(std::string(filename), policy_in, dimod_policy, dimod_notify)
{}

//-------------------------------------------------------------------------------------------------
MdlMol::MdlMol(const ChemicalFeatures *chemfe, const CoordinateFrameReader cfr,
               const int molecule_index) :
    MdlMol(chemfe, cfr.xcrd, cfr.ycrd, cfr.zcrd, molecule_index)
{}

//-------------------------------------------------------------------------------------------------
MdlMol::MdlMol(const ChemicalFeatures *chemfe, const CoordinateFrame *cf,
               const int molecule_index) :
    MdlMol(chemfe, cf->data(), molecule_index)
{}

//-------------------------------------------------------------------------------------------------
MdlMol::MdlMol(const ChemicalFeatures &chemfe, const CoordinateFrame &cf,
               const int molecule_index) :
    MdlMol(chemfe.getSelfPointer(), cf.data(), molecule_index)
{}

//-------------------------------------------------------------------------------------------------
MdlMol::MdlMol(const ChemicalFeatures *chemfe, const PhaseSpaceReader psr,
               const int molecule_index) :
    MdlMol(chemfe, psr.xcrd, psr.ycrd, psr.zcrd, molecule_index)
{}

//-------------------------------------------------------------------------------------------------
MdlMol::MdlMol(const ChemicalFeatures *chemfe, const PhaseSpace *ps,
               const int molecule_index) :
    MdlMol(chemfe, ps->data(), molecule_index)
{}

//-------------------------------------------------------------------------------------------------
MdlMol::MdlMol(const ChemicalFeatures &chemfe, const PhaseSpace &ps,
               const int molecule_index) :
    MdlMol(chemfe.getSelfPointer(), ps.data(), molecule_index)
{}

//-------------------------------------------------------------------------------------------------
MdlMol::MdlMol(const ChemicalFeatures *chemfe, const PsSynthesisReader poly_psr,
               const int system_index, const int molecule_index) :
    MdlMol(ExceptionResponse::SILENT)
{
  // Check that the system index and system size are valid
  if (system_index < 0 || system_index >= poly_psr.system_count) {
    rtErr("System index " + std::to_string(system_index) + " is invalid for a synthesis of " +
          std::to_string(poly_psr.system_count) + " systems.");
  }
  const AtomGraph *ag = chemfe->getTopologyPointer();
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  if (cdk.natom != poly_psr.atom_counts[system_index]) {
    rtErr("System " + std::to_string(system_index) + " contains " +
          std::to_string(poly_psr.atom_counts[system_index]) + " atoms, whereas the topology "
          "expects " + std::to_string(cdk.natom) + ".", "MdlMol");
  }
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  allocate(cdk, nbk, molecule_index);

  // Fill in details from the topology
  transferTopologicalDetails(chemfe, molecule_index);
  
  // Fill in the coordinates
  const size_t init_index = poly_psr.atom_starts[system_index];
  impartCoordinates(&poly_psr.xcrd[init_index], &poly_psr.ycrd[init_index],
                    &poly_psr.zcrd[init_index], &poly_psr.xcrd_ovrf[init_index],
                    &poly_psr.ycrd_ovrf[init_index], &poly_psr.zcrd_ovrf[init_index],
                    poly_psr.inv_gpos_scale, cdk, molecule_index);
}

//-------------------------------------------------------------------------------------------------
MdlMol::MdlMol(const ChemicalFeatures *chemfe, const PhaseSpaceSynthesis *poly_ps,
               const int system_index, const int molecule_index) :
    MdlMol(chemfe, poly_ps->data(), system_index, molecule_index)
{}

//-------------------------------------------------------------------------------------------------
MdlMol::MdlMol(const ChemicalFeatures &chemfe, const PhaseSpaceSynthesis &poly_ps,
               const int system_index, const int molecule_index) :
    MdlMol(chemfe.getSelfPointer(), poly_ps.data(), system_index, molecule_index)
{}

//-------------------------------------------------------------------------------------------------
const std::string& MdlMol::getTitle() const {
  return title;
}

//-------------------------------------------------------------------------------------------------
int MdlMol::getAtomCount() const {
  return atom_count;
}

//-------------------------------------------------------------------------------------------------
int MdlMol::getBondCount() const {
  return bond_count;
}

//-------------------------------------------------------------------------------------------------
int MdlMol::getPropertiesCount() const {
  return properties_count;
}

//-------------------------------------------------------------------------------------------------
int MdlMol::getDataItemCount() const {
  return data_item_count;
}

//-------------------------------------------------------------------------------------------------
double3 MdlMol::getCoordinates(const int index) const {
  return coordinates[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<double3>& MdlMol::getCoordinates() const {
  return coordinates;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> MdlMol::getCoordinates(const CartesianDimension dim) const {
  std::vector<double> result(atom_count);
  for (int i = 0; i < atom_count; i++) {
    switch (dim) {
    case CartesianDimension::X:
      result[i] = coordinates[i].x;
      break;
    case CartesianDimension::Y:
      result[i] = coordinates[i].y;
      break;
    case CartesianDimension::Z:
      result[i] = coordinates[i].z;
      break;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const MdlMol* MdlMol::getSelfPointer() const {
  return this;
}

//-------------------------------------------------------------------------------------------------
PhaseSpace MdlMol::exportPhaseSpace() const {
  PhaseSpace result(atom_count);
  PhaseSpaceWriter rw = result.data();
  for (int i = 0; i < atom_count; i++) {
    rw.xcrd[i] = coordinates[i].x;
    rw.ycrd[i] = coordinates[i].y;
    rw.zcrd[i] = coordinates[i].z;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame MdlMol::exportCoordinateFrame() const {
  CoordinateFrame result(atom_count);
  CoordinateFrameWriter rw = result.data();
  for (int i = 0; i < atom_count; i++) {
    rw.xcrd[i] = coordinates[i].x;
    rw.ycrd[i] = coordinates[i].y;
    rw.zcrd[i] = coordinates[i].z;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
char4 MdlMol::getAtomSymbol(const int index) const {
  validateAtomIndex(index, "getAtomSymbol");
  return atomic_symbols[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<char4>& MdlMol::getAtomSymbols() const {
  return atomic_symbols;
}

//-------------------------------------------------------------------------------------------------
int MdlMol::getAtomicNumber(const int index) const {
  validateAtomIndex(index, "getAtomicNumber");
  return atomic_numbers[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& MdlMol::getAtomicNumbers() const {
  return atomic_numbers;
}

//-------------------------------------------------------------------------------------------------
int MdlMol::getFormalCharge(const int index) const {
  validateAtomIndex(index, "getFormalCharge");
  return formal_charges[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<int>& MdlMol::getFormalCharges() const {
  return formal_charges;
}

//-------------------------------------------------------------------------------------------------
void MdlMol::impartCoordinates(const llint* xcrd, const llint* ycrd, const llint* zcrd,
                               const int* xcrd_ovrf, const int* ycrd_ovrf, const int* zcrd_ovrf,
                               const double scale_factor) {
  for (int i = 0; i < atom_count; i++) {
    coordinates[i].x = hostInt95ToDouble(xcrd[i], xcrd_ovrf[i]) * scale_factor;
    coordinates[i].y = hostInt95ToDouble(ycrd[i], ycrd_ovrf[i]) * scale_factor;
    coordinates[i].z = hostInt95ToDouble(zcrd[i], zcrd_ovrf[i]) * scale_factor;
  }
}

//-------------------------------------------------------------------------------------------------
void MdlMol::impartCoordinates(const llint* xcrd, const llint* ycrd, const llint* zcrd,
                               const int* xcrd_ovrf, const int* ycrd_ovrf, const int* zcrd_ovrf,
                               const double scale_factor, const ChemicalDetailsKit &cdk,
                               const int molecule_index) {
  if (molecule_index < 0 || molecule_index >= cdk.nmol) {
    rtErr("Molecule index " + std::to_string(molecule_index) + " is invalid for a topology with " +
          std::to_string(cdk.nmol) + " molecules.", "MdlMol", "impartCoordinates");
  }
  const int llim = cdk.mol_limits[molecule_index];
  const int hlim = cdk.mol_limits[molecule_index + 1];
  if (hlim - llim != atom_count) {
    rtErr("Molecule " + std::to_string(molecule_index) + " with " + std::to_string(hlim - llim) +
          " atoms is incompatible with an MDL MOL object of " + std::to_string(atom_count) +
          " atoms.", "MdlMol", "impartCoordinates");
  }
  int atomcon = 0;
  for (int i = llim; i < hlim; i++) {
    const int iatom = cdk.mol_contents[i];
    if (cdk.z_numbers[iatom] > 0) {
      if (atomcon >= atom_count) {
        rtErr("Data from integral types is for fixed-precision representations and must be "
              "accompanied by a scaling factor to take the values back into real, internal units.",
              "MdlMol");
      }
      coordinates[atomcon].x = hostInt95ToDouble(xcrd[iatom], xcrd_ovrf[iatom]) * scale_factor;
      coordinates[atomcon].y = hostInt95ToDouble(ycrd[iatom], ycrd_ovrf[iatom]) * scale_factor;
      coordinates[atomcon].z = hostInt95ToDouble(zcrd[iatom], zcrd_ovrf[iatom]) * scale_factor;
      atomcon++;
    }
  }
}

//-------------------------------------------------------------------------------------------------
void MdlMol::impartCoordinates(const PhaseSpaceReader &psr) {
  checkAtomCount(psr.natom);
  impartCoordinates(psr.xcrd, psr.ycrd, psr.zcrd, 1.0);
}

//-------------------------------------------------------------------------------------------------
void MdlMol::impartCoordinates(const PhaseSpace *ps, const CoordinateCycle orientation,
                               const HybridTargetLevel tier) {
  impartCoordinates(ps->data(orientation, tier));
}

//-------------------------------------------------------------------------------------------------
void MdlMol::impartCoordinates(const PhaseSpace *ps, const HybridTargetLevel tier) {
  impartCoordinates(ps, ps->getCyclePosition(), tier);
}

//-------------------------------------------------------------------------------------------------
void MdlMol::impartCoordinates(const PhaseSpace &ps, const CoordinateCycle orientation,
                               const HybridTargetLevel tier) {
  impartCoordinates(ps.data(orientation, tier));
}

//-------------------------------------------------------------------------------------------------
void MdlMol::impartCoordinates(const PhaseSpace &ps, const HybridTargetLevel tier) {
  impartCoordinates(ps, ps.getCyclePosition(), tier);
}

//-------------------------------------------------------------------------------------------------
void MdlMol::impartCoordinates(const CoordinateFrameReader &cfr) {
  checkAtomCount(cfr.natom);
  impartCoordinates(cfr.xcrd, cfr.ycrd, cfr.zcrd, 1.0);
}

//-------------------------------------------------------------------------------------------------
void MdlMol::impartCoordinates(const CoordinateFrameWriter &cfw) {
  checkAtomCount(cfw.natom);
  impartCoordinates(cfw.xcrd, cfw.ycrd, cfw.zcrd, 1.0);
}

//-------------------------------------------------------------------------------------------------
void MdlMol::impartCoordinates(const CoordinateFrame *cf, const HybridTargetLevel tier) {
  impartCoordinates(cf->data());
}

//-------------------------------------------------------------------------------------------------
void MdlMol::impartCoordinates(const CoordinateFrame &cf, const HybridTargetLevel tier) {
  impartCoordinates(cf.data());
}

//-------------------------------------------------------------------------------------------------
void MdlMol::addProperty(const MdlMolProperty &mmprop) {
  bool update_needed = false;
  properties.push_back(mmprop);
  switch (mmprop.getKind()) {
  case MdlMolPropertyKind::CHARGE:
    property_formal_charges = true;
    update_needed = true;
    break;
  case MdlMolPropertyKind::RADICAL:
    property_radicals = true;
    update_needed = true;
    break;
  case MdlMolPropertyKind::ISOTOPE:
    property_isotopes = true;
    update_needed = true;
    break;
  case MdlMolPropertyKind::ATOM_LIST:
    property_element_lists = true;
    update_needed = true;
    break;
  case MdlMolPropertyKind::ATOM_ALIAS:
  case MdlMolPropertyKind::ATOM_VALUE:
  case MdlMolPropertyKind::GROUP_ABBREVIATION:
  case MdlMolPropertyKind::SGROUP_SUBSCRIPT:
  case MdlMolPropertyKind::SGROUP_BOND_VECTOR:
  case MdlMolPropertyKind::SGROUP_FIELD:
  case MdlMolPropertyKind::SGROUP_DISPLAY:
  case MdlMolPropertyKind::SGROUP_DATA:
  case MdlMolPropertyKind::SPATIAL_FEATURE:
  case MdlMolPropertyKind::PHANTOM_ATOM:
  case MdlMolPropertyKind::SGROUP_CLASS:
  case MdlMolPropertyKind::LARGE_REGNO:
  case MdlMolPropertyKind::SGROUP_EXPANSION:
  case MdlMolPropertyKind::SGROUP_ATOM_LIST:
  case MdlMolPropertyKind::SGROUP_BOND_LIST:
  case MdlMolPropertyKind::MG_PARENT_ATOM_LIST:
  case MdlMolPropertyKind::RING_BOND_COUNT:
  case MdlMolPropertyKind::SUBSTITUTION_COUNT:
  case MdlMolPropertyKind::UNSATURATED_COUNT:
  case MdlMolPropertyKind::RGROUP_LABEL_LOCATION:
  case MdlMolPropertyKind::SGROUP_TYPE:
  case MdlMolPropertyKind::SGROUP_SUBTYPE:
  case MdlMolPropertyKind::SGROUP_LABELS:
  case MdlMolPropertyKind::SGROUP_CONNECTIVITY:
  case MdlMolPropertyKind::SGROUP_HIERARCHY:
  case MdlMolPropertyKind::SGROUP_COMP_NUMBER:
  case MdlMolPropertyKind::SGROUP_BRACKET_STYLE:
  case MdlMolPropertyKind::SGROUP_CORRESPONDENCE:
  case MdlMolPropertyKind::SGROUP_ATTACH_POINT:
  case MdlMolPropertyKind::LINK_ATOM:
  case MdlMolPropertyKind::SGROUP_DISPLAY_INFO:
  case MdlMolPropertyKind::ATTACHMENT_POINT:
  case MdlMolPropertyKind::ATTACHMENT_ORDER:
  case MdlMolPropertyKind::RGROUP_LOGIC:
  case MdlMolPropertyKind::SKIP:
  case MdlMolPropertyKind::NONE:
    break;
  }
  if (update_needed) {
    updateV2kAtomAttributes();
  }
  properties_count += 1;
}

//-------------------------------------------------------------------------------------------------
void MdlMol::addProperty(const MdlMolProperty *mmprop) {
  addProperty(*mmprop);
}

//-------------------------------------------------------------------------------------------------
void MdlMol::addDataItem(const MdlMolDataRequest &ask, const AtomGraph &ag,
                         const RestraintApparatus &ra) {

  // If the text of the data item can be prepared in advance, do so.
  std::vector<std::string> di_body;
  switch (ask.getKind()) {
  case DataRequestKind::STATE_VARIABLE:
    
    // The text will be a single number, supplied by an external function or re-evaluated.
    break;
  case DataRequestKind::ATOM_INFLUENCES:

    // The text will be constructed during the energy re-evaluation when the data item is printed.
    break;
  case DataRequestKind::TOPOLOGY_PARAMETER:
    {
      // Search for all instances of a particular parameter type affecting the stated atom types
      // in the given order (or the reverse order).  Each instance of the parameter in the topology
      // will become one data line of the data item.  If no instances are found, a single data line
      // is entered to indicate this fact.  The first string in the body will be reserved to hold
      // a message indicating however many instances were found.
      int nfound = 0;
      di_body.push_back("");
      char buffer_data[256];
      const ValenceKit<double> vk  = ag.getDoublePrecisionValenceKit();
      const ChemicalDetailsKit cdk = ag.getChemicalDetailsKit();
      const std::vector<char4> atom_types = ask.getAtomTypes();
      switch (ask.getValenceParameter()) {
      case StateVariable::BOND:
        {
          const char4 i_tref = atom_types[0];
          const char4 j_tref = atom_types[1];
          for (int pos = 0; pos < vk.nbond; pos++) {
            const int i_atom = vk.bond_i_atoms[pos];
            const int j_atom = vk.bond_j_atoms[pos];
            const char4 i_type = cdk.atom_types[i_atom];
            const char4 j_type = cdk.atom_types[j_atom];
            if ((i_type == i_tref && j_type == j_tref) || (i_type == j_tref && j_type == i_tref)) {
              const bool forward = (i_type == i_tref);
              const char4 i_name = (forward) ? cdk.atom_names[i_atom] : cdk.atom_names[j_atom];
              const char4 j_name = (forward) ? cdk.atom_names[j_atom] : cdk.atom_names[i_atom];
              const int param_idx = vk.bond_param_idx[pos];
              if (nfound == 0) {
                snprintf(buffer_data, 256, " Atom Names   Indices     Eq. L0    Stiffness");
                di_body.push_back(buffer_data);
                snprintf(buffer_data, 256, "  ---- ----  ---- ----  ---------- ----------");
                di_body.push_back(buffer_data);
              }
              const int i_index = (forward) ? i_atom + 1 : j_atom + 1;
              const int j_index = (forward) ? j_atom + 1 : i_atom + 1;
              snprintf(buffer_data, 256, "  %c%c%c%c %c%c%c%c  %4d %4d  %10.4lf %10.4lf", i_name.x,
                       i_name.y, i_name.z, i_name.w, j_name.x, j_name.y, j_name.z, j_name.w,
                       i_index, j_index, vk.bond_leq[param_idx], vk.bond_keq[param_idx]);
              di_body.push_back(buffer_data);
              nfound++;
            }
          }
        }
        break;
      case StateVariable::ANGLE:
        {
          const char4 i_tref = atom_types[0];
          const char4 j_tref = atom_types[1];
          const char4 k_tref = atom_types[2];
          for (int pos = 0; pos < vk.nangl; pos++) {
            const int i_atom = vk.angl_i_atoms[pos];
            const int j_atom = vk.angl_j_atoms[pos];
            const int k_atom = vk.angl_k_atoms[pos];
            const char4 i_type = cdk.atom_types[i_atom];
            const char4 j_type = cdk.atom_types[j_atom];
            const char4 k_type = cdk.atom_types[k_atom];
            if ((i_type == i_tref && j_type == j_tref && k_type == k_tref) ||
                (i_type == k_tref && j_type == j_tref && k_type == i_tref)) {
              const bool forward = (i_type == i_tref);
              const char4 i_name = (forward) ? cdk.atom_names[i_atom] : cdk.atom_names[k_atom];
              const char4 j_name = cdk.atom_names[j_atom];
              const char4 k_name = (forward) ? cdk.atom_names[k_atom] : cdk.atom_names[i_atom];
              const int param_idx = vk.angl_param_idx[pos];
              if (nfound == 0) {
                snprintf(buffer_data, 256,
                         "    Atom Names        Indices       Eq. L0    Stiffness");
                di_body.push_back(buffer_data);
                snprintf(buffer_data, 256,
                         "  ---- ---- ----  ---- ---- ----  ---------- ----------");
                di_body.push_back(buffer_data);
              }
              const int i_index = (forward) ? i_atom + 1 : k_atom + 1;
              const int j_index = j_atom + 1;
              const int k_index = (forward) ? k_atom + 1 : i_atom + 1;
              snprintf(buffer_data, 256,
                       "  %c%c%c%c %c%c%c%c %c%c%c%c  %4d %4d %4d  %10.4lf %10.4lf",
                       i_name.x, i_name.y, i_name.z, i_name.w, j_name.x, j_name.y, j_name.z,
                       j_name.w, k_name.x, k_name.y, k_name.z, k_name.w, i_index, j_index, k_index,
                       vk.angl_theta[param_idx] * 180.0 / symbols::pi, vk.angl_keq[param_idx]);
              di_body.push_back(buffer_data);
              nfound++;
            }
          }
        }
        break;
      case StateVariable::PROPER_DIHEDRAL:
      case StateVariable::IMPROPER_DIHEDRAL:
        {
          const char4 i_tref = atom_types[0];
          const char4 j_tref = atom_types[1];
          const char4 k_tref = atom_types[2];
          const char4 l_tref = atom_types[3];
          const bool is_proper = (ask.getValenceParameter() == StateVariable::PROPER_DIHEDRAL);
          for (int pos = 0; pos < vk.ndihe; pos++) {
            switch (static_cast<TorsionKind>(vk.dihe_modifiers[pos].w)) {
            case TorsionKind::PROPER:
            case TorsionKind::PROPER_NO_14:
              if (is_proper == false) {
                continue;
              }
              break;
            case TorsionKind::IMPROPER:
            case TorsionKind::IMPROPER_NO_14:
              if (is_proper) {
                continue;
              }
            }
            const int i_atom = vk.dihe_i_atoms[pos];
            const int j_atom = vk.dihe_j_atoms[pos];
            const int k_atom = vk.dihe_k_atoms[pos];
            const int l_atom = vk.dihe_l_atoms[pos];
            const char4 i_type = cdk.atom_types[i_atom];
            const char4 j_type = cdk.atom_types[j_atom];
            const char4 k_type = cdk.atom_types[k_atom];
            const char4 l_type = cdk.atom_types[l_atom];
            if ((i_type == i_tref && j_type == j_tref && k_type == k_tref && l_type == l_tref) ||
                (i_type == l_tref && j_type == k_tref && k_type == j_tref && l_type == i_tref)) {
              const bool forward = (i_type == i_tref);
              const char4 i_name = (forward) ? cdk.atom_names[i_atom] : cdk.atom_names[l_atom];
              const char4 j_name = (forward) ? cdk.atom_names[j_atom] : cdk.atom_names[k_atom];
              const char4 k_name = (forward) ? cdk.atom_names[k_atom] : cdk.atom_names[j_atom];
              const char4 l_name = (forward) ? cdk.atom_names[l_atom] : cdk.atom_names[i_atom];
              const int param_idx = vk.dihe_param_idx[pos];
              if (nfound == 0) {
                snprintf(buffer_data, 256, "       Atom Names            Indices         "
                         "Amplitude    Phase    N");
                di_body.push_back(buffer_data);
                snprintf(buffer_data, 256, "  ---- ---- ---- ----  ---- ---- ---- ----  "
                         "---------- ---------- --");
                di_body.push_back(buffer_data);
              }
              const int i_index = (forward) ? i_atom + 1 : l_atom + 1;
              const int j_index = (forward) ? j_atom + 1 : k_atom + 1;
              const int k_index = (forward) ? k_atom + 1 : j_atom + 1;
              const int l_index = (forward) ? l_atom + 1 : i_atom + 1;
              snprintf(buffer_data, 256, "  %c%c%c%c %c%c%c%c %c%c%c%c %c%c%c%c  %4d %4d %4d %4d  "
                       "%10.4lf %10.4lf %2d", i_name.x, i_name.y, i_name.z, i_name.w, j_name.x,
                       j_name.y, j_name.z, j_name.w, k_name.x, k_name.y, k_name.z, k_name.w,
                       l_name.x, l_name.y, l_name.z, l_name.w, i_index, j_index, k_index, l_index,
                       vk.dihe_amp[param_idx], vk.dihe_phi[param_idx] * 180.0 / symbols::pi,
                       static_cast<int>(std::round(vk.dihe_freq[param_idx])));
              di_body.push_back(buffer_data);
              nfound++;
            }
          }
        }
        break;
      case StateVariable::UREY_BRADLEY:
        {
          const char4 i_tref = atom_types[0];
          const char4 k_tref = atom_types[1];
          for (int pos = 0; pos < vk.nubrd; pos++) {
            const int i_atom = vk.ubrd_i_atoms[pos];
            const int k_atom = vk.ubrd_k_atoms[pos];
            const char4 i_type = cdk.atom_types[i_atom];
            const char4 k_type = cdk.atom_types[k_atom];
            if ((i_type == i_tref && k_type == k_tref) || (i_type == k_tref && k_type == i_tref)) {
              const bool forward = (i_type == i_tref);
              const char4 i_name = (forward) ? cdk.atom_names[i_atom] : cdk.atom_names[k_atom];
              const char4 k_name = (forward) ? cdk.atom_names[k_atom] : cdk.atom_names[i_atom];
              const int param_idx = vk.ubrd_param_idx[pos];
              if (nfound == 0) {
                snprintf(buffer_data, 256, " Atom Names   Indices     Eq. L0    Stiffness");
                di_body.push_back(buffer_data);
                snprintf(buffer_data, 256, "  ---- ----  ---- ----  ---------- ----------");
                di_body.push_back(buffer_data);
              }
              const int i_index = (forward) ? i_atom + 1 : k_atom + 1;
              const int k_index = (forward) ? k_atom + 1 : i_atom + 1;
              snprintf(buffer_data, 256, "  %c%c%c%c %c%c%c%c  %4d %4d  %10.4lf %10.4lf", i_name.x,
                       i_name.y, i_name.z, i_name.w, k_name.x, k_name.y, k_name.z, k_name.w,
                       i_index, k_index, vk.ubrd_leq[param_idx], vk.ubrd_keq[param_idx]);
              di_body.push_back(buffer_data);
              nfound++;
            }
          }
        }
        break;
      case StateVariable::CHARMM_IMPROPER:
        {
          const char4 i_tref = atom_types[0];
          const char4 j_tref = atom_types[1];
          const char4 k_tref = atom_types[2];
          const char4 l_tref = atom_types[3];
          for (int pos = 0; pos < vk.ndihe; pos++) {
            const int i_atom = vk.cimp_i_atoms[pos];
            const int j_atom = vk.cimp_j_atoms[pos];
            const int k_atom = vk.cimp_k_atoms[pos];
            const int l_atom = vk.cimp_l_atoms[pos];
            const char4 i_type = cdk.atom_types[i_atom];
            const char4 j_type = cdk.atom_types[j_atom];
            const char4 k_type = cdk.atom_types[k_atom];
            const char4 l_type = cdk.atom_types[l_atom];
            if ((i_type == i_tref && j_type == j_tref && k_type == k_tref && l_type == l_tref) ||
                (i_type == l_tref && j_type == k_tref && k_type == j_tref && l_type == i_tref)) {
              const bool forward = (i_type == i_tref);
              const char4 i_name = (forward) ? cdk.atom_names[i_atom] : cdk.atom_names[l_atom];
              const char4 j_name = (forward) ? cdk.atom_names[j_atom] : cdk.atom_names[k_atom];
              const char4 k_name = (forward) ? cdk.atom_names[k_atom] : cdk.atom_names[j_atom];
              const char4 l_name = (forward) ? cdk.atom_names[l_atom] : cdk.atom_names[i_atom];
              const int param_idx = vk.cimp_param_idx[pos];
              if (nfound == 0) {
                snprintf(buffer_data, 256, "       Atom Names            Indices         "
                         "Stiffness    Phase  ");
                di_body.push_back(buffer_data);
                snprintf(buffer_data, 256, "  ---- ---- ---- ----  ---- ---- ---- ----  "
                         "---------- ----------");
                di_body.push_back(buffer_data);
              }
              const int i_index = (forward) ? i_atom + 1 : l_atom + 1;
              const int j_index = (forward) ? j_atom + 1 : k_atom + 1;
              const int k_index = (forward) ? k_atom + 1 : j_atom + 1;
              const int l_index = (forward) ? l_atom + 1 : i_atom + 1;
              snprintf(buffer_data, 256, "  %c%c%c%c %c%c%c%c %c%c%c%c %c%c%c%c  %4d %4d %4d %4d  "
                       "%10.4lf %10.4lf", i_name.x, i_name.y, i_name.z, i_name.w, j_name.x,
                       j_name.y, j_name.z, j_name.w, k_name.x, k_name.y, k_name.z, k_name.w,
                       l_name.x, l_name.y, l_name.z, l_name.w, i_index, j_index, k_index, l_index,
                       vk.cimp_keq[param_idx], vk.cimp_phi[param_idx] * 180.0 / symbols::pi);
              di_body.push_back(buffer_data);
              nfound++;
            }
          }
        }
        break;
      case StateVariable::CMAP:
        {
          const char4 i_tref = atom_types[0];
          const char4 j_tref = atom_types[1];
          const char4 k_tref = atom_types[2];
          const char4 l_tref = atom_types[3];
          const char4 m_tref = atom_types[4];
          for (int pos = 0; pos < vk.ndihe; pos++) {
            const int i_atom = vk.cmap_i_atoms[pos];
            const int j_atom = vk.cmap_j_atoms[pos];
            const int k_atom = vk.cmap_k_atoms[pos];
            const int l_atom = vk.cmap_l_atoms[pos];
            const int m_atom = vk.cmap_m_atoms[pos];
            const char4 i_type = cdk.atom_types[i_atom];
            const char4 j_type = cdk.atom_types[j_atom];
            const char4 k_type = cdk.atom_types[k_atom];
            const char4 l_type = cdk.atom_types[l_atom];
            const char4 m_type = cdk.atom_types[m_atom];
            if ((i_type == i_tref && j_type == j_tref && k_type == k_tref && l_type == l_tref &&
                 m_type == m_tref) ||
                (i_type == m_tref && j_type == l_tref && k_type == k_tref && l_type == j_tref &&
                 m_type == i_tref)) {
              const bool forward = (i_type == i_tref);
              const char4 i_name = (forward) ? cdk.atom_names[i_atom] : cdk.atom_names[m_atom];
              const char4 j_name = (forward) ? cdk.atom_names[j_atom] : cdk.atom_names[l_atom];
              const char4 k_name = (forward) ? cdk.atom_names[k_atom] : cdk.atom_names[k_atom];
              const char4 l_name = (forward) ? cdk.atom_names[l_atom] : cdk.atom_names[j_atom];
              const char4 m_name = (forward) ? cdk.atom_names[m_atom] : cdk.atom_names[i_atom];
              const int param_idx = vk.cimp_param_idx[pos];
              if (nfound == 0) {
                snprintf(buffer_data, 256, "          Atom Names                 Indices          "
                         " Map");
                di_body.push_back(buffer_data);
                snprintf(buffer_data, 256, "  ---- ---- ---- ---- ----  ---- ---- ---- ---- ----  "
                        "----");
                di_body.push_back(buffer_data);
              }
              const int i_index = (forward) ? i_atom + 1 : m_atom + 1;
              const int j_index = (forward) ? j_atom + 1 : l_atom + 1;
              const int k_index = (forward) ? k_atom + 1 : k_atom + 1;
              const int l_index = (forward) ? l_atom + 1 : j_atom + 1;
              const int m_index = (forward) ? m_atom + 1 : i_atom + 1;
              snprintf(buffer_data, 256, "  %c%c%c%c %c%c%c%c %c%c%c%c %c%c%c%c %c%c%c%c  %4d %4d "
                       "%4d %4d %4d  %4d", i_name.x, i_name.y, i_name.z, i_name.w, j_name.x,
                       j_name.y, j_name.z, j_name.w, k_name.x, k_name.y, k_name.z, k_name.w,
                       l_name.x, l_name.y, l_name.z, l_name.w, m_name.x, m_name.y, m_name.z,
                       m_name.w, i_index, j_index, k_index, l_index, m_index,
                       vk.cmap_surf_idx[param_idx]);
              di_body.push_back(buffer_data);
              nfound++;
            }
          }
        }
        break;
      case StateVariable::RESTRAINT:
      case StateVariable::VDW:
      case StateVariable::VDW_ONE_FOUR:
      case StateVariable::ELECTROSTATIC:
      case StateVariable::ELEC_ONE_FOUR:
      case StateVariable::GENERALIZED_BORN:
      case StateVariable::KINETIC:
      case StateVariable::PRESSURE:
      case StateVariable::VIRIAL_11:
      case StateVariable::VIRIAL_12:
      case StateVariable::VIRIAL_22:
      case StateVariable::VIRIAL_13:
      case StateVariable::VIRIAL_23:
      case StateVariable::VIRIAL_33:
      case StateVariable::VOLUME:
      case StateVariable::TEMPERATURE_ALL:
      case StateVariable::TEMPERATURE_PROTEIN:
      case StateVariable::TEMPERATURE_LIGAND:
      case StateVariable::TEMPERATURE_SOLVENT:
      case StateVariable::DU_DLAMBDA:
      case StateVariable::POTENTIAL_ENERGY:
      case StateVariable::TOTAL_ENERGY:
      case StateVariable::ALL_STATES:
        break;
      }

      // Mark the number of findings on the first line
      if (nfound == 0) {
        di_body[0] = std::string("No");
      }
      else {
        di_body[0] = std::string("A total of ") + std::to_string(nfound);
      }
      di_body[0].append(" instances of a " + getEnumerationName(ask.getValenceParameter()) +
                        " involving atom type");
      const int natypes = atom_types.size();
      if (natypes > 1) {
        di_body[0] += 's';
      }
      for (int i = 0; i < natypes; i++) {
        di_body[0] += ' ';
        di_body[0] += '[';
        if (atom_types[i].x != ' ') di_body[0] += atom_types[i].x;
        if (atom_types[i].y != ' ') di_body[0] += atom_types[i].y;
        if (atom_types[i].z != ' ') di_body[0] += atom_types[i].z;
        if (atom_types[i].w != ' ') di_body[0] += atom_types[i].w;
        di_body[0] += ']';        
      }
      di_body[0].append(" were found.");
    }
    break;
  case DataRequestKind::STRING:

    // The text is the user-supplied message
    di_body.push_back(ask.getMessage());
    break;
  case DataRequestKind::ALL_KINDS:
    break;
  }
  data_items.emplace_back(ask, di_body);
  data_item_count = data_items.size();
  compareExternalRegistryNumbers(data_items.back().getExternalRegistryNumber());
}

//-------------------------------------------------------------------------------------------------
MdlMolDataItemKind MdlMol::getDataItemKind(const int item_index) const {
  checkDataItemIndex(item_index, "getDataItemKind");
  return data_items[item_index].getKind();
}

//-------------------------------------------------------------------------------------------------
StateVariable MdlMol::getTrackedState(const int item_index) const {
  checkDataItemIndex(item_index, "getTrackedState");
  return data_items[item_index].getTrackedState();
}

//-------------------------------------------------------------------------------------------------
void MdlMol::addLineToDataItem(const std::string &text, int item_index) {
  checkDataItemIndex(item_index, "addLineToDataItem");
  data_items[item_index].addDataLine(text);
}

//-------------------------------------------------------------------------------------------------
const std::string& MdlMol::getDataItemName(const int item_index) const {
  checkDataItemIndex(item_index, "getDataItemName");
  return data_items[item_index].getItemName();
}

//-------------------------------------------------------------------------------------------------
const std::string& MdlMol::getDataItemOutputName(const int item_index) const {
  checkDataItemIndex(item_index, "getDataItemOutputName");
  return data_items[item_index].getOutputItemName();
}

//-------------------------------------------------------------------------------------------------
int MdlMol::getDataItemIndex(const std::string &item_name) const {
  for (int i = 0; i < data_item_count; i++) {
    if (item_name == data_items[i].getItemName()) {
      return i;
    }
  }
  return data_item_count;
}

//-------------------------------------------------------------------------------------------------
const std::vector<std::string>& MdlMol::getDataItemContent(const int item_index) const {
  checkDataItemIndex(item_index, "getDataItemContent");
  return data_items[item_index].getBody();
}

//-------------------------------------------------------------------------------------------------
std::vector<std::string> MdlMol::getDataItemContent(const std::string &item_name) const {
  std::vector<std::string> result;
  for (int i = 0; i < data_item_count; i++) {
    if (item_name == data_items[i].getItemName()) {
      result = data_items[i].getBody();
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void MdlMol::writeMdl(std::ofstream *foutp, const MdlMolVersion vformat) const {
  const TextFile result(writeMdl(vformat), TextOrigin::RAM);
  writeFrame(foutp, result);
}

//-------------------------------------------------------------------------------------------------
void MdlMol::writeMdl(const std::string &fname, const MdlMolVersion vformat,
                      const PrintSituation expectation) const {
  const std::string activity = (data_item_count > 0) ?
    "Open an output file for writing an MDL MOL format structure" :
    "Open an SDF archive for writing MDL MOL format output with additional data items";
  std::ofstream foutp = openOutputFile(fname, expectation, activity);
  writeMdl(&foutp, vformat);
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
std::string MdlMol::writeMdl(const MdlMolVersion vformat) const {

  // Build the result based on the MDL MOL leading three lines, which are common to both V2000 and
  // V3000 formats.
  std::string result(title + '\n' + software_details + '\n' + general_comment + '\n');
  std::string buffer(512, ' ');
  char* buffer_data = buffer.data();
  switch (vformat) {
  case MdlMolVersion::V2000:

    // Write out the counts line.
    snprintf(buffer.data(), 512, "%3d%3d%3d   %3d%3d            999 V2000\n", atom_count,
             bond_count, list_count, static_cast<int>(chirality), stext_entry_count);
    result.append(buffer_data, 40);
    
    // Write out the atom block.
    for (int i = 0; i < atom_count; i++) {
      snprintf(buffer_data, 512, "%10.4lf%10.4lf%10.4lf %c%c%c%2d%3d%3d%3d%3d%3d%3d  0  "
               "0%3d%3d%3d\n", coordinates[i].x, coordinates[i].y, coordinates[i].z,
               atomic_symbols[i].x, atomic_symbols[i].y, atomic_symbols[i].z,
               getIsotopicShiftCode(i), getFormalChargeCode(i), static_cast<int>(parities[i]),
               getImplicitHydrogenCode(i), static_cast<int>(stereo_considerations[i]),
               static_cast<int>(valence_connections[i]),
               static_cast<int>(hydrogenation_protocol[i]), atom_atom_mapping_count[i],
               static_cast<int>(orientation_stability[i]),
               static_cast<int>(exact_change_enforced[i]));
      result.append(buffer_data, 70);
    }

    // Write out the bond block.
    for (int i = 0; i < bond_count; i++) {

      // Add 1 to the bond atom indices to get back into the file format.  This adjustment is
      // automated for the properties.
      snprintf(buffer.data(), 256, "%3d%3d%3d%3d  0%3d%3d\n", bonds[i].getFirstAtom() + 1,
               bonds[i].getSecondAtom() + 1, static_cast<int>(bonds[i].getOrder()),
               static_cast<int>(bonds[i].getStereochemistry()),
               static_cast<int>(bonds[i].getRingStatus()),
               static_cast<int>(bonds[i].getReactivity()));
      result.append(buffer_data, 22);
    }

    // Write out the atom list block, if appropriate.
    if (property_element_lists == false) {
      for (int i = 0; i < list_count; i++) {
        const int n_entry = element_lists[i].getEntryCount();
        snprintf(buffer_data, 256, "%3d %c    %d", element_lists[i].getAttachmentPoint(),
                 element_lists[i].getExclusionCode(), n_entry);
        int nchar = 10;
        for (int j = 0; j < n_entry; j++) {
          snprintf(&buffer_data[nchar], 512 - nchar, " %3d", element_lists[i].getEntry(j));
          nchar += 4;
        }
        buffer_data[nchar] = '\n';
        nchar++;
        buffer_data[nchar] = '\0';
        result.append(buffer_data, nchar);
      }
    }

    // Write out the S-text block

    // Write out the properties block.
    for (int i = 0; i < properties_count; i++) {
      result.append(properties[i].getMdlText());
    }

    // Write out the terminal line (not considered a property, even in the V2000 format)
    result.append("M  END\n");
    
    break;
  case MdlMolVersion::V3000:
  case MdlMolVersion::UNKNOWN:

    // Make the default writing method the V3000 format.
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void MdlMol::writeDataItems(std::ofstream *foutp, const int mol_index) const {
  const TextFile result(writeDataItems(mol_index), TextOrigin::RAM);
  writeFrame(foutp, result);
}

//-------------------------------------------------------------------------------------------------
void MdlMol::writeDataItems(const std::string &fname, const PrintSituation expectation,
                            const int mol_index) const {
  const std::string activity("Open an SDF archive for writing MDL MOL format output with "
                             "additional data items");
  std::ofstream foutp = openOutputFile(fname, expectation, activity);
  writeDataItems(&foutp, mol_index);
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
std::string MdlMol::writeDataItems(const int mol_index) const {
  std::string result;
  for (int i = 0; i < data_item_count; i++) {

    // Compose the header line
    std::string header_line("> ");
    bool required_id = false;
    if (data_items[i].placeInternalRegnoInHeader()) {
      header_line += std::to_string(mol_index);
    }
    if (data_items[i].placeExternalRegnoInHeader()) {
      header_line += '(' + data_items[i].getExternalRegistryNumber() + ')';
    }
    if (data_items[i].placeTitleInHeader()) {
      header_line += '<' + data_items[i].getOutputItemName() + '>';
      required_id = true;
    }
    if (data_items[i].placeMaccsIIFieldInHeader()) {
      header_line += "DT" + std::to_string(data_items[i].getMaccsFieldNumber());
      required_id = true;
    }
    if (data_items[i].noteArchivesInHeader()) {
      header_line += " FROM ARCHIVES";
    }
    header_line += '\n';
    result.append(header_line);
    const int ndl = data_items[i].getDataLineCount();
    for (int j = 0; j < ndl; j++) {
      result.append(data_items[i].getDataLine(j));
      result += '\n';
    }

    // There must be a blank line between any data item and the next, but not between the last
    // data item and the SD file's terminating "$$$$" card.
    if (i < data_item_count - 1) {
      result += '\n';
    }
  }
  result.append("$$$$\n");
  return result;
}

//-------------------------------------------------------------------------------------------------
void MdlMol::validateAtomIndex(const int index, const char* caller) const {
  if (index < 0 || index >= atom_count) {
    rtErr("Atom index " + std::to_string(index) + " is out of bounds for an MDL MOL entry with " +
          std::to_string(atom_count) + " atoms.", "MdlMol", caller);
  }
}

//-------------------------------------------------------------------------------------------------
void MdlMol::allocate() {

  // Atom property fields are resized and then set as part of a loop in the parent MdlMol
  // constructor.
  coordinates.resize(atom_count);
  atomic_symbols.resize(atom_count, default_mdl_atomic_symbol);
  atomic_numbers.resize(atom_count, default_mdl_atomic_number);
  formal_charges.resize(atom_count, default_mdl_formal_charge);
  radicals.resize(atom_count, default_mdl_radical_state);
  isotopic_shifts.resize(atom_count, default_mdl_isotopic_shift);
  parities.resize(atom_count, default_mdl_stereo_parity);
  implicit_hydrogens.resize(atom_count, default_mdl_implicit_hydrogen);
  stereo_considerations.resize(atom_count, default_mdl_stereo_considerations);
  valence_connections.resize(atom_count, default_mdl_valence_connections);
  atom_atom_mapping_count.resize(atom_count, default_mdl_map_count);
  exact_change_enforced.resize(atom_count, default_mdl_exact_change);
  hydrogenation_protocol.resize(atom_count, default_hydrogenation);
  orientation_stability.resize(atom_count, default_mdl_stereo_retention);

  // Other arrays are reserved and built with emplace_back().  The MDL MOL properties array is
  // not reserved to any specific length, however, as the number of properties is not known from
  // the counts line, where it is set to 999 by default.
  bonds.reserve(bond_count);
  element_lists.reserve(list_count);
  stext_entries.reserve(stext_entry_count);
}

//-------------------------------------------------------------------------------------------------
void MdlMol::allocate(const ChemicalDetailsKit &cdk, const NonbondedKit<double> &nbk,
                      const int molecule_index) {

  // Determine the number of atoms and bonds in the molecule
  if (molecule_index < 0 || molecule_index >= cdk.nmol) {
    rtErr("Molecule index " + std::to_string(molecule_index) + " is invalid for a topology with " +
          std::to_string(cdk.nmol) + " molecules.", "MdlMol", "allocate");
  }
  const int hlim = cdk.mol_limits[molecule_index + 1];
  const int llim = cdk.mol_limits[molecule_index];
  int tmp_atoms = 0;
  for (int i = llim; i < hlim; i++) {
    tmp_atoms += (cdk.z_numbers[i] > 0);
  }
  atom_count = tmp_atoms;
  coordinates.resize(atom_count);
  atomic_symbols.resize(atom_count);
  atomic_numbers.resize(atom_count);
  formal_charges.resize(atom_count);
  radicals.resize(atom_count);
  isotopic_shifts.resize(atom_count);
  parities.resize(atom_count);
  implicit_hydrogens.resize(atom_count);
  stereo_considerations.resize(atom_count);
  valence_connections.resize(atom_count);
  atom_atom_mapping_count.resize(atom_count);
  exact_change_enforced.resize(atom_count);
  hydrogenation_protocol.resize(atom_count);
  orientation_stability.resize(atom_count);
  
  // Determine the number of bonds in the molecule
  int tmp_bonds = 0;
  for (int i = llim; i < hlim; i++) {
    if (cdk.z_numbers[i] <= 0) {
      continue;
    }
    const int atom_idx = cdk.mol_contents[i];
    int nconn = 0;
    for (int j = nbk.nb12_bounds[atom_idx]; j < nbk.nb12_bounds[atom_idx + 1]; j++) {
      nconn += (cdk.z_numbers[nbk.nb12x[j]] > 0);
    }
    valence_connections[i - llim] = nconn;
    tmp_bonds += nconn;
  }

  // Divide the number of bonds by two here, to avoid unecessary arithmetic and also prevent
  // integer roundoff error from affecting the tally.
  bond_count = tmp_bonds / 2;
  bonds.reserve(bond_count);
  element_lists.resize(0);
  stext_entries.resize(0);
  element_lists.shrink_to_fit();
  stext_entries.shrink_to_fit();
}

//-------------------------------------------------------------------------------------------------
int MdlMol::getIsotopicShiftCode(const int atom_index) const {
  if (property_isotopes) {
    return 0;
  }
  else {
    return isotopic_shifts[atom_index];
  }
  __builtin_unreachable();
}
 
//-------------------------------------------------------------------------------------------------
void MdlMol::interpretFormalCharge(const int charge_in, const int atom_index) {
  switch (charge_in) {
  case 0:
    formal_charges[atom_index] = 0;
    break;
  case 1:
    formal_charges[atom_index] = 3;
    break;
  case 2:
    formal_charges[atom_index] = 2;
    break;
  case 3:
    formal_charges[atom_index] = 1;
    break;
  case 4:
    radicals[atom_index] = RadicalState::DOUBLET;
    break;
  case 5:
    formal_charges[atom_index] = -1;
    break;
  case 6:
    formal_charges[atom_index] = -2;
    break;
  case 7:
    formal_charges[atom_index] = -3;
    break;
  default:
    rtErr("A formal charge code of " + std::to_string(charge_in) + " is invalid for an MDL MOL "
          "entry.  Title of entry: \"" + title + "\".", "MdlMol", "interpretFormalCharge");    
  }
}

//-------------------------------------------------------------------------------------------------
int MdlMol::getFormalChargeCode(const int atom_index) const {
  if (property_formal_charges) {
    return 0;
  }
  switch (formal_charges[atom_index]) {
  case 0:
    return (radicals[atom_index] == RadicalState::DOUBLET) ? 4 : 0;
  case 1:
    return 3;
  case 2:
    return 2;
  case 3:
    return 1;
  case -1:
    return 5;
  case -2:
    return 6;
  case -3:
    return 7;
  default:

    // Formal charges outside the V2000 range must be handled by properties, and any such case
    // implies that all formal charges are handled this way.
    rtErr("A formal charge that cannot be expressed in the atoms block of a V2000 format MDL MOL "
          "file exists, but there are no properties to account for it.", "MdlMol",
          "getFormalChargeCode");
  }
  __builtin_unreachable();
}
 
//-------------------------------------------------------------------------------------------------
MolObjAtomStereo MdlMol::interpretStereoParity(const int setting_in) {
  switch (setting_in) {
  case 0:
    return MolObjAtomStereo::NOT_STEREO;
  case 1:
    return MolObjAtomStereo::ODD;
  case 2:
    return MolObjAtomStereo::EVEN;
  case 3:
    return MolObjAtomStereo::UNMARKED;
  default:
    rtErr("A stereochemical parity setting of " + std::to_string(setting_in) + " is invalid.  "
          "Title of entry: \"" + title + "\".", "MdlMol", "interpretStereoParity");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void MdlMol::interpretImplicitHydrogenContent(const int nh_in, const int atom_index) {
  if (nh_in > 5 || nh_in < 0) {
    rtErr("An implicit hydrogen content of " + std::to_string(nh_in) + " would imply " +
          std::to_string(nh_in - 1) + " hydrogens can be inferred around an atom in MDL MOL "
          "entry \"" + title + "\".", "MdlMol", "interpretImplicitHydrogenContent");
  }
  if (nh_in > 0) {
    implicit_hydrogens[atom_index] = nh_in - 1;
    hydrogenation_protocol[atom_index] = (nh_in == 1) ? HydrogenAssignment::DO_NOT_HYDROGENATE :
                                                        HydrogenAssignment::VALENCE_SHELL;
  }
  else {

    // An implicit hydrogen indicator of 0 does not correspond to H0, H1, H2, H3, or H4, but it is
    // very common.  This final possibility implies "free hydrogen content." While the actual
    // number of hydrogens will read 0, the flag will be set to apply as many as are needed to
    // fill the valence shell given the bonding considerations.
    implicit_hydrogens[atom_index] = 0;
    hydrogenation_protocol[atom_index] = HydrogenAssignment::VALENCE_SHELL;
  }
}

//-------------------------------------------------------------------------------------------------
int MdlMol::getImplicitHydrogenCode(const int atom_index) const {
  switch (hydrogenation_protocol[atom_index]) {
  case HydrogenAssignment::VALENCE_SHELL:
    return (implicit_hydrogens[atom_index] == 0) ? 0 : implicit_hydrogens[atom_index] + 1;
  case HydrogenAssignment::DO_NOT_HYDROGENATE:
    return 1;
  }
  __builtin_unreachable();
}
 
//-------------------------------------------------------------------------------------------------
bool MdlMol::interpretBooleanValue(const int value_in, const std::string &desc) {
  if (value_in != 0 && value_in != 1) {
    rtErr("A directive of " + std::to_string(value_in) + " is invalid when " + desc + ".  Title "
          "of entry: \"" + title + "\".", "MdlMol", "interpretBooleanValue");
  }
  return (value_in == 1);
}

//-------------------------------------------------------------------------------------------------
int MdlMol::interpretValenceNumber(const int count_in) {
  if (count_in == 15) {
    return 0;
  }
  else if (count_in >= 0 && count_in < 15) {
    return count_in;
  }
  else {
    rtErr("An atom cannot have " + std::to_string(count_in) + " valence connections, as is the "
          "case for one atom in MDL MOL entry \"" + title + "\".", "MdlMol",
          "interpretValenceNumber");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
StereoRetention MdlMol::interpretStereoStability(const int code_in) {
  switch (code_in) {
  case 0:
    return StereoRetention::NOT_APPLIED;
  case 1:
    return StereoRetention::INVERTED;
  case 2:
    return StereoRetention::RETAINED;
  default:
    rtErr("A stereochemistry retention code of " + std::to_string(code_in) + " is invalid in MDL "
          "MOL entry \"" + title + "\".", "MdlMol", "interpretStereoStability");
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
void MdlMol::updateV2kAtomAttributes() {

  // Return immediately if the version is not V2000
  switch (version_no) {
  case MdlMolVersion::V2000:
    break;
  case MdlMolVersion::V3000:
  case MdlMolVersion::UNKNOWN:
    return;
  }
  
  // Scan all properties for items that would invalidate atom-block information.  Wipe the
  // relevant arrays.
  for (int i = 0; i < properties_count; i++) {
    if (properties[i].getCode() == char4({ 'C', 'H', 'G', 'M' })) {
      for (int j = 0; j < atom_count; j++) {
        formal_charges[j] = 0;
      }
      property_formal_charges = true;
    }
    if (properties[i].getCode() == char4({ 'R', 'A', 'D', 'M' })) {
      for (int j = 0; j < atom_count; j++) {
        formal_charges[j] = 0;
        radicals[j] = RadicalState::NONE;
      }
      property_radicals = true;
    }
    if (properties[i].getCode() == char4({ 'I', 'S', 'O', 'M' })) {
      for (int j = 0; j < atom_count; j++) {
        isotopic_shifts[j] = 0;
      }
      property_isotopes = true;
    }
    if (properties[i].getCode() == char4({ 'A', 'L', 'S', 'M' })) {
      property_element_lists = true;
    }
  }

  // Scan the properties again and add details to the atoms.
  for (int i = 0; i < properties_count; i++) {
    const int n_entry = properties[i].getEntryCount();
    if (properties[i].getCode() == char4({ 'C', 'H', 'G', 'M' })) {
      for (int j = 0; j < n_entry; j++) {
        formal_charges[properties[i].getIntegerValue(j, 0)] = properties[i].getIntegerValue(j, 1);
      }
    }
    if (properties[i].getCode() == char4({ 'R', 'A', 'D', 'M' })) {
      for (int j = 0; j < n_entry; j++) {
        const int atom_idx = properties[i].getIntegerValue(j, 0);
        switch (properties[i].getIntegerValue(j, 1)) {
        case 0:
          radicals[atom_idx] = RadicalState::NONE;
          break;
        case 1:
          radicals[atom_idx] = RadicalState::SINGLET;
          break;
        case 2:
          radicals[atom_idx] = RadicalState::DOUBLET;
          break;
        case 3:
          radicals[atom_idx] = RadicalState::TRIPLET;
          break;
        }
      }
    }
    if (properties[i].getCode() == char4({ 'I', 'S', 'O', 'M' })) {
      for (int j = 0; j < n_entry; j++) {
        isotopic_shifts[properties[i].getIntegerValue(j, 0)] = properties[i].getIntegerValue(j, 1);
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void MdlMol::hydrogenate() {

}

//-------------------------------------------------------------------------------------------------
void MdlMol::compareExternalRegistryNumbers(const std::string &regno_in) {
  if (regno_in.size() > 0LLU) {
    if (external_regno.size() == 0LLU) {
      external_regno = regno_in;
    }
    else if (regno_in == external_regno) {
      return;
    }
    else {
      switch (policy) {
      case ExceptionResponse::DIE:
        rtErr("The external registry number associated with a data item to add (" + regno_in +
              ") is inconsistent with the external registry number already associated with the "
              "MDL MOL entry (" + external_regno + ").", "MdlMol",
              "compareExternalRegistryNumbers");
      case ExceptionResponse::WARN:
        rtWarn("The external registry number associated with a data item to add (" + regno_in +
               ") is inconsistent with the external registry number already associated with the "
               "MDL MOL entry (" + external_regno + ").  The existing registry number will take "
               "precedence.", "MdlMol", "compareExternalRegistryNumbers");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
    }
  }
}

//-------------------------------------------------------------------------------------------------
void MdlMol::checkAtomCount(const int ext_atom_count) const {
  if (ext_atom_count != atom_count) {
    rtErr("The number of atoms coming in (" + std::to_string(ext_atom_count) + ") is not "
          "consistent with the number of atoms in the molecule (" + std::to_string(atom_count) +
          ").", "MdlMol", "checkAtomCount");
  }
}

//-------------------------------------------------------------------------------------------------
void MdlMol::checkDataItemIndex(const int item_index, const char* caller) const {
  if (item_index < 0 || item_index >= data_item_count) {
    rtErr("An index of " + std::to_string(item_index) + " is invalid for an MDL MOL entry with " +
          std::to_string(data_item_count) + " entries.", "MdlMol", caller);
  }
}

//-------------------------------------------------------------------------------------------------
void MdlMol::addAtomIntProperty(const std::vector<int2> &notable_data, const char4 prcode) {
  const size_t nd_count = notable_data.size();
  if (nd_count > 0) {
    std::vector<int> int_data_in(2 * nd_count);
    for (size_t i = 0; i < nd_count; i++) {
      int_data_in[(2 * i)    ] = notable_data[i].x;
      int_data_in[(2 * i) + 1] = notable_data[i].y;
    }
    const MdlMolProperty qprop(prcode, -1, nd_count, 2, false, -1, 6, 9, { 0, 4, 8 },
                               { MolObjPropField::INTEGER, MolObjPropField::INTEGER },
                               { MolObjIndexKind::ATOM, MolObjIndexKind::OTHER }, int_data_in);
    addProperty(qprop);
  }
}

//-------------------------------------------------------------------------------------------------
void MdlMol::transferTopologicalDetails(const ChemicalFeatures *chemfe, const int molecule_index) {
  const AtomGraph *ag = chemfe->getTopologyPointer();
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  const ValenceKit<double> vk = ag->getDoublePrecisionValenceKit();
  const ChemicalFeaturesReader chemfer = chemfe->data();
  const int llim = cdk.mol_limits[molecule_index];
  const int hlim = cdk.mol_limits[molecule_index + 1];

  // To ensure that the correspondence of atoms is maintained even if the order of the molecule's
  // atoms in the topology is not consecutive, make a map of the atoms as they will appear in this
  // MDL MOL object and their indices as they are found in the topology, as well as a map of the
  // indices in the topology (x member) and in the MDL mol object, taking the molecule's minimum
  // index as an offset to avoid having to allocate huge amounts of memory in order to transfer,
  // say, a single ligand out of a very large simulation.
  int min_top_idx = cdk.natom;
  int max_top_idx = 0;
  for (int i = llim; i < hlim; i++) {
    min_top_idx = std::min(cdk.mol_contents[i], min_top_idx);
    max_top_idx = std::max(cdk.mol_contents[i], max_top_idx);
  }
  std::vector<int> top2mol_map(max_top_idx - min_top_idx + 1, -1);
  std::vector<int> mol2top_map(atom_count);
  int mol_atom_idx = 0;
  for (int i = llim; i < hlim; i++) {
    if (cdk.z_numbers[cdk.mol_contents[i]] <= 0) {
      continue;
    }
    mol2top_map[mol_atom_idx] = cdk.mol_contents[i];
    top2mol_map[cdk.mol_contents[i] - min_top_idx] = mol_atom_idx;
    mol_atom_idx++;
  }

  // Add formal charges
  const std::vector<int>& fc_viewer = chemfe->getZeroKelvinFormalCharges();
  const std::vector<int>& fe_viewer = chemfe->getZeroKelvinFreeElectrons();
  std::vector<int2> nonzero_formal_charges;
  std::vector<int2> noteworthy_radicals;
  std::vector<int2> nonzero_isotopic_shifts;
  for (int i = 0; i < atom_count; i++) {
    const int top_idx = mol2top_map[i];
    atomic_numbers[i] = cdk.z_numbers[top_idx];
    const char2 isymb = zNumberToSymbol(cdk.z_numbers[top_idx]);
    atomic_symbols[i].x = isymb.x;
    atomic_symbols[i].y = isymb.y;
    atomic_symbols[i].z = ' ';
    atomic_symbols[i].w = ' ';
    formal_charges[i] = fc_viewer[top_idx];
    if (formal_charges[i] != 0) {
      nonzero_formal_charges.push_back({ i, formal_charges[i] });
    }
    radicals[i] = (fabs(chemfer.free_electrons[top_idx] - 1.0) < 0.1) ? RadicalState::SINGLET :
                                                                        RadicalState::NONE;
    if (radicals[i] != RadicalState::NONE) {
      noteworthy_radicals.push_back({ i, static_cast<int>(radicals[i]) });
    }
    isotopic_shifts[i] = round(cdk.masses[top_idx] - zNumberToNaturalMass(cdk.z_numbers[top_idx]));
    if (isotopic_shifts[i]) {
      nonzero_isotopic_shifts.push_back({ i, isotopic_shifts[i] });
    }
    switch (chemfe->getAtomChirality(top_idx)) {
    case ChiralOrientation::NONE:
      if (chemfe->chiralitiesComputed()) {
        parities[i] = MolObjAtomStereo::NOT_STEREO;
        stereo_considerations[i] = false;
      }
      else {

        // If the chiralities of atoms are not known, say that each atom has stereo considerations
        // if it has four bonds to it.
        if (valence_connections[i] == 4) {
          parities[i] = MolObjAtomStereo::UNMARKED;
          stereo_considerations[i] = true;
        }
        else {
          parities[i] = MolObjAtomStereo::NOT_STEREO;
          stereo_considerations[i] = false;
        }
      }
      break;
    case ChiralOrientation::RECTUS:
      parities[i] = MolObjAtomStereo::EVEN;
      stereo_considerations[i] = true;
      break;
    case ChiralOrientation::SINISTER:
      parities[i] = MolObjAtomStereo::ODD;
      stereo_considerations[i] = true;
      break;
    }
    implicit_hydrogens[i] = 0;
    atom_atom_mapping_count[i] = 0;
    exact_change_enforced[i] = false;
    hydrogenation_protocol[i] = HydrogenAssignment::DO_NOT_HYDROGENATE;
    orientation_stability[i] = StereoRetention::NOT_APPLIED;
  }

  // Create properties for significant formal charges, radicals, and isotopic shifts.
  addAtomIntProperty(nonzero_formal_charges, { 'C', 'H', 'G', 'M' });
  addAtomIntProperty(noteworthy_radicals, { 'R', 'A', 'D', 'M' });
  addAtomIntProperty(nonzero_isotopic_shifts, { 'I', 'S', 'O', 'M' });
  
  // Loop over all atoms of the molecule, adding bonds which are relevant to the MDL MOL entry.
  const std::vector<int>& bo_viewer = chemfe->getZeroKelvinBondOrders();
  for (int i = llim; i < hlim; i++) {
    const int itop_idx = cdk.mol_contents[i];
    if (cdk.z_numbers[itop_idx] <= 0) {
      continue;
    }
    const int imol_idx = top2mol_map[itop_idx - min_top_idx];
    for (int j = vk.bond_asgn_bounds[itop_idx]; j < vk.bond_asgn_bounds[itop_idx + 1]; j++) {
      const int jtop_idx = vk.bond_asgn_atoms[j];
      if (cdk.z_numbers[jtop_idx] <= 0) {
        continue;
      }
      const int jmol_idx = top2mol_map[jtop_idx - min_top_idx];
      MdlMolBondOrder ij_order;
      switch(bo_viewer[vk.bond_asgn_terms[j]]) {
      case 1:
        ij_order = MdlMolBondOrder::SINGLE;
        break;

      case 2:
        ij_order = MdlMolBondOrder::DOUBLE;
        break;
      case 3:
        ij_order = MdlMolBondOrder::TRIPLE;
        break;
      default:
        rtErr("A bond of order " + std::to_string(bo_viewer[vk.bond_asgn_terms[j]]) +
              " (between atoms " + char4ToString(cdk.atom_names[itop_idx]) + " and " +
              char4ToString(cdk.atom_names[jtop_idx]) + ") has no cognate in the MDL MOL "
              "designations.", "MdlMol", "transferTopologicalDetails");
      }
      MolObjRingState ij_ring_state =  (chemfe->bondIsInRing(vk.bond_asgn_terms[j])) ?
                                       MolObjRingState::RING : MolObjRingState::CHAIN;
      bonds.emplace_back(imol_idx, jmol_idx, ij_order, MdlMolBondStereo::NOT_STEREO,
                         ij_ring_state, MolObjReactionCenter::NON_CENTER);
    }
  }
}

//-------------------------------------------------------------------------------------------------
std::vector<MdlMolBond> operator+(const std::vector<MdlMolBond> &lhs,
				  const std::vector<MdlMolBond> &rhs) {
  std::vector<MdlMolBond> result(lhs.begin(), lhs.end());
  result.insert(result.end(), rhs.begin(), rhs.end());
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<MdlMol> readStructureDataFile(const TextFile &tf, const int low_frame_limit,
                                          const int high_frame_limit,
                                          const CaseSensitivity capitalization,
                                          const ExceptionResponse policy) {
  std::vector<MdlMol> result;
  
  // Find the limits for different MDL MOL entries
  const std::vector<int2> mol_entry_limits = findSdfMolEntryLimits(tf);
  const int nsection = mol_entry_limits.size();
  int actual_low_limit, actual_high_limit;
  if (low_frame_limit >= nsection) {
    rtErr("An SD file with " + std::to_string(nsection) + " frames cannot be read starting at "
          "frame index " + std::to_string(low_frame_limit) + ".", "readStructureDataFile");
  }
  else if (low_frame_limit < 0 || high_frame_limit < 0 || high_frame_limit >= nsection) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("The frame range " + std::to_string(low_frame_limit) + " to " +
            std::to_string(high_frame_limit) + " is invalid for a file with " +
            std::to_string(nsection) + " frames.", "readStructureDataFile");
    case ExceptionResponse::WARN:
      rtWarn("The frame range " + std::to_string(low_frame_limit) + " to " +
             std::to_string(high_frame_limit) + " is invalid for a file with " +
             std::to_string(nsection) + " frames.  Only the valid range will be taken",
             "readStructureDataFile");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    actual_low_limit = std::max(0, low_frame_limit);
    actual_high_limit = (high_frame_limit < low_frame_limit || high_frame_limit >= nsection) ?
                        nsection - 1 : high_frame_limit;
  }
  else {
    actual_low_limit = low_frame_limit;
    actual_high_limit = (high_frame_limit < low_frame_limit) ? nsection - 1: high_frame_limit;
  }
  result.reserve(actual_high_limit - actual_low_limit + 1);
  for (int i = actual_low_limit; i <= actual_high_limit; i++) {
    result.emplace_back(tf, mol_entry_limits[i].x, mol_entry_limits[i].y, capitalization, policy);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<MdlMol> readStructureDataFile(const TextFile &tf, const CaseSensitivity capitalization,
                                          const ExceptionResponse policy) {
  std::vector<MdlMol> result;
  
  // Find the limits for different MDL MOL entries
  const std::vector<int2> mol_entry_limits = findSdfMolEntryLimits(tf);

  // Parse each MDL MOL entry
  const int nsection = mol_entry_limits.size();
  result.reserve(nsection);
  for (int i = 0; i < nsection; i++) {
    result.emplace_back(tf, mol_entry_limits[i].x, mol_entry_limits[i].y, capitalization, policy);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<MdlMol> readStructureDataFile(const std::string &file_name,
                                          const CaseSensitivity capitalization,
                                          const ExceptionResponse policy) {
  const TextFile tf(file_name);
  return readStructureDataFile(tf, capitalization, policy);
}
  
} // namespace structure
} // namespace stormm
