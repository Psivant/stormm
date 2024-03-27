#include "copyright.h"
#include "FileManagement/file_listing.h"
#include "Math/vector_ops.h"
#include "MoleculeFormat/mdlmol.h"
#include "Parsing/textfile.h"
#include "Namelists/nml_files.h"
#include "Trajectory/trajectory_enumerators.h"
#include "test_system_manager.h"
#include "test_environment.h"

namespace stormm {
namespace testing {

using constants::getEnumerationName;
using diskutil::DrivePathType;
using diskutil::getDefaultFileExtension;
using diskutil::getDrivePathType;
using diskutil::getBaseName;
using diskutil::osSeparator;
using diskutil::splitPath;
using namelist::FilesControls;
using parse::TextFile;
using parse::TextOrigin;
using stmath::minValue;
using stmath::maxValue;
using structure::MdlMol;
using trajectory::CoordinateFileKind;
using trajectory::translateCoordinateFileKind;
using trajectory::detectCoordinateFileKind;

//-------------------------------------------------------------------------------------------------
TestSystemManager::TestSystemManager() :
    topology_base{std::string("")}, topology_extn{std::string("")},
    coordinate_base{std::string("")}, coordinate_extn{std::string("")}, 
    fault_response{ExceptionResponse::WARN}, fault_found{false},
    topology_names{}, coordinate_names{}, topology_exist{}, topology_success{}, coordinate_exist{},
    coordinate_success{}, compatibility{}, all_topologies{}, all_coordinates{}
{}
  
//-------------------------------------------------------------------------------------------------
TestSystemManager::TestSystemManager(const std::string &topology_base_in,
                                     const std::string &topology_extn_in,
                                     const std::vector<std::string> &topology_names_in,
                                     const std::string &coordinate_base_in,
                                     const std::string &coordinate_extn_in,
                                     const std::vector<std::string> &coordinate_names_in,
                                     const ExceptionResponse policy,
                                     const TestPriority fault_response_in,
                                     const TestPriority all_go_response_in) :
    TestSystemManager()
{
  topology_base = topology_base_in;
  topology_extn = topology_extn_in;
  coordinate_base = coordinate_base_in;
  coordinate_extn = coordinate_extn_in;
  const size_t n_tops = topology_names_in.size();
  const size_t n_crds = coordinate_names_in.size();
  if (n_tops != n_crds) {
    rtErr("The number of topology (" + std::to_string(n_tops) + ") and coordinate (" +
          std::to_string(n_crds) + ") file names must be equal.", "TestSystemManager");
  }
  system_count = n_tops;
  topology_names.reserve(system_count);
  coordinate_names.reserve(system_count);
  const char osc = osSeparator();
  for (int i = 0; i < system_count; i++) {
    topology_names.emplace_back(topology_base_in + osc + topology_names_in[i] + '.' +
                                topology_extn_in);
    coordinate_names.emplace_back(coordinate_base_in + osc + coordinate_names_in[i] + '.' +
                                  coordinate_extn_in);
  }
  fault_response = fault_response_in;
  all_go_response = all_go_response_in;
  topology_exist.resize(system_count);
  topology_success.resize(system_count);
  coordinate_exist.resize(system_count);
  coordinate_success.resize(system_count);
  compatibility.resize(system_count, false);
  all_topologies.reserve(system_count);
  for (int i = 0; i < system_count; i++) {
    topology_exist[i] = (getDrivePathType(topology_names[i]) == DrivePathType::FILE);
    coordinate_exist[i] = (getDrivePathType(coordinate_names[i]) == DrivePathType::FILE);
    if (topology_exist[i]) {
      try {
        all_topologies.emplace_back(topology_names[i], policy);
        topology_success[i] = true;
      }
      catch (std::runtime_error) {
        all_topologies.emplace_back();
        topology_success[i] = false;
        fault_found = true;
      }
    }
    else {
      all_topologies.emplace_back();
      topology_success[i] = false;
      fault_found = true;
    }
    bool coordinates_valid = coordinate_exist[i];
    if (coordinate_exist[i]) {
      try {
        all_coordinates.emplace_back(coordinate_names[i]);
      }
      catch (std::runtime_error) {

        // If the coordinate reading was unsuccessful but the topology was read, try again using
        // the topology to help guide the interpretation.
        if (topology_success[i]) {
          try {
            all_coordinates.emplace_back(all_topologies[i].getAtomCount(),
                                         all_topologies[i].getUnitCellType());
            all_coordinates[i].buildFromFile(coordinate_names[i]);
          }
          catch (std::runtime_error) {

            // If the format is MDL MOL or a BIOVIA SD file, create an MdlMol object and export the
            // PhaseSpace from that.  As with trajectory files, only the first system will be read.
            if (detectCoordinateFileKind(coordinate_names[i]) == CoordinateFileKind::SDF) {
              try {
                const MdlMol tmp_mdl(coordinate_names[i]);
                all_coordinates.push_back(tmp_mdl.exportPhaseSpace());
              }
              catch (std::runtime_error) {
                coordinates_valid = false;
              }
            }
            else {
              coordinates_valid = false;
            }
          }
        }
      }
    }
    if (coordinates_valid) {
      coordinate_success[i] = true;
    }
    else {
      all_coordinates.emplace_back();
      coordinate_success[i] = false;
      fault_found = true;
    }
    if (topology_success[i] && coordinate_success[i]) {
      if (all_topologies[i].getAtomCount() == all_coordinates[i].getAtomCount()) {
        compatibility[i] = true;
      }
      else {
        fault_found = true;
      }
    }
  }

  // Report errors
  if (fault_found) {
    std::string file_errors;
    for (int i = 0; i < system_count; i++) {
      if (topology_exist[i] == false) {
        file_errors += "  - " + getBaseName(topology_names[i]) + " (does not exist)\n";
      }
      else if (topology_success[i] == false) {
        file_errors += "  - " + getBaseName(topology_names[i]) + " (unreadable)\n";
      }
      if (coordinate_exist[i] == false) {
        file_errors += "  - " + getBaseName(coordinate_names[i]) + " (does not exist)\n";
      }
      else if (coordinate_success[i] == false) {
        file_errors += "  - " + getBaseName(coordinate_names[i]) + " (unreadable)\n";
      }
      if (compatibility[i] == false) {
        file_errors += "  - " + getBaseName(coordinate_names[i]) + " does not appear to be "
                       "compatible with " + getBaseName(topology_names[i]) + "\n";
      }
    }
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("Some critical files were not found or were unreadable.  These include:\n" +
            file_errors + "This is unrecoverable.", "TestSystemManager");
    case ExceptionResponse::WARN:
      rtWarn("Some critical files were not found or were unreadable.  These include:\n" +
             file_errors + "Subsequent tests will be handled with priority " +
             getEnumerationName(getTestingStatus()) + ".", "TestSystemManager");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
  }
}

//-------------------------------------------------------------------------------------------------
TestSystemManager::TestSystemManager(const std::string &topology_base_in,
                                     const std::vector<std::string> &topology_names_in,
                                     const std::string &coordinate_base_in,
                                     const std::vector<std::string> &coordinate_names_in,
                                     const ExceptionResponse policy,
                                     const TestPriority fault_response_in,
                                     const TestPriority all_go_response_in) :
    TestSystemManager(topology_base_in, std::string(""), topology_names_in, coordinate_base_in,
                      std::string(""), coordinate_names_in, policy, fault_response_in,
                      all_go_response_in)
{}

//-------------------------------------------------------------------------------------------------
TestSystemManager::TestSystemManager(const std::vector<std::string> &topology_names_in,
                                     const std::vector<std::string> &coordinate_names_in,
                                     const ExceptionResponse policy,
                                     const TestPriority fault_response_in,
                                     const TestPriority all_go_response_in) :
    TestSystemManager(std::string(""), std::string(""), topology_names_in, std::string(""),
                      std::string(""), coordinate_names_in, policy, fault_response_in,
                      all_go_response_in)
{}

//-------------------------------------------------------------------------------------------------
void TestSystemManager::checkIndexing(const int index, const char* caller) const {
  if (index < 0 || index >= system_count) {
    rtErr("Index " + std::to_string(index) + " is invalid for a collection of " +
          std::to_string(system_count) + " systems.", "TestSystemManager", caller);
  }
}

//-------------------------------------------------------------------------------------------------
int TestSystemManager::getSystemCount() const {
  return system_count;
}

//-------------------------------------------------------------------------------------------------
int TestSystemManager::getSystemCount(const UnitCellType query_bc) const {
  int result = 0;
  for (int i = 0; i < system_count; i++) {
    result += (all_topologies[i].getUnitCellType() == query_bc);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
int TestSystemManager::getSystemCount(const std::vector<UnitCellType> &query_bc) const {
  int result = 0;
  std::vector<bool> included(system_count, false);
  const size_t nq = query_bc.size();
  for (size_t i = 0; i < nq; i++) {
    for (int j = 0; j < system_count; j++) {
      if (included[j]) {
        continue;
      }
      if (all_topologies[j].getUnitCellType() == query_bc[i]) {
        result++;
        included[j] = true;
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int> TestSystemManager::getQualifyingSystems(const UnitCellType uc_choice) const {
  int nqual = 0;
  for (int i = 0; i < system_count; i++) {
    nqual += (all_topologies[i].getUnitCellType() == uc_choice);
  }
  std::vector<int> result(nqual);
  nqual = 0;
  for (int i = 0; i < system_count; i++) {
    if (all_topologies[i].getUnitCellType() == uc_choice) {
      result[nqual] = i;
      nqual++;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<int>
TestSystemManager::getQualifyingSystems(const std::vector<UnitCellType> &uc_choice) const {
  std::vector<bool> included(system_count, false);
  int nqual = 0;
  const int nchoice = uc_choice.size();
  for (int i = 0; i < nchoice; i++) {
    for (int j = 0; j < system_count; j++) {
      if (included[j] == false && all_topologies[j].getUnitCellType() == uc_choice[i]) {
        nqual++;
        included[j] = true;
      }
    }
  }
  std::vector<int> result(nqual);
  nqual = 0;
  for (int i = 0; i < system_count; i++) {
    if (included[i]) {
      result[nqual] = i;
      nqual++;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::string TestSystemManager::getTopologyBasePath() const {
  return topology_base;
}

//-------------------------------------------------------------------------------------------------
std::string TestSystemManager::getCoordinateBasePath() const {
  return coordinate_base;
}

//-------------------------------------------------------------------------------------------------
std::string TestSystemManager::getTopologyExtension() const {
  return topology_extn;
}

//-------------------------------------------------------------------------------------------------
std::string TestSystemManager::getCoordinateExtension() const {
  return coordinate_extn;
}

//-------------------------------------------------------------------------------------------------
std::string TestSystemManager::getTopologyFile(const int index) const {
  checkIndexing(index, "getTopologyFile");
  return topology_names[index];
}

//-------------------------------------------------------------------------------------------------
std::string TestSystemManager::getCoordinateFile(const int index) const {
  checkIndexing(index, "getCoordinateFile");
  return coordinate_names[index];
}

//-------------------------------------------------------------------------------------------------
TestPriority TestSystemManager::getTestingStatus() const {
  return (fault_found) ? fault_response : all_go_response;
}

//-------------------------------------------------------------------------------------------------
TestPriority TestSystemManager::getTestingStatus(const int index) const {
  return (topology_success[index] && coordinate_success[index] && compatibility[index]) ?
    all_go_response : fault_response;
}

//-------------------------------------------------------------------------------------------------
TestPriority TestSystemManager::getTestingStatus(const std::vector<int> &indices) const {
  bool all_pass = true;
  const size_t n_idx = indices.size();
  for (size_t i = 0; i < n_idx; i++) {
    all_pass = (all_pass && topology_success[indices[i]] && coordinate_success[indices[i]] &&
                compatibility[indices[i]]);
  }
  return (all_pass) ? all_go_response : fault_response;
}

//-------------------------------------------------------------------------------------------------
CoordinateFrame TestSystemManager::exportCoordinateFrame(const int index) const {
  checkIndexing(index, "exportCoordinateFrame");
  return CoordinateFrame(all_coordinates[index]);
}

//-------------------------------------------------------------------------------------------------
PhaseSpace TestSystemManager::exportPhaseSpace(const int index) const {
  checkIndexing(index, "exportPhaseSpace");
  return all_coordinates[index];
}

//-------------------------------------------------------------------------------------------------
AtomGraph TestSystemManager::exportAtomGraph(const int index) const {
  checkIndexing(index, "exportAtomGraph");
  return all_topologies[index];
}

//-------------------------------------------------------------------------------------------------
PhaseSpace& TestSystemManager::viewCoordinates(const int index) {
  checkIndexing(index, "viewCoordinates");
  return all_coordinates[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<PhaseSpace>& TestSystemManager::viewCoordinates() {
  return all_coordinates;
}

//-------------------------------------------------------------------------------------------------
const AtomGraph& TestSystemManager::getTopologyReference(const int index) const {
  checkIndexing(index, "viewTopology");
  return all_topologies[index];
}

//-------------------------------------------------------------------------------------------------
const std::vector<AtomGraph>& TestSystemManager::getTopologyReference() const {
  return all_topologies;
}

//-------------------------------------------------------------------------------------------------
AtomGraph* TestSystemManager::getTopologyPointer(const int index) {
  checkIndexing(index, "viewTopology");
  return &all_topologies[index];
}

//-------------------------------------------------------------------------------------------------
const AtomGraph* TestSystemManager::getTopologyPointer(const int index) const {
  checkIndexing(index, "viewTopology");
  return &all_topologies[index];
}

//-------------------------------------------------------------------------------------------------
std::vector<AtomGraph*> TestSystemManager::getTopologyPointer() {
  std::vector<AtomGraph*> result(system_count);
  for (int i = 0; i < system_count; i++) {
    result[i] = &all_topologies[i];
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
const std::vector<AtomGraph*> TestSystemManager::getTopologyPointer() const {
  std::vector<AtomGraph*> result(system_count);
  for (int i = 0; i < system_count; i++) {
    result[i] = const_cast<AtomGraph*>(&all_topologies[i]);
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis TestSystemManager::exportPhaseSpaceSynthesis(const std::vector<int> &index_key,
                                                                 const double perturbation_sigma,
                                                                 const int xrs_seed,
                                                                 const int scale_bits,
                                                                 const int vel_bits,
                                                                 const int frc_bits) const {
  std::vector<AtomGraph*> ag_pointers(system_count);
  for (size_t i = 0; i < system_count; i++) {
    ag_pointers[i] = const_cast<AtomGraph*>(&all_topologies[i]);
  }
  PhaseSpaceSynthesis result(all_coordinates, ag_pointers, index_key, scale_bits, 24, vel_bits,
                             frc_bits);
  PsSynthesisWriter resr = result.data();
  Xoshiro256ppGenerator xrs(xrs_seed);
  if (fabs(perturbation_sigma) > constants::tiny) {
    const double rscale = perturbation_sigma * resr.gpos_scale;
    for (int i = 0; i < resr.system_count; i++) {
      const int llim = resr.atom_starts[i];
      const int hlim = resr.atom_starts[i] + resr.atom_counts[i];
      for (int j = llim; j < hlim; j++) {
        const int95_t xpert = hostDoubleToInt95(rscale * xrs.gaussianRandomNumber());
        const int95_t ypert = hostDoubleToInt95(rscale * xrs.gaussianRandomNumber());
        const int95_t zpert = hostDoubleToInt95(rscale * xrs.gaussianRandomNumber());
        const int95_t xnew = hostSplitFPSum(xpert, resr.xcrd[j], resr.xcrd_ovrf[j]);
        const int95_t ynew = hostSplitFPSum(ypert, resr.ycrd[j], resr.ycrd_ovrf[j]);
        const int95_t znew = hostSplitFPSum(zpert, resr.zcrd[j], resr.zcrd_ovrf[j]);
        resr.xcrd[j]      = xnew.x;
        resr.xcrd_ovrf[j] = xnew.y;
        resr.ycrd[j]      = ynew.x;
        resr.ycrd_ovrf[j] = ynew.y;
        resr.zcrd[j]      = znew.x;
        resr.zcrd_ovrf[j] = znew.y;
      }
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis
TestSystemManager::exportPhaseSpaceSynthesis(const UnitCellType uc_choice) const {
  return exportPhaseSpaceSynthesis(getQualifyingSystems(uc_choice));
}

//-------------------------------------------------------------------------------------------------
PhaseSpaceSynthesis
TestSystemManager::exportPhaseSpaceSynthesis(const std::vector<UnitCellType> &uc_choice) const {
  std::vector<int> all_quals;
  for (size_t i = 0; i < uc_choice.size(); i++) {
    const std::vector<int> matches = getQualifyingSystems(uc_choice[i]);
    all_quals.insert(all_quals.end(), matches.begin(), matches.end());
  }
  return exportPhaseSpaceSynthesis(all_quals);
}

//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis
TestSystemManager::exportAtomGraphSynthesis(const std::vector<int> &index_key,
                                            const ExceptionResponse policy) const {
  std::vector<AtomGraph*> ag_pointers(system_count);
  for (size_t i = 0; i < system_count; i++) {
    ag_pointers[i] = const_cast<AtomGraph*>(&all_topologies[i]);
  }
  AtomGraphSynthesis result(ag_pointers, index_key, policy);
  return result;
}

//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis
TestSystemManager::exportAtomGraphSynthesis(const UnitCellType uc_choice,
                                            const ExceptionResponse policy) const {
  return exportAtomGraphSynthesis(getQualifyingSystems(uc_choice));
}

//-------------------------------------------------------------------------------------------------
AtomGraphSynthesis
TestSystemManager::exportAtomGraphSynthesis(const std::vector<UnitCellType> &uc_choice,
                                            const ExceptionResponse policy) const {
  std::vector<int> all_quals;
  for (size_t i = 0; i < uc_choice.size(); i++) {
    const std::vector<int> matches = getQualifyingSystems(uc_choice[i]);
    all_quals.insert(all_quals.end(), matches.begin(), matches.end());
  }
  return exportAtomGraphSynthesis(all_quals);
}

//-------------------------------------------------------------------------------------------------
SystemCache
TestSystemManager::exportSystemCache(const std::vector<int> &index_key, const TestEnvironment &oe,
                                     const std::vector<std::string> &x_name_list,
                                     const std::vector<CoordinateFileKind> &x_type_list,
                                     const std::vector<std::string> &r_name_list,
                                     const std::vector<CoordinateFileKind> &r_type_list,
                                     const std::vector<std::string> &label_list) {
  const char osc = osSeparator();
  const size_t nitem = index_key.size();
  if (minValue(index_key) < 0 || minValue(index_key) >= system_count ||
      maxValue(index_key) < 0 || maxValue(index_key) >= system_count) {
    rtErr("The index range " + std::to_string(minValue(index_key)) + " - " +
          std::to_string(maxValue(index_key)) + " is invalid for a collection of " +
          std::to_string(system_count) + " systems.", "TestSystemManager", "exportSystemCache");
  }
  const std::vector<CoordinateFileKind> actual_x_types =
    (x_type_list.size() > 0) ? x_type_list :
                               std::vector<CoordinateFileKind>(nitem,
                                                               CoordinateFileKind::AMBER_CRD);
  std::vector<std::string> actual_x_names = x_name_list;
  if (x_name_list.size() == 0) {
    actual_x_names.resize(nitem);
    for (size_t i = 0; i < nitem; i++) {
      const std::string base = getBaseName(coordinate_names[index_key[i]]);
      std::string before, after;
      splitPath(base, &before, &after);
      const std::string fext = getDefaultFileExtension(actual_x_types[i]);
      actual_x_names[i] = oe.getTemporaryDirectoryPath() + osc + before + "." + fext;
    }
  }
  else {
    for (size_t i = 0; i < nitem; i++) {
      actual_x_names[i] = oe.getTemporaryDirectoryPath() + osc + actual_x_names[i];
    }
  }
  const std::vector<CoordinateFileKind> actual_r_types =
    (r_type_list.size() > 0) ?
    r_type_list : std::vector<CoordinateFileKind>(nitem, CoordinateFileKind::AMBER_ASCII_RST);
  std::vector<std::string> actual_r_names = r_name_list;
  if (r_name_list.size() == 0) {
    actual_r_names.resize(nitem);
    for (size_t i = 0; i < nitem; i++) {
      const std::string base = getBaseName(coordinate_names[index_key[i]]);
      std::string before, after;
      splitPath(base, &before, &after);
      const std::string fext = getDefaultFileExtension(actual_r_types[i]);
      actual_r_names[i] = oe.getTemporaryDirectoryPath() + osc + before + "." + fext;
    }
  }
  else {
    for (size_t i = 0; i < nitem; i++) {
      actual_r_names[i] = oe.getTemporaryDirectoryPath() + osc + actual_r_names[i];
    }
  }
  std::vector<std::string> actual_label_list = label_list;  
  if (label_list.size() == 0) {
    actual_label_list.resize(nitem);
    for (int i = 0; i < nitem; i++) {
      actual_label_list[i] = "system_" + std::to_string(i);
    }
  }
  if (nitem != actual_x_names.size() || nitem != actual_x_types.size() ||
      nitem != actual_r_names.size() || nitem != actual_r_types.size() ||
      nitem != actual_label_list.size()) {
    rtErr("If provided, the number of trajectory and restart file names and types, as well as the "
          "number of labels, must equal the number of items to be placed in the cache.",
          "TestSystemManager", "exportSystemCache");
  }
  std::string fcon_str("&files\n");
  for (int i = 0; i < nitem; i++) {
    fcon_str += "  -sys { -p " + topology_names[index_key[i]] + "\n         -c " +
                coordinate_names[index_key[i]] + "\n         -x " + actual_x_names[index_key[i]] +
                "\n         -r " + actual_r_names[index_key[i]] + "}\n";
  }
  fcon_str += "&end\n";
  return exportSystemCache(fcon_str);
}

//-------------------------------------------------------------------------------------------------
SystemCache TestSystemManager::exportSystemCache(const std::string &fnml) {
  TextFile ftf(fnml, TextOrigin::RAM);
  int start_line = 0;
  FilesControls fcon(ftf, &start_line);
  return SystemCache(fcon);
}

} // namespace testing
} // namespace stormm
