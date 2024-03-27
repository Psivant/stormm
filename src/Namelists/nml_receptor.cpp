#include "copyright.h"
#include "Constants/symbol_values.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "Structure/mesh_parameters.h"
#include "input.h"
#include "nml_receptor.h"

namespace stormm {
namespace namelist {

using constants::getEnumerationName;
using parse::minimalRealFormat;
using structure::default_mesh_scaling_bits;
using structure::getEnumerationName;
using structure::translateBoundaryCondition;
using structure::translateGridDetail;
using structure::translateMeshPosition;
using structure::translateRMSDMethod;

//-------------------------------------------------------------------------------------------------
ReceptorControls::ReceptorControls(const ExceptionResponse policy_in) :
    policy{policy_in},
    label_group{std::string("")}, align_mask{std::string("")},
    align_method{std::string(default_receptor_alignment_method)},
    mesh_position{std::string(default_receptor_grid_alignment)},
    nml_transcript{"receptor"}
{}

//-------------------------------------------------------------------------------------------------
ReceptorControls::ReceptorControls(const TextFile &tf, int *start_line, bool *found_nml,
                                   const ExceptionResponse policy_in, const WrapTextSearch wrap) :
    ReceptorControls(policy_in)
{
  bool nml_exist;
  const NamelistEmulator t_nml = receptorInput(tf, start_line, &nml_exist, policy, wrap);
  if (found_nml != nullptr) {
    *found_nml = nml_exist;
  }

  // The namelist's existence must be detected in order to proceed.  If the namelist has been
  // provided, it must have a label group.  Otherwise, take the default values for the
  // ReceptorControls and return.
  if (nml_exist == false) {
    return;
  }
  nml_transcript = t_nml;
  if (t_nml.getKeywordStatus("label_group") == InputStatus::MISSING) {

    // Always fatal: the receptor grid must have a structure to operate on, and a distinct label
    // is the only way to specify that structure.
    rtErr("A label group must be provided in order to construct a mesh based on a rigid or "
          "semi-rigid receptor molecule.", "ReceptorControls");
  }
  else {
    label_group = t_nml.getStringValue("label_group");
  }
  t_nml.assignVariable(&align_mask, "alignment_mask");
  align_method = t_nml.getStringValue("alignment_method");
  setMeshPosition(t_nml.getStringValue("mesh_position"));
}

//-------------------------------------------------------------------------------------------------
const std::string& ReceptorControls::getLabelGroup() const {
  return label_group;
}

//-------------------------------------------------------------------------------------------------
const std::string& ReceptorControls::getAlignmentMask() const {
  return align_mask;
}

//-------------------------------------------------------------------------------------------------
RMSDMethod ReceptorControls::getAlignmentMethod() const {
  return translateRMSDMethod(align_method);
}
  
//-------------------------------------------------------------------------------------------------
MeshPosition ReceptorControls::getMeshPosition() const {
  return translateMeshPosition(mesh_position);
}

//-------------------------------------------------------------------------------------------------
const NamelistEmulator& ReceptorControls::getTranscript() const {
  return nml_transcript;
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setLabelGroup(const std::string &label_group_in) {
  label_group = label_group_in;
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setMeshPosition(const std::string &mesh_position_in) {
  mesh_position = mesh_position_in;
  try {
    const MeshPosition interp = translateMeshPosition(mesh_position);
  }
  catch (std::runtime_error) {
    switch (policy) {
    case ExceptionResponse::DIE:
      rtErr("An invalid mesh positioning " + mesh_position + " was specified.", "ReceptorControls",
            "setMeshPosition");
    case ExceptionResponse::WARN:
      rtWarn("An invalid mesh positioning " + mesh_position + " was specified.  The default of " +
             std::string(default_receptor_grid_alignment) + " will be restored.",
             "ReceptorControls", "setMeshPosition");
      break;
    case ExceptionResponse::SILENT:
      break;
    }
    mesh_position = std::string(default_receptor_grid_alignment);
  }
}

//-------------------------------------------------------------------------------------------------
void ReceptorControls::setMeshPosition(const MeshPosition mesh_position_in) {
  mesh_position = getEnumerationName(mesh_position_in);
}

//-------------------------------------------------------------------------------------------------
NamelistEmulator receptorInput(const TextFile &tf, int *start_line, bool *found,
                               const ExceptionResponse policy, const WrapTextSearch wrap) {
  NamelistEmulator t_nml("receptor", CaseSensitivity::AUTOMATIC, policy, "Collects details of a "
                         "receptor density field for docking and conformer selection in STORMM.");

  // Keyword: label_group
  t_nml.addKeyword("label_group", NamelistType::STRING, std::string(""));
  t_nml.addHelp("label_group", "Files management is confined to the &files namelist by specifying "
                "a label group to define the coordinates and topology (or topologies) composing "
                "the receptor.");

  // Keyword: alignment_mask
  t_nml.addKeyword("alignment_mask", NamelistType::STRING, std::string(""));
  t_nml.addHelp("alignment_mask", "When multiple structures are provided to compose the mesh, a "
                "method for aligning them is critical.  This atom mask indicates which atoms and "
                "residues of each molecule shall be compared in order to bring all structures of "
                "the receptor into the most relevant alignment.");

  // Keyword: alignment_method
  t_nml.addKeyword("alignment_method", NamelistType::STRING,
                   std::string(default_receptor_alignment_method));
  t_nml.addHelp("alignment_method", "Multiple receptor structures should be aligned to one "
                "another in order to provide the best focus on a particular binding site.  Use "
                "this keyword to indicate how to position each structure relative to an average "
                "set of positions for atoms described by alignment_mask.");
  
  // Keyword: mesh_position
  t_nml.addKeyword("mesh_position", NamelistType::STRING,
                   std::string(default_receptor_grid_alignment));
  t_nml.addHelp("mesh_position", "The manner in which the receptor's potential mesh shall be "
                "aligned to the receptor molecule itself.  Aligning to the molecule centers the "
                "mesh on the receptor with equal space (or overhang) between the molecule's atoms "
                "and the planar faces of the mesh.  Options include MOLECULE or STRUCTURE (center "
                "on the molecule), ORIGIN (set the origin of the mesh to the origin of the "
                "receptor's Cartesian coordinate system), and ARBITRARY (set the origin of the "
                "mesh to some user-specified point in space).");

  // Search the input file, read the namelist if it can be found, and update the current line
  // for subsequent calls to this function or other namelists.
  *start_line = readNamelist(tf, &t_nml, *start_line, wrap, tf.getLineCount(), found);
  return t_nml;
}
  
} // namespace namelist
} // namespace stormm
