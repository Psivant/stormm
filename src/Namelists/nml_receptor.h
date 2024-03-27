// -*-c++-*-
#ifndef STORMM_NML_RECEPTOR_H
#define STORMM_NML_RECEPTOR_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Constants/symbol_values.h"
#include "Parsing/textfile.h"
#include "Structure/structure_enumerators.h"
#include "namelist_element.h"
#include "namelist_emulator.h"

namespace stormm {
namespace namelist {

using constants::CartesianDimension;
using constants::ExceptionResponse;
using constants::UnitCellAxis;
using parse::TextFile;
using parse::WrapTextSearch;
using structure::BoundaryCondition;
using structure::GridDetail;
using structure::MeshPosition;
using structure::RMSDMethod;

/// \brief Default parameters for the &receptor namelist
constexpr char default_receptor_grid_alignment[] = "molecule";
constexpr char default_receptor_alignment_method[] = "align_mass";

/// \brief Encapsulate the data extracted from a &receptor namelist to define a grid-mapped
///        representation of a rigid macromolecular structure.
///
class ReceptorControls {
public:
  
  /// \brief The constructor can prepare an object with default settings or read the corresponding
  ///        namelist to accept user input.
  ///
  /// \param tf          Input file translated into RAM
  /// \param start_line  Line of the input file to begin searching for the &solvent namelist
  /// \param found_nml   Indication of whether the namelist was found in the input file
  /// \param policy_in   Requested error handling behavior
  /// \param wrap        Indicate that the search for a &conformer namelist should carry on from
  ///                    the beginning of an input file if no such namelist is found starting
  ///                    from the original starting point
  /// \{
  ReceptorControls(ExceptionResponse policy_in = ExceptionResponse::DIE);

  ReceptorControls(const TextFile &tf, int *start_line, bool *found_nml,
                   ExceptionResponse policy_in = ExceptionResponse::DIE,
                   WrapTextSearch wrap = WrapTextSearch::NO);
  /// \}

  /// \brief As with other control objects, copy and move constructors, plus copy and move
  ///        assignment operators, can all take their default forms.
  /// \{
  ReceptorControls(const ReceptorControls &original) = default;
  ReceptorControls(ReceptorControls &&original) = default;
  ReceptorControls& operator=(const ReceptorControls &original) = default;
  ReceptorControls& operator=(ReceptorControls &&original) = default;
  /// \}

  /// \brief Get the label group in which to find all systems composing the receptor grid.  This
  ///        will correspond to a label from the &files namelist.
  const std::string& getLabelGroup() const;

  /// \brief Get the alignment mask showing which atoms and residues are most critical to the
  ///        alignment of multiple structures if more than one receptor structure is at hand.
  const std::string& getAlignmentMask() const;

  /// \brief Get the alignment method.
  RMSDMethod getAlignmentMethod() const;
  
  /// \brief Get the manner in which the mesh is aligned to the rigid molecule it represents.
  MeshPosition getMeshPosition() const;

  /// \brief Get the original namelist emulator object as a transcript of the user input.
  const NamelistEmulator& getTranscript() const;

  /// \brief Set the label group for structures to be used in composing the mesh.
  void setLabelGroup(const std::string &label_group_in);

  /// \brief Set the positioning of the mesh relative to the receptor molecule is describes.
  ///        Overloading in this function follows from setMeshPotential() above.
  ///
  /// \param alignment_in  The chosen mesh positioning
  /// \{
  void setMeshPosition(const std::string &alignment_in);
  void setMeshPosition(MeshPosition alignment_in);
  /// \}
  
private:
  ExceptionResponse policy;   ///< The course to take when encountering bad input
  std::string label_group;    ///< The label group from which to draw structures for the mesh.
                              ///<   Meshes expressing OCCLUSION potentials can draw from more
                              ///<   than one structure, but others may require single structures
                              ///<   or more strenuous approximations in order to draw meaningful
                              ///<   maps of receptor ensembles.
  std::string align_mask;     ///< Atom mask for aligning multiple structures
  std::string align_method;   ///< Chosen method for aligning multiple receptor structures
                              ///<   according to the atom mask given in align_mask
  std::string mesh_position;  ///< The way in which the mesh aligns to the molecule it represents

  /// Store a deep copy of the original namelist emulator as read from the input file.
  NamelistEmulator nml_transcript;
};

/// \brief Free function to read the &receptor namelist.
///
/// \param tf          Text of file containing the input deck, read into RAM
/// \param start_line  Line of the input file at which to begin the scan
/// \param found       Indicator that the namelist was found in the input file
/// \param policy      Response to bad inputs
/// \param wrap        Indicate that the search for a &conformer namelist should carry on from the
///                    beginning of an input file if no such namelist is found starting from the
///                    original starting point
NamelistEmulator receptorInput(const TextFile &tf, int *start_line, bool *found,
                               ExceptionResponse policy = ExceptionResponse::DIE,
                               WrapTextSearch wrap = WrapTextSearch::NO);

} // namespace namelist
} // namespace stormm

#endif
