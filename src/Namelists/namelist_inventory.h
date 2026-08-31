// -*-c++-*-
#ifndef STORMM_NAMELIST_INVENTORY_H
#define STORMM_NAMELIST_INVENTORY_H

#include <string>
#include <vector>
#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/parsing_enumerators.h"
#include "Parsing/textfile.h"
#include "namelist_emulator.h"
#include "nml_analysis.h"
#include "nml_conformer.h"
#include "nml_debug.h"
#include "nml_dynamics.h"
#include "nml_emulate.h"
#include "nml_ffmorph.h"
#include "nml_files.h"
#include "nml_mesh.h"
#include "nml_minimize.h"
#include "nml_nice.h"
#include "nml_pppm.h"
#include "nml_precision.h"
#include "nml_random.h"
#include "nml_receptor.h"
#include "nml_remd.h"
#include "nml_report.h"
#include "nml_restraint.h"
#include "nml_scene.h"
#include "nml_solvent.h"

namespace stormm {
namespace namelist {

using constants::ExceptionResponse;
using parse::WrapTextSearch;
using parse::TextFile;
using NmlFuncPtr = NamelistEmulator(*)(const TextFile&, int*, bool*, const ExceptionResponse,
                                       WrapTextSearch);
  
/// \brief Link a string containing the human-readable title of a namelist to a function pointer
///        which produces a complete, but probably uninitialized, namelist emulator containing all
///        of the relevant keywords.
class NamelistToken {
public:

  /// \brief The constructor will take both the namelist title and a pointer to the function that
  ///        produces its NamelistEmulator.
  NamelistToken(const std::string &title_in, NmlFuncPtr producer_in);

  /// \brief Get the namelist title
  const std::string& getTitle() const;

  /// \brief Invoke the corresponding function to produce a namelist emulator, perhaps loaded with
  ///        input data.
  NamelistEmulator invoke(const TextFile &tf, int *start_line, bool *found,
                          ExceptionResponse policy = ExceptionResponse::DIE,
                          WrapTextSearch wrap = WrapTextSearch::YES) const;
  
private:
  std::string title;    ///< Title of the namelist, e.g. &files
  NmlFuncPtr producer;  ///< Function producing a relevant example of the namelist loaded with
                        ///<   defaults or even actual user data, e.g. filesInput()
};

const std::vector<NamelistToken> namelist_inventory = {
  NamelistToken(std::string("&analysis"), analysisInput),
  NamelistToken(std::string("&conformer"), conformerInput),
  NamelistToken(std::string("&debug"), debugInput),
  NamelistToken(std::string("&dynamics"), dynamicsInput),
  NamelistToken(std::string("&emulator"), emulatorInput),
  NamelistToken(std::string("&ffmorph"), ffmorphInput),
  NamelistToken(std::string("&files"), filesInput),
  NamelistToken(std::string("&mesh"), meshInput),
  NamelistToken(std::string("&minimize"), minimizeInput),
  NamelistToken(std::string("&nice"), niceInput),
  NamelistToken(std::string("&pppm"), pppmInput),
  NamelistToken(std::string("&precision"), precisionInput),
  NamelistToken(std::string("&random"), randomInput),
  NamelistToken(std::string("&receptor"), receptorInput),
  NamelistToken(std::string("&remd"), remdInput),
  NamelistToken(std::string("&report"), reportInput),
  NamelistToken(std::string("&restraint"), restraintInput),
  NamelistToken(std::string("&scene"), sceneInput),
  NamelistToken(std::string("&solvent"), solventInput)
};

} // namespace namelist
} // namespace stormm

#endif
