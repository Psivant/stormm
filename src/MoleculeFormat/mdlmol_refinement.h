// -*-c++-*-
#ifndef STORMM_MDLMOL_REFINEMENT_H
#define STORMM_MDLMOL_REFINEMENT_H

#include <vector>
#include <string>
#include "copyright.h"
#include "Namelists/nml_report.h"
#include "Potential/scorecard.h"
#include "Restraints/restraint_apparatus.h"
#include "Synthesis/systemcache.h"
#include "Topology/atomgraph.h"
#include "mdlmol.h"

namespace stormm {
namespace structure {

using energy::ScoreCard;
using namelist::ReportControls;
using restraints::RestraintApparatus;
using synthesis::SystemCache;
using topology::AtomGraph;

/// \brief Add a series of user-defined data items to one or more MDL MOL entries.
///
/// Overloaded:
///   - Customize the data items of a single MDL MOL object
///   - Customize the data items of a list of MDL MOL objects
///
/// \param mol_entry    The molecule to customize
/// \param mol_entries  List of molecules to customize
/// \param label        The system label, passed down from a SystemCache's array, probably
///                     originating in a &files namelist -sys keyword
/// \param ag           Topology describing the molecule 
/// \param ra           Restraint apparatus serving the molecule
/// \param sysc         The cache of systems, to which the list of molecules corresponds in order
/// \param repcon       Diagnostics output control instructions (from a &report namelist in the
///                     input deck)
/// \{
void customizeDataItems(MdlMol *mol_entry, const std::string &label, const AtomGraph &ag,
                        const RestraintApparatus &ra, const ReportControls &repcon);

void customizeDataItems(std::vector<MdlMol> *mol_entries, const SystemCache &sysc,
                        const ReportControls &repcon);
/// \}

/// \brief Update the data items of one of more MDL MOL entries based on prior calculations.
///
/// Overloaded:
///   - Update the data items of a single MDL MOL object
///   - Update the data items of a list of MDL MOL objects
///
/// \param mol_entry     The molecule to update
/// \param mol_entries   List of molecules to update
/// \param sysc          The cache of systems, to which the list of molecules corresponds in order
/// \param nrg           Energy tracking for all systems (the order must match that of the
///                      systems cache)
/// \param system_index  Index of the system of interested (for accessing the proper entries of the
///                      energy tracker, if it contains information on more than one system)
/// \{
void updateDataItemReadouts(MdlMol *mol_entry, const SystemCache &sysc, const ScoreCard &nrg,
                            int system_index = 0);

void updateDataItemReadouts(std::vector<MdlMol> *mol_entries, const SystemCache &sysc,
                            const ScoreCard &nrg);
/// \}

} // namespace structure
} // namespace stormm

#endif
