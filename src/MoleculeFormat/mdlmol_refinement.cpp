#include "copyright.h"
#include "Constants/behavior.h"
#include "Parsing/parse.h"
#include "Parsing/polynumeric.h"
#include "mdlmol_refinement.h"

namespace stormm {
namespace structure {

using constants::CaseSensitivity;
using parse::NumberFormat;
using parse::realToString;
using parse::strcmpCased;
  
//-------------------------------------------------------------------------------------------------
void customizeDataItems(MdlMol *mol_entry, const std::string &label, const AtomGraph &ag,
                        const RestraintApparatus &ra, const ReportControls &repcon) {
  const int nreq = repcon.getSDFileDataRequestCount();
  const std::vector<MdlMolDataRequest>& reqs = repcon.getSDFileDataRequests();
  for (int i = 0; i < nreq; i++) {
    if (reqs[i].getSystemLabel() == label ||
        strcmpCased(reqs[i].getSystemLabel(), "ALL", CaseSensitivity::NO)) {
      mol_entry->addDataItem(reqs[i], ag, ra);
    }
  }
}

//-------------------------------------------------------------------------------------------------
void customizeDataItems(std::vector<MdlMol> *mol_entries, const SystemCache &sysc,
                        const ReportControls &repcon) {
  const int nmol = mol_entries->size();
  MdlMol* ment_data = mol_entries->data();
  for (int i = 0; i < nmol; i++) {
    customizeDataItems(&ment_data[i], sysc.getSystemLabel(i), sysc.getSystemTopology(i),
                       sysc.getRestraints(i), repcon);
  }
}

//-------------------------------------------------------------------------------------------------
void updateDataItemReadouts(MdlMol *mol_entry, const SystemCache &sysc, const ScoreCard &nrg,
                            const int system_index) {
  const int nitems = mol_entry->getDataItemCount();
  for (int i = 0; i < nitems; i++) {
    switch (mol_entry->getDataItemKind(i)) {
    case MdlMolDataItemKind::STATE_VARIABLE:
      {
        const StateVariable detail = mol_entry->getTrackedState(i);
        double dt_val;
        if (detail == StateVariable::POTENTIAL_ENERGY) {
          dt_val = nrg.reportPotentialEnergy(system_index);
        }
        else if (detail == StateVariable::TOTAL_ENERGY) {
          dt_val = nrg.reportTotalEnergy(system_index);
        }
        else {
          dt_val = nrg.reportInstantaneousStates(detail, system_index);
        }
        mol_entry->addLineToDataItem(realToString(dt_val, 14, 4, NumberFormat::STANDARD_REAL), i);
      }
      break;
    case MdlMolDataItemKind::ATOM_INFLUENCES:
      break;
    case MdlMolDataItemKind::TOPOLOGY_PARAMETER:
    case MdlMolDataItemKind::STRING:
    case MdlMolDataItemKind::NATIVE:
    case MdlMolDataItemKind::NONE:
      break;
    }
  }
}
  
//-------------------------------------------------------------------------------------------------
void updateDataItemReadouts(std::vector<MdlMol> *mol_entries, const SystemCache &sysc,
                            const ScoreCard &nrg) {
  const int nmol = mol_entries->size();
  MdlMol* ment_data = mol_entries->data();
  for (int i = 0; i < nmol; i++) {
    updateDataItemReadouts(&ment_data[i], sysc, nrg, i);
  }
}
  
} // namespace structure
} // namespace stormm
