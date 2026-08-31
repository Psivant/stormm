#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "DataTypes/stormm_vector_types.h"
#include "Debug/watcher.h"
#include "FileManagement/file_listing.h"
#include "Math/vector_ops.h"
#include "Namelists/nml_debug.h"
#include "Namelists/nml_files.h"
#include "Namelists/nml_report.h"
#include "Parsing/parse.h"
#include "Reporting/ordered_list.h"
#include "Reporting/report_table.h"
#include "Reporting/reporting_enumerators.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Synthesis/synthesis_abstracts.h"
#include "Topology/atomgraph.h"
#include "present_debug.h"
#include "section_contents.h"

namespace stormm {
namespace review {

using card::HybridTargetLevel;
using debug::Watcher;
using diskutil::getBaseName;
using namelist::DebugControls;
using namelist::FilesControls;
using namelist::ReportControls;
using parse::char4ToString;
using parse::realToString;
using review::commentSymbol;
using review::ListEnumeration;
using review::OrderedList;
using review::OutputSyntax;
using review::ReportTable;
using stmath::findBin;
using synthesis::PhaseSpaceSynthesis;
using synthesis::SyNonbondedKit;
using topology::AtomGraph;
  
//-------------------------------------------------------------------------------------------------
void createDebugReport(const UserSettings &ui, const DynamicsIntervention &mm_intv) {
  std::vector<SectionContents> output_blocks(1);

  // The first section will provide some overview of what was debugged.
  output_blocks[0].setTitle("Results of debugging activity");
  output_blocks[0].addNarration("Users can add watch requests and other features to the "
                                "simulation which will report apparent anomalies.");
  const DebugControls dbgcon = ui.getDebugNamelistInfo();

  // Subsequent sections will be written on an ad-hoc basis.
  if (mm_intv.hasAnomalyReporting()) {
    
    // Extract a pointer to the simulations' anomaly tracker.
    const Watcher *bugrep = mm_intv.getAnomalyReportingPointer();
    const int nfrc_anomaly = bugrep->getLargeForceCount();
    const std::vector<float4> frc_anomalies = bugrep->getLargeForces();
    const std::vector<int> frc_anom_steps = bugrep->getLargeForceSteps();
    const SyNonbondedKit<double, double2> poly_nbk = mm_intv.getDoublePrecisionNonbondedKit();
    std::vector<int> system_bounds(2 * poly_nbk.nsys);
    for (int i = 0; i < poly_nbk.nsys; i++) {
      system_bounds[(2 * i)    ] = poly_nbk.atom_offsets[i];
      system_bounds[(2 * i) + 1] = poly_nbk.atom_offsets[i] + poly_nbk.atom_counts[i];
    }

    // Prepare a table of all large forces and the steps on which they occurred.  The columns will
    // be as follows:
    //
    // 1.)  Step number
    // 2.)  System Index
    // 3.)  Atom index within the system
    // 4.)  Cartesian X force
    // 5.)  Cartesian Y force
    // 6.)  Cartesian Z force
    // 7.)  Notes
    // 8.)  Atom name
    // 9.)  Residue name containing the atom
    // 10.) Topology file name
    std::vector<std::string> force_data(10 * nfrc_anomaly);
    const PhaseSpaceSynthesis *poly_ps = mm_intv.getPhaseSpaceSynthesisPointer();
    for (int i = 0; i < nfrc_anomaly; i++) {
      const Ecumenical4 idx_intrp = { .f = frc_anomalies[i].w };
      const int twice_sys_idx = findBin(system_bounds, idx_intrp.i);
      const int sys_idx = twice_sys_idx / 2;
      const AtomGraph *i_ag = poly_ps->getSystemTopologyPointer(sys_idx);
      if (twice_sys_idx & 0x1) {
        force_data[(6 * nfrc_anomaly) + i] = std::string("% [! PR !]");
      }
      else {
        force_data[(6 * nfrc_anomaly) + i] = std::string("%");
      }
      force_data[(9 * nfrc_anomaly) + i] = getBaseName(i_ag->getFileName());
      force_data[                     i] = std::to_string(frc_anom_steps[i]);
      force_data[     nfrc_anomaly  + i] = std::to_string(twice_sys_idx);
      const int anom_atom_idx = idx_intrp.i - poly_nbk.atom_offsets[sys_idx];
      force_data[(2 * nfrc_anomaly) + i] = std::to_string(anom_atom_idx);
      if (twice_sys_idx & 0x1) {
        force_data[(7 * nfrc_anomaly) + i] = std::string("    ");
        force_data[(8 * nfrc_anomaly) + i] = std::string("    ");
      }
      else {
        force_data[(7 * nfrc_anomaly) + i] = char4ToString(i_ag->getAtomName(anom_atom_idx));
        const int res_idx = i_ag->getResidueIndex(anom_atom_idx);
        force_data[(8 * nfrc_anomaly) + i] = char4ToString(i_ag->getResidueName(res_idx));
      }
      force_data[(3 * nfrc_anomaly) + i] = realToString(frc_anomalies[i].x, 1);
      force_data[(4 * nfrc_anomaly) + i] = realToString(frc_anomalies[i].y, 1);
      force_data[(5 * nfrc_anomaly) + i] = realToString(frc_anomalies[i].z, 1);
    }
    std::vector<std::string> force_headings = {
      "Step", "System", "Atom Index", "Atom Name", "Residue", "Force in X", "Force in Y",
      "Force in Z", "Notes", "Topology File"
    };
    std::vector<JustifyText> force_justifications(10, JustifyText::RIGHT);
    for (int i = 6; i < 10; i++) {
      force_justifications[i] = JustifyText::LEFT;
    }
    const std::string base_varname = dbgcon.getVariableBase();
    ReportTable force_list(force_data, force_headings, base_varname + "_anomalous_forces", 200,
                           force_justifications);
    force_list.unprotectContent();
    SectionContents result;
    result.setTitle("Anomalous forces observed during the simulation");
    if (nfrc_anomaly > 0) {
      result.addNarration("Forces with magnitude above the threshold of " +
                          realToString(dbgcon.getLargeForceThreshold(), 1) + " kcal/mol-A were "
                          "found in the simulation.  Columns of the following table show:");
      OrderedList force_desc(ListEnumeration::NUMBERED);
      force_desc.addItem("Step number at which the large force was observed");
      force_desc.addItem("System index in which the large force occurred");
      force_desc.addItem("Atom index within the system which was subjected to the large force");
      force_desc.addItem("Name of the atom");
      force_desc.addItem("Name of the residue containing the atom");
      force_desc.addItem("Cartesian X component of the force");
      force_desc.addItem("Cartesian Y component of the force");
      force_desc.addItem("Cartesian Z component of the force");
      force_desc.addItem("Notes:");
      force_desc.addNestedItem("PR: The anomalous force was found in the padding region between "
                               "systems, not on an actual atom.  This may not have caused any "
                               "spurious movement, but typical force routines do not write data "
                               "in this region.");
      force_desc.addItem("Topology file name (directory paths are clipped)");
      result.addList(force_desc);
      result.addTable(force_list);
    }
    else {
      result.addNarration("No anomalous forces exceeding a threshold of " +
                          realToString(dbgcon.getLargeForceThreshold(), 1) + " kcal/mol-A were "
                          "found in the simulation.");
    }
    output_blocks.push_back(result);
  }

  // Produce the report file
  const FilesControls& ficon = ui.getFilesNamelistInfo();
  const ReportControls& repcon = ui.getReportNamelistInfo();
  printAllSections(ficon.getDebugFile(), ui.getPrintingPolicy(), output_blocks,
                   repcon.getOutputSyntax(), ListEnumeration::NUMBERED,
                   ListEnumeration::ALPHABETIC, repcon.getReportFileWidth());
}

} // namespace review
} // namespace stormm
