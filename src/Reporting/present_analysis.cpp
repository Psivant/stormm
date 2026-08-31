#include "copyright.h"
#include "Analysis/hydrogen_bond_analysis.h"
#include "Namelists/nml_files.h"
#include "Namelists/nml_report.h"
#include "present_analysis.h"
#include "section_contents.h"

namespace stormm {
namespace review {

using analysis::HydrogenBondAnalysis;
using namelist::FilesControls;
using namelist::ReportControls;

//-------------------------------------------------------------------------------------------------
void createAnalysisReport(const UserSettings &ui, const DynamicsIntervention &mm_intv) {
  std::vector<SectionContents> output_blocks;
  const int n_hba = mm_intv.getHydrogenBondAnalysisCount();
  output_blocks.reserve(n_hba);
  
  // Create the writeup of each hydrogen bonding analysis
  const ReportControls& repcon = ui.getReportNamelistInfo();
  for (int i = 0; i < n_hba; i++) {
    const HydrogenBondAnalysis *hba_report = dyna_tk.getHydrogenBondAnalysisPointer(i);
    output_blocks.push_back(hba_report->reportResults(repcon));
  }

  // Produce the report file
  const FilesControls& ficon = ui.getFilesNamelistInfo();
  printAllSections(ficon.getAnalysisFile(), ui.getPrintingPolicy(), output_blocks,
                   repcon.getOutputSyntax(), ListEnumeration::NUMBERED,
                   ListEnumeration::ALPHABETIC, repcon.getReportFileWidth());
}

} // namespace review
} // namespace stormm
