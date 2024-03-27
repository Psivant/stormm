#include "copyright.h"
#include "Constants/scaling.h"
#include "FileManagement/file_listing.h"
#include "FileManagement/file_util.h"
#include "Math/summation.h"
#include "Math/statistical_enumerators.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "Topology/atomgraph_enumerators.h"
#include "present_energy.h"
#include "section_contents.h"
#include "summary_file.h"

namespace stormm {
namespace review {

using constants::verytiny;
using diskutil::extractCommonPaths;
using diskutil::listCommonPaths;
using diskutil::openOutputFile;
using energy::getEnumerationName;
using parse::lowercase;
using stmath::sum;
using stmath::colorVectorMask;
using stmath::elementwiseMultiply;
using stmath::mean;
using stmath::variance;
using stmath::VarianceMethod;
using synthesis::SystemCache;
using synthesis::SynthesisMapReader;
using topology::AtomGraph;
using topology::ImplicitSolventModel;
using topology::TorsionKind;
using topology::NonbondedKit;
using topology::UnitCellType;
using topology::ValenceKit;

//-------------------------------------------------------------------------------------------------
std::vector<StateVariable> assessEnergyDetails(const ScoreCard &nrg,
                                               const SynthesisCacheMap &scmap,
                                               const std::vector<StateVariable> &quantities) {
  if (quantities.size() > 0) {
    return quantities;
  }
  else {

    // Determine which cache systems and cache topologies are in use
    const SynthesisMapReader scmapr = scmap.data();
    std::vector<bool> cache_system_in_use(scmapr.ncache, false);
    std::vector<bool> cache_topol_in_use(scmapr.ntopol, false);
    for (int i = 0; i < scmapr.nsynth; i++) {
      cache_system_in_use[scmapr.cache_origins[i]] = true;
      cache_topol_in_use[scmapr.topology_origins[i]] = true;
    }

    // If the systems exist in isolate dboundary conditions, pressure and volume have no meaning.
    bool has_prss, has_volm;
    switch (scmap.getCoordinateSynthesisPointer()->getUnitCellType()) {
    case UnitCellType::NONE:
      has_prss = false;
      has_volm = false;
      break;
    case UnitCellType::ORTHORHOMBIC:
    case UnitCellType::TRICLINIC:
      has_prss = true;
      has_volm = true;
      break;
    }
    
    // If no details of the desired energies are indicated, assume that all relevant energies
    // should be detected and reported.
    const SystemCache *sc = scmap.getCachePointer();
    const std::vector<AtomGraph*> topologies = sc->getTopologyPointer();
    const int ntop = scmapr.ntopol;
    bool has_bond = false;
    bool has_angl = false;
    bool has_dihe = false;
    bool has_impr = false;
    bool has_ubrd = false;
    bool has_cimp = false;
    bool has_cmap = false;
    bool has_qq14 = false;
    bool has_lj14 = false;
    bool has_qqnb = false;
    const bool has_ljnb = true;
    bool has_gb = false;
    std::vector<StateVariable> result;
    for (int i = 0; i < scmapr.ntopol; i++) {
      if (cache_topol_in_use[i] == false) {
        continue;
      }

      // Check for basic valence energy terms
      const ValenceKit<double> vk = topologies[i]->getDoublePrecisionValenceKit();
      has_bond = (has_bond || vk.nbond > 0);
      has_angl = (has_angl || vk.nangl > 0);
      bool found_14 = false;
      for (int j = 0; j < vk.ndihe; j++) {
        switch (static_cast<TorsionKind>(vk.dihe_modifiers[j].w)) {
        case TorsionKind::PROPER:
          has_dihe = true;
          found_14 = true;
          break;
        case TorsionKind::PROPER_NO_14:
          has_dihe = true;
          break;
        case TorsionKind::IMPROPER:
        case TorsionKind::IMPROPER_NO_14:
          has_impr = true;
          break;
        }
      }

      // Check for quantities relevant to CHARMM forcefields.
      has_ubrd = (has_ubrd || vk.nubrd > 0);
      has_cimp = (has_cimp || vk.ncimp > 0);
      has_cmap = (has_cmap || vk.ncmap > 0);
      
      // Check for non-bonded properties
      const NonbondedKit<double> nbk = topologies[i]->getDoublePrecisionNonbondedKit();
      const bool has_charge = (variance(nbk.charge, nbk.natom,
                                        VarianceMethod::VARIANCE) > verytiny);
      has_gb = (has_gb || topologies[i]->getImplicitSolventModel() != ImplicitSolventModel::NONE);
      has_qqnb = (has_qqnb || has_charge);

      // The presence of 1:4 terms and presence of nonzero charges will be taken to imply the
      // relevance of 1:4 electrostatic energies.  The mere presence of 1:4 terms will imply the
      // relevance of 1:4 van-der Waals energies (all systems are assumed to have some van-der
      // Waals properties).
      has_qq14 = (has_qq14 || found_14 && has_charge);
      has_lj14 = (has_lj14 || found_14);
    }

    // Check for restraints
    bool has_rstr = false;
    for (int i = 0; i < scmapr.ncache; i++) {
      if (cache_system_in_use[i] == false) {
        continue;
      }
      has_rstr = (has_rstr || (sc->getRestraintPointer(i)->getTotalRestraintCount() > 0));
    }

    // Assemble the list of applicable energy terms
    if (has_bond) result.push_back(StateVariable::BOND);
    if (has_angl) result.push_back(StateVariable::ANGLE);
    if (has_dihe) result.push_back(StateVariable::PROPER_DIHEDRAL);
    if (has_impr) result.push_back(StateVariable::IMPROPER_DIHEDRAL);
    if (has_ubrd) result.push_back(StateVariable::UREY_BRADLEY);
    if (has_cimp) result.push_back(StateVariable::CHARMM_IMPROPER);
    if (has_cmap) result.push_back(StateVariable::CMAP);
    if (has_qq14) result.push_back(StateVariable::ELEC_ONE_FOUR);
    if (has_lj14) result.push_back(StateVariable::VDW_ONE_FOUR);
    if (has_qqnb) result.push_back(StateVariable::ELECTROSTATIC);
    if (has_ljnb) result.push_back(StateVariable::VDW);
    if (has_gb)   result.push_back(StateVariable::GENERALIZED_BORN);
    if (has_rstr) result.push_back(StateVariable::RESTRAINT);
    if (has_prss) result.push_back(StateVariable::PRESSURE);
    if (has_volm) result.push_back(StateVariable::VOLUME);

    // Check all systems for the presence of nonzero kientic energy
    if (mean(nrg.reportAverageStates(StateVariable::KINETIC)) > verytiny) {
      result.push_back(StateVariable::KINETIC);
      result.push_back(StateVariable::TEMPERATURE_ALL);
    }
    return result;
  }
  __builtin_unreachable();
}

//-------------------------------------------------------------------------------------------------
std::vector<ReportTable> tabulateAverageEnergy(const ScoreCard &nrg, const EnergySample measure,
                                               const std::string &varname,
                                               const SynthesisCacheMap &scmap,
                                               const std::vector<StateVariable> &quantities,
                                               const ReportControls &repcon) {
  const std::vector<StateVariable> used_quantities = assessEnergyDetails(nrg, scmap, quantities);
  const int nquant = used_quantities.size();
  std::vector<ReportTable> result;
  result.reserve(nquant);
  const VarianceMethod stdev = VarianceMethod::STANDARD_DEVIATION;
  for (int i = 0; i < nquant; i++) {
    switch (measure) {
    case EnergySample::TIME_SERIES:
      {
        const int npts = nrg.getSampleSize();
        const int nsys = nrg.getSystemCount();
        std::vector<double2> iresult(npts, { 0.0, 0.0 });
        std::vector<double> table_data(static_cast<size_t>(npts) *
                                       static_cast<size_t>((2 * nsys) + 1));
        for (int j = 0; j < npts; j++) {
          table_data[j] = nrg.getTimeStep(j);
        }
        size_t adcon = npts;
        for (int j = 0; j < nsys; j++) {
          const std::vector<double> ehist = nrg.reportHistory(used_quantities[i], j);
          for (int k = 0; k < npts; k++) {
            iresult[k].x += ehist[k];
            iresult[k].y += ehist[k] * ehist[k];
          }
        }
        const double dnsys = nsys;
        for (int j = 0; j < npts; j++) {
          iresult[j].y = sqrt((dnsys * iresult[i].y) - (iresult[i].x * iresult[i].x)) /
                         sqrt(dnsys * (dnsys - 1.0));
          iresult[j].x /= dnsys;
        }
        for (int j = 0; j < npts; j++) {
          table_data[adcon] = iresult[j].x;
          adcon++;
        }
        for (int j = 0; j < npts; j++) {
          table_data[adcon] = iresult[j].y;
          adcon++;
        }
        std::vector<int> table_decimals((2 * nsys) + 1, repcon.getEnergyDecimalPlaces());
        std::vector<std::string> table_headings;
        table_headings.reserve((2 * nsys) + 1);
        table_headings.push_back("Step Number");
        table_decimals[0] = 0;
        result.emplace_back(table_data, table_headings, table_decimals,
                            varname + "_" + lowercase(getEnumerationName(used_quantities[i])),
                            repcon.getReportFileWidth());
      }
      break;
    case EnergySample::FINAL:
    case EnergySample::TIME_AVERAGE:
      {
        const std::vector<double> equant = (measure == EnergySample::FINAL) ?
                                           nrg.reportInstantaneousStates(used_quantities[i]) :
                                           nrg.reportAverageStates(used_quantities[i]);
        const std::vector<double> table_data = { mean(equant), variance(equant, stdev) };
        const std::vector<std::string> table_headings = { "Mean", "Std. Dev." };
        const std::vector<int> table_decimals = { repcon.getEnergyDecimalPlaces(),
                                                  repcon.getEnergyDecimalPlaces() };
        result.emplace_back(table_data, table_headings, table_decimals,
                            varname + "_" + lowercase(getEnumerationName(used_quantities[i])),
                            repcon.getReportFileWidth());
      }
      break;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
ReportTable groupSystemDesignations(const SynthesisCacheMap &scmap,
                                    const SystemGrouping organization, std::string *shortcut_key,
                                    const ReportControls &repcon) {
  const SynthesisMapReader scmapr = scmap.data();
  const SystemCache *sc = scmap.getCachePointer();
  std::vector<std::string> table_headings;
  table_headings.push_back("Group Index");
  std::vector<std::string> table_data;
  const std::string common_path_preamble("Common roots paths were identified which can help "
                                         "condense the table that follows:\n");
  switch (organization) {
  case SystemGrouping::SOURCE:
    {
      table_headings.push_back("Coordinate Source");
      table_headings.push_back("Synthesis Count");
      std::vector<std::string> all_sources;
      all_sources.reserve(scmapr.ncache);
      std::vector<int> all_counts(scmapr.ncache);
      table_data.reserve(3 * scmapr.ncache);
      for (int i = 0; i < scmapr.ncache; i++) {
        all_sources.push_back(sc->getInputCoordinatesName(i));
        all_counts[i] = scmapr.csystem_bounds[i + 1] - scmapr.csystem_bounds[i];
        table_data.push_back(std::to_string(i + 1));
      }
      *shortcut_key = common_path_preamble +
                      listCommonPaths(extractCommonPaths(&all_sources, repcon.getCommonPathLimit(),
                                                         repcon.getCommonPathThreshold()));
      table_data.insert(table_data.end(), all_sources.begin(), all_sources.end());
      for (int i = 0; i < scmapr.ncache; i++) {
        table_data.push_back(std::to_string(all_counts[i]));
      }
    }
    break;
  case SystemGrouping::TOPOLOGY:
    {
      table_headings.push_back("Topology File");
      table_headings.push_back("Source Count");
      table_headings.push_back("Synthesis Count");
      std::vector<std::string> all_topologies;
      all_topologies.reserve(scmapr.ntopol);
      std::vector<int> src_counts(scmapr.ntopol);
      std::vector<int> syn_counts(scmapr.ntopol);
      table_data.reserve(4 * scmapr.ntopol);
      for (int i = 0; i < scmapr.ntopol; i++) {
        const int llim = scmapr.ctopol_bounds[i];
        const int hlim = scmapr.ctopol_bounds[i + 1];
        syn_counts[i] = hlim - llim;
        std::vector<bool> uses_topology(scmapr.ncache, false);
        for (int j = llim; j < hlim; j++) {
          uses_topology[scmapr.cache_origins[scmapr.ctopol_proj[j]]] = true;
        }
        int nsrc = 0;
        for (int j = 0; j < scmapr.ncache; j++) {
          nsrc += uses_topology[j];
        }
        src_counts[i] = nsrc;
        all_topologies.push_back(sc->getTopologyPointer(i)->getFileName());
        table_data.push_back(std::to_string(i + 1));
      }
      *shortcut_key = common_path_preamble +
                      listCommonPaths(extractCommonPaths(&all_topologies,
                                                         repcon.getCommonPathLimit(),
                                                         repcon.getCommonPathThreshold()));
      table_data.insert(table_data.end(), all_topologies.begin(), all_topologies.end());
      for (int i = 0; i < scmapr.ntopol; i++) {
        table_data.push_back(std::to_string(src_counts[i]));
      }
      for (int i = 0; i < scmapr.ntopol; i++) {
        table_data.push_back(std::to_string(syn_counts[i]));
      }
    }
    break;
  case SystemGrouping::LABEL:
    {
      table_headings.push_back("Label");
      table_headings.push_back("Source Count");
      table_headings.push_back("Topology Count");
      table_headings.push_back("Synthesis Count");
      std::vector<std::string> all_labels;
      all_labels.reserve(scmapr.nlabel);
      std::vector<int> src_counts(scmapr.nlabel);
      std::vector<int> top_counts(scmapr.nlabel);
      std::vector<int> syn_counts(scmapr.nlabel);
      table_data.reserve(5 * scmapr.nlabel);
      for (int i = 0; i < scmapr.nlabel; i++) {
        const int llim = scmapr.clabel_bounds[i];
        const int hlim = scmapr.clabel_bounds[i + 1];
        syn_counts[i] = hlim - llim;
        std::vector<bool> src_in_label(scmapr.ncache, false);
        std::vector<bool> top_in_label(scmapr.ntopol, false);
        for (int j = llim; j < hlim; j++) {
          src_in_label[scmapr.cache_origins[scmapr.clabel_proj[j]]] = true;
          top_in_label[scmapr.topology_origins[scmapr.clabel_proj[j]]] = true;
        }
        int nsrc = 0;
        int ntop = 0;
        for (int j = llim; j < hlim; j++) {
          nsrc += src_in_label[j];
          ntop += top_in_label[j];
        }
        src_counts[i] = nsrc;
        top_counts[i] = ntop;
        all_labels.push_back(sc->getLabel(i));
        table_data.push_back(std::to_string(i + 1));
      }
      table_data.insert(table_data.end(), all_labels.begin(), all_labels.end());
      for (int i = 0; i < scmapr.nlabel; i++) {
        table_data.push_back(std::to_string(src_counts[i]));
      }
      for (int i = 0; i < scmapr.nlabel; i++) {
        table_data.push_back(std::to_string(top_counts[i]));
      }
      for (int i = 0; i < scmapr.nlabel; i++) {
        table_data.push_back(std::to_string(syn_counts[i]));
      }      
    }
    break;
  }

  // The leftmost column with group indices will be right-justified, but the second column with
  // source, topology, or label names will be left-justified.  All subsequent columns with system
  // counts will be right-justified.
  return ReportTable(table_data, table_headings, std::string(""), repcon.getReportFileWidth(),
                     { JustifyText::RIGHT, JustifyText::LEFT });
}
  
//-------------------------------------------------------------------------------------------------
ReportTable groupEnergyAccumulation(const ScoreCard &nrg, const EnergySample measure,
                                    const int* system_indices, const int* group_bounds,
                                    const int group_count, const StateVariable quantity,
                                    const std::string &varname, const ReportControls &repcon,
                                    const std::vector<bool> &report_mask) {
  size_t total_data = (group_count * 2) + 1;
  size_t npts;
  switch (measure) {
  case EnergySample::TIME_SERIES:
    npts = nrg.getSampleSize();
    total_data *= npts;
    break;
  case EnergySample::FINAL:
  case EnergySample::TIME_AVERAGE:
    npts = 1;
    break;
  }
  std::vector<double> table_data(total_data);
  std::vector<std::string> table_headings;
  table_headings.reserve((group_count * 2) + 1);
  table_headings.push_back("Step Number");
  for (size_t i = 0; i < npts; i++) {
    table_data[i] = nrg.getTimeStep(i);
  }
  size_t adcon = npts;
  for (int grp = 0; grp < group_count; grp++) {
    const int low_limit  = group_bounds[grp];
    const int high_limit = group_bounds[grp + 1];
    std::vector<double2> iresult;
    switch (measure) {
    case EnergySample::TIME_SERIES:
      iresult.resize(nrg.getSampleSize(), { 0.0, 0.0 });
      break;
    case EnergySample::FINAL:
    case EnergySample::TIME_AVERAGE:
      iresult.resize(1, { 0.0, 0.0 });
      break;
    }
    for (int i = low_limit; i < high_limit; i++) {
      switch (measure) {
      case EnergySample::TIME_SERIES:
        {
          const int npts = nrg.getSampleSize();
          const std::vector<double> contrib = nrg.reportHistory(quantity, i);
          for (int j = 0; j < npts; j++) {
            iresult[j].x += contrib[j];
            iresult[j].y += contrib[j] * contrib[j];
          }
        }
        break;
      case EnergySample::FINAL:
        {
          const double contrib = nrg.reportInstantaneousStates(quantity, i);
          iresult[0].x += contrib;
          iresult[0].y += contrib * contrib;
        }
        break;
      case EnergySample::TIME_AVERAGE:
        {
          const double contrib = nrg.reportAverageStates(quantity, i);
          iresult[0].x += contrib;
          iresult[0].y += contrib * contrib;
        }
        break;
      }
    }
    const int nsys = high_limit - low_limit;
    const double dnsys = nsys;
    for (size_t i = 0; i < npts; i++) {
      iresult[i].y = sqrt((dnsys * iresult[i].y) - (iresult[i].x * iresult[i].x)) /
                     sqrt(dnsys * (dnsys - 1.0));
      iresult[i].x /= dnsys;
    }
    for (size_t i = 0; i < npts; i++) {
      table_data[adcon] = iresult[i].x;
      adcon++;
    }
    for (size_t i = 0; i < npts; i++) {
      table_data[adcon] = iresult[i].y;
      adcon++;
    }
    table_headings.push_back("Mean " + std::to_string(grp + 1));
    table_headings.push_back("Std. Dev. " + std::to_string(grp + 1));
  }
  std::vector<int> table_decimals((group_count * 2) + 1, repcon.getEnergyDecimalPlaces());
  table_decimals[0] = 0;
  return ReportTable(table_data, table_headings, table_decimals,
                     varname + "_" + lowercase(getEnumerationName(quantity)),
                     repcon.getReportFileWidth());
}

//-------------------------------------------------------------------------------------------------
std::vector<ReportTable> tabulateGroupedEnergy(const ScoreCard &nrg,
                                               const SystemGrouping organization,
                                               const EnergySample measure,
                                               const std::string &varname,
                                               std::string *shortcut_key,
                                               const SynthesisCacheMap &scmap,
                                               const std::vector<StateVariable> &quantities,
                                               const ReportControls &repcon,
                                               const std::vector<bool> &report_mask) {
  const SynthesisMapReader scmapr = scmap.data();

  // Translate the mask of reported systems from an order corresponding to a straight reading of
  // the systems in the systems cache into one in which systems for each of the relevant groups
  // appear in order.
  std::vector<bool> group_ordered_mask(scmap.getSynthesisSystemCount(), false);
  int ijcon = 0;
  switch (organization) {
  case SystemGrouping::SOURCE:
    for (int i = 0; i < scmapr.ncache; i++) {
      for (int j = scmapr.csystem_bounds[i]; j < scmapr.csystem_bounds[i + 1]; j++) {
        group_ordered_mask[ijcon] = report_mask[scmapr.csystem_proj[j]];
        ijcon++;
      }
    }
    break;
  case SystemGrouping::TOPOLOGY:
    for (int i = 0; i < scmapr.ntopol; i++) {
      for (int j = scmapr.ctopol_bounds[i]; j < scmapr.ctopol_bounds[i + 1]; j++) {
        group_ordered_mask[ijcon] = report_mask[scmapr.ctopol_proj[j]];
        ijcon++;
      }
    }
    break;
  case SystemGrouping::LABEL:
    for (int i = 0; i < scmapr.nlabel; i++) {
      for (int j = scmapr.clabel_bounds[i]; j < scmapr.clabel_bounds[i + 1]; j++) {
        group_ordered_mask[ijcon] = report_mask[scmapr.clabel_proj[j]];
        ijcon++;
      }
    }
    break;
  }
  const std::vector<StateVariable> used_quantities = assessEnergyDetails(nrg, scmap, quantities);
  const int nquant = used_quantities.size();
  std::vector<ReportTable> result;
  result.reserve(1 + nquant);
  result.push_back(groupSystemDesignations(scmap, organization, shortcut_key, repcon));
  for (int i = 0; i < nquant; i++) {
    switch (organization) {
    case SystemGrouping::SOURCE:
      result.push_back(groupEnergyAccumulation(nrg, measure, scmapr.csystem_proj,
                                               scmapr.csystem_bounds, scmapr.ncache,
                                               used_quantities[i], varname, repcon, report_mask));
      break;
    case SystemGrouping::TOPOLOGY:
      result.push_back(groupEnergyAccumulation(nrg, measure, scmapr.ctopol_proj,
                                               scmapr.ctopol_bounds, scmapr.ntopol,
                                               used_quantities[i], varname, repcon, report_mask));
      break;
    case SystemGrouping::LABEL:
      result.push_back(groupEnergyAccumulation(nrg, measure, scmapr.clabel_proj,
                                               scmapr.clabel_bounds, scmapr.nlabel,
                                               used_quantities[i], varname, repcon, report_mask));
      break;
    }
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<ReportTable> tabulateGroupedEnergy(const ScoreCard &nrg,
                                               const SystemGrouping organization,
                                               const EnergySample measure,
                                               const std::string &varname,
                                               std::string *shortcut_key,
                                               const SynthesisCacheMap &scmap,
                                               const std::vector<StateVariable> &quantities,
                                               const ReportControls &repcon,
                                               const std::vector<int> &report_mask) {
  const std::vector<bool> brep_mask = colorVectorMask(report_mask,
                                                      scmap.getSynthesisSystemCount());
  return tabulateGroupedEnergy(nrg, organization, measure, varname, shortcut_key, scmap,
                               quantities, repcon, brep_mask);
}

//-------------------------------------------------------------------------------------------------
std::vector<ReportTable> tabulateGroupedEnergy(const ScoreCard &nrg,
                                               const SystemGrouping organization,
                                               const EnergySample measure,
                                               const std::string &varname,
                                               std::string *shortcut_key,
                                               const SynthesisCacheMap &scmap,
                                               const std::vector<StateVariable> &quantities,
                                               const ReportControls &repcon) {
  const std::vector<bool> report_mask(scmap.getSynthesisSystemCount(), true);
  return tabulateGroupedEnergy(nrg, organization, measure, varname, shortcut_key, scmap,
                               quantities, repcon, report_mask);
}

//-------------------------------------------------------------------------------------------------
ReportTable allSystemDesignations(std::string *shortcut_key, const SynthesisCacheMap &scmap,
                                  const ReportControls &repcon) {
  std::string explain_paths("In the following energetics table(s), each system is given by a "
                            "numerical identifier as follows");
  const SynthesisMapReader scmapr = scmap.data();
  std::vector<std::string> topology_names;
  topology_names.reserve(scmapr.nsynth);
  const SystemCache *sc = scmap.getCachePointer();
  for (int i = 0; i < scmapr.nsynth; i++) {
    const AtomGraph *iag_ptr = sc->getTopologyPointer(scmapr.topology_origins[i]);
    topology_names.push_back(iag_ptr->getFileName());
  }
  const std::vector<std::string> tn_shortcuts = extractCommonPaths(&topology_names, 8, 3);
  if (tn_shortcuts.size() > 0) {
    const std::string pluralize = (tn_shortcuts.size() > 1) ? "path shortcuts" : "a path shortcut";
    explain_paths += ".  This table of keys uses " + pluralize + " to condense its contents:\n";
    explain_paths += listCommonPaths(tn_shortcuts);
  }
  else {
    explain_paths += ":";
  }
  *shortcut_key = explain_paths;
  const std::vector<std::string> all_system_name_headings = {
    "System Number", "&files Entry",  "Topology Origin", "Label"
  };
  std::vector<std::string> all_system_names_data;
  all_system_names_data.reserve(scmapr.nsynth * 4);
  for (int i = 0; i < scmapr.nsynth; i++) {
    all_system_names_data.push_back(std::to_string(i + 1));
  }
  for (int i = 0; i < scmapr.nsynth; i++) {
    all_system_names_data.push_back(std::to_string(scmapr.cache_origins[i] + 1));
  }
  for (int i = 0; i < scmapr.nsynth; i++) {
    all_system_names_data.push_back(topology_names[i]);
  }
  for (int i = 0; i < scmapr.nsynth; i++) {
    all_system_names_data.push_back(sc->getSystemLabel(scmapr.cache_origins[i]));
  }
  return ReportTable(all_system_names_data, all_system_name_headings, std::string(""),
                     repcon.getReportFileWidth(), { JustifyText::RIGHT, JustifyText::RIGHT,
                                                    JustifyText::LEFT, JustifyText::RIGHT });
}
  
//-------------------------------------------------------------------------------------------------
std::vector<ReportTable> tabulateFullEnergy(const ScoreCard &nrg, const EnergySample measure,
                                            const std::string &varname, std::string *shortcut_key,
                                            std::string *post_script,
                                            const SynthesisCacheMap &scmap,
                                            const std::vector<StateVariable> &quantities,
                                            const ReportControls &repcon,
                                            const std::vector<bool> &report_mask) {
  const std::vector<StateVariable> used_quantities = assessEnergyDetails(nrg, scmap, quantities);
  const int nquant = used_quantities.size();
  std::vector<ReportTable> result;
  result.reserve(1 + nquant);
  result.push_back(allSystemDesignations(shortcut_key, scmap, repcon));
  const int nsys = scmap.getSynthesisSystemCount();
  if (report_mask.size() < nsys) {
    rtErr("A mask of " + std::to_string(report_mask.size()) + " is insufficient to cover a "
          "collection of " + std::to_string(nsys) + " systems.", "tabulateFullEnergy");
  }
  int nrep = 0;
  for (int i = 0; i < nsys; i++) {
    nrep += static_cast<int>(report_mask[i]);
  }
  std::vector<double> table_data;
  std::vector<std::string> table_headings;
  const int npts = nrg.getSampleSize();
  post_script->resize(0);
  switch (measure) {
  case EnergySample::TIME_SERIES:
    table_data.resize(static_cast<size_t>(nrep + 1) * static_cast<size_t>(npts));
    table_headings.reserve(nrep + 1);
    table_headings.push_back("Step Number");
    for (int i = 0; i < nsys; i++) {
      if (report_mask[i]) {
        table_headings.push_back("System " + std::to_string(i + 1));
      }
    }
    break;
  case EnergySample::FINAL:
    table_data.resize(nrep * nquant);
    table_headings.reserve(nrep);
    for (int i = 0; i < nsys; i++) {
      if (report_mask[i]) {
        table_headings.push_back("System " + std::to_string(i + 1));
      }
    }
    break;
  case EnergySample::TIME_AVERAGE:
    table_data.resize(2LLU * static_cast<size_t>(nrep * nquant));
    table_headings.reserve(2LLU * static_cast<size_t>(nrep));
    for (int i = 0; i < nsys; i++) {
      if (report_mask[i]) {
        table_headings.push_back("System " + std::to_string(i + 1) + " Mean");
        table_headings.push_back("System " + std::to_string(i + 1) + " Std. Dev.");
      }
    }
    break;
  }
  switch (measure) {
  case EnergySample::TIME_SERIES:
    for (int i = 0; i < nquant; i++) {
      for (int j = 0; j < npts; j++) {
        table_data[j] = nrg.getTimeStep(j);
      }
      size_t adcon = npts;
      for (int j = 0; j < nsys; j++) {
        if (report_mask[j]) {
          const std::vector<double> jhist = nrg.reportHistory(used_quantities[i], j);
          for (int k = 0; k < npts; k++) {
            table_data[adcon] = jhist[k];
            adcon++;
          }
        }
      }
      std::vector<int> table_decimals(nrep + 1, repcon.getEnergyDecimalPlaces());
      table_decimals[0] = 0;
      result.emplace_back(table_data, table_headings, table_decimals,
                          varname + "_" + lowercase(getEnumerationName(used_quantities[i])));
    }
    break;
  case EnergySample::FINAL:
    {
      // Make a temporary table of all energy components, then have the script break it into
      // separate variables.  This is a clean and legible approach that will also be performant in
      // posterior analysis.
      size_t adcon = 0;
      for (int i = 0; i < nsys; i++) {
        if (report_mask[i]) {
          for (int j = 0; j < nquant; j++) {
            table_data[adcon] = nrg.reportInstantaneousStates(used_quantities[j], i);
            adcon++;
          }
        }
      }
      result.emplace_back(table_data, table_headings,
                          std::vector<int>(nrep, repcon.getEnergyDecimalPlaces()),
                          varname + "_aggregate");
    }
    break;
  case EnergySample::TIME_AVERAGE:

    // Make a temporary table of all time-averaged energy components and their error bars, then
    // have the script break it into separate variables.
    for (int i = 0; i < nsys; i++) {
      if (report_mask[i]) {
        for (int j = 0; j < nquant; j++) {
          std::vector<double> ehist = nrg.reportHistory(used_quantities[j], i);
          table_data[  (2 * i * nrep)       + j] = mean(ehist);
          table_data[(((2 * i) + 1) * nrep) + j] = variance(ehist,
                                                            VarianceMethod::STANDARD_DEVIATION);
        }
      }
    }
    result.emplace_back(table_data, table_headings,
                        std::vector<int>(nrep, repcon.getEnergyDecimalPlaces()),
                        varname + "_aggregate");
    break;
  }
  switch (measure) {
  case EnergySample::TIME_SERIES:
    break;
  case EnergySample::FINAL:
  case EnergySample::TIME_AVERAGE:
    for (int i = 0; i < nquant; i++) {
      switch (repcon.getOutputSyntax()) {
      case OutputSyntax::MATPLOTLIB:
        post_script->append(varname + "_" + lowercase(getEnumerationName(used_quantities[i])) +
                            " = " + varname + "_aggregate[" + std::to_string(i) + "]");
        break;
      case OutputSyntax::MATRIX_PKG:
        post_script->append(varname + "_" + lowercase(getEnumerationName(used_quantities[i])) +
                            " = " + varname + "_aggregate(" + std::to_string(i + 1) + ", :);\n");
        break;
      case OutputSyntax::STANDALONE:
        post_script->append(std::string("| Row ") + ((i < 9) ? " " : "") + std::to_string(i + 1) +
                            " : " + lowercase(getEnumerationName(used_quantities[i])));
        break;
      }
      switch (repcon.getOutputSyntax()) {
      case OutputSyntax::MATPLOTLIB:
      case OutputSyntax::STANDALONE:
        if (i < nquant - 1) {
          post_script->append("\n");
        }
        break;
      case OutputSyntax::MATRIX_PKG:
        break;
      }
    }
    switch (repcon.getOutputSyntax()) {
    case OutputSyntax::MATPLOTLIB:
      break;
    case OutputSyntax::MATRIX_PKG:
      post_script->append(std::string("clear ") + varname + "_aggregate;");
      break;
    case OutputSyntax::STANDALONE:
      break;
    }
    break;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<ReportTable> tabulateFullEnergy(const ScoreCard &nrg, const EnergySample measure,
                                            const std::string &varname, std::string *shortcut_key,
                                            std::string *post_script,
                                            const SynthesisCacheMap &scmap,
                                            const std::vector<StateVariable> &quantities,
                                            const ReportControls &repcon,
                                            const std::vector<int> &report_mask) {
  const std::vector<bool> brep_mask = colorVectorMask(report_mask,
                                                      scmap.getSynthesisSystemCount());
  return tabulateFullEnergy(nrg, measure, varname, shortcut_key, post_script, scmap, quantities,
                            repcon, brep_mask);
}

//-------------------------------------------------------------------------------------------------
std::vector<ReportTable> tabulateFullEnergy(const ScoreCard &nrg, const EnergySample measure,
                                            const std::string &varname, std::string *shortcut_key,
                                            std::string *post_script,
                                            const SynthesisCacheMap &scmap,
                                            const std::vector<StateVariable> &quantities,
                                            const ReportControls &repcon) {
  const std::vector<bool> report_mask(scmap.getSynthesisSystemCount(), true);
  return tabulateFullEnergy(nrg, measure, varname, shortcut_key, post_script, scmap, quantities,
                            repcon, report_mask);
}

//-------------------------------------------------------------------------------------------------
std::vector<ReportTable> tabulateOutlierEnergy(const ScoreCard &nrg, const EnergySample measure,
                                               std::string *shortcut_key, std::string *post_script,
                                               const SynthesisCacheMap &scmap,
                                               const SystemGrouping organization,
                                               const std::vector<StateVariable> &quantities,
                                               const ReportControls &repcon,
                                               const std::vector<bool> &report_mask) {
  std::vector<ReportTable> result;
  
  // Determine the outliers based on the final energies, whether across the entire synthesis or
  // in terms of specific groups of systems.
  const SynthesisMapReader scmapr = scmap.data();
  const std::vector<double> final_e = nrg.reportInstantaneousStates();
  std::vector<double> outlier_scores;
  std::vector<int> outlier_indices;
  switch (repcon.getOutputScope()) {
  case OutputScope::AVERAGES:
  case OutputScope::CLUSTER_AVERAGES:
  case OutputScope::FULL:

    // These cases will have redirected to tabulateFullEnergy() or tabulateAverageEnerg(), above.
    // The output sccope is evaluated a second time internally to differentiate between cluster
    // outliers and cache-wide outliers.
    break;
  case OutputScope::OUTLIERS:
    {
      // Analyze the final energies of each system and compare them based on the mean final energy
      // for all structures sharing the same topology.  Report outliers based on deviations from
      // their respective mean energies.
      int maxsys = 0;
      for (int i = 0; i < scmapr.ntopol; i++) {
        maxsys = std::max(maxsys, scmapr.ctopol_bounds[i + 1] - scmapr.ctopol_bounds[i]);
      }
      std::vector<double> common_e(maxsys);
      outlier_scores.resize(repcon.getOutlierCount());
      outlier_indices.resize(repcon.getOutlierCount());
      double min_outlier_score = 0.0;
      int noutlier = 0;
      for (int i = 0; i < scmapr.ntopol; i++) {
        const int llim = scmapr.ctopol_bounds[i];
        const int hlim = scmapr.ctopol_bounds[i + 1];

        // Pass to the next topology if there are no applicable systems.
        if (hlim == llim) {
          continue;
        }
        for (int j = llim; j < hlim; j++) {
          common_e[j - llim] = final_e[scmapr.ctopol_proj[j]];
        }
        const int nt = hlim - llim;
        const double mean_e = mean(common_e.data(), nt);
        const double stdev_e = variance(common_e.data(), nt, VarianceMethod::STANDARD_DEVIATION);
        const double othresh = repcon.getOutlierSigmaFactor() * stdev_e;
        for (int j = 0; j < nt; j++) {
          const double oscore = fabs(common_e[j] - mean_e);
          if (oscore > othresh && oscore > min_outlier_score) {
            int ocon = noutlier - 1;
            while (ocon >= 0 && oscore > outlier_scores[ocon]) {
              ocon--;
            }
            ocon++;
            for (int i = std::min(noutlier, repcon.getOutlierCount()) - 1; i > ocon; i--) {
              outlier_scores[i] = outlier_scores[i - 1];
              outlier_indices[i] = outlier_indices[i - 1];
            }
            outlier_scores[ocon]  = oscore;
            outlier_indices[ocon] = scmapr.ctopol_proj[scmapr.ctopol_bounds[i] + j];
            noutlier += (noutlier < repcon.getOutlierCount());
            min_outlier_score = outlier_scores[noutlier - 1];
          }
        }
      }
    }
    break;
  case OutputScope::CLUSTER_OUTLIERS:
    {
      int maxsys = 0;
      int ngrp;
      switch (organization) {
      case SystemGrouping::SOURCE:
        for (int i = 0; i < scmapr.ncache; i++) {
          maxsys = std::max(maxsys, scmapr.csystem_bounds[i + 1] - scmapr.csystem_bounds[i]);
        }
        ngrp = scmapr.ncache;
        break;
      case SystemGrouping::TOPOLOGY:
        for (int i = 0; i < scmapr.ntopol; i++) {
          maxsys = std::max(maxsys, scmapr.ctopol_bounds[i + 1] - scmapr.ctopol_bounds[i]);
        }
        ngrp = scmapr.ntopol;
        break;
      case SystemGrouping::LABEL:
        for (int i = 0; i < scmapr.nlabel; i++) {
          maxsys = std::max(maxsys, scmapr.clabel_bounds[i + 1] - scmapr.clabel_bounds[i]);
        }
        ngrp = scmapr.nlabel;
        break;
      }
      std::vector<double> common_e(maxsys);
      std::vector<int> common_idx(maxsys);
      for (int i = 0; i < ngrp; i++) {
        int llim, hlim;
        switch (organization) {
        case SystemGrouping::SOURCE:
          llim = scmapr.csystem_bounds[i];
          hlim = scmapr.csystem_bounds[i + 1];
          for (int j = llim; j < hlim; j++) {
            common_e[j - llim] = final_e[scmapr.csystem_proj[j]];
            common_idx[j - llim] = scmapr.csystem_proj[j];
          }
          break;
        case SystemGrouping::TOPOLOGY:
          llim = scmapr.ctopol_bounds[i];
          hlim = scmapr.ctopol_bounds[i + 1];
          for (int j = llim; j < hlim; j++) {
            common_e[j - llim] = final_e[scmapr.ctopol_proj[j]];
            common_idx[j - llim] = scmapr.ctopol_proj[j];
          }
          break;
        case SystemGrouping::LABEL:
          llim = scmapr.clabel_bounds[i];
          hlim = scmapr.clabel_bounds[i + 1];
          for (int j = llim; j < hlim; j++) {
            common_e[j - llim] = final_e[scmapr.clabel_proj[j]];
            common_idx[j - llim] = scmapr.clabel_proj[j];
          }
          break;
        }
        const int nt = hlim - llim;
        const double mean_e = mean(common_e.data(), nt);
        const double stdev_e = variance(common_e.data(), nt, VarianceMethod::STANDARD_DEVIATION);
        const double othresh = repcon.getOutlierSigmaFactor() * stdev_e;
        double min_outlier_score = othresh;
        std::vector<double> sysi_scores(repcon.getOutlierCount());
        std::vector<double> sysi_indices(repcon.getOutlierCount());
        int noutlier = 0;
        for (int j = 0; j < nt; j++) {
          const double oscore = fabs(common_e[j] - mean_e);
          if (oscore > min_outlier_score) {
            int ocon = noutlier - 1;
            while (ocon >= 0 && oscore > sysi_scores[ocon]) {
              ocon--;
            }
            ocon++;
            for (int i = std::min(noutlier, repcon.getOutlierCount()) - 1; i > ocon; i--) {
              sysi_scores[i] = sysi_scores[i - 1];
              sysi_indices[i] = sysi_indices[i - 1];
            }
            sysi_scores[ocon]  = oscore;
            sysi_indices[ocon] = common_idx[j];
            noutlier += (noutlier < repcon.getOutlierCount());
            min_outlier_score = sysi_scores[noutlier - 1];
          }
        }
        for (int j = 0; j < noutlier; j++) {
          outlier_scores.push_back(sysi_scores[j]);
          outlier_indices.push_back(sysi_indices[j]);
        }
      }
    }
    break;
  }

  return result;
}

//-------------------------------------------------------------------------------------------------
std::vector<ReportTable> tabulateOutlierEnergy(const ScoreCard &nrg, const EnergySample measure,
                                               std::string *shortcut_key, std::string *post_script,
                                               const SynthesisCacheMap &scmap,
                                               const SystemGrouping organization,
                                               const std::vector<StateVariable> &quantities,
                                               const ReportControls &repcon,
                                               const std::vector<int> &report_mask) {
  const std::vector<bool> brep_mask = colorVectorMask(report_mask,
                                                      scmap.getSynthesisSystemCount());
  return tabulateOutlierEnergy(nrg, measure, shortcut_key, post_script, scmap, organization,
                               quantities, repcon, brep_mask);
}

//-------------------------------------------------------------------------------------------------
std::vector<ReportTable> tabulateOutlierEnergy(const ScoreCard &nrg, const EnergySample measure,
                                               std::string *shortcut_key, std::string *post_script,
                                               const SynthesisCacheMap &scmap,
                                               const SystemGrouping organization,
                                               const std::vector<StateVariable> &quantities,
                                               const ReportControls &repcon) {
  const std::vector<bool> report_mask(scmap.getSynthesisSystemCount(), true);
  return tabulateOutlierEnergy(nrg, measure, shortcut_key, post_script, scmap, organization,
                               quantities, repcon, report_mask);
}

//-------------------------------------------------------------------------------------------------
void createDiagnosticReport(const ScoreCard &nrg, const SynthesisCacheMap &scmap,
                            const ReportControls &repcon, std::ofstream *foutp) {
  
  const std::string end_line("\n");
  foutp->write(end_line.data(), end_line.size());
}

//-------------------------------------------------------------------------------------------------
void createDiagnosticReport(const ScoreCard &nrg, const SynthesisCacheMap &scmap,
                            const UserSettings &ui) {
  const ReportControls repcon = ui.getReportNamelistInfo();
  const FilesControls ficon = ui.getFilesNamelistInfo();
  const std::vector<StateVariable> quantities = repcon.getReportedQuantities();
  std::vector<ReportTable> diag_contents;
  std::string shortcut_key(""), post_script("");
  switch (repcon.getOutputScope()) {
  case OutputScope::AVERAGES:
  case OutputScope::CLUSTER_AVERAGES:
  case OutputScope::FULL:
    diag_contents = tabulateFullEnergy(nrg, repcon.getEnergySamplingMethod(),
                                       repcon.getReportVariable(), &shortcut_key, &post_script,
                                       scmap, repcon.getReportedQuantities(), repcon);
    break;
  case OutputScope::OUTLIERS:
  case OutputScope::CLUSTER_OUTLIERS:
    //tabulateOutlierEnergy();
    break;
  }

  // Various sections of the input will be ordered within a Standard Template Library vector.
  std::vector<SectionContents> all_sect;
  
  // Create a summary of the user input
  if (ui.getInputFileName().size() > 0) {
    all_sect.emplace_back("User Input", ficon.getReportFile(), repcon.getReportFileWidth(),
                          repcon.getOutputSyntax(), 1.0e-4);
    all_sect.back().addNarration("This is the input file:");
    const TextFile input_tf(ui.getInputFileName());
    all_sect.back().addNarration(input_tf);
  }
  all_sect.emplace_back("Energy and State Variables", ficon.getReportFile(),
                        repcon.getReportFileWidth(), repcon.getOutputSyntax(),
                        pow(10.0, -repcon.getEnergyDecimalPlaces()));
  all_sect.back().addNarration(shortcut_key);
  for (size_t i = 0; i < diag_contents.size(); i++) {
    all_sect.back().addTable(diag_contents[i]);
  }
  printAllSections(ficon.getReportFile(), ui.getPrintingPolicy(), all_sect,
                   repcon.getOutputSyntax());
}
  
} // namespace review
} // namespace stormm
