#include <string>
#include <vector>
#include "copyright.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/series_ops.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Namelists/nml_conformer.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Synthesis/atomgraph_synthesis.h"
#include "../../src/Synthesis/phasespace_synthesis.h"
#include "../../src/Synthesis/synthesis_cache_map.h"
#include "../../src/Synthesis/synthesis_enumerators.h"
#include "../../src/Synthesis/synthesis_permutor.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::diskutil;
using namespace stormm::namelist;
using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::stmath;
using namespace stormm::symbols;
using namespace stormm::synthesis;
using namespace stormm::testing;

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Lay out the test materials
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }
  StopWatch timer;
  const char osc = osSeparator();
  const std::string base_top_name = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_name = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  const std::vector<std::string> molmix = { "symmetry_L1", "med_1", "med_4", "stereo_L1",
                                            "med_3", "symmetry_L1_vs", "bromobenzene_vs_iso",
                                            "med_2", "bromobenzene_iso", "symmetry_C2",
                                            "symmetry_C3", "drug_example_iso", "symmetry_C6" };
  TestSystemManager tsm(base_top_name, "top", molmix, base_crd_name, "inpcrd", molmix);

  // Create a synthesis with permutor
  const std::vector<int> inc_indices = incrementingSeries<int>(0, molmix.size());
  AtomGraphSynthesis poly_ag = tsm.exportAtomGraphSynthesis(inc_indices);
  PhaseSpaceSynthesis poly_ps = tsm.exportPhaseSpaceSynthesis(inc_indices);
  SystemCache sc = tsm.exportSystemCache(inc_indices, oe);
  const SynthesisCacheMap scmap(inc_indices, sc, poly_ag, poly_ps);
  SynthesisPermutor syper(sc.getFeaturesPointer(), true, &timer);
  const std::vector<double> rot_values = { pi / 3.0, pi, -pi / 3.0 };
  const std::vector<double> ctx_values = { 0.0, pi };
  syper.setVariableRanges(rot_values, ctx_values);
  syper.applySynthesis(poly_ps, VariableTorsionAdjustment::DO_NOT_CHANGE);
  Xoshiro256ppGenerator xrs;
  ClashReport clrep;
#ifndef STORMM_USE_HPC
  PhaseSpaceSynthesis minimal_synth = syper.buildSynthesis(SamplingIntensity::MINIMAL, &xrs,
                                                           10000, 28, 26, 40, 28,
                                                           PrecisionModel::DOUBLE);
  PhaseSpaceSynthesis light_synth = syper.buildSynthesis(SamplingIntensity::LIGHT, &xrs, 10000,
                                                         28, 26, 40, 28, PrecisionModel::DOUBLE);
  PhaseSpaceSynthesis heavy_synth = syper.buildSynthesis(SamplingIntensity::HEAVY, &xrs, 10000,
                                                         28, 26, 40, 28, PrecisionModel::DOUBLE);
  PhaseSpaceSynthesis exhaustive_synth = syper.buildSynthesis(SamplingIntensity::EXHAUSTIVE, &xrs,
                                                              100000, 28, 26, 40, 28,
                                                              PrecisionModel::DOUBLE, &clrep);
#endif
  
  //syper.applySynthesis(poly_ps, VariableTorsionAdjustment::ADJUST_NEARBY_VALUES);

  // Summary evaluation
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
