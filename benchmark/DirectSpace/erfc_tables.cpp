#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>
#include <fstream>
#include <iostream>
#include "copyright.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/log_scale_spline.h"
#include "../../src/Math/matrix_ops.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Namelists/command_line_parser.h"
#include "../../src/Parsing/parse.h"
#include "../../src/Potential/pme_util.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Reporting/report_table.h"
#include "../../src/Reporting/reporting_enumerators.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Topology/atomgraph.h"
#include "../../src/Trajectory/coordinateframe.h"
#include "../../src/Trajectory/phasespace.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/test_system_manager.h"

#ifndef STORMM_USE_HPC
using stormm::data_types::double4;
using stormm::data_types::float4;
#endif
using namespace stormm::constants;
using namespace stormm::diskutil;
using namespace stormm::energy;
using namespace stormm::errors;
using namespace stormm::namelist;
using namespace stormm::parse;
using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::stmath;
using namespace stormm::testing;
using namespace stormm::topology;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Place a pair of points and separate them by a set distance.  Return the result as a six-element
// vector orders point 1 {x, y, z}, point 2 {x, y, z}.
//
// Arguments:
//   xrs:       The random number generator used to create the initial locations of the points.
//              Modified as it churns through random numbers.
//   r_target:  The target distance between the two points
//-------------------------------------------------------------------------------------------------
std::vector<double> placePointPair(Xoshiro256ppGenerator *xrs, const double r_target) {
  std::vector<double> result = uniformRand(xrs, 6, 10.0);
  addScalarToVector(&result, -5.0);
  double dx = result[3] - result[0];
  double dy = result[4] - result[1];
  double dz = result[5] - result[2];
  double r2 = (dx * dx) + (dy * dy) + (dz * dz);
  bool refresh;
  do {
    refresh = (r2 < 1.0e-4);
    if (refresh) {
      for (int i = 0; i < 6; i++) {
        result[i] = uniformRand(xrs, 10.0) - 5.0;
        dx = result[3] - result[0];
        dy = result[4] - result[1];
        dz = result[5] - result[2];
        r2 = (dx * dx) + (dy * dy) + (dz * dz);
      }
    }
  } while (refresh);

  // Determine the current separation and expand or contract the distance between the points as
  // needed.
  const double midx = 0.5 * (result[0] + result[3]);
  const double midy = 0.5 * (result[1] + result[4]);
  const double midz = 0.5 * (result[2] + result[5]);
  const double r_init = sqrt(r2);
  const double rscale = 0.5 * r_target / r_init;
  result[0] = midx - (dx * rscale);
  result[1] = midy - (dy * rscale);
  result[2] = midz - (dz * rscale);
  result[3] = midx + (dx * rscale);
  result[4] = midy + (dy * rscale);
  result[5] = midz + (dz * rscale);
  return result;
}

//-------------------------------------------------------------------------------------------------
// Open the output file and add some comments to indicate to any reader the content of the
// variables that follow.
//
// Arguments:
//   output_file_name:  Name of the output file to open
//   lgsp:              The logarithmic spline table which will become the source of results,
//                      containing critical parameters to be conveyed as annotation in the output
//-------------------------------------------------------------------------------------------------
template <typename T4>
std::ofstream prepOutputFile(const std::string &output_file_name, const LogScaleSpline<T4> &lgsp) {

  // Prepare to write the summary
  std::ofstream foutp = openOutputFile(output_file_name, PrintSituation::APPEND, "open the output "
                                       "for additional test recording");
  const char comment_guard = commentSymbol(OutputSyntax::MATRIX_PKG);
  printProtectedText("Potential form:  " + getEnumerationName(lgsp.getForm()) +
                     "\nIndexing method: " + getEnumerationName(lgsp.getIndexingMethod()) +
                     "\nBasis functions: " + getEnumerationName(lgsp.getBasisSet()) + "\n\n",
                     comment_guard, &foutp, 120);
  return foutp;
}

//-------------------------------------------------------------------------------------------------
// Set a display variable name based on a spline parameters and any custom extensions.  This is
// abstracted from the function below for accessibility throughout the benchmarking program.
//
// Arguments:
//   target_form_in:      The form of the function being approximated
//   indexing_method_in:  The means of deriving a spline table index out of the function argument
//   basis_set_in:        The basis set for spline elements
//   custom_extension:    Additional details to add to the variable name
//-------------------------------------------------------------------------------------------------
std::string setDisplayVariable(const LogSplineForm target_form_in,
                               const TableIndexing indexing_method_in,
                               const BasisFunctions basis_set_in,
                               const std::string &custom_extension = std::string("")) {
  std::string result;
  switch (target_form_in) {
  case LogSplineForm::ELEC_PME_DIRECT:
    result = "u_";
    break;
  case LogSplineForm::ELEC_PME_DIRECT_EXCL:
    result = "ux_";
    break;
  case LogSplineForm::DELEC_PME_DIRECT:
    result = "du_";
    break;
  case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    result = "dux_";
    break;
  case LogSplineForm::CUSTOM:
    break;
  }
  switch (indexing_method_in) {
  case TableIndexing::ARG:
    result += "r_";
    break;
  case TableIndexing::SQUARED_ARG:
    result += "r2_";
    break;
  case TableIndexing::ARG_OFFSET:
    result += "ro_";
    break;
  case TableIndexing::SQ_ARG_OFFSET:
    result += "r2o_";
    break;
  }
  switch (basis_set_in) {
  case BasisFunctions::MIXED_FRACTIONS:
    result += "frac";
    break;
  case BasisFunctions::POLYNOMIAL:
    result += "poly";
    break;
  }
  result += custom_extension;
  return result;
}

//-------------------------------------------------------------------------------------------------
// Set variable names and display text based on critical spline parameters.
//
// Arguments:
//   lgsp:       The logarithmic spline which will be plotted
//   var_name:   Name of the variable to store results in the display script (filled and
//               returned)
//   unit_str:   The units of relevant spline values and comparisons (filled and returned)
//   title_str:  The title of any plot to display (filled and returned)
//-------------------------------------------------------------------------------------------------
template <typename T4>
void setDisplayStrings(const LogScaleSpline<T4> &lgsp, std::string *var_name,
                       std::string *unit_str, std::string *title_str) {
  *var_name = setDisplayVariable(lgsp.getForm(), lgsp.getIndexingMethod(), lgsp.getBasisSet());
  switch (lgsp.getForm()) {
  case LogSplineForm::ELEC_PME_DIRECT:
    *unit_str = "kcal/mol";
    *title_str = "PME Energy";
    break;
  case LogSplineForm::ELEC_PME_DIRECT_EXCL:
    *unit_str = "kcal/mol";
    *title_str = "PME Energy (with Exclusion)";
    break;
  case LogSplineForm::DELEC_PME_DIRECT:
    *unit_str = "kcal/mol-A";
    *title_str = "PME Force";
    break;
  case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    *unit_str = "kcal/mol-A";
    *title_str = "PME Force (with Exclusion)";
    break;
  case LogSplineForm::CUSTOM:
    break;
  }
  switch (lgsp.getIndexingMethod()) {
  case TableIndexing::ARG:
    title_str->append(", Arg");
    break;
  case TableIndexing::SQUARED_ARG:
    title_str->append(", Sq. Arg");
    break;
  case TableIndexing::ARG_OFFSET:
    title_str->append(", Arg Offset");
    break;
  case TableIndexing::SQ_ARG_OFFSET:
    title_str->append(", Sq. Arg Offset");
    break;
  }
  switch (lgsp.getBasisSet()) {
  case BasisFunctions::MIXED_FRACTIONS:
    title_str->append(", Fraction Series");
    break;
  case BasisFunctions::POLYNOMIAL:
    title_str->append(", Polynomial");
    break;
  }
}

//-------------------------------------------------------------------------------------------------
// Evaluate a spline table at many points over its relevant range.
//
// Arguments:
//-------------------------------------------------------------------------------------------------
template <typename T4>
void computeSplineErrors(const LogScaleSpline<T4> &lgsp, const double min_analysis_range,
                         const double max_analysis_range, const bool test_natural,
                         std::vector<double> *rpts, std::vector<double> *flt_eval_mue,
                         std::vector<double> *spl_eval_mue, Xoshiro256ppGenerator *xrs) {
  const double dscr = 0.00390625;
  const int npts = (max_analysis_range - min_analysis_range) / dscr;
  std::vector<double> dbl_eval(npts);
  std::vector<double> flt_eval_means(npts), flt_eval_stdev(npts);
  std::vector<double> spl_eval_means(npts), spl_eval_stdev(npts);
  flt_eval_mue->resize(npts);
  spl_eval_mue->resize(npts);
  rpts->resize(npts);
  double* flt_eval_mue_ptr = flt_eval_mue->data();
  double* spl_eval_mue_ptr = spl_eval_mue->data();
  for (int i = 0; i < npts; i++) {
    flt_eval_mue_ptr[i] = 0.0;
    spl_eval_mue_ptr[i] = 0.0;
  }
  double* rpts_ptr = rpts->data();
  const int ntrials = (test_natural) ? 32 : 1;
  std::vector<double> flt_eval(ntrials), spl_eval(ntrials);
  const double bfac = 2.0 * lgsp.getEwaldCoefficient() / sqrt(stormm::symbols::pi);
  const double kcoul = lgsp.getCoulombConstant();
  const float kcoulf = kcoul;
  const float bfacf = bfac;
  const double ew_coeff = lgsp.getEwaldCoefficient();
  const float ew_coeff_f = ew_coeff;
  for (int i = 0; i < npts; i++) {
    const double r = min_analysis_range + (static_cast<double>(i) * dscr);
    rpts_ptr[i] = r;
    const double ewr = ew_coeff * r;
    switch (lgsp.getForm()) {
    case LogSplineForm::ELEC_PME_DIRECT:
      dbl_eval[i] = kcoul * erfc(ew_coeff * r) / r;
      break;
    case LogSplineForm::ELEC_PME_DIRECT_EXCL:
      dbl_eval[i] = kcoul * (erfc(ew_coeff * r) - 1.0) / r;
      break;
    case LogSplineForm::DELEC_PME_DIRECT:
      dbl_eval[i] = -kcoul * ((bfac * exp(-ewr * ewr)) + (erfc(ew_coeff * r) / r)) / (r * r);
      break;
    case LogSplineForm::DELEC_PME_DIRECT_EXCL:
      dbl_eval[i] = -kcoul * ((bfac * exp(-ewr * ewr)) + ((erfc(ew_coeff * r) - 1.0) / r)) /
                    (r * r);
      break;
    case LogSplineForm::CUSTOM:
      break;
    }

    // For any given value of r, there are many orientations of two particles that could create
    // the proper distance, and do so with different combinations of displacments along the
    // Cartesian axes.  Find a series of these points.
    for (int j = 0; j < ntrials; j++) {

      // Create a pair of points in a random orientation with a distance r between them.  Compute
      // displacements in float in order evaluate the spline as well as the single-precision
      // analytic result.
      float sp_r2;
      if (test_natural) {
        const std::vector<double> coords = placePointPair(xrs, r);
        const float sp_dx = static_cast<float>(coords[3]) - static_cast<float>(coords[0]);
        const float sp_dy = static_cast<float>(coords[4]) - static_cast<float>(coords[1]);
        const float sp_dz = static_cast<float>(coords[5]) - static_cast<float>(coords[2]);    
        sp_r2 = (sp_dx * sp_dx) + (sp_dy * sp_dy) + (sp_dz * sp_dz);
      }
      else {
        sp_r2 = r * r;
      }
      const float sp_r  = sqrtf(sp_r2);
      const float sp_ewr = ew_coeff_f * sp_r;
      switch (lgsp.getForm()) {
      case LogSplineForm::ELEC_PME_DIRECT:
        flt_eval[j] = kcoulf * erfcf(sp_ewr) / sp_r;
        break;
      case LogSplineForm::ELEC_PME_DIRECT_EXCL:
        flt_eval[j] = kcoulf * (erfcf(sp_ewr) - 1.0f) / sp_r;
        break;
      case LogSplineForm::DELEC_PME_DIRECT:
        flt_eval[j] = -kcoulf * ((bfacf * expf(-sp_ewr * sp_ewr)) + (erfcf(sp_ewr) / sp_r)) /
                       sp_r2;
        break;
      case LogSplineForm::DELEC_PME_DIRECT_EXCL:
        flt_eval[j] = -kcoulf * ((bfacf * expf(-sp_ewr * sp_ewr)) +
                                 ((erfcf(sp_ewr) - 1.0f) / sp_r)) / sp_r2;
        break;
      case LogSplineForm::CUSTOM:
        break;
      }
      flt_eval_mue_ptr[i] += fabs(flt_eval[j] - dbl_eval[i]);

      // Evaluate the spline based on the way it is indexed.
      spl_eval[j] = lgsp.evaluateByRealArg(sp_r);
      spl_eval_mue_ptr[i] += fabs(spl_eval[j] - dbl_eval[i]);
    }
    if (test_natural) {
      flt_eval_means[i] = mean(flt_eval);
      flt_eval_stdev[i] = variance(flt_eval, VarianceMethod::STANDARD_DEVIATION);
      spl_eval_means[i] = mean(spl_eval);
      spl_eval_stdev[i] = variance(spl_eval, VarianceMethod::STANDARD_DEVIATION);
      flt_eval_mue_ptr[i] /= static_cast<double>(ntrials);
      spl_eval_mue_ptr[i] /= static_cast<double>(ntrials);
    }
    else {
      flt_eval_means[i] = flt_eval[0];
      flt_eval_stdev[i] = 0.0;
      spl_eval_means[i] = spl_eval[0];
      spl_eval_stdev[i] = 0.0;
    }
  }
}

//-------------------------------------------------------------------------------------------------
// Test a logarithmic spline table in terms of accuracy over a particular range.  Print the result
// for display and processing in popular matrix algebra packages.
//
// Arguments:
//   target_form_in:      Form of the function to spline
//   ew_coeff:            Ewald coefficient defining the splined function
//   min_analysis_range:  Minimum relevant range of the spline, for statistical error computation
//   max_analysis_range:  Maximum relevant range of the spline, for statistical error computation
//   mantissa_bits_in:    Number of bits after the exponent which become part of the table index
//   indexing_method_in:  The indexing method for pulling out spline coefficients
//   basis_set_in:        Basis functions for the spline table that will be constructed
//   indexing_offset_in:  The real-valued offset for the spline argument
//   test_environment:    Not a TestEnvironment object, but rather an indication of whether the
//                        function argument is to be computed in float_32t or float_64t precision
//   print_figure:        Flag to have the figure printed.  Set to FALSE to suppress display of the
//                        results.
//   xrs:                 Source of random numbers for range sampling
//   output_file_name:    Name of the output file to write results into
//   timer:               Tracks the time required to construct various splines
//-------------------------------------------------------------------------------------------------
template <typename T4>
void analyzeTable(const LogSplineForm target_form_in, const double ew_coeff,
                  const double min_analysis_range, const double max_analysis_range,
                  const int mantissa_bits_in, const TableIndexing indexing_method_in,
                  const BasisFunctions basis_set_in, const float indexing_offset_in,
                  const std::string &test_environment, const bool print_figure,
                  Xoshiro256ppGenerator *xrs, const std::string &output_file_name,
                  StopWatch *timer) {
  
  // Compute the spline tables.
  const int table_construction = timer->addCategory("Log Spline Table Building");
  timer->assignTime(0);
  LogScaleSpline<T4> lgsp(target_form_in, ew_coeff, stormm::symbols::charmm_gromacs_bioq,
                          mantissa_bits_in, 4096.0, 0.015625, indexing_method_in, basis_set_in,
                          0, indexing_offset_in);
  timer->assignTime(table_construction);

  // Begin building the report file.  The variable name will be useful for reporting timings of
  // each logarithmic spline table construction.
  std::string var_name, unit_str, title_str;
  setDisplayStrings(lgsp, &var_name, &unit_str, &title_str);
  std::ofstream foutp = prepOutputFile(output_file_name, lgsp);
  std::vector<double> rpts, flt_eval_mue, spl_eval_mue;
  const bool test_natural = strcmpCased(test_environment, "natural", CaseSensitivity::NO);
  computeSplineErrors(lgsp, min_analysis_range, max_analysis_range, test_natural, &rpts,
                      &flt_eval_mue, &spl_eval_mue, xrs);

  // Smooth the data and contribute to the output.
  const int ave_window = 8;
  const int npts = flt_eval_mue.size();
  const int plot_pts = npts / ave_window;
  const double ave_wt = 1.0 / static_cast<double>(ave_window);
  std::vector<double> test_data(plot_pts * 3);
  for (int i = 0; i < plot_pts; i++) {
    double r_ave = 0.0;
    double spl_mue = 0.0;
    double flt_mue = 0.0;
    for (int j = 0; j < ave_window; j++) {
      r_ave += rpts[(i * ave_window) + j];
      spl_mue += spl_eval_mue[i];
      flt_mue += flt_eval_mue[i];
    }
    r_ave *= ave_wt;
    spl_mue = log10(spl_mue * ave_wt);
    flt_mue = log10(flt_mue * ave_wt);
    test_data[i] = r_ave;
    test_data[i +       plot_pts] = spl_mue;
    test_data[i + (2 * plot_pts)] = flt_mue;
  }

  // One final edit to the variable name, irrelevant to the table of timings
  if (test_natural) {
    var_name += "_natural";
    title_str += ", Natural Process";
  }
  ReportTable test_tab(test_data, { "Distance, A", "log(10) Mean Unsigned Error, Spline",
                                    "log(10) Mean Unsigned Error, FP32 Analytic" },
                       std::vector<int>(3, 12), var_name, 120);
  test_tab.printTable(&foutp, OutputSyntax::MATRIX_PKG);
  if (print_figure) {
    std::string plot_command("figure;\nhold on;\n");
    plot_command += "plot(" + var_name + "(:,1), " + var_name + "(:,3), " +
                    encodePsivantColor(PsivantColor::BLUE, OutputSyntax::MATRIX_PKG) +
                    ", 'linewidth', 8);\n";
    plot_command += "plot(" + var_name + "(:,1), " + var_name +
                    "(:,2), 'color', [ 0.1 0.1 0.1 ], 'linewidth', 8);\n";
    plot_command += "axis([ 0 12 -8 -3 ]);\n";
    plot_command += "daspect([ 2.4 1 1 ]);\n";
    plot_command += "xlabel('Distance, A');\nylabel('log10 Mean Unsigned Error, " + unit_str +
                    "');\n";
    plot_command += "title('" + title_str + "');\n";
    plot_command += "legend('Analytic FP32', 'Spline');\n";
    plot_command += "set(gca, 'fontsize', 36);\n";
    foutp.write(plot_command.data(), plot_command.size());
  }
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
// Compute the effect of optimizations with different units of least place.  Descriptions of input
// arguments follow from analyzeTable(), above, in addition to:
//
// Arguments:
//   ulp_opt:   The number of units of least place to let each spline coefficient wander
//-------------------------------------------------------------------------------------------------
template <typename T4>
void testUlpOptimization(const LogSplineForm target_form_in, const double ew_coeff,
                         const double min_analysis_range, const double max_analysis_range,
                         const int mantissa_bits_in, const TableIndexing indexing_method_in,
                         const BasisFunctions basis_set_in, const float indexing_offset_in,
                         const int max_ulp, Xoshiro256ppGenerator *xrs,
                         const std::string &output_file_name, StopWatch *timer) {
  const int table_construction = timer->addCategory("Log Spline Table Building");
  timer->assignTime(0);
  LogScaleSpline<T4> lgsp(target_form_in, ew_coeff, stormm::symbols::charmm_gromacs_bioq,
                          mantissa_bits_in, 4096.0, 0.015625, indexing_method_in, basis_set_in,
                          0, indexing_offset_in);
  timer->assignTime(table_construction);

  // Prepare to write the summary
  std::string var_name, unit_str, title_str;
  setDisplayStrings(lgsp, &var_name, &unit_str, &title_str);
  std::ofstream foutp;
  std::vector<double> rpts, flt_eval_mue, spl_eval_mue;
  std::vector<double> test_data(2 * (max_ulp + 1));
  var_name += "_idxb";
  for (int i = 0; i <= max_ulp; i++) {
    timer->assignTime(0);
    LogScaleSpline<T4> lgsp(target_form_in, ew_coeff, stormm::symbols::charmm_gromacs_bioq,
                            mantissa_bits_in, 4096.0, 0.015625, indexing_method_in, basis_set_in,
                            i, indexing_offset_in);
    timer->assignTime(table_construction);
    computeSplineErrors(lgsp, min_analysis_range, max_analysis_range, false, &rpts,
                        &flt_eval_mue, &spl_eval_mue, xrs);
    if (i == 0) {
      foutp = prepOutputFile(output_file_name, lgsp);
    }
    const int npts = spl_eval_mue.size();
    double weighted_error = 0.0;
    for (int j = 0; j < npts; j++) {
      weighted_error += spl_eval_mue[j] * rpts[j] * rpts[j];
    }
    test_data[i] = i;
    test_data[i + (max_ulp + 1)] = weighted_error / static_cast<double>(npts);
  }
  const std::vector<int> dec_places = { 0, 12 };
  ReportTable test_tab(test_data, { "ULP Perturbations", "Weighted Mean Unsigned Error" },
                       dec_places, var_name, 120);
  test_tab.printTable(&foutp, OutputSyntax::MATRIX_PKG);
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
// Tabulate the mean unsigned error inherent in logarithmic splines for different numbers of bits
// taken out of the mantissa.
//
// Arguments:
//   ew_coeff:          The Ewald coefficient selected for this round of tests
//   xrs:               Source of random numbers for spline evaluations
//   timer:             Tracks the time required to make various splines
//   output_file_name:  Name of the output file into which to funnel results
//-------------------------------------------------------------------------------------------------
void testMantissaBits(const double ew_coeff, Xoshiro256ppGenerator *xrs, StopWatch *timer,
                      const std::string &output_file_name) {
  const int table_construction = timer->addCategory("Log Spline Table Building");
  const std::vector<BasisFunctions> tbasis = { BasisFunctions::POLYNOMIAL,
                                               BasisFunctions::MIXED_FRACTIONS };
  const std::vector<TableIndexing> tindex = { TableIndexing::ARG, TableIndexing::SQUARED_ARG };
  std::vector<double> test_data;
  for (int i = 4; i <= 7; i++) {
    test_data.push_back(i);
  }
  std::vector<double> rpts, flt_eval_mue, spl_eval_mue;
  std::vector<std::string> legend;
  legend.push_back(std::string("Mantissa Bits"));
  for (size_t bidx = 0; bidx < tbasis.size(); bidx++) {
    for (size_t tidx = 0; tidx < tindex.size(); tidx++) {
      legend.push_back(getEnumerationName(tbasis[bidx]) + " / " +
                       getEnumerationName(tindex[tidx]));
      for (int i = 4; i <= 7; i++) {
        timer->assignTime(0);
        LogScaleSpline<float4> lgsp(LogSplineForm::DELEC_PME_DIRECT, ew_coeff,
                                    stormm::symbols::charmm_gromacs_bioq, i, 4096.0, 0.015625,
                                    tindex[tidx], tbasis[bidx], 0, 0.0);
        timer->assignTime(table_construction);
        computeSplineErrors(lgsp, 0.5, 12.0, false, &rpts, &flt_eval_mue, &spl_eval_mue, xrs);
        const int npts = spl_eval_mue.size();
        double weighted_error = 0.0;
        for (int j = 0; j < npts; j++) {
          weighted_error += spl_eval_mue[j] * rpts[j] * rpts[j];
        }
        test_data.push_back(weighted_error / static_cast<double>(npts));
      }
    }
  }
  std::ofstream foutp = openOutputFile(output_file_name, PrintSituation::APPEND, "open the output "
                                       "for additional test recording");
  std::string buffer("\n");
  foutp.write(buffer.data(), buffer.size());
  const std::vector<int> dec_places = { 0, 12, 12, 12, 12 };
  const std::string var_name = "mantissa_dependence";
  
  ReportTable test_tab(test_data, legend, dec_places, var_name, 120);
  test_tab.printTable(&foutp, OutputSyntax::MATRIX_PKG);
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Baseline variables
  StopWatch timer;
  
  // Take in additional inputs
  const std::vector<TableIndexing> tidx_methods = { TableIndexing::SQUARED_ARG, TableIndexing::ARG,
                                                    TableIndexing::SQ_ARG_OFFSET,
                                                    TableIndexing::ARG_OFFSET };
  const std::vector<LogSplineForm> lsfrm_methods = { LogSplineForm::DELEC_PME_DIRECT,
                                                     LogSplineForm::DELEC_PME_DIRECT_EXCL };
  std::vector<bool> disp_tidx(tidx_methods.size(), false);
  std::vector<bool> disp_lsfrm(lsfrm_methods.size(), false);

  // Take in command line inputs
  CommandLineParser clip("erfc_tables", "A benchmarking program which analyzes error levels in "
                         "spline tables for interpolating the electrostatic particle-particle "
                         "interaction in PME calculations.", { "-timings" });
  clip.addStandardAmberInputs("-ig_seed");
  clip.suppressHelpOnNoArgs();
  NamelistEmulator *t_nml= clip.getNamelistPointer();
  t_nml->addKeyword("-o", NamelistType::STRING, std::string("erfc_table_results.m"));
  t_nml->addKeyword("-cut", NamelistType::REAL, std::to_string(10.0));
  t_nml->addHelp("-cut", "The cutoff to apply to the spline tables, in units of Angstroms.  All "
                 "particle-particle displacments less than the cutoff will have "
                 "spline-interpolated values.");
  t_nml->addKeyword("-dsum_tol", NamelistType::REAL, std::to_string(default_dsum_tol));
  t_nml->addHelp("-dsum_tol", "Direct sum tolerance for the PME calculation ,stipulating how much "
                 "of the electrostatic interaction will be discarded.");
  t_nml->addKeyword("-mnbits", NamelistType::INTEGER, std::to_string(5));
  t_nml->addHelp("-mnbits", "The number of bits of the mantissa to use (in conjunction with the "
                 "exponent bits) when producing an index key into the spline tables.");
  t_nml->addKeyword("-ulp", NamelistType::INTEGER, std::to_string(4));
  t_nml->addHelp("-ulp", "The number of bits of least precision to optimize such that spline "
                 "values computed in float32_t will better match those computed in float64_t.");
  t_nml->addKeyword("-func", NamelistType::STRING,
                    getEnumerationName(LogSplineForm::DELEC_PME_DIRECT), DefaultIsObligatory::NO,
                    InputRepeats::YES);
  t_nml->addHelp("-func", "Form of the function to express in splines.  Acceptable values include "
                 "ELEC_PME (electrostatic particle-particle potential used in PME), DELEC_PME "
                 "(derivative of the electrostatic potential), ELEC_PME_EXCL (electrostatic PME "
                 "potential, with exclusion between particles in the primary image), and "
                 "DELEC_PME_EXCL. (These specifications are case-insensitive, though the keyword "
                 "is case-sensitive.) This may be specified multiple times to display multiple "
                 "functional forms.");
  t_nml->addKeyword("-index", NamelistType::STRING,
                    getEnumerationName(TableIndexing::ARG), DefaultIsObligatory::NO,
                    InputRepeats::YES);
  t_nml->addHelp("-index", "The table indexing method used to prepare logarithmic spline tables.  "
                 "The acceptable forms include ARG (the table will be indexed by the actual value "
                 "of the inter-particle displacement), SQUARED_ARG (the table will be indexed by "
                 "the square of the inter-particle displacement), ARG_OFFSET (the table index "
                 "will include a flat offset added to the inter-particle displacement), or "
                 "SQ_ARG_OFFSET.  These values are case insensitive.  The keyword may be "
                 "specified multiple times in order to display results for mutliple table "
                 "indexing methods.");

  // Load the testing environment and have it cooperate with the=is program's own CommandLineParser
  // to read user input.
  TestEnvironment oe(argc, argv, &clip, TmpdirStatus::NOT_REQUIRED, ExceptionResponse::SILENT);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Create a Hybrid object to engage the GPU and absorb any bootup time into "miscellaneous"
  if (oe.getDisplayTimingsOrder()) {
    Hybrid<int> gpu_trigger(1);
    timer.assignTime(0);
  }

  // Parse user input
  clip.parseUserInput(argc, argv);
  const double cutoff = t_nml->getRealValue("-cut");
  const double dsum_tol = t_nml->getRealValue("-dsum_tol");
  const int mantissa_bits = t_nml->getIntValue("-mnbits");
  const int ig_seed = t_nml->getIntValue("-ig_seed");
  const int max_ulp = t_nml->getIntValue("-ulp");
  const std::string output_file_name = t_nml->getStringValue("-o");
  for (int i = 0; i < t_nml->getKeywordEntries("-func"); i++) {
    const std::string lgform = t_nml->getStringValue("-func", i);
    try {
      const LogSplineForm test_lsfrm = translateLogSplineForm(lgform);
      for (size_t j = 0; j < lsfrm_methods.size(); j++) {
        if (test_lsfrm == lsfrm_methods[j]) {
          disp_lsfrm[j] = true;
        }
      }
    }
    catch (std::runtime_error) {
      rtWarn("\"" + lgform + "\" was not recognized as any of the LogScaleSpline functions "
             "covered by these benchmarks.", "main");
    }
  }
  for (int i = 0; i < t_nml->getKeywordEntries("-index"); i++) {
    const std::string tbmeth = t_nml->getStringValue("-index", i);
    try {
      const TableIndexing test_tidx = translateTableIndexing(tbmeth);
      for (size_t j = 0; j < tidx_methods.size(); j++) {
        if (test_tidx == tidx_methods[j]) {
          disp_tidx[j] = true;
        }
      }
    }
    catch (std::runtime_error) {
      rtWarn("\"" + tbmeth + "\" was not recognized as any of the table indexing methods  "
             "covered by these benchmarks.", "main");
    }
  }

  // Initialize the random number generator.
  Xoshiro256ppGenerator xrs(ig_seed);

  // Initialize the output.
  std::ofstream foutp = openOutputFile(output_file_name, PrintSituation::OVERWRITE, "prime the "
                                       "output file for printing");
  std::string primer("%% Clear the decks\nclear all\nclose all\n\n");
  foutp.write(primer.data(), primer.size());
  foutp.close();
  
  // Compute the Ewald coefficient.
  const double ew_coeff = ewaldCoefficient(cutoff, dsum_tol);
  for (size_t i = 0; i < tidx_methods.size(); i++) {
    if (disp_tidx[i] == false) {
      continue;
    }
    float idx_offset;
    switch (tidx_methods[i]) {
    case TableIndexing::ARG:
    case TableIndexing::SQUARED_ARG:
      idx_offset = 0.0;
      break;
    case TableIndexing::ARG_OFFSET:
    case TableIndexing::SQ_ARG_OFFSET:
      idx_offset = 0.0625;
      break;
    }
    for (size_t j = 0; j < lsfrm_methods.size(); j++) {
      if (disp_lsfrm[j] == false) {
        continue;
      }
      analyzeTable<float4>(lsfrm_methods[j], ew_coeff, 0.5, 12.0, mantissa_bits, tidx_methods[i],
                           BasisFunctions::POLYNOMIAL, idx_offset, "natural",
                           disp_tidx[i] && disp_lsfrm[j], &xrs, output_file_name, &timer);
      analyzeTable<float4>(lsfrm_methods[j], ew_coeff, 0.5, 12.0, mantissa_bits, tidx_methods[i],
                           BasisFunctions::POLYNOMIAL, idx_offset, "lab",
                           disp_tidx[i] && disp_lsfrm[j], &xrs, output_file_name, &timer);
      analyzeTable<float4>(lsfrm_methods[j], ew_coeff, 0.5, 12.0, mantissa_bits, tidx_methods[i],
                           BasisFunctions::MIXED_FRACTIONS, idx_offset, "natural",
                           disp_tidx[i] && disp_lsfrm[j], &xrs, output_file_name, &timer);
      analyzeTable<float4>(lsfrm_methods[j], ew_coeff, 0.5, 12.0, mantissa_bits, tidx_methods[i],
                           BasisFunctions::MIXED_FRACTIONS, idx_offset, "lab",
                           disp_tidx[i] && disp_lsfrm[j], &xrs, output_file_name, &timer);

      // Only perform ULP optimization for the non-excluded force, ARG or SQUARED_ARG indexing
      switch (lsfrm_methods[j]) {
      case LogSplineForm::ELEC_PME_DIRECT:
      case LogSplineForm::ELEC_PME_DIRECT_EXCL:
      case LogSplineForm::DELEC_PME_DIRECT_EXCL:
      case LogSplineForm::CUSTOM:
        continue;
      case LogSplineForm::DELEC_PME_DIRECT:
        break;
      }
      switch (tidx_methods[i]) {
      case TableIndexing::ARG_OFFSET:
      case TableIndexing::SQ_ARG_OFFSET:
        continue;
      case TableIndexing::ARG:
      case TableIndexing::SQUARED_ARG:
        break;
      }
      testUlpOptimization<float4>(lsfrm_methods[j], ew_coeff, 0.5, 12.0, mantissa_bits,
                                  tidx_methods[i], BasisFunctions::POLYNOMIAL, idx_offset, max_ulp,
                                  &xrs, output_file_name, &timer);
      testUlpOptimization<float4>(lsfrm_methods[j], ew_coeff, 0.5, 12.0, mantissa_bits,
                                  tidx_methods[i], BasisFunctions::MIXED_FRACTIONS, idx_offset,
                                  max_ulp, &xrs, output_file_name, &timer);
    }
  }

  // Aggregate the results of ULP optimization tests
  foutp = openOutputFile(output_file_name, PrintSituation::APPEND, "prime the output file for "
                         "printing");
  std::string buffer;
  const std::vector<BasisFunctions> spline_basis = { BasisFunctions::POLYNOMIAL,
                                                     BasisFunctions::MIXED_FRACTIONS };
  const std::vector<PsivantColor> spline_colors = { PsivantColor::BLACK, PsivantColor::RED,
                                                    PsivantColor::YELLOW, PsivantColor::PURPLE };
  for (size_t i = 0; i < lsfrm_methods.size(); i++) {
    switch (lsfrm_methods[i]) {
    case LogSplineForm::ELEC_PME_DIRECT:
    case LogSplineForm::ELEC_PME_DIRECT_EXCL:
    case LogSplineForm::DELEC_PME_DIRECT_EXCL:
    case LogSplineForm::CUSTOM:
      continue;
    case LogSplineForm::DELEC_PME_DIRECT:
      break;
    }
    buffer = "\nfigure;\nhold on;\n";
    foutp.write(buffer.data(), buffer.size());
    std::string legend("legend(");
    int nlegend_entry = 0;
    for (size_t bss_con = 0; bss_con < spline_basis.size(); bss_con++) {
      for (size_t j = 0; j < tidx_methods.size(); j++) {
        switch (tidx_methods[j]) {
        case TableIndexing::ARG:
        case TableIndexing::SQUARED_ARG:
          break;
        case TableIndexing::ARG_OFFSET:
        case TableIndexing::SQ_ARG_OFFSET:
          continue;
        }
        const std::string var_name = setDisplayVariable(lsfrm_methods[i], tidx_methods[j],
                                                        spline_basis[bss_con], "_idxb");
        buffer = "plot(" + var_name + "(:,1), " + var_name + "(:,2), " +
                 encodePsivantColor(spline_colors[nlegend_entry], OutputSyntax::MATRIX_PKG) +
                 ", 'linewidth', 8);\n"; 
        buffer += "plot(" + var_name + "(:,1), " + var_name + "(:,2), '.', 'markersize', 48, " +
                  encodePsivantColor(spline_colors[nlegend_entry], OutputSyntax::MATRIX_PKG) +
                  ", 3);\n"; 
        buffer += "plot(" + var_name + "(:,1), " + var_name + "(:,2), 'w.', 'markersize', 32);\n"; 
        foutp.write(buffer.data(), buffer.size());
        if (nlegend_entry == 0) {
          legend += "'" + getEnumerationName(tidx_methods[j]) + ", " +
                    getEnumerationName(spline_basis[bss_con]) + "'";
        }
        else {
          legend += ", '" + getEnumerationName(tidx_methods[j]) + ", " +
                    getEnumerationName(spline_basis[bss_con]) + "'";
        }
        nlegend_entry += 1;
      }
    }
    buffer = "xlabel('ULP Optimization Range');\nylabel('Error, kcal/mol-A');\n";
    buffer += "axis([ 0 " + std::to_string(max_ulp) + " 1.0e-6 4.0e-6 ]);\n";
    buffer += "daspect([ 2.0 1.0e-6 1.0 ]);\n";
    legend += ");\n";
    foutp.write(legend.data(), legend.size());
    buffer += "set(gca, 'fontsize', 36);\nlegend('Location', 'southeastoutside');\n";
    foutp.write(buffer.data(), buffer.size());
  }
  foutp.close();
  testMantissaBits(ew_coeff, &xrs, &timer, output_file_name);
  
  // Create some cell grids out of periodic systems and compute direct-space interactions.  Test
  // the spline tables in "real-world" applications.
  const double half_cutoff = 0.5 * cutoff;
  const std::vector<std::string> systems = { "trpcage_in_water", "drug_example", "ubiquitin" };
  const char osc = osSeparator();
  const std::string base_top_path = oe.getStormmSourcePath() + osc + "test" + osc + "Topology";
  const std::string base_crd_path = oe.getStormmSourcePath() + osc + "test" + osc + "Trajectory";
  TestSystemManager tsm(base_top_path, "top", systems, base_crd_path, "inpcrd", systems);
  for (int i = 0; i < tsm.getSystemCount(); i++) {
    PhaseSpace ps = tsm.exportPhaseSpace(i);
    PhaseSpaceWriter psw = ps.data();
    double cell_a, cell_b, cell_c;
    hessianNormalWidths(psw.invu, &cell_a, &cell_b, &cell_c);
    const int na_cell = cell_a / half_cutoff;
    const int nb_cell = cell_b / half_cutoff;
    const int nc_cell = cell_c / half_cutoff;
  }
  
  // Print results
  if (oe.getDisplayTimingsOrder()) {
    timer.assignTime(0);
    timer.printResults();
  }
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
}
