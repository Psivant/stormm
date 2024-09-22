#include <cmath>
#include "copyright.h"
#include "DataTypes/common_types.h"
#include "Math/vector_ops.h"
#include "Parsing/parse.h"
#include "vector_report.h"
#include "approx.h"

namespace stormm {
namespace testing {

using stmath::maxValue;
using stmath::mean;
using stmath::meanUnsignedError;
using stmath::pearson;
using stmath::variance;
using stmath::VarianceMethod;
using parse::realDecimalPlaces;
using parse::realToString;

//-------------------------------------------------------------------------------------------------
std::string vectorAlignmentReport(const std::vector<PolyNumeric> &va,
                                  const std::vector<PolyNumeric> &vb, NumberFormat data_format,
                                  const double tol) {

  // Initialize the output
  std::string result;

  // Bail out immediately if the vectors are not the same length
  if (va.size() != vb.size()) {
    result += "The vectors are of different sizes (" + std::to_string(va.size()) + " and " +
              std::to_string(vb.size()) + ").";
    return result;
  }
  
  // Different numerical formats (plus char4)
  switch (data_format) {
  case NumberFormat::SCIENTIFIC:
  case NumberFormat::STANDARD_REAL:
  case NumberFormat::INTEGER:
  case NumberFormat::LONG_LONG_INTEGER:
  case NumberFormat::UNSIGNED_INTEGER:
  case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
    {
      // Convert both vectors to double precision real for analysis
      std::vector<double> dva, dvb;
      switch (data_format) {
      case NumberFormat::SCIENTIFIC:
      case NumberFormat::STANDARD_REAL:
        dva = doubleFromPolyNumeric(va);
        dvb = doubleFromPolyNumeric(vb);
        break;
      case NumberFormat::INTEGER:
        {
          const std::vector<int> iva = intFromPolyNumeric(va);
          const std::vector<int> ivb = intFromPolyNumeric(vb);
          dva = std::vector<double>(iva.begin(), iva.end());
          dvb = std::vector<double>(ivb.begin(), ivb.end());
        }
        break;
      case NumberFormat::LONG_LONG_INTEGER:
        {
          const std::vector<llint> lliva = llintFromPolyNumeric(va);
          const std::vector<llint> llivb = llintFromPolyNumeric(vb);
          dva = std::vector<double>(lliva.begin(), lliva.end());
          dvb = std::vector<double>(llivb.begin(), llivb.end());
        }
        break;
      case NumberFormat::UNSIGNED_INTEGER:
        {
          const std::vector<uint> uiva = uintFromPolyNumeric(va);
          const std::vector<uint> uivb = uintFromPolyNumeric(vb);
          dva = std::vector<double>(uiva.begin(), uiva.end());
          dvb = std::vector<double>(uivb.begin(), uivb.end());
        }
        break;
      case NumberFormat::UNSIGNED_LONG_LONG_INTEGER:
        {
          const std::vector<ullint> ulliva = ullintFromPolyNumeric(va);
          const std::vector<ullint> ullivb = ullintFromPolyNumeric(vb);
          dva = std::vector<double>(ulliva.begin(), ulliva.end());
          dvb = std::vector<double>(ullivb.begin(), ullivb.end());
        }
        break;
      case NumberFormat::CHAR4:
        break;
      }      
      std::vector<double> work(std::max(va.size(), vb.size()), 0.0);
      const int n_va = va.size();
      const int n_vb = vb.size();
      
      // Look for a common multiple if the vectors are of the same size
      if (n_va == n_vb) {
        for (int i = 0; i < n_va; i++) {
          work[i] = dva[i] / dvb[i];
        }
        if (variance(work, VarianceMethod::NORMALIZED_RMSD) < tol && n_va > 2) {
          result += "The vectors appear to be multiples of one another: first vector = " +
                    realToString(mean(work)) + " x second vector.";
        }
        else {
          int n_match = 0;
          int n_mismatch = 0;
          std::vector<double> vec_differences;
          for (int i = 0; i < n_va; i++) {
            if (dva[i] == Approx(dvb[i], ComparisonType::ABSOLUTE, tol)) {
              n_match++;
            }
            else {
              n_mismatch++;
              vec_differences.push_back(fabs(dva[i] - dvb[i]));
            }
          }
          const NumberFormat scifm = NumberFormat::SCIENTIFIC;
          if (n_va > 1) {
            if (n_mismatch > 0 && n_match > 0) {
              result += "The vectors fail to match in " + std::to_string(n_mismatch) +
                        " indices out of " + std::to_string(n_va) + ", with a mean unsigned error "
                        "(among the deviating entries) of " +
                        realToString(mean(vec_differences), 11, 4, scifm) +
                        " and maximum unsigned deviation " +
                        realToString(maxValue(vec_differences), 11, 4, scifm) + ".  ";
            }
            else {
              result += "Deviations occur throughout the data.  ";
            }
            if (n_va > 2) {
              if (dva.size() > 1LLU) {
                result += "Pearson correlation between the vectors is " +
                          realToString(pearson(dva, dvb), 4) + ".  ";
              }
            }
          }
          if (n_mismatch > 0) {
            std::vector<double2> d_vabn(n_va);
            for (int j = 0; j < n_va; j++) {
              d_vabn[j] = { fabs(dva[j] - dvb[j]), static_cast<double>(j) };
            }
            std::sort(d_vabn.begin(), d_vabn.end(),
                      []( double2 a, double2 b ) { return a.x > b.x; });
            std::vector<int> discrepancy_ids(n_va);
            for (int j = 0; j < n_va; j++) {
              discrepancy_ids[j] = d_vabn[j].y;
            }
            result += "Mismatched entries:\n";
            int n_reported = 0;
            const int ndec = realDecimalPlaces(tol);
            int i = 0;
            while (n_reported < 16 && n_reported < n_mismatch && i < n_va) {
              if (dva[discrepancy_ids[i]] != Approx(dvb[discrepancy_ids[i]],
                                                    ComparisonType::ABSOLUTE, tol)) {
                result += "    " + realToString(dva[discrepancy_ids[i]], ndec + 7, ndec, scifm) +
                          " != " + realToString(dvb[discrepancy_ids[i]], ndec + 7, ndec, scifm) +
                          " (error " + realToString(fabs(dvb[discrepancy_ids[i]] -
                                                         dva[discrepancy_ids[i]]), ndec + 7, ndec,
                                                    scifm) +
                          ", entry " + std::to_string(discrepancy_ids[i]) + ")\n";
                n_reported++;
              }
              i++;
            }
          }
	}
      }
    }
    break;
  case NumberFormat::CHAR4:
    rtErr("Vector analysis method not yet implemented.", "vectorAlignmentReport");
    break;
  }
  return result;
}

} // namespace testing
} // namespace stormm
