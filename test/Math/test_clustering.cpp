#include <cmath>
#include "copyright.h"
#include "../../src/Constants/symbol_values.h"
#include "../../src/Math/clustering.h"
#include "../../src/Math/cluster_manager.h"
#include "../../src/Math/vector_ops.h"
#include "../../src/Random/random.h"
#include "../../src/Reporting/summary_file.h"
#include "../../src/Structure/structure_enumerators.h"
#include "../../src/UnitTesting/approx.h"
#include "../../src/UnitTesting/test_environment.h"
#include "../../src/UnitTesting/test_system_manager.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::stmath;
using namespace stormm::structure;
using namespace stormm::testing;
using namespace stormm::trajectory;

//-------------------------------------------------------------------------------------------------
// Prepare a data set and then cluster it.  Check that the points fall into the correct clusters
// based on the locations of the chosen centroids.
//
// Arguments:
//   npoints:       The number of points to place in the data set
//   nclusters:     The number of clusters to form
//   xrs:           Random number generator with which to generate the data set
//   boundary:      Boundary conditions on the data set
//   range:         The range over which to distribute data
//   beats:         Make beating structures in the data which clustering should pick up on
//   scale_factor:  The scaling factor applied to all data points, if in fixed-precision format
//-------------------------------------------------------------------------------------------------
template<typename Tdata, typename Tcalc>
void clusterOneDimData(const int npoints, const int nclusters, Xoshiro256ppGenerator *xrs,
                       const BoundaryCondition boundary, const double range,
                       const bool beats = false, const Tcalc scale_factor = 1.0) {
  std::vector<double> r_distro = uniformRand(xrs, npoints, 1.0);
  elementwiseMultiply(&r_distro, range);
  ClusterManager<double, double> cls_obj(npoints, nclusters);
  kMeans(&cls_obj, r_distro, nclusters, 1.0, boundary, range);
  check(cls_obj.getClusterCount(), RelationalOperator::EQUAL, nclusters, "The number of clusters "
        "in a ClusterManager does not meet expectations.");
  int total_points = 0;
  for (int i = 0; i < cls_obj.getClusterCount(); i++) {
    total_points += cls_obj.getClusterSize(i);
  }
  check(npoints, RelationalOperator::EQUAL, total_points, "The total number of points spread "
        "across all clusters does not match the size of the data set.");
  bool clusters_violated = false;
  std::vector<size_t> manual_primes(nclusters), obj_primes(nclusters);
  std::vector<Tcalc> manual_centroids(nclusters), obj_centroids(nclusters);
  const Tcalc zero = 0.0;
  for (int i = 0; i < cls_obj.getClusterCount(); i++) {
    const std::vector<size_t> clmemi = cls_obj.getClusterMembers(i);
    size_t iprime;
    Tcalc best_prime_distance = 2.0 * range;
    Tcalc icent = zero;
    const Tcalc obj_icent = cls_obj.getClusterCentroid(i);
    for (size_t j = 0; j < clmemi.size(); j++) {
      Tcalc cdist;
      const Tdata djval = r_distro[clmemi[j]]; 
      switch (boundary) {
      case BoundaryCondition::ISOLATED:
        cdist = std::abs(djval - obj_icent);
        icent += djval;
        break;
      case BoundaryCondition::PERIODIC:
        cdist = std::abs(imageValue(djval - obj_icent, range,
                                    ImagingMethod::MINIMUM_IMAGE));
        icent += imageValue(djval - obj_icent, range, ImagingMethod::MINIMUM_IMAGE, scale_factor) +
                 obj_icent;
        break;
      }
      if (cdist < best_prime_distance) {
        best_prime_distance = cdist;
        iprime = clmemi[j];
      }
      int best_centroid;
      Tcalc best_distance = 2.0 * range;
      for (int k = 0; k < nclusters; k++) {
        switch (boundary) {
        case BoundaryCondition::ISOLATED:
          cdist = std::abs(djval - cls_obj.getClusterCentroid(k));
          break;
        case BoundaryCondition::PERIODIC:
          cdist = std::abs(imageValue(djval - cls_obj.getClusterCentroid(k), range,
                                      ImagingMethod::MINIMUM_IMAGE));
          break;
        }
        if (cdist < best_distance) {
          best_distance = cdist;
          best_centroid = k;
        }
      }
      clusters_violated = (clusters_violated || best_centroid != i);
    }
    manual_primes[i] = iprime;
    obj_primes[i] = cls_obj.getClusterPrime(i);
    manual_centroids[i] = icent / static_cast<Tcalc>(cls_obj.getClusterSize(i));
    obj_centroids[i] = obj_icent;
  }
  check(clusters_violated == false, "The assignment of data points to various clusters does not "
        "minimize the local sum of squares (one or more points would be better assigned to a "
        "cluster with a nearer centroid).  Data type: " + getStormmScalarTypeName<Tdata>() + ".");
  check(manual_primes, RelationalOperator::EQUAL, obj_primes, "The prime members of one or more "
        "clusters are incorrect.");
  check(manual_centroids, RelationalOperator::EQUAL, obj_centroids, "The centroids of one or more "
        "clusters are incorrect.");
}

//-------------------------------------------------------------------------------------------------
// main
//-------------------------------------------------------------------------------------------------
int main(const int argc, const char* argv[]) {

  // Some baseline initialization
  TestEnvironment oe(argc, argv);
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmSplash();
  }

  // Section 1: cluster one-dimensonal data
  section("Clustering in one-dimension");
  
  // Check the clustering of one-dimensional data.
  section(1);
  Xoshiro256ppGenerator xrs(7183529);
  clusterOneDimData<double, double>(200, 5, &xrs, BoundaryCondition::PERIODIC,
                                    stormm::symbols::twopi, false);
  clusterOneDimData<double, double>(2000, 10, &xrs, BoundaryCondition::PERIODIC,
                                    stormm::symbols::twopi, false);
  clusterOneDimData<double, double>(200, 5, &xrs, BoundaryCondition::ISOLATED,
                                    stormm::symbols::twopi, false);
  clusterOneDimData<double, double>(100, 2, &xrs, BoundaryCondition::ISOLATED,
                                    stormm::symbols::twopi, false);
  clusterOneDimData<double, double>(200, 2, &xrs, BoundaryCondition::PERIODIC,
                                    stormm::symbols::twopi, false);
  
  // Print results
  printTestSummary(oe.getVerbosity());
  if (oe.getVerbosity() == TestVerbosity::FULL) {
    stormmWatermark();
  }
  return countGlobalTestFailures();
}
