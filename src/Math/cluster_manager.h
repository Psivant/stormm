// -*-c++-*-
#ifndef STORMM_CLUSTER_MANAGER_H
#define STORMM_CLUSTER_MANAGER_H

#include <vector>
#include "copyright.h"
#include "Accelerator/hybrid.h"
#include "Constants/hpc_bounds.h"
#include "Synthesis/condensate.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Trajectory/coordinateframe.h"
#include "Trajectory/coordinate_series.h"
#include "Trajectory/phasespace.h"
#include "rounding.h"

namespace stormm {
namespace stmath {

using card::Hybrid;
using card::HybridKind;
using card::HybridTargetLevel;
using synthesis::Condensate;
using synthesis::PhaseSpaceSynthesis;
using trajectory::CoordinateFrame;
using trajectory::PhaseSpace;
using trajectory::CoordinateSeries;

/// \brief An abstract of the cluster manager suitable for rapid access and modification in the C++
///        layer, or submission as an argument to an HPC kernel.  As with any collection of C-style
///        pointers, the abstract is only valid so long as the object it references does not get
///        resized.
template <typename Tcalc>
struct ClusterManagerKit {
public:

  /// \brief The constructor takes arguments for every member variable.
  ClusterManagerKit(size_t data_points_in, int clusters_in, int* cls_homes_in,
                    size_t* cls_membership_in, size_t* cls_bounds_in, size_t* cls_primes_in,
                    Tcalc* cls_rmsd_in, Tcalc* centroid_dist_in);

  /// \brief The const sizing elements mean that the default copy and move assignment operators
  ///        will be implicitly deleted.
  /// \{
  ClusterManagerKit(const ClusterManagerKit &original) = default;
  ClusterManagerKit(ClusterManagerKit &&original) = default;
  /// \}

  const size_t data_points;  ///< The number of data points in the entire set
  const int clusters;        ///< The number of clusters into which data points are grouped
  int* cls_homes;            ///< The clusters to which each data point belongs (in the current
                             ///<   clustering arrangement)
  size_t* cls_membership;    ///< Concatenated lists of the data points each cluster comprises
  size_t* cls_bounds;        ///< Bounds array on cls_membership
  size_t* cls_primes;        ///< Members of each cluster nearest the cluster's centroid
  Tcalc* cls_rmsd;           ///< Root mean squared deviation in the spread of points in each
                             ///<   cluster
  Tcalc* centroid_dist;      ///< Distance of each data point to its cluster's centroid
};

/// \brief An object to track the clusters into which a set of data points can be sorted.
template <typename Tdata, typename Tcalc>
class ClusterManager {
public:

  /// \brief The constructor accepts a number of data points and a number of clusters, as well as
  ///        the number of seeding attempts, if available.
  explicit ClusterManager(size_t data_point_count_in = 0, int cluster_count_in = 1,
                          int attempt_count_in = 0);

  /// \brief The copy and move constructors, as well as copy and move assignment operators, must
  ///        be explicitly instantiated in order to repair POINTER-kind Hybrid objects among the
  ///        member variables.
  ///
  /// \param original  The object to copy or move
  /// \param other     A pre-existing object to place on the right-hand side of the assignment
  ///                  statement
  /// \{
  ClusterManager(const ClusterManager &original);
  ClusterManager(ClusterManager &&original);
  ClusterManager& operator=(const ClusterManager &original);
  ClusterManager& operator=(ClusterManager &&original);
  /// \}

  /// \brief Get the number of data points.
  size_t getDataPointCount() const;

  /// \brief Get the number of clusters.
  int getClusterCount() const;
  
  /// \brief Get the number of attempts made to cluster the data resulting in the reported
  ///        arrangement.
  int getAttemptCount() const;
  
  /// \brief Get the cluster to which a specific data point belongs.
  ///
  /// \param point_index  Index of the data point of interest
  int getClusterHome(size_t point_index) const;

  /// \brief Get the size of a cluster.
  ///
  /// \param cluster_index  Index of the cluster of interest
  size_t getClusterSize(int cluster_index) const;

  /// \brief Get the cluster member nearest the cluster centroid.
  ///
  /// \param cluster_index  Index of the cluster of interest
  size_t getClusterPrime(int cluster_index) const;

  /// \brief Get all members of a cluster.
  ///
  /// \param cluster_index  Index of the cluster of interest
  std::vector<size_t> getClusterMembers(int cluster_index) const;

  /// \brief Get the centroid for a particular cluster.
  ///
  /// \param cluster_index  Index of the cluster of interest
  const Tdata& getClusterCentroid(int cluster_index) const;

  /// \brief Get a const references to the vector of centroids for all clusters.
  const std::vector<Tdata>& getClusterCentroids() const;
  
  /// \brief Get the abstract.
  ///
  /// Overloaded:
  ///   - Get a read-only abstract from a const object.
  ///   - Get a non-const, writeable abstract from a non-const object.
  ///
  /// \param tier  Set points to data on the GPU device or CPU host
  /// \{
  const ClusterManagerKit<Tcalc> data(HybridTargetLevel tier = HybridTargetLevel::HOST) const;
  ClusterManagerKit<Tcalc> data(HybridTargetLevel tier = HybridTargetLevel::HOST);
  /// \}
  
  /// \brief Set the number of clusters and data points.
  ///
  /// \param data_point_count_in  The new number of data points
  /// \param cluster_count_in     The new number of clusters to partition data points into
  void resize(size_t data_point_count_in, int cluster_count_in);
  
  /// \brief Set the number of attempts.
  ///
  /// \param attempt_count_in  The number of attempted clusterings to log
  void setAttemptCount(int attempt_count_in);

  /// \brief Set the value of one cluster's centroid.
  ///
  /// \param cluster_index  Index of the cluster to set
  /// \param value          Value of the cluster's new centroid
  void setCentroid(int cluster_index, const Tdata value);
  
private:
  size_t data_point_count;            ///< The total number of data points
  int cluster_count;                  ///< Total number of clusters
  int attempt_count;                  ///< The number of attempts used to arrive at the clustering
                                      ///<   result stored in the object
  Hybrid<int> cluster_homes;          ///< The clusters to which each data point belongs
  Hybrid<size_t> cluster_membership;  ///< The membership of each cluster, which can be queried to
                                      ///<   determine which data points a given cluster contains.
  Hybrid<size_t> cluster_bounds;      ///< Bounds array for cluster_membership
  Hybrid<size_t> cluster_primes;      ///< Indices of the data points nearest to each cluster
                                      ///<   centroid
  Hybrid<Tcalc> cluster_rmsd;         ///< Root-mean-squared deviations among the members of each
                                      ///<   cluster, resulting from the distance metric
                                      ///<   implemented in the kmeans() function (see
                                      ///<   series_ops.h)
  Hybrid<Tcalc> centroid_distances;   ///< Distances of each cluster member to its respective
                                      ///<   cluster centroid.  This array follows the indexing of
                                      ///<   the original data items.
  Hybrid<size_t> sizet_data;          ///< ARRAY-kind Hybrid object targeted by the preceding
                                      ///<   size_t type Hybrid member variables
  Hybrid<Tcalc> real_data;            ///< ARRAY-kind Hybrid object targeted by the preceding
                                      ///<   Tcalc type (real) Hybrid member variables
  std::vector<Tdata> centroids;       ///< A vector of centroids for each cluster of data points.
                                      ///<   The centroids have the same type as the original data
                                      ///<   points, but are not part of the original data even if
                                      ///<   they coincide perfectly with specific data points.
                                      ///<   Free functions in the implementation file
                                      ///<   cluster_manager.cpp handle specific data types.  While
                                      ///<   this member variable is a std::vector<T>, its elements
                                      ///<   are not restricted from having Hybrid memory of their
                                      ///<   own, which would then be accessible to the GPU.  In
                                      ///<   selected cases, the templated type may, itself, imply
                                      ///<   an array of structures (e.g. a PhaseSpaceSynthesis).
                                      ///<   The associated free functions can then call the first
                                      ///<   element of the centroids array to open up the internal
                                      ///<   array underneath.

  /// \brief Allocate the necessary data and set POINTER-kind Hybrid objects for clusters and data
  ///        of the stated sizes.
  void allocate();

  /// \brief Validate the index of the requested cluster.
  ///
  /// \param cluster_index  The cluster of interest
  void validateClusterIndex(int cluster_index) const;
};

} // namespace stmath
} // namespace stormm

#include "cluster_manager.tpp"

#endif
