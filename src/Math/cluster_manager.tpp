// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
template <typename Tcalc>
ClusterManagerKit<Tcalc>::ClusterManagerKit(const size_t data_points_in, const int clusters_in,
                                            int* cls_homes_in, size_t* cls_membership_in,
                                            size_t* cls_bounds_in, size_t* cls_primes_in,
                                            Tcalc* cls_rmsd_in, Tcalc* centroid_dist_in) :
    data_points{data_points_in}, clusters{clusters_in}, cls_homes{cls_homes_in},
    cls_membership{cls_membership_in}, cls_bounds{cls_bounds_in}, cls_primes{cls_primes_in},
    cls_rmsd{cls_rmsd_in}, centroid_dist{centroid_dist_in}
{}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
ClusterManager<Tdata, Tcalc>::ClusterManager(const size_t data_point_count_in,
                                             const int cluster_count_in,
                                             const int attempt_count_in) :
    data_point_count{data_point_count_in}, cluster_count{cluster_count_in},
    attempt_count{attempt_count_in},
    cluster_homes{data_point_count_in, "clsrep_homes"},
    cluster_membership{HybridKind::POINTER, "clsrep_membership"},
    cluster_bounds{HybridKind::POINTER, "clsrep_cluster_bounds"},
    cluster_primes{HybridKind::POINTER, "clsrep_cluster_primes"},
    cluster_rmsd{HybridKind::POINTER, "clsrep_cluster_rmsd"},
    centroid_distances{HybridKind::POINTER, "clsrep_distances"},
    sizet_data{HybridKind::ARRAY, "clsrep_sizet_data"},
    real_data{HybridKind::ARRAY, "clsrep_real_data"},
    centroids{}
{
  allocate();
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
ClusterManager<Tdata, Tcalc>::ClusterManager(const ClusterManager<Tdata, Tcalc> &original) :
    data_point_count{original.data_point_count},
    cluster_count{original.cluster_count},
    attempt_count{original.attempt_count},
    cluster_homes{original.cluster_homes},
    cluster_membership{original.cluster_membership},
    cluster_bounds{original.cluster_bounds},
    cluster_primes{original.cluster_primes},
    cluster_rmsd{original.cluster_rmsd},
    centroid_distances{original.centroid_distances},
    sizet_data{original.sizet_data},
    real_data{original.real_data},
    centroids{original.centroids}
{
  // Rebase the POINTER-kind Hybrid objects using the allocator.
  allocate();
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
ClusterManager<Tdata, Tcalc>::ClusterManager(ClusterManager<Tdata, Tcalc> &&original) :
    data_point_count{std::move(original.data_point_count)},
    cluster_count{std::move(original.cluster_count)},
    attempt_count{std::move(original.attempt_count)},
    cluster_homes{std::move(original.cluster_homes)},
    cluster_membership{std::move(original.cluster_membership)},
    cluster_bounds{std::move(original.cluster_bounds)},
    cluster_primes{std::move(original.cluster_primes)},
    cluster_rmsd{std::move(original.cluster_rmsd)},
    centroid_distances{std::move(original.centroid_distances)},
    sizet_data{std::move(original.sizet_data)},
    real_data{std::move(original.real_data)},
    centroids{std::move(original.centroids)}
{}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
ClusterManager<Tdata, Tcalc>&
ClusterManager<Tdata, Tcalc>::operator=(const ClusterManager<Tdata, Tcalc> &other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }
  data_point_count = other.data_point_count;
  cluster_count = other.cluster_count;
  attempt_count = other.attempt_count;
  cluster_homes = other.cluster_homes;
  cluster_membership = other.cluster_membership;
  cluster_bounds = other.cluster_bounds;
  cluster_primes = other.cluster_primes;
  cluster_rmsd = other.cluster_rmsd;
  centroid_distances = other.centroid_distances;
  sizet_data = other.sizet_data;
  real_data = other.real_data;
  centroids = other.centroids;

  // Rebase the POINTER-kind Hybrid objects using the allocator.
  allocate();
  return *this;
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
ClusterManager<Tdata, Tcalc>&
ClusterManager<Tdata, Tcalc>::operator=(ClusterManager<Tdata, Tcalc> &&other) {

  // Guard against self-assignment
  if (this == &other) {
    return *this;
  }
  data_point_count = other.data_point_count;
  cluster_count = other.cluster_count;
  attempt_count = other.attempt_count;
  cluster_homes = std::move(other.cluster_homes);
  cluster_membership = std::move(other.cluster_membership);
  cluster_bounds = std::move(other.cluster_bounds);
  cluster_primes = std::move(other.cluster_primes);
  cluster_rmsd = std::move(other.cluster_rmsd);
  centroid_distances = std::move(other.centroid_distances);
  sizet_data = std::move(other.sizet_data);
  real_data = std::move(other.real_data);
  centroids = std::move(other.centroids);

  // As with other move assignment functions, pointer repair is not needed.
  return *this;
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
size_t ClusterManager<Tdata, Tcalc>::getDataPointCount() const {
  return data_point_count;
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
int ClusterManager<Tdata, Tcalc>::getClusterCount() const {
  return cluster_count;
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
int ClusterManager<Tdata, Tcalc>::getAttemptCount() const {
  return attempt_count;
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
int ClusterManager<Tdata, Tcalc>::getClusterHome(const size_t point_index) const {
  return cluster_homes.readHost(point_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
size_t ClusterManager<Tdata, Tcalc>::getClusterSize(const int cluster_index) const {
  return cluster_bounds.readHost(cluster_index + 1) - cluster_bounds.readHost(cluster_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
size_t ClusterManager<Tdata, Tcalc>::getClusterPrime(const int cluster_index) const {
  return cluster_primes.readHost(cluster_index);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
std::vector<size_t>
ClusterManager<Tdata, Tcalc>::getClusterMembers(const int cluster_index) const {
  const size_t llim = cluster_bounds.readHost(cluster_index);
  const size_t hlim = cluster_bounds.readHost(cluster_index + 1);
  return cluster_membership.readHost(llim, hlim - llim);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
const Tdata& ClusterManager<Tdata, Tcalc>::getClusterCentroid(const int cluster_index) const {
  return centroids[cluster_index];
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
const std::vector<Tdata>& ClusterManager<Tdata, Tcalc>::getClusterCentroids() const {
  return centroids;
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
const ClusterManagerKit<Tcalc>
ClusterManager<Tdata, Tcalc>::data(const HybridTargetLevel tier) const {
  return ClusterManagerKit<Tcalc>(data_point_count, cluster_count, cluster_homes.data(tier),
                                  cluster_membership.data(tier), cluster_bounds.data(tier),
                                  cluster_primes.data(tier), cluster_rmsd.data(tier),
                                  centroid_distances.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
ClusterManagerKit<Tcalc> ClusterManager<Tdata, Tcalc>::data(const HybridTargetLevel tier) {
  return ClusterManagerKit<Tcalc>(data_point_count, cluster_count, cluster_homes.data(tier),
                                  cluster_membership.data(tier), cluster_bounds.data(tier),
                                  cluster_primes.data(tier), cluster_rmsd.data(tier),
                                  centroid_distances.data(tier));
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
void ClusterManager<Tdata, Tcalc>::resize(const size_t data_point_count_in,
                                          const int cluster_count_in) {
  data_point_count = data_point_count_in;
  cluster_count = cluster_count_in;
  allocate();
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
void ClusterManager<Tdata, Tcalc>::setAttemptCount(const int attempt_count_in) {
  attempt_count = attempt_count_in;
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
void ClusterManager<Tdata, Tcalc>::setCentroid(const int cluster_index, const Tdata value) {
  validateClusterIndex(cluster_index);
  centroids[cluster_index] = value;
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc> void ClusterManager<Tdata, Tcalc>::allocate() {
  cluster_homes.resize(data_point_count);

  // Check whether the templated data type is one of those which serves as its own vector.
  const size_t ct = std::type_index(typeid(Tdata)).hash_code();
  if (ct == std::type_index(typeid(CoordinateSeries<Tdata>)).hash_code() ||
      ct == condensate_type_index ||
      ct == phasespace_synthesis_type_index) {
    centroids.resize(1);
  }
  else {
    centroids.resize(cluster_count);
  }
  const size_t padded_data_point_count = roundUp(data_point_count, warp_size_zu);
  const size_t padded_bounds_count = roundUp<size_t>(cluster_count + 1, warp_size_zu);
  const size_t padded_cluster_count = roundUp<size_t>(cluster_count, warp_size_zu);
  sizet_data.resize(padded_data_point_count + padded_bounds_count + padded_cluster_count);
  real_data.resize(padded_data_point_count + padded_cluster_count);
  cluster_membership.setPointer(&sizet_data, 0, data_point_count);
  cluster_bounds.setPointer(&sizet_data, padded_data_point_count, cluster_count + 1);
  cluster_primes.setPointer(&sizet_data, padded_data_point_count + padded_bounds_count,
                            cluster_count);
  cluster_rmsd.setPointer(&real_data, 0, cluster_count);
  centroid_distances.setPointer(&real_data, padded_cluster_count, data_point_count);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
void ClusterManager<Tdata, Tcalc>::validateClusterIndex(const int cluster_index) const {
  if (cluster_index < 0 || cluster_index >= cluster_count) {
    rtErr("Cluster index " + std::to_string(cluster_index) + " is invalid for a set of " +
          std::to_string(cluster_count) + " clusters.", "ClusterManager", "validateClusterIndex");
  }
}

} // namespace stmath
} // namespace stormm
