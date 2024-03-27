// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace stmath {

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
Tcalc clusteringDistance(const Tdata* data_points, const size_t pt, const Tdata centroid,
                         const Tcalc scale_factor, const BoundaryCondition boundaries,
                         const Tcalc range) {
  Tcalc dv = centroid - data_points[pt];
  switch (boundaries) {
  case BoundaryCondition::PERIODIC:
    dv = imageValue(dv, range, ImagingMethod::MINIMUM_IMAGE, scale_factor);
    break;
  case BoundaryCondition::ISOLATED:
    break;
  }

  // The distance is returned in the unit system of the data itself.
  return dv;
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
void seedKMeans(ClusterManager<Tdata, Tcalc> *cls, const Tdata* data_points, const size_t count,
                const int clusters, Xoshiro256ppGenerator *xrs, const Tcalc scale_factor,
                const BoundaryCondition boundaries, const Tcalc range) {

  // Prepare the cluster manager to accept this data set.
  cls->resize(count, clusters);
  ClusterManagerKit<Tcalc> clskit = cls->data();
  
  // Choose one centroid from a random, uniform distribution of the data points.
  const size_t first_choice = xrs->uniformRandomNumber() * static_cast<double>(count);
  cls->setCentroid(0, data_points[first_choice]);
  clskit.cls_primes[0] = first_choice;
  for (int i = 1; i <= clusters; i++) {

    // Find the nearest cluster for each point.  Catalog the results.  Track the maximum distance
    // that any point lies from the closest cluster centroid.
    Tcalc max_distance = 0.0;
    size_t max_loc = first_choice;
    for (size_t j = 0; j < count; j++) {
      Tcalc min_distance = clusteringDistance(data_points, j, cls->getCentroid(0), scale_factor,
                                              boundaries, range);
      int min_loc = 0;
      for (int k = 1; k < i; k++) {
        const Tcalc trial_distance = clusteringDistance(data_points, j, cls->getCentroid(k),
                                                        scale_factor, boundaries, range);
        if (trial_distance < min_distance) {
          min_distance = trial_distance;
          min_loc = k;
        }
      }
      clskit.cls_homes[j] = min_loc;
      clskit.centroid_dist[j] = min_distance;
      if (min_distance > max_distance) {
        max_distance = min_distance;
        max_loc = j;
      }
    }

    // The point with the maximum minimum distance to any cluster centroid is chosen as the
    // location of the next cluster centroid.  The final loop sets all data points into their
    // appropriate seed clusters, but will not create another new cluster. 
    if (i < clusters) {
      cls->setCentroid(i, data_points[max_loc]);
      clskit.cls_primes[i] = max_loc;
    }
  }
  
  // Compute the home clusters of each data point.
  indexingArray(clskit.cls_homes, clskit.cls_membership, clskit.cls_bounds, count, clusters);

  // Calculate the mean values and update each centroid.  The centroids are calculated in the
  // calculation type and then recast back to the data set's native type.  This is critical as
  // fixed-precision data sets may not be suitable for computing a mean, which would involve a
  // running sum that could exceed the fixed precision model bounds even if each member does not.
  std::vector<Tcalc> next_centroids(clusters, 0.0);
  for (size_t i = 0; i < count; i++) {
    next_centroids[clskit.cls_homes[i]] += data_points[i];
  }
  for (int i = 0; i < clusters; i++) {
    cls->setCentroid(i, next_centroids[i] / static_cast<Tcalc>(clskit.cls_bounds[i + 1] -
                                                               clskit.cls_bounds[i]));
  }
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
void kMeans(ClusterManager<Tdata, Tcalc> *cls, const Tdata* data_points, const size_t count,
            const int clusters, const Tcalc scale_factor, const BoundaryCondition boundaries,
            const Tcalc range) {
  const bool tcalc_is_double = (std::type_index(typeid(Tcalc)).hash_code() == double_type_index);
  if (clusters > count) {
    rtErr("The number of clusters (" + std::to_string(clusters) + ") cannot exceed the number of "
          "data points (" + std::to_string(count) + ").", "kMeans");
  }

  // For one-dimensional data, whether in two or three dimensions, the problem is trivial in
  // comparison to the multi-dimensional case: every dominance relationship chains.  Sort the data
  // points, after re-imaging if periodic boundary conditions exist, to expedite partitioning.
  struct TSort {
    Tdata dpt;
    size_t loc;
  };
  std::vector<TSort> tagged_data(count);
  const Tcalc scaled_range = range * scale_factor;
  for (size_t i = 0; i < count; i++) {
    const Tcalc reim_dp = imageValue(data_points[i], range, ImagingMethod::PRIMARY_UNIT_CELL,
                                     scale_factor);
    tagged_data[i].dpt = reim_dp;
    tagged_data[i].loc = i;
  }
  std::sort(tagged_data.begin(), tagged_data.end(),
            [](TSort x, TSort y) { return x.dpt < y.dpt; });
  std::vector<size_t> partitions(clusters);
  const double dpt_count = count;
  size_t part_loc = 0;
  const double portion = 1.0 / static_cast<double>(clusters);
  for (int i = 0; i < clusters; i++) {
    partitions[i] = part_loc;
    part_loc = std::max(part_loc + 1,
                        static_cast<size_t>(static_cast<double>(i + 1) * portion * dpt_count));
  }
  const Tcalc zero = 0.0;
  std::vector<Tdata> tmp_centroids(clusters);
  bool unfinished;
  do {
    unfinished = false;
    switch (boundaries) {
    case BoundaryCondition::ISOLATED:

      // Compute the centroids of each cluster.
      for (int i = 0; i < clusters; i++) {
        const size_t llim = partitions[i];
        const size_t hlim = (i + 1 == clusters) ? count : partitions[i + 1];
        Tcalc centsum = zero;
        size_t pos = llim;
        for (size_t pos = llim; pos < hlim; pos++) {
          centsum += static_cast<Tcalc>(tagged_data[pos].dpt);        
        }
        tmp_centroids[i] = centsum / static_cast<Tcalc>(hlim - llim);
      }

      // Calculate the new best positions for the partitions.
      for (int i = 1; i < clusters; i++) {
        size_t pos = partitions[i];
        const Tdata low_centroid = tmp_centroids[i - 1];
        const Tdata cur_centroid = tmp_centroids[i];

        // While the partition between clusters i - 1 and i is closer to the centroid of cluster
        // i - 1, shift it to the right.
        while (pos < count && tagged_data[pos].dpt > low_centroid &&
               (i == clusters - 1 || pos < partitions[i + 1]) &&
               tagged_data[pos].dpt - low_centroid < cur_centroid - tagged_data[pos].dpt) {
          pos++;
        }

        // The process above would be blind to a partition that begins too far to the right.
        bool move_point;
        do {
          move_point = false;
          if (pos > partitions[i - 1]) {
            const size_t pos_m = pos - 1;
            if (tagged_data[pos_m].dpt > low_centroid &&
                cur_centroid - tagged_data[pos_m].dpt < tagged_data[pos_m].dpt - low_centroid) {
              pos--;
              move_point = true;
            }
          }
        } while (move_point);
        unfinished = (unfinished || partitions[i] != pos);
        partitions[i] = pos;
      }
      break;
    case BoundaryCondition::PERIODIC:
      {
        // Compute the centroids of each cluster.
        Tcalc offset = zero;
        for (int i = 0; i < clusters; i++) {
          const size_t llim = partitions[i];
          const size_t hlim = (i + 1 == clusters) ? partitions[0] : partitions[i + 1];
          Tcalc centsum = zero;
          size_t pos = llim;
          while (pos != hlim) {
            centsum += static_cast<Tcalc>(tagged_data[pos].dpt + offset);
            pos++;

            // If the end is reached, start back from the beginning but offset all values to
            // keep them all imaged with respect to one another.
            if (pos == count) {
              pos = 0;
              offset = range;
            }
          }
          if (hlim > llim) {
            tmp_centroids[i] = centsum / static_cast<Tcalc>(hlim - llim);
          }
          else {
            tmp_centroids[i] = centsum / static_cast<Tcalc>(hlim + (count - llim));
          }
        }

        // Calculate the new best positions for the partitions.
        for (int i = 0; i < clusters; i++) {
          size_t pos = partitions[i];
          const Tdata low_centroid = (i == 0) ? tmp_centroids[clusters - 1] : tmp_centroids[i - 1];
          const Tdata cur_centroid = tmp_centroids[i];

          // While the point at the partition between clusters i - 1 and i is closer to the
          // centroid of cluster i - 1, shift the partition to the left.
          bool move_point;
          do {
            move_point = false;

            // Check that the partitions do not collide
            const int next_partition = (i == clusters - 1) ? 0 : i + 1;
            const size_t ilimit = (partitions[next_partition] - 1 > count) ?
                                  count - 1 : partitions[next_partition] - 1;
            if (pos == ilimit) {
              continue;
            }

            // Compute distances for the point to the nearest centroids and check whether the
            // partition can be shifted to the right.
            const Tcalc d_to_low = imageValue<Tdata, Tcalc>(tagged_data[pos].dpt - low_centroid,
                                                            range, ImagingMethod::MINIMUM_IMAGE,
                                                            scale_factor);
            const Tcalc d_to_cur = imageValue<Tdata, Tcalc>(tagged_data[pos].dpt - cur_centroid,
                                                            range, ImagingMethod::MINIMUM_IMAGE,
                                                            scale_factor);
            const Tcalc abs_to_low = std::abs(d_to_low);
            const Tcalc abs_to_cur = std::abs(d_to_cur);

            // If the value at the partition (representing the low end of the ith cluster) is
            // closer to the (i - 1)th cluster's centroid than to the ith cluster's centroid,
            // advance the partition location by one to include the partition's current location
            // in the (i - 1)th cluster.
            if (abs_to_low < abs_to_cur) {
              pos = (pos == count - 1) ? 0 : pos + 1;
              move_point = true;
            }
          } while (move_point);

          // The above analysis might be blind to a situation in which the partition is closer to
          // the ith cluster's centroid than the (i - 1)th cluster's centroid, but the preceding
          // point would also be closer to the ith cluster's centroid, and thus should not be part
          // of the (i - 1)th cluster.
          do {
            move_point = false;
            const size_t pos_minus_one = (pos == 0) ? count - 1 : pos - 1;

            // Check that the partitions do not collide.
            const int next_partition = (i == 0) ? clusters - 1 : i - 1;
            const size_t ilimit = (partitions[next_partition] + 1 >= count) ?
                                  0 : partitions[next_partition] + 1;
            if (pos_minus_one == ilimit) {
              continue;
            }

            // Check whether the partition can be shifted to the left.
            const Tcalc dm_to_low = imageValue<Tdata, Tcalc>(tagged_data[pos_minus_one].dpt -
                                                             low_centroid, range,
                                                             ImagingMethod::MINIMUM_IMAGE,
                                                             scale_factor);
            const Tcalc dm_to_cur = imageValue<Tdata, Tcalc>(tagged_data[pos_minus_one].dpt -
                                                             cur_centroid, range,
                                                             ImagingMethod::MINIMUM_IMAGE,
                                                             scale_factor);
            const Tcalc mabs_to_low = std::abs(dm_to_low);
            const Tcalc mabs_to_cur = std::abs(dm_to_cur);
            if (mabs_to_cur < mabs_to_low) {
              pos = pos_minus_one;
              move_point = true;
            }
          } while (move_point);
          unfinished = (unfinished || partitions[i] != pos);
          partitions[i] = pos;
        }
      }
      break;
    }
  } while (unfinished);

  // The partitions are now set.  Trace back to the original data point indices and compose each
  // cluster.  Accumulate the results in the input cluster manager.
  ClusterManagerKit clskit = cls->data();
  size_t cls_bound_tracker = 0;
  for (int i = 0; i < clusters; i++) {
    cls->setCentroid(i, tmp_centroids[i]);
    const size_t llim = partitions[i];
    const size_t hlim = (i == clskit.clusters - 1) ? partitions[0] : partitions[i + 1];
    clskit.cls_bounds[i] = cls_bound_tracker;
    size_t pos = llim;
    Tcalc prime_centroid_proximity;
    Tcalc tmp_rmsd = zero;
    size_t prime_pos;
    while (pos != hlim) {
      const size_t orig_dpt_idx = tagged_data[pos].loc;
      clskit.cls_membership[cls_bound_tracker] = orig_dpt_idx;
      clskit.cls_homes[orig_dpt_idx] = i;
      clskit.centroid_dist[orig_dpt_idx] = std::abs(data_points[orig_dpt_idx] - tmp_centroids[i]);
      if (pos == llim || clskit.centroid_dist[orig_dpt_idx] < prime_centroid_proximity) {
        prime_centroid_proximity = clskit.centroid_dist[orig_dpt_idx];
        prime_pos = pos;
      }
      tmp_rmsd += clskit.centroid_dist[orig_dpt_idx] * clskit.centroid_dist[orig_dpt_idx];
      cls_bound_tracker++;
      pos++;
      pos = (pos < count) * pos;
    }
    clskit.cls_primes[i] = tagged_data[prime_pos].loc;
    tmp_rmsd /= static_cast<Tcalc>(cls_bound_tracker - clskit.cls_bounds[i]);
    clskit.cls_rmsd[i] = (tcalc_is_double) ? sqrt(tmp_rmsd) : sqrtf(tmp_rmsd);
  }
  clskit.cls_bounds[clusters] = cls_bound_tracker;
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
void kMeans(ClusterManager<Tdata, Tcalc> *cls, const std::vector<Tdata> &data_points,
            const int clusters, const Tcalc scale_factor, const BoundaryCondition boundaries,
            const Tcalc range) {
  kMeans(cls, data_points.data(), data_points.size(), clusters, scale_factor, boundaries, range);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
void kMeans(ClusterManager<Tdata, Tcalc> *cls, const Hybrid<Tdata> &data_points,
            const int clusters, const Tcalc scale_factor, const BoundaryCondition boundaries,
            const Tcalc range) {
  kMeans(cls, data_points.data(), data_points.size(), clusters, scale_factor, boundaries, range);
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
ClusterManager<Tdata, Tcalc> kMeans(const Tdata* data_points, const size_t count,
                                    const int clusters, const Tcalc scale_factor,
                                    const BoundaryCondition boundaries, const Tcalc range) {
  ClusterManager result(data_points->size(), clusters);
  kMeans(&result, data_points, count, clusters, scale_factor, boundaries, range);
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
ClusterManager<Tdata, Tcalc> kMeans(const std::vector<Tdata> &data_points, const int clusters,
                                    const Tcalc scale_factor, const BoundaryCondition boundaries,
                                    const Tcalc range) {
  ClusterManager result(data_points.size(), clusters);
  kMeans(&result, data_points.data(), data_points.size(), clusters, scale_factor, boundaries,
         range);
  return result;
}

//-------------------------------------------------------------------------------------------------
template <typename Tdata, typename Tcalc>
ClusterManager<Tdata, Tcalc> kMeans(const Hybrid<Tdata> &data_points, const int clusters,
                                    const Tcalc scale_factor, const BoundaryCondition boundaries,
                                    const Tcalc range) {
  ClusterManager result(data_points.size(), clusters);
  kMeans(&result, data_points.data(), data_points.size(), clusters, scale_factor, boundaries,
         range);
  return result;
}

} // namespace stmath
} // namespace stormm
