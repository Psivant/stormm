// -*-c++-*-
#ifndef STORMM_CLUSTERING_H
#define STORMM_CLUSTERING_H

#include "copyright.h"
#include "Random/random.h"
#include "Structure/local_arrangement.h"
#include "Structure/structure_enumerators.h"
#include "cluster_manager.h"
#include "series_ops.h"

namespace stormm {
namespace stmath {

using random::Xoshiro256ppGenerator;
using structure::BoundaryCondition;
using structure::imageCoordinates;
using structure::imageValue;
using structure::ImagingMethod;

/// \brief Compute the distance between two data points for clustering purposes.  The weighted
///        distance will be computed in all dimensions, clustering according to the unit cell
///        dimensions.
template <typename Tdata, typename Tcalc>
Tcalc clusteringDistance(const Tdata* data_points, size_t pt, Tdata centroid,
                         Tcalc scale_factor = 1.0,
                         BoundaryCondition boundaries = BoundaryCondition::ISOLATED,
                         Tcalc range = 0.0);
  
/// \brief Select initial centroids from points of the original data set for a K-means clustering,
///        using the K-means++ algorithm of Arthur and Vassilvitskii:
///
/// David Arthur and Sergei Vassilvitskii. "K-means++: the advantages of careful seeding." (2007)
/// Proceedings of the eighteenth annual ACM-SIAM symposium on Discrete Algorithms.  Society for
/// Industrial and Applied Mathematics, Philadelphia, PA (United States of America) pp. 1027-1035.
///
/// Overloaded:
///   - Operate on C-style arrays with trusted length
///   - Operate on Standard Template Library vectors
///   - Operate on Hybrid objects
///
/// \param data_points           The series of M data points, arranged in order with N dimensions
///                              as { pt0_1, pt1_1, ..., pt(M-1)_1, pt0_2, pt1_2, ..., pt(M-1)_2,
///                              ..., pt0_(N-1), ..., pt(M-1)_(N-1) }.  The distance between each
///                              point is calculated using the pythagorean theorem in all N
///                              dimensions.
/// \param count                 The number of points, a trusted length if the data is provided in
///                              a C-style array.  The array data_points is expected to be count *
///                              dimensions in length.
/// \param clusters              The number of clusters to sort all data points into.  Converged
///                              clustering will minimize the variance, as judged by mutual
///                              point-to-point distances, within each cluster.
/// \param dimensions            The number of dimensions to the data.  The array data points is
///                              arranged with all points' values along a particular dimension in
///                              contiguous memory.
/// \param weights               Weights assigned to each dimension of a given data point (this
///                              array must have a length equal to the dimensions parameter)
/// \param boundaries            Boundary condition to apply when computing distance between data
///                              points
/// \param umat                  Transformation matrix to take data points into the unit cell space
///                              (only valid if periodic boundary conditions are in effect)
/// \param invu                  Transformation matrix to take data points from the unit cell space
///                              back into real space (only valid if periodic boundary conditions
///                              are in effect)
/// \param unit_cell_dimensions  The number of dimensions in the unit cell space (typically 1 or 3)
/// \param report                The clustering management object, used to report results and also
///                              track the progress of a particular clustering procedure.  This
///                              will be returned with the seed data points copied into its
///                              centroids array.
template <typename Tdata, typename Tcalc>
void seedKMeans(const Tdata* data_points, size_t count, int clusters, Xoshiro256ppGenerator *xrs,
                Tcalc scale_factor = 1.0,
                BoundaryCondition boundaries = BoundaryCondition::ISOLATED, Tcalc range = 0.0);

/// \brief Cluster a series of values into a preset number of groups using the K-means algorithm
///        with k-means++ initialization.
///
/// Overloaded:
///   - Operate on C-style arrays with trusted length
///   - Operate on Standard Template Library vectors
///   - Operate on Hybrid objects
///   - Output the results as a new ClusterReport, or fill in an existing ClusterReport object
///
/// Parameter descriptions follow from the seeding function seedKMeansPP, above.
/// \{
template <typename Tdata, typename Tcalc>
void kMeans(ClusterManager<Tdata, Tcalc> *cls, const Tdata* data_points, size_t count,
            int clusters, Tcalc scale_factor = 1.0,
            BoundaryCondition boundaries = BoundaryCondition::ISOLATED, Tcalc range = 0.0);

template <typename Tdata, typename Tcalc>
void kMeans(ClusterManager<Tdata, Tcalc> *cls, const std::vector<Tdata> &data_points, int clusters,
            Tcalc scale_factor = 1.0, BoundaryCondition boundaries = BoundaryCondition::ISOLATED,
            Tcalc range = 0.0);

template <typename Tdata, typename Tcalc>
void kMeans(ClusterManager<Tdata, Tcalc> *cls, const Hybrid<Tdata> &data_points, int clusters,
            Tcalc scale_factor = 1.0, BoundaryCondition boundaries = BoundaryCondition::ISOLATED,
            Tcalc range = 0.0);

template <typename Tdata, typename Tcalc>
ClusterManager<Tdata, Tcalc> kMeans(const Tdata* data_points, size_t count, int clusters,
                                    Tcalc scale_factor = 1.0,
                                    BoundaryCondition boundaries = BoundaryCondition::ISOLATED,
                                    Tcalc range = 0.0);

template <typename Tdata, typename Tcalc>
ClusterManager<Tdata, Tcalc> kMeans(const std::vector<Tdata> &data_points, int clusters,
                                    Tcalc scale_factor = 1.0,
                                    BoundaryCondition boundaries = BoundaryCondition::ISOLATED,
                                    Tcalc range = 0.0);

template <typename Tdata, typename Tcalc>
ClusterManager<Tdata, Tcalc> kMeans(const Hybrid<Tdata> &data_points, int clusters,
                                    Tcalc scale_factor = 1.0,
                                    BoundaryCondition boundaries = BoundaryCondition::ISOLATED,
                                    Tcalc range = 0.0);
/// \}
  
} // namespace stmath
} // namespace stormm

#include "clustering.tpp"

#endif
