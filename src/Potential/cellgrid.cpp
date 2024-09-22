#include "copyright.h"
#include "cellgrid.h"

namespace stormm {
namespace energy {

using stmath::mean;

//-------------------------------------------------------------------------------------------------
CellOriginsWriter::CellOriginsWriter(const int stride_in, llint* ax_in, int* ax_ovrf_in,
                                     llint* bx_in, int* bx_ovrf_in, llint* by_in, int* by_ovrf_in,
                                     llint* cx_in, int* cx_ovrf_in, llint* cy_in, int* cy_ovrf_in,
                                     llint* cz_in, int* cz_ovrf_in) :
  stride{stride_in}, ax{ax_in}, ax_ovrf{ax_ovrf_in}, bx{bx_in}, bx_ovrf{bx_ovrf_in}, by{by_in},
  by_ovrf{by_ovrf_in}, cx{cx_in}, cx_ovrf{cx_ovrf_in}, cy{cy_in}, cy_ovrf{cy_ovrf_in}, cz{cz_in},
  cz_ovrf{cz_ovrf_in}
{}

//-------------------------------------------------------------------------------------------------
CellOriginsReader::CellOriginsReader(const int stride_in, const llint* ax_in,
                                     const int* ax_ovrf_in, const llint* bx_in,
                                     const int* bx_ovrf_in, const llint* by_in,
                                     const int* by_ovrf_in, const llint* cx_in,
                                     const int* cx_ovrf_in, const llint* cy_in,
                                     const int* cy_ovrf_in, const llint* cz_in,
                                     const int* cz_ovrf_in) :
  stride{stride_in}, ax{ax_in}, ax_ovrf{ax_ovrf_in}, bx{bx_in}, bx_ovrf{bx_ovrf_in}, by{by_in},
  by_ovrf{by_ovrf_in}, cx{cx_in}, cx_ovrf{cx_ovrf_in}, cy{cy_in}, cy_ovrf{cy_ovrf_in}, cz{cz_in},
  cz_ovrf{cz_ovrf_in}
{}

//-------------------------------------------------------------------------------------------------
CellOriginsReader::CellOriginsReader(const CellOriginsWriter &w) :
  stride{w.stride}, ax{w.ax}, ax_ovrf{w.ax_ovrf}, bx{w.bx}, bx_ovrf{w.bx_ovrf}, by{w.by},
  by_ovrf{w.by_ovrf}, cx{w.cx}, cx_ovrf{w.cx_ovrf}, cy{w.cy}, cy_ovrf{w.cy_ovrf}, cz{w.cz},
  cz_ovrf{w.cz_ovrf}
{}

//-------------------------------------------------------------------------------------------------
CellOriginsReader::CellOriginsReader(const CellOriginsWriter *w) :
  stride{w->stride}, ax{w->ax}, ax_ovrf{w->ax_ovrf}, bx{w->bx}, bx_ovrf{w->bx_ovrf}, by{w->by},
  by_ovrf{w->by_ovrf}, cx{w->cx}, cx_ovrf{w->cx_ovrf}, cy{w->cy}, cy_ovrf{w->cy_ovrf}, cz{w->cz},
  cz_ovrf{w->cz_ovrf}
{}

//-------------------------------------------------------------------------------------------------
double computeMigrationRate(const double effective_cutoff, const double sigma) {

  // Compute a grid plotting the Gaussian propensity for movement.
  const double invs2 = 1.0 / (3.0 * sigma * sigma);
  const int ngrid = 12;
  const int ngrid3 = ngrid * ngrid * ngrid;
  std::vector<double> density(ngrid3);
  double gss = 0.0;
  for (int i = 0; i < ngrid; i++) {
    const double di = (0.3125 + (static_cast<double>(i) * 0.625)) * sigma;
    for (int j = 0; j < ngrid; j++) {
      const double dj = (0.3125 + (static_cast<double>(j) * 0.625)) * sigma;
      for (int k = 0; k < ngrid; k++) {
        const double dk = (0.3125 + (static_cast<double>(k) * 0.625)) * sigma;
        const double r2 = (di * di) + (dj * dj) + (dk * dk);
        const int cubelet_idx = (((k * ngrid) + j) * ngrid) + i;
        density[cubelet_idx] = exp(-0.5 * r2 * invs2);
        gss += density[cubelet_idx];
      }
    }
  }
  gss = 0.125 / gss;
  for (int i = 0; i < ngrid3; i++) {
    density[i] *= gss;
  }

  // The probability of an atom diffusing out of the cell from some point in the cell can now be
  // computed as the integral of the integral of the probability density of movement which would
  // exit the cell times the volume element.  Foremost, there are six slabs within the cell, which
  // here is approximated as a cube of the effective cutoff in width, from which particles could
  // diffuse outward through one of the faces.
  double one_way_out = 0.0;
  for (int k = 0; k < ngrid; k++) {
    const double kfac = static_cast<double>(4 * (k + 1));
    double ij_sum = 0.0;
    for (int j = 0; j < ngrid; j++) {
      for (int i = 0; i < ngrid; i++) {
        ij_sum += density[(((k * ngrid) + j) * ngrid) + i];
      }
    }
    one_way_out += kfac * ij_sum;
  }
  const double d_ngrid = ngrid;
  const double depth = d_ngrid * sigma;
  const double slab_width = effective_cutoff - (2.0 * depth);
  const double slab_volume = 6.0 * (slab_width * slab_width * depth);
  const double slab_propensity = (one_way_out / d_ngrid) * slab_volume;

  // The probability of an atom near one of the edges can be computed based on the original
  // density map as well.
  double two_ways_out = 0.0;
  for (int k = 0; k < ngrid; k++) {
    for (int j = 0; j <= k; j++) {
      const int jk_cross_sect = (2 * ngrid * (j + 1)) - ((j + 1) * (j + 2) / 2);
      const double jk_fac = 3.0 * static_cast<double>((jk_cross_sect) * (1 + (j != k)));
      double i_sum = 0.0;
      const int jk_idx = ((k * ngrid) + j) * ngrid;
      for (int i = 0; i < ngrid; i++) {
        i_sum += density[jk_idx + i];
      }
      two_ways_out += jk_fac * i_sum;
    }
  }
  const double pencil_volume = 12.0 * slab_width * depth * depth;
  const double edge_propensity = (two_ways_out / (d_ngrid * d_ngrid)) * pencil_volume;
  
  // The probability of an atom near one of the corners can finally be computed.
  double three_ways_out = 0.0;
  for (int k = 0; k < ngrid; k++) {
    for (int j = 0; j <= k; j++) {
      for (int i = 0; i <= j; i++) {
        const int ijk_cross_sect = (3 * ngrid * ngrid * (i + 1)) - ((i + 1) * (i + 1) * (i + 1)) -
                                   (3 * (i + 1) * (i + 1) * (ngrid - (i + 1)));
        int cs_mult;
        if (i == j && j == k) {
          cs_mult = 1;
        }
        else if (i == j || j == k || i == k) {
          cs_mult = 3;
        }
        else {
          cs_mult = 6;
        }
        const double ijk_fac = 7.0 * static_cast<double>(ijk_cross_sect * cs_mult);
        three_ways_out += ijk_fac * density[(((k * ngrid) + j) * ngrid) + i];
      }
    }
  }
  const double corner_volume = 8.0 * depth * depth * depth;
  const double corner_propensity = three_ways_out / (d_ngrid * d_ngrid * d_ngrid) * corner_volume;
  
  // The probability of a particle exiting the cell upon movement is the sum of propensities from
  // one of the six slab regions, one of the twelve edge regions, or one of the eight corner
  // regions, all divided by the total volume of the original cell.
  const double cell_volume = effective_cutoff * effective_cutoff * effective_cutoff;
  return (slab_propensity + edge_propensity + corner_propensity) / cell_volume;
}

//-------------------------------------------------------------------------------------------------
int3 optimizeCellConfiguration(const int cell_na, const int cell_nb, const int cell_nc,
                               const int subdivisions) {
  const std::vector<uint> prime_factors = { 2, 3, 5, 7, 11 };
  
  // While there is no restriction on the FFT dimension in any particular direction, it is
  // important to try and avoid combining "bad" radices in the FFT: no more than one factor of 11
  // or 7, and avoid combinations of 7 and 11.  If the subdivision is alreayd 7 or 11, take that
  // into account.
  const ullint big_product = (subdivisions == 7 || subdivisions == 11) ?
                             ipowl(2, 14) * ipowl(3, 10) * ipowl(5, 6) :
                             ipowl(2, 14) * ipowl(3, 10) * ipowl(5, 6) * 7LL * 11LL;
  int3 result;
  result.x = nearestFactor(big_product, cell_na, prime_factors, LimitApproach::BELOW);
  result.y = nearestFactor(big_product, cell_nb, prime_factors, LimitApproach::BELOW);
  result.z = nearestFactor(big_product, cell_nc, prime_factors, LimitApproach::BELOW);
  if ((result.x % 77) == 0) {
    result.x = (result.x / 77) * 75;
  }
  if ((result.y % 77) == 0) {
    result.y = (result.y / 77) * 75;
  }
  if ((result.z % 77) == 0) {
    result.z = (result.z / 77) * 75;
  }
  return result;
}

//-------------------------------------------------------------------------------------------------
void loadRuler(const int cell_count, const int95_t cell_len, const int95_t box_len,
               const int remainder, llint* ruler, int* ruler_ovrf) {
  int95_t running_sum = { 0LL, 0 };
  for (int i = 0; i < cell_count; i++) {
    const int i_bump = (remainder * i) / cell_count;
    const int95_t adj_tick = hostSplitFPSum(running_sum, i_bump, 0);
    ruler[i] = adj_tick.x;
    ruler_ovrf[i] = adj_tick.y;
    running_sum = hostSplitFPSum(running_sum, cell_len);
  }
  ruler[cell_count] = box_len.x;
  ruler_ovrf[cell_count] = box_len.y;
}

} // namespace energy
} // namespace stormm
