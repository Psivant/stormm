#include <cmath>
#include <vector>
#include "copyright.h"
#include "Namelists/nml_remd.h"
#include "Parsing/parsing_enumerators.h"
#include "Parsing/parse.h"
#include "Reporting/error_format.h"
#include "Topology/atomgraph.h"
#include "Topology/atomgraph_abstracts.h"
#include "temp_distributions.h"

namespace stormm {
namespace sampling {

using constants::PrecisionModel;
using constants::ExceptionResponse;
using symbols::boltzmann_constant;
using topology::AtomGraph;

const double a0 = -59.2194;
const double a1 = 0.07594;
const double b0 = -22.8396;
const double b1 = 0.01347;
const double d0 = 1.1677;
const double d1 = 0.002976;
const int max_iter = 100;

//-------------------------------------------------------------------------------------------------
double calc_mu(int n_w, int n_p, double temp, double f_ener) {
  return ((a0 + a1 * temp) * n_w + (b0 + b1 * temp) * n_p - temp * f_ener);
}

//-------------------------------------------------------------------------------------------------
double myeval(double m12, double s12, double cc, double u) {
  double argument = -cc * u - (u - m12) * (u - m12) / (2 * s12 * s12);
  return std::exp(argument);
}

//-------------------------------------------------------------------------------------------------
double myintegral(double m12, double s12, double cc) {
  double integral = 0.0;
  double umax = m12 + 5 * s12;
  double du = umax / 100;
  for (double u = 0; u < umax; u += du) {
    double di = myeval(m12, s12, cc, u + du / 2);
    integral += di;
  }
  const double pi = 3.14159265358979;
  return du * integral / (s12 * std::sqrt(2 * pi));
}

//-------------------------------------------------------------------------------------------------
std::vector<double> vanDerSpoel(double p_des, double t_low, double t_high, int n_w, int n_p,
                                double tol, int pc, int wc, int hff, int vs, int alg,
                                ExceptionResponse policy) {
  std::vector<double> temperatures;

  double kB = 0.008314;
  int n_pp = 0;
  int n_prot = 0;
  int nh = 0;
  int nc = 0;
  int vc = 0;
  double flex_ener = 0.5 * kB * (nc + vc + wc * n_w);

  if (hff == 0) {
    nh = static_cast<int>(std::round(n_p * 0.5134));
    if (vs == 1) {
      vc = static_cast<int>(std::round(1.91 * nh));
    }
    n_prot = n_p;
  } else {
    n_pp = static_cast<int>(std::round(n_p / 0.65957));
    nh = static_cast<int>(std::round(n_p * 0.22));
    if (vs == 1) {
      vc = static_cast<int>(std::round(n_p + 1.91 * nh));
    }
    n_prot = n_pp;
  }

  if (pc == 1) {
    nc = nh;
  } else if (pc == 2) {
    nc = n_p;
  }

  int ndf = (9 - wc) * n_w + 3 * n_p - nc - vc;

  // Iterate through temperatures
  int index = 1;
  double T = t_low;
  while (T < t_high) {
    double piter = 0;
    bool forward = true;
    int iter = 0;
    double T1 = T;
    double T2 = T1 + 1;
    if (T2 >= t_high) {
      T2 = t_high;
    }
    double low = T1;
    double high = t_high;

    while (std::abs(p_des - piter) > tol && iter < max_iter) {
      iter++;
      double mu12 = (T2 - T1) * ((a1 * n_w) + (b1 * n_prot) - flex_ener);
      double cc = (1 / kB) * ((1 / T1) - (1 / T2));
      double delta = cc * mu12;
      double var = ndf * (d1 * d1 * (T1 * T1 + T2 * T2) +
                          2 * d1 * d0 * (T1 + T2) +
                          2 * d0 * d0);

      double sig12 = std::sqrt(var);

      if (sig12 == 0) {
        switch(policy) {
        case ExceptionResponse::DIE:
          rtErr("Sigma = 0.", "vanDerSpoel");
          break;
        case ExceptionResponse::WARN:
          rtWarn("Sigma=0.", "vanDerSpoel");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
      }

      double erfarg1 = mu12 / (sig12 * std::sqrt(2));
      double I1 = 0.5 * std::erfc(erfarg1);

      double I2 = myintegral(mu12, sig12, cc);
      piter = (I1 + I2);

      if (piter > p_des) {
        if (forward) {
          T2 += 1.0;
        } else {
          low = T2;
          T2 = low + ((high - low) / 2);
        }
        if (T2 >= t_high) {
          T2 = t_high;
        }
      } else {
        if (forward) {
          forward = false;
          low = T2 - 1.0;
        }
        high = T2;
        T2 = low + ((high - low) / 2);
      }
    }
    temperatures.push_back(std::floor(T * 100.0) / 100.0);
    index++;
    T = T2;
  }

  temperatures.push_back(std::floor(T * 100.0) / 100.0);

  return temperatures;
}

} // namespace sampling
} // namespace stormm
