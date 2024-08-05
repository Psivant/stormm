// -*-c++-*-
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

using namelist::default_exchange_probability;
using parse::NumberFormat;
using parse::realToString;
using topology::AtomGraph;

// Constants
const double a0 = -59.2194;
const double a1 = 0.07594;
const double b0 = -22.8396;
const double b1 = 0.01347;
const double d0 = 1.1677;
const double d1 = 0.002976;
const int max_iter = 100;
const double kb = 0.008314;
const double pi = 3.14159265358979;

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
  return du * integral / (s12 * std::sqrt(2 * pi));
}

//-------------------------------------------------------------------------------------------------
std::vector<double> vanDerSpoel(double p_des, double t_low, double t_high, int n_w, int n_p,
																double tol, int pc, int wc, int hff, int vs, int alg,
																ExceptionResponse policy) {
	std::vector<double> t;
	std::vector<double> p;
	std::vector<double> siigma;
	std::vector<double> muu;
	std::vector<double> mm;
	std::vector<double> ss;

	int n_pp = 0;
	int n_prot = 0;
	int nh = 0;
	int vc = 0;
	int nc = 0;

  // Input variables error checking
  if(!std::isfinite(p_des) || !std::isfinite(t_low) || !std::isfinite(t_high) || 
          !std::isfinite(n_p) || !std::isfinite(n_w) || !std::isfinite(tol)) {
      switch(policy){
      case ExceptionResponse::DIE:
        rtErr("Some of the input variables are not numbers", "vanDerSpoel");
        break;
      case ExceptionResponse::WARN:
        rtWarn("Some of your input variables are not numbers", "vanDerSpoel");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
  }
  if((p_des > 1) || (p_des < 0)) {
        switch(policy){
        case ExceptionResponse::DIE:
          rtErr("You have to give a probability between 0 and 1.", "vanDerSpoel");
          break;
        case ExceptionResponse::WARN:
          rtWarn("You have to give a probability between 0 and 1.", "vanDerSpoel");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
  }
  if ((t_low <= 0) || (t_high <= 0)){
      switch(policy) {
      case ExceptionResponse::DIE:
        rtErr("The temperatures have to be above 0.", "vanDerSpoel");
        break;
      case ExceptionResponse::WARN:
        rtWarn("The temperatures have to be above 0.", "vanDerSpoel");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
  }
  if(n_p == 0){
      switch(policy) {
      case ExceptionResponse::DIE:
        rtErr("Cannot have 0 protein atoms.", "vanDerSpoel");
        break;
      case ExceptionResponse::WARN:
        rtWarn("Cannot have 0 protein atoms.", "vanDerSpoel");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
  }
  if(alg != 0) {
      switch(policy) {
      case ExceptionResponse::DIE:
        rtErr("Cannot do constant volume yet.", "vanDerSpoel");
        break;
      case ExceptionResponse::WARN:
        rtWarn("Cannot do constant volume yet.", "vanDerSpoel");
        break;
      case ExceptionResponse::SILENT:
        break;
      }
  }
  //-----------------------------------------------------------------------------------------------
	if (hff == 0) {
			nh = round(n_p * 0.5134);
			if (vs == 1) {
					vc = round(1.91 * nh);
			}
			n_prot = n_p;
	} else {
			n_pp = round(n_p / 0.65957);
			nh = round(n_p * 0.22);
			if (vs == 1) {
					vc = round(n_p + 1.91 * nh);
			}
			n_prot = n_pp;
	}

	if (pc == 1) {
			nc = nh;
	} else if (pc == 2) {
			nc = n_p;
	}

	int Ndf = (9 - wc) * n_w + 3 * n_p - nc - vc;
	double FlexEner = 0.5 * kb * (nc + vc + wc * n_w);

	int index = 0;
	t.push_back(t_low);

	while (t[index] < t_high) {
			double piter = 0.0;
			int forward = 1;
			int iter = 0;
			double T1 = t[index];
			double T2 = T1 + 1;
			if (T2 >= t_high) {
					T2 = t_high;
			}
			double low = T1;
			double high = t_high;

			while (abs(p_des - piter) > tol && iter < max_iter) {
					iter++;
					double mu12 = (T2 - T1) * ((a1 * n_w) + (b1 * n_prot) - FlexEner);
					mm.push_back(mu12);

					double cc = (1 / kb) * ((1 / T1) - (1 / T2));
					double Delta = cc * mu12;

					double var = Ndf * (d1 * d1 * (T1 * T1 + T2 * T2) +
							2 * d1 * d0 * (T1 + T2) +
							2 * d0 * d0);

					double sig12 = sqrt(var);
					ss.push_back(sig12);

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

					double erfarg1 = mu12 / (sig12 * sqrt(2));
					double I1 = 0.5 * (std::erfc(erfarg1));

					double I2 = myintegral(mu12, sig12, cc);
					piter = (I1 + I2);

					if (piter > p_des) {
							if (forward == 1) {
									T2 = T2 + 1.0;
							} else if (forward == 0) {
									low = T2;
									T2 = low + ((high - low) / 2);
							}
							if (T2 >= t_high) {
									T2 = t_high;
							}
					} else if (piter < p_des) {
							if (forward == 1) {
									forward = 0;
									low = T2 - 1.0;
							}
							high = T2;
							T2 = low + ((high - low) / 2);
					}
			}

			p.push_back(piter);
			siigma.push_back(sqrt(Ndf) * (d0 + d1 * T1));
			muu.push_back(calc_mu(n_w, n_prot, T1, FlexEner));

			index++;
			t.push_back(T2);
	}

	siigma.push_back(sqrt(Ndf) * (d0 + d1 * t[index]));
	muu.push_back(calc_mu(n_w, n_prot, t[index], FlexEner));

	return t;
	}

}	// namespace sampling
}	// namespace stormm
