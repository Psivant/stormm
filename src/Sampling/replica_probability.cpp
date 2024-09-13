#include <cmath>
#include "copyright.h"
#include "replica_probability.h"

namespace stormm {
namespace sampling {

//-------------------------------------------------------------------------------------------------
double getTotalHamiltonian(std::vector<double> hamiltonian){
    double totalHamiltonian = 0;
    for (double h : hamiltonian) {
        totalHamiltonian += h;
    }
    return totalHamiltonian;
}

//-------------------------------------------------------------------------------------------------
std::vector<double> getProbabilityAtState(const std::vector<double> &temperatures, 
                                          const std::vector<double> &hamiltonian) {
  std::vector<double> probabilities;
  for (size_t i = 0; i < temperatures.size(); ++i) {
    double totalHamiltonian = getTotalHamiltonian(hamiltonian);
    double beta = 1.0 / (temperatures[i] * boltzmann_constant);
    probabilities.push_back(exp((-beta) * totalHamiltonian));
  }
  return probabilities;
}

} // namespace sampling
} // namespace stormm