// -*-c++-*-

using stormm::stmath::sum;
using stormm::stmath::mean;
using stormm::stmath::maxAbsValue;
using stormm::constants::tiny;
using stormm::topology::BondTerm;

//-------------------------------------------------------------------------------------------------
template <typename T>
std::vector<double> realAGProp(const std::vector<AtomGraph*> &topols, const bool all_top_exist,
                               const RealInfoCode analysis_code) {

  // Return a vector of blank data if no topologies have been filled out
  const int ntop = topols.size();
  if (all_top_exist == false) {
    return std::vector<double>(ntop, 0);
  }

  // Loop over all topologies and apply the requested analysis
  std::vector<double> result(ntop);
  for (int i = 0; i < ntop; i++) {
    switch (analysis_code) {
    case RealInfoCode::FIRST_MOLECULE_MASS:
      {
        const std::vector<int> mol_atoms = topols[i]->getMoleculeContents(0);
        const size_t n_mol_atoms = mol_atoms.size();
        double msum = 0.0;
        for (size_t j = 0; j < n_mol_atoms; j++) {
          msum += topols[i]->getAtomicMass<T>(j);
	}
        result[i] = msum;
      }
      break;
    case RealInfoCode::TOTAL_CHARGE:
      result[i] = sum<double>(topols[i]->getPartialCharge<T>());
      break;
    case RealInfoCode::MAX_ABS_CHARGE:
      result[i] = maxAbsValue<T>(topols[i]->getPartialCharge<T>());
      break;
    case RealInfoCode::AVERAGE_LJ_SIGMA:
      {
        const int natom = topols[i]->getAtomCount();
        std::vector<int> lj_idx = topols[i]->getLennardJonesIndex();
        double meansig = 0.0;
        int n_nonzero_sig = 0;
        for (int j = 0; j < natom; j++) {
          double tmpsig = topols[i]->getLennardJonesSigma<T>(lj_idx[j]);
          if (tmpsig > tiny) {
            meansig += tmpsig;
            n_nonzero_sig++;
          }
        }
        result[i] = meansig / static_cast<double>(n_nonzero_sig);
      }
      break;
    case RealInfoCode::AVERAGE_LJ_EPS:
      {
        const int natom = topols[i]->getAtomCount();
        std::vector<int> lj_idx = topols[i]->getLennardJonesIndex();
        double meaneps = 0.0;
        int n_nonzero_eps = 0;
        for (int j = 0; j < natom; j++) {
          double tmpeps = topols[i]->getLennardJonesEpsilon<T>(lj_idx[j]);
          if (tmpeps > tiny) {
            meaneps += tmpeps;
            n_nonzero_eps++;
          }
        }
        result[i] = meaneps / static_cast<double>(n_nonzero_eps);
      }
      break;
    case RealInfoCode::AVERAGE_BOND_STIFFNESS:
      {
        double mean_keq = 0.0;
        const int nbonds = topols[i]->getBondTermCount();
        for (int j = 0; j < nbonds; j++) {
          const BondTerm<T> btj = topols[i]->getBondTerm<T>(j);
          mean_keq += btj.keq;
        }
        result[i] = mean_keq / static_cast<double>(nbonds);
      }
      break;
    case RealInfoCode::AVERAGE_ANGL_THETA:
      {
        double mean_teq = 0.0;
        const int nangls = topols[i]->getAngleTermCount();
        for (int j = 0; j < nangls; j++) {
          const AngleTerm<T> atj = topols[i]->getAngleTerm<T>(j);
          mean_teq += atj.theta_eq;
        }
        result[i] = (nangls > 0) ? mean_teq / static_cast<double>(nangls) : 0.0;
      }
      break;
    case RealInfoCode::AVERAGE_DIHE_COMP:
      {
        double mean_comp = 0.0;
        const int ndihes = topols[i]->getDihedralTermCount();
        for (int j = 0; j < ndihes; j++) {
          const DihedralTerm<T> htj = topols[i]->getDihedralTerm<T>(j);
          mean_comp += htj.amplitude * htj.periodicity * htj.phase;
        }
        result[i] = (ndihes > 0) ? mean_comp / static_cast<double>(ndihes) : 0.0;
      }
      break;
    case RealInfoCode::AVERAGE_UBRD_COMP:
      {
        double mean_comp = 0.0;
        const int nubrds = topols[i]->getUreyBradleyTermCount();
        for (int j = 0; j < nubrds; j++) {
          const UreyBradleyTerm<T> utj = topols[i]->getUreyBradleyTerm<T>(j);
          mean_comp += utj.keq * utj.leq;
        }
        result[i] = (nubrds > 0) ? mean_comp / static_cast<double>(nubrds) : 0.0;
      }
      break;
    case RealInfoCode::AVERAGE_CIMP_COMP:
      {
        double mean_comp = 0.0;
        const int ncimps = topols[i]->getCharmmImprTermCount();
        for (int j = 0; j < ncimps; j++) {
          const CharmmImprTerm<T> citj = topols[i]->getCharmmImprTerm<T>(j);
          mean_comp += citj.keq + citj.phi_eq;
        }
        result[i] = (ncimps > 0) ? mean_comp / static_cast<double>(ncimps) : 0.0;
      }
      break;
    case RealInfoCode::AVERAGE_CMAP_VALUE:
      {
        double mean_value = 0.0;
        int nsample = 0;
        const int ncmaps = topols[i]->getCmapTermCount();
        for (int j = 0; j < ncmaps; j++) {
          const CmapTerm<T> ctj = topols[i]->getCmapTerm<T>(j);
          for (int k = 0; k < ctj.surf_dim * ctj.surf_dim; k++) {
            mean_value += ctj.surf[k];
          }
          nsample += ctj.surf_dim * ctj.surf_dim;
        }
        result[i] = (ncmaps > 0) ? mean_value / static_cast<double>(nsample) : 0.0;
      }
      break;
    }
  }
  return result;
}

