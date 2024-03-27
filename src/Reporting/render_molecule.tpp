// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace review {

//-------------------------------------------------------------------------------------------------
template <typename T>
void renderMolecule(std::ofstream *foutp, const T* atom_x, const T* atom_y, const T* atom_z,
                    const NonbondedKit<double> &nbk, const ChemicalDetailsKit &cdk,
                    const RenderOptions &ropt, const GridFileSyntax syntax, const size_t nframe) {

  // Begin at atom zero, tracing a list of all possible connections and keeping a history for
  // back-tracking purposes.  This is similar to the ring-tracing algorithm.
  std::vector<uint> atom_coverage((nbk.natom + (sizeof(uint) * 8) - 1) / sizeof(uint), 0U);
  std::vector<uint> bond_coverage((nbk.nb12_bounds[nbk.natom] +
                                   (sizeof(uint) * 8) - 1) / sizeof(uint), 0U);
  std::vector<int> history;
  
  // Loop over all molecules in the topology.  Plot their lines for all structures.
  int n_water_oxygens = 0;
  int n_water_hydrogens = 0;
  std::vector<bool> is_water(cdk.nmol, false);
  int trace_count = 0;
  const size_t frame_natom = nbk.natom;
  const size_t padded_natom = roundUp(cdk.natom, warp_size_int);
  for (int i = 0; i < cdk.nmol; i++) {

    // Skip water molecules--these will be displayed with special markers but would be too
    // taxing on the matrix package plotter to display each molecule as an individual plot.
    int n_real_atom = 0;
    int n_oxygen = 0;
    int n_hydrogen = 0;
    for (int j = cdk.mol_limits[i]; j < cdk.mol_limits[i + 1]; j++) {
      const int atom_idx = cdk.mol_contents[j];
      const int atom_zn  = cdk.z_numbers[atom_idx]; 
      n_real_atom += (atom_zn > 0);
      n_oxygen += (atom_zn == 8);
      n_hydrogen += (atom_zn == 1);
    }
    if (n_real_atom == 3 && n_hydrogen == 2 && n_oxygen == 1) {
      n_water_oxygens++;
      n_water_hydrogens += 2;
      is_water[i] = true;
    }
    else {
      trace_count++;
    }
  }
  std::vector<std::vector<int>> atom_traces(trace_count);
  std::vector<int> n_traced_elements(element_maximum_count, 0);
  trace_count = 0;
  for (int i = 0; i < cdk.nmol; i++) {
    if (is_water[i]) {
      continue;
    }
    atom_traces[trace_count] = traceBondLines(nbk, cdk, i, &atom_coverage, &bond_coverage);
    for (int j = cdk.mol_limits[i]; j < cdk.mol_limits[i + 1]; j++) {
      const int atom_zn = cdk.z_numbers[cdk.mol_contents[j]];
      n_traced_elements[atom_zn] += 1;
    }
    trace_count++;
  }
  
  // Multiply the numbers of traced elements and water atoms.  Allocate arrays to hold the
  // matrix package arrays for printing.
  std::vector<int> water_oxygens, water_hydrogens;
  water_oxygens.reserve(n_water_oxygens);
  water_hydrogens.reserve(n_water_hydrogens);
  std::vector<std::vector<int>> traced_elements(element_maximum_count);
  for (int i = 0; i < element_maximum_count; i++) {
    traced_elements[i].reserve(n_traced_elements[i]);
  }
  for (int i = 0; i < cdk.nmol; i++) {
    if (is_water[i]) {
      continue;
    }
    for (int j = cdk.mol_limits[i]; j < cdk.mol_limits[i + 1]; j++) {
      const int j_atom = cdk.mol_contents[j];
      traced_elements[cdk.z_numbers[j_atom]].push_back(j_atom);
    }
  }
  
  // Fill out the water oxygen and hydrogen arrays
  for (int i = 0; i < cdk.nmol; i++) {
    if (is_water[i]) {
      for (int j = cdk.mol_limits[i]; j < cdk.mol_limits[i + 1]; j++) {
        const int j_atom = cdk.mol_contents[j];
        const int j_znum = cdk.z_numbers[j_atom];
        if (j_znum == 1) {
          water_hydrogens.push_back(j_atom);
        }
        else if (j_znum == 8) {
          water_oxygens.push_back(j_atom);
        }
      }
    }
  }

  // Write the result, beginning with the complete array of all atom positions.
  const NumberFormat real_fmt = NumberFormat::STANDARD_REAL;
  std::string result;
  switch (syntax) {
  case GridFileSyntax::MATPLOTLIB:
    break;
  case GridFileSyntax::MATRIX_PKG:
    for (size_t fpos = 0; fpos < nframe; fpos++) {
      result.reserve((30LLU * static_cast<size_t>(nbk.natom)) + 64LLU);
      result += "frame_count = " + std::to_string(nframe) + ";\n";
      result += "atomic_coordinates = zeros(" + std::to_string(nbk.natom) + ", 3, frame_count);\n";
      result += "atomic_coordinates(:,:," + std::to_string(fpos + 1) + ") = [\n";
      for (size_t i = 0; i < frame_natom; i++) {
        result += realToString(atom_x[(padded_natom * fpos) + i], 9, 4, real_fmt) + " " +
                  realToString(atom_y[(padded_natom * fpos) + i], 9, 4, real_fmt) + " " +
                  realToString(atom_z[(padded_natom * fpos) + i], 9, 4, real_fmt) + "\n";
      }
      result += "];\n";
      foutp->write(result.data(), result.size());
    }
    break;    
  case GridFileSyntax::OPEN_DX:
    break;
  case GridFileSyntax::CUBEGEN:
    rtErr("Gaussian cube files are not paired with any visualization package capable of rendering "
          "molecules in the scene.");
  }
  
  // List the water atom positions.
  switch (syntax) {
  case GridFileSyntax::MATPLOTLIB:
    break;
  case GridFileSyntax::MATRIX_PKG:
    if (n_water_oxygens > 0) {
      result.reserve((7LLU * static_cast<size_t>(n_water_oxygens)) + 64LLU);
      result += "water_ox_idx = [\n";
      size_t line_length = 0;
      for (int i = 0; i < n_water_oxygens; i++) {
        const std::string next_number = std::to_string(water_oxygens[i] + 1) + " ";
        if (line_length + next_number.size() >= 116) {
          result += "...\n";
          line_length = 0;
        }
        result += next_number;
        line_length += next_number.size();
      }
      result += "water_hd_idx = [\n";
      line_length = 0;
      for (int i = 0; i < n_water_hydrogens; i++) {
        const std::string next_number = std::to_string(water_hydrogens[i] + 1) + " ";
        if (line_length + next_number.size() >= 116) {
          result += "...\n";
          line_length = 0;
        }
        result += next_number;
        line_length += next_number.size();
      }
      result += "\n];\n";
      result += "per_frame_waters = size(water_ox_idx, 2);\n";
      result += "water_ox = zeros(" + std::to_string(n_water_oxygens) + ",3);\n";
      result += "water_hd = zeros(" + std::to_string(n_water_hydrogens) + ",3);\n";
      result += "h = 0;\n";
      result += "for int i = 1:1:" + std::to_string(nframe) + "\n";
      result += "  for int j = 1:1:per_frame_waters" + std::to_string(nframe) + "\n";
      result += "    water_ox(h + 1,:) = atomic_coordinates(water_ox_idx(j),:,i);\n";
      result += "    water_hd((2 * h) + 1,:) = "
                "atomic_coordinates(water_hd_idx((2 * j) + 1),:,i);\n";
      result += "    water_hd((2 * h) + 2,:) = "
                "atomic_coordinates(water_hd_idx((2 * j) + 2),:,i);\n";
      result += "    h = h + 1;\n";
      result += "  end\n";
      result += "end\n";
      result += "];\n";
      foutp->write(result.data(), result.size());
      result.resize(0);
      result += "plot3(water_ox(:,1), water_ox(:,2), water_ox(:,3), '.', 'markersize', 32, "
                "'color', [ 0.878 0.000 0.204 ]);\n";
      result += "plot3(water_hd(:,1), water_hd(:,2), water_hd(:,3), '.', 'markersize', 32, "
                "'color', [ 0.950 0.950 0.950 ]);\n";
      foutp->write(result.data(), result.size());
    }

    // Provide the atom indices for each trace and write script code to construct them based on
    // the general cache of atomic coordinates.
    for (int i = 0; i < trace_count; i++) {
      result.resize(0);
      const size_t tr_length = atom_traces[i].size();
      result.reserve((7LLU * tr_length) + 64LLU);
      result += "moltrace_" + std::to_string(i) + "_idx = [\n";
      size_t line_length = 0;
      for (size_t j = 0; j < tr_length; j++) {
        const std::string next_number = std::to_string(atom_traces[i][j] + 1) + " ";
        if (line_length + next_number.size() >= 116) {
          result += "...\n";
          line_length = 0;
        }
        result += next_number;
        line_length += next_number.size();
      }
      result += "\n];\n";
      foutp->write(result.data(), result.size());
    }
    
    // Add script to assemble each molecular trace based on the coordinates given to the matrix
    // software package.
    for (size_t fpos = 0; fpos < nframe; fpos++) {
      const size_t fr_offset = padded_natom * fpos;
      for (int i = 0; i < trace_count; i++) {
        result.resize(0);
        const size_t tr_length = atom_traces[i].size();
        result += "moltrace_" + std::to_string(i) + " = zeros(" + std::to_string(tr_length) +
                  ", 3, frame_count);\n";
        result += "for frm = 1:1:" + std::to_string(nframe) + "\n";
        result += "  for i = 1:1:" + std::to_string(tr_length) + "\n";
        result += "    moltrace_" + std::to_string(i) + "(i,:,frm) = "
                  "atomic_coordinates(moltrace_" + std::to_string(i) + "_idx(i),:,frm);\n";
        result += "  end\n";
        result += "end\n";
        foutp->write(result.data(), result.size());
      }
    }
    for (size_t fpos = 0; fpos < nframe; fpos++) {
      for (int i = 0; i < trace_count; i++) {
        result = "plot3(moltrace_" + intToString(i) + "(:,1," + intToString(fpos + 1) + "), " +
                 "moltrace_" + intToString(i) + "(:,2," + intToString(fpos + 1) + "), " +
                 "moltrace_" + intToString(i) + "(:,3," + intToString(fpos + 1) + "), " +
                 ropt.getFormattedReceptorLineColor(syntax) + ", 'linewidth', " +
                 realToString(ropt.getReceptorLineWidth(), 5, 2, real_fmt) + ");\n";
        foutp->write(result.data(), result.size());
      }
    }

    // Create scatter plots of all elements in the traced systems.
    for (int i = 0; i < element_maximum_count; i++) {
      const size_t n_elemental_atoms = traced_elements[i].size();
      if (n_elemental_atoms > 0) {
        result.resize(0);
        result.reserve((7LLU * n_elemental_atoms) + 64LLU);
        const std::string atm_prfx = removeTailingWhiteSpace(char2ToString(elemental_symbols[i])) +
                                     "_atom";
        result = atm_prfx + "_idx = [\n";
        size_t line_length = 0;
        for (size_t j = 0; j < n_elemental_atoms; j++) {
          const std::string next_number = std::to_string(traced_elements[i][j] + 1) + " ";
          if (line_length + next_number.size() >= 116) {
            result += "...\n";
            line_length = 0;
          }
          result += next_number;
          line_length += next_number.size();
        }
        result += "\n];\n";
        result += atm_prfx + "s = zeros(" + std::to_string(n_elemental_atoms) +
                  ", 3, frame_count);\n";
        result += "for frm = 1:1:" + std::to_string(nframe) + "\n";
        result += "  for i = 1:1:" + std::to_string(n_elemental_atoms) + "\n";
        result += "    " + atm_prfx + "s(i,:,frm) = atomic_coordinates(" + atm_prfx +
                  "_idx(i),:,frm);\n";
        result += "  end\n";
        result += "end\n";
        foutp->write(result.data(), result.size());
      }
    }
    for (int i = 0; i < element_maximum_count; i++) {
      const size_t n_atoms = traced_elements[i].size();
      if (n_atoms > 0) {
        const std::string varname = removeTailingWhiteSpace(char2ToString(elemental_symbols[i])) +
                                    "_atoms";
        const int marker_size = (i <= 4) ? ropt.getReceptorLightAtomSize() :
                                           ropt.getReceptorHeavyAtomSize();
        result = "plot3(" + varname + "(:,1), " + varname + "(:,2), " + varname + "(:,3), " +
                 "'.', 'markersize', " + std::to_string(marker_size) + ", " +
          ropt.getFormattedElementColor(i, syntax) + ");\n";
        foutp->write(result.data(), result.size());
      }
    }
    break;
  case GridFileSyntax::OPEN_DX:
  case GridFileSyntax::CUBEGEN:
    break;
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void renderMolecule(std::ofstream *foutp, const CoordinateSeries<T> *cs, const AtomGraph *ag,
                    const RenderOptions *ropt, const GridFileSyntax syntax) {
  const CoordinateSeriesReader csr = cs->data();
  const NonbondedKit<double> nbk = ag->getDoublePrecisionNonbondedKit();
  const ChemicalDetailsKit cdk = ag->getChemicalDetailsKit();
  renderMolecule(foutp, csr.xcrd, csr.ycrd, csr.zcrd, nbk, cdk, *ropt, syntax, csr.nframe);
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void renderMolecule(std::ofstream *foutp, const CoordinateSeries<T> &cs, const AtomGraph &ag,
                    const RenderOptions &ropt, const GridFileSyntax syntax) {
  renderMolecule(foutp, cs.getSelfPointer(), ag.getSelfPointer(), ropt.getSelfPointer(), syntax);
}

} // namespace review
} // namespace stormm
