// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace review {

//-------------------------------------------------------------------------------------------------
template <typename T>
void printToFile(const BackgroundMesh<T> &bgm, const MeshParameters &mps,
                 const std::string &file_name, const PrintSituation expectation,
                 const GridFileSyntax syntax, const RenderOptions &ropt, const int decimal_places,
                 const ExceptionResponse policy) {
  
  // Check that the dimensions of the mesh to print are compatible with the mesh provided.
  const MeshParamKit input_dims = bgm.getDimensions().data();
  const MeshParamKit output_dims = mps.data();
  const std::string input_mesh_desc  = bgm.getDimensions().printDimensions();
  const std::string output_mesh_desc = mps.printDimensions();
  const NumberFormat real_fmt = NumberFormat::STANDARD_REAL;
  const double output_orig_x = hostInt95ToDouble(output_dims.orig_x) * output_dims.inv_scale;
  const double output_orig_y = hostInt95ToDouble(output_dims.orig_y) * output_dims.inv_scale;
  const double output_orig_z = hostInt95ToDouble(output_dims.orig_z) * output_dims.inv_scale;
  switch (input_dims.bounds) {
  case BoundaryCondition::ISOLATED:
    {
      // A periodic mesh can only be created based on an existing periodic mesh
      switch (output_dims.bounds) {
      case BoundaryCondition::ISOLATED:
        break;
      case BoundaryCondition::PERIODIC:
        rtErr("An output mesh with " + getEnumerationName(output_dims.bounds) + " boundary "
              "conditions cannot be output from a mesh with " +
              getEnumerationName(input_dims.bounds) + " boundary conditions.", "printToFile");
      }
    
      // Check the minimum and maximum dimensions.  If the dimensions called for are very near to
      // validity but violate it for something that looks like roundoff error, adjust the requested
      // output dimensions.
      const double input_orig_x = hostInt95ToDouble(input_dims.orig_x) * input_dims.inv_scale;
      const double input_orig_y = hostInt95ToDouble(input_dims.orig_y) * input_dims.inv_scale;
      const double input_orig_z = hostInt95ToDouble(input_dims.orig_z) * input_dims.inv_scale;
      const double dx_orig = output_orig_x - input_orig_x;
      const double dy_orig = output_orig_y - input_orig_y;
      const double dz_orig = output_orig_z - input_orig_z;
      for (int i = 0; i < 2; i++) {
        const double ldi = static_cast<double>(output_dims.na * i);
        for (int j = 0; j < 2; j++) {
          const double ldj = static_cast<double>(output_dims.nb * j);
          for (int k = 0; k < 2; k++) {
            const double ldk = static_cast<double>(output_dims.nc * k);
            const double rel_xpt = (output_dims.invu[0] * ldi) + (output_dims.invu[3] * ldj) +
                                   (output_dims.invu[6] * ldk) + dx_orig;
            const double rel_ypt = (output_dims.invu[4] * ldj) + (output_dims.invu[7] * ldk) +
                                   dy_orig;
            const double rel_zpt = (output_dims.invu[8] * ldk) + dz_orig;
            const double disp_na = (input_dims.umat[0] * rel_xpt) +
                                   (input_dims.umat[3] * rel_ypt) +
                                   (input_dims.umat[6] * rel_zpt);
            const double disp_nb = (input_dims.umat[4] * rel_ypt) + (input_dims.umat[7] * rel_zpt);
            const double disp_nc = (input_dims.umat[8] * rel_zpt);
            if (disp_na < -constants::small || disp_nb < -constants::small ||
                disp_nc < -constants::small ||
                disp_na > static_cast<double>(input_dims.na + 1) + constants::small ||
                disp_nb > static_cast<double>(input_dims.nb + 1) + constants::small ||
                disp_nc > static_cast<double>(input_dims.nc + 1) + constants::small) {
              rtErr("A mesh with " + getEnumerationName(input_dims.bounds) + " boundary "
                    "conditions can only produce data for a mesh within its boundaries.  A "
                    "request for data at (" + realToString(output_orig_x, 9, 4, real_fmt) + ", " +
                    realToString(output_orig_y, 9, 4, real_fmt) + ", " +
                    realToString(output_orig_z, 9, 4, real_fmt) + ") is invalid.  [Input] " +
                    input_mesh_desc + ".  Requested point: corner {" + std::to_string(i) + "," +
                    std::to_string(j) + "," + std::to_string(k) + "}.  [Output] " +
                    output_mesh_desc + ".", "printToFile");
            }
          }
        }
      }
    }
    break;
  case BoundaryCondition::PERIODIC:
    break;
  }

  // Do not allow printing of excessively large meshes
  if (static_cast<size_t>(output_dims.na) * static_cast<size_t>(output_dims.nb) *
      static_cast<size_t>(output_dims.nc) > static_cast<size_t>(INT_MAX)) {
    rtErr("Output dimensions " + std::to_string(output_dims.na) + " x " +
          std::to_string(output_dims.nb) + " x " + std::to_string(output_dims.nc) +
          " would create too large a mesh, not least in terms of file size.  Choose a number of "
          "points <= " + std::to_string(INT_MAX) + ".", "printToFile");
  }
  const int npts_out = output_dims.na * output_dims.nb * output_dims.nc;
  
  // Fill out the preamble and closing text for whatever file type.
  std::string preamble, post_script;
  const MeshFoundation& bgm_bss = bgm.getMolecularBasis();
  const int n_comments = bgm_bss.getCommentCount();
  const std::vector<double> torig = {
    hostInt95ToDouble(output_dims.orig_x) * output_dims.inv_scale,
    hostInt95ToDouble(output_dims.orig_y) * output_dims.inv_scale,
    hostInt95ToDouble(output_dims.orig_z) * output_dims.inv_scale
  };
  const int orig_format_length = findAlignmentWidth(torig, 9);
  const int xform_format_length = findAlignmentWidth(output_dims.invu, 9, 9);
  const char comment_symbol = commentSymbol(syntax);

  // Determine the content format length based on the requested number of decimal places.
  const BackgroundMeshReader<T> bgmr = bgm.data();
  const MeshFFKit<double> mnbk = bgm.getReferenceNonbondedKit();
  std::vector<double> prop_a, prop_b;
  std::vector<int> index_zero(output_dims.nc, 0);
  switch (bgm.getNonbondedPotential()) {
  case NonbondedPotential::VAN_DER_WAALS:
    {
      const double pr_six = pow(bgm.getProbeRadius(), 6.0);
      const double wd_four = 4.0 * bgm.getWellDepth();
      const double lj_bcoef = wd_four * pr_six;
      const double lj_acoef = lj_bcoef * pr_six;
      prop_a.resize(1, lj_acoef);
      prop_b.resize(1, lj_bcoef);
    }
    break;
  case NonbondedPotential::ELECTROSTATIC:
    prop_a.resize(output_dims.nc, 1.0);
    break;
  case NonbondedPotential::CLASH:
    break;
  }
  int content_format_length = 0;
  std::vector<double> xpt_out(output_dims.nc), ypt_out(output_dims.nc);
  std::vector<double> zpt_out(output_dims.nc), u_intrp(output_dims.nc);
  for (int i = 0; i < output_dims.na; i++) {
    const double di = i;
    for (int j = 0; j < output_dims.nb; j++) {
      const double dj = j;
      for (int k = 0; k < output_dims.nc; k++) {
        const double dk = k;
        xpt_out[k] = (di * output_dims.invu[0]) + (dj * output_dims.invu[3]) +
                     (dk * output_dims.invu[6]) + output_orig_x;
        ypt_out[k] = (dj * output_dims.invu[4]) + (dk * output_dims.invu[7]) + output_orig_y;
        zpt_out[k] = (dk * output_dims.invu[8]) + output_orig_z;
      }
      interpolate<T, double, double, double>(bgmr, mnbk, prop_a.data(), prop_b.data(),
                                             index_zero.data(), xpt_out.data(), ypt_out.data(),
                                             zpt_out.data(), nullptr, nullptr, nullptr,
                                             output_dims.nc, u_intrp.data(), nullptr, nullptr,
                                             nullptr, nullptr, nullptr, nullptr,
                                             OffMeshProtocol::EXTRAPOLATE);
      content_format_length = std::max(content_format_length,
                                       findAlignmentWidth(u_intrp, decimal_places));
    }
  }
  
  // In some formats, data may go to an auxiliary file for faster loading.  Set that file's name.
  const std::string data_file_name = getSceneDataFileName(file_name);
  switch (syntax) {
  case GridFileSyntax::MATPLOTLIB:
    break;
  case GridFileSyntax::MATRIX_PKG:
    {
      for (int i = 0; i < n_comments; i++) {
        preamble += protectText(bgm_bss.getComment(i), comment_symbol);
      }
      preamble += "clear all\nclose all\n\n";
      preamble += "na = " + std::to_string(output_dims.na) + ";\n";
      preamble += "nb = " + std::to_string(output_dims.nb) + ";\n";
      preamble += "nc = " + std::to_string(output_dims.nc) + ";\n";
      preamble += "x = zeros(na, nb, nc);\n";
      preamble += "y = zeros(na, nb, nc);\n";
      preamble += "z = zeros(na, nb, nc);\n";
      for (int i = 0; i < 3; i++) {
        preamble += (i == 0) ? "invU = [ " : "         ";
        preamble += realToString(output_dims.invu[i    ], xform_format_length, 9) + " " +
                    realToString(output_dims.invu[i + 3], xform_format_length, 9) + " " +
                    realToString(output_dims.invu[i + 6], xform_format_length, 9);
        preamble += (i == 2) ? " ];\n" : "\n";
      }
      preamble += "mg_origin = [ " +
                  realToString(output_orig_x, orig_format_length, 9, real_fmt) + " " +
                  realToString(output_orig_y, orig_format_length, 9, real_fmt) + " " +
                  realToString(output_orig_z, orig_format_length, 9, real_fmt) + " ];\n";
      preamble += "for i = 1:1:na\n  for j = 1:1:nb\n    for k = 1:1:nc\n";
      preamble += "      pt = invU * [ i - 1; j - 1; k - 1 ];\n";
      preamble += "      x(i,j,k) = pt(1) + mg_origin(1);\n";
      preamble += "      y(i,j,k) = pt(2) + mg_origin(2);\n";
      preamble += "      z(i,j,k) = pt(3) + mg_origin(3);\n";
      preamble += "    end\n  end\nend\n\n";
      preamble += protectText("Load the field data from an auxiliary file.", comment_symbol);
      preamble += "dock_field = load(\"" + data_file_name + "\");\n\n";
      preamble += "dock_field = reshape(dock_field, na, nb, nc);\n";
      preamble += protectText("Check if the program running the script is GNU Octave.  As of "
                              "version 6.4.0, Octave did not support 'facealpha' modifications in "
                              "the \"patch\" command.", comment_symbol);
      preamble += "is_octave = (exist('OCTAVE_VERSION', 'builtin') ~= 0);\n";
      preamble += "figure;\nhold on;\n";
      const int n_isosurf = ropt.getIsosurfaceCount();
      bool octave_warn = false;
      for (int i = 0; i < n_isosurf; i++) {
        post_script += "[ f, vert ] = isosurface(x, y, z, dock_field, " +
                       realToString(ropt.getIsosurfaceValue(i), content_format_length,
                                    decimal_places, real_fmt) + ");\n";
        const uchar4 fcolor = ropt.getIsosurfaceColor(i);
        const double f_red = static_cast<double>(fcolor.x) / 255.0;
        const double f_grn = static_cast<double>(fcolor.y) / 255.0;
        const double f_blu = static_cast<double>(fcolor.z) / 255.0;
        const double f_alf = static_cast<double>(fcolor.w) / 255.0;
        switch (ropt.getIsosurfaceDrawingMethod(i)) {
        case SurfaceRender::SOLID:
          post_script += "patch('faces', f, 'vertices', vert, 'facecolor', [ " +
                         realToString(f_red, 6, 4, real_fmt) + ", " +
                         realToString(f_grn, 6, 4, real_fmt) + ", " +
                         realToString(f_blu, 6, 4, real_fmt) + " ], "
                         "'facealpha', " + realToString(f_alf, 6, 4, real_fmt) +
                         ", 'edgecolor', 'none');\n";
          break;
        case SurfaceRender::WIRE:
          post_script += "patch('faces', f, 'vertices', vert, 'facecolor', 'none', "
                         "'edgecolor', [ " +
                         realToString(f_red, 6, 4, real_fmt) + ", " +
                         realToString(f_grn, 6, 4, real_fmt) + ", " +
                         realToString(f_blu, 6, 4, real_fmt) + " ]);\n";
          break;
        case SurfaceRender::SCAFFOLD:
          post_script += "patch('faces', f, 'vertices', vert, 'facecolor', [ " +
                         realToString(f_red, 6, 4, real_fmt) + ", " +
                         realToString(f_grn, 6, 4, real_fmt) + ", " +
                         realToString(f_blu, 6, 4, real_fmt) + " ], "
                         "'facealpha', " + realToString(f_alf, 6, 4, real_fmt) +
                         ", 'edgecolor', [ 0.1, 0.1, 0.1 ]);\n";
          break;
        }
        switch (ropt.getIsosurfaceDrawingMethod(i)) {
        case SurfaceRender::SOLID:
        case SurfaceRender::SCAFFOLD:
          if (fcolor.w != 255 && octave_warn == false) {
            octave_warn = true;
            post_script += "if (is_octave == 1)\n";
            post_script += "  fprintf('GNU Octave does not support a \"facealpha\" property in "
                           "its \"patch\" command.\n');\n";
            post_script += "  fprintf('Plot the field isosurfaces using a wire mesh outline "
                           "instead, or use a licensed\n');\n";
            post_script += "  fprintf('copy of Matlab which does support alpha channels in the "
                           "\"patch\" command.\n');\n";
            post_script += "end\n";
          }
          break;
        case SurfaceRender::WIRE:
          break;
        }
      }
      post_script += "daspect([ 1, 1, 1 ]);\n";
    }
    break;
  case GridFileSyntax::OPEN_DX:
    for (int i = 0; i < n_comments; i++) {
      preamble += protectText(bgm_bss.getComment(i), comment_symbol);
    }
    preamble += "object 1 class gridpositions counts " + std::to_string(output_dims.na) + " " +
                std::to_string(output_dims.nb) + " " + std::to_string(output_dims.nc) + "\n";
    preamble += "origin " +
                realToString(output_orig_x, orig_format_length, 9, real_fmt) + " " +
                realToString(output_orig_y, orig_format_length, 9, real_fmt) + " " +
                realToString(output_orig_z, orig_format_length, 9, real_fmt) + "\n";
    for (int i = 0; i < 3; i++) {
      preamble += "delta " +
                  realToString(output_dims.invu[i    ], xform_format_length, 9, real_fmt) + " " +
                  realToString(output_dims.invu[i + 3], xform_format_length, 9, real_fmt) + " " +
                  realToString(output_dims.invu[i + 6], xform_format_length, 9, real_fmt) + "\n";
    }
    preamble += "object 2 class gridconnections counts " + std::to_string(output_dims.na) + " " +
                std::to_string(output_dims.nb) + " " + std::to_string(output_dims.nc) + "\n";
    preamble += "object 3 class array type double rank 0 items " + std::to_string(npts_out) +
                " data follows\n";
    post_script += "attribute \"dep\" string \"positions\"\n";
    post_script += "object \"regular positions regular connections\" class field\n";
    post_script += "component \"positions\" value 1\n";
    post_script += "component \"connections\" value 2\n";
    post_script += "component \"data\" value 3\n";
    break;
  case GridFileSyntax::CUBEGEN:
    {
      for (int i = 0; i < std::min(n_comments, 2); i++) {
        preamble += protectText(bgm_bss.getComment(i), comment_symbol,
                                bgm_bss.getComment(i).size() + 8);
      }
      if (n_comments > 2) {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("The Gaussian cubegen format accepts at most two comment lines.  A total of " +
                std::to_string(n_comments) + " comments were logged in the BackgroundMesh object.",
                "printToFile");
        case ExceptionResponse::WARN:
          rtWarn("The Gaussian cubegen format accepts at most two comment lines, but a total of " +
                 std::to_string(n_comments) + " comments were logged in the BackgroundMesh "
                 "object.  Excess comments will be ignored.", "printToFile");
          break;
        case ExceptionResponse::SILENT:
          break;
        }
      }
      else {
        for (int i = n_comments; i < 2; i++) {
          if (i == 0) {
            const std::string first_comment("Generated by STORMM");
            preamble += protectText(first_comment, comment_symbol, first_comment.size() + 8);
          }
          else {
            const std::string second_comment("OUTER LOOP: X, MIDDLE LOOP: Y, INNER LOOP: Z");
            preamble += protectText(second_comment, comment_symbol, second_comment.size() + 8);
          }
        }
      }
      const AtomGraph *ag_ptr = bgm_bss.getTopologyPointer();
      
      if (ag_ptr == nullptr) {
        switch (policy) {
        case ExceptionResponse::DIE:
          rtErr("No topology was referenced by the BackgroundMesh object.", "printToFile");
        case ExceptionResponse::WARN:
          rtWarn("No topology was referenced by the BackgroundMesh object.  No Cube format file "
                 "can be printed without molecular information.", "printToFile");
          return;
        case ExceptionResponse::SILENT:
          return;
        }
      }
      preamble += std::to_string(ag_ptr->getAtomCount()) +
                  realToString(output_orig_x, orig_format_length, 9, real_fmt) + " " +
                  realToString(output_orig_y, orig_format_length, 9, real_fmt) + " " +
                  realToString(output_orig_z, orig_format_length, 9, real_fmt) + "\n";
      const std::vector<int> n_voxels = { output_dims.na, output_dims.nb, output_dims.nc };
      const int voxel_format_length = findAlignmentWidth<int>(n_voxels, 0);

      // Present the geometry in units of Bohr, standard for Gaussian Cube files.  While the format
      // will imply Angstrom units by specifying negative voxel (mesh element) counts along each
      // axis, other pgorams may not understand the negative counts.
      std::vector<double> bohr_invu(9);
      for (int i = 0; i < 9; i++) {
        bohr_invu[i] = angstrom_to_bohr * output_dims.invu[i];
      }
      for (int i = 0; i < 3; i++) {
        preamble += intToString(n_voxels[i], voxel_format_length) + " " +
                    realToString(bohr_invu[i    ], xform_format_length, 9, real_fmt) + " " +
                    realToString(bohr_invu[i + 3], xform_format_length, 9, real_fmt) + " " +
                    realToString(bohr_invu[i + 6], xform_format_length, 9, real_fmt) + "\n";
      }

      // Print the atom information.  This could originate from a coordinate frame or series.
      CoordinateFrame cf(ag_ptr->getAtomCount());
      if (bgm_bss.getCoordinatePointer() != nullptr) {
        coordCopy(&cf, *(bgm_bss.getCoordinatePointer()));
      }
      else {

        // The BackgroundMesh must reference a coordinate series if not a single frame.
        const size_t ens_code = bgm_bss.getEnsembleTypeCode();
        if (ens_code == double_type_index) {
          coordCopy(&cf, *(bgm_bss.getEnsemblePointer<double>()), 0);
        }
        else if (ens_code == float_type_index) {
          coordCopy(&cf, *(bgm_bss.getEnsemblePointer<float>()), 0);
        }
        else if (ens_code == llint_type_index) {
          coordCopy(&cf, *(bgm_bss.getEnsemblePointer<llint>()), 0);
        }
        else if (ens_code == int_type_index) {
          coordCopy(&cf, *(bgm_bss.getEnsemblePointer<int>()), 0);
        }
        else if (ens_code == short_type_index) {
          coordCopy(&cf, *(bgm_bss.getEnsemblePointer<short int>()), 0);
        }
        else {
          switch (policy) {
          case ExceptionResponse::DIE:
            rtErr("The coordinate ensemble type in the BackgroundMesh object is not sanctioned, "
                  "or the BackgroundMesh may contain no coordinate structure.", "printToFile");
          case ExceptionResponse::WARN:
            rtWarn("The coordinate ensemble type in the BackgroundMesh object is not sanctioned, "
                   "or the BackgroundMesh may contain no coordinate structure.  Coordinates for "
                   "the molecule will be printed as zeros in the Cube file.", "printToFile");
            break;
          case ExceptionResponse::SILENT:
            break;
          }
        }
      }
      const CoordinateFrameReader cfr = cf.data();
      const ChemicalDetailsKit cdk = ag_ptr->getChemicalDetailsKit();
      const NonbondedKit<double> nbk = ag_ptr->getDoublePrecisionNonbondedKit();
      for (int i = 0; i < cdk.natom; i++) {
        preamble += intToString(cdk.z_numbers[i], 3) + " " +
                    realToString(nbk.charge[i], 8, 5, real_fmt) + " " +
                    realToString(angstrom_to_bohr * cfr.xcrd[i], 10, 5, real_fmt) + " " +
                    realToString(angstrom_to_bohr * cfr.ycrd[i], 10, 5, real_fmt) + " " +
                    realToString(angstrom_to_bohr * cfr.zcrd[i], 10, 5, real_fmt) + "\n";
      }
      post_script = "";
    }
    break;
  }

  // Write the contents to disk, beginning with the preamble
  std::ofstream foutp = openOutputFile(file_name, expectation, "Print a mesh to " +
                                       getEnumerationName(syntax) + " format.");
  foutp.write(preamble.data(), preamble.size());

  // Reprint the mesh data, transposed as necessary for the format
  std::vector<PolyNumeric> out_buffer;
  switch (syntax) {
  case GridFileSyntax::MATPLOTLIB:
    break;
  case GridFileSyntax::MATRIX_PKG:
    {
      // Seek an output file width of about 100 characters.
      const std::vector<uint> primes = { 2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31 };
      const int n_primes = primes.size();
      std::vector<uint> factors = primeFactors(npts_out, primes);
      factors.push_back(output_dims.nb);
      const int target_per_line_count = (120 + content_format_length - 1) / content_format_length;
      const int per_line_count = nearestFactor(npts_out, target_per_line_count, factors,
                                               LimitApproach::ABOVE);
      const int data_nrow = npts_out / per_line_count;
      std::ofstream data_foutp = openOutputFile(data_file_name, expectation, "Print data for a "
                                                "mesh to " + getEnumerationName(syntax) +
                                                " format.");
      
      // Resize and populate supporting vectors to handle up to 100 rows at a time.
      u_intrp.resize(100 * per_line_count);
      xpt_out.resize(100 * per_line_count);
      ypt_out.resize(100 * per_line_count);
      zpt_out.resize(100 * per_line_count);
      switch (bgm.getNonbondedPotential()) {
      case NonbondedPotential::VAN_DER_WAALS:
        index_zero.resize(100 * per_line_count, 0);
        break;
      case NonbondedPotential::CLASH:
      case NonbondedPotential::ELECTROSTATIC:
        prop_a.resize(100 * per_line_count, 1.0);
        break;
      }
      int lines_buffered = 0;
      const int output_ab_slab = output_dims.na * output_dims.nb;
      int h = 0;
      for (int i = 0; i < data_nrow; i++) {
        for (int j = 0; j < per_line_count; j++) {
          const int native_idx = (j * data_nrow) + i;
          const int native_k = native_idx / output_ab_slab;
          const int native_j = (native_idx - (native_k * output_ab_slab)) / output_dims.na;
          const int native_i = native_idx - (native_k * output_ab_slab) -
                               (native_j * output_dims.na);
          const double di = native_i;
          const double dj = native_j;
          const double dk = native_k;
          xpt_out[h] = (output_dims.invu[0] * di) + (output_dims.invu[3] * dj) +
                       (output_dims.invu[6] * dk) + output_orig_x;
          ypt_out[h] = (output_dims.invu[4] * dj) + (output_dims.invu[7] * dk) + output_orig_y;
          zpt_out[h] = (output_dims.invu[8] * dk) + output_orig_z;
          h++;
        }
        lines_buffered++;
        if (lines_buffered == 100 || i == data_nrow - 1) {
          interpolate<T, double, double, double>(bgmr, mnbk, prop_a.data(), prop_b.data(),
                                                 index_zero.data(), xpt_out.data(), ypt_out.data(),
                                                 zpt_out.data(), nullptr, nullptr, nullptr, h,
                                                 u_intrp.data(), nullptr, nullptr, nullptr,
                                                 nullptr, nullptr, nullptr,
                                                 OffMeshProtocol::EXTRAPOLATE, 0, 0, 0);

          // Resize the output buffer to meet the demand.  It will likely be re-allocated once on
          // the first pass and may be re-allocated once more on the final pass if the remaining
          // number of rows is significantly different from 100.
          out_buffer.resize(h);
          for (int j = 0; j < h; j++) {
            out_buffer[j].d = u_intrp[j];
          }
          printNumberSeries(&data_foutp, out_buffer, per_line_count, content_format_length + 1,
                            decimal_places, NumberFormat::STANDARD_REAL, "printToFile",
                            "Print a " + getEnumerationName(bgmr.field) + " mesh in " +
                            getEnumerationName(syntax) + " format.");
          lines_buffered = 0;
          h = 0;
        }
      }

      // Close the output data file.  The named output file will reference this data in order to
      // display the mesh in the context of the underlying structure.
      data_foutp.close();

      // Data is also written to the matrix file in the form of the molecular structure.  Draw the
      // molecule as a 3D line snaking through all of its various atoms.  The line will backtrack
      // on itself.
      renderMolecule(&foutp, bgm.getCoordinatePointer(), bgm.getTopologyPointer(),
                     ropt.getSelfPointer(), syntax);
      if (ropt.displayFieldBorders()) {
        foutp.write("\n", 1);
        drawFieldBorders(&foutp, mps, ropt, syntax);
      }
    }
    break;
  case GridFileSyntax::OPEN_DX:
  case GridFileSyntax::CUBEGEN:
    for (int i = 0; i < output_dims.na; i++) {
      const double di = i;
      for (int j = 0; j < output_dims.nb; j++) {
        const double dj = j;
        for (int k = 0; k < output_dims.nc; k++) {
          const double dk = k;
          xpt_out[k] = (output_dims.invu[0] * di) + (output_dims.invu[3] * dj) +
                       (output_dims.invu[6] * dk) + output_orig_x;
          ypt_out[k] = (output_dims.invu[4] * dj) + (output_dims.invu[7] * dk) + output_orig_y;
          zpt_out[k] = (output_dims.invu[8] * dk) + output_orig_z;
        }
        interpolate<T, double, double, double>(bgmr, mnbk, prop_a.data(), prop_b.data(),
                                               index_zero.data(), xpt_out.data(), ypt_out.data(),
                                               zpt_out.data(), nullptr, nullptr, nullptr,
                                               output_dims.nc, u_intrp.data(), nullptr, nullptr,
                                               nullptr, nullptr, nullptr, nullptr,
                                               OffMeshProtocol::EXTRAPOLATE, 0, 0, 0);
        out_buffer.resize(output_dims.nc);
        for (int k = 0; k < output_dims.nc; k++) {
          out_buffer[k].d = u_intrp[k];
        }
        if (syntax == GridFileSyntax::OPEN_DX) {
          printNumberSeries(&foutp, out_buffer, 3, content_format_length + 1, decimal_places,
                            NumberFormat::STANDARD_REAL, "printToFile", "Print a " +
                            getEnumerationName(bgmr.field) + " mesh in " +
                            getEnumerationName(syntax) + " format.");
        }
        else {
          printNumberSeries(&foutp, out_buffer, 6, decimal_places + 8, decimal_places,
                            NumberFormat::SCIENTIFIC, "printToFile", "Print a " +
                            getEnumerationName(bgmr.field) + " mesh in " +
                            getEnumerationName(syntax) + " format.");
        }
      }
    }
    break;
  }
  foutp.write(post_script.data(), post_script.size());

  // Close the output file
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void printToFile(const BackgroundMesh<T> &bgm, const std::string &file_name,
                 const PrintSituation expectation, const GridFileSyntax syntax,
                 const RenderOptions &ropt, const int decimal_places,
                 const ExceptionResponse policy) {
  printToFile(bgm, bgm.getDimensions(), file_name, expectation, syntax, ropt, decimal_places,
              policy);
}

} // namespace review
} // namespace stormm
