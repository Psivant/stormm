// -*-c++-*-
#include "copyright.h"
#include "Accelerator/gpu_details.h"
#include "Math/rounding.h"
#include "MolecularMechanics/dynamics_intervention.h"
#include "Parsing/parse.h"
#include "Parsing/parsing_enumerators.h"
#include "Potential/cellgrid.h"
#include "Reporting/error_format.h"
#include "Synthesis/atomgraph_synthesis.h"
#include "Synthesis/phasespace_synthesis.h"
#include "Trajectory/coordinate_copy.h"
#include "Trajectory/trajectory_enumerators.h"
#include "debug_nl_util.h"

namespace stormm {
namespace mm {

using debug::cellAssignmentFailure;
using energy::CellGrid;
using energy::CellGridReader;
using energy::CellGridWriter;
using parse::NumberFormat;
using parse::realToString;
using synthesis::AtomGraphSynthesis;
using synthesis::PhaseSpaceSynthesis;
using trajectory::coordCopy;
using trajectory::CoordinateCycle;

//-------------------------------------------------------------------------------------------------
template <typename Tcoord, typename Tacc, typename Tcalc, typename Tcoord4>
void debugNLComposition(const int step, const NonbondedTheme theme, const GpuDetails &gpu) {

  // By this point, the layout of the neighbor list has been verified, with one neighbor list
  // theme chosen.  The neighbor list workspace and download will be chosen depending on the theme.
  CellGrid<Tcoord, Tacc,
           Tcalc, Tcoord4> *test_cg = dyna_tk.getNLWorkspacePointer<Tcoord, Tacc,
                                                                    Tcalc, Tcoord4>(theme);
  const CellGrid<Tcoord, Tacc,
                 Tcalc, Tcoord4> *base_cg = dyna_tk.getNeighborListPointer<Tcoord, Tacc,
                                                                           Tcalc, Tcoord4>(theme);
  test_cg->updateCyclePosition(base_cg->getCyclePosition());

  const PhaseSpaceSynthesis *poly_ps = dyna_tk.getPhaseSpaceSynthesisPointer();
  PhaseSpaceSynthesis *crd_wkspc = dyna_tk.getWorkspacePointer(0);
  crd_wkspc->updateCyclePosition(poly_ps->getCyclePosition());
  const int nsys = poly_ps->getSystemCount();
  Hybrid<int2> system_pairs(nsys, "synth_xfer_pairs", HybridFormat::HOST_MOUNTED);
  for (int i = 0; i < nsys; i++) {
    system_pairs.putHost({i, i}, i);
  }
  coordCopy(crd_wkspc, *poly_ps, system_pairs, HybridTargetLevel::HOST, HybridTargetLevel::DEVICE,
            gpu);
  test_cg->populateImage(CoordinateCycle::WHITE);
  test_cg->populateImage(CoordinateCycle::BLACK);
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  const CellGridReader<Tcoord, Tacc,
                       Tcalc, Tcoord4> base_cgr = base_cg->data(devc_tier);
  CellGridWriter<Tcoord, Tacc, Tcalc, Tcoord4> test_cgw = test_cg->data(HybridTargetLevel::HOST);
  
  // Confirm the consistency of each neighbor list
  if (base_cgr.system_count != test_cgw.system_count) {
    rtErr("The number of systems in the workspace (" + std::to_string(test_cgw.system_count) +
          ") is inconsistent with that in the original neighbor list (" +
          std::to_string(base_cgr.system_count) + ").", "debugNLComposition");
  }
  if (base_cgr.total_cell_count != test_cgw.total_cell_count) {
    rtErr("The number of decomposition cells in the workspace (" +
          std::to_string(test_cgw.total_cell_count) + ") is inconsistent with that in the "
          "original neighbor list (" + std::to_string(base_cgr.total_cell_count) + ").",
          "debugNLComposition");
  }
  if (base_cgr.cell_base_capacity != test_cgw.cell_base_capacity) {
    rtErr("The base cell capacity in the workspace (" +
          std::to_string(test_cgw.cell_base_capacity) + ") is inconsistent with that in the "
          "original neighbor list (" + std::to_string(base_cgr.cell_base_capacity) + ").",
          "debugNLComposition");
  }
  
  // Check the populations of each cell.
  const size_t total_image_volume = base_cgr.total_cell_count * base_cgr.cell_base_capacity;
  const std::vector<ullint> sim_cell_grids = deepCopy<ullint>(base_cgr.system_cell_grids,
                                                              base_cgr.system_count, devc_tier,
                                                              "copy neighbor list cell grid "
                                                              "dimensions for NL debugging");
  const std::vector<uint2> sim_cell_limits = deepCopy<uint2>(base_cgr.cell_limits,
                                                             base_cgr.total_cell_count, devc_tier,
                                                             "copy cell limits for NL debugging");
  const std::vector<Tcoord4> sim_img = deepCopy<Tcoord4>(base_cgr.image, total_image_volume,
                                                         devc_tier,
                                                         "copy image contents for NL debugging");
  const std::vector<Tcalc> sim_cell_umat = deepCopy<Tcalc>(base_cgr.system_cell_umat,
                                                           base_cgr.system_count * warp_size_int,
                                                           devc_tier, "copy neighbor list "
                                                           "cell-specific transform matrices");
  const std::vector<int> sim_nonimg_atom_idx = deepCopy<int>(base_cgr.nonimg_atom_idx,
                                                             total_image_volume, devc_tier,
                                                             "copy non-image atom indices for NL "
                                                             "debugging");
  
  struct NLNote {
    uint test_image_index;  // Index of the noteworthy atom within the test neighbor list created
                            //   by initialization from the current state of the simulation.  This
                            //   refers to coordinates downloaded into the staging
                            //   PhaseSpaceSynthesis, attached to dyna_tk and accessed here by the
                            //   pointer crd_wrkspc) 
    uint base_image_index;  // Index of the noteworthy atom within the original neighbor list as
                            //   utilized by the simulation
    int synthesis_index;    // Synthesis index of the atom, which will be consistent among both
                            //   neighbor lists.  Neighbor list notes are constructed by searching
                            //   for atoms of the same synthesis index in one list or the other.
    int system_index;       // Index of the atom within one of the systems in the synthesis
    int test_nl_cell;       // Cell index of the noteworthy atom in the test neighbor list
    int base_nl_cell;       // Cell index of the noteworthy atom in the original neighbor list
    int sys_test_cell_a;    // The position, along the system's unit cell A axis, at which the
                            //   noteworthy atom is found in the test neighbor list
    int sys_test_cell_b;    // Cell position along the unit cell B axis in the test neighbor list
    int sys_test_cell_c;    // Cell position along the unit cell C axis in the test neighbor list
    int sys_base_cell_a;    // Cell position along the A axis in the original neighbor list
    int sys_base_cell_b;    // Cell position along the B axis in the original neighbor list
    int sys_base_cell_c;    // Cell position along the C axis in the original neighbor list
    double base_fa;         // Fractional coordinates of the atom, within the neighbor list cell
                            //   of the original neighbor list, along the unit cell A axis
    double base_fb;         // Fractional coordinates of the atom along the unit cell B axis
    double base_fc;         // Fractional coordinates of the atom along the unit cell C axis
    double test_fa;         // Fractional coordinates of the atom, within the neighbor list cell
                            //   of the original neighbor list, along the unit cell A axis
    double test_fb;         // Fractional coordinates of the atom along the unit cell B axis
    double test_fc;         // Fractional coordinates of the atom along the unit cell C axis
  };
  
  // Check for possible problems within each cell, based on inconsistencies between the CPU-based
  // initialization of the neighbor list cells from the current coordinates and the simulation's
  // state of the cell grid.  Because the CPU initialization is done in double-precision while the
  // GPU may be calculating in some other mode, false positives may arise in which the two cell
  // grids display inconsistent populations from cell to cell but the discrepancy is merely an
  // atom near the cell boundary.  Filter the false positives but report anything which is not so
  // easily dismissed.
  std::vector<NLNote> population_discrepancies;
  for (int sysidx = 0; sysidx < test_cgw.system_count; sysidx++) {
    const ullint icg_test_dims = test_cgw.system_cell_grids[sysidx];
    const int sys_st = (icg_test_dims & 0xfffffff);
    const int sys_na = ((icg_test_dims >> 28) & 0xfff);
    const int sys_nb = ((icg_test_dims >> 40) & 0xfff);
    const int sys_nc = ((icg_test_dims >> 52) & 0xfff);
    const ullint icg_base_dims = sim_cell_grids[sysidx];
    const int chk_sys_st = (icg_test_dims & 0xfffffff);
    const int chk_sys_na = ((icg_test_dims >> 28) & 0xfff);
    const int chk_sys_nb = ((icg_test_dims >> 40) & 0xfff);
    const int chk_sys_nc = ((icg_test_dims >> 52) & 0xfff);
    if (sys_st != chk_sys_st ||
        sys_na != chk_sys_na || sys_nb != chk_sys_nb || sys_nc != chk_sys_nc) {
      rtErr("The cell grid layout for system " + std::to_string(sysidx) + " is inconsistent in "
            "the simulation and workspace neighbor lists.\n  Simulation: [ init " +
            std::to_string(chk_sys_st) + "   grid " + std::to_string(chk_sys_na) + " x " +
            std::to_string(chk_sys_nb) + " x " + std::to_string(chk_sys_nc) + " ]\n  "
            "Workspace:  [ init " + std::to_string(sys_st) + "   grid " + std::to_string(sys_na) +
            " x " + std::to_string(sys_nb) + " x " + std::to_string(sys_nc) + " ]",
            "debugNLComposition");
    }
    const size_t xfrm_start = sysidx * roundUp(9, warp_size_int);
    for (size_t i = 0; i < 9; i++) {
      const Tcalc dmt = sim_cell_umat[xfrm_start + i] - test_cgw.system_cell_umat[xfrm_start + i];
      if (fabs(dmt) > 1.0e-5) {
        rtErr("Transform element " + std::to_string(i) + " of system " + std::to_string(sysidx) +
              " differs by " + realToString(dmt, 3, NumberFormat::SCIENTIFIC) + " in the "
              "original and workspace neighbor lists.", "debugNLComposition");
      }
    }
    for (int j = 0; j < sys_nb; j++) {
      for (int k = 0; k < sys_nc; k++) {

        // Check the chain.  If one cell has an errant atom count due to roundoff error placing
        // the atom in an adjacent cell, that discrepancy will propagate down the chain.  Such a
        // result must be anticipated as it is benign, so that true errors in the cell population
        // can rise to the developers' attention.  Each benign error will end up being verified
        // twice, to reduce the overhead of tracking but they will be relatively rare.
        for (int i = 0; i < sys_na; i++) {
          const int cell_ijk = sys_st + (((k * sys_nb) + j) * sys_na) + i;
          const uint2 icell_gpu_lims = sim_cell_limits[cell_ijk];

          // Find the list synthesis atom indices which are in the GPU's neighbor list cell but
          // not in the CPU neighbor list cell, and vice-versa.
          const uint base_llim = sim_cell_limits[cell_ijk].x;
          const uint base_hlim = base_llim + (sim_cell_limits[cell_ijk].y >> 16);
          for (int mb = base_llim; mb < base_hlim; mb++) {
            const int synth_idx = sim_nonimg_atom_idx[mb];
            const uint mx = test_cgw.img_atom_idx[synth_idx];
            NLNote tnote;
            tnote.test_image_index = mx;
            tnote.base_image_index = mb;
            tnote.synthesis_index = synth_idx;
            tnote.system_index = (icell_gpu_lims.y & 0xffff);
            tnote.base_nl_cell = cell_ijk;

            // Check first whether the atom is in the same cell of the test neighbor list.  If not,
            // look around.  If that still doesn't work, do an exhaustive search to find the actual
            // cell that contains image atom mx.
            if (mx >= test_cgw.cell_limits[cell_ijk].x &&
                mx < test_cgw.cell_limits[cell_ijk].x + (test_cgw.cell_limits[cell_ijk].y >> 16)) {
              tnote.test_nl_cell = cell_ijk;
            }
            else {
              const uint sys_chain_len = sys_na * test_cgw.cell_base_capacity;
              const uint chn_bc = (mx - (sys_st * test_cgw.cell_base_capacity)) / sys_chain_len;
              bool found = false;
              for (int aidx = 0; aidx < sys_na; aidx++) {
                const int cta = sys_st + (chn_bc * sys_na) + aidx;
                if (mx >= test_cgw.cell_limits[cta].x &&
                    mx < test_cgw.cell_limits[cta].x + (test_cgw.cell_limits[cta].y >> 16)) {
                  tnote.test_nl_cell = cta;
                  found = true;
                }
              }
              if (found == false) {
                rtErr("Unable to find synthesis atom " + std::to_string(synth_idx) + " in the "
                      "rebuilt neighbor list.", "debugNLComposition");
              }
            }
            if (tnote.test_nl_cell != tnote.base_nl_cell) {

              // If an inconsistency in the cell locations is detected, finish filling out the
              // entry for later inspection.
              const int test_cidx = tnote.test_nl_cell - sys_st;
              tnote.sys_test_cell_c = test_cidx / (sys_na * sys_nb);
              tnote.sys_test_cell_b = (test_cidx - (tnote.sys_test_cell_c * sys_na * sys_nb)) /
                                      sys_na;
              tnote.sys_test_cell_a = test_cidx - (tnote.sys_test_cell_c * sys_na * sys_nb) -
                                      (tnote.sys_test_cell_b * sys_na);
              const int base_cidx = tnote.base_nl_cell - sys_st;
              tnote.sys_base_cell_c = base_cidx / (sys_na * sys_nb);
              tnote.sys_base_cell_b = (base_cidx - (tnote.sys_base_cell_c * sys_na * sys_nb)) /
                                      sys_na;
              tnote.sys_base_cell_a = base_cidx - (tnote.sys_base_cell_c * sys_na * sys_nb) -
                                      (tnote.sys_base_cell_b * sys_na);
              const Tcoord4 base_crdq = sim_img[mb];
              const Tcoord4 test_crdq = test_cgw.image[mx];
              tnote.base_fa = (sim_cell_umat[xfrm_start    ] * base_crdq.x) +
                              (sim_cell_umat[xfrm_start + 3] * base_crdq.y) +
                              (sim_cell_umat[xfrm_start + 6] * base_crdq.z);
              tnote.base_fb = (sim_cell_umat[xfrm_start + 4] * base_crdq.y) +
                              (sim_cell_umat[xfrm_start + 7] * base_crdq.z);
              tnote.base_fc = (sim_cell_umat[xfrm_start + 8] * base_crdq.z);
              tnote.test_fa = (test_cgw.system_cell_umat[xfrm_start    ] * test_crdq.x) +
                              (test_cgw.system_cell_umat[xfrm_start + 3] * test_crdq.y) +
                              (test_cgw.system_cell_umat[xfrm_start + 6] * test_crdq.z);
              tnote.test_fb = (test_cgw.system_cell_umat[xfrm_start + 4] * test_crdq.y) +
                              (test_cgw.system_cell_umat[xfrm_start + 7] * test_crdq.z);
              tnote.test_fc = (test_cgw.system_cell_umat[xfrm_start + 8] * test_crdq.z);
              if (cellAssignmentFailure(tnote.sys_test_cell_a, tnote.sys_base_cell_a, sys_na,
                                        tnote.test_fa, tnote.base_fa) ||
                  cellAssignmentFailure(tnote.sys_test_cell_b, tnote.sys_base_cell_b, sys_nb,
                                        tnote.test_fb, tnote.base_fb) ||
                  cellAssignmentFailure(tnote.sys_test_cell_c, tnote.sys_base_cell_c, sys_nc,
                                        tnote.test_fc, tnote.base_fc)) {
                printf("  Atom %6d  System %d  Orig %2d %2d %2d  %9.6lf %9.6lf %9.6lf :: Test %2d "
                       "%2d %2d  %9.6lf %9.6lf %9.6lf\n", tnote.synthesis_index,
                       tnote.system_index, tnote.sys_base_cell_a, tnote.sys_base_cell_b,
                       tnote.sys_base_cell_c, tnote.base_fa, tnote.base_fb, tnote.base_fc,
                       tnote.sys_test_cell_a, tnote.sys_test_cell_b, tnote.sys_test_cell_c,
                       tnote.test_fa, tnote.test_fb, tnote.test_fc);
                population_discrepancies.push_back(tnote);
              }
            }
          }

          // Examine the outliers.  The list of atoms found in any given neighbor list cell of the
          // CPU-based test apparatus but not found in the same neighbor list cell of the
          // simulation should contain the same synthesis
        }
      }
    }
  }
  if (population_discrepancies.size() > 0) {
    printf("A total of %zu discrepancies were found on step %4d.\n",
           population_discrepancies.size(), step);
  }
}

} // namespace mm
} // namespace stormm
