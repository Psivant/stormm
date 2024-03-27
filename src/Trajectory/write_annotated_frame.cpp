#include "copyright.h"
#include "write_annotated_frame.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
void writeFrame(std::ofstream *foutp, const TextFile &tf) {
  tf.write(foutp);
}

//-------------------------------------------------------------------------------------------------
void writeFrame(const std::string &filename, const PrintSituation expectation,
                const TextFile &tf) {
  tf.write(filename, expectation);
}

//-------------------------------------------------------------------------------------------------
void writeFrame(std::ofstream *foutp, const TextFile &tf, const double* xcrd, const double* ycrd,
                const double* zcrd, const CoordinateSwapPlan &excision) {

}

//-------------------------------------------------------------------------------------------------
void writeFrame(const std::string &filename, const PrintSituation expectation, const TextFile &tf,
                const double* xcrd, const double* ycrd, const double* zcrd,
                const CoordinateSwapPlan &excision) {
  std::ofstream foutp = openOutputFile(filename, expectation, "Open a coordinate file write "
                                       "annotated coordinates.");
  writeFrame(&foutp, tf, xcrd, ycrd, zcrd, excision);
  foutp.close();
}

//-------------------------------------------------------------------------------------------------
void writeFrame(const std::string &filename, const PrintSituation expectation, const TextFile &tf,
                const PhaseSpaceReader &psr, const CoordinateSwapPlan &excision) {
  writeFrame(filename, expectation, tf, psr.xcrd, psr.ycrd, psr.zcrd, excision);
}

//-------------------------------------------------------------------------------------------------
void writeFrame(const std::string &filename, const PrintSituation expectation, const TextFile &tf,
                const PhaseSpace &ps, const CoordinateSwapPlan &excision) {
  const PhaseSpaceReader psr = ps.data(ps.getCyclePosition());
  writeFrame(filename, expectation, tf, psr, excision);
}

//-------------------------------------------------------------------------------------------------
void writeFrame(const std::string &filename, const PrintSituation expectation, const TextFile &tf,
                const PhaseSpace &ps, const CoordinateSwapPlan &excision,
                const CoordinateCycle time_point) {
  const PhaseSpaceReader psr = ps.data(time_point);
  writeFrame(filename, expectation, tf, psr, excision);
}

//-------------------------------------------------------------------------------------------------
void writeFrame(const std::string &filename, const PrintSituation expectation, const TextFile &tf,
                const CoordinateFrameReader &cfr, const CoordinateSwapPlan &excision) {
  writeFrame(filename, expectation, tf, cfr.xcrd, cfr.ycrd, cfr.zcrd, excision);
}

//-------------------------------------------------------------------------------------------------
void writeFrame(const std::string &filename, const PrintSituation expectation, const TextFile &tf,
                const CoordinateFrame &cf, const CoordinateSwapPlan &excision) {
  const CoordinateFrameReader cfr = cf.data();
  writeFrame(filename, expectation, tf, cfr, excision);
}

//-------------------------------------------------------------------------------------------------
void writeFrame(std::ofstream *foutp, const std::vector<TextFile> &tf_list,
                const PsSynthesisReader &poly_psr,
                const std::vector<CoordinateSwapPlan> &excision_list,
                const std::vector<int> &plan_indices) {
  for (int i = 0; i < poly_psr.system_count; i++) {
    const int natom  = poly_psr.atom_counts[i];
    const int offset = poly_psr.atom_starts[i];
    std::vector<double> xcrd(natom), ycrd(natom), zcrd(natom);
    for (int i = 0; i < natom; i++) {
      xcrd[i] = static_cast<double>(poly_psr.xcrd[offset + i]) * poly_psr.inv_gpos_scale;
      ycrd[i] = static_cast<double>(poly_psr.ycrd[offset + i]) * poly_psr.inv_gpos_scale;
      zcrd[i] = static_cast<double>(poly_psr.zcrd[offset + i]) * poly_psr.inv_gpos_scale;
    }
    writeFrame(foutp, tf_list[plan_indices[i]], xcrd.data(), ycrd.data(),
               zcrd.data(), excision_list[plan_indices[i]]);
  }
}

//-------------------------------------------------------------------------------------------------
void writeFrame(const std::string &filename, const PrintSituation expectation,
                const std::vector<TextFile> &tf_list, const PhaseSpaceSynthesis &poly_ps,
                const std::vector<CoordinateSwapPlan> &excision_list,
                const std::vector<int> &plan_indices) {
  std::ofstream foutp = openOutputFile(filename, expectation, "Open an archive for writing "
                                       "PhaseSpaceSynthesis coordinates from a variety of "
                                       "systems.");
  writeFrame(&foutp, tf_list, poly_ps.data(), excision_list, plan_indices);
}

//-------------------------------------------------------------------------------------------------
void writeFrame(const std::string &filename, const PrintSituation expectation,
                const std::vector<TextFile> &tf_list, const PhaseSpaceSynthesis &poly_ps,
                const std::vector<CoordinateSwapPlan> &excision_list,
                const std::vector<int> &plan_indices, const CoordinateCycle time_point) {
  std::ofstream foutp = openOutputFile(filename, expectation, "Open an SD file to archive "
                                       "PhaseSpaceSynthesis coordinates from a variety of "
                                       "systems.");
  writeFrame(&foutp, tf_list, poly_ps.data(time_point), excision_list, plan_indices);
}

} // namespace trajectory
} // namespace stormm
