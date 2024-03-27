// -*-c++-*-
#include "copyright.h"

namespace stormm {
namespace trajectory {

//-------------------------------------------------------------------------------------------------
template <typename T>
void writeFrame(std::ofstream *foutp, const TextFile &tf, const CoordinateSeries<T> &cs,
                const CoordinateSwapPlan &excision) {
  const int nframe = cs.getFrameCount();
  CoordinateFrame stage(cs.getAtomCount());
  CoordinateFrameWriter stagew = stage.data();
  for (int i = 0; i < nframe; i++) {
    cs.extractFrame(&stage, i);
    writeFrame(foutp, tf, stagew.xcrd, stagew.ycrd, stagew.zcrd, excision);
  }
}

//-------------------------------------------------------------------------------------------------
template <typename T>
void writeFrame(const std::string &filename, PrintSituation expectation, const TextFile &tf,
                const CoordinateSeries<T> &cs, const CoordinateSwapPlan &excision) {
  std::ofstream foutp = openOutputFile(filename, expectation, "Open an archive for writing all "
                                       "frames of a CoordinateSeries in an annotated format.");
  writeFrame(&foutp, tf, cs, excision);
  foutp.close();
}

} // namespace trajectory
} // namespace stormm
