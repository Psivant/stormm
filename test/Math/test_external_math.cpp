// -*-c++-*-
#include <cmath>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if STORMM_INCLUDE_NETCDF
#  include <netcdf.h>
#endif
#  include "pocketfft_hdronly.h"
#ifdef STORMM_USE_HPC
#  ifdef STORMM_USE_CUDA
#    include <cufft.h>
#  endif
#endif
#include <string.h>
#include "copyright.h"
#include "../../src/Constants/behavior.h"
#include "../../src/Constants/hpc_bounds.h"
#include "../../src/DataTypes/common_types.h"
#include "../../src/DataTypes/stormm_vector_types.h"
#include "../../src/Accelerator/hpc_config.h"
#include "../../src/Accelerator/hybrid.h"
#include "../../src/FileManagement/file_listing.h"
#include "../../src/Math/bspline.h"
#include "../../src/Math/fft_stage.h"
#include "../../src/Math/summation.h"
#include "../../src/Numerics/split_fixed_precision.h"
#include "../../src/Reporting/error_format.h"
#include "../../src/Random/random.h"
#include "../../src/Random/hpc_random.h"
#include "../../src/UnitTesting/unit_test.h"

using namespace stormm::data_types;
using namespace stormm::random;
using namespace stormm::review;
using namespace stormm::stmath;
using namespace stormm::testing;
using namespace pocketfft;

const std::vector<double> five_six_eight = {
    0.23561365,  -1.48268197,   0.22822286,   1.71344362,   0.33405646,  -2.01013931,
    0.35689393,  -1.14162142,  -1.44353881,  -0.70749066,  -1.28968792,   0.01851056,
    0.60666665,  -2.16329268,  -0.62773473,   1.02803268,   0.09074048,   0.62099517,
   -1.61119606,  -1.09926623,   0.42645758,  -1.00952006,   1.68298606,  -0.15159400,
    0.46751418,   0.41377707,  -0.44422088,   1.64234361,   0.10727689,  -0.34424157,
    1.36103807,   0.77701904,   0.03351500,  -0.43374849,  -0.32534136,   0.72404286,
   -0.16036023,   0.41833060,  -0.47646626,   1.45103153,   0.32335401,   0.16063269,
   -0.83933848,  -0.32857699,   1.28150914,  -1.49647110,  -0.65922042,   0.36632144,
   -0.65587288,   1.85606646,  -0.90628305,   0.01471202,  -0.57748024,  -0.12842493,
    1.13404761,  -0.95068390,  -0.09574151,   0.26724517,   0.09267206,  -0.60516867,
    1.91574366,   1.08291777,   0.09192416,   0.14352491,  -0.50932370,  -1.83282802,
   -0.65296841,  -0.56653663,   0.21985596,  -0.08819616,   0.95827239,  -0.18231796,
   -1.33640166,   0.80625718,   0.53020014,  -1.64453823,  -0.17422112,  -1.92345358,
   -0.24957657,   1.16931836,  -1.00701365,  -0.26617480,   1.35801174,  -0.55614727,
    0.65291326,   0.72109507,  -1.32034236,   1.52397440,   0.24137837,  -0.06796653,
   -0.69681188,  -1.30040716,  -0.50684943,   0.73373434,   0.40005207,   1.27517710,
    0.44934973,  -0.27877955,   0.21305959,  -0.15175912,  -1.23950211,   0.37515122,
   -0.52296297,   1.22406115,  -0.78918458,   0.93116711,   0.39296802,   1.71805306,
   -0.11809381,  -0.81139642,  -1.05623996,  -0.69043243,   0.54746314,  -0.20811051,
    1.16362095,  -0.29200480,  -0.07729086,   0.87994897,  -0.40483997,  -0.66788911,
   -1.48633626,   0.03692348,   0.59140726,   0.78573582,  -0.41206623,   0.59786897,
   -0.83734444,   0.66895850,   0.55654466,  -1.09368130,  -0.54147987,  -1.32213020,
   -1.19319530,   1.34880462,  -0.23731254,  -0.11757148,   1.82389741,   0.00331725,
   -0.32555577,   2.09311919,   1.44359775,   1.11218297,   2.06375284,  -2.30544081,
   -0.24425427,  -2.65119112,   1.29463193,  -0.28213662,  -0.48051671,   2.10417064,
    0.50621488,   0.94443114,  -0.97810234,  -0.88811855,  -2.06844300,  -0.19114342,
   -0.16810124,  -1.57142322,   0.82480635,  -1.41920836,  -0.53484718,  -0.86952860,
    0.37496199,   0.24968398,  -0.31957573,  -0.46835642,   0.23541842,  -0.41525381,
   -0.63075030,   1.41821369,  -0.32664792,  -0.48291617,  -0.08603563,   0.90363597,
    0.51925941,   1.43802401,   1.42389384,   0.05597802,   1.41618107,  -0.09471747,
   -1.76181323,  -0.92103928,  -0.65157751,  -0.03881014,  -0.76551793,  -0.76492448,
   -0.78990745,   0.11829086,  -0.36993883,  -1.45492164,   0.18712974,   0.99856928,
    1.36749131,   1.32316703,  -1.25382806,   1.26913056,   0.43826269,  -1.34703531,
   -1.07537209,  -0.03493886,  -1.04767856,   0.80075685,  -0.89459365,   1.21728151,
   -0.92472402,   0.48265555,   0.92885820,  -1.66164597,  -0.03359298,   0.87621395,
    1.06758372,   0.46138843,  -1.69091690,   0.39565109,   0.02670228,  -0.14086059,
   -1.26511792,   0.45415082,  -0.33273476,   1.86455964,  -1.61003522,   0.98634721,
   -0.37366074,  -0.58006470,  -0.05685893,   1.14383989,  -1.89326922,  -0.11051784,
    0.25556975,  -1.40434275,  -0.52397432,   1.58208146,   1.87196696,   2.83300654,
    1.03893543,  -0.17444100,  -1.86539154,   1.05696610,  -1.25215311,   0.80210980
};  

const std::vector<double> sixteen_four_seven = {
    0.72447985,   1.75201878,  -0.28372798,   0.15994284,  -1.63646208,   0.38657712,
    1.53166195,  -0.24417209,  -0.10943985,  -0.01502616,  -1.49830529,   0.27024711,
   -1.71656603,  -1.52624478,  -1.45445336,  -0.15037633,   1.10651483,  -0.25324536,
    0.79475156,  -0.52273076,   2.03474419,   1.90662492,  -0.67911477,  -0.38971166,
   -0.07012029,   0.97608999,  -0.63658922,   0.59607114,   2.30339013,   0.95820941,
    0.01612747,   0.53424655,   1.69728927,   0.21652037,   0.97245801,  -0.84793122,
    1.20035205,   0.30417433,  -0.13270383,   0.64001461,  -0.81793756,  -1.77866677,
    1.42672641,  -2.18398261,  -0.37492039,   0.05432246,  -0.86920731,  -0.05209435,
    0.98295108,  -1.26913829,  -1.10281990,  -1.84312506,  -0.47332580,  -0.41602871,
    0.27561692,  -0.22870741,  -0.20712625,   0.12192829,  -0.26890303,  -0.27414726,
    2.02594129,  -1.73249840,   0.19429491,  -0.37541359,   2.01826211,  -0.56217066,
   -0.18206370,   0.45947102,   0.58768167,   0.64818320,   0.70460887,   1.35105322,
   -1.39780496,  -0.91574346,  -0.07563850,   0.98175196,  -0.16785577,  -0.39397784,
   -2.17311330,  -0.12688891,   0.30843970,  -0.81593206,   1.24123646,   0.07037778,
   -0.72370251,   1.16178660,   0.40170143,   0.59887245,  -0.38396047,   0.20518420,
    0.67673330,   1.51132230,  -0.17922098,  -0.05195773,   1.44268265,   1.91543827,
    1.52273433,   0.93424092,   0.09244797,  -0.39644335,   0.74640335,  -0.66154297,
   -1.05565615,  -0.66076470,  -3.18065231,   0.75115428,   1.22586634,  -0.08563510,
    0.13733619,   0.39219665,   1.22982309,   0.83005928,   1.97046692,   0.20397678,
    0.57333862,   0.10873303,   0.57919968,  -1.72462149,  -1.60115213,   1.92016957,
    1.89439511,   2.20925000,  -0.60809349,  -0.69287241,  -0.27139947,  -1.10912926,
    0.22944764,   0.89110937,  -0.39001964,   0.77264610,  -1.21091252,   0.33812431,
    0.66584047,  -0.82729876,  -0.29608153,   1.81010517,  -1.28204858,   0.65939148,
    0.26622133,  -0.25956688,   0.50648723,   0.09533395,   1.32591719,   0.02163427,
   -1.19246873,   1.22251327,  -0.32528721,   0.89350946,  -1.10728892,  -1.55959689,
   -0.31165747,   0.26085593,   0.71968584,  -0.31877247,  -1.30050486,  -0.08340267,
    1.41931630,   0.23020373,   0.31089505,  -1.99671231,   0.33219541,  -1.46393963,
   -0.58839891,   0.09578559,   0.27072256,   0.91284770,  -0.79320375,   0.65049550,
    0.53837662,   0.37574620,   1.68984657,  -1.82814879,  -0.75652914,   0.00981733,
   -0.40404023,   2.21449929,  -0.41822246,   0.08578655,  -0.87427102,  -0.62533607,
   -0.21980852,   1.29298122,   0.38764863,  -0.82527197,   0.01683894,   0.02876320,
    0.25535238,  -1.81301094,   0.55787262,   1.24042320,   1.04001436,  -1.50149265,
    0.67467053,   0.25155575,  -0.25157681,   0.16060238,   0.84405965,   0.65209515,
    0.67534233,  -0.86067756,  -1.21556516,   0.60529846,  -0.32894262,   0.20766376,
   -1.06693550,  -0.75505938,   2.32793832,  -0.38091714,  -0.79153448,  -0.84853009,
   -0.87273136,   0.30826863,   0.26594882,  -0.83534240,   0.32009972,   1.29490407,
   -0.27843445,   0.37601008,  -0.68400667,  -0.31511294,   0.58478386,  -0.50255175,
   -0.02914289,   1.43252934,   0.71570464,   1.14021715,  -0.17668914,   2.06216699,
   -1.61268802,  -0.61144504,   1.68459322,  -0.22269785,  -0.48889279,   0.18281786,
    0.09590751,   0.25469778,   1.58956194,   0.46988526,   1.74265469,  -0.13451233,
    0.18551447,   0.23530844,   1.42119214,  -0.21092768,   0.79747941,   0.82961449,
    0.24872532,  -0.18357345,  -0.67349570,  -0.07903986,   1.22736061,   1.64925820,
    0.57169080,  -0.57044717,   0.17983231,  -0.83859338,   0.11606877,  -0.90705341,
    0.94977519,   1.81966463,  -0.85208865,  -2.29715927,   1.32084593,   0.11374786,
   -0.06340263,  -0.02009809,   0.52486086,  -0.12983146,  -1.48712625,  -1.65616626,
   -0.64284505,  -1.03260698,   1.68527564,  -0.55692328,  -0.53176539,   1.08795502,
    0.09725703,  -2.17056662,   0.56072450,  -0.65863482,  -0.31919207,   0.37850798,
   -1.86876186,   0.09723926,   1.46643377,  -0.17937838,   1.55555207,   0.20827132,
    0.86484809,  -0.54879777,   0.63744192,   0.85596358,   0.72625306,  -0.00734712,
   -1.57306259,  -0.93689201,  -0.30005537,  -0.49804151,  -1.15148895,   0.33947881,
   -1.56724312,   1.90640016,  -0.28802564,   0.05563947,  -0.28629281,   0.05147013,
    0.32937343,   1.51031386,   0.26018635,  -0.28619596,  -0.52502187,   0.76570069,
   -0.57562597,  -2.26309055,  -0.91196715,  -0.33245477,  -0.79605619,  -0.01429554,
    0.34396104,   1.47077855,  -1.21636687,  -1.07224661,  -0.37091417,   0.19720375,
   -2.19180173,  -0.48715055,   1.28719204,   1.20111651,   0.78034550,  -1.26324548,
   -0.07011629,  -0.97826528,  -0.29551355,  -0.93399998,  -0.04121046,  -0.05094014,
    0.69192986,  -0.35210787,   2.46866494,  -1.76176187,   0.68726709,   0.26861531,
    0.09257712,   0.72509152,   1.94459145,   1.81921095,   1.67807604,  -0.52817126,
   -2.05024290,   1.17582932,   0.43195629,   0.33256868,   0.97680663,  -0.61975080,
   -0.52184655,   0.94608706,  -0.00546799,  -0.21750957,   2.08163374,  -1.23517192,
    0.50937603,   0.44259025,   2.19166319,  -0.74781347,  -0.19302186,  -0.34754287,
   -1.39236684,  -0.78550648,  -0.12355717,   1.04595997,  -0.12559170,  -1.83556788,
   -1.01453831,   1.33347592,  -0.81909963,   0.38424993,   0.03150246,   0.14699707,
    0.95304884,  -0.23892603,  -0.67945977,  -0.13942425,  -0.97366722,   0.34726833,
   -0.27404220,   1.79888914,   0.15013219,  -0.22463619,   1.85304396,  -1.12150106,
    1.07914011,   1.58174016,   1.06782882,  -0.44398164,  -0.68090334,  -0.20673989,
   -0.59442310,  -0.09620844,  -1.28877873,  -1.51578952,  -0.84649154,   0.93020353,
   -0.11075728,  -1.10553778,  -2.09540315,  -1.09147399,   0.30670773,  -0.87657730,
   -0.08371890,  -0.45817912,  -1.13226487,   0.98482639,  -0.47147283,  -0.06108671,
   -0.44259872,  -0.64860974,   0.43865628,  -0.32901696,   0.70009298,  -1.26368077,
    0.80404024,   0.22886170,  -0.17554467,  -1.09231060,  -2.15242264,  -1.87788989,
    0.04289880,  -0.99149711,  -0.54452874,  -0.81007972,  -1.45122202,  -0.07128785,
    0.51605788,   0.21831101,   1.08973995,   2.00893823,   0.88118737,  -0.59756638,
    0.56282122,  -1.12301670,   0.29834822,   1.94680233,  -1.33242359,   1.81700713,
   -2.20463280,  -0.57930862,  -0.16361596,   0.38315203
};

//-------------------------------------------------------------------------------------------------
// Tests the real FFT (RFFT) functionality by generating random data, performing forward and
// backward FFT operations, and comparing the result with the original data.
//
// Arguments:
//   rng:      Random number generator to use for generating input data.
//   maxlen:   The maximum length of the real FFT to test (several values in the range will be
//             tested);
//   epsilon:  Any given point in the FFT grid must return to its oiginal value after the forward
//             FFT, reverse FFT, and normalization cycle.  This is the tolerance for judging
//             whether that condition is met.
//-------------------------------------------------------------------------------------------------
void testRealFFT(Xoshiro256ppGenerator *rng, const int maxlen = 8192,
                 const double epsilon = 5.0e-10) {
  
  // Perform FFT test for different lengths of the data
  std::vector<int> test_lengths;
  int tl = 1;
  while (tl < maxlen) {
    test_lengths.push_back(tl);
    if (tl < 8) {
      tl++;
    }
    else if (tl < maxlen / 8) {
      tl = maxlen / 8;
    }
    else {
      tl += maxlen / 8;
    }
  }
  for (size_t l_idx = 0; l_idx < test_lengths.size(); l_idx++) {
    
    // Create a copy for the FFT operation
    std::vector<double> real_data = gaussianRand(rng, test_lengths[l_idx], 1.0);
    std::vector<double> original_grid_values = real_data;
    
    // Define shape and stride for the FFT
    std::vector<size_t> shape = { static_cast<size_t>(test_lengths[l_idx]) };
    std::vector<size_t> axes = { 0 };
    std::vector<ptrdiff_t> stride_in = { static_cast<ptrdiff_t>(sizeof(double)) };
    std::vector<ptrdiff_t> stride_out = { static_cast<ptrdiff_t>(sizeof(std::complex<double>)) };
    
    // Forward FFT (real to complex)
    std::vector<std::complex<double>> complex_data((test_lengths[l_idx] / 2) + 1);
    r2c(shape, stride_in, stride_out, axes, FORWARD, real_data.data(), complex_data.data(), 1.0);
    
    // Backward FFT (complex to real)
    c2r(shape, stride_out, stride_in, axes, BACKWARD, complex_data.data(), real_data.data(),
        1.0 / static_cast<double>(test_lengths[l_idx]));
    
    // Calculate error between transformed data and original data
    check(real_data, RelationalOperator::EQUAL, Approx(original_grid_values).margin(epsilon),
          "The root mean-squared error in real FFT transformation of dimension " +
          std::to_string(test_lengths[l_idx]) + " exceeds the permissible threshold.");
  }
}

//-------------------------------------------------------------------------------------------------
// Tests the complex FFT (CFFT) functionality by generating random complex data, performing forward
// and backward FFT operations, and comparing the result with the original complex data.
// Descriptions of input arguments follow from testRealFFT(), above.
//-------------------------------------------------------------------------------------------------
void testComplexFFT(Xoshiro256ppGenerator* rng, const int maxlen = 8192,
                    const double epsilon = 5.0e-10) {

  // Perform FFT test for different lengths of the data
  std::vector<int> test_lengths;
  int tl = 1;
  while (tl < maxlen) {
    test_lengths.push_back(tl);
    if (tl < 4) {
      tl++;
    }
    else if (tl < maxlen / 8) {
      tl = maxlen / 8;
    }
    else {
      tl += maxlen / 8;
    }
  }
  for (size_t l_idx = 0; l_idx < test_lengths.size(); l_idx++) {
    std::vector<double> grid_values = gaussianRand(rng, 2 * test_lengths[l_idx], 1.0);
    std::vector<double> original_grid_values = grid_values;
    
    // Create a copy for the FFT operation
    std::vector<std::complex<double>> complex_data(test_lengths[l_idx]);
    for (size_t i = 0; i < test_lengths[l_idx]; ++i) {
      complex_data[i] = std::complex<double>(grid_values[2*i], grid_values[2*i+1]);
    }
    
    // Define shape and stride for the FFT
    std::vector<size_t> shape = { static_cast<size_t>(test_lengths[l_idx]) };
    std::vector<ptrdiff_t> stride_in = { static_cast<ptrdiff_t>(sizeof(std::complex<double>)) };
    std::vector<ptrdiff_t> stride_out = { static_cast<ptrdiff_t>(sizeof(std::complex<double>)) };
    
    // Define axes for the FFT
    std::vector<size_t> axes;
    for (size_t i = 0; i < shape.size(); i++) {
      axes.push_back(i);
    }
    
    // Forward FFT
    c2c(shape, stride_in, stride_out, axes, true, complex_data.data(), complex_data.data(), 1.0);
    
    // Backward FFT
    c2c(shape, stride_in, stride_out, axes, false, complex_data.data(), complex_data.data(),
        1.0 / static_cast<double>(test_lengths[l_idx]));
    
    // Convert back to real array
    std::vector<double> fft_data(2 * test_lengths[l_idx]);
    for (size_t i = 0; i < test_lengths[l_idx]; ++i) {
      fft_data[ 2 * i     ] = complex_data[i].real();
      fft_data[(2 * i) + 1] = complex_data[i].imag();
    }
    
    // Calculate error between transformed data and original data
    check(fft_data, RelationalOperator::EQUAL, Approx(original_grid_values).margin(epsilon),
          "The root mean-squared error in complex FFT transformation of dimension " +
          std::to_string(test_lengths[l_idx]) + " exceeds the permissible threshold.");
  }
}

//-------------------------------------------------------------------------------------------------
// Test the real FFT (RFFT) functionality in a three-dimensional array by generating random data,
// performing forward and backward FFT operations, and comparing the result with the original data.
//
// Arguments:
//   rng:      Random number generator to use for generating input data.
//   shape:    Dimensions of the transform to set up and perform
//   epsilon:  Tolerance for deviations in the input and output arguments
//-------------------------------------------------------------------------------------------------
void testReal2DFFT(Xoshiro256ppGenerator *rng, const std::vector<size_t> &shape,
                   const double epsilon = 5.0e-10) {
  if (shape.size() != 2) {
    rtErr("A two-dimensional array is expected (" + std::to_string(shape.size()) +
          " dimensions were provided).", "testReal2DFFT");
  }
  if (shape[0] <= 0 || shape[1] <= 0) {
    rtErr("Positive dimensions for the grid must be specified (" + std::to_string(shape[0]) +
          " and " + std::to_string(shape[1]) + " were provided).", "testReal2DFFT");
  }
  std::vector<double> real_grid = gaussianRand(rng, shape[0] * shape[1], 1.0);
  std::vector<double> original_grid = real_grid;
  std::vector<size_t> cshape = { shape[0], (shape[1] / 2) + 1 };

  // Compute the strides
  std::vector<size_t> axes = { 0, 1 };
  std::vector<ptrdiff_t> stride_in(2);
  std::vector<ptrdiff_t> stride_out(2);
  size_t tmp_in = sizeof(double);
  size_t tmp_out = sizeof(std::complex<double>);
  for (int i = 1; i >= 0; i--) {
    stride_in[i] = tmp_in;
    tmp_in *= shape[i];
    stride_out[i] = tmp_out;
    tmp_out *= cshape[i];
  }
    
  // Forward FFT (real to complex)
  std::vector<std::complex<double>> complex_grid(cshape[0] * cshape[1]);
  for (int i = 0; i < cshape[0] * cshape[1]; i++) {
    complex_grid[i] = { 0.0, 0.0 };
  }
  r2c(shape, stride_in, stride_out, axes, FORWARD, real_grid.data(), complex_grid.data(), 1.0);
  
  // Backward FFT (complex to real)
  c2r(shape, stride_out, stride_in, axes, BACKWARD, complex_grid.data(), real_grid.data(),
      1.0 / static_cast<double>(shape[0] * shape[1]));
  
  // Calculate error between transformed data and original data
  check(real_grid, RelationalOperator::EQUAL, Approx(original_grid).margin(epsilon), "The root "
        "mean-squared error in real-to-complex FFT transformation of dimension " +
        std::to_string(shape[0]) + " x " + std::to_string(shape[1]) +
        " exceeds the permissible threshold.");
}

//-------------------------------------------------------------------------------------------------
// Test the real FFT (RFFT) functionality in a three-dimensional array by generating random data,
// performing forward and backward FFT operations, and comparing the result with the original data.
//
// Arguments:
//   rng:      Random number generator to use for generating input data.
//   shape:    Dimensions of the transform to set up and perform
//   epsilon:  Tolerance for deviations in the input and output arguments
//-------------------------------------------------------------------------------------------------
void testComplex2DFFT(Xoshiro256ppGenerator *rng, const std::vector<size_t> &shape,
                      const double epsilon = 5.0e-10) {
  if (shape.size() != 2) {
    rtErr("A two-dimensional array is expected (" + std::to_string(shape.size()) +
          " dimensions were provided).", "testComplex2DFFT");
  }
  if (shape[0] <= 0 || shape[1] <= 0) {
    rtErr("Positive dimensions for the grid must be specified (" + std::to_string(shape[0]) +
          " and " + std::to_string(shape[1]) + " were provided).", "testComplex2DFFT");
  }
  std::vector<size_t> actual_shape = { shape[1], shape[0] };
  std::vector<std::complex<double>> complex_grid(shape[0] * shape[1]);
  for (int i = 0; i < complex_grid.size(); i++) {
    complex_grid[i] = { rng->gaussianRandomNumber(), rng->gaussianRandomNumber() };
  }
  std::vector<std::complex<double>> original_grid = complex_grid;
  std::vector<size_t> cshape = actual_shape;

  // Compute the strides
  std::vector<size_t> axes = { 0, 1 };
  std::vector<ptrdiff_t> stride_in(2);
  std::vector<ptrdiff_t> stride_out(2);
  size_t tmp_in = sizeof(std::complex<double>);
  size_t tmp_out = sizeof(std::complex<double>);
  for (int i = 1; i >= 0; i--) {
    stride_in[i] = tmp_in;
    tmp_in *= actual_shape[i];
    stride_out[i] = tmp_out;
    tmp_out *= cshape[i];
  }
  
  // Forward FFT (real to complex)
  c2c(actual_shape, stride_in, stride_out, axes, FORWARD, complex_grid.data(),
      complex_grid.data(), 1.0);
  
  // Backward FFT (complex to real)
  c2c(cshape, stride_out, stride_in, axes, BACKWARD, complex_grid.data(),
      complex_grid.data(), 1.0 / static_cast<double>(shape[0] * shape[1]));
  
  // Calculate error between transformed data and original data
  std::vector<double> real_part(complex_grid.size());
  std::vector<double> imag_part(complex_grid.size());
  std::vector<double> real_init(complex_grid.size());
  std::vector<double> imag_init(complex_grid.size());
  for (int i = 0; i < complex_grid.size(); i++) {
    real_part[i] = complex_grid[i].real();
    imag_part[i] = complex_grid[i].imag();
    real_init[i] = original_grid[i].real();
    imag_init[i] = original_grid[i].imag();
  }
  check(real_part, RelationalOperator::EQUAL, Approx(real_init).margin(epsilon), "The root "
        "mean-squared error of the real part in complex FFT transformation of dimension " +
        std::to_string(shape[0]) + " x " + std::to_string(shape[1]) +
        " exceeds the permissible threshold.");
  check(imag_part, RelationalOperator::EQUAL, Approx(imag_init).margin(epsilon), "The root "
        "mean-squared error of the imaginary part in complex FFT transformation of dimension " +
        std::to_string(shape[0]) + " x " + std::to_string(shape[1]) +
        " exceeds the permissible threshold.");
}

//-------------------------------------------------------------------------------------------------
// Test the real FFT (RFFT) functionality in a three-dimensional array by generating random data,
// performing forward and backward FFT operations, and comparing the result with the original data.
//
// Arguments:
//   rng:      Random number generator to use for generating input data.
//   shape:    Dimensions of the transform to set up and perform
//   epsilon:  Tolerance for deviations in the input and output arguments
//-------------------------------------------------------------------------------------------------
void testReal3DFFT(Xoshiro256ppGenerator *rng, const std::vector<size_t> &shape,
                   const double epsilon = 5.0e-10) {
  if (shape.size() != 3) {
    rtErr("A three-dimensional array is expected (" + std::to_string(shape.size()) +
          " dimensions were provided).", "testReal3DFFT");
  }
  if (shape[0] <= 0 || shape[1] <= 0 || shape[2] <= 0) {
    rtErr("Positive dimensions for the grid must be specified (" + std::to_string(shape[0]) +
          ", " + std::to_string(shape[1]) + ", and " + std::to_string(shape[2]) +
          " were provided).", "testReal3DFFT");
  }
  std::vector<double> real_grid = gaussianRand(rng, shape[0] * shape[1] * shape[2], 1.0);
  const size_t nreal = real_grid.size();
  for (size_t i = 0; i < nreal; i++) {
    real_grid[i] = round(real_grid[i] * 1.0e8) / 1.0e8;
  }
  std::vector<double> original_grid = real_grid;
  std::vector<size_t> actual_shape = { shape[2], shape[1], shape[0] };
  std::vector<size_t> cshape = { actual_shape[0], actual_shape[1], (actual_shape[2] / 2) + 1 };

  // Compute the strides
  std::vector<size_t> axes = { 0, 1, 2 };
  std::vector<ptrdiff_t> stride_in(3);
  std::vector<ptrdiff_t> stride_out(3);
  size_t tmp_in = sizeof(double);
  size_t tmp_out = sizeof(std::complex<double>);
  for (int i = 2; i >= 0; i--) {
    stride_in[i] = tmp_in;
    tmp_in *= actual_shape[i];
    stride_out[i] = tmp_out;
    tmp_out *= cshape[i];
  }
    
  // Forward FFT (real to complex)
  std::vector<std::complex<double>> complex_grid(cshape[0] * cshape[1] * cshape[2]);
  r2c(actual_shape, stride_in, stride_out, axes, FORWARD, real_grid.data(), complex_grid.data(),
      1.0);
  
  // Backward FFT (complex to real)
  c2r(actual_shape, stride_out, stride_in, axes, BACKWARD, complex_grid.data(), real_grid.data(),
      1.0 / static_cast<double>(shape[0] * shape[1] * shape[2]));
  
  // Calculate error between transformed data and original data
  check(real_grid, RelationalOperator::EQUAL, Approx(original_grid).margin(epsilon), "The root "
        "mean-squared error in real-to-complex FFT transformation of dimension " +
        std::to_string(shape[0]) + " x " + std::to_string(shape[1]) + " x " +
        std::to_string(shape[2]) + " exceeds the permissible threshold.");
}

//-------------------------------------------------------------------------------------------------
// Test the real FFT (RFFT) functionality in a three-dimensional array by generating random data,
// performing forward and backward FFT operations, and comparing the result with the original data.
//
// Arguments:
//   rng:      Random number generator to use for generating input data.
//   shape:    Dimensions of the transform to set up and perform
//   epsilon:  Tolerance for deviations in the input and output arguments
//-------------------------------------------------------------------------------------------------
void testComplex3DFFT(Xoshiro256ppGenerator *rng, const std::vector<size_t> &shape,
                      const double epsilon = 5.0e-10) {
  if (shape.size() != 3) {
    rtErr("A two-dimensional array is expected (" + std::to_string(shape.size()) +
          " dimensions were provided).", "testComplex2DFFT");
  }
  if (shape[0] <= 0 || shape[1] <= 0 || shape[2] <= 0) {
    rtErr("Positive dimensions for the grid must be specified (" + std::to_string(shape[0]) +
          " and " + std::to_string(shape[1]) + " were provided).", "testComplex2DFFT");
  }
  std::vector<std::complex<double>> complex_grid(shape[0] * shape[1] * shape[2]);
  for (int i = 0; i < complex_grid.size(); i++) {
    complex_grid[i] = { round(1.0e8 * rng->gaussianRandomNumber()) / 1.0e8,
                        round(1.0e8 * rng->gaussianRandomNumber()) / 1.0e8 };
  }
  std::vector<std::complex<double>> original_grid = complex_grid;
  std::vector<size_t> actual_shape = { shape[2], shape[1], shape[0] };
  std::vector<size_t> cshape = actual_shape;

  // Compute the strides
  std::vector<size_t> axes = { 0, 1, 2 };
  std::vector<ptrdiff_t> stride_in(3);
  std::vector<ptrdiff_t> stride_out(3);
  size_t tmp_in = sizeof(std::complex<double>);
  size_t tmp_out = sizeof(std::complex<double>);
  for (int i = 2; i >= 0; i--) {
    stride_in[i] = tmp_in;
    tmp_in *= actual_shape[i];
    stride_out[i] = tmp_out;
    tmp_out *= cshape[i];
  }
  
  // Forward FFT (real to complex)
  c2c(actual_shape, stride_in, stride_out, axes, FORWARD, complex_grid.data(), complex_grid.data(),
      1.0);
  
  // Backward FFT (complex to real)
  c2c(cshape, stride_out, stride_in, axes, BACKWARD, complex_grid.data(), complex_grid.data(),
      1.0 / static_cast<double>(shape[0] * shape[1] * shape[2]));
  
  // Calculate error between transformed data and original data
  std::vector<double> real_part(complex_grid.size());
  std::vector<double> imag_part(complex_grid.size());
  std::vector<double> real_init(complex_grid.size());
  std::vector<double> imag_init(complex_grid.size());
  for (int i = 0; i < complex_grid.size(); i++) {
    real_part[i] = complex_grid[i].real();
    imag_part[i] = complex_grid[i].imag();
    real_init[i] = original_grid[i].real();
    imag_init[i] = original_grid[i].imag();
  }
  check(real_part, RelationalOperator::EQUAL, Approx(real_init).margin(epsilon), "The root "
        "mean-squared error of the real part in complex FFT transformation of dimension " +
        std::to_string(shape[0]) + " x " + std::to_string(shape[1]) + " x " +
        std::to_string(shape[2]) + " exceeds the permissible threshold.");
  check(imag_part, RelationalOperator::EQUAL, Approx(imag_init).margin(epsilon), "The root "
        "mean-squared error of the imaginary part in complex FFT transformation of dimension " +
        std::to_string(shape[0]) + " x " + std::to_string(shape[1]) + " x " +
        std::to_string(shape[2]) + " exceeds the permissible threshold.");
}

//-------------------------------------------------------------------------------------------------
// Test one of the frequencies in the complex data output from a forward FFT against a known
// result.
//
// Arguments:
//   cmplx_result:  The complex data output from an FFT (the host side is expected to be updated,
//                  perhaps with a download from the GPU memory where the transform may have been
//                  performed)
//   row:           The row at which to query the result
//   column:        The column at which to query the result
//   slab:          The slab at which to query the result (relevant for three-dimensional FFTs)
//   orientation:   Specify whether the row, column, and slab indices are expected to be transposed
//                  upon completion of the FFT operation
//   target:        Real and complex parts of the expected result
//   tol:           The margin of error in FFT results for a successful test
//-------------------------------------------------------------------------------------------------
template <typename T, typename T2>
void checkFFTComplexPoint(const Hybrid<T> &cmplx_result, const int row, const int column,
                          const int slab, const int row_count, const int column_count,
                          const int slab_count, const TransposeState orientation,
                          const T2 target, const double tol) {
  int index;
  std::string transpose_detail;
  switch (orientation) {
  case TransposeState::AS_IS:
    if (slab_count > 0) {
      index = (((column_count * slab) + column) * row_count) + row;
    }
    else {
      index = (row_count * column) + row;
    }
    transpose_detail = std::string("");
    break;
  case TransposeState::TRANSPOSE:
    if (slab_count > 0) {
      index = (((column_count * row) + column) * row_count) + slab;
    }
    else {
      index = (column_count * row) + column;
    }
    transpose_detail = std::string(" transposed");
    break;
  }
  std::string slab_str;
  if (slab_count > 0) {
    slab_str = std::string(" x ") + std::to_string(slab_count);
  }
  const T2* host_ptr = reinterpret_cast<const T2*>(cmplx_result.data());
  const std::string dimensionality = (slab_count > 0) ? std::string("three") : std::string("two");
  check(host_ptr[index].x, RelationalOperator::EQUAL, Approx(target.x).margin(tol), "The real "
        "part of the frequency at index (" + std::to_string(row) + ", " + std::to_string(column) +
        transpose_detail + ") was not computed as expected in a " + dimensionality +
        "-dimensional problem of size " + std::to_string(row_count) + " x " +
        std::to_string(column_count) + slab_str + ".");
  check(host_ptr[index].y, RelationalOperator::EQUAL, Approx(target.y).margin(tol),
        "The imaginary part of the frequency at index (" + std::to_string(row) + ", " +
        std::to_string(column) + transpose_detail + ") was not computed as expected in a " +
        dimensionality + "-dimensional problem of size " + std::to_string(row_count) + " x " +
        std::to_string(column_count) + slab_str + ".");
}

//-------------------------------------------------------------------------------------------------
// Run tests of various HPC packages for GPU-based Fast Fourier Transforms.
//
// Arguments:
//   rng:      The source of random numbers for very large arrays
//-------------------------------------------------------------------------------------------------
#ifdef STORMM_USE_HPC
template <typename T, typename T2> void testHpcFFT(Xoshiro256ppGenerator *rng) {

  // Try one-dimensional arrays of 8, 16, 24, 25, 32, 48, 64, 128, 256, 512, and larger sizes.
  const std::vector<int> sd_fft_lengths = { 8, 16, 24, 25, 32, 48, 64, 128, 256, 512, 1024,
                                            327680 };
  Hybrid<T> g_field(maxValue(sd_fft_lengths), "fft_real_buffer");
  Hybrid<T> g_complex(maxValue(sd_fft_lengths) + 2, "fft_complex_buffer");
  const size_t ct = std::type_index(typeid(T)).hash_code();
  const bool tcalc_is_double = (ct == double_type_index);
  const HybridTargetLevel devc_tier = HybridTargetLevel::DEVICE;
  double tol;
  if (ct == double_type_index) {
    tol = 2.0e-8;
  }
  else if (ct == float_type_index) {
    tol = 6.0e-6;      
  }
  else {
    std::string inv_type("");
    if (isScalarType<T>()) {
      inv_type = getStormmScalarTypeName<T>() + " ";
    }
    else if (isHpcVectorType<T>()) {
      inv_type = getHpcVectorTypeName<T>() + " ";
    }
    rtErr("An invalid data type " + inv_type + "was specified.", "testHpcFFT");
  }
  for (size_t i = 0; i < sd_fft_lengths.size(); i++) {
    const int npts = sd_fft_lengths[i];
    std::vector<T> x;
    if (npts == 8) {
      x = {  3.17694902, -2.65208103,  1.05017029,  4.25793448,  5.31215437, -9.86488889,
             2.27428943, -7.00268227 };
    }
    else if (npts == 16) {
      x = {  4.74111119,  0.81511486, -2.13973013,  0.52978928, -1.79337729, -2.91455126,
             1.01172005, -2.01859928,  0.87512968, -2.00326782, -2.24594094,  1.47920875,
            -2.60376280, -1.18260077,  0.51767092,  2.36817165 };
    }
    else if (npts == 24) {
      x = { -4.21365442,  1.43155025, -3.45461600, -7.74975828,  9.90276153, -1.56531276,
             1.33016668,  0.32688376, -1.96036372,  2.03212827,  1.58277438, -2.91534328,
             1.11695635, 12.98817405,  4.34801988,  3.68996465, -2.32321051, -5.83896880,
             8.45186116, -7.70571218,  3.17082051,  2.27453828,  0.00715778,  6.49104959 };
    }
    else if (npts == 25) {
      x = {  1.83705088, -2.03380454, -5.28353594, -2.62746028,  0.74534075, -3.87908843,
            -0.52614981, -2.68795284,  1.93858433,  2.63784353, -7.98046265,  4.28247492,
            -6.51106568, -0.53860660, -0.54042969, -1.54137652, -7.22990675, -2.21069468,
             4.26841557, -2.97178506,  3.91905109,  6.93504517,  5.88603285, -1.87730830,
             2.48063541 };
    }
    else if (npts == 32) {
      x = {  1.02529432, -2.05148706, -4.67411700,  3.05841986,  3.00151914,  5.68637002,
            -5.00326272, -6.68233314,  1.18487037, -2.17423502, -8.47796520, -4.92390517,
            -8.39267598, -3.74424208, -5.53718829,  3.38483913, -1.67039108, -2.84660492,
             0.76504200, -0.22941044,  3.23196408,  8.25001148,  2.44066217,  2.06605699,
            -7.19599978,  4.74648156, -6.50755875, -1.08087607,  7.48658759, -6.07368744,
             4.10370011,  2.45931335 };
    }
    else if (npts == 64) {
      x = { -4.93594761, -0.20253486, -2.50931319,  4.38201731,  1.56055113, -9.02763516,
             3.22710670,  3.73649728,  8.55665492,  6.67341529,  7.69461118, 10.87441097,
            -2.92115477,  1.49442042, -7.10826686,  2.87653163, -0.59716435, -0.63629142,
            -1.71445559,  6.26916474,  0.52334479, -0.54750472, -4.60258735, -2.17784809,
             1.61992188,  7.97814129, -3.97719935,  4.41320997,  4.05793457,  1.04244760,
             2.49813091,  1.78843231, -1.17885417, -1.01018380, -5.52695861,  3.20837466,
            -4.62910431, -3.35058143, -2.55509548, -5.89897419, -5.54722964,  2.19763116,
            -3.97311475, -7.52060406, -4.78342637, -6.39556963,  2.75808337,  1.08719980,
            -1.82049628,  5.09841578,  2.34304133,  4.29005070,  2.93348144, -0.23195780,
             2.01237258, -3.91734698, -9.99758991,  2.65495001,  3.63422497, -7.96062746,
             2.01683224,  6.85325902,  9.17878839,  1.51070916 };
    }
    else {
      x.resize(npts);
      for (int j = 0; j < npts; j++) {
        x[j] = rng->gaussianRandomNumber();
      }
    }
    g_field.putHost(x, 0, npts);

    // Plan the forward and backward transforms
#  ifdef STORMM_USE_CUDA
    const std::string pkg_name("cuFFT");
    cufftHandle fwd_plan, bkwd_plan;
    if (ct == double_type_index) {
      cufftPlan1d(&fwd_plan, npts, CUFFT_D2Z, 1);
      cufftPlan1d(&bkwd_plan, npts, CUFFT_Z2D, 1);
    }
    else {
      cufftPlan1d(&fwd_plan, npts, CUFFT_R2C, 1);
      cufftPlan1d(&bkwd_plan, npts, CUFFT_C2R, 1);
    }
#  endif
    
    // Upload the data
    g_field.upload(0, npts);

    // Pre-computed answers to various forward FFT operations
    std::vector<double> fftx_ans;
    if (npts == 8) {
      fftx_ans = {  -3.44815460,   0.00000000,  -4.99743847, -11.83856467,   5.16464367,
                     9.77222213,   0.72702777, -14.28680295,  27.07528082,   0.00000000 };
    }
    else if (npts == 16) {
      fftx_ans = {  -4.56391391,   0.00000000,  10.54789401,   1.84262564,  10.89696279,
                     2.68469793,   6.17455314,   0.80882358,   4.07538088,   7.64387539,
                     2.10589609,   0.03694831,   9.12979913,  -9.14542615,  -3.36441720,
                     4.31229241,   1.28955527,   0.00000000 };
    }
    else if (npts == 24) {
      fftx_ans = {  21.41786717,   0.00000000, -15.74138225,   5.89015441,  12.35498705,
                    -1.95540846, -13.63254212,  30.89076569,   2.40794524,  -3.29422954,
                    15.20009031,  25.18258199,  -6.57205414, -19.18502503,   7.73945420,
                    -7.57635921,  -0.31062955, -15.53224671, -31.74299196,  32.68011609,
                   -44.41911064,   7.86976397,   6.19370719,  -2.29172088,  14.49948007,
                     0.00000000 };
    }
    else if (npts == 25) {
      fftx_ans = { -13.50915327,   0.00000000,  16.24316241,  17.02680539,  -8.71421112,
                    29.44674899,  -8.78249438,   6.00732458,   0.85449066,  -8.41182705,
                     1.32486062,   7.91500862,  15.86976950,  -7.66790926,  28.11309276,
                    -1.45380743,  -8.59684465,  18.92837839,   9.02128508,  -8.02997445,
                   -13.68234806, -10.84186454,   5.77494098,  30.20565873,  -7.70799118,
                   -10.75875568 };
    }
    else if (npts == 32) {
      fftx_ans = { -24.37480797,   0.00000000,  16.57153219,  34.71751535,  11.04931494,
                   -37.41907339,  -2.23970188,  11.85780646, -13.42657794,  22.56916609,
                   -17.31698373,  -9.54753370,  11.66441879,  29.14253784,  14.27223999,
                    14.94365480,  21.56185634,  -3.74050203,   9.67536841,  -7.97270191,
                   -18.19362385,  28.93105292, -12.75202479,   9.33789784, -10.54066406,
                    -7.22785435,  -1.17009109, -48.34537912,  16.94402072,  -9.07227187,
                    14.52514409,  -0.24049728, -24.06423007,   0.00000000 };
    }
    else if (npts == 64) {
      fftx_ans = {  15.78874131,   0.00000000,  40.77554505, -48.54247231,   0.71405460,
                    17.16401876, -43.77914765, -48.77194472,   8.26847719,  22.28174193,
                    19.31674613,  51.87697120,   2.52438123,  72.32849109,  20.86483901,
                     9.15342735,  -2.78380304, -27.23813717, -28.67909645,  -0.72828063,
                   -24.42039244,   1.29987372,  -5.25615981,  -2.01397691, -34.22498373,
                    46.37722790, -19.05976789,   0.61935536, -32.72706661,  -5.45171397,
                   -55.64149975,  34.70831307, -16.52161469,   4.37077600,   3.23423900,
                   -31.77305004,  22.60180121,  -4.87111761,  21.81667482, -31.46956952,
                    30.07622758,  14.01273271,  16.96751335, -15.90066880,  -0.34972425,
                     1.55687666, -16.95367543,  39.71048117, -22.53452472, -46.11352971,
                    -5.09888935, -17.47778037, -13.06015052, -30.39332084, -43.65056563,
                     1.38806198, -16.77659968,  -2.01040373,   9.52291080,  18.25244787,
                    15.13996759,   5.86254673,  25.50683877, -26.80495927, -43.31449769,
                     0.00000000 };
    }

    // Perform the transform
    if (ct == double_type_index) {
      double* gfield_ptr = reinterpret_cast<double*>(g_field.data(devc_tier));
      double2* gcmplx_ptr = reinterpret_cast<double2*>(g_complex.data(devc_tier));
#  ifdef STORMM_USE_CUDA
      cufftExecD2Z(fwd_plan, gfield_ptr, gcmplx_ptr);
#  endif
      if (npts <= 64 && npts != 48) {
        g_complex.download();
        check(g_complex.readHost(0, 2 * ((npts / 2) + 1)), RelationalOperator::EQUAL,
              Approx(fftx_ans).margin(1.0e-7), "The complex transform of a data set with " +
              std::to_string(npts) + " real entries does not evaluate by " + pkg_name + " to the "
              "expected result.");
      }
#  ifdef STORMM_USE_CUDA
      cufftExecZ2D(bkwd_plan, gcmplx_ptr, gfield_ptr);
#  endif
    }
    else {
      float* gfield_ptr = reinterpret_cast<float*>(g_field.data(devc_tier));
      float2* gcmplx_ptr = reinterpret_cast<float2*>(g_complex.data(devc_tier));
#  ifdef STORMM_USE_CUDA
      cufftExecR2C(fwd_plan, gfield_ptr, gcmplx_ptr);
#  endif
      if (npts <= 64 && npts != 48) {
        g_complex.download();
        check(g_complex.readHost(0, 2 * ((npts / 2) + 1)), RelationalOperator::EQUAL,
              Approx(fftx_ans).margin(5.0e-5), "The complex transform of a data set with " +
              std::to_string(npts) + " real entries does not evaluate by " + pkg_name + " to the "
              "expected result.");
      }
#  ifdef STORMM_USE_CUDA
      cufftExecC2R(bkwd_plan, gcmplx_ptr, gfield_ptr);
#  endif
    }
    
    // Download the result and normalize
    g_field.download(0, npts);
    const double dinv_npts = 1.0 / static_cast<double>(npts);
    for (int i = 0; i < npts; i++) {
      g_field.putHost(g_field.readHost(i) * dinv_npts, i);
    }

    // Check the outcome
    check(g_field.readHost(0, npts), RelationalOperator::EQUAL, Approx(x).margin(tol),
          "A one-dimensional array of " + std::to_string(npts) + " " +
          getStormmScalarTypeName<T>() + "s did not make a round-trip through forward and "
          "backward FFTs intact.");
  }

  // Try two-dimensional real-to-complex transformations
  const std::vector<int2> rectangle_fft_sizes = { { 5, 5 }, { 6, 9 }, { 16, 4 }, { 12, 9 } };
  for (int i = 0; i < rectangle_fft_sizes.size(); i++) {
    const int2 npts = rectangle_fft_sizes[i];
    if (npts.x * npts.y > g_field.size()) {
      rtErr("A rectangle of " + std::to_string(npts.x) + " x " + std::to_string(npts.y) +
            " points is too large for the pre-allocated real array.", "testHpcFFT");
    }
    if (((npts.x / 2) + 1) * npts.y > g_complex.size()) {
      rtErr("A rectangle of " + std::to_string(npts.x) + " x " + std::to_string(npts.y) +
            " points is too large for the pre-allocated complex array.", "testHpcFFT");
    }
    std::vector<T> x;
    if (npts.x == 5 && npts.y == 5) {
      x = {  -3.27754176,  -0.47440397,   1.53624322,  -5.14275674,   3.76596582,  -2.32511370,
              2.35402898,  -5.27320945,  -3.00371158,   4.18858300,  -1.50040821,  -6.78514855,
             -3.75402925,  -0.53374168,   8.32293817, -13.51371302,   7.44328740,  -5.19651117,
             -2.07922717,   8.59473250,  -2.38735179,  -2.54920395,   5.13080051,  -1.37141558, 
              5.21069997 };
    }
    else if (npts.x == 6 && npts.y == 9) {
      x = {  -2.32692176,  -3.17158209,  -4.61037865,  -3.79507679,  -2.06538315,   4.18155992, 
              4.02168703,   0.02405967,  -9.87909411,  -2.56341878,   3.00922334,   2.19068500, 
              5.66042327,  -2.09602974,   2.69083684,   4.30740378,   4.04392720,   1.45378851, 
              1.50089367,  -0.75098882,  -0.46193515,  -0.00933502, -11.53959375,   0.34932809, 
              8.07396796,   7.86537186,  -4.09815649,   3.32990557,  -0.31250363,   5.82272797, 
              3.99003801,  -3.87024846,  -0.05673339,  14.35624835,  -3.18260638,  -5.89975462, 
             -5.26221512,  -8.62116398,   9.00719579,  -5.29126952,   0.14560620,   3.73183172, 
            -11.74025138,  -7.17313068,  -1.46888136,   2.94001304,  -2.53583777,   4.78613316, 
              6.25479576,  -6.73146984,  -0.36084164,   6.74213452,   1.47331757, -13.18483423 };
    }
    else if (npts.x == 16 && npts.y == 4) {
      x = {   0.22957476,   0.32853126,   0.36787771,   0.49044577,   0.85450582,   0.66989332,
              0.91872269,   0.52863291,   0.37780102,   0.23688771,   0.80429305,   0.05240436,
              0.12539040,   0.35385460,   0.25804278,   0.82130340,   0.42964504,   0.96651358,
              0.74281294,   0.02580511,   0.16842366,   0.25273676,   0.29048074,   0.59656608,
              0.13011002,   0.72130976,   0.56621904,   0.69726983,   0.03986420,   0.28816825,
              0.61441204,   0.52329818,   0.56819675,   0.88456553,   0.23093443,   0.59543566,
              0.15914994,   0.79031135,   0.01609413,   0.48744888,   0.88607393,   0.20187446,
              0.47818122,   0.18023682,   0.35631707,   0.07248257,   0.22502734,   0.72925102,
              0.46384574,   0.25311277,   0.00435210,   0.44414388,   0.32602516,   0.82073284,
              0.73543170,   0.71186790,   0.23895895,   0.36151463,   0.11134394,   0.56212906,
              0.16673364,   0.08248443,   0.10354043,   0.88930170 };
    }
    else if (npts.x == 12 && npts.y == 9) {
      x = {  -1.15414360,   2.34857216,   1.24978327,   0.91783576,   0.47791522,   0.87757142,
              1.15234820,  -1.36877466,   0.27095391,  -0.15903985,  -1.47475186,   1.24664710,
             -1.11119068,   0.03574012,  -0.43371790,  -0.60413390,   0.48559955,   0.93431797,
             -1.80862172,  -0.21270329,   1.05479336,  -0.22524861,  -1.46604875,  -2.17749624,
              0.34056939,  -1.57958428,  -0.78828890,  -0.22503347,   0.20597271,  -0.58062807,
             -0.08086487,  -0.33712036,   1.10958460,   0.04583221,  -0.51421385,   1.27503867,
              1.23098951,  -1.10058325,  -1.28419620,  -0.37376524,   0.80758595,  -0.19223241,
             -0.83319934,   0.91388857,   2.30534617,  -0.43596816,   2.29885699,  -1.32228684,
             -0.42146045,   0.50869319,  -1.15942148,  -0.15614589,  -0.42084715,  -1.39447000,
             -1.81487769,   1.87159709,  -0.11553168,  -1.16432473,   0.96381812,   1.67797343,
             -0.77558695,  -1.14208526,   0.46790966,  -0.66403821,  -1.23879395,  -0.39330126,
              0.52967961,   0.09684297,  -0.94860163,   0.55067579,   0.80645139,  -0.22163221,
              1.13834526,   0.61686222,   0.68040667,  -0.99348972,  -0.33918065,  -0.00450279,
             -1.11174721,  -0.20255913,   0.24853527,  -0.06002839,   0.15799193,   1.76875839,
             -1.04542658,  -1.57282139,  -1.90756102,  -1.36702979,   0.47865183,   0.35347184,
             -1.20960817,  -0.38704038,  -1.26047014,   1.57554457,  -1.23085815,   2.22709534,
              0.39525485,   0.80404599,   1.20948821,  -1.79575658,   0.94496571,  -0.77887987,
             -0.20647136,  -0.94579827,  -2.26717075,  -0.57115710,   0.76430001,   1.10699990 };
    }
    g_field.putHost(x, 0, npts.x * npts.y);

    // Plan the forward and backward transforms
    cufftHandle fwd_plan, bkwd_plan;
    if (ct == double_type_index) {
#  ifdef STORMM_USE_CUDA
      cufftPlan2d(&fwd_plan, npts.y, npts.x, CUFFT_D2Z);
      cufftPlan2d(&bkwd_plan, npts.y, npts.x, CUFFT_Z2D);
#  endif
    }
    else {
#  ifdef STORMM_USE_CUDA
      cufftPlan2d(&fwd_plan, npts.y, npts.x, CUFFT_R2C);
      cufftPlan2d(&bkwd_plan, npts.y, npts.x, CUFFT_C2R);
#  endif
    }
    
    // Upload the data
    g_field.upload(0, npts.x * npts.y);

    // Correct answers produced by a third-party matrix algebra program
    T2 zeroth_freq, first_freq, second_freq;
    zeroth_freq.x = static_cast<T>(sum<double, T>(x));
    zeroth_freq.y = 0.0;
    if (npts.x == 5 && npts.y == 5) {
      first_freq = { 8.77531005, -26.83360653 };
      second_freq = { -20.88624495, 29.57157352 };
    }
    else if (npts.x == 6 && npts.y == 9) {
      first_freq = { 0.06737186, -23.89487939 };
      second_freq = { 32.08395597, -8.61275927 };
    }
    else if (npts.x == 16 && npts.y == 4) {
      first_freq = { 0.58537092, -1.94205728 };
      second_freq = { -0.58678388, 0.17166108 };
    }
    else if (npts.x == 12 && npts.y == 9) {
      first_freq = { 0.61400869, -5.24894025 };
      second_freq = { -5.48162413, -6.01775245 };
    }
    
    // Perform the transform
    if (ct == double_type_index) {
      double* gfield_ptr = reinterpret_cast<double*>(g_field.data(devc_tier));
      double2* gcmplx_ptr = reinterpret_cast<double2*>(g_complex.data(devc_tier));
#  ifdef STORMM_USE_CUDA
      cufftExecD2Z(fwd_plan, gfield_ptr, gcmplx_ptr);
#  endif
      g_complex.download();
      checkFFTComplexPoint<T, T2>(g_complex, 0, 0, -1, (npts.x / 2) + 1, npts.y, 0,
                                  TransposeState::AS_IS, zeroth_freq, tol);
      checkFFTComplexPoint<T, T2>(g_complex, 1, 1, -1, (npts.x / 2) + 1, npts.y, 0,
                                  TransposeState::AS_IS, first_freq, tol);
      checkFFTComplexPoint<T, T2>(g_complex, 2, 2, -1, (npts.x / 2) + 1, npts.y, 0,
                                  TransposeState::AS_IS, second_freq, tol);
#  ifdef STORMM_USE_CUDA
      cufftExecZ2D(bkwd_plan, gcmplx_ptr, gfield_ptr);
#  endif
    }
    else {
      float* gfield_ptr = reinterpret_cast<float*>(g_field.data(devc_tier));
      float2* gcmplx_ptr = reinterpret_cast<float2*>(g_complex.data(devc_tier));
#  ifdef STORMM_USE_CUDA
      cufftExecR2C(fwd_plan, gfield_ptr, gcmplx_ptr);
#  endif
      g_complex.download();
      checkFFTComplexPoint<T, T2>(g_complex, 0, 0, -1, (npts.x / 2) + 1, npts.y, 0,
                                  TransposeState::AS_IS, zeroth_freq, tol);
      checkFFTComplexPoint<T, T2>(g_complex, 1, 1, -1, (npts.x / 2) + 1, npts.y, 0,
                                  TransposeState::AS_IS, first_freq, tol);
      checkFFTComplexPoint<T, T2>(g_complex, 2, 2, -1, (npts.x / 2) + 1, npts.y, 0,
                                  TransposeState::AS_IS, second_freq, tol);
#  ifdef STORMM_USE_CUDA
      cufftExecC2R(bkwd_plan, gcmplx_ptr, gfield_ptr);
#  endif
    }
    
    // Download the result and normalize
    g_field.download(0, npts.x * npts.y);
    const double dinv_npts = 1.0 / static_cast<double>(npts.x * npts.y);
    for (int i = 0; i < npts.x * npts.y; i++) {
      g_field.putHost(g_field.readHost(i) * dinv_npts, i);
    }

    // Check the outcome
    check(g_field.readHost(0, npts.x * npts.y), RelationalOperator::EQUAL, Approx(x).margin(tol),
          "A two-dimensional array of " + std::to_string(npts.x) + " x " + std::to_string(npts.y) +
          " " + getStormmScalarTypeName<T>() + "s did not make a round-trip through forward and "
          "backward FFTs intact.");
  }

  // Try three-dimensional real-to-complex transformations
  const std::vector<int3> block_fft_sizes = { {  4,  4,  4 }, {  5,  5,  5 }, {  5,  6,  8 },
                                              { 16,  4,  7 } };
  for (int i = 0; i < block_fft_sizes.size(); i++) {
    const int3 npts = block_fft_sizes[i];
    if (npts.x * npts.y * npts.z > g_field.size()) {
      rtErr("A block of " + std::to_string(npts.x) + " x " + std::to_string(npts.y) + " x " +
            std::to_string(npts.z) + " points is too large for the pre-allocated real array.",
            "testHpcFFT");
    }
    if (((npts.x / 2) + 1) * npts.y * npts.z > g_complex.size()) {
      rtErr("A block of " + std::to_string(npts.x) + " x " + std::to_string(npts.y) + " x " +
            std::to_string(npts.z) + " points is too large for the pre-allocated complex array.",
            "testHpcFFT");
    }
    std::vector<T> x;
    if (npts.x == 4 && npts.y == 4 && npts.z == 4) {
      x = {  -0.24182799,  -1.10068872,  -0.31872197,   1.92209237,  -0.39652697,  -0.58492346,
             -0.30576802,   1.54740771,   0.10916699,  -0.25552662,   0.01456756,  -0.07569496,
              0.66439004,   0.76021431,   0.38605576,   1.39529652,  -1.18703048,  -0.12954200,
             -0.54828318,  -0.36819622,  -0.08705386,   2.59850990,  -0.08511539,  -0.26862134,
             -0.02751446,  -1.52674993,   0.89406158,   1.04298375,  -0.71901925,   1.95192262,
              0.08321914,  -1.42182847,   0.80632141,   1.27734419,   0.16322959,  -0.59216798,
             -0.35399023,   1.09402464,   0.14469383,   0.17366149,   0.03662405,  -1.75575704,
             -0.62766204,  -1.18376229,  -0.04312667,  -1.29821333,  -2.44974055,  -0.12305022,
             -0.17406734,   0.96150724,  -1.10778502,   1.95584298,  -1.19133676,  -0.49627837,
             -0.32982965,  -1.09291909,   1.71806197,   0.44009324,  -0.54220322,  -0.09029556,
              0.73021468,   0.15647556,   0.07151922,   0.14636218 };
    }
    else if (npts.x == 5 && npts.y == 5 && npts.z == 5) {
      x = {  -0.87609096,  -0.77957073,  -0.63420409,  -0.08526856,   2.46449347,   1.58082636,
              0.35619115,  -0.93068242,  -1.13405480,   1.33434307,   0.32031741,  -0.07192847,
              0.79668848,  -1.25116873,  -0.34284036,   0.83551973,  -0.73949967,   1.02236044,
             -0.53736742,  -1.12471200,  -0.54511479,   0.87183470,  -0.26975722,   0.45462716,
             -0.90470034,   0.43155649,   0.30811858,   1.34588637,   0.50558002,  -0.40204743,
              2.03846792,   1.56621079,   0.72524273,  -0.70108889,  -0.38178494,  -1.59009185,
             -0.28780127,  -0.04162096,  -0.03146920,  -0.32771009,   0.20173189,  -0.38174555,
             -0.45765539,   0.04874962,  -0.76310332,   1.03211354,   0.14895746,   1.81868755,
              0.03261233,   0.36858538,   1.49037898,   0.05296502,  -1.32649099,  -1.57330374,
             -1.58053665,   0.24419772,  -0.72059230,  -1.02553541,  -2.08543464,   0.27658976,
              1.29799848,   0.92407494,   2.12988604,   0.08007394,  -1.35900169,   1.11381200,
             -2.24436575,  -0.05981749,  -1.27245440,  -0.75776356,   2.31840908,   0.51685653,
             -0.39417538,   0.49915238,  -0.43487834,   1.34983842,  -0.80292207,   0.45321704,
              0.59484347,   0.74029429,  -1.29287070,  -0.52058155,  -0.34945023,  -0.13051290,
              0.73431025,  -0.98897610,  -0.88122308,   0.67310982,   0.06736756,   0.15115284,
             -0.44869075,  -0.45354876,   0.96895351,  -0.39067758,  -1.42237541,   0.40311752,
             -0.56648948,   0.07770063,  -0.07211839,   0.19157903,  -1.36202727,  -0.33965758,
              0.19071050,  -0.80741863,   0.28108408,   1.25018052,   0.27484534,  -1.12150942,
             -0.77926643,   2.55474155,   2.03014093,  -0.42395556,   2.23481587,   0.77201007,
              1.20732182,  -0.11967842,  -0.87707608,   0.05868741,  -2.14638752,   0.15701399,
             -1.13362310,   0.70199726,  -0.93275222,   0.34620291,   1.87221895 };
    }
    else if (npts.x == 5 && npts.y == 6 && npts.z == 8) {
      x = std::vector<T>(five_six_eight.begin(), five_six_eight.end());
    }
    else if (npts.x == 16 && npts.y == 4 && npts.z == 7) {
      x = std::vector<T>(sixteen_four_seven.begin(), sixteen_four_seven.end());
    }
    g_field.putHost(x, 0, npts.x * npts.y * npts.z);

    // Plan the forward and backward transforms
    cufftHandle fwd_plan, bkwd_plan;
    if (ct == double_type_index) {
#  ifdef STORMM_USE_CUDA
      cufftPlan3d(&fwd_plan, npts.z, npts.y, npts.x, CUFFT_D2Z);
      cufftPlan3d(&bkwd_plan, npts.z, npts.y, npts.x, CUFFT_Z2D);
#  endif
    }
    else {
#  ifdef STORMM_USE_CUDA
      cufftPlan3d(&fwd_plan, npts.z, npts.y, npts.x, CUFFT_R2C);
      cufftPlan3d(&bkwd_plan, npts.z, npts.y, npts.x, CUFFT_C2R);
#  endif
    }
    
    // Upload the data
    g_field.upload(0, npts.x * npts.y * npts.z);

    // Correct answers produced by a third-party matrix algebra program
    T2 zeroth_freq, first_freq, second_freq;
    zeroth_freq.x = static_cast<T>(sum<double, T>(x));
    zeroth_freq.y = 0.0;
    if (npts.x == 4 && npts.y == 4 && npts.z == 4) {
      first_freq = { -3.05735089, 0.04572820 };
      second_freq = { 7.18366791, 0.00000000 };
    }
    else if (npts.x == 5 && npts.y == 5 && npts.z == 5) {
      first_freq = { -7.52964025, -2.51476524 };
      second_freq = { -3.77320356, 7.92055736 };
    }
    else if (npts.x == 5 && npts.y == 6 && npts.z == 8) {
      first_freq = { -15.05667162, -10.78948826 };
      second_freq = { 11.71704074, -1.52016957 };
    }
    else if (npts.x == 16 && npts.y == 4 && npts.z == 7) {
      first_freq = { -16.09308673, -6.53780183 };
      second_freq = { 0.67388598, 1.58951931 };
    }
    
    // Perform the transform
    if (ct == double_type_index) {
      double* gfield_ptr = reinterpret_cast<double*>(g_field.data(devc_tier));
      double2* gcmplx_ptr = reinterpret_cast<double2*>(g_complex.data(devc_tier));
#  ifdef STORMM_USE_CUDA
      cufftExecD2Z(fwd_plan, gfield_ptr, gcmplx_ptr);
#  endif
      g_complex.download();
      checkFFTComplexPoint<T, T2>(g_complex, 0, 0, 0, (npts.x / 2) + 1, npts.y, npts.z,
                                  TransposeState::AS_IS, zeroth_freq, tol);
      checkFFTComplexPoint<T, T2>(g_complex, 1, 1, 1, (npts.x / 2) + 1, npts.y, npts.z,
                                  TransposeState::AS_IS, first_freq, tol);
      checkFFTComplexPoint<T, T2>(g_complex, 2, 2, 2, (npts.x / 2) + 1, npts.y, npts.z,
                                  TransposeState::AS_IS, second_freq, tol);
#  ifdef STORMM_USE_CUDA
      cufftExecZ2D(bkwd_plan, gcmplx_ptr, gfield_ptr);
#  endif
    }
    else {
      float* gfield_ptr = reinterpret_cast<float*>(g_field.data(devc_tier));
      float2* gcmplx_ptr = reinterpret_cast<float2*>(g_complex.data(devc_tier));
#  ifdef STORMM_USE_CUDA
      cufftExecR2C(fwd_plan, gfield_ptr, gcmplx_ptr);
#  endif
      g_complex.download();
      checkFFTComplexPoint<T, T2>(g_complex, 0, 0, 0, (npts.x / 2) + 1, npts.y, npts.z,
                                  TransposeState::AS_IS, zeroth_freq, tol);
      checkFFTComplexPoint<T, T2>(g_complex, 1, 1, 1, (npts.x / 2) + 1, npts.y, npts.z,
                                  TransposeState::AS_IS, first_freq, tol);
      checkFFTComplexPoint<T, T2>(g_complex, 2, 2, 2, (npts.x / 2) + 1, npts.y, npts.z,
                                  TransposeState::AS_IS, second_freq, tol);
#  ifdef STORMM_USE_CUDA
      cufftExecC2R(bkwd_plan, gcmplx_ptr, gfield_ptr);
#  endif
    }
    
    // Download the result and normalize
    g_field.download(0, npts.x * npts.y * npts.z);
    const double dinv_npts = 1.0 / static_cast<double>(npts.x * npts.y * npts.z);
    for (int i = 0; i < npts.x * npts.y * npts.z; i++) {
      g_field.putHost(g_field.readHost(i) * dinv_npts, i);
    }

    // Check the outcome
    check(g_field.readHost(0, npts.x * npts.y * npts.z), RelationalOperator::EQUAL,
          Approx(x).margin(tol), "A three-dimensional array of " + std::to_string(npts.x) + " x " +
          std::to_string(npts.y) + " x " + std::to_string(npts.z) + " " +
          getStormmScalarTypeName<T>() + "s did not make a round-trip through forward and "
          "backward FFTs intact.");
  }
}
#endif

//-------------------------------------------------------------------------------------------------
// Check a specific frequency with an FFT management object.
//
// Arguments:
//   ffts:  The FFT management object
//   ans:   The expected frequency.  Its real and imaginary parts are found in the tuple's "x" and
//          "y" members, respectively.
//   tol:   Tolerance for getting the answer correct
//   i:     Position along the problem's first dimension, in frequency space
//   j:     Position along the problem's second dimension, in frequency space (if applicable--the
//          default value of -1 indicates that the problem has only a single dimesion)
//   k:     Position along the problem's third dimension, in frequency space (if applicable)
//   k:     Position along the problem's fourth dimension, in frequency space (if applicable)
//-------------------------------------------------------------------------------------------------
void checkSpecificFrequency(const FFTStage &ffts, double2 ans, double tol, const int i,
                            const int j = -1, const int k = -1, const int m = -1) {
  if (j < 0) {
    check(ffts.getFrequencyData(i).x, RelationalOperator::EQUAL, Approx(ans.x).margin(tol),
          "The real part of the frequency at position { " + std::to_string(i) + " } did not meet "
          "expectations.");
    check(ffts.getFrequencyData(i).y, RelationalOperator::EQUAL, Approx(ans.y).margin(tol),
          "The imaginary part of the frequency at position { " + std::to_string(i) + " } did not "
          "meet expectations.");
  }
  else if (k < 0) {
    check(ffts.getFrequencyData(i, j).x, RelationalOperator::EQUAL, Approx(ans.x).margin(tol),
          "The real part of the frequency at position { " + std::to_string(i) + ", " +
          std::to_string(j) + " } did not meet expectations.");
    check(ffts.getFrequencyData(i, j).y, RelationalOperator::EQUAL, Approx(ans.y).margin(tol),
          "The imaginary part of the frequency at position { " + std::to_string(i) + ", " +
          std::to_string(j) + " } did not meet expectations.");
  }
  else if (m < 0) {
    check(ffts.getFrequencyData(i, j, k).x, RelationalOperator::EQUAL, Approx(ans.x).margin(tol),
          "The real part of the frequency at position { " + std::to_string(i) + ", " +
          std::to_string(j) + ", " + std::to_string(k) + " } did not meet expectations.");
    check(ffts.getFrequencyData(i, j, k).y, RelationalOperator::EQUAL, Approx(ans.y).margin(tol),
          "The imaginary part of the frequency at position { " + std::to_string(i) + ", " +
          std::to_string(j) + ", " + std::to_string(k) + " } did not meet expectations.");
  }
  else {
    check(ffts.getFrequencyData(i, j, k, m).x, RelationalOperator::EQUAL,
          Approx(ans.x).margin(tol), "The real part of the frequency at position { " +
          std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ", " +
          std::to_string(m) + " } did not meet expectations.");
    check(ffts.getFrequencyData(i, j, k, m).y, RelationalOperator::EQUAL,
          Approx(ans.y).margin(tol), "The imaginary part of the frequency at position { " +
          std::to_string(i) + ", " + std::to_string(j) + ", " + std::to_string(k) + ", " +
          std::to_string(m) + " } did not meet expectations.");
  }
}

//-------------------------------------------------------------------------------------------------
// Run a battery of tests on the FFTStage class, for encapsulating the PocketFFT and cuFFT APIs.
//
// Arguments:
//   rng:  The source of random numbers for very large arrays
//   tol:  Tolerance for deviations between the initial data, round-trip results, and intermediate
//         transformed (complex) data
//-------------------------------------------------------------------------------------------------
template <typename T, typename T2> void testFFTEncapsulation(Xoshiro256ppGenerator *rng,
                                                             const HybridTargetLevel mount,
                                                             const double tol) {
  const size_t ct = std::type_index(typeid(T)).hash_code();
  if (ct != float_type_index && ct != double_type_index) {
    rtErr("This function must be invoked with float or double type, plus the two-tuple of the "
          "same.", "testFFTEncapsulation");
  }
  
  // Test one-dimensional FFTs
  const std::vector<int> mono_npts = { 10, 13, 19 };
  const std::vector<double2> third_freq = { { -0.71109066, -1.17186700 },
                                            {  0.84572661,  0.84330843 },
                                            {  1.07159447,  0.61331392 } };
  Hybrid<T> g_field(HybridKind::ARRAY, "fft_substrate");
  Hybrid<T2> g_complex(HybridKind::ARRAY, "fft_freq_space");
  std::vector<T> gfld_ref;
  for (size_t i = 0; i < mono_npts.size(); i++) {
    g_field.resize(mono_npts[i]);
    g_complex.resize(mono_npts[i]);
    gfld_ref.resize(mono_npts[i]);
    const double dn = mono_npts[i];
    double sgn = 1.0;
    for (int j = 0; j < mono_npts[i]; j++) {
      gfld_ref[j] = sgn * round(0.5e8 * pow(static_cast<double>(j) + 0.5, 1.5) / dn) / 1.0e8;
      sgn *= -1.0;
    }
    g_field.putHost(gfld_ref);
#ifdef STORMM_USE_HPC
    if (mount == HybridTargetLevel::DEVICE) {
      g_field.upload();
    }
#endif
    FFTStage mono_ffts(g_field.data(mount), g_complex.data(mount), mount, Normalization::YES,
                       FFTMode::OUT_OF_PLACE, mono_npts[i]);
    mono_ffts.forwardFFT();
#ifdef STORMM_USE_HPC
    if (mount == HybridTargetLevel::DEVICE) {
      g_complex.download();
    }
#endif
    checkSpecificFrequency(mono_ffts, { sum<double>(g_field), 0.0 }, tol, 0);
    checkSpecificFrequency(mono_ffts, third_freq[i], tol, 3);
    for (int j = 0; j < mono_npts[i]; j++) {
      g_field.putHost(0.0, j);
    }
    mono_ffts.backwardFFT();
#ifdef STORMM_USE_HPC
    if (mount == HybridTargetLevel::DEVICE) {
      g_field.download();
    }
#endif
    const std::vector<double> dgfld_ref(gfld_ref.begin(), gfld_ref.end());
    check(g_field.readHost(), RelationalOperator::EQUAL, dgfld_ref, "An FFTStage object "
          "managing " + getEnumerationName(mono_ffts.getPrecision()) + "-precision transforms of "
          "one dimension on the " + getEnumerationName(mount) + " did not guide the data on a "
          "proper round-trip.");
  }

  // Test two-dimensional FFTs
  std::vector<int2> dual_npts = { { 6, 6 }, { 7, 8 } };
  for (size_t i = 0; i < dual_npts.size(); i++) {
    const int total_ipts = dual_npts[i].x * dual_npts[i].y;
    g_field.resize(total_ipts);
    g_complex.resize(total_ipts);
    gfld_ref.resize(total_ipts);
    double sgn = 1.0;
    for (int j = 0; j < dual_npts[i].x; j++) {
      for (int k = 0; k < dual_npts[i].y; k++) {
        const int jk = (k * dual_npts[i].x) + j;
        gfld_ref[jk] = sgn * round(1.0e8 * sin(5.0 / static_cast<double>(jk + 1))) / 1.0e8;
        sgn *= -1.0;
      }
    }
    g_field.putHost(gfld_ref);
#ifdef STORMM_USE_HPC
    if (mount == HybridTargetLevel::DEVICE) {
      g_field.upload();
    }
#endif
    FFTStage dual_ffts(g_field.data(mount), g_complex.data(mount), mount, Normalization::YES,
                       FFTMode::OUT_OF_PLACE, dual_npts[i].x, dual_npts[i].y);
    dual_ffts.forwardFFT();
#ifdef STORMM_USE_HPC
    if (mount == HybridTargetLevel::DEVICE) {
      g_complex.download();
    }
#endif
    double2 first_freq, second_freq;
    if (dual_npts[i].x == 6 && dual_npts[i].y == 6) {
      first_freq = { -2.14207129, 0.23639995 };
      second_freq = { -1.50480871, 0.41246737 };
    }
    else if (dual_npts[i].x == 7 && dual_npts[i].y == 8) {
      first_freq = { -2.18957473, -0.00126432 };
      second_freq = { -1.64761168, 0.40640957 };
    }
    else {
      rtErr("An invalid problem size was submitted for testing.", "testFFTEncapsulation");
    }
    checkSpecificFrequency(dual_ffts, { sum<double>(g_field), 0.0 }, tol, 0, 0);
    checkSpecificFrequency(dual_ffts, first_freq, tol, 1, 1);
    checkSpecificFrequency(dual_ffts, second_freq, tol, 2, 2);
    for (int j = 0; j < total_ipts; j++) {
      g_field.putHost(0.0, j);
    }
    dual_ffts.backwardFFT();
#ifdef STORMM_USE_HPC
    if (mount == HybridTargetLevel::DEVICE) {
      g_field.download();
    }
#endif
    const std::vector<double> dgfld_ref(gfld_ref.begin(), gfld_ref.end());
    check(g_field.readHost(), RelationalOperator::EQUAL, dgfld_ref, "An FFTStage object "
          "managing " + getEnumerationName(dual_ffts.getPrecision()) + "-precision transforms of "
          "two dimensions on the " + getEnumerationName(mount) + " did not guide the data on a "
          "proper round-trip.");
  }
  
  // Test three-dimensional FFTs
  const std::vector<int3> tri_npts = { { 5, 6, 8 } };
  for (size_t i = 0; i < tri_npts.size(); i++) {
    const int total_ipts = tri_npts[i].x * tri_npts[i].y * tri_npts[i].z;
    g_field.resize(total_ipts);
    g_complex.resize(total_ipts);
    std::vector<T> local_fse(five_six_eight.begin(), five_six_eight.end());
    g_field.putHost(local_fse);
#ifdef STORMM_USE_HPC
    if (mount == HybridTargetLevel::DEVICE) {
      g_field.upload();
    }
#endif
    FFTStage tri_ffts(g_field.data(mount), g_complex.data(mount), mount, Normalization::YES,
                      FFTMode::OUT_OF_PLACE, tri_npts[i].x, tri_npts[i].y, tri_npts[i].z);
    tri_ffts.forwardFFT();
#ifdef STORMM_USE_HPC
    if (mount == HybridTargetLevel::DEVICE) {
      g_complex.download();
    }
#endif
    const double2 first_freq = { -15.05667162, -10.78948826 };
    const double2 second_freq = { 11.71704074, -1.52016957 };
    checkSpecificFrequency(tri_ffts, { sum<double>(g_field), 0.0 }, tol, 0, 0, 0);
    checkSpecificFrequency(tri_ffts, first_freq, tol, 1, 1, 1);
    checkSpecificFrequency(tri_ffts, second_freq, tol, 2, 2, 2);
  }

  // Test batched FFTs with PocketFFT
  if (mount == HybridTargetLevel::HOST) {
    std::vector<double> dsix_dozen = gaussianRand(rng, 72, 5.0);
    std::vector<T> six_dozen(dsix_dozen.begin(), dsix_dozen.end());
    std::vector<T> six_dozen_copy = six_dozen;
    std::vector<T2> fsix_dozen(42), fsix_dozen_copy(42);
    CHECK_THROWS(FFTStage bad_batch(six_dozen.data(), fsix_dozen.data(), HybridTargetLevel::HOST,
                                    12, 0, 0, 0, 6),
                 "A batched FFT was staged with no batch strides provided.");
    FFTStage six_batch(six_dozen.data(), fsix_dozen.data(), HybridTargetLevel::HOST,
                       Normalization::YES, FFTMode::OUT_OF_PLACE, 12, 0, 0, 0, 6, 12, 7);
    six_batch.forwardFFT();
    std::vector<double> real_batch(42), real_unitary(42), imag_batch(42), imag_unitary(42);
    for (int i = 0; i < 6; i++) {
      std::vector<T> one_dozen(12);
      std::vector<T2> fone_dozen(7);
      for (int j = 0; j < 12; j++) {
        one_dozen[j] = six_dozen_copy[(12 * i) + j];
      }
      FFTStage one_batch(one_dozen.data(), fone_dozen.data(), HybridTargetLevel::HOST,
                         Normalization::YES, FFTMode::OUT_OF_PLACE, 12);
      one_batch.forwardFFT();
      for (int j = 0; j < 7; j++) {
        real_batch[(7 * i) + j] = fsix_dozen[(7 * i) + j].x;
        real_unitary[(7 * i) + j] = fone_dozen[j].x;
        imag_batch[(7 * i) + j] = fsix_dozen[(7 * i) + j].y;
        imag_unitary[(7 * i) + j] = fone_dozen[j].y;
      }
    }
    check(real_batch, RelationalOperator::EQUAL, Approx(real_unitary).margin(tol),
          "The real-values of the transformed result of batched real-to-complex operations do not "
          "agree with the results of individual transforms performed on the same data.");
    check(imag_batch, RelationalOperator::EQUAL, Approx(imag_unitary).margin(tol), "The imaginary "
          "values of the transformed result of batched real-to-complex operations do not agree "
          "with the results of individual transforms performed on the same data.");
    std::vector<double> dsix_patches = gaussianRand(rng, 216, 5.0);
    std::vector<T> six_patches(dsix_patches.begin(), dsix_patches.end());
    std::vector<T> six_patches_copy = six_patches;
    std::vector<T2> fsix_patches(216), fsix_patches_copy(216);
    FFTStage six_fields(six_patches.data(), fsix_patches.data(), HybridTargetLevel::HOST,
                        Normalization::YES, FFTMode::OUT_OF_PLACE, 4, 9, 0, 0, 6, 36, 27);
    CHECK_THROWS(FFTStage bad_fields(six_patches.data(), fsix_patches.data(),
                                     HybridTargetLevel::HOST, Normalization::YES,
                                     FFTMode::OUT_OF_PLACE, 4, 9, 0, 0, 6, 36, 20),
                 "A batched FFT was staged with an insufficient frequency batch stride.");
    six_fields.forwardFFT();
    real_batch.resize(162);
    real_unitary.resize(162);
    imag_batch.resize(162);
    imag_unitary.resize(162);
    for (int i = 0; i < 6; i++) {
      std::vector<T> one_patch(36);
      std::vector<T2> fone_patch(27);
      for (int j = 0; j < 36; j++) {
        one_patch[j] = six_patches_copy[(36 * i) + j];
      }
      FFTStage one_field(one_patch.data(), fone_patch.data(), HybridTargetLevel::HOST,
                         Normalization::YES, FFTMode::OUT_OF_PLACE, 4, 9);
      one_field.forwardFFT();
      for (int j = 0; j < 27; j++) {
        real_batch[(27 * i) + j] = fsix_patches[(27 * i) + j].x;
        real_unitary[(27 * i) + j] = fone_patch[j].x;
        imag_batch[(27 * i) + j] = fsix_patches[(27 * i) + j].y;
        imag_unitary[(27 * i) + j] = fone_patch[j].y;
      }
    }
    check(real_batch, RelationalOperator::EQUAL, Approx(real_unitary).margin(tol),
          "The real-values of the transformed result of batched, two-dimensional real-to-complex "
          "operations do not agree with the results of individual transforms performed on the "
          "same data.");
    check(imag_batch, RelationalOperator::EQUAL, Approx(imag_unitary).margin(tol), "The imaginary "
          "values of the transformed result of batched. two-dimensional real-to-complex "
          "operations do not agree with the results of individual transforms performed on the "
          "same data.");
  }
}

//-------------------------------------------------------------------------------------------------
// Run a battery of tests for NetCDF output and then input.
//
// Arguments:
//   
//-------------------------------------------------------------------------------------------------

//-------------------------------------------------------------------------------------------------
// Main function to run the FFT tests (both real and complex) and report results.
//
// Arguments:
//   argc: The number of command-line arguments.
//   argv: The array of command-line arguments (unused in this case).
//-------------------------------------------------------------------------------------------------
int main(int argc, const char* argv[]) {

  // Initialize STORMM test environment
  TestEnvironment oe(argc, argv);
  Xoshiro256ppGenerator rng(7183529);  // Using STORMM's random number generator

  // Section 1
  section("Running a real FFT test in PocketFFT");

  // Section 2
  section("Running a complex FFT in PocketFFT");

  // Section 3
  section("Test real-to-(half) complex 2D FFTs in PocketFFT");

  // Section 4
  section("Test real-to-(half) complex 3D FFTs in PocketFFT");

  // Section 5
  section("Test complex 2D FFTs in PocketFFT");
  
  // Section 6
  section("Test complex 3D FFTs in PocketFFT");
  
  // Section 7
  section("Test the FFT management system");

  // Section 8
  section("Test Basic NetCDF Functionality");

#ifdef STORMM_USE_HPC
  // Section 9
  section("Test cuFFT, double-precision");
  
  // Section 10
  section("Test cuFFT, single-precision");
#endif  
  
  // Begin testing
  section(1);
  testRealFFT(&rng);
  section(2);
  testComplexFFT(&rng);
  section(3);
  testReal2DFFT(&rng, { 8,  8 });
  testReal2DFFT(&rng, { 5,  6 });
  testReal2DFFT(&rng, { 4,  7 });
  section(4);
  testReal3DFFT(&rng, { 5,  6,  8 });
  testReal3DFFT(&rng, { 8,  8,  8 });
  testReal3DFFT(&rng, { 9,  9, 12 });
  testReal3DFFT(&rng, { 9, 12, 15 });
  section(5);
  testComplex2DFFT(&rng, { 8,  8 });
  testComplex2DFFT(&rng, { 5,  6 });
  section(6);
  testComplex3DFFT(&rng, { 8,  8,  8 });
  testComplex3DFFT(&rng, { 5,  6,  7 });
  testFFTEncapsulation<double, double2>(&rng, HybridTargetLevel::HOST, 1.0e-8);
  testFFTEncapsulation<float, float2>(&rng, HybridTargetLevel::HOST, 1.0e-5);
#ifdef STORMM_USE_HPC
  testFFTEncapsulation<double, double2>(&rng, HybridTargetLevel::DEVICE, 1.0e-8);
  testFFTEncapsulation<float, float2>(&rng, HybridTargetLevel::DEVICE, 1.0e-5);
#endif
#ifdef STORMM_USE_HPC
  section(9);
  testHpcFFT<double, double2>(&rng);
  section(10);
  testHpcFFT<float, float2>(&rng);
#endif
  
  // Print results
  printTestSummary(oe.getVerbosity());
  return countGlobalTestFailures();
}
