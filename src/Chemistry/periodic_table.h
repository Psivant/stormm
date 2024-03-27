// -*-c++-*-
#ifndef STORMM_PERIODIC_TABLE_H
#define STORMM_PERIODIC_TABLE_H

#include <array>
#include "copyright.h"
#include "DataTypes/stormm_vector_types.h"

namespace stormm {
namespace chemistry {

constexpr int element_maximum_count = 96;

/// \brief Atomic masses, weighted averages of the natural abundance, for the first 96 elements
///        plus a virtual site at the front, to take the 0th position.  The Z-number is thus
///        implied by the position in the array.
constexpr std::array<double, element_maximum_count + 1> elemental_masses = {
      0.0000,
      1.0080,   4.0026,   6.9410,   9.0122,  10.8100,  12.0110,  14.0067,  15.9994,
     18.9984,  20.1790,  22.9898,  24.3050,  26.9815,  28.0855,  30.9738,  32.0600,
     35.4530,  39.9480,  39.0983,  40.0800,  44.9559,  47.9000,  50.9415,  51.9960,
     54.9380,  55.8470,  58.9332,  58.7000,  63.5460,  65.3800,  69.7200,  72.5900,
     74.9216,  78.9600,  79.9040,  83.8000,  85.4678,  87.6200,  88.9059,  91.2200,
     92.9064,  95.9400,  98.0000, 101.0700, 102.9055, 106.4000, 107.8680, 112.4100,
    114.8200, 118.6900, 121.7500, 127.6000, 126.9045, 131.3000, 132.9054, 137.3300,
    138.9055, 140.1200, 140.9077, 144.2400, 145.0000, 150.4000, 151.9600, 157.2500,
    158.9254, 162.5000, 164.9304, 167.2600, 168.9342, 173.0400, 174.9670, 178.4900,
    180.9479, 183.8500, 186.2070, 190.2000, 192.2200, 195.0900, 196.9665, 200.5900,
    204.3700, 207.2000, 208.9804, 209.0000, 210.0000, 222.0000, 223.0000, 226.0254,
    227.0278, 232.0381, 231.0359, 238.0290, 237.0482, 242.0000, 243.0000, 247.0000
};

/// \brief Atomic symbols as they appear in the periodic table.  The ith element of the array is
///        the symbol of the ith element, with VS being a virtual site and a placeholder in the
///        0th position (VS is deliberately two capital letters to break from the capital -
///        lowercase convention for all natural elements).
constexpr char2 element_VS = {'V', 'S'};
constexpr char2 element_H  = {'H', ' '};
constexpr char2 element_He = {'H', 'e'};
constexpr char2 element_Li = {'L', 'i'};
constexpr char2 element_Be = {'B', 'e'};
constexpr char2 element_B  = {'B', ' '};
constexpr char2 element_C  = {'C', ' '};
constexpr char2 element_N  = {'N', ' '};
constexpr char2 element_O  = {'O', ' '};
constexpr char2 element_F  = {'F', ' '};
constexpr char2 element_Ne = {'N', 'e'};
constexpr char2 element_Na = {'N', 'a'};
constexpr char2 element_Mg = {'M', 'g'};
constexpr char2 element_Al = {'A', 'l'};
constexpr char2 element_Si = {'S', 'i'};
constexpr char2 element_P  = {'P', ' '};
constexpr char2 element_S  = {'S', ' '};
constexpr char2 element_Cl = {'C', 'l'};
constexpr char2 element_Ar = {'A', 'r'};
constexpr char2 element_K  = {'K', ' '};
constexpr char2 element_Ca = {'C', 'a'};
constexpr char2 element_Sc = {'S', 'c'};
constexpr char2 element_Ti = {'T', 'i'};
constexpr char2 element_V  = {'V', ' '};
constexpr char2 element_Cr = {'C', 'r'};
constexpr char2 element_Mn = {'M', 'n'};
constexpr char2 element_Fe = {'F', 'e'};
constexpr char2 element_Co = {'C', 'o'};
constexpr char2 element_Ni = {'N', 'i'};
constexpr char2 element_Cu = {'C', 'u'};
constexpr char2 element_Zn = {'Z', 'n'};
constexpr char2 element_Ga = {'G', 'a'};
constexpr char2 element_Ge = {'G', 'e'};
constexpr char2 element_As = {'A', 's'};
constexpr char2 element_Se = {'S', 'e'};
constexpr char2 element_Br = {'B', 'r'};
constexpr char2 element_Kr = {'K', 'r'};
constexpr char2 element_Rb = {'R', 'b'};
constexpr char2 element_Sr = {'S', 'r'};
constexpr char2 element_Y  = {'Y', ' '};
constexpr char2 element_Zr = {'Z', 'r'};
constexpr char2 element_Nb = {'N', 'b'};
constexpr char2 element_Mo = {'M', 'o'};
constexpr char2 element_Tc = {'T', 'c'};
constexpr char2 element_Ru = {'R', 'u'};
constexpr char2 element_Rh = {'R', 'h'};
constexpr char2 element_Pd = {'P', 'd'};
constexpr char2 element_Ag = {'A', 'g'};
constexpr char2 element_Cd = {'C', 'd'};
constexpr char2 element_In = {'I', 'n'};
constexpr char2 element_Sn = {'S', 'n'};
constexpr char2 element_Sb = {'S', 'b'};
constexpr char2 element_Te = {'T', 'e'};
constexpr char2 element_I  = {'I', ' '};
constexpr char2 element_Xe = {'X', 'e'};
constexpr char2 element_Cs = {'C', 's'};
constexpr char2 element_Ba = {'B', 'a'};
constexpr char2 element_La = {'L', 'a'};
constexpr char2 element_Ce = {'C', 'e'};
constexpr char2 element_Pr = {'P', 'r'};
constexpr char2 element_Nd = {'N', 'd'};
constexpr char2 element_Pm = {'P', 'm'};
constexpr char2 element_Sm = {'S', 'm'};
constexpr char2 element_Eu = {'E', 'u'};
constexpr char2 element_Gd = {'G', 'd'};
constexpr char2 element_Tb = {'T', 'b'};
constexpr char2 element_Dy = {'D', 'y'};
constexpr char2 element_Ho = {'H', 'o'};
constexpr char2 element_Er = {'E', 'r'};
constexpr char2 element_Tm = {'T', 'm'};
constexpr char2 element_Yb = {'Y', 'b'};
constexpr char2 element_Lu = {'L', 'u'};
constexpr char2 element_Hf = {'H', 'f'};
constexpr char2 element_Ta = {'T', 'a'};
constexpr char2 element_W  = {'W', ' '};
constexpr char2 element_Re = {'R', 'e'};
constexpr char2 element_Os = {'O', 's'};
constexpr char2 element_Ir = {'I', 'r'};
constexpr char2 element_Pt = {'P', 't'};
constexpr char2 element_Au = {'A', 'u'};
constexpr char2 element_Hg = {'H', 'g'};
constexpr char2 element_Tl = {'T', 'l'};
constexpr char2 element_Pb = {'P', 'b'};
constexpr char2 element_Bi = {'B', 'i'};
constexpr char2 element_Po = {'P', 'o'};
constexpr char2 element_At = {'A', 't'};
constexpr char2 element_Rn = {'R', 'n'};
constexpr char2 element_Fr = {'F', 'r'};
constexpr char2 element_Ra = {'R', 'a'};
constexpr char2 element_Ac = {'A', 'c'};
constexpr char2 element_Th = {'T', 'h'};
constexpr char2 element_Pa = {'P', 'a'};
constexpr char2 element_U  = {'U', ' '};
constexpr char2 element_Np = {'N', 'p'};
constexpr char2 element_Pu = {'P', 'u'};
constexpr char2 element_Am = {'A', 'm'};
constexpr char2 element_Cm = {'C', 'm'};
constexpr std::array<char2, element_maximum_count + 1> elemental_symbols = {
  element_VS,
  element_H , element_He, element_Li, element_Be, element_B , element_C , element_N , element_O ,
  element_F , element_Ne, element_Na, element_Mg, element_Al, element_Si, element_P , element_S ,
  element_Cl, element_Ar, element_K , element_Ca, element_Sc, element_Ti, element_V , element_Cr,
  element_Mn, element_Fe, element_Co, element_Ni, element_Cu, element_Zn, element_Ga, element_Ge,
  element_As, element_Se, element_Br, element_Kr, element_Rb, element_Sr, element_Y , element_Zr,
  element_Nb, element_Mo, element_Tc, element_Ru, element_Rh, element_Pd, element_Ag, element_Cd,
  element_In, element_Sn, element_Sb, element_Te, element_I , element_Xe, element_Cs, element_Ba,
  element_La, element_Ce, element_Pr, element_Nd, element_Pm, element_Sm, element_Eu, element_Gd,
  element_Tb, element_Dy, element_Ho, element_Er, element_Tm, element_Yb, element_Lu, element_Hf,
  element_Ta, element_W , element_Re, element_Os, element_Ir, element_Pt, element_Au, element_Hg,
  element_Tl, element_Pb, element_Bi, element_Po, element_At, element_Rn, element_Fr, element_Ra,
  element_Ac, element_Th, element_Pa, element_U , element_Np, element_Pu, element_Am, element_Cm
};

/// \brief Atomic periods (rows in the periodic table).  The virtual site again takes row 0.  The
///        Lanthanides are completely represented, but only the first eight Actinides are present
///        before half-lives become so short that they are irrelevant.
constexpr std::array<int, element_maximum_count + 1> elemental_periods = {
  0,
  1,                                                                                           1,
  2, 2,                                                                         2, 2, 2, 2, 2, 2,
  3, 3,                                                                         3, 3, 3, 3, 3, 3,
  4, 4,                                           4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
  5, 5,                                           5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5, 5,  
  6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6, 6,
  7, 7, 7, 7, 7, 7, 7, 7, 7, 7
};

/// \brief Atomic groups (columns in the periodic table).  The virtual site takes group 0, and the
///        groups are represented with integers rather than roman numerals.
constexpr std::array<int, element_maximum_count + 1> elemental_groups = {
   0,  1,                                                                 18,
       1,  2,                                         13, 14, 15, 16, 17, 18,
       1,  2,                                         13, 14, 15, 16, 17, 18,
       1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
       1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
       1,  2,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  3,  4,
                       5,  6,  7,  8,  9, 10, 11, 12, 13, 14, 15, 16, 17, 18,
       1,  2,  3,  3,  3,  3,  3,  3,  3,  3
};

} // namespace chemistry
} // namespace stormm

#endif
