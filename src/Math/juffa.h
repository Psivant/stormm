// -*-c++-*-
#ifndef STORMM_JUFFA_H
#define STORMM_JUFFA_H

#include "copyright.h"

// Copyright (c) 2015-2021 Norbert Juffa
// All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without modification, are permitted
//  provided that the following conditions are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//
//  2. Redistributions in binary form must reproduce the above copyright notice, this list of
//     conditions and the following disclaimer in the documentation and/or other materials provided
//     with the distribution.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR
//  IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
//  AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
//  CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
//  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
//  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
//  THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
//  OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
//  POSSIBILITY OF SUCH DAMAGE.

namespace stormm {
namespace stmath {

/// \brief Implement Norbert Juffa's fast erfc(alpha * x) approximation for 32-bit floating point
///        computations.  In the most prominent case for molecular dynamics, alpha corresponds to
///        the "Ewald coefficient."
class JfErfc {
public:

  /// \brief The constructor requires a pre-factor for the argument that the object will eventually
  ///        take in order to return erfc(alpha * x).
  JfErfc(float alpha_in = 1.0);

  /// \brief Evaluate the approximate erfc(alpha * x) using all 32-bit floating-point arithmetic.
  ///
  /// \param x  The function argument to evaluate
  float erfcJf(float x);
    
private:
  float alpha;  ///< The original prefactor applied to every argument processed by this object.
  float c0;     ///< Highest-order polynomial coefficient, containing a factor of alpha^10
  float c1;     ///< Second highest-order polynomial coefficient
  float c2;     ///< Third polynomial coefficient
  float c3;     ///< Fourth polynomial coefficient
  float c4;     ///< Fifth polynomial coefficient
  float c5;     ///< Sixth polynomial coefficient
  float c6;     ///< Seventh polynomial coefficient
  float c7;     ///< Eighth polynomial coefficient
  float c8;     ///< Ninth polynomial coefficient
  float c9;     ///< Tenth and lowest-order polynomial coefficient, containing a factor of alpha
};

/// \brief Implement Norbert Juffa's performant, more accurate expf().  This stands as a free
///        function, with no pre-computations necessary.
///
/// \param x  The exponential argument to evaluate
float expJf(float x);
  
} // namespace stmath
} // namespace stormm
#endif
