
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2014 Anton Bikineev
//  Copyright 2014 Christopher Kormanyos
//  Copyright 2014 John Maddock
//  Copyright 2014 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <vector>

template <class T,
          class forward_iterator_a_type,
          class forward_iterator_b_type>
inline T hypergeometric_pfq_imp(forward_iterator_a_type coefficients_a_begin,
                                forward_iterator_a_type coefficients_a_end,
                                forward_iterator_b_type coefficients_b_begin,
                                forward_iterator_b_type coefficients_b_end,
                                T x)
{
  // TBD: Chris: Implement the small-argument Taylor series for
  // hypergeometric_pfq(a0, a1, ... an ; b0, b1, ... bm, x).

  static_cast<void>(coefficients_a_begin);
  static_cast<void>(coefficients_a_end);
  static_cast<void>(coefficients_b_begin);
  static_cast<void>(coefficients_b_end);
  static_cast<void>(x);

  return T();
}

namespace
{
  typedef double float_type;
}

int main()
{
  const std::vector<float_type> an( { float_type(1) / 2, float_type(1) / 3, float_type(1) / 4, float_type(1) / 5, } );
  const std::vector<float_type> bm( { float_type(2) / 3, float_type(2) / 4, float_type(2) / 5, float_type(2) / 6, float_type(2) / 7 } );

  const float_type x = float_type(1) / 7;

  const float_type h = hypergeometric_pfq_imp(an.begin(), an.end(), bm.begin(), bm.end(), x);

  std::cout << std::setprecision(std::numeric_limits<float_type>::digits10)
            << h
            << std::endl;
}
