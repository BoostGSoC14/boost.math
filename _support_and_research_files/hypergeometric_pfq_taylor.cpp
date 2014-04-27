
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2014 Anton Bikineev
//  Copyright 2014 Christopher Kormanyos
//  Copyright 2014 John Maddock
//  Copyright 2014 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <cstdint>
#include <functional>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include <boost/multiprecision/cpp_dec_float.hpp>

template <class T,
          class forward_iterator_a_type,
          class forward_iterator_b_type>
inline T hypergeometric_pfq_imp(forward_iterator_a_type coefficients_a_begin,
                                forward_iterator_a_type coefficients_a_end,
                                forward_iterator_b_type coefficients_b_begin,
                                forward_iterator_b_type coefficients_b_end,
                                T x)
{
  // Compute the Taylor series expansion of
  // hypergeometric_pfq[{a0, a1, a2, ... an}; {b0, b1, b2, ... bm}, x].

  // TBD: There are no checks on input range or parameter boundaries.
  // But maybe we don't need any since this is a *detail* file.

  T x_pow_n_div_n_fact(x);

  // The pochhammer symbols for the multiplications in the series expansion
  // will be stored in STL-containers.
  std::vector<T> ap(coefficients_a_begin, coefficients_a_end);
  std::vector<T> bp(coefficients_b_begin, coefficients_b_end);

  const T my_one(1);

  // Initialize the pochhammer product terms with the products of the form:
  // [(a0)_1 * (a1)_1 * (a2)_1 * ...], or [(b0)_1 * (b1)_1 * (b2)_1 * ...].
  T pochhammer_sequence_a = std::accumulate(ap.begin(), ap.end(), my_one, std::multiplies<T>());
  T pochhammer_sequence_b = std::accumulate(bp.begin(), bp.end(), my_one, std::multiplies<T>());

  T hypergeometric_pfq_result = my_one + ((pochhammer_sequence_a / pochhammer_sequence_b) * x_pow_n_div_n_fact);

  int n;

  // This is the maximum number of iterations allowed.
  // TBD: Should this upper limit be scaled according to the number
  // of digits in the type T?
  const int max_iteration = 10000;

  for(n = 2; n < max_iteration; ++n)
  {
    x_pow_n_div_n_fact *= x;
    x_pow_n_div_n_fact /= n;

    // Increment each of the pochhammer elements in a and b.
    std::transform(ap.begin(), ap.end(), ap.begin(), std::bind1st(std::plus<T>(), my_one));
    std::transform(bp.begin(), bp.end(), bp.begin(), std::bind1st(std::plus<T>(), my_one));

    // Multiply the pochhammer product terms with the products of the incremented
    // pochhammer elements. These are products of the form:
    // [(a0)_k * (a1)_k * (a2)_k * ...], or [(b0)_k * (b1)_k * (b2)_k * ...].
    pochhammer_sequence_a *= std::accumulate(ap.begin(), ap.end(), my_one, std::multiplies<T>());
    pochhammer_sequence_b *= std::accumulate(bp.begin(), bp.end(), my_one, std::multiplies<T>());

    // TBD: Which algebraic order is better here for overflow / accuracy concerns?
    // Should we use?
    // (an * (x^n / n!)) / bn
    // Or rather use?
    // (an / bn) * (x^n / n!)
    const T next_term = (pochhammer_sequence_a * x_pow_n_div_n_fact) / pochhammer_sequence_b;

    using std::fabs;

    // TBD: Should we use more clever limiting here based on the digits in the type T?
    // TBD: Use boost:tools::max_value or something like that here instead of
    // numeric_limits (real-concept has no limits).
    if((n > 10) && fabs(next_term) < std::numeric_limits<T>::epsilon())
    {
      break;
    }

    hypergeometric_pfq_result += next_term;
  }

  // TBD: Use boost:tools::max_value or something like that here instead of
  // numeric_limits (real-concept has no limits).
  return ((n < max_iteration) ? hypergeometric_pfq_result : std::numeric_limits<T>::quiet_NaN());
}

namespace
{
//  typedef double float_type;
  typedef boost::multiprecision::cpp_dec_float_100 float_type;
}

int main()
{
  const std::vector<float_type> an( { float_type(1) / 2, float_type(1) / 3, float_type(1) / 4, float_type(1) / 5, } );
  const std::vector<float_type> bm( { float_type(2) / 3, float_type(2) / 4, float_type(2) / 5, float_type(2) / 6, float_type(2) / 7 } );

  const float_type x = float_type(1) / 7;

  const float_type h = hypergeometric_pfq_imp(an.begin(), an.end(), bm.begin(), bm.end(), x);

  // Control result from Wolfram's Alpha or Mathermaica(R).
  // N[HypergeometricPFQ[{1/2, 1/3, 1/4, 1/5}, {2/3, 2/4, 2/5, 2/6, 2/7}, 1/7], 100]
  // 1.097152657057488780105930458245903202876825351898026852284800644895088124044818091921532757454387577
  std::cout << std::setprecision(std::numeric_limits<float_type>::digits10)
            << h
            << std::endl;
}
