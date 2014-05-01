
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2014 Anton Bikineev
//  Copyright 2014 Christopher Kormanyos
//  Copyright 2014 John Maddock
//  Copyright 2014 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#include <algorithm>
#include <array>
#include <functional>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <vector>

#include <boost/math/constants/constants.hpp>
#include <boost/math/tools/precision.hpp>
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

  // Calculate the first term in the Taylor series expansion.
  // Use either: (an * (x^n / n!)) / bn
  // or else use: (an / bn) * (x^n / n!)
  // based on whether or not (x^n / n!) > 1.
  const T first_term = ((x_pow_n_div_n_fact > 1)
                         ? T((pochhammer_sequence_a / pochhammer_sequence_b) * x_pow_n_div_n_fact)
                         : T((pochhammer_sequence_a * x_pow_n_div_n_fact) / pochhammer_sequence_b));

  T hypergeometric_pfq_result = my_one + first_term;

  int n;

  // Calculate the maximum number of iterations allowed.
  // Here we use an expression similar to (std::numeric_limits<T>::digits10 * 10).
  const int max_iteration = static_cast<int>(static_cast<float>(boost::math::tools::digits<T>()) * 3.01F);

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

    // Calculate the next term in the Taylor series expansion.
    // Use either: (an * (x^n / n!)) / bn
    // or else use: (an / bn) * (x^n / n!)
    // based on whether or not (x^n / n!) > 1.
    const T next_term = ((x_pow_n_div_n_fact > 1)
                          ? T((pochhammer_sequence_a / pochhammer_sequence_b) * x_pow_n_div_n_fact)
                          : T((pochhammer_sequence_a * x_pow_n_div_n_fact) / pochhammer_sequence_b));

    using std::fabs;

    // TBD: Should we use a more clever limiting here (other than n > 10)
    // based on the digits in the type T?
    if((n > 10) && (fabs(next_term) < boost::math::tools::epsilon<T>()))
    {
      break;
    }

    hypergeometric_pfq_result += next_term;
  }

  if(n < max_iteration)
  {
    return hypergeometric_pfq_result;
  }
  else
  {
    // TBD: Do we need to throw an exception here.
    // TBD: Is T() the right return value here?
    return T();
  }
}

namespace
{
//  typedef double float_type;
  typedef boost::multiprecision::cpp_dec_float_100 float_type;
}

int main()
{
  const std::array<float_type, 4U> an = { float_type(1) / 2, float_type(1) / 3, float_type(1) / 4, float_type(1) / 5, };
  const std::array<float_type, 5U> bm = { float_type(2) / 3, float_type(2) / 4, float_type(2) / 5, float_type(2) / 6, float_type(2) / 7 };

  const float_type h = hypergeometric_pfq_imp(an.begin(),
                                              an.end(),
                                              bm.begin(),
                                              bm.end(),
                                              boost::math::constants::euler<float_type>());

  // Here is a control value from Wolfram's Alpha or Mathematica(R).
  // N[HypergeometricPFQ[{1/2, 1/3, 1/4, 1/5}, {2/3, 2/4, 2/5, 2/6, 2/7}, EulerGamma], 100]
  // 1.437152091623117098817180937046270756251132185487659323159061684966332133966272470711486705986290248

  std::cout << std::setprecision(std::numeric_limits<float_type>::digits10)
            << h
            << std::endl;
}
