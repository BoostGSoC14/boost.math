
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2014 Anton Bikineev
//  Copyright 2014 Christopher Kormanyos
//  Copyright 2014 John Maddock
//  Copyright 2014 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
#ifndef BOOST_MATH_HYPERGEOMETRIC_0F1_BESSEL_HPP
  #define BOOST_MATH_HYPERGEOMETRIC_0F1_BESSEL_HPP

  namespace boost { namespace math { namespace detail {

  template <class T, class Policy>
  inline T hypergeometric_0f1_bessel(const T& b, const T& z, const Policy& pol)
  {
    BOOST_MATH_STD_USING

    const bool is_z_nonpositive = z <= 0;

    const T sqrt_z = is_z_nonpositive ? T(sqrt(-z)) : T(sqrt(z));
    const T sqrt_z_pow_b = pow(sqrt_z, b);
    const T bessel_mult = is_z_nonpositive ?
      boost::math::cyl_bessel_j(b - 1, 2 * sqrt_z, pol) :
      boost::math::cyl_bessel_i(b - 1, 2 * sqrt_z, pol) ;

    return ((boost::math::tgamma(b, pol) * sqrt_z) / sqrt_z_pow_b) * bessel_mult;
  }

  } } } // namespaces

#endif // BOOST_MATH_HYPERGEOMETRIC_0F1_BESSEL_HPP
