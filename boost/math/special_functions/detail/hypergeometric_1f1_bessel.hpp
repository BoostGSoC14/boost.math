
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2014 Anton Bikineev
//  Copyright 2014 Christopher Kormanyos
//  Copyright 2014 John Maddock
//  Copyright 2014 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
#ifndef BOOST_MATH_HYPERGEOMETRIC_1F1_BESSEL_HPP
  #define BOOST_MATH_HYPERGEOMETRIC_1F1_BESSEL_HPP

  namespace boost { namespace math { namespace detail {

  // helper definitions:
  template <class T>
  inline T hypergeometric_bessel_j_recurrence_next(const T& jvm2, const T& jvm1, const T& vm1, const T& z)
  {
    return (((2 * vm1) / z) * jvm1) - jvm2;
  }

  template <class T>
  inline void hypergeometric_bessel_j_recurrence_iterate(T& jvm2, T& jvm1, T& jv, const T& vm1, const T& z)
  {
    using std::swap;

    swap(jvm2, jvm1);
    swap(jvm1, jv);
    jv = detail::hypergeometric_bessel_j_recurrence_next(jvm2, jvm1, vm1, z);
  }

  template <class T>
  inline T hypergeometric_bessel_j_coefficient_next(const T& cnm3, const T& cnm2, const T& cnm1, const T& a, const T& b, const T& h, const unsigned n)
  {
    const T one_minus_two_h = 1 - (2 * h);
    const T h_minus_one = h - 1;

    const T term_cnm1  = ((one_minus_two_h * n) - (b * h)) * cnm1;
    const T term_cnm2 = ((one_minus_two_h * a) - ((h * h_minus_one) * (b + (n - 1)))) * cnm2;
    const T term_cnm3 = ((-h * h_minus_one) * a) * cnm3;

    return ((term_cnm1 + term_cnm2) + term_cnm3) / (n + 1);
  }

  template <class T>
  inline T hypergeometric_bessel_j_coefficient_iterate(T& cnm3, T& cnm2, T& cnm1, T& cn, const T& a, const T& b, const T& h, const unsigned n)
  {
    using std::swap;

    swap(cnm3, cnm2);
    swap(cnm2, cnm1);
    swap(cnm1, cn);
    cn = detail::hypergeometric_bessel_j_coefficient_next(cnm3, cnm2, cnm1, a, b, h, n);
  }

  // term class of Abramowitz & Stegun 13_3_8 formula
  template <class T, class Policy>
  struct hypergeometric_1f1_bessel_j_series_term
  {
    typedef T result_type;

    static const T h;

    hypergeometric_1f1_bessel_j_series_term(const T& a, const T& b, const T& z):
      a(a), b(b), z(z), n(0)
    {
      BOOST_MATH_STD_USING

      sqrt_minus_az = sqrt(-a * z);
      const T double_sqrt_minus_az = 2 * sqrt_minus_az;

      v_current = b - 1;
      z_pow_n = z; sqrt_minus_az_pow_n = sqrt_minus_az;

      cnm3 = 1;
      cnm2 = -b * h;
      cnm1 = (((1 - (2 * h)) * a) + ((b * (b + 1)) * (h * h))) / 2;
      cn = detail::hypergeometric_bessel_j_coefficient_next(cnm3, cnm2, cnm1, a, b, h, 3u);

      jvm2 = boost::math::cyl_bessel_j(v_current, double_sqrt_minus_az);
      jvm1 = boost::math::cyl_bessel_j(b, double_sqrt_minus_az);
      jv = detail::hypergeometric_bessel_j_recurrence_next(jvm2, jvm1, b, double_sqrt_minus_az);
      ++v_current;
    }

    T operator()()
    {
      ++n;
      switch (n - 1)
      {
        case 0u:
          return jvm2;
        case 1u:
          return ((cnm2 * z) / sqrt_minus_az) * jvm1;
      }

      ++v_current;
      z_pow_n *= z;
      sqrt_minus_az_pow_n *= sqrt_minus_az;
      const T result = ((cnm1 * z_pow_n) / sqrt_minus_az_pow_n) * jv;

      detail::hypergeometric_bessel_j_coefficient_iterate(cnm3, cnm2, cnm1, cn, a, b, h, n + 1);
      detail::hypergeometric_bessel_j_recurrence_iterate(jvm2, jvm1, jv, v_current, 2 * sqrt_minus_az);
      return result;
    }

  private:
    const T a, b, z;
    unsigned n;
    T sqrt_minus_az;
    T z_pow_n, sqrt_minus_az_pow_n;
    T v_current;
    T cnm3, cnm2, cnm1, cn;
    T jvm2, jvm1, jv;
  };

  template <class T, class Policy>
  const T hypergeometric_1f1_bessel_j_series_term<T, Policy>::h = -boost::math::constants::pi<T>() / T(10.);

  // "higher level" function
  template <class T, class Policy>
  inline T hypergeometric_1f1_bessel_j_series(const T& a, const T& b, const T& z, const Policy& pol)
  {
    BOOST_MATH_STD_USING

    const T sqrt_minus_az_pow_b_minus_one = pow(sqrt(-a * z), b - 1);
    const T prefix = (boost::math::tgamma(b) / sqrt_minus_az_pow_b_minus_one) *
                      exp(detail::hypergeometric_1f1_bessel_j_series_term<T, Policy>::h * z);

    detail::hypergeometric_1f1_bessel_j_series_term<T, Policy> s(a, b, z);
    boost::uintmax_t max_iter = boost::math::policies::get_max_series_iterations<Policy>();
#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
    T zero = 0;
    T result = boost::math::tools::sum_series(s, boost::math::policies::get_epsilon<T, Policy>(), max_iter, zero);
#else
    T result = boost::math::tools::sum_series(s, boost::math::policies::get_epsilon<T, Policy>(), max_iter);
#endif
    boost::math::policies::check_series_iterations<T>("boost::math::hypergeometric_1f1_bessel_j_series<%1%>(%1%,%1%,%1%)", max_iter, pol);
    return prefix * result;
  }

  } } } // namespaces

#endif // BOOST_MATH_HYPERGEOMETRIC_1F1_BESSEL_HPP
