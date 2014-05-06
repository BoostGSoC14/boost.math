
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2014 Anton Bikineev
//  Copyright 2014 Christopher Kormanyos
//  Copyright 2014 John Maddock
//  Copyright 2014 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
#ifndef BOOST_MATH_HYPERGEOMETRIC_ASYM_HPP
  #define BOOST_MATH_HYPERGEOMETRIC_ASYM_HPP

  namespace boost { namespace math { namespace detail {

  template <class T, class Policy>
  struct hypergeometric_1f1_asym_series_term_a
  {
    typedef T result_type;

    hypergeometric_1f1_asym_series_term_a(const T& a, const T& b, const T& z):
      n(0), b_sub_a(b - a), one_sub_a(1 - a), z(z), term(1)
    {
    }

    result_type operator()()
    {
      const T result = term;
      ++n;
      term *= (b_sub_a * one_sub_a) / (n * z);
      ++b_sub_a, ++one_sub_a;
      return result;
    }

  private:
    unsigned n;
    T b_sub_a;
    T one_sub_a;
    T z;
    T term;
  };

  template <class T, class Policy>
  struct hypergeometric_1f1_asym_series_term_b
  {
    typedef T result_type;

    hypergeometric_1f1_asym_series_term_b(const T& a, const T& b, const T& z):
      n(0), a(a), one_plus_a_sub_b(1 + a - b), z(z), term(1)
    {
    }

    result_type operator()()
    {
      const T result = term;
      ++n;
      term *= -(a * one_plus_a_sub_b) / (n * z);
      ++a; ++one_plus_a_sub_b;
      return result;
    }

  private:
    unsigned n;
    T a;
    T one_plus_a_sub_b;
    T z;
    T term;
  };

  template <class T, class Policy>
  inline T hypergeometric_1f1_asym_series(const T& a, const T& b, const T& z, const Policy& pol)
  {
    BOOST_MATH_STD_USING
    static const char* const function = "boost::math::hypergeometric_1f1_asym_series<%1%>(%1%,%1%,%1%)";

    const T prefix_a = (exp(z) * boost::math::tgamma_ratio(b, a, pol)) * pow(z, (a - b));
    hypergeometric_1f1_asym_series_term_a<T, Policy> s_a(a, b, z);
    boost::uintmax_t max_iter = policies::get_max_series_iterations<Policy>();
#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
    const T zero = 0;
    T result_a = boost::math::tools::sum_series(s_a, boost::math::policies::get_epsilon<T, Policy>(), max_iter, zero);
#else
    T result_a = boost::math::tools::sum_series(s_a, boost::math::policies::get_epsilon<T, Policy>(), max_iter);
#endif
    policies::check_series_iterations<T>(function, max_iter, pol);
    result_a *= prefix_a;

    const T prefix_b = boost::math::cos_pi(a) /
        (boost::math::tgamma((b - a), pol) * pow(z, a));
    hypergeometric_1f1_asym_series_term_b<T, Policy> s_b(a, b, z);
    max_iter = policies::get_max_series_iterations<Policy>();
#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
    const T zero = 0;
    T result_b  = boost::math::tools::sum_series(s_b, boost::math::policies::get_epsilon<T, Policy>(), max_iter, zero);
#else
    T result_b = boost::math::tools::sum_series(s_b, boost::math::policies::get_epsilon<T, Policy>(), max_iter);
#endif
    policies::check_series_iterations<T>(function, max_iter, pol);
    result_b *= prefix_b;

    return result_a + result_b;
  }

  } } } // namespaces

#endif // BOOST_MATH_HYPERGEOMETRIC_ASYM_HPP


