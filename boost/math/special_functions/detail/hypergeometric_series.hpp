
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2014 Anton Bikineev
//  Copyright 2014 Christopher Kormanyos
//  Copyright 2014 John Maddock
//  Copyright 2014 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
#ifndef BOOST_MATH_HYPERGEOMETRIC_SERIES_HPP
  #define BOOST_MATH_HYPERGEOMETRIC_SERIES_HPP

  namespace boost { namespace math { namespace detail {

  // primary template for term of Taylor series
  template <class T, unsigned p, unsigned q, class Policy>
  struct hypergeometric_pfq_generic_series_term;

  // partial specialization for 0F1
  template <class T, class Policy>
  struct hypergeometric_pfq_generic_series_term<T, 0u, 1u, Policy>
  {
    typedef T result_type;

    hypergeometric_pfq_generic_series_term(T b, T z)
       : N(0), term(1), b(b), z(z)
    {
    }

    T operator()()
    {
      BOOST_MATH_STD_USING
      const T r = term;
      term *= ((1 / ((b + N) * (N + 1))) * z);
      ++N;
      return r;
    }

  private:
    unsigned N;
    T term;
    const T b, z;
  };

  // partial specialization for 1F0
  template <class T, class Policy>
  struct hypergeometric_pfq_generic_series_term<T, 1u, 0u, Policy>
  {
    typedef T result_type;

    hypergeometric_pfq_generic_series_term(T a, T z)
       : N(0), term(1), a(a), z(z)
    {
    }

    T operator()()
    {
      BOOST_MATH_STD_USING
      const T r = term;
      term *= (((a + N) / (N + 1)) * z);
      ++N;
      return r;
    }

  private:
    unsigned N;
    T term;
    const T a, z;
  };

  // partial specialization for 1F1
  template <class T, class Policy>
  struct hypergeometric_pfq_generic_series_term<T, 1u, 1u, Policy>
  {
    typedef T result_type;

    hypergeometric_pfq_generic_series_term(T a, T b, T z)
       : N(0), term(1), a(a), b(b), z(z)
    {
    }

    T operator()()
    {
      BOOST_MATH_STD_USING
      const T r = term;
      term *= (((a + N) / ((b + N) * (N + 1))) * z);
      ++N;
      return r;
    }

  private:
    unsigned N;
    T term;
    const T a, b, z;
  };

  // partial specialization for 1F2
  template <class T, class Policy>
  struct hypergeometric_pfq_generic_series_term<T, 1u, 2u, Policy>
  {
    typedef T result_type;

    hypergeometric_pfq_generic_series_term(T a, T b1, T b2, T z)
       : N(0), term(1), a(a), b1(b1), b2(b2), z(z)
    {
    }

    T operator()()
    {
      BOOST_MATH_STD_USING
      const T r = term;
      term *= (((a + N) / ((b1 + N) * (b2 + N) * (N + 1))) * z);
      ++N;
      return r;
    }

  private:
    unsigned N;
    T term;
    const T a, b1, b2, z;
  };

  // partial specialization for 2F1
  template <class T, class Policy>
  struct hypergeometric_pfq_generic_series_term<T, 2u, 1u, Policy>
  {
    typedef T result_type;

    hypergeometric_pfq_generic_series_term(T a1, T a2, T b, T z)
       : N(0), term(1), a1(a1), a2(a2), b(b), z(z)
    {
    }

    T operator()()
    {
      BOOST_MATH_STD_USING
      const T r = term;
      term *= (((a1 + N) * (a2 + N) / ((b + N) * (N + 1))) * z);
      ++N;
      return r;
    }

  private:
    unsigned N;
    T term;
    const T a1, a2, b, z;
  };

  template <class T, class Policy, unsigned p, unsigned q>
  inline T sum_pfq_series(detail::hypergeometric_pfq_generic_series_term<T, p, q, Policy>& term, const Policy& pol)
  {
    BOOST_MATH_STD_USING
    boost::uintmax_t max_iter = policies::get_max_series_iterations<Policy>();
#if BOOST_WORKAROUND(__BORLANDC__, BOOST_TESTED_AT(0x582))
    const T zero = 0;
    const T result = boost::math::tools::sum_series(term, boost::math::policies::get_epsilon<T, Policy>(), max_iter, zero);
#else
    const T result = boost::math::tools::sum_series(term, boost::math::policies::get_epsilon<T, Policy>(), max_iter);
#endif
    policies::check_series_iterations<T>("boost::math::hypergeometric_pfq_generic_series<%1%>(%1%,%1%,%1%)", max_iter, pol);
    return result;
  }

  template <class T, class Policy>
  inline T hypergeometric_0f1_generic_series(T b, T z, const Policy& pol)
  {
    detail::hypergeometric_pfq_generic_series_term<T, 0u, 1u, Policy> s(b, z);
    return detail::sum_pfq_series(s, pol);
  }

  template <class T, class Policy>
  inline T hypergeometric_1f0_generic_series(T a, T z, const Policy& pol)
  {
    detail::hypergeometric_pfq_generic_series_term<T, 1u, 0u, Policy> s(a, z);
    return detail::sum_pfq_series(s, pol);
  }

  template <class T, class Policy>
  inline T hypergeometric_1f1_generic_series(T a, T b, T z, const Policy& pol)
  {
    detail::hypergeometric_pfq_generic_series_term<T, 1u, 1u, Policy> s(a, b, z);
    return detail::sum_pfq_series(s, pol);
  }

  template <class T, class Policy>
  inline T hypergeometric_1f2_generic_series(T a, T b1, T b2, T z, const Policy& pol)
  {
    detail::hypergeometric_pfq_generic_series_term<T, 1u, 2u, Policy> s(a, b1, b2, z);
    return detail::sum_pfq_series(s, pol);
  }

  template <class T, class Policy>
  inline T hypergeometric_2f1_generic_series(T a1, T a2, T b, T z, const Policy& pol)
  {
    detail::hypergeometric_pfq_generic_series_term<T, 2u, 1u, Policy> s(a1, a2, b, z);
    return detail::sum_pfq_series(s, pol);
  }

  } } } // namespaces

#endif // BOOST_MATH_HYPERGEOMETRIC_SERIES_HPP
