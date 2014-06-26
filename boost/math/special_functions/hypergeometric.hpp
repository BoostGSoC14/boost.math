
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2014 Anton Bikineev
//  Copyright 2014 Christopher Kormanyos
//  Copyright 2014 John Maddock
//  Copyright 2014 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef _BOOST_HYPERGEOMETRIC_2014_04_07_HPP_
  #define _BOOST_HYPERGEOMETRIC_2014_04_07_HPP_

  #include <boost/math/special_functions/detail/hypergeometric_series.hpp>
  #include <boost/math/special_functions/detail/hypergeometric_0f1_bessel.hpp>
  #include <boost/math/special_functions/detail/hypergeometric_asym.hpp>

  namespace boost { namespace math { namespace detail {

  // check when 1f1 series can't decay to polynom
  template <class T>
  inline bool check_hypergeometric_1f1_parameters(const T& a, const T& b)
  {
    BOOST_MATH_STD_USING

    if (b <= 0 && b == floor(b))
    {
      if (a >= 0)
        return false;

      if (a < b && a == floor(a))
        return false;
    }

    return true;
  }

  template <class T, class Policy>
  inline T hypergeometric_0f1_imp(T b, T z, const Policy& pol)
  {
    BOOST_MATH_STD_USING

    // some special cases
    if (z == 0)
      return T(1);

    if (b <= 0 && b == floor(b))
      return policies::raise_pole_error<T>(
        "boost::math::hypergeometric_0f1<%1%,%1%>(%1%, %1%)",
        "Evaluation of 0f1 with nonpositive integer b = %1%.", b, pol);

    // evaluation through Taylor series looks
    // more precisious than Bessel relation:
    // detail::hypergeometric_0f1_bessel(b, z, pol);
    return detail::hypergeometric_0f1_generic_series(b, z, pol);
  }

  template <class T, class Policy>
  inline T hypergeometric_1f0_imp(T a, T z, const Policy& pol)
  {
    BOOST_MATH_STD_USING // pow

    if (z == 1)
      return policies::raise_pole_error<T>(
        "boost::math::hypergeometric_1f0<%1%,%1%>(%1%, %1%)",
        "Evaluation of 1f0 with z = %1%.",
        z,
        pol);

    // more naive and convergent method than series
    return pow(1 - z, -a);
  }

  template <class T, class Policy>
  inline T hypergeometric_1f1_imp(T a, T b, T z, const Policy& pol)
  {
    BOOST_MATH_STD_USING // exp, fabs, sqrt

    static const char* const function = "boost::math::hypergeometric_1f1<%1%,%1%,%1%>(%1%,%1%,%1%)";

    if (z == 0 || a == 0)
      return T(1);

    // undefined result:
    if (!detail::check_hypergeometric_1f1_parameters(a, b))
      return policies::raise_domain_error<T>(
        function,
        "Function is indeterminate for negative integer b = %1%.",
        b,
        pol);

    // other checks:
    if (a == -1)
      return 1 - (z / b);

    const T b_minus_a = b - a;

    // 0f0 (exp) case;
    if (b_minus_a == 0)
      return exp(z);

    if (b_minus_a == -1)
      return (1 + (z / b)) * exp(z);

    // check for poles in gamma
    if ((a > 0 || a != floor(a)) &&
        (b > 0 || b != floor(b)))
    {
      // asymp expansion
      if (detail::hypergeometric_1f1_asym_region(a, b, z))
        return boost::math::detail::hypergeometric_1f1_asym_series(a, b, z, pol);
    }

    // kummer's transformation
    if (z < 0)
      return exp(z) * detail::hypergeometric_1f1_imp<T>(b_minus_a, b, -z, pol);

    return detail::hypergeometric_1f1_generic_series(a, b, z, pol);
  }

  template <class T, class Policy>
  inline T hypergeometric_1f2_imp(T a, T b1, T b2, T z, const Policy& pol)
  {
    static const char* const function = "boost::math::hypergeometric_1f2<%1%,%1%,%1%,%1%>(%1%,%1%,%1%,%1%)";

    if (z == 0 || a == 0)
      return T(1);

    // check for parameter equality
    if (a == b1 || a == b2)
    {
      const T b = (a == b1) ? b2 : b1; 
      return detail::hypergeometric_0f1_imp(b, z, pol);
    }

    // undefined result:
    if (!detail::check_hypergeometric_1f1_parameters(a, b1) ||
        !detail::check_hypergeometric_1f1_parameters(a, b2))
    {
      return policies::raise_domain_error<T>(
        function,
        "Function is indeterminate for negative integer b = %1%.",
        b1 <= b2 ? b1 : b2,
        pol);
    }

    return detail::hypergeometric_1f2_generic_series(a, b1, b2, z, pol);
  }

  template <class T, class Policy>
  inline T hypergeometric_2f1_imp(T a1, T a2, T b, T z, const Policy& pol)
  {
    static const char* const function = "boost::math::hypergeometric_2f1<%1%,%1%,%1%,%1%>(%1%,%1%,%1%,%1%)";

    // undefined result:
    if (!detail::check_hypergeometric_1f1_parameters(a1, b) ||
        !detail::check_hypergeometric_1f1_parameters(a2, b))
    {
      return policies::raise_domain_error<T>(
        function,
        "Function is indeterminate for negative integer b = %1%.",
        b,
        pol);
    }

    return detail::hypergeometric_2f1_generic_series(a1, a2, b, z, pol);
  }

  } // namespace detail

  template <class T1, class T2, class Policy>
  inline typename tools::promote_args<T1, T2>::type hypergeometric_0f1(T1 b, T2 z, const Policy& /* pol */)
  {
    BOOST_FPU_EXCEPTION_GUARD
    typedef typename tools::promote_args<T1, T2>::type result_type;
    typedef typename policies::evaluation<result_type, Policy>::type value_type;
    typedef typename policies::normalise<
       Policy,
       policies::promote_float<false>,
       policies::promote_double<false>,
       policies::discrete_quantile<>,
       policies::assert_undefined<> >::type forwarding_policy;
    return policies::checked_narrowing_cast<result_type, Policy>(
          detail::hypergeometric_0f1_imp<value_type>(
                static_cast<value_type>(b),
                static_cast<value_type>(z),
                forwarding_policy()),
          "boost::math::hypergeometric_0f1<%1%>(%1%,%1%)");
  }

  template <class T1, class T2>
  inline typename tools::promote_args<T1, T2>::type hypergeometric_0f1(T1 b, T2 z)
  {
    return hypergeometric_0f1(b, z, policies::policy<>());
  }

  template <class T1, class T2, class Policy>
  inline typename tools::promote_args<T1, T2>::type hypergeometric_1f0(T1 a, T2 z, const Policy& /* pol */)
  {
    BOOST_FPU_EXCEPTION_GUARD
    typedef typename tools::promote_args<T1, T2>::type result_type;
    typedef typename policies::evaluation<result_type, Policy>::type value_type;
    typedef typename policies::normalise<
       Policy,
       policies::promote_float<false>,
       policies::promote_double<false>,
       policies::discrete_quantile<>,
       policies::assert_undefined<> >::type forwarding_policy;
    return policies::checked_narrowing_cast<result_type, Policy>(
          detail::hypergeometric_1f0_imp<value_type>(
                static_cast<value_type>(a),
                static_cast<value_type>(z),
                forwarding_policy()),
          "boost::math::hypergeometric_1f0<%1%>(%1%,%1%)");
  }

  template <class T1, class T2>
  inline typename tools::promote_args<T1, T2>::type hypergeometric_1f0(T1 a, T2 z)
  {
    return hypergeometric_1f0(a, z, policies::policy<>());
  }

  template <class T1, class T2, class T3, class Policy>
  inline typename tools::promote_args<T1, T2, T3>::type hypergeometric_1f1(T1 a, T2 b, T3 z, const Policy& /* pol */)
  {
    BOOST_FPU_EXCEPTION_GUARD
    typedef typename tools::promote_args<T1, T2, T3>::type result_type;
    typedef typename policies::evaluation<result_type, Policy>::type value_type;
    typedef typename policies::normalise<
       Policy,
       policies::promote_float<false>,
       policies::promote_double<false>,
       policies::discrete_quantile<>,
       policies::assert_undefined<> >::type forwarding_policy;
    return policies::checked_narrowing_cast<result_type, Policy>(
          detail::hypergeometric_1f1_imp<value_type>(
                static_cast<value_type>(a),
                static_cast<value_type>(b),
                static_cast<value_type>(z),
                forwarding_policy()),
          "boost::math::hypergeometric_1f1<%1%>(%1%,%1%,%1%)");
  }

  template <class T1, class T2, class T3>
  inline typename tools::promote_args<T1, T2, T3>::type hypergeometric_1f1(T1 a, T2 b, T3 z)
  {
    return hypergeometric_1f1(a, b, z, policies::policy<>());
  }

  template <class T1, class T2, class T3, class T4, class Policy>
  inline typename tools::promote_args<T1, T2, T3, T4>::type hypergeometric_1f2(T1 a, T2 b1, T3 b2, T4 z, const Policy& /* pol */)
  {
    BOOST_FPU_EXCEPTION_GUARD
    typedef typename tools::promote_args<T1, T2, T3, T4>::type result_type;
    typedef typename policies::evaluation<result_type, Policy>::type value_type;
    typedef typename policies::normalise<
       Policy,
       policies::promote_float<false>,
       policies::promote_double<false>,
       policies::discrete_quantile<>,
       policies::assert_undefined<> >::type forwarding_policy;
    return policies::checked_narrowing_cast<result_type, Policy>(
          detail::hypergeometric_1f2_imp<value_type>(
                static_cast<value_type>(a),
                static_cast<value_type>(b1),
                static_cast<value_type>(b2),
                static_cast<value_type>(z),
                forwarding_policy()),
          "boost::math::hypergeometric_1f2<%1%>(%1%,%1%,%1%,%1%)");
  }

  template <class T1, class T2, class T3, class T4>
  inline typename tools::promote_args<T1, T2, T3, T4>::type hypergeometric_1f2(T1 a, T2 b1, T3 b2, T4 z)
  {
    return hypergeometric_1f2(a, b1, b2, z, policies::policy<>());
  }

  template <class T1, class T2, class T3, class T4, class Policy>
  inline typename tools::promote_args<T1, T2, T3, T4>::type hypergeometric_2f1(T1 a1, T2 a2, T3 b, T4 z, const Policy& /* pol */)
  {
    BOOST_FPU_EXCEPTION_GUARD
    typedef typename tools::promote_args<T1, T2, T3, T4>::type result_type;
    typedef typename policies::evaluation<result_type, Policy>::type value_type;
    typedef typename policies::normalise<
       Policy,
       policies::promote_float<false>,
       policies::promote_double<false>,
       policies::discrete_quantile<>,
       policies::assert_undefined<> >::type forwarding_policy;
    return policies::checked_narrowing_cast<result_type, Policy>(
          detail::hypergeometric_2f1_imp<value_type>(
                static_cast<value_type>(a1),
                static_cast<value_type>(a2),
                static_cast<value_type>(b),
                static_cast<value_type>(z),
                forwarding_policy()),
          "boost::math::hypergeometric_2f1<%1%>(%1%,%1%,%1%,%1%)");
  }

  template <class T1, class T2, class T3, class T4>
  inline typename tools::promote_args<T1, T2, T3, T4>::type hypergeometric_2f1(T1 a1, T2 a2, T3 b, T4 z)
  {
    return hypergeometric_2f1(a1, a2, b, z, policies::policy<>());
  }

  } } // namespace boost::math

#endif // _BOOST_HYPERGEOMETRIC_2014_04_07_HPP_
