
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

  namespace boost { namespace math { namespace detail {

  template <class T, class Policy>
  inline T hypergeometric_0f1_imp(T b, T z, const Policy& pol)
  {
    BOOST_MATH_STD_USING

    // some special cases
    if (z == 0)
      return 1;

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
    BOOST_MATH_STD_USING

    // some special cases
    //if (fabs(z - 1) >= 1)
    //{
      //std::cout << "special case\n";
      //return pow(1 - z, -a);
    //}

    return detail::hypergeometric_1f0_generic_series(a, z, pol);
  }

  template <class T, class Policy>
  inline T hypergeometric_1f1_imp(T a, T b, T z, const Policy& pol)
  {
    // some special cases
    // ...

    return detail::hypergeometric_1f1_generic_series(a, b, z, pol);
  }

  template <class T, class Policy>
  inline T hypergeometric_1f2_imp(T a, T b1, T b2, T z, const Policy& pol)
  {
    // some special cases
    // ...

    return detail::hypergeometric_1f2_generic_series(a, b1, b2, z, pol);
  }

  template <class T, class Policy>
  inline T hypergeometric_2f1_imp(T a1, T a2, T b, T z, const Policy& pol)
  {
    // some special cases
    // ...

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
