
///////////////////////////////////////////////////////////////////////////////
//  Copyright 2014 Anton Bikineev
//  Copyright 2014 Christopher Kormanyos
//  Copyright 2014 John Maddock
//  Copyright 2014 Paul Bristow
//  Distributed under the Boost
//  Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_HYPERGEOMETRIC_1F1_RECURRENCE_HPP_
  #define BOOST_HYPERGEOMETRIC_1F1_RECURRENCE_HPP_

  #include <boost/math/special_functions/modf.hpp>

  namespace boost { namespace math { namespace detail {

  // forward declaration
  template <class T, class Policy>
  inline T hypergeometric_1f1_imp(const T& a, const T& b, const T& z, const Policy& pol);

  template<class T>
  inline T hypergeometric_1f1_recurrence_a(T&a, const T& b, const T& z, const T& last, T& first, T& second)
  {
    T third = 0;

    // probably we will need to change type of last and k from T to boost::uintmax_t
    for (T k = 0; k < last; ++k)
    {
      const T an = (a - b) / a;
      const T bn = (((2 * a) - b) + z) / a;

      third = ((bn * second) - first) / an;

      --a;
      swap(first, second);
      swap(second, third);
    }

    return first;
  }

  template <class T, class Policy>
  inline T hypergeometric_1f1_backward_recurrence_for_negative_a(const T& a, const T& b, const T& z, const Policy& pol)
  {
    BOOST_MATH_STD_USING // modf, fabs

    T integer_part = 0;
    T ak = modf(a, &integer_part);

    // we will never get to infinite recursion here:
    T first = detail::hypergeometric_1f1_imp(ak, b, z, pol);
    --ak;
    T second = detail::hypergeometric_1f1_imp(ak, b, z, pol);

    // probably we will need to change type of last and k from T to boost::uintmax_t
    const T last = fabs(integer_part);
    return detail::hypergeometric_1f1_recurrence_a(ak, b, z, last, first, second);
  }

  } } } // namespaces

#endif // BOOST_HYPERGEOMETRIC_1F1_RECURRENCE_HPP_
