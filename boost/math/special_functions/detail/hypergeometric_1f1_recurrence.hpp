
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

  template <class T>
  struct hypergeometric_1f1_recurrence_a_next_coefficients
  {
    typedef std::pair<T, T> result_type;

    hypergeometric_1f1_recurrence_a_next_coefficients(const T& a, const T& b, const T& z):
      b(b), z(z), a(a)
    {
    }

    result_type operator()()
    {
      const result_type result = std::make_pair<T, T>((a - b) / a, (((2 * a) - b) + z) / a);
      --a;
      return result;
    }

  private:
    const T b, z;
    T a;
  };

  template <class T>
  struct hypergeometric_1f1_recurrence_b_next_coefficients
  {
    typedef std::pair<T, T> result_type;

    hypergeometric_1f1_recurrence_b_next_coefficients(const T& a, const T& b, const T& z):
      a(a), z(z), b(b)
    {
    }

    result_type operator()()
    {
      const T cn = z * (b - a);
      const result_type result = std::make_pair<T, T>((b * (b - 1)) / cn, (b * ((b - 1) + z)) / cn);
      --b;
      return result;
    }

  private:
    const T a, z;
    T b;
  };

  template <class T, class NextCoefs>
  inline T hypergeometric_1f1_recurrence(NextCoefs& get_coefs, const T& last_index, T& first, T& second)
  {
    using std::swap;

    T third = 0;

    // probably we will need to change type of last and k from T to boost::uintmax_t
    for (T k = 0; k < last_index; ++k)
    {
      std::pair<T, T> next = get_coefs();

      third = ((next.second * second) - first) / next.first;

      swap(first, second);
      swap(second, third);
    }

    return first;
  }

  // forward declaration
  template <class T, class Policy>
  inline T hypergeometric_1f1_imp(const T& a, const T& b, const T& z, const Policy& pol);

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
    const T last_index = fabs(integer_part);
    detail::hypergeometric_1f1_recurrence_a_next_coefficients<T> s(ak, b, z);

    return detail::hypergeometric_1f1_recurrence(s, last_index, first, second);
  }

  template <class T, class Policy>
  inline T hypergeometric_1f1_backward_recurrence_for_negative_b(const T& a, const T& b, const T& z, const Policy& pol)
  {
    BOOST_MATH_STD_USING // modf, fabs

    T integer_part = 0;
    T bk = modf(b, &integer_part);

    // we will never get to infinite recursion here:
    T first = detail::hypergeometric_1f1_imp(a, bk, z, pol);
    --bk;
    T second = detail::hypergeometric_1f1_imp(a, bk, z, pol);

    // probably we will need to change type of last and k from T to boost::uintmax_t
    const T last_index = fabs(integer_part);
    detail::hypergeometric_1f1_recurrence_b_next_coefficients<T> s(a, bk, z);

    return detail::hypergeometric_1f1_recurrence(s, last_index, first, second);
  }

  } } } // namespaces

#endif // BOOST_HYPERGEOMETRIC_1F1_RECURRENCE_HPP_
