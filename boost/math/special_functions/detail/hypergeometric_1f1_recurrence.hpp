
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

  #include <deque>

  #include <boost/math/special_functions/modf.hpp>
  #include <boost/math/special_functions/next.hpp>

  namespace boost { namespace math { namespace detail {

  template <class T>
  struct hypergeometric_1f1_backward_recurrence_a_next_coefficients
  {
    typedef std::pair<T, T> result_type;

    hypergeometric_1f1_backward_recurrence_a_next_coefficients(const T& a, const T& b, const T& z):
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
  struct hypergeometric_1f1_forward_recurrence_a_next_coefficients
  {
    typedef std::pair<T, T> result_type;

    hypergeometric_1f1_forward_recurrence_a_next_coefficients(const T& a, const T& b, const T& z):
      b(b), z(z), a(a)
    {
    }

    result_type operator()()
    {
      const T cn = a - b; 
      const result_type result = std::make_pair<T, T>(a / cn, (((2 * a) - b) + z) / cn);
      ++a;
      return result;
    }

  private:
    const T b, z;
    T a;
  };

  template <class T>
  struct hypergeometric_1f1_backward_recurrence_b_next_coefficients
  {
    typedef std::pair<T, T> result_type;

    hypergeometric_1f1_backward_recurrence_b_next_coefficients(const T& a, const T& b, const T& z):
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

  template <class T>
  struct hypergeometric_1f1_backward_recurrence_a_and_b_next_coefficients
  {
    typedef std::pair<T, T> result_type;

    hypergeometric_1f1_backward_recurrence_a_and_b_next_coefficients(const T& a, const T& b, const T& z):
      z(z), a(a), b(b)
    {
    }

    result_type operator()()
    {
      const T cn = a * z;
      const result_type result = std::make_pair<T, T>((b * (1 - b)) / cn, (b * ((1 - b) + z)) / cn);
      --a; --b;
      return result;
    }

  private:
    const T z;
    T a, b;
  };

  template <class T, class NextCoefs>
  inline T hypergeometric_1f1_recurrence(NextCoefs& get_coefs, boost::uintmax_t last_index, T& first, T& second, std::deque<T>* hyp_results = 0)
  {
    using std::swap;

    T third = 0;

    if (hyp_results)
    {
      hyp_results->clear();
      hyp_results->resize(last_index);
    }

    for (boost::uintmax_t k = 0; k < last_index; ++k)
    {
      std::pair<T, T> next = get_coefs();

      third = ((next.second * second) - first) / next.first;
      if (hyp_results) hyp_results->push_back(third);

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
    BOOST_MATH_STD_USING // modf, frexp, fabs, pow

    boost::intmax_t integer_part = 0;
    const T bk = modf(b, &integer_part);
    T ak = modf(a, &integer_part);

    int exp_of_a = 0; frexp(a, &exp_of_a);
    int exp_of_b = 0; frexp(b, &exp_of_b);

    const bool are_fractional_parts_close_enough =
      fabs(boost::math::float_distance(ak, bk)) <= pow(2, ((std::max)(exp_of_a, exp_of_b)));

    if ((a < b) && (b < 0) && (are_fractional_parts_close_enough)) // TODO: to be researched in detail
    {
      ak = b - 1;
      integer_part -= ceil(b) - 1;
    }
    
    T first = detail::hypergeometric_1f1_imp(ak, b, z, pol);
    --ak;
    T second = detail::hypergeometric_1f1_imp(ak, b, z, pol);

    const boost::uintmax_t last_index = fabs(integer_part);
    detail::hypergeometric_1f1_backward_recurrence_a_next_coefficients<T> s(ak, b, z);

    return detail::hypergeometric_1f1_recurrence(s, last_index, first, second);
  }

  template <class T, class Policy>
  inline T hypergeometric_1f1_forward_recurrence_for_positive_a(const T& a, const T& b, const T& z, const Policy& pol)
  {
    BOOST_MATH_STD_USING // modf, fabs

    boost::intmax_t integer_part = 0;
    T ak = modf(a, &integer_part);

    T first = detail::hypergeometric_1f1_imp(ak, b, z, pol);
    ++ak;
    T second = detail::hypergeometric_1f1_imp(ak, b, z, pol);

    const boost::uintmax_t last_index = fabs(integer_part);
    detail::hypergeometric_1f1_forward_recurrence_a_next_coefficients<T> s(ak, b, z);

    return detail::hypergeometric_1f1_recurrence(s, last_index, first, second);
  }

  template <class T, class Policy>
  inline T hypergeometric_1f1_backward_recurrence_for_negative_b(const T& a, const T& b, const T& z, const Policy& pol)
  {
    BOOST_MATH_STD_USING // modf, fabs

    boost::intmax_t integer_part = 0;
    T bk = modf(b, &integer_part);

    T first = detail::hypergeometric_1f1_imp(a, bk, z, pol);
    --bk;
    T second = detail::hypergeometric_1f1_imp(a, bk, z, pol);

    const boost::uintmax_t last_index = fabs(integer_part);
    detail::hypergeometric_1f1_backward_recurrence_b_next_coefficients<T> s(a, bk, z);

    return detail::hypergeometric_1f1_recurrence(s, last_index, first, second);
  }

  // this method works provided that integer part of a is the same as integer part of b
  // we don't make this check inside it.
  template <class T, class Policy>
  inline T hypergeometric_1f1_backward_recurrence_for_negative_a_and_b(const T& a, const T& b, const T& z, const Policy& pol)
  {
    BOOST_MATH_STD_USING // modf, fabs

    boost::intmax_t integer_part = 0;
    T ak = modf(a, &integer_part);
    T bk = modf(b, &integer_part);

    T first = detail::hypergeometric_1f1_imp(ak, bk, z, pol);
    --ak; --bk;
    T second = detail::hypergeometric_1f1_imp(ak, bk, z, pol);

    const boost::uintmax_t last_index = fabs(integer_part);
    detail::hypergeometric_1f1_backward_recurrence_a_and_b_next_coefficients<T> s(ak, bk, z);

    return detail::hypergeometric_1f1_recurrence(s, last_index, first, second);
  }

  // ranges
  template <class T>
  inline bool hypergeometric_1f1_is_a_small_enough(const T& a)
  {
    return a < -50; // TODO: make dependent on precision
  }

  } } } // namespaces

#endif // BOOST_HYPERGEOMETRIC_1F1_RECURRENCE_HPP_
