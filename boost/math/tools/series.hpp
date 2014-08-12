//  (C) Copyright John Maddock 2005-2006.
//  (C) Copyright Anton Bikineev 2014.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_SERIES_INCLUDED
#define BOOST_MATH_TOOLS_SERIES_INCLUDED

#ifdef _MSC_VER
#pragma once
#endif

#include <utility>

#include <boost/config/no_tr1/cmath.hpp>
#include <boost/cstdint.hpp>
#include <boost/limits.hpp>
#include <boost/math/tools/config.hpp>
#include <boost/math/tools/precision.hpp>

namespace boost{ namespace math{ namespace tools{

//
// Simple series summation come first:
//
template <class Functor, class U, class V>
inline typename Functor::result_type sum_series(Functor& func, const U& factor, boost::uintmax_t& max_terms, const V& init_value)
{
   BOOST_MATH_STD_USING

   typedef typename Functor::result_type result_type;

   boost::uintmax_t counter = max_terms;

   result_type result = init_value;
   result_type next_term;
   do{
      next_term = func();
      result += next_term;
   }
   while((fabs(factor) < fabs(next_term / result)) && --counter);

   // set max_terms to the actual number of terms of the series evaluated:
   max_terms = max_terms - counter;

   return result;
}

template <class Functor, class U>
inline typename Functor::result_type sum_series(Functor& func, const U& factor, boost::uintmax_t& max_terms)
{
   typename Functor::result_type init_value = 0;
   return sum_series(func, factor, max_terms, init_value);
}

//
// This experimental function works as previous one,
// but also returns an approximate amount of exact digits,
// quantifying possible cancellations on the fly.
//
// According to Henrici P. "Essentials of numerical analysis
// with pocket calculator demonstrations",
// d(z) ~ log10(mu(z)) - log10(f(z)),
// where d(z) is number of lost digits,
//       mu(z) is max term,
//       f(z) is function to be found.
//
// To accelerate evaluation for multiprecision numbers, we may assume that
// d(z) ~ digits10(mu(z)) + exponent(mu(z)) - digits10(f(z)) - exponent(f(z)).
// Since we use a fixed precision:
// d(z) ~ exponent(mu(z)) - exponent(f(z)).
//
template <class Functor, class U, class V>
inline std::pair<typename Functor::result_type, unsigned> sum_cancelled_series(Functor& func, const U& factor, boost::uintmax_t& max_terms, const V& init_value)
{
   BOOST_MATH_STD_USING

   typedef typename Functor::result_type result_type;

   boost::uintmax_t counter = max_terms;

   result_type result = init_value;
   result_type next_term;
   result_type max_term = 0;
   do{
      next_term = func();
      const result_type abs_next_term = fabs(next_term);
      if (abs_next_term > max_term)
         max_term = abs_next_term;
      result += next_term;
   }
   while((fabs(factor) < fabs(next_term / result)) && --counter);

   // set max_terms to the actual number of terms of the series evaluated:
   max_terms = max_terms - counter;

   int max_term_exp; frexp(max_term, &max_term_exp);
   int result_exp; frexp(result, &result_exp);

   int approximate_lost_digits = max_term_exp - result_exp;
   unsigned exact_digits = (std::max)(boost::math::tools::digits<result_type>() - approximate_lost_digits, 0);

   return std::make_pair(result, (exact_digits * 30103ul) / 100000ul);
}

template <class Functor, class U>
inline std::pair<typename Functor::result_type, unsigned> sum_cancelled_series(Functor& func, const U& factor, boost::uintmax_t& max_terms)
{
   typename Functor::result_type init_value = 0;
   return sum_cancelled_series(func, factor, max_terms, init_value);
}

template <class Functor, class U>
inline typename Functor::result_type sum_series(Functor& func, int bits, boost::uintmax_t& max_terms, const U& init_value)
{
   BOOST_MATH_STD_USING
   typedef typename Functor::result_type result_type;
   result_type factor = ldexp(result_type(1), 1 - bits);
   return sum_series(func, factor, max_terms, init_value);
}

template <class Functor>
inline typename Functor::result_type sum_series(Functor& func, int bits)
{
   BOOST_MATH_STD_USING
   typedef typename Functor::result_type result_type;
   boost::uintmax_t iters = (std::numeric_limits<boost::uintmax_t>::max)();
   result_type init_val = 0;
   return sum_series(func, bits, iters, init_val);
}

template <class Functor>
inline typename Functor::result_type sum_series(Functor& func, int bits, boost::uintmax_t& max_terms)
{
   BOOST_MATH_STD_USING
   typedef typename Functor::result_type result_type;
   result_type init_val = 0;
   return sum_series(func, bits, max_terms, init_val);
}

template <class Functor, class U>
inline typename Functor::result_type sum_series(Functor& func, int bits, const U& init_value)
{
   BOOST_MATH_STD_USING
   boost::uintmax_t iters = (std::numeric_limits<boost::uintmax_t>::max)();
   return sum_series(func, bits, iters, init_value);
}

//
// Algorithm kahan_sum_series invokes Functor func until the N'th
// term is too small to have any effect on the total, the terms
// are added using the Kahan summation method.
//
// CAUTION: Optimizing compilers combined with extended-precision
// machine registers conspire to render this algorithm partly broken:
// double rounding of intermediate terms (first to a long double machine
// register, and then to a double result) cause the rounding error computed
// by the algorithm to be off by up to 1ulp.  However this occurs rarely, and
// in any case the result is still much better than a naive summation.
//
template <class Functor>
inline typename Functor::result_type kahan_sum_series(Functor& func, int bits)
{
   BOOST_MATH_STD_USING

   typedef typename Functor::result_type result_type;

   result_type factor = pow(result_type(2), bits);
   result_type result = func();
   result_type next_term, y, t;
   result_type carry = 0;
   do{
      next_term = func();
      y = next_term - carry;
      t = result + y;
      carry = t - result;
      carry -= y;
      result = t;
   }
   while(fabs(result) < fabs(factor * next_term));
   return result;
}

template <class Functor>
inline typename Functor::result_type kahan_sum_series(Functor& func, int bits, boost::uintmax_t& max_terms)
{
   BOOST_MATH_STD_USING

   typedef typename Functor::result_type result_type;

   boost::uintmax_t counter = max_terms;

   result_type factor = ldexp(result_type(1), bits);
   result_type result = func();
   result_type next_term, y, t;
   result_type carry = 0;
   do{
      next_term = func();
      y = next_term - carry;
      t = result + y;
      carry = t - result;
      carry -= y;
      result = t;
   }
   while((fabs(result) < fabs(factor * next_term)) && --counter);

   // set max_terms to the actual number of terms of the series evaluated:
   max_terms = max_terms - counter;

   return result;
}

} // namespace tools
} // namespace math
} // namespace boost

#endif // BOOST_MATH_TOOLS_SERIES_INCLUDED

