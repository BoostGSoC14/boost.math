//  (C) Copyright Anton Bikineev 2014
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_RECURRENCE_HPP_
  #define BOOST_MATH_TOOLS_RECURRENCE_HPP_

  #include <vector>

  #include <boost/utility/enable_if.hpp>
  #include <boost/mpl/if.hpp>
  #include <boost/math/tools/config.hpp>
  #include <boost/math/tools/precision.hpp>
  #include <boost/math/tools/tuple.hpp>

  namespace boost { namespace math { namespace tools {

  namespace detail {

  template <class T>
  struct is_homogeneous:
      boost::mpl::if_c<
        boost::math::tuple_size<typename T::result_type>::value == 3u,
        boost::mpl::true_,
        boost::mpl::false_>::type {};

  // solves recurrence relation defined by homogeneous difference equation
  template <class Functor, class U, class T>
  inline T solve_recurrence_relation_by_oliver_imp(Functor& get_next_coefs, const U& factor, unsigned index, const T& init_value, boost::mpl::true_)
  {
    BOOST_MATH_STD_USING // fabs
    typedef typename Functor::result_type coef_tuple;

    std::vector<T> p, e;
    p.reserve(index); e.reserve(index);

    // initialization
    get_next_coefs(); //idle (because p[1] == 1)
    coef_tuple coefs = get_next_coefs();
    T an = boost::math::get<0>(coefs),
      bn = boost::math::get<1>(coefs),
      cn = boost::math::get<2>(coefs);

    p.push_back(0); p.push_back(1);
    e.push_back(init_value); e.push_back((cn * init_value) / an);

    T check_n = 0;
    T min_check_n = boost::math::tools::max_value<T>();
    unsigned i = 2;

    // forward recurrence
    do {
      const T next_p = ((bn * p[i-1]) - (cn * p[i-2])) / an;

      coefs = get_next_coefs();
      an = boost::math::get<0>(coefs),
      bn = boost::math::get<1>(coefs),
      cn = boost::math::get<2>(coefs);

      const T next_e = (cn * e[i-1]) / an;

      p.push_back(next_p); e.push_back(next_e);

      check_n = fabs(e[i-1] / (p[i-1] * p[i]));

      if ((i <= index) && (check_n < min_check_n))
        min_check_n = check_n;

      ++i;
    } while (check_n > fabs(factor * min_check_n));

    std::vector<T> w;
    w.resize(p.size());

    unsigned k = w.size();
    w[--k] = 0;

    // backward recurrence
    for (; k >= index; --k)
      w[k-1] = (p[k-1] * w[k] + e[k-1]) / p[k];

    return w[index];
  }

  // TODO: after some optimization of previous function
  // there is a need to add function for solving inhomogeneous
  // difference equation here. It's just a dummy function now.
  template <class Functor, class U, class T>
  inline T solve_recurrence_relation_by_oliver_imp(Functor& get_next_coefs, const U& factor, unsigned index, const T& init_value, boost::mpl::false_)
  {
    return 0;
  }

  } // namespace detail

  // solves difference equations of the following form:
  // a(n)w(n+1) - b(n)w(n) + c(n)w(n-1) = d(n) - inhomogeneous case
  // a(n)w(n+1) - b(n)w(n) + c(n)w(n-1) = 0    - homogeneous case
  //
  // This implementation uses Olver's algorithm because of some
  // valuable advantages as opposed to usual Miller's algorithm
  template <class Coefficients, class U, class T>
  inline T solve_recurrence_relation_by_oliver(Coefficients& coefs, const U& factor, unsigned index, const T& init_value)
  {
    typedef typename detail::is_homogeneous<Coefficients>::type is_homogeneous;

    return detail::solve_recurrence_relation_by_oliver_imp(coefs, factor, index, init_value, is_homogeneous());
  }

  } } } // namespaces

#endif // BOOST_MATH_TOOLS_RECURRENCE_HPP_
