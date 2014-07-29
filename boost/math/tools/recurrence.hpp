//  (C) Copyright Anton Bikineev 2014
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#ifndef BOOST_MATH_TOOLS_RECURRENCE_HPP_
  #define BOOST_MATH_TOOLS_RECURRENCE_HPP_

  #include <vector>
  #include <algorithm>
  #include <functional>

  #include <boost/utility/enable_if.hpp>
  #include <boost/mpl/if.hpp>
  #include <boost/bind.hpp>

  #include <boost/math/tools/config.hpp>
  #include <boost/math/tools/precision.hpp>
  #include <boost/math/tools/tuple.hpp>
  #include <boost/math/special_functions/next.hpp>

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
  inline std::pair<T, unsigned> olver_checked_recurrence_imp(Functor& get_coefs, const U& factor, const T& init_value, unsigned init_pos, unsigned index, boost::mpl::true_)
  {
    BOOST_MATH_STD_USING // fabs
    typedef typename Functor::result_type coef_tuple;

    std::vector<T> p, e, check_ns;
    p.reserve(index); e.reserve(index); check_ns.reserve(index);

    // initialization
    coef_tuple coefs = get_coefs(init_pos + 1);
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

      coefs = get_coefs(init_pos + i);
      an = boost::math::get<0>(coefs),
      bn = boost::math::get<1>(coefs),
      cn = boost::math::get<2>(coefs);

      const T next_e = (cn * e[i-1]) / an;

      p.push_back(next_p); e.push_back(next_e);

      // check for close overflow
      if (log((bn * p[i]) - (cn * p[i-1])) >= log(boost::math::tools::max_value<T>()) + log(an)) // TODO: this check takes some time
      {
        typename std::vector<T>::iterator min_check_it =
            std::find_if(check_ns.begin(), check_ns.end(), boost::bind(std::less<T>(), _1, check_ns.back() / factor)); // TODO: find better algorithm with no boost::bind
        index = std::distance(check_ns.begin(), min_check_it);

        break;
      }

      check_n = fabs(e[i-1] / (p[i-1] * p[i]));
      if ((i <= index) && (check_n < min_check_n))
        min_check_n = check_n;

      check_ns.push_back(check_n);

      ++i;
    } while (check_n > fabs(factor * min_check_n));

    std::vector<T> w;
    w.resize(p.size());

    unsigned k = w.size();
    w[--k] = 0;

    // backward recurrence
    for (; k >= index; --k)
      w[k-1] = (p[k-1] * w[k] + e[k-1]) / p[k];

    return std::make_pair(w[index], index);
  }

  // TODO: after some optimization of previous function
  // there is a need to add function for solving inhomogeneous
  // difference equation here. It's just a dummy function now.
  template <class Functor, class U, class T>
  inline T olver_checked_recurrence_imp(Functor& get_coefs, const U& factor, unsigned index, const T& init_value, boost::mpl::false_)
  {
    return 0;
  }

  // this wrapper-implementation protects us from possible overflow
  template <class Functor, class U, class T, class IsHomogeneous>
  inline T solve_recurrence_relation_by_olver_imp(Functor& get_coefs, const U& factor, unsigned index, T init_value, const IsHomogeneous& is_homogeneous)
  {
    unsigned init_pos = 0u;
    unsigned new_index = index;

    do
    {
      std::pair<T, unsigned> result = detail::olver_checked_recurrence_imp(get_coefs, factor, init_value, init_pos, new_index, is_homogeneous);

      init_pos += result.second;
      new_index -= result.second;
      init_value = result.first;

    } while (init_pos < index);

    return init_value;
  }

  } // namespace detail

  // solves difference equations of the following form:
  // a(n)w(n+1) - b(n)w(n) + c(n)w(n-1) = d(n) - inhomogeneous case
  // a(n)w(n+1) - b(n)w(n) + c(n)w(n-1) = 0    - homogeneous case
  //
  // This implementation uses Olver's algorithm because of some
  // valuable advantages as opposed to usual Miller's algorithm
  template <class Coefficients, class U, class T>
  inline T solve_recurrence_relation_by_olver(Coefficients& coefs, const U& factor, unsigned index, const T& init_value)
  {
    typedef typename detail::is_homogeneous<Coefficients>::type is_homogeneous;

    return detail::solve_recurrence_relation_by_olver_imp(coefs, factor, index, init_value, is_homogeneous());
  }

  } } } // namespaces

#endif // BOOST_MATH_TOOLS_RECURRENCE_HPP_
