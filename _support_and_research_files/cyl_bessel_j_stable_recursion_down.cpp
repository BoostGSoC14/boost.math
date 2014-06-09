
#include <algorithm>
#include <deque>
#include <functional>
#include <iomanip>
#include <iostream>

#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/gamma.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>

namespace boost { namespace math { namespace hypergeometric_detail {
template<typename T>
T cyl_bessel_j_stable_recursion_down(const T& v,
                                     const T& x,
                                     std::deque<T>* const cyl_bessel_j_results = static_cast<std::deque<T>* const>(0U))
{
  using std::fabs;
  using std::modf;
  using std::pow;

  // Use recursion and normalization with a Neumann sum.

  T n;
  const T v_frac = modf(v, &n);

  // TBD: Calculate the upper bound of the recursion dynamically.
  const int upper_bound_of_recursion = 14;

  int k                   = upper_bound_of_recursion / 2;
  T   one_over_k_fact     = T(1) / boost::math::unchecked_factorial<T>(static_cast<unsigned int>(k));
  T   k_plus_v_frac       = v_frac + k;
  T   gamma_k_plus_v_frac = boost::math::tgamma(k_plus_v_frac);

  T Jv_p2(0);
  T Jv_p1(1);
  T Jv;
  T Jvn;

  const T v_frac_plus_two_k_over_k_fact = (v_frac + (k * 2)) * one_over_k_fact;

  T sum_v_frac = v_frac_plus_two_k_over_k_fact * gamma_k_plus_v_frac;

  const bool want_cyl_bessel_j_results = (cyl_bessel_j_results != static_cast<std::deque<T>* const>(0U));

  if(want_cyl_bessel_j_results)
  {
    cyl_bessel_j_results->resize(upper_bound_of_recursion);
  }

  // Do the downward recursion of Jv:
  //
  //                  Jv+1
  //   Jv = [ 2 (v+1) ---- ] - Jv+2
  //                   x 

  const T two_over_x = T(2) / x;

  for(int m = (upper_bound_of_recursion - 1); m >= 0; --m)
  {
    const T n_plus_v_frac_plus_one = v_frac + (m + 1);

    // Downward recursion for Jv_frac.
    Jv    = ((n_plus_v_frac_plus_one * Jv_p1) * two_over_x) - Jv_p2;
    Jv_p2 = Jv_p1;
    Jv_p1 = Jv;

    if(want_cyl_bessel_j_results)
    {
      cyl_bessel_j_results->at(m) = Jv;
    }

    // Store the value of the Bessel function which has the sought order.
    if(n == m)
    {
      Jvn = Jv;
    }

    // Do the normalization using a Neumann sum which is:
    // (x/2)^v = Sum_k { [((v + 2k) gamma(v + k)) / k!] * J_v+2k }
    if((m % 2) == 0)
    {
      one_over_k_fact *= k;

      --k;

      --k_plus_v_frac;

      gamma_k_plus_v_frac /= k_plus_v_frac;

      // Increment the normalization sum for Jv_frac.
      sum_v_frac += (((v_frac + m) * gamma_k_plus_v_frac) * Jv) * one_over_k_fact;
    }
  }

  const T x_half_pow_v_frac = pow(x / 2, T(v_frac));

  // Compute the normalizations using the Neumann sums in combination with (x/2)^v.
  const T norm =  x_half_pow_v_frac / sum_v_frac;

  if(want_cyl_bessel_j_results)
  {
    std::transform(cyl_bessel_j_results->begin(),
                   cyl_bessel_j_results->end(),
                   cyl_bessel_j_results->begin(),
                   std::bind2nd(std::multiplies<T>(), norm));
  }

  // Normalize In_v_frac.
  return Jvn * norm;
}
} } } // namespace boost::math::hypergeometric_detail

int main()
{
//  typedef boost::multiprecision::cpp_dec_float_100 float_type;
  typedef double float_type;

  // Mathematica(R) or Wolfram's Alpha:
  // Table[N[BesselJ[n + (1/3), EulerGamma], 100], {n, 0, 56, 1}]

  const float_type v = float_type(5) + (float_type(1) / 3);
  const float_type x = boost::math::constants::euler<float_type>();

  const float_type jv = boost::math::hypergeometric_detail::cyl_bessel_j_stable_recursion_down(v, x);

  std::cout << std::setprecision(std::numeric_limits<float_type>::digits10)
            << std::scientific
            << jv
            << std::endl;
}
