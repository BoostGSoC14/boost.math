#include "instantiated.hpp"

#include <boost/math/special_functions.hpp>
#include <boost/math/special_functions/hypergeometric.hpp>
#include <boost/math/special_functions/detail/hypergeometric_rational.hpp>
#include <boost/math/special_functions/detail/hypergeometric_1f1_bessel.hpp>

#include <boost/math/policies/policy.hpp>

namespace boost{ namespace math {

namespace
{
using namespace policies;

const policy<
  overflow_error<ignore_error>,
  underflow_error<ignore_error>,
  denorm_error<ignore_error>
      > pol;

}

// only function boilerplates are here

// interface functions
float hypergeometric_1f1_interface_f(float a, float b, float z)
{
  return boost::math::detail::hypergeometric_1f1_imp(a, b, z, pol);
}
double hypergeometric_1f1_interface_d(double a, double b, double z)
{
  return boost::math::detail::hypergeometric_1f1_imp(a, b, z, pol);
}
long double hypergeometric_1f1_interface_l(long double a, long double b, long double z)
{
  return boost::math::detail::hypergeometric_1f1_imp(a, b, z, pol);
}

// generic Taylor series functions
float hypergeometric_1f1_series_f(float a, float b, float z)
{
  return boost::math::detail::hypergeometric_1f1_generic_series(a, b, z, pol);
}
double hypergeometric_1f1_series_d(double a, double b, double z)
{
  return boost::math::detail::hypergeometric_1f1_generic_series(a, b, z, pol);
}
long double hypergeometric_1f1_series_l(long double a, long double b, long double z)
{
  return boost::math::detail::hypergeometric_1f1_generic_series(a, b, z, pol);
}

// asymptotic approx series functions
float hypergeometric_1f1_asym_f(float a, float b, float z)
{
  return boost::math::detail::hypergeometric_1f1_asym_series(a, b, z, pol);
}
double hypergeometric_1f1_asym_d(double a, double b, double z)
{
  return boost::math::detail::hypergeometric_1f1_asym_series(a, b, z, pol);
}
long double hypergeometric_1f1_asym_l(long double a, long double b, long double z)
{
  return boost::math::detail::hypergeometric_1f1_asym_series(a, b, z, pol);
}

// Luke's rational functions
float hypergeometric_1f1_rational_f(float a, float b, float z)
{
  return boost::math::detail::hypergeometric_1f1_rational(a, b, z, pol);
}
double hypergeometric_1f1_rational_d(double a, double b, double z)
{
  return boost::math::detail::hypergeometric_1f1_rational(a, b, z, pol);
}
long double hypergeometric_1f1_rational_l(long double a, long double b, long double z)
{
  return boost::math::detail::hypergeometric_1f1_rational(a, b, z, pol);
}

// A&S 13_3_7 formula through Bessel functions
float hypergeometric_1f1_13_3_7_f(float a, float b, float z)
{
  return boost::math::detail::hypergeometric_1f1_13_3_7_series(a, b, z, pol);
}
double hypergeometric_1f1_13_3_7_d(double a, double b, double z)
{
  return boost::math::detail::hypergeometric_1f1_13_3_7_series(a, b, z, pol);
}
long double hypergeometric_1f1_13_3_7_l(long double a, long double b, long double z)
{
  return boost::math::detail::hypergeometric_1f1_13_3_7_series(a, b, z, pol);
}

}} // namespaces
