#ifndef BOOST_MATH_HYPERGEOMETRIC_INSTANTIATED_HPP
  #define BOOST_MATH_HYPERGEOMETRIC_INSTANTIATED_HPP

#ifdef __cplusplus
  extern "C" {

  namespace boost{ namespace math{
#endif

  float hypergeometric_1f1_interface_f(float a, float b, float z);
  double hypergeometric_1f1_interface_d(double a, double b, double z);
  long double hypergeometric_1f1_interface_l(long double a, long double b, long double z);

  float hypergeometric_1f1_series_f(float a, float b, float z);
  double hypergeometric_1f1_series_d(double a, double b, double z);
  long double hypergeometric_1f1_series_l(long double a, long double b, long double z);

  float hypergeometric_1f1_asym_f(float a, float b, float z);
  double hypergeometric_1f1_asym_d(double a, double b, double z);
  long double hypergeometric_1f1_asym_l(long double a, long double b, long double z);

  float hypergeometric_1f1_rational_f(float a, float b, float z);
  double hypergeometric_1f1_rational_d(double a, double b, double z);
  long double hypergeometric_1f1_rational_l(long double a, long double b, long double z);

  float hypergeometric_1f1_13_3_7_f(float a, float b, float z);
  double hypergeometric_1f1_13_3_7_d(double a, double b, double z);
  long double hypergeometric_1f1_13_3_7_l(long double a, long double b, long double z);

#ifdef __cplusplus
  }} // namespaces

  } // extern "C"
#endif

#endif // BOOST_MATH_HYPERGEOMETRIC_INSTANTIATED_HPP
