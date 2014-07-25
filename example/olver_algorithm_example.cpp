#include <iostream>

#include <boost/math/tools/config.hpp>
#include <boost/math/special_functions.hpp>
#include <boost/math/tools/recurrence.hpp>
#include <boost/math/special_functions/hypergeometric.hpp>

template <class T>
struct hypergeometric_1f1_recurrence_b_next_coefficients
{
  typedef boost::math::tuple<T, T, T> result_type;

  hypergeometric_1f1_recurrence_b_next_coefficients(const T& a, const T& b, const T& z):
  a(a), z(z), b(b)
  {
  }

  result_type operator()()
  {
    const T an = z * (b - a);
    const T bn = b * ((z + b) - 1);
    const T cn = b * (b - 1);
    const result_type result = boost::math::tuple<T, T, T>(an, bn, cn);
    ++b;
    return result;
  }

private:
  const T a, z;
  T b;
};

template <class T>
struct hypergeometric_1f1_recurrence_a_and_b_next_coefficients
{
  typedef boost::math::tuple<T, T, T> result_type;

  hypergeometric_1f1_recurrence_a_and_b_next_coefficients(const T& a, const T& b, const T& z):
  z(z), a(a), b(b)
  {
  }

  result_type operator()()
  {
    const T an = a * z;
    const T bn = b * ((1 - b) + z);
    const T cn = b * (1 - b);
    const result_type result = boost::math::tuple<T, T, T>(an, bn, cn);
    ++a; ++b;
    return result;
  }

private:
  const T z;
  T a, b;
};

int main()
{
  std::cout.precision(16);
  using float_type = long double;

  float_type a = 0.1, b = 0.2, z = 5;
  unsigned m = 10u;

  hypergeometric_1f1_recurrence_a_and_b_next_coefficients<float_type> s(a, b, z);
  float_type initial = boost::math::hypergeometric_1f1(a, b, z);

  std::cout << boost::math::tools::solve_recurrence_relation_by_oliver(
      s,
      boost::math::tools::epsilon<float_type>(),
      m,
      initial
      )
    << std::endl;
}
