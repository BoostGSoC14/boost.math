#include <iostream>
#include <algorithm>
#include <iterator>

#include <boost/math/tools/config.hpp>
#include <boost/math/tools/recurrence.hpp>
#include <boost/math/special_functions/hypergeometric.hpp>

#include <boost/multiprecision/cpp_dec_float.hpp>

template <class T>
struct hypergeometric_1f1_recurrence_b_coefficients
{
  typedef boost::math::tuple<T, T, T> result_type;

  hypergeometric_1f1_recurrence_b_coefficients(const T& a, const T& b, const T& z):
  a(a), b(b), z(z)
  {
  }

  result_type operator()(int i) const
  {
    const T bi = b + i;

    const T an = z * (bi - a);
    const T bn = bi * ((z + bi) - 1);
    const T cn = bi * (bi - 1);

    return boost::math::tuple<T, T, T>(an, bn, cn);
  }

private:
  const T a, b, z;
};

template <class T>
struct hypergeometric_1f1_recurrence_a_and_b_coefficients
{
  typedef boost::math::tuple<T, T, T> result_type;

  hypergeometric_1f1_recurrence_a_and_b_coefficients(const T& a, const T& b, const T& z):
  a(a), b(b), z(z)
  {
  }

  result_type operator()(int i) const
  {
    const T ai = a + i;
    const T bi = b + i;

    const T an = ai * z;
    const T bn = bi * ((1 - bi) + z);
    const T cn = bi * (1 - bi);

    return boost::math::tuple<T, T, T>(an, bn, cn);
  }

private:
  const T a, b, z;
};

int main()
{
  using float_type = double;

  const float_type a = 0.1, b = 0.2, z = 5;
  unsigned m = 1000000u;

  hypergeometric_1f1_recurrence_b_coefficients<float_type> s(a, b, z);
  const float_type initial = boost::math::hypergeometric_1f1(a, b, z);

  std::cout.precision(std::numeric_limits<float_type>::digits10);

  // outputs 1f1(0.1, 1000000.2, 5)
  std::cout << boost::math::tools::solve_recurrence_relation_by_olver(
      s,
      boost::math::tools::epsilon<float_type>(),
      m,
      initial
      )
    << std::endl;
}
