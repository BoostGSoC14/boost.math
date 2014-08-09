//  (C) Copyright John Maddock 2007.
//  (C) Copyright Anton Bikineev 2014.
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

#define BOOST_MATH_OVERFLOW_ERROR_POLICY ignore_error
#include <boost/math/concepts/real_concept.hpp>
#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
#include <boost/math/special_functions/math_fwd.hpp>
#include <boost/type_traits/is_floating_point.hpp>
#include <boost/array.hpp>
#include "functor.hpp"

#include "handle_test_result.hpp"
#include "table_type.hpp"

#ifndef SC_
#  define SC_(x) static_cast<typename table_type<T>::type>(BOOST_JOIN(x, L))
#endif

template <class T>
T hypergeometric_1f1_int_wrapper(T a, T b, T z)
{
  return static_cast<T>(
    boost::math::hypergeometric_1f1(
      boost::math::itrunc(a),
      boost::math::itrunc(b),
      z)
  );
}

template <class Real, class T>
void do_test_hypergeometric_1f1(const T& data, const char* type_name, const char* test_name)
{
  typedef Real                   value_type;

  typedef value_type (*pg)(value_type, value_type, value_type);
#if defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
  pg funcp = boost::math::hypergeometric_1f1<value_type, value_type, value_type>;
#else
  pg funcp = boost::math::hypergeometric_1f1;
#endif

  boost::math::tools::test_result<value_type> result;

  std::cout << "Testing " << test_name << " with type " << type_name
    << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

  //
  // test hypergeometric_1f1 against data:
  //
  result = boost::math::tools::test_hetero<Real>(
    data, 
    bind_func<Real>(funcp, 0, 1, 2), 
    extract_result<Real>(3));
  handle_test_result(result, data[result.worst()], result.worst(), type_name, "boost::math::hypergeometric_1f1", test_name);
  std::cout << std::endl;
}

template <class Real, class T>
void do_test_hypergeometric_1f1_int(const T& data, const char* type_name, const char* test_name)
{
  typedef Real                   value_type;

  typedef value_type (*pg)(value_type, value_type, value_type);
#if defined(BOOST_MATH_NO_DEDUCED_FUNCTION_POINTERS)
  pg funcp = hypergeometric_1f1_int_wrapper<value_type>;
#else
   pg funcp = hypergeometric_1f1_int_wrapper;
#endif

  boost::math::tools::test_result<value_type> result;

  std::cout << "Testing " << test_name << " with type " << type_name
    << "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n";

  //
  // test hypergeometric_1f1 against data:
  //
  result = boost::math::tools::test_hetero<Real>(
    data,
    bind_func<Real>(funcp, 0, 1),
    extract_result<Real>(2));
  handle_test_result(result, data[result.worst()], result.worst(), type_name, "boost::math::hypergeometric_1f1", test_name);
  std::cout << std::endl;
}

template <class T>
void test_hypergeometric(T, const char* name)
{
  // function values calculated on http://wolframalpha.com/
  static const boost::array<boost::array<T, 4>, 35> pearson_1f1_data = {{
    {{ SC_(0.1), SC_(0.2), SC_(0.5), SC_(1.31762717827850999771120412953367530104071548286484701180675) }},    // basic case, positive a, with b = 2*a
    {{ SC_(-0.1), SC_(0.2), SC_(0.5), SC_(0.695536565102261062224048114156662645634235812108270459351843) }},   // basic case, negative a
    {{ SC_(0.1), SC_(0.2), SC_(-0.5), SC_(0.799181281696560321539893636914544570179421736784530881768694) }},    // basic case, negative z, with b = 2*a
    {{ SC_(1), SC_(1), SC_(1), SC_(2.71828182845904523536028747135266249775724709369995957496697) }},    // exp case
    {{ SC_(10e-8), SC_(10e-8), SC_(10e-10), SC_(1.00000000100000000050000000016666666670833333334166666666806) }},    // small parameters and variable
    {{ SC_(10e-8), SC_(10e-12), SC_(-10e-10), SC_(0.999990000000005000000498283333078775491746764290939543303969) }},  // small parameters and negative variable
    {{ SC_(1), SC_(1), SC_(10), SC_(22026.4657948067165169579006452842443663535126185567810742354) }},    // a == b
    {{ SC_(1), SC_(3), SC_(10), SC_(440.309315896134330339158012905684887327070252371135621484709) }},    // b > a > 0
    {{ SC_(500), SC_(511), SC_(10), SC_(17796.6855333739325171845968940947655234783824325983251918552) }},    // large b > a > 0
    {{ SC_(8.1), SC_(10.1), SC_(100), SC_(1.72413107599268832161436460524695239617795929232918418998156e41) }},    // larger z, with b > a > 0
    {{ SC_(1), SC_(2), SC_(600), SC_(6.28836716821656637233571865580522516197129772710204948904503e257) }},    // very large z, with b = 2*a
    {{ SC_(100), SC_(1.5), SC_(2.5), SC_(2.748892975858683147847923899606490482485018872830787291358434e12) }},    // large positive a
    {{ SC_(-60), SC_(1), SC_(10), SC_(-10.0489541129649484585795209536988461964429445871176041239949) }},    // large negative a, positive z
    {{ SC_(60), SC_(1), SC_(10), SC_(1.81808688761894542863565474130721710227820680443373130533771e22) }},    // large positive a, positive z
    {{ SC_(60), SC_(1), SC_(-10), SC_(-0.000671306684545906746425986256966243643593235570842157239963030) }},    // large positive a, negative z
    {{ SC_(-60), SC_(1), SC_(-10), SC_(1.23314254099858891104546255150303465981995132515458313535635e18) }},    // large negative a, negative z
    {{ SC_(1000), SC_(1), SC_(10e-3), SC_(90.7978833631836412934772519836415738922522669606760796835796) }},    // very large a > 0, small z > 0, z = 1/a
    {{ SC_(10e-3), SC_(1), SC_(700), SC_(1.55801242207709612448488212978843007187046609293269937340995e299) }},    // very large z > 0, small a > 0
    {{ SC_(500), SC_(1), SC_(-5), SC_(0.00105389594336545171918394318988443915771905392258994659221871) }},    // very large a > 0, z < 0
    {{ SC_(-500), SC_(1), SC_(5), SC_(0.251406264291805126115940947726022616964993619603008772953284) }},    // very large a < 0, z > 0
    {{ SC_(20), SC_(-boost::math::constants::pi<T>() * 3), SC_(-2.5), SC_(-459.7875667063203123347166727428484134851591729699138974283514) }},
    {{ SC_(20), SC_(boost::math::constants::pi<T>() * 3), SC_(2.5), SC_(122.18665421816610412294237478730436434907793724942950) }},
    {{ SC_(-20), SC_(-boost::math::constants::pi<T>() * 3), SC_(2.5), SC_(269.066377517348825230568727303741204558261168571043416124554) }},
    {{ SC_(50), SC_(10), SC_(200), SC_(1.70970417612603891522605058021030241438424161666807852280059e125) }},
    {{ SC_(-5), SC_(-boost::math::constants::pi<T>() * 2), SC_(-1), SC_(0.444913437706720747615165269249706557611011297777085582480832) }},
    {{ SC_(4), SC_(80), SC_(200), SC_(3.44855150621665366519268801773620305050820893622757983240575e27) }},                // large b and larger z
    {{ SC_(-4), SC_(500), SC_(300), SC_(0.0249062013158541889868414106813659161997265162513954135719595) }},              // large z and larger b
    {{ SC_(5), SC_(0.1), SC_(-2), SC_(2.69810607783535476510057484623671017756932808501877574485808) }},
    {{ SC_(-5), SC_(0.1), SC_(2), SC_(19.1331126256381960551905477161181331126256381960551905477161) }},
    {{ SC_(5), SC_(2), SC_(100), SC_(1.25851372306500557378632832198689216128543415504279140966030e48) }},
    {{ SC_(-5), SC_(2), SC_(-100), SC_(1.84891398888888888888888888888888888888888888888888888888889e7) }},
    {{ SC_(1), SC_(10e-12), SC_(1), SC_(2.71828182844739141320726036451388968000401036004750004644223e11) }},       // very small b
    {{ SC_(10), SC_(10e-12), SC_(10), SC_(1.332534440738693105395841520859877375366639491043385439173436e22) }},       // very small b with larger a and z
    {{ SC_(-1000), SC_(1), SC_(1000), SC_(-2.59382078336200571793976408157920288042214518531712673813455e215) }},       // very large negative real a and very large positive real z
    {{ SC_(20), SC_(10), SC_(-5), SC_(-6.25682720117872755719988465952783439826492132332740217148550e-6) }},
  }};

  do_test_hypergeometric_1f1<T>(pearson_1f1_data, name, "John Pearson dissertation: WolframAlpha Data");

#include "hypergeometric_1f1_luke_rational_data.ipp"
  do_test_hypergeometric_1f1<T>(hypergeometric_1f1_luke_rational_data, name, "Random positive data with b >= (10 * a) (rational case)");
#include "hypergeometric_1f1_luke_pade_moderate_data.ipp"
  do_test_hypergeometric_1f1<T>(hypergeometric_1f1_luke_pade_moderate_data, name, "Random data with a == 1 (pade case)");
#include "hypergeometric_1f1_moderate_data.ipp"
  do_test_hypergeometric_1f1<T>(hypergeometric_1f1_moderate_data, name, "Random moderate data");
}

