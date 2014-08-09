//  Copyright (c) 2007 John Maddock
//  Copyright (c) 2014 Anton Bikineev
//  Use, modification and distribution are subject to the
//  Boost Software License, Version 1.0. (See accompanying file
//  LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
//
// This program computes hypergeometric 1f1 function using
// naive Taylor series at quite large internal precision -
// 1000 decimal digits.
//
//
#include <fstream>

#include <boost/multiprecision/mpfr.hpp>
#include <boost/math/tools/test_data.hpp>
#include <boost/test/included/prg_exec_monitor.hpp>
#include <boost/math/special_functions/hypergeometric.hpp>

using namespace boost::math;
using namespace boost::math::tools;
using namespace boost::multiprecision;
using namespace std;

template <class T>
T hypergeometric_1f1_taylor_bare(const T& a, const T& b, const T& z)
{
   const boost::math::policies::policy<> pol;
   return boost::math::detail::hypergeometric_1f1_generic_series(a, b, z, pol);
}

int cpp_main(int argc, char* argv[])
{
   typedef number<mpfr_float_backend<1000u> > float_type;

   parameter_info<float_type> arg1, arg2, arg3;
   test_data<float_type> data;

   bool cont;
   std::string line;

   std::cout << "Welcome.\n"
      "This program will generate spot tests for confluent hypergeometric 1f1 function\n\n";
   do{
      get_user_parameter_info(arg1, "a");
      get_user_parameter_info(arg2, "b");
      get_user_parameter_info(arg3, "z");
      float_type (*fp)(const float_type&, const float_type&, const float_type&) = hypergeometric_1f1_taylor_bare;
      data.insert(fp, arg1, arg2, arg3);

      std::cout << "Any more data [y/n]?";
      std::getline(std::cin, line);
      boost::algorithm::trim(line);
      cont = (line == "y");
   }while(cont);

   std::cout << "Enter name of test data file [default=hypergeometric_1f1.ipp]";
   std::getline(std::cin, line);
   boost::algorithm::trim(line);
   if(line == "")
      line = "hypergeometric_1f1.ipp";
   std::ofstream ofs(line.c_str());
   line.erase(line.find('.'));
   ofs << std::scientific << std::setprecision(60u);
   write_code(ofs, data, line.c_str());

   return 0;
}
