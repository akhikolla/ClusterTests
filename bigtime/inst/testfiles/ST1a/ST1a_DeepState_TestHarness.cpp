#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double ST1a(double z, double gam);

TEST(bigtime_deepstate_test,ST1a_test){
  std::ofstream z_stream;
  std::ofstream gam_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  double z  = RcppDeepState_double();
  z_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bigtime/inst/testfiles/ST1a/inputs/z");
  z_stream << z;
  std::cout << "z values: "<< z << std::endl;
  z_stream.close();
  double gam  = RcppDeepState_double();
  gam_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bigtime/inst/testfiles/ST1a/inputs/gam");
  gam_stream << gam;
  std::cout << "gam values: "<< gam << std::endl;
  gam_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    ST1a(z,gam);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
