#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

void outlierCpp(const int K, const int R, NumericMatrix xy, NumericMatrix ratio, NumericMatrix imat, NumericVector rmin);

TEST(Benchmarking_deepstate_test,outlierCpp_test){
  std::ofstream K_stream;
  std::ofstream R_stream;
  std::ofstream xy_stream;
  std::ofstream ratio_stream;
  std::ofstream imat_stream;
  std::ofstream rmin_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  int K  = RcppDeepState_int();
  K_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Benchmarking/inst/testfiles/outlierCpp/inputs/K");
  K_stream << K;
  std::cout << "K values: "<< K << std::endl;
  K_stream.close();
  int R  = RcppDeepState_int();
  R_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Benchmarking/inst/testfiles/outlierCpp/inputs/R");
  R_stream << R;
  std::cout << "R values: "<< R << std::endl;
  R_stream.close();
  NumericMatrix xy  = RcppDeepState_NumericMatrix();
  xy_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Benchmarking/inst/testfiles/outlierCpp/inputs/xy");
  xy_stream << xy;
  std::cout << "xy values: "<< xy << std::endl;
  xy_stream.close();
  NumericMatrix ratio  = RcppDeepState_NumericMatrix();
  ratio_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Benchmarking/inst/testfiles/outlierCpp/inputs/ratio");
  ratio_stream << ratio;
  std::cout << "ratio values: "<< ratio << std::endl;
  ratio_stream.close();
  NumericMatrix imat  = RcppDeepState_NumericMatrix();
  imat_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Benchmarking/inst/testfiles/outlierCpp/inputs/imat");
  imat_stream << imat;
  std::cout << "imat values: "<< imat << std::endl;
  imat_stream.close();
  NumericVector rmin  = RcppDeepState_NumericVector();
  rmin_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Benchmarking/inst/testfiles/outlierCpp/inputs/rmin");
  rmin_stream << rmin;
  std::cout << "rmin values: "<< rmin << std::endl;
  rmin_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    outlierCpp(K,R,xy,ratio,imat,rmin);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
