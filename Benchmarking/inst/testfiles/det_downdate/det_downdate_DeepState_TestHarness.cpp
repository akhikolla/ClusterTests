#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double det_downdate(const NumericMatrix A, const NumericVector v, const double det);

TEST(Benchmarking_deepstate_test,det_downdate_test){
  std::ofstream A_stream;
  std::ofstream v_stream;
  std::ofstream det_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix A  = RcppDeepState_NumericMatrix();
  A_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Benchmarking/inst/testfiles/det_downdate/inputs/A");
  A_stream << A;
  std::cout << "A values: "<< A << std::endl;
  A_stream.close();
  NumericVector v  = RcppDeepState_NumericVector();
  v_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Benchmarking/inst/testfiles/det_downdate/inputs/v");
  v_stream << v;
  std::cout << "v values: "<< v << std::endl;
  v_stream.close();
  double det  = RcppDeepState_double();
  det_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/Benchmarking/inst/testfiles/det_downdate/inputs/det");
  det_stream << det;
  std::cout << "det values: "<< det << std::endl;
  det_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    det_downdate(A,v,det);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
