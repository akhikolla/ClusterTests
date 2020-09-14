#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

List vusC_full_core(NumericVector tt1, NumericVector tt2, NumericVector tt3);

TEST(bcROCsurface_deepstate_test,vusC_full_core_test){
  std::ofstream tt1_stream;
  std::ofstream tt2_stream;
  std::ofstream tt3_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector tt1  = RcppDeepState_NumericVector();
  tt1_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bcROCsurface/inst/testfiles/vusC_full_core/inputs/tt1");
  tt1_stream << tt1;
  std::cout << "tt1 values: "<< tt1 << std::endl;
  tt1_stream.close();
  NumericVector tt2  = RcppDeepState_NumericVector();
  tt2_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bcROCsurface/inst/testfiles/vusC_full_core/inputs/tt2");
  tt2_stream << tt2;
  std::cout << "tt2 values: "<< tt2 << std::endl;
  tt2_stream.close();
  NumericVector tt3  = RcppDeepState_NumericVector();
  tt3_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bcROCsurface/inst/testfiles/vusC_full_core/inputs/tt3");
  tt3_stream << tt3;
  std::cout << "tt3 values: "<< tt3 << std::endl;
  tt3_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    vusC_full_core(tt1,tt2,tt3);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
