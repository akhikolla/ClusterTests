#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericVector densityMixture(NumericVector weights, NumericMatrix densities);

TEST(beyondWhittle_deepstate_test,densityMixture_test){
  std::ofstream weights_stream;
  std::ofstream densities_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector weights  = RcppDeepState_NumericVector();
  weights_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/densityMixture/inputs/weights");
  weights_stream << weights;
  std::cout << "weights values: "<< weights << std::endl;
  weights_stream.close();
  NumericMatrix densities  = RcppDeepState_NumericMatrix();
  densities_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/beyondWhittle/inst/testfiles/densityMixture/inputs/densities");
  densities_stream << densities;
  std::cout << "densities values: "<< densities << std::endl;
  densities_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    densityMixture(weights,densities);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
