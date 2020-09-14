#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

IntegerVector PELT_meanvar_norm2(NumericVector resp, double pen);

TEST(bartBMA_deepstate_test,PELT_meanvar_norm2_test){
  std::ofstream resp_stream;
  std::ofstream pen_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericVector resp  = RcppDeepState_NumericVector();
  resp_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/PELT_meanvar_norm2/inputs/resp");
  resp_stream << resp;
  std::cout << "resp values: "<< resp << std::endl;
  resp_stream.close();
  double pen  = RcppDeepState_double();
  pen_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bartBMA/inst/testfiles/PELT_meanvar_norm2/inputs/pen");
  pen_stream << pen;
  std::cout << "pen values: "<< pen << std::endl;
  pen_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    PELT_meanvar_norm2(resp,pen);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
