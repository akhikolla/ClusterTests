#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

double Window_helper(int b, NumericVector xb, NumericVector y, NumericVector t, NumericVector Lfactors, NumericVector J1_b, NumericVector J2_b);

TEST(binnednp_deepstate_test,Window_helper_test){
  std::ofstream b_stream;
  std::ofstream xb_stream;
  std::ofstream y_stream;
  std::ofstream t_stream;
  std::ofstream Lfactors_stream;
  std::ofstream J1_b_stream;
  std::ofstream J2_b_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  int b  = RcppDeepState_int();
  b_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/Window_helper/inputs/b");
  b_stream << b;
  std::cout << "b values: "<< b << std::endl;
  b_stream.close();
  NumericVector xb  = RcppDeepState_NumericVector();
  xb_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/Window_helper/inputs/xb");
  xb_stream << xb;
  std::cout << "xb values: "<< xb << std::endl;
  xb_stream.close();
  NumericVector y  = RcppDeepState_NumericVector();
  y_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/Window_helper/inputs/y");
  y_stream << y;
  std::cout << "y values: "<< y << std::endl;
  y_stream.close();
  NumericVector t  = RcppDeepState_NumericVector();
  t_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/Window_helper/inputs/t");
  t_stream << t;
  std::cout << "t values: "<< t << std::endl;
  t_stream.close();
  NumericVector Lfactors  = RcppDeepState_NumericVector();
  Lfactors_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/Window_helper/inputs/Lfactors");
  Lfactors_stream << Lfactors;
  std::cout << "Lfactors values: "<< Lfactors << std::endl;
  Lfactors_stream.close();
  NumericVector J1_b  = RcppDeepState_NumericVector();
  J1_b_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/Window_helper/inputs/J1_b");
  J1_b_stream << J1_b;
  std::cout << "J1_b values: "<< J1_b << std::endl;
  J1_b_stream.close();
  NumericVector J2_b  = RcppDeepState_NumericVector();
  J2_b_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/binnednp/inst/testfiles/Window_helper/inputs/J2_b");
  J2_b_stream << J2_b;
  std::cout << "J2_b values: "<< J2_b << std::endl;
  J2_b_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    Window_helper(b,xb,y,t,Lfactors,J1_b,J2_b);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
