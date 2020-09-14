#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

IntegerVector mySeq(int start, int end);

TEST(bigReg_deepstate_test,mySeq_test){
  std::ofstream start_stream;
  std::ofstream end_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  int start  = RcppDeepState_int();
  start_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bigReg/inst/testfiles/mySeq/inputs/start");
  start_stream << start;
  std::cout << "start values: "<< start << std::endl;
  start_stream.close();
  int end  = RcppDeepState_int();
  end_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/bigReg/inst/testfiles/mySeq/inputs/end");
  end_stream << end;
  std::cout << "end values: "<< end << std::endl;
  end_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    mySeq(start,end);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
