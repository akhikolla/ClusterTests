#include <fstream>
#include <RInside.h>
#include <iostream>
#include <RcppDeepState.h>
#include <DeepState.hpp>

NumericMatrix updateSubmat(NumericMatrix submat, NumericMatrix sublb_red);

TEST(BacArena_deepstate_test,updateSubmat_test){
  std::ofstream submat_stream;
  std::ofstream sublb_red_stream;
  RInside();
  std::cout << "input starts" << std::endl;
  NumericMatrix submat  = RcppDeepState_NumericMatrix();
  submat_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BacArena/inst/testfiles/updateSubmat/inputs/submat");
  submat_stream << submat;
  std::cout << "submat values: "<< submat << std::endl;
  submat_stream.close();
  NumericMatrix sublb_red  = RcppDeepState_NumericMatrix();
  sublb_red_stream.open("/home/akolla/R/x86_64-pc-linux-gnu-library/3.6/RcppDeepState/extdata/compileAttributes/BacArena/inst/testfiles/updateSubmat/inputs/sublb_red");
  sublb_red_stream << sublb_red;
  std::cout << "sublb_red values: "<< sublb_red << std::endl;
  sublb_red_stream.close();
  std::cout << "input ends" << std::endl;
  try{
    updateSubmat(submat,sublb_red);
  }
  catch(Rcpp::exception& e){
    std::cout<<"Exception Handled"<<std::endl;
  }
}
