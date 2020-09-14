
#include <RcppArmadillo.h>
#include <R.h>
#include <Rmath.h>
#include <Rdefines.h>
#include "nor_generator.h"
extern "C"
{
#include "gig_generator.h"
}


using namespace Rcpp;



// [[Rcpp::export]]


List predictive_gen(arma::mat Xrep, arma::mat Zrep,
               arma::mat beta, arma::mat u, arma::vec sigma2,
               int nsamp, int verbose){


  int nrep = Xrep.n_rows;
  arma::mat yrep = arma::zeros<arma::mat>(nsamp, nrep);
  yrep.row(0) = arma::trans(mvrnorm_Ind(Xrep * arma::trans(beta.row(0)) + Zrep * arma::trans(u.row(0)), sigma2(0) ));



  //loop
  for(int i = 1; i < nsamp; i ++){
    //check user interruption every 1000 iterations
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    //generation posterior
    yrep.row(i) = trans(mvrnorm_Ind(Xrep * trans(beta.row(i)) + Zrep * trans(u.row(i)), sigma2(i)));


    //Progress bar
    double Perc=100*(i+1)/(nsamp);
    if(((i+1)%(nsamp/100))==0 && (verbose==1)){
      Rprintf("-");
      R_FlushConsole();
      //R_ProcessEvents(); for windows
      if(((i+1)%(nsamp/10))==0){
        Rprintf(":%.1f%%\n",Perc);
        R_FlushConsole();
        //R_ProcessEvents(); for windows
      }
    }
  }

  return List::create(Named("yrep") = yrep);

}




















