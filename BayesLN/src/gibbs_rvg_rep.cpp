//C++ implementation of a Gibbs sampler for a normal hierarchical model with only one variance component as random effect.


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


List gibbs_rep(arma::vec y, arma::mat X, arma::mat Z,arma::mat Xrep, arma::mat Zrep,
           double l_s,double l_t,
           double d_s,double d_t,
           double g_s,double g_t,
           int nsamp, int verbose,
           arma::mat beta, arma::vec sigma2, arma::vec tau2){

  //dimensions
  int n = X.n_rows;
  int p = X.n_cols;
  int m = Z.n_cols;

  //statistics
  arma::mat XtX = trans(X) * X;
  arma::mat ZtZ = trans(Z) * Z;
  arma::mat pinvZtZ = arma::pinv(ZtZ);
  arma::vec beta_hat = inv(XtX)*trans(X)*y;
  arma::mat Pz = Z * pinvZtZ * trans(Z);

  //parameters samples
  arma::mat u = arma::zeros<arma::mat>(nsamp, m);
  arma::mat I = arma::eye(m,m);



  //inits parameters
  u.row(0) = arma::zeros<arma::mat>(1, m);


    int nrep = Xrep.n_rows;
    arma::mat yrep = arma::zeros<arma::mat>(nsamp, nrep);
    yrep.row(0) = arma::trans(mvrnorm_Ind(Xrep * arma::trans(beta.row(0)) + Zrep * arma::trans(u.row(0)), sigma2(0) ));

  //initialization of useful components
  arma::mat Vz = arma::zeros<arma::mat>(n, n);
  arma::mat Vbeta = arma::zeros<arma::mat>(p, p);
  arma::vec beta_bar = arma::zeros<arma::vec>(p);
  arma::mat Vu = arma::zeros<arma::mat>(m, m);
  arma::mat u_bar = arma::zeros<arma::vec>(m);
  arma::vec beta_tilde = arma::zeros<arma::vec>(p);
  double uu = 0.0;
  double SS = 0.0;

  //loop
  for(int i = 1; i < nsamp; i ++){
    //check user interruption every 1000 iterations
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    //beta
    Vz = Z * pinvZtZ * inv(pinvZtZ + I * (tau2(i-1) / sigma2(i-1))) * pinvZtZ * trans(Z);
    Vbeta = inv((1.0 / sigma2(i-1)) * XtX + trans(X) * ((1.0 / sigma2(i-1)) * (Vz-Pz)) * X );
    beta_tilde = inv(trans(X) * ((1.0 / sigma2(i-1)) * (Vz - Pz)) * X) * trans(X) * ((1.0 / sigma2(i-1)) * (Vz - Pz)) * y;
    beta_bar = Vbeta * ((1.0 / sigma2(i-1)) * XtX * beta_hat + trans(X) * ((1.0 / sigma2(i-1)) * (Vz - Pz)) * X * beta_tilde );
    beta.row(i) = trans(mvrnorm(beta_bar, Vbeta));

    //u
    u_bar = inv(ZtZ + (sigma2(i-1) / tau2(i-1)) * I) * trans(Z) * (y - X * trans(beta.row(i)));
    Vu = sigma2(i-1) * inv(ZtZ + (sigma2(i-1) / tau2(i-1)) * I);
    u.row(i) =trans(mvrnorm(u_bar, Vu));

    //tau2
    uu = arma::as_scalar(u.row(i) * trans(u.row(i)));
    tau2(i) = rgig(l_t - m / 2.0, d_t *d_t + uu, g_t * g_t);

    //sigma2
    SS = arma::as_scalar(trans(y - X * trans(beta.row(i)) - Z * trans(u.row(i))) * (y - X * trans(beta.row(i)) - Z * trans(u.row(i))));
    sigma2(i) = rgig(l_s - n / 2.0, d_s * d_s+ SS, g_s * g_s);

    // Include the the deviance


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

    return List::create(Named("beta") = beta, Named("u") = u, Named("tau2") = tau2, Named("sigma2") = sigma2, Named("yrep") = yrep);

}











// [[Rcpp::export]]


List gibbs_rep_r(arma::vec y, arma::mat X, arma::mat Z,arma::mat Xrep, arma::mat Zrep,
                 double l_s,arma::vec l_t,
                 double d_s,arma::vec d_t,
                 double g_s,arma::vec g_t,
                 int r, arma::vec m_s_cum, arma::vec m_s,
               int nsamp, int verbose,
               arma::mat beta, arma::vec sigma2, arma::mat tau2){

  //dimensions
  int n = X.n_rows;
  int p = X.n_cols;
  int m = Z.n_cols;

  //statistics
  arma::mat XtX = trans(X) * X;
  arma::mat ZtZ = trans(Z) * Z;
  arma::mat pinvZtZ = arma::pinv(ZtZ);
  arma::vec beta_hat = inv(XtX)*trans(X)*y;
  arma::mat Pz = Z * pinvZtZ * trans(Z);

  //parameters samples
  arma::mat u = arma::zeros<arma::mat>(nsamp, m);
  arma::mat D = arma::eye(m,m);//inits D

  //inits parameters
  u.row(0) = arma::zeros<arma::mat>(1, m);

  int nrep = Xrep.n_rows;
  arma::mat yrep = arma::zeros<arma::mat>(nsamp, nrep);
  yrep.row(0) = arma::trans(mvrnorm_Ind(Xrep * arma::trans(beta.row(0)) + Zrep * arma::trans(u.row(0)), sigma2(0) ));


  //initialization of useful components
  arma::mat Vz = arma::zeros<arma::mat>(n, n);
  arma::mat Vbeta = arma::zeros<arma::mat>(p, p);
  arma::vec beta_bar = arma::zeros<arma::vec>(p);
  arma::mat Vu = arma::zeros<arma::mat>(m, m);
  arma::mat u_bar = arma::zeros<arma::vec>(m);
  arma::vec beta_tilde = arma::zeros<arma::vec>(p);
  arma::vec D_d = arma::ones<arma::vec>(m);
  double uus ;
  double SS ;
  arma::vec us = arma::zeros<arma::vec>(m);


  //loop
  for(int i = 1; i < nsamp; i ++){
    //check user interruption every 1000 iterations
    if (i % 1000 == 0)
      Rcpp::checkUserInterrupt();

    //beta
    Vz = Z * pinvZtZ * inv(pinvZtZ + D / sigma2(i-1)) * pinvZtZ * trans(Z);
    Vbeta = inv((1.0 / sigma2(i-1)) * XtX + trans(X) * ((1.0 / sigma2(i-1)) * (Vz-Pz)) * X );
    beta_tilde = inv(trans(X) * ((1.0 / sigma2(i-1)) * (Vz - Pz)) * X) * trans(X) * ((1.0 / sigma2(i-1)) * (Vz - Pz)) * y;
    beta_bar = Vbeta * ((1.0 / sigma2(i-1)) * XtX * beta_hat + trans(X) * ((1.0 / sigma2(i-1)) * (Vz - Pz)) * X * beta_tilde );
    beta.row(i) = trans(mvrnorm(beta_bar, Vbeta));

    //u
    u_bar = inv(ZtZ + (sigma2(i-1)) * inv(D)) * trans(Z) * (y - X * trans(beta.row(i)));
    Vu = sigma2(i-1) * inv(ZtZ + sigma2(i-1) * inv(D));
    u.row(i) =trans(mvrnorm(u_bar, Vu));



    us.resize(m_s(0));
    for(int k =0; k<m_s(0); k++){
      us(k)= u(i,k);
    }

    uus = arma::as_scalar(trans(us)*us);
    tau2(i,0) = rgig(l_t(0) - m_s(0) / 2.0, d_t(0) * d_t(0) + uus, g_t(0) * g_t(0));
    //tau2
    /////////partizionare u
    for(int j = 1; j < r; j++){
      us.resize(m_s(j));
      for(int k = m_s_cum(j-1); k<m_s_cum(j); k++){
        us(k-m_s_cum(j-1))= u(i,k);
      }
      uus = arma::as_scalar(trans(us)*us);
      tau2(i,j) = rgig(l_t(j) - m_s(j) / 2.0, d_t(j) * d_t(j) + uus, g_t(j) * g_t(j));
    }


    for(int j=0; j<m; j++){
      if(j<m_s_cum(0)){
        D_d(j)=tau2(i,0);
      }else{
        for(int s=1;s<r; s++){
          if(j>=m_s_cum(s-1)&&j<m_s_cum(s)){
            D_d(j)=tau2(i,s);
          }else{
            ;
          }
        }
      }
    }

    D.diag() = D_d;

    //sigma2
    SS = arma::as_scalar(trans(y - X * trans(beta.row(i)) - Z * trans(u.row(i))) * (y - X * trans(beta.row(i)) - Z * trans(u.row(i))));
    sigma2(i) = rgig(l_s - n / 2.0, d_s * d_s + SS, g_s * g_s);

    // Include the the deviance

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
  return List::create(Named("beta") = beta, Named("u") = u, Named("tau2") = tau2, Named("sigma2") = sigma2, Named("yrep") = yrep);

}









