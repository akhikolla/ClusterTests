#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export()]]
Rcpp::List Sinkhorn(arma::mat a, arma::mat b, arma::mat costm, double lambda = 1, int maxIter = 10000, double tolerance = 10^(-8)) {

  lambda = 1/lambda;
  // Transforming Input
  arma::mat K = exp(-lambda*costm);
  arma::mat U = K % costm;

  //Checking the type of computation: 1-vs-N points or many pairs
  bool ONE_VS_N;
  if(a.n_cols==1){
    ONE_VS_N = true;
  }
  else{
    ONE_VS_N = false;
  }

  //Checking dimensionality
  bool BIGN;
  if(b.n_cols>b.n_rows){
    BIGN = true;
  }
  else{
    BIGN = false;
  }

  //Checking for the case that a transportplan can be returned, i.e. a and b only represent one marginal, respectively
  bool TRANSPORT = false;
  if((ONE_VS_N == true) & (b.n_cols == 1)){
    TRANSPORT = true;
  }

  //Initializing certain matrices for the Sinkhorn algorithm
  arma::mat ainvK;
  arma::mat v;

  //Need to update some vectors and matrices if a does not have full support
  bool ZeroValues = false;
  arma::uvec ids;
  arma::mat Ktemp;
  if(ONE_VS_N == true){
    arma::vec a1 = a;
    if(all(a1)==false){
      ZeroValues = true;
      Ktemp = K;  //save original matrix K to compute the transportplan at the end
      ids = find(a1 > 0); // Find indices with support
      K = K.rows(ids);
      a = a.rows(ids);
      U = U.rows(ids);
    }
    ainvK = diagmat(1/a) * K;//Precomputation of this matrix saves a d1 x N Schur product at each iteration
  }

  //Fixed point counter
  int iteration=0;

  //Initialization of left scaling Factors and the matrix v for the Sinkhorn algorithm
  arma::mat u(a.n_rows,b.n_cols);
  u = u.ones()/a.n_rows;

  // main iteration of Sinkhorn's algorithm
  while(iteration < maxIter){
    if(ONE_VS_N==true){
      if(BIGN==true){
        u = 1 / (ainvK * (b / (K.t()*u)));
      }
      else{
        u = 1 / (ainvK * (b / (u.t()*K).t() ));
      }
    }
    else{
      if(BIGN==true){
        u = a / (K*(b/(u.t()*K).t()));
      }
      else{
        u = a / (K*(b / (K.t()*u)));
      }
    }
    iteration++;
    if(iteration % 20 == 0 || iteration == maxIter){
      if(BIGN==true){
        v = b / (K.t()*u);
      }
      else{
        v = b / ((u.t()*K).t());
      }
      if(ONE_VS_N==true){
        u = 1 / (ainvK*v);
      }
      else{
        u = a / (K*v);
      }
      //Checking the stopping criterion with marginal differences
      double Criterion = norm(sum(abs(v % (K.t()*u)-b)));
      if(Criterion < tolerance){
        break;
      }
    }


  }

  arma::mat Distance = sum(u % (U * v)) ; //computing Sinkhorn Distance

  // compute the transportplan if possible, in case ZeroValues is TRUE we have to add artificial zero values
  arma::mat Transportplan;
  arma::mat uu(b.n_rows,1);
  if(TRANSPORT == true){
    if(ZeroValues == false){
      Transportplan = diagmat(u)*K*diagmat(v);
    }
    else{
      uu.rows(ids) = u;
      Transportplan = diagmat(uu)*Ktemp*diagmat(v);
    }
    return Rcpp::List::create(
      Rcpp::Named("Transportplan") = Transportplan,
      Rcpp::Named("Distance") = Distance,
      Rcpp::Named("alpha") = lambda*log(u),
      Rcpp::Named("beta") = lambda*log(v)
    );
  }
  else{
    return Rcpp::List::create(
      Rcpp::Named("Distance") = Distance
    );
  }


  //arma::mat Transportplan = diagmat(u)*K*diagmat(v);

  //arma::mat alt(u.n_rows,1);
  //alt = alt.ones();

  //return Distance ;
  //return Rcpp::List::create(
   // Rcpp::Named("Distance") = Distance,
  //  Rcpp::Named("Transportplan") = Transportplan
  //);

  //return Rcpp::List::create(
    //Rcpp::Named("Distance") = Distance,
    //Rcpp::Named("u") = u,
   // Rcpp::Named("v") = v
   // Rcpp::Named("Dual_Solution") = -log(u)/lambda + accu(log(u))/(lambda*u.n_rows) * alt
  //)
   // ;
}
