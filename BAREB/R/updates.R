#' update_w_beta
#'
#' This function updates the weights for each patient level cluster
#'
#' @param S The number of patient level clusters
#' @param E A vector that records the current clustering membership.
#' @param hyper_delta The hyper-parameter with default value being 1
#' @return The updated weights for each patient level cluster
#' @seealso \link{update_RJ} for a complete example for all functions in this package.
#' @examples
#' #Suppose we know the number of patient level cluster is 4
#' #Suppose the current clustering membership indicates 3 patients in cluster 1,
#' #2 patients in cluster 2, 3 patinets in cluster 3, 1 patient in cluster 4
#' #Use the default value, 1, for hyper-parameter
#' update_w_beta(S=4,E=c(1,1,1,2,2,3,3,3,4))
#'
#'
#' #To change the hyper-parameter to, for example 2
#' update_w_beta(S=4,E=c(1,1,1,2,2,3,3,3,4),hyper_delta = 2)
#' @export
update_w_beta <- function(S,  E, hyper_delta = 1){
  w_star = rep(NA, S)
  for (l in 1:S){
    Dir_par = hyper_delta + sum((E == l))
    w_star[l] = rgamma(1, shape = Dir_par, rate = 1)
  }
  w_star = w_star/sum(w_star)
  return(w_star)
}





#' update_w
#'
#' This function updates the weights for site level clusters
#'
#' @param K The number of site level cluster for each patient level cluster. Should be a vector of length \code{S}
#' @param R The current site level clustering membership. Should be a matrix with \code{S} rows.
#' @param S The number of patient level clusters.
#' @param hyper_delta The hyper-parameter with default value being 1
#' @return The updated weights for site level clusters. Should be a matrix with \code{S} rows.
#' @seealso \link{update_RJ} for a complete example for all functions in this package.
#' @examples
#' #Suppose we know the number of patient level cluster is 2,
#' #one has 2 site level clusters and one has 3.
#' #Use the default value, 1, for hyper-parameter
#' update_w(K=c(2,3),R=matrix(c(1,2,2,1,1,2,3,2),nrow=2,byrow=TRUE), S=2)
#'
#' #To change the hyper-parameter to, for example 2
#' update_w(K=c(2,3),R=matrix(c(1,2,2,1,1,2,3,2),nrow=2,byrow=TRUE), S=2, hyper_delta = 2)
#' @details
#' It returns a matrix with 10 columns.
#' For example, first patient cluster has 2 site level clusters. The first row's first 2 elements give the weights for site level clusters in patient cluster 1. Last 8 elements are NA's
#' @export
update_w <- function(K, R, S, hyper_delta = 1){
  #S is the number of patient level clusters
  #K should be a vector of length S
  ans <- matrix(NA, nrow = S, ncol = 10)
  for(i in 1:S){

    w_star = rep(NA, K[i])
    for (l in 1:K[i]){
      Dir_par = hyper_delta + sum((R[i,] == l))
      w_star[l] = rgamma(1, shape = Dir_par, rate = 1)
    }
    w_star = w_star/sum(w_star)
    ans[i, 1:K[i]]<- w_star
  }
  return(ans)
}


#' update_RJ
#'
#' Update the number of site level clusters for each patient level cluster via RJMCMC
#'
#' @param w The weights for site level clusters. Should be a matrix with S rows.
#' @param K The number of site level cluster for each patient level cluster. Should be a vector of length S
#' @param Gamma The linear coefficients for site level covariates. Should be a 3-dimensional array.
#' @param Beta The linear coefficients for patient level covariates. Should be a matrix with S rows
#' @param E A vector that records the current clustering membership.
#' @param Z The design matrix for site level clusters.
#' @param X The design matrix for patient level clusters.
#' @param R The current site level clustering membership
#' @param mu The current CAL mean matrix
#' @param mu_star The latent value matrix for missingness model
#' @param Y The observed CAL value matrix
#' @param delta The missing indicator matrix
#' @param c The linear coefficients for missingness model
#' @param sigma_square The current noise variance
#' @param C The DPP related kernel matrices. Should be an array of 3 dimensions
#' @param S The number of patient level clusters.
#' @param theta The DPP hyper-parameter for site level
#' @param tau A fixed DPP hyper-parameter, which we suggest high value, say 10^5
#' @param q The number of site level covariates
#' @param m The number of sites
#' @param T0 The number of teeth
#' @param hyper_delta The hyper-parameter with default value being 1
#' @return A list with following updated parameters:
#' @return \item{K}{The numbers of site level clusters}
#' @return \item{w}{The weights for site level clusters}
#' @return \item{Gamma}{Linear coefficients for site level covariates}
#' @return \item{R}{The site level clusteirng membership}
#' @return \item{C}{The DPP related kernel matrice}
#' @examples
#' library(BAREB)
#' data("obs")
#' X <- obs$X
#' Y <- obs$Y
#' Z <- obs$Z
#' delta <- obs$delta
#' data("truth")
#'
#' set.seed(1)
#'
#' n<-80
#' m<-168
#' T0<-28
#' q<-3
#' p<-3
#' S<-3
#' theta1 <- theta2 <- 5
#' tau <- 100000
#' D<-matrix(0, nrow = T0, ncol = m)
#' for(i in 1:T0){
#'   indi<-1:6
#'   indi<-indi+6*(i-1)
#'   D[i,indi]<-rep(1/6,6)
#' }
#' nu_gamma<-0.05
#' nu_beta <- 0.05
#' data("Init")
#' Beta0 <- Init$Beta
#' Gamma0 <- Init$Gamma
#'
#' Niter<-10
#' record <- NULL
#' record$E<-matrix(NA,nrow = Niter, ncol = n)
#' record$R<-array(NA, dim = c(Niter, S, m))
#' record$Gamma <-array(NA,dim = c(Niter, 10, q, S))
#' record$Beta <- array(NA,dim = c(Niter, S, p))
#' record$K <- matrix(NA,nrow = Niter, ncol = S)
#' record$sigma_square <-rep(NA,Niter)
#' record$theta1<-rep(0,Niter)
#' record$theta2<-rep(0,Niter)
#' record$c<-matrix(0,nrow =  Niter,ncol = T0)
#' record$mu<-array(NA,dim = c(Niter,n,m))
#' record$w_beta <- array(NA, dim = c(Niter, S))
#' record$w <- array(NA, dim = c(Niter, S, 10))
#' set.seed(1)
#' E<-sample.int(S,n,TRUE)
#' Beta <- matrix(NA,nrow = S, ncol = p)
#' Beta[1,] <- Beta[2,] <- Beta[3,]  <- Beta0
#' Beta<- Beta + matrix(rnorm(S*p, 0, 1) , nrow = S, ncol = p)
#' Gamma <- array(NA,dim = c(10,q,S))
#' Gamma[1,,1] <- Gamma[2,,1] <- Gamma[3,,1] <- Gamma0
#' Gamma[,,2] <- Gamma[,,3]   <- Gamma[,,1]
#' Gamma <- Gamma + array(rnorm(10*q*S, 0, 5), dim = c(10, q, S))
#' K <- rep(3,S)
#' R <- matrix(NA,nrow = S, ncol = m)
#' R[1,]<- R[2,]<- R[3,]  <- sample.int(3,m,TRUE)
#' mu<-updatemu(R,Z,X,Gamma,K,Beta,E,m,n,p,q)
#' mu_star<-updatemustar(mu,rep(0.01,2),n,T0,D)
#' z_star<-updateZstar(mu_star,delta,n,T0)
#' sigma_square <- 10
#' C<-array(NA,dim = c(10,10,S))
#' w<-matrix(NA,nrow = S, ncol = 10)
#' w_beta<-rep(1/S,S)
#' for(i in 1:S){
#'  C[1:K[i],1:K[i],i]<-updateC(Gamma[1:K[i],,i],theta2,tau)
#'   w[i, 1:K[i]]<-rep(1/K[i],K[i])
#' }
#' c<-c(0,0.01)
#' start <- Sys.time()
#' for(iter in 1:Niter){
#'   c<-updatec(z_star, mu,D, 100,sigma_square, n, T0)
#'   mu_star<-updatemustar(mu,c,n,T0,D)
#'   z_star<-updateZstar(mu_star,delta,n, T0)
#'   w <- update_w(K, R, S)
#'   R <- updateR(w, Gamma, Beta,
#'                Y, Z, delta, mu, mu_star, c[2], S,
#'                sigma_square, K, E, X,
#'                m, n, q, p, T0)
#'   for(i in unique(E)){
#'     ind<-sort(unique(R[i,]))
#'     KK<-length(ind)
#'     Gamma_temp<-Gamma[ind,,i]
#'     Gamma[,,i]<-NA
#'     Gamma[1:KK,,i]<-Gamma_temp
#'     w_temp<-w[i,ind]
#'     w_temp<-w_temp/sum(w_temp)
#'     w[i,]<-NA
#'     w[i,1:KK]<-w_temp
#'     for(k in 1:KK){
#'       R[i,which(R[i,]==ind[k])]<-k
#'     }
#'     K[i]<-KK
#'   }
#'   mu<-updatemu(R,Z,X,Gamma,K,Beta,E,m,n,p,q)
#'   mu_star<-updatemustar(mu,c,n,T0,D)
#'   step <- array(rnorm(max(K) * S *q, 0, nu_gamma),dim=c( max(K), q,S))
#'   run<- array(runif(max(K) * S *q, 0, 1),dim=c( max(K), q,S))
#'   A<-updateGamma(X,Y, Z, delta, Beta, Gamma, E, R, S, K , mu, mu_star, sigma_square, rep(c,T0),
#'                  step,  run,  n,  m, T0,  p,  q,  D,theta2,  tau)
#'   Gamma<-A$Gamma
#'   mu<-A$mu
#'   mu_star<-A$mustar
#'   for(i in 1:S){
#'     if(K[i]==1){
#'       Gammai = t(as.matrix(Gamma[1:K[i],,i]))
#'       C[1:K[i],1:K[i],i]<-updateC(Gammai,theta2,tau)
#'     }
#'     else{
#'       C[1:K[i],1:K[i],i]<-updateC(Gamma[1:K[i],,i],theta2,tau)
#'     }
#'   }
#'   A<-update_RJ(w, K, Gamma,Beta, E,
#'                Z, X, R, mu, mu_star, Y, delta, c,sigma_square, C,
#'                S, theta2, tau, q, m, T0)
#'   K<-A$K
#'   w<-A$w
#'   Gamma<-A$Gamma
#'   R<-A$R
#'   C<-A$C
#'   mu<-updatemu(R,Z,X,Gamma,K,Beta,E,m,n,p,q)
#'   mu_star<-updatemustar(mu,rep(c,T0),n,T0,D)
#'   C_beta<-updateC(Beta,theta1,tau)
#'   step<-matrix(rnorm(S*p,0,nu_beta),nrow = S)
#'   runif<-matrix(runif(S*p,0,1),nrow = S)
#'   A<- updateBeta( X,Y, Z, delta,
#'                   Beta, Gamma, E,R,S,K,mu_star, mu,sigma_square,rep(c,T0), C_beta, step,runif,
#'                   n, m, T0,  p,  q,  D,
#'                   theta1, tau)
#'   Beta<-A$Beta
#'   mu<-A$mu
#'   mu_star<-A$mustar
#'   record$Beta[iter,,]<-Beta

#'   A<-updateE( Beta,Gamma, w_beta, X,Y,Z,delta,E, R, S,K, mu, mu_star,sigma_square,rep(c,T0),
#'               n,  m,  T0,  p,  q,  D)
#'   E<-A$E
#'   mu<-A$mu
#'   mu_star<-A$mustar
#'   K<-A$Ds
#'   w_beta<- update_w_beta(S, w_beta, E)
#'   theta1<-update_theta_beta(theta1,tau,Beta)
#'   theta2<-update_theta_gamma(theta2,tau,Gamma,S,K)
#'   sigma_square <- update_sigma_squre(Y,mu)
#' }
#'
#' @export
update_RJ <- function(w, K, Gamma,Beta, E,
                      Z, X, R, mu, mu_star, Y, delta, c, sigma_square, C,
                      S, theta,tau, q,m, T0,
                      hyper_delta = 1){
  for(i in 1:S){
    wi<- w[i, 1:K[i]]
    Gammai<- Gamma[1:K[i],,i]
    Betai<- Beta[i,]
    indi<- which(E == i)
    Ri<-R[i,]
    Ki<-K[i]
    Ci<- C[1:Ki,1:Ki, i]
    if(Ki==1){
      Gammai = t(as.matrix(Gammai))
      Ci = as.matrix(Ci)
      wi = as.vector(wi)
    }
    if(length(indi) == 0){
      A<-RJi_empty(wi, Ki, Gammai,  Ri,  Ci,theta, tau, m,  q, hyper_delta)
    }
    else{
      Yi<-Y[indi,]
      deltai<-delta[indi,]
      mui<-mu[indi,]
      mu_stari<-mu_star[indi,]
      Xi<-X[indi,]
      if(length(indi)==1){
        Yi <- t(as.matrix(Yi))
        deltai <- t(as.matrix(deltai))
        mui <- t(as.matrix(mui))
        mu_stari <- t(as.matrix(mu_stari))
        Xi <- t(as.matrix(Xi))
      }
      A<-RJi(wi, Ki, as.matrix(Gammai), Betai, Xi, Yi, Z, Ri, deltai, mui[1,], mu_stari, 0.01, sigma_square, Ci,
             theta, tau, m, length(indi),q,T0, hyper_delta=1)
    }

    K[i]<-A$K
    w[i, 1:K[i]]<- A$w
    Gamma[,,i]<-NA
    Gamma[1:K[i],,i]<-A$Gamma
    R[i,]<-A$R
    C[1:K[i],1:K[i], i]<-A$C

  }
  return(list(K = K, w = w, Gamma = Gamma, R = R, C = C))
}


#' update_sigma_squre
#'
#' Update the noise variance
#'
#' @param Y The CAL observation matrix, with missing values
#' @param mu Current estimated mean matrix for CAL
#' @param a The hyper-parameter with default value being 1
#' @param b The hyper-parameter with default value being 1
#' @return Updated noise variance
#' @seealso \link{update_RJ} for a complete example for all functions in this package.
#' @export
update_sigma_squre<-function(Y, mu, a =1,b=1){
  temp <- na.omit(as.vector(Y-mu))
  s <- sum(temp^2)
  sigma2<-1/rgamma(1,a+length(temp)/2,b+s/2)
  return (sigma2)
}

#' const
#'
#' The function to compute the normalizing constant for Determinantal Point Process (DPP)
#'
#' @param theta The parameter that controls how repulsive the DPP is
#' @param tau The parameter we fixed at a large number
#' @return The computed normalizing constant
#' @export
const <- function(theta, tau){
  a = 1/(4*tau^2)
  b = 1/(theta^2)
  c = sqrt(a^2 + 2*a*b)
  sum = a + b + c
  lambda = rep(0, 100)
  for (i in 1:100){
    eig1 = sqrt(2*a/sum)*(b/sum)^(i-1)
    if (eig1 > 0.000001) {lambda[i] = eig1}
    else{
      break
    }
  }
  return(prod(lambda+1))
}

#' likeli_theta
#'
#' The function computes the likelihood of the hyperparameter
#'
#' @param theta The parameter that controls how repulsive the DPP is
#' @param tau The parameter we fixed at a large number
#' @param Beta The matrix that records the samples
#' @return The likelihood
#' @export
likeli_theta<-function(theta,tau, Beta){
  #function to minimize in order to update theta
  if(is.matrix(Beta)){
    C<-updateC(Beta,theta,tau)
    Likeli<-log(det(C))-log(const(theta,tau))
  }
  else{
    C<-kernelC(Beta,Beta,theta,tau)
    Likeli<-log(C)-log(const(theta,tau))
  }

  return(Likeli)
}

#' update_theta_beta
#'
#' Update the DPP hyper-parameter for patient level
#'
#' @param theta The DPP hyper-parameter for patient level
#' @param tau A fixed DPP hyper-parameter, which we suggest high value, say 10^5
#' @param Beta The linear coefficients for patient level covariates. Should be a matrix with S rows
#' @param sig The hyper-parameter with default value being 10
#' @return updated DPP hyper-parameter for patient level
#' @seealso \link{update_RJ} for a complete example for all functions in this package.
#' @export
update_theta_beta<-function(theta, tau, Beta,sig = 10){
  theta.new<- theta + (rnorm(1,0,0.1))
  ratio <- (likeli_theta(theta.new,tau,Beta) + dnorm(theta.new,0,sig,log=TRUE) -
              likeli_theta(theta,tau,Beta) - dnorm(theta,0,sig,log=TRUE))
  if(ratio>log(runif(1))){
    return(theta.new)
  }
  else{
    return(theta)
  }
}

#' update_theta_gamma
#'
#' Update the DPP hyper-parameter for site level
#'
#' @param theta The DPP hyper-parameter for patient level
#' @param tau A fixed DPP hyper-parameter, which we suggest high value, say 10^5
#' @param Gamma The linear coefficients for site level covariates. Should be a 3-dimensional array.
#' @param S The number of patient level clusters
#' @param Ds The number for site level clusters for each patient level cluster
#' @param sig The hyper-parameter with default value being 10
#' @return updated DPP hyper-parameter for site level
#' @seealso \link{update_RJ} for a complete example for all functions in this package.
#' @export
update_theta_gamma<-function(theta,tau,Gamma, S, Ds, sig = 10 ){
  theta.new<- theta + (rnorm(1,0,0.1))
  ratio <- 0
  for(i in 1:S){
    Gammatemp<-Gamma[1:Ds[i],,i]
    ratio <- ratio + (likeli_theta(theta.new,tau,Gammatemp) + dnorm(theta.new,0,sig,log=TRUE) -
                        likeli_theta(theta,tau,Gammatemp) - dnorm(theta,0,sig,log=TRUE))
  }
  if(ratio>log(runif(1))){
    return(theta.new)
  }
  else{
    return(theta)
  }
}

#tested






