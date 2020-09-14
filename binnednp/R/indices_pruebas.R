#' Nonparametric estimates of indices measuring the global slope and curvature of the density function for binned data.
#'
#' @param n Positive integer. Size of the complete sample.
#' @param y Vector. Observed values. They define the extremes of the sequence of intervals in which data is binned.
#' @param w Vector. Proportion of observations within each interval.
#' @param ni Vector. Number of observations within each interval.
#' @param hseq Vector. Grid of bandwidths for which MSE is approximated through bootstrap. If missing, a grid of size hn is considered.
#' @param hn Positive integer. Size of the grid of bandwidths for which MSE will be approximated. Defaults to 200.
#' @param nmix Positive integer. Maximum number of components for the normal mixture model. Defaults to 4.
#' @param B Positive integer. Number of bootstrap resamples used to find the bandwidth that minimizes MSE. Defaults to 500.
#' @param method Character. If method="np", resamples are taken from kernel density estimator with pilot bandwidth. If method="mix", a normal mixture pilot model is considered. If method="plugin", plug-in estimates are returned. Defaults to "np".
#' @param confint Logical. If TRUE, bootstrap confidence intervals are constructed for the indices. Defaults to FALSE.
#' @param B.conf Positive integer. Number of bootstrap resamples considered to construct the confidence intervals.
#' @param alpha Real number between 0 and 1. Significance level considered for the confidence intervals.
#' @param print Logical. If TRUE, status of the script and results are printed. Defaults to TRUE.
#' @param last.iter.np Logical. If FALSE, normality is assumed at the last step when calculating the plug-in bandwidth. Otherwise, rule-of-thumb selector is used. Defaults to FALSE.
#' @param parallel Logical. If TRUE, confidence bands are estimated using parallel computing with sockets.
#' @param pars Environment. Needed for the well functioning of the script. DO NOT modify this argument.
#'
#' @return Nonparametric estimates of the indices and (optional) confidence intervals.
#'
#' @examples
#' set.seed(1)
#' n <- 200 #complete sample size
#' k <- 30 #number of intervals
#' x <- rnorm(n,6,1) #complete sample
#' y <- seq(min(x)-0.2,max(x)+0.2,len=k+1) #intervals
#' w <- c(sapply(2:k,function(i)sum( x<y[i]&x>=y[i-1] )), sum(x<=y[k+1]&x>=y[k]) )/n #proportions
#' emergence.indices(n,y,w)
#'
#' @references
#' \insertRef{Indices2011}{binnednp}
#'
#' \insertRef{Test2017}{binnednp}
#'
#' @useDynLib binnednp
#' @importFrom Rcpp sourceCpp
#'
#' @export
emergence.indices = function(n, y, w, ni, hseq, hn=200, nmix=4, B=500, method="np", last.iter.np=F, confint=FALSE, B.conf=1000, alpha=0.05, print=TRUE, parallel=FALSE, pars=new.env()) {

  main <- !exists("k", where = pars, inherits = FALSE)

  if(main == TRUE){
    t <- y[-length(y)]+diff(y)/2
    k <- length(t)

    if(missing(w)){

      if(!missing(ni)){
        if(missing(n)) n <- sum(ni)
        w <- ni/n
      } else {
        stop("Arguments w or ni must be provided.")
      }

    }

    if(anyDuplicated(y) != 0){
      dup <- which(duplicated(y))
      comby <- y[-dup]
      lcy <- length(comby)
      combw <- numeric(lcy-1)
      for(i in 2:lcy){
        combw[i-1] <- sum(w[which(y==comby[i])-1])
      }
      combt <- comby[-lcy]+diff(comby)/2

      t <- combt
      y <- comby
      w <- combw
      k <- lcy-1
    }

  } else {
    t <- pars$t
    k <- pars$k
    y <- pars$y
  }

  mu_hat <- sum(w*t)
  sigma_hat <- sqrt( sum(w*(t-mu_hat)^2) )

  if(main == TRUE){
    I1 <- sigma_hat/mu_hat
    I2 <- sum(w*(t-mu_hat)^4)/sigma_hat^4

    l <- mean(diff(y))
    r <- max(y)-min(y)


    K2 <- function(x) (x^2-1)*stats::dnorm(x)

    K4 <- function(x) ((x^2-3)*(x^2-1)-2*x^2)*stats::dnorm(x)




    L1 <- function(w,h) -sum( sapply( 1:k,function(i) sum(K2((t[i]-t)/h)*w[i]*w) ) )/h^3

    L2 <- function(w,h) sum( sapply( 1:k,function(i) sum(K4((t[i]-t)/h)*w[i]*w) ) )/h^5

    L3 <- function(w,h) -sum( sapply( 1:k,function(i) sum(kedd::kernel.fun((t[i]-t)/h,deriv.order=6,kernel="gaussian")$kx*w[i]*w) ) )/h^7

    L4 <- function(w,h) sum( sapply( 1:k,function(i) sum(kedd::kernel.fun((t[i]-t)/h,deriv.order=8,kernel="gaussian")$kx*w[i]*w) ) )/h^9

    L5 <- function(w,h) -sum( sapply( 1:k,function(i) sum(kedd::kernel.fun((t[i]-t)/h,deriv.order=10,kernel="gaussian")$kx*w[i]*w) ) )/h^11



    if(method=="plugin") {



      K <- function(x)stats::dnorm(x)
      mu2K <- stats::integrate(function(x)x^2*K(x),-Inf,Inf)$value
      hAf <- 1.06*sigma_hat*n^(-1/5)
      f_hAf <- Vectorize(function(x)1/hAf*sum(w*K((x-t)/hAf)))
      try(Af <- stats::integrate(function(x)f_hAf(x)^2,-Inf,Inf)$value, silent=T)
      if(!exists("Af")) Af <- stats::integrate(function(x)f_hAf(x)^2,y[1],y[k+1])$value

      if(last.iter.np==T){
        eta5 <- hAf
        Af5 <- -L5(w,eta5)
      } else {
        Af <- 1/(2*sqrt(pi)*sigma_hat)
        Af5 <- -30240/((2*sigma_hat)^11*sqrt(pi))
      }

      K8_0 <- kedd::kernel.fun(0,deriv.order=8,kernel="gaussian")$kx

      eta4 <- (-2*K8_0*Af*l/(mu2K*Af5))^(1/11)
      Af4 <- L4(w,eta4)
      K6_0 <- kedd::kernel.fun(0,deriv.order=6,kernel="gaussian")$kx

      eta3 <- (-2*K6_0*Af*l/(mu2K*Af4))^(1/9)
      Af3 <- -L3(w,eta3)

      eta2 <- (-2*K4(0)*Af*l/(mu2K*Af3))^(1/7)
      Af2 <- L2(w,eta2)

      eta1 <- (-2*K2(0)*Af*l/(mu2K*Af2))^(1/5)
      Af1 <- L1(w,eta1)

      J1_plugin <- sigma_hat^3*Af1
      J2_plugin <- sigma_hat^5*Af2



    }

  } else {
    I1 <- pars$I1
    I2 <- pars$I2
    l <- pars$l
    r <- pars$r
    K2 <- pars$K2
    K4 <- pars$K4
    L1 <- pars$L1
    L2 <- pars$L2
    L3 <- pars$L3
    L4 <- pars$L4
    L5 <- pars$L5
  }





  if(method=="np") {

    if(main == TRUE & print == TRUE) cat("Searching for bandwidth minimizing MSE...")




    h_boot <- bw.dens.binned(n,y,w,hn=100,pilot.type=3,plot=FALSE,print=FALSE)$h_boot

    J1_np <- sigma_hat^3*L1(w,h_boot)
    J2_np <- sigma_hat^5*L2(w,h_boot)




    if(main == TRUE){

      if(missing(hseq)){
        hseq <- seq(min(diff(t))/2,max(t)-min(t),len=hn)
      }
      else{
        hn <- length(hseq)
      }

    } else {
      hseq <- pars$hseq
    }

    MSE_J1 <- numeric(hn)
    MSE_J2 <- numeric(hn)

    J1b <- numeric(B)
    J2b <- numeric(B)


    xb <- matrix(0,nrow=n,ncol=B)
    for(b in 1:B){

      xb[,b] <- sample(t,size=n,replace=TRUE,prob=w) + h_boot*stats::rnorm(n)

    }

    main_method_np(hn, B, hseq, xb, y, t, MSE_J1, MSE_J2, J1_np, J2_np, J1b, J2b)

    sd1 <- stats::sd(log(J1b))
    sd2 <- stats::sd(log(J2b))

    if(main == TRUE & print == TRUE) cat(" Done!\n")

  }





  if(method=="mix") {




    # fits <- list()
    # rept <- rep(t,n*w)
    # for(r in 1:nmix){
    #
    #   fits[[r]] <- nor1mix::norMixEM(rept,m=stats::kmeans(rept,centers=r)$cluster,maxiter=10^4,trace=FALSE)
    #
    # }
    #
    #
    # loglik <- unlist(lapply(fits,function(i)attr(i,"loglik")))
    #
    #
    #
    # test <- 2*diff(loglik)
    # nparams <- 3*1:nmix-1
    # df <- diff(nparams)
    # pvalue <- stats::pchisq(test,df,lower.tail=FALSE)
    # alpha <- 0.05
    #
    # j <- 1
    # while(j < nmix){
    #   if(pvalue[j] > alpha) break
    #   else j <- j+1
    # }
    #
    #
    #
    #
    # mixfit <- fits[[j]]
    # nmix <- j
    # mumix <- mixfit[,1]
    # sigmix <- mixfit[,2]
    # alfamix <- mixfit[,3]


    clustering <- mclust::Mclust(rep(t,n*w),G=1:5,verbose=FALSE,modelNames="V")
    mix_params <- clustering$parameters
    mumix <- mix_params$mean
    sigmix <- sqrt(mix_params$variance$sigmasq)
    alfamix <- mix_params$pro
    mixfit <- nor1mix::norMix(mu=mumix,sigma=sigmix,w=alfamix)




    fmix <- Vectorize( function(x) sum(alfamix/sigmix*stats::dnorm((x-mumix)/sigmix)) )

    fmix1 <- Vectorize( function(x) -sum(alfamix*(x-mumix)/sigmix^3*stats::dnorm((x-mumix)/sigmix)) )

    fmix2 <- Vectorize( function(x) -sum( alfamix*stats::dnorm((x-mumix)/sigmix)*( 1/sigmix^3-(x-mumix)^2/sigmix^5 ) ) )


    mu_pooled <- sum(alfamix*mumix)
    sigma_pooled <- sqrt(sum(alfamix*(sigmix^2+(mumix-mu_pooled)^2)))
    J1_mix <- sigma_pooled^3*stats::integrate(function(x)fmix1(x)^2,-Inf,Inf)$value
    J2_mix <- sigma_pooled^5*stats::integrate(function(x)fmix2(x)^2,-Inf,Inf)$value




    if(main == TRUE & print == TRUE) cat("Searching for bandwidth minimizing MSE...")


    if(main == TRUE){

      if(missing(hseq)){
        hseq <- seq(min(diff(t))/2,max(t)-min(t),len=hn)
      }
      else{
        hn <- length(hseq)
      }

    } else {
      hseq <- pars$hseq
    }

    MSE_J1 <- numeric(hn)
    MSE_J2 <- numeric(hn)

    J1b <- numeric(B)
    J2b <- numeric(B)

    xb <- matrix(0,nrow=n,ncol=B)
    for(b in 1:B){

      # ib <- sample(1:nmix, size=n, replace=T, prob=alfamix)
      # xb[,b] <- stats::rnorm(n,mumix[ib],sigmix[ib])
      xb[,b] <- nor1mix::rnorMix(n,mixfit)

    }

    J1_np <- J1_mix
    J2_np <- J2_mix

    main_method_np(hn, B, hseq, xb, y, t, MSE_J1, MSE_J2, J1_np, J2_np, J1b, J2b)

    sd1 <- stats::sd(log(J1b))
    sd2 <- stats::sd(log(J2b))

    if(main == TRUE & print == TRUE) cat(" Done!\n")



  }




  if(method=="np" | method=="mix") {

    h_J1 <- hseq[which.min(MSE_J1)]
    h_J2 <- hseq[which.min(MSE_J2)]


    J1_opt <- sigma_hat^3*L1(w,h_J1)
    J2_opt <- sigma_hat^5*L2(w,h_J2)




    if(confint == TRUE){

      if(print == TRUE) cat("Constructing bootstrap confidence intervals...\n")


      hseq.res <- hseq


      pars$t <- t
      pars$k <- k
      pars$y <- y
      pars$I1 <- I1
      pars$I2 <- I2
      pars$l <- l
      pars$r <- r
      pars$K2 <- K2
      pars$K4 <- K4
      pars$L1 <- L1
      pars$L2 <- L2
      pars$L3 <- L3
      pars$L4 <- L4
      pars$L5 <- L5
      pars$hseq <- hseq.res



      est1 <- numeric(B.conf)
      est2 <- numeric(B.conf)

      if(!parallel){

      for(b in 1:B.conf){

        if(method=="np"){
          res <- sample(t,size=n,replace=TRUE,prob=w) + h_boot * stats::rnorm(n)
        }
        else{
          # ib <- sample(1:nmix,size=n,replace=T,prob=alfamix)
          # res <- stats::rnorm(n,mumix[ib],sigmix[ib])
          res <- nor1mix::rnorMix(n,mixfit)
        }


        wb <- calcw_cpp(res,y)

        obj <- emergence.indices(n,y,wb,B=100,hseq=hseq.res,confint=FALSE,pars=pars)



        est1[b] <- log(obj$J1/J1_np)/obj$sd1
        est2[b] <- log(obj$J2/J2_np)/obj$sd2


        if(print == TRUE) cat("\rProgress:",floor(100*b/B.conf),"%")

      }

      } else {

        parfun <- function(b){
          if(method=="np"){
            res <- sample(t,size=n,replace=TRUE,prob=w) + h_boot * stats::rnorm(n)
          }
          else{
            # ib <- sample(1:nmix,size=n,replace=T,prob=alfamix)
            # res <- stats::rnorm(n,mumix[ib],sigmix[ib])
            res <- nor1mix::rnorMix(n,mixfit)
          }


          wb <- calcw_cpp(res,y)

          obj <- emergence.indices(n,y,wb,B=100,hseq=hseq.res,confint=FALSE,pars=pars)

          return(c(log(obj$J1/J1_np)/obj$sd1,log(obj$J2/J2_np)/obj$sd2))

        }

        ncores <- parallel::detectCores()
        cl <- parallel::makeCluster(ncores)
        parallel::clusterEvalQ(cl, {
          library(binnednp)
          library(nor1mix)
          })
        parallel::clusterExport(cl, 'calcw_cpp')
        paroutput <- parallel::parSapply(cl, 1:B.conf, parfun)
        est1 <- paroutput[1,]
        est2 <- paroutput[2,]
        parallel::stopCluster(cl)

      }

      if(print == TRUE) cat("\nDone!\n")



      bias_J1_opt <- mean(est1)-ifelse(method=="np",J1_np,J1_mix)
      bias_J2_opt <- mean(est2)-ifelse(method=="np",J2_np,J2_mix)

      var_J1_opt <- stats::var(est1)
      var_J2_opt <- stats::var(est2)

      MSE_J1_opt <- var_J1_opt+bias_J1_opt^2
      MSE_J2_opt <- var_J2_opt+bias_J2_opt^2


      int1_basic <- exp( log(J1_opt)-stats::quantile(est1,c(1-alpha/2,alpha/2))*sd1 )
      int2_basic <- exp( log(J2_opt)-stats::quantile(est2,c(1-alpha/2,alpha/2))*sd2 )




    }



  }

  if(method == "plugin") return(list(I1=I1,I2=I2,J1=J1_plugin,J2=J2_plugin))

  if(confint == FALSE) return(list(I1=I1,I2=I2,h_J1=h_J1,h_J2=h_J2,J1=J1_opt,J2=J2_opt,
                                   sd1=sd1,sd2=sd2))
  if(confint == TRUE) return(list(I1=I1,I2=I2,h_J1=h_J1,h_J2=h_J2,J1=J1_opt,J2=J2_opt,
                                  int1=int1_basic,int2=int2_basic,est1=est1,est2=est2,
                                  J1_np=J1_np,J2_np=J2_np,sd1=sd1,sd2=sd2))




}
