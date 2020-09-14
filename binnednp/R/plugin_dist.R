#' Plug-in bandwidth selector for kernel distribution estimation and binned data.
#'
#' @param n Positive integer. Size of the complete sample.
#' @param y Vector. Observed values. They define the extremes of the sequence of intervals in which data is binned.
#' @param w Vector. Proportion of observations within each interval.
#' @param ni Vector. Number of observations within each interval.
#' @param gplugin Positive real number. Pilot bandwidth. If missing, rule-of-thumb bandwidth is considered.
#' @param confband Logical. If TRUE, bootstrap confidence bands for the distribution are constructed. Defaults to FALSE.
#' @param alpha Significance level for the bootstrap confidence bands. Defaults to 0.05.
#' @param B Number of bootstrap resamples. Defaults to 1000.
#' @param type Character. If \code{type} = "N", normality is assumed at the last step when calculating the plug-in bandwidth. If \code{type} = "A", parameter at last step is estimated nonparametrically using \code{gplugin} as bandwidth. Otherwise, the unknown parameter is estimated fitting a normal mixture. Defaults to \code{type} = "N".
#' @param plot Logical. If TRUE, results are plotted. Defaults to TRUE.
#' @param print Logical. If TRUE, script current status is printed. Defaults to TRUE.
#' @param model Character. Name of the parametric family of distributions to be fitted for the grouped sample. Parameters are estimated by maximum likelihood.
#' @param parallel Logical. If TRUE, confidence bands are estimated using parallel computing with sockets.
#' @param pars Environment. Needed for the well functioning of the script. DO NOT modify this argument.
#'
#' @return A list with components
#' \item{h}{Plug-in bandwidth.}
#' \item{Fh}{Function. Kernel distribution estimator with bandwidth \code{h}.}
#' \item{confband}{(Optional) Bootstrap confidence bands for the distribution function.}
#'
#' @examples
#' set.seed(1)
#' n <- 200 #complete sample size
#' k <- 30 #number of intervals
#' x <- rnorm(n,6,1) #complete sample
#' y <- seq(min(x)-0.2,max(x)+0.2,len=k+1) #intervals
#' w <- c(sapply(2:k,function(i)sum( x<y[i]&x>=y[i-1] )), sum(x<=y[k+1]&x>=y[k]) )/n #proportions
#' bw.dist.binned(n,y,w,plot=FALSE)
#'
#' @references
#' \insertRef{TesisMiguel2015}{binnednp}
#'
#' \insertRef{Phalaris2016}{binnednp}
#'
#' @useDynLib binnednp
#' @importFrom Rcpp sourceCpp
#'
#' @export
bw.dist.binned <- function(n,y,w,ni,gplugin,type="N",confband=F,B=1000,alpha=0.05,plot=TRUE,print=TRUE,model,parallel=FALSE,pars=new.env()){

  main <- !exists("k", where = pars, inherits = FALSE)

  if(main==TRUE){

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



    K <- function(x) stats::dnorm(x)
    pK <- Vectorize(function(x) stats::pnorm(x))

    C0 <- 2*stats::integrate(function(x)x*K(x)*pK(x),-Inf,Inf)$value

    l <- mean(diff(y))


    K2 <- function(x) (x^2-1)*stats::dnorm(x)
    K4 <- function(x) ((x^2-3)*(x^2-1)-2*x^2)*stats::dnorm(x)
    K6_0 <- kedd::kernel.fun(0,deriv.order=6,kernel="gaussian")$kx
    K8_0 <- kedd::kernel.fun(0,deriv.order=8,kernel="gaussian")$kx
    L1 <- function(w,h) -sum( sapply( 1:k,function(i) sum(K2((t[i]-t)/h)*w[i]*w) ) )/h^3 # Af1=-psi2
    L2 <- function(w,h) sum( sapply( 1:k,function(i) sum(K4((t[i]-t)/h)*w[i]*w) ) )/h^5 # Af2=psi4
    L3 <- function(w,h) -sum( sapply( 1:k,function(i) sum(kedd::kernel.fun((t[i]-t)/h,deriv.order=6,kernel="gaussian")$kx*w[i]*w) ) )/h^7
    L4 <- function(w,h) sum( sapply( 1:k,function(i) sum(kedd::kernel.fun((t[i]-t)/h,deriv.order=8,kernel="gaussian")$kx*w[i]*w) ) )/h^9
    L5 <- function(w,h) -sum( sapply( 1:k,function(i) sum(kedd::kernel.fun((t[i]-t)/h,deriv.order=10,kernel="gaussian")$kx*w[i]*w) ) )/h^11

  } else {
    t <- pars$t
    k <- pars$k
    y <- pars$y
    K <- pars$K
    pK <- pars$pK
    C0 <- pars$C0
    l <- pars$l
    K2 <- pars$K2
    K4 <- pars$K4
    K6_0 <- pars$K6_0
    K8_0 <- pars$K8_0
    L1 <- pars$L1
    L2 <- pars$L2
    L3 <- pars$L3
    L4 <- pars$L4
    L5 <- pars$L5
  }

  mu_hat <- sum(w*t)
  sigma_hat <- sqrt( sum(w*(t-mu_hat)^2) )

  if(type == "N" || type == "A"){

  if(missing(gplugin)) gplugin <- 1.59*sigma_hat*n^(-1/3)
  if(type == "A"){
    f_gplugin <- Vectorize(function(x)1/gplugin*sum(w*K((x-t)/gplugin)))
    llim <- mu_hat-3*sigma_hat
    rlim <- mu_hat+3*sigma_hat
    Af <- stats::integrate(function(x)f_gplugin(x)^2,llim,rlim)$value
    if(!exists("Af")) Af <- stats::integrate(function(x)f_gplugin(x)^2,y[1],y[k+1])$value
    psi10 <- -L5(w,gplugin)
  } else {
    Af <- 1/(2*sqrt(pi)*sigma_hat)
    psi10 <- -30240/((2*sigma_hat)^11*sqrt(pi))
  }
  eta8 <- (-2*K8_0*Af*l/psi10)^(1/11)
  psi8 <- L4(w,eta8)
  eta6 <- (-2*K6_0*Af*l/psi8)^(1/9)
  psi6 <- -L3(w,eta6)
  eta4 <- (-2*K4(0)*Af*l/psi6)^(1/7)
  psi4 <- L2(w,eta4)
  eta2 <- (-2*K2(0)*Af*l/psi4)^(1/5)
  Af1 <- L1(w,eta2)

  } else {

    clustering <- mclust::Mclust(rep(t,n*w),G=1:5,verbose=FALSE,modelNames="V")
    mix_params <- clustering$parameters
    mu_mixt <- mix_params$mean
    sigma_mixt <- sqrt(mix_params$variance$sigmasq)
    alfa_mixt <- mix_params$pro
    normal_mixt <- nor1mix::norMix(mu=mu_mixt,sigma=sigma_mixt,w=alfa_mixt)
    f1_mixt <- Vectorize( function(x) -sum(alfa_mixt*(x-mu_mixt)/sigma_mixt^3*stats::dnorm((x-mu_mixt)/sigma_mixt)) )
    mixlim1 <- nor1mix::qnorMix(0.001,normal_mixt)
    mixlim2 <- nor1mix::qnorMix(0.999,normal_mixt)
    Af1 <- stats::integrate(function(x)f1_mixt(x)^2,mixlim1,mixlim2)$value

  }

  h_AMISE <- ( C0/(n*Af1) )^(1/3)

  if(main == TRUE){
    x <- NULL
    Fh <- Vectorize( function(x)sum(w*stats::pnorm((x-t)/h_AMISE)) )
  }

  if(confband == T){


    pars$t <- t
    pars$k <- k
    pars$y <- y
    pars$K <- K
    pars$pK <- pK
    pars$C0 <- C0
    pars$l <- l
    pars$K2 <- K2
    pars$K4 <- K4
    pars$K6_0 <- K6_0
    pars$K8_0 <- K8_0
    pars$L1 <- L1
    pars$L2 <- L2
    pars$L3 <- L3
    pars$L4 <- L4
    pars$L5 <- L5


    gplugin <- bw.dist.binned.boot(n,y,w,pilot.type=2,print=FALSE)$h

    Fhy <- sapply(y, function(x)sum( w*stats::pnorm((x-t)/h_AMISE) ))

    Fg <- Vectorize(function(x)sum( w*stats::pnorm((x-t)/gplugin) ))
    Fgy <- Fg(y)

    Dn <- matrix(0,nrow=length(y),ncol=B)

    if(!parallel){

    for(b in 1:B){
      rx <- sort(sample(t,replace=T,size=n,prob=w)+gplugin*stats::rnorm(n))
      wb <- calcw_cpp(rx,y)
      h_AMISE_boot <- bw.dist.binned(n,y,wb,plot=FALSE,print=FALSE,type=type,pars=pars)
      Fhby <- sapply(y,function(x)sum( wb*stats::pnorm( (x-t)/h_AMISE_boot ) ))
      Dn[,b] <- Fhby
      if(print == TRUE) cat("\rConstructing bootstrap confidence bands. Progress:",floor(100*b/B),"%")
    }

    } else {

      parfun <- function(b){
        rx <- sort(sample(t,replace=T,size=n,prob=w)+gplugin*stats::rnorm(n))
        wb <- calcw_cpp(rx,y)
        h_AMISE_boot <- bw.dist.binned(n,y,wb,plot=FALSE,print=FALSE,type=type,pars=pars)
        Fhby <- sapply(y,function(x)sum( wb*stats::pnorm( (x-t)/h_AMISE_boot ) ))
        return(Fhby)
      }

      ncores <- parallel::detectCores()
      cl <- parallel::makeCluster(ncores)
      parallel::clusterEvalQ(cl, library(binnednp))
      parallel::clusterExport(cl, 'calcw_cpp')
      paroutput <- parallel::parSapply(cl, 1:B, parfun)
      Dn <- paroutput
      parallel::stopCluster(cl)

    }

    alfa <- alpha/length(y)

    count <- 0
    maxcount <- 100
    low.alpha <- alfa
    high.alpha <- alpha
    if(print == TRUE) cat("\n")
    while(count < maxcount){
      if(print == TRUE) cat("\rAdjusting significance level. Progress:",floor(100*(count+1)/maxcount),"%")

      mean.alpha <- 0.5*(low.alpha+high.alpha)

      q1 <- apply(Dn,1,function(i)stats::quantile(i,low.alpha/2))
      q2 <- apply(Dn,1,function(i)stats::quantile(i,1-low.alpha/2))
      banda1 <- pmin(pmax(Fhy+Fgy-q1,0),1)
      banda2 <- pmin(pmax(Fhy+Fgy-q2,0),1)
      p_low.alpha <- sapply(1:B,function(b){
        aux <- ifelse(Dn[,b] < banda1 && Dn[,b] > banda2,1,0)
        return(aux)
      })
      p_low.alpha <- sum(p_low.alpha)/B

      q1 <- apply(Dn,1,function(i)stats::quantile(i,high.alpha/2))
      q2 <- apply(Dn,1,function(i)stats::quantile(i,1-high.alpha/2))
      banda1 <- pmin(pmax(Fhy+Fgy-q1,0),1)
      banda2 <- pmin(pmax(Fhy+Fgy-q2,0),1)
      p_high.alpha <- sapply(1:B,function(b){
        aux <- ifelse(Dn[,b] < banda1 && Dn[,b] > banda2,1,0)
        return(aux)
      })
      p_high.alpha <- sum(p_high.alpha)/B

      q1 <- apply(Dn,1,function(i)stats::quantile(i,mean.alpha/2))
      q2 <- apply(Dn,1,function(i)stats::quantile(i,1-mean.alpha/2))
      banda1 <- pmin(pmax(Fhy+Fgy-q1,0),1)
      banda2 <- pmin(pmax(Fhy+Fgy-q2,0),1)
      p_mean.alpha <- sapply(1:B,function(b){
        aux <- ifelse(Dn[,b] < banda1 && Dn[,b] > banda2,1,0)
        return(aux)
      })
      p_mean.alpha <- sum(p_mean.alpha)/B

      if(p_mean.alpha >= 1-alpha){
        low.alpha <- mean.alpha
      } else {
        high.alpha <- mean.alpha
      }

      count <- count+1
    }


    q1 <- apply(Dn,1,function(i)stats::quantile(i,mean.alpha/2))
    q2 <- apply(Dn,1,function(i)stats::quantile(i,1-mean.alpha/2))

    banda1 <- pmin(pmax(Fhy+Fgy-q1,0),1)
    banda2 <- pmin(pmax(Fhy+Fgy-q2,0),1)

    confband <- cbind(banda1,banda2)


    if(plot == TRUE){

      gu <- Vectorize( function(t){
        j <- max(which(y<t))
        a <- y[j]
        d <- confband[j,1]
        fa <- Fh(a)
        return(Fh(t)+d-fa)
      } )

      hu <- Vectorize( function(t){
        j <- max(which(y<t))
        a <- y[j]
        b <- y[j+1]
        f <- confband[j+1,1]
        return(gu(t)+(t-a)/(b-a)*(f-gu(b)))
      } )

      gl <- Vectorize( function(t){
        j <- max(which(y<t))
        a <- y[j]
        d <- confband[j,2]
        fa <- Fh(a)
        return(Fh(t)+d-fa)
      } )

      hl <- Vectorize( function(t){
        j <- max(which(y<t))
        a <- y[j]
        b <- y[j+1]
        f <- confband[j+1,2]
        return(gl(t)+(t-a)/(b-a)*(f-gl(b)))
      } )

      grid <- seq(y[1],y[k+1],len=500)[-c(1,500)]

      grDevices::dev.new(noRStudioGD=TRUE)
      graphics::curve(Fh,min(y),max(y),type="n",lwd=2,ylab="Cumulative probability",main="Kernel distribution")
      graphics::polygon(c(grid,rev(grid)),c(hl(grid),rev(hu(grid))),col="grey",border=NA)
      graphics::curve(Fh,lwd=2,add=T)
      graphics::lines(grid,hu(grid),col=2,lty=2)
      graphics::lines(grid,hl(grid),col=2,lty=2)
      if(!missing(model)){
        pmodel <- paste0("p",model)
        try(parfit <- fitdistrplus::fitdist(rep(t,ceiling(n*w)),model),silent=T)
        if(exists("parfit")){
          distr <- get(pmodel)
          params <- as.list(parfit$estimate)
          pargrid <- seq(min(y),max(y),len=101)
          clist <- c(list(q=pargrid),params)
          graphics::lines(pargrid, do.call( distr, clist ), col=3, lwd=2 )

        }
      }
    }

    return(list(h=h_AMISE,confband=confband,Fh=Fh))

  } else {

    if(plot == TRUE){
      grDevices::dev.new(noRStudioGD=TRUE)
      graphics::curve(Fh,min(y),max(y),lwd=2,ylab="Cumulative probability",main="Kernel distribution")
      if(!missing(model)){
        pmodel <- paste0("p",model)
        try(parfit <- fitdistrplus::fitdist(rep(t,ceiling(n*w)),model),silent=T)
        if(exists("parfit")){
          distr <- get(pmodel)
          params <- as.list(parfit$estimate)
          pargrid <- seq(min(y),max(y),len=101)
          clist <- c(list(q=pargrid),params)
          graphics::lines(pargrid, do.call( distr, clist ), col=3, lwd=2 )

        }
      }
    }

    if(main == TRUE){
      return(list(h=h_AMISE, Fh=Fh))
    } else {
      return(h=h_AMISE)
    }

  }

}
