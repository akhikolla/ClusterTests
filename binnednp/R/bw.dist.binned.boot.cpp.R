#' Bootstrap bandwidth selector for kernel distribution estimation and binned data.
#'
#' @param n Positive integer. Size of the complete sample.
#' @param y Vector. Observed values. They define the extremes of the sequence of intervals in which data is binned.
#' @param w Vector. Proportion of observations within each interval.
#' @param ni Vector. Number of observations within each interval.
#' @param g Positive real number. Pilot bandwidth. If missing, plug-in N bandwidth for the distribution is considered.
#' @param pilot.type 1 or 2. If \code{g} is missing, pilot bandwidth for the bootstrap bandwidth selector is automatically selected using methods 1 or 2. Defaults to 1. See details for more information.
#' @param nit Positive integer. Number of iterations in the dichotomy algorithm for the estimation of the bootstrap bandwidth.
#' @param confband Logical. If TRUE, bootstrap confidence bands are constructed for the estimator. Defaults to FALSE.
#' @param B Positive integer. Number of bootstrap resamples used for the construction of the confidence bands. Defaults to 1000.
#' @param alpha Real number between 0 and 1. Significance level for the confidence bands. Defaults to 0.05
#' @param print Logical. If TRUE, script current status is printed. Defaults to TRUE.
#' @param plot Logical. If TRUE, results are plotted. Defaults to FALSE.
#' @param parallel Logical. If TRUE, confidence bands are estimated using parallel computing with sockets.
#' @param pars Environment. Needed for the well functioning of the script. DO NOT modify this argument.
#'
#' @details
#' If \code{pilot.type} = 1, plug-in bandwidth for the distribution is considered as pilot bandwidth for the bootstrap selector.
#'
#' If \code{pilot.type} = 2, the pilot bandwidth is such that the kernel distribution estimator with bandwidth \code{g} approximates the empirical distribution of the grouped sample minimizing the residual sum of squares. Also, a penalty is imposed on the global slope of the kernel density estimator with bandwidth \code{g}. The penalty parameter is selected as to best approximate the global slope of the true density.
#'
#' @return A list with components:
#' \item{h}{Bootstrap bandwidth for the distribution function.}
#' \item{Fh}{Function. Kernel distribution estimator with bandwidth \code{h}.}
#' \item{confband (optional)}{Matrix. Its columns contain the bootstrap confidence bands for the estimator.}
#'
#' @examples
#' set.seed(1)
#' n <- 200 #complete sample size
#' k <- 30 #number of intervals
#' x <- rnorm(n,6,1) #complete sample
#' y <- seq(min(x)-0.2,max(x)+0.2,len=k+1) #intervals
#' w <- c(sapply(2:k,function(i)sum( x<y[i]&x>=y[i-1] )), sum(x<=y[k+1]&x>=y[k]) )/n #proportions
#' bw.dist.binned.boot(n,y,w,plot=FALSE)
#'
#' @references
#' \insertRef{TesisMiguel2015}{binnednp}
#'
#' @useDynLib binnednp
#' @importFrom Rcpp sourceCpp
#'
#' @export
bw.dist.binned.boot <- function(n,y,w,ni,g,pilot.type=2,nit=10,confband=FALSE,B=1000,alpha=0.05,print=TRUE,plot=TRUE,parallel=FALSE,pars=new.env()){

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

  if(missing(g)){

    if(pilot.type == 1) g <- bw.dist.binned(n,y,w,plot=F,print=F)$h

    if(pilot.type == 2){

      clustering <- mclust::Mclust(rep(t,n*w),G=1:5,verbose=FALSE,modelNames="V")
      mix_params <- clustering$parameters
      mu_mixt <- mix_params$mean
      sigma_mixt <- sqrt(mix_params$variance$sigmasq)
      alfa_mixt <- mix_params$pro
      normal_mixt <- nor1mix::norMix(mu=mu_mixt,sigma=sigma_mixt,w=alfa_mixt)
      f1_mixt <- Vectorize( function(x) -sum(alfa_mixt*(x-mu_mixt)/sigma_mixt^3*stats::dnorm((x-mu_mixt)/sigma_mixt)) )
      mixlim1 <- nor1mix::qnorMix(0.001,normal_mixt)
      mixlim2 <- nor1mix::qnorMix(0.999,normal_mixt)
      Af1_mixt <- stats::integrate(function(x)f1_mixt(x)^2,mixlim1,mixlim2)$value

      Fht <- function(h)sapply(y,function(x)sum(w*stats::pnorm((x-t)/h)))
      emp <- c(0,cumsum(w))
      h0 <- min(diff(y))/4
      h1 <- max(y)-min(y)
      rho <- exp((log(h1)-log(h0))/4)
      l0 <- 1e-3
      l1 <- 10
      lrho <- exp((log(l1)-log(l0))/4)

      g <- zeta_hist_p_dist(emp,t,y,w,Af1_mixt,l0,l1,h0,h1,lrho,rho,10,10,mixlim1,mixlim2)

    }

  }

  x <- NULL
  Fg <- Vectorize( function(x)sum(w*stats::pnorm((x-t)/g)) )

  p <- Fg(y[2:(k+1)])-Fg(y[1:k])

  mu <- sum(w*t)
  sgm <- sqrt(sum(w*(t-mu)^2))
  lim1 <- mu-3*sgm
  lim2 <- mu+3*sgm

  if(main == TRUE){

  h0 <- min(diff(unique(y)))/2
  h1 <- max(y)-min(y)
  rho <- exp( (log(h1)-log(h0))/4 )

  } else {

    h0 <- pars$h0
    h1 <- pars$h1
    rho <- pars$rho

  }

  hboot <- boot_bw_dist(nit,h0,h1,rho,n,t,w,p,g,lgrid=100,lim1,lim2)

  Fh <- Vectorize(function(x)sum(w*stats::pnorm((x-t)/hboot)))


  if(confband == TRUE){

    pars$t <- t
    pars$k <- k
    pars$y <- y
    pars$h0 <- h0
    pars$h1 <- h1
    pars$rho <- rho

    Dn <- matrix(0,nrow=(k+1),ncol=B)
if(!parallel){
    for(b in 1:B){
      rx <- sort(sample(t,size=n,replace=TRUE,prob=w)+hboot*stats::rnorm(n))
      wb <- calcw_cpp(rx,y)
      Fhboot_wb <- bw.dist.binned.boot(n,y,wb,plot=FALSE,print=FALSE,pilot.type=pilot.type,pars=pars)$Fh
      Dn[,b] <- Fhboot_wb(y)
      if(print == TRUE) cat("\rConstructing bootstrap confidence bands. Progress:",floor(100*b/B),"%")
    }
} else {
    parfun <- function(b)
    {
      rx <- sort(sample(t,size=n,replace=TRUE,prob=w)+hboot*stats::rnorm(n))
      wb <- calcw_cpp(rx,y)
      Fhboot_wb <- bw.dist.binned.boot(n,y,wb,plot=FALSE,print=FALSE,pilot.type=pilot.type,pars=pars)$Fh
      return(Fhboot_wb(y))
    }
    ncores <- parallel::detectCores()
    cl <- parallel::makeCluster(ncores)
    parallel::clusterEvalQ(cl, library(binnednp))
    parallel::clusterExport(cl, 'calcw_cpp')
    paroutput <- parallel::parSapply(cl, 1:B, parfun)
    Dn <- paroutput
    parallel::stopCluster(cl)
}


    alfa <- alpha/(k+1)

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
      banda1 <- pmin(pmax(2*Fh(y)-q1,0),1)
      banda2 <- pmin(pmax(2*Fh(y)-q2,0),1)
      p_low.alpha <- sapply(1:B,function(b){
        aux <- ifelse(Dn[,b] < banda1 && Dn[,b] > banda2,1,0)
        return(aux)
      })
      p_low.alpha <- sum(p_low.alpha)/B

      q1 <- apply(Dn,1,function(i)stats::quantile(i,high.alpha/2))
      q2 <- apply(Dn,1,function(i)stats::quantile(i,1-high.alpha/2))
      banda1 <- pmin(pmax(2*Fh(y)-q1,0),1)
      banda2 <- pmin(pmax(2*Fh(y)-q2,0),1)
      p_high.alpha <- sapply(1:B,function(b){
        aux <- ifelse(Dn[,b] < banda1 && Dn[,b] > banda2,1,0)
        return(aux)
      })
      p_high.alpha <- sum(p_high.alpha)/B

      q1 <- apply(Dn,1,function(i)stats::quantile(i,mean.alpha/2))
      q2 <- apply(Dn,1,function(i)stats::quantile(i,1-mean.alpha/2))
      banda1 <- pmin(pmax(2*Fh(y)-q1,0),1)
      banda2 <- pmin(pmax(2*Fh(y)-q2,0),1)
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
    band1 <- pmin(pmax(2*Fh(y)-q1,0),1)
    band2 <- pmin(pmax(2*Fh(y)-q2,0),1)
    confband <- cbind(band1,band2)

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
      grDevices::dev.new(noRStudioGD = TRUE)
      graphics::curve(Fh,min(y),max(y),type="n")
      graphics::polygon(c(grid,rev(grid)),c(hl(grid),rev(hu(grid))),col="grey",border=NA)
      graphics::curve(Fh,lwd=2,add=T)
      graphics::lines(grid,hu(grid),col=2,lty=2)
      graphics::lines(grid,hl(grid),col=2,lty=2)


    }

    return(list(h=hboot,Fh=Fh,confband=confband))

  } else {

    return(list(h=hboot,Fh=Fh))

  }

}
