#' ANOVA in kernel distribution estimation with binned data using bootstrap.
#' 
#' @param n Vector of positive integers. Sizes of the complete samples corresponding to each treatment.
#' @param y Vector. Observed values. They define the extremes of the sequence of intervals in which data is binned.
#' @param trt.w Matrix. Proportion of observations within each interval. Each column corresponds to a different treatment.
#' @param abs.values Logical. Indicates if the values of trt.w are given in absolute (TRUE) or relative (FALSE) format.
#' @param B Positive integer. Number of bootstrap replicates used to compute the confidence bands.
#' 
#' @details 
#' ANOVA for interval-grouped data.
#' 
#' @return p-value of the test.
#' 
#' @references 
#' \insertRef{TesisMiguel2015}{binnednp}
#' 
#' @useDynLib binnednp
#' @importFrom Rcpp sourceCpp
#' 
#' @export
anv.binned <- function(n,y,trt.w,abs.values=FALSE,B=500)
{
  
  ntrt <- ncol(trt.w)
  lwp <- nrow(trt.w)
  npooled <- sum(n)
  if(abs.values)
  {
    for(i in 1:ntrt)
    {
      trt.w[,i] <- trt.w[,i]/n[i]
    }
  }
  
  wpooled <- numeric(lwp)
  for(i in 1:lwp)
  {
    wpooled[i] <- sum(trt.w[i,]*n)
  }
  wpooled <- wpooled/npooled
  
  t <- y[-length(y)]+diff(y)/2
  lim2 <- t[length(t)]
  lim1 <- t[1]
  obj <- bw.dist.binned(npooled,y,wpooled,print=FALSE,plot=FALSE)
  hpooled <- obj$h
  Fh <- obj$Fh
  D <- 0
  fh <- Vectorize(function(x)1/hpooled*sum(wpooled*stats::dnorm((x-t)/hpooled)))
  for(i in 1:ntrt)
  {
    Fhi <- bw.dist.binned(n[i],y,trt.w[,i],plot=FALSE,print=FALSE)$Fh
    D <- D+stats::integrate(function(x)(Fhi(x)-Fh(x))^2*fh(x),lim1+hpooled*stats::qnorm(0.99),lim2+hpooled*stats::qnorm(0.01))$value*n[i]
  }
  D <- D/ntrt
  
  
  Db <- numeric(B)
  for(b in 1:B)
  {
    cat("\r",b)
    wbi <- matrix(NA,nrow=lwp,ncol=ntrt)
    wbpooled <- numeric(lwp)
    for(i in 1:ntrt)
    {
      xbi <- sort(sample(t,size=n[i],prob=wpooled,replace=TRUE)+hpooled*stats::rnorm(n[i]))
      wbi[,i] <- calcw_cpp(xbi,y)
    }
    for(i in 1:lwp)
    {
      wbpooled[i] <- sum(wbi[i,]*n)
    }
    wbpooled <- wbpooled/npooled
    objb <- bw.dist.binned(npooled,y,wbpooled,print=FALSE,plot=FALSE)
    hb <- objb$h
    Fhb <- objb$Fh
    fhb <- Vectorize(function(x)1/hb*sum(wbpooled*stats::dnorm((x-t)/hb)))
    for(i in 1:ntrt)
    {
      Fhbi <- bw.dist.binned(n[i],y,wbi[,i],print=FALSE,plot=FALSE)$Fh
      Db[b] <- Db[b]+stats::integrate(function(x)(Fhbi(x)-Fhb(x))^2*fhb(x),lim1+hb*stats::qnorm(0.99),lim2+hb*stats::qnorm(0.01))$value*n[i]
    }
    
  }
  
  Db <- Db/ntrt
  
  pvalue <- mean(Db>D)
  
  return(pvalue)
  
}