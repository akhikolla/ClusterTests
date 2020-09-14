#' Bootstrap and plug-in bandwidth selectors for kernel density estimation with binned data.
#'
#' @param n Positive integer. Size of the complete sample.
#' @param y Vector. Observed values. They define the extremes of the sequence of intervals in which data is binned.
#' @param w Vector. Proportion of observations within each interval.
#' @param ni Vector. Number of observations within each interval.
#' @param gboot Positive real number. Pilot bandwidth for the bootstrap bandwidth selector.
#' @param pilot.type 1, 2 or 3. If \code{gboot} is missing, pilot bandwidth for the bootstrap bandwidth selector is automatically selected using methods 1, 2 or 3. Defaults to 3. See details for more information.
#' @param hn Positive integer. Size of the grid of bandwidths for which MISE will be approximated. Defaults to 100.
#' @param plugin.type Character. If \code{plugin.type} = "N", normality is assumed at the last step when calculating the plug-in bandwidth. If \code{plugin.type} = "A", parameter at last step is estimated nonparametrically using \code{gplugin} as bandwidth. Otherwise, the unknown parameter is estimated fitting a normal mixture. Defaults to \code{type} = "N".
#' @param confband Logical. If TRUE, bootstrap confidence bands are constructed for the density function. Defaults to FALSE.
#' @param alpha Real number between 0 and 1. Significance level for the bootstrap confidence bands. Defaults to 0.05.
#' @param B Positive integer. Number of bootstrap resamples used when constructing confidence bands. Defaults to 1000.
#' @param plot Logical. If TRUE, kernel density estimators are plotted along with (optional) bootstrap confidence bands. Defaults to TRUE.
#' @param print Logical. If TRUE and confband is TRUE, the percentage of bootstrap resamples already evaluated is printed. Defaults to TRUE.
#' @param model Character. Name of the parametric family of distributions to be fitted for the grouped sample. Parameters are estimated by maximum likelihood.
#' @param parallel Logical. If TRUE, confidence bands are estimated using parallel computing with sockets.
#' @param pars Environment. Needed for the well functioning of the script. DO NOT modify this argument.
#'
#' @details
#' If \code{pilot.type} = 1, an heuristic rule is used for calculating the pilot bandwidth. It's not recommended when population's density function is suspected to be highly multimodal.
#'
#' If \code{pilot.type} = 2, the pilot bandwidth is such that the kernel density estimator with bandwidth \code{gboot} approximates the histogram of the grouped sample minimizing the residual sum of squares. If \code{pilot.type} = 3, a penalty is imposed on the curvature of the kernel density estimator with bandwidth \code{gboot}. The penalty parameter is selected as to best approximate the curvature of the true density.
#'
#' @return A list with components
#' \item{h_boot}{Bootstrap bandwidth selector.}
#' \item{h_plugin}{Plug-in bandwidth selector.}
#'
#' @examples
#' set.seed(1)
#' n <- 200 #complete sample size
#' k <- 30 #number of intervals
#' x <- rnorm(n,6,1) #complete sample
#' y <- seq(min(x)-0.2,max(x)+0.2,len=k+1) #intervals
#' w <- c(sapply(2:k,function(i)sum( x<y[i]&x>=y[i-1] )), sum(x<=y[k+1]&x>=y[k]) )/n #proportions
#' bw.dens.binned(n,y,w,plot=FALSE)
#'
#' @references
#' \insertRef{TesisMiguel2015}{binnednp}
#'
#' \insertRef{JNS2016}{binnednp}
#'
#' \insertRef{Test2017}{binnednp}
#'
#' @useDynLib binnednp
#' @importFrom Rcpp sourceCpp
#'
#' @export
bw.dens.binned <- function(n,y,w,ni,gboot,pilot.type=3,hn=100,plugin.type="N",confband=FALSE,alpha=0.05,B=1000,plot=TRUE,print=TRUE,model,parallel=FALSE,pars=new.env()){


  main <- !exists("k", where = pars, inherits = FALSE)

  if(main == TRUE){
    t <- y[-length(y)]+diff(y)/2
    k <- length(t)

    if(missing(w)){

      if(!missing(ni)){
        if(missing(n)) n <- floor(sum(ni))
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


  } else{
    t <- pars$t
    k <- pars$k
    y <- pars$y
  }
  mu_hat <- sum(w*t)
  sigma_hat <- sqrt(sum(w*(t-mu_hat)^2))

  if(main == TRUE){
    l <- mean(diff(y))
    r <- max(t)-min(t)


    K <- function(x) stats::dnorm(x)
    AK <- 1/(2*sqrt(pi))

    K4_0 <- kedd::kernel.fun(0,deriv.order=4)$kx
    K6_0 <- kedd::kernel.fun(0,deriv.order=6)$kx
    K8_0 <- kedd::kernel.fun(0,deriv.order=8)$kx
    K10_0 <- kedd::kernel.fun(0,deriv.order=10)$kx

    psi4 <- function(w,h) sum( sapply( 1:k, function(i){
      sum( kedd::kernel.fun((t[i]-t)/h, deriv.order=4)$kx*w[i]*w )
    } ) )/h^5

    psi6 <- function(w,h) sum( sapply( 1:k, function(i){
      sum( kedd::kernel.fun((t[i]-t)/h, deriv.order=6)$kx*w[i]*w )
    } ) )/h^7

    psi8 <- function(w,h) sum( sapply( 1:k, function(i){
      sum( kedd::kernel.fun((t[i]-t)/h, deriv.order=8)$kx*w[i]*w )
    } ) )/h^9

    psi10 <- function(w,h) sum( sapply( 1:k, function(i){
      sum( kedd::kernel.fun((t[i]-t)/h, deriv.order=10)$kx*w[i]*w )
    } ) )/h^11


  } else{
    l <- pars$l
    r <- pars$r
    K <- pars$K
    AK <- pars$AK
    K4_0 <- pars$K4_0
    K6_0 <- pars$K6_0
    K8_0 <- pars$K8_0
    psi4 <- pars$psi4
    psi6 <- pars$psi6
    psi8 <- pars$psi8
    psi10 <- pars$psi10
  }

  if(plugin.type == "N" || plugin.type == "A"){


  if(plugin.type == "A"){
    hAf <- 1.06*sigma_hat*n^(-1/5)
    fhAf <- Vectorize(function(x)1/hAf*sum(w*stats::dnorm((x-t)/hAf)))
    try(Af <- stats::integrate(function(x)fhAf(x)^2,-Inf,Inf)$value, silent=TRUE)
    if(!exists("Af")) Af <- stats::integrate(function(x)fhAf(x)^2,mu_hat-3*sigma_hat,mu_hat+3*sigma_hat)$value
    psi10 <- psi10(w,hAf)
  } else {
    Af <- 1/(2*sqrt(pi)*sigma_hat)
    psi10 <- -30240/((2*sigma_hat)^11*sqrt(pi))
  }

  eta8 <- (-2*K8_0*Af*l/psi10)^(1/11)

  eta6 <- (-2*K6_0*Af*l/psi8(w,eta8))^(1/9)

  eta4 <- (-2*K4_0*Af*l/psi6(w,eta6))^(1/7)

  Af2 <- psi4(w,eta4)

  } else {

    clustering <- mclust::Mclust(rep(t,n*w),G=1:5,verbose=FALSE,modelNames="V")
    mix_params <- clustering$parameters
    mu_mixt <- mix_params$mean
    sigma_mixt <- sqrt(mix_params$variance$sigmasq)
    alfa_mixt <- mix_params$pro
    normal_mixt <- nor1mix::norMix(mu=mu_mixt,sigma=sigma_mixt,w=alfa_mixt)
    f2_mixt <- Vectorize( function(x) -sum( alfa_mixt*stats::dnorm((x-mu_mixt)/sigma_mixt)*( 1/sigma_mixt^3-(x-mu_mixt)^2/sigma_mixt^5 ) ) )
    mixlim1 <- nor1mix::qnorMix(0.001,normal_mixt)
    mixlim2 <- nor1mix::qnorMix(0.999,normal_mixt)
    Af2 <- stats::integrate(function(x)f2_mixt(x)^2,mixlim1,mixlim2)$value

  }

  h_plugin <- (AK/(Af2*n))^(1/5)


  if(main == TRUE){
    aux <- NULL
    Fplugin <- Vectorize(function(aux)sum(w*stats::dnorm((aux-t)/h_plugin))/h_plugin)

    if(!missing(model)){
      dmodel <- paste0("d",model)
      try(parfit <- fitdistrplus::fitdist(rep(t,ceiling(n*w)),model),silent=T)
      if(exists("parfit")){
        distr <- get(dmodel)
        params <- as.list(parfit$estimate)
        pargrid <- seq(min(y),max(y),len=101)
        clist <- c(list(x=pargrid),params)
      }
    }
  }




  if(main == TRUE){
    omega <- l/r
  } else{
    omega <- pars$omega
  }

  if(missing(gboot)){

    if(pilot.type == 1){
      eta.norm <- 0.78*sigma_hat*n^(-2/13)

      if(n <= 150) gboot <- ifelse(omega <= 0.15,0.78*eta.norm,eta.norm*(4*omega+0.4))
      if(n > 150) gboot <- ifelse(omega <= 0.075,eta.norm,eta.norm*(7.2*omega+0.46))
    }

    if(pilot.type == 2){

      hist <- hist(rep(t,n*w),breaks=y,plot=FALSE)$density
      fht <- function(h)sapply(t,function(x)sum(w*stats::dnorm((x-t)/h))/h)
      h0 <- min(diff(y))/4
      h1 <- max(y)-min(y)
      rho <- exp((log(h1)-log(h0))/4)

      it <- 0
      while(it < 10){
        hseq <- h0*rho^(0:4)
        dist <- sapply(hseq,function(h)sum((hist-fht(h))^2))
        minind <- which.min(dist)
        gboot <- hseq[minind]
        if(minind == 1){
          h1 <- hseq[2]
          h0 <- gboot/rho
          rho <- exp((log(h1)-log(h0))/4)
        } else {
          if(minind == 5){
            h0 <- hseq[4]
            h1 <- gboot*rho
            rho <- exp((log(h1)-log(h0))/4)
          } else {
            h0 <- hseq[minind-1]
            h1 <- hseq[minind+1]
            rho <- exp((log(h1)-log(h0))/4)
          }
        }
        it <- it+1
      }


    }



    if(pilot.type == 3)
    {

      clustering <- mclust::Mclust(rep(t,n*w),G=1:5,verbose=FALSE,modelNames="V")
      mix_params <- clustering$parameters
      mu_mixt <- mix_params$mean
      sigma_mixt <- sqrt(mix_params$variance$sigmasq)
      alfa_mixt <- mix_params$pro
      normal_mixt <- nor1mix::norMix(mu=mu_mixt,sigma=sigma_mixt,w=alfa_mixt)
      f2_mixt <- Vectorize( function(x) -sum( alfa_mixt*stats::dnorm((x-mu_mixt)/sigma_mixt)*( 1/sigma_mixt^3-(x-mu_mixt)^2/sigma_mixt^5 ) ) )
      mixlim1 <- nor1mix::qnorMix(0.001,normal_mixt)
      mixlim2 <- nor1mix::qnorMix(0.999,normal_mixt)
      Af2_mixt <- stats::integrate(function(x)f2_mixt(x)^2,mixlim1,mixlim2)$value


      hist <- hist(rep(t,n*w),breaks=y,plot=FALSE)$density
      fht <- function(h)sapply(t,function(x)sum(w*stats::dnorm((x-t)/h))/h)
      d2fht <- Vectorize( function(x,h)sum(w*kedd::kernel.fun((x-t)/h,2)$kx)/h^3, "x" )

      lim_int_1 <- mixlim1
      lim_int_2 <- mixlim2

      l0 <- 1e-3
      l1 <- 10
      lrho <- exp((log(l1)-log(l0))/4)

      h0 <- min(diff(y))/4
      h1 <- max(y)-min(y)
      rho <- exp((log(h1)-log(h0))/4)

      gboot <- zeta_hist_p(hist,t,y,w,Af2_mixt,l0,l1,h0,h1,lrho,rho,10,10,lim_int_1,lim_int_2)


    }

  }

  wg <- sapply( 2:(k+1), function(i) sum( w*( stats::pnorm((y[i]-t)/gboot) - stats::pnorm((y[i-1]-t)/gboot) ) ) )



  gaussian_mise_initialize(n, k, wg, w, t, gboot, AK)
  h_boot <- gaussian_dichotomy(hn, t)

  if(main == TRUE){
    aux <- NULL
    Fboot <- Vectorize(function(aux)sum(w*stats::dnorm((aux-t)/h_boot))/h_boot)
  }


  if(confband == T){

    pars$t <- t
    pars$k <- k
    pars$y <- y
    pars$l <- l
    pars$r <- r
    pars$omega <- omega
    pars$K <- K
    pars$AK <- AK
    pars$K4_0 <- K4_0
    pars$K6_0 <- K6_0
    pars$K8_0 <- K8_0
    pars$psi4 <- psi4
    pars$psi6 <- psi6
    pars$psi8 <- psi8
    pars$psi10 <- psi10

    fgbooty <- sapply(y, function(x)sum( w*stats::dnorm((x-t)/gboot) )/gboot)
    fhbooty <- sapply(y, function(x)sum( w*stats::dnorm((x-t)/h_boot) )/h_boot)
    fhpluginy <- sapply(y, function(x)sum( w*stats::dnorm((x-t)/h_plugin) )/h_plugin)
    Dnboot <- matrix(0,nrow=k+1,ncol=B)
    Dnplugin <- matrix(0,nrow=k+1,ncol=B)

    if(!parallel){

    for(b in 1:B){
      rx <- sort(sample(t,replace=TRUE,size=n,prob=w)+gboot*stats::rnorm(n))
      wb <- calcw_cpp(rx,y)
      obj <- bw.dens.binned(n,y,wb,hn=10,plot=FALSE,pilot.type=pilot.type,plugin.type=plugin.type,pars=pars)
      h_boot.boot <- obj$h_boot
      h_plugin.boot <- obj$h_plugin
      Dnboot[,b] <- sapply(y,function(x)sum( wb*stats::dnorm( (x-t)/h_boot.boot ) )/h_boot.boot)
      Dnplugin[,b] <- sapply(y,function(x)sum( wb*stats::dnorm( (x-t)/h_plugin.boot ) )/h_plugin.boot)
      if(print == TRUE) cat("\rConstructing confidence bands... Progress:",floor(100*b/B),"%")
    }

    } else {

      parfun <- function(b){
        rx <- sort(sample(t,replace=TRUE,size=n,prob=w)+gboot*stats::rnorm(n))
        wb <- calcw_cpp(rx,y)
        obj <- bw.dens.binned(n,y,wb,hn=10,plot=FALSE,pilot.type=pilot.type,plugin.type=plugin.type,pars=pars)
        h_boot.boot <- obj$h_boot
        h_plugin.boot <- obj$h_plugin
        Dnbootcolb <- sapply(y,function(x)sum( wb*stats::dnorm( (x-t)/h_boot.boot ) )/h_boot.boot)
        Dnplugincolb <- sapply(y,function(x)sum( wb*stats::dnorm( (x-t)/h_plugin.boot ) )/h_plugin.boot)
        return(c(boot=Dnbootcolb,pi=Dnplugincolb))
        }

      ncores <- parallel::detectCores()
      cl <- parallel::makeCluster(ncores)
      parallel::clusterEvalQ(cl, library(binnednp))
      parallel::clusterExport(cl, 'calcw_cpp')
      paroutput <- parallel::parSapply(cl, 1:B, parfun)
      Dn <- paroutput
      idx <- seq(1,2*(k+1),by=2)
      Dnboot <- Dn[1:(k+1),]
      Dnplugin <- Dn[(k+2):(2*(k+1)),]
      parallel::stopCluster(cl)

    }

    alfa <- alpha/(k+1)

    count <- 0
    maxcount <- 100
    low.alpha_boot <- alfa
    high.alpha_boot <- alpha
    low.alpha_plugin <- alfa
    high.alpha_plugin <- alpha
    if(print == TRUE) cat("\n")
    while(count < maxcount){
      if(print == TRUE) cat("\rAdjusting significance level. Progress:",floor(100*(count+1)/maxcount),"%")

      mean.alpha_boot <- 0.5*(low.alpha_boot+high.alpha_boot)

      q1 <- apply(Dnboot,1,stats::quantile,probs=low.alpha_boot/2)
      q2 <- apply(Dnboot,1,stats::quantile,probs=1-low.alpha_boot/2)
      banda1 <- pmin(pmax(fgbooty+fhbooty-q1,0),1)
      banda2 <- pmin(pmax(fgbooty+fhbooty-q2,0),1)
      p_low.alpha_boot <- sapply(1:B,function(b){
        aux <- ifelse(Dnboot[,b] < banda1 && Dnboot [,b] > banda2,1,0)
        return(aux)
      })
      p_low.alpha_boot <- sum(p_low.alpha_boot)/B

      q1 <- apply(Dnboot,1,stats::quantile,probs=high.alpha_boot/2)
      q2 <- apply(Dnboot,1,stats::quantile,probs=1-high.alpha_boot/2)
      banda1 <- pmin(pmax(fgbooty+fhbooty-q1,0),1)
      banda2 <- pmin(pmax(fgbooty+fhbooty-q2,0),1)
      p_high.alpha_boot <- sapply(1:B,function(b){
        aux <- ifelse(Dnboot[,b] < banda1 && Dnboot[,b] > banda2,1,0)
        return(aux)
      })
      p_high.alpha_boot <- sum(p_high.alpha_boot)/B

      q1 <- apply(Dnboot,1,stats::quantile,probs=mean.alpha_boot/2)
      q2 <- apply(Dnboot,1,stats::quantile,probs=1-mean.alpha_boot/2)
      banda1 <- pmin(pmax(fgbooty+fhbooty-q1,0),1)
      banda2 <- pmin(pmax(fgbooty+fhbooty-q2,0),1)
      p_mean.alpha_boot <- sapply(1:B,function(b){
        aux <- ifelse(Dnboot[,b] < banda1 && Dnboot[,b] > banda2,1,0)
        return(aux)
      })
      p_mean.alpha_boot <- sum(p_mean.alpha_boot)/B

      if(p_mean.alpha_boot >= 1-alpha){
        low.alpha_boot <- mean.alpha_boot
      } else {
        high.alpha_boot <- mean.alpha_boot
      }




      mean.alpha_plugin <- 0.5*(low.alpha_plugin+high.alpha_plugin)

      q1 <- apply(Dnplugin,1,stats::quantile,probs=low.alpha_plugin/2)
      q2 <- apply(Dnplugin,1,stats::quantile,probs=1-low.alpha_plugin/2)
      banda1 <- pmin(pmax(fhpluginy+fgbooty-q1,0),1)
      banda2 <- pmin(pmax(fhpluginy+fgbooty-q2,0),1)
      p_low.alpha_plugin <- sapply(1:B,function(b){
        aux <- ifelse(Dnplugin[,b] < banda1 && Dnplugin [,b] > banda2,1,0)
        return(aux)
      })
      p_low.alpha_plugin <- sum(p_low.alpha_plugin)/B

      q1 <- apply(Dnplugin,1,stats::quantile,probs=high.alpha_plugin/2)
      q2 <- apply(Dnplugin,1,stats::quantile,probs=1-high.alpha_plugin/2)
      banda1 <- pmin(pmax(fhpluginy+fgbooty-q1,0),1)
      banda2 <- pmin(pmax(fhpluginy+fgbooty-q2,0),1)
      p_high.alpha_plugin <- sapply(1:B,function(b){
        aux <- ifelse(Dnplugin[,b] < banda1 && Dnplugin[,b] > banda2,1,0)
        return(aux)
      })
      p_high.alpha_plugin <- sum(p_high.alpha_plugin)/B

      q1 <- apply(Dnplugin,1,stats::quantile,probs=mean.alpha_plugin/2)
      q2 <- apply(Dnplugin,1,stats::quantile,probs=1-mean.alpha_plugin/2)
      banda1 <- pmin(pmax(fhpluginy+fgbooty-q1,0),1)
      banda2 <- pmin(pmax(fhpluginy+fgbooty-q2,0),1)
      p_mean.alpha_plugin <- sapply(1:B,function(b){
        aux <- ifelse(Dnplugin[,b] < banda1 && Dnplugin[,b] > banda2,1,0)
        return(aux)
      })
      p_mean.alpha_plugin <- sum(p_mean.alpha_plugin)/B

      if(p_mean.alpha_plugin >= 1-alpha){
        low.alpha_plugin <- mean.alpha_plugin
      } else {
        high.alpha_plugin <- mean.alpha_plugin
      }




      count <- count+1
    }

    q1 <- apply(Dnboot,1,stats::quantile,probs=mean.alpha_boot/2)
    q2 <- apply(Dnboot,1,stats::quantile,probs=1-mean.alpha_boot/2)

    banda1boot <- pmin(pmax(fgbooty+fhbooty-q1,0),1)
    banda2boot <- pmin(pmax(fgbooty+fhbooty-q2,0),1)

    confbands.boot <- cbind(banda1boot,banda2boot)


    q1 <- apply(Dnplugin,1,stats::quantile,probs=alfa/2)
    q2 <- apply(Dnplugin,1,stats::quantile,probs=1-alfa/2)

    banda1plugin <- pmin(pmax(fhpluginy+fgbooty-q1,0),1)
    banda2plugin <- pmin(pmax(fhpluginy+fgbooty-q2,0),1)

    confbands.plugin <- cbind(banda1plugin,banda2plugin)


    if(plot == TRUE){

      fh1 <- Vectorize( function(x)sum( w*stats::dnorm((x-t)/h_boot) )/h_boot )
      fh2 <- Vectorize( function(x)sum( w*stats::dnorm((x-t)/h_plugin) )/h_plugin )

      gub <- Vectorize( function(t){
        j <- max(which(y<t))
        a <- y[j]
        d <- confbands.boot[j,1]
        fa <- fh1(a)
        return(fh1(t)+d-fa)
      } )

      hub <- Vectorize( function(t){
        j <- max(which(y<t))
        a <- y[j]
        b <- y[j+1]
        f <- confbands.boot[j+1,1]
        return(gub(t)+(t-a)/(b-a)*(f-gub(b)))
      } )

      glb <- Vectorize( function(t){
        j <- max(which(y<t))
        a <- y[j]
        d <- confbands.boot[j,2]
        fa <- fh1(a)
        return(fh1(t)+d-fa)
      } )

      hlb <- Vectorize( function(t){
        j <- max(which(y<t))
        a <- y[j]
        b <- y[j+1]
        f <- confbands.boot[j+1,2]
        return(glb(t)+(t-a)/(b-a)*(f-glb(b)))
      } )


      gup <- Vectorize( function(t){
        j <- max(which(y<t))
        a <- y[j]
        d <- confbands.plugin[j,1]
        fa <- fh2(a)
        return(fh2(t)+d-fa)
      } )

      hup <- Vectorize( function(t){
        j <- max(which(y<t))
        a <- y[j]
        b <- y[j+1]
        f <- confbands.plugin[j+1,1]
        return(gup(t)+(t-a)/(b-a)*(f-gup(b)))
      } )

      glp <- Vectorize( function(t){
        j <- max(which(y<t))
        a <- y[j]
        d <- confbands.plugin[j,2]
        fa <- fh2(a)
        return(fh2(t)+d-fa)
      } )

      hlp <- Vectorize( function(t){
        j <- max(which(y<t))
        a <- y[j]
        b <- y[j+1]
        f <- confbands.plugin[j+1,2]
        return(glp(t)+(t-a)/(b-a)*(f-glp(b)))
      } )


      grid <- seq(y[1],y[k+1],len=500)[-c(1,500)]


      grDevices::dev.new(noRStudioGD=TRUE)
      graphics::curve(fh2,min(y),max(y),type="n",ylim=c(0,max(banda1plugin)),lwd=2,ylab="Density",main="Kernel density (plug-in selector)")
      graphics::polygon(c(grid,rev(grid)),c(hlp(grid),rev(hup(grid))),col="grey",border=NA)
      graphics::curve(fh2,add=TRUE,lwd=2,ylab="Density",main="Kernel density (plug-in selector)")
      graphics::lines(grid,hup(grid),col=2,lty=2)
      graphics::lines(grid,hlp(grid),col=2,lty=2)

      if(exists("parfit")) graphics::lines(pargrid, do.call( distr, clist ), col=3, lwd=2 )

      grDevices::dev.new(noRStudioGD=TRUE)
      graphics::curve(fh1,min(y),max(y),type="n",ylim=c(0,max(banda1boot)),lwd=2,ylab="Density",main="Kernel density (bootstrap selector)")
      graphics::polygon(c(grid,rev(grid)),c(hlb(grid),rev(hub(grid))),col="grey",border=NA)
      graphics::curve(fh1,add=TRUE,lwd=2,ylab="Density",main="Kernel density (bootstrap selector)")
      graphics::lines(grid,hub(grid),col=2,lty=2)
      graphics::lines(grid,hlb(grid),col=2,lty=2)

      if(exists("parfit"))graphics::lines(pargrid, do.call( distr, clist ), col=3, lwd=2 )
    }

    return(list(h_plugin=h_plugin,h_boot=h_boot,
                confbands.boot=confbands.boot,confbands.plugin=confbands.plugin))
  } else {


    if(plot == TRUE){
      fh1 <- Vectorize(function(x)sum(w*stats::dnorm((x-t)/h_boot)))
      fh2 <- Vectorize(function(x)sum(w*stats::dnorm((x-t)/h_plugin)))
      grDevices::dev.new(noRStudioGD=TRUE)
      graphics::curve(fh2,min(y),max(y),lwd=2,ylab="Density",main="Kernel density (plug-in selector)")

      if(exists("parfit")) graphics::lines(pargrid, do.call( distr, clist ), col=3, lwd=2 )

      grDevices::dev.new(noRStudioGD=TRUE)
      graphics::curve(fh1,min(y),max(y),lwd=2,ylab="Density",main="Kernel density (bootstrap selector)")

      if(exists("parfit")) graphics::lines(pargrid, do.call( distr, clist ), col=3, lwd=2 )
    }


    return(list(h_plugin=h_plugin,h_boot=h_boot))

  }



}
