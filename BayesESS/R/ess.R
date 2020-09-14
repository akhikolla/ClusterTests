#' @export
#' @importFrom Rcpp evalCpp
#' @importFrom stats runif approx rnorm rexp
#' @importFrom LaplacesDemon rinvchisq
#' @importFrom MCMCpack rinvgamma
#' @importFrom dfcrm titesim crmsim
#' @useDynLib BayesESS, .registration = TRUE

#####################################################################
##
## A wrapper function for ESS calculation written by Jaejoon Song
##
#####################################################################

ess <- function(model=NULL,label=NULL,
                prior=NULL,
                m=NULL,nsim=NULL,
                ncov=NULL,svec1=NULL,svec2=NULL,
                PI=NULL,betaSD=NULL,target=NULL,
                obswin=NULL,rate=NULL,accrual=NULL,
                shapeParam=NULL,scaleParam=NULL,
                fast=TRUE){

  Prior_0 <- Prior_1 <- Prior_2 <- Prior_3 <- Prior_4 <- Prior_5 <- NULL
  Prior_6 <- Prior_7 <- Prior_8 <- Prior_9 <- Prior_10 <- Prior_11 <- NULL

  modelStatement <- switch(EXPR = model,
                           'betaBin' = 'beta-binomial',
                           'Betabin' = 'beta-binomial',
                           'betabin' = 'beta-binomial',
                           'BetaBin' = 'beta-binomial',
                           'gammaEx' = 'gamma-exponential',
                           'Gammaex' = 'gamma-exponential',
                           'gammaex' = 'gamma-exponential',
                           'GammaEx' = 'gamma-exponential',
                           'dirMult' = 'dirichlet-multinomial',
                           'Dirmult' = 'dirichlet-multinomial',
                           'dirmult' = 'dirichlet-multinomial',
                           'DirMult' = 'dirichlet-multinomial',
                           'gammaPois' = 'gamma-poisson',
                           'GammaPois' = 'gamma-poisson',
                           'Gammapois' = 'gamma-poisson',
                           'gammapois' = 'gamma-poisson',
                           'normNorm' = 'normal-normal',
                           'NormNorm' = 'normal-normal',
                           'Normnorm' = 'normal-normal',
                           'normnorm' = 'normal-normal',
                           'tite.crm' = 'TITE CRM',
                           'crm' = 'CRM',
                           'surv' = 'time to event',
                           'invChisqNorm' = 'scaled inverse-chi-squared-normal',
                           'InvChisqNorm' = 'scaled inverse-chi-squared-normal',
                           'invchisqnorm' = 'scaled inverse-chi-squared-normal',
                           'invChiSqNorm' = 'scaled inverse-chi-squared-normal',
                           'invGammaNorm' = 'inverse-gamma-normal',
                           'InvGammaNorm' = 'inverse-gamma-normal',
                           'invgammanorm' = 'inverse-gamma-normal')

  ##
  ## ESS for conjugate models
  ##
  if(model %in% c('betaBin','Betabin','betabin','BetaBin')){

    if( as.numeric(prior[2]) <= 0 ) stop('For Beta(alpha,beta) distribution, alpha should be alpha > 0')
    if( as.numeric(prior[3]) <= 0 ) stop('For Beta(alpha,beta) distribution, beta should be beta > 0')

    label <- c('alpha','beta')
    overall_out_label <- paste("(",paste(label,collapse=","),")",sep="")
    myESS <- as.numeric(prior[2]) + as.numeric(prior[3])
    cat("\n")
    cat(noquote(paste("ESS was calculated for a ",modelStatement," model", sep="")))
    cat("\n")
    cat("\n")
    cat(noquote(paste("ESS for the beta",overall_out_label," prior is: ",myESS, sep="")))
    cat("\n")
    cat("\n")
    cat("\n")
    ESS <- myESS
  }


  if(model %in% c('gammaEx','Gammaex','gammaex','GammaEx')){

    if( as.numeric(prior[2]) <= 0 ) stop('For Gamma(alpha,beta) distribution, alpha should be alpha > 0')
    if( as.numeric(prior[3]) <= 0 ) stop('For Gamma(alpha,beta) distribution, beta should be beta > 0')

    label <- c('alpha','beta')
    overall_out_label <- paste("(",paste(label,collapse=","),")",sep="")
    myESS <- as.numeric(prior[3])
    cat("\n")
    cat(noquote(paste("ESS was calculated for a ",modelStatement," model", sep="")))
    cat("\n")
    cat("\n")
    cat(noquote(paste("ESS for the gamma",overall_out_label," prior is: ",myESS, sep="")))
    cat("\n")
    cat("\n")
    cat("\n")
    ESS <- myESS
  }


  if(model %in% c('gammaPois', 'GammaPois', 'gammapois', 'Gammapois')){

    if( as.numeric(prior[2]) <= 0 ) stop('For Gamma(alpha,beta) distribution, alpha should be alpha > 0')
    if( as.numeric(prior[3]) <= 0 ) stop('For Gamma(alpha,beta) distribution, beta should be beta > 0')

    label <- c('alpha','beta')
    overall_out_label <- paste("(",paste(label,collapse=","),")",sep="")
    myESS <- as.numeric(prior[3])
    cat("\n")
    cat(noquote(paste("ESS was calculated for a ",modelStatement," model", sep="")))
    cat("\n")
    cat("\n")
    cat(noquote(paste("ESS for the gamma",overall_out_label," prior is: ",myESS, sep="")))
    cat("\n")
    cat("\n")
    cat("\n")
    ESS <- myESS
  }


  if(model %in% c('normNorm', 'NormNorm', 'normnorm', 'Normnorm')){

    label <- c('mu_0','sigma^2/n_0')
    overall_out_label <- paste("(",paste(label,collapse=","),")",sep="")

    prior <- suppressWarnings(as.numeric(gsub("[^0-9\\.]", "", prior)))
    prior <- prior[!is.na(prior)]

    itIsInteger <- round(as.numeric(prior[3])) == as.numeric(prior[3])

    if( as.numeric(prior[2]) <= 0 ) stop(paste('For Normal ',overall_out_label,' distribution, sigma^2 should be sigma^2 > 0',sep=""))
    if( as.numeric(prior[3]) <= 0 | itIsInteger == FALSE) stop(paste('For Normal ',overall_out_label,' distribution, n_0 should be a positive integer',sep=""))

    myESS <- as.numeric(prior[3])
    cat("\n")
    cat(noquote(paste("ESS was calculated for a ",modelStatement," model", sep="")))
    cat("\n")
    cat("\n")
    cat(noquote(paste("ESS for the normal",overall_out_label," prior is: ",myESS, sep="")))
    cat("\n")
    cat("\n")
    cat("\n")
    ESS <- myESS
  }



  if(model %in% c('dirMult','Dirmult','dirmult','DirMult')){
    label <- 'alpha1'
    for(i in 2:(length(prior)-1)){
      label <- c(label,paste('alpha',i,sep=""))
    }

    overall_out_label <- paste("(",paste(label,collapse=","),")",sep="")

    aProblem1 <- sum( as.numeric(prior[2:length(prior)]) <= 0 ) > 0
    aProblem2 <- sum( round(as.numeric(prior[2:length(prior)])) ==  as.numeric(prior[2:length(prior)]) ) < ( length(prior) - 1 )

    if( aProblem1 == TRUE | aProblem2 == TRUE) stop(paste('For dirichlet ',overall_out_label,' distribution, all alpha_i should be a positive integers',sep=""))

    myESS <- sum(as.numeric(prior[2:length(prior)]))

    cat("\n")
    cat(noquote(paste("ESS was calculated for a ",modelStatement," model", sep="")))
    cat("\n")
    cat("\n")
    cat(noquote(paste("ESS for the dirichlet",overall_out_label," prior is: ",myESS, sep="")))
    cat("\n")
    cat("\n")
    cat("\n")
    ESS <- myESS
  }


  ##
  ## ESS for linear regression or logistic regression models
  ##
  if(model %in% c('linreg','Linreg','logistic','Logistic')){
    default_prior <- list(c(1,0,1000),c(1,0,1000),c(1,0,1000),
                          c(1,0,1000),c(1,0,1000),c(1,0,1000),
                          c(1,0,1000),c(1,0,1000),c(1,0,1000),
                          c(1,0,1000),c(1,0,1000),c(1,0,1000))

    if(length(svec1)<12){svec1 <- c(svec1,rep(0,(12-length(svec1))))}
    if(length(svec2)<12){svec2 <- c(svec2,rep(0,(12-length(svec2))))}

    if(model=='linreg'){model <- 1}
    if(model=='logistic'){model <- 2}

    for(i in 1:length(prior)){
      prior[[i]][prior[[i]] %in% c('norm','normal','N','Norm','Normal')] <- 1
      prior[[i]][prior[[i]] %in% c('gamma','Gamma','Gam','gam')] <- 2
      prior[[i]] <- as.numeric(prior[[i]])
    }

    for(i in 1:length(default_prior)){
      assign(paste('Prior_',(i-1),sep=''),default_prior[[i]])
    }

    if(length(prior)>0){
      for(i in 1:length(prior)){
        assign(paste('Prior_',(i-1),sep=''),prior[[i]])
      }
    }

    if(is.null(label)){
      svec <- svec1 + svec2
      svec[1] <- 1
      label <- 'theta1'
      for(i in 2:sum(svec)){
        label <- c(label,paste('theta',i,sep=""))
      }
    }

    svec1_label <- svec1[1:length(label)]
    svec2_label <- svec2[1:length(label)]

    svec1_label <- label[which(svec1_label==1)]
    svec2_label <- label[which(svec2_label==1)]

    overall_out_label <- paste("(",paste(label,collapse=","),")",sep="")
    svec1_out_label <- paste("(",paste(svec1_label,collapse=","),")",sep="")
    svec2_out_label <- paste("(",paste(svec2_label,collapse=","),")",sep="")

    if(fast==FALSE){
      myESS <- ESS_RegressionCalc (Reg_model = model, Num_cov = ncov,
                                   Prior_0 = Prior_0, Prior_1 = Prior_1, Prior_2 = Prior_2,
                                   Prior_3 = Prior_3, Prior_4 = Prior_4, Prior_5 = Prior_5,
                                   Prior_6 = Prior_6, Prior_7 = Prior_7, Prior_8 = Prior_8,
                                   Prior_9 = Prior_9, Prior_10 = Prior_10, Prior_11 = Prior_11,
                                   M = m, NumSims = nsim,
                                   theta_sub1=svec1, theta_sub2=svec2 )
    }

    if(fast==TRUE){
      myESS <- ESS_RegressionCalcFast1 (Reg_model = model, Num_cov = ncov,
                                        Prior_0 = Prior_0, Prior_1 = Prior_1, Prior_2 = Prior_2,
                                        Prior_3 = Prior_3, Prior_4 = Prior_4, Prior_5 = Prior_5,
                                        Prior_6 = Prior_6, Prior_7 = Prior_7, Prior_8 = Prior_8,
                                        Prior_9 = Prior_9, Prior_10 = Prior_10, Prior_11 = Prior_11,
                                        M = m, NumSims = nsim,
                                        theta_sub1=svec1, theta_sub2=svec2 )
    }

    ESSoutput <- list(#ESSoverall = myESS$ESSwholetheta,
                      ESSsubvec1 = myESS$ESSsubvector1,
                      ESSsubvec2 = myESS$ESSsubvector2)

    modelStatement <- switch(EXPR = model, 'linreg' = 'linear regression',
                             'logistic' = 'logistic regression')
    cat("\n")
    cat(noquote(paste("ESS was calculated for a ",modelStatement," model", sep="")))
    cat("\n")
    cat("\n")
    #cat(noquote(paste("ESSoverall: Overall ESS for the whole vector ",overall_out_label, sep="")))
    #cat("\n")
    cat(noquote(paste("ESSsubvector1: ESS for the first sub-vector ",svec1_out_label, sep="")))
    cat("\n")
    cat(noquote(paste("ESSsubvector1: ESS for the second sub-vector ",svec2_out_label, sep="")))
    cat("\n")
    cat("\n")

    ESS <- ESSoutput
  }


  ##
  ## ESS for CRM models
  ## Uses the essCRM function in 'internal'
  ##

  # For CRM
  if(model=='crm'){

    myEssCRM <- essCRM(type=model,
                       prior=prior,
                       m=m,nsim=nsim,
                       PI=PI,
                       betaSD=betaSD,
                       target=target)
    label <- c('beta')

    cat("\n")
    cat(noquote(paste("ESS was calculated for a ",modelStatement," model", sep="")))
    cat("\n")
    cat("\n")
    cat(noquote(paste("ESS for the N(0,",round(betaSD^2,10),") ",label," prior"," is: ",round(myEssCRM,10), sep="")))
    cat("\n")
    cat("\n")
    cat("\n")
    ESS <- myEssCRM
  }

  # For TITE CRM
  if(model=='tite.crm'){

    myEssCRM <- essCRM(type=model,
                       prior=prior,
                       m=m,nsim=nsim,
                       PI=PI,
                       betaSD=betaSD,
                       target=target,
                       obswin=obswin,rate=rate,accrual=accrual)
    label <- c('beta')

    cat("\n")
    cat(noquote(paste("ESS was calculated for a ",modelStatement," model", sep="")))
    cat("\n")
    cat("\n")
    cat(noquote(paste("ESS for the N(0,",round(betaSD^2,10),") ",label," prior"," is: ",round(myEssCRM,10), sep="")))
    cat("\n")
    cat("\n")
    cat("\n")
    ESS <- myEssCRM
  }



  ##
  ## ESS for time to event models
  ## Uses the essSurv function in 'internal'
  ##

  if(model=='surv'){

    myEssSurv <- essSurv(shapeParam=shapeParam,
                         scaleParam=scaleParam,
                         m=m,nsim=nsim)

    label <- c('alpha','beta')
    overall_out_label <- paste("(",paste(label,collapse=","),")",sep="")

    cat("\n")
    cat(noquote(paste("ESS was calculated for a ",modelStatement," model", sep="")))
    cat("\n")
    cat("\n")
    cat(noquote(paste("ESS for the inverse-gamma ",overall_out_label," prior is: ",round(myEssSurv,2), sep="")))
    cat("\n")
    cat("\n")
    cat("\n")
    ESS <- myEssSurv
  }



  ##
  ## ESS for scaled-inverse-chi-squared-normal
  ## Uses the essNormal function in 'internal'
  ##

  if(model %in% c('invChisqNorm','InvChisqNorm','invchisqnorm','invChiSqNorm')){

    #if(is.character(prior) == TRUE){ prior <- prior[2:length(prior)] }

    prior <- suppressWarnings(as.numeric(gsub("[^0-9\\.]", "", prior)))
    prior <- prior[!is.na(prior)]

    if(length(prior)==2){

      aProblem1 <- sum( as.numeric(prior) <= 0 ) > 0
      if( aProblem1 == TRUE ) stop(paste("Both parameters for the scaled-inv-chi-squared prior should be positive", sep=""))

      myESSNorm <- essNormal(nu=prior[1],sigma0=prior[2],knownMean=TRUE,m=m,nsim=nsim)

      cat("\n")
      cat(noquote(paste("ESS was calculated for a ",modelStatement," model (with known mean and unknown variance)", sep="")))
      cat("\n")
      cat("\n")
      cat(noquote(paste("ESS for the scaled-inv-chi-squared(",prior[1],",",prior[2],") prior for variance"," is: ",myESSNorm, sep="")))
      cat("\n")
      cat("\n")
      cat("\n")
      ESS <- myESSNorm

    }

    if(length(prior)==4){

      if( length(prior) != 4 ) stop(paste("Please specify two parameters (nu_0,sigma^2) for scaled-inverse-chi-squared(nu_0,sigma^2) prior for variance,
                                     and two parameters (mu_0 and phi) for normal(mu_0,sigma^2/phi) prior for mean", sep=""))

      aProblem1 <- sum( as.numeric(prior[1:2]) <= 0 ) > 0
      if( aProblem1 == TRUE ) stop(paste("Both parameters for the scaled-inv-chi-squared prior should be positive", sep=""))

      itIsInteger <- round(as.numeric(prior[4])) == as.numeric(prior[4])
      if( as.numeric(prior[4]) <= 0 | itIsInteger == FALSE) stop(paste('For Normal (mu_0,sigma^2/phi) distribution, phi should be a positive integer',sep=""))

      myESSNorm <- essNormal(nu=prior[1],sigma0=prior[2],mu0=prior[3],phi=prior[4],knownMean=FALSE,m=m,nsim=nsim)

      cat("\n")
      cat(noquote(paste("ESS was calculated for a ",modelStatement," model (with unknown mean and unknown variance)", sep="")))
      cat("\n")
      cat("\n")
      #cat(noquote(paste("ESS_overall: Overall ESS for the priors specified for the mean and variance is ",myESSNorm$ESS_overall, sep="")))
      #cat("\n")
      cat(noquote(paste("ESS_sigmasq: ESS for the scaled-inv-chi-squared(",prior[1],",",prior[2],") prior for variance is ",round(myESSNorm$ESS_sigmasq,2), sep="")))
      cat("\n")
      cat(noquote(paste("ESS_mu: ESS for the normal prior for mean is ",round(myESSNorm$ESS_mu,2), sep="")))
      cat("\n")
      cat("\n")

      ESS <- myESSNorm
    }




  }



  ##
  ## ESS for inverse-gamma-normal
  ## Uses the essNormal function in 'internal'
  ##

  if(model %in% c('invGammaNorm','InvGammaNorm','inverse-gamma-normal','invgammanorm')){

    #if(is.character(prior) == TRUE){ prior <- prior[2:length(prior)] }
    prior <- suppressWarnings(as.numeric(gsub("[^0-9\\.]", "", prior)))
    prior <- prior[!is.na(prior)]

    if(length(prior)==2){

      aProblem1 <- sum( as.numeric(prior) <= 0 ) > 0
      if( aProblem1 == TRUE ) stop(paste("Both parameters for the inverse-gamma prior should be positive", sep=""))

      myESSNorm <- essNormal(nu=(prior[1]*2),sigma0=(prior[2]/prior[1]),knownMean=TRUE,m=m,nsim=nsim)

      cat("\n")
      cat(noquote(paste("ESS was calculated for a ",modelStatement," model (with known mean and unknown variance)", sep="")))
      cat("\n")
      cat("\n")
      cat(noquote(paste("ESS for the inv-gamma(",prior[1],",",prior[2],") prior for variance"," is: ",myESSNorm, sep="")))
      cat("\n")
      cat("\n")
      cat("\n")
      ESS <- myESSNorm
    }

    if(length(prior)==4){

      if( length(prior) != 4 ) stop(paste("Please specify two parameters (alpha,beta) for inverse-gamma (alpha,beta) prior for variance,
                                     and two parameters (mu_0 and phi) for normal(mu_0,sigma^2/phi) prior for mean", sep=""))

      aProblem1 <- sum( as.numeric(prior[1:2]) <= 0 ) > 0
      if( aProblem1 == TRUE ) stop(paste("Both parameters for the inverse-gamma prior should be positive", sep=""))

      itIsInteger <- round(as.numeric(prior[4])) == as.numeric(prior[4])
      if( as.numeric(prior[4]) <= 0 | itIsInteger == FALSE) stop(paste('For Normal (mu_0,sigma^2/phi) distribution, phi should be a positive integer',sep=""))


      myESSNorm <- essNormal(nu=(prior[1]*2),sigma0=(prior[2]/prior[1]),mu0=prior[3],phi=prior[4],knownMean=FALSE,m=m,nsim=nsim)

      cat("\n")
      cat(noquote(paste("ESS was calculated for a ",modelStatement," model (with unknown mean and unknown variance)", sep="")))
      cat("\n")
      cat("\n")
      #cat(noquote(paste("ESS_overall: Overall ESS for the priors specified for the mean and variance is ",myESSNorm$ESS_overall, sep="")))
      #cat("\n")
      cat(noquote(paste("ESS_sigmasq: ESS for the inv-gamma(",prior[1],",",prior[2],") prior for variance is ",round(myESSNorm$ESS_sigmasq,2), sep="")))
      cat("\n")
      cat(noquote(paste("ESS_mu: ESS for the normal prior for mean is ",round(myESSNorm$ESS_mu,2), sep="")))
      cat("\n")
      cat("\n")

      ESS <- myESSNorm
    }

  }



  ESS
}
