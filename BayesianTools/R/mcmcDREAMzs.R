### DREAMzs algorithm

##' DREAMzs
##' @author Stefan Paul
##' @param bayesianSetup Object of class 'bayesianSetup' or 'bayesianOuput'.
##' @param settings list with parameter values
##' @param iterations Number of model evaluations
##' @param nCR parameter determining the number of cross-over proposals.  If nCR = 1 all parameters are updated jointly.
##' @param updateInterval determining the intervall for the pCR (crossover probabilities) update
##' @param gamma Kurtosis parameter Bayesian Inference Scheme. 
##' @param eps Ergodicity term
##' @param e Ergodicity term
##' @param pCRupdate Update of crossover probabilities
##' @param burnin number of iterations treated as burn-in. These iterations are not recorded in the chain.
##' @param thin thin thinning parameter. Determines the interval in which values are recorded.
##' @param adaptation Number or percentage of samples that are used for the adaptation in DREAM (see Details)
##' @param DEpairs Number of pairs used to generate proposal
##' @param ZupdateFrequency frequency to update Z matrix
##' @param pSnooker probability of snooker update
##' @param Z starting matrix for Z
##' @param startValue eiter a matrix containing the start values (see details), an integer to define the number of chains that are run, a function to sample the start values or NUll, in which case the values are sampled from the prior.
##' @param consoleUpdates Intervall in which the sampling progress is printed to the console
##' @param message logical determines whether the sampler's progress should be printed
##' @return mcmc.object containing the following elements: chains, X, pCR, Z
##' @references Vrugt, Jasper A., et al. "Accelerating Markov chain Monte Carlo simulation by differential evolution with self-adaptive randomized subspace sampling." International Journal of Nonlinear Sciences and Numerical Simulation 10.3 (2009): 273-290.
#' @references ter  Braak C. J. F., and Vrugt J. A. (2008). Differential Evolution Markov Chain with snooker updater and fewer chains. Statistics and Computing http://dx.doi.org/10.1007/s11222-008-9104-9 
##' @details Insted of a bayesianSetup, the function can take the output of a previous run to restart the sampler 
##' from the last iteration. Due to the sampler's internal structure you can only use the output
##' of DREAMzs.
##' If you provide a matrix with start values the number of rows detemines the number of chains that are run.
##' The number of coloumns must be equivalent to the number of parameters in your bayesianSetup. \cr\cr
##' There are several small differences in the algorithm presented here compared to the original paper by Vrugt et al. (2009). Mainly
##' the algorithm implemented here does not have an automatic stopping criterion. Hence, it will
##' always run the number of iterations specified by the user. Also, convergence is not
##' monitored and left to the user. This can easily be done with coda::gelman.diag(chain).
##' Further the proposed delayed rejectio step in Vrugt et al. (2009) is not implemented here.\cr\cr
##' During the adaptation phase DREAM is running two mechanisms to enhance the sampler's efficiency.
##' First the disribution of crossover values is tuned to favor large jumps in the parameter space.
##' The crossover probabilities determine how many parameters are updated simultaneously. 
##' Second outlier chains are replanced as they can largely deteriorate the sampler's performance.
##' However, these steps destroy the detailed balance of the chain. Consequently these parts of the chain
##' should be discarded when summarizing posterior moments. This can be done automatically during the
##' sampling process (i.e. burnin > adaptation) or subsequently by the user. We chose to distinguish between
##' the burnin and adaptation phase to allow the user more flexibility in the sampler's settings.
##' 
##' @seealso \code{\link{DREAM}}
##' @export

DREAMzs <- function(bayesianSetup,
                    settings = list(iterations = 10000,
                                    nCR = 3,
                                    gamma = NULL, 
                                    eps = 0,
                                    e = 5e-2,
                                    pCRupdate = FALSE, 
                                    updateInterval = 10,
                                    burnin = 0,
                                    thin = 1,
                                    adaptation = 0.2,
                                    parallel = NULL,
                                    
                                    Z = NULL,
                                    ZupdateFrequency = 10,
                                    pSnooker = 0.1,
                                    
                                    
                                    DEpairs = 2, 
                                    consoleUpdates = 10, 
                                    startValue = NULL,
                                    currentChain = 1,
                                    message = FALSE)) {
  
  
  
  if("bayesianOutput" %in% class(bayesianSetup)){
    restart <- TRUE
  } else restart <- FALSE
  
  
  if(restart){
    if(is.null(settings)) settings <- bayesianSetup$settings
    else  settings <- applySettingsDefault(settings = settings, sampler = "DREAMzs")
    
    settings$adaptation <- 0 # set burnIn to 0 if restart because it has already been
    # applied in chain that is restarted and destroys detailed balance.
    
  }else{
    # If nothing provided use default settings
    settings <- applySettingsDefault(settings = settings, sampler = "DREAMzs")
  }
  
  if(!restart){ 
    setup <- bayesianSetup
  } else setup <- bayesianSetup$setup
  
  setup <- checkBayesianSetup(setup, parallel = settings$parallel) # calling parallel will check if requested parallelization in settings is provided by the BayesianSetup
  if(is.null(settings$parallel)) settings$parallel = setup$parallel # checking back - if no parallelization is provided, we use the parallelization in the BayesianSetup. We could also set parallel = F, but I feel it makes more sense to use the Bayesiansetup as default
  
  if(!restart){
    if(is.null(settings$startValue)){
      parLen = length(bayesianSetup$prior$sampler(1))
      X = bayesianSetup$prior$sampler(3)
    }
    if(is.function(settings$startValue)){
      X = settings$startValue()
    }
    if(class(settings$startValue)[1] == "numeric"){
      X = bayesianSetup$prior$sampler(settings$startValue)
    }
    
    if(is.matrix(settings$startValue)) X <- settings$startValue
    
    if(is.null(settings$Z)){
      parLen = length(bayesianSetup$prior$sampler(1))
      Z = bayesianSetup$prior$sampler(parLen * 10)
    } 
    if(is.function(settings$Z)){
      Z = settings$Z()
    }
    
    if(class(settings$Z)[1] == "numeric"){
      Z = bayesianSetup$prior$sampler(settings$Z)
    }
    if(is.matrix(settings$Z)) Z <- settings$Z
    
  }else{
    X <- bayesianSetup$X
    Z <- bayesianSetup$Z
    if(is.vector(Z)) Z = as.matrix(Z)
  }
  
  
  if (! is.matrix(X)) stop("wrong starting values")
  if (! is.matrix(Z)) stop("wrong Z values")
  
  
  FUN = setup$posterior$density
  
  pCRupdate <- settings$pCRupdate
  nCR <- settings$nCR
  Npar <- ncol(X)
  
  Npar12 <- (Npar - 1)/2  # factor for Metropolis ratio DE Snooker update
  
  parallel <- settings$parallel
  if(!is.null(parallel)){
  if(is.numeric(parallel) | parallel == "external") parallel <- TRUE
  }else parallel <- FALSE
  
  pCRupdate <- settings$pCRupdate
  nCR <- settings$nCR
  Npar <- ncol(X)
  Npop <- nrow(X)
  
  
  # Set adaptation if percentage is supplied
  if(settings$adaptation <1) settings$adaptation <- settings$adaptation*settings$iterations
  
  # Determine number of iterations and initialize chain
  n.iter <- ceiling(settings$iterations/Npop)
  if (n.iter < 2) stop ("The total number of iterations must be greater than 3")
  settings$burnin <- settings$burnin/Npop
  lChain <- ceiling((n.iter - settings$burnin)/settings$thin)+1
  pChain <- array(NA, dim=c(lChain, Npar+3, Npop))

  
  # assign memory for Z and write first values in Z
  M <- nrow(Z[complete.cases(Z),,drop = FALSE])
  Zold <- Z[complete.cases(Z),,drop = FALSE]
  Z <- matrix(NA, nrow= M + floor((n.iter) /settings$ZupdateFrequency) * Npop, ncol=Npar)
  Z[1:M,] <- Zold
  
  
  colnames(pChain) <- c(setup$names, "LP", "LL", "LPr")
  
  
  # Evaluate start values and write them in the chain
  logfitness_X <- FUN(X, returnAll = T)
  pChain[1,,] <- t(cbind(X,logfitness_X))
  
  
  # Set counter
  counter <- 1
  iseq <- 1:Npop
  
  
  #### gamma, initialization. However gamma is calculated every iteration (see below).
  gamma <- 2.38/sqrt(settings$DEpairs*Npar)
  
  
  ## delta initialization
  delta <- rep(0, settings$nCR)
  
  funevals <- 0
  #### pCR update
  # Initialization 
  if(!restart){
    pCR = rep(1/nCR, nCR)
    lCR <- rep(0,nCR)

    CR <- matrix(1/nCR, nrow = Npop, ncol = settings$updateInterval)
  }else{
    pCR <- bayesianSetup$pCR
    CR <- generateCRvalues(pCR, settings, Npop)
    
  }
  
  
  # helper counter for CR value index
  counter_update <- 0
  
  ##  Omega initialization
  omega <- numeric()
  
  ## eps and e
  eps <- settings$eps
  e <- settings$e
  
  
  ##################### Start iterations ##############################
  for(iter in 2:n.iter){
    
    xOld <- X
    
    
    if(parallel == TRUE){
      x_prop <- matrix(NA, nrow= Npop, ncol=Npar)
      r_extra <- numeric(Npop)
      
      for(i in 1:Npop){
        
        if(runif(1)>settings$pSnooker){
          selectedChains1 <- sample((1:M), settings$DEpairs, replace = FALSE)
          selectedChains2 <- numeric(settings$DEpairs)
          
          # Avoid that selected chains are identical
          for(k in 1:settings$DEpairs){
            selectedChains2[k] <- sample((1:M)[-c(selectedChains1[k],selectedChains2[1:k]) ],1)
          }
          
          # Get indices of parameters that are updated = indX
          rn <- runif(Npar)
          indX <- which(rn>(1-CR[i]))
          
          # Make sure at least one dimension is updated
          if(length(indX) == 0) indX <- sample(1:Npar, 1)
          
          # First update proposal
          x_prop[i,] <- X[i,]
          
          # Calculate gamma based on DEpairs and number of dimensions 
          # that are updated simulateously.
          # To jump between modes gamma is set to 1 every fifth iteration.
          if(runif(1)>4/5){ 
            gamma <- 1
          }else{
            gamma <-2.38/sqrt(settings$DEpairs* length(indX))
          }
          
          # No snooker update
          # Replace with new proposal for indX
          x_prop[i,indX] <- X[i,indX] + (1+e)*gamma*(apply(as.matrix(Z[selectedChains1,indX]),2,sum)- 
                                                     apply(as.matrix(Z[selectedChains2,indX]),2,sum)) + eps*rnorm(length(indX),0,1)
          r_extra[i] <- 0
          
        }else{ # Make proposal using snooker update
          selectSnooker <-  sample((1:M),replace = FALSE, 3)
          
          z <- Z[selectSnooker[1],]
          x_z <- X[i,] - z  
          D2 <- max(sum(x_z*x_z), 1.0e-300)
          projdiff <- sum((Z[selectSnooker[1],] -Z[selectSnooker[2],]) * x_z)/D2 # inner_product of difference with x_z / squared norm x_z
          gamma_snooker <- runif(1, min=1.2,max=2.2)
          
          x_prop[i,] <-  X[i,] + gamma_snooker * projdiff * x_z
          x_z <- x_prop[i,] - z
          D2prop <- max(sum(x_z*x_z), 1.0e-300)
          r_extra[i] <- Npar12 * (log(D2prop) - log(D2))
          
        }
      } # Npop
        
        
        # run proposals
        logfitness_x_prop <- FUN(x_prop, returnAll = T)
        
        # evaluate acceptance
        for(i in 1:Npop){
          if(!is.na(logfitness_x_prop[i,1] - logfitness_X[i,1])){
            if ((logfitness_x_prop[i,1] - logfitness_X[i,1] + r_extra[i]) > log(runif(1))){
              # accept <- accept + 1 
              X[i,] <- x_prop[i,]
              logfitness_X[i,] <- logfitness_x_prop[i,]
            }
          }
        }
        
    
      }else{   ## If not parallel
    for(i in 1:Npop){
      
      if(runif(1)>settings$pSnooker){
        selectedChains1 <- sample((1:M), settings$DEpairs, replace = FALSE)
        selectedChains2 <- numeric(settings$DEpairs)
        
        # Avoid that selected chains are identical
        for(k in 1:settings$DEpairs){
          selectedChains2[k] <- sample((1:M)[-c(selectedChains1[k],selectedChains2[1:k]) ],1)
        }
        
        # Get indices of parameters that are updated = indX
        rn <- runif(Npar)
        indX <- which(rn>(1-CR[i]))
        
        # Make sure at least one dimension is updated
        if(length(indX) == 0) indX <- sample(1:Npar, 1)
        
        # First update proposal
        x_prop <- X[i,]
      
        # Calculate gamma based on DEpairs and number of dimensions 
        # that are updated simulateously.
        # To jump between modes gamma is set to 1 every fifth iteration.
        if(runif(1)>4/5){ 
          gamma <- 1
        }else{
          gamma <-2.38/sqrt(settings$DEpairs* length(indX))
        }
        
        # No snooker update
        # Replace with new proposal for indX
        x_prop[indX] <- X[i,indX] + (1+e)*gamma*(apply(as.matrix(Z[selectedChains1,indX]),2,sum)- 
                                                   apply(as.matrix(Z[selectedChains2,indX]),2,sum)) + eps*rnorm(length(indX),0,1)
        r_extra <- 0
        
      }else{ # Make proposal using snooker update
        selectSnooker <-  sample((1:M),replace = FALSE, 3)
        
        z <- Z[selectSnooker[1],]
        x_z <- X[i,] - z  
        D2 <- max(sum(x_z*x_z), 1.0e-300)
        projdiff <- sum((Z[selectSnooker[1],] -Z[selectSnooker[2],]) * x_z)/D2 # inner_product of difference with x_z / squared norm x_z
        gamma_snooker <- runif(1, min=1.2,max=2.2)
        
        x_prop <-  X[i,] + gamma_snooker * projdiff * x_z
        x_z <- x_prop - z
        D2prop <- max(sum(x_z*x_z), 1.0e-300)
        r_extra <- Npar12 * (log(D2prop) - log(D2))
        
      }
      
      
      logfitness_x_prop <- FUN(x_prop, returnAll = T)
      
      if(!is.na(logfitness_x_prop[1] - logfitness_X[i,1])){ # To catch possible error
        if ((logfitness_x_prop[1] - logfitness_X[i,1] + r_extra) > log(runif(1))){
          X[i,] <- x_prop
          logfitness_X[i,] <- logfitness_x_prop
        }
        
      }
      
    } #Npop
    
      } # not parallel
    
    
    ## Write values in chain
    
    if((iter > settings$burnin) && (iter %% settings$thin == 0)){
      counter <- counter+1
      pChain[counter,,] <- t(cbind(X,logfitness_X))
    }
    
    # Update Z
    if(counter%%settings$ZupdateFrequency == 0){
      Z[(M+1):(M+Npop),] <- X
      M <- M+Npop
      
    }
    
    ###################################
    
    if(iter < settings$adaptation){
      
      if(pCRupdate){  ## Calculate delta, this is (unlike the update) done every iteration
        ## Calculate delta
        
        ## Calculate standard deviation of each dimension of X
        sdX <- apply(X[,1:Npar,drop=FALSE],2,sd)
        
        ## Compute Euclidean distance between old and new X values
        delta_Norm <- rowSums(((xOld-X[,1:Npar,drop=FALSE])/sdX)^2)
        
        ## Now delta can be calculated
        for (k in 1:settings$nCR){ # Loop over CR values
          
          # Find updated chains         
          ind <- which(abs(CR[,k]-(k/nCR)) < 1e-5)
          
          ## Add normalized squared distance to the current delta
          delta[k] <- delta[k]+sum(delta_Norm[ind])
          #delta[k] <- delta[k]+sum(delta_Norm)
          
        } 
        
      }
      
      
      if(iter%%settings$updateInterval == 0){
        
        
        if(pCRupdate){
          # Update CR values
          tmp <- AdaptpCR(CR, delta, lCR, settings, Npop)
          pCR <- tmp$pCR
          lCR <- tmp$lCR
        }
        
        ## remove outliers
        ## TODO include if(remOutliers = TRUE) ??
        for(out in 1:Npop){
          omega[out] <- mean(pChain[((counter/2):counter),Npar+1, out])
        }
        
        if(NaN %in% omega){
          outlierChain <- NULL # Prevent possible error
        }else{
          # Inter quantile range
          IQR <- quantile(omega, probs = c(0.25, 0.75))
          
          # Determine outlier chains
          outlierChain <- which(omega< IQR[1] - 2*(IQR[2]-IQR[1]))
        }
        
        
        # Replace with best chain
        if(length(outlierChain) > 0){
          best <- which.max(pChain[counter,Npar+1,])
          pChain[counter,,outlierChain] <- pChain[counter,,best]
          
        } # Remove outliers
        
      }
    }
    
    if(iter%%settings$updateInterval == 0){
      counter_update <- 0 # set counter back to zero
      CR <- generateCRvalues(pCR, settings, Npop)  
      
    }
    ###############################
    
    
    
    if(settings$message){  
      if( (iter %% settings$consoleUpdates == 0) | (iter == n.iter)) cat("\r","Running DREAM-MCMC, chain ", settings$currentChain, 
                                                                         "iteration" ,iter*Npop,"of",n.iter*Npop,". Current logp ",
                                                                         logfitness_X[,1],
                                                                         "Please wait!","\r")
      flush.console()
    }
    
    
  } # niter
  
  iterationsOld <- 0
  
  pChain <- pChain[1:counter,,]
  
  if(restart){ # Combine chains
    newchains <- array(NA, dim = c((counter+nrow(bayesianSetup$chain[[1]])), (Npar+3), Npop))
    
    for(i in 1:Npop){
      for(k in 1:(Npar+3)){
        newchains[,k,i] <- c(bayesianSetup$chain[[i]][,k],pChain[,k,i])
      }
    }
    pChain <- newchains
  }
  
  
  
  pChain<- coda::as.mcmc.list(lapply(1:Npop,function(i) coda::as.mcmc(pChain[,1:(Npar+3),i])))
  
  
  list(chains = pChain, X = as.matrix(X[,1:Npar]), Z = Z, pCR = pCR)

}
