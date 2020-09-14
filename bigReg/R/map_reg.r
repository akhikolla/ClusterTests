
#------------------------------------------------------------------------------------------------
# Function specification for blm, bglm and supporting functions
#------------------------------------------------------------------------------------------------

# Function for control parameters
#' @title Function for creating control parameters for the GLM fit
#' @description Function for creating control parameters for the GLM fit
#' @param epsilon defaults to 1E-8
#' @param maxit defaults 25 maximum number of iterations
#' @param trace defaults to TRUE
#' @export
.control <- function (epsilon = 1E-8, maxit = 25, trace = TRUE) 
{
  if (!is.numeric(epsilon) || epsilon <= 0) 
    stop("value of 'epsilon' must be > 0")
  if (!is.numeric(maxit) || maxit <= 0) 
    stop("maximum number of iterations must be > 0")
  list(epsilon = epsilon, maxit = maxit, trace = trace)
}

#------------------------------------------------------------------------------

#' @title Function to carry out linear regression on a data_frame data object
#' @description Function to carry out linear regression on a data_frame data object
#' @param formula formula that defines your regression model
#' @param data data_frame object containing data for linear regression
#' @param weights weights for the model
#' @param offset offsets for the model
#' @param control list of parameters for control() function
#' @export
blm <- function(formula = stop("formula: not supplied"), data = stop("data: data not supplied"), 
					control = list(), weights = NULL, offset = NULL)
{
	mCall <- match.call()
	weights <- substitute(weights)
	offset <- substitute(offset)
  
	if(class(weights) == "character")
	{
		weights <- as.call(parse(text = weights))
		weights <- as.list(weights)[[1]]
	}
  
	if(class(offset) == "character")
	{
		offset <- as.call(parse(text = offset))
		offset <- as.list(offset)[[1]]
	}
  
	# The model terms
	mterms <- terms(formula)
	dataCols <- as.list(attr(mterms, "variables"))[-1]
	dataCols <- sapply(dataCols, function(x){
		if(is.call(x))
		{
			x <- as.list(x)[-1]
			x <- x[[which(sapply(x, class) == "name")]]
		}
		return(deparse(x))
	})
  
	xnames <- dimnames(model.matrix(formula, 
					data = data$Reader(data$block_names[1])[1,]))[[2]]

	# Offsets
	if(!is.null(offset)){
		if(is.call(offset)){
			osfCol <- deparse(as.list(offset)[[length(offset)]])
		}else{
			osfCol <- deparse(offset)
		}
		dataCols <- c(dataCols, osfCol)
	}else{
		offset <- NULL
	}
  
	# Weights
	if(!is.null(weights)){
		if(is.call(weights)){
			weightCol <- deparse(as.list(weights)[[length(weights)]])
		}else{
			weightCol <- deparse(weights)
		}
		dataCols <- c(dataCols, weightCol)
	}else{
		weights <- NULL
	}
  
	control <- do.call(.control, control)
  
	# For loop to process function
	# block_names <- data$block_names
	variables <- attr(mterms, "variables")
	varNames <- as.character(variables[-1])
  
	mmCall <- as.list(mCall)
	mmCall <- mmCall[!names(mmCall) %in% c("control", "start", "family")]
	mmCall[[1L]] <- as.name("model.frame")
	mmCall[["offset"]] <- offset
	mmCall[["weights"]] <- weights
	mmCall <- as.call(mmCall)
	mmCall[["data"]] <- quote(.mf)

	# Function to be mapped over the blocks
	process_blm_block <- function(.mf, .call, .formula){
		.mf <- eval(.call)
		X <- model.matrix(.formula, .mf)
		y <- matrix(model.response(.mf))
		if(!is.null(offset))
		{
			y <- y - model.offset(.mf)
		}
	  
		if(is.null(weights))
		{
			out <- XWXMatrix(X, y)
		}else{
			W <- model.weights(.mf)
			out <- XWXMatrixW(X, y, W)
		}

		return(out)
	}

	sum_blm_block <- function(x1, x2){
		x1[["XWX"]] = x1[["XWX"]] + x2[["XWX"]]
	    x1[["XWy"]] = x1[["XWy"]] + x2[["XWy"]]
	    x1[["yWy"]] = x1[["yWy"]] + x2[["yWy"]]
	    x1[["nobs"]] = x1[["nobs"]] + x2[["nobs"]]
	    x1[["wt"]] = x1[["wt"]] + x2[["wt"]]
	    x1[["ySum"]] = x1[["ySum"]] + x2[["ySum"]]
	    return(x1)
	}

	fit <- Reduce(sum_blm_block, data$map(process_blm_block, .call = mmCall, .formula = formula))

	# Coefficient Calculation
	fit <- SVD(fit, control$epsilon)
	fit$call <- mCall
	fitNames <- c("call", "coefficients", "rank", "var", "cov", "deviance", "nobs",
					"aliased", "cov.unscaled", "XWX", "XWy", "yWy", "wt", "ySum")
	fit <- fit[fitNames]
	dimnames(fit$coefficients)[[1]] <- xnames
	dimnames(fit$cov) <- list(xnames, xnames)
	dimnames(fit$aliased)[[1]] <- xnames
	dimnames(fit$cov.unscaled) <- list(xnames, xnames)
	dimnames(fit$XWX) <- list(xnames, xnames)
	dimnames(fit$XWy)[[1]] <- xnames
	fit$NMissing <- data$nrows - fit$nobs
	class(fit) <- "blm"
	return(fit)
}

#------------------------------------------------------------------------------

#' @title print function for the blm object
#' @description print function for the blm object
#' @param x blm object to be displayed
#' @param digits number of significant digits to use
#' @param ... not yet used
#' @export
print.blm <- function(x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n",
      paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if(length(drop(x$coefficients))) {
    cat("Coefficients:\n")
    print.default(format(drop(x$coefficients), digits = digits),
                  print.gap = 2L, quote = FALSE)
  } else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}

#------------------------------------------------------------------------------

#' @title summary function for the blm object
#' @description summary function for the blm object
#' @param object blm object to be summarized
#' @param ... not used
#' @export
summary.blm <- function(object, ...){
  ans <- list()
  rdf <- object$nobs - object$rank
  N <- object$nobs
  coeffs <- object$coefficients
  coeffs <- cbind(coeffs, sqrt(diag(object$cov)))
  coeffs <- cbind(coeffs, coeffs[,1]/coeffs[,2])
  coeffs <- cbind(coeffs, 2 * pt(abs(coeffs[,3]), rdf, lower.tail = FALSE))
  dimnames(coeffs)[[2]] <- c("Estimate", "Std.Error", "t.value", "P(>|t|)")
  ans$call <- object$call
  ans$coefficients <- coeffs
  ans$logLik <- list("logLik" = .5* (object$wt - N * (log(2 * pi) + 1 - 
                  log(N) + log(object$deviance))), "df" = object$rank + 1)
  ans$AIC <- -2*ans$logLik$logLik + 2*(object$rank)
  ans$BIC <- -2*ans$logLik$logLik + 2*(object$rank)*log(N)
  N.y2 <- ((object$ySum)^2/N)
  S0 <- drop(object$yWy - N.y2)
  S1 <- drop(object$yWy - t(object$coefficients) %*% object$XWy)
  ans$r.squared <- (S0 - S1)/S0
  ans$f.statistic <- matrix(c(((S0 - S1)/(object$rank - 1))/(S1/(N - object$rank)), 
                       object$rank - 1, N - object$rank), nrow = 1)
  fpval <- pf(ans$f.statistic[,1], ans$f.statistic[,2], ans$f.statistic[,3], 
              lower.tail = FALSE)
  ans$f.statistic <- cbind(ans$f.statistic, fpval)
  dimnames(ans$f.statistic)[[2]] <- c("value", "numdf", "dendf", "pval")
  ans$NMissing <- object$NMissing
  class(ans) <- "summary.blm"
  ans
}

#------------------------------------------------------------------------------

#' @title Function to print the summary object from the blm object
#' @description Function to print the summary object from the blm object
#' @param x summary blm object
#' @param digits - the digits to be displayed
#' @param signif.stars passed to printCoefmat
#' @param ... arguments passed to \code{printCoefmat()} function
#' @export
print.summary.blm <- function(x, digits = max(3L, getOption("digits") - 3L),
                        signif.stars = getOption("show.signif.stars"), ...)
{
  cat("\nCall:\n", paste(deparse(x$call), sep="\n", collapse = "\n"), "\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
               na.print = "NA", ...)
  cat("\nlogLik:", x$logLik$logLik, ", df:", x$logLik$df)
  cat("\nAIC:", x$AIC, ", BIC:", x$BIC)
  cat("\nR-Squared:", x$r.squared)
  cat("\nF-Statistic:", x$f.statistic[1L], "on", x$f.statistic[2L], "and",
      x$f.statistic[3L], "DF, pvalue:", x$f.statistic[4L], "\n\n")
  invisible(x)
}

#------------------------------------------------------------------------------
#' @title Function to print the summary object from the blm object
#' @description Function to print the summary object from the blm object
#' @param mf the data block to be processed
#' @param formula the formula of for the model
#' @param mmCall the call object of the model
#' @param family the family object for the model
#' @param offset the model offset
#' @param weights the model weights
#' @param start the starting coefficient estimates
#' @param niter the current number of iterations
#' @param etastart the start for eta
#' @param mustart the start for mu
process_bglm_block <- function(mf, formula, mmCall, family, offset, weights, start, 
					niter, etastart, mustart){
  
  mf <- eval(mmCall)
  X <- model.matrix(formula, mf)
  y <- matrix(model.response(mf))
  if(!is.null(offset))
  {
    ofst <- model.offset(mf)
  }else{
    ofst <- 0
  }
  
  if(is.null(weights))
  {
    wghts <- 1
  }else{
    wghts <- drop(model.weights(mf))
  }
  
  if(niter != 0)
  {
    eta <- X %*% start + ofst
    mu <- family$linkinv(eta)
    
  }else{
    
    eval(family$initialize)
    
    if(!is.null(mustart))
    {
      mu <- mustart
      eta <- family$linkfun(mu)
    }else{
      eta <- X %*% start + ofst
      mu <- family$linkinv(eta)
    }
    
    if(!is.null(etastart))
    {
      eta <- etastart
    }
  }
  
  mu.eta.val <- family$mu.eta(eta)
  
  W <- drop((wghts * mu.eta.val^2)/family$variance(mu))
  residuals <- (y - mu)/mu.eta.val
  z <- (eta - ofst) + residuals
  
  out <- XWXMatrixW(X, z, W)
  out$dev <- sum(family$dev.resids(y, mu, wghts))
  out$aic <- family$partialAIC(y, out$nobs, mu, wghts, out$dev)
  out$resid.sqr <- sum(W*residuals^2)
  out$col_sums <- colSums(X*wghts)
  return(out)
}

#------------------------------------------------------------------------------

#' @title The reduction function for the algorithm
#' @description The reduction function for the algorithm
#' @param x1 the first list object to be reduced
#' @param x2 the second list object to be reduced
sum_bglm_block <- function(x1, x2){
    x1[["XWX"]] = x1[["XWX"]] + x2[["XWX"]]
    x1[["XWy"]] = x1[["XWy"]] + x2[["XWy"]]
    x1[["yWy"]] = x1[["yWy"]] + x2[["yWy"]]
    x1[["nobs"]] = x1[["nobs"]] + x2[["nobs"]]
    x1[["dev"]] = x1[["dev"]] + x2[["dev"]]
    x1[["resid.sqr"]] = x1[["resid.sqr"]] + x2[["resid.sqr"]]
    x1[["ySum"]] = x1[["ySum"]] + x2[["ySum"]]
    x1[["aic"]] = x1[["aic"]] + x2[["aic"]]
    x1[["col_sums"]] = x1[["col_sums"]] + x2[["col_sums"]]
    return(x1)
}

#------------------------------------------------------------------------------
# Function to carry out the iterative fit
#------------------------------------------------------------------------------
.loopFit <- function(formula, data, mmCall, family, offset, weights, start, 
					niter, etastart, mustart){

	fit <- Reduce(sum_bglm_block, data$map(process_bglm_block, formula = formula, mmCall = mmCall, 
        family = family, offset = offset, weights = weights,
        start = start, niter = niter, etastart = etastart,
        mustart = mustart))

	return(fit)
}

#------------------------------------------------------------------------------
#' @title Function to carry out generalized linear regression on a data_frame data object
#' @description Function to carry out generalized linear regression on a data_frame data object
#' @param formula formula that defines your regression model
#' @param family family object from activeReg, e.g. .gaussian(), .binomial(), .poisson(),
#'                .quasipoisson(), .quasibinomial(), .Gamma(), .inverse.gaussian(),
#'                .quasi()
#' @param data data_frame object containing data for linear regression
#' @param weights weights for the model
#' @param offset offsets for the model
#' @param start starting values for the linear predictor
#' @param control list of parameters for .control() function
#' @param etastart starting values for the linear predictor
#' @param mustart starting values for vector of means
#' @examples
#' require(parallel)
#' data("plasma", package = "bigReg")
#' plasma1 <- plasma
#' plasma1 <- data_frame(plasma1, 10, path = "outputs", nCores = 1)
#' plasma_glm <- bglm(ESR ~ fibrinogen + globulin, data = plasma1, family = binomial_("logit"))
#' summary(plasma_glm)
#' @export
bglm <- function(formula, family = gaussian_(), data, weights = NULL, offset = NULL, 
                 start = NULL, control = list(), etastart = NULL, mustart = NULL)
{
  
  mCall <- match.call()
  weights <- substitute(weights)
  offset <- substitute(offset)
  
  if(class(weights) == "character")
  {
    weights <- as.call(parse(text = weights))
    weights <- as.list(weights)[[1]]
  }
  
  if(class(offset) == "character")
  {
    offset <- as.call(parse(text = offset))
    offset <- as.list(offset)[[1]]
  }
  
  # The model terms
  mterms <- terms(formula)
  dataCols <- as.list(attr(mterms, "variables"))[-1]
  dataCols <- sapply(dataCols, function(x){
    if(is.call(x))
    {
      x <- as.list(x)[-1]
      x <- x[[which(sapply(x, class) == "name")]]
    }
    return(deparse(x))
  })
  
  xnames <- dimnames(model.matrix(formula, data = data$Reader(data$block_names[1])[1,]))[[2]]
  
  nvars <- length(xnames)
  
  if(is.null(start))
    start <- rep(0, nvars)
  
  coefold <- start
  # coefold[] <- 0
  
  conv <- FALSE
  niter <- 0
  
  env <- environment()
  
  # Offsets
  if(!is.null(offset)){
    if(is.call(offset)){
      osfCol <- deparse(as.list(offset)[[length(offset)]])
    }else{
      osfCol <- deparse(offset)
    }
    dataCols <- c(dataCols, osfCol)
  }else{
    offset <- NULL
  }
  
  # Weights
  if(!is.null(weights)){
    if(is.call(weights)){
      weightCol <- deparse(as.list(weights)[[length(weights)]])
    }else{
      weightCol <- deparse(weights)
    }
    dataCols <- c(dataCols, weightCol)
  }else{
    weights <- NULL
  }
  
  control <- do.call(.control, control)
  
  # For loop to process function
  variables <- attr(mterms, "variables")
  
  # Call for constructing model matrix
  mmCall <- as.list(mCall)
  mmCall <- mmCall[!names(mmCall) %in% c("control", "start", "family")]
  mmCall[[1L]] <- as.name("model.frame")
  mmCall[["offset"]] <- offset
  mmCall[["weights"]] <- weights
  mmCall <- as.call(mmCall)
  mmCall[["data"]] <- quote(mf)
  
  coeffDiff <- 1
  
  timeA <- proc.time()
  # Convergence criterion
  while(all(start == 0) | coeffDiff > control$epsilon)
  {
    if(niter > control$maxit)
    {
      cat("The process has not converged, number of iterations is greater than maxit\n")
      break
    }
    
    fit <- .loopFit(formula = formula, data = data, mmCall = mmCall, family = family, 
                offset = offset, weights = weights, start = start, niter = niter, 
                etastart = etastart, mustart = mustart)
    # Coeff Calculation
    fit <- SVD(fit, control$epsilon)
    coefold <- start
    start <- fit$coefficients
    coeffDiff <- sum(abs(coefold - start))
    # A deviance approximation
    ndev <- fit$deviance
    
    if (!is.finite(ndev))
    {
      if (is.null(coefold))
        stop("no valid set of coefficients has been found: please supply starting values", 
             call. = FALSE)
      warning("step size truncated due to divergence", call. = FALSE)
      
      ii <- 1
      while (!is.finite(ndev)) {
        if (ii > control$maxit) 
          stop("inner loop 1; cannot correct step size", 
               call. = FALSE)
        if (control$trace) 
          cat("Step halved \n")
        ii <- ii + 1
        start <- (start + coefold)/2
        fit <- .loopFit(formula = formula, data = data, mmCall = mmCall, family = family, 
                    offset = offset, weights = weights, start = start, niter = niter, 
                    etastart = etastart, mustart = mustart)
        fit <- SVD(fit, control$epsilon)
        ndev <- fit$deviance
      }
    }
    
    niter <- niter + 1
  }
  timeB <- proc.time()
  
  fit$call <- mCall
  fitNames <- c("call", "coefficients", "col_sums", "rank", "var", "cov", "dev", "nobs", "aliased",
                "aic", "cov.unscaled", "resid.sqr", "XWX", "XWy", "yWy", "wt", "ySum")
  fit <- fit[fitNames]
  fit$aic <- family$finalAIC(fit$aic, fit$nobs, fit$dev) + 2*fit$rank
  
  if(family$family %in% c("binomial", "poisson"))
  {
    est.disp <- FALSE
    fit$dispersion <- 1
  }else{
    est.disp <- TRUE
    fit$dispersion <- fit$resid.sqr/(fit$nobs - fit$rank)
  }
  
  fit$formula <- formula
  fit$terms <- mterms
  fit$cov <- fit$dispersion * fit$cov.unscaled
  dimnames(fit$coefficients)[[1]] <- xnames
  dimnames(fit$cov) <- list(xnames, xnames)
  dimnames(fit$aliased)[[1]] <- xnames
  dimnames(fit$cov.unscaled) <- list(xnames, xnames)
  dimnames(fit$XWX) <- list(xnames, xnames)
  dimnames(fit$XWy)[[1]] <- xnames
  fit$est.disp <- est.disp
  fit$NMissing <- data$nrows - fit$nobs
  fit$col_sums <- matrix(fit$col_sums, nrow = 1)
  dimnames(fit$col_sums)[[2]] <- xnames
  fit$col_means <- fit$col_sums/fit$nobs
  fit$niter <- niter
  fit$solver_time <- (timeB - timeA)["elapsed"]
  names(fit$solver_time) <- NULL
  class(fit) <- "bglm"
  return(fit)
}

#------------------------------------------------------------------------------

#' @title print function for the bglm object
#' @description print function for the bglm object
#' @param x bglm object to be displayed
#' @param digits number of significant digits to use
#' @param ... not yet used
#' @export
print.bglm <- print.blm

#------------------------------------------------------------------------------

#' @title summary function for the bglm object
#' @description summary function for the bglm object
#' @param object bglm object to be summarized
#' @param ... not used
#' @export
summary.bglm <- function(object, ...){
  ans <- list()
  rdf <- object$nobs - object$rank
  N <- object$nobs
  coeffs <- object$coefficients
  coeffs <- cbind(coeffs, sqrt(diag(object$cov)))
  coeffs <- cbind(coeffs, coeffs[,1]/coeffs[,2])
  if(!object$est.disp)
  {
    coeffs <- cbind(coeffs, 2 * pnorm(-abs(coeffs[,3])))
    dimnames(coeffs)[[2]] <- c("Estimate", "Std.Error", "z.value", "P(>|z|)")
  }else{
    coeffs <- cbind(coeffs, 2 * pt(-abs(coeffs[,3]), rdf))
    dimnames(coeffs)[[2]] <- c("Estimate", "Std.Error", "t.value", "P(>|t|)")
  }
  
  ans$call <- object$call
  ans$coefficients <- coeffs
  ans$logLik <- list("logLik" = .5* (object$wt - 
      N * (log(2 * pi) + 1 - log(N) + log(object$dev))), "df" = object$rank + 1)
  ans$AIC <- ifelse(is.na(object$aic), -2*ans$logLik$logLik + 2*(object$rank), object$aic)
  ans$BIC <- -2*ans$logLik$logLik + 2*(object$rank)*log(N)
  ans$NMissing <- object$NMissing
  class(ans) <- "summary.bglm"
  ans
}

#------------------------------------------------------------------------------

#' @title Function to print the summary object from the bglm object
#' @description Function to print the summary object from the bglm object
#' @param x summary blm object
#' @param digits - the digits to be displayed
#' @param signif.stars passed to printCoefmat
#' @param ... arguments passed to \code{printCoefmat()} function
#' @export
print.summary.bglm <- function(x, digits = max(3L, getOption("digits") - 3L),
                               signif.stars = getOption("show.signif.stars"), ...)
{
  cat("\nCall:\n", paste(deparse(x$call), sep="\n", collapse = "\n"), "\n")
  printCoefmat(x$coefficients, digits = digits, signif.stars = signif.stars,
               na.print = "NA", ...)
  cat("\nlogLik:", x$logLik$logLik, ", df:", x$logLik$df)
  cat("\nAIC:", x$AIC, ", BIC:", x$BIC, "\n\n")
  invisible(x)
}

# Predict bglm
#' @title predict function for bglm object
#' @description predict function for bglm object
#' @param mf model frame
#' @param object a bglm object
#' @param type one of c("link", "response", "terms")
#' @export
bglm_predict <- function(mf = stop("mf: model frame must be supplied"), 
                          object = stop("object: bglm object must be supplied"),
                          type = stop("type: either \"link\", \"response\", \"terms\"")){
  mf <- model.frame(object$formula, data = mf)
  mm <- model.matrix(object$terms, data = mf)
  betas <- object$coefficients
  if(type == "terms"){
    aa <- attr(mm, "assign")
    tt <- object$terms
    ll <- attr(tt, "term.labels")
    mf <- mf[,-attr(tt, "response")]
    hasintercept <- attr(tt, "intercept") > 0L
    ll <- c(rep("(intercept)", hasintercept), ll)
    mm <- sweep(mm, 2, object$col_means)
    asgn <- split(order(aa), aa)
    names(asgn) <- ll
    mfDim <- dim(mf)
    predictors <- matrix(NA, nrow = mfDim[1], ncol = length(ll))
    dimnames(predictors)[[2]] <- ll #names(mf)
    if(hasintercept){
      asgn["(intercept)"] <- NULL
      predictors <- predictors[,-1, drop = FALSE]
      ll <- ll[-1]
      termsconst <- object$col_means %*% object$coefficients
    }
    for(i in ll){
      predictors[,i] <- mm[,asgn[[i]], drop = FALSE] %*% betas[asgn[[i]], , drop = FALSE]
    }
    attr(predictors, "constant") <- if (hasintercept) 
        termsconst
    else 0
  }else if(type %in% c("link", "response")){
    predictors <- mm %*% betas
    if(type == "response")
      predictors <- eval(object$call$family)$linkinv(predictors)
  }

  return(predictors)
}

