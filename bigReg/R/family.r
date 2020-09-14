#------------------------------------------------------------------------------------------------
# Family functions for blm and bglm big data functions
#------------------------------------------------------------------------------------------------

# The binomial family
#' @title binomial family function
#' @description binomial family function
#' @param link function character
#' @export
binomial_ <- function (link = "logit") 
{
  linktemp <- substitute(link)
  if (!is.character(linktemp)) 
    linktemp <- deparse(linktemp)
  okLinks <- c("logit", "probit", "cloglog", "cauchit", "log")
  if (linktemp %in% okLinks) 
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  }
  else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name)) 
        linktemp <- stats$name
    }
    else {
      stop(gettextf("link \"%s\" not available for binomial family; available links are %s", 
                    linktemp, paste(sQuote(okLinks), collapse = ", ")), 
           domain = NA)
    }
  }
  variance <- function(mu) mu * (1 - mu)
  validmu <- function(mu) all(mu > 0) && all(mu < 1)
  dev.resids <- binomial()$dev.resids # function(y, mu, wt) .Call(stats:::C_binomial_dev_resids, y, mu, wt)
  partialAIC <- function(y, n, mu, wt, dev)NA
  finalAIC <- function(paic, n, dev)NA
  initialize <- expression({
    if (NCOL(y) == 1) {
      if (is.factor(y)) y <- y != levels(y)[1L]
      #n <- rep.int(1, nrow(y))
      y[wghts == 0] <- 0
      if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
      mustart <- (wghts * y + 0.5)/(wghts + 1)
      m <- wghts * y
      if (any(abs(m - round(m)) > 0.001)) warning("non-integer #successes in a binomial glm!")
    } 
  })
  simfun <- function(object, nsim) {
    ftd <- fitted(object)
    n <- length(ftd)
    ntot <- n * nsim
    wts <- object$prior.weights
    if (any(wts%%1 != 0)) 
      stop("cannot simulate from non-integer prior.weights")
    if (!is.null(m <- object$model)) {
      y <- model.response(m)
      if (is.factor(y)) {
        yy <- factor(1 + rbinom(ntot, size = 1, prob = ftd), 
                     labels = levels(y))
        split(yy, rep(seq_len(nsim), each = n))
      }
      else if (is.matrix(y) && ncol(y) == 2) {
        yy <- vector("list", nsim)
        for (i in seq_len(nsim)) {
          Y <- rbinom(n, size = wts, prob = ftd)
          YY <- cbind(Y, wts - Y)
          colnames(YY) <- colnames(y)
          yy[[i]] <- YY
        }
        yy
      }
      else rbinom(ntot, size = wts, prob = ftd)/wts
    }
    else rbinom(ntot, size = wts, prob = ftd)/wts
  }
  structure(list(family = "binomial", link = linktemp, linkfun = stats$linkfun, 
                 linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, 
                 partialAIC = partialAIC, finalAIC = finalAIC, mu.eta = stats$mu.eta, 
                 initialize = initialize, validmu = validmu, valideta = stats$valideta, 
                 simulate = simfun), class = "family")
}

# The poisson family
#' @title poisson family function
#' @description poisson family function
#' @param link function character
#' @export
poisson_ <- function (link = "log")
{
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("log", "identity", "sqrt")
  if (linktemp %in% okLinks)
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    ## what else shall we allow?  At least objects of class link-glm.
    if(inherits(link, "link-glm")) {
      stats <- link
      if(!is.null(stats$name)) linktemp <- stats$name
    } else {
      stop(gettextf('link "%s" not available for poisson family; available links are %s',
                    linktemp, paste(sQuote(okLinks), collapse =", ")),
           domain = NA)
    }
  }
  variance <- function(mu) mu
  validmu <- function(mu) all(mu>0)
  dev.resids <- function(y, mu, wt)
    2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))
  partialAIC <- function(y, n, mu, wt, dev) -2*sum(dpois(y, mu, log=TRUE)*wt)
  finalAIC <- function(paic, n, dev){
    paic
  }
  initialize <- expression({
    if (any(y < 0))
      stop("negative values not allowed for the 'Poisson' family")
    mustart <- y + 0.1
  })
  simfun <- function(object, nsim) {
    ## A Poisson GLM has dispersion fixed at 1, so prior weights
    ## do not have a simple unambiguous interpretation:
    ## they might be frequency weights or indicate averages.
    wts <- object$prior.weights
    if (any(wts != 1)) warning("ignoring prior weights")
    ftd <- fitted(object)
    rpois(nsim*length(ftd), ftd)
  }
  structure(list(family = "poisson",
                 link = linktemp,
                 linkfun = stats$linkfun,
                 linkinv = stats$linkinv,
                 variance = variance,
                 dev.resids = dev.resids,
                 partialAIC = partialAIC, finalAIC = finalAIC,
                 mu.eta = stats$mu.eta,
                 initialize = initialize,
                 validmu = validmu,
                 valideta = stats$valideta,
                 simulate = simfun),
            class = "family")
}

# Gaussian family
#' @title gaussian family function
#' @description gaussian family function
#' @param link function character
#' @export
gaussian_ <- function (link = "identity")
{
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("inverse", "log", "identity")
  if (linktemp %in% okLinks)
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    ## what else shall we allow?  At least objects of class link-glm.
    if(inherits(link, "link-glm")) {
      stats <- link
      if(!is.null(stats$name)) linktemp <- stats$name
    } else {
      stop(gettextf('link "%s" not available for gaussian family; available links are %s',
                    linktemp, paste(sQuote(okLinks), collapse =", ")),
           domain = NA)
    }
  }
  partialAIC <-  function(y, n, mu, wt, dev) {
    sum(log(wt))
  }
  finalAIC <- function(paic, n, dev){
    n*(log(dev/n*2*pi) + 1) + 2 - paic
  }
  structure(list(family = "gaussian",
                 link = linktemp,
                 linkfun = stats$linkfun,
                 linkinv = stats$linkinv,
                 variance = function(mu) rep.int(1, length(mu)),
                 dev.resids = function(y, mu, wt) wt * ((y - mu)^2),
                 partialAIC = partialAIC, finalAIC = finalAIC,
                 mu.eta = stats$mu.eta,
                 initialize = expression({
                   n <- rep.int(1, nrow(y))
                   if(is.null(etastart) && is.null(start) &&
                        is.null(mustart) &&
                        ((family$link == "inverse" && any(y == 0)) ||
                           (family$link == "log" && any(y <= 0))))
                     stop("cannot find valid starting values: please specify some")
                   
                   mustart <- y 
                   }),
                 validmu = function(mu) TRUE,
                 valideta = stats$valideta
  ),
            class = "family")
}


#' @title quasipoisson family function
#' @description quasipoisson family function
#' @param link function character
#' @export
quasipoisson_ <- function (link = "log")
{
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("log", "identity", "sqrt")
  if (linktemp %in% okLinks)
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    ## what else shall we allow?  At least objects of class link-glm.
    if(inherits(link, "link-glm")) {
      stats <- link
      if(!is.null(stats$name)) linktemp <- stats$name
    } else {
      stop(gettextf('link "%s" not available for quasipoisson family; available links are %s',
                    linktemp, paste(sQuote(okLinks), collapse =", ")),
           domain = NA)
    }
  }
  variance <- function(mu) mu
  validmu <- function(mu) all(mu>0)
  dev.resids <- function(y, mu, wt)
    2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))
  partialAIC <- function(y, n, mu, wt, dev) NA
  finalAIC <- function(paic, n, dev) NA
  initialize <- expression({
    if (any(y < 0))
      stop("negative values not allowed for the 'Poisson' family")
    mustart <- y + 0.1
  })
  structure(list(family = "quasipoisson",
                 link = linktemp,
                 linkfun = stats$linkfun,
                 linkinv = stats$linkinv,
                 variance = variance,
                 dev.resids = dev.resids,
                 partialAIC = partialAIC, finalAIC = finalAIC,
                 mu.eta = stats$mu.eta,
                 initialize = initialize,
                 validmu = validmu,
                 valideta = stats$valideta),
            class = "family")
}

#' @title quasibinomial family function
#' @description quasibinomial family function
#' @param link function character
#' @export
quasibinomial_ <- function (link = "logit")
{
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("logit", "probit", "cloglog", "cauchit", "log")
  if (linktemp %in% okLinks)
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    ## what else shall we allow?  At least objects of class link-glm.
    if(inherits(link, "link-glm")) {
      stats <- link
      if(!is.null(stats$name)) linktemp <- stats$name
    } else {
      stop(gettextf('link "%s" not available for quasibinomial family; available links are %s',
                    linktemp, paste(sQuote(okLinks), collapse =", ")),
           domain = NA)
    }
  }
  variance <- function(mu) mu * (1 - mu)
  validmu <- function(mu) all(mu>0) && all(mu<1)
  dev.resids <- function(y, mu, wt)
    2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) +
                (1 - y) * log(ifelse(y == 1, 1, (1 - y)/(1 - mu))))
  partialAIC <- function(y, n, mu, wt, dev) NA
  finalAIC <- function(paic, n, dev) NA
  initialize <- expression({
    if (NCOL(y) == 1) {
      if (is.factor(y)) y <- y != levels(y)[1L]
      y[wghts == 0] <- 0
      if (any(y < 0 | y > 1)) stop("y values must be 0 <= y <= 1")
      mustart <- (wghts * y + 0.5)/(wghts + 1)
    }
  })
  structure(list(family = "quasibinomial",
                 link = linktemp,
                 linkfun = stats$linkfun,
                 linkinv = stats$linkinv,
                 variance = variance,
                 dev.resids = dev.resids,
                 partialAIC = partialAIC, finalAIC = finalAIC,
                 mu.eta = stats$mu.eta,
                 initialize = initialize,
                 validmu = validmu,
                 valideta = stats$valideta),
            class = "family")
}

#' @title Gamma family function
#' @description Gamma family function
#' @param link function character
#' @export
Gamma_ <- function (link = "inverse")
{
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("inverse", "log", "identity")
  if (linktemp %in% okLinks)
    stats <- make.link(linktemp)
  else if(is.character(link)) stats <- make.link(link)
  else {
    ## what else shall we allow?  At least objects of class link-glm.
    if(inherits(link, "link-glm")) {
      stats <- link
      if(!is.null(stats$name)) linktemp <- stats$name
    } else {
      stop(gettextf('link "%s" not available for gamma family; available links are %s',
                    linktemp, paste(sQuote(okLinks), collapse =", ")),
           domain = NA)
    }
  }
  variance <- function(mu) mu^2
  validmu <- function(mu) all(mu>0)
  dev.resids <- function(y, mu, wt)
    -2 * wt * (log(ifelse(y == 0, 1, y/mu)) - (y - mu)/mu)
  aic <- function(y, n, mu, wt, dev){
    n <- sum(wt)
    disp <- dev/n
    -2*sum(dgamma(y, 1/disp, scale=mu*disp, log=TRUE)*wt) + 2
  }
  partialAIC <- function(y, n, mu, wt, dev) sum(wt)
  finalAIC <- function(paic, n, dev) NA
  initialize <- expression({
    if (any(y <= 0))
      stop("non-positive values not allowed for the 'gamma' family")
    mustart <- y
  })
  simfun <- function(object, nsim) {
    wts <- object$prior.weights
    if (any(wts != 1)) message("using weights as shape parameters")
    ftd <- fitted(object)
    shape <- MASS::gamma.shape(object)$alpha * wts
    rgamma(nsim*length(ftd), shape = shape, rate = shape/ftd)
  }
  structure(list(family = "Gamma",
                 link = linktemp,
                 linkfun = stats$linkfun,
                 linkinv = stats$linkinv,
                 variance = variance,
                 dev.resids = dev.resids,
                 partialAIC = partialAIC,
                 finalAIC = finalAIC,
                 mu.eta = stats$mu.eta,
                 initialize = initialize,
                 validmu = validmu,
                 valideta = stats$valideta,
                 simulate = simfun),
            class = "family")
}

#' @title inverse.gaussian family function
#' @description inverse.gaussian family function
#' @param link function character
#' @export
inverse.gaussian_ <- function(link = "1/mu^2")
{
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  okLinks <- c("inverse", "log", "identity", "1/mu^2")
  if (linktemp %in% okLinks)
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    ## what else shall we allow?  At least objects of class link-glm.
    if(inherits(link, "link-glm")) {
      stats <- link
      if(!is.null(stats$name)) linktemp <- stats$name
    } else {
      stop(gettextf('link "%s" not available for inverse.gaussian family; available links are %s',
                    linktemp, paste(sQuote(okLinks), collapse =", ")),
           domain = NA)
    }
  }
  variance <- function(mu) mu^3
  dev.resids <- function(y, mu, wt)  wt*((y - mu)^2)/(y*mu^2)
  aic <- function(y, n, mu, wt, dev)
    sum(wt)*(log(dev/sum(wt)*2*pi)+1)+3*sum(log(y)*wt)+2
  partialAIC <- function(y, n, mu, wt, dev) NA
  finalAIC <- function(paic, n, dev) NA
  initialize <- expression({
    if(any(y <= 0))
      stop("positive values only are allowed for the 'inverse.gaussian' family")
    #n <- rep.int(1, nobs)
    mustart <- y
  })
  validmu <- function(mu) TRUE
  simfun <- function(object, nsim) NA
  
  structure(list(family = "inverse.gaussian",
                 link = linktemp,
                 linkfun = stats$linkfun,
                 linkinv = stats$linkinv,
                 variance = variance,
                 dev.resids = dev.resids,
                 partialAIC = partialAIC, finalAIC = finalAIC,
                 mu.eta = stats$mu.eta,
                 initialize = initialize,
                 validmu = validmu,
                 valideta = stats$valideta,
                 simulate = simfun),
            class = "family")
}


#' @title quasi family function
#' @description quasi family function
#' @param link function character
#' @param variance choice character
#' @export
quasi_ <- function(link = "identity", variance = "constant")
{
  linktemp <- substitute(link)
  if (!is.character(linktemp)) linktemp <- deparse(linktemp)
  if (linktemp %in% c("logit", "probit", "cloglog", "identity",
                      "inverse", "log", "1/mu^2", "sqrt"))
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  } else {
    stats <- link
    linktemp <- if(!is.null(stats$name)) stats$name else deparse(linktemp)
  }
  vtemp <- substitute(variance)
  if (!is.character(vtemp)) vtemp <- deparse(vtemp)
  variance_nm <- vtemp
  switch(vtemp,
         "constant" = {
           varfun <- function(mu) rep.int(1, length(mu))
           dev.resids <- function(y, mu, wt) wt * ((y - mu)^2)
           validmu <- function(mu) TRUE
           initialize <- expression({
             mustart <- y
             })
         },
         "mu(1-mu)" = {
           varfun <- function(mu) mu * (1 - mu)
           validmu <- function(mu) all(mu>0) && all(mu<1)
           dev.resids <- function(y, mu, wt)
             2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) +
                         (1 - y) * log(ifelse(y == 1, 1, (1 - y)/(1 - mu))))
           initialize <- expression({
             #n <- rep.int(1, nobs)
             mustart <- pmax(0.001, pmin(0.999, y))
             })
         },
         "mu" = {
           varfun <- function(mu) mu
           validmu <- function(mu) all(mu>0)
           dev.resids <- function(y, mu, wt)
             2 * wt * (y * log(ifelse(y == 0, 1, y/mu)) - (y - mu))
           ## 0.1 fudge here matches poisson: S has 1/6.
           initialize <- expression({
              mustart <- y + 0.1 * (y == 0)
             })
         },
         "mu^2" = {
           varfun <- function(mu) mu^2
           validmu <- function(mu) all(mu>0)
           dev.resids <- function(y, mu, wt)
             pmax(-2 * wt * (log(ifelse(y == 0, 1, y)/mu) - (y - mu)/mu), 0)
           initialize <- expression({
              mustart <- y + 0.1 * (y == 0)
             })
         },
         "mu^3" = {
           varfun <- function(mu) mu^3
           validmu <- function(mu) all(mu>0)
           dev.resids <- function(y, mu, wt)
             wt * ((y - mu)^2)/(y * mu^2)
           initialize <- expression({
             # n <- rep.int(1, nobs)
              mustart <- y + 0.1 * (y == 0)
             })
         },
         variance_nm <- NA
  )# end switch(.)
  
  if(is.na(variance_nm)) {
    if(is.character(variance))
      stop(gettextf('\'variance\' "%s" is invalid: possible values are "mu(1-mu)", "mu", "mu^2", "mu^3" and "constant"', variance_nm), domain = NA)
    ## so we really meant the object.
    varfun <- variance$varfun
    validmu <- variance$validmu
    dev.resids <- variance$dev.resids
    initialize <- variance$initialize
    variance_nm <- variance$name
  }
  aic <- function(y, n, mu, wt, dev) NA
  partialAIC <- function(y, n, mu, wt, dev) NA
  finalAIC <- function(paic, n, dev) NA
  structure(list(family = "quasi",
                 link = linktemp,
                 linkfun = stats$linkfun,
                 linkinv = stats$linkinv,
                 variance = varfun,
                 dev.resids = dev.resids,
                 partialAIC = partialAIC, finalAIC = finalAIC,
                 mu.eta = stats$mu.eta,
                 initialize = initialize,
                 validmu = validmu,
                 valideta = stats$valideta,
                 ## character form of the var fun is needed for gee
                 varfun = variance_nm),
            class = "family")
}


# The family function
#' @title family function
#' @description family function
#' @param distr distr character one of "binomial", "poisson", "gaussian", "quasipoisson", 
#'              "quasibinomial", "Gamma", "inverse.gaussian", "quasi"
#' @param link function character
#' @export
family_ <- function(distr, link)
{
  switch(distr,
    "binomial" = {family <- binomial_(link = link)},
    "poisson" = {family <- poisson_(link = link)},
    "gaussian" = {family <- gaussian_(link = link)},
    "quasipoisson" = {family <- quasipoisson_(link = link)},
    "quasibinomial" = {family <- quasibinomial_(link = link)},
    "Gamma" = {family <- Gamma_(link = link)},
    "inverse.gaussian" = {family <- inverse.gaussian_(link = link)},
    "quasi" = {family <- quasi_(link = "identity", variance = "constant")}
  )
  
  return(family)
}
