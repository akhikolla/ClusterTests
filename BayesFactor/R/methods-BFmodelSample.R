setMethod('show', signature = c("BFmcmc"),
  function(object){
    show(S3Part(object))
    show(object@model)
  }
)

setAs("BFmcmc" , "mcmc",
      function ( from , to ){
        as.mcmc(from)
      })

setAs("BFmcmc" , "matrix",
      function ( from , to ){
        as.matrix(from)
      })

setAs("BFmcmc" , "data.frame",
      function ( from , to ){
        as.data.frame(from)
      })


#' @rdname recompute-methods
#' @aliases recompute,BFmcmc-method
setMethod('recompute', signature(x = "BFmcmc", progress="ANY"),
  function(x, progress, ...){
    posterior(model=x@model, data = x@data, progress = progress, ...)
  }
)

setMethod('compare', signature(numerator = "BFmcmc", denominator = "BFmcmc"),
  function(numerator, denominator, ...){
    compare(numerator = numerator@model, data = numerator@data, ...) /
      compare(numerator = denominator@model, data = denominator@data, ...)
  }
)

setMethod('compare', signature(numerator = "BFmcmc", denominator = "missing"),
          function(numerator, denominator, ...){
            compare(numerator = numerator@model, data = numerator@data, ...)
          }
)

#' @rdname posterior-methods
#' @aliases posterior,BFmodel,missing,data.frame,missing-method
setMethod("posterior", signature(model="BFmodel", index="missing", data="data.frame", iterations="missing"),
  function(model, index, data, iterations, ...)
    stop("Iterations must be specified for posterior sampling.")
  )

#' @rdname posterior-methods
#' @aliases posterior,BFBayesFactor,missing,missing,missing-method
setMethod("posterior", signature(model="BFBayesFactor", index="missing", data="missing", iterations="missing"),
          function(model, index, data, iterations, ...)
            stop("Iterations must be specified for posterior sampling.")
)

#' @rdname posterior-methods
#' @aliases posterior,BFBayesFactor,numeric,missing,numeric-method
setMethod('posterior', signature(model = "BFBayesFactor", index = "numeric", data = "missing", iterations = "numeric"),
  function(model, index, data, iterations, ...){
    if(length(model[index])>1) stop("Index must specify single element.")
    posterior(model = model[index], iterations = iterations, ...)
  }
)

#' @rdname posterior-methods
#' @aliases posterior,BFBayesFactor,missing,missing,numeric-method
setMethod('posterior', signature(model = "BFBayesFactor", index = "missing", data = "missing", iterations = "numeric"),
  function(model, index=NULL, data, iterations, ...){
    if(length(model)>1) stop("Index argument required for posterior with multiple numerators.")
     posterior(model = model@numerator[[1]], data = model@data, iterations = iterations, ...)
  }
)

#' @rdname posterior-methods
#' @aliases posterior,BFlinearModel,missing,data.frame,numeric-method
setMethod('posterior', signature(model = "BFlinearModel", index = "missing", data = "data.frame", iterations = "numeric"),
  function(model, index = NULL, data, iterations, ...){

    rscaleFixed = rpriorValues("allNways","fixed",model@prior$rscale[['fixed']])
    rscaleRandom = rpriorValues("allNways","random",model@prior$rscale[['random']])
    rscaleCont = rpriorValues("regression",,model@prior$rscale[['continuous']])
    rscaleEffects = model@prior$rscale[['effects']]

    formula = formula(model@identifier$formula)
    checkFormula(formula, data, analysis = "lm")

    factors = fmlaFactors(formula, data)[-1]
    nFactors = length(factors)
    dataTypes = model@dataTypes
    relevantDataTypes = dataTypes[names(dataTypes) %in% factors]

    dv = stringFromFormula(formula[[2]])
    dv = composeTerm(dv)
    if(model@type != "JZS") stop("Unknown model type.")

    if( nFactors == 0 ){
      stop("Sampling from intercept-only model not implemented.")
    }else if(all(relevantDataTypes == "continuous")){
      ## Regression
      X = fullDesignMatrix(formula, data, dataTypes)
      chains = linearReg.Gibbs(y = data[[dv]],covariates = X,iterations = iterations, rscale = rscaleCont, ...)
    }else if(all(relevantDataTypes != "continuous")){
      # ANOVA or t test
      chains = nWayFormula(formula=formula, data = data,
                       dataTypes = dataTypes,
                       rscaleFixed = rscaleFixed,
                       rscaleRandom = rscaleRandom,
                       rscaleEffects = rscaleEffects,
                       iterations = iterations,
                       posterior = TRUE, ...)
    }else{
      # GLM
      chains = nWayFormula(formula=formula, data = data,
                       dataTypes = dataTypes,
                       rscaleFixed = rscaleFixed,
                       rscaleRandom = rscaleRandom,
                       rscaleCont = rscaleCont,
                       rscaleEffects = rscaleEffects,
                       iterations = iterations,
                       posterior = TRUE, ...)
    }

    return(new("BFmcmc",chains, model = model, data = data))

  }
)

#' @rdname posterior-methods
#' @aliases posterior,BFindepSample,missing,data.frame,numeric-method
setMethod('posterior', signature(model = "BFindepSample", index = "missing", data = "data.frame", iterations = "numeric"),
  function(model, index = NULL, data, iterations, ...){
    formula = formula(model@identifier$formula)
    rscale = model@prior$rscale
    interval = model@prior$nullInterval
    nullModel = ( formula[[3]] == 1 )
    chains = ttestIndepSample.Gibbs(formula, data, nullModel, iterations,rscale, interval,...)
    new("BFmcmc",chains,model = model, data = data)
})

#' @rdname posterior-methods
#' @aliases posterior,BFcontingencyTable,missing,data.frame,numeric-method
setMethod('posterior', signature(model = "BFcontingencyTable", index = "missing", data = "data.frame", iterations = "numeric"),
          function(model, index = NULL, data, iterations, ...){
            mod = formula(model@identifier)
            type = model@type
            prior = model@prior$a
            marg = model@prior$fixedMargin
            chains = sampleContingency(mod, type, marg, prior, data = data, iterations = iterations, ...)
            new("BFmcmc",chains,model = model, data = data)
          })

#' @rdname posterior-methods
#' @aliases posterior,BFoneSample,missing,data.frame,numeric-method
setMethod('posterior', signature(model = "BFoneSample", index = "missing", data = "data.frame", iterations = "numeric"),
  function(model, index = NULL, data, iterations, ...){
     mu = model@prior$mu
     rscale = model@prior$rscale
     interval = model@prior$nullInterval
     nullModel = ( model@identifier$formula == "y ~ 0" )
     chains = ttestOneSample.Gibbs(y = data$y, nullModel, iterations = iterations, rscale = rscale,
                          nullInterval = interval, ...)
     new("BFmcmc",chains,model = model, data = data)
})

#' @rdname posterior-methods
#' @aliases posterior,BFmetat,missing,data.frame,numeric-method
setMethod('posterior', signature(model = "BFmetat", index = "missing", data = "data.frame", iterations = "numeric"),
          function(model, index = NULL, data, iterations, ...){
            rscale = model@prior$rscale
            interval = model@prior$nullInterval
            nullModel = ( model@identifier$formula == "d = 0" )
            chains = meta.t.Metrop(t = data$t, n1 = data$n1, n2 = data$n2, nullModel, iterations = iterations,
                                   rscale = rscale, nullInterval = interval,  ...)
            new("BFmcmc",chains, model = model, data = data)
          })


#' @rdname posterior-methods
#' @aliases posterior,BFproportion,missing,data.frame,numeric-method
setMethod('posterior', signature(model = "BFproportion", index = "missing", data = "data.frame", iterations = "numeric"),
          function(model, index = NULL, data, iterations, ...){
            rscale = model@prior$rscale
            p = model@prior$p0
            interval = model@prior$nullInterval
            nullModel = ( model@identifier$formula == "p = p0" )

            chains = proportion.Metrop(y = data$y, N = data$N, nullModel, iterations = iterations,
                                      nullInterval = interval, p = p, rscale = rscale, ...)
            new("BFmcmc",chains, model = model, data = data)
          })

#' @rdname posterior-methods
#' @aliases posterior,BFcorrelation,missing,data.frame,numeric-method
setMethod('posterior', signature(model = "BFcorrelation", index = "missing", data = "data.frame", iterations = "numeric"),
          function(model, index = NULL, data, iterations, ...){
            rscale = model@prior$rscale
            interval = model@prior$nullInterval
            nullModel = ( model@identifier$formula == "rho = 0" )

            chains = correlation.Metrop(y = data$y, x = data$x, nullModel, iterations = iterations,
                                       nullInterval = interval, rscale = rscale, ...)
            new("BFmcmc",chains, model = model, data = data)
          })



###########
## S3
###########

as.mcmc.BFmcmc <- function(x, ...){
  return(S3Part(x))
}

as.matrix.BFmcmc <- function(x,...){
  return(as.matrix(S3Part(x)))
}

as.data.frame.BFmcmc <- function(x, row.names=NULL,optional=FALSE,...){
  return(as.data.frame(S3Part(x)))
}


