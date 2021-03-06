#' Stepwise backward regression
#'
#' @description
#' Build regression model from a set of candidate predictor variables by
#' removing predictors based on p values, in a stepwise manner until there is
#' no variable left to remove any more.
#'
#' @param model An object of class \code{lm}; the model should include all
#'   candidate predictor variables.
#' @param prem p value; variables with p more than \code{prem} will be removed
#'   from the model.
#' @param details Logical; if \code{TRUE}, will print the regression result at
#'   each step.
#' @param x An object of class \code{blr_step_p_backward}.
#' @param print_plot logical; if \code{TRUE}, prints the plot else returns a plot object.
#' @param ... Other inputs.
#'
#' @return \code{blr_step_p_backward} returns an object of class \code{"blr_step_p_backward"}.
#' An object of class \code{"blr_step_p_backward"} is a list containing the
#' following components:
#'
#' \item{model}{model with the least AIC; an object of class \code{glm}}
#' \item{steps}{total number of steps}
#' \item{removed}{variables removed from the model}
#' \item{aic}{akaike information criteria}
#' \item{bic}{bayesian information criteria}
#' \item{dev}{deviance}
#' \item{indvar}{predictors}
#'
#' @references
#' Chatterjee, Samprit and Hadi, Ali. Regression Analysis by Example. 5th ed. N.p.: John Wiley & Sons, 2012. Print.
#'
#' @examples
#' \dontrun{
#' # stepwise backward regression
#' model <- glm(honcomp ~ female + read + science + math + prog + socst,
#'   data = hsb2, family = binomial(link = 'logit'))
#' blr_step_p_backward(model)
#'
#' # stepwise backward regression plot
#' model <- glm(honcomp ~ female + read + science + math + prog + socst,
#'   data = hsb2, family = binomial(link = 'logit'))
#' k <- blr_step_p_backward(model)
#' plot(k)
#'
#' # final model
#' k$model
#'
#' }
#'
#'
#' @family variable selection procedures
#'
#' @export
#'
blr_step_p_backward <- function(model, ...) UseMethod("blr_step_p_backward")

#' @export
#' @rdname blr_step_p_backward
#'
blr_step_p_backward.default <- function(model, prem = 0.3, details = FALSE, ...) {

  blr_check_model(model)
  blr_check_logic(details)
  blr_check_npredictors(model, 3)
  blr_check_values(prem, 0, 1)


  l        <- model$data
  nam      <- colnames(attr(model$terms, "factors"))
  response <- names(model$model)[1]
  preds    <- nam
  cterms   <- preds
  ilp      <- length(preds)
  end      <- FALSE
  step     <- 0
  rpred    <- c()
  aic      <- c()
  bic      <- c()
  dev      <- c()

  cat(format("Backward Elimination Method", justify = "left", width = 27), "\n")
  cat(rep("-", 27), sep = "", "\n\n")
  cat(format("Candidate Terms:", justify = "left", width = 16), "\n\n")
  for (i in seq_len(length(nam))) {
    cat(paste(i, ".", nam[i]), "\n")
  }
  cat("\n")

  cat("We are eliminating variables based on p value...")
  cat("\n")

  cat("\n")
  if (!details) {
    cat("Variables Removed:", "\n\n")
  }

  while (!end) {
    m <- glm(paste(response, "~", paste(preds, collapse = " + ")), l, family = binomial(link = 'logit'))
    m_sum <- Anova(m, test.statistic = "Wald")
    pvals <- m_sum$`Pr(>Chisq)`
    maxp  <- which(pvals == max(pvals))

    suppressWarnings(
      if (pvals[maxp] > prem) {

        step   <- step + 1
        rpred  <- c(rpred, preds[maxp])
        preds  <- preds[-maxp]
        lp     <- length(rpred)
        fr     <- glm(paste(response, "~",
                                  paste(preds, collapse = " + ")), l, family = binomial(link = 'logit'))
        mfs    <- blr_model_fit_stats(fr)
        aic    <- c(aic, mfs$m_aic)
        bic    <- c(bic, mfs$m_bic)
        dev    <- c(dev, mfs$m_deviance)

        if (interactive()) {
          cat(paste("-", rev(rpred)[1]), "\n")
        } else {
          cat(paste("-", rev(rpred)[1]), "\n")
        }

        if (details) {
          cat("\n")
          cat(paste("Backward Elimination: Step", step, "\n\n"), paste("Variable", rpred[lp], "Removed"), "\n\n")
          m <- blr_regress(paste(response, "~", paste(preds, collapse = " + ")), l)
          print(m)
          cat("\n\n")
        }
      } else {
        end <- TRUE
        cat("\n")
        cat("No more variables satisfy the condition of p value = ", prem)
        cat("\n")
      }
    )
  }

  if (details) {
    cat("\n\n")
    cat("Variables Removed:", "\n\n")
    for (i in seq_len(length(rpred))) {
      if (interactive()) {
        cat(paste("-", rpred[i]), "\n")
      } else {
        cat(paste("-", rpred[i]), "\n")
      }
    }
  }

  cat("\n\n")
  cat("Final Model Output", "\n")
  cat(rep("-", 18), sep = "", "\n\n")

  fi <- blr_regress(
    paste(response, "~", paste(preds, collapse = " + ")),
    data = l
  )
  print(fi)

  final_model <- glm(paste(response, "~", paste(preds, collapse = " + ")),
    data = l, family = binomial(link = 'logit'))

  out <- list(removed    = rpred,
              indvar     = cterms,
              steps      = step,
              bic        = bic,
              aic        = aic,
              dev        = dev,
              model      = final_model)

  class(out) <- "blr_step_p_backward"

  return(out)
}

#' @export
#'
print.blr_step_p_backward <- function(x, ...) {
  if (x$steps > 0) {
    print_step_backward(x)
  } else {
    print("No variables have been removed from the model.")
  }
}



#' @export
#' @rdname blr_step_p_backward
#'
plot.blr_step_p_backward <- function(x, model = NA, print_plot = TRUE, ...) {

  a <- NULL
  b <- NULL

  y <- seq_len(x$steps)

  d4 <- data.frame(a = y, b = x$aic)
  d5 <- data.frame(a = y, b = x$bic)
  d6 <- data.frame(a = y, b = x$dev)

  p4 <- plot_stepwise(d4, "AIC")
  p5 <- plot_stepwise(d5, "BIC")
  p6 <- plot_stepwise(d6, "Deviance")

  myplots <- list(aic = p4, bic = p5, deviance = p6)
  
  if (print_plot) {
    gridExtra::marrangeGrob(myplots, nrow = 2, ncol = 2)
  }
  
  invisible(myplots)
}


