#' Simulate samples of failure time using parameters \eqn{\theta}
#' @description
#' Simulate n observations of failure time from
#'  \eqn{\theta = (\mu_a, \sigma_a, \mu_b, \sigma_b, \mu_c, \sigma_c, \mu_n, \sigma_n, \mu_{\sigma_0} , \sigma_{\sigma_0})}{\theta = (\mu_a, \sigma_a, \mu_b, \sigma_b, \mu_c, \sigma_c, \mu_n, \sigma_n, \mu_\sigma0 , \sigma_\sigma0)}
#'  for constant load test or ramp load test.
#'
#'
#' @param paras either a vector or a matrix where each row is \eqn{\theta}
#' @param n number of samples
#' @param tau_c constant load level (can be set to infinity for ramp load test)
#' @param t_c ending time of the constant load test in hours (can be set to infinity for ramp load test)
#' @param verbose print information messages to console if TRUE
#'
#'
#' @return
#' Return a matrix of failure times. Each column contains samples from one \eqn{\theta}.
#' The name of the output matrix must be in the format "ID_\eqn{\tau_c}_t_c(in years)", e.g. simuConstLoad_4500_1Y.
#'
#'
#' @details
#' * The default \eqn{\theta} is provided and can be invoked via \code{\link{default_param}}.
#'
#' * When sample pieces are subject to the load profile
#'
#' \eqn{\tau(t) = kt} if \eqn{t \le T_0}{t <= T_0}
#'
#' \eqn{\tau(t) = \tau_c} if \eqn{t > T_0}
#'
#' where \eqn{\tau_c} is the selected constant-load level under the ramp-loading rate k, and \eqn{T_0} is the
#' time required for the load to reach \eqn{\tau_c} under the ramp-loading rate k.
#'
#' * The constant load level is assumed to be reached at the ramp-loading rate (k). The ramp-loading rate is 388,440 psi/hour.
#'
#' Utilizes routines from the BRENT C++ root-finding library.
#'
#' @examples
#' # simulate data with constant load 4500 psi and test of duration 1Y (8760 hours)
#' simuConstLoad_4500_1Y = simu_failTime(default_param, 30, 4500, 8760, TRUE)
#'
#' # simulate data for ramp load
#' simuRampLoad_Inf_1Y = simu_failTime(default_param, 30, Inf, 8760)
#'
#'
#' @references
#'
#' Brent, R. (2002). Algorithms for Minimization without Derivatives. Dover. ISBN: 0-486-41998-3.
#'
#' Yang, C. H., Zidek, J. V., & Wong, S. W. (2019). Bayesian analysis of accumulated damage models in lumber reliability. Technometrics, 61(2), 233-245.
#'
#'
#' @export

simu_failTime = function(paras, n, tau_c, t_c, verbose = FALSE) {

  terminatemsg = "Terminating program..."
  origmsg = "Here is the original error message: "

  N = tryCatch(
    {
      as.integer(n)
    },
    warning = function(warning_msg) {
      message("n provided is invalid, please try again")
      message(origmsg)
      message(paste(warning_msg, "\n"))
      return(NULL)
    }
  )

  w0 = tryCatch(
    {
      as.double(tau_c)
    },
    warning = function(warning_msg) {
      message("tau_c provided is invalid, please try again")
      message(origmsg)
      message(paste(warning_msg, "\n"))
      return(NULL)
    }
  )

  censored = tryCatch(
    {
      as.double(t_c)
    },
    warning = function(warning_msg) {
      message("t_c provided is invalid, please try again")
      message(origmsg)
      message(paste(warning_msg, "\n"))
      return(NULL)
    }
  )

  if (is.null(N) || is.na(N) || is.infinite(N)) {
    stop(terminatemsg)
  }

  if (is.null(w0) || is.na(w0)) {
    stop(terminatemsg)
  }

  if (is.null(censored) || is.na(censored)) {
    stop(terminatemsg)
  }

  k = 388440
  T0 = w0 / k

  if(!is.numeric(paras)) {
    paramsg = "paras provided is not valid, pleas try again. "
    stop(paste(paramsg, terminatemsg))
  }

  #Check is paras is a vector, if yes, take the transpose
  if (is.vector(paras)) {
    paras = t(paras)
  }

  parasMat = as.matrix(paras)

  #NEW
  if (dim(parasMat)[2] != 10) {
    paramsg = paste("Expecting 10 parameters, but", dim(parasMat)[2], "is supplied. ")
    stop(paste(paramsg, terminatemsg))
  }

  if (verbose) {
    finReadingmsg = paste("Finish reading", dim(parasMat)[1], "parameters\n")
    cat(finReadingmsg)

    startGen = "Start generating data from the parameters\n"
    cat(startGen)
  }

  mat = use_CANL(parasMat, N, w0, k, 0.5)

  if (verbose) {
    fingenratemsg = "Finish generating data\n"
    cat(fingenratemsg)
  }

  return(t(mat))
}


#' Default \eqn{\theta}
#'
#' @description
#' Contains the values of \eqn{\theta = (\mu_a, \sigma_a, \mu_b, \sigma_b, \mu_c, \sigma_c, \mu_n, \sigma_n, \mu_{\sigma_0} , \sigma_{\sigma_0})}{\theta = (\mu_a, \sigma_a, \mu_b, \sigma_b, \mu_c, \sigma_c, \mu_n, \sigma_n, \mu_\sigma0 , \sigma_\sigma0)}.
#'
#'
#' @details
#' Example values of \eqn{\theta} from preliminary data analysis of Hemlock load duration data.
#'
#'
#' @examples
#' default_param
#' #modify one of the values
#' default_param[10] = 0.049
#'
#'
#' @references
#' Yang, C. H., Zidek, J. V., & Wong, S. W. (2019). Bayesian analysis of accumulated damage models in lumber reliability. Technometrics, 61(2), 233-245.
#'
#'
#' @export

default_param = c(-7.5, 0.5, 3.2, 0.2, -22, 0.3, -1, 0.2, 0.15, 0.05)
