#' Simulate the time-to-failures given \eqn{\theta} (posterior samples) under load profile
#' @description
#' Give the estimated probability of failure after 30 years
#' (or can be set via \code{\link{set_returnPeriod}})
#' based on a large number of replications with given load profile.
#' The \eqn{\theta}'s, where \deqn{\theta = (\mu_a, \sigma_a, \mu_b, \sigma_b, \mu_c, \sigma_c, \mu_n, \sigma_n, \mu_{\sigma_0}, \sigma_{\sigma_0}),}{\theta = (\mu_a, \sigma_a, \mu_b, \sigma_b, \mu_c, \sigma_c, \mu_n, \sigma_n, \mu_\sigma0, \sigma_\sigma0),}
#' will then be used to solve the Canadian ADM for time-to-failure T_f and
#' probability of failure P_f with and without duration of load (DOL) effect.
#'
#'
#' @param inputTheta a matrix where each row is \eqn{\theta}
#' @param inputPhi performance factor \eqn{\phi} of the load
#' @param verbose displays the index, the number of failure of the current \eqn{\theta}, and the number of failure of the current \eqn{\theta} with no DOL effect to console if TRUE
#'
#'
#' @return
#' Returns a list of three elements. The first element of the list is the time-to-failure T_f.
#' The second element is the probability of failture P_f.
#' The third element is the probability of failure P_f without DOL effect.
#' The index of \eqn{\theta}, the number of failure of the current \eqn{\theta}, and the number of failure of the current \eqn{\theta} with no DOL effect will be printed in this order along the way if \code{verbose} is TRUE.
#'
#'
#' @details
#' The equation of the Canadian model is:
#' \deqn{\frac{d}{dt} \alpha(t) = [(a\tau_s)(\tau(t)/\tau_s - \sigma_0)_+]^b + [(c\tau_s)(\tau(t)/\tau_s - \sigma_0)_+]^n\alpha(t),}{d/dt \alpha(t) = [(a * \tau_s) * (\tau(t) / \tau_s - \sigma_0)_+]^b + [(c * \tau_s) * (\tau(t) / \tau_s - \sigma_0)_+]^n * \alpha(t),} where
#'
#' * (a, b, c, n, \eqn{\eta}) are piece-specific random effects, and
#'
#'       - \eqn{a|\mu_a, \sigma_a \sim Log-Normal(\mu_a, \sigma_a)}{a|\mu_a, \sigma_a ~ Log-Normal(\mu_a, \sigma_a)};
#'
#'       - \eqn{b|\mu_b, \sigma_b \sim Log-Normal(\mu_b, \sigma_b)}{b|\mu_b, \sigma_b ~ Log-Normal(\mu_b, \sigma_b)};
#'
#'       - \eqn{c|\mu_c, \sigma_c \sim Log-Normal(\mu_c, \sigma_c)}{c|\mu_c, \sigma_c ~ Log-Normal(\mu_c, \sigma_c)};
#'
#'       - \eqn{n|\mu_n, \sigma_n \sim Log-Normal(\mu_n, \sigma_n)}{n|\mu_n, \sigma_n ~ Log-Normal(\mu_n, \sigma_n)};
#'
#'       - \eqn{\eta|\mu_{\sigma_0}, \sigma_{\sigma_0} \sim Log-Normal(\mu_{\sigma_0}, \sigma_{\sigma_0})}{\eta|\mu_\sigma0, \sigma_\sigma0 ~ Log-Normal(\mu_\sigma0, \sigma_\sigma0)} and set \eqn{\sigma_0 = \frac{\eta}{1+\eta}}{\sigma_0 = \eta / (1 + \eta)}.
#'
#' * \eqn{(x)_+ = max(x, 0)}.
#'
#' * \eqn{\sigma_0} serves as the stress ratio threshold in that damage starts to accumulate only when \eqn{\tau(t)/\tau_s > \sigma_0}.
#'
#' * The performance factor \eqn{\phi} comes from the load \eqn{\tau(t)=\phi R_o\frac{\gamma{D}_d + {D}_s(t) + {D}_e(t)}{\gamma\alpha_d  + \alpha_l}}{\tau(t) = \phi*R_o*(\gamma*D_d+D_s(t)+D_e(t)) / (\gamma*\alpha_d+\alpha_l)}.
#'
#' * The default time step when solving the Canadian model numerically is 100 hours. It can be set via \code{\link{set_timeStep}}.
#'
#' * The probability is calculated as (number of failed samples / total number of simulation samples). Total number of simulation samples can be set via \code{\link{set_simuSampleSize}}.
#'
#'
#' @examples
#' # This is a small demo with 50 simulated failure times 
#' set_simuSampleSize(50)
#' # simulate the time-to-failures given theta under the residential loads
#' resList = solveCanADM(default_param, 1, TRUE)
#'
#' # Below is a more realistic example, using 1000 simulated failure times
#' # for each posterior draw of theta
#' \donttest{set_simuSampleSize(1000)
#' resTheta = abcMCMC(100, 100, 10, 0.4, c("constLoad_4500_1Y"))
#' resList = solveCanADM(resTheta, 1, TRUE)}
#' 
#'
#' @references
#' Yang, C. H., Zidek, J. V., & Wong, S. W. (2019). Bayesian analysis of accumulated damage models in lumber reliability. Technometrics, 61(2), 233-245.
#'
#' Wong, S. W., & Zidek, J. V. (2018). Dimensional and statistical foundations for accumulated damage models. Wood science and technology, 52(1), 45-65.
#'
#'
#' @export

solveCanADM = function(inputTheta, inputPhi, verbose = FALSE) {

  terminatepg = "Terminating program ... \n"

  phi = tryCatch(
    {
      as.double(inputPhi)
    },
    warning = function(warning_msg) {
      message("inputPhi provided is invalid, please try again")
      message("Here is the original message:")
      message(paste(warning_msg, "\n"))
      return(NULL)
    }
  )
  if (is.null(phi) || is.na(phi) || is.infinite(phi)) {
    stop(terminatepg)
  }

  #Check is paras is a vector, if yes, take the transpose
  if (is.vector(inputTheta)) {
    inputTheta = t(inputTheta)
  }

  theta = as.matrix(inputTheta)

  if (!is.numeric(theta)) {
    inputThetamsg = "inputTheta provided is not valid, please try again. "
    stop(paste(inputTheta, terminatepg))
  }

  if (dim(theta)[2] != 10) {
    odeparamsg = paste("The Canadian ODE requires 10 parameters, only", dim(theta)[2], "is supplied. ")
    stop(paste(odeparamsg, terminatepg))
  }

  res = solveODE(theta, phi, as.integer(verbose))
  #res is a list
  #First element of the list is time, second element is the prob
  return(res)
}
