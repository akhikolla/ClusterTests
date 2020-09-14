#' Set the number of simulated samples
#' @description
#' Modifies the sample size of reliability simulations for \code{\link{solveCanADM}}.
#'
#'
#' @param n the number of samples to simulate for each \eqn{\theta} to estimate failure probability
#' @param verbose prints information message to the console if TRUE
#'
#'
#' @details
#' The default number of samples is 1000
#'
#'
#' @examples
#' set_simuSampleSize(10000)
#'
#'
#' @export

set_simuSampleSize = function(n, verbose = FALSE) {
  res = tryCatch(
    {
      as.integer(n)
    },
    warning = function(warning_msg) {
      errmsg = "Input is invalid, please try again"
      message(errmsg)
      message("Here is the original message:")
      message(paste(warning_msg,'\n'))
      return(NULL)
    }
  )
  if (is.null(res) || is.na(res) || is.infinite(res)) {
    terminatemsg = "Terminating program..."
    stop(terminatemsg)
  }
  modifynPerTheta(res, as.integer(verbose))
}


#' Set end time of reliability assessment
#' @description
#' Modifies the ending time of the simulation period (in hours) for the reliability assessment.
#' In \code{\link{solveCanADM}}, the load profile is simulated for a return period of \code{time_end} to determine whether failure has occurred.
#'
#'
#' @param time_end end of simulation period (in hours)
#' @param verbose prints information message to the console if TRUE
#'
#'
#' @details
#' The default is a reliability assessment period of 262980 hours (30 years).
#'
#'
#' @examples
#' set_returnPeriod(131400)
#'
#'
#' @export

set_returnPeriod = function(time_end, verbose = FALSE) {
  res = tryCatch(
    {
      as.double(time_end)
    },
    warning = function(warning_msg) {
      errmsg = "Input is invalid, please try again"
      message(errmsg)
      message("Here is the original message:")
      message(paste(warning_msg,'\n'))
      return(NULL)
    }
  )
  if (is.null(res) || is.na(res) || is.infinite(res)) {
    terminatemsg = "Terminating program..."
    stop(terminatemsg)
  }
  modifyTEnd(res, as.integer(verbose))
}


#' Set time step
#' @description
#' Modifies time step (in hours).
#' The step size to use when solving the Canadian model numerically with \code{\link{solveCanADM}} for reliability assessments.
#'
#'
#' @param time_step step of time
#' @param verbose prints information message to the console if TRUE
#'
#'
#' @details
#' The default step of time is 100 hours.
#'
#'
#' @examples
#' set_timeStep(50)
#'
#'
#' @export

set_timeStep= function(time_step, verbose = FALSE) {
  res = tryCatch(
    {
      as.double(time_step)
    },
    warning = function(warning_msg) {
      errmsg = "Input is invalid, please try again"
      message(errmsg)
      message("Here is the original message:")
      message(paste(warning_msg,'\n'))
      return(NULL)
    }
  )
  if (is.null(res) || is.na(res) || is.infinite(res)) {
    terminatemsg = "Terminating program..."
    stop(terminatemsg)
  }
  modifyTimeStep(res, as.integer(verbose))
}
