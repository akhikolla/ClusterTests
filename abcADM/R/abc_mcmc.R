#' ABC-MCMC for load duration experiment data
#' @description
#' Use modified ABC-MCMC algorithm to obtain posterior samples of
#' \eqn{\theta = (\mu_a, \sigma_a, \mu_b, \sigma_b, \mu_c, \sigma_c, \mu_n, \sigma_n, \mu_{\sigma_0} , \sigma_{\sigma_0})}{\theta = (\mu_a, \sigma_a, \mu_b, \sigma_b, \mu_c, \sigma_c, \mu_n, \sigma_n, \mu_\sigma0 , \sigma_\sigma0)},
#' given ramp and constant load failure time data.
#'
#'
#' @param n number of posterior samples
#' @param numBurning number of burn-in iterations
#' @param numThining number of thining iterations
#' @param inputD bandwidth \eqn{\delta} for ABC approximation
#' @param dataNames a vector of strings of the names of the datasets, which
#' must be in the format "ID_(\eqn{\tau_c})_(t_c)Y" (see Details and Example)
#' @param verbose displays information messages to console if TRUE
#'
#'
#' @return
#' Returns a matrix of posterior samples where each row is one \eqn{\theta}, and if verbose is TRUE, prints the acceptance rate.
#'
#' @details
#' The generated posterior samples are the parameters associated with (a, b, c, n, \eqn{\eta}), which are the
#' random effects in the Canadian Model for load duration,
#' \deqn{\frac{d}{dt} \alpha(t) = [(a\tau_s)(\tau(t)/\tau_s - \sigma_0)_+]^b + [(c\tau_s)(\tau(t)/\tau_s - \sigma_0)_+]^n\alpha(t),}{d/dt \alpha(t) = [(a * \tau_s) * (\tau(t) / \tau_s - \sigma_0)_+]^b + [(c * \tau_s) * (\tau(t) / \tau_s - \sigma_0)_+]^n * \alpha(t),}
#' where
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
#' * When sample pieces are subject to the load profile
#'
#' \eqn{\tau(t) = kt} if \eqn{t \le T_0}{t <= T_0}
#'
#' \eqn{\tau(t) = \tau_c} if \eqn{t > T_0}
#'
#' where \eqn{\tau_c} is the selected constant-load level under the ramp-loading rate k, and \eqn{T_0} is the
#' time required for the load to reach \eqn{\tau_c} under the ramp-loading rate k.
#'
#'  * The constant load level is assumed to be reached at the ramp-loading rate (k). The ramp-loading rate is 388,440 psi/hour.
#'
#'  * The constant load test ends at time \eqn{t_c} (in years).
#'
#'  * To achieve a ramp-load test, set \eqn{\tau_c} to \code{Inf}.
#'
#'
#'
#' @examples
#' # run the abc-mcmc algorithm to obtain 10 posterior samples
#' # example only, more iterations needed for convergence
#' resTheta = abcMCMC(10, 100, 10, 0.3, c("constLoad_4500_1Y"), TRUE)
#'
#'
#' @references
#'
#' Foschi, R. O., Folz, B., and Yao, F. (1989), Reliability-Based Design of Wood Structures (Vol. 34), Vancouver, BC: Department of Civil Engineering, University of British Columbia.
#'
#' Wong, S. W., & Zidek, J. V. (2018). Dimensional and statistical foundations for accumulated damage models. Wood science and technology, 52(1), 45-65.
#'
#' Yang, C. H., Zidek, J. V., & Wong, S. W. (2019). Bayesian analysis of accumulated damage models in lumber reliability. Technometrics, 61(2), 233-245.
#'
#'
#' @export

abcMCMC = function(n, numBurning, numThining, inputD, dataNames, verbose = FALSE) {

  terminatepg = "Terminating program ... \n"
  origmsg = "Here is the original message: "

  m = tryCatch(
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

  burning = tryCatch(
    {
      as.integer(numBurning)
    },
    warning = function(warning_msg) {
      message("numBurning provided is invalid, please try again")
      message(origmsg)
      message(paste(warning_msg, "\n"))
      return(NULL)
    }
  )

  thining = tryCatch(
    {
      as.integer(numThining)
    },
    warning = function(warning_msg) {
      message("numThining provided is invalid, please try again")
      message(origmsg)
      message(paste(warning_msg, "\n"))
      return(NULL)
    }
  )

  delta = tryCatch(
    {
      as.double(inputD)
    },
    warning = function(warning_msg) {
      message("inputD provided is invalid, please try again")
      message(origmsg)
      message(paste(warning_msg, "\n"))
      return(NULL)
    }
  )

  if (is.null(m) || is.na(m) || is.infinite(m)) {
    stop(terminatepg)
  }
  if (is.null(burning) || is.na(burning) || is.infinite(burning)) {
    stop(terminatepg)
  }
  if (is.null(thining) || is.na(thining) || is.infinite(thining)) {
    stop(terminatepg)
  }
  if (is.null(delta) || is.na(delta) || is.infinite(delta)) {
    stop(terminatepg)
  }

  nameMsg = "Dataframe name format not supported! "

  if (length(dataNames) < 1) {
    lenmsg = "Not enough data frames. "
    stop(paste(lenmsg, terminatepg))
  }

  k = 388440
  numData = length(dataNames)

  w0 = c()
  censored = c()
  T0 = c()

  load = list()
  Ldata = list()

  for (i in 1:numData) {

    dname = dataNames[i]
    underScore_pos = unlist(gregexpr(pattern ='_', dname))

    if (underScore_pos[1] == -1) {
      stop(paste(nameMsg, terminatepg))
    }

    pos1 = underScore_pos[1]
    pos2 = underScore_pos[2]

    # Can add a bunch of condition checking here to ensure the file names are in the correct format
    w0[i] = as.double(substr(dname, pos1 + 1, pos2 - 1)) # get the MOR PSI
    if (is.na(w0[i])) {
      stop(paste(nameMsg, terminatepg))
    }
    T0[i] = w0[i] / k   # get the hour
    duration = substr(dname, pos2 + 1, nchar(dname))
    durLen = nchar(duration)
    if (substr(duration, durLen, durLen) == 'Y') {
      censored[i] = as.double(substr(duration, 1, durLen - 1)) * 8760

      if (is.na(censored[i])) {
        stop(paste(nameMsg, terminatepg))
      }

    } else {
      if (substr(duration, durLen, durLen) == 'M') {
        censored[i] = as.double(substr(duration, 1, durLen - 1)) * 720

        if (is.na(censored[i])) {
          stop(paste(nameMsg, terminatepg))
        }

      } else {
        stop(paste(nameMsg, terminatepg))
      }
    }

    if (verbose) {
      readingmsg = paste("Reading", dname, "\n")
      cat(readingmsg)
    }

    curdat = tryCatch(
      {
        eval(parse(text = dname))    #convert the name of the dataframe to the actual dataframe
      },
      error = function(error_msg) {
        datanfmsg = "Dataframe not found, please check dataframe's name"
        message(datanfmsg)
        message("Here is the original error message: ")
        message(paste(error_msg, "\n"))
        return(NULL)
      }
    )

    if (is.null(curdat)) {
      stop(terminatepg)
    }

    tmpData = as.matrix(curdat)

    if (dim(tmpData)[2] < 1) {

      dataframemsg = "Incorrect format of dataframe. "
      stop(paste(dataframemsg, terminatepg))

    } else if (dim(tmpData)[2] == 1) {

      load[[i]] = as.vector(tmpData[ , 1])
      Ldata[[i]] = as.vector(tmpData[ , 1])

    } else { #number of column is greater than 1

      load[[i]] = as.vector(tmpData[ , 1])
      Ldata[[i]] = as.vector(tmpData[ , 2])
    }

  }
  resTheta = MCMC_ABC(m, burning, thining, delta, k, w0, censored, Ldata)
  acceptRate = resTheta[1]

  if (verbose) {
    acceptMsg = paste("The acceptance rate is: ", acceptRate, "\n")
    cat(acceptMsg)
  }
  return(resTheta[[2]])
}
