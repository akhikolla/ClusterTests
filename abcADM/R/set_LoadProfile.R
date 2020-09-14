#' Set load profile
#'
#' @description
#' Modifies the current load profile for reliability assessments.
#'
#'
#' @param profile desired load profile
#' @param verbose prints the result of the new load profile to console if TRUE
#'
#' @details
#'
#' - The load is defined as \deqn{\tau(t)=\phi R_o\frac{\gamma{D}_d + {D}_s(t) + {D}_e(t)}{\gamma\alpha_d  + \alpha_l}}{\tau(t) = \phi*R_o*(\gamma*D_d+D_s(t)+D_e(t)) / (\gamma*\alpha_d+\alpha_l)}
#'
#' - Following the National Building Code of Canada (NBCC) standards document CAN/CSA-O86, assume that \eqn{\gamma} = 0.25, \eqn{\alpha_d} = 1.25, \eqn{\alpha_l} = 1.5.
#'
#' - \eqn{R_o} is the characteristic value depending on the lumber population. By default, \eqn{R_o} = 2722psi.
#'
#' - \eqn{\phi} is the performance factor.
#'
#' - According to Foschi, Folz, and Yao(1989)
#'
#'     1. \eqn{D_d} is the normalized dead load for the weight of the structure, and \eqn{D_d} ~ N(load_d_mean, load_d_sd). By default, \code{load_d_mean = 1}, \code{load_d_sd = 0.01}.
#'
#'     2. \eqn{D_s(t)} is the sustained load. \eqn{D_e(t)} is the extraordinary load. \eqn{D_s(t)} and \eqn{D_e(t)} are two independent processes.
#'
#'     3. The sizes of the loads are modelled using gamma distributions \eqn{G(k,\theta)} where k and \eqn{\theta} represent the shape and scale parameters.
#'     The random times between and during live load events are modeled using exponential distributions \eqn{Exp(\lambda)} with mean \eqn{\lambda^{-1}}{\lambda^(-1)}.
#'     Parameters for these models were previously fitted using survey data.
#'
#'     4. The process \eqn{D_s(t)} consists of a sequence of successive periods of sustained occupancy each with iid duration \eqn{T_s} ~ Exp(1 / mean_Ts).
#'     During these periods of occupancy \eqn{D_{ls}}{D_ls} ~ G(load_s_shape, load_s_scale) iid. By default, \code{mean_Ts = 10}, \code{load_s_shape = 3.122}, and \code{load_s_scale = 0.0481}.
#'
#'     5. The process \eqn{D_e(t)} consists of brief periods of extraordinary loads, separated by longer periods with no load \eqn{T_e} ~ Exp(mean_Te) of expected duration 1 year.
#'     When extraordinary loads occur, they last for iid periods of random duration \eqn{T_p} ~ Exp(1 / mean_Tp). The normalized loads \eqn{D_{le}}{D_le} during these brief periods
#'     are iid with gamma distribution \eqn{D_{le}}{D_le} ~ G(load_p_shape, load_p_scale). By default, \code{mean_Te = 1}, \code{mean_Tp = 0.03835}, \code{load_p_shape = 0.826}, and \code{load_p_scale = 0.1023}.
#'
#' - The default load profile is the residential load profile, which is included in the package and
#' can be invoked via \code{\link{sample_load_profile}}. It can also be modified with \code{sample_load_profile["parameterName"]}.
#'
#'
#' @examples
#' p1 = sample_load_profile
#' set_LoadProfile(p1)
#'
#'
#' @references
#' Foschi, R. O., Folz, B., and Yao, F. (1989), Reliability-Based Design of Wood Structures (Vol. 34), Vancouver, BC: Department of Civil Engineering, University of British Columbia.
#'
#' Corotis, R. B., and Doshi, V. A. (1977), “Probability Models for Live-Load Survey Results,” Journal of the Structural Division, 103, 1257–1274.
#'
#' Chalk, P. L., and Corotis, R. B. (1980), “Probability Model for Design Live Loads,” Journal of the Structural Division, 106, 2017–2033.
#'
#' Harris, M. E., Bova, C. J., and Corotis, R. B. (1981), “Area-Dependent Pro-cesses for Structural Live Loads,” Journal of the Structural Division, 107,857–872.
#'
#' Yang, C. H., Zidek, J. V., & Wong, S. W. (2019). Bayesian analysis of accumulated damage models in lumber reliability. Technometrics, 61(2), 233-245.
#'
#'
#' @export

set_LoadProfile = function(profile, verbose = FALSE) {

  profilemsg = "Profile format not supported. "
  terminatemsg = "Terminating program..."

  if(is.null(profile)) {
    stop(paste(profilemsg, terminatemsg))
  }

  if (dim(as.matrix(profile))[1] > 1) {
    stop(paste(profilemsg, terminatemsg))
  }

  params_vec = c("load_s_shape", "load_s_scale", "load_p_shape",
                 "load_p_scale", "R0", "alpha_d", "alpha_l", "load_d_mean",
                 "load_d_sd", "gamma", "mean_Ts", "mean_Te", "mean_Tp")

  for (param in params_vec) {

    cur_param = tryCatch(
      {
        as.double(profile[param])
      },
      error = function(error_msg) {
        datanfmsg = paste("Dataframe provided does not contain:", param)
        message(datanfmsg)
        return(NULL)
      }
    )

    if (is.null(cur_param)) {
      next
    }

    if (!is.na(cur_param) && !is.infinite(cur_param)) {
      modifyLoadP(param, cur_param, as.integer(verbose))
    } else {
      paramerrormsg = paste("The parameter", param, "is not available")
      message(paramerrormsg)
    }
  }
  if (verbose) {
    cat("Finish setting load profile parameters")
  }
}


#' Residential load profile
#'
#' @description
#' \code{sample_load_profile} contains the default residential load profile.
#' It can be invoked via \code{sample_load_profile}.
#' Its parameters can be modified with \code{sample_load_profile["parameterName"]}.
#'
#'
#' @details
#'
#' - The load is defined as \deqn{\tau(t)=\phi R_o\frac{\gamma{D}_d + {D}_s(t) + {D}_e(t)}{\gamma\alpha_d  + \alpha_l}}{\tau(t) = \phi*R_o*(\gamma*D_d+D_s(t)+D_e(t)) / (\gamma*\alpha_d+\alpha_l)}
#'
#' - Following the National Building Code of Canada (NBCC) standards document CAN/CSA-O86, assume that \eqn{\gamma} = 0.25, \eqn{\alpha_d} = 1.25, \eqn{\alpha_l} = 1.5.
#'
#' - \eqn{R_o} is the characteristic value depending on the lumber population. By default, \eqn{R_o} = 2722psi.
#'
#' - \eqn{\phi} is the performance factor.
#'
#' - According to Foschi, Folz, and Yao(1989)
#'
#'     1. \eqn{D_d} is the normalized dead load for the weight of the structure, and \eqn{D_d} ~ N(load_d_mean, load_d_sd). By default, \code{load_d_mean = 1}, \code{load_d_sd = 0.01}.
#'
#'     2. \eqn{D_s(t)} is the sustained load. \eqn{D_e(t)} is the extraordinary load. \eqn{D_s(t)} and \eqn{D_e(t)} are two independent processes.
#'
#'     3. The sizes of the loads are modelled using gamma distributions \eqn{G(k,\theta)} where k and \eqn{\theta} represent the shape and scale parameters.
#'     The random times between and during live load events are modeled using exponential distributions \eqn{Exp(\lambda)} with mean \eqn{\lambda^{-1}}{\lambda^(-1)}.
#'     Parameters for these models were previously fitted using survey data.
#'
#'     4. The process \eqn{D_s(t)} consists of a sequence of successive periods of sustained occupancy each with iid duration \eqn{T_s} ~ Exp(1 / mean_Ts).
#'     During these periods of occupancy \eqn{D_{ls}}{D_ls} ~ G(load_s_shape, load_s_scale) iid. By default, \code{mean_Ts = 10}, \code{load_s_shape = 3.122}, and \code{load_s_scale = 0.0481}.
#'
#'     5. The process \eqn{D_e(t)} consists of brief periods of extraordinary loads, separated by longer periods with no load \eqn{T_e} ~ Exp(mean_Te) of expected duration 1 year.
#'     When extraordinary loads occur, they last for iid periods of random duration \eqn{T_p} ~ Exp(1 / mean_Tp). The normalized loads \eqn{D_{le}}{D_le} during these brief periods
#'     are iid with gamma distribution \eqn{D_{le}}{D_le} ~ G(load_p_shape, load_p_scale). By default, \code{mean_Te = 1}, \code{mean_Tp = 0.03835}, \code{load_p_shape = 0.826}, and \code{load_p_scale = 0.1023}.
#'
#'
#' @example
#' lp = sample_load_profile
#' # modify the value of load_s_shape
#' lp["load_s_shape"] = 4.387
#'
#'
#' @references
#' Foschi, R. O., Folz, B., and Yao, F. (1989), Reliability-Based Design of Wood Structures (Vol. 34), Vancouver, BC: Department of Civil Engineering, University of British Columbia.
#'
#' Corotis, R. B., and Doshi, V. A. (1977), “Probability Models for Live-Load Survey Results,” Journal of the Structural Division, 103, 1257–1274.
#'
#' Chalk, P. L., and Corotis, R. B. (1980), “Probability Model for Design Live Loads,” Journal of the Structural Division, 106, 2017–2033.
#'
#' Harris, M. E., Bova, C. J., and Corotis, R. B. (1981), “Area-Dependent Pro-cesses for Structural Live Loads,” Journal of the Structural Division, 107,857–872.
#'
#' Yang, C. H., Zidek, J. V., & Wong, S. W. (2019). Bayesian analysis of accumulated damage models in lumber reliability. Technometrics, 61(2), 233-245.
#'
#'
#' @export

sample_load_profile = data.frame("load_s_shape" = 3.122,
                                 "load_s_scale" = 0.0481,
                                 "load_p_shape" = 0.826,
                                 "load_p_scale" = 0.1023,
                                 "R0" = 2722,
                                 "alpha_d" = 1.25,
                                 "alpha_l" = 1.5,
                                 "load_d_mean" = 1,
                                 "load_d_sd" = 0.1,
                                 "gamma" = 0.25,
                                 "mean_Ts" = 10,
                                 "mean_Te" = 1,
                                 "mean_Tp" = 0.03835)
