\name{RJi}
\alias{RJi}
\title{Function to update the number of site level clusters for one patient level cluster via RJMCMC}
\description{
  This function updates the number of site level clusters for one patient level cluster via RJMCMC. Only used for \code{updateRJ} function.
}
\usage{
RJi(w, K, Gamma, Beta, X, Y, Z, R, delta, mu, mu_star, c, sigma_square, C,
             theta, tau, m, n,q,T0, hyper_delta)
}
\arguments{
\item{w}{The weights for site level clusters. Should be a vector.}
\item{K}{The number of site level cluster for the patient level cluster. }
\item{Gamma}{The linear coefficients for site level covariates. Should be a matrix.}
\item{Beta}{The linear coefficients for patient level covariates. Should be a vector.}
\item{X}{The design matrix for patient level clusters.}
\item{Y}{The observed CAL value matrix}
\item{Z}{The design matrix for site level clusters.}
\item{R}{The current site level clustering membership}
\item{delta}{The missing indicator matrix}
\item{mu}{The current CAL mean vector}
\item{mu_star}{The latent value vector for missingness model}
\item{c}{The linear coefficients for missingness model}
\item{sigma_square}{The current noise variance}
\item{C}{The DPP related kernel matrices. Should be an array of 3 dimensions}
\item{theta}{The DPP hyper-parameter for site level}
\item{tau}{A fixed DPP hyper-parameter, which we suggest high value, say 10^5}
\item{m}{The number of sites}
\item{n}{The number of patients}
\item{q}{The number of site level covariates}
\item{T0}{The number of teeth}
\item{hyper_delta}{The hyper-parameter with default value being 1}
}
\value{
  \code{RJi(w, K, Gamma, Beta, X, Y, Z, R, delta, mu, mu_star, c, sigma_square, C, theta, tau, m, n,q,T0, hyper_delta)} returns a list with following variables:
  \item{K}{The numbers of site level clusters}
  \item{w}{The weights for site level clusters}
  \item{Gamma}{Linear coefficients for site level covariates}
  \item{R}{The site level clusteirng membership}
  \item{C}{The DPP related kernel matrice}
}
\author{Yuliang Li}

