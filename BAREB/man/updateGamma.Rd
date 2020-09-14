\name{updateGamma}
\alias{updateGamma}
\title{Function to update site level linear coefficients in the BAREB model}
\description{
  This function takes current parameters and observed data, gives an updated site level linear coefficients.
}
\usage{
updateGamma(X, Y, Z, delta,
    Beta, Gamma, E, R,
    S, Ds, mu, mustar,
    sigma, c, step, runif,
    n, m, T0, p, q, D, theta,tau)
}

\arguments{
  \item{X}{the patient level covariate matrix}
  \item{Y}{the CAL observation matrix, with missing values}
  \item{Z}{the site level covariate matrix}
  \item{delta}{the missing indicator matrix, with 1 means missing}
  \item{Beta}{current patient level linear coefficients matrix}
  \item{Gamma}{current site level linear coefficients array}
  \item{E}{current patient level clustering vector}
  \item{R}{current site level clustering matrix}
  \item{S}{number of patient level clusters}
  \item{Ds}{a vector recording numbers of site level clusters}
  \item{mu}{current estimated mean matrix for CAL}
  \item{mustar}{current matrix of latent value for missingness model}
  \item{sigma}{current estimated noise variance}
  \item{c}{current c for missingness model. It is a vector}
  \item{step}{an array of steps for M-H}
  \item{runif}{an array of uniform random variables for deciding whether to accept new proposed point in M-H}
  \item{n}{number of patients}
  \item{m}{number of sites}
  \item{T0}{number of teeth}
  \item{p}{dimension of patient level covariates}
  \item{q}{dimension of site level covariates}
  \item{D}{the D matrix in the paper}
  \item{theta}{parameter for DPP}
  \item{tau}{parameter for DPP}
}
\value{
  \code{updateGamma(X, Y, Z, delta, Beta, Gamma, E, R, S, Ds, mu, mustar, sigma, c, step, runif, n, m, T0, p, q, D, theta,tau)} returns a list with following variables:
  \item{Gamma}{the updated site level linear coefficients}
  \item{mu}{the updated mu computed by updated Gamma}
  \item{mustar}{the updated mustar computed by updated Gamma}
}
\seealso{
\link{update_RJ} for a complete example for all functions in this package.
}
\author{Yuliang Li}
