\name{updateE}
\alias{updateE}
\title{Function to update patient level clustering in the BAREB model}
\description{
  This function takes current parameters and observed data, gives an updated patient level clustering.
}
\usage{
updateE( Beta, Gamma,w,
    X, Y, Z, delta,
    E, R, S, Ds,
    mu, mustar,
    sigma, c,
    n, m, T0, p, q, D)
}


\arguments{
  \item{Beta}{current patient level linear coefficients matrix}
  \item{Gamma}{current site level linear coefficients array}
  \item{w}{current patient level clustering prior prob, a vector}
  \item{X}{the patient level covariate matrix}
  \item{Y}{the CAL observation matrix, with missing values}
  \item{Z}{the site level covariate matrix}
  \item{delta}{the missing indicator matrix, with 1 means missing}
  \item{E}{current patient level clustering vector}
  \item{R}{current site level clustering matrix}
  \item{S}{number of patient level clusters}
  \item{Ds}{a vector recording numbers of site level clusters}
  \item{mu}{current estimated mean matrix for CAL}
  \item{mustar}{current matrix of latent value for missingness model}
  \item{sigma}{current estimated noise variance}
  \item{c}{current c for missingness model. It is a vector}
  \item{n}{number of patients}
  \item{m}{number of sites}
  \item{T0}{number of teeth}
  \item{p}{dimension of patient level covariates}
  \item{q}{dimension of site level covariates}
  \item{D}{the D matrix in the paper}
}
\value{
  \code{updateE( Beta, Gamma,w, X, Y, Z, delta,E, R, S, Ds, mu, mustar, sigma, c, n, m, T0, p, q, D)} returns a list with following variables:
  \item{E}{the updated patient level clustering}
  \item{Ds}{new vector recording the numbers of site level clusters}
  \item{mu}{the updated mu computed by updated E}
  \item{mustar}{the updated mustar computed by updated E}
}
\seealso{
\link{update_RJ} for a complete example for all functions in this package.
}
\author{Yuliang Li}
