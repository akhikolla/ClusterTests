\name{updateR}
\alias{updateR}
\title{Function to update site level clustering in the BAREB model}
\description{
  This function takes current parameters and observed data, gives an updated site level clustering.
}
\usage{
updateR( w ,Gamma,Beta,
    Y, Z, delta,
    mu, mu_star,
    c, S, sigma_square,
    K, E, X,
    m,n, q, p, T0)
}




\arguments{
  \item{w}{current site level clustering prior prob, a matrix}
  \item{Gamma}{current site level linear coefficients array}
  \item{Beta}{current patient level linear coefficients matrix}
  \item{Y}{the CAL observation matrix, with missing values}
  \item{Z}{the site level covariate matrix}
  \item{delta}{the missing indicator matrix, with 1 means missing}
  \item{mu}{current estimated mean matrix for CAL}
  \item{mu_star}{current matrix of latent value for missingness model}
  \item{c}{current c for missingness model}
  \item{S}{number of patient level clusters}
  \item{sigma_square}{current estimated noise variance}
  \item{K}{a vector recording numbers of site level clusters}
  \item{E}{current patient level clustering vector}
  \item{X}{the patient level covariate matrix}
  \item{m}{number of sites}
  \item{n}{number of patients}
  \item{p}{dimension of patient level covariates}
  \item{q}{dimension of site level covariates}
  \item{T0}{number of teeth}
}
\value{
  \code{updateR( w ,Gamma,Beta,  Y, Z, delta, mu, mu_star, c, S, sigma_square, K, E, X, m,n,  p, q, T0)} returns the updated site level clustering.
}
\seealso{
\link{update_RJ} for a complete example for all functions in this package.
}
\author{Yuliang Li}
