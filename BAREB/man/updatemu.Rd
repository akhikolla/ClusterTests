\name{updatemu}
\alias{updatemu}
\title{Function to update estimated mean CAL values based on current parameters}
\description{
  This function takes current parameters, gives an estimated mean CAL values.
}
\usage{
updatemu(R, Z, X, Gamma, K, Beta, E,  m,n,p, q)
}



\arguments{
  \item{R}{current site level clustering matrix}
  \item{Z}{the site level covariate matrix}
  \item{X}{the patient level covariate matrix}
  \item{Gamma}{current site level linear coefficients array}
  \item{K}{a vector recording numbers of site level clusters}
  \item{Beta}{current patient level linear coefficients matrix}
  \item{E}{current patient level clustering vector}
  \item{m}{number of sites}
  \item{n}{number of patients}
  \item{p}{dimension of patient level covariates}
  \item{q}{dimension of site level covariates}
  
}
\value{
  \code{updatemu(R, Z, X, Gamma, K, Beta, E,  m,n,p, q)} returns the updated estimated mean CAL matrix.
}
\seealso{
\link{update_RJ} for a complete example for all functions in this package.
}
\author{Yuliang Li}