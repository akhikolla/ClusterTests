\name{RJi_empty}
\alias{RJi_empty}
\title{Function to update the number of site level clusters for one patient level cluster via RJMCMC, when there is no patient in this cluster}
\description{
  This function updates the number of site level clusters for one patient level cluster via RJMCMC, when there is no patient in this cluster. Only used for \code{updateRJ} function.
}
\usage{
RJi_empty(w, K, Gamma,  R,  C,theta, tau, m,  q, hyper_delta)
}
\arguments{
\item{w}{The weights for site level clusters. Should be a vector.}
\item{K}{The number of site level cluster for the patient level cluster. }
\item{Gamma}{The linear coefficients for site level covariates. Should be a matrix.}
\item{R}{The current site level clustering membership}
\item{C}{The DPP related kernel matrices. Should be an array of 3 dimensions}
\item{theta}{The DPP hyper-parameter for site level}
\item{tau}{A fixed DPP hyper-parameter, which we suggest high value, say 10^5}
\item{m}{The number of sites}
\item{q}{The number of site level covariates}
\item{hyper_delta}{The hyper-parameter with default value being 1}
}
\value{
  \code{RJi_empty(w, K, Gamma,  R,  C,theta, tau, m,  q, hyper_delta)} returns a list with following variables:
  \item{K}{The numbers of site level clusters}
  \item{w}{The weights for site level clusters}
  \item{Gamma}{Linear coefficients for site level covariates}
  \item{R}{The site level clusteirng membership}
  \item{C}{The DPP related kernel matrice}
}
\author{Yuliang Li}

