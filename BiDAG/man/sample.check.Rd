% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performanceassess.R
\name{sample.check}
\alias{sample.check}
\title{Performance assessment of sampling algorithms against a known Bayesian network}
\usage{
sample.check(
  MCMCchain,
  truedag,
  pbarrier = c(0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2),
  pdag = TRUE,
  burnin = 0.2,
  trans = TRUE
)
}
\arguments{
\item{MCMCchain}{an object of class MCMCres, representing the output of structure sampling function \code{\link{partitionMCMC}} or \code{\link{orderMCMC}} (the latter when parameter \code{chainout}=TRUE)}

\item{truedag}{ground truth DAG which generated the data used in the search procedure; represented by an object of class  \code{\link[graph]{graphNEL}}}

\item{pbarrier}{(optional) a vector of numeric values between 0 and 1, defining posterior probabilities according to which the edges of assessed structures are drawn, please note very low barriers can lead to very dense structures; by default 
\eqn{pbarrier=c(0.99, 0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2)}}

\item{pdag}{logical, if TRUE (default) all DAGs in the MCMCchain are first converted to equivalence class (CPDAG) before the averaging}

\item{burnin}{(optional)  number between \code{0} and \code{1}, indicates the percentage of the samples which will be  the discarded as `burn-in' of the MCMC chain; the rest  of the samples will be used to calculate the posterior probabilities; 0.2 by default}

\item{trans}{logical, for DBNs indicates if model comparions are performed for transition structure; when \code{trans} equals FALSE the comparison is performed for initial structures of estimated models and the ground truth DBN; for usual BNs the parameter is disregarded}
}
\value{
A matrix with the number of rows equal to the number of posterior thresholds tested, and 4 columns reporting for each  thresholded directed graphs the number of true positive edges (`TP'), the number of false positive edges (`FP'), the structural Hamming distance (`SHD') and the posterior threshold
}
\description{
This function calculates the number of true and false positives and the structural Hamming distance between a ground truth DAG and a directed graph summarising a sample of DAGs obtained from an MCMC scheme, as the posterior probability threshold is varied
}
\examples{
gsim.score<-scoreparameters(100, "bge", gsim)
\dontrun{
mapest<-iterativeMCMC(gsim.score)
ordersample<-orderMCMC(gsim.score, MAP=FALSE, startspace=mapest$endspace)
sample.check(ordersample, gsimmat)
}
}
