% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/performanceassess.R
\name{edges.posterior}
\alias{edges.posterior}
\title{Estimating posterior probabilities of single edges}
\usage{
edges.posterior(MCMCchain, pdag = FALSE, burnin = 0.2, endstep = 1)
}
\arguments{
\item{MCMCchain}{list of square matrices with elements in \code{\{0,1\}} and representing adjacency matrices of a sample of DAGs obtained via an MCMC scheme
(objects of classes 'MCMCtrace' or 'MCMCres' are also valid data types)}

\item{pdag}{logical, if TRUE (FALSE by default) all DAGs in the MCMCchain are first converted to equivalence class (CPDAG) before the averaging}

\item{burnin}{(optional) number between \code{0} and \code{1}, indicates the percentage of the samples which will be discarded as `burn-in' of the MCMC chain; the rest  of the samples will be used to calculate the posterior probabilities; 0.2 by default}

\item{endstep}{(optional) number between \code{0} and \code{1}; 1 by default}
}
\value{
a square matrix with dimensions equal to the number of variables; each entry \code{[i,j]} is an estimate of the posterior probability of the edge from node \code{i} to node \code{j}
}
\description{
This function estimates the posterior probabilities of edges by averaging over a sample of DAGs
obtained via an MCMC scheme.
}
\examples{
Bostonscore<-scoreparameters(14, "bge", Boston)
\dontrun{
samplefit<-orderMCMC(Bostonscore, iterations=25000,chainout=TRUE)
edgesposterior<-edges.posterior(samplefit, pdag=TRUE, burnin=0.2)
}
}
