\name{Init}
\docType{data}
\alias{Init}

\title{The initial value of patient- and site-level covariates for simulation}


\usage{
data("Init")
}

\description{
These data record the initial value of patient- and site-level covariates in simulation in the paper "BAREB: A Bayesian Repulsive Biclustering Model for Periodontal Data". It is obtained by simple linear regression.

The variables are:
\tabular{ll}{
	\code{Init}          \tab  It has two parts: the initial value of patient-level covariates, \code{Beta};  the initial value of site-level covariates, \code{Gamma}
}
}

\examples{
# output patient level covariates
data("Init")
Init$Beta
Init$Gamma
}
\keyword{datasets}
