% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotfns.R
\name{plotDBN}
\alias{plotDBN}
\title{Plotting a DBN}
\usage{
plotDBN(DBN, struct = c("init", "trans"), n.dynamic, n.static)
}
\arguments{
\item{DBN}{binary matrix (or a graph object) representing a 2-step DBN (compact or unrolled)}

\item{struct}{option used to determine if the initial or the transition structure should be plotted; accaptable values are init or trans}

\item{n.dynamic}{number of dynamic variables in one time slice of a DBN}

\item{n.static}{number of static variables in one time slice of a DBN; note that for function to work correctly all static variables have to be in the first n.static columns of the matrix}
}
\description{
This function can be used for plotting initial and transition structures of a dynamic Bayesian network.
}
\examples{
\dontrun{
plotDBN(DBNmat, "trans", n.dynamic=12,n.static=3)
plotDBN(DBNmat, "init", n.dynamic=12,n.static=3)
}
}
