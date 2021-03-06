% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/output.R
\name{residuals.bcp}
\alias{residuals.bcp}
\title{Extract model residuals}
\usage{
\method{residuals}{bcp}(object, ...)
}
\arguments{
\item{object}{the result of a call to \code{bcp()}.}

\item{...}{(optional) additional arguments, ignored.}
}
\value{
Residuals extracted from the \code{bcp} object.
}
\description{
residuals method for class \code{bcp}.
}
\examples{
##### A random sample from a few normal distributions #####
testdata <- c(rnorm(50), rnorm(50, 5, 1), rnorm(50))
bcp.0 <- bcp(testdata)
residuals(bcp.0)

}
\seealso{
\code{\link{bcp}} and \code{\link{plot.bcp}}
}
\author{
Xiaofei Wang, Chandra Erdman, and John W. Emerson
}
\keyword{datasets}
