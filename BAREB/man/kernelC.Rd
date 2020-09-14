\name{kernelC}
\alias{kernelC}
\title{The function to get the kernel function value}
\description{
  This function take two configurations x and y, two parameters of the kernel function and returns its kernel function value.
}
\usage{
kernelC(x,y,theta,tau)
}
\arguments{
  \item{x}{a numeric vector, representing one configuration}
  \item{y}{a numeric vector, representing one configuration}
  \item{theta}{a parameter of the DPP's kernel function}
  \item{tau}{a parameter of the DPP's kernel function}
}
\value{
  \code{kernelC(x,y,theta,tau)} returns the value of the kernel function
}
\examples{
  x <- rnorm(5)
  y <- rnorm(5) + 1
  kernelC(x,y,1,1)
}
\seealso{
\link{update_RJ} for a complete example for all functions in this package.
}
\author{Yuliang Li}

