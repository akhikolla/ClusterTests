\name{updateC}
\alias{updateC}
\title{Function to obtain the kernel matrix of the determinantal point process}
\description{
  This function takes a matrix and two parameters of kernel function for the determinantal point process (DPP) and gives the kernel matrix for that DPP.
}
\usage{
updateC(Z,theta,tau)
}
\arguments{
  \item{Z}{a matrix, whose rows stand for configurations of the DPP.}
  \item{theta}{a parameter of the DPP's kernel function}
  \item{tau}{a parameter of the DPP's kernel function}
}
\value{
  \code{updateC(Z,theta,tau)} returns the kernel matrix
}
\examples{
  Z <- matrix(rnorm(15), nrow = 5)
  updateC(Z,1,1)
}
\seealso{
\link{update_RJ} for a complete example for all functions in this package.
}
\author{Yuliang Li}

