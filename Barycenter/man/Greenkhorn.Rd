\name{Greenkhorn}
\alias{Greenkhorn}
\title{Greenkhorn Distances (approximation to EMD)}
\description{The \code{Greenkhorn} algorithm to approximate the earth movers distance (EMD), a.k.a. Wasserstein distance, between two probability vectors \code{r} and \code{c} with specified cost-matrix \code{costm}.}
\usage{
Greenkhorn(r, c, costm, lambda = 1, maxIter = 10000, tolerance=10^(-8))
}
\arguments{
\item{r}{(n x 1) row vector in the probability simplex (nonnegative summing to one).}
\item{c}{(1 x m) row vector in the probability simplex (nonnegative summing to one).}
\item{costm}{(n x m) matrix of pairwise distances/costs between bins with mass described by \code{r} and bins with mass described by \code{c}.}
\item{lambda}{Non-negative regularization parameter (for small \code{lambda} the Sinkhorn Distance is close to the EMD).}
\item{maxIter}{Maximum number of iterations.}
\item{tolerance}{A threshold for the integrated stopping criterion based on marginal differences.}
}
\value{Returns a list containing the regularized transport plan represented as a \eqn{n x m} matrix as well as the Sinkhorn distance between the given marginals \code{r} and \code{c}.}
\examples{
#Sinkhorn Distances between the first image to the second image in the dataset eight.
#We creat costm simply using a distance matrix on the grid [0,1]x[0,1].
n <- seq(0,1,length.out = dim(eight[[1]])[2])
costm <- as.matrix(dist(expand.grid(n,rev(n)), diag=TRUE, upper=TRUE))
r <- matrix(eight[[1]],28*28,1)
c <- matrix(eight[[2]],1,28*28)
Greenkhorn(r, c, costm)$Distance
}
\author{
Marcel Klatt
}
\references{
Altschuler, J., Weed, J. and Rigollet, P.: \code{Near-linear time approximation algorithms for optimal transport via Sinkhorn iteration}, Advances in Neural Information Processing Systems 30 (NIPS 2017)
}
