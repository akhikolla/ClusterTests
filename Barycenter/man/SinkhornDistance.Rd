\name{Sinkhorn}
\alias{Sinkhorn}
\title{Sinkhorn Distances (upper bound to EMD)}
\description{The \code{Sinkhorn} algorithm to compute \eqn{N} dual-Sinkhorn divergences, i.e. upper bounds on the earth movers distance (EMD), a.k.a. Wasserstein distance.}
\usage{
Sinkhorn(a, b, costm, lambda = 1, maxIter = 10000, tolerance=10^(-8))
}
\arguments{
\item{a}{Either a \eqn{d_1 x 1} column vector in the probability simplex (nonnegative summing to one) or a \eqn{d_1 x N} matrix, where each column is in the probability simplex.}
\item{b}{A \eqn{d_1 x N} matrix of \eqn{N} vectors in the probability simplex.}
\item{costm}{A matrix of pairwise distances/costs between bins described in \code{a} and bins in the \eqn{b_1},...,\eqn{b_N} histograms.}
\item{lambda}{Non-negative regularization parameter (for small \code{lambda} the Sinkhorn Distance is close to the EMD).}
\item{maxIter}{Maximum number of iterations.}
\item{tolerance}{A threshold for the integrated stopping criterion based on marginal differences.}
}
\value{Returns the Sinkhorn Distances between the given bins. If \code{a} is given by a \eqn{d_1 x 1} column vector the function returns the distances \deqn{[D(a,b_1),...,D(a,b_N)].} If \code{a} is given by a \eqn{d_1 x N} matrix the function returns the distances \deqn{[D(a_1,b_1),...,D(a_N,b_N)].} If \code{a} and \code{b} are simply given by two \eqn{d_1 x 1} and \eqn{d_2 x 1} vectors each in the probability simplex, respectively, then the functions returns a list containing the Sinkhorn Distance as well as the corresponding regularized transport plan.}
\examples{
#Sinkhorn Distances between the first image to the remaining four images in the dataset eight.
#We creat costm simply using a distance matrix on the grid [0,1]x[0,1].
n <- seq(0,1,length.out = dim(eight[[1]])[2])
costm <- as.matrix(dist(expand.grid(n,rev(n)), diag=TRUE, upper=TRUE))
a <- matrix(eight[[1]],28*28,1)
b <- matrix(c(eight[[2]],eight[[3]],eight[[4]],eight[[5]]),28*28,4)
Sinkhorn(a, b, costm)
}
\author{
Marcel Klatt
}
\references{
Cuturi, M.: \code{Sinkhorn Distances: Lightspeed Computation of Optimal Transport}, Advances in Neural Information Processing Systems 26, 2013
}
