% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bdm_merge.R
\name{bdm.optk.s2nr}
\alias{bdm.optk.s2nr}
\title{Find optimal number of clusters based on signal-to-noise-ratio.}
\usage{
bdm.optk.s2nr(bdm, info = T, plot.optk = T, ret.optk = F,
  layer = 1)
}
\arguments{
\item{bdm}{A clustered \var{bdm} instance (\var{i.e.} all up-stream steps performed: \code{bdm.ptse(), bdm.pakde() and bdm.wtt()}.}

\item{info}{Logical value. If TRUE, all merging steps are shown (default value is \code{info = FALSE}).}

\item{plot.optk}{Logical value. If TRUE, this function plots the heuristic measure versus the number of clusters (default value is \code{plot.optk = TRUE})}

\item{ret.optk}{Logical value. For large datasets this computation can take a while and it might be interesting to save it. If TRUE, the function returns a copy of the \var{bdm} instance with the values of S2NR attached as \var{bdm$optk} (default value is \code{ret.optk = FALSE}).}

\item{layer}{The \var{bdm$ptsne} layer to be used (default value is \code{layer = 1}).}
}
\value{
None if \code{ret.optk = FALSE}. Else, a copy of the input \var{bdm} instance with new element \var{bdm$optk} (a matrix).
}
\description{
Performs a recursive merging of clusters based on minimum loss of signal-to-noise-ratio (S2NR). The S2NR is the explained/unexplained variance ratio measured in the high dimensional space based on the given low dimensional clustering. Merging is applied recursively until reaching a configuration of only 2 clusters and the S2NR is measured at each step.
}
\details{
The logic under this heuristic is that neigbouring clusters in the embedding correspond to close clusters in the high dimensional space, \var{i.e.} it is a merging heuristic based on the spatial distribution of clusters. For each cluster (child cluster) we choose the neighboring cluster with steepest gradient along their common border (father cluster). Thus, we get a set of pairs of clusters (child/father) as potential mergings. Given this set of candidates, the merging is performed recursively choosing, at each step, the pair of child/father clusters that results in a minimum loss of S2NR.
A typical situation is that some clusters dominate over all of their neighboring clusters. This clusters have no \var{father}. Thus, once all candidate mergings have been performed we reach a \var{blocked} state where only the dominant clusters remain. This situation identifies a hierarchy level in the clustering. When this situation is reached, the algorithm starts a new merging round, identifying the child/father relations at that level of hierarchy. The process stops when only two clusters remain.
Usually, the clustering hierarchy is clearly depicted by singular points in the S2NR function. This is a hint that the low dimensional clustering configuration is an image of a hierarchycal spatial configuration in the high dimensional space. See \code{bdm.optk.plot()}.
}
\examples{

# --- load mapped dataset
bdm.example()
# --- compute optimal number of clusters and attach the computation
bdm.optk.s2nr(exMap, plot.optk = TRUE, ret.optk = FALSE)
}
