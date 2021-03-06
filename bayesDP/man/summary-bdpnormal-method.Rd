% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/summary.R
\name{summary,bdpnormal-method}
\alias{summary,bdpnormal-method}
\title{bdpnormal Object Summary}
\usage{
\S4method{summary}{bdpnormal}(object)
}
\arguments{
\item{object}{object of class \code{bdpnormal}. The result of a call to the
\code{\link{bdpnormal}} function.}
}
\description{
\code{summary} method for class \code{bdpnormal}.
}
\details{
Displays a summary of the \code{bdpnormal} fit including the
  input data, the stochastic comparison between current and historical
  data, and the resulting historical data weight (alpha). If historical
  data is missing then no stochastic comparison nor weight are displayed.

  In the case of a one-arm analysis, the displayed 95 percent
  CI contains the lower and upper limits of the
  augmented mean value of the current data. The displayed
  \code{mean of treatment group} is the mean of the current data
  augmented by the historical data.

  When a control arm is present, a two-arm analysis is carried out.
  Now, the displayed 95 percent CI contains the
  lower and upper limits of the difference between the treatment and
  control arms with the historical data augmented to current data, if
  present. The displayed posterior sample estimates are the
  mean of the treatment and control arms, each of
  which are augmented when historical data are present.
}
