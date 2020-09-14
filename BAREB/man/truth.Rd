\name{truth}
\docType{data}
\alias{truth}

\title{The simulation truth}


\usage{
data("truth")
}

\description{
These data record the simulation truth in the paper "BAREB: A Bayesian Repulsive Biclustering Model for Periodontal Data".
It includes the true simulated parameters.

The variables are:
\tabular{ll}{
	\code{S}              \tab  the number of patient level clusters \cr
	\code{E}              \tab  the clustering membership of patient level\cr
	\code{K}              \tab  the numbers of site level clusters\cr
	\code{R}              \tab  the site level clustering membership\cr
	\code{Beta}           \tab  the patient level linear coefficients\cr
	\code{Gamma}          \tab  the site level linear coefficients\cr
	\code{mu}             \tab  the underlying mean for CAL values \cr
	\code{sigma_square}   \tab  the variance of noise for CAL values \cr
	\code{noise}          \tab  the noise for CAL values\cr
	\code{c}              \tab  the parameter for missingness model\cr
	\code{mu.star}        \tab  the mean of latent values for missingness model\cr
	\code{z.star}         \tab  the latent values for missingness model

}
}

\examples{
# output true patient level clustering membership
data("truth")
truth$E

#get the details of the list
str(truth)
}
\keyword{datasets}
