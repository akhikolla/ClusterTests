\name{updatemustar}
\alias{updatemustar}
\title{Function to update mean latent values for missingness model}
\description{
  This function takes current parameters, gives updated mean latent values for missingness model
}
\usage{
updatemustar(mu, c, n, T0, D)
}



\arguments{
  \item{mu}{current estimated mean matrix for CAL}
  \item{c}{current c for missingness model}
  \item{n}{number of patients}
  \item{T0}{number of teeth}
  \item{D}{the D matrix in the paper}
}
\value{
  \code{updatemustar(mu, c, n, T0, D)} returns the updated mean latent values for missingness model.
}
\seealso{
\link{update_RJ} for a complete example for all functions in this package.
}
\author{Yuliang Li}