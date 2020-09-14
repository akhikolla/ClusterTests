\name{updateZstar}
\alias{updateZstar}
\title{Function to generate new latent values for missingness model}
\description{
  This function takes current parameters, gives updated c in missingess model. Note a double type of value is returned
}
\usage{
updateZstar(mu_star, delta, n, T0)
}


\arguments{
  \item{mu_star}{current matrix of latent value for missingness model}
  \item{delta}{the missing indicator matrix, with 1 means missing}
  \item{n}{number of patients}
  \item{T0}{number of teeth}
  
}
\value{
  \code{updateZstar(mu_star, delta, n, T0)} returns a matrix of new generated latent values.
}
\seealso{
\link{update_RJ} for a complete example for all functions in this package.
}
\author{Yuliang Li}