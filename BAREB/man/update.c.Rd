\name{updatec}
\alias{updatec}
\title{Function to update c in missingess model}
\description{
  This function takes current parameters, gives updated c in missingess model. Note a double type of value is returned
}
\usage{
updatec(Zstar, mu,  D, sigmac,sigma_square, n, T0)
}



\arguments{
  \item{Zstar}{generated latent value for missingness model}
  \item{mu}{current estimated mean matrix for CAL}
  \item{D}{the D matrix in the paper}
  \item{sigmac}{hyperparamter for c (variance)}
  \item{sigma_square}{current estimated noise variance}
  \item{n}{number of patients}
  \item{T0}{number of teeth}
  
}
\value{
  \code{updatec(Zstar, mu,  D, sigmac, n, T0)} returns the updated c in missingess model.
}
\seealso{
\link{update_RJ} for a complete example for all functions in this package.
}
\author{Yuliang Li}