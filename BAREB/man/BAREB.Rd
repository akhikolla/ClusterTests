\name{BAREB}
\alias{BAREB}
\title{
The BAREB package: summary information
}
\description{
This package implements BAREB (A Bayesian Repulsive Biclustering Model for Periodontal Data)
}
\section{Main Features}{
The functions in the pacakge implement BAREB. BAREB can simultaneously cluster periodontal disease patients and their tooth site after taking the patient- and site-level covariates into consideration. BAREB uses the determinantal point process prior to induce diversity among different biclusters to facilitate parsimony and interpretability. In addition, since periodontal diseases are the leading cause for tooth loss, the missing data mechanism is non-ignorable. Such nonrandom missingness is incorporated into BAREB.
}
\section{Functions}{
The main functions are \code{updateBeta}, \code{update.theta.beta}, \code{update.theta.gamma}, \code{updatec}, \code{updateE}, \code{updateGamma}, \code{updatemu}, \code{updatemustar}, \code{updateR}, \code{updateZstar}, \code{update_RJ}, \code{update_sigma_square}, \code{update_w}, and \code{update_w_beta}; other
functions intended for direct access by the user are:
\code{kernelC} and \code{updateC}.  
There are undocumented functions which are called by these.
}
\section{Requirements}{
R version >= 3.4.3. Packages \code{Rcpp} and \code{RcppArmadillo} are used so that complicated functions are implemented in C++ to speed up.
}
\section{Version}{
This is version 1.0.
}
\author{
Yuliang Li (Dept Applied Mathematics and Statistics, Johns Hopkins University, USA) and
Yanxun Xu (Dept Applied Mathematics and Statistics, Johns Hopkins University, USA) and 
Dipankar Bandyopadhyay (Dept Biostatistics, Virginia Commonwealth University, USA)
.
Please send comments, error reports, etc. to the maintainer (Yuliang) via email.
}
\section{Licence}{
This package and its documentation are usable under the terms of the "GNU
General Public License", a copy of which is distributed with the package.
}
