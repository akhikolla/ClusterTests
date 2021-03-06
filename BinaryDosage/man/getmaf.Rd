% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/WriteBinaryDosage.R
\name{getmaf}
\alias{getmaf}
\title{Calculate minor allele frequency}
\usage{
getmaf(dosage, p0, p1, p2)
}
\arguments{
\item{dosage}{Dosage values}

\item{p0}{Pr(g=0) - unused}

\item{p1}{Pr(g=1) - unused}

\item{p2}{Pr(g=2) - unused}
}
\value{
Minor allele frequency
}
\description{
Routine to calculate the minor allele frequency given the dosages.
Missing values for dosage ignored. This function is used internally and
is exported for use in examples. Note: The minor allele in one data set
may be different from another data set. This can make comparing minor
allele frequencies between data sets nonsensical.
}
\examples{
# Get information about binary dosage file
bdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
bdinfo <- getbdinfo(bdfiles = bdfile)
snp1 <- getsnp(bdinfo = bdinfo, 1)
maf <- getmaf(snp1$dosage)
}
