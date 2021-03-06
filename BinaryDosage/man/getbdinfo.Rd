% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getinfo.R
\name{getbdinfo}
\alias{getbdinfo}
\title{Get information about a binary dosage file}
\usage{
getbdinfo(bdfiles)
}
\arguments{
\item{bdfiles}{Vector of file names. The first is the
binary dosage data containing the dosages and genetic
probabilities. The second file name is the family information
file. The third file name is the SNP information file.
The family and SNP information files are not used if the
binary dosage file is in format 4. For this format the
family and SNP information are in the file with the dosages
and genetic probabilities.}
}
\value{
List with information about the binary dosage file.
This includes family and subject IDs along with
a list of the SNPs in the file. Other information needed
to read the file is also included.
}
\description{
Routine to return information about a binary dosage file.
This information is used by other routines to
allow for quicker extraction of values from the
file.
}
\examples{
vcf1abdfile <- system.file("extdata", "vcf1a.bdose", package = "BinaryDosage")
bdinfo <- getbdinfo(bdfiles = vcf1abdfile)
}
