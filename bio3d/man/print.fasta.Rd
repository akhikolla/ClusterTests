\name{print.fasta}
\alias{print.fasta}
\alias{.print.fasta.ali}
\title{ Printing Sequence Alignments }
\description{
    Print method for fasta and pdbs sequence alignment objects.
}
\usage{
\method{print}{fasta}(x, alignment=TRUE, ...)
.print.fasta.ali(x, width = NULL, col.inds = NULL, numbers = TRUE,
         conservation=TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ a sequence alignment object as obtained from the functions  
    \code{\link{read.fasta}}, \code{\link{read.fasta.pdb}},
    \code{\link{pdbaln}}, \code{\link{seqaln}}, etc. }
  \item{alignment}{ logical, if TRUE the sequence alignment will be
    printed to screen. }
  \item{width}{ a single numeric value giving the number of residues per 
    printed sequence block. By default this is determined from considering 
    alignment identifier widths given a standard 85 column terminal
    window. }
  \item{col.inds}{ an optional numeric vector that can be used to select 
    subsets of alignment positions/columns for printing. } 
  \item{numbers}{ logical, if TRUE position numbers and a tick-mark every 
    10 positions are printed above and below sequence blocks. }
  \item{conservation}{ logical, if TRUE conserved and semi-conserved
    columns in the alignment are marked with an \sQuote{*} and 
    \sQuote{^}, respectively. }
  \item{\dots}{ additional arguments to \sQuote{.print.fasta.ali}. }
}
\value{
  Called mostly for its effect but also silently returns block
  divided concatenated sequence strings as a matrix.  
}
\references{
  Grant, B.J. et al. (2006) \emph{Bioinformatics} \bold{22}, 2695--2696.}
\author{ Barry Grant }
\seealso{
  \code{\link{read.fasta}}, \code{\link{read.fasta.pdb}},
  \code{\link{pdbaln}}, \code{\link{seqaln}}
}
\examples{

     file <- system.file("examples/kif1a.fa",package="bio3d")
     aln  <- read.fasta(file)
     print(aln)

     # print(aln, col.inds=30:100, numbers=FALSE)
}
\keyword{ utilities }

