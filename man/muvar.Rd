\name{muvar}
\alias{muvar}
\alias{muvar1}
\alias{muvar2}
\title{Means and variances under 1- and 2- locus (biallelic) QTL model}
\usage{
muvar1(y1,p1) 

or 

muvar2(y12,p1,p2}
\description{
Function muvar1() gives means and variances under 1-locus QTL model (simple), while muvar2() gives
means and variances under 2-locus QTL model -- it gives results from different avenues.

The results can be kept via R's sink() command or via modifying the C/R codes
}
\arguments{
  \item{y1}{the genotypic means of aa, Aa and AA}
  \item{p1}{the frequency of the lower allele, or the that for the first locus under a 2-locus model}
  \item{y12}{the genotypic means of aa, Aa and AA at the first locus and bb, Bb and BB at the second locus}
  \item{p2}{the frequency of the lower allele at the second locus}
}

\value{
It does not return any value except screent output.
}

\section{References}{
Sham P (1998). Statistics in Human Genetics. Arnold
}

\examples{
\dontrun{
# the default 1-locus model
muvar1()
# the default 2-locus model
muvar2()
}
}
\author{Jing hua Zhao}
\note{Adapted from an earlier C program written for the above book}
\keyword{}
