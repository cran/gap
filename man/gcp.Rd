% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gcp.R
\name{gcp}
\alias{gcp}
\title{Permutation tests using GENECOUNTING}
\usage{
gcp(
  y,
  cc,
  g,
  handle.miss = 1,
  miss.val = 0,
  n.sim = 0,
  locus.label = NULL,
  quietly = FALSE
)
}
\arguments{
\item{y}{A column of 0/1 indicating cases and controls.}

\item{cc}{analysis indicator, 0 = marker-marker, 1 = case-control.}

\item{g}{the multilocus genotype data.}

\item{handle.miss}{a flag with value 1 indicating missing data are allowed.}

\item{miss.val}{missing value.}

\item{n.sim}{the number of permutations.}

\item{locus.label}{label of each locus.}

\item{quietly}{a flag if TRUE will suppress the screen output.}
}
\value{
The returned value is a list containing (p.sim and ph when n.sim > 0):
\itemize{
\item x2obs the observed chi-squared statistic.
\item pobs the associated p value.
\item zobs the observed z value for individual haplotypes.
\item p.sim simulated p value for the global chi-squared statistic.
\item ph simulated p values for individual haplotypes.
}
}
\description{
Permutation tests using GENECOUNTING
}
\details{
This function is a R port of the GENECOUNTING/PERMUTE program which
generates EHPLUS-type statistics including z-tests for individual haplotypes
}
\note{
Built on gcp.c.
}
\examples{
\dontrun{
data(fsnps)
y<-fsnps$y
cc<-1
g<-fsnps[,3:10]

gcp(y,cc,g,miss.val="Z",n.sim=5)
hap.score(y,g,method="hap",miss.val="Z")
}

}
\references{
\insertRef{zhao00}{gap}

\insertRef{zhao04}{gap}

Zhao JH, Qian WD (2003) Association analysis of unrelated individuals using
polymorphic genetic markers -- methods, implementation and application, Royal
Statistical Society, Hassallt-Diepenbeek, Belgium.
}
\seealso{
\code{\link{genecounting}}
}
\author{
Jing Hua Zhao
}
\keyword{htest}
\keyword{models}
