\name{tbyt}
\alias{tbyt}
\title{LD statistics for two SNPs}
\usage{tbyt(h,n)}
\description{
LD statistics for two SNPs
}

\arguments{
  \item{h}{a vector of haplotype frequencies}
  \item{n}{number of haplotypes}
}

\value{
The returned value is a list containing:

\item{h}{the original haplotype frequency vector}
\item{n}{the number of haplotypes}
\item{D}{the linkage disequilibrium parameter}
\item{VarD}{the variance of D}
\item{Dmax}{the maximum of D}
\item{VarDmax}{the variance of Dmax}
\item{Dprime}{the scaled disequilibrium parameter}
\item{VarDprime}{the variance of Dprime}
\item{x2}{the Chi-squared statistic}
}

\section{References}{
Zapata C, Alvarez G, Carollo C (1997) Approximate variance of the standardized
measure of gametic disequilibrium D'. Am. J. Hum. Genet. 61:771-774
}

\seealso{
\code{\link[gap]{kbyl}}
}
\examples{
\dontrun{
h <- c(0.442356,0.291532,0.245794,0.020319)
n <- 481*2
tbyt(h,n)
}
}
\author{Jing hua Zhao}
\note{extracted from 2ld.c}
\keyword{models}
