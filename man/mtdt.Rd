\name{mtdt}
\alias{mtdt}
\title{Transmission/disequilibrium test of a multiallelic marker}
\usage{mtdt(x,n.sim=0)}
\description{
This function calculates transmission-disequilibrium statistics involving
multiallelic marker.

Inside the function are tril and triu used to obtain lower and upper triangular
matrices.
}
\arguments{
\item{x}{the data table}
\item{n.sim}{the number of simulations}
}

\value{
It returned list contains the following components:
\item{SE}{Spielman-Ewens Chi-square from the observed data}
\item{ST}{Stuart or score Statistic from the observed data}
\item{pSE}{the simulated p value}
\item{sSE}{standard error of the simulated p value}
\item{pST}{the simulated p value}
\item{sST}{standard error of the simulated p value}
}

\section{References}{
Miller MB (1997) Genomic scanning and the transmission/disequilibrium test: 
analysis of error rates. Genet. Epidemiol. 14:851-856

Sham PC (1997) Transmission/disequilibrium tests for multiallelic loci. 
Am. J. Hum. Genet. 61:774-778

Spielman RS, Ewens WJ (1996) The TDT and other family-based tests for
linkage disequilibrium and association. Am. J. Hum. Genet. 59:983-989

Zhao JH, Sham PC, Curtis D (1999) A program for the Monte Carlo evaluation 
of significance of the extended transmission/disequilibrium test. 
Am. J. Hum. Genet. 64:1484-1485

}
\seealso{
\code{\link[gap]{bt}}
}

\examples{
\dontrun{
# Copeman et al (1995) Nat Genet 9: 80-5

x <- matrix(c(0,0, 0, 2, 0,0, 0, 0, 0, 0, 0, 0,
              0,0, 1, 3, 0,0, 0, 2, 3, 0, 0, 0,
              2,3,26,35, 7,0, 2,10,11, 3, 4, 1,
              2,3,22,26, 6,2, 4, 4,10, 2, 2, 0,
              0,1, 7,10, 2,0, 0, 2, 2, 1, 1, 0,
              0,0, 1, 4, 0,1, 0, 1, 0, 0, 0, 0,
              0,2, 5, 4, 1,1, 0, 0, 0, 2, 0, 0,
              0,0, 2, 6, 1,0, 2, 0, 2, 0, 0, 0,
              0,3, 6,19, 6,0, 0, 2, 5, 3, 0, 0,
              0,0, 3, 1, 1,0, 0, 0, 1, 0, 0, 0,
              0,0, 0, 2, 0,0, 0, 0, 0, 0, 0, 0,
              0,0, 1, 0, 0,0, 0, 0, 0, 0, 0, 0),nrow=12)

# See note to bt for the score test obtained by SAS

mtdt(x)
}
}
\author{Mike Miller, Jing Hua Zhao}
\keyword{models}
