\name{kbyl}
\alias{kbyl}
\title{LD statistics for two multiallelic loci}
\usage{kbyl(n1,n2,h,n,optrho=2)}
\description{
LD statistics for two multiallelic loci}
}
\arguments{
  \item{n1}{number of alleles at marker 1}
  \item{n2}{number of alleles at marker 2}
  \item{h}{a vector of haplotype frequencies}
  \item{n}{number of haplotypes}
  \item{optrho}{type of contingency table association,
0=Pearson, 1=Tschuprow, 2=Cramer (default)}
}

\value{
The returned value is a list containing:

\item{n1}{the number of alleles at marker 1}
\item{n2}{the number of alleles at marker 2}
\item{h}{the haplotype frequency vector}
\item{n}{the number of haplotypes}
\item{VarDp}{variance of D'}
\item{Dijtable}{table of Dij}
\item{Dmaxtable}{table of Dmax}
\item{Dijptable}{table of Dij'}
\item{VarDijtable}{table of variances for Dij}
\item{VarDijptable}{table of variances for Dij'}
\item{x2}{the Chi-squared statistic}
\item{seX2}{the variance of x2}
\item{rho}{the measure of association}
\item{seR}{the standard error of rho}
\item{optrho}{the method for calculating rho}
\item{klinfo}{the Kullback-Leibler information}
}

\section{References}{
Bishop YMM, Fienberg SE, Holland PW (1975) Discrete Multivariate Analysis
-- Theory and Practice, The MIT press

Cramer H (1946) Mathematical Methods of Statistics. Princeton Univ. Press

Zapata C, Carollo C, Rodriquez S (2001) Sampleing variance and distribution
of the D' measure of overall gametic disequilibrium between multiallelic loci.
Ann. Hum. Genet. 65: 395-406

Zhao, J. H. (2004) 2LD, GENECOUNTING and HAP: Computer programs for
linkage disequilibrium analysis. Bioinformatics, 20, 1325-1326

}
\seealso{
\code{\link[gap]{tbyt}}
}

\examples{
\dontrun{
# example of two SNPs
h <- c(0.442356,0.291532,0.245794,0.020319)
n <- 481*2
kbyl(2,2,h,n)
}
}
\author{Jing hua Zhao}
\note{adapted from 2ld.c}
\keyword{models}
