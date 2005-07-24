\name{kbyl}
\alias{kbyl}
\title{LD statistics for two multiallelic loci}
\usage{kbyl(n1=2,n2=2,h,n,optrho=2,verbose=FALSE)}
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
  \item{verbose}{detailed output of individual statistics}
}

\value{
The returned value is a list containing:

\item{n1}{the number of alleles at marker 1}
\item{n2}{the number of alleles at marker 2}
\item{h}{the haplotype frequency vector}
\item{n}{the number of haplotypes}
\item{Dp}{D'}
\item{VarDp}{variance of D'}
\item{Dijtable}{table of Dij}
\item{VarDijtable}{table of variances for Dij}
\item{Dmaxtable}{table of Dmax}
\item{Dijptable}{table of Dij'}
\item{VarDijptable}{table of variances for Dij'}
\item{X2table}{table of Chi-squares (based on Dij)}
\item{ptable}{table of p values}
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

Zhao, JH (2004). 2LD, GENECOUNTING and HAP: Computer programs for
linkage disequilibrium analysis. Bioinformatics 20:1325-1326

}
\seealso{
\code{\link[gap]{tbyt}}
}

\examples{
\dontrun{

# two examples in the C program 2LD:
# two SNPs as in 2by2.dat
# this can be compared with output from tbyt

h <- c(0.442356,0.291532,0.245794,0.020319)
n <- 481*2
t <- kbyl(2,2,h,n)
t

# two multiallelic markers as in kbyl.dat
# the two-locus haplotype vector is in the file "h"

h <- scan("h")
t <- kbyl(9,5,h,213*2,verbose=T)
}
}
\author{Jing hua Zhao}
\note{adapted from 2ld.c}
\keyword{models}
