\name{hwe.hardy}
\alias{hwe.hardy}
\title{Hardy-Weinberg equlibrium test}
\usage{hwe.hardy(a,alleles,seed,sample)}
\description{
Hardy-Weinberg equilibrium test by MCMC
}
\arguments{
  \item{a}{a trangular array containing the genotype}
  \item{alleles}{number of allele at the locus}
  \item{seed}{pseudo-random number seed, with default value 3000}
  \item{sample}{optional, parameters for MCMC containing \# of chunks, 
                size of chunk, burn-in steps, with default values to be 1000, 1000, 5000}
}

\source{
http://www.stat.washington.edu/thompson/Genepi/pangaea.shtml,
}

\value{
The returned value is a list containing:

\item{p}{Monte Carlo p value}
\item{se}{standard error of Monte Carlo p value}
\item{swp}{percentage of switches}
}

\section{References}{

Guo, S.-W. and E. A. Thompson (1992) Performing the Exact Test of
Hardy-Weinberg Proportion for Multiple Alleles. Biometrics. 48:361--372.
}
\seealso{
\code{\link[gap]{genecounting}} 
}

\examples{
\dontrun{
# example 2 from hwe.doc:
a<-c(
3,
4, 2,
2, 2, 2,
3, 3, 2, 1,
0, 1, 0, 0, 0,
0, 0, 0, 0, 0, 1,
0, 0, 1, 0, 0, 0, 0,
0, 0, 0, 2, 1, 0, 0, 0)
ex2<-hwe.hardy(a,alleles=8,seed=3000,sample=c(1000,1000,5000))
}
}
\author{Sun-Wei Guo, Jing hua Zhao}
\note{adapted from HARDY, testable with -Dexecutable as standalone program, bug to be fixed under Unix}
\keyword{}
