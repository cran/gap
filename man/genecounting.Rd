\name{genecounting}
\alias{genecounting}
\title{Gene counting for haplotype analysis}
\usage{genecounting(data,weight=NULL,convll,handle.miss,eps,maxit,pl)}
\description{
Gene counting for haplotype analysis with missing data
}
\arguments{
  \item{data}{genotype table}
  \item{weight}{a column of frequencies}
  \item{convll}{set convergence criteria according to log-likelihood, if its value set to 1}
  \item{handle.miss}{to handle missing data, if its value set to 1}
  \item{eps}{the actual convergence criteria, with default value 1e-5}
  \item{maxit}{maximum number of iterations, with default value 50}
  \item{pl}{criteria for trimming haplotypes according to posterior probabilities}
}

\value{
The returned value is a list containing:

\item{h}{haplotype frequency estimates under linkage disequilibrium}
\item{h0}{haplotype frequency estimates under linkage equilibrium}
\item{prob}{genotype probability estimates}
\item{lnl0}{log-likelihood under linkage equilibrium}
\item{lnl1}{log-likelihood under linkage disequilibrium}
\item{hapid}{unique haplotype identifier (defunct, see gc.em)}
\item{npusr}{number of parameters according user-given alleles}
\item{npdat}{number of parameters according to observed}
\item{htrtable}{design matrix for haplotype trend regression (defunct, see gc.em)}
\item{iter}{number of iterations used in gene counting}
\item{converge}{a flag indicating convergence status of gene counting}
\item{Dprime}{a multilocus version of Hedrick's D'}
}

\section{References}{

Zhao, J. H., Lissarrague, S., Essioux, L. and P. C. Sham (2002).
Gene-counting for haplotype analysis with missing genotypes.
Bioinformatics 18(12):1694-1695

Zhao, J. H. and P. C. Sham (2003). Generic number systems and haplotype
analysis. Comp Meth Prog Biomed 70: 1-9

}
\seealso{
\code{\link[gap]{gc.em}}, \code{\link[gap]{kbyl}}
}

\examples{
\dontrun{
# Now we use the HLA data for testing
data(hla)
hla.gc<-genecounting(hla[,3:8])
summary(hla.gc)
hla.gc$lnl0
hla.gc$lnl1

# Now we use ALDH2 data
data(aldh2)
aldh2.gc<-genecounting(aldh2[,3:6],handle.miss=1)
summary(aldh2.gc)
aldh2.gc$lnl0
aldh2.gc$lnl1
}
}
\author{Jing hua Zhao}
\note{adapted from GENECOUNTING}
\keyword{}
