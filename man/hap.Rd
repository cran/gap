\name{hap}
\alias{hap}
\title{Haplotype reconstruction}
\usage{hap(id,data,nloci,loci,names,mb,pr,po,to,th,maxit,n,ss,rs,rp,ro,rv,sd,mm,
           mi,mc,ds,de,q)}
\description{
Haplotype reconstruction using sorting and trimming algorithms
}
\arguments{
\item{id}{a column of subject id}
\item{data}{genotype table}
\item{nloci}{number of loci}
\item{loci}{number of alleles at all loci}
\item{names}{locus names}
\item{mb}{Maximum dynamic storage to be allocated, in Mb}
\item{pr}{Prior (ie population) probability threshold}
\item{po}{Posterior probability threshold}
\item{to}{Log-likelihood convergence tolerance}
\item{th}{Posterior probability threshold for output}
\item{maxit}{Maximum EM iteration}
\item{n}{Force numeric allele coding (1/2) on output (off)}
\item{ss}{Tab-delimited speadsheet file output (off)}
\item{rs}{Random starting points for each EM iteration (off)}
\item{rp}{Restart from random prior probabilities}
\item{ro}{Loci added in random order (off)}
\item{rv}{Loci added in reverse order (off)}
\item{sd}{Set seed for random number generator (use date+time)}
\item{mm}{Repeat final maximization multiple times}
\item{mi}{Create multiple imputed datasets. If set >0}
\item{mc}{ Number of MCMC steps between samples}
\item{ds}{ Starting value of Dirichlet prior parameter}
\item{de}{ Finishing value of Dirichlet prior parameter}
\item{q}{Quiet operation (off)}
}
\details{
The package can hanlde much larger number of multiallelic loci. 
For large sample size with relatively small number of multiallelic
loci, genecounting should be used.

}

\value{
The returned value is a list containing:

\item{l1}{log-likelihood assuming linkage disequilibrium}
\item{converge}{convergence status, 0=failed, 1=succeeded}
\item{niter}{number of iterations}
}

\section{References}{

Clayton DG (2001) SNPHAP. http://www-gene.cimr.cam.ac.uk/clayton/software

Zhao JH and W Qian (2003) Association analysis of unrelated individuals
using polymorphic genetic markers. RSS 2003, Hassalt, Belgium

}
\seealso{
\code{\link[gap]{genecounting}}
}

\examples{
\dontrun{
# 4 SNP example, to generate hap.out and assign.out alone
data(fsnps)
hap(id=fsnps[,1],data=fsnps[,3:10],nloci=4)

# to generate results of imputations
hap(id=fsnps[,1],data=fsnps[,3:10],nloci=4,ss=1,mi=5)
}
}
\note{adapted from hap}
\keyword{}
