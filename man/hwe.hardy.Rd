\name{hwe.hardy}
\alias{hwe.hardy}
\title{Hardy-Weinberg equlibrium test using MCMC}
\usage{hwe.hardy(a,alleles=3,seed=3000,sample=c(1000,1000,5000))}

\description{
Hardy-Weinberg equilibrium test by MCMC
}
\arguments{
  \item{a}{a trangular array containing the genotype}
  \item{alleles}{number of allele at the locus, greater than or equal to 3}
  \item{seed}{pseudo-random number seed}
  \item{sample}{optional, parameters for MCMC containing \# of chunks, 
                size of chunk and burn-in steps}
}

\source{
http://www.stat.washington.edu/thompson/Genepi/pangaea.shtml,
}

\value{
The returned value is a list containing:

\item{p}{Monte Carlo p value}
\item{se}{standard error of Monte Carlo p value}
\item{swp}{percentage of switches (partial, full and altogether)}
}

\section{References}{

Guo, S.-W. and E. A. Thompson (1992) Performing the exact test of
Hardy-Weinberg proportion for multiple alleles. Biometrics. 48:361--372.
}
\seealso{
\code{\link[gap]{hwe}} 
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
ex2<-hwe.hardy(a,alleles=8)

# example using HLA
data(hla)
x<-hla[,3:4]
y<-pgc(x,handle.miss=0,with.id=1)
n.alleles<-max(x,na.rm=TRUE)
z<-vector("numeric",n.alleles*(n.alleles+1)/2)
z[y$idsave]<-y$wt
hwe.hardy(z,alleles=n.alleles)
## comparison with hwe
hwe(z,data.type="count")
## to create input file for HARDY

print.tri<-function (xx,n)
{
  cat(n,"\n")
  for(i in 1:n) 
  {
     for(j in 1:i)
     {
        cat(xx[i,j]," ")
     }
     cat("\n")
  }
  cat("100 170 1000\n")
}
xx<-matrix(0,n.alleles,n.alleles)
xxx<-lower.tri(xx,diag=TRUE)
xx[xxx]<-z
sink("z.dat")
print.tri(xx,n.alleles)
sink()
## now call as: hwe z.dat z.out
}
}
\author{Sun-Wei Guo, Jing hua Zhao}
\note{Adapted from HARDY, testable with -Dexecutable as standalone program}
\keyword{htest}
