\name{htr}
\alias{htr}
\title{Haplotype trend regression}
\usage{htr(y,x,n.sim=0)}
\description{
Haplotype trend regression (with permutation)
}
\arguments{
  \item{y}{a vector of phenotype}
  \item{x}{a haplotype table}
  \item{n.sim}{the number of permutations}
}

\value{
The returned value is a list containing:

\item{f}{the F statistic for overall association}
\item{p}{the p value for overall association}
\item{fv}{the F statistics for individual haplotypes}
\item{pi}{the p values for individual haplotypes}
}

\section{References}{
Zaykin DV, Westfall PH, Young SS, Karnoub MA, Wagner MJ, Ehm MG (2002) Testing association of statistically inferred haplotypes with 
discrete and continuous traits in samples of unrelated individuals. Hum. Hered. 53:79-91
}
\seealso{
\code{\link[gap]{htr}}
}

\examples{
\dontrun{
# 26-10-03
test2<-read.table("test2.dat")
y<-test2[,1]
x<-test2[,-1]
y<-as.matrix(y)
x<-as.matrix(x)
htr.test2<-htr(y,x)
htr.test2
htr.test2<-htr(y,x,n.sim=10)
htr.test2

# 13-11-2003
data(apoeapoc)
apoeapoc.gc<-gc.em(apoeapoc[,5:8])
y<-apoeapoc$y
for(i in 1:length(y)) if(y[i]==2) y[i]<-1
htr(y,apoeapoc.gc$htrtable)
}
}

\author{Dimitri Zaykin, Jing hua Zhao}
\note{adapted from emgi.cpp, a pseudorandom number seed will be added on}
\keyword{regression}
