\name{gcontrol}
\alias{gcontrol}
\title{genomic control}
\usage{gcontrol(data,zeta,kappa,tau2,epsilon,ngib,burn,idum)}

\description{
the genomic control statistics of Devlin and Roeder, the list of
parameters after data are optional (with default value)
}

\arguments{
\item{data}{the data matrix}
\item{zeta}{with default value 1000}
\item{kappa}{with default value 4}
\item{tau2}{with default value 1}
\item{epsilon}{with default value 0.01}
\item{ngib}{number of Gibbs steps, with default value 500}
\item{burn}{number of burn-ins with default value 50}
\item{idum}{seed for pseudorandom number sequence}
}

\source{
\url{http://www.stat.cmu.edu}
}

\value{
The returned value is a list containing:

\item{deltot}{the probability}
\item{x2}{the statistic}
\item{A}{the A vector}
}

\section{References}{
Devlin B, Roeder K (1999) Genomic control for association studies. 
Biometrics 55:997-1004
}

\examples{
\dontrun{
test<-c(1,2,3,4,5,6,  1,2,1,23,1,2, 100,1,2,12,1,1, 
        1,2,3,4,5,61, 1,2,11,23,1,2, 10,11,2,12,1,11)
test<-t(matrix(test,nrow=6))
gcontrol(test)
}
}

\author{Bobby Jones, Jing hua Zhao}

\note{Adapted from gcontrol.c Bobby Jones and Kathryn Roeder, 
use -Dexecutable for standalone program, function getnum in the original 
code needs \%*s to skip id string}

\keyword{}

