\name{kin.morgan}
\alias{kin.morgan}
\title{kinship matrix for simple pedigree}
\usage{kin.morgan(ped)}
\description{
kinship matrix according to Morgan v2.1
}

\arguments{
\item{ped}{pedigree id and family trio (id, father id, morther id}
}

\source{
CRAN \url{http://cran.r-project.org}
}

\value{
The returned value is a list containing:

\item{kin}{the kinship matrix}
}

\section{References}{
Morgan V2.1 http://www.stat.washington.edu/thompson/Genepi/MORGAN/morgan.shtml
}

\seealso{
\code{\link[gap]{gif}}
}
\examples{
\dontrun{
# Werner syndrome pedigree
werner<-c(
 1, 0,  0,  1,
 2, 0,  0,  2,
 3, 0,  0,  2,
 4, 1,  2,  1,
 5, 0,  0,  1,
 6, 1,  2,  2,
 7, 1,  2,  2,
 8, 0,  0,  1,
 9, 4,  3,  2,
10, 5,  6,  1,
11, 5,  6,  2,
12, 8,  7,  1,
13,10,  9,  2,
14,12, 11,  1,
15,14, 13,  1)
werner<-t(matrix(werner,nrow=4))
kin.morgan(werner[,1:3])
}
}
\author{Morgan development team, Jing hua Zhao}
\note{The input data is required to be sorted so that parents preceed their children}
\keyword{}
