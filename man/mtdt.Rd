% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mtdt.R
\name{mtdt}
\alias{mtdt}
\title{Transmission/disequilibrium test of a multiallelic marker}
\usage{
mtdt(x, n.sim = 0)
}
\arguments{
\item{x}{the data table.}

\item{n.sim}{the number of simulations.}
}
\value{
It returned list contains the following components:
\itemize{
\item SE Spielman-Ewens Chi-square from the observed data.
\item ST Stuart or score Statistic from the observed data.
\item pSE the simulated p value.
\item sSE standard error of the simulated p value.
\item pST the simulated p value.
\item sST standard error of the simulated p value.
}
}
\description{
Transmission/disequilibrium test of a multiallelic marker
}
\details{
This function calculates transmission-disequilibrium statistics involving multiallelic marker.
Inside the function are tril and triu used to obtain lower and upper triangular matrices.
}
\examples{
\dontrun{
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
\references{
\insertRef{miller97}{gap}

\insertRef{sham95}{gap}

\insertRef{spielman96}{gap}

\insertRef{zhao99}{gap}
}
\seealso{
\code{\link{bt}}
}
\author{
Mike Miller, Jing Hua Zhao
}
\keyword{models}
