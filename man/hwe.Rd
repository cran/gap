\name{hwe}
\alias{hwe}
\title{Hardy-Weinberg equlibrium test}
\synopsis{hwe(data, is.count=FALSE, is.genotype=FALSE, yates.correct=FALSE, miss.val=0)}
\usage{
hwe(data,yates.correct=FALSE, miss.val=0)
hwe(data,is.genotype=FALSE, yates.correct=FALSE, miss.val=0)
hwe(data,is.count=FALSE, yates.correct=FALSE, miss.val=0)
}

\description{
Hardy-Weinberg equilibrium test
}
\arguments{
  \item{data}{a rectangular data containing the genotype, or an array of genotype counts}
  \item{is.genotype}{A flag indicating if the data is an array of genotypes}
  \item{is.count}{A flag indicating if the data is an array of genotypes count}
  \item{yates.correct}{A flag indicating if Yates' correction is used for Pearson \eqn{\chi^2}{chi-squared} statistic}
  \item{miss.val}{A list of missing values}
}

\details{
This function obtains Hardy-Weinberg equilibrium test statistics. It can
handle data coded as allele numbers (default), genotype identifiers (by
setting is.genotype=TRUE) and counts corresponding to individual genotypes
 (by setting is.count=TRUE) ; the latter does not need is.genotype to be
specified but requires that genotype counts for all possible genotypes,
i.e. n(n+1)/2, where n is the number of alleles.
}

\value{
The returned value is a list containing:

\item{x2}{Pearson \eqn{\chi^2}{chi-square}}
\item{p.x2}{p value for \eqn{\chi^2}{chi-square}}
\item{lrt}{Log-likelihood ratio test statistic}}
\item{p.lrt}{p value for lrt}}
\item{df}{Degree(s) of freedom}
\item{rho}{\eqn{\chi^2/N}{chi-square/N} the effect size}
}

\seealso{
\code{\link[gap]{hwe.hardy}} 
}

\examples{
\dontrun{
a <- c(3,2,2)
a.out <- hwe(a,is.genotype=TRUE)
a.out
a.out <- hwe(a,is.count=TRUE)
a.out
}
}
\author{Jing hua Zhao}
\keyword{}

