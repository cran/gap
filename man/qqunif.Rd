\name{qqunif}
\alias{qqunif}
\title{Q-Q plot for uniformly distributed random variable}
\usage{qqunif(u,logscale=TRUE,base=10,...)}
\description{
This function produces Q-Q plot for a random variable following uniform distribution using log-scale.
}
\arguments{
  \item{u}{a vector of uniformly distributed random variables}
  \item{logscale}{to use logscale}
  \item{base}{the base of the log function}
  \item{...}{other options as appropriae for the qqplot function}
}

\value{
The returned value is a list with components of a qqplot:

\item{x}{-log(expected value for uniform order statistics)}
\item{y}{-log(observed value)}
}

\section{References}{
Davison AC. Statistical Models. Cambridge University Press 2003.
}
\seealso{
\code{\link[SNPassoc]{qqpval}}
}
\examples{
\dontrun{
p <- runif(100)
r <- qqunif(p)
}
}
\author{Jing Hua Zhao}
\keyword{models}
