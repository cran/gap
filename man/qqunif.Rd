\name{qqunif}
\alias{qqunif}
\title{Q-Q plot for uniformly distributed random variable}
\usage{qqunif(u,logscale=TRUE,base=10,col=palette()[4],lcol=palette()[2],...)}
\description{
This function produces Q-Q plot for a random variable following uniform distribution with or
without using log-scale.
}
\arguments{
  \item{u}{a vector of uniformly distributed random variables}
  \item{logscale}{to use logscale}
  \item{base}{the base of the log function}
  \item{col}{color for points}
  \item{lcol}{color for the diagonal line}
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
\code{\link[gap]{mhtplot}}
}
\examples{
\dontrun{
# Q-Q Plot for 1000 U(0,1) r.v., marking those <= 1e-5
u_obs <- runif(1000)
r <- qqunif(u_obs,pch=21,bg="blue",bty="n")
u_exp <- r$y
hits <- u_exp >= 2.30103
points(r$x[hits],u_exp[hits],pch=21,bg="green")
}
}
\author{Jing Hua Zhao}
\keyword{hplot}
