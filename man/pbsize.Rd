\name{pbsize}
\alias{pbsize}
\title{Power for population-based association design}
\usage{pbsize(gamma=4.5, p=0.15, kp, x2alpha=29.72, zalpha=5.45, z1beta=-0.84)}
\description{
This function implements Long et al. (1997) statistics for population-based association
design
}

\arguments{
  \item{gamma}{genotype relative risk assuming multiplicative model}
  \item{p}{frequency of disease allele}
  \item{kp}{population disease prevalence}
  \item{x2alpha}{normal z-deviate}
  \item{zalpha}{normal z-deviate}
  \item{z1beta}{normal z-deviate}
}

\value{
The returned value is scaler containing the required sample size
}

\section{References}{

Scott, W. K., M. A. Pericak-Vance, et al. (1997). Genetic analysis of complex 
diseases. Science 275: 1327.
		
Long, A. D. Grote, M. N. and C. H. Langley (1997). Genetic analysis of complex traits.
Science 275: 1328.
	
}

\seealso{
\code{\link[gap]{fbsize}}
}
\examples{
options(echo=FALSE)
models <- matrix(c(
    4.0, 0.01,
    4.0, 0.10,
    4.0, 0.50, 
    4.0, 0.80,
    2.0, 0.01,
    2.0, 0.10,
    2.0, 0.50,
    2.0, 0.80,
    1.5, 0.01,    
    1.5, 0.10,
    1.5, 0.50,
    1.5, 0.80), ncol=2, byrow=TRUE)
    
g <- 4.5
p <- 0.15
cat("\nAlzheimer's:\n\n")
 
zalpha <- 5.45   # 5.4513104
z1beta <- -0.84

q <- 1-p
pi <- 0.065      # 0.07 generates 163, equivalent to ASP
k <- pi*(g*p+q)^2
s <- (1-pi*g^2)*p^2+(1-pi*g)*2*p*q+(1-pi)*q^2
# LGL formula
lambda <- pi*(g^2*p+q-(g*p+q)^2)/(1-pi*(g*p+q)^2)
# my own
lambda <- pi*p*q*(g-1)^2/(1-pi*(g*p+q)^2)
# not sure about +/-!
n <- (z1beta+zalpha)^2/lambda

# may be used to correct for population prevalence
cat("\nThe population-based result: Kp=",k, "Kq=",s, "n=",ceiling(n),"\n")

# population-based sample size
strlen <- function(x) length(unlist(strsplit(as.character(x),split="")))
kp <- c(0.01,0.05,0.10)
cat("\nRandom ascertainment with disease prevalence\n")
cat("\n          1\%          5\%         10\%\n\n")
for(i in 1:12) {
  g <- models[i,1]
  p <- models[i,2]
  q <- 1-p
  for(j in 1:3) {
    n <- pbsize(g,p,kp[j])
    cat(rep("",12-strlen(ceiling(n))),format(ceiling(n)))
  }
  cat("\n")
  if(i\%\%4==0) cat("\n")
} 
cat("This is only an approximation, a more accurate result\n")
cat("can be obtained by Fisher's exact test\n")
options(echo=TRUE)
}
\author{Jing Hua Zhao}
\note{extracted from rm.c}
\keyword{models}
