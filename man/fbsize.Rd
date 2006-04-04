\name{fbsize}
\alias{fbsize}
\title{Sample size for family-based linkage and association design}
\usage{fbsize(gamma,p,debug=0,error=0)}
\description{
This function implements Risch and Merikangas (1996) statistics 
evaluating power for family-based linkage and association design.
They are potentially useful in the prospect of genome-wide 
association studies.

The function calls auxiliary functions sn() and strlen; sn() 
contains the necessary thresholds for power calculation while
strlen() evaluates length of a string (generic).
}

\arguments{
  \item{gamma}{genotype relative risk assuming multiplicative model}
  \item{p}{frequency of disease allele}
  \item{debug}{verbose output}
  \item{error}{0=use the correct formula,1=the original paper}
}

\value{
The returned value is a list containing:

  \item{gamma}{input gamma}
  \item{p}{input p}
  \item{n1}{sample size for ASP}
  \item{n2}{sample size for TDT}
  \item{n3}{sample size for ASP-TDT}
  \item{lambdao}{lambda o}
  \item{lambdas}{lambda s}
}

\section{References}{

Risch, N. and K. Merikangas (1996). The future of genetic studies of
complex human diseases. Science 273(September): 1516-1517.

Risch, N. and K. Merikangas (1997). Reply to Scott el al. Science
275(February): 1329-1330.

Scott, W. K., M. A. Pericak-Vance, et al. (1997). Genetic analysis of 
complex diseases. Science 275: 1327.

}

\seealso{
\code{\link[gap]{pbsize}}
}
\examples{
\dontrun{
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
    
cat("\nThe family-based result: \n")
cat("\ngamma   p     Y     N_asp   P_A    Het    N_tdt  Het N_asp/tdt  L_o  L_s\n\n")
for(i in 1:12) \{
  g <- models[i,1]
  p <- models[i,2]
  fbsize(g,p)
  if(i\%\%4==0) cat("\n")
\}

# APOE-4, Scott WK, Pericak-Vance, MA & Haines JL
# Genetic analysis of complex diseases 1327
g <- 4.5
p <- 0.15
cat("\nAlzheimer's:\n\n")
fbsize(g,p)
}
}
\author{Jing Hua Zhao}
\note{extracted from rm.c}
\keyword{models}
