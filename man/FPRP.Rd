\name{FPRP}
\alias{FPRP}
\title{False-positive report probability}
\usage{FPRP(a,b,pi0,ORlist,logscale=FALSE)}
\description{
The function calculates the false positive report probability (FPRP), the probability of no true
association beteween a genetic variant and disease given a statistically significant finding,
which depends not only on the observed P value but also on both the prior probability that the
assocition is real and the statistical power of the test. See example for the recommended steps.
}
\arguments{
  \item{a}{parameter value at which the power is to be evaluated}
  \item{b}{the variance for a, or the uppoer point of a 95\%CI if logscale=FALSE}
  \item{pi0}{the prior probabiility that \eqn{H_0}{H0} is true}
  \item{ORlist}{a vector of ORs that is most likely}
  \item{logscale}{FALSE=a,b in orginal scale, TRUE=a, b in log scale}
}

\value{
The returned value is a list with compoents,

\item{p}{p value corresponding to a,b}
\item{power}{the power corresponding to the vector of ORs}
\item{FPRP}{False-positive report probability}
}
\section{References}{
Wacholder S, Chanock S, Garcia-Closas M, El ghomli L, Rothman N. (2004) Assessing the probability that a positive
report is false: an approach for molecular epidemiology studies. J Natl Cancer Inst 96:434-442
}
\seealso{
\code{\link[gap]{BFDP}}
}
\examples{
\dontrun{
# Example by Laure El ghormli & Sholom Wacholder on 25-Feb-2004
# Step 1 - Pre-set an FPRP-level criterion for noteworthiness

T <- 0.2

# Step 2 - Enter values for the prior that there is an association

pi0 <- c(0.25,0.1,0.01,0.001,0.0001,0.00001)

# Step 3 - Enter values of odds ratios (OR) that are most likely, assuming that there is a non-null association

ORlist <- c(1.2,1.5,2.0)

# Step 4 - Enter OR estimate and 95% confidence interval (CI) to obtain FPRP 												

OR <- 1.316
ORlo <- 1.08
ORhi <- 1.60

logOR <- log(OR)
selogOR <- abs(logOR-log(ORhi))/1.96
p <- ifelse(logOR>0,2*(1-pnorm(logOR/selogOR)),2*pnorm(logOR/selogOR))
p
q <- qnorm(1-p/2)
POWER <- ifelse(log(ORlist)>0,1-pnorm((q*selogOR-log(ORlist))/selogOR),pnorm((-q*selogOR-log(ORlist))/selogOR))
POWER
FPRPex <- t(p*(1-pi0)/(p*(1-pi0)+POWER\%o\%pi0))
row.names(FPRPex) <- pi0
colnames(FPRPex) <- ORlist
FPRPex
FPRPex>T

## now turn to FPRP
OR <- 1.316
ORhi <- 1.60
ORlist <- c(1.2,1.5,2.0)
pi0 <- c(0.25,0.1,0.01,0.001,0.0001,0.00001)
z <- FPRP(OR,ORhi,pi0,ORlist,logscale=FALSE)
z
}
}
\author{Jing Hua Zhao}
\keyword{models}
