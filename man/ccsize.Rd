\name{ccsize}
\alias{ccsize}
\title{Power and sample size for case-cohort design}
\usage{ccsize(n,q,pD,p1,alpha=0.05,theta,power=NULL)}
\description{
The power of the test is according to 
\deqn{\Phi\left(Z_\alpha+m^{1/2}\theta\sqrt{\frac{p_1p_2p_D}{q+(1-q)p_D}}\right)}{Phi(Z_alpha+m^0.5*theta*sqrt(p_1p_2p_D/q+(1-q)p_D))}
where \eqn{\alpha}{alpha} is the significance level, \eqn{\theta}{theta} is the log-hazard ratio for two groups, \eqn{p_j}{p_j}, 
j=1, 2, are the proportion of the two groups in the population. \eqn{m} is the total number of subjects in the subcohort, 
\eqn{p_D} is the proportion of the failures in the full cohort, and \eqn{q} is the sampling fraction of the subcohort.

Alternatively, the sample size required for the subcohort is \deqn{m=nBp_D/(n-B(1-p_D))}{m=nBp_D/(n-B(1-p_D))}
where \eqn{B=(Z_{1-\alpha}+Z_\beta)^2/(\theta^2p_1p_2p_D)}{B=(Z_{1-alpha}+Z_beta)^2/(theta^2p_1p_2p_D))}, and \eqn{n} is the size of cohort.

}
\arguments{
  \item{n}{the total number of subjects in the subcohort}
  \item{q}{the sampling fraction of the subcohort}
  \item{pD}{the proportion of the failures in the full cohort}
  \item{p1}{proportions of the two groups (p2=1-p1)}
  \item{alpha}{significant level}
  \item{theta}{log-hazard ratio for two groups}
  \item{power}{if specified, the power for which sample size is calculated}
}

\value{
The returned value is a value indicating the power or required sample size.

}

\section{References}{
Cai J, Zeng D. Sample size/power calculation for case-cohort studies.
Biometrics 2004, 60:1015-1024

}
\seealso{
\code{\link[gap]{pbsize}}
}

\examples{
\dontrun{
# Table 1 of Cai & Zeng (2004).

cat("n\tpD\tp1\ttheta\tq\tpower\n")
alpha <- 0.05
n <- 1000
for(pD in c(0.10,0.05))
{
   for(p1 in c(0.3,0.5))
   {
      for(theta in c(0.5,1.0))
      {
         for(q in c(0.1,0.2))
         {
            power <- ccsize(n,q,pD,p1,alpha,theta)
            cat(n,"\t",pD,"\t",p1,"\t",theta,"\t",q,"\t",signif(power,digits=3),"\n")
         }
      }
   }
}

n <- 5000
for(pD in c(0.05,0.01))
{
   for(p1 in c(0.3,0.5))
   {
      for(theta in c(0.5,1.0))
      {
         for(q in c(0.01,0.02))
         {
            power <- ccsize(n,q,pD,p1,alpha,theta)
            cat(n,"\t",pD,"\t",p1,"\t",theta,"\t",q,"\t",signif(power,digits=3),"\n")
         }
      }
   }
}

# ARIC study

n <- 15792
pD <- 0.03
p1 <- 0.25
alpha <- 0.05
theta <- c(1.35,1.40,1.45)
power <- 0.8

s_nb <- c(1463,722,468)

for(i in 1:3)
{
  q <- s_nb[i]/n
  power <- ccsize(n,q,pD,p1,alpha,log(theta[i]))
  ssize <- ccsize(n,q,pD,p1,alpha,log(theta[i]),power)
  cat(n,"\t",pD,"\t",p1,"\t",theta[i],"\t",q,"\t",signif(power,digits=3),"\t",ceiling(ssize),"\n")
}

# EPIC study?

n <- 25000
alpha <- 0.00000001
power <- 0.8
s_pD <- c(0.3,0.2,0.1,0.05)
s_p1 <- seq(0.1,0.5,by=0.1)
s_theta <- seq(1.2,1.8,by=0.2)
s_q <- seq(0.01,0.5,by=0.01)

# direct calculation
for(pD in s_pD)
{
   for(p1 in s_p1)
   {
      for(theta in s_theta)
      {
         ssize <- ccsize(n,q,pD,p1,alpha,log(theta),power)
         if(ssize>0) cat(n,"\t",pD,"\t",p1,"\t",theta,"\t",ssize,"\n")
      }
   }
}

# exhaustive search
nrows <- length(s_pD) * length(s_p1) * length(s_theta) * length(s_q)
powtable <- matrix(rep(0,nrows * 5),ncol=5,byrow=T)
ijkl <- 0
for(pD in s_pD)
{
   for(p1 in s_p1)
   {
      for(theta in s_theta)
      {
         for(q in s_q)
         {
            ijkl <- ijkl + 1
            power <- ccsize(n,q,pD,p1,alpha,log(theta))
            powtable[ijkl,] <- c(pD,p1,theta,q*n,power)
            cat(n,"\t",pD,"\t",p1,"\t",theta,"\t",q*n,"\t",signif(power,digits=3),"\n")
         }
      }
   }
}
}
}
\author{Jing hua Zhao}
\note{Programmed for EPIC study}
\keyword{misc}
