pbsize <- function (gamma=4.5, p=0.15, kp, x2alpha=29.72, zalpha=5.45, z1beta=-0.84)
# population-based sample size
# alpha=5e-8,beta=0.8, alpha would give 5% genome-wide significance level
# x2alpha = 29.72 (Q=29.7168)
#
# lambda is the NCP from the marginal table
# pi is the pr(Affected|aa)
{
  q <- 1-p
  pi <- kp/(gamma*p+q)^2
  lambda <- pi*p*q*(gamma-1)^2/(1-pi*(gamma*p+q)^2)
  n <- (z1beta-zalpha)^2/lambda
  n
}
