qqunif <- function(u,logscale=TRUE,base=10,...)
{
  u <- u[!is.na(u)]
  n <- length(u)
  r <- (1:n)/(n+1)
  y <- -log(u,base)
  x <- -log(r,base)
  if (logscale) z <- qqplot (y,x,xlab=paste("-log",base,"(expected value)",sep=""), ylab=paste("-log",base,"(observed value)",sep=""),...)
  else z <- qqplot (u,r,xlab="expected value", ylab="observed value",...)
  abline(0,1)
  invisible(z)
}

# 16-8-2007

