qqunif <- function(u,logscale=TRUE,base=10,col=palette()[4],lcol=palette()[2],...)
{
  u <- u[!is.na(u)]
  n <- length(u)
  r <- (1:n)/(n+1)
  y <- -log(u,base)
  x <- -log(r,base)
  if (logscale) z <- qqplot (x,y,xlab=paste("-log",base,"(expected value)",sep=""),
     ylab=paste("-log",base,"(observed value)",sep=""),col=col,...)
  else z <- qqplot (r,u,xlab="expected value", ylab="observed value",...)
  abline(0,1,col=lcol)
  invisible(z)
}

# 9-4-2008
