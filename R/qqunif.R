qqunif <- function(u,base=10,...)
{
  u <- u[!is.na(u)]
  n <- length(u)
  r <- (1:n)/(n+1)
  y <- -log(u,base)
  x <- -log(r,base) 
  z <- qqplot (u,r,xlab="-log(expected P value)", ylab="-log(observed P value",...)
  abline(0,1)
  invisible(z)
}

# 28-7-2007

