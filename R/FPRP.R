FPRP <- function(a,b,pi0,ORlist,logscale=FALSE)
{
  if (logscale) {
     thetahat <- a
     V <- b
  } else {
     thetahat <- log(a)
     V <- ((log(b)-log(a))/1.96)^2
  }
  if (thetahat>0) {
     p <- 2*(1-pnorm(thetahat/sqrt(V)))
     q <- qnorm(1-p/2)
     b <- 1-pnorm((q*sqrt(V)-log(ORlist))/sqrt(V))
  } else {
     p <- 2*pnorm(thetahat/sqrt(V))
     q <- qnorm(1-p/2)
     b <- pnorm((-q*sqrt(V)-log(ORlist))/sqrt(V))
  }
  FPRP <- t(p*(1-pi0)/(p*(1-pi0)+b%o%pi0))
  row.names(FPRP) <- pi0
  colnames(FPRP) <- ORlist
  invisible(list(p=p,power=b,FPRP=FPRP))
}

# 28-7-2007
