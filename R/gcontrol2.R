gcontrol2 <- function(p,col=palette()[4],lcol=palette()[2],...)
{
   p <- p[!is.na(p)]
   n <- length(p)
   x2obs <- qchisq(p,1,lower.tail=FALSE)
   x2exp <- qchisq((1:n)/n,1,lower.tail=FALSE)
   lambda <- median(x2obs)/median(x2exp)
   qqplot(x2exp,x2obs,xlab="expected value",ylab="observed value",col=col,...)
   abline(0,1,col=lcol)
   invisible(list(x=x2exp,y=x2obs,lambda=lambda))
}
