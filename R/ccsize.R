ccsize <- function(n,q,pD,p1,alpha=0.05,theta,power=NULL)
{
   if(p1>1) stop("p1 is in [0,1]")
   if(pD>1) stop("pD is in [0,1]")
   if(q>1) stop("q is in [0,1]")
   p2 <- 1 - p1
   if(is.null(power))
# formula (5)
   {
     s <- p1 * p2 * pD/(q + (1 - q) * pD)
     z_alpha <- qnorm(p=alpha)
     nb <- n * q
     z <- z_alpha + theta * (nb * s)^0.5
     power <- pnorm(q=z)
     invisible(power)
   }
   else
# formula (6)
   {
     bad_nb <- -999
     z1_alpha <- qnorm(p=1 - alpha)
     z_beta <- qnorm(p=power)
     s <- (z1_alpha + z_beta)/(p1 * p2 * pD)^0.5
     if(theta <= s * ((1-pD) / n)^0.5) 
     {
        nb <- bad_nb
        cat("bad exp(theta)=",exp(theta),"\n")
     }
     else 
     {
       B <- s^2 / theta^2
       nb <- ceiling(n * B * pD / (n - B * (1 - pD)))
       if(nb > n)
       {
          cat("bad subcohort size", nb, "exp(theta)=",exp(theta),"\n")
          nb <- bad_nb
       }
     }
     invisible(nb)
   }
}
