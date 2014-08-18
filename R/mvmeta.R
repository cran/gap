#19-8-2009 MRC-Epid JHZ

mvmeta <- function(b,V)
{
   library(magic)
   n1 <- dim(b)[1]
   n2 <- dim(b)[2]
   d <- as.vector(t(b))
   d <- d[!is.na(d)]
   X <- vector()
   Vi <- matrix(NA,n2,n2)
   for (i in 1:n1)
   {
      dz <- !is.na(b[i,])
      dx <- diag(ifelse(dz,1,NA))
      X <- rbind(X,dx[dz,])
      Vi[upper.tri(Vi,diag=TRUE)] <- V[i,]
      if (i==1) Psi <- magic::adiag(Vi[dz,dz])
      else Psi <- magic::adiag(Psi,Vi[dz,dz])
   }
   dl <- length(d)
   for (i in 1:dl) Psi[i:dl,i] <- Psi[i,i:dl]
#  Psi <- replace(Psi,is.na(Psi),0)
   require(MASS)
   cpd <- t(X) %*% MASS::ginv(Psi) %*% X
   beta <- MASS::ginv(cpd) %*% t(X) %*% MASS::ginv(Psi) %*% d
   cov.beta <- MASS::ginv(cpd)
   X2 <- t(d) %*% MASS::ginv(Psi) %*% d - t(beta) %*% cpd %*% beta
   df <- length(d)-n2
   p <- 1-pchisq(X2,df)

   invisible (list(d=d,Psi=Psi,X=X,beta=beta,cov.beta=cov.beta,X2=X2,df=df,p=p))
}
