kin.morgan<-function(ped)
{
   vector.to.kin <- function (pars) 
   # 5/6/2004
   # this function is in spirit similar to g2a and make.del(mvnmle)
   {
       k <- floor((-1 + sqrt(1 + 8 * length(pars)))/2)
       mymatrix <- diag(1:k)
       if (k > 1) {
           for (i in 1:k) {
               mymatrix[i, 1:i] <- mymatrix[1:i,i] <- pars[1:i]
               pars <- pars[-(1:i)]
           }
       }
       mymatrix
   }
   pedsize<-dim(ped)[1]
   kin<-rep(0,pedsize*(pedsize+1)/2)
   z<-.C("kin_morgan",data=as.integer(t(ped)),
          pedsize=as.integer(pedsize),kin=as.double(array(kin)),PACKAGE="gap")

   list(kin=z$kin,kin.matrix=vector.to.kin(z$kin))
}
