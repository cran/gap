kin.morgan<-function(ped)
{
   pedsize<-dim(ped)[1]
   kin<-rep(0,pedsize*(pedsize+1)/2)
   z<-.C("kin_morgan",data=as.integer(t(ped)),
          pedsize=as.integer(pedsize),kin=as.double(array(kin)),PACKAGE="gap")

    list(kin=z$kin)
}
