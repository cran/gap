# worked 28/6/03
tbyt<-function(h,n) {
   D<-VarD<-Dmax<-VarDmax<-Dprime<-VarDprime<-x2<-0
   z<-.C("tbyt",h=as.vector(h), haplotypes=as.double(n),
          D=as.double(D), VarD=as.double(VarD),
          Dmax=as.double(Dmax), VarDmax=as.double(VarDmax),
          Dprime=as.double(Dprime), VarDprime=as.double(VarDprime),
          x2=as.double(x2))

    list(h=h,n=n,D=z$D,VarD=z$VarD,
         Dmax=z$Dmax,VarDmax=z$VarDmax,Dprime=z$Dprime,
         VarDprime=z$VarDprime,x2=z$x2)
}

kbyl<-function(n1,n2,h,n,optrho=2)
{
   x2<-seX2<-rho<-seR<-klinfo<-0
   VarDp<-0
   Dijtable<-Dmaxtable<-Dijptable<-VarDijtable<-VarDijptable<-matrix(rep(0,n1*n2),nrow=n1)
   z<-.C("kbyl",nalleles1=as.integer(n1), nalleles2=as.integer(n2),
          h=as.double(h), haplotypes=as.double(n),
          VarDp=as.double(VarDp),Dijtable=matrix(Dijtable,nrow=n1),
          Dmaxtable=matrix(Dmaxtable,nrow=n1),
          Dijptable=matrix(Dijptable,nrow=n1),
          VarDijtable=matrix(VarDijtable,nrow=n1),
          VarDijptable=matrix(VarDijptable,nrow=n1),
          x2=as.double(x2), seX2=as.double(seX2),
          rho=as.double(rho), seR=as.double(seR), optrho=as.integer(optrho),
          klinfo=as.double(klinfo))

   list(n1=z$nalleles1, n2=z$nalleles2, h=z$h, n=z$haplotypes,
   VarDp=z$VarDp,Dijtable=z$Dijtable, Dmaxtable=z$Dmaxtable,
   Dijptable=z$Dijptable, VarDijtable=z$VarDijtable, VarDijptable=z$VarDijptable,
   x2=z$x2, seX2=z$seX2,
   rho=z$rho, seR=z$seR, optrho=z$optrho, klinfo=z$klinfo)
}
