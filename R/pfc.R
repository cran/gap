pfc <- function(famdata,enum=0)
{
  famsize<-dim(famdata)[1]
  p<-tailp<-sump<-1.0
  stat<-rep(0,20)
  nenum<-0
  z<-.Fortran("family",famdata=as.integer(matrix(famdata,ncol=3)),
               famsize=as.integer(famsize),p=as.double(p),,stat=as.double(stat),toenum=as.integer(enum),
               tailp=as.double(tailp),sump=as.double(sump),nenum=as.double(nenum))

  if(enum==0) list(p=z$p,stat=z$stat[1:5])
  else list(p=z$p,stat=z$stat[1:5],tailp=z$tailp,sump=z$sump,nenum=z$nenum)
}
