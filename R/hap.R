hap<-function(id,data,nloci,loci=rep(2,nloci),names=paste("loci",1:nloci,sep=""),
              mb=0,pr=0,po=0.001,to=0.001,th=1,maxit=100,n=0,
              ss=0,rs=0,rp=0,ro=0,rv=0,sd=0,mm=0,mi=0,mc=50,ds=0.1,de=0,q=0)
{
  if (rv & ro) stop("rv and ro flags cannot both be set\n");
# if (mi==0 & (mc | ds | de)) stop("mc, ds, de parameters are only legal if mi is set\n");
  if (rp & mm==0) stop("rp option only relevant with mm # option\n");
  nobs<-dim(data)[1]
  data<-as.matrix(data)
  if(length(id)!=dim(data)[1]) stop("id and data should have the same length")
  l1<-niter<-converge<-0

  z<-.C("hap",nobs=as.integer(nobs),idstr=as.character(id),data=as.character(t(data)),
        nloci=as.integer(nloci),loci=as.integer(loci),names=as.character(names),mb=as.double(mb),
        pr=as.double(pr),po=as.double(po),to=as.double(to),th=as.double(th),
        maxitt=as.double(maxit),n=as.integer(n),sst=as.integer(ss),rst=as.integer(rs),
        rp=as.integer(rp),ro=as.integer(ro),rv=as.integer(rv),sd=as.double(sd),
        mm=as.integer(mm),mi=as.integer(mi),mc=as.integer(mc),ds=as.double(ds),
        de=as.double(de),q=as.integer(q),l1=as.double(l1),niter=as.integer(niter),
        converged=as.integer(converge))

  list(l1=z$l1,converge=z$converged,niter=z$niter)
}
