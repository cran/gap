# 15-10-03 WAH
# 16-10-03 UCL office to add sort=F in merge
hap.em<-function(id,data,locus.label=NA,converge.eps=0.000001,maxiter=500)
{
  data<-as.matrix(data)
  nloci<-dim(data)[2]/2
  loci<-rep(0,nloci)
  for (i in 1:nloci)
  {
     loci[i]<-max(data[,c(2*i-1,2*i)])
  }
  if(all(is.na(locus.label))) {
     locus.label<- paste("loc-",1:nloci,sep="")
  }
##
  z<- hap(id=id,data=data,nloci,loci=loci,ss=1)
#
  tmp1<-read.table("hap.out",header=T)
# unlink("hap.out")
  haplotype<-as.matrix(tmp1[,1:nloci])
  uhap<-hapid<-1:(dim(tmp1)[1])
  tmp1<-data.frame(tmp1,hapid)
  hap.prob<-tmp1[,nloci+1]
#
  tmp2<-read.table("assign.out",header=T)
# unlink("assign.out")
  nrow<-dim(tmp2)[1]/2
  indx1<-2*1:nrow-1
  indx2<-2*1:nrow
  indx.subj<-tmp2[indx1,1]
  post<-tmp2[indx1,nloci+3]

  tmp<-merge(tmp1[,-(nloci+1)],tmp2[,-c(1,2)],sort=F)
  hap1<-tmp[indx1,(nloci+1)]
  hap2<-tmp[indx2,(nloci+1)]
  nreps<-tapply(indx.subj,indx.subj,length)

  list (lnlike=z$l1,hap.prob=hap.prob,indx.subj=indx.subj,post=post,
        hap1code=hap1,hap2code=hap2,haplotype=haplotype,nreps=nreps,
        converge=z$converge,niter=z$niter,uhap=uhap)
}
