hwe.hardy<-function(a,alleles=3,seed=3000,sample=c(1000,1000,5000))
{
  if (alleles<3) stop("number of alleles should be at least 3")
  p<-1.0
  se<-0.0
  swp<-rep(0,3)
  z<-.C("hwe_hardy",a=as.integer(a), alleles=as.integer(alleles), 
        seed=as.integer(seed), gss=as.integer(sample),
        p=as.double(p), se=as.double(se), swp=as.double(swp))

  list(p=z$p, se=z$se, swp=z$swp)
}
