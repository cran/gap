hap.score<-function(y, geno, trait.type="gaussian",
                    offset = NA, x.adj = NA, skip.haplo=.005,
                    locus.label=NA, miss.val=0, n.sim=0, method="gc", id=NA, handle.miss=0, n.miss.loci=NA, sexid=NA)
{
  require(haplo.score)
  trait.int <- charmatch(trait.type, c("gaussian", "binomial", "poisson", "ordinal"))
  if(is.na(trait.int)) stop("Invalid trait type")
  if(trait.int == 0)   stop("Ambiguous trait type")
  if(length(y)!=nrow(geno)) stop("Dims of y and geno are not compatible")
  n.loci <- ncol(geno)/2
  if(n.loci != (floor(ncol(geno)/2)) )stop("Odd number of cols of geno")
  if(handle.miss==0)
  {
    miss <- apply(is.na(geno),1,any)
    if(!all(is.na(miss.val))) {
       for(mval in miss.val){
          miss <- miss | apply(geno==mval, 1, any)
       }
    }
  }
  else
  {
    if(is.na(n.miss.loci)) stop("Maximum number of missing loci (n.miss.loci) not specified")
    nmiss <- apply(is.na(geno),1,sum)
    if(!all(is.na(miss.val))) {
       for(mval in miss.val) {
          nmiss <- nmiss + apply(geno==mval, 1, sum)
       }
    }
    if(n.miss.loci<0 | n.miss.loci >= n.loci) stop("Invalid control for number of missing loci")
    miss <- rep(F, length(y))
    for(i in 1:length(y)) if(nmiss[i] > n.miss.loci*2) miss[i] <- T
  }
  adjusted <- T
  if( all(is.na(x.adj)) ) adjusted <- F
  if(adjusted){
    x.adj <- as.matrix(x.adj)
    if(nrow(x.adj)!=length(y)) stop("Dims of y and x.adj are not compatible")
  }
  miss <- miss | is.na(y) 
  if(adjusted) miss <- miss| apply(is.na(x.adj),1,any)
  if(trait.int==3) {
     if(all(is.na(offset))) stop("Missing offset")
     miss <- miss | is.na(offset)
     offset <- offset[!miss]
  }
  y <- as.numeric(y[!miss])
  geno <- geno[!miss,]
  if(adjusted) x.adj <- x.adj[!miss,,drop=F]
  if(trait.int==2) {
    if(!all(y==1|y==0)) stop("Invalid y values")
    if(all(y==1) | all(y==0)) stop("No variation in y values")
  }
  if(trait.int==4){
     y <- factor(y)
     y.lev <- levels(y)
     y <- as.numeric(y)
     if(max(y) < 3) stop("Less than 3 levels for y values")
  }
  n.subj <- length(y)
  if(all(is.na(id))) id <- 1:n.subj
  method.id<-charmatch(method, c("gc", "hap"))
  if(is.na(method.id)) stop("Invalid selection of method")
  if(method.id == 0)   stop("Ambiguous method")
  else if(method.id==1) haplo <- gc.em(data=geno, locus.label, converge.eps=0.00001, maxiter=5000, handle.miss=handle.miss)
  else haplo <- hap.em(id, data=geno, locus.label, converge.eps=0.00001, maxiter=5000) 
  if(!haplo$converge) stop("EM for haplo failed to converge")
  hap1 <- haplo$hap1code
  hap2 <- haplo$hap2code
  indx <- haplo$indx.subj
  post <- haplo$post
  nreps <- as.vector(haplo$nreps)
  uhap<-haplo$uhap
  which.haplo<-haplo$hap.prob>=skip.haplo
  uhap<-uhap[which.haplo]
  x <- outer(hap1,uhap,"==") + outer(hap2,uhap,"==")
  n.x <- ncol(x)
  x.post<-matrix(rep(NA, n.subj * n.x), ncol=n.x)
  for(j in 1:n.x){
     x.post[,j] <- tapply(x[,j]*post, indx, sum)
  }
  if(trait.int <= 3){ 
    if(!adjusted){
       mu <- switch(trait.int, mean(y), mean(y), sum(y)/sum(offset) )
       a  <- switch(trait.int, var(y), 1, 1)
       x.adj <- matrix(rep(1,n.subj),ncol=1)
     }
     if(adjusted){
        reg.out <- glm(y ~ x.adj, family=trait.type)
        x.adj <- cbind(rep(1,n.subj),x.adj)
        mu <- reg.out$fitted.values
        a  <- switch(trait.int,sum(reg.out$residuals^2)/reg.out$df.residual,1, 1)
      }
     v <- switch(trait.int, 1/a, mu*(1-mu), mu )
     tmp <- haplo.score.glm(y, mu, a, v, x.adj, nreps, x.post, post, x)
     u.score <- tmp$u.score
     v.score <- tmp$v.score
   }
   if(trait.int ==4) {
      if(adjusted){
         library("Design")
         library("Hmisc")
         reg.out <- lrm(y ~ x.adj)
         K <- max(y)
         n.xadj <- ncol(x.adj)
         alpha <- reg.out$coef[1:(K-1)]
         beta <- reg.out$coeff[K:(K-1 + n.xadj)]

         tmp <- haplo.score.podds(y, alpha, beta, x.adj, nreps, x.post,
                              post, x)
       }
      if(!adjusted){
         tbl <- table(y)
         s <- 1- (cumsum(tbl)-tbl)/n.subj
         alpha <-  - log((1-s[-1])/s[-1])
         tmp <- haplo.score.podds(y, alpha, beta=NA, x.adj=NA, nreps, x.post, post, x)
       }
      u.score <- tmp$u.score
      v.score <- tmp$v.score
    }
   tmp<-Ginv(v.score)
   df <- tmp$rank
   g.inv <- tmp$Ginv
   score.global <- u.score%*% g.inv %*%u.score
   score.haplo <- u.score / sqrt(diag(v.score))
   score.max <-  max(score.haplo^2)
   if(n.sim==0){
      score.global.p.sim <- NA
      score.haplo.p.sim <- rep(NA,length(score.haplo))
      score.max.p.sim <- NA
      n.val.global <- NA
      n.val.haplo <- NA
   }
   if(n.sim > 0){
      score.global.rej <- 0
      score.haplo.rej  <- rep(0,length(score.haplo))
      score.max.rej    <- 0
      n.val.global <- 0
      n.val.haplo <- 0
      if(trait.int<=3){
         mu.rand <- mu
         v.rand <- v  
       }
      for(i in 1:n.sim){
         rand.ord <- order(runif(n.subj))
         if(trait.int <=3){ 
           if(adjusted){
              mu.rand <- mu[rand.ord]
              v.rand <- switch(trait.int, v, v[rand.ord], v[rand.ord])
            }
           tmp <- haplo.score.glm(y[rand.ord], mu.rand, a, v.rand, 
                                x.adj[rand.ord,], nreps, x.post, post, x)
         }
         if(trait.int ==4){
            if(adjusted){             
               tmp <- haplo.score.podds(y[rand.ord], alpha, beta, 
                               x.adj[rand.ord,,drop=F],nreps, x.post, post, x)
            }
            if(!adjusted) {
               tmp <- haplo.score.podds(y[rand.ord], alpha, beta=NA, 
                               x.adj=NA,nreps, x.post, post, x)
             }

          }
         u.score <- tmp$u.score
         v.score <- tmp$v.score
         tmp <- Ginv(v.score)
         g.inv <- tmp$Ginv  
         score.global.sim <- u.score %*% g.inv %*% u.score
         score.haplo.sim  <- (u.score / sqrt(diag(v.score)))^2
         score.max.sim <- max(score.haplo.sim)
         if(!is.na(score.global.sim)) {
            n.val.global <- n.val.global +1
            if(score.global.sim >= score.global) score.global.rej <- score.global.rej +1
          }
         if(!any(is.na(score.haplo.sim))){
            n.val.haplo <- n.val.haplo + 1
            score.haplo.rej <- score.haplo.rej +
                               ifelse(score.haplo.sim >= score.haplo^2, 1, 0)
            if(score.max.sim >= score.max) score.max.rej <- score.max.rej +1
          }
      }
      score.global.p.sim <- score.global.rej /  n.val.global
      score.haplo.p.sim <- score.haplo.rej / n.val.haplo
      score.max.p.sim <- score.max.rej / n.val.haplo
    }
   score.global.p <- 1 - pchisq(score.global,df)
   score.haplo.p <- 1-pchisq(score.haplo^2,1)
   if(all(is.na(locus.label))) {
      locus.label<- paste("loc-",1:n.loci,sep="")
    }
   obj <- (list(score.global=score.global, df=df,score.global.p=score.global.p,
       score.global.p.sim=score.global.p.sim,
       score.haplo=score.haplo,score.haplo.p=score.haplo.p,
       score.haplo.p.sim=score.haplo.p.sim,
       score.max.p.sim=score.max.p.sim,
       haplotype=haplo$haplotype[which.haplo,],
       hap.prob=haplo$hap.prob[which.haplo],
       locus.label=locus.label,
       n.sim=n.sim, n.val.global=n.val.global, n.val.haplo=n.val.haplo))
   class(obj) <- "haplo.score"
   return(obj)
}
 
# 13-9-2003 start to implement
# 14-9-2003 in shape
# 21-9-2003 start extensive checking
# 23-9-2003 rewrite interface to genecounting
# 26-9-2003 done with successful use of byand order
# 17-10-2003 start to implement missing genotype code
  

