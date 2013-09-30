\name{l51}
\docType{data}
\alias{l51}
\title{An example pedigree data}
\description{
The data contains data on 51 individuals in a pedigree. Below it is used for comparing results from various packages.
}
\usage{data(l51)}
\format{A data frame}
\source{Morgan v3.}
\references{
Morgan v3. \url{http://www.stat.washington.edu/thompson/Genepi/MORGAN/Morgan.shtml}
} 
\examples{\dontrun{
library(kinship2)
library(regress)
library(coxme)
library(pls)
library(POET)
library(pedigreemm)
N <- dim(l51)[1]
ped51 <- within(l51, {
     hist(qt)
     pid <- "51"
     w <- quantile(qt,probs=0.75,na.rm=TRUE)
     r <- rnorm(N)
     bt=ifelse(qt<=w,0,1)
     }
)
ped <- with(ped51,kinship2::pedigree(id,fid,mid,sex,aff=bt))
plot(ped)
kmat <- with(ped51,kinship(id, fid, mid))
k2 <- as.matrix(2*kmat)
I <- diag(N)
dimnames(I) <- dimnames(k2)
k3 <- k2+I

regress(qt ~ 1, ~I, data=ped51)
regress(qt ~ 1, ~k2, data=ped51)
regress(qt ~ 1, ~k2+I, data=ped51)
regress(qt ~ 1, ~k3, data=ped51)

lmekin(qt ~ 1 + (1|id), data=ped51)
lmekin(qt ~ 1 + (1|id), data=ped51, varlist=list(k2))
lmekin(qt ~ 1 + (1|id), data=ped51, varlist=list(k2,I))
lmekin(qt ~ 1 + (1|id), data=ped51, varlist=list(k3))

coxme(Surv(qt,bt)~1+(1|id),data=ped51)
coxme(Surv(qt,bt)~1+(1|id),data=ped51,varlist=list(k2))
coxme(Surv(qt,bt)~1+(1|id),data=ped51,varlist=list(k3))
ped <- with(ped51,pedigreemm::pedigree(fid,mid,id))
U <- relfactor(ped)
A <- as.matrix(t(U) %*% U)
sum(A-k2)
F <- inbreeding(ped)
pedigreemm(bt~(1|id),data=ped51, family="binomial"(link="logit"),
           pedigree=list(id=ped))
pedigreemm(bt~(1|id),data=ped51, family="binomial"(link="logit"),
           REML=FALSE, pedigree=list(id=ped))

with(ped51,summary(qt))
m0 <- lm(qt ~ 1, data=ped51)
summary(m0)
m1 <- lm(qt ~ r, data=ped51)
summary(m1)
rawdata <- function(){
  pca <- prcomp(k2)
  screeplot(pca,25)
  biplot(pca)
  pred <- predict(pca)
  l1 <- lm(qt~1+pred[,1:6],data=ped51)
  summary(l1)
  l2 <- lm(qt~1+pred[,1:35],data=ped51)
  summary(l2)
  
  p0 <- pcr(qt ~ pred, 25, data=ped51, validation="CV")
  summary(p0)
  coef(p0)
  coefplot(p0)
  p1 <- pcr(qt ~ k2, 35, data=ped51, validation="CV")
  summary(p1)
  coef(p0)
  coefplot(p1)
  scores(p1)
  poet <- POET(k2)
  pred <- k2%*%poet$loadings
  t0 <- lm(qt~pred,data=ped51)
  summary(t0)
  poet <- POET(k2,4)
  pred <- k2%*%poet$loadings
  t1 <- lm(qt~pred,data=ped51)
  summary(t1)
}

# GCTA

k2l <- matrix(0,N*(N+1)/2,4)
p <- matrix(0,N,4)
k <- 1
for(i in 1:N)
{
  p[i,] <- with(ped51[i,],c(i,i,qt,bt))
  for(j in 1:i) 
  {
    k2l[k,] <- c(i,j,51,k2[i,j])
    k <- k+1
  }
}
write(t(p),file="51.txt",4,sep="\t")
write(t(p[,-(3:4)]),file="51.grm.id",2,sep="\t")
write(t(k2l),file="51.grm",4,sep="\t")

# gzip -f 51.grm
# gcta64 --grm 51 --pca 51 --out 51
# gcta64 --reml --grm 51 --pheno 51.txt --out 51_qt
# gcta64 --reml --grm 51 --pheno 51.txt --mpheno 2 --prevalence 0.25 --out 51_bt_25
# gcta64 --reml --grm 51 --pheno 51.txt --mpheno 2 --out 51_bt

# SOLAR

solar_ped <- ped51[,c("id","fid","mid","sex")]
write.table(solar_ped,file="51.ped",col.names=c("id","fa","mo","sex"),
          quote=FALSE, row.names=FALSE,sep=",",na="")
solar_pheno <- ped51[,c("id","qt","bt")]
write.table(solar_pheno,file="51.phen",col.names=c("id","qt","bt"),
          quote=FALSE, row.names=FALSE,sep=",",na="")

# load pedigree 51.ped
# load phenotypes 51.phen
# model new
# trait qt
# polygenic -screen
# trait qt
# covariate sex
# polygenic -screen -fix sex

# The following is from experiments with GCTA on heritability estimates
# involving a binary trait
# 13-9-2013 MRC-Epid JHZ

VG <- 0.017974; VVG <- 0.003988^2
VGE <- 0.002451; VVGE <- 0.005247^2
Ve <- 0.198894; VVe <- 0.005764^2

cVGVGE <- -7.93348e-06
cVGVe <- -5.54006e-06
cVGEVe <- -1.95297e-05

Vp <- VG + VGE + Ve
s <- VVG + VVGE + VVe
s12 <- 2*cVGVGE
s13 <- 2*cVGVe
s23 <- 2*cVGEVe
VVp <- s + s12 + s13 + s23
cVpVG <- VVG + cVGVGE + cVGVe
cVpVGE <- cVGVGE + VVGE + cVGEVe

Vp
sqrt(VVp)
VG/Vp
sqrt(VR(VG, VVG, Vp, VVp, cVpVG))
VGE/Vp
sqrt(VR(VGE,VVGE,Vp,VVp, cVpVGE))

K <- 0.05
x <- qnorm(1-K)
z <- dnorm(x)
1/sqrt(2*pi)*exp(-x^2/2)
P <- 0.496404
fK <- (K*(1-K)/z)^2
fP <- P*(1-P)
f <- fK/fP
ho <- 0.274553
vo <- 0.067531
f*ho
f*vo
hl <- 0.232958
vl <- 0.057300
r1 <- hl/ho
r2 <- vl/vo
r1==r2
z2 <- K^2*(1-K)^2/(f*fP)
x2 <- -log(2*pi*z2)
sqrt(x2)

V <- c(0.017974, 0.002451, 0.198894)
VCOV <- matrix(0,3,3)
diag(VCOV) <- c(0.003988, 0.005247, 0.005764)^2
VCOV[2,1] <- -7.93348e-06
VCOV[3,1] <- -5.54006e-06
VCOV[3,2] <- -1.95297e-05
z <- h2GE(V,VCOV)

h2 <- 0.274553
se <- 0.067531
P <- 0.496404
z <- h2l(P=P,h2=h2,se=se)

hl <- 0.232958
vl <- 0.057300

R <- 50
kk <- h2all <- seall <- rep(0,R)
for(i in 1:R)
{
  kk[i] <- 0.4*i/R
  z <- h2l(kk[i],P=P,h2=h2,se=se)
  h2all[i] <- z$h2l
  seall[i] <- z$se
}

tiff("f.tiff",width=0.03937*189,height=0.03937*189/1.5,units="in",res=1200,compress="lzw")
plot(kk,h2all,type="l",ylab="Adjusted heritability",xlab="Prevalence")
lines(kk,h2all-seall,lty="dashed")
lines(kk,h2all+seall,lty="dashed")
title("Adjusted heritability for observed value of .23 and cases% = 0.495")
dev.off()
}
}
\keyword{datasets}
