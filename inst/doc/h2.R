### R code from vignette source 'h2.Rnw'

###################################################
### code chunk number 1: h2.Rnw:42-49
###################################################
library(gap)
head(l51,10)
library(kinship2)
ped <- with(l51,pedigree(id,fid,mid,sex))
pdf("figures/l51.pdf")
plot(ped)
dev.off()


###################################################
### code chunk number 2: h2.Rnw:57-64
###################################################
library(gap)
k2 <- kin.morgan(l51)$kin.matrix*2
k2[1:10,1:10]
library(regress)
r <- regress(qt ~ 1, ~k2, data=l51)
r$sigma
r$sigma.cov


###################################################
### code chunk number 3: h2.Rnw:70-77
###################################################
N <- dim(l51)[1]
w <- with(l51,quantile(qt,probs=0.75,na.rm=TRUE))
l51 <- within(l51, bt <- ifelse(qt<=w,0,1))
with(l51,table(bt))
d <- regress(bt ~ 1, ~k2, data=l51)
d$sigma
d$sigma.cov


###################################################
### code chunk number 4: h2.Rnw:88-101
###################################################
library(gap)
# qt
sigma <- c(0.2817099, 0.4444962)
sigma.cov <- matrix(
c(0.07163300, -0.03991478,
-0.03991478, 0.04042731), 2, 2) 
h2G(sigma,sigma.cov)
# bt
sigma <- c(0.0307703, 0.1678370)
sigma.cov <- matrix(
c(0.003615481, -0.002525622,
 -0.002525622,  0.003492826), 2, 2)
h2G(sigma,sigma.cov)


###################################################
### code chunk number 5: h2.Rnw:108-109
###################################################
h2l(K=0.25, P=11/51, h2=0.1549304, se=0.2904298)


###################################################
### code chunk number 6: h2.Rnw:121-129
###################################################
library(gap)
V <- c(0.017974, 0.002451, 0.198894)
VCOV <- matrix(0,3,3)
diag(VCOV) <- c(0.003988, 0.005247, 0.005764)^2
VCOV[2,1] <- -7.93348e-06
VCOV[3,1] <- -5.54006e-06
VCOV[3,2] <- -1.95297e-05
z <- h2GE(V,VCOV)


###################################################
### code chunk number 7: h2.Rnw:133-173
###################################################
library(gap)
h2 <- 0.274553
se <- 0.067531
P <- 0.496404

z <- h2l(P=P,h2=h2,se=se)

R <- 50
kk <- h2all <- seall <- rep(0,R)
for(i in 1:R)
{
  kk[i] <- 0.4*i/R
  z <- h2l(kk[i],P=P,h2=h2,se=se,verbose=FALSE)
  h2all[i] <- z$h2l
  seall[i] <- z$se
}

h2 <- 0.044
se <- 0.061
z <- h2l(P=P,h2=h2,se=se)

h2alls <- sealls <- rep(0,R)
for(i in 1:R)
{
  z <- h2l(kk[i],P=P,h2=h2,se=se,verbose=FALSE)
  h2alls[i] <- z$h2l
  sealls[i] <- z$se 
}

pdf("figures/h2l.pdf")
par(mfrow=c(1,2))
plot(kk,h2all,type="l",ylab="Adjusted heritability",xlab="Prevalence")
lines(kk,h2all-seall,lty="dashed")
lines(kk,h2all+seall,lty="dashed")
title("(a) h2 = .274 and cases% = 50")
plot(kk,h2alls,type="l",ylab="Adjusted heritability",xlab="Prevalence",ylim=c(0,0.15))
lines(kk,h2alls-sealls,lty="dashed")
lines(kk,h2alls+sealls,lty="dashed")
title("(b) h2 = .044 and cases% = 50")
dev.off()


###################################################
### code chunk number 8: h2.Rnw:188-198 (eval = FALSE)
###################################################
## p <- matrix(0,N,4)
## for(i in 1:N) p[i,] <- with(l51[i,],c(i,i,qt,bt))
## write(t(p),file="51.txt",4,sep="\t")
## NN <- rep(51, N * (N + 1)/2)
## WriteGRM(51,p[,1:2],NN,k2)
## one <- ReadGRM(51)
## grm <- one$grm
## WriteGRMBin(51,grm,NN,p[,1:2])
## two <- ReadGRMBin(51,TRUE)
## sum(one$GRM-two$GRM)


###################################################
### code chunk number 9: h2.Rnw:209-211
###################################################
library(foreign)
write.dta(l51, "l51.dta")


