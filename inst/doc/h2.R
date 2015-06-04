### R code from vignette source 'h2.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: h2.Rnw:45-52
###################################################
library(gap)
head(l51,10)
library(kinship2)
ped <- with(l51,pedigree(id,fid,mid,sex))
pdf("figures/l51.pdf")
plot(ped)
dev.off()


###################################################
### code chunk number 2: h2.Rnw:60-67
###################################################
library(gap)
k2 <- kin.morgan(l51)$kin.matrix*2
k2[1:10,1:10]
library(regress)
r <- regress(qt ~ 1, ~k2, data=l51)
r$sigma
r$sigma.cov


###################################################
### code chunk number 3: h2.Rnw:72-79
###################################################
N <- dim(l51)[1]
w <- with(l51,quantile(qt,probs=0.75,na.rm=TRUE))
l51 <- within(l51, bt <- ifelse(qt<=w,0,1))
with(l51,table(bt))
d <- regress(bt ~ 1, ~k2, data=l51)
d$sigma
d$sigma.cov


###################################################
### code chunk number 4: h2.Rnw:90-103
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
### code chunk number 5: h2.Rnw:110-111
###################################################
h2l(K=0.25, P=11/51, h2=0.1549304, se=0.2904298)


###################################################
### code chunk number 6: h2.Rnw:123-131
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
### code chunk number 7: h2.Rnw:138-168
###################################################
library(gap)
P <- 0.496404
R <- 50
kk <- h2all <- seall <- h2alls <- sealls <- rep(0,R)
for(i in 1:R)
{
  kk[i] <- i/R
  h2 <- 0.274553
  se <- 0.067531
  z <- h2l(kk[i],P=P,h2=h2,se=se,verbose=FALSE)
  h2all[i] <- z$h2l
  seall[i] <- z$se
  h2 <- 0.044
  se <- 0.061
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
### code chunk number 8: h2.Rnw:180-190 (eval = FALSE)
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
### code chunk number 9: h2.Rnw:205-207
###################################################
library(foreign)
write.dta(l51, "l51.dta")


###################################################
### code chunk number 10: h2.Rnw:226-234
###################################################
library(gap)
library(MCMCglmm)
prior<-list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))
m <- MCMCglmm(qt~sex,random=~id,data=l51,prior=prior,burnin=10000,nitt=100000,verbose=FALSE)
summary(m)
pdf("figures/MCMCglmm1.pdf")
plot(m)
dev.off()


###################################################
### code chunk number 11: h2.Rnw:240-263 (eval = FALSE)
###################################################
## library(gap)
## km <- kin.morgan(l51)
## k2 <- km$kin.matrix*2
## 
## N <- 51
## i <- rep(1:N,rep(N,N))
## j <- rep(1:N,N)
## 
## library(Matrix)
## s <-spMatrix(N,N,i,j,as.vector(k2))
## Ginv<-solve(s)
## class(Ginv) <- "dgCMatrix"
## rownames(Ginv) <- Ginv@Dimnames[[1]] <- with(l51,id)
## 
## library(MCMCglmm)
## prior<-list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))
## m <- MCMCglmm(qt~1, random=~id, ginverse=list(id=Ginv), data=l51, prior=prior, 
##               burnin=10000, nitt=100000, verbose=FALSE)
## summary(m)
## save(m,file="MCMCglmm.fit")
## pdf("MCMCglmm2.pdf")
## plot(m$VCV)
## dev.off()


###################################################
### code chunk number 12: h2.Rnw:300-306 (eval = FALSE)
###################################################
## s <- kin.morgan(l51)
## K <- with(s,kin.matrix*2)
## prior <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))
## m <- MCMCgrm(qt~1,prior,l51,K,n.burnin=10000, n.iter=100000)
## save(m,file="l51.m")
## plot(m)


###################################################
### code chunk number 13: h2.Rnw:317-348
###################################################
meyer <- within(meyer,{
   g1 <- ifelse(generation==1,1,0)
   g2 <- ifelse(generation==2,1,0)
})
# library(kinship)
# A <- with(meyer,kinship(animal,sire,dam))*2
# Here we convert NAs to 0s to be compatible with kin.morgan
meyer0 <- within(meyer,{
   id <- animal
   animal <- ifelse(!is.na(animal),animal,0)
   dam <- ifelse(!is.na(dam),dam,0)
   sire <- ifelse(!is.na(sire),sire,0)
   g1 <- ifelse(generation==1,1,0)
   g2 <- ifelse(generation==2,1,0)
})
A <- kin.morgan(meyer0)$kin.matrix*2

library(regress)
r <- regress(y~-1+g1+g2,~A,data=meyer0)
summary(r)
with(r,h2G(sigma,sigma.cov))

library(MCMCglmm)
m <-MCMCglmm(y~-1+g1+g2,random=animal~1,pedigree=meyer[,1:3],data=meyer,verbose=FALSE)
summary(m)
plot(m)

prior <- list(R=list(V=1, nu=0.002), G=list(G1=list(V=1, nu=0.002)))
m2 <- MCMCgrm(y~-1+g1+g2,prior,meyer0,A,singular.ok=TRUE,verbose=FALSE)
summary(m2)
plot(m2)


###################################################
### code chunk number 14: h2.Rnw:439-488
###################################################
library(gap)
set.seed(1234567)
km <- kin.morgan(l51)
k2 <- km$kin.matrix*2
l51 <-
within(l51, {qt[is.na(qt)] <- rnorm(length(qt[is.na(qt)]),
             mean(qt,na.rm=TRUE),sd(qt,na.rm=TRUE))})
N <- dim(l51)[1]
data=with(l51,list(N=N,qt=qt,sex=sex,GI=solve(k2),u=rep(0,N)))
library(regress)
r <- regress(qt ~ sex, ~k2, data=data)
r
with(r,{
  print(sqrt(sigma+1.96*sqrt(diag(sigma.cov))))
  h2G(sigma,sigma.cov)
})
inits=function()list(b1=0,b2=0,sigma.p=0.001,sigma.r=0.001)
modelfile=function() {
    b1 ~ dnorm(0, 0.001)
    b2 ~ dnorm(0, 0.001)
    sigma.p ~ dunif(0,0.9)
    sigma.r ~ dunif(0,0.9)
    p <- pow(sigma.p, 2)
    r <- pow(sigma.r, 2)
    h2 <- p / (p + r)
    tau <- pow(sigma.r, -2)
    xi ~ dnorm(0,tau.xi)
    tau.xi <- pow(0.9,-2)
    g[1:N] ~ dmnorm(u[],GI[,]/p)   
    for(i in 1:N) {qt[i] ~ dnorm(b1 + b2 * sex[i] + xi*g[i],tau)}
}
library(R2jags)
jagsfit <- jags(data,
                inits,
                parameters.to.save=c("b1","b2","p","r","h2"),
                model.file=modelfile,
                n.chains=3,
                n.burnin=1000,
                n.iter=10000)
save(jagsfit,file="jags.fit")
print(jagsfit)
pdf("figures/jags.pdf")
plot(jagsfit)
library(lattice)
jagsfit.mcmc <- as.mcmc(jagsfit)
traceplot(jagsfit.mcmc)
xyplot(jagsfit.mcmc)
densityplot(jagsfit.mcmc)
dev.off()


###################################################
### code chunk number 15: h2.Rnw:508-540
###################################################
library(gap)
set.seed(1234567)
ped51 <- 
within(l51, {
  qt[is.na(qt)] <- rnorm(length(qt[is.na(qt)]),mean(qt,na.rm=TRUE),sd(qt,na.rm=TRUE))
})
km <- kin.morgan(l51)
k2 <- km$kin.matrix*2
N <- dim(l51)[1]
data=with(ped51,list(N=N,qt=qt,sex=sex,G=k2,I=diag(N),alpha=1,gamma=1))
inits=function()list(b1=0,b2=0,h2=0.4)
modelfile=function() {
    h2 ~ dunif(0,1)
    Omega <- inverse((h2*G[,] + (1-h2)*I[,])*gamma/alpha)
    qt[1:N] ~ dmt(mu[],Omega[,],2*alpha)
    mu[1:N] <- b1 + b2 * sex[]
    b1 ~ dnorm(0, 0.001)
    b2 ~ dnorm(0, 0.001)
}
library(R2jags)
jagsfit <- jags(data,inits,parameters.to.save=c("b1","b2","h2"),
                model.file=modelfile, n.chains=3, n.burnin=1000, n.iter=10000)
save(jagsfit,file="h2.fit")
print(jagsfit)
pdf("figures/h2.pdf")
plot(jagsfit)
library(lattice)
jagsfit.mcmc <- as.mcmc(jagsfit)
traceplot(jagsfit.mcmc)
xyplot(jagsfit.mcmc)
densityplot(jagsfit.mcmc)
dev.off()


