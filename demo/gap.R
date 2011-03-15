## Created on 9-1-2004 (GMT)
## Change from cat() to matrix following example from gif()
##
## fisher.test
#
# This will be suitable for analysis of population-based case-control data with single multiallelic marker
# library(ctest)
# ? fisher.test

## GENECOUNTING
#
# This is the orginal example used in EHPLUS (Zhao et al. 2000 Hum Hered), based on Dr T Li's COMT data
#
ex<-c(
1003.0, 1,  1,2,  1,2,
1003.1, 0,  1,2,  1,2,
1003.2, 0,  1,2,  1,1,
1005.0, 1,  1,1,  1,2,
1005.1, 0,  1,1,  1,2,
1005.2, 0,  1,2,  1,1,
1006.0, 1,  2,2,  2,2,
1006.1, 0,  1,2,  2,2,
1006.2, 0,  2,2,  1,2,
1007.0, 1,  1,1,  1,2,
1007.1, 0,  1,1,  1,2,
1007.2, 0,  1,2,  1,1)

ex<-matrix(ex,ncol=6,byrow=TRUE)
ex
ex.gc<-genecounting(ex[,3:6])
summary(ex.gc)

## HAPLOTYPE TREND REGRESSION
#
# This is an example extracted from HTR program by Dr D Zaykin
#
filespec <- file.path(.path.package("gap"),"tests/htr.R")
source(filespec)
htr.test2<-htr(y,x)
htr.test2
htr.test2<-htr(y,x,n.sim=10)
htr.test2

#
# Generate .dot file
#
cat(
"1 a x x m y",
"1 b x x f y",
"1 c a b m n",
"1 d a b m n",
"1 e a b f y",sep="\n", file="1.ped")
ped <- read.table("1.ped", as.is=T)
pedtodot(ped,sink=F)
unlink("1.ped")

#
# To produce a pedigree diagram
#
filespec <- file.path(.path.package("gap"),"tests/ped.1.3.pre")
pre <- read.table(filespec,as.is=TRUE)
pre
pedtodot(pre,dir="forward")

#
# Q-Q Plot for 1000 U(0,1) r.v., marking those <= 1e-5
#
u_obs <- runif(1000)
r <- qqunif(u_obs,pch=21,bg="blue",bty="n")
u_exp <- r$y
hits <- u_exp >= 2.30103
points(r$x[hits],u_exp[hits],pch=21,bg="green")

#
# To produce a Manhattan plot
#
affy <-c(40220, 41400, 33801, 32334, 32056, 31470, 25835, 27457, 22864, 28501, 26273, 
         24954, 19188, 15721, 14356, 15309, 11281, 14881, 6399, 12400, 7125, 6207)
CM <- cumsum(affy)
n.markers <- sum(affy)
n.chr <- length(affy)
test <- data.frame(chr=rep(1:n.chr,affy),pos=1:n.markers,p=runif(n.markers))
par(las="2",cex=0.6)
colors <- rep(c("blue","green"),11)
mhtplot(test,color=colors,labels=paste(1:n.chr,sep=""),pch=21,bg=colors)
title("A simulated example according to EPIC-Norfolk QCed SNPs")

#
# To produce a regional association plot
#
data(CDKN)
asplot("rs10811661", "CDKN2A/CDKN2B region", "9", CDKNlocus, 9, 5.4e-8, 23, CDKNmap, CDKNgenes)
