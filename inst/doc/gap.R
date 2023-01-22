## ----setup, include=FALSE-----------------------------------------------------
set.seed(0)
options(rmarkdown.html_vignette.check_title = FALSE)
knitr::opts_chunk$set(
  out.extra = 'style="display:block; margin: auto"',
  fig.align = "center",
  fig.path = "./",
  collapse = TRUE,
  comment = "#>",
  dev = "png")

## ---- echo=FALSE--------------------------------------------------------------
p1 <- "
 1  2   3  2  2  7/7  7/10
 2  0   0  1  1  -/-  -/-
 3  0   0  2  2  7/9  3/10
 4  2   3  2  2  7/9  3/7
 5  2   3  2  1  7/7  7/10
 6  2   3  1  1  7/7  7/10
 7  2   3  2  1  7/7  7/10
 8  0   0  1  1  -/-  -/-
 9  8   4  1  1  7/9  3/10
10  0   0  2  1  -/-  -/-
11  2  10  2  1  7/7  7/7
12  2  10  2  2  6/7  7/7
13  0   0  1  1  -/-  -/-
14 13  11  1  1  7/8  7/8
15  0   0  1  1  -/-  -/-
16 15  12  2  1  6/6  7/7
"

p2 <- as.data.frame(scan(file=textConnection(p1),what=list(0,0,0,0,0,"","")))
names(p2) <-c("id","fid","mid","sex","aff","GABRB1","D4S1645")
p3 <- data.frame(pid=10081,p2)

## ----pedtodot, fig.cap="An example pedigree", fig.height=8, fig.width=7-------
library(gap)
knitr::kable(p3,caption="An example pedigree")
library(DOT)
# one can see the diagram in RStudio
pedtodot_verbatim(p3,run=TRUE,toDOT=TRUE,return="verbatim")
library(DiagrammeR)
pedtodot_verbatim(p3)
grViz(readr::read_file("10081.dot"))

## ----lukas, fig.cap="A pedigree diagram", fig.height=8, fig.width=7-----------
# pedigree diagram
data(lukas, package="gap.datasets")
library(kinship2)
ped <- with(lukas,pedigree(id,father,mother,sex))
plot(ped,cex=0.4)

## -----------------------------------------------------------------------------
# unordered individuals
gk1 <- kin.morgan(lukas)
write.table(gk1$kin.matrix,"gap_1.txt",quote=FALSE)

library(kinship2)
kk1 <- kinship(lukas[,1],lukas[,2],lukas[,3])
write.table(kk1,"kinship_1.txt",quote=FALSE)

d <- gk1$kin.matrix-kk1
sum(abs(d))

# order individuals so that parents precede their children
library(pedigree)
op <- orderPed(lukas)
olukas <- lukas[order(op),]
gk2 <- kin.morgan(olukas)

write.table(olukas,"olukas.csv",quote=FALSE)
write.table(gk2$kin.matrix,"gap_2.txt",quote=FALSE)

kk2 <- kinship(olukas[,1],olukas[,2],olukas[,3])
write.table(kk2,"kinship_2.txt",quote=FALSE)

z <- gk2$kin.matrix-kk2
sum(abs(z))

## ----fb---------------------------------------------------------------------------------------------------------------------------------------------
options(width=150)
library(gap)
models <- matrix(c(
         4.0, 0.01,
         4.0, 0.10,
         4.0, 0.50, 
         4.0, 0.80,
         2.0, 0.01,
         2.0, 0.10,
         2.0, 0.50,
         2.0, 0.80,
         1.5, 0.01,    
         1.5, 0.10,
         1.5, 0.50,
         1.5, 0.80), ncol=2, byrow=TRUE)
outfile <- "fbsize.txt"
cat("gamma","p","Y","N_asp","P_A","H1","N_tdt","H2","N_asp/tdt",
    "L_o","L_s\n",file=outfile,sep="\t")
for(i in 1:12) {
    g <- models[i,1]
    p <- models[i,2]
    z <- fbsize(g,p)
    cat(z$gamma,z$p,z$y,z$n1,z$pA,z$h1,z$n2,z$h2,z$n3,
        z$lambdao,z$lambdas,file=outfile,append=TRUE,sep="\t")
    cat("\n",file=outfile,append=TRUE)
}
table1 <- read.table(outfile,header=TRUE,sep="\t")
nc <- c(4,7,9)
table1[,nc] <- ceiling(table1[,nc])
dc <- c(3,5,6,8,10,11)
table1[,dc] <- round(table1[,dc],2)
unlink(outfile)
knitr::kable(table1,caption="Power/Sample size of family-based designs")

## ----alz--------------------------------------------------------------------------------------------------------------------------------------------
g <- 4.5
p <- 0.15
alz <- data.frame(fbsize(g,p))
knitr::kable(alz,caption="Power/Sample size of study on Alzheimer's disease")

## ----pb---------------------------------------------------------------------------------------------------------------------------------------------
library(gap)
kp <- c(0.01,0.05,0.10,0.2)
models <- matrix(c(
          4.0, 0.01,
          4.0, 0.10,
          4.0, 0.50, 
          4.0, 0.80,
          2.0, 0.01,
          2.0, 0.10,
          2.0, 0.50,
          2.0, 0.80,
          1.5, 0.01,    
          1.5, 0.10,
          1.5, 0.50,
          1.5, 0.80), ncol=2, byrow=TRUE)
outfile <- "pbsize.txt"
cat("gamma","p","p1","p5","p10","p20\n",sep="\t",file=outfile)
for(i in 1:dim(models)[1])
{
   g <- models[i,1]
   p <- models[i,2]
   n <- vector()
   for(k in kp) n <- c(n,ceiling(pbsize(k,g,p)))
   cat(models[i,1:2],n,sep="\t",file=outfile,append=TRUE)
   cat("\n",file=outfile,append=TRUE)
} 
table5 <- read.table(outfile,header=TRUE,sep="\t")
knitr::kable(table5,caption="Sample size of population-based design")

## ---------------------------------------------------------------------------------------------------------------------------------------------------
library(gap)
# ARIC study
outfile <- "aric.txt"
n <- 15792
pD <- 0.03
p1 <- 0.25
alpha <- 0.05
theta <- c(1.35,1.40,1.45)
beta <- 0.2
s_nb <- c(1463,722,468)
cat("n","pD","p1","hr","q","power","ssize\n",file=outfile,sep="\t")
for(i in 1:3)
{
  q <- s_nb[i]/n
  power <- ccsize(n,q,pD,p1,log(theta[i]),alpha,beta,TRUE)
  ssize <- ccsize(n,q,pD,p1,log(theta[i]),alpha,beta)
  cat(n,"\t",pD,"\t",p1,"\t",theta[i],"\t",q,"\t",
      signif(power,3),"\t",ssize,"\n",
      file=outfile,append=TRUE)
}
read.table(outfile,header=TRUE,sep="\t")
unlink(outfile)
# EPIC study
outfile <- "epic.txt"
n <- 25000
alpha <- 0.00000005
beta <- 0.2
s_pD <- c(0.3,0.2,0.1,0.05)
s_p1 <- seq(0.1,0.5,by=0.1)
s_hr <- seq(1.1,1.4,by=0.1)
cat("n","pD","p1","hr","alpha","ssize\n",file=outfile,sep="\t")
# direct calculation
for(pD in s_pD)
{
   for(p1 in s_p1)
   {
      for(hr in s_hr)
      {
         ssize <- ccsize(n,q,pD,p1,log(hr),alpha,beta)
         if (ssize>0) cat(n,"\t",pD,"\t",p1,"\t",hr,"\t",alpha,"\t",
                          ssize,"\n",
                          file=outfile,append=TRUE)
      }
   }
}
knitr::kable(read.table(outfile,header=TRUE,sep="\t"),caption="Sample size of case-cohort designs")
unlink(outfile)

## ----qq, fig.cap="A Q-Q plot", fig.height=7, fig.width=7--------------------------------------------------------------------------------------------
library(gap)
u_obs <- runif(1000)
r <- qqunif(u_obs,pch=21,bg="blue",bty="n")
u_exp <- r$y
hits <- u_exp >= 2.30103
points(r$x[hits],u_exp[hits],pch=21,bg="green")
legend("topleft",sprintf("GC.lambda=%.4f",gc.lambda(u_obs)))

## ----chicken, fig.cap="A genome-wide association study on chickens", fig.height=7, fig.width=7, results="hide"--------------------------------------
ord <- with(w4,order(chr,pos))
w4 <- w4[ord,]
oldpar <- par()
par(cex=0.6)
colors <- c(rep(c("blue","red"),15),"red")
mhtplot(w4,control=mht.control(colors=colors,gap=1000),pch=19,srt=0)
axis(2,cex.axis=2)
suggestiveline <- -log10(3.60036E-05)
genomewideline <- -log10(1.8E-06)
abline(h=suggestiveline, col="blue")
abline(h=genomewideline, col="green")
abline(h=0)

## ----mhtplot, fig.cap="A Manhattan plot with gene annotation", fig.height=7, fig.width=7, messages=FALSE, results="hide"----------------------------
data <- with(mhtdata,cbind(chr,pos,p))
glist <- c("IRS1","SPRY2","FTO","GRIK3","SNED1","HTR1A","MARCH3","WISP3",
           "PPP1R3B","RP1L1","FDFT1","SLC39A14","GFRA1","MC4R")
hdata <- subset(mhtdata,gene%in%glist)[c("chr","pos","p","gene")]
color <- rep(c("lightgray","gray"),11)
glen <- length(glist)
hcolor <- rep("red",glen)  
par(las=2, xpd=TRUE, cex.axis=1.8, cex=0.4)
ops <- mht.control(colors=color,yline=1.5,xline=3)
hops <- hmht.control(data=hdata,colors=hcolor)
mhtplot(data,ops,hops,pch=19)
axis(2,pos=2,at=1:16,cex.axis=0.5)

## ----circos, fig.cap="A circos Manhattan plot", fig.height=7, fig.width=8---------------------------------------------------------------------------
circos.mhtplot(mhtdata, glist)

## ----circos2, fig.cap="Another circos Manhattan plot", fig.height=7, fig.width=8--------------------------------------------------------------------
require(gap.datasets)
library(dplyr)
testdat <- mhtdata[c("chr","pos","p","gene","start","end")] %>%
           rename(log10p=p) %>%
           mutate(chr=paste0("chr",chr),log10p=-log10(log10p))
dat <- mutate(testdat,start=pos,end=pos) %>%
       select(chr,start,end,log10p)
labs <- subset(testdat,gene %in% glist) %>%
        group_by(gene,chr,start,end) %>%
        summarize() %>%
        mutate(cols="blue") %>%
        select(chr,start,end,gene,cols)
circos.mhtplot2(dat,labs,ticks=0:2*10)

## ----miami, fig.cap="Miami plots", fig.height=7, fig.width=7----------------------------------------------------------------------------------------
mhtdata <- within(mhtdata,{pr=p})
miamiplot(mhtdata,chr="chr",bp="pos",p="p",pr="pr",snp="rsn")
# An alternative implementation
gwas <- select(mhtdata,chr,pos,p) %>%
        mutate(z=qnorm(p/2))
chrmaxpos <- miamiplot2(gwas,gwas,name1="Batch 2",name2="Batch 1",z1="z",z2="z")
labelManhattan(chr=c(2,16),pos=c(226814165,52373776),name=c("AnonymousGene","FTO"),gwas,gwasZLab="z",chrmaxpos=chrmaxpos)

## ----il12b, echo=FALSE, fig.align="left", fig.cap="Association of IL-12B", fig.height=6, fig.width=6------------------------------------------------
knitr::include_graphics("IL-12B.png")

## ----asplot, fig.cap="A regional association plot", fig.height=7, fig.width=7-----------------------------------------------------------------------
asplot(CDKNlocus, CDKNmap, CDKNgenes, best.pval=5.4e-8, sf=c(3,6))
title("CDKN2A/CDKN2B Region")

## ----cnv, fig.cap="A CNV plot", fig.height=7, fig.width=7-------------------------------------------------------------------------------------------
cnvplot(gap.datasets::cnv)

## ----esplot, fig.cap="rs12075 and traits", fig.height=7, fig.width=7--------------------------------------------------------------------------------
library(gap)
rs12075 <- data.frame(id=c("CCL2","CCL7","CCL8","CCL11","CCL13","CXCL6","Monocytes"),
                      b=c(0.1694,-0.0899,-0.0973,0.0749,0.189,0.0816,0.0338387),
                      se=c(0.0113,0.013,0.0116,0.0114,0.0114,0.0115,0.00713386))
ESplot(rs12075)

## ----forest, fig.cap="Forest plots", fig.height=6, fig.width=9, results="hide"----------------------------------------------------------------------
data(OPG,package="gap.datasets")
METAL_forestplot(OPGtbl[5:6,],OPGall,OPGrsid,width=6.75,height=5,digits.TE=2,digits.se=2)

## ----normal, echo=FALSE, fig.cap="Normal(0,1) distribution", fig.height=5, fig.width=9--------------------------------------------------------------
library(lattice)
z0 <- 1.96
z <- 5
a <- seq(-z, z, length = 10000)
b <- dnorm(a, 0, 1)
xyplot(b ~ a,
       type = "l",
       panel = function(x,y, ...)
               {
                   panel.xyplot(x,y, ...)
                   panel.abline(v = 0, lty = 2)
                   xx <- c(-z, x[x>=-z & x<=-z0], -z0)
                   yy <- c(0, y[x>=-z & x<=-z0], 0)
                   panel.polygon(xx,yy, ..., col='red')
                   xx <- c(z0, x[x>=z0 & x<=z], z)
                   yy <- c(0, y[x>=z0 & x<=z], 0)
                   panel.polygon(xx,yy, ..., col='red')
               },
       xlab="z",
       ylab=expression(1/sqrt(2*pi) * exp(-z^2/2))
)

## ---------------------------------------------------------------------------------------------------------------------------------------------------
require(gap)
v <- data.frame()
for (z in c(5,10,30,40,50,100,500,1000,2000,3000,5000))
{
  vi <- c(z,pvalue(z),logp(z),log10p(z))
  v <- rbind(v,vi)
}
names(v) <- c("z","P","log(P)","log10(P)")
knitr::kable(v,caption="z,P,log(P) and log10(P)")

## ----pve, echo=FALSE, fig.cap="Comparisons of Var($R^2$) estimates", fig.height=6, fig.width=14, messages=FALSE, warning=FALSE----------------------
oldpar <- par()
par(mfrow=c(1,2))
N <- 2:500
R2 <- 2*(N-2)/(N-1)^2/(N+1)
R2LR <- 2/(N+1)^2
R2t <- 2/(N-1)^2
plot(N,R2,cex=0.6,xaxt="n",xlab="Sample size",ylab=expression(Var(R^2)),col="black",pch=20)
points(N,R2LR,cex=0.6,pch=15,col="red")
points(N,R2t,cex=0.6,pch=17,col="blue")
axis(1,at=c(2,(1:5)*100))
legend(400,0.03,c("Asymptotic","LR","t-statistic"),col=c("black","red","blue"),pch=c(20,15,17))
sLR <- (N-2)*(N+1)/(N-1)^2
st <- (N-2)/(N+1)
plot(N,sLR,cex=0.6,xaxt="n",xlab="Sample size",ylab="Asymptotic/approximation estimator ratio",col="red",pch=20)
points(N,st,cex=0.6,pch=15,col="blue")
abline(h=1,col="black")
axis(1,at=c(2,(1:5)*100))
legend(400,0.2,c("Asymptotic","LR","t-statistic"),col=c("black","red","blue"),pch=c(20,15,17))
par(oldpar)

## ---------------------------------------------------------------------------------------------------------------------------------------------------
get_b_se(0.6396966,23991,4.7245)

## ---------------------------------------------------------------------------------------------------------------------------------------------------
get_sdy(0.6396966,23991,0.04490488,0.009504684)

## ----mrdat, echo=FALSE------------------------------------------------------------------------------------------------------------------------------
mrdat <- '
 rs188743906  0.6804   0.1104  0.00177 0.01660        NA        NA
   rs2289779 -0.0788   0.0134  0.00104 0.00261 -0.007543 0.0092258
 rs117804300 -0.2281   0.0390 -0.00392 0.00855  0.109372 0.0362219
   rs7033492 -0.0968   0.0147 -0.00585 0.00269  0.022793 0.0119903
  rs10793962  0.2098   0.0212  0.00378 0.00536 -0.014567 0.0138196
    rs635634 -0.2885   0.0153 -0.02040 0.00334  0.077157 0.0117123
    rs176690 -0.0973   0.0142  0.00293 0.00306 -0.000007 0.0107781
 rs147278971 -0.2336   0.0378 -0.01240 0.00792  0.079873 0.0397491
  rs11562629  0.1155   0.0181  0.00960 0.00378 -0.010040 0.0151460
'
mrdat <- setNames(as.data.frame(scan(file=textConnection(mrdat), what=list("",0,0,0,0,0,0))), 
         c("SNP", "b.LIF.R", "SE.LIF.R", "b.FEV1", "SE.FEV1", "b.CAD", "SE.CAD"))
knitr::kable(mrdat,caption="LIF-R and CAD/FEV1")

## ----mr, echo=FALSE, fig.cap="Mendelian randomization", fig.height=7, fig.width=7-------------------------------------------------------------------
res <- mr(mrdat, "LIF.R", c("CAD","FEV1"), other_plots=TRUE)
s <- res$r[-1,]
colnames(s) <- res$r[1,]
r <- matrix(as.numeric(s[,-1]),nrow(s),dimnames=list(res$r[-1,1],res$r[1,-1]))
p <- sapply(c("IVW","EGGER","WM","PWM"), function(x)
     format(2*pnorm(-abs(r[,paste0("b",x)]/r[,paste0("seb",x)])),digits=3,scientific=TRUE))
rp <- t(data.frame(round(r,3),p))
knitr::kable(rp, align="r", caption="LIFR variant rs635634 and CAD/FEV1")

## ----mr2, fig.cap="Combined forest plots for LIF.R, FEV1 and CAD", fig.height=7, fig.width=7--------------------------------------------------------
mr_names <- names(mrdat)
LIF.R <- cbind(mrdat[grepl("SNP|LIF.R",mr_names)],trait="LIF.R"); names(LIF.R) <- c("SNP","b","se","trait")
FEV1 <- cbind(mrdat[grepl("SNP|FEV1",mr_names)],trait="FEV1"); names(FEV1) <- c("SNP","b","se","trait")
CAD <- cbind(mrdat[grepl("SNP|CAD",mr_names)],trait="CAD"); names(CAD) <- c("SNP","b","se","trait")
mrdat2 <- within(rbind(LIF.R,FEV1,CAD),{y=b})
library(ggplot2)
p <- ggplot2::ggplot(mrdat2,aes(y = SNP, x = y))+
     ggplot2::theme_bw()+
     ggplot2::geom_point()+
     ggplot2::facet_wrap(~ trait, ncol=3, scales="free_x")+
     ggplot2::geom_segment(aes(x = b-1.96*se, xend = b+1.96*se, yend = SNP))+
     ggplot2::geom_vline(lty=2, ggplot2::aes(xintercept=0), colour = 'red')+
     ggplot2::xlab("Effect size")+
     ggplot2::ylab("")
p

## ----tnfb, echo=FALSE-------------------------------------------------------------------------------------------------------------------------------
tnfb <- '
              "multiple sclerosis"  0.69058600 0.059270400
    "systemic lupus erythematosus"  0.76687500 0.079000500
          "sclerosing cholangitis"  0.62671500 0.075954700
   "juvenile idiopathic arthritis" -1.17577000 0.160293000
                       "psoriasis"  0.00582586 0.000800016
            "rheumatoid arthritis" -0.00378072 0.000625160
      "inflammatory bowel disease" -0.14334200 0.025272500
          "ankylosing spondylitis" -0.00316852 0.000626225
                  "hypothyroidism" -0.00432054 0.000987324
               "allergic rhinitis"  0.00393075 0.000926002
          "IgA glomerulonephritis" -0.32696600 0.105262000
                   "atopic eczema" -0.00204018 0.000678061
'

tnfb <- as.data.frame(scan(file=textConnection(tnfb),what=list("",0,0))) %>%
        setNames(c("outcome","Effect","StdErr")) %>%
        mutate(outcome=gsub("\\b(^[a-z])","\\U\\1",outcome,perl=TRUE))

## ---- fig.cap="Forest plots for MR results on TNFB", fig.align="left", fig.height=6, fig.width=9, results="hide"------------------------------------
mr_forestplot(tnfb, colgap.forest.left="0.05cm", fontsize=14, leftlabs=c("Outcome","b","SE"),
              common=FALSE, random=FALSE, print.I2=FALSE, print.pval.Q=FALSE, print.tau2=FALSE,
              spacing=1.6,digits.TE=2,digits.se=2)

## ---- fig.cap="Forest plots for MR results on TNFB (no summary statistics)", fig.align="left", fig.height=6, fig.width=9, results="hide"------------
mr_forestplot(tnfb, colgap.forest.left="0.05cm", fontsize=14,
              leftcols="studlab", leftlabs="Outcome", plotwidth="3inch", sm="OR", rightlabs="ci",
              sortvar=tnfb[["Effect"]],
              common=FALSE, random=FALSE, print.I2=FALSE, print.pval.Q=FALSE, print.tau2=FALSE,
              backtransf=TRUE, spacing=1.6)

## ---- fig.cap="Forest plots for MR results on TNFB (with P values)", fig.align="left", fig.height=6, fig.width=9, results="hide"--------------------
mr_forestplot(tnfb,colgap.forest.left="0.05cm", fontsize=14,
              leftcols=c("studlab"), leftlabs=c("Outcome"),
              plotwidth="3inch", sm="OR", sortvar=tnfb[["Effect"]],
              rightcols=c("effect","ci","pval"), rightlabs=c("OR","95%CI","P"),
              digits=3, digits.pval=2, scientific.pval=TRUE,
              common=FALSE, random=FALSE, print.I2=FALSE, print.pval.Q=FALSE, print.tau2=FALSE,
              addrow=TRUE, backtransf=TRUE, spacing=1.6)

## ---------------------------------------------------------------------------------------------------------------------------------------------------
require(gap)
s <- chr_pos_a1_a2(1,c(123,321),letters[1:2],letters[2:1])
s
inv_chr_pos_a1_a2(s)
inv_chr_pos_a1_a2("chr1:123-A_B",seps=c(":","-","_"))

## ---------------------------------------------------------------------------------------------------------------------------------------------------
gc.lambda <- function(p) {
  p <- p[!is.na(p)]
  n <- length(p)

  obs <- qchisq(p,1,lower.tail=FALSE)
  exp <- qchisq(1:n/n,1,lower.tail=FALSE)

  lambda <- median(obs)/median(exp)
  return(lambda)
}

# A simplified version is as follows,
# obs <- median(chisq)
# exp <- qchisq(0.5, 1) # 0.4549364
# lambda <- obs/exp
# see also estlambda from GenABEL and qq.chisq from snpStats

# A related function

lambda1000 <- function(lambda, ncases, ncontrols)
  1 + (lambda - 1) * (1 / ncases + 1 / ncontrols)/( 1 / 1000 + 1 / 1000)

## ---------------------------------------------------------------------------------------------------------------------------------------------------
invnormal <- function(x) qnorm((rank(x,na.last="keep")-0.5)/sum(!is.na(x)))

## ---------------------------------------------------------------------------------------------------------------------------------------------------
set.seed(12345)
Ni <- rpois(50, lambda = 4); table(factor(Ni, 0:max(Ni)))
y <- invnormal(Ni)
sd(y)
mean(y)
Ni <- 1:50
y <- invnormal(Ni)
mean(y)
sd(y)

## ---------------------------------------------------------------------------------------------------------------------------------------------------
alleles <- c("a","c","G","t")
revStrand(alleles)

## ---- echo=FALSE------------------------------------------------------------------------------------------------------------------------------------
library(gap)
search()
lsf.str("package:gap")
data(package="gap")$results

