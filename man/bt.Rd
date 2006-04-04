\name{bt}
\alias{bt}
\title{Bradley-Terry model for contingency table}
\usage{bt(x)}
\description{
This function calculates statistics under Bradley-Terry model. 
}
\arguments{
\item{x}{the data table}
}

\value{
The returned value is a list containing:
\item{y}{A column of 1}
\item{count}{the frequency count/weight}
\item{allele}{the design matrix}
\item{bt.glm}{a glm.fit object}
\item{etdt.dat}{a data table that can be used by ETDT}
}

\section{References}{
Bradley RA, Terry ME (1952) Rank analysis of incomplete block designs I. 
the method of paired comparisons. Biometrika 39:324--345

Sham PC, Curtis D (1995) An extended transmission/disequilibrium 
test ({TDT}) for multi-allelic marker loci. Ann. Hum. Genet. 59:323-336
}
\seealso{
\code{\link[gap]{mtdt}}
}

\examples{
\dontrun{
# Copeman JB, Cucca F, Hearne CM, Cornall RJ, Reed PW, 
# Ronningen KS, Undlien DE, Nistico L, Buzzetti R, Tosi R, et al.
# (1995) Linkage disequilibrium mapping of a type 1 
# diabetes susceptibility gene (IDDM7) to chromosome 2q31-q33. 
# Nat Genet 9: 80-5

x <- matrix(c(0,0, 0, 2, 0,0, 0, 0, 0, 0, 0, 0,
              0,0, 1, 3, 0,0, 0, 2, 3, 0, 0, 0,
              2,3,26,35, 7,0, 2,10,11, 3, 4, 1,
              2,3,22,26, 6,2, 4, 4,10, 2, 2, 0,
              0,1, 7,10, 2,0, 0, 2, 2, 1, 1, 0,
              0,0, 1, 4, 0,1, 0, 1, 0, 0, 0, 0,
              0,2, 5, 4, 1,1, 0, 0, 0, 2, 0, 0,
              0,0, 2, 6, 1,0, 2, 0, 2, 0, 0, 0,
              0,3, 6,19, 6,0, 0, 2, 5, 3, 0, 0,
              0,0, 3, 1, 1,0, 0, 0, 1, 0, 0, 0,
              0,0, 0, 2, 0,0, 0, 0, 0, 0, 0, 0,
              0,0, 1, 0, 0,0, 0, 0, 0, 0, 0, 0),nrow=12)

# Bradley-Terry model, only deviance is available in glm
# (SAS gives score and Wald statistics as well)
bt.ex<-bt(x)
anova(bt.ex$bt.glm)
summary(bt.ex$bt.glm)
}
}

\author{Jing Hua Zhao}
\keyword{models}
\note{
\preformatted{
/*Adapted from the SAS macro below for data in the example section*/
\%macro mtdt(data,n);
data _bt_;
  set &data;
  array x {&n} x1-x&n;
  array allele {&n} y1-y&n;
  do i=1 to &n; allele{i}=0; end;
  y=1;
  do i=1 to &n;
     allele{_n_}=1;
     allele{i}=-1;
     count=x{i};
     if _n_ ne i then output;
     allele{i}=0;
  end;
  keep y count y1-y&n;
run;
/*Bradly-Terry model*/
proc logistic data=_bt_;
  freq count;
  model y=y1-y&n / noint;
  output out=out p=p;
  run;
/*Bowker's test of symmetry*/   
data b;
  array x x1-x&n;
  do i=1 to &n;
     set &data;
     do j=1 to &n; w=x[j]; output; end; 
  end;
  drop x1-x&n;
run;
proc freq;
   weight w;   
   table i*j / agree noprint;
   run;
\%mend;
data a;
  input x1-x12;
cards;
0 0  0  2  0 0  0  0  0  0  0  0
0 0  1  3  0 0  0  2  3  0  0  0
2 3 26 35  7 0  2 10 11  3  4  1
2 3 22 26  6 2  4  4 10  2  2  0
0 1  7 10  2 0  0  2  2  1  1  0
0 0  1  4  0 1  0  1  0  0  0  0
0 2  5  4  1 1  0  0  0  2  0  0
0 0  2  6  1 0  2  0  2  0  0  0
0 3  6 19  6 0  0  2  5  3  0  0
0 0  3  1  1 0  0  0  1  0  0  0
0 0  0  2  0 0  0  0  0  0  0  0
0 0  1  0  0 0  0  0  0  0  0  0
;
\%mtdt(a,12);
}}
