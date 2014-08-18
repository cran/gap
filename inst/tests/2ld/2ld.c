#include <stdio.h>
#include <math.h>
#include <ctype.h>
#include <time.h>

#define name "A naive 2-locus LD calculator"
#define version 1.0
#define author "JH Zhao"
#define st_date "12/05/2000"
#define ct_date "15/02/2002"
#define maxalleles 50
#define maxgenotypes maxalleles*(maxalleles+1)/2
#define minval(a,b) ((a<b)?a:b)

enum {EHOUTPUT,RAWDATA,CONTINGTABLE};
enum {PEARSON,TSCHUPROW,CRAMER};
int filetype;
int alleles1,alleles2,dfobs;
double obs[maxgenotypes][maxgenotypes];
double p[maxalleles],q[maxalleles],h[maxalleles*maxalleles],haplotypes;
double sample_size,z1,z2,x2obs,x2lrt;

int getobs(char *);
int getdat(char*);
int tbyt();
int kbyl();
int kbylem(double *,double *);
int kbylse(short,double*,double*,double*);
int kbylinfo(double *);
double probchi(double,int);
double probnorm(double);

int main(int argc, char *argv[])
{
time_t run_time;
double l0,l1,x2;
int df;

printf("%s, version %.2f %s %s-%s IoP\n",name,version,author,st_date,ct_date);
time(&run_time);
if(run_time!=-1) printf("Running on %s",asctime(localtime(&run_time)));
printf("\nMaximum number of alleles = %d\n\n",maxalleles);
if(argc>1) {
#ifdef USEOLDCODE
  getdat(argv[1]);
#else
  getobs(argv[1]);
  haplotypes=2*sample_size;
  if(filetype==RAWDATA||filetype==CONTINGTABLE) {
    kbylem(&l0,&l1);
    x2=2*(l1-l0);
    df=z1*z2-1;
    df-=(z1-1)+(z2-1);
    printf("\nPearson Chi-squared statistic of genotype table");
    printf(" = %.2f df=%d p=%.4lf\n",x2obs,dfobs,probchi(x2obs,dfobs));
    printf("Likelihood ratio Chi-squared statistic of genotype table");
    printf(" = %.2f df=%d p=%.4lf\n",x2lrt,dfobs,probchi(x2lrt,dfobs));
    printf("Chi-squared statistic due to allelelic association = %.2f",x2);
    printf(" df=%d p=%.4lf\n",df,probchi(x2,df));
  }
#endif
  if(alleles1==2&&alleles2==2) tbyt();
  kbyl();
} else {
  printf("Usage: 2ld <ld file> [>outfile]\n\n");
  printf("ld file contains the following information\n\n");
  printf("line 1: <allles at locus 1> <alleles at locus 2> <sample size>\n");
  printf("line 2- haplotype frequencies\n\n");
  printf("Otherwise uses internal example:\n\n");
  sample_size=481;
  h[0]=0.442356;
  h[1]=0.291532;
  h[2]=0.245794;
  h[3]=0.020319;
  haplotypes=2*sample_size;
  printf("2 2 481\n");
  printf("%lf\n%lf\n%lf\n%lf\n",h[0],h[1],h[2],h[3]);
  printf("\nPress <Return> to continue ...\n");
  getc(stdin);
  tbyt();
}
return 0;
}

/*

Weir, BS (1996) Genetic Data Analysis II. Sinauer, p113

Zapata C, Alvarez G, Carollo C (1997) Approximate variance of the standardized
measure of gametic disequilibrium D'. Am. J. Hum. Genet. 61:771-774

notations for a sample of 2n (haplotypes)

     B1  B2
   +--------+
 A1| h0  h1 | p
 A2| h2  h3 | q
   +--------+
     u   v

 */
int tbyt(void)
{
double p,q,u,v,t;
double a,b,xi;
double x2;
double D,ED,VarD;
double Dmax,EDmax,VarDmax;
double Dprime,VarDprime;

p=h[0]+h[1];
q=h[2]+h[3];
u=h[0]+h[2];
v=h[1]+h[3];

D=h[0]-p*u;
ED=D*(haplotypes-1)/haplotypes;
VarD=(p*q*u*v+D*(q-p)*(v-u)-D*D)/haplotypes;
if(D<0) {
  if(p*u<q*v) {
    Dmax=p*u;
    xi=h[0];
  } else {
    Dmax=q*v;
    xi=h[3];
  }
} else {
  if(p*v<q*u) {
    Dmax=p*v;
    xi=h[1];
  } else {
    Dmax=q*u;
    xi=h[2];
  }
}
EDmax=Dmax*(haplotypes-1)/haplotypes;
Dprime=D/Dmax;
if(Dprime<0) {
   a=v;
   b=u;
} else {
   a=u;
   b=v;
}
VarDmax=Dmax*(p*a+q*b-2*fabs(D))/haplotypes;
if(fabs(Dprime)==1) VarDprime=0;
else {
  t=haplotypes*VarD-fabs(Dprime)*Dmax*(p*a+q*b-2*fabs(D));
  VarDprime=((1-fabs(Dprime))*t+fabs(Dprime)*xi*(1-xi))/haplotypes/Dmax/Dmax;
}
x2=haplotypes*D*D/p/q/u/v;

printf("\nDisequilibria, expectations and variances\n\n");
printf("   D \t%9.6f %9.6f %9.6f\n",D,ED,VarD);
printf("Dmax \t%9.6f %9.6f %9.6f\n",Dmax,EDmax,VarDmax);
printf("  D' \t%9.6f           %9.6f\n\n",Dprime,VarDprime);
printf("Chi-squared statistic for D = %.2lf, df = 1, p = %.4lf\n",x2,probchi(x2,1));

return 0;
}

/*

Klitz W, Stephen JC, Grote M, Carrington M (1995) Discordant patterns of
linkage disequilibrium of the peptide transporter loci within the HLA class II
region. Am. J. Hum. Genet. 57:1436-1444

Long JC, Williams RC, Unbanek M (1995) An E-M  algorithm and testing strategy
for multiple-locus haplotypes. Am. J. Hum. Genet. 56:799-810

Weir BS (1996) Genetic Data Analysis II. Sinaur

Zapata C, Carollo C, Rodriquez S (2001) Sampleing variance and distribution
of the D' measure of overall gametic disequilibrium between multiallelic loci.
Ann. Hum. Genet. 65: 395-406

*/

int Dmaxtable[maxalleles][maxalleles];
void abp(int,int,double*,double*,double*,double*);

int kbyl(void)
{
double Dij,VarDij,Dmax,Xij,Dijp,VarDijp,a,b,Eijtable[maxalleles][maxalleles];
double Dijtable[maxalleles][maxalleles],Dijptable[maxalleles][maxalleles];
double Dp,VarDp;
double VarDijtable[maxalleles][maxalleles];
double VarDijptable[maxalleles][maxalleles];
double ai,aip,bj,bjp,ak,akp,bk,bkp,al,alp,bl,blp;
double AI,AIP,BJ,BJP;
double W,Wn,x2,t,tt[4],sgn;
double seX2,rho,seW,klinfo;
double Pij,Pil,Pkj,Dkj,Dkl,Dil,Eij,Eil,Ekj,Ekl;
int i,j,k,l,ij;

printf("\nAllele frequencies\n\n");
printf("locus/allele");
k=(alleles1>alleles2)?alleles1:alleles2;
for(j=0;j<k;j++) printf(" %6d",j+1);
printf("\n");
printf("%6d      ",1);
z1=0;
for(i=0;i<alleles1;i++) {
  if(p[i]==0) {
    printf("       ");
    continue;
  } else {
    printf(" %6.4f",p[i]);
    z1++;
  }
}
printf("\n");
z2=0;
printf("%6d      ",2);
for(j=0;j<alleles2;j++) {
  if(q[j]==0) {
    printf("       ");
    continue;
  } else {
    printf(" %6.4f",q[j]);
    z2++;
  }
}
printf("\n\n * alleles present in the data are %.f, %.f\n",z1,z2);

/*pi*qj*/
printf("\nHaplotype frequencies at linkage equilibrium\n\n");
printf("locus/allele");
for(j=0;j<alleles2;j++) printf(" %6d",j+1);printf("\n");
for(i=0;i<alleles1;i++) {
  printf("%6d      ",i+1);
  if(p[i]>0)
  for(j=0;j<alleles2;j++) {
     if(q[j]==0) {
       printf("       ");
       continue;
     }
     else printf(" %6.4lf",p[i]*q[j]);
   }
   printf("\n");
}
printf("\n");

/*h*/
printf("\nEstimated haplotype frequencies\n\n");
printf("locus/allele");
for(j=0;j<alleles2;j++) printf(" %6d",j+1);printf("\n");
for(i=0;i<alleles1;i++) {
  printf("%6d      ",i+1);
  if(p[i]>0)
  for(j=0;j<alleles2;j++) {
     if(q[j]==0) {
       printf("       ");
       continue;
     }
     else printf(" %6.4lf",h[i*alleles2+j]);
   }
   printf("\n");
}
printf("\n");

Dp=0;
for(i=0;i<alleles1;i++) {
  for(j=0;j<alleles2;j++) {
     ij=i*alleles2+j;
     Dijtable[i][j]=Dij=h[ij]-p[i]*q[j];
     t=p[i]*(1-p[i])*q[j]*(1-q[j]); /*+Dij*(1-2*p[i])*(1-2*q[j])-Dij*Dij;*/
     VarDijtable[i][j]=VarDij=t/haplotypes;
     tt[0]=p[i]*q[j];
     tt[1]=(1-p[i])*(1-q[j]);
     tt[2]=p[i]*(1-q[j]);
     tt[3]=q[j]*(1-p[i]);
     if(Dij<0) {
       a=1-q[j];
       b=q[j];
       if(tt[0]<tt[1]) k=0;
       else k=1;
     } else {
       a=q[j];
       b=1-q[j];
       if(tt[2]<tt[3]) k=2;
       else k=3;
     }
     Dmax=tt[k];
     Dmaxtable[i][j]=k;
     t=h[ij];
     switch(k) {
       case 0:Xij=t;break;
       case 1:Xij=1-p[i]-q[j]+t;break;
       case 2:Xij=p[i]-t;break;
       case 3:Xij=q[j]-t;break;
     }
     abp(i,j,&ai,&bj,&aip,&bjp);
     if(Dmax!=0)
     {
       Dijp=Dij/Dmax;
       if(h[ij]==0||1-p[i]-q[j]+h[ij]==0) Dijp=-1;
       if(p[i]-h[ij]==0||q[j]-h[ij]==0) Dijp=1;
     }
     if(fabs(Dijp)==1) VarDijp=0;
     else if(Dmax!=0) {
       t=haplotypes*VarDij-fabs(Dijp)*Dmax*(a*p[i]+b*(1-p[i])-2*fabs(Dij));
       VarDijp=((1-fabs(Dijp))*t+fabs(Dijp)*Xij*(1-Xij))/haplotypes/Dmax/Dmax;
     }
     Dijptable[i][j]=Dijp;
     VarDijptable[i][j]=VarDijp;
     Dp+=p[i]*q[j]*fabs(Dijp);
     Eijtable[i][j]=(p[i]*aip*bj+q[j]*ai*bjp)*Dij+ai*bj*(Dij-p[i]*q[j]);
   }
}

/*D*/
printf("\nTable of D\n\n");
printf("locus/allele");
for(j=0;j<alleles2;j++) printf(" %6d",j+1);printf("\n");
for(i=0;i<alleles1;i++) {
  printf("%6d      ",i+1);
  if(p[i]>0)
  for(j=0;j<alleles2;j++) {
     if(q[j]==0) {
       printf("       ");
       continue;
     }
     else printf(" %6.3lf",Dijtable[i][j]);
   }
   printf("\n");
}
printf("\n");

/*Dmax*/
printf("\nTable of Dmax\n\n");
printf("locus/allele");
for(j=0;j<alleles2;j++) printf(" %6d",j+1);printf("\n");
for(i=0;i<alleles1;i++) {
  printf("%6d      ",i+1);
  if(p[i]>0)
  for(j=0;j<alleles2;j++) {
     if(q[j]==0) {
       printf("       ");
       continue;
     } else {
       tt[0]=p[i]*q[j];
       tt[1]=(1-p[i])*(1-q[j]);
       tt[2]=p[i]*(1-q[j]);
       tt[3]=q[j]*(1-p[i]);
       printf(" %6.3lf",tt[Dmaxtable[i][j]]);
     }
   }
   printf("\n");
}
printf("\n");

/*D'*/
printf("\nTable of Dij'\n\n");
printf("locus/allele");
for(j=0;j<alleles2;j++) printf(" %6d",j+1);printf("\n");
for(i=0;i<alleles1;i++) {
  printf("%6d      ",i+1);
  if(p[i]>0)
  for(j=0;j<alleles2;j++) {
     if(q[j]==0) {
       printf("       ");
       continue;
     }
     else printf(" %6.3lf",Dijptable[i][j]);
   }
   printf("\n");
}
printf("\n");

/*SE(Dij')*/
printf("\nTable of SE(Dij')\n\n");
printf("locus/allele");
for(j=0;j<alleles2;j++) printf(" %6d",j+1);printf("\n");
for(i=0;i<alleles1;i++) {
  printf("%6d      ",i+1);
  if(p[i]>0)
  for(j=0;j<alleles2;j++) {
     if(q[j]==0) {
       printf("       ");
       continue;
     }
     else printf(" %6.3lf",sqrt(VarDijptable[i][j]));
   }
   printf("\n");
}
printf("\n");

/*Var(D')*/
VarDp=0;
for(i=0;i<alleles1;i++)
for(j=0;j<alleles2;j++)
for(k=0;k<alleles1;k++)
for(l=0;l<alleles2;l++) {
  Pij=h[i*alleles2+j];
  Pil=h[i*alleles2+l];
  Pkj=h[k*alleles2+j];
  Dij=Dijtable[i][j];
  Dil=Dijtable[i][l];
  Dkj=Dijtable[k][j];
  Dkl=Dijtable[k][l];
  Eij=Eijtable[i][j];
  Eil=Eijtable[i][l];
  Ekj=Eijtable[k][j];
  Ekl=Eijtable[k][l];
  abp(i,j,&ai,&bj,&aip,&bjp);
  if(i==k&&j==l) {
    t=Pij*pow((aip*bj+ai*bjp)*Dij+ai*bj*(1-p[i]-q[j]),2)
     +(p[i]-Pij)*pow(aip*bj*Dij-ai*bj*q[j],2)
     +(q[j]-Pij)*pow(ai*bjp*Dij-ai*bj*p[i],2)-Eij*Eij;
     sgn=1;
  }
  else if(i==k&&j!=l) {
    abp(i,l,&AI,&bl,&AIP,&blp);
    t=Pij*((aip*bj+ai*bjp)*Dij+ai*bj*(1-p[i]-q[j]))*(AIP*bl*Dil-AI*bl*q[l])
     +Pil*(aip*bj*Dij-ai*bj*q[j])*((AIP*bl+AI*blp)*Dil+AI*bl*(1-p[i]-q[l]))
     +(p[i]-Pij-Pil)*(aip*bj*Dij-ai*bj*q[j])*(AIP*bl*Dil-AI*bl*q[l])-Eij*Eil;
     if((Dij>0&&Dil>0)||(Dij<0&&Dil<0)) sgn=1;
     else sgn=-1;
  }
  else if(i!=k&&j==l) {
    abp(k,j,&ak,&BJ,&akp,&BJP);
    t=Pij*((aip*bj+ai*bjp)*Dij+ai*bj*(1-p[i]-q[j]))*(ak*BJP*Dkj-ak*BJ*p[k])
     +Pkj*(ai*bjp*Dij-ai*bj*p[i])*((akp*BJ+ak*BJP)*Dkj+ak*BJ*(1-p[k]-q[j]))
     +(q[j]-Pij-Pkj)*(ai*bjp*Dij-ai*bj*p[i])*(ak*BJP*Dkj-ak*BJ*p[k])-Eij*Ekj;
     if((Dij>0&&Dkj>0)||(Dij<0&&Dkj<0)) sgn=1;
     else sgn=-1;
  }
  else if(i!=k&&j!=l) {
    abp(k,l,&ak,&bl,&akp,&blp);
    t=Pkj*(akp*bl*Dkl-ak*bl*q[l])*(ai*bjp*Dij-ai*bj*p[i])
     +Pil*(aip*bj*Dij-ai*bj*q[j])*(ak*blp*Dkl-ak*bl*p[k])-Eij*Ekl;
     if((Dij>0&&Dkl>0)||(Dij<0&&Dkl<0)) sgn=1;
     else sgn=-1;
  }
  VarDp+=sgn*t;
}
VarDp/=haplotypes;

/*individual coefficients tested separately with chi-square*/
/*doubled for the overall hypothesis of all nonzeros*/
printf("\nTable of chi-square [based on D]\n\n");
printf("locus/allele");
for(j=0;j<alleles2;j++) printf(" %6d",j+1);printf("\n");
W=0;
for(i=0;i<alleles1;i++) {
  printf("%6d      ",i+1);
  for(j=0;j<alleles2;j++) {
     if(p[i]==0||q[j]==0) printf("       ");
     else {
       Dij=Dijtable[i][j];
       printf(" %6.2lf",Dij*Dij/VarDijtable[i][j]);
       W+=Dij*Dij/p[i]/q[j];
     }
  } printf("\n");
}
printf("\n");

/*P value*/
printf("\nTable of p value [based on D]\n\n");
printf("locus/allele");
for(j=0;j<alleles2;j++) printf(" %6d",j+1);printf("\n");
for(i=0;i<alleles1;i++) {
  printf("%6d      ",i+1);
  for(j=0;j<alleles2;j++) {
     if(p[i]==0||q[j]==0) printf("       ");
     else {
       Dij=Dijtable[i][j];
       printf(" %6.4lf",probchi(Dij*Dij/VarDijtable[i][j],1));
     }
  } printf("\n");
}
printf("\n");

x2=haplotypes*W;
W=sqrt(W);
t=(z1<z2)?z1:z2;
Wn=W/sqrt(t-1);

kbylse(CRAMER, &seX2, &rho, &seW);

k=(z1-1)*(z2-1);
printf("\nChi-squared statistic = %.2f, df = %d, p = %.4lf\n",x2,k,probchi(x2,k));
printf("(Pearson Chi-squared statistic of haplotype table)\n");
printf("\n\nGlobal disequilibrium statistics and their standard errors:\n");
printf("\nPhi coefficient W = (%.2f / %.0f)^0.5 = %.4f, SD = %.4f\n",x2,haplotypes,W,seX2);
printf("Cramer's V (W/[min(%.lf,%.lf)-1]^0.5) = %.4f, SD = %.4f\n",z1,z2,Wn,seW);
printf("\n");
printf("D' coefficient = %f, SD = %.4f (Var = %.6f)\n",Dp,sqrt(VarDp),VarDp);
kbylinfo(&klinfo);
printf("\nKullback-Leibler information = %lf\n",klinfo);

return 0;
}

/*
 alpha, beta and their derivatives as in Zapata et al. (2001)
 */
void abp(int i, int j, double *a, double *b, double *ap, double *bp)
{
double pi;
double qj;

pi=p[i];
qj=q[j];
switch (Dmaxtable[i][j]) {
case 0:
     *a=*b=1;
     *ap=*bp=0;
     break;
case 1:
     *a=pi/(1-pi);*ap=1/(1-pi)/(1-pi);
     *b=qj/(1-qj);*bp=1/(1-qj)/(1-qj);
     break;
case 2:
     *a=1;*ap=0;
     *b=qj/(1-qj);*bp=1/(1-qj)/(1-qj);
     break;
case 3:
     *a=pi/(1-pi);*ap=1/(1-pi)/(1-pi);
     *b=1;*bp=0;
     break;
default: printf("definitely something went wrong\n");
}
}

/*
 Cramer H (1946) Mathematical Methods of Statistics. Princeton Univ. Press
 Bishop YMM, Fienberg SE, Holland PW (1975) Discrete Multivariate Analysis
       -- Theory and Practice, The MIT press
 */

int kbylse(short opt, double *seX2, double *rho, double *se)
{
double phi2=0,p2[maxalleles],q2[maxalleles];
double s0,s1,s2,s3,s4,sk,sl,t;
int i,j,k,l;

for(i=0;i<alleles1;i++) p2[i]=0;
for(j=0;j<alleles2;j++) q2[j]=0;

s1=s2=s3=s4=0;

phi2=0;
for(i=0;i<alleles1;i++)
  for(j=0;j<alleles2;j++) {
  if(p[i]==0||q[j]==0) s0=0;
  else s0=pow(h[i*alleles2+j],2)/p[i]/q[j];
  p2[i]+=s0;
  q2[j]+=s0;
  phi2+=s0;
  if(p[i]==0||q[j]==0) s0=0;
  else s0=pow(h[i*alleles2+j],3)/pow(p[i]*q[j],2);
  s1+=s0;
}
phi2-=1;

for(i=0;i<alleles1;i++) if(p[i]!=0) s2+=p2[i]*p2[i]/p[i];
for(j=0;j<alleles2;j++) if(q[j]!=0) s3+=q2[j]*q2[j]/q[j];

#ifdef FOOL
fprintf(stderr,"\ts2=%lf\ts3=%lf",s2,s3);
s2=0;
for(i=0;i<alleles1;i++) {
 sk=0;
 for(j=0;j<alleles2;j++) sk+=pow(h[i*alleles2+j],2)/p[i]/q[j];
 s2+=sk*sk/p[i];
}
s3=0;
for(j=0;j<alleles2;j++) {
 sl=0;
 for(i=0;i<alleles1;i++) sl+=pow(h[i*alleles2+j],2)/p[i]/q[j];
 s3+=sl*sl/q[j];
}
fprintf(stderr,"\ts2=%lf\ts3=%lf",s2,s3);
#endif

for(i=0;i<alleles1;i++)
  for(j=0;j<alleles2;j++) {

#ifdef FOOL
    sk=0;
    for(k=0;k<alleles1;k++) sk+=pow(h[k*alleles2+j],2)/p[k]/q[j];
    sl=0;
    for(l=0;l<alleles2;l++) sl+=pow(h[i*alleles2+l],2)/p[i]/q[l];
    s4+=p[i][j]/p[i]/q[j]*sk*sl;
#endif
    if(p[i]!=0&&q[j]!=0) s4+=h[i*alleles2+j]/p[i]/q[j]*p2[i]*q2[j];
}
t=sqrt((4*s1-3*s2-3*s3+2*s4)/haplotypes);
*se=t;
*seX2=t;
switch(opt)
{
case PEARSON:
     *rho=sqrt(phi2/(phi2+1));
     *se*=0.5/sqrt(phi2)/pow(1+phi2,1.5);
     break;
case TSCHUPROW:
     s0=(z1-1)*(z2-1);
     *rho=sqrt(phi2/sqrt(s0));
     *se*=0.5/s0/(*rho);
     break;
case CRAMER:
     s0=(z1<z2)?(z1-1):(z2-1);
     *rho=sqrt(phi2/s0);
     *se*=0.5/sqrt(s0)/(*rho);
     break;
default:
     fprintf(stderr,"Invalid option not supported");
}
return 0;
}

/*
  Kullback-Leibler information measure
 */

double cellinfo(int,int,int,int,int);

int kbylinfo(double *t)
{
int i,j,k,l,ik,jl,il,jk;
int nhet;
double po,pe;

*t=0;
for(i=0;i<alleles1;i++) {
  for(j=0;j<=i;j++) {
    for(k=0;k<alleles2;k++) {
      for(l=0;l<=k;l++) {
#ifdef USEOLDCODE
        if(i==j&&k==l) nhet=0;
        if((i==j&&k!=l)||(i!=j&&k==l)) nhet=1;
        if(i!=j&&k!=l) nhet=2;
        *t+=cellinfo(nhet,i,j,k,l);
#else
        ik=i*alleles2+k;
        jl=j*alleles2+l;
        if((i!=j)&&(k!=l)) {
          il = i * alleles2 + l;
          jk = j * alleles2 + k;
          po = 2.0 * (h[ik] * h[jl] + h[il] * h[jk]);
          pe=2.0*(p[i]*q[k]*p[j]*q[l]+p[i]*q[l]*p[j]*q[k]);
        } else {
          if((i==j)&&(k==l)) {
            po = h[ik]*h[ik];
            pe = p[i]*q[k]*p[i]*q[k];
          }
          else {
            po = 2.0*h[ik]*h[jl];
            pe = 2.0*p[i]*q[k]*p[j]*q[l];
          }
        }
        if(po!=0&&pe!=0) *t+=po*log(po/pe);
#endif
      }
    }
  }
}

return 0;
}

/*
  information from [ij][kl] based on genotype frequencies
 */

double cellinfo(int nhet,int i,int j,int k,int l)
{
double po,pe,ci;
int ik,il,jk,jl;

ik=i*alleles2+k;
il=i*alleles2+l;
jk=j*alleles2+k;
jl=j*alleles2+l;

ci=0;
switch(nhet) {
case 0:
      po=h[ik]*h[ik];
      pe=p[i]*q[k]*p[i]*q[k];
      break;
case 1:
      po=2*h[ik]*h[jl];
      pe=2*p[i]*q[k]*p[j]*q[l];
      break;
case 2:
      po=2*(h[ik]*h[jl]+h[il]*h[jk]);
      pe=2*(p[i]*q[k]*p[j]*q[l]+p[i]*q[l]*p[j]*q[k]);
      break;
}

if(po!=0&&pe!=0) ci=po*log(po/pe);

return ci;
}

/*
 obsolete
*/

int getdat(char *infile)
{
FILE *fin;
char line[200];
double t;
int i,j;

fin=fopen(infile,"r");
if(!fin) {
  fprintf(stderr,"I cannot open a file called %s\n",infile);
  exit(1);
} else {
  if(fgets(line,200,fin)&&sscanf(line,"%d %d %lf",&alleles1,&alleles2,&sample_size)<2) {
    fprintf(stderr,"I need read two numbers representing # alleles at two loci and # subjects\n");
    exit(1);
  }
  i=1;
  while (fgets(line,200,fin)&&sscanf(line,"%lf",&h[i-1])&&(i<alleles1*alleles2)) i++;
  if(i<alleles1*alleles2) fprintf(stderr,"I need %d haplotype frequencies\n",alleles1*alleles2);
  printf("%d haplotype frequencies based on a sample of %.0f diploid individuals\n",i,sample_size);
  for(i=0;i<alleles1;i++) {
    t=0;
    for(j=0;j<alleles2;j++) {
      t+=h[i*alleles2+j];
    }
    p[i]=t;
  }
  for(j=0;j<alleles2;j++) {
    t=0;
    for(i=0;i<alleles1;i++) {
      t+=h[i*alleles2+j];
    }
    q[j]=t;
  }
}
fclose(fin);
return 0;
}

/*
 obtain haplotype frequencies and log-likelihood for two loci
 changed from ASSOCIAT.PAS but trade efficiency for clarification

 19/02/2001 fix bug on k1,k2
 */

int kbylem(double *l0,double *l1)
{
int iter,i,j,k,l,k1,k2,ik,jl,il,jk;
double pobs,e1,e2,r1,r2;
double hc[maxalleles*maxalleles];

k=0;
for(i=0;i<alleles1;i++) for(j=0;j<alleles2;j++) {
  h[i*alleles2+j]=p[i]*q[j];
  hc[k++]=0;
}
iter=0;
do {
  *l1=0;
  k1=0;
  for(i=0;i<alleles1;i++) {
    for(j=0;j<=i;j++) {
      k2 = 0;
      for(k=0;k<alleles2;k++) {
        for(l=0;l<=k;l++) {
          ik = i * alleles2 + k;
          jl = j * alleles2 + l;
          if((i!=j)&&(k!=l)) {
            il = i * alleles2 + l;
            jk = j * alleles2 + k;
            r1 = 2.0 * h[ik] * h[jl];
            r2 = 2.0 * h[il] * h[jk];
            pobs=r1+r2;
            if(obs[k1][k2]>0) {
               e1=r1/pobs;
               e2=r2/pobs;
               hc[ik] += e1*obs[k1][k2];
               hc[il] += e2*obs[k1][k2];
               hc[jk] += e2*obs[k1][k2];
               hc[jl] += e1*obs[k1][k2];
            }
          } else {
            if((i==j)&&(k==l)) {
              pobs = h[ik]*h[ik];
              hc[ik] += 2*obs[k1][k2];
            } else {
              pobs = 2.0*h[ik]*h[jl];
              hc[ik] += obs[k1][k2];
              hc[jl] += obs[k1][k2];
            }
          }
          if (obs[k1][k2]>0) *l1+= obs[k1][k2] * log(pobs);
          ++k2;
        }
      }
      ++k1;
    }
  }
  for(i=0;i<alleles1*alleles2;i++) {
    h[i]=hc[i]/sample_size/2;
    hc[i]=0;
  }
  if(iter==0) *l0=*l1;
} while(iter++<15);

return 0;
}

/*
  now support both marker genotypes and EH output

  JH Zhao 19/02/2001

 */

int getobs(char *obsfile)
{
FILE *fp;
char line[301],id[20];
int a1,a2,b1,b2;
int i,j,k,l,l1,l2,u1,u2,nmiss;
double t,rt[maxgenotypes],ct[maxgenotypes];

for(i=0;i<maxalleles;i++) p[i]=0;
for(j=0;j<maxalleles;j++) q[j]=0;
for(i=0;i<maxgenotypes;i++) rt[i]=ct[i]=0;
fp=fopen(obsfile,"r");
if(!fp) {
  printf("Sorry, I cannot open file %s\n",obsfile);
  exit(1);
}
if(fgets(line,301,fp)&&sscanf(line,"%s %d %d %d %d",id,&a1,&a2,&b1,&b2)>4)
{
  filetype=RAWDATA;
  rewind(fp);
  alleles1=alleles2=0;
  sample_size=0;
  nmiss=0;
  for(i=0;i<maxgenotypes;i++) for(j=0;j<maxgenotypes;j++) obs[i][j]=0;
  goto rawcodes;
}
rewind(fp);
if(fgets(line,301,fp)&&sscanf(line,"%d%d%lf",&alleles1,&alleles2,&sample_size)>2)
filetype=EHOUTPUT;
else {
  rewind(fp);
  if(fgets(line,301,fp)&&sscanf(line,"%d%d",&alleles1,&alleles2)) {
    filetype=CONTINGTABLE;
    sample_size=0;
    l1=0;
    for(i=0;i<alleles1;i++) {
      for(j=0;j<=i;j++) {
        l2=0;
        for(k=0;k<alleles2;k++) {
          for(l=0;l<=k;l++) {
            fscanf(fp,"%lf",&obs[l1][l2]);
            sample_size+=obs[l1][l2];
            p[i]+=obs[l1][l2];
            p[j]+=obs[l1][l2];
            q[k]+=obs[l1][l2];
            q[l]+=obs[l1][l2];
            rt[l1]+=obs[l1][l2];
            ct[l2]+=obs[l1][l2];
            l2++;
          }
        }
        l1++;
      }
    }
    goto ok;
  }
  printf("Sorry, but what type of file is this ?\n");
  exit(1);
}
if(filetype==EHOUTPUT)
{
  i=1;
  while(fgets(line,301,fp)&&sscanf(line,"%lf",&h[i-1])&&(i<alleles1*alleles2)) i++;
  if(i<alleles1*alleles2) fprintf(stderr,"I need %d haplotype frequencies\n",alleles1*alleles2);
  printf("%d haplotype frequencies based on a sample of %.0f diploid individuals\n",i,sample_size);
  for(i=0;i<alleles1;i++) {
    t=0;
    for(j=0;j<alleles2;j++) {
      t+=h[i*alleles2+j];
    }
    p[i]=t;
  }
  for(j=0;j<alleles2;j++) {
    t=0;
    for(i=0;i<alleles1;i++) {
      t+=h[i*alleles2+j];
    }
    q[j]=t;
  }
  return 0;
}
rawcodes:
while(fgets(line,301,fp)&&sscanf(line,"%s %d %d %d %d",id,&a1,&a2,&b1,&b2)>3) {
  if((a1!=0)&&(a2!=0)&&(b1!=0)&&(b2!=0)) {
    l1=(a1<a2)?a1:a2;
    u1=(a1>=a2)?a1:a2;
    p[a1-1]++;
    p[a2-1]++;
    if(u1>alleles1) alleles1=u1;
    l2=(b1<b2)?b1:b2;
    u2=(b1>=b2)?b1:b2;
    q[b1-1]++;
    q[b2-1]++;
    if(u2>alleles2) alleles2=u2;
    obs[l1+u1*(u1-1)/2-1][l2+u2*(u2-1)/2-1]++;
    sample_size++;
  } else nmiss++;
}
fclose(fp);
printf("%.0f individuals with full genotypic information and %d not.\n",sample_size,nmiss);
l1=0;
for(i=0;i<alleles1;i++) {
  for(j=0;j<=i;j++) {
    l2=0;
    for(k=0;k<alleles2;k++) {
      for(l=0;l<=k;l++) {
        rt[l1]+=obs[l1][l2];
        ct[l2]+=obs[l1][l2];
        l2++;
      }
    }
    l1++;
  }
}
ok:
x2obs=x2lrt=0;
l1=0;
for(i=0;i<alleles1;i++) {
  for(j=0;j<=i;j++) {
    l2=0;
    for(k=0;k<alleles2;k++) {
      for(l=0;l<=k;l++) {
        if(rt[l1]>0&&ct[l2]>0) {
          t=rt[l1]*ct[l2]/sample_size;
          x2obs+=pow(obs[l1][l2]-t,2)/t;
          if(obs[l1][l2]>0)
          x2lrt+=obs[l1][l2]*log(obs[l1][l2]/rt[l1]/ct[l2]*sample_size);
        }
        l2++;
      }
    }
    l1++;
  }
}
x2lrt*=2;
z1=0;
z2=0;
for(i=0;i<alleles1;i++) {
  p[i]/=sample_size*2;
  if(p[i]!=0) z1++;
}
for(j=0;j<alleles2;j++) {
  q[j]/=sample_size*2;
  if(q[j]!=0) z2++;
}
dfobs=(z1*(z1+1)/2-1)*(z2*(z2+1)/2-1);

return 0;
}

/*
 Linkage Utility programs.

 Error probability, probchi = 1-f, where f is the chi-square distribution
 function with ndf degrees of freedom evaluated at x2.

 Function required: probnorm.

 Abramowitz M and Stegun IA (1968) Handbook of mathematical functions. New York

 26.4.4 and 26.4.21
 */

double probchi(double x2, int ndf)
{
double z,p,x,sum,re,ch,chp;
short i,n1,n2;

if(x2<=0.0&&ndf<=0) return 1.0;
/*ndf=1*/
if(ndf==1)
{
  x=sqrt(x2);return(2.0*probnorm(x));
}
/*ndf>1*/
if(x2>169.0){                   /*formula 26.4.14, p.941*/
  ch=2.0/(9.0*ndf);
  x=(exp(log(x2/ndf)/3.0)-1.0+ch)/sqrt(ch);
  return(probnorm(x));
}
/*ndf=2*/
if(ndf==2) return(exp(-0.5*x2));
/*ndf>2*/
n1=(ndf-1)/2;
n2=(ndf-2)/2;
if(n1==n2) /*ndf is even and >2*/
{
  sum=1.0;
  re=0.5*x2;
  ch=1.0;
  for(i=1;i<=n2;i++)
  {
    ch=ch*re/i;
    sum+=ch;
  }
  return(exp(-re)*sum);
}
else       /*ndf is odd and >1*/
{
  ch=sqrt(x2);
  z=0.39894228*exp(-0.5*x2);
  p=probnorm(ch);
  if(ndf==3) return(2.0*(p+z*ch)); /*ndf=3*/
  else     /*ndf odd and >3*/
  {
    chp=ch;
    re=1.0;
    for(i=2;i<=n1;i++)
    {
      re+=2.0;
      chp=chp*x2/re;
      ch+=chp;
    }
    return(2.0*(p+z*ch));
  }
}
}

/*
 upper one sided tail probability of the normal distribution
 for a given normal deviate, x, by 26.2.16 in Abramowitz and Stegun.
*/

double probnorm(double x)
{
double z,t,p,xa;

xa=fabs(x);
if(xa>12.0) p=0.0;
else
{
  z=0.39894228*exp(-0.5*xa*xa);
  t=1.0/(1.0+0.33267*xa);
  p=z*t*(0.4361836+t*(0.937298*t-0.1201676));
}
if(x>=0.0) return p;else return(1.0-p);
}

/*
 2ld - two-locus LD calculator JH Zhao 2000-2002

 History.

 12-05-2000 add s.e. to David Collier's Excel macro (2x2, sphadco@iop)
 11-08-2000 add s.e. to Cramer's V statistic, etc.
 08-11-2000 start experiment on multiallelic S.E. (D')
 08-02-2001 fix ``||'' bug and handle missing alleles by Zhicheng Lin's data
 10-02-2001 simplify implementation
 15-02-2001 careful check (2x2,kxl)
 16-02-2001 add Kullback-Leibler information measure
 19-02-2001 accommodate raw genotype data, fix bug in K-L information
 21-02-2001 accommodate contingency table
 30-03-2001 add contingency table Chi-squares, fix bug for phenotype data
 16-07-2001 fix bug in obtaining maximum allele at the second locus
 28-09-2001 change variance formula for Dij
 10-10-2001 fix table layout when the first locus has missing alleles
 11-10-2001 update reference
 15-02-2002 name/author/version/date & time stamp, modify documentation
 17-05-2002 change label for Chi-squared statistic, move this block to the end
 10-09-2002 fix error in Var(Dij')

 */
