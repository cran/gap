#ifndef FASTPM_H
#define FASTPM_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <assert.h>

#define MAX_LOC 30
#define maxalleles    50         /*25*/

/*EM code*/

#define maxloops     100
#define eps          0.1e-10
#define exp2(h)      (pow(2,h))

#define maxloci   MAX_LOC
#define maxalle   maxalleles

typedef enum         BOOL {false,true} BOOL;
BOOL delta=true;
BOOL cc,err,onemark,wd,addum;

char ch;
FILE *finp,*fout;
int sample_size,obscom,totalloci,haplnum,maxall,loop,loop1;
int locihn[maxloci+1],locigeno[maxloci+1],loci[maxloci+1];
long int obshap;
short locip[maxloci+1];
double *indivN,*caseN,*contrlN,*dg11,*dg12,*dg22,p[maxloci][maxalle];
double *oldhaplo,*newhaplo,*inihaplo,*indhaplo;
double oldloglike,loglike,iniloglike,diff,indloglike,oldlog;
double totalall,totalind,ff[3],ss[3],casep[3],controlp[3],s2,pp,qq;

#define useid
#ifdef useid
double *id,*idsave;
#endif

typedef struct
{
  short l[maxloci],u[maxloci];
} phenotype;
phenotype *alist;

struct DBT {
  BOOL dma;
  int time,nloci,fl,ll;
  long int l,u;
};

void eh(int);
void dirichlet(long int);

/*
 * correct_df=true would correct for any missing alleles
 */
BOOL correct_df=true,correct_den=true;

/*PM code*/

enum {WHOLE,BLOCK1,BLOCK2};
#define _swap_(a,b) __swap__(a,b)
void __swap__(int *a,int *b) {int t;t=*a;*a=*b;*b=t;}
#define LL 1500
typedef struct
{
  char id[10];
  BOOL affection;
  int no;
  int locus[MAX_LOC][2];
  int gtype[MAX_LOC];
} individual;
individual p_t, *person, *person_s;
typedef struct
{
  int nloci,nloci_s,selected,selectp,selectn;
  BOOL sel[MAX_LOC],selp[MAX_LOC];
  BOOL selidx[MAX_LOC],selpidx[MAX_LOC],selnpdx[MAX_LOC];
} pm_model;
pm_model orig_model,any_model;

int getloci(char*),getdat(char*);
int a2g(int,int);
int g2a(int,int*,int*,int*);
int pgtype(int),ehout(char *),plabel(void),ind2eht(void);
int muvar(int, double, double*, double*);
double xi(int, double, double, double, double);
double varxi(int, double, double);

int nloci,nloci_s,alleles[MAX_LOC],permute,npermute;
double nall[MAX_LOC],np[MAX_LOC],nnp[MAX_LOC],position(int,int*,int);
BOOL sel[MAX_LOC],selp[MAX_LOC],isgenotype,iogenotype;
BOOL selidx[MAX_LOC],selpidx[MAX_LOC],selnpdx[MAX_LOC];
int selected,selectn,selectp,rows,cols;
float freq,pen0,pen1,pen2;
int cases;
int seed=3000;
int iter;

/*BST code: insert,search,delete,print,remove,confer*/

typedef struct node_type
{
  double genid;           /*unique genotype identifier*/
  int nca;                  /*number of cases*/
  int nco;                  /*number of controls*/
  short l[MAX_LOC],u[MAX_LOC];  /*actual phenotypes*/
  struct node_type *left, *right;
} node;

node *rt=0, *rt1=0, *rt2=0;

node *itree(node*,double,int);
node *stree(node*,double);
node *dtree(node*,double);
void inorder(node*),preorder(node*),postorder(node*);
void rtree(node*),etree(int),ptree(node*,int,FILE*,int);
void ctree(node*,int);
int outehpdat(char*,node*,int),getehpdat(void);

void *xmalloc(long len)
{
void *mem;

mem=(void *)malloc(len);
if(mem) return (void *)mem;
else
{
fprintf(stderr,"Sorry, but I cannot allocate memory\n");
return NULL;
}
}

/*EH+ procedures*/

double likelihood(double *,int,BOOL,BOOL);
void alistout(void);
void haplotype(double *);
void probgeno(double *, struct DBT *, double *);
void probgen(double *, struct DBT *, short *, short *, double *);
void getden(double *,int,BOOL,double *,double *);
int genecount(int,int,long int *,double *,struct DBT *,double *);
int genec(int, int, double *, struct DBT *, short *, short *, double *);
double getdiff(void);

void hilo(short *a, short *b)
{
short temp;
temp = *a;*a = *b;*b = temp;
}

void initP(void)
/*initialize allele frequencies across loci*/
{
int i,j;
for(i=0;i<totalloci;i++) for(j=0;j<maxall;j++) p[i][j]=0.0;
totalind=0.0;
}

void calcP(int nloci)
/*calculate allele frequencies based on allele counts*/
{
int i,j;
for(i=0;i<nloci;i++) for(j=0;j<maxall;j++) p[i][j]/=totalall;
}

void initV(int nloci, int *locii)
/*initialize numbers of haplotypes and genotypes*/
{
int i;
haplnum=1;
for(i=0;i<nloci;i++){
  haplnum*=locii[i];
  locigeno[i]=locii[i]*(locii[i]+1)/2;
}
}

double getdiff()
/*get absolute difference between haplotype frequency estimates*/
{
  int i;
  double temp=0.0;
  for(i=0;i<haplnum;i++) temp+=fabs(oldhaplo[i]-newhaplo[i]);
  return temp;
}

int linenum(int *locia, short *locip, int sp)
/*the genotype/haplotype/allele identifier*/
{
int sum,j;
sum=0;
for(j=sp;j<=totalloci;j++)
  if(j==totalloci) sum+=locip[j-1]; else sum=(sum+locip[j-1]-1)*locia[j];
return(sum);
}

void getgeno(int totalloci)
/*get number of possible genotypes of each locus*/
{
int i;
for(i=0;i<totalloci;i++) locigeno[i]=loci[i]*(loci[i]+1)/2;
}

void getfirst(int *first, int *last)
/*get the first heterozygous locus*/
{
int j;
BOOL ini;
*first=*last=1;
ini=true;
for (j=1;j<=totalloci;j++)
  if (locihn[j-1]==1) {
    if (ini) if(*first<=j) {
        *first=j;
        ini=false;
      }
    if(j>*last) *last=j;
  }
}

double nget(int n,int nloci)
/*gene counting without recursion*/
{
int i,j,j1,j2;
double tempn,all=0.0;
for (i=0;i<n;i++) {
  tempn=indivN[i];
  for (j=0;j<nloci;j++) {
    j1=alist[i].l[j]-1;
    j2=alist[i].u[j]-1;
    p[j][j1]+=tempn;
    p[j][j2]+=tempn;
  }
  all+=tempn;
}
return all;
}

void getH(int nloci)
/*get haplotype frequencies from allele frequencies*/
{
int i;
long int line, times, n;

n=times=1;
for(i=0;i<nloci;i++){
  locip[i]=1;
  n*=loci[i];
}
do {
  line=linenum(loci,locip,1);
  newhaplo[line-1]=1.0;
  for(i=0;i<nloci;i++) newhaplo[line-1]*=p[i][locip[i]-1];
  locip[nloci-1]++;
  for(i=nloci-1;i>0;i--) {
    if(locip[i]>loci[i]){
      locip[i]=1;
      locip[i-1]++;
    }
  }
  times++;
} while(times<=n);
for(i=0;i<n;i++) inihaplo[i]=newhaplo[i];
#ifdef DIRICHLET
if(iter>0) dirichlet(n);
#endif
}

int chklimit(void)
/*Check default upper bounds*/
{
int i;

if(totalloci>maxloci) fprintf(stderr,"set maxloci=%d\n",totalloci);
if(maxall>maxalle) fprintf(stderr,"set maxalle=%d\n",maxall);
/*maxhap is obsolete*/
/*
lh=lg=1;
for(i=0;i<nloci;i++) if(sel[i]) {
  lh*=alleles[i];
  lg*=alleles[i]*(alleles[i]+1)/2;
}
if(cc) lh*=2;
*/

return 0;
}

int nph0=0,nph1=0,nph2=0,malleles[MAX_LOC],cnph0=0,cnph1=0,cnph2=0;

void getdf(void)
/*obtain correct degrees of freedom*/
{
int i,j,n;
if(cc) {
  n=0;
  for(i=1;i<totalloci;i++) n+=loci[i]-1;
  nph0=n;
  n=1;
  for(i=1;i<totalloci;i++) n*=loci[i];
  n--;
  nph1=n;
  nph2=2*n;
  if(onemark) nph1=2*n;
} else { /*todo - reduce df by data*/
  n=0;
  for(i=0;i<totalloci;i++) n+=loci[i]-1;
  nph0=n;
  n=1;
  for(i=0;i<totalloci;i++) n*=loci[i];
  n--;
  nph1=n;
  if(addum) nph0=nph1=loci[0]-1;
}
/*correct missing alleles here*/
for(i=0;i<totalloci;i++) {
  n=0;
  for(j=0;j<loci[i];j++) if(p[i][j]<eps) n++;
  malleles[i]=n;
}
if(cc) {
  n=0;
  for(i=1;i<totalloci;i++) n+=(loci[i]-malleles[i])-1;
  cnph0=n;
  n=1;
  for(i=1;i<totalloci;i++) n*=(loci[i]-malleles[i]);
  n--;
  cnph1=n;
  cnph2=2*n;
  if(onemark) nph1=2*n;
} else {
  n=0;
  for(i=0;i<totalloci;i++) n+=(loci[i]-malleles[i])-1;
  cnph0=n;
  n=1;
  for(i=0;i<totalloci;i++) n*=(loci[i]-malleles[i]);
  n--;
  cnph1=n;
  if(addum) cnph0=cnph1=loci[0]-malleles[0]-1;
}
}

void outH(int *locia, double *a)
/*output haplotype frequencies and likelihoods*/
{
int times,i,line,n;
double oi,ei,ftdev;
times=n=1;
for(i=0;i<totalloci;i++){
  n*=locia[i];locip[i]=1;
}
fprintf(fout, "There are %d Possible Haplotypes of These %d Loci.\n",n, totalloci);
fprintf(fout, "They are Listed Below, with their Estimated Frequencies:\n");putc('-', fout);
for(i=0;i<totalloci;i++) fprintf(fout,"---------");
if(cc&&!onemark)
{
for(i=0;i<41;i++) putc('-',fout);fprintf(fout, "\n|");
for(i=0;i<totalloci;i++) fprintf(fout, " Allele  ");
fprintf(fout, "|         Haplotype   Frequency         |\n");putc('|',fout);
for(i=0;i<totalloci;i++) fprintf(fout, "   at    ");
fprintf(fout, "|                                       |\n");putc('|',fout);
}
else
{
for(i=0;i<33;i++) putc('-',fout);fprintf(fout, "\n|");
for(i=0;i<totalloci;i++) fprintf(fout, " Allele  ");
fprintf(fout, "|      Haplotype Frequency      |\n");putc('|', fout);
for(i=0;i<totalloci;i++) fprintf(fout, "   at    ");
fprintf(fout, "|                               |\n");putc('|', fout);
}
if(cc) {
  fprintf(fout, " Disease ");
  for(i=1;i<totalloci;i++)
    if(totalloci>9) fprintf(fout, " Marker%2d", i);else fprintf(fout, " Marker%d ", i);
} else for(i=1;i<=totalloci;i++)
  if(totalloci>9) fprintf(fout, " Locus %2d", i);else fprintf(fout, " Locus %d ", i);
if(cc&&!onemark) fprintf(fout, "|  Independent   Ind-Disease   w/Asso.  |\n");
else fprintf(fout, "|  Independent   w/Association  |\n");putc('-', fout);
for(i=0;i<totalloci;i++) fprintf(fout, "---------");
if(cc&&!onemark) for(i=0;i<41;i++) putc('-', fout);
else for(i=0;i<33;i++) putc('-', fout);
if(!delta) putc('\n', fout);
else fprintf(fout,"     Observed     Expected    Freeman-Tukey Zi\n");
do {
  line=linenum(locia,locip,1);
  for(i=1;i<=totalloci;i++) {
    if(cc&&i==1) {
      if(locip[i-1]==1) fprintf(fout, "%4c     ", '+');else fprintf(fout, "%4c     ", 'D');
    } else fprintf(fout, "%4d     ", locip[i-1]);
  }
  if(cc&&!onemark) fprintf(fout, "      %8.6f      %8.6f     %8.6f",
            inihaplo[line-1],indhaplo[line-1],a[line-1]);
  else fprintf(fout, "      %8.6f     %8.6f",inihaplo[line-1],a[line-1]);
  if(!delta) fprintf(fout,"\n");
  else {
    if(cc&&!onemark) {
      ei=indhaplo[line-1]*s2;
      oi=a[line-1]*s2;
      ftdev=sqrt(oi)+sqrt(oi+1)-sqrt(4.0*ei+1);
      fprintf(fout, "     %8.2f      %8.2f     %8.2f",oi,ei,ftdev);
    }
    else {
      ei=inihaplo[line-1]*s2;
      oi=a[line-1]*s2;
      ftdev=sqrt(oi)+sqrt(oi+1)-sqrt(4.0*ei+1);
      fprintf(fout, "           %8.2f     %8.2f     %8.2f",oi,ei,ftdev);
    }
    if(fabs(ftdev)<5.0) fprintf(fout,"\n");else fprintf(fout," *\n");
  }
  locip[totalloci-1]++;
  for(i=totalloci;i>=2;i--)
    if(locip[i-1]>loci[i-1]){
      locip[i-1]=1;
      locip[i-2]++;
    }
  times++;
} while (times<=n);putc('-', fout);
for(i=0;i<totalloci;i++) fprintf(fout, "---------");
if(cc && !onemark) fprintf(fout, "-----------------------------------------\n");
else fprintf(fout, "---------------------------------\n");
fprintf(fout, "# of Iterations = %d\n\n", loop);
fprintf(fout,"                                       np   Ln(L)     Chi-square");
fprintf(fout,"\n");
fprintf(fout,"-------------------------------------------------------------------\n");
if(cc) {
  n=0;
  for(i=1;i<totalloci;i++) n+=loci[i]-1;
  fprintf(fout,"H0: No Association                     %2d  %8.2f      0.00\n",n,iniloglike);
  n=1;
  for(i=1;i<totalloci;i++) n*=loci[i];
  n--;
  if(!onemark) {
    if(fabs(indloglike-iniloglike)<eps)
      fprintf(fout,"H1: Markers Asso., Indep. of Disease   %2d  %8.2f      0.00",n,indloglike);
    else fprintf(fout,"H1: Markers Asso., Indep. of Disease   %2d  %8.2f %8.2f\n",
              n,indloglike,2*(indloglike-iniloglike));
    if(fabs(loglike-iniloglike)<eps)
      fprintf(fout,"H2: Markers and Disease Associated     %2d  %8.2f      0.00",n*2,loglike);
    else fprintf(fout,"H2: Markers and Disease Associated     %2d  %8.2f %8.2f\n",
              n*2,loglike,2*(loglike-iniloglike));
  } else {
    if(fabs(loglike-iniloglike)<eps)
      fprintf(fout,"H1: Markers and Disease Associated     %2d  %8.2f      0.00",n*2,loglike);
    else fprintf(fout,"H1: Markers and Disease Associated     %2d  %8.2f %8.2f\n",
              n*2,loglike,2*(loglike-iniloglike));
  }
} else {
  n=0;
  for(i=0;i<totalloci;i++) n+=loci[i]-1;
  fprintf(fout,"H0: No Association                     %2d  %8.2f      0.00\n",n,iniloglike);
  n=1;
  for(i=0;i<totalloci;i++) n*=loci[i];
  n--;
  fprintf(fout, "H1: Allelic Associations Allowed       %2d  %8.2f",n, loglike);
  if(fabs(loglike-iniloglike)<eps) fprintf(fout, "      0.00\n");
  else fprintf(fout, "%10.2f\n", 2*(loglike-iniloglike));
} putc('\n', fout);
}

void outP(int totalloci)
/*output allele frequency estimates*/
{
int i,j;
fprintf(fout, "Estimates of Gene Frequencies (Assuming Independence)\n");
if(cc) fprintf(fout, "(Disease gene frequencies are user specified)\n");
fprintf(fout, "----\\------------");
for(i=1; i <=maxall; i++) fprintf(fout, "--------");
fprintf(fout, "\nlocus \\ allele ");
for(i=1; i <=maxall; i++) fprintf(fout, "%3c%2d%3c", ' ', i, ' ');
fprintf(fout, "\n--------\\--------");
for(i=1; i <=maxall; i++) fprintf(fout, "--------");
putc('\n', fout);
if(cc) {
  fprintf(fout, "Disease |%7c", ' ');
  for(j=0;j<maxall; j++)
    if(p[0][j]<eps) fprintf(fout, "%8c", ' '); else fprintf(fout, "%8.4f", p[0][j]);
  putc('\n', fout);
  for(i=1;i<totalloci;i++) {
    fprintf(fout, "%4d    |%7c", i, ' ');
    for(j=0;j<maxall;j++)
      if(p[i][j]<eps) fprintf(fout, "%8c", ' ');else fprintf(fout, "%8.4f", p[i][j]);
    putc('\n', fout);
  }
} else {
  for(i=0;i<totalloci;i++) {
    fprintf(fout, "%4d    |%7c", i+1, ' ');
    for(j=0;j<maxall;j++)
      if(p[i][j]<eps) fprintf(fout, "%8c", ' ');else fprintf(fout, "%8.4f", p[i][j]);
    putc('\n', fout);
  }
}
fprintf(fout, "-----------------");
for(i=0;i<maxall;i++) fprintf(fout, "--------");
fprintf(fout, "\n# of Typed Individuals: %2ld\n\n",(long)floor(totalind + 0.5));
}

double probnorm(double),probchi(double, int);

double normapprox(double x2, int df, int op)
/*
 * Manoukian E.B. (1986)
 * Moder Concepts and Theorems of Mathematical Statistics
 * Springer p137
*/
{
double x,m,v,t;

if(df>30) switch (op) {
case 1: /*moderate*/
     x=2.0*x2;
     m=2.0*df-1.0;
     v=1;
     break;
case 2: /*better*/
     x=pow(x2/df,1/3);
     t=(2/9/df);
     m=1-t;
     v=t;
     break;
case 3: /*by chi-square, worst*/
     x=x2;
     m=df;
     v=2*df;
     break;
default:break;
     fprintf(stderr,"no normal approximation for such option\n");
}
return (x-m)/sqrt(v);
}


/*LD code*/

#define minval(a,b) ((a<b)?a:b)

short locus1,locus2;
short locik[maxloci], lociq[maxloci];

void ehtest(void)
/*
  to check for the genotype/allele indices
  corrected 17/8/1999, 23/8/1999, 13/5/2000 JH Zhao
*/
{
int genonum,i,time;

genonum = 1;
totalloci = maxloci;
for(i=0;i<totalloci;i++) {
    locip[i] = 1;
    locik[i] = 1;
    lociq[i] = 1;
    locigeno[i] = loci[i] * (loci[i] + 1) / 2;
    genonum = genonum * locigeno[i];
}
for(time=0;time<genonum;time++) {
    for(i=0;i<totalloci;i++) {
        printf(" %1d/%1d[%2d]",locik[i],lociq[i],locip[i]);
        printf(" %2d",locip[i]);
    } printf("\n");
    locip[totalloci-1]++;
    for(i=totalloci-1;i>0;i--) {
        if(locip[i]>locigeno[i]) {
           locip[i]=1;
           locip[i-1]++;
        }
    }
    locik[totalloci-1]++;
    for(i=totalloci-1;i>0;i--) {
        if(locik[i]>lociq[i]) {
           lociq[i]++;
           locik[i]=1;
           if(lociq[i]>loci[i]) {
              lociq[i]=1;
              locik[i-1]++;
              if(i==1) if(locik[0]>lociq[0]) {
                lociq[0]++;
                locik[0]=1;
              }
           }
        }
    }
}
}

int linenum2(int *locia, short *locip, int sp, int subset)
/*the genotype/haplotype/allele identifier*/
{
int sum,j;
sum=0;
for(j=sp;j<=subset;j++)
  if(j==subset) sum+=locip[j-1]; else sum=(sum+locip[j-1]-1)*locia[j];
return(sum);
}

void getHij(double *hij, int locusi, int locusj)
/*get 2-locus haplotype frequencies from all haplotype frequencies*/
{
int line, line2, i, k, n, times, loci2[maxloci], l1, l2;

k=0;
n=1;
for(i=0;i<nloci;i++) {
  if(!sel[i]) continue;
  loci2[k]=alleles[i];
  locip[k]=1;
  n*=alleles[i];
  if(i==locusi) l1=k;
  if(i==locusj) l2=k;
  k++;
}
for(i=0;i<alleles[locusi]*alleles[locusj];i++) hij[i]=0;
times=0;
do {
  line=linenum2(loci2,locip,1,k)-1;
  line2=(locip[l1]-1)*loci2[l2]+locip[l2]-1;
  hij[line2]+=newhaplo[line];
/*
  for(i=0;i<k;i++) printf(" %d",locip[i]);
  printf(" %d %f -> h[%d]\n",line,newhaplo[times],line2);
*/
  locip[k-1]++;
  for(i=k-1;i>0;i--) {
    if(locip[i]>loci2[i]) {
      locip[i]=1;
      locip[i-1]++;
    }
  }
  times++;
} while(times<n);
}

int getp(double *hij)
{
double t,pi[maxalleles],qj[maxalleles];
int i,j;

for(i=0;i<alleles[locus1];i++) {
  t=0;
  for(j=0;j<alleles[locus2];j++) {
    t+=hij[i*alleles[locus2]+j];
  }
  pi[i]=t;
}
for(j=0;j<alleles[locus2];j++) {
  t=0;
  for(i=0;i<alleles[locus1];i++) {
    t+=hij[i*alleles[locus2]+j];
  }
  qj[j]=t;
}

return 0;
}
/*

Zapata C, Alvarez G, Carollo C (1997) Approximate variance of the standardized
measure of gametic disequilibrium D'. Am. J. Hum. Genet. 61:771-774

Weir, BS (1996) Genetic Data Analysis II. Sinauer

for a sample of n haplotypes

     B1  B2
   +--------+
 A1| x0  x1 | p
 A2| x2  x3 | q
   +--------+
     u   v

 */

double D,Dmax,Dprime,VarD,VarDmax,VarDprime,x22;

int twobytwo(double *hij)
{
int i;
double xx[4];
double p,q,u,v;
double a,b,xxi,t;
double ED,EDmax;
double haplotypes;

haplotypes=sample_size;
haplotypes*=2;/*is actually # of haplotypes*/
for(i=0;i<4;i++) xx[i]=hij[i];
p=xx[0]+xx[1];
q=xx[2]+xx[3];
u=xx[0]+xx[2];
v=xx[1]+xx[3];

D=xx[0]-p*u;
ED=D*(haplotypes-1)/haplotypes;
VarD=(p*q*u*v+D*(q-p)*(v-u)-D*D)/haplotypes;
if(D<0) {
  if(p*u<q*v) {
    Dmax=p*u;
    xxi=xx[0];
  } else {
    Dmax=q*v;
    xxi=xx[3];
  }
} else {
  if(p*v<q*u) {
    Dmax=p*v;
    xxi=xx[1];
  } else {
    Dmax=q*u;
    xxi=xx[2];
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
  VarDprime=((1-fabs(Dprime))*t+fabs(Dprime)*xxi*(1-xxi))/haplotypes/Dmax/Dmax;
}
x22=haplotypes*D*D/p/q/u/v;

return 0;
}

/*

Klitz W, Stephen JC, Grote M, Carrington M (1995) Discordant patterns of
linkage disequilibrium of the peptide transporter loci within the HLA class II
region. Am. J. Hum. Genet. 57:1436-1444

Long JC, Williams RC, Unbanek M (1995) An E-M  algorithm and testing strategy
for multiple-locus haplotypes. Am. J. Hum. Genet. 56:799-810

*/

double W,Wn,x2D,pD;
int dfD;

void kbyl(double *hij)
/*with adjustments to d.o.f*/
{
double Dij,Dmax,t,n;
int i,j,k,l,l1,l2,s1,s2;

n=sample_size;
l1=locus1;
l2=locus2;
k=alleles[l1];
l=alleles[l2];
j=0;
for(i=0;i<nloci;i++) {
  if(!sel[i]) continue;
  if(i==locus1) l1=j;
  if(i==locus2) l2=j;
  j++;
}
s1=s2=0;
for(i=0;i<k;i++) if(p[l1][i]>eps) s1++;
for(j=0;j<l;j++) if(p[l2][j]>eps) s2++;
t=0;
for(i=0;i<k;i++) {
   if(p[l1][i]==0) continue;
   for(j=0;j<l;j++) {
      Dij=hij[i*l+j]-p[l1][i]*p[l2][j];
      if(p[l2][j]==0) continue;
      t+=Dij*Dij/p[l1][i]/p[l2][j];
      if(Dij<=0) Dmax=minval(p[l1][i]*p[l2][j],(1-p[l1][i])*(1-p[l2][j]));
      else Dmax=minval(p[l1][i]*(1-p[l2][j]),p[l2][j]*(1-p[l1][i]));
   }
}
x2D=t*n*2;
W=sqrt(t);
t=(s1<s2)?s1:s2;
Wn=W/sqrt(t-1);
dfD=(s1-1)*(s2-1);
pD=probchi(x2D,dfD);

}

void ld2(FILE *fp)
/*pairwise LD measures*/
{
double hij[maxalleles*maxalleles];
int i,j,k;

fprintf(fp,"\nPairwise disequilibrium measures from a sample of %d individuals\n",sample_size);
for(k=0;k<72;k++) fprintf(fp,"-");fprintf(fp,"\n");
fprintf(fp,"  loci       W       Wn     X^2   df       p    D/SE   Dmax/SE    D'/SE\n");
for(k=0;k<72;k++) fprintf(fp,"-");fprintf(fp,"\n");
for(i=0;i<=selected;i++) {
  for(j=i+1;j<=selected;j++) {
    locus1=selidx[i];
    locus2=selidx[j];
    if(nloci>2) getHij(hij,locus1,locus2);
    else for(k=0;k<alleles[locus1]*alleles[locus2];k++) hij[k]=newhaplo[k];
/*
    getp(hij);
*/
    kbyl(hij);
    fprintf(fp,"%3d%3d %f %f %7.2f %4d %.5f",locus1+1,locus2+1,W,Wn,x2D,dfD,pD);
    if(alleles[locus1]==2&&alleles[locus2]==2) {
      twobytwo(hij);
      fprintf(fp," %.5f %.5f %.5f\n\t\t\t\t\t      ",D,Dmax,Dprime);
      fprintf(fp," %.5f  %.5f  %.5f\n",sqrt(VarD),sqrt(VarDmax),sqrt(VarDprime));
    } else fprintf(fp,"\n");
  }
}
for(k=0;k<72;k++) fprintf(fp,"-");fprintf(fp,"\n");
fprintf(fp,"W = Global disequilibrium test statistic (D)\n");
fprintf(fp,"Wn = Global normalized disequilibrium, or Cramer's V \n\n");
fprintf(fp,"Note haplotype frequencies are based on joint analysis of all selected markers\n");
fprintf(fp,"Selecting only marker pairs would normally be more informative\n");
fprintf(fp,"\n");
fprintf(fp,"Dmax, S.E.(Dmax), and D', S.E.(D') only given for two biallelic loci\n");
fprintf(fp,"Degrees of freedom have been adjusted based on alleles observed\n\n");
}


/*CP code*/

/*

to generate all subset in natural order, JH Zhao 14/5/2000 IoP

Xiru Chen, Songgui Wang (1987). Modern Regression Analysis. Anhui
Educational Publishing Co, pp203-206

Furnival GM, Wilson RW Jr (1974) Regression by leaps and bounds.
Technometrics 16: 499-511

*/

void cp_natural(FILE *fp)
{
int ind[MAX_LOC], name[MAX_LOC];
int iz, is, ia, ib, im, ip, n0, n1, j;
double l0, l1, x2;
double hl0[3], hl1[3], hl2[3];
int hn0[3], hn1[3], hn2[3];

is = 0;
for (iz = 0; iz < nloci; iz++) {
    if (!sel[iz]) continue;
    name[is] = iz;
    is++;
}
fprintf(fp,"\nAll-subset disequilibrium statistics from a sample of %d individuals\n",sample_size);
for(iz = 0; iz < 72; iz++) fprintf(fp,"-");fprintf(fp,"\n");
fprintf(fp,"subsets           H0  n0           H1  n1         x^2   df      p\n");
for(iz = 0; iz < 72; iz++) fprintf(fp,"-");fprintf(fp,"\n");
ip=im=selected+1;
for (iz = 0; iz < MAX_LOC; ++iz) ind[iz] = sel[iz] = selp[iz] = 0;
do {
   for (iz = im - 1; iz < ip; ++iz) {
       if (ind[iz] < iz + 1) continue;
       ++ind[iz - 1];
       ind[iz] = ind[iz - 1];
   }
   do {
      ++ind[ip - 1];
      for (iz = im - 1; iz < ip; ++iz) sel[name[ind[iz] - 1]] = true;
      fprintf(fp,"%d", name[ind[im - 1] - 1] + 1);
      for (iz = im; iz < ip; ++iz) fprintf(fp,"-%d", name[ind[iz] - 1] + 1);
      fprintf(fp,"\n");
      selected = selectn = ip - im;
      selectp = -1;
      is = ia = ib = 0;
      for (iz = 0; iz < nloci; iz++) {
          if (!sel[iz]) continue;
          selidx[is] = iz;
          if (selp[iz]) {
             selpidx[ia]= iz;
             ia++;
          } else {
             selnpdx[ib] = iz;
             ib++;
          }
        is++;
      }
      ind2eht();
/*
      for (iz = im - 1; iz < ip; ++iz) sel[ind[iz] - 1] = 1;
      fprintf(fp,"%d", ind[im - 1]);
      for (iz = im; iz < ip; ++iz) fprintf(fp,"-%d", ind[iz]);
      fprintf(fp,"\n");
      fprintf(fp,"%d", sel[0]);
      for (iz = 1; iz < nloci; ++iz) fprintf(fp,"-%d", sel[iz]);
      fprintf(fp,"\n");
      printf("[%d %d %d]\n",selected,selectp,selectn);
      ptree(rt,0,stdout,WHOLE);
      ptree(rt2,0,stdout,BLOCK2);
*/
      obscom = 0;
      ctree(rt, WHOLE);
/*
      fout=stdout;
      outP(totalloci);outH(loci,newhaplo);
*/
      if(!cc) {
/*all data*/
        eh(WHOLE);
        l0=iniloglike;
        l1=loglike;
        n0=(correct_df)?cnph0:nph0;
        n1=(correct_df)?cnph1:nph1;
        x2=2*(l1-l0);
        fprintf(fp,"\t");
        fprintf(fp,"%12.2f %3d %12.2f %3d %12.2f %3d %.5f\n",l0,n0,l1,n1,x2,n1-n0,probchi(x2,n1-n0));
      } else {
/*all data*/
        cc=false;
        eh(WHOLE);
        l0=hl0[0]=iniloglike;
        l1=hl0[1]=loglike;
        n0=hn0[0]=(correct_df)?cnph0:nph0;
        n1=hn0[1]=(correct_df)?cnph1:nph1;
        x2=2*(l1-l0);
        fprintf(fp,"\t");
        fprintf(fp,"%12.2f %3d %12.2f %3d %12.2f %3d %.5f\n",l0,n0,l1,n1,x2,n1-n0,probchi(x2,n1-n0));
/*cases*/
        totalind=0;
        for(j=0;j<obscom;j++) {
          indivN[j]=caseN[j];
          totalind+=caseN[j];
        }
        eh(WHOLE);
        l0=hl1[0]=iniloglike;
        l1=hl1[1]=loglike;
        n0=hn1[0]=(correct_df)?cnph0:nph0;
        n1=hn1[1]=(correct_df)?cnph1:nph1;
        x2=2*(l1-l0);
        fprintf(fp,"\t");
        fprintf(fp,"%12.2f %3d %12.2f %3d %12.2f %3d %.5f\n",l0,n0,l1,n1,x2,n1-n0,probchi(x2,n1-n0));
/*controls*/
        totalind=0;
        for(j=0;j<obscom;j++) {
          indivN[j]=contrlN[j];
          totalind+=contrlN[j];
        }
        eh(WHOLE);
        l0=hl2[0]=iniloglike;
        l1=hl2[1]=loglike;
        n0=hn2[0]=(correct_df)?cnph0:nph0;
        n1=hn2[1]=(correct_df)?cnph1:nph1;
        x2=2*(l1-l0);
        fprintf(fp,"\t");
        fprintf(fp,"%12.2f %3d %12.2f %3d %12.2f %3d %.5f\n",l0,n0,l1,n1,x2,n1-n0,probchi(x2,n1-n0));
/*heterogeneity*/
        l0=l1=0;
        n0=n1=0;
        x2=-2*(hl0[1]-hl2[1]-hl1[1]);
        j=(hn1[1]+hn2[1]-hn0[1]);
        fprintf(fp,"\t");
        fprintf(fp,"%12.2f %3d %12.2f %3d %12.2f %3d %.5f\n",l0,n0,l1,n1,x2,j,probchi(x2,j));
        cc=true;
      }
      for (iz = 0; iz < orig_model.nloci; ++iz) sel[iz] = 0;
   }  while (ind[ip - 1] < ip);
   if (ind[im - 1] == im) --im;
}  while (im > 0);
for(iz = 0; iz < 72; iz++) fprintf(fp,"-");fprintf(fp,"\n");
fprintf(fp,"Single marker results are from adding dummy marker with no effect\n");
if(cc)
fprintf(fp,"Statistics from pooled data, case only, control only and heterogeneity test\n\n");
fprintf(fp,"H0 -- assuming no association\n");
fprintf(fp,"H1 -- association\n");
fprintf(fp,"n0, n1 -- number of parameters\n");
fprintf(fp,"x^2, df, p -- LRT chi-square, degrees of freedom and p-value\n\n");

/*reset original model, selidx, etc won't work !*/
nloci_s = nloci = orig_model.nloci_s;
selected = orig_model.selected;
selectp = orig_model.selectp;
selectn = orig_model.selectn;
for (iz = 0; iz < nloci; iz++) {
    sel[iz] = orig_model.sel[iz];
    selp[iz] = orig_model.selp[iz];
}
is = ia = ib = 0;
for (iz = 0; iz < nloci; iz++) if (sel[iz]) {
    selidx[is] = iz;
    if (selp[iz]) {
       selpidx[ia] = iz;
       ia++;
    } else {
       selnpdx[ib] = iz;
       ib++;
    }
    is++;
}
}


/*code for dynamic allocation*/

/*
 * pass various tests after s/t initialized
 * JH Zhao 16-17/5/2000
 */

struct obs_type {
  individual anybody;
  struct obs_type *next;
};
typedef struct obs_type obs;
typedef obs *link;

link ilink(individual new_ind,link t)
/*insert*/
{
/*
 * list from t to s, head h is maintained
 */
  link s,p1,p2,h;

  s=(obs *)malloc(sizeof(obs));
  s->anybody=new_ind;
  s->next=NULL;
  p1=t;
  if(p1==NULL) h=s;
  else {
    h=p1;
    do {
      p2=p1;
      p1=p1->next;
    } while (p1!=NULL);
    p2->next=s;
  }
  return h;
}

void clink(link h)
/*copy*/
{
  link t;
  int valids;

  valids=0;
  if(h==NULL) {
    fprintf(stderr,"Empty list.\n");
  } else {
    while(h!=NULL)
    {
      person_s[valids]=person[valids]=h->anybody;
      t=h;
      h=h->next;
      valids++;
      free(t);
    }
  }
}

void dirichlet(long int n)
/*9-10-2001 (IoP)*/
/*randu and 1-randu are both U(0,1)*/
{
time_t t;
double s,randu;
long int i;

/*srand((unsigned)time(&t));*/
srand(seed+iter);
s=0;
for(i=0;i<n;i++)
{
  randu=(double)rand()/RAND_MAX;
  if(randu<eps) randu=eps;
  inihaplo[i]=-log(randu);
  s+=inihaplo[i];
}
for(i=0;i<n;i++) inihaplo[i]/=s;
}

#endif
