/*                                                           */
/* GENECOUNTING - haplotype analysis with missing genotypes  */
/* (C) Copyright JH Zhao 2001-2003 University College London */
/*                                                           */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define version 2.2
double pl=0.001;
double epsh=0.000001;

/*global definitions*/
#ifdef DEBUG
int debug=1;
#else
int debug=0;
#endif
enum {NO, YES};

#ifdef NOHM
  #define handle_missing NO
#else
  #define handle_missing YES
#endif

enum {by_l1_norm,by_ll};
#ifndef L1
  #define conv_ll by_ll
  #define eps 0.01
#else
  #define conv_ll by_l1_norm
  #define eps 0.01
#endif

#ifdef Y
  #define read_y YES
#else
  #define read_y NO
#endif

#ifdef HTR
  #define htr YES
#else
  #define htr NO
#endif

#define mxloc 60
#define mxalleles 50
#define maxit 50

double precis=1.0,tol;

FILE *fo;
int loci[mxloc];
int nloci,nalleles;
int obscom;
long int hapall,nhaplotypes;
double d_sample,d_msample;
double p[mxloc][mxalleles];
double *c,*h,*hs,*h0;
void *xmalloc(long);

typedef struct
{
  long int id;
  double le,ld;
  short l[mxloc];
} hap;

hap *haplotypes,haplotype_t;

typedef struct
{
  long gid;
  short sex;
  double count,prob,y;
  union {
  short h[mxloc];
  short d[mxloc][2];
  } hd;
} pat;

typedef struct patp_type
{
  pat anyline;
  struct patp_type *next;
} patp, *patlink;

patlink pilink(pat,patlink);
void pclink(patlink,pat*);

typedef struct hnode_type
{
  long int id;
  int n;
  struct hnode_type *left, *right;
  short l[mxloc];
} hnode;

hnode *hrt=0;
hnode *hitree(hnode*,long int,short*);
void hrtree(hnode*),hptree(hnode*,long int*);

/*ordinary data*/

void ehm(pat *),inith(long),hilo(short*,short*);
void getp(pat *),geth(pat *),ch(pat *),counting(pat *, long int);
void genp(pat *,long int);
void digit2(int,short*,int),digitm(short*,short*,int);
int linenum(int*,short*);
double phasep(short*,short*);
double ll(pat *);

/*missing data*/

int linenums(int*,short*);

int nlocim,locim[mxloc+1];
int idm[mxloc];
long int hapallm;
double *hm;

/*HTR of Zaykin et al. (2002) Hum Hered 53:79-91*/

double **htr_table;
double  *htr_beta;

/*X chromosome data*/
int nhaploid,ndiploid,h_sample,h_msample;
double pm[mxloc][mxalleles];

#ifdef X
int xdata=1;
#else
int xdata=0;
#endif
enum {MALE,FEMALE};

int main(int argc,char *argv[])
{
int gcmain(char*);
time_t t;

printf("\nGENECOUNTING, %s (%s) version %.2f \n\t\tJH Zhao 03/01--07/05\n\n",
(!handle_missing)?"ordinary":"missing-value handling",
(!xdata)?"autosome data":"X chromosome data",version);
time(&t);
printf("%s\n\n",ctime(&t));
printf("max. loci=%d, max. alleles=%d\n\n",mxloc,mxalleles);

/*     PRECIS    Machine precision ( changed on different machine )       */
/*               But it is automatically ascertained in revised GEMINI    */
/*               according to Forsythe,G.E., Malcolm, M. A. and Moler,    */
/*               C. B. (1977). Computer mathods for mathematical          */
/*                computations. Prentice-Hall.                            */

do
{
  precis/=2.0;
  tol=1.0+precis;
} while(tol>1.0);
fo=stdout;
if(argc>1)
{
#ifdef DEBUG_HAPLOID
  xehm(argv[1]);
  return 0;
#endif
  if(argc>2)
  {
    fo=fopen(argv[2],"w");
    if(!fo)
    {
      fprintf(stderr,"I can not open output file %s,\n",argv[2]);
      fprintf(stderr,"so I will write to the screen\n");
    }
    else
    {
      fprintf(fo,"The input/output files are: %s/%s\n\n",argv[1],argv[2]);
      if(argc>3)
      {
       pl=atof(argv[3]);
       if(pl<precis) pl=precis;
      }
    }
  }
  gcmain(argv[1]);
  fclose(fo);
}
else
{
  printf("Usage: %s <input file> [outfile] [threshold] [>screen file]\n\n",argv[0]);
  printf("<input file> contains the following information:\n\n");
  printf("line 1: <#alleles at locus 1> <#alleles at locus 2> ...\n");
  printf("line 2- genotype identifier, instances, marker phenotypes\n\n");
  printf("[threshold] is cutoff value for posterior probability (%.6f)\n",pl);
}

return 0;
}

int gcmain(char *gcdat)
{
FILE *fp;
long int i,j,k,l;
char line[201],rest[201];
short loci1[mxloc],d[mxloc];
short sex,la,ua;
double v0,v1,sn,ft,s0,s1;
pat *table,tt;
patlink s,t;

void qsorts(long int,long int);

fp=fopen(gcdat,"r");
if(!fp)
{
  fprintf(stderr,"Sorry, is there file called %s ?\n",gcdat);
  exit(1);
}
fgets(line,200,fp);
nloci=0;
while(sscanf(line,"%d %[^\n]",&loci[nloci++],rest)>1) strcpy(line,rest);
hapall=1;
nalleles=1;
for(i=0;i<nloci;i++)
{
  hapall*=loci[i];
  if(loci[i]>nalleles) nalleles=loci[i];
}
inith(hapall);
tt.gid=0;
tt.count=0;
tt.sex=0;
for(i=0;i<mxloc;i++) tt.hd.h[i]=tt.hd.d[i][0]=tt.hd.d[i][1]=0;
s=t=NULL;
nhaploid=ndiploid=0;
i=l=0;
while(fgets(line,200,fp)&&
sscanf(line,"%ld %lf %[^\n]",&tt.gid,&tt.count,rest)>2)
{
  strcpy(line,rest);
  if(xdata)
  {
     sscanf(line,"%hd %[^\n]",&tt.sex,rest);
     strcpy(line,rest);
  }
  k=0;
  for(j=0;j<nloci;j++)
  {
    if (!xdata)
    {
       if((sscanf(line,"%hd/%hd %[^\n]",&la,&ua,rest)<2)&&
          (sscanf(line,"%hd %hd %[^\n]",&la,&ua,rest)<2))
       {
         fprintf(stderr,"There isn't enough data at line %d\n",i+1);
         return 1;
       }
       if ((la>loci[j]||ua>loci[j])||(la<=0||ua<=0)) ++k;
    }
    else
    {
       if ((tt.sex==FEMALE&&sscanf(line,"%hd/%hd %[^\n]",&la,&ua,rest)<2 &&
                            sscanf(line,"%hd %hd %[^\n]",&la,&ua,rest)<2)||
           (tt.sex==MALE&&sscanf(line,"%hd %[^\n]",&la,rest)<1))
       {
         fprintf(stderr,"There isn't enough x data at line %d\n",i+1);
         return 1;
       }
       if (tt.sex==FEMALE&&((la>loci[j]||ua>loci[j])||(la<=0||ua<=0))) ++k;
       if (tt.sex==MALE&&(la>loci[j]||la<=0)) ++k;
    }
    strcpy(line,rest);
    *rest='\0';
    if(xdata&&tt.sex==MALE)
    tt.hd.h[j]=la;
    else
    {
      if (la>ua) hilo(&la,&ua);
      tt.hd.d[j][0]=la;
      tt.hd.d[j][1]=ua;
    }
  }
  if(read_y)
  {
    if(sscanf(line,"%lf %[^\n]",&tt.y,rest)<1)
    {
      fprintf(stderr,"There isn't response variable at line %d\n",i+1);
      return 1;
    }
  }
  else tt.y=0;
  if(k<nloci)
  {
    if(xdata)
    {
      if(tt.sex==MALE) nhaploid++;
      else ndiploid++;
    }
    t=pilink(tt,s);
    s=t;
    ++i;
  }
  else ++l;
}
/*if(l>0) fprintf(stderr,"%d data lines have been skipped\n",l);*/
fclose(fp);
obscom=i;
table=(pat*)xmalloc(obscom*sizeof(pat));
pclink(s,table);
for(i=0;i<obscom;i++) if(debug)
{
  printf("%5d %.0f %hd",table[i].gid,table[i].count,table[i].sex);
  for(j=0;j<nloci;j++) if(xdata&&table[i].sex==MALE) printf(" %2d",table[i].hd.h[j]);
  else printf(" %2d/%2d",table[i].hd.d[j][0],table[i].hd.d[j][1]);
  printf("\n");
}
ehm(table);
free(hs);
free(c);
if(handle_missing) free(hm);
haplotypes=(hap*)xmalloc(nhaplotypes*sizeof(hap));
for(i=0;i<=nloci;i++) d[i]=0;
for(i=0;i<nloci;i++) loci1[i]=loci[nloci-i-1]-1;
l=0;
for(i=0;i<hapall;i++)
{
  if(h[i]>epsh)
  {
    haplotypes[l].id=i+1;
    haplotypes[l].le=h0[i];
    haplotypes[l].ld=h[i];
    for(j=nloci-1;j>=0;j--) haplotypes[l].l[j]=d[j]+1;
    ++l;
  }
  digitm(loci1,d,0);
}
fprintf(fo,"Haplotype frequency estimates by decreasing order\n\n");
fprintf(fo,"    h1    accumul.    h0    accumul.     FT haplotypes/id\n\n");
qsorts(0,nhaplotypes-1);
sn=2.0*(d_sample+d_msample)+(h_sample+h_msample);
s0=s1=0;
for(i=0;i<nhaplotypes;i++)
{
  fprintf(fo,"");
  v0=haplotypes[i].le;
  v1=haplotypes[i].ld;
  s0+=v0;
  s1+=v1;
  ft=sqrt(v1*sn)+sqrt(v1*sn+1)-sqrt(4*v0*sn+1);
  fprintf(fo," %.6lf %.6lf %.6lf %.6lf %6.2lf",v1,s1,v0,s0,ft);
  for(j=nloci-1;j>=0;j--) fprintf(fo,"%3d",haplotypes[i].l[j]);
  fprintf(fo," %d\n",haplotypes[i].id);
}
free(haplotypes);
htr_beta=(double*)xmalloc((1+hapall)*sizeof(double));
htr_table=(double**)xmalloc(obscom*sizeof(double));
for(i=0;i<obscom;i++)
{
  htr_table[i]=(double*)xmalloc((1+hapall)*sizeof(double));
  htr_table[i][0]=1.0;
  for(j=1;j<=hapall;j++) htr_table[i][j]=0;
}
ch(table);
if(debug)
{
  fprintf(fo,"\nFrequencies of possible haplotypes\n\n");
  l=0;
  hptree(hrt,&l);
  fprintf(fo,"\nPossible haplotypes excluding missing loci = %ld\n",l);
}
hrtree(hrt);
free(h0);
if(htr)
{
  /*call up regression routines here, use SPlus at the moment*/
  fprintf(fo,"\n\nHTR working matrices\n\n");
  fprintf(fo,"   ID    w");
  if(read_y) fprintf(fo,"      y");
  for(j=0;j<hapall;j++) if(h[j]>epsh) fprintf(fo," h%-d",j+1);
  fprintf(fo," sum\n");
  for(i=0;i<obscom;i++)
  {
    fprintf(fo,"%5ld %4.0lf",table[i].gid,table[i].count);
    if(read_y) fprintf(fo," %6.2lf",table[i].y);
    sn=0;
    for(j=1;j<=hapall;j++)
    if (h[j-1]>epsh){
      htr_table[i][j]/=2.0;
      sn+=htr_table[i][j];
      fprintf(fo," %lf",htr_table[i][j]);
    }
    fprintf(fo," %lf\n",sn);
  }
}
free(h);
free(table);
free(htr_beta);
for(i=0;i<nhaplotypes;i++) free(htr_table[i]);
free(htr_table);

return 0;
}

void inith(long obshap)
/*allocate haplotype arrays*/
{
h=(double*)xmalloc(obshap*sizeof(double));
h0=(double*)xmalloc(obshap*sizeof(double));
hs=(double*)xmalloc(obshap*sizeof(double));
c=(double*)xmalloc(obshap*sizeof(double));
if(handle_missing) hm=(double*)xmalloc(obshap*sizeof(double));
}

void ehm(pat *table)
/*estimating haplotypes with missing data*/
{
int iter,i,j,k,l,dfusr[3],dfdat[3];
long int io;
short loci1[mxloc],d[mxloc];
double s,lnl0,lnl1,lnls,x2mm,ft;

for(i=0;i<2;i++) dfusr[i]=dfdat[i]=0;
hapall=1;
k=1;
for(i=0;i<nloci;i++)
{
  hapall*=loci[i];
  dfusr[0]+=loci[i]-1;
  k*=loci[i];
}
dfusr[1]=k-1;
getp(table);
k=1;
for(i=0;i<nloci;i++)
{
  loci1[i]=loci[i];
  for(j=0;j<loci[i];j++) if(p[i][j]<precis) --loci1[i];
  dfdat[0]+=loci1[i]-1;
  k*=loci1[i];
}
dfdat[1]=k-1;
for(io=0;io<obscom;io++) genp(table,io);
lnl0=ll(table);
lnls=lnl0;
fprintf(fo,"\nlog-likelihood assuming linkage equilibrium %.2lf\n",lnl0);
fprintf(fo,"# parameters %d (%d, by # of specified alleles)\n\n",dfdat[0],dfusr[0]);
k=1;
iter=1;
do
{
  for(i=0;i<hapall;i++) hs[i]=h[i];
  geth(table);
  if(!conv_ll)
  {
    s=0;
    for(i=0;i<hapall;i++) s+=fabs(hs[i]-h[i]);
  }
  else
  {
    fprintf(fo,"Iteration %3d, ",k++);
    for(io=0;io<obscom;io++) genp(table,io);
    lnl1=ll(table);
    fprintf(fo,"log-likelihood=%.2lf\n",lnl1);
    s=lnl1-lnls;
    lnls=lnl1;
  }
} while(s>eps&&iter++<maxit);
if(!conv_ll)
{
  for(io=0;io<obscom;io++) genp(table,io);
  lnl1=ll(table);
}
fprintf(fo,"\nlog-likelihood assuming linkage disequilibrium =%.2lf\n",lnl1);
fprintf(fo,"# parameters %d (%d, by # of specified alleles)\n",dfdat[1],dfusr[1]);
x2mm=2*(lnl1-lnl0);
fprintf(fo,"\nTest of marker-marker disequilibrium, chi-square=%.2lf\n",x2mm);
fprintf(fo,"degree(s) of freedom = %d (%d, by # of specified alleles) \n",dfdat[1]-dfdat[0],dfusr[1]-dfusr[0]);
s=2.0*(d_sample+d_msample)+(h_sample+h_msample);
fprintf(fo,"\nHaplotype frequency estimates and Freeman-Tukey statistic (FT)\n\n");
fprintf(fo," n = number of haplotypes (%.0lf) in the sample\n",s);
fprintf(fo,"h1 = frequencies under linkage disequilibrium (>%lf)\n",epsh);
fprintf(fo,"h0 = frequencies under linkage equilibrium\n");
fprintf(fo,"FT = sqrt(h1*n)+sqrt(h1*n+1)-sqrt(4*h0+1)\n");
fprintf(fo,"\n       h1             h1       h0     FT  haplotype/id\n\n");
for(i=0;i<=nloci;i++) d[i]=0;
for(i=0;i<nloci;i++) loci1[i]=loci[nloci-i-1]-1;
l=0;
for(i=0;i<hapall;i++)
{
  if(h[i]>epsh)
  {
    ft=sqrt(h[i]*s)+sqrt(h[i]*s+1)-sqrt(4*h0[i]*s+1);
    fprintf(fo," %.6lf %.12lf %.6lf %6.2lf",h[i],h[i],h0[i],ft);
    for(j=nloci-1;j>=0;j--) fprintf(fo,"%3d",d[j]+1);
    fprintf(fo," %d\n",i+1);
    l++;
  }
  digitm(loci1,d,0);
}
nhaplotypes=l;
fprintf(fo,"\nNumber of nonzero haplotypes = %ld\n\n",l);
}

void getp(pat *table)
/*allele frequencies*/
{
int i,j,k,l,u,iter;
short d[mxloc+1],loci1[mxloc],d1[mxloc+1];
double s[mxloc],n[mxloc],ss;

for(i=0;i<nloci;i++)
{
  s[i]=n[i]=0;
  for(j=0;j<nalleles;j++) p[i][j]=pm[i][j]=0;
}
ss=d_msample=0;
for(i=0;i<obscom;i++)
{
  if(xdata&&table[i].sex==MALE) goto x;
  k=0;
  for(j=0;j<nloci;j++)
  {
    l=table[i].hd.d[j][0]-1;
    u=table[i].hd.d[j][1]-1;
    if((l>=loci[j]||u>=loci[j])||(l<0||u<0))
    {
      ++k;
      continue;
    }
    p[j][l]+=table[i].count;
    p[j][u]+=table[i].count;
    s[j]+=table[i].count;
  }
  if(k==0) ss+=table[i].count;
  else d_msample+=table[i].count;
  continue;
x:
  k=0;
  for(j=0;j<nloci;j++)
  {
    l=table[i].hd.h[j]-1;
    if(l>=0&&l<=loci[j])
    {
      ++k;
      pm[j][l]+=table[i].count;
      n[j]+=table[i].count;
    }
  }
  if(k==nloci) h_sample+=table[i].count;
  else h_msample+=table[i].count;
}
for(i=0;i<nloci;i++) fprintf(fo,"locus %3d # individuals = %.0lf\n",i+1,s[i]+n[i]);
fprintf(fo,"\n%.lf/%.lf individuals with F/P information\n\n",ss+h_sample,d_msample+h_msample);
d_sample=ss;
fprintf(fo,"\nMarker allele frequencies:\n\nMarker\\Allele number\n\n");
fprintf(fo,"   ");
for(i=0;i<nalleles;i++) fprintf(fo,"%9d",i+1);
fprintf(fo,"\n\n");
for(i=0;i<nloci;i++)
{
  fprintf(fo,"%2d ",i+1);
  for(j=0;j<loci[i];j++)
  {
    if(!xdata) p[i][j]/=2*s[i];
    else p[i][j]=(p[i][j]+pm[i][j])/(2*s[i]+n[i]);
    if(p[i][j]<precis) fprintf(fo,"         ");
    else fprintf(fo," %.6lf",p[i][j]);
  }
  fprintf(fo,"\n");
  loci1[i]=loci[i]-1;
}
for(i=0;i<=nloci;i++) d[i]=0;
for(i=0;i<hapall;i++)
{
  ss=1;
  for(j=0;j<nloci;j++) ss*=p[j][d[j]];
  for(j=0;j<=nloci;j++) d1[j]=d[j]+1;
  l=linenum(loci,d1)-1;
  h[l]=h0[l]=ss;
  digitm(loci1,d,0);
}
if(debug)
{
  fprintf(fo,"\nHaplotype frequencies under linkage equilibrium\n\n");
  for(i=0;i<=nloci;i++) d[i]=0;
  for(i=0;i<nloci;i++) loci1[i]=loci[nloci-i-1]-1;
  for(i=0;i<hapall;i++)
  {
    fprintf(fo,"%.6lf [%.12lf]",h[i],h[i]);
    for(j=nloci-1;j>=0;j--) fprintf(fo," %3d",d[j]+1);
    fprintf(fo," %d\n",i+1);
    digitm(loci1,d,0);
  }
  }
}

void counting(pat *table, long int io)
{
long int i,j,k,l,k1,k2,nhet2;
short la[mxloc],ua[mxloc],hetid[mxloc],nhet,d[mxloc+1];
double s,ej;

l=0;
for(j=0;j<nloci;j++) hetid[j]=0;
for(j=0;j<nloci;j++)
{
  if(table[io].hd.d[j][0]!=table[io].hd.d[j][1])
  {
    hetid[l]=j;
    ++l;
  }
}
nhet=l;
if(nhet>0)
{
  nhet2=(int)pow(2,nhet-1);
  s=0;
  for(j=0;j<=nloci;j++) d[j]=0;
  for(j=0;j<nhet2;j++)
  {
    for(k=nloci-1;k>=0;k--)
    {
      la[k]=table[io].hd.d[k][0];
      ua[k]=table[io].hd.d[k][1];
    }
    for(l=0;l<nhet;l++) if(d[l]==1)
    {
      k=la[hetid[l]];
      la[hetid[l]]=ua[hetid[l]];
      ua[hetid[l]]=k;
    }
    k1=linenum(loci,la)-1;
    k2=linenum(loci,ua)-1;
    s+=2*h[k1]*h[k2];
    digit2(1,d,0);
  }
  for(j=0;j<=nloci;j++) d[j]=0;
  for(j=0;j<nhet2;j++)
  {
    for(k=0;k<nloci;k++)
    {
      la[k]=table[io].hd.d[k][0];
      ua[k]=table[io].hd.d[k][1];
    }
    for(l=0;l<nhet;l++) if(d[l]==1)
    {
      k=la[hetid[l]];
      la[hetid[l]]=ua[hetid[l]];
      ua[hetid[l]]=k;
    }
    k1=linenum(loci,la)-1;
    k2=linenum(loci,ua)-1;
    ej=2*h[k1]*h[k2]/s;
    c[k1]+=ej*table[io].count;
    c[k2]+=ej*table[io].count;
    digit2(1,d,0);
  }
}
else
{
  for(j=0;j<nloci;j++) la[j]=table[io].hd.d[j][0];
  k=linenum(loci,la)-1;
  c[k]+=2*table[io].count;
}
}

void geth(pat *table)
/*haplotype frequencies with missing data (conditioned)*/
{
long int io,i,j,k,l,k1,k2,cycle,ncycle,nhet2;
short loci1[mxloc],l0[mxloc],l1[mxloc],lk[mxloc],lq[mxloc];
short la[mxloc],ua[mxloc],hetid[mxloc],nhet,d[mxloc+1];
double s,ss,tc,ej;

if(!handle_missing)
{
  for(i=0;i<hapall;i++) c[i]=0;
  for(i=0;i<obscom;i++)
  {
     if((!xdata)||(xdata&&table[i].sex==FEMALE)) counting(table,i);
     else
     {
       k=linenum(loci,table[i].hd.h)-1;
       c[k]+=table[i].count;
     }
  }
  for(i=0;i<hapall;i++) h[i]=c[i]/(2*d_sample+h_sample);
}
else
{
  for(i=0;i<hapall;i++) c[i]=hm[i]=0;
  for(io=0;io<obscom;io++)
  {
    if(xdata&&table[io].sex==MALE)
    {
      j=0;
      cycle=1;
      for(i=0;i<nloci;i++)
      {
         idm[i]=0;
         locim[i]=0;
         k=table[io].hd.h[i];
         if(k<1||k>loci[i])
         {
           idm[i]=1;
           locim[j]=loci[i];
           cycle*=loci[i];
           j++;
         }
      }
      ncycle=cycle;
      nlocim=j;
      if(nlocim==0||!handle_missing)
      {
        k=linenum(loci,table[io].hd.h)-1;
        c[k]+=table[io].count;
        continue;
      }
      tc=0;
      for(i=0;i<=nloci;i++) d[i]=0;
      for(i=0;i<nloci;i++) lk[i]=loci1[i]=0;
      for(i=0;i<nlocim;i++) loci1[i]=locim[nlocim-i-1]-1;
      do
      {
        for(i=nlocim-1;i>=0;i--) lk[nlocim-i-1]=d[i]+1;
        j=0;
        for(i=0;i<nloci;i++)
        {
          lq[i]=table[io].hd.h[i];
          if(idm[i]==1)
          {
            lq[i]=lk[j];
            j++;
          }
        }
        k=linenum(loci,lq)-1;
        tc+=h[k];
        digitm(loci1,d,0);
      } while(--cycle>0);
      for(i=0;i<=nloci;i++) d[i]=0;
      for(i=0;i<nloci;i++) lk[i]=loci1[i]=0;
      for(i=0;i<nlocim;i++) loci1[i]=locim[nlocim-i-1]-1;
      do
      {
        for(i=nlocim-1;i>=0;i--) lk[nlocim-i-1]=d[i]+1;
        j=0;
        for(i=0;i<nloci;i++)
        {
          lq[i]=table[io].hd.h[i];
          if(idm[i]==1)
          {
            lq[i]=lk[j];
            j++;
          }
        }
        if(tc==0) s=0;
        else s=table[io].count/tc;
        k=linenum(loci,lq)-1;
        hm[k]+=h[k]*s;
        digitm(loci1,d,0);
      } while(--ncycle>0);
      continue;
    }
    l=0;
    cycle=1;
    for(j=0;j<nloci;j++)
    {
      idm[j]=0;
      locim[j]=0;
      lk[j]=lq[j]=0;
      k1=table[io].hd.d[j][0];
      k2=table[io].hd.d[j][1];
      if((k1<1||k1>loci[j])||(k2<1||k2>loci[j]))
      {
        cycle*=loci[j]*(loci[j]+1)/2;
        idm[j]=1;
        locim[l]=loci[j];
        lk[l]=lq[l]=1;
        ++l;
      }
    }
    ncycle=cycle;
    nlocim=l;
    if(nlocim==0)
    {
      counting(table,io);
      continue;
    }
    /*obtain total weight*/
    tc=0;
    do
    {
      l=0;
      for(j=0;j<nloci;j++)
      {
        l0[j]=table[io].hd.d[j][0];
        l1[j]=table[io].hd.d[j][1];
        if(idm[j]==1)
        {
          l0[j]=lk[l];
          l1[j]=lq[l];
          ++l;
        }
      }
      tc+=phasep(l0,l1);
      if(nlocim==1)
      {
        ++lk[0];
        if(lk[0]>lq[0])
        {
          lk[0]=1;
          ++lq[0];
        }
      }
      else if(nlocim>1)
      {
        ++lk[nlocim-1];
        for(i=nlocim-1;i>0;i--) if(lk[i]>lq[i])
        {
          ++lq[i];
          lk[i]=1;
          if(lq[i]>locim[i])
          {
            lq[i]=1;
            ++lk[i-1];
            if((i==1)&&(lk[i-1]>lq[i-1]))
            {
              ++lq[i-1];
              lk[i-1]=1;
            }
          }
        }
      }
    } while(--cycle>0);
    l=0;
    for(j=0;j<nloci;j++)
    {
      lk[j]=lq[j]=0;
      if(idm[j]==1)
      {
        lk[l]=lq[l]=1;
        ++l;
      }
    }
    /*update haplotype counts*/
    do
    {
      k1=k2=0;
      for(j=0;j<nloci;j++)
      {
        l0[j]=table[io].hd.d[j][0];
        l1[j]=table[io].hd.d[j][1];
        if(idm[j]==1)
        {
          l0[j]=lk[k1];
          l1[j]=lq[k1];
          ++k1;
        }
        hetid[j]=0;
        if(l0[j]!=l1[j])
        {
          hetid[k2]=j;
          ++k2;
        }
      }
      nhet=k2;
      if(tc==0) s=0;
      else s=table[io].count/tc;
      if(nhet>0)
      {
        nhet2=(int)pow(2,nhet-1);
        for(j=0;j<=nloci;j++) d[j]=0;
        for(j=0;j<nhet2;j++)
        {
          for(k=0;k<nloci;k++)
          {
            la[k]=l0[k];
            ua[k]=l1[k];
          }
          for(l=0;l<nhet;l++) if(d[l]==1)
          {
            k=la[hetid[l]];
            la[hetid[l]]=ua[hetid[l]];
            ua[hetid[l]]=k;
          }
          k1=linenum(loci,la)-1;
          k2=linenum(loci,ua)-1;
          ej=2*h[k1]*h[k2];
          hm[k1]+=ej*s;
          hm[k2]+=ej*s;
          digit2(1,d,0);
        }
      }
      else
      {
        k=linenum(loci,l0)-1;
        hm[k]+=2*h[k]*h[k]*s;
      }
      if(nlocim==1)
      {
        ++lk[0];
        if(lk[0]>lq[0])
        {
          lk[0]=1;
          ++lq[0];
        }
      }
      else if(nlocim>1)
      {
        ++lk[nlocim-1];
        for(i=nlocim-1;i>0;i--) if(lk[i]>lq[i])
        {
          ++lq[i];
          lk[i]=1;
          if(lq[i]>locim[i])
          {
            lq[i]=1;
            ++lk[i-1];
            if((i==1)&&(lk[i-1]>lq[i-1]))
            {
              ++lq[i-1];
              lk[i-1]=1;
            }
          }
        }
      }
    } while(--ncycle>0);
  }
  tc=0;
  for(i=0;i<hapall;i++) tc+=hm[i];
  for(i=0;i<hapall;i++) h[i]=(c[i]+hm[i])/(d_sample*2+tc+h_sample);
}
/* if calling geth() for observed data first
tc*=0.5;
for(i=0;i<hapall;i++) h[i]=(h[i]*d_sample+hm[i]/2)/(d_sample+tc);
*/
}

double phasep(short *l0,short *l1)
/*phase probability*/
{
long int k,k1,k2,nhet2;
short i,j,l,nhet,hetid[mxloc],d[mxloc+1],la[mxloc],ua[mxloc];
double s;

for(j=0;j<nloci;j++) hetid[j]=0;
l=0;
for(i=0;i<nloci;i++) if(l0[i]!=l1[i])
{
  hetid[l]=i;
  ++l;
}
nhet=l;
if(nhet>0)
{
  nhet2=(int)pow(2,nhet-1);
  s=0;
  for(j=0;j<=nloci;j++) d[j]=0;
  for(j=0;j<nhet2;j++)
  {
    for(k=0;k<nloci;k++)
    {
      la[k]=l0[k];
      ua[k]=l1[k];
    }
    for(l=0;l<nhet;l++) if(d[l]==1)
    {
      k=la[hetid[l]];
      la[hetid[l]]=ua[hetid[l]];
      ua[hetid[l]]=k;
    }
    k1=linenum(loci,la)-1;
    k2=linenum(loci,ua)-1;
    s+=2*h[k1]*h[k2];
    digit2(1,d,0);
  }
}
else
{
  k=linenum(loci,l0)-1;
  s=h[k]*h[k];
}

return s;
}

void genp(pat *table,long int io)
/*get genotype probability */
{
long int i,j,k,l,k1,k2;
long int cycle,nhet2;
short loci1[mxloc],l0[mxloc],l1[mxloc],lk[mxloc],lq[mxloc];
short la[mxloc],ua[mxloc],hetid[mxloc],nhet,d[mxloc+1];
double tc,s;

if(xdata&&table[io].sex==MALE)
{
  j=0;
  cycle=1;
  for(i=0;i<nloci;i++)
  {
     idm[i]=0;
     locim[i]=0;
     k=table[io].hd.h[i];
     if(k<1||k>loci[i])
     {
       idm[i]=1;
       locim[j]=loci[i];
       cycle*=loci[i];
       j++;
     }
  }
  nlocim=j;
  if(nlocim==0)
  {
    k=linenum(loci,table[io].hd.h)-1;
    tc=h[k];
  }
  else
  {
    tc=0;
    for(i=0;i<nloci;i++) lk[i]=loci1[i]=0;
    for(i=0;i<=nloci;i++) d[i]=0;
    for(i=0;i<nlocim;i++) loci1[i]=locim[nlocim-i-1]-1;
    do
    {
      for(i=nlocim-1;i>=0;i--) lk[nlocim-i-1]=d[i]+1;
      j=0;
      for(i=0;i<nloci;i++)
      {
        lq[i]=table[io].hd.h[i];
        if(idm[i]==1)
        {
          lq[i]=lk[j];
          j++;
        }
      }
      k=linenum(loci,lq)-1;
      tc+=h[k];
      digitm(loci1,d,0);
    } while(--cycle>0);
  }
  table[io].prob=tc;
  return;
}
l=0;
cycle=1;
for(j=0;j<nloci;j++)
{
  idm[j]=0;
  locim[j]=0;
  lk[j]=lq[j]=0;
  k1=table[io].hd.d[j][0];
  k2=table[io].hd.d[j][1];
  if((k1<1||k1>loci[j])||(k2<1||k2>loci[j]))
  {
    cycle*=loci[j]*(loci[j]+1)/2;
    idm[j]=1;
    locim[l]=loci[j];
    lk[l]=lq[l]=1;
    ++l;
  }
}
nlocim=l;
if(nlocim>0)
{
  s=0;
  do
  {
    l=0;
    for(j=0;j<nloci;j++)
    {
      l0[j]=table[io].hd.d[j][0];
      l1[j]=table[io].hd.d[j][1];
      if(idm[j]==1)
      {
        l0[j]=lk[l];
        l1[j]=lq[l];
        ++l;
      }
    }
    s+=phasep(l0,l1);
    if(nlocim==1)
    {
      ++lk[0];
      if(lk[0]>lq[0])
      {
        lk[0]=1;
        ++lq[0];
      }
    }
    else if(nlocim>1)
    {
      ++lk[nlocim-1];
      for(i=nlocim-1;i>0;i--) if(lk[i]>lq[i])
      {
        ++lq[i];
        lk[i]=1;
        if(lq[i]>locim[i])
        {
          lq[i]=1;
          ++lk[i-1];
          if((i==1)&&(lk[i-1]>lq[i-1]))
          {
            ++lq[i-1];
            lk[i-1]=1;
          }
        }
      }
    }
  } while(--cycle>0);
}
else
{
  l=0;
  for(j=0;j<nloci;j++) hetid[j]=0;
  for(j=0;j<nloci;j++)
  {
    if(table[io].hd.d[j][0]!=table[io].hd.d[j][1])
    {
      hetid[l]=j;
      ++l;
    }
  }
  nhet=l;
  if(nhet>0)
  {
    nhet2=(int)pow(2,nhet-1);
    s=0;
    for(j=0;j<=nloci;j++) d[j]=0;
    for(j=0;j<nhet2;j++)
    {
      for(k=nloci-1;k>=0;k--)
      {
        la[k]=table[io].hd.d[k][0];
        ua[k]=table[io].hd.d[k][1];
      }
      for(l=0;l<nhet;l++) if(d[l]==1)
      {
        k=la[hetid[l]];
        la[hetid[l]]=ua[hetid[l]];
        ua[hetid[l]]=k;
      }
      k1=linenum(loci,la)-1;
      k2=linenum(loci,ua)-1;
      s+=2*h[k1]*h[k2];
      digit2(1,d,0);
    }
  }
  else
  {
    for(j=0;j<nloci;j++) la[j]=table[io].hd.d[j][0];
    k=linenum(loci,la)-1;
    s=h[k]*h[k];
  }
}
table[io].prob=s;
}

double ll(pat *table)
/*log-likelihood*/
{
int i,j,l;
long int k,k1,k2;
double t,lnl,xlnl;

lnl=xlnl=0;
for(i=0;i<obscom;i++)
{
  if(xdata&&table[i].sex==MALE)
  {
    l=0;
    for(j=0;j<nloci;j++)
    {
      k=table[i].hd.h[j];
      if(k<1||k>loci[j]) ++l;
    }
    if(l>0&&!handle_missing) continue;
    t=table[i].count;
    if(t!=0&&table[i].prob>0) xlnl+=t*log(table[i].prob);
    continue;
  }
  l=0;
  for(j=0;j<nloci;j++)
  {
    k1=table[i].hd.d[j][0];
    k2=table[i].hd.d[j][1];
    if((k1<1||k1>loci[j])||(k2<1||k2>loci[j])) ++l;
  }
  if(l>0&&!handle_missing) continue;
  t=table[i].count;
  if(t!=0&&table[i].prob>0) lnl+=t*log(table[i].prob);
}
return (lnl+xlnl);
}

void ch(pat *table)
/*trace observed haplotypes*/
{
long int i,j,l,a1,a2;
long int k,k1,k2,nhet2;
short k0,la[mxloc],ua[mxloc];
short hetid[mxloc],nhet;
short d[mxloc+1],loci1[mxloc];
long int cycle, ncycle;
short l0[mxloc],l1[mxloc],lk[mxloc],lq[mxloc];
double tc,ej;

fprintf(fo,"\nAssignment of haplotypes by posterior probabilities\n");
fprintf(fo,"(ID[chromosome] haplotype, probability[cutoff=%lf] id)\n\n",pl);
for(i=0;i<obscom;i++)
{
  if(xdata&&table[i].sex==MALE)
  {
    l=0;
    cycle=1;
    for(j=0;j<nloci;j++)
    {
       idm[j]=0;
       locim[j]=0;
       k=table[i].hd.h[j];
       if(k<1||k>loci[j])
       {
         idm[j]=1;
         locim[l]=loci[j];
         cycle*=loci[j];
         l++;
       }
    }
    nlocim=l;
    ncycle=cycle;
    if(nlocim==0)
    {
      k=linenum(loci,table[i].hd.h)-1;
      fprintf(fo,"%5d [1]",table[i].gid);
      for(j=0;j<nloci;j++) fprintf(fo," %2d",table[i].hd.h[j]);
      fprintf(fo," %lf %ld\n",1.0,k+1);
      if(!hrt) hrt=hitree(hrt,k+1,table[i].hd.h);
      htr_table[i][k]=1.0;
    }
    else
    {
      tc=0;
      for(j=0;j<nloci;j++) lk[j]=loci1[j]=0;
      for(j=0;j<=nloci;j++) d[j]=0;
      for(j=0;j<nlocim;j++) loci1[j]=locim[nlocim-j-1]-1;
      do
      {
        for(j=nlocim-1;j>=0;j--) lk[nlocim-j-1]=d[j]+1;
        l=0;
        for(j=0;j<nloci;j++)
        {
          lq[j]=table[i].hd.h[j];
          if(idm[j]==1)
          {
            lq[j]=lk[l];
            l++;
          }
        }
        k=linenum(loci,lq)-1;
        tc+=h[k];
        digitm(loci1,d,0);
      } while(--cycle>0);
      for(j=0;j<=nloci;j++) d[j]=0;
      for(j=0;j<nloci;j++) lk[j]=loci1[j]=0;
      for(j=0;j<nlocim;j++) loci1[j]=locim[nlocim-j-1]-1;
      do
      {
        for(j=nlocim-1;j>=0;j--) lk[nlocim-j-1]=d[j]+1;
        l=0;
        for(j=0;j<nloci;j++)
        {
          lq[j]=table[i].hd.h[j];
          if(idm[j]==1)
          {
            lq[j]=lk[l];
            l++;
          }
        }
        k=linenum(loci,lq)-1;
        if(tc==0) ej=0;
        else ej=h[k]/tc;
        if(ej>pl)
        {
          fprintf(fo,"%5d [1]",table[i].gid);
          for(a1=0;a1<nloci;a1++) fprintf(fo," %2d",lq[a1]);
          fprintf(fo," %lf %ld\n",ej,k+1);
          htr_table[i][k+1]+=ej;
        }
        digitm(loci1,d,0);
      } while(--ncycle>0);
    }
    continue;
  }
  l=0;
  cycle=1;
  for(j=0;j<nloci;j++)
  {
    idm[j]=0;
    locim[j]=0;
    lk[j]=lq[j]=0;
    k1=table[i].hd.d[j][0];
    k2=table[i].hd.d[j][1];
    if((k1<1||k1>loci[j])||(k2<1||k2>loci[j]))
    {
      cycle*=loci[j]*(loci[j]+1)/2;
      idm[j]=1;
      locim[l]=loci[j];
      lk[l]=lq[l]=1;
      ++l;
    }
  }
  nlocim=l;
  tc=table[i].prob;
  if(nlocim>0) do
  {
    k1=k2=0;
    for(j=0;j<nloci;j++)
    {
      l0[j]=table[i].hd.d[j][0];
      l1[j]=table[i].hd.d[j][1];
      if(idm[j]==1)
      {
        l0[j]=lk[k1];
        l1[j]=lq[k1];
        ++k1;
      }
      hetid[j]=0;
      if(l0[j]!=l1[j])
      {
        hetid[k2]=j;
        ++k2;
      }
    }
    nhet=k2;
    if(nhet>0)
    {
      nhet2=(int)pow(2,nhet-1);
      for(j=0;j<=nloci;j++) d[j]=0;
      for(j=0;j<nhet2;j++)
      {
        for(k=0;k<nloci;k++)
        {
          la[k]=l0[k];
          ua[k]=l1[k];
        }
        for(l=0;l<nhet;l++) if(d[l]==1)
        {
          k=la[hetid[l]];
          la[hetid[l]]=ua[hetid[l]];
          ua[hetid[l]]=k;
        }
        k1=linenum(loci,la)-1;
        k2=linenum(loci,ua)-1;
        if(tc==0) ej=0;
        else ej=2*h[k1]*h[k2]/tc;
        if(ej>pl)
        {
          fprintf(fo,"%5d [1]",table[i].gid);
          for(a1=0;a1<nloci;a1++) fprintf(fo," %2d",la[a1]);
          fprintf(fo," %lf %ld\n",ej,k1+1);
          fprintf(fo,"%5d [2]",table[i].gid);
          for(a2=0;a2<nloci;a2++) fprintf(fo," %2d",ua[a2]);
          fprintf(fo," %lf %ld\n",ej,k2+1);
          htr_table[i][k1+1]+=ej;
          htr_table[i][k2+1]+=ej;
        }
        digit2(1,d,0);
      }
    }
    else
    {
      k=linenum(loci,l0)-1;
      if(tc==0) ej=0;
      else ej=h[k]*h[k]/tc;
      if(ej>pl)
      {
        fprintf(fo,"%5d [1]",table[i].gid);
        for(a1=0;a1<nloci;a1++) fprintf(fo," %2d",l0[a1]);
        fprintf(fo," %lf %ld\n",ej,k+1);
        fprintf(fo,"%5d [2]",table[i].gid);
        for(a2=0;a2<nloci;a2++) fprintf(fo," %2d",l0[a2]);
        fprintf(fo," %lf %ld\n",ej,k+1);
        htr_table[i][k+1]+=2.0*ej;
      }
    }
    if(nlocim==1)
    {
      ++lk[0];
      if(lk[0]>lq[0])
      {
        lk[0]=1;
        ++lq[0];
      }
    }
    else if(nlocim>1)
    {
      ++lk[nlocim-1];
      for(k=nlocim-1;k>0;k--) if(lk[k]>lq[k])
      {
        ++lq[k];
        lk[k]=1;
        if(lq[k]>locim[k])
        {
          lq[k]=1;
          ++lk[k-1];
          if((k==1)&&(lk[k-1]>lq[k-1]))
          {
            ++lq[k-1];
            lk[k-1]=1;
          }
        }
      }
    }
  } while(--cycle>0);
  else
  {
    l=0;
    for(j=0;j<nloci;j++)
    {
      if(table[i].hd.d[j][0]!=table[i].hd.d[j][1])
      {
        hetid[l]=j;
        ++l;
      }
    }
    nhet=l;
    if(nhet>0)
    {
      nhet2=(int)pow(2,nhet-1);
      for(j=0;j<=nloci;j++) d[j]=0;
      for(j=0;j<nhet2;j++)
      {
        for(k=nloci-1;k>=0;k--)
        {
          la[k]=table[i].hd.d[k][0];
          ua[k]=table[i].hd.d[k][1];
        }
        for(l=0;l<nhet;l++) if(d[l]==1)
        {
          k0=la[hetid[l]];
          la[hetid[l]]=ua[hetid[l]];
          ua[hetid[l]]=k0;
        }
        k1=linenums(loci,la)-1;
        k2=linenums(loci,ua)-1;
        if(tc==0) ej=0;
        else ej=2.0*h[k1]*h[k2]/tc;
        if(ej>pl)
        {
          fprintf(fo,"%5d [1]",table[i].gid);
          for(a1=0;a1<nloci;a1++) fprintf(fo," %2d",la[a1]);
          fprintf(fo," %lf %ld\n",ej,k1+1);
          fprintf(fo,"%5d [2]",table[i].gid);
          for(a2=0;a2<nloci;a2++) fprintf(fo," %2d",ua[a2]);
          fprintf(fo," %lf %ld\n",ej,k2+1);
          htr_table[i][k1+1]+=ej;
          htr_table[i][k2+1]+=ej;
        }
        if(!hrt) hrt=hitree(hrt,k1+1,la);
        else hitree(hrt,k1+1,la);
        if(!hrt) hrt=hitree(hrt,k2+1,ua);
        else hitree(hrt,k2+1,ua);
        digit2(1,d,0);
      }
    } else {
      for(j=0;j<nloci;j++) la[j]=table[i].hd.d[j][0];
      k=linenums(loci,la)-1;
      fprintf(fo,"%5d [1]",table[i].gid);
      for(a1=0;a1<nloci;a1++) fprintf(fo," %2d",la[a1]);
      fprintf(fo," %lf %ld\n",1.0,k+1);
      fprintf(fo,"%5d [2]",table[i].gid);
      for(a2=0;a2<nloci;a2++) fprintf(fo," %2d",la[a2]);
      fprintf(fo," %lf %ld\n",1.0,k+1);
      htr_table[i][k+1]=2.0;
      if(!hrt) hrt=hitree(hrt,k+1,la);
      else hitree(hrt,k+1,la);
      if(!hrt) hrt=hitree(hrt,k+1,ua);
      else hitree(hrt,k+1,ua);
    }
  }
}
}

hnode *hitree(hnode *r,long int id,short l[mxloc])
/*insert*/
{
int i;
if (r==NULL)
{
  r=malloc(sizeof(hnode));
  r->left=r->right=NULL;
  r->id=id;
  r->n=0;
  for(i=0;i<nloci;i++) r->l[i]=l[i];
}
else
if (id<r->id) r->left=hitree(r->left,id,l);
else if (id>r->id) r->right=hitree(r->right,id,l);
else r->n++;
return r;
}

void hptree(hnode *r,long int *l)
/*inorder*/
{
int i;

if (!r) return;
hptree(r->left,l);
++*l;
fprintf(fo," %.6lf [%.12lf]",h0[r->id-1],h0[r->id-1]);
fprintf(fo," %.6lf [%.12lf]",h[r->id-1],h[r->id-1]);
for(i=0;i<nloci;i++) fprintf(fo," %2hd",r->l[i]);
fprintf(fo," %d\n",r->id);
hptree(r->right,l);
}

void hrtree(hnode *t)
/*postorder*/
{
if (!t) return;
hrtree(t->left);
hrtree(t->right);
free(t);
}

void digit2(int radix, short d[], int i)
/*binary routine*/
{
if(d[i]<radix)
{
  ++d[i];
  return;
}
else
{
  d[i]=0;
  ++d[i+1];
  if(d[i+1]<=radix) return;
}
digit2(radix,d,i+1);
}

void digitm(short radix[], short d[], int i)
/*mixed-radix routine*/
{
if(d[i]<radix[i])
{
  ++d[i];
  return;
}
else
{
  d[i]=0;
  ++d[i+1];
  if(d[i+1]<=radix[i+1]) return;
}
digitm(radix,d,i+1);
}

int linenum(int *loci, short *ai)
/*haplotype array index*/
{
int loc,j;
loc=0;
for(j=1;j<=nloci;j++)
  if(j==nloci) loc+=ai[j-1];
  else loc=(loc+ai[j-1]-1)*loci[j];
return(loc);
}

int linenums(int *loci, short *ai)
/*array index of existing haplotypes*/
{
int loc,j;
loc=0;
for(j=1;j<=nloci;j++)
  if(j==nloci) loc+=ai[j-1];
  else loc=(loc+ai[j-1]-1)*loci[j];
return(loc);
}

void hilo(short *a, short *b)
{
short temp;
temp = *a;*a = *b;*b = temp;
}

/*code for dynamic allocation*/

/*
 * pass various tests after s/t initialized
 * JH Zhao 16-17/5/2000, 22/03/2001
 */

void *xmalloc(long len)
{
void *mem;

mem=malloc(len);
if(mem) return mem;
fprintf(stderr,"Sorry, but I cannot allocate memory\n");
exit(-1);

}

patlink pilink(pat new_pat,patlink t)
/*insert*/
{
/*
 * list from t to s, head h is maintained
 */
  patlink s,p1,p2,h;

  s=(patp *)xmalloc(sizeof(patp));
  s->anyline=new_pat;
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

void pclink(patlink h,pat *table)
/*copy*/
{
  patlink t;
  int i;

  i=0;
  if(h==NULL) {
    fprintf(stderr,"Empty list.\n");
  } else {
    while(h!=NULL)
    {
      table[i]=h->anyline;
      t=h;
      h=h->next;
      i++;
      free(t);
    }
  }
}

void qsorts(long int from,long int to)
/************************************************/
/* David Brunskill and John Turner (1996)       */
/* Understanding Algorithms and Data Structures */
/* McGraw-Hill                                  */
/************************************************/
{
long int i,pivot;

if(to>from)
{
  pivot=from;
  for(i=from+1;i<=to;++i)
  {
    haplotype_t=haplotypes[i];
    if(haplotype_t.ld>haplotypes[pivot].ld)
    {
      haplotypes[i]=haplotypes[pivot+1];
      haplotypes[pivot+1]=haplotypes[pivot];
      haplotypes[pivot]=haplotype_t;
      pivot++;
    }
  }
  qsorts(from,pivot-1);
  qsorts(pivot+1,to);
}
}

/* History:                                                                */
/*                                                                         */
/* -2001-                                                                  */
/*                                                                         */
/* 08-01 Start implementation                                              */
/* 09-01 Work for 2x2 case                                                 */
/* 10-01 Add linenum, digit2, etc                                          */
/* 11-01 Fix allele and initial haplotype frequencies                      */
/* 12-01 Fix EM and likelihood, add ALDH example                           */
/* 13-01 Integrate three tests into one routine, use eps, and fix digitm   */
/* 14-03 Work through missing data problem                                 */
/* 16-03 Check out a number of pecular scenarios                           */
/* 19-03 Check out c[i],hm[i],ss in gethm                                  */
/* 20-03 Check out error in geths and use new algorithms                   */
/* 21-03 Tune new algorithms in geths, change program name                 */
/* 22-03 Use dynamic allocation, add xmalloc, remove internal data         */
/* 12-06 Add I/O filenames to output, add lln()                            */
/* 14-06 Check out k1,k2 (-1) in gethm() and llm()                         */
/* 15-06 Check out dm, dm1 in geths()                                      */
/* 05-07 change default as missing value option and debug dynamics         */
/* 06-07 change geth(), ll(), add genc(), gethc(), maxit, remove lln()     */
/* 07-07 extensive check on gethc()                                        */
/* 08-07 check out sample in getp(), ncycle in gethc()                     */
/* 20-09 Resume work                                                       */
/* 24-09 delete linenumm(),gethm(),geths(),genp(),llm(), add lln()         */
/* 25-09 rearrange statement in lln() to obtain correct lln()              */
/* 01-11 absorb geth() into gethc() and add counting()                     */
/*                                                                         */
/* -2002-                                                                  */
/*                                                                         */
/* 19-04 starting code optimisation                                        */
/* 08-07 add posterior probability, adjusting for int/short                */
/* 09-07 rewrap work in 2002                                               */
/* 11-07 fix bug in posterior probabilty, add getprob()                    */
/* 12-07 add code to deal with all data missing                            */
/*       add Freeman-Tukey statistic                                       */
/*       fix bug in convergence by loglikelihood                           */
/*       merge geth() and gethc(), merge ll() and lln()                    */
/*       move digit2() output of if(ej>pl) in ch()                         */
/*       rename getprob to genp, genp to phasep                            */
/*       subroutines rearrangement                                         */
/* 15-07 add qsorts                                                        */
/* 16-07 tidy qsorts                                                       */
/* 20-07 add linenum for haplotype assignments                             */
/*       start implementing Haplotype Trend Regression                     */
/* 21-07 fix nhaplotypes->hapall                                           */
/*       output only nonzero haplotypes                                    */
/*       reset epsh=0.000001                                               */
/*                                                                         */
/* Compiling options:                                                      */
/*                                                                         */
/* -DDEBUG to view all possible haplotypes                                 */
/* -DNOHM  to turn off missing value handling mechanism                    */
/* -L1     to use L1 norm of haplotype frequencies as convergence criteria */

/* -2003- */
/* 14-08  X chromosome version in shape, start to incorporate              */

/* -2005- */
/* 13-04  change digit2, digitm to void as in R and add code to geth       */
/*        to allow for x/x genotypes in all cases                          */
/*        add documentation                                                */

typedef struct {
  int id,count;
  short l[mxloc];
  double p;
} haploid;
typedef struct {
  int id,count;
  short l[mxloc][2];
} diploid;

int xgeth(haploid *);
void xgenp(haploid *);
double xll(haploid *);

#define checkid() {\
    printf("[%5d]",men[io].id); \
    for(i=0;i<nlocim;i++) printf(" %2d",lk[i]);printf("\n");}
#define lli() {\
  l=0;\
  for(j=0;j<nloci;j++)\
  {\
    k=htable[i].l[j];\
    if(k<1||k>loci[j]) ++l;\
  }\
  if(l>0&&!handle_missing) continue;\
  t=htable[i].count;\
  if(t!=0&&htable[i].p>epsh) xlnl+=t*log(htable[i].p);}
#define xgethi() {\
for(io=0;io<nhaploid;io++)\
{\
  j=0;\
  cycle=1;\
  for(i=0;i<nloci;i++)\
  {\
     idm[i]=0;\
     locim[i]=0;\
     k=htable[io].l[i];\
     if(k<1||k>loci[i])\
     {\
       idm[i]=1;\
       locim[j]=loci[i];\
       cycle*=loci[i];\
       j++;\
     }\
  }\
  ncycle=cycle;\
  nlocim=j;\
  if(nlocim==0||!handle_missing)\
  {\
    k=linenum(loci,htable[io].l)-1;\
    c[k]+=htable[io].count;\
    continue;\
  }\
  tc=0;\
  for(i=0;i<=nloci;i++) d[i]=0;\
  for(i=0;i<nloci;i++) lk[i]=loci1[i]=0;\
  for(i=0;i<nlocim;i++) loci1[i]=locim[nlocim-i-1]-1;\
  do\
  {\
    for(i=nlocim-1;i>=0;i--) lk[nlocim-i-1]=d[i]+1;\
    j=0;\
    for(i=0;i<nloci;i++)\
    {\
      l[i]=htable[io].l[i];\
      if(idm[i]==1)\
      {\
        l[i]=lk[j];\
        j++;\
      }\
    }\
    k=linenum(loci,l)-1;\
    tc+=h[k];\
    digitm(loci1,d,0);\
  } while(--cycle>0);\
  for(i=0;i<=nloci;i++) d[i]=0;\
  for(i=0;i<nloci;i++) lk[i]=loci1[i]=0;\
  for(i=0;i<nlocim;i++) loci1[i]=locim[nlocim-i-1]-1;\
  do\
  {\
    for(i=nlocim-1;i>=0;i--) lk[nlocim-i-1]=d[i]+1;\
    j=0;\
    for(i=0;i<nloci;i++)\
    {\
      l[i]=htable[io].l[i];\
      if(idm[i]==1)\
      {\
        l[i]=lk[j];\
        j++;\
      }\
    }\
    if(tc==0) s=0;\
    else s=htable[io].count/tc;\
    k=linenum(loci,l)-1;\
    hm[k]+=h[k]*s;\
    digitm(loci1,d,0);\
  } while(--ncycle>0);\
}\
tc=0;\
for(i=0;i<hapall;i++) tc+=hm[i];\
for(i=0;i<hapall;i++) h[i]=(c[i]+hm[i])/(h_sample+tc);}

#define xgenpi() {\
  j=0;\
  cycle=1;\
  for(i=0;i<nloci;i++)\
  {\
     idm[i]=0;\
     locim[i]=0;\
     k=htable[io].l[i];\
     if(k<1||k>loci[i])\
     {\
       idm[i]=1;\
       locim[j]=loci[i];\
       cycle*=loci[i];\
       j++;\
     }\
  }\
  nlocim=j;\
  if(nlocim==0)\
  {\
    k=linenum(loci,htable[io].l)-1;\
    tc=h[k];\
  }\
  else\
  {\
    tc=0;\
    for(i=0;i<nloci;i++) lk[i]=loci1[i]=0;\
    for(i=0;i<=nloci;i++) d[i]=0;\
    for(i=0;i<nlocim;i++) loci1[i]=locim[nlocim-i-1]-1;\
    do\
    {\
      for(i=nlocim-1;i>=0;i--) lk[nlocim-i-1]=d[i]+1;\
      j=0;\
      for(i=0;i<nloci;i++)\
      {\
        l[i]=htable[io].l[i];\
        if(idm[i]==1)\
        {\
          l[i]=lk[j];\
          j++;\
        }\
      }\
      k=linenum(loci,l)-1;\
      tc+=h[k];\
      digitm(loci1,d,0);\
    } while(--cycle>0);\
  }\
  htable[io].p=tc;}

int xehm(char *xgcdat)
{
FILE *fi;
char record[301],eol[301];
int idn,count,sex,i,j,k,a1,a2,iter;
double s,ss,ft,lnl0,lnl1,lnls,n[mxloc];
short loci1[mxloc],d[mxloc+1],d1[mxloc+1];
haploid *htable;
diploid *dtable;

fi=fopen(xgcdat,"r");
if(!fi)
{
  perror("I cannot open the file");
  return 1;
}
fgets(record,300,fi);
nloci=0;
while(sscanf(record,"%d %[^\n]",&loci[nloci++],eol)>1) strcpy(record,eol);
nhaploid=ndiploid=0;
while(fgets(record,300,fi)
      &&sscanf(record,"%d %d %d %[^\n]",&idn,&count,&sex,eol)>3)
if(sex==0) nhaploid++;
else ndiploid++;
rewind(fi);
fgets(record,300,fi);
htable=(haploid*)malloc(nhaploid*sizeof(haploid));
dtable=(diploid*)malloc(ndiploid*sizeof(diploid));
if(!htable||!dtable)
{
  perror("can not allocate memory");
  return 1;
}
for(i=0;i<nloci;i++) n[i]=0;
for(i=0;i<mxloc;i++) for(j=0;j<loci[i];j++) pm[i][j]=0;
i=j=0;
h_sample=h_msample=0;
while(fgets(record,300,fi)
      &&sscanf(record,"%d %d %d %[^\n]",&idn,&count,&sex,eol)>3)
{
  if(debug) printf("%3d %d %2d",idn,count,sex);
  if(sex==0)
  {
    htable[i].id=idn;
    htable[i].count=1;
  }
  else
  {
    dtable[j].id=idn;
  }
  strcpy(record,eol);
  iter=0;
  for(k=0;k<nloci;k++)
  {
    if(sex==0)
    {
      sscanf(record,"%d %[^\n]",&a1,eol);
      htable[i].l[k]=a1;
      if(a1>=1&&a1<=loci[k])
      {
        ++iter;
        pm[k][a1-1]+=htable[i].count;
        n[k]+=htable[i].count;
      }
      if(debug) printf(" %5d",a1);
    }
    else
    {
      sscanf(record,"%d/%d %[^\n]",&a1,&a2,eol);
      dtable[j].l[k][0]=a1;
      dtable[j].l[k][1]=a2;
      if(debug) printf(" %2d/%2d",a1,a2);
    }
    strcpy(record,eol);
  }
  if(iter==nloci) h_sample+=htable[i].count;
  else h_msample+=htable[i].count;
  if(sex==0) i++;
  else j++;
  if(debug) printf("\n");
}
fclose(fi);
printf("\nNo of nonmissing individuals at each locus\n\n");
for(i=0;i<nloci;i++) printf("%3d %.0lf\n",i+1,n[i]);
printf("\n%d/%d individuals with F/P information\n",h_sample,h_msample);
for(i=0;i<nloci;i++) for(j=0;j<loci[i];j++) pm[i][j]/=n[i];
printf("\nAllele frequencies\n\nLocus/Frequecies\n\n");
for(i=0;i<nloci;i++)
{
  printf("%3d",i+1);
  for(j=0;j<loci[i];j++) printf(" %.4lf",pm[i][j]);
  printf("\n");
}
hapall=1;
for(i=0;i<nloci;i++) hapall*=loci[i];
inith(hapall);
for(i=0;i<hapall;i++) h[i]=h0[i]=0;
for(i=0;i<nloci;++i) loci1[i]=loci[i]-1;
for(i=0;i<=nloci;i++) d[i]=0;
for(i=0;i<hapall;i++)
{
  ss=1;
  for(j=0;j<nloci;j++) ss*=pm[j][d[j]];
  for(j=0;j<=nloci;j++) d1[j]=d[j]+1;
  j=linenum(loci,d1)-1;
  h[j]=h0[j]=ss;
  digitm(loci1,d,0);
}
if(debug) printf("\nHaplotype frequencies under linkage equilibrium\n\n");
for(i=0;i<=nloci;i++) d[i]=0;
for(i=0;i<nloci;i++) loci1[i]=loci[nloci-i-1]-1;
k=0;
for(i=0;i<hapall;i++) if(debug)
{
  if(h[i]>epsh)
  {
    printf("%6d %.6lf [%.12lf]",k+1,h[i],h[i]);
    for(j=nloci-1;j>=0;j--) printf(" %3d",d[j]+1);
    printf(" %d\n",i+1);
    ++k;
  }
  digitm(loci1,d,0);
}
if(debug) printf("\n%d haplotypes >%lf\n\n",k,epsh);
lnl0=xll(htable);
printf("\nlog-likelihood assuming linkage equilibrium = %.2lf\n\n",lnl0);
k=1;
iter=1;
do
{
  for(i=0;i<hapall;i++) hs[i]=h[i];
  xgeth(htable);
  if(!conv_ll)
  {
    s=0;
    for(i=0;i<hapall;i++) s+=fabs(hs[i]-h[i]);
  }
  else
  {
    printf("Iteration %3d, ",k++);
    lnl1=xll(htable);
    printf("log-likelihood=%.2lf\n",lnl1);
    s=lnl1-lnls;
    lnls=lnl1;
  }
} while(s>eps&&iter++<maxit);
if(!conv_ll) lnl1=xll(htable);
printf("\nHaplotype frequencies under linkage disequilibrium\n\n");
s=h_sample+h_msample;
for(i=0;i<=nloci;i++) d[i]=0;
for(i=0;i<nloci;i++) loci1[i]=loci[nloci-i-1]-1;
k=0;
for(i=0;i<hapall;i++)
{
  if(h[i]>epsh)
  {
    ft=sqrt(h[i]*s)+sqrt(h[i]*s+1)-sqrt(4*h0[i]*s+1);
    printf("%6d %.6lf %.12lf %.6lf %6.2lf",k+1,h[i],h[i],h0[i],ft);
    for(j=nloci-1;j>=0;j--) printf("%3d",d[j]+1);
    printf(" %d\n",i+1);
    k++;
  }
  digitm(loci1,d,0);
}
printf("\n%d haplotypes >%lf\n\n",k,epsh);
printf("log-likelihood assuming linkage disequilibrium = %.2lf\n",lnl1);

return 0;
}

int xgeth(haploid *htable)
/*handle missing genotype*/
{
int idm[mxloc],locim[mxloc],lk[mxloc];
int io,i,j,k,cycle,ncycle;
double tc,s;
short loci1[mxloc],l[mxloc],d[mxloc+1];

for(i=0;i<hapall;i++) c[i]=hm[i]=0;
xgethi();

return 0;
}

void xgenp(haploid *htable)
/*obtain genotype probability*/
{
int idm[mxloc],locim[mxloc],lk[mxloc];
int io,i,j,k,cycle;
double tc,s;
short loci1[mxloc],l[mxloc],d[mxloc+1];

for(io=0;io<nhaploid;io++) xgenpi();

}

double xll(haploid *htable)
/*log-likelihood including missing genotype*/
{
double xlnl=0,t;
int i,j,k,l;

xgenp(htable);
for(i=0;i<nhaploid;i++) lli();

return xlnl;
}
/*history:*/
/*7-8-3 start*/
/*10-8-3 use rewind to read data*/
/*11-8-3 rewrite loop in digitm*/
/*12-8-3 in shape with xgenp() and xll() working*/
/*14-8-3 remove xcounting, check for xgeth/xgenp, add argc/argv*/
