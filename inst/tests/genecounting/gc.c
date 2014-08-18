/*                                                           */
/* GENECOUNTING - haplotype analysis with missing genotypes  */
/* (C) Copyright JH Zhao 2001 Institute of Psychiatry, KCL   */
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
enum {NO, YES};
#ifdef NOHM
  #define handle_miss NO
#else
  #define handle_miss YES
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
double sample,msample;
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
  double count,prob,y;
  short l[mxloc][2];
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
int digit2(int,short*,int),digitm(short*,short*,int);
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

int main(int argc,char *argv[])
{
int gcmain(char*);
time_t t;

printf("\nGENECOUNTING, %s %.2f %s",(!handle_miss)?"ordinary":
"missing-value handling",version, " JH Zhao 03/01--07/05\n",version);
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
short la,ua;
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
for(i=0;i<mxloc;i++) tt.l[i][0]=tt.l[i][1]=0;
s=t=NULL;
i=0;
l=0;
while(fgets(line,200,fp)&&
sscanf(line,"%ld %lf %[^\n]",&tt.gid,&tt.count,rest)>2)
{
  k=0;
  strcpy(line,rest);
  for(j=0;j<nloci;j++)
  {
    if (sscanf(line,"%hd/%hd %[^\n]",&la,&ua,rest)<2&&
        sscanf(line,"%hd %hd %[^\n]",&la,&ua,rest)<2)
    {
      fprintf(stderr,"There isn't enough data at line %d\n",i+1);
      return 1;
    }
    strcpy(line,rest);
    *rest='\0';
    if (la>ua) hilo(&la,&ua);
    tt.l[j][0]=la;
    tt.l[j][1]=ua;
    if((la>loci[j]||ua>loci[j])||(la<=0||ua<=0)) ++k;
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
ehm(table);
free(hs);
free(c);
if(handle_miss) free(hm);
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
sn=2.0*(sample+msample);
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
#ifdef DEBUG
fprintf(fo,"\nFrequencies of possible haplotypes\n\n");
l=0;
hptree(hrt,&l);
fprintf(fo,"\nPossible haplotypes excluding missing loci = %ld\n",l);
#endif
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
if(handle_miss) hm=(double*)xmalloc(obshap*sizeof(double));
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
s=2.0*(sample+msample);
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
int i,j,k,l,u;
short d[mxloc+1],loci1[mxloc],d1[mxloc+1];
double s[mxloc],ss;

for(i=0;i<nloci;i++)
{
  s[i]=0;
  for(j=0;j<nalleles;j++) p[i][j]=0;
}
ss=msample=0;
for(i=0;i<obscom;i++)
{
  k=0;
  for(j=0;j<nloci;j++)
  {
    l=table[i].l[j][0]-1;
    u=table[i].l[j][1]-1;
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
  else msample+=table[i].count;
}
for(i=0;i<nloci;i++) fprintf(fo,"locus %3d # individuals = %.0lf\n",i+1,s[i]);
fprintf(fo,"\n%.lf/%.lf individuals with F/P information\n\n",ss,msample);
sample=ss;
fprintf(fo,"\nMarker allele frequencies:\n\nMarker\\Allele number\n\n");
fprintf(fo,"   ");
for(i=0;i<nalleles;i++) fprintf(fo,"%9d",i+1);
fprintf(fo,"\n\n");
for(i=0;i<nloci;i++)
{
  fprintf(fo,"%2d ",i+1);
  for(j=0;j<loci[i];j++)
  {
    p[i][j]/=2*s[i];
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
#ifdef DEBUG
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
#endif
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
  if(table[io].l[j][0]!=table[io].l[j][1])
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
      la[k]=table[io].l[k][0];
      ua[k]=table[io].l[k][1];
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
      la[k]=table[io].l[k][0];
      ua[k]=table[io].l[k][1];
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
  for(j=0;j<nloci;j++) la[j]=table[io].l[j][0];
  k=linenum(loci,la)-1;
  c[k]+=2*table[io].count;
}
}

void geth(pat *table)
/*haplotype frequencies with missing data (conditioned)*/
{
long int io,i,j,k,l,k1,k2,cycle,ncycle,nhet2;
short l0[mxloc],l1[mxloc],lk[mxloc],lq[mxloc];
short la[mxloc],ua[mxloc],hetid[mxloc],nhet,d[mxloc+1];
double s,ss,tc,ej;

if(!handle_miss)
{
  for(i=0;i<hapall;i++) c[i]=0;
  for(i=0;i<obscom;i++) counting(table,i);
  for(i=0;i<hapall;i++) h[i]=c[i]/2/sample;
}
else
{
  for(i=0;i<hapall;i++) c[i]=hm[i]=0;
  for(io=0;io<obscom;io++)
  {
    l=0;
    cycle=1;
    for(j=0;j<nloci;j++)
    {
      idm[j]=0;
      locim[j]=0;
      lk[j]=lq[j]=0;
      k1=table[io].l[j][0];
      k2=table[io].l[j][1];
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
        l0[j]=table[io].l[j][0];
        l1[j]=table[io].l[j][1];
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
        l0[j]=table[io].l[j][0];
        l1[j]=table[io].l[j][1];
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
  for(i=0;i<hapall;i++) h[i]=(c[i]+hm[i])/(sample*2+tc);
}
/* if calling geth() for observed data first
tc*=0.5;
for(i=0;i<hapall;i++) h[i]=(h[i]*sample+hm[i]/2)/(sample+tc);
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
long int cycle, nhet2;
short l0[mxloc],l1[mxloc],lk[mxloc],lq[mxloc];
short la[mxloc],ua[mxloc],hetid[mxloc],nhet,d[mxloc+1];
double s;

l=0;
cycle=1;
for(j=0;j<nloci;j++)
{
  idm[j]=0;
  locim[j]=0;
  lk[j]=lq[j]=0;
  k1=table[io].l[j][0];
  k2=table[io].l[j][1];
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
      l0[j]=table[io].l[j][0];
      l1[j]=table[io].l[j][1];
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
    if(table[io].l[j][0]!=table[io].l[j][1])
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
        la[k]=table[io].l[k][0];
        ua[k]=table[io].l[k][1];
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
    for(j=0;j<nloci;j++) la[j]=table[io].l[j][0];
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
long int k,k1,k2,nhet2;
short la[mxloc],ua[mxloc],d[mxloc+1],hetid[mxloc],nhet;
double t,s,lnl;

lnl=0;
for(i=0;i<obscom;i++)
{
  l=0;
  for(j=0;j<nloci;j++)
  {
    k1=table[i].l[j][0];
    k2=table[i].l[j][1];
    if((k1<1||k1>loci[j])||(k2<1||k2>loci[j])) ++l;
  }
  if(l>0&&!handle_miss) continue;
  t=table[i].count;
  if(t!=0&&table[i].prob>0) lnl+=t*log(table[i].prob);
}
return (lnl);
}

void ch(pat *table)
/*trace observed haplotypes*/
{
long int i,j,l,a1,a2;
long int k,k1,k2,nhet2;
short k0,la[mxloc],ua[mxloc];
short hetid[mxloc],nhet;
short d[mxloc+1];
long int cycle;
short l0[mxloc],l1[mxloc],lk[mxloc],lq[mxloc];
double tc,ej;

fprintf(fo,"\nAssignment of haplotypes by posterior probabilities\n");
fprintf(fo,"(ID[chromosome] haplotype, probability[cutoff=%lf] id)\n\n",pl);
for(i=0;i<obscom;i++)
{
  l=0;
  cycle=1;
  for(j=0;j<nloci;j++)
  {
    idm[j]=0;
    locim[j]=0;
    lk[j]=lq[j]=0;
    k1=table[i].l[j][0];
    k2=table[i].l[j][1];
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
      l0[j]=table[i].l[j][0];
      l1[j]=table[i].l[j][1];
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
      if(table[i].l[j][0]!=table[i].l[j][1])
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
          la[k]=table[i].l[k][0];
          ua[k]=table[i].l[k][1];
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
      for(j=0;j<nloci;j++) la[j]=table[i].l[j][0];
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

int digit2(int radix, short d[], int i)
/*binary routine*/
{
if(d[i]<radix)
{
  ++d[i];
  goto ok;
}
else
{
  d[i]=0;
  ++d[i+1];
  if(d[i+1]<=radix) goto ok;
}
digit2(radix,d,i+1);

ok:
return 0;
}

int digitm(short radix[], short d[], int i)
/*mixed-radix routine*/
{
if(d[i]<radix[i])
{
  ++d[i];
  goto ok;
}
else
{
  d[i]=0;
  ++d[i+1];
  if(d[i+1]<=radix[i+1]) goto ok;
}
digitm(radix,d,i+1);

ok:
return 0;
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
  int valids;

  valids=0;
  if(h==NULL) {
    fprintf(stderr,"Empty list.\n");
  } else {
    while(h!=NULL)
    {
      table[valids]=h->anyline;
      t=h;
      h=h->next;
      valids++;
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


