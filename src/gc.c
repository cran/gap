#include <R.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define version 1.30
#define mxloc 50
#define mxalleles 100

static double pl;
static int handlemissing;
static int nloci,nloci2,nalleles,loci[mxloc],obscom,nlocim,locim[mxloc+1];
static double sample,msample,pp[mxloc][mxalleles];
static long int hapall;
static double *c,*h,*hs,*h0,*hm;
static int idm[mxloc];

void *xmalloc(long),hilo(int*,int*),getp(int*,int *),geth(int*,int *),
     ch(int*,int*,double*,double*),counting(int*,int *, long int),
     genp(int*,long int,double*);
int digit2(int,int*,int),digitm(int*,int*,int);
int linenum(int*,int*),linenums(int*,int*);
double phasep(int*,int*),ll(int*,int*,double*);

void gc(
int *Rhandlemissing,
int *convll,
double *eps,
int *maxit,
double *Rpl,
double *precis,
int *gid,
int *Rnloci,
int *Rloci,
int *Robscom,
int *Rhapall,
int *genotype,
int *count,
int *hapid,
double *prob,
double *Rh0,
double *Rh1,
double *lnl0,
double *lnl1,
int *npusr,
int *npdat,
double *htrtable,
int *iter,
int *converge
)
{
long int i,j,k,io,it;
int loci1[mxloc];
double s,lnls;
/*
time_t t;
printf("\nGENECOUNTING, %s %.2f %s",(!handlemissing)?"ordinary":"missing-value handling",version, " JH Zhao 03/01--05/03\n");
time(&t);
printf("%s\n\n",ctime(&t));
printf("max. loci=%d, max. alleles=%d\n\n",mxloc,mxalleles);
*/
handlemissing=*Rhandlemissing;
pl=*Rpl;
nloci=*Rnloci;
for(i=0;i<nloci;i++) loci[i]=Rloci[i];
hapall=*Rhapall;
h0=Rh0;
h=Rh1;
obscom=*Robscom;
nloci2=2*nloci;
nalleles=1;
for(i=0;i<nloci;i++)
{
  if(loci[i]>nalleles) nalleles=loci[i];
}

hs=(double*)xmalloc(hapall*sizeof(double));
c=(double*)xmalloc(hapall*sizeof(double));
if(handlemissing) hm=(double*)xmalloc(hapall*sizeof(double));

getp(genotype,count);

for(i=0;i<2;i++) npusr[i]=npdat[i]=0;
k=1;
for(i=0;i<nloci;i++)
{
  npusr[0]+=loci[i]-1;
  k*=loci[i];
}
npusr[1]=k-1;
k=1;
for(i=0;i<nloci;i++)
{
  loci1[i]=loci[i];
  for(j=0;j<loci[i];j++) if(pp[i][j]<*precis) --loci1[i];
  npdat[0]+=loci1[i]-1;
  k*=loci1[i];
}
npdat[1]=k-1;

for(io=0;io<obscom;io++) genp(genotype,io,&prob[io]);

*lnl0=ll(genotype,count,prob);
lnls=*lnl0;
it=1;
do
{
  for(i=0;i<hapall;i++) hs[i]=h[i];
  geth(genotype,count);
  if(!*convll)
  {
    s=0;
    for(i=0;i<hapall;i++) s+=fabs(hs[i]-h[i]);
  }
  else
  {
/*
    printf("Iteration %3ld, ",it);
*/
    for(io=0;io<obscom;io++) genp(genotype,io,&prob[io]);
    *lnl1=ll(genotype,count,prob);
/*
    printf("log-likelihood=%.2f\n",*lnl1);
*/
    s=*lnl1-lnls;
    lnls=*lnl1;
  }
} while((s>*eps)&&(*maxit)>it++);
*iter=it;
if(!*convll)
{
  for(io=0;io<obscom;io++) genp(genotype,io,&prob[io]);
  *lnl1=ll(genotype,count,prob);
}
if(s>*eps) *converge=0;
else *converge=1;
free(hs);
free(c);
if(handlemissing) free(hm);

ch(gid,genotype,prob,htrtable);

}

void getp(int *genotype,int *count)
{
int i,j,k,l,u;
int d[mxloc+1],loci1[mxloc],d1[mxloc+1];
double s[mxloc],ss;
for(i=0;i<nloci;i++)
{
  s[i]=0;
  for(j=0;j<nalleles;j++) pp[i][j]=0;
}
ss=msample=0;
for(i=0;i<obscom;i++)
{
  k=0;
  for(j=0;j<nloci;j++)
  {
    l=genotype[i*nloci2+2*j]-1;
    u=genotype[i*nloci2+2*j+1]-1;
    if((l>=loci[j]||u>=loci[j])||(l<0||u<0))
    {
      ++k;
      continue;
    }
    pp[j][l]+=count[i];
    pp[j][u]+=count[i];
    s[j]+=count[i];
  }
  if(k==0) ss+=count[i];
  else msample+=count[i];
}
sample=ss;
for(i=0;i<nloci;i++)
{
  for(j=0;j<loci[i];j++)
  {
    pp[i][j]/=2*s[i];
  }
  loci1[i]=loci[i]-1;
}
for(i=0;i<=nloci;i++) d[i]=0;
for(i=0;i<hapall;i++)
{
  ss=1;
  for(j=0;j<nloci;j++) ss*=pp[j][d[j]];
  for(j=0;j<=nloci;j++) d1[j]=d[j]+1;
  l=linenum(loci,d1)-1;
  h[l]=h0[l]=ss;
  digitm(loci1,d,0);
}
}

void counting(int *genotype, int *count, long int io)
{
long int j,k,l,k1,k2;
int la[mxloc],ua[mxloc],hetid[mxloc],nhet,nhet2,d[mxloc+1];
double s,ej;

l=0;
for(j=0;j<nloci;j++) hetid[j]=0;
for(j=0;j<nloci;j++)
{
  if(genotype[io*nloci2+2*j]!=genotype[io*nloci2+2*j+1])
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
      la[k]=genotype[io*nloci2+2*k];
      ua[k]=genotype[io*nloci2+2*k+1];
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
      la[k]=genotype[io*nloci2+2*k];
      ua[k]=genotype[io*nloci2+2*k+1];
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
    c[k1]+=ej*count[io];
    c[k2]+=ej*count[io];
    digit2(1,d,0);
  }
}
else
{
  for(j=0;j<nloci;j++) la[j]=genotype[io*nloci2+2*j];
  k=linenum(loci,la)-1;
  c[k]+=2*count[io];
}
}

void geth(int *genotype,int *count)
{
long int io,i,j,k,l,k1,k2,cycle,ncycle;
int l0[mxloc],l1[mxloc],lk[mxloc],lq[mxloc];
int la[mxloc],ua[mxloc],hetid[mxloc],nhet,nhet2,d[mxloc+1];
double s,tc=0,ej;

if(!handlemissing)
{
  for(i=0;i<hapall;i++) c[i]=0;
  for(i=0;i<obscom;i++) counting(genotype,count,i);
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
      k1=genotype[io*nloci2+2*j];
      k2=genotype[io*nloci2+2*j+1];
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
      counting(genotype,count,io);
      continue;
    }
    tc=0;
    do
    {
      l=0;
      for(j=0;j<nloci;j++)
      {
        l0[j]=genotype[io*nloci2+2*j];
        l1[j]=genotype[io*nloci2+2*j+1];
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
    do
    {
      k1=k2=0;
      for(j=0;j<nloci;j++)
      {
        l0[j]=genotype[io*nloci2+2*j];
        l1[j]=genotype[io*nloci2+2*j+1];
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
      else s=count[io]/tc;
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
}

double phasep(int *l0,int *l1)
{
long int k,k1,k2;
int i,j,l,nhet,nhet2,hetid[mxloc],d[mxloc+1],la[mxloc],ua[mxloc];
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

void genp(int *genotype,long int io,double *prob)
{
long int i,j,k,l,k1,k2;
long int cycle;
int l0[mxloc],l1[mxloc],lk[mxloc],lq[mxloc];
int la[mxloc],ua[mxloc],hetid[mxloc],nhet,nhet2,d[mxloc+1];
double s;

l=0;
cycle=1;
for(j=0;j<nloci;j++)
{
  idm[j]=0;
  locim[j]=0;
  lk[j]=lq[j]=0;
  k1=genotype[io*nloci2+2*j];
  k2=genotype[io*nloci2+2*j+1];
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
      l0[j]=genotype[io*nloci2+2*j];
      l1[j]=genotype[io*nloci2+2*j+1];
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
    if(genotype[io*nloci2+2*j]!=genotype[io*nloci2+2*j+1])
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
        la[k]=genotype[io*nloci2+2*k];
        ua[k]=genotype[io*nloci2+2*k+1];
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
    for(j=0;j<nloci;j++) la[j]=genotype[io*nloci2+2*j];
    k=linenum(loci,la)-1;
    s=h[k]*h[k];
  }
}
*prob=s;
}

double ll(int *genotype,int *count, double *prob)
{
int i,j,l;
long int k1,k2;
double t,lnl;

lnl=0;
for(i=0;i<obscom;i++)
{
  l=0;
  for(j=0;j<nloci;j++)
  {
    k1=genotype[i*nloci2+2*j];
    k2=genotype[i*nloci2+2*j+1];
    if((k1<1||k1>loci[j])||(k2<1||k2>loci[j])) ++l;
  }
  if(l>0&&!handlemissing) continue;
  t=count[i];
  if(t!=0) lnl+=t*log(prob[i]);
}
return (lnl);
}

typedef struct hnode_type
{
  long int id;
  int n;
  struct hnode_type *left, *right;
  int l[mxloc];
} hnode;

hnode *hrt=0;
hnode *hitree(hnode*,long int,int*);
void hrtree(hnode*),hptree(hnode*,long int*);

void ch(int *gid, int *genotype,double *prob,double *htrtable)
{
FILE *fo;
long int i,j,l,a1,a2;
long int k,k1,k2;
int k0,la[mxloc],ua[mxloc];
int hetid[mxloc],nhet,nhet2;
int d[mxloc+1];
long int cycle;
int l0[mxloc],l1[mxloc],lk[mxloc],lq[mxloc];
double tc,ej;

fo=fopen("assign.dat","w");
if(!fo) 
{
  Rprintf("error openning assign.dat\n");
  return;
}
for(i=0;i<obscom;i++)
{
  l=0;
  cycle=1;
  for(j=0;j<nloci;j++)
  {
    idm[j]=0;
    locim[j]=0;
    lk[j]=lq[j]=0;
    k1=genotype[i*nloci2+2*j];
    k2=genotype[i*nloci2+2*j+1];
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
  tc=prob[i];
  if(nlocim>0) do
  {
    k1=k2=0;
    for(j=0;j<nloci;j++)
    {
      l0[j]=genotype[i*nloci2+2*j];
      l1[j]=genotype[i*nloci2+2*j+1];
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
        if(ej>=pl)
        {
          fprintf(fo,"%5d 1",gid[i]);
          for(a1=0;a1<nloci;a1++) fprintf(fo," %2d",la[a1]);
          fprintf(fo," %e %ld\n",ej,k1+1);
          fprintf(fo,"%5d 2",gid[i]);
          for(a2=0;a2<nloci;a2++) fprintf(fo," %2d",ua[a2]);
          fprintf(fo," %e %ld\n",ej,k2+1);
/*        htrtable[i*hapall+k1]+=ej;
          htrtable[i*hapall+k2]+=ej;
*/      }
        digit2(1,d,0);
      }
    }
    else
    {
      k=linenum(loci,l0)-1;
      if(tc==0) ej=0;
      else ej=h[k]*h[k]/tc;
      if(ej>=pl)
      {
        fprintf(fo,"%5d 1",gid[i]);
        for(a1=0;a1<nloci;a1++) fprintf(fo," %2d",l0[a1]);
        fprintf(fo," %e %ld\n",ej,k+1);
        fprintf(fo,"%5d 2",gid[i]);
        for(a2=0;a2<nloci;a2++) fprintf(fo," %2d",l0[a2]);
        fprintf(fo," %e %ld\n",ej,k+1);
/*      htrtable[i*hapall+k]+=2.0*ej;*/
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
      if(genotype[i*nloci2+2*j]!=genotype[i*nloci2+2*j+1])
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
          la[k]=genotype[i*nloci2+2*k];
          ua[k]=genotype[i*nloci2+2*k+1];
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
        if(ej>=pl)
        {
          fprintf(fo,"%5d 1",gid[i]);
          for(a1=0;a1<nloci;a1++) fprintf(fo," %2d",la[a1]);
          fprintf(fo," %e %ld\n",ej,k1+1);
          fprintf(fo,"%5d 2",gid[i]);
          for(a2=0;a2<nloci;a2++) fprintf(fo," %2d",ua[a2]);
          fprintf(fo," %e %ld\n",ej,k2+1);
/*        htrtable[i*hapall+k1]+=ej;
          htrtable[i*hapall+k2]+=ej;
*/      }
        if(!hrt) hrt=hitree(hrt,k1+1,la);
        else hitree(hrt,k1+1,la);
        if(!hrt) hrt=hitree(hrt,k2+1,ua);
        else hitree(hrt,k2+1,ua);
        digit2(1,d,0);
      }
    } else {
      for(j=0;j<nloci;j++) la[j]=genotype[i*nloci2+2*j];
      k=linenums(loci,la)-1;
      fprintf(fo,"%5d 1",gid[i]);
      for(a1=0;a1<nloci;a1++) fprintf(fo," %2d",la[a1]);
      fprintf(fo," %e %ld\n",1.0,k+1);
      fprintf(fo,"%5d 2",gid[i]);
      for(a2=0;a2<nloci;a2++) fprintf(fo," %2d",la[a2]);
      fprintf(fo," %e %ld\n",1.0,k+1);
/*    htrtable[i*hapall+k]=2.0;*/
      if(!hrt) hrt=hitree(hrt,k+1,la);
      else hitree(hrt,k+1,la);
      if(!hrt) hrt=hitree(hrt,k+1,ua);
      else hitree(hrt,k+1,ua);
    }
  }
}
fclose(fo);
}

hnode *hitree(hnode *r,long int id,int l[mxloc])
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
{

if (!r) return;
hptree(r->left,l);
++*l;
hptree(r->right,l);
}

void hrtree(hnode *t)
{
if (!t) return;
hrtree(t->left);
hrtree(t->right);
free(t);
}

int digit2(int radix, int d[], int i)
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

int digitm(int radix[], int d[], int i)
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

int linenum(int *loci, int *ai)
{
int loc,j;
loc=0;
for(j=1;j<=nloci;j++)
  if(j==nloci) loc+=ai[j-1];
  else loc=(loc+ai[j-1]-1)*loci[j];
return(loc);
}

int linenums(int *loci, int *ai)
{
int loc,j;
loc=0;
for(j=1;j<=nloci;j++)
  if(j==nloci) loc+=ai[j-1];
  else loc=(loc+ai[j-1]-1)*loci[j];
return(loc);
}

void hilo(int *a, int *b)
{
int temp;
temp = *a;*a = *b;*b = temp;
}

void *xmalloc(long len)
{
void *mem;

mem=malloc(len);
if(mem) return mem;
fprintf(stderr,"Sorry, but I cannot allocate memory\n");
exit(-1);

}

