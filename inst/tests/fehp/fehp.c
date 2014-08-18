#include "fehp.h"

#define vector(a,b) ((b*)xmalloc((a)*sizeof(b)))
void haplotype(double *);
void getden(double *,int, boolean, double *, double *);
void probgen(double *, struct DBT *, short *, short *, double *);
int genec(int, int,  double *, struct DBT *, short *, short *, double *);
double likelihood(double *,int,boolean,boolean);
int genecount(int,int,int *,double *,struct DBT *,double *);
void probgeno(double *, struct DBT *, double *);
void getden(double *,int,boolean,double *,double *);

int main(int argc, char **argv)
{

#ifdef useid
double *id,*idsave;
#endif
char line[LL+1],rest[LL+1];
int i,j,n,nloci;
double tid;
short l,u;
double sumca,sumco,sump,sumq,asump,asumq;
double *hap,x2,x2df,lrt,lrtdf;
fout=finp=NULL;

getoptions();
indivN=vector(3*maxobscom,double);
if(cc) {
  caseN=vector(maxobscom,double);
  contrlN=vector(maxobscom,double);
}
#ifdef useid
id=vector(3*maxobscom,double);
idsave=vector(maxobscom,double);
#endif
if(cc) alist=(phenotype *)malloc(3*maxobscom*sizeof(phenotype));
else alist=(phenotype *)malloc(maxobscom*sizeof(phenotype));
if(!alist) {
  perror("Error to allocate the big caches !");
  return 1;
}
finp=fopen(infile_name,"r");
fout=fopen(outfile_name,"w");
if(!finp||!fout) {
  fprintf(stderr,"Could not open file %s %s\n",infile_name,outfile_name);
  return 1;
}
fgets(line,LL,finp);
maxall=nloci=0;
while(sscanf(line,"%d %[^\n]",&loci[nloci++],rest)>1) strcpy(line,rest);
for(i=0;i<nloci;i++) if(loci[i]>maxall) maxall=loci[i];
totalloci=nloci;
initV(nloci,loci);
if(cc) obshap=2*haplnum;
else obshap=haplnum;
totalind=i=0;
while(fgets(line,LL,finp)&&(cc)?
sscanf(line,"%lf %lf %lf %lf %[^\n]",&tid,&indivN[i],&caseN[i],&contrlN[i],rest)>4
:sscanf(line,"%lf %lf %[^\n]",&tid,&indivN[i],rest)>2)
{
#ifdef useid
  idsave[i]=tid;
  id[i]=idsave[i];
#endif
  strcpy(line,rest);
  for (j=0;j<nloci;j++)
  {
      if (sscanf(line,"%hd %hd %[^\n]",&l,&u,rest)<2) return 1;
      strcpy(line,rest);
      *rest='\0';
      if (l>u) hilo(&l,&u);
      alist[i].l[j]=l;
      alist[i].u[j]=u;
  }
  totalind+=indivN[i];
  ++i;
}
fclose(finp);
sample_size=obscom=i;
chklimit();
hap=vector(obshap,double);
oldhaplo=vector(obshap,double);
newhaplo=vector(obshap,double);
inihaplo=vector(obshap,double);
indhaplo=vector(obshap,double);
if(cc)
{
  qq=1-pp;
  ss[0]=1-ff[0];
  ss[1]=1-ff[1];
  ss[2]=1-ff[2];
  sumca=ff[0]*qq*qq+2*ff[1]*pp*qq+ff[2]*pp*pp;
  sumco=ss[0]*qq*qq+2*ss[1]*pp*qq+ss[2]*pp*pp;
  casep[0]  =   ff[0]*qq*qq/sumca;
  casep[1]  = 2*ff[1]*pp*qq/sumca;
  casep[2]  =   ff[2]*pp*pp/sumca;
  controlp[0]=  ss[0]*qq*qq/sumco;
  controlp[1]=2*ss[1]*pp*qq/sumco;
  controlp[2]=  ss[2]*pp*pp/sumco;

  dg11=vector(obscom,double);
  dg12=vector(obscom,double);
  dg22=vector(obscom,double);
  for(i=0;i<obscom;i++)
  {
    dg11[i]=casep[0]*caseN[i]+controlp[0]*contrlN[i];
    dg12[i]=casep[1]*caseN[i]+controlp[1]*contrlN[i];
    dg22[i]=casep[2]*caseN[i]+controlp[2]*contrlN[i];
  }
  getgeno(totalloci);
  initP();
  i=j=0;
  totalind=nget(obscom,totalloci);
  totalall=s2=2*totalind;
  calcP(totalloci);
  getH(totalloci);
  onemark=(boolean)(totalloci==1);
  if(!onemark) {
    loop=0;
    do {
      loop++;
      wd=false;
      haplotype(hap);
      for(i=0;i<haplnum;i++) newhaplo[i]/=totalall;
      loglike=likelihood(newhaplo,totalloci,false,false);
      diff=getdiff();
      err=(loop>maxloops||diff<0.0001);
    } while(!err);
    if(loop>maxloops) printf("Warning: Exceeds iteration limit\007\n");
  }
  for(i=0;i<haplnum;i++) {
    indhaplo[i]=newhaplo[i]*(1-pp);
    indhaplo[i+haplnum]=newhaplo[i]*pp;
  }
  n=1;
  for(i=0;i<totalloci;i++) n*=locigeno[i];
  totalloci++;nloci++;
  for(i=totalloci;i>0;i--) loci[i]=loci[i-1];loci[0]=2;
  for(j=0;j<obscom;j++)
  {
      indivN[j]=dg11[j];
      indivN[obscom+j]=dg12[j];
      indivN[obscom*2+j]=dg22[j];
#ifdef useid
      id[j]=idsave[j];
      id[obscom+j]=n+idsave[j];
      id[obscom*2+j]=2*n+idsave[j];
#endif
      for (i=nloci-2;i>=0;i--)
      {
         l=alist[j].l[i];
         u=alist[j].u[i];
         alist[j].l[i+1]=l;
         alist[j].u[i+1]=u;
         alist[obscom+j].l[i+1]=l;
         alist[obscom+j].u[i+1]=u;
         alist[obscom*2+j].l[i+1]=l;
         alist[obscom*2+j].u[i+1]=u;
      }
      alist[j].l[0]=1;
      alist[j].u[0]=1;
      alist[obscom+j].l[0]=1;
      alist[obscom+j].u[0]=2;
      alist[obscom*2+j].l[0]=2;
      alist[obscom*2+j].u[0]=2;
  }
  initV(nloci,loci);
  loop=0;
  do {
    getgeno(totalloci);
    initP();
    i=j=0;
    totalind=nget(obscom*3,totalloci);
    totalall=s2=2*totalind;
    calcP(totalloci);
    if(loop==0) {
      p[0][0]=1-pp;
      p[0][1]=pp;
      outP(totalloci);
      getH(totalloci);
      iniloglike=likelihood(newhaplo,totalloci,true,false);
      indloglike=likelihood(indhaplo,totalloci,true,false);
      loglike=iniloglike;
    }
    oldlog=loglike;
    loop1=0;
    do {
      wd=true;
      haplotype(hap);
      for(i=0;i<haplnum;i++) newhaplo[i]/=totalall;
      sump=sumq=asump=asumq=0.0;
      n=1;
      for(i=0;i<totalloci;i++) n*=loci[i];
      for(i=0;i<n/2;i++){
        sump+=inihaplo[i];asump+=newhaplo[i];
      }
      for(i=n/2;i<n;i++){
        sumq+=inihaplo[i];asumq+=newhaplo[i];
      }
      for(i=0;i<n/2;i++){
        inihaplo[i]*=(1-pp)/sump;
        newhaplo[i]*=(1-pp)/asump;
      }
      for(i=n/2;i<n;i++){
        inihaplo[i]*=pp/sumq;
        newhaplo[i]*=pp/asumq;
      }
      diff=getdiff();
      loop1++;
    } while(diff>=0.0001&&loop1<=maxloops);
    loop++;
    loglike=likelihood(newhaplo,totalloci,true,true);
    diff=loglike-oldlog;
  } while(diff>=0.0001&&loop<=15);
} else {
  getgeno(totalloci);
  initP();
  i=j=0;
  totalind=nget(obscom,totalloci);
  totalall=s2=2*totalind;
  calcP(totalloci);
  outP(totalloci);
  getH(totalloci);
  loop=0;
  iniloglike=likelihood(newhaplo,totalloci,false,false);
  do {
    loop++;
    wd=false;
    haplotype(hap);
    for(i=0;i<haplnum;i++) newhaplo[i]/=totalall;
    loglike=likelihood(newhaplo,totalloci,false,false);
    diff=getdiff();
    err=(loop>maxloops||diff<0.0001);
  } while(!err);
  x2=0;
  x2df=0;
  lrt=0;
  lrtdf=0;
  for(i=0;i<haplnum;i++) if(inihaplo[i]>0) {
    if(newhaplo[i]>0) {
      lrt+=newhaplo[i]*log(newhaplo[i]/inihaplo[i]);
      lrtdf++;
    }
    x2+=pow(newhaplo[i]-inihaplo[i],2)/inihaplo[i];
    x2df++;
  }
}
outH(loci,newhaplo);
if(!cc) {
#ifdef DEBUG
  fprintf(fout,"(%.2lf, %.0lf %.lf %.2lf %.0lf)",x2*s2,x2df,s2,lrt*s2,lrtdf);
#endif
}
if(loop>maxloops) printf("Warning: bad data set\007\n");
fclose(fout);
free(indivN);
#ifdef useid
free(id);
free(idsave);
#endif
free(oldhaplo);
free(newhaplo);
free(inihaplo);
free(indhaplo);
free(hap);
free(alist);
if(cc) {
free(caseN);
free(contrlN);
free(dg11);
free(dg12);
free(dg22);
}
return 0;
}

void haplotype(double *hap)
/*counting haplotypes*/
{
int i,h,subn,lines,gindex,loops;
double totalg;
struct DBT V;

for(i=0;i<haplnum;i++) {
  oldhaplo[i]=newhaplo[i];
  newhaplo[i]=0.0;
}
loops=(wd)?3*obscom:obscom;
gindex=0;
do {
  V.time=gindex;
  for(i=0;i<totalloci;i++)
     if(alist[V.time].l[i]!=alist[V.time].u[i]) locihn[i]=1;
     else locihn[i]=0;
  h=0;
  for(i=0;i<totalloci;i++) h+=locihn[i];
  if(h<2) {
    V.l=linenum(loci,alist[V.time].l,1)-1;
    V.u=linenum(loci,alist[V.time].u,1)-1;
    newhaplo[V.l]+=indivN[V.time];
    newhaplo[V.u]+=indivN[V.time];
  } else if(h>1) {
    getfirst(&V.fl,&V.ll);
    lines=0;
    totalg=0.0;
    genecount(0,V.fl,&lines,&totalg,&V,hap);
    subn=exp2(h-1);
    if(totalg<eps) for(i=0;i<subn;i++) hap[i]=0.0;
    else for(i=0;i<subn;i++) hap[i]/=totalg;
    lines=0;
    genecount(1,V.fl,&lines,&totalg,&V,hap);
  }
  gindex++;
} while(gindex<loops);
}

int genecount(int opt, int i, int *line, double *totalg, struct DBT *LINK, double *hap)
/*gene counting*/
{
int j;
if(locihn[i-1]!=1) {
  genecount(opt,i+1,line,totalg,LINK,hap);
  return 0;
}
if(i==LINK->fl) {
  if(i!=LINK->ll) {
    genecount(opt,i+1,line,totalg,LINK,hap);
    return 0;
  }
  LINK->l=linenum(loci,alist[LINK->time].l,1)-1;
  LINK->u=linenum(loci,alist[LINK->time].u,1)-1;
  switch(opt) {
  case 0:
    hap[*line]=2*oldhaplo[LINK->l]*oldhaplo[LINK->u];
    *totalg +=hap[*line];
    break;
  case 1:
    newhaplo[LINK->l]+=hap[*line]*indivN[LINK->time];
    newhaplo[LINK->u]+=hap[*line]*indivN[LINK->time];
    break;
  case 2:
    *totalg+=2*hap[LINK->l]*hap[LINK->u];
  }
  (*line)++;
  return 0;
}
for(j=1;j<=2;j++){
  if(j==2) hilo(&alist[LINK->time].l[i-1],&alist[LINK->time].u[i-1]);
  if(i==LINK->ll){
    LINK->l=linenum(loci,alist[LINK->time].l,1)-1;
    LINK->u=linenum(loci,alist[LINK->time].u,1)-1;
    switch(opt) {
    case 0:
      hap[*line]=2*oldhaplo[LINK->l]*oldhaplo[LINK->u];
      *totalg+=hap[*line];
      break;
    case 1:
      newhaplo[LINK->l]+=hap[*line]*indivN[LINK->time];
      newhaplo[LINK->u]+=hap[*line]*indivN[LINK->time];
      break;
    case 2:
      *totalg+=2*hap[LINK->l]*hap[LINK->u];
      break;
    }
    (*line)++;
  } else genecount(opt,i+1,line,totalg,LINK,hap);
}
hilo(&alist[LINK->time].l[i-1],&alist[LINK->time].u[i-1]);return 0;
}

double likelihood(double *hap, int nloci, boolean dma, boolean dg)
/*yes, the likelihoods, P(A|g)P(g)+P(U|g)P(g)*/
{
int i,gindex;
double like, l11, l12, l22, ta, tu, sa, su, ssa, ssu;
struct DBT V;

V.nloci=nloci;
V.dma=dma;
sa=su=ssa=ssu=0.0;
gindex=0;
if(V.dma) do {
   V.time=gindex;
   alist[V.time].l[0]=1;alist[V.time].u[0]=1;
   probgeno(&l11, &V, hap);
   V.time=obscom+gindex;
   alist[V.time].l[0]=1;alist[V.time].u[0]=2;
   probgeno(&l12, &V, hap);
   V.time=obscom+obscom+gindex;
   alist[V.time].l[0]=2;alist[V.time].u[0]=2;
   probgeno(&l22, &V, hap);
   ta=ff[0]*l11+ff[1]*l12+ff[2]*l22;
   tu=ss[0]*l11+ss[1]*l12+ss[2]*l22;
   sa+=ta;
   su+=tu;
   gindex++;
}  while(gindex<obscom);
if(V.dma&&correct_den) getden(hap,nloci,dma,&sa,&su);
gindex=0;
like=0.0;
do {
   V.time=gindex;
   if (V.dma) {
      alist[V.time].l[0]=1;alist[V.time].u[0]=1;
      probgeno(&l11, &V, hap);
      V.time=obscom+gindex;
      alist[V.time].l[0]=1;alist[V.time].u[0]=2;
      probgeno(&l12, &V, hap);
      V.time=obscom+obscom+gindex;
      alist[V.time].l[0]=2;alist[V.time].u[0]=2;
      probgeno(&l22, &V, hap);
      ta=ff[0]*l11+ff[1]*l12+ff[2]*l22;
      tu=ss[0]*l11+ss[1]*l12+ss[2]*l22;
      if(caseN[gindex]>0.01) {
        ssa+=caseN[gindex];
        like+=caseN[gindex]*log(ta);
      }
      if(contrlN[gindex]>0.01) {
        ssu+=contrlN[gindex];
        like+=contrlN[gindex]*log(tu);
      }
      dg11[gindex]=dg12[gindex]=dg22[gindex]=0.0;
      if(caseN[gindex]>0.01){
        dg11[gindex]+=l11*caseN[gindex]*ff[0]/ta;
        dg12[gindex]+=l12*caseN[gindex]*ff[1]/ta;
        dg22[gindex]+=l22*caseN[gindex]*ff[2]/ta;
      }
      if(contrlN[gindex]>0.01){
        dg11[gindex]+=l11*contrlN[gindex]*ss[0]/tu;
        dg12[gindex]+=l12*contrlN[gindex]*ss[1]/tu;
        dg22[gindex]+=l22*contrlN[gindex]*ss[2]/tu;
      }
   } else {
      probgeno(&l11, &V, hap);
      like += l11;
   }
   gindex++;
}  while(gindex<obscom);
if(V.dma&&correct_den) cf=ssa*log(sa)+ssu*log(su);
if(V.dma&&!correct_den) cf=ssa*log(sumcase)+ssu*log(sumcontrol);
if(!dg) return like;
gindex=1;
for(i=1;i<totalloci;i++) gindex*=locigeno[i];
for(i=0;i<obscom;i++)
{
  indivN[i]=dg11[i];
  indivN[obscom+i]=dg12[i];
  indivN[obscom*2+i]=dg22[i];
#ifdef useid
  id[i]=idsave[i];
  id[obscom+i]=gindex+idsave[i];
  id[obscom*2+i]=2*gindex+idsave[i];
#endif
}
return like;
}

void probgeno(double *prob, struct DBT *LINK, double *hap)
/*P(g), g=multilocus genotype*/
{
int i,lines,h;
double totalg,temp;

h=0;
for(i=0;i<LINK->nloci;i++) {
  if(alist[LINK->time].l[i]!=alist[LINK->time].u[i]) locihn[i]=1;
  else locihn[i]=0;
  h+=locihn[i];
}
*prob=0.0;
if(h<2) {
  LINK->l=linenum(loci,alist[LINK->time].l,1)-1;
  LINK->u=linenum(loci,alist[LINK->time].u,1)-1;
  temp=hap[LINK->l]*hap[LINK->u];
  if(temp>eps)
  if(LINK->dma) {
    if(h==0) *prob+=temp;
    else *prob+=2*temp;
  } else {
    if(h==0) *prob+=log(temp)*indivN[LINK->time];
    else *prob+=log(2*temp)*indivN[LINK->time];
  }
} else {
  getfirst(&LINK->fl,&LINK->ll);
  lines=0;
  totalg=0.0;
  genecount(2,LINK->fl,&lines,&totalg,LINK,hap);
  if(totalg>eps) {
    if(LINK->dma) *prob+=totalg;
    else *prob+=log(totalg)*indivN[LINK->time];
  }
}
}

void probgen(double *prob, struct DBT *LINK, short *l, short *u, double *hap)
/*P(g) for P(A) and P(U)*/
{
int i,h;
double totalg,temp;

h=0;
for(i=0;i<LINK->nloci;i++) {
  if(l[i]!=u[i]) locihn[i]=1;
  else locihn[i]=0;
  h+=locihn[i];
}
*prob=0.0;
if(h<2) {
  LINK->l=linenum(loci,l,1)-1;
  LINK->u=linenum(loci,u,1)-1;
  temp=hap[LINK->l]*hap[LINK->u];
  if(temp>eps) {
    if(h==0) *prob+=temp;
    else *prob+=2*temp;
  }
} else {
  getfirst(&LINK->fl,&LINK->ll);
  totalg=0.0;
  genec(2,LINK->fl,&totalg,LINK,l,u,hap);
  if(totalg>eps) if(LINK->dma) *prob+=totalg;
}
}

int genec(int opt, int i, double *totalg, struct DBT *LINK, short *l, short *u, double *hap)
/*genecounting for P(A), P(U)*/
{
int j;
if(locihn[i-1]!=1) {
  genec(opt,i+1,totalg,LINK,l,u,hap);
  return 0;
}
if(i==LINK->fl) {
  if(i!=LINK->ll) {
    genec(opt,i+1,totalg,LINK,l,u,hap);
    return 0;
  }
  LINK->l=linenum(loci,l,1)-1;
  LINK->u=linenum(loci,u,1)-1;
  if(opt==2) *totalg+=2*hap[LINK->l]*hap[LINK->u];
  return 0;
}
for(j=1;j<=2;j++){
  if(j==2) hilo(&l[i-1],&u[i-1]);
  if(i==LINK->ll){
    LINK->l=linenum(loci,l,1)-1;
    LINK->u=linenum(loci,u,1)-1;
    if(opt==2) *totalg+=2*hap[LINK->l]*hap[LINK->u];
  } else genec(opt,i+1,totalg,LINK,l,u,hap);
}
hilo(&l[i-1],&u[i-1]);return 0;
}

void getden(double *hap,int nloci,boolean dma,double *ssa,double *ssu)
/*obtain correct denominator, slow*/
{
struct DBT V;
int i,time;
long int genonum;
double l11, l12, l22, sa, su;
short locik[maxloci+1],lociq[maxloci+1];

V.nloci=nloci;
V.dma=dma;
genonum=1;
for(i=0;i<nloci;i++) {
  lociq[i]=locik[i]=1;
  locihn[i]=0;
  genonum*=locigeno[i];
}
time=1;
sa=su=0.0;
do {
   lociq[0]=locik[0]=1;
   probgen(&l11, &V, locik, lociq, hap);
   lociq[0]=2;
   locik[0]=1;
   probgen(&l12, &V, locik, lociq, hap);
   locik[0]=lociq[0]=2;
   probgen(&l22, &V, locik, lociq, hap);
   sa+=l11*ff[0]+l12*ff[1]+l22*ff[2];
   su+=l11*ss[0]+l12*ss[1]+l22*ss[2];
   locik[V.nloci-1]++;
   for(i=V.nloci-1;i>0;i--) {
     if(locik[i]>lociq[i]) {
       lociq[i]++;
       locik[i]=1;
       if(lociq[i]>loci[i]) {
         lociq[i]=1;
         locik[i-1]++;
         if(i==1&&(locik[i-1]>lociq[i-1])) {
           lociq[i-1]++;
           locik[i-1]=1;
         }
       }
     }
   }
   time+=3;
} while(time<=genonum);
*ssa=sa;
*ssu=su;
}
