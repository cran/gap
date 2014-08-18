#include "fpmp.h"

#undef  DEBUG
#define PROGNAME "FASTPMPLUS"
#define VERSION "1.0"
#define DATE "1-MAY-2000"
#define vector(a,b) ((b*)xmalloc((a)*sizeof(b)))

int initrun(void),tidyup(void),runobs(void),runperm(void);

int np0[3],np1[3],np2[3],d[3],np0p[3],np1p[3];
double l0[3],l1[3],l2[3],l0p[3],l1p[3];
double x[3],y[3],n[3];
double pmmu[3],pmvar[3],xihat[3],xihatv[3],ep[5],se[5];
char *mmmodel[3]={"One block association:            ",
                  "One versus two block association: ",
                  "Permuted block association:       "};
char *mm0="chi-squared statistic for marker-association";
char *mm1="Chi-squared statistic for one block association ";
char *mm12="Chi-squared statistic for one versus two block association ";
char *mmp="Chi-squared statistic for permuted block association ";
char *cch="Chi-squared statistic for heterogeneity model  ";
char *ccmodel="T5 - Heterogeneity model:      ";
FILE *fp;
char *outfile;

int main(int argc, char *argv[])
{
time_t t;
int i,j;

printf("Program for marker-marker association analysis");
printf(" %s %s %8s\n\n",PROGNAME,VERSION,DATE);
printf("Maximum number of loci = %d\n", MAX_LOC);
printf("Maximum number of alleles = %d\n", maxalleles);
time(&t);
printf("\nStarting: %s\n",ctime(&t));

if(argc>3) {
  if (!getloci(argv[1])) getdat(argv[2]);
  outfile=argv[3];
  if(argc>4) seed=atoi(argv[4]);
} else {
  fprintf(stderr,"Usage: %s parfile datfile outfile [seed]\n",argv[0]);
  exit(1);
}
initrun();
for(i=0;i<3;i++) n[i]=0;
runobs();
for(i=1;i<=npermute;i++)
{
  iter=i;
  if(iter==1)
  {
    fprintf(fp,"\n");
    fprintf(fp,"Random number seed = %d\n",seed);
    fprintf(fp,"Number of replicates = %d\n\n",npermute);
  } else {
    for(j=0;j<60;j++) printf("\b");
    printf("Running permutation %6d out of %6d ... ",iter,npermute);
  }
  runperm();
}
printf("done\n");
tidyup();

return 0;
}

int initrun(void)
{
int i,j,k1,k2;
long int l;

chklimit();
orig_model.nloci=nloci;
orig_model.nloci_s=nloci_s;
orig_model.selected=selected;
orig_model.selectp=selectp;
orig_model.selectn=selectn;
for(i=0;i<nloci;i++) {
  orig_model.sel[i]=sel[i];
  orig_model.selp[i]=selp[i];
}
j=0;
l=1;
k1=k2=0;
for(i=0;i<nloci;i++) {
  if(!sel[i]) continue;
  l*=alleles[i];
  selidx[j]=orig_model.selidx[j]=i;
  if(selp[i]) {
    selpidx[k1]=orig_model.selpidx[k1]=i;
    k1++;
  } else {
    selnpdx[k2]=orig_model.selnpdx[k2]=i;
    k2++;
  }
  j++;
}
obshap=l;
/***modify***/
if(cc) obshap*=1/*2*/;
if(selected==0) obshap*=2;
#ifdef useid
  id=vector(3*sample_size,double);
  idsave=vector(sample_size,double);
#endif
/***modify***/
/*3*sample_size*/
alist=vector((cc)?sample_size:sample_size,phenotype);
newhaplo=vector(obshap,double);
oldhaplo=vector(obshap,double);
inihaplo=vector(obshap,double);
indhaplo=vector(obshap,double);
indivN=vector(3*sample_size,double);
if(cc) {
  caseN=vector(sample_size,double);
  contrlN=vector(sample_size,double);
}
fp=fopen(outfile,"w");
if(!fp) {
  fprintf(stderr,"I can not open output file\n");
  return -1;
}
printf("Random number seed = %d \n",seed);
srand(seed);

return 0;
}

int runobs(void)
{
int j;

printf("\nAnalysing observed data...");
/*marker-marker and disease-marker*/
if(!cc) {
  ind2eht();
  outehpdat("ehplus.sav",rt,WHOLE);
  if(selectp>=0) outehpdat("ehplus1.sav",rt1,BLOCK1);
  if(selectn>=0) outehpdat("ehplus2.sav",rt2,BLOCK2);
/*whole block*/
  obscom=0;
  ctree(rt,WHOLE);
  eh(WHOLE);
  l0[0]=iniloglike;
  l0[1]=loglike;
  np0[0]=(correct_df)?cnph0:nph0;
  np0[1]=(correct_df)?cnph1:nph1;
#ifdef DEBUG
  ehout("ehplus.out");
#endif
/*all or none selection ends here*/
  if(selectp<0||selectn<0)
  {
      x[0]=2*(l0[1]-l0[0]);d[0]=np0[1]-np0[0];
      if(selected>0&&npermute==0) {
        ld2(fp);
        cp_natural(fp);
      }
      fprintf(fp,"%s = %.2f, df= %d, p= %.4f\n",mm0,x[0],d[0],probchi(x[0],d[0]));
  } else {
/*block permutation*/
 /*intact within each block*/
    ind2eht();
 /*recalculate whole block*/
    obscom=0;
    ctree(rt,WHOLE);
    eh(WHOLE);
    l0p[0]=iniloglike;
    l0p[1]=loglike;
    np0p[0]=(correct_df)?cnph0:nph0;
    np0p[1]=(correct_df)?cnph1:nph1;
 /*only need for the observed*/
    obscom=0;
    ctree(rt1,BLOCK1);
    eh(BLOCK1);
    l1[0]=iniloglike;
    l1[1]=loglike;
    np1[0]=(correct_df)?cnph0:nph0;
    np1[1]=(correct_df)?cnph1:nph1;
#ifdef DEBUG
    ehout("ehplus1.out");
#endif
    obscom=0;
    ctree(rt2,BLOCK2);
    eh(BLOCK2);
    l2[0]=iniloglike;
    l2[1]=loglike;
    np2[0]=(correct_df)?cnph0:nph0;
    np2[1]=(correct_df)?cnph1:nph1;
#ifdef DEBUG
    ehout("ehplus2.out");
#endif
 /*only when multiple markers in the permuted block*/
    x[0]=2*(l0[1]-l0[0]);d[0]=np0[1]-np0[0];
    x[1]=2*(l0[1]-l1[1]-l2[1]);d[1]=np0[1]-np1[1]-np2[1];
    x[2]=2*(l1[1]-l1[0]);d[2]=np1[1]-np1[0];
    fprintf(fp,"%s = %.2f, df= %d, p= %.4f\n",mm1,x[0],d[0],probchi(x[0],d[0]));
    fprintf(fp,"%s = %.2f, df= %d, p= %.4f\n",mm12,x[1],d[1],probchi(x[1],d[1]));
 /*omit only one marker*/
    if(selectp>0)
    fprintf(fp,"%s = %.2f, df= %d, p= %.4f\n",mmp,x[2],d[2],probchi(x[2],d[2]));
  }
} else {
  if(npermute==0) cp_natural(fp);
  ind2eht();
  obscom=0;
  ctree(rt,WHOLE);
/***modify***/
/*cc=true*/
  cc=false;
/***modify***/
/*eh(WHOLE)*/
  outehpdat("ehplus.sav",rt,WHOLE);
  eh(WHOLE);
  l0[0]=iniloglike;
  l0[1]=loglike;
  np0[0]=(correct_df)?cnph0:nph0;
  np0[1]=(correct_df)?cnph1:nph1;
#ifdef DEBUG
  ehout("ehplus.out");
#endif
  totalind=0;
  for(j=0;j<obscom;j++) {
    indivN[j]=caseN[j];
    totalind+=caseN[j];
  }
  eh(WHOLE);
  l1[0]=iniloglike;
  l1[1]=loglike;
  np1[0]=(correct_df)?cnph0:nph0;
  np1[1]=(correct_df)?cnph1:nph1;
#ifdef DEBUG
  ehout("ehplus1.out");
#endif
  totalind=0;
  for(j=0;j<obscom;j++) {
    indivN[j]=contrlN[j];
    totalind+=contrlN[j];
  }
  eh(WHOLE);
  l2[0]=iniloglike;
  l2[1]=loglike;
  np2[0]=(correct_df)?cnph0:nph0;
  np2[1]=(correct_df)?cnph1:nph1;
#ifdef DEBUG
  ehout("ehplus2.out");
#endif
  x[0]=-2*(l0[1]-l1[1]-l2[1]);
  d[0]=(np1[1]+np2[1]-np0[1]);
  fprintf(fp,"%s  = %.2f, df=%d, p=%.4f\n",cch,x[0],d[0],probchi(x[0],d[0]));
  cc=true;
}
printf("done\n");

return 0;
}

int runperm(void)
{
int j;

/*marker-marker and disease-marker*/
if(!cc) {
  pgtype(1);
  ind2eht();
/*whole block*/
  obscom=0;
  ctree(rt,WHOLE);
  eh(WHOLE);
  l0[0]=iniloglike;
  l0[1]=loglike;
  np0[0]=(correct_df)?cnph0:nph0;
  np0[1]=(correct_df)?cnph1:nph1;
/*all or none selection ends here*/
  if(selectp<0||selectn<0)
  {
    y[0]=2*(l0[1]-l0[0]);
    if(y[0]>=x[0]) ++n[0];
    muvar(iter,y[0],&pmmu[0],&pmvar[0]);
  } else {
/*block permutation*/
 /*intact within each block*/
    pgtype(2);
    ind2eht();
 /*recalculate whole block*/
    obscom=0;
    ctree(rt,WHOLE);
    eh(WHOLE);
    l0p[0]=iniloglike;
    l0p[1]=loglike;
    np0p[0]=(correct_df)?cnph0:nph0;
    np0p[1]=(correct_df)?cnph1:nph1;
 /*only when multiple markers in the permuted block*/
    if(selectp>0) {
      pgtype(3);
      ind2eht();
      obscom=0;
      ctree(rt1,BLOCK1);
      eh(BLOCK1);
      l1p[0]=iniloglike;
      l1p[1]=loglike;
      np1p[0]=(correct_df)?cnph0:nph0;
      np1p[1]=(correct_df)?cnph1:nph1;
    }
    y[0]=2*(l0[1]-l0[0]);
    y[1]=2*(l0p[1]-l1[1]-l2[1]);
    if(y[0]>=x[0]) ++n[0];
    if(y[1]>=x[1]) ++n[1];
    muvar(iter,y[0],&pmmu[0],&pmvar[0]);
    muvar(iter,y[1],&pmmu[1],&pmvar[1]);
    if(selectp>0) {
      y[2]=2*(l1p[1]-l1p[0]);
      if(y[2]>=x[2]) ++n[2];
      muvar(iter,y[2],&pmmu[2],&pmvar[2]);
    }
  }
} else {
  plabel();
  ind2eht();
  obscom=0;
  ctree(rt,WHOLE);
/***modify***/
/*cc=true*/
  cc=false;
/***modify***/
/*eh(WHOLE)*/
  totalind=0;
  for(j=0;j<obscom;j++) {
    indivN[j]=caseN[j];
    totalind+=caseN[j];
  }
  eh(WHOLE);
  l1[0]=iniloglike;
  l1[1]=loglike;
  totalind=0;
  for(j=0;j<obscom;j++) {
    indivN[j]=contrlN[j];
    totalind+=contrlN[j];
  }
  eh(WHOLE);
  l2[0]=iniloglike;
  l2[1]=loglike;
  y[0]=-2*(l0[1]-l1[1]-l2[1]);
  d[0]=(np1[1]+np2[1]-np0[1]);
  if(y[0]>=x[0]) ++n[0];
  muvar(iter,y[0],&pmmu[0],&pmvar[0]);
  cc=true;
}

return 0;
}

int tidyup(void)
{
int i,l;

if(npermute>0) if(!cc)
{
  if(selectp<0||selectn<0)
  {
    ep[0]=n[0]/npermute;
    se[0]=sqrt(ep[0]*(1-ep[0])/npermute);
    fprintf(fp,"%s (%.2f) was reached %.0f times\n\n",mm0,x[0],n[0]);
    fprintf(fp,"The empirical p-value is as follows:\n");
    fprintf(fp,"%s ",mmmodel[0]);
    fprintf(fp,"P-value = %.4f ",ep[0]);
  /*if (ep[0]>0) fprintf(fp,"S.E. = %.4f ",se[0]);*/
    fprintf(fp,"\n");
    xihat[0]=xi(sample_size,d[0],pmmu[0],sqrt(pmvar[0]),x[0]);
    xihatv[0]=varxi(sample_size,d[0],xihat[0]);
    fprintf(fp,"xihat=%.4f std(xihat)=%.4f\n\n",xihat[0],sqrt(xihatv[0]));
  }
  else {
    l=(selectp>0)?3:2;
    for (i=0;i<l;i++) {
      ep[i]=n[i]/npermute;
      se[i]=sqrt(ep[i]*(1-ep[i])/npermute);
      xihat[i]=xi(sample_size,d[i],pmmu[i],sqrt(pmvar[i]),x[i]);
      xihatv[i]=varxi(sample_size,d[i],xihat[i]);
    }
    fprintf(fp,"%s (%.2f) was reached %.0f times\n",mm1,x[0],n[0]);
    fprintf(fp,"%s (%.2f) was reached %.0f times\n",mm12,x[1],n[1]);
    if(selectp>0)
    fprintf(fp,"%s (%.2f) was reached %.0f times\n\n",mmp,x[2],n[2]);
    fprintf(fp,"Empirical p-values for these statistics are as follows:\n");
    for (i=0;i<l;i++) {
      fprintf(fp,"%s ",mmmodel[i]);
      fprintf(fp,"P-value = %.4f ",ep[i]);
    /*if (ep[i]>0) fprintf(fp,"S.E. = %.4f ",se[i]);*/
      fprintf(fp,"\n");
      fprintf(fp,"xihat=%.4f std(xihat)=%.4f\n\n",xihat[i],sqrt(xihatv[i]));
    }
  }
} else {
  ep[0]=n[0]/npermute;
  se[0]=sqrt(ep[0]*(1-ep[0])/npermute);
  xihat[0]=xi(sample_size,d[0],pmmu[0],sqrt(pmvar[0]),x[0]);
  xihatv[0]=varxi(sample_size,d[0],xihat[0]);
  fprintf(fp,"%s (%.2f) was reached %.0f times\n\n",cch,x[0],n[0]);
  fprintf(fp,"The empirical p-value is as follows:\n");
  fprintf(fp,"%s ",ccmodel);
  fprintf(fp,"P-value = %.4f ",ep[0]);
  /*if (ep[i]>0) fprintf(fp,"S.E. = %.4f ",se[0]);*/
  fprintf(fp,"\n");
  fprintf(fp,"xihat=%.4f std(xihat)=%.4f\n\n",xihat[0],sqrt(xihatv[0]));
}
fclose(fp);
free(person);
free(person_s);
free(indivN);
if(cc) {
  free(caseN);
  free(contrlN);
}
#ifdef useid
  free(id);
  free(idsave);
#endif
free(alist);
free(newhaplo);
free(oldhaplo);
free(inihaplo);
free(indhaplo);

return 0;
}

void eh(int op)
/*estimating haplotypes, op=0,1,2 for whole block, 1 and 2*/
{
int i,j,n;
short l,u;
double sumca,sumco,sump,sumq,asump,asumq;
double *hap;

hap=vector(obshap,double);
j=0;
for(i=0;i<nloci;i++)
  if(sel[i]) switch(op) {
  case WHOLE:
       loci[j]=alleles[i];
       j++;
       break;
  case BLOCK1:
       if(selp[i]) {
         loci[j]=alleles[i];
         j++;
       }
       break;
  case BLOCK2:
       if(!selp[i]) {
         loci[j]=alleles[i];
         j++;
       }
       break;
  default:
       fprintf(stderr,"Wrong selection in EH\n");
       break;
  }
nloci=j;
maxall=0;
for(i=0;i<nloci;i++) if(loci[i]>maxall) maxall=loci[i];
totalloci=nloci;
if(!cc) if(nloci>1) addum=false;
else {/*add dummy locus*/
  addum=true;
  loci[1]=2;
  nloci++;
  totalloci++;
  for(i=0;i<obscom;i++) alist[i].l[1]=alist[i].u[1]=1;
}
initV(nloci,loci);
if(cc)
{
  pp=freq;
  ff[0]=pen0;
  ff[1]=pen1;
  ff[2]=pen2;
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
  onemark=(BOOL)(totalloci==1);
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
  free(dg11);
  free(dg12);
  free(dg22);
} else {
  getgeno(totalloci);
  initP();
  i=j=0;
  totalind=nget(obscom,totalloci);
  totalall=s2=2*totalind;
  calcP(totalloci);
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
}
if(loop>maxloops) {
#ifdef DEBUG
  fprintf(stderr,"loop over %d for option %d\n",maxloops,op);
  switch(op) {
  case WHOLE:
       outehpdat("ehplus.err",rt,op);
       break;
  case BLOCK1:
       outehpdat("ehplus1.err",rt,op);
       break;
  case BLOCK2:
       outehpdat("ehplus2.err",rt,op);
       break;
  default:
       fprintf(stderr,"This shouldn't happen in EH\n");
       break;
  }
#endif
}
getdf();
nloci=nloci_s=orig_model.nloci_s;
free(hap);
}

int ehout(char *outfile)
{
fout=fopen(outfile,"w");
if(!fout) {
  fprintf(stderr,"Could not open file\n");
  exit(1);
}
outP(totalloci);
outH(loci,newhaplo);
if(fout!=NULL) fclose(fout);
return 0;
}

void alistout(void)
{
int i,j;

for(i=0;i<nloci;i++) printf(" %d",loci[i]);
printf("\n");
printf("\nobscom=%d\n",obscom);
for(i=0;i<obscom;i++) {
  printf("%6.0f",indivN[i]);
  for(j=0;j<nloci;j++) printf(" %d %d",alist[i].l[j],alist[i].u[j]);
  printf("\n");
}
}

void haplotype(double *hap)
/*counting haplotypes*/
{
int i,h,subn,gindex,loops;
long int lines;
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
    subn=(int)exp2(h-1);
    if(totalg<eps) for(i=0;i<subn;i++) hap[i]=0.0;
    else for(i=0;i<subn;i++) hap[i]/=totalg;
    lines=0;
    genecount(1,V.fl,&lines,&totalg,&V,hap);
  }
  gindex++;
} while(gindex<loops);
}

int genecount(int opt, int i, long int *line, double *totalg, struct DBT *LINK, double *hap)
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
hilo(&alist[LINK->time].l[i-1],&alist[LINK->time].u[i-1]);
return 0;
}

double likelihood(double *hap, int nloci, BOOL dma, BOOL dg)
/*yes, the likelihoods, P(A|g)P(g)+P(U|g)P(g)*/
{
int i,gindex;
double like, l11, l12, l22, ta, tu, sa, su, ssa, ssu, cf;
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
        like+=caseN[gindex]*log(ta/sa);
      }
      if(contrlN[gindex]>0.01) {
        ssu+=contrlN[gindex];
        like+=contrlN[gindex]*log(tu/su);
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
if(V.dma) cf=ssa*log(sa)+ssu*log(su);
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
int i,h;
long int lines;
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
  for(i=0;i<nloci;i++)
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

void probgen(double *prob, struct DBT *LINK, short *l, short *u,double *hap)
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

void getden(double *hap,int nloci,BOOL dma,double *ssa,double *ssu)
/*obtain correct denominator, slow*/
{
struct DBT V;
int i;
double time,genonum;
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
   probgen(&l11, &V, locik, lociq,hap);
   lociq[0]=2;
   locik[0]=1;
   probgen(&l12, &V, locik, lociq,hap);
   locik[0]=lociq[0]=2;
   probgen(&l22, &V, locik, lociq,hap);
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

/*Binary Search Trees*/

node *itree(node *r,double genid,int op)
/*insert and sort*/
{
int i,j;
if (r==NULL)
{
  r=(node*)malloc(sizeof(node));
  if (!r)
  {
    fprintf(stderr,"out of memory in itree\n");exit(1);
  }
  r->left=r->right=NULL;
  r->genid=genid;
  r->nca=0;
  r->nco=0;
  if (p_t.affection) r->nca++;
  else r->nco++;
  switch(op) {
  case WHOLE:
       for(i=0;i<=selected;i++) {
         j=selidx[i];
         r->l[i]=p_t.locus[j][0];
         r->u[i]=p_t.locus[j][1];
       }
       break;
  case BLOCK1:
       for(i=0;i<=selectp;i++) {
         j=selpidx[i];
         r->l[i]=p_t.locus[j][0];
         r->u[i]=p_t.locus[j][1];
       }
       break;
  case BLOCK2:
       for(i=0;i<=selectn;i++) {
         j=selnpdx[i];
         r->l[i]=p_t.locus[j][0];
         r->u[i]=p_t.locus[j][1];
       }
       break;
  default:
       fprintf(stderr,"we'll need to make sure of option in itree\n");
       break;
  }
} else
if (genid<r->genid) r->left=itree(r->left,genid,op);
else if (genid>r->genid) r->right=itree(r->right,genid,op);
else
{
  if (p_t.affection) r->nca++;
  else r->nco++;
}
return r;
}

node *stree(node *t,double key)
/*search*/
{
if (!t) return t;
while (t->genid!=key)
{
  if (key<t->genid) t=t->left;
  else t=t->right;
  if (t==NULL) break;
}
return t;
}

node *dtree(node *t,double key)
/*delete*/
{
node *p,*p2;

if (t->genid==key)
{
  if(t->left==t->right)
  {
    free(t);
    return 0;
  }
  else if (t->left==0)
  {
    p=t->right;
    free(t);
    return p;
  }
  else if (t->right==0)
  {
    p=t->left;
    free(t);
    return p;
  }
  else
  {
    p2=t->right;
    p=t->right;
    while (p->left) p=p->left;
    p->left=t->left;
    free(t);
    return p2;
  }
}
if (t->genid<key) t->right=dtree(t->right,key);
else t->left=dtree(t->left,key);
return t;
}

void inorder(node *t)
/*left subtree->t->right subtree*/
{
  if (t) {
     inorder(t->left);
     printf("%d ",t->genid);
     inorder(t->right);
  }
}

void preorder(node *t)
/*t->left subtree->right subtree*/
{
  if (t) {
     printf("%d ",t->genid);
     preorder(t->left);
     preorder(t->right);
  }
}

void postorder(node *t)
/*left subtree->right subtree->t*/
{
  if (t) {
     postorder(t->left);
     postorder(t->right);
     printf("%d ",t->genid);
  }
}

void ptree(node *r,int l,FILE *gdat,int op)
/*print tree inorder*/
{
int i;
if (r) {
   ptree(r->left,l+1,gdat,op);
   fprintf(gdat,"%20.0lf %4d",r->genid,r->nca+r->nco);
   if (cc) fprintf(gdat," %4d %4d",r->nca,r->nco);
   switch(op) {
   case WHOLE:
        for(i=0;i<=selected;i++)
        {
          fprintf(gdat,"%3d%3d",r->l[i],r->u[i]);
        }
        break;
   case BLOCK1:
        for(i=0;i<=selectp;i++)
        {
          fprintf(gdat,"%3d%3d",r->l[i],r->u[i]);
        }
        break;
   case BLOCK2:
        for(i=0;i<=selectn;i++)
        {
          fprintf(gdat,"%3d%3d",r->l[i],r->u[i]);
        }
        break;
   default:
        fprintf(stderr,"what a choice for ptree !\n");
        break;
   }
   fprintf(gdat,"\n");
   ptree(r->right,l+1,gdat,op);
}
}

void ctree(node *r,int op)
/*confer tree to arrays*/
{
int i;
if (r) {
   ctree(r->left,op);
#ifdef usedid
   idsave[obscom]=r->genid;
#endif
   indivN[obscom]=r->nca+r->nco;
   if (cc) {
      caseN[obscom]=r->nca;
      contrlN[obscom]=r->nco;
   }
   switch(op) {
     case WHOLE:
          for(i=0;i<=selected;i++)
          {
            alist[obscom].l[i]=r->l[i];
            alist[obscom].u[i]=r->u[i];
          }
          break;
     case BLOCK1:
          for(i=0;i<=selectp;i++) {
            alist[obscom].l[i]=r->l[i];
            alist[obscom].u[i]=r->u[i];
          }
          break;
     case BLOCK2:
          for(i=0;i<=selectn;i++) {
            alist[obscom].l[i]=r->l[i];
            alist[obscom].u[i]=r->u[i];
          }
          break;
     default:
          fprintf(stderr,"What does this mean in ctree ?\n");
          break;
   }
   obscom++;
   ctree(r->right,op);
}
}

void rtree(node *t)
/*left subtree->right subtree->t*/
{
if (!t) return;
rtree(t->left);
rtree(t->right);
free(t);
}

void etree(int op)
{
switch(op) {
case WHOLE:
     rtree(rt);
     rt=NULL;
     break;
case BLOCK1:
     rtree(rt1);
     rt1=NULL;
     break;
case BLOCK2:
     rtree(rt2);
     rt2=NULL;
     break;
default:;
}
}

int getloci(char *locfile)
/*
retrieves locus information
*/
{
FILE *fp;
int i,j,l,l1,l2,m,n[MAX_LOC],ngtype[MAX_LOC],pg[MAX_LOC],npg[MAX_LOC];
double ll,k1,k2;
char line[LL+1],rest[LL+1];
double kp;

fp=fopen(locfile,"r");
if(!fp)
{
  fprintf(stderr,"Error opening %s\n",locfile);exit(1);
}
fgets(line,LL,fp);
sscanf(line,"%d %d %d %d",&nloci,&cc,&permute,&npermute);
if(nloci>=MAX_LOC)
{
  fprintf(stderr,"Error: maximum number of loci exceeded\n");exit(1);
}
nloci_s=nloci;
fgets(line,LL,fp);
for(i=0;i<nloci;++i,strcpy(line,rest),*rest='\0')
   if(sscanf(line,"%d %[^\n]",&alleles[i],rest)<1) return 0;
fgets(line,LL,fp);
sscanf(line,"%d %d",&isgenotype,&iogenotype);
fgets(line,LL,fp);
for(i=0;i<nloci;++i,strcpy(line,rest),*rest='\0')
   if(sscanf(line,"%d %[^\n]",&sel[i],rest)<1) return 0;
fgets(line,LL,fp);
for(i=0;i<nloci;++i,strcpy(line,rest),*rest='\0')
   if(sscanf(line,"%d %[^\n]",&selp[i],rest)<1) return 0;
if(fgets(line,LL,fp)&&sscanf(line,"%f %f %f %f",&freq,&pen0,&pen1,&pen2)==4)
  if(pen0>pen2) fprintf(stderr,"The order of penetrances is wrong !\n");
fclose(fp);

printf("Number of loci in this analysis = %d \n",nloci);
i=(permute||(!cc&&(npermute>0)));
printf("Permutation procedure %s invoked ",(i)?"will be":"will not be");
if(i) printf("%d times\n",npermute);else printf("\n");
printf("Number of alleles at these loci and their\n");
printf("  selection/permutation statuses [1=yes,0=no]:\n");
for(i=0;i<nloci;++i)
printf("  locus %3d: alleles=%2d selection= %d permutation= %d\n",i+1,alleles[i],sel[i],selp[i]);
if(cc)
{
  printf("The disease model (q,f0,f1,f2) specified: %.4f %.4f %.4f %.4f \n",freq,pen0,pen1,pen2);
  kp=pow(1-freq,2)*pen0+2*freq*(1-freq)*pen1+pow(freq,2)*pen2;
  printf("  amounts to a population disease prevalence of %.4f\n",kp);
}
m=l=0;
l1=l2=0;
for(i=0;i<nloci;++i)
{
  if(sel[i])
  {
    n[l]=alleles[i];
    ngtype[l]=n[l]*(n[l]+1)/2;
    if(selp[i])
    {
      pg[l1]=ngtype[l];
      ++l1;
    }
    else
    {
      npg[l2]=ngtype[l];
      ++l2;
    }
    ++l;++m;
  }
}
selected=m-1;
selectp=l1-1;
selectn=l2-1;
ll=1;
for(i=m;i>0;--i)
{
  j=i-1;
  ll*=ngtype[j];nall[j]=ll;
}
k1=1;
for(i=l1;i>0;--i)
{
  j=i-1;
  k1*=pg[j];np[j]=k1;
}
k2=1;
for(i=l2;i>0;--i)
{
  j=i-1;
  k2*=npg[j];nnp[j]=k2;
}
if(!cc) printf("Blocks 1 and 2 have %d, %d loci\n",l1,l2);
return 0;
}

int getdat(char *datfile)
/*
retrieves data file
*/
{
FILE *fp;
int i,j,k,l,n,a1,a2,gid;
char line[LL+1],rest[LL+1];
link s,t;

fp=fopen(datfile,"r");
if(!fp) {
  fprintf(stderr,"Error opening %s\n",datfile);exit(1);
}
i=n=0;
cases=0;
s=t=NULL;
if(iogenotype) printf("\n   ID  label locus genotype \n\n");
while(fgets(line,LL,fp)
   &&sscanf(line,"%s %d %[^\n]",p_t.id,&p_t.affection,rest)==3)
{
   strcpy(line,rest);
   k=0;
   for(j=0;j<nloci;++j,strcpy(line,rest),*rest='\0')
   {
      if(isgenotype)
      {
        sscanf(line,"%d %[^\n]",&p_t.gtype[j],rest);
        g2a(p_t.gtype[j],&a1,&a2,&gid);
      }
      else
      {
        sscanf(line,"%d %d %[^\n]",&a1,&a2,rest);
        if(a1>a2) _swap_(&a1,&a2);
        p_t.gtype[j]=a2g(a1,a2);
      }
      p_t.locus[j][0]=a1;
      p_t.locus[j][1]=a2;
      if(sel[j]&&p_t.gtype[j]==0) ++k;
   }
   if(iogenotype)
   {
      printf("%5s %3d",p_t.id,p_t.affection);
      l=0;
      for(j=0;j<nloci;++j)
      {
        if(!sel[j]) continue;
        printf(" %6d",p_t.gtype[l]);
        l++;
      }
      printf("\n");
   }
   if(k!=0)
   {
     ++n;continue;
   }
   if(!cc) p_t.affection=false;
   else if(p_t.affection==true) ++cases;
        else p_t.affection=false;
   p_t.no=i;
// person_s[i]=person[i]=p_t;
   t=ilink(p_t,s);
   s=t;
   ++i;
}
fclose(fp);
sample_size=i;
person=(individual *)malloc(sample_size*sizeof(individual));
person_s=(individual *)malloc(sample_size*sizeof(individual));
if(!person||!person_s) {
  fprintf(stderr,"I can not allocate room for individual data\n");
  return -1;
}
clink(s);
printf("There are %d cases out of %d individuals\n",cases,sample_size);
if(n>0) printf("%d records with partial information have been left out \n",n);
return 0;
}

int a2g(int l1,int l2)
/*
converts alleles to genotype
*/
{
int lo,up;
lo=l1; up=l2;
if (l1>l2)
{
   lo=l2;
   up=l1;
}
if(lo==0) return 0;
else return (up*(up-1)/2+lo);
}

int g2a(int s,int *l,int *u,int *t)
/*
  Recover alleles from genotype identifier
  18/01/02 JH Zhao
*/
{
double d;

if(s==0)
{
  *l=*u=*t=0;
  return 1;
}
d = 1 + 8 * (s - 1);
*u = 1 + (int)((1 + (sqrt(d) - 1) - 1) / 2);
*l = s - *u * (*u - 1) / 2;
*t = *l + *u * (*u - 1) / 2;

return 0;

}

double position(int m,int *genotype,int op)

/*(loc1-1)*PROD n[2,...,m]+(loc2-1)*PROD n[3,...,m]+...+locm */

{
int l;
double pos,sum;
pos=sum=0;
for(l=0;l<m+1;++l) if(genotype[l]==0) return(pos);
switch(op)
{
case 0:
for(l=0;l<m;++l) sum+=(genotype[l]-1)*nall[l+1];break;
case 1:
for(l=0;l<m;++l) sum+=(genotype[l]-1)*np[l+1];break;
case 2:
for(l=0;l<m;++l) sum+=(genotype[l]-1)*nnp[l+1];break;
default:break;
}
pos=sum+genotype[m];
return(pos);
}

int plabel(void)
/*
permutes case-control labels
*/
{
int n,i,id,l,*x;
float r;

x=vector(sample_size,int);
n=sample_size;
for(i=0;i<n;i++) x[i]=person[i].no;
while(n>0)
{
  n-=1;
  r=rand()/(float)RAND_MAX;
  id=(int)floor(n*r+1);
  l=x[id];
  x[id]=x[n];
  x[n]=l;
  p_t=person[id];
  person[id]=person[n];
  person[n]=p_t;
}
for(i=0;i<sample_size;++i) person[i].affection=(i<cases)?1:0;
free(x);

return 0;
}

int pgtype(int op)
/*
permutes genotypes
*/
{
int n,i,j,id,k,l,*x,*gtype,*l1,*l2;
int a1,a2;
float r;

x=vector(sample_size,int);
gtype=vector(sample_size,int);
l1=vector(sample_size,int);
l2=vector(sample_size,int);

for(i=0;i<sample_size;++i) person[i]=person_s[i];

switch(op)
{
case 1: /*everything*/
     k=0;
     for(i=0;i<nloci;++i)
     {
        n=sample_size;
        if(!sel[i]) continue;
        ++k;
        if(k==1) continue;
        while(n>0)
        {
          n--;
          r=rand()/(float)RAND_MAX;
          id=(int)floor(n*r+1);
          l=person[id].gtype[i];
          a1=person[id].locus[i][0];
          a2=person[id].locus[i][1];
          person[id].gtype[i]=person[n].gtype[i];
          person[id].locus[i][0]=person[n].locus[i][0];
          person[id].locus[i][1]=person[n].locus[i][1];
          person[n].gtype[i]=l;
          person[n].locus[i][0]=a1;
          person[n].locus[i][1]=a2;
        }
     } break;
case 2:/*blockwise, so only change for two blocks together*/
     n=sample_size;
     for(l=0;l<n;l++) x[l]=person_s[l].no;
     while(n>0)
     {
        n--;
        r=rand()/(float)RAND_MAX;
        id=(int)floor(n*r+1);
        l=x[id];
        x[id]=x[n];
        x[n]=l;
     }
     for(i=0;i<nloci;++i) if(sel[i]) if(selp[i])
     {
       for(l=0;l<sample_size;++l) {
         gtype[l]=person[l].gtype[i];
         l1[l]=person[l].locus[i][0];
         l2[l]=person[l].locus[i][1];
       }
       for(l=0;l<sample_size;++l) {
         j=x[l];
         person[l].gtype[i]=gtype[j];
         person[l].locus[i][0]=l1[j];
         person[l].locus[i][1]=l2[j];
       }
     } break;
case 3:/*everything in one block*/
     for(i=0;i<nloci;++i) if(sel[i]) if(selp[i])
     {
       n=sample_size;
       for(l=0;l<n;l++) x[l]=person_s[l].no;
       while(n>0)
       {
          n--;
          r=rand()/(float)RAND_MAX;
          id=(int)floor(n*r+1);
          l=x[id];
          x[id]=x[n];
          x[n]=l;
       }
       for(l=0;l<sample_size;++l) {
         gtype[l]=person[l].gtype[i];
         l1[l]=person[l].locus[i][0];
         l2[l]=person[l].locus[i][1];
       }
       for(l=0;l<sample_size;++l) {
         j=x[l];
         person[l].gtype[i]=gtype[j];
         person[l].locus[i][0]=l1[j];
         person[l].locus[i][1]=l2[j];
       }
     } break;
default:
     fprintf(stderr,"This option hasn't been implemented in pgtype yet\n");
     break;
}
free(x);
free(gtype);
free(l1);
free(l2);

return 0;
}

int ind2eht(void)
{
int i,j,l0,l1,l2;
int genotype[3][MAX_LOC];
double at0,at1,at2;

etree(WHOLE);
etree(BLOCK1);
etree(BLOCK2);
for(i=0;i<sample_size;++i)
{
   p_t=person[i];
   l0=l1=l2=0;
   for(j=0;j<nloci;++j) if(sel[j])
   {
     genotype[0][l0]=p_t.gtype[j];
     ++l0;
     if(selp[j])
     {
       genotype[1][l1]=p_t.gtype[j];
       ++l1;
     }
     if(!selp[j])
     {
       genotype[2][l2]=p_t.gtype[j];
       ++l2;
     }
   }
/* for(j=0;j<=selected;++j) printf("%d ",genotype[0][j]);printf("\n");*/
   at0=position(selected,genotype[0],WHOLE);
   if(at0>0)
   {
     rt=itree(rt,at0,WHOLE);
     if(!cc)
     {
       at1=position(selectp,genotype[1],BLOCK1);
       rt1=itree(rt1,at1,BLOCK1);
       at2=position(selectn,genotype[2],BLOCK2);
       rt2=itree(rt2,at2,BLOCK2);
     }
   }
}
/*
ptree(rt,0,stdout,WHOLE);
ptree(rt1,0,stdout,BLOCK1);
ptree(rt2,0,stdout,BLOCK2);
*/
return 0;
}

int outehpdat(char *outfile, node *r, int op)
{
int j;

fout=fopen(outfile,"w");
if (!fout) {
   fprintf(stderr,"Can't open file for output !!!\n");
   return 1;
}
for(j=0;j<nloci;++j)
  if(sel[j]) switch(op) {
  case WHOLE:
       fprintf(fout,"%2d ",alleles[j]);
       break;
  case BLOCK1:
       if(selp[j]) fprintf(fout,"%2d ",alleles[j]);
       break;
  case BLOCK2:
       if(!selp[j]) fprintf(fout,"%2d ",alleles[j]);
       break;
  default:
       fprintf(stderr,"what a choice in outehpdat ? \n");
       break;
  }
fprintf(fout,"\n");
ptree(r,0,fout,op);
fclose(fout);

return 0;
}

int getehpdat(void)
/*read ehplus.dat*/
{
char line[LL+1],rest[LL+1];
int i,j;
short l,u;
double tid;

finp=fopen("ehplus.dat","r");
if(!finp) {
  fprintf(stderr,"Could not open file in getehpdat\n");
  return 1;
}
fgets(line,LL,finp);
nloci=0;
while(sscanf(line,"%d %[^\n]",&loci[nloci++],rest)>1) strcpy(line,rest);
i=0;
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
obscom=i;

return 0;
}

double probnorm(double x)
/*
Computes the upper one sided tail probability of the normal distribution
for a given normal deviate, x. After formula 26.2.16 in Abramowitz and
Stegun.
*/
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

double probchi(double x2, int ndf)
/*
INCLUDE file for Linkage Utility programs.
Computes the error probability, probchi = 1-f, where f is the
chi-square distribution function with ndf degrees of freedom
evaluated at x2. Function required: probnorm.
Follows formulas 26.4.4 and 26.4.21 in Abramowitz and Stegun,
"Handbook of mathematical functions", Dover, New York 1968
*/
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

int muvar(int n, double x, double *m1, double *v1)
/*28-4-2000, recursive scheme for mean and variance*/
{

if (n==1) {
   *m1=x;
   *v1=0;
}
else {
   *v1=(*v1*(n-2)+x*x-pow(*m1*(n-1)+x,2)/n)/(n-1)+*m1**m1;
   *m1=(*m1*(n-1)+x)/n;
}
return 0;
}

/*LD measure and notations based on

Zhao H, Pakstis AJ, Kidd JR, Kidd KK (1999) Assessing linkage disequilibrium in
a complex genetic system. I. Overall deviation from random association. Ann.
Hum. Genet. 63:167-179

n = number of individuals
t = sample loglikelihood ratio test statistic
nu = degrees of freedom
mu = empirical mean
sigma = empirical standard deviation

*/
double xi(int n, double nu, double mu, double sigma, double t)
{

if(n>0) {
  if(sigma>0) return sqrt(2.0*nu)/n*(t-mu)/sigma;
  else {
    fprintf(stderr,"Sample variance is nonpositive in xi\n");
    return 1;
  }
}
else {
  fprintf(stderr,"The # n is wrong in xi\n");
  return 1;
}
}

double varxi(int n, double nu, double xihat)
{
double s;

if(n>0) {
  s=(2.0*nu+4.0*n*xihat)/n/n;
  if(s>0) return s;
  else {
    fprintf(stderr,"xi has nonpositive variance\n");
    return 1;
  }
}
else {
  fprintf(stderr,"The # n is wrong in varxi\n");
  return 1;
}
}
/*
EHPLUS/FAST - Faster program for allelic association analysis

JH Zhao (IOP/UK) j.zhao@iop.kcl.ac.uk

History.

21-09-1999 fasteh with allele caching
28-04-2000 add permutation of genotypes
08-05-2000 fix pgtype
13-05-2000 disequilibrium measures
14-05-2000 all-subset analysis
18-05-2000 machine-dependent sample size
19-05-2000 clarify block/segment statistics
30-05-2000 fix P(A) and P(U) for cc analysis
31-05-2000 get rid of hap[] in DBT structure and use obshap
01-06-2000 add all-subset option for heterogeneity analysis
02-06-2000 fix problem generating large win386.swp and example ok on Alpha
03-06-2000 change genotype identifier from long int to double
04-06-2000 fix a2g when there is partial genotyping
08-06-2000 fix haplnum/obshap for cc==true and add LL
20-02-2001 change twobytwo and kbyl according to 2LD
09-10-2001 suppress repeated calls to eh(WHOLE), initial values by Dirichlet
18-02-2002 add g2a
02-04-2002 change itree definition
03-04-2002 disable nrutil.c/nrutil.h, add initrun/runobs/runperm/tidyup
04-04-2002 check out i=j=0 after initP(), comparison of vector with N.R. (ok!)
16-04-2002 check out error in dimension of haplotype when cc&selected==0

Things to do.

. uncomment ***modify*** for obshap and alist of cc when necessary
. genetic distance by haplotypes ?
. add HWE/Workman-Niswander test ?

*/
