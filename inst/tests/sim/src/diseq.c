/***J.H. Zhao & Sham P.C. 03MAR97 I.O.P.            ***/
/******************************************************
 D I S E Q -- disequilibruim optimization problem
*******************************************************/
#include "sim.h"

char *disfile="samples/diseq.dat";
char *locfile="samples/diseq.loc";
char *selped="diseq.out";

int ld_sim(int, char**);

#define NB_MU 5
#define NB_P 0.629
#define TP_MU 3
int getsibs(int,NUC_FAMILY*),dtnb(float,float),dtp(float);
struct trunc_neg_binom {float n,q,pr,p[MAXSIBS],cump[MAXSIBS];} tnb;
struct trunc_poisson {float mu,p[MAXSIBS],cump[MAXSIBS];} tp;

int nfam=150,pedigree;

#define n_ranf() (((double)rand() + 0.5)/((double)RAND_MAX + 1.0))

main(int argc, char *argv[])
{
ld_sim(argc, argv);
}

int ld_sim(int argc, char **argv)
{
FILE *fi,*fo;
NUC_FAMILY nuc_fam[MAXFAM];
float kp,c[4],f0,f1,f2;
float mu,p,q,s,s2,pden,ax,delta,f;
int nreps;
int pid,id,fid,mid,nid,pbnd,sex,aff,liab;
int i,j,k,k1,k2,l,l1,l2,m,m1,m2,skiplines=0;
int op=2;
char line[2000],rest[2000];
long seed=3000,seed1,seed2;

person=(PERSON *)malloc(MAXIND*sizeof(PERSON));
locus=(LOCUS *)malloc(MAXLOCI*sizeof(LOCUS));
typed=imatrix(0,MAXIND,0,MAXLOCI);
locus_order=ivector(0,MAXLOCI);
mg=ivector(0,NMG);
r=vector(0,MAXLOCI-1);
beta=matrix(0,NTRAIT,0,NMG+NPG+NCE+NUE);
kce=matrix(0,NTRAIT,0,NCE);

fi=fopen(locfile,"r");
if (fi==NULL) perror("error opening locus file");
locprep(fi);
fclose(fi);

nobody.pid=0;
nobody.id=0;
nobody.fid=0;
nobody.mid=0;
nobody.sex=UNKNOWN;
nobody.simmed=false;
for(i=0;i<MAXIND;i++) {
  person[i]=nobody;
  for(j=0;j<MAXLOCI;j++) typed[i][j]=1;
}

setall(123,321);
if(argc>1) skiplines=atoi(argv[1]);
setsd(123456780+skiplines*10,987654321-skiplines*10);
fi=fopen(disfile,"r");
for (i=1;i<skiplines;++i) for (j=0;j<3;++j) fgets(line,240,fi);
printf("allele frequency, penetrance and haplotypes\n");
if (fi==NULL)
{
  perror("error reading disequilibrium parameter file\n");exit(0);
}
fgets(line,240,fi);
sscanf(line,"%f",&q);
fgets(line,240,fi);
sscanf(line,"%f %f %f",&f0,&f1,&f2);
fgets(line,240,fi);
sscanf(line,"%f %f %f %f",&c[0],&c[1],&c[2],&c[3]);
fclose(fi);
p=1.0-q;
kp=p*p*f0+2.0*p*q*f1+q*q*f2;
ax=f0-2.0*f1+f2;
delta=f1*f1-f0*f2+ax*kp;
if(ax==0.0) p=(kp-f2)/(2.0*(f1-f2));
else p=(-f1+f2-sqrt(delta))/ax;
q=1.0-p;
#ifdef DEBUG
  printf("%f\n",q);
  printf("%f %f %f\n",f0,f1,f2);
  printf("%f %f %f %f\n",c[0],c[1],c[2],c[3]);
  printf("The disease model: k=%f q=%f f=%.4f %.4f %.4f\n",kp,q,f0,f1,f2);
#endif
diseqloci[0].cumf[0][0]=c[0]/(c[0]+c[1]);
diseqloci[0].cumf[0][1]=1;
diseqloci[0].cumf[1][0]=c[2]/(c[2]+c[3]);
diseqloci[0].cumf[1][1]=1;

getsd(&seed1,&seed2);
#ifdef DEBUG
printf("The current seeds: %d %d\n",seed1,seed2);
#endif

if(op==1) {
  mu=NB_MU;
  p=NB_P;
  dtnb(mu,p);
} else {
  mu=TP_MU;
  dtp(mu);
}
fo=fopen(selped,"w");
if(fo==NULL) perror("error openning output");
selected=0;
pid=0;
do {
  pedigree=selected+1;
  getsibs(op,nuc_fam);
  simped();
  id=0;
  k=0;
  for (l=0;l<pedsize;++l){
    ++id;
    for (j=0;j<nloci;++j) {
      l1=person[id].gtype[j][0];
      l2=person[id].gtype[j][1];
      m=l1+l2;
      switch(locus[j].type) {
      case AFFECTION:
        f=ranf();
        switch(m) {
        case 2:aff=1;if(f<f2) aff=2;break;
        case 3:aff=1;if(f<f1) aff=2;break;
        case 4:aff=1;if(f<f0) aff=2;break;
        default: break;
        }
        person[id].ptype[j].aff.aff=aff;break;
      default: break;
      }
    }
    if (l>1) if(aff==2) ++k;
  }
  if (k>1){
    ++selected;
    outped(fo);
  }
} while (selected<SELECTED_FAMILY);
fclose(fo);

free(person);
free(locus);
free_imatrix(typed,0,MAXIND,0,MAXLOCI);
free_ivector(locus_order,0,MAXLOCI);
free_ivector(mg,0,NMG);
free_vector(r,0,MAXLOCI-1);
free_matrix(beta,0,NTRAIT,0,NMG+NPG+NCE+NUE);
free_matrix(kce,0,NTRAIT,0,NCE);

return 0;
}

/*
History -
draft on 28-12-98
simplified 29-12-98, but problem arises as to ranf()! so n_ranf()
pending for marry-in based on demographic models (POPGEN)
*/

int getsibs(int option,NUC_FAMILY *nuc_fam)
/*get sibships, a case for C++*/
{
int i,j,pid,id,fid,mid,sex;
float s,f;

pid=pedigree;
id=0;
id++;
fid=id;
nuc_fam[pid].pid=pid;
person[id].pid=pid;
person[id].id=id;
person[id].fid=0;
person[id].mid=0;
person[id].sex=MALE;
nuc_fam[pid].father=fid;
id++;
mid=id;
person[id].pid=pid;
person[id].id=id;
person[id].fid=0;
person[id].mid=0;
person[id].sex=FEMALE;
nuc_fam[pid].mother=mid;
f=ranf();
for(j=0;j<MAXSIBS;j++) {
  s=(option==1)?tnb.cump[j]:tp.cump[j];
  if(f>=s) i=j;
}
nuc_fam[pid].nsibs=pedsize=i;
for (j=0;j<pedsize;++j) {
  id++;
  sex=ranf()<SEXRATIO/(SEXRATIO+1.0)?MALE:FEMALE;
  nuc_fam[pid].sib[j]=id;
  person[id].pid=pid;
  person[id].id=id;
  person[id].fid=fid;
  person[id].mid=mid;
  person[id].sex=sex;
}
pedsize+=2;
return 0;
}

int dtnb(float n,float p)
/*tuncated N.B. distribution function*/
{
float p0,pr,s,f;
int i;

tnb.n=n;
tnb.pr=p;
tnb.q=1+p;
pr=pow(tnb.q,-n);
p0=1-pr;
s=0.0;
for (i=1;i<MAXSIBS;i++) {
  pr*=((tnb.n+i-1)/i)*(tnb.pr/tnb.q);
  tnb.p[i]=pr/p0;
  s+=tnb.p[i];
} p0=s;
s=tnb.cump[0]=0.0;
for (i=1;i<MAXSIBS;i++) {
  tnb.p[i]/=p0;
  s+=tnb.p[i];
  tnb.cump[i]=s;
}
#ifdef DEBUG /*trunc_neg_binom*/
printf("Information about truncated nb distribution:\n");
for (i=1;i<MAXSIBS;i++) printf("%d %f %f\n",i,tnb.p[i],tnb.cump[i]);
printf("\n");
#endif
return 0;
}

int dtp(float m)
/*truncated Poisson distribution function*/
{
float p,p0,s,f;
int i;

tp.mu=m;
p=exp(-tp.mu);
p0=1-p;
s=0.0;
for (i=1;i<MAXSIBS;i++) {
  p*=(tp.mu/i);
  tp.p[i]=p/p0;
  s+=tp.p[i];
} p0=s;
s=tp.cump[0]=0.0;
for (i=1;i<MAXSIBS;i++) {
  tp.p[i]/=p0;
  s+=tp.p[i];
  tp.cump[i]=s;
}
#ifdef DEBUG /*trunc_poisson*/
printf("Information about truncated Poisson distribution:\n");
for (i=1;i<MAXSIBS;i++) printf("%d %f %f\n",i,tp.p[i],tp.cump[i]);
printf("\n");
#endif
return 0;
}

/*end of diseq.c*/
