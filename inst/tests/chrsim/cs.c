/*CHRSIM (C version) JH Zhao 1/4/1999 IoP*/
/*genome scanning features 11/4/1999*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define debug
#define mped   250 /*maximum pedigrees per dataset, obsolete*/
#define mind    30 /*maximum number of individuals per pedigree*/
#define mmark  300 /*maximum number of markers to be simulated*/
#define mallel  35 /*maximum number of alleles per marker*/
#define mcross   5 /*maximum number of crossovers per chromosome*/
#define mline  400 /*length of line*/
#define nl      20 /*length of name*/
#define swap(a,b)((tc=a,a=b,b=tc))

FILE *input, *pedin, *pedout, *datout;
char pedinn[30], inputn[30], pedoutn[30], datoutn[30];

typedef struct {char name[12];int nall;float q[mallel],s[mallel];} marker;
typedef struct {int l,o;} selection;
struct cross_process {double l,p,d,c,t[mcross];} poisson;
void isort(selection *,int,int,int);
int setmark(int),cmpf(const float *,const float *);
marker *keep,markers[mmark];
float *pos;
int tc, nfam=0, nloop, nsim=0, seed1, seed2=28007, seed3=24001, avail=0;
int fsize, fsizes[mped], markpos[mmark], gens[mind][mmark*2];
int nfamily, nreplicates, peds[mind][mmark+8], simmed[mind];
int chrom, nmark, mat[mind][mcross], pat[mind][mcross];

int n[23]={461,452,353,280,312,311,272,249,189,281,273,
    249,164,162,145,180,186,136,121,144,61,67,209}; // x was 216
double uran,
morton[]={
 294,266,213,238,211,210,178,172,147,174,150,
 179,130,118,151,151,202,143,147,151,92,81},
browman[]={
 290,269,228,212,198,193,182,168,169,173,148,
 171,115,138,122,134,126,126,105,101,58,62},
dib[]={
 292.7,277.0,233.0,212.0,197.6,201.1,184.0,166.4,166.5,181.7,156.1,
 169.1,117.5,128.6,110.2,130.8,128.7,123.8,109.9,96.5,59.6,58.1};
char errors[][80]={
 "error: seed value must be between 0 and 30,000",
 "error: sum of allele frequencies does not equal 1.0",
 "error: open/read file"};


int getsize()
/*decide family size*/
{
int i=0,j=0,k=0,eof=0;
char ll[101],name[nl],id[nl],l1[nl],l2[nl];
typedef struct {char name[nl],id[mind][nl];int n;} ids;
ids *ff;

ff=(ids*)malloc(mped*sizeof(ids));
if (!ff) perror("error allocating memory in getsize()");
pedin=fopen(pedinn,"r");
while (1) {
  strcpy(ll,"\0");
  if (fgets(ll,100,pedin)==NULL) eof=1;
  sscanf(ll,"%s%s",&name,&id);
  strcpy(l2,name);
  if (eof) strcpy(l2,"\0");
  if (j==0) {
    strcpy(l1,l2);
    goto next;
  }
  if (k==0) strcpy(ff[i].name,name);
  if (!strcmp(l1,l2)) k++;
  else {
    strcpy(l1,l2);
    ff[i].n=k+1;
    ++i;
    k=0;
  }
next:
  strcpy(ff[i].id[k],id);
  j++;
  if (eof) break;
}
fclose(pedin);
#ifdef debug
printf("\nThere are %d families with sizes:\n",i);
for (j=0;j<i;j++) {
  fsizes[j]=ff[j].n;
  printf("%s %d",ff[j].name,ff[j].n);
  for (k=0; k<ff[j].n; k++) printf(" %s",ff[j].id[k]);printf("\n");
}
#endif
free(ff);
return i;
}

int setmark(int nsel)
/*sample markers from Genethon map and place their positions*/
{
FILE *f;
char ll[mline],after[mline],name[12];
time_t seed;
int i,j,k,l;
float q[mallel];
marker *select;
selection *sel;

pos=(float*)malloc(nsel*sizeof(float));
keep=(marker*)malloc(nsel*sizeof(marker));
sel=(selection*)malloc(nsel*sizeof(selection));
select=(marker*)malloc(nsel*sizeof(marker));
if (!pos||!keep||!sel||!select) perror("too many markers perhaps");
/*get a random sequence and positions*/
srand((unsigned) time(&seed));
for (j=0;j<nsel;j++) {
  pos[j]=rand()/(float)RAND_MAX;
  sel[j].l=rand()%n[chrom-1];
  sel[j].o=j;
}
/*an ordered uniform(0,1) positions*/
qsort(pos,nsel,sizeof(float),(int(*)(const void*,const void*))cmpf);
/*assign the order systematically*/
isort(sel,0,nsel-1,0);
sprintf(name,"frq%d",chrom);
f=fopen(name,"r");
if (!f) perror("cannot read marker information");
/*i=selected index, j=alleles, l=record*/
/*this is purely memory efficient*/
i=j=l=0;
while (i<nsel) {
  fgets(ll,mline,f);
  sscanf(ll,"%s %d %[^\n]",&name,&j,after);
  for (k=0;k<j;k++) {
    strcpy(ll,after);
    *after='\0';
    sscanf(ll,"%f %[^\n]",&q[k],after);
  }
  while (sel[i].l==l) {
    strcpy(select[i].name,name);
    select[i].nall=j;
    for (k=0;k<j;k++) select[i].q[k]=q[k];
    i++;
  }
  l++;
}
fclose(f);
/*reverse the order but with newly numbered subset*/
for (i=0;i<nsel;i++) sel[i].l=i;
isort(sel,0,nsel-1,1);
for (k=0;k<nsel;k++) {
  i=sel[k].l;
  keep[k]=select[i];
}
#ifdef debug
printf("\nThe following markers are imaged:\n");
for (k=0;k<nsel;k++) {
  printf("%12s %2d",keep[k].name,keep[k].nall);
  for (j=0;j<keep[k].nall;j++) printf(" %.3f",keep[k].q[j]);
  printf("\n");
}
#endif
free(sel);
free(select);
}

int cmpf(const float *a,const float *b)
{
return((*a>*b)?1:((*a==*b)?0:-1));
}

void isort(selection *a,int from,int to,int type)
/*sort records*/
{
int i,pivot,new_val,new_o,new_l;

if (to>from)
{
  pivot=from;
  for(i=from+1;i<=to;++i) switch(type)
  {
  case 0:
    new_val=a[i].l;
    new_o=a[i].o;
    if(new_val<a[pivot].l)
    {
      a[i]=a[pivot+1];
      a[pivot+1]=a[pivot];
      a[pivot].l=new_val;
      a[pivot].o=new_o;
      pivot++;
    } break;
  case 1:
    new_val=a[i].o;
    new_l=a[i].l;
    if(new_val<a[pivot].o)
    {
      a[i]=a[pivot+1];
      a[pivot+1]=a[pivot];
      a[pivot].o=new_val;
      a[pivot].l=new_l;
      pivot++;
    } break;
  }
  isort(a,from,pivot-1,type);
  isort(a,pivot+1,to,type);
}
}

void initialize()
{
  int i, j, mapfun, mdist;
  float q, sumq, alpha;
  char markch;

  input=(input!=NULL)?freopen(inputn,"r",input):fopen(inputn,"r");
  if (input==NULL) perror(errors[2]);
  fscanf(input, "%d%*[^\n]", &seed1);
  if (seed1==0 || seed1>30000) perror(errors[0]);
  fscanf(input, "%d%f%*[^\n]", &nreplicates, &alpha);
  fscanf(input, "%d%*[^\n]", &mapfun);
  nfamily=getsize();
#ifdef genome_scan
  fscanf(input, "%d%d%*[^\n]", &chrom, &mdist);
  nmark=morton[chrom]/mdist;
  setmark(nmark);
  for (i=0; i<nmark; i++) markpos[i]=i*mdist;
  for (i=0;i<nmark;i++) {
    markers[i]=keep[i];
    printf("%d \n",markpos[i]);
    sumq=0;
    for (j=0; j<markers[i].nall; j++) {
      q=markers[i].q[j];
      markers[i].s[j]=sumq+=q;
    }
    if (fabs(sumq-1.0)>0.01) perror(errors[1]);
  }
#else
  fscanf(input, "%d%d%*[^\n]", &chrom, &nmark);
#ifdef oldcode
  for (i=0; i<nmark; i++) fscanf(input, "%d", &markpos[i]);
  fscanf(input, "%*[^\n]");
  for (i=0; i<nmark; i++) {
    do markch=getc(input); while (markch!=' ');
    fscanf(input, "%d%*[^\n]", &markers[i].nall);
    sumq=0;
    for (j=0; j<markers[i].nall; j++) {
      fscanf(input, "%f", &q);
      markers[i].q[j]=q;
      markers[i].s[j]=sumq+=q;
    }
    if (fabs(sumq-1.0)>0.01) perror(errors[1]);
    fscanf(input, "%*[^\n]");
  }
#else
  setmark(nmark);
  printf("\nThe positions of markers are as follows:\n");
  for (i=0;i<nmark;i++) {
    markpos[i]=pos[i]*morton[chrom-1];
    markers[i]=keep[i];
    printf("%d ",markpos[i]);
    sumq=0;
    for (j=0; j<markers[i].nall; j++) {
      q=markers[i].q[j];
      markers[i].s[j]=sumq+=q;
    }
    if (fabs(sumq-1.0)>0.01) perror(errors[1]);
  } printf("\n");
#endif
#endif
  for (i=0;i<nmark;i++) {
    fprintf(datout,"%s %2d\n",keep[i].name,keep[i].nall);
    for (j=0;j<keep[i].nall;j++) fprintf(datout," %.3f",keep[i].q[j]);
    fprintf(datout,"\n");
  }
  for (i=0;i<nmark;i++) fprintf(datout," %d",markpos[i]);
  fprintf(datout,"\n");
  free(keep);
  fclose(input);
  pedin=(pedin!=NULL)?freopen(pedinn,"r",pedin):fopen(pedinn,"r");
  if (pedin==NULL) perror(errors[2]);
}

void read_oneped()
/*read simped.dat (post-MAKEPED) with availability codes*/
{
  int i, j, k, pid;
  double doubleds=0.0;

  nloop=0;
  for (j=0;j<fsize;j++) {
    fscanf(pedin, "%d", &pid);
    for (i=0;i<8;i++) fscanf(pedin, "%d", &peds[j][i]);
    if (peds[j][7]>1) doubleds += 0.5;
    for (k=8; k<nmark+8; k++)
    if (avail) fscanf(pedin, "%d", &peds[j][k]);
    else peds[j][k]=1;
    fscanf(pedin, "%*[^\n]");
  }
  if (doubleds>0.0) nloop=(long)floor(doubleds+0.5);
  for (i=0; i<fsize; i++) simmed[i]=0;
}

double ranuni(int *seed1)
/*platform-independent uniform pseudorandom generator*/
{
  double r;

  *seed1=*seed1 % 177 * 171 - *seed1 / 177 * 2;
  seed2=seed2 % 176 * 172 - seed2 / 176 * 35;
  seed3=seed3 % 178 * 170 - seed3 / 178 * 63;
  if (*seed1<0) *seed1 += 30269;
  if (seed2<0) seed2 += 30367;
  if (seed3<0) seed3 += 30323;
  r=*seed1 / 30269.0 + seed2 / 30307.0 + seed3 / 30323.0;
  return (r - (long)r);
}

void crossover(int id, int avec, int chrno, int fm[][mcross])
/*place founder cross-overs*/
{
  int ncm, i, j, k, n, minc, t;

  ncm=(long)floor(morton[chrno-1]+0.5);
  for (i=0; i<avec; i++) {
L1: /*to guarantee nonidentical crossovers*/
    k=1;
    uran=ranuni(&seed1);
    for (j=1; j<=ncm; j++) {
      if (uran>=(double)j/ncm) continue;
        k=j;
        break;
    }
    fm[id][i]=k;
    if (i>0) for (n=0;n<i;n++) if (fm[id][i]==fm[fsize-1][n]) goto L1;
  }
  if (avec>1) for (i=0; i<avec-1; i++) {
    minc=i;
    for (j=i+1; j<avec; j++) if(fm[id][j]<fm[id][minc]) minc=j;
    swap(fm[id][i],fm[id][minc]);
  }
  if (avec<5) for (j=avec; j<mcross; j++) fm[id][j]=0;
}

void founders()
/*simulate founder genotypes*/
{
  int i1, i, j, k, l;

  for (i1=0; i1<fsize; i1++) {
    k=0;
    if (peds[i1][1]!=0) continue;
    simmed[i1]=1;
    nsim++;
    for (i=0; i<nmark; i++) {
      uran=ranuni(&seed1);
      l=1;
      for (j=0;j<markers[i].nall;j++) if (uran>=markers[i].s[j]) l++;
      gens[i1][i+k]=l;
      uran=ranuni(&seed1);
      l=1;
      for (j=0;j<markers[i].nall;j++) if (uran>=markers[i].s[j]) l++;
      gens[i1][i+k+1]=l;
      k++;
    }
    for (l=0; l<mcross; l++) mat[i1][l]=pat[i1][l]=0;
  }
}

void nonfounders()
/*simulate nonfounder genotypes*/
{
  int i, j, k, npos, chrp, chrm, ncross, t2[2], fid, mid;

  for (i=0; i<fsize; i++) {
    fid=peds[i][1]-1;
    mid=peds[i][2]-1;
    if (!peds[i][1]||simmed[i]||!simmed[fid]||!simmed[mid]) continue;
    nsim++;
    simmed[i]=1;
    chrp=(ranuni(&seed1)>0.50)?2:1;
    chrm=(ranuni(&seed1)>0.50)?2:1;
    ncross=1;
    do uran=ranuni(&seed1); while(uran>=poisson.t[mcross-1]);
    while (ncross<=mcross) {/*make crossovers*/
      if (uran<poisson.t[ncross-1]) break;
      ncross++;
    }
    crossover(i, ncross, chrom, pat);
    ncross=1;
    do uran=ranuni(&seed1); while(uran>=poisson.t[mcross-1]);
    while (ncross<=mcross) {/*make crossovers*/
      if (uran<poisson.t[ncross-1]) break;
      ncross++;
    }
    crossover(i, ncross, chrom, mat);
    for (j=0; j<nmark; j++) {
      t2[0]=gens[fid][j*2];
      t2[1]=gens[fid][j*2+1];
      npos=1;
      for (k=0; k<mcross; k++)
        if ((pat[i][k]!=0) && (markpos[j]>pat[i][k])) npos++;
      if ((npos & 1)==0) swap(t2[0],t2[1]); /*even intervals*/
      gens[i][j*2]=(chrp==1)?t2[0]:t2[1];
      t2[0]=gens[mid][j*2];
      t2[1]=gens[mid][j*2+1];
      npos=1;
      for (k=0; k<mcross; k++)
        if ((mat[i][k]!=0) && (markpos[j]>mat[i][k])) npos++;
      if ((npos & 1)==0) swap(t2[0],t2[1]); /*even intervals*/
      gens[i][j*2+1]=(chrm==1)?t2[0]:t2[1];
    }
  }
}

main(int argc, char *argv[])
{
  int i, j, i1, j1, k, factden, t2[mmark*2];

  if (argc>1) avail=atoi(argv[1]);
  if (argc>2) strcpy(inputn,argv[2]);
  else strcpy(inputn, "cs.par");
  if (argc>3) strcpy(inputn,argv[3]);
  else strcpy(pedinn, "cs.dat");
  printf("The pedigree data is in post-makeped format named ");
  printf("\"simped.dat\".\n");
  printf("The parameter file is named \"input.dat\".\n");
  printf("The output files are \"datafile.dat\" and ");
  printf("\"pedfile.dat\".\n");
  strcpy(datoutn, "datafile.dat");
  datout=(datout!=NULL)?freopen(datoutn,"w",datout):fopen(datoutn,"w");
  if (datout==NULL) perror(errors[2]);
  initialize();
  printf("\nBeginning simulations . . .\n");
  strcpy(pedoutn, "pedfile.dat");
  pedout=(pedout!=NULL)?freopen(pedoutn,"w",pedout):fopen(pedoutn,"w");
  if (pedout==NULL) perror(errors[2]);
  poisson.l=morton[chrom-1]/100.0;
  for (i=0; i<mcross; i++) poisson.t[i]=0.0;
  poisson.d=0.0;
  for (i=1; i<=mcross; i++) {
    factden=1;
    for (j=i; j>=1; j--) factden *= j;
    poisson.p=pow(poisson.l,(double)i)*exp(-poisson.l)/factden;
    poisson.c=poisson.p/(1-exp(-poisson.l));
    poisson.d+=poisson.c;
    poisson.t[i-1]=poisson.d;
  }
  for (i1=0; i1<nreplicates; i1++) {
  for (j1=0; j1<nfamily; j1++) {
    nfam++;
    fsize=fsizes[j1];
    read_oneped();
    founders();
    do nonfounders(); while (nsim != fsize); /*get all genotypes*/
    nsim=0;
    if (nloop>0) for (i=0; i<nloop;i++) {
      for (j=0; j<fsize; j++) if (peds[j][7]==i && peds[j][1])
        for (k=0; k<nmark*2; k++) t2[k]=gens[j][k];
      for (j=0; j<fsize; j++) if (peds[j][7]==i && !peds[j][1])
        for (k=0; k<nmark*2; k++) gens[j][k]=t2[k];
    }
    for (i=0; i<fsize; i++) {
      fprintf(pedout, "%6d", nfam);
      for(j=0;j<8;j++) fprintf(pedout, (j<6)?"%5d":"%3d", peds[i][j]);
      for (k=0; k<nmark; k++) if (peds[i][k+8]==1)
      fprintf(pedout,"%5d%3d", gens[i][k*2], gens[i][k*2+1]);
      else fprintf(pedout, "    0  0");
      putc('\n', pedout);
    }
  }
  pedin=(pedin!=NULL)?freopen(pedinn, "r", pedin):fopen(pedinn, "r");
  if (pedin==NULL) perror(errors[2]);
  putchar('.');
  }
  printf("\nsimulations done.\n");
  fclose(pedin);
  fclose(datout);
  fclose(pedout);
}

//Broman KW et al (1998). Comphensive human genetic maps: individual
//and sex-specific variation in recombination. AJHG 63:861-869

//David Brunskill and John Turner (1996). Understanding Algorithms
//and Data Structures. The McGraw-Hill Co.

//Dib C et al.(1996). A comp-rehensive genetic map of the human
//genome based on 5,264 microsatellites. Nature 380:152-154

//Morton NE (1991). Parameters of the human genome. PNAS 88: 7474

//Terwilliger JD, Speer M and Ott (1993). Chromosome-based method
//for rapid computer simulations in human genetic linkage analysis.
//Genet Epidemiol 10: 217-224

//Wichman BA, Hill ID (1982). Algorithm 193:  an efficient and
//portable pseudo-random number generator.  Appl Stat 31:188-190

//Full features, other types of markers (QTL,AFF) and disequilibrium
