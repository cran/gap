#include <sys/stat.h>
#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include <time.h>

#define PROGNAME "PM+"
#define VERSION "1.2"
#define DATE "14-APR-2002"
#define MAX_LOC 30
#define MAX_IND 800

#define MAXEVALS 5
#define MAXGRIDLINES 20
#define QTOL 0.001
#define _swap_(a,b) __swap__(a,b)
void __swap__(int *a,int *b) {int t;t=*a;*a=*b;*b=t;}
typedef enum boolean{false,true} boolean;
struct item_type {
  int id,all,cases,controls;
  short l[MAX_LOC],u[MAX_LOC];
  struct item_type *next;
};
typedef struct item_type item;
typedef item *link;

/*internal procedures*/

void set_all(link,int);
int pm(FILE *),getloci(char*),getdat(char*),plabel(),pgtype(int);
int ind2eh(int,int),ind2ehs(int,int),copyfile(char*,char*);
int a2g(int,int),position(int,int*,int),addummy(char*),adjustn(int,int);
int calleh(),mfeh(float*,int*);
double probnorm(double),probchi(double, int);
int get_ln(float*,int*,float*,int*),get_x2(float*,int*,int);

/*internal variables and arrays*/

int nloci,alleles[MAX_LOC],permute,npermute,valids;
int nall[MAX_LOC],np[MAX_LOC],nnp[MAX_LOC];
boolean sel[MAX_LOC],selp[MAX_LOC],isgenotype,iogenotype,iscase;
int selected,selectn,selectp,rows,cols;
float freq,pen0,pen1,pen2;
int sample_size,cases;
int seed=3000;
typedef struct
{
  char id[10];
  int no;
  boolean affection;
  int locus[MAX_LOC][2];
  int gtype[MAX_LOC];
} individual;
individual p_t,person[MAX_IND],person_s[MAX_IND];

boolean useehplus=true,cc,usefehp=true;
link insert(int,link,int,individual);
item *used;
char cmd1[240],cmd2[240],cmd3[240];

main(int argc, char *argv[])
{
FILE *fo;
char *outfile;
time_t t;

printf("Permutation & Model-free analysis");
printf(" %s %s %25s\n\n",PROGNAME,VERSION,DATE);
time(&t);
printf("Starting: %s\n",ctime(&t));
printf("Maximum number of loci = %d\n", MAX_LOC);
printf("Maximum number of individuals = %d\n\n", MAX_IND);

if(argc>3) {
  getloci(argv[1]);
  getdat(argv[2]);
  outfile=argv[3];
  if(argc>4) seed=atoi(argv[4]);
} else {
  fprintf(stderr,"Usage: %s parfile datfile outfile [seed]\n",argv[0]);
  exit(1);
}
printf("Random number seed = %d \n",seed);
srand(seed);

fo=fopen(outfile,"w+");
if(!fo) fprintf(stderr,"Error opening %s \n",outfile);
pm(fo);
fclose(fo);
printf("\nOutput has been written to %s\n",outfile);
return 0;
}

int pm(FILE *fp)
{
FILE *fo;
int i,j,np0,np1,np2,d[5],op,n0,n1,n2;
float x[5],y[5],n[5],xlnl0,xlnl1,xlnl2,ylnl0,ylnl1,ylnl2_a,ylnl2_b;
float t0,t1,t2,s1,s2;
float p[5],se[5];
char *mmmodel[3]={"One block association:            ",
                  "One versus two block association: ",
                  "Two block versus no association:  "};
char *ccmodel[5]={"T1 - User specified model:     ",
                  "T2 - Mendelian recessive model:",
                  "T3 - Mendelian dominant model: ",
                  "T4 - Model-free analysis:      ",
                  "T5 - Heterogeneity model:      "};

sprintf(cmd1,"eh<scratch >pmplus.tmp");
sprintf(cmd2,"ehplus<scratch >pmplus.tmp");
sprintf(cmd3,"fehp<scratch >fehp.tmp");
printf("\nAnalysing observed data...");
if(!cc) /*to allow for all-or-none selection and 1-marker*/
{
  if(useehplus) {
    ind2ehs(selected,0);
    copyfile("ehplus.dat","ehplus.sav");
  }
  else {
    ind2eh(selected,0);
    copyfile("eh.dat","eh.sav");
  }
  get_ln(&xlnl0,&np0,&xlnl2,&np2);
  if(selected==0) np0=np2=np0-1;
  if(selectp<0||selectn<0)
  {
    x[0]=2*(xlnl2-xlnl0);d[0]=np2-np0;
    fprintf(fp,"The chi-squared statistic for marker-association = %.2f, df= %d, p= %.4f\n",x[0],d[0],probchi(x[0],d[0]));
  }
  else
  {
    if(useehplus) {
      ind2ehs(selectp,1);
      copyfile("ehplus.dat","ehplus1.sav");
    }
    else {
      ind2eh(selectp,1);
      copyfile("control.dat","control1.sav");
    }
    get_ln(&s1,&n0,&t1,&n1);
    if(selectp==0) n1=n0-1;
    if(useehplus) {
      ind2ehs(selectn,2);
      copyfile("ehplus.dat","ehplus2.sav");
    }
    else {
      ind2eh(selectn,2);
      copyfile("control.dat","control2.sav");
    }
    get_ln(&s2,&n0,&t2,&n2);
    if(selectn==0) n2=n0-1;
    xlnl1=t1+t2;
    np1=n1+n2;
    x[0]=2*(xlnl2-xlnl0);d[0]=np2-np0;
    x[1]=2*(xlnl2-xlnl1);d[1]=np2-np1;
    x[2]=2*(xlnl1-xlnl0);d[2]=np1-np0;
#ifdef DEBUG
    printf("\nln(L0), ln(L1)[ln(L0)], lnl(L2) for observed data: \n\n");
    printf("%f %f[%f] %f \n",xlnl0,xlnl1,s1+s2,xlnl2);
#endif
    fprintf(fp,"Chi-squared statistic for one block association            = %.2f, df= %d, p= %.4f\n",x[0],d[0],probchi(x[0],d[0]));
    fprintf(fp,"Chi-squared statistic for one versus two block association = %.2f, df= %d, p= %.4f\n",x[1],d[1],probchi(x[1],d[1]));
    fprintf(fp,"Chi-squared statistic for two block versus no association  = %.2f, df= %d, p= %.4f\n",x[2],d[2],probchi(x[2],d[2]));
  }
}
else /*to be looked after within mfeh*/
{
  if(useehplus) {
    ind2ehs(selected,0);
    copyfile("ehplus.dat","ehplus.sav");
  }
  else {
    ind2eh(selected,0);
    copyfile("case.dat","case.sav");
    copyfile("control.dat","control.sav");
    copyfile("eh.dat","eh.sav");
  }
  mfeh(x,d);
  copyfile("mfeh.out","mfeh.sav");
  free(used);
  fprintf(fp,"Chi-squared statistic for user-specified model = %.2f, df=%d, p=%.4f\n",x[0],d[0],probchi(x[0],d[0]));
  fprintf(fp,"Chi-squared statistic for recessive model      = %.2f, df=%d, p=%.4f\n",x[1],d[1],probchi(x[1],d[1]));
  fprintf(fp,"Chi-squared statistic for dominant model       = %.2f, df=%d, p=%.4f\n",x[2],d[2],probchi(x[2],d[2]));
  fprintf(fp,"Chi-squared statistic for model-free analysis  = %.2f, df=%d, p=%.4f\n",x[3],d[3],probchi(x[3],d[3]));
  fprintf(fp,"Chi-squared statistic for heterogeneity model  = %.2f, df=%d, p=%.4f\n",x[4],d[4],probchi(x[4],d[4]));
}
printf("done\n");
if((cc&&!permute)||npermute<=0) return 0;
#ifdef KEEP_X2
fo=fopen("chisq.out","w");
if(!fo)
{
  perror("Could not open chisq.out");exit(1);
}
fprintf(fo,"\nChi-square statistic(s) from the permuted data:\n\n");
#endif
for(i=0;i<5;++i) n[i]=0;
for(i=0;i<npermute;i++)
{
   for(j=0;j<60;j++) printf("\b");
   printf("Running permutation %6d out of %6d ... ",i+1,npermute);
   if(!cc)
   {
/* Permute everything to get lnl1 and lnl2; lnl0 as observed */
     pgtype(1);
     if(useehplus) ind2ehs(selected,0);else ind2eh(selected,0);
     get_ln(&ylnl0,&n0,&ylnl2_a,&n2);
     if(selectp<0||selectn<0)
     {
       y[0]=2*(ylnl2_a-ylnl0);
#ifdef KEEP_X2
       fprintf(fo,"%5d  %.2f \n",i+1,y[0]);
#endif
       if(y[0]>=x[0]) ++n[0];
       continue;
     }
     if(useehplus) ind2ehs(selectp,1);else ind2eh(selectp,1);
     get_ln(&s1,&n0,&t1,&n1);
     if(useehplus) ind2ehs(selectn,2);else ind2eh(selectn,2);
     get_ln(&s2,&n0,&t2,&n2);
     ylnl1=t1+t2;
/* only one block to get lnl2; lnl0 and lnl1 as observed */
     pgtype(2);
     if(useehplus) ind2ehs(selected,0);else ind2eh(selected,0);
     get_ln(&ylnl0,&n0,&ylnl2_b,&n2);
     y[0]=2*(ylnl1-xlnl0);
     y[1]=2*(ylnl2_b-xlnl1);
     y[2]=2*(ylnl2_a-xlnl0);
     if(y[0]>=x[0]) ++n[0];
     if(y[1]>=x[1]) ++n[1];
     if(y[2]>=x[2]) ++n[2];
#ifdef DEBUG
     if(i==0) printf("\nln(L0), ln(L1)[Ln(L0)], ln(L2)[all] ln(L2)[block]\n\n");
     printf("%f %f[%f] %f %f\n",ylnl0,ylnl1,s1+s2,ylnl2_a,ylnl2_b);
#endif
#ifdef KEEP_X2
     fprintf(fo,"%5d  %.2f %.2f %.2f \n",i+1,y[0],y[1],y[2]);
#endif
   }
   else
   {
     plabel();
     if(useehplus) ind2ehs(selected,0);else ind2eh(selected,0);
     mfeh(y,d);
     free(used);
#ifdef KEEP_X2
     fprintf(fo,"%5d  %.2f %.2f %.2f %.2f %.2f \n",i+1,y[0],y[1],y[2],y[3],y[4]);
#endif
     for(j=0;j<5;++j) if(y[j]>=x[j]) ++n[j];
   }
}
#ifdef KEEP_X2
fclose(fo);
#endif
printf("done\n");

fprintf(fp,"\n");
fprintf(fp,"Random number seed = %d\n",seed);
fprintf(fp,"Number of replicates = %d\n\n",npermute);

if(!cc)
{
  if(selectp<0||selectn<0)
  {
    p[0]=n[0]/npermute;
    se[0]=sqrt(p[0]*(1-p[0])/npermute);
    fprintf(fp,"Chi-square statistic for marker-association (%.2f) was reached %.0f times\n\n",x[0],n[0]);
    fprintf(fp,"The empirical p-value is as follows:\n");
    fprintf(fp,"%s ",mmmodel[0]);
    fprintf(fp,"P-value = %.4f ",p[0]);
  /*if (p[0]>0) fprintf(fp,"S.E. = %.4f ",se[0]);*/
    fprintf(fp,"\n");
  }
  else {
    for (i=0;i<3;i++) {
      p[i]=n[i]/npermute;se[i]=sqrt(p[i]*(1-p[i])/npermute);
    }
    fprintf(fp,"One block association chi-squared statistic (%.2f) was reached %.0f times\n",x[0],n[0]);
    fprintf(fp,"One block versus two block association chi-squared statistic (%.2f) was reached %.0f times\n",x[1],n[1]);
    fprintf(fp,"Two block versus no association chi-squared statistic (%.2f) was reached %.0f times\n\n",x[2],n[2]);
    fprintf(fp,"The empirical p-values for these statistics are as follows:\n");
    for (i=0;i<3;i++) {
      fprintf(fp,"%s ",mmmodel[i]);
      fprintf(fp,"P-value = %.4f ",p[i]);
    /*if (p[i]>0) fprintf(fp,"S.E. = %.4f ",se[i]);*/
      fprintf(fp,"\n");
    }
 }
}
else {
  for (i=0;i<5;i++) {
    p[i]=n[i]/npermute;se[i]=sqrt(p[i]*(1-p[i])/npermute);
  }
  fprintf(fp,"User-specified model chi-squared statistic (%.2f) was reached %.0f times\n",x[0],n[0]);
  fprintf(fp,"Recessive model chi-squared statistic (%.2f) was reached %.0f times\n",x[1],n[1]);
  fprintf(fp,"Dominant model chi-squared statistic (%.2f) was reached %.0f times\n",x[2],n[2]);
  fprintf(fp,"Model-free chi-squared statistic (%.2f) was reached %.0f times\n",x[3],n[3]);
  fprintf(fp,"Heterogeneity model chi-squared statistic (%.2f) was reached %.0f times\n\n",x[4],n[4]);
  fprintf(fp,"Empirical p-values for these statistics are as follows:\n");
  for (i=0;i<5;i++) {
    fprintf(fp,"%s ",ccmodel[i]);
    fprintf(fp,"P-value = %.4f ",p[i]);
  /*if (p[i]>0) fprintf(fp,"S.E. = %.4f ",se[i]);*/
    fprintf(fp,"\n");
  }
}
return 0;
}

int getloci(char *locfile)
/*
retrieves locus information
*/
{
FILE *fp;
int i,j,l,l1,l2,k1,k2,m,n[MAX_LOC],ngtype[MAX_LOC],pg[MAX_LOC],npg[MAX_LOC];
char line[241],rest[241];
float kp;

fp=fopen(locfile,"r");
if(!fp)
{
  fprintf(stderr,"Error opening %s",locfile);exit(1);
}
fgets(line,240,fp);
sscanf(line,"%d %d %d %d",&nloci,&cc,&permute,&npermute);
if(nloci>=MAX_LOC)
{
  perror("Error: maximum number of loci exceeded");
}
fgets(line,240,fp);
for(i=0;i<nloci;++i,strcpy(line,rest),*rest='\0')
   if(sscanf(line,"%d %[^\n]",&alleles[i],rest)<1) return 0;
fgets(line,240,fp);
sscanf(line,"%d %d",&isgenotype,&iogenotype);
fgets(line,240,fp);
for(i=0;i<nloci;++i,strcpy(line,rest),*rest='\0')
   if(sscanf(line,"%d %[^\n]",&sel[i],rest)<1) return 0;
fgets(line,240,fp);
for(i=0;i<nloci;++i,strcpy(line,rest),*rest='\0')
   if(sscanf(line,"%d %[^\n]",&selp[i],rest)<1) return 0;
if(fgets(line,240,fp)&&sscanf(line,"%f %f %f %f",&freq,&pen0,&pen1,&pen2)==4)
  if(pen0>pen2) perror("The order of penetrances is wrong !");
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
#ifdef DEBUG
for(l=0;l<m;++l)
{
   printf("alleles of selected marker: %d\n",n[l]);
}
#endif
l=1;
for(i=m;i>0;--i)
{
   j=i-1;
   l*=ngtype[j];nall[j]=l;
#ifdef DEBUG
   printf("alleles %2d nall %d\n",n[j],nall[j]);
#endif
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
if(!cc) printf("blocks 1 and 2 have  %d, %d loci\n",l1,l2);
return 0;
}

int getdat(char *datfile)
/*
retrieves data
*/
{
FILE *fp;
int i,j,k,n,a1,a2,gid;
char line[1000],rest[1000];

if((fp=fopen(datfile,"r"))==NULL) fprintf(stderr,"Error opening %s",datfile);
i=0;n=0;
cases=0;
if(iogenotype) printf("\n   ID  label locus genotype \n\n");
while(fgets(line,1000,fp)
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
      if(sel[j]&&(a1==0||a2==0)) ++k;
   }
   if(iogenotype)
   {
      printf("%5s %3d",p_t.id,p_t.affection);
      for(j=0;j<nloci;++j)
      printf(" %6d",p_t.gtype[j]);
      printf("\n");
   }
   if(k!=0)
   {
     ++n;continue;
   }
   if(!cc) p_t.affection=false;
   else if(p_t.affection==true) ++cases; else p_t.affection=false;/*not 0,1*/
   p_t.no=i;
   person_s[i]=person[i]=p_t;
   if(i>=MAX_IND) perror("Error: maximum number of individuals exceeded");
   ++i;
}
fclose(fp);
sample_size=i;
printf("There are %d cases out of %d individuals\n",cases,sample_size);
if(n>0) printf("%d records with partial information have been left out \n",n);
return 0;
}

int plabel()
/*
permutes case-control labels
*/
{
int n,i,id,l,x[MAX_IND];
float r;

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
return 0;
}

int pgtype(int op)
/*
permutes genotypes
*/
{
int n,i,id,k,l,x[MAX_IND],gtype[MAX_IND];
float r;

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
     n-=1;
     r=rand()/(float)RAND_MAX;
     id=(int)floor(n*r+1);
     l=person[id].gtype[i];
     person[id].gtype[i]=person[n].gtype[i];
     person[n].gtype[i]=l;
   }
} break;
case 2:/*blockwise*/
  n=sample_size;
  for(l=0;l<n;l++) x[l]=person[l].no;
  while(n>0)
  {
     n-=1;
     r=rand()/(float)RAND_MAX;
     id=(int)floor(n*r+1);
     l=x[id];
     x[id]=x[n];
     x[n]=l;
  }
  for(i=0;i<nloci;++i)
  {
     if(sel[i]&&selp[i])
     {
       for(l=0;l<sample_size;++l) gtype[l]=person_s[l].gtype[i];
       for(l=0;l<sample_size;++l) person[l].gtype[i]=gtype[x[l]];
     }
  } break;
}
return 0;
}

int ind2eh(int m,int op)
/*
for PM, superseded in PM+
*/
{
int i,j,l,lo,lu,offset,attable,row,col;
int genotype[MAX_LOC];
struct {double lall,lcase,lctrl;} obs;

FILE *fp,*ef,*ca,*co;
char *ehfile="eh.dat";
char *cafile="case.dat";
char *cofile="control.dat";
char *scratch,line[2500];

switch(op)
{
case 0:lo=nall[m];lu=nall[0];break;
case 1:lo=np[m];lu=np[0];break;
case 2:lo=nnp[m];lu=nnp[0];break;
default:break;
}
if((scratch = tmpnam (NULL)) == NULL)
{
   perror ("Unable to create temporary filename");
   exit(EXIT_FAILURE);
}
/* printf ("Temporary filename \"%s\" created.\n", scratch); */
obs.lall=obs.lcase=obs.lctrl=0;
fp = fopen(scratch, "w+");
for (i=0;i<lu;++i) fwrite(&obs,sizeof(obs),1,fp);
for(i=0;i<sample_size;++i)
{
   p_t=person[i];
   l=0;
   switch(op)
   {
   case 0:
     for(j=0;j<nloci;++j)
     {
       if(!sel[j]) continue;
       genotype[l]=p_t.gtype[j];++l;
     } break;
   case 1:
     for(j=0;j<nloci;++j)
     {
        if(!sel[j]) continue;
        if(selp[j])
        {
          genotype[l]=p_t.gtype[j];++l;
        }
     } break;
   case 2:
     for(j=0;j<nloci;++j)
     {
        if(!sel[j]) continue;
        if(!selp[j])
        {
          genotype[l]=p_t.gtype[j];++l;
        }
     } break;
   }
 /*for(l=0;l<=m;++l) printf("%d ",genotype[l]);printf("\n");*/
   l=0;
   attable=position(m,genotype,op);
   if(attable>0)
   {
     offset=(attable-1)*(sizeof(obs));
     fseek(fp,offset,SEEK_SET);
     fread(&obs,sizeof(obs),1,fp);
     obs.lall++;
     if(p_t.affection) obs.lcase++;else obs.lctrl++;
     fseek(fp,offset,SEEK_SET);
     fwrite(&obs,sizeof(obs),1,fp);
   }
}
col=lo; row=lu/lo;

ef=fopen(ehfile,"w+");
if(ef==NULL) printf("Error opening %s \n",ehfile);
ca=fopen(cafile,"w+");
if(ca==NULL) printf("Error opening %s \n",cafile);
co=fopen(cofile,"w+");
if(co==NULL) printf("Error opening %s \n",cofile);
assert(ca && co && ef);
switch(op)
{
case 0:
   for(i=0;i<nloci;++i)
   if(sel[i])
   {
     fprintf(ef,"%3d ",alleles[i]);
     fprintf(ca,"%3d ",alleles[i]);
     fprintf(co,"%3d ",alleles[i]);
   } break;
case 1:
   for(i=0;i<nloci;++i)
   if(sel[i]&&selp[i])
   {
     fprintf(ef,"%3d ",alleles[i]);
     fprintf(ca,"%3d ",alleles[i]);
     fprintf(co,"%3d ",alleles[i]);
   } break;
case 2:
   for(i=0;i<nloci;++i)
   if(sel[i]&&!selp[i])
   {
     fprintf(ef,"%3d ",alleles[i]);
     fprintf(ca,"%3d ",alleles[i]);
     fprintf(co,"%3d ",alleles[i]);
   } break;
default:break;
}
fprintf(ef,"\n");
fprintf(ca,"\n");
fprintf(co,"\n");
for(i=0;i<row;++i){
   for(j=0;j<col;++j){
      attable=i*col+j;
      offset=attable*(sizeof(obs));
      fseek(fp,offset,SEEK_SET);
      fread(&obs,sizeof(obs),1,fp);
      fprintf(ef,"%5.0f",obs.lall);
      fprintf(ca,"%5.0f",obs.lcase);
      fprintf(co,"%5.0f",obs.lctrl);
   }
   fprintf(ef,"\n");
   fprintf(ca,"\n");
   fprintf(co,"\n");
}
fclose(fp);
fclose(ef);
fclose(ca);
fclose(co);
remove(scratch);
rows=row;cols=col;
if(cc||row>1) return 0;
addummy(ehfile);
return 0;
}

int a2g(int a1,int a2)
{
int lo,hi;
lo=a1; hi=a2;
if(a1>a2)
{
  lo=a2; hi=a1;
}
return (hi*(hi-1)/2+lo);
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

int position(int m,int *genotype,int op)

/*(loc1-1)*PROD n[2,...,m]+(loc2-1)*PROD n[3,...,m]+...+locm */

{
int l,pos,sum;
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

ind2ehs(int m,int op)
{
int i,j,k,l,lo,lu,attable,row,col,genotype[MAX_LOC];
link p,q;

FILE *ef,*ca,*co,*ehp;
char *ehfile="eh.dat",*ehpfile="ehplus.dat";
char *cafile="case.dat";
char *cofile="control.dat";
char line[2500];

switch(op)
{
case 0:lo=nall[m];lu=nall[0];break;
case 1:lo=np[m];lu=np[0];break;
case 2:lo=nnp[m];lu=nnp[0];break;
default:break;
}
valids=0;
p=NULL;
for(i=0;i<sample_size;++i)
{
   p_t=person[i];
   l=0;
   switch(op)
   {
   case 0:
     for(j=0;j<nloci;++j)
     {
       if(!sel[j]) continue;
       genotype[l]=p_t.gtype[j];++l;
     } break;
   case 1:
     for(j=0;j<nloci;++j)
     {
        if(!sel[j]) continue;
        if(selp[j])
        {
          genotype[l]=p_t.gtype[j];++l;
        }
     } break;
   case 2:
     for(j=0;j<nloci;++j)
     {
        if(!sel[j]) continue;
        if(!selp[j])
        {
          genotype[l]=p_t.gtype[j];++l;
        }
     } break;
   }
 /*for(l=0;l<=m;++l) printf("%d ",genotype[l]);printf("\n");*/
   l=0;
   attable=position(m,genotype,op);
   if(attable>0)
   {
     iscase=(boolean)(person[i].affection==true);
     q=insert(attable,p,op,p_t);
     p=q;
   }
}
used=(item*)calloc(valids,sizeof(item));
if(used==NULL) perror("memory allocation error");
set_all(p,op);
ehp=fopen(ehpfile,"w");
if(!ehp) {
  fprintf(stderr,"I can not open file %s\n",ehpfile);
  exit(1);
}
l=0;
for(i=0;i<nloci;i++) if(sel[i]) l++;
if(cc) {
  for(i=0;i<nloci;i++) if(sel[i]) fprintf(ehp,"%2d ",alleles[i]);
} else
switch(op) {
case 0: if(l==1) fprintf(ehp,"2 ");
   for(i=0;i<nloci;i++) {
     if(!sel[i]) continue;
     fprintf(ehp,"%2d ",alleles[i]);
   } break;
case 1: if(selectp==0) fprintf(ehp,"2 ");
   for(i=0;i<nloci;i++) {
     if(sel[i]&&selp[i]) fprintf(ehp,"%2d ",alleles[i]);
   } break;
case 2: if(selectn==0) fprintf(ehp,"2 ");
   for(i=0;i<nloci;i++) {
     if(sel[i]&&!selp[i]) fprintf(ehp,"%2d ",alleles[i]);
   } break;
default:break;
}
fprintf(ehp,"\n");
for(l=0;l<valids;++l)
{
  if(cc)
  fprintf(ehp,"%6d  %3d %3d %3d",used[l].id,
         used[l].cases+used[l].controls,used[l].cases,used[l].controls);
  else fprintf(ehp,"%6d  %3d",used[l].id,used[l].cases+used[l].controls);
  j=0;
  for(i=0;i<nloci;i++) switch(op) {
  case 0:
       if(!sel[i]) continue;
       fprintf(ehp," %2hd %2hd",used[l].l[j],used[l].u[j]);
       j++;
       break;
  case 1:
       if(!(sel[i]&&selp[i])) continue;
       if(selectp==0) fprintf(ehp," 1 1");
       fprintf(ehp," %2hd %2hd",used[l].l[j],used[l].u[j]);
       j++;
       break;
  case 2:
       if(!(sel[i]&&!selp[i])) continue;
       if(selectn==0) fprintf(ehp," 1 1");
       fprintf(ehp," %2hd %2hd",used[l].l[j],used[l].u[j]);
       j++;
       break;
  }
  fprintf(ehp,"\n");
}
fclose(ehp);
#ifdef useoldcode
col=lo;row=lu/lo;

ef=fopen(ehfile,"w+");
if(ef==NULL) printf("Error opening %s \n",ehfile);
ca=fopen(cafile,"w+");
if(ca==NULL) printf("Error opening %s \n",cafile);
co=fopen(cofile,"w+");
if(co==NULL) printf("Error opening %s \n",cofile);
assert(ca && co && ef);
switch(op)
{
case 0:
   for(i=0;i<nloci;++i)
   if(sel[i])
   {
     fprintf(ef,"%3d ",alleles[i]);
     fprintf(ca,"%3d ",alleles[i]);
     fprintf(co,"%3d ",alleles[i]);
   } break;
case 1:
   for(i=0;i<nloci;++i)
   if(sel[i]&&selp[i])
   {
     fprintf(ef,"%3d ",alleles[i]);
     fprintf(ca,"%3d ",alleles[i]);
     fprintf(co,"%3d ",alleles[i]);
   } break;
case 2:
   for(i=0;i<nloci;++i)
   if(sel[i]&&!selp[i])
   {
     fprintf(ef,"%3d ",alleles[i]);
     fprintf(ca,"%3d ",alleles[i]);
     fprintf(co,"%3d ",alleles[i]);
   } break;
default:break;
}
fprintf(ef,"\n");
fprintf(ca,"\n");
fprintf(co,"\n");
k=0;
for(i=0;i<row;++i){
   for(j=0;j<col;++j){
      attable=i*col+j+1;
      if(attable==used[k].id)
      {
        fprintf(ef,"%5d",used[k].cases+used[k].controls);
        fprintf(ca,"%5d",used[k].cases);
        fprintf(co,"%5d",used[k].controls);
        ++k;
      }
      else
      {
        fprintf(ef,"%5d",0);
        fprintf(ca,"%5d",0);
        fprintf(co,"%5d",0);
      }
   }
   fprintf(ef,"\n");
   fprintf(ca,"\n");
   fprintf(co,"\n");
}
fclose(ef);
fclose(ca);
fclose(co);
if(!cc) free(used);
rows=row;cols=col;
if(cc||row>1) return 0;
addummy(ehfile);
#endif
return 0;
}

link insert(int new_id,link p,int op,individual p_tt)
{
  link h,s,p1,p2;
  int i,k=0;
  p1=p;
  s=malloc(sizeof(item));
  s->id=new_id;
  if((p1==NULL)||(p1->id>s->id))
  {
    s->all=1;
    if(iscase) {s->cases=1;s->controls=0;}
    else {s->cases=0;s->controls=1;}
    for(i=0;i<nloci;i++)
    {
      s->l[i]=p_tt.locus[i][0];
      s->u[i]=p_tt.locus[i][1];
    }
    s->next=p1;
    p1=s;
    h=p1;
  }
  else
  {
    h=p1;
    while((p1!=NULL)&&(p1->id<=s->id))
    {
      if(p1->id==s->id)
      {
        p1->all++;k=1;
        if(iscase) p1->cases++;else p1->controls++;
      }
      p2=p1;
      p1=p1->next;
    }
    if(k==0)
    {
      s->next=p1;
      s->all=1;
      if(iscase) {s->cases=1; s->controls=0;}
      else {s->cases=0; s->controls=1;}
      for(i=0;i<nloci;i++)
      {
        s->l[i]=p_tt.locus[i][0];
        s->u[i]=p_tt.locus[i][1];
      }
      p2->next=s;
    }
    else free(s);
  }
  if(k==0) valids++;
  return h;
}

void set_all(link h,int op)
{
  int k,j,i=0;
  link p;
  while(h!=NULL)
  {
    used[i].id=h->id;
    used[i].cases=h->cases;
    used[i].controls=h->controls;
    j=0;
    for(k=0;k<nloci;k++) switch(op) {
    case 0:
         if(!sel[k]) continue;
         used[i].l[j]=h->l[k];
         used[i].u[j]=h->u[k];
         j++;
         break;
    case 1:
         if(!(sel[k]&&selp[k])) continue;
         used[i].l[j]=h->l[k];
         used[i].u[j]=h->u[k];
         j++;
         break;
    case 2:
         if(!(sel[k]&&!selp[k])) continue;
         used[i].l[j]=h->l[k];
         used[i].u[j]=h->u[k];
         j++;
         break;
    }
    ++i;
    p=h;
    h=h->next;
    if(p!=NULL)free(p);
  }
}

int addummy(char *filename)
/*
add dummy marker
*/
{
FILE *fp;
char line[2500];
int i,j,l;

fp=fopen(filename,"r+");
if(fp==NULL) printf("Error opening dummy file %s \n",filename);
fgets(line,2500,fp);sscanf(line,"%d",&l);
fgets(line,2500,fp);
fseek(fp,0,SEEK_SET);
fprintf(fp,"%3d %3d \n",2,l);
fprintf(fp,"%s",line);
for(i=0;i<2;++i)
{
  for(j=0;j<cols;++j) fprintf(fp,"%5.0f",0.0);
  fprintf(fp,"\n");
}
fclose(fp);
return 0;
}

int adjustn(int op,int op2)
/*
prepares appropriate counts and add dummy marker for EHPLUS
*/
{
FILE *ehp;
char *ehpfile="ehplus.dat";
int i,j,l;
ehp=fopen(ehpfile,"w");
if(!ehp) {
  fprintf(stderr,"I can not open file %s\n",ehpfile);
  exit(1);
}
l=0;
for(i=0;i<nloci;i++) if(sel[i]) l++;
if(l==1) fprintf(ehp,"2 ");
for(i=0;i<nloci;i++)
  if(sel[i]) fprintf(ehp,"%2d ",alleles[i]);fprintf(ehp,"\n");
for(i=0;i<valids;++i)
{
  switch(op) {
  case 0:
    fprintf(ehp,"%6d  %3d",used[i].id,used[i].controls);
    break;
  case 1:
    fprintf(ehp,"%6d  %3d",used[i].id,used[i].cases);
    break;
  case 2:
    fprintf(ehp,"%6d  %3d",used[i].id,used[i].cases+used[i].controls);
    break;
  }
  j=0;
  for(l=0;l<nloci;l++) switch(op2)
  {
    case 0:
         if(!sel[l]) continue;
         if(selected==0) fprintf(ehp," 1 1");
         fprintf(ehp," %2hd %2hd",used[i].l[j],used[i].u[j]);
         j++;
         break;
    case 1:
         if(!(sel[l]&&selp[l])) continue;
         if(selectp==0) fprintf(ehp," 1 1");
         fprintf(ehp," %2hd %2hd",used[i].l[j],used[i].u[j]);
         j++;
         break;
    case 2:
         if(!(sel[l]&&!selp[l])) continue;
         if(selectn==0) fprintf(ehp," 1 1");
         fprintf(ehp," %2hd %2hd",used[i].l[j],used[i].u[j]);
         j++;
         break;
  }
  fprintf(ehp,"\n");
}
fclose(ehp);
return 0;
}

int copyfile(char *from,char *to)
/*
keep input files from observed data
*/
{
FILE *fs,*ft;
char ch;

fs=fopen(from,"r");
if (!fs) {
  fprintf(stderr,"Cannot open %s",from);return 1;
}
ft=fopen(to,"w");
if (!ft) {
  fprintf(stderr,"Cannot open %s",to); return 1;
}
while ((ch=getc(fs))!=EOF) putc(ch,ft);
fclose(fs);
fclose(ft);

return 0;
}

int get_ln(float *l0, int *d0, float *l1, int *d1)
/*
gets loglikelihoods
*/
{
FILE *fs;
char line[241],*ptr;
char *h0 ="H0: No Association";
char *h1 ="H1: Allelic Associations Allowed";
char *scratch="scratch";
int exec,np0,np1;
float lnl0,lnl1;

fs=fopen(scratch,"w");
if(fs==NULL) perror("error openning scratch file in get_in");
if(useehplus) fprintf(fs,"n \nehplus.dat \nscratch.out \n");
else fprintf(fs,"n \neh.dat \nscratch.out \n");
fclose(fs);

remove("scratch.out");
if(calleh()) perror("Error executing EH/EH+ in get_ln");
fs=fopen("scratch.out","r");
if(fs==NULL) perror("Error openning file for get_ln");

np0=np1=0;
lnl0=lnl1=0;
while(fgets(line,240,fs))
{
   ptr=strstr(line,h0);
   if(ptr) sscanf(ptr+strlen(h0),"%d %f",&np0,&lnl0);
   *l0=lnl0;*d0=np0;
   ptr=strstr(line,h1);
   if(ptr) sscanf(ptr+strlen(h1),"%d %f",&np1,&lnl1);
   *l1=lnl1;*d1=np1;
}fclose(fs);
return 0;
}

int mfeh(float *x2, int *d)
/*
calculates model-free and other statistics using EH
*/
{
FILE *fo,*fs,*fgrid;
char cmdln[50],line[240];
int exec,mf;
int df,np0,np1,np2,mxstep;
float chisq,lnl0,lnl1,lnl2,chisqmx;
float chi1,chi2,chi3,chi4,chi5;
int df1,df2,df3,df4,df5;
int i,j,k,l,m,p;
float Kp,q,f0,f1,f2,ax,bx,cx;
int nlines;
float startp[MAXGRIDLINES][3],endp[MAXGRIDLINES][3],nevals[MAXGRIDLINES];
char *title="  #       q     f0     f1     f2      K  Chi-square  DF     P";

char *outfile="mfeh.out";
char *scratch="scratch";

fo=fopen(outfile,"w");
if(fo==NULL) perror("error openning output file");

q=freq;f0=pen0;f1=pen1;f2=pen2;
Kp=pow(1-q,2)*f0+2*q*(1-q)*f1+pow(q,2)*f2;
fs=fopen(scratch,"w");
if(fs==NULL) perror("error openning file for true model");
if(useehplus) fprintf(fs,"y \nehplus.dat \nscratch.out \n");
else fprintf(fs,"y \ncontrol.dat \ncase.dat \nscratch.out \n");
fprintf(fs,"%f \n%f %f %f \n",q,f0,f1,f2);
fclose(fs);
get_x2(&chi1,&df1,0);
/* under recessive model */
fs=fopen(scratch,"w");
if(fs==NULL) perror("error openning file for recessive model");
if(useehplus) fprintf(fs,"y \nehplus.dat \nscratch.out \n");
else fprintf(fs,"y \ncontrol.dat \ncase.dat \nscratch.out \n");
fprintf(fs,"%f \n%f %f %f \n",sqrt(Kp),0.0,0.0,1.0);
fclose(fs);
get_x2(&chi2,&df2,0);
/* under dominant model */
fs=fopen(scratch,"w");
if(fs==NULL) perror("error openning file for dominant model");
if(useehplus) fprintf(fs, "y \nehplus.dat \nscratch.out \n");
else fprintf(fs,"y \ncontrol.dat \ncase.dat \nscratch.out \n");
fprintf(fs,"%f \n%f %f %f \n",1-sqrt(1-Kp),0.0,1.0,1.0);
fclose(fs);
get_x2(&chi3,&df3,0);
fprintf(fo,"%s\n",title);
fgrid=fopen("grid.dat","r");
if(fgrid) {
   nlines=0;
   while(fgets(line,240,fgrid)){
         if(sscanf(line,"%f",&startp[nlines][0])!=1) continue;
         if(sscanf(line,"%f,%f,%f %f,%f,%f %f",
           &startp[nlines][0],&startp[nlines][1],&startp[nlines][2],
           &endp[nlines][0],  &endp[nlines][1],  &endp[nlines][2],
           &nevals[nlines])!=7)
           {printf("Bad grid line format: %s \n",line);return 1;}
           ++nlines;
   } fclose(fgrid);
} else {
        nlines=3;
        nevals[0]=nevals[2]=MAXEVALS;nevals[1]=1;
        startp[0][0]=0.000; startp[0][1]=0.000; startp[0][2]=1.000;
        startp[1][0]=startp[1][1]=startp[1][2]=Kp;
        startp[2][0]=0.000; startp[2][1]=1.000; startp[2][2]=1.000;
        endp[0][0]=endp[0][1]=endp[0][2]=Kp;
        endp[1][0]=endp[1][1]=endp[1][2]=Kp;
        endp[2][0]=endp[2][1]=endp[2][2]=Kp;
}
mxstep=0;
chisqmx=-1e20;
for(i=0,l=0;l<nlines;++l)
  for(p=0;p<nevals[l];++p)
  {
      f0=startp[l][0]+p/nevals[l]*(endp[l][0]-startp[l][0]);
      f1=startp[l][1]+p/nevals[l]*(endp[l][1]-startp[l][1]);
      f2=startp[l][2]+p/nevals[l]*(endp[l][2]-startp[l][2]);
      if(f0==f1 && f0==f2) q=0.5;
      else if(fabs((f0-Kp)/Kp)<QTOL) q=1;
      else if(fabs((f2-Kp)/Kp)<QTOL) q=0;
      else if(fabs(f0+f2-2*f1)<QTOL) q=0.5*(Kp-f0)/(f1-f0);
      else {
        ax=f0-2*f1+f2;bx=2*(f1-f0);cx=f0-Kp;
        q=(-bx+sqrt(bx*bx-4*ax*cx))/(2*ax);
      }
      fs=fopen(scratch,"w");
      if (fs==NULL) perror("error openning file for model-free");
      if (useehplus) fprintf(fs,"y \nehplus.dat \nscratch.out \n");
      else fprintf(fs,"y \ncontrol.dat \ncase.dat \nscratch.out\n");
      fprintf(fs,"%f \n%f %f %f\n",q,f0,f1,f2);
      fclose(fs);
      get_x2(&chisq,&df4,0);
      df4++;
      Kp=pow(1-q,2)*f0+2*q*(1-q)*f1+pow(q,2)*f2;
      fprintf(fo,"%3d %.5f %.4f %.4f %.4f %.4f %11.2f %3d %.4f\n",
                 i,q,f0,f1,f2,Kp,chisq,df4,probchi(chisq,df4));
      i++;
      if (chisq>chisqmx) {chisqmx=chisq;mxstep=i;}
  }
chi4=chisqmx;
/* under heterogeneity model*/
fs=fopen(scratch,"w");
if(fs==NULL) perror("error openning file for control");
if(useehplus) {
   fprintf(fs,"n \nehplus.dat \nscratch.out \n");
   adjustn(0,0);
}
else {
   fprintf(fs,"n \ncontrol.dat \nscratch.out \n");
   if(rows==1) addummy("control.dat");
}
fclose(fs);
get_x2(&lnl0,&np0,1);
fs=fopen(scratch,"w");
if(fs==NULL) perror("error openning file for case");
if(useehplus) {
   fprintf(fs,"n \nehplus.dat \nscratch.out \n");
   adjustn(1,0);
}
else {
   fprintf(fs,"n \ncase.dat \nscratch.out \n");
   if(rows==1) addummy("case.dat");
}
fclose(fs);
get_x2(&lnl1,&np1,1);
fs=fopen(scratch,"w");
if(fs==NULL) perror("error openning file for case/control");
if(useehplus) {
   fprintf(fs,"n \nehplus.dat \nscratch.out \n");
   adjustn(2,0);
}
else {
   fprintf(fs,"n \neh.dat \nscratch.out \n");
   if(rows==1) addummy("eh.dat");
}
fclose(fs);
get_x2(&lnl2,&np2,1);
chi5=-2*(lnl0+lnl1-lnl2);
if(usefehp) df5=np0+np1-np2;
else {
  if(rows==1) df5=(sqrt(8*cols+1)-3)/2;
  else df5=np0+np1-np2;
}
#ifdef DEBUG
fprintf(fo,"chi-square %9.2f %9.2f %9.2f %9.2f %9.2f\n",
        chi1,chi2,chi3,chi4,chi5);
#endif
fclose(fo);
x2[0]=chi1;d[0]=df1;
x2[1]=chi2;d[1]=df2;
x2[2]=chi3;d[2]=df3;
x2[3]=chi4;d[3]=df4;
x2[4]=chi5;d[4]=df5;
return 0;
}

int get_x2(float *x2, int *df, int mf)
/*
gets chi-square statistics from EH output
*/
{
FILE *fs;
char line[241],*ptr;
char *h0 ="H0: No Association";
char *h1,*h2;
int np0,np1,np2;
float lnl0,lnl1,lnl2;

remove("scratch.out");
if(calleh()) perror("Error executing EH/EH+");
fs=fopen("scratch.out","r");
if(fs==NULL) perror("error openning file for get_x2");

np0=np1=np2=0;
lnl0=lnl1=lnl2=0;
while(fgets(line,240,fs)) {
      ptr=strstr(line,h0);
      if(ptr)
      {
        if(usefehp) sscanf(ptr+strlen(h0),"%d %*f %*f %f",&np0,&lnl0);
        else sscanf(ptr+strlen(h0),"%d %f",&np0,&lnl0);
      }
      switch(mf){
      case 0:
           h1="H1: Markers Asso., Indep. of Disease";
           h2="H2: Markers and Disease Associated";
           ptr=strstr(line,h1);
           if(ptr)
           {
             if(usefehp) sscanf(ptr+strlen(h1),"%d %*f %*f %f",&np1,&lnl1);
             else sscanf(ptr+strlen(h1),"%d %f",&np1,&lnl1);
           }
           h1="H1: Markers and Disease Associated";
           ptr=strstr(line,h1);
           if(ptr)
           {
             if(usefehp) sscanf(ptr+strlen(h1),"%d %*f %*f %f",&np2,&lnl2);
             else sscanf(ptr+strlen(h1),"%d %f",&np2,&lnl2);
             np1=np0;lnl1=lnl0;
           }
           ptr=strstr(line,h2);
           if(ptr)
           {
             if(usefehp) sscanf(ptr+strlen(h2),"%d %*f %*f %f",&np2,&lnl2);
             else sscanf(ptr+strlen(h2),"%d %f",&np2,&lnl2);
           }
           *x2=-2*(lnl1-lnl2);*df=(np2-np1);
           break;
      case 1:
           h1="H1: Allelic Associations Allowed";
           ptr=strstr(line,h1);
           if(ptr) sscanf(ptr+strlen(h1),"%d %f",&np1,&lnl1);
           *x2=-lnl1;*df=np1;
           break;
      default: break;
      }
}fclose(fs);
return 0;
}

int calleh()
/*run EH*/
{
return((useehplus)?((usefehp)?system(cmd3):system(cmd2)):system(cmd1));
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

/*
Program for permutation and model-free association analysis

Questions and problems please refer to JH Zhao (IOP)
j.zhao@iop.kcl.ac.uk

History.

  7-8-1996 implement program to generate EH input from LINKAGE files
11-10-1996 model-free analysis for case-control data
 15-1-1998 PED2EH.CPP run under UNIX with better linkfile.x
 26-5-1998 multiple group permutation codes done
  3-6-1998 start to draft the combined program
 14-6-1998 change from unistd.h to stdio.h
 26-6-1998 fix various errors under Unix
 13-7-1998 avoid temporary file using linked list
16-10-1998 make program able to handle both EH and EHPLUS
 16-6-1999 add verbose explanation
  1-7-1999 order statistics according to paper, save EH inputs
 15-4-2002 change haplotype analyzer to be fastEH+
*/
