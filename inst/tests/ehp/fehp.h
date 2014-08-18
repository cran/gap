#ifndef FASTEH_H
#define FASTEH_H

#include <stdio.h>
#include <string.h>
#include <math.h>

#define maxalleles    50
#define maxloci       30
#define maxobscom   2000
#define maxloops     100
#define eps          0.1e-10
#define exp2(h)      (pow(2,h))

#define maxalle    maxalleles

#define v               "1.11"
#define FNL             120
#define LL              1500
char infile_name[FNL],wf_NAME[FNL],outfile_name[FNL];

typedef enum         boolean {false,true} boolean;
boolean delta=true,correct_den=false;
boolean cc,err,onemark,wd;

char ch;
FILE *finp,*fout;
int sample_size,obscom,totalloci,haplnum,maxall,loop,loop1;
int locihn[maxloci+1],locigeno[maxloci+1],loci[maxloci+1];
short locip[maxloci+1];
double *indivN,*caseN,*contrlN,*dg11,*dg12,*dg22,p[maxloci][maxalle];
double *oldhaplo,*newhaplo,*inihaplo,*indhaplo;
double oldloglike,loglike,iniloglike,diff,indloglike,oldlog,cf;
double totalall,totalind,ff[3],ss[3],casep[3],controlp[3],s2,pp,qq;
long int obshap;

typedef struct
{
  short l[maxloci],u[maxloci];
} phenotype;
phenotype *alist;

struct DBT {
  boolean dma;
  int time,nloci,fl,ll,l,u;
};

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
boolean ini;
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
int line, i, n, times;
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
}

int chklimit(void)
/*Check default upper bounds*/
{
int i,lh,lg;
if(totalloci>maxloci) fprintf(stderr,"set maxloci=%d\n",totalloci);
lh=lg=1;
for(i=0;i<totalloci;i++) {
  lh*=loci[i];lg*=loci[i]*(loci[i]+1)/2;
}
if(cc) lh*=2;
if(maxall>maxalle) fprintf(stderr,"set maxalle=%d\n",maxall);
/*maxhap is obsolete now*/
if(obscom>maxobscom) fprintf(stderr,"set maxobscom=%d\n",obscom);
return 0;
}

void outH(int *locia, double *a)
/*output haplotype frequencies and likelihoods*/
{
int times,i,j,line,n,malleles[maxloci],cn0=0,cn1=0,cn2=0;
double oi,ei,ftdev,x2ft,x2ftdf,x2,x2df,lrt,lrtdf;

/*correct missing alleles here*/
for(i=0;i<totalloci;i++) {
  n=0;
  for(j=1;j<loci[i];j++) if(p[i][j]<eps) n++;
  malleles[i]=n;
}
if(cc) {
  n=0;
  for(i=1;i<totalloci;i++) n+=(loci[i]-malleles[i])-1;
  cn0=n;
  n=1;
  for(i=1;i<totalloci;i++) n*=(loci[i]-malleles[i]);
  n--;
  cn1=n;
  cn2=2*n;
  if(onemark) cn1=2*n;
} else {
  n=0;
  for(i=0;i<totalloci;i++) n+=(loci[i]-malleles[i])-1;
  cn0=n;
  n=1;
  for(i=0;i<totalloci;i++) n*=(loci[i]-malleles[i]);
  n--;
  cn1=n;
}
times=n=1;
for(i=0;i<totalloci;i++){
  n*=locia[i];
  locip[i]=1;
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
x2ft=0;
x2ftdf=0;
x2=0;
x2df=0;
lrt=0;
lrtdf=0;
do {
  line=linenum(locia,locip,1);
  for(i=1;i<=totalloci;i++) {
    if(cc&&i==1) {
      if(locip[i-1]==1) fprintf(fout, "%4c     ", '+');else fprintf(fout, "%4c     ", 'D');
    } else fprintf(fout, "%4d     ", locip[i-1]);
  }
  if(cc&&!onemark) {
    fprintf(fout, "      %8.6f      %8.6f     %8.6f",
            inihaplo[line-1],indhaplo[line-1],a[line-1]);
  } else {
    fprintf(fout, "      %8.6f     %8.6f",inihaplo[line-1],a[line-1]);
    if(inihaplo[line-1]>0) {
      if(a[line-1]>0) {
        lrt+=a[line-1]*log(a[line-1]/inihaplo[line-1]);
        lrtdf++;
      }
      x2+=pow(a[line-1]-inihaplo[line-1],2)/inihaplo[line-1];
      x2df++;
    }
  }
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
      x2ft+=ftdev*ftdev;
      if(inihaplo[line-1]>0) x2ftdf++;
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
x2*=s2;
lrt*=s2;
for(i=0;i<totalloci;i++) fprintf(fout, "---------");
if(cc && !onemark) fprintf(fout, "-----------------------------------------\n");
else fprintf(fout, "---------------------------------\n");
fprintf(fout, "# of Iterations = %d\n\n", loop);
fprintf(fout,"                                    #parm   Ln(L)     Chi-square");
if (cc) fprintf(fout,"    Conditional Ln(L)\n");else fprintf(fout,"\n");
fprintf(fout,"-------------------------------------------------------------------\n");
if(cc) {
  n=0;
  for(i=1;i<totalloci;i++) n+=loci[i]-1;
  fprintf(fout,"H0: No Association                     %2d  %8.2f      0.00",cn0,iniloglike);
  fprintf(fout,"          %8.2f\n",iniloglike-cf);
  n=1;
  for(i=1;i<totalloci;i++) n*=loci[i];
  n--;
  if(!onemark) {
    if(fabs(indloglike-iniloglike)<eps)
      fprintf(fout,"H1: Markers Asso., Indep. of Disease   %2d  %8.2f      0.00",cn1,indloglike);
    else fprintf(fout,"H1: Markers Asso., Indep. of Disease   %2d  %8.2f %8.2f",
              cn1,indloglike,2*(indloglike-iniloglike));
    fprintf(fout,"          %8.2f\n",indloglike-cf);
    if(fabs(loglike-iniloglike)<eps)
      fprintf(fout,"H2: Markers and Disease Associated     %2d  %8.2f      0.00",cn2,loglike);
    else fprintf(fout,"H2: Markers and Disease Associated     %2d  %8.2f %8.2f",
              cn2,loglike,2*(loglike-iniloglike));
    fprintf(fout,"          %8.2f\n",loglike-cf);
  } else {
    if(fabs(loglike-iniloglike)<eps)
      fprintf(fout,"H1: Markers and Disease Associated     %2d  %8.2f      0.00",cn1,loglike);
    else fprintf(fout,"H1: Markers and Disease Associated     %2d  %8.2f %8.2f",
              cn1,loglike,2*(loglike-iniloglike));
    fprintf(fout,"          %8.2f\n",loglike-cf);
  }
} else {
  n=0;
  for(i=0;i<totalloci;i++) n+=loci[i]-1;
  fprintf(fout,"H0: No Association                     %2d  %8.2f      0.00\n",cn0,iniloglike);
  n=1;
  for(i=0;i<totalloci;i++) n*=loci[i];
  n--;
  fprintf(fout, "H1: Allelic Associations Allowed       %2d  %8.2f",cn1, loglike);
  if(fabs(loglike-iniloglike)<eps) fprintf(fout, "      0.00\n");
  else fprintf(fout, "%10.2f\n", 2*(loglike-iniloglike));
  fprintf(fout,"\nNonzero haplotypes h[indep] (excluding missing alleles) = %.0lf\n",x2df);
  fprintf(fout,"Nonzero haplotypes h[w/assoc] (excluding missing alleles/haplotypes) = %.0lf\n",lrtdf);
  fprintf(fout,"\nChi-squares based directly on haplotype frequencies (not recommended)\n\n");
  fprintf(fout,"1. Pearson                %.2lf\n",x2);
  fprintf(fout,"2. Log-likelihood ratio   %.2lf\n",2*lrt);
  if(delta) fprintf(fout,"3. Freeman-Tukey          %.2lf\n",x2ft);
  fprintf(fout,"\nChi-square 1 = %.0lf * sum_{i=1}^%d x_i\n\n",s2,n+1);
  fprintf(fout,"1. x_i=(h_i[w/assoc]-h_i[indep])^2/h_i[indep], h_i[indep]>0\n");
  fprintf(fout,"\nChi-square 2 = %.0lf * sum_{i=1}^%d chisq_i\n\n",2*s2,n+1);
  fprintf(fout,"2. x_i=h_i[w/assoc]*ln(h_i[w/assoc]/h_i[indep]), h_i[indep], h_i[w/assoc]>0\n");
  if(delta) {
    fprintf(fout,"\nChi-square 3 = sum_{i=1}^%d x_i*x_i\n\n",n+1);
    fprintf(fout,"3. x_i=sqrt(h_i[w/assoc]*%.0lf)",s2);
    fprintf(fout,"+sqrt(h_i[w/assoc]*%.0lf+1)\n",s2);
    fprintf(fout,"          -sqrt(4*h_i[indep]*%.0lf+1)\n",s2);
  }
}
putc('\n', fout);
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

void logo()
{
int i;
printf("\033[2J\n");for(i=0;i<5;i++) putchar('\n');
printf("          +------------------------------------------+\n");
printf("          |                                          \n");
printf("          |   EH with allele caching version %s    |\n",v);
printf("          |                                          |\n");
printf("          +------------------------------------------+\n\n");
printf("          Based on EH by J. Ott and Xiaoli Xie 1990.5\n\n");
}

char P_peek(FILE *fp)
{
char s;
s=getc(fp);ungetc(s,fp);return(s);
}

boolean P_eoln(FILE *fp)
{
char s;
s=getc(fp);ungetc(s,fp);
if(feof(fp)) return(feof(fp));return((boolean)(s=='\n'));
}

boolean seekeof(FILE **f)
{
while((P_peek(*f)==' '||P_peek(*f)== '\t')&&!feof(*f)) getc(*f);
return feof(*f);
}

boolean seekeoln(FILE **f)
{
while((P_peek(*f)==' '||P_peek(*f)=='\t')&&!P_eoln(*f)) if(!feof(*f)) getc(*f);
return (P_eoln(*f)||feof(*f));
}
void programsuccess()
{
logo();printf("        Running program is completed successfully\n");
}

void programaborted()
{
logo();printf("        Running program is aborted\n");
}

double sumcase, sumcontrol;

void getfreq(FILE *fi)
{
double TEMP;

printf("Enter frequency of disease allele\n");
fscanf(fi,"%lg%*[^\n]", &pp);
printf("Enter penetrances for each genotype in the following order\n");
printf(" Genotype   +/+   +/D    D/D\n");
fscanf(fi,"%lg%lg%lg%*[^\n]",&ff[0],&ff[1],&ff[2]);
TEMP=1-pp;
sumcase=ff[0]*(TEMP*TEMP)+2*ff[1]*pp*(1-pp)+ff[2]*pp*pp;
sumcontrol=(1-ff[0])*(TEMP*TEMP)+2*(1-ff[1])*pp*(1-pp)+(1-ff[2])*pp*pp;
casep[0]=ff[0]*(TEMP*TEMP)/sumcase;
casep[1]=2*ff[1]*pp*(1-pp)/sumcase;
casep[2]=ff[2]*pp*pp/sumcase;
controlp[0]=(1-ff[0])*(TEMP*TEMP)/sumcontrol;
controlp[1]=2*(1-ff[1])*pp*(1-pp)/sumcontrol;
controlp[2]=(1-ff[2])*pp*pp/sumcontrol;
}

void getoptions()
{
FILE *fi,*infile,*wf;
char s[FNL];
boolean file_exist=false;

wf=infile=NULL;
logo();
ch='n';
fi=stdin;
printf("Do you wish to use the case-control sampling option?  [N]\n");
if(seekeoln(&fi)) {
  fscanf(fi,"%*[^\n]");
  fgetc(fi);
} else {
  fscanf(fi,"%c%*[^\n]", &ch);
  fgetc(fi);
}
cc=(boolean)(ch=='Y'||ch=='y');
do {
  err=false;
  strcpy(outfile_name, "ehplus.dat");
  printf("Enter name of data file  [EHPLUS.DAT]\n");
  if(seekeoln(&fi)) fgetc(fi);
  else {fgets(s,FNL,fi);
  sscanf(s,"%s",&outfile_name);}
  printf("you entered: %s\n", outfile_name);
  strcpy(infile_name, outfile_name);
  infile=fopen(infile_name, "r");
  if(!infile) err=true;
} while(err);
printf("Enter name of output file.  [EHPLUS.OUT]\n");
strcpy(outfile_name, "ehplus.out");
if(seekeoln(&fi)) fgetc(fi);
else {fgets(s,FNL,fi);
sscanf(s,"%s",&outfile_name);}
strcpy(wf_NAME, outfile_name);
wf=fopen(wf_NAME,"r");
if(wf) {
  file_exist=true;
  fclose(wf);
}
if(file_exist) {
  printf("File %s exist. Over write it? [y/n]\n", outfile_name);
  if(seekeoln(&fi)) {
    fscanf(fi,"%*[^\n]");
    fgetc(fi);
    strcpy(wf_NAME, outfile_name);
    wf=fopen(wf_NAME,"w");
  } else {
    fscanf(fi,"%c%*[^\n]", &ch);
    fgetc(fi);
    if(ch!='Y'&&ch!='y') {
      programaborted();
      exit(1);
    }
    strcpy(wf_NAME, outfile_name);
    wf=fopen(wf_NAME,"w");
  }
} else {
  strcpy(wf_NAME, outfile_name);
  wf=fopen(wf_NAME,"w");
}
if(!wf) {
  fprintf(stderr,"I can not open file %s\n",wf_NAME);
}
if(cc) getfreq(fi);
logo();
printf("        Program is running\n\n");

}
#endif
