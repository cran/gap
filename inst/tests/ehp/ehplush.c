#include <stdio.h>
#include <string.h>
#include <math.h>

#define maxalleles       15 /*25*/
#define maxloci           6 /*6*/
#define maxhap          286 /*2*25*10*15*/
#define maxposcom       400 /*85995L*/
#define maxobscom       400 /*1500*/
#define maxalle         maxalleles
#define v               "1.11"
#define bs              '\b'
#define verysmall       0.1e-10
#define FNL             120
#define exp2(h) (pow(2,h))
typedef enum boolean {false, true} boolean;
typedef int larr[maxloci+1];
typedef double linetype[maxhap];

boolean ehplus=true ,delta=true;
/* The 10-10000000th prime number table
 10  29  100  541      1000      7919
 20  71  150  863      5000     48611
 30 113  200 1223     10000    104729
 40 173  250 1583     50000    611953
 50 229  300 1987    100000   1299709
 60 281  400 2471    500000   7368787
 70 349  500 3571   1000000  15485863
 80 409  600 4409  10000000 179424673
 90 463  700 5279
         800 6133
         900 6997
Case-control data has maxbinsize=3 maxbinsize(obscom)
Change hash(.) for alternative hash function
*/
#define PRIME 541
#define hash(x)(x%PRIME)
#define maxbinsize 12
struct {int n,obs[maxbinsize];long id[maxbinsize];} bin[PRIME];
int bindex,mmid;
long id[3*maxobscom],idsave[maxobscom];
double p[maxloci][maxalle],*indivN,*caseN,*contrlN,*dg11,*dg12,*dg22,*dg33;
double totalg,oldloglike,loglike,iniloglike,diff,indloglike,oldlog;
double pp,totalall,totalind,genefreq,ttemp,ff[3],casep[3],controlp[3],s2;
FILE *controlf,*casef,*wf,*tmp;
char fname[FNL],ch;
char controlf_NAME[FNL],casef_NAME[FNL],wf_NAME[FNL],tmp_NAME[FNL];
int sample_size,obscom,i,totalloci,cnloci,haplnum,genonum,maxall,loop,loop1;
boolean casestudy,err,onemark,genotype;
larr locihn,locigeno,locik,loci,caseloci,locip,lociq;
linetype oldhaplo,newhaplo,inihaplo,indephaplo;

void logo()
{
int i;
printf("\033[2J\n");for(i=0;i<7;i++) putchar('\n');
printf("              Program  EH  version %s      \n\n\n\n", v);
printf("        Developed under a grant from NIMH to J. Ott\n");
printf("        Programmed by Xiaoli Xie        August 1995\n\n");
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

void initial(double (*p)[maxalle])
{
int i,j;
for(i=0;i<totalloci;i++) for(j=0;j<maxall;j++) p[i][j] = 0.0;
totalind=0.0;
}

void programsuccess()
{
logo();printf("        Running program EH is completed successfully\n");
}

void programaborted()
{
logo();printf("        Running  program EH is aborted\007\n");
}

void getfreq(FILE *fi)
{
double sumcase, sumcontrol, TEMP;
printf("Enter gene frequency of disease allele\n");
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

void outind(FILE **ff, int *locia, double *a, int sp)
{
int i,j,n;
n=1;
for(i=sp-1;i<totalloci;i++) n*=locia[i];
j=locia[totalloci-1];
for(i=1;i<=n;i++) {
  fprintf(*ff, "%5.2f ", a[i-1]);if(i%j==0) putc('\n',*ff);
}
}

void casecontrol(f1, f2, l, k, totalloci, line)
FILE **f1, **f2;
int *l, k, totalloci, *line;
{
int i,j;
double ca,co;
for (i=1;i<=l[k-1];i++){
  for (j=1;j<=i;j++){
    if(k==totalloci){
      if(!seekeoln(f1))fscanf(*f1,"%lg",&co);
      if(!seekeoln(f2))fscanf(*f2,"%lg",&ca);
      contrlN[*line]=co;
      caseN[*line]=ca;
      dg33[*line]=co+ca;
      dg11[*line]=casep[0]*ca+controlp[0]*co;
      dg12[*line]=casep[1]*ca+controlp[1]*co;
      dg22[*line]=casep[2]*ca+controlp[2]*co;
      (*line)++;
    } else casecontrol(f1,f2,l,k+1,totalloci,line);
  }
}
if(k!=totalloci) return;
fscanf(*f1, "%*[^\n]");getc(*f1);
fscanf(*f2, "%*[^\n]");getc(*f2);
}

void initdg()
{
FILE *tmp;
int i;
char tmp_NAME[FNL];
strcpy(tmp_NAME,"temp.dat");
tmp=fopen(tmp_NAME,"w");
fprintf(tmp, "2 ");
for(i=0;i<totalloci;i++) fprintf(tmp," %2d",loci[i]);putc('\n', tmp);
outind(&tmp,locigeno,dg11,1);
outind(&tmp,locigeno,dg12,1);
outind(&tmp,locigeno,dg22,1);
if(tmp!=NULL) fclose(tmp);
}

double getN(f, l, k, totalloci, p, line)
FILE **f;
int *l, k, totalloci;
double (*p)[maxalle];
int *line;
{
int i,j;
double n,tempn,temp;
tempn=0.0;
for(i=0;i<l[k-1];i++){
  for(j=0;j<i+1;j++){
    if(k==totalloci){
      if(!seekeoln(f))fscanf(*f,"%lg",&n);
      indivN[*line]=n;
      p[k-1][i]+=n;
      p[k-1][j]+=n;
      tempn+=n;
      (*line)++;
    } else{
      temp=getN(f,l,k+1,totalloci,p,line);
      p[k-1][i]+=temp;
      p[k-1][j]+=temp;
      tempn+=temp;
    }
  }
}
if(k==totalloci) {
  fscanf(*f, "%*[^\n]");getc(*f);
}
return tempn;
}

void calcula(int totalloci)
{
int i,j;
for(i=0;i<totalloci;i++) for(j=0;j<maxall;j++) p[i][j]/=totalall;
}

int linenum(int *locia, int *locip, int sp)
{
int sum,j;
sum=0;
for(j=sp;j<=totalloci;j++)
  if(j==totalloci) sum+=locip[j-1]; else sum=(sum+locip[j-1]-1)*locia[j];
bindex=hash(sum);
if(ehplus&&genotype) {
  j=bsch(sum,bin[bindex].n)+1;
  /*
  if(j!=-1) printf("j=%d id=%ld sum=%ld case=%lf\n",j,bin[bindex].id[mmid],sum,indivN[j-1]);
  */
  return(j);
}
else return(sum);
}

int bsch(int k, int n)
{
int j,lo,hi,mid;
lo=0;hi=n;
while(lo<=hi)
{
  mmid=mid=(lo+hi)/2;
  if(bin[bindex].id[mid]<k) lo=mid+1;
  else if(bin[bindex].id[mid]>k) hi=mid-1; else return(bin[bindex].obs[mid]);
}
return -2;
}

void normalizeH(int *locia, double *a)
{
double sump,sumq,asump,asumq;
int i,n;
sump=sumq=asump=asumq=0.0;
n=1;
for(i=0;i<totalloci;i++) n*=locia[i];
for(i=0;i<n/2;i++){
  sump+=inihaplo[i];asump+=a[i];
}
for(i=n/2;i<n;i++){
  sumq+=inihaplo[i];asumq+=a[i];
}
for(i=0;i<n/2;i++){
  inihaplo[i]=inihaplo[i]*(1-pp)/sump;
  a[i]*=(1-pp)/asump;
}
for(i=n/2;i<n;i++){
  inihaplo[i]=inihaplo[i]*pp/sumq;
  a[i]*=pp/asumq;
}
}

void outvec(int *locia, double *a)
{
int times,i,line,n;
double oi,ei,ftdev;
times=n=1;
for(i=0;i<totalloci;i++){
  n*=locia[i];locip[i]=1;
}
fprintf(wf, "There are %d Possible Haplotypes of These %d Loci.\n",n, totalloci);
fprintf(wf, "They are Listed Below, with their Estimated Frequencies:\n");putc('-', wf);
for(i=0;i<totalloci;i++) fprintf(wf,"---------");
if(casestudy&&!onemark)
{
for(i=0;i<41;i++) putc('-',wf);fprintf(wf, "\n|");
for(i=0;i<totalloci;i++) fprintf(wf, " Allele  ");
fprintf(wf, "|         Haplotype   Frequency         |\n");putc('|',wf);
for(i=0;i<totalloci;i++) fprintf(wf, "   at    ");
fprintf(wf, "|                                       |\n");putc('|',wf);
}
else
{
for(i=0;i<33;i++) putc('-',wf);fprintf(wf, "\n|");
for(i=0;i<totalloci;i++) fprintf(wf, " Allele  ");
fprintf(wf, "|      Haplotype Frequency      |\n");putc('|', wf);
for(i=0;i<totalloci;i++) fprintf(wf, "   at    ");
fprintf(wf, "|                               |\n");putc('|', wf);
}
if(casestudy) {
  fprintf(wf, " Disease ");
  for(i=1;i<totalloci;i++)
    if(totalloci>9) fprintf(wf, " Marker%2d", i);else fprintf(wf, " Marker%d ", i);
} else for(i=1;i<=totalloci;i++)
  if(totalloci>9) fprintf(wf, " Locus %2d", i);else fprintf(wf, " Locus %d ", i);
if(casestudy&&!onemark) fprintf(wf, "|  Independent   Ind-Disease   w/Asso.  |\n");
else fprintf(wf, "|  Independent   w/Association  |\n");putc('-', wf);
for(i=0;i<totalloci;i++) fprintf(wf, "---------");
if(casestudy&&!onemark) for(i=0;i<41;i++) putc('-', wf);
else for(i=0;i<33;i++) putc('-', wf);
if(!delta) putc('\n', wf);
else fprintf(wf,"     Observed     Expected    Freeman-Tukey Zi\n");
do {
  line=linenum(locia,locip,1);
  for(i=1;i<=totalloci;i++) {
    if(casestudy&&i==1) {
      if(locip[i-1]==1) fprintf(wf, "%4c     ", '+');else fprintf(wf, "%4c     ", 'D');
    } else fprintf(wf, "%4d     ", locip[i-1]);
  }
  if(casestudy&&!onemark) fprintf(wf, "      %8.6f      %8.6f     %8.6f",
            inihaplo[line-1],indephaplo[line-1],a[line-1]);
  else fprintf(wf, "      %8.6f     %8.6f",inihaplo[line-1],a[line-1]);
  if(!delta) fprintf(wf,"\n");
  else {
    if(casestudy&&!onemark) {
      ei=indephaplo[line-1]*s2;
      oi=a[line-1]*s2;
      ftdev=sqrt(oi)+sqrt(oi+1)-sqrt(4.0*ei+1);
      fprintf(wf, "     %8.2f      %8.2f     %8.2f",oi,ei,ftdev);
    }
    else {
      ei=inihaplo[line-1]*s2;
      oi=a[line-1]*s2;
      ftdev=sqrt(oi)+sqrt(oi+1)-sqrt(4.0*ei+1);
      fprintf(wf, "           %8.2f     %8.2f     %8.2f",oi,ei,ftdev);
    }
    if(fabs(ftdev)<5.0) fprintf(wf,"\n");else fprintf(wf," *\n");
  }
  locip[totalloci-1]++;
  for(i=totalloci;i>=2;i--)
    if(locip[i-1]>loci[i-1]){
      locip[i-1]=1;
      locip[i-2]++;
    }
  times++;
} while (times<=n);putc('-', wf);
for(i=0;i<totalloci;i++) fprintf(wf, "---------");
if(casestudy && !onemark) fprintf(wf, "-----------------------------------------\n");
else fprintf(wf, "---------------------------------\n");
fprintf(wf, "# of Iterations = %d\n\n", loop);
fprintf(wf,"                                       df   Ln(L)     Chi-square\n");
fprintf(wf,"-------------------------------------------------------------------\n");
if(casestudy) {
  n=0;
  for(i=1;i<totalloci;i++) n+=loci[i]-1;
  fprintf(wf,"H0: No Association                     %2d  %8.2f      0.00\n",n,iniloglike);
  n=1;
  for(i=1;i<totalloci;i++) n*=loci[i];
  n--;
  if(!onemark) {
    if(fabs(indloglike-iniloglike)<verysmall)
      fprintf(wf,"H1: Markers Asso., Indep. of Disease   %2d  %8.2f      0.00\n",n,indloglike);
    else fprintf(wf,"H1: Markers Asso., Indep. of Disease   %2d  %8.2f %8.2f\n",
              n,indloglike,2*(indloglike-iniloglike));
    if(fabs(loglike-iniloglike)<verysmall)
      fprintf(wf,"H2: Markers and Disease Associated     %2d  %8.2f      0.00\n",n*2,loglike);
    else fprintf(wf,"H2: Markers and Disease Associated     %2d  %8.2f %8.2f\n",
              n*2,loglike,2*(loglike-iniloglike));
  } else {
    if(fabs(loglike-iniloglike)<verysmall)
      fprintf(wf,"H1: Markers and Disease Associated     %2d  %8.2f      0.00\n",n*2,loglike);
    else fprintf(wf,"H1: Markers and Disease Associated     %2d  %8.2f %8.2f\n",
              n*2,loglike,2*(loglike-iniloglike));
  }
} else {
  n=0;
  for(i=0;i<totalloci;i++) n+=loci[i]-1;
  fprintf(wf,"H0: No Association                     %2d  %8.2f      0.00\n",n,iniloglike);
  n=1;
  for(i=0;i<totalloci;i++) n*=loci[i];
  n--;
  fprintf(wf, "H1: Allelic Associations Allowed       %2d  %8.2f",n, loglike);
  if(fabs(loglike-iniloglike)<verysmall) fprintf(wf, "      0.00\n");
  else fprintf(wf, "%10.2f\n", 2*(loglike-iniloglike));
} putc('\n', wf);
}

void outindP(int totalloci)
{
int i,j;
fprintf(wf, "Estimates of Gene Frequencies (Assuming Independence)\n");
if(casestudy) fprintf(wf, "(Disease gene frequencies are user specified)\n");
fprintf(wf, "----\\------------");
for(i=1; i <=maxall; i++) fprintf(wf, "--------");
fprintf(wf, "\nlocus \\ allele ");
for(i=1; i <=maxall; i++) fprintf(wf, "%3c%2d%3c", ' ', i, ' ');
fprintf(wf, "\n--------\\--------");
for(i=1; i <=maxall; i++) fprintf(wf, "--------");
putc('\n', wf);
if(casestudy) {
  fprintf(wf, "Disease |%7c", ' ');
  for(j=0;j<maxall; j++)
    if(p[0][j]<verysmall) fprintf(wf, "%8c", ' '); else fprintf(wf, "%8.4f", p[0][j]);
  putc('\n', wf);
  for(i=1;i<totalloci;i++) {
    fprintf(wf, "%4d    |%7c", i, ' ');
    for(j=0;j<maxall;j++)
      if(p[i][j]<verysmall) fprintf(wf, "%8c", ' ');else fprintf(wf, "%8.4f", p[i][j]);
    putc('\n', wf);
  }
} else {
  for(i=0;i<totalloci;i++) {
    fprintf(wf, "%4d    |%7c", i+1, ' ');
    for(j=0;j<maxall;j++)
      if(p[i][j]<verysmall) fprintf(wf, "%8c", ' ');else fprintf(wf, "%8.4f", p[i][j]);
    putc('\n', wf);
  }
}
fprintf(wf, "-----------------");
for(i=0;i<maxall;i++) fprintf(wf, "--------");
fprintf(wf, "\n# of Typed Individuals: %2ld\n\n",(long)floor(totalind + 0.5));
}

void getlocih(int totalloci)
{
int i;
for(i=0;i<totalloci;i++) locigeno[i]=loci[i]*(loci[i]+1)/2;
}

void gethapP(int totalloci)
{
int line, i, n, times;
n=times=1;
for(i=0;i<totalloci;i++){
  locip[i]=1;
  n*=loci[i];
}
do {
  line=linenum(loci,locip,1);
  newhaplo[line-1]=1.0;
  for(i=0;i<totalloci;i++) newhaplo[line-1]*=p[i][locip[i]-1];
  locip[totalloci-1]++;
  for(i=totalloci;i>=2;i--) {
    if(locip[i-1]>loci[i-1]){
      locip[i-1]=1;
      locip[i-2]++;
    }
  }
  times++;
} while(times<=n);
for(i=0;i<n;i++) inihaplo[i]=newhaplo[i];
}

void getfirst(int *first, int *last)
{
int j;
boolean ini;
*first=*last=1;
ini=true;
for(j=1;j<=totalloci;j++)
  if(locihn[j-1]==1) {
    if(ini) if(*first<=j) {
        *first=j;
        ini=false;
      }
    if(j>*last) *last=j;
  }
}

void swap__(int *a, int *b)
{
int temp;
temp = *a;*a = *b;*b = temp;
}

void newtoold()
{
int i;
for (i = 0; i < haplnum; i++) {
  oldhaplo[i] = newhaplo[i];newhaplo[i] = 0.0;
}
}

void iniindex()
{
int i;
for(i=0;i<totalloci;i++) {
  locip[i]=lociq[i]=locik[i]=1; locihn[i]=0;
}
}

void inivar(int nloci, int *locii)
{
int i;
haplnum=genonum=1;
for(i=0;i<nloci;i++){
  haplnum*=locii[i];
  locigeno[i]=locii[i]*(locii[i]+1)/2;
  genonum*=locigeno[i];
}
}

boolean inputdataok(FILE **inf,FILE **outf,int *nloci,int *locii)
{
double data;
boolean b;
*nloci=0;
maxall=0;
b=false;
while(!seekeof(inf)) {
  if(!b) {
    while(!seekeoln(inf)) {
      (*nloci)++;
      fscanf(*inf, "%hd", &locii[*nloci-1]);
      if(locii[*nloci-1]>maxall) maxall=locii[*nloci-1];
    }
    fscanf(*inf, "%*[^\n]");getc(*inf);
    inivar(*nloci, locii);
    putc('\n', *outf);
    b=true;
    i=0;
    continue;
  }
  while(!seekeoln(inf)) {
    fscanf(*inf, "%lg", &data);
    i++;
  }
  fscanf(*inf, "%*[^\n]");getc(*inf);
}
return(i==genonum);
}

boolean dataconsist(int *loci1,int *loci2, int n1, int n2)
{
int i;
boolean equal;
i=0;
equal=(n1==n2);
while(equal&&i<n1){
  equal=(loci1[i]==loci2[i]);
  i++;
}
return equal;
}

struct LOC_likelihood {
  linetype hap;
  int totalloci;
  boolean dma;
  int first, last, h, linek, lineq, linep;
  double multilog, temp;
} ;

int multilike(int i, double *logtemp, struct LOC_likelihood *LINK)
{
int j;
if(locihn[i-1]!=1) {
  multilike(i+1,logtemp,LINK);return;
}
if(i==LINK->first) {
  if(i!=LINK->last) {
    multilike(i+1,logtemp,LINK);return;
  }
  LINK->linek=linenum(loci,locik,1);
  LINK->lineq=linenum(loci,lociq,1);
  *logtemp+=2*LINK->hap[LINK->linek-1]*LINK->hap[LINK->lineq-1];return;
}
for(j=1;j<=2;j++) {
  if(j==2) swap__(&locik[i-1],&lociq[i-1]);
  if(i==LINK->last) {
    LINK->linek=linenum(loci, locik, 1);
    LINK->lineq=linenum(loci, lociq, 1);
    *logtemp+=2*LINK->hap[LINK->linek-1]*LINK->hap[LINK->lineq-1];
  } else multilike(i+1, logtemp, LINK);
}
swap__(&locik[i-1], &lociq[i-1]);
}

int probpcell(double *prob, struct LOC_likelihood *LINK)
{
int i;
*prob=0.0;
for(i=0;i<LINK->totalloci;i++) locihn[i]=0;
for(i=0;i<LINK->totalloci;i++) if(locik[i]!=lociq[i]) locihn[i]=1;
LINK->h=0;
for(i=0;i<LINK->totalloci;i++) LINK->h+=locihn[i];
genotype=true;
if(LINK->dma) LINK->linep=linenum(locigeno,locip,2);
else LINK->linep=linenum(locigeno,locip,1);
genotype=false;
if(LINK->h<2) {
  LINK->linek=linenum(loci,locik,1);
  LINK->lineq=linenum(loci,lociq,1);
  LINK->temp=LINK->hap[LINK->linek-1]*LINK->hap[LINK->lineq-1];
  if(LINK->temp>verysmall)
  if(LINK->dma) {
    if(LINK->h==0) *prob+=LINK->temp;else *prob+=2*LINK->temp;
  } else {
    if(LINK->linep!=-1)
    if(LINK->h==0) *prob+=log(LINK->temp)*indivN[LINK->linep-1];
    else *prob+=log(2*LINK->temp)*indivN[LINK->linep-1];
  }
} else {
  getfirst(&LINK->first,&LINK->last);
  LINK->multilog=0.0;
  multilike(LINK->first,&LINK->multilog,LINK);
  if(LINK->multilog>verysmall) {
    if(LINK->dma) *prob+=LINK->multilog;
    else if(LINK->linep!=-1) *prob+=log(LINK->multilog)*indivN[LINK->linep-1];
  }
}
return;
}

double likelihood(double *hap_, int totalloci_, boolean dma_, boolean dg)
{
struct LOC_likelihood V;
double Result;
int i,j,time;
double like, like11, like12, like22, temp1, sa, su;
memcpy(V.hap, hap_, sizeof(linetype));
V.totalloci=totalloci_;
V.dma=dma_;
if(V.dma) {
  iniindex();
  time=1;
  sa=su=0.0;
  do {
    lociq[0]=locik[0]=1;
    probpcell(&like11, &V);
    lociq[0]=2;
    locik[0]=1;
    probpcell(&like12, &V);
    locik[0]=lociq[0]=2;
    probpcell(&like22, &V);
    sa+=like11*ff[0]+like12*ff[1]+like22*ff[2];
    su+=like11*(1-ff[0])+like12*(1-ff[1])+like22*(1-ff[2]);
    locip[V.totalloci-1]++;
    for(i=V.totalloci;i>=2;i--)
      if(locip[i-1]>locigeno[i-1]){
        locip[i-1]=1;
        locip[i-2]++;
      }
    locik[V.totalloci-1]++;
    for(i=V.totalloci;i>=2;i--){
      if(locik[i-1]>lociq[i-1]){
        lociq[i-1]++;
        locik[i-1]=1;
        if(lociq[i-1]>loci[i-1]){
          lociq[i-1]=1;
          locik[i-2]++;
          if(i==2){
            if(locik[i-2]>lociq[i-2]){
              lociq[i-2]++;
              locik[i-2]=1;
            }
          }
        }
      }
    }
    time+=3;
  } while(time<=genonum);
}
iniindex();
time=1;
like=0.0;
do {
  if(V.dma) {
    lociq[0]=locik[0]=1;
    probpcell(&like11, &V);
    lociq[0]=2;
    locik[0]=1;
    probpcell(&like12, &V);
    locik[0]=lociq[0]=2;
    probpcell(&like22, &V);
    if(V.linep!=-1) {
      if(caseN[V.linep-1]>0.01) like+=caseN[V.linep-1]*
        log((like11*ff[0]+like12*ff[1]+like22*ff[2])/sa);
      if(contrlN[V.linep-1]>0.01)like+=contrlN[V.linep-1]*
        log((like11*(1-ff[0])+like12*(1-ff[1])+like22*(1-ff[2]))/su);
    }
    V.temp=ff[0]*like11+ff[1]*like12+ff[2]*like22;
    temp1=(1-ff[0])*like11+(1-ff[1])*like12+(1-ff[2])*like22;
    if(V.linep!=-1) dg11[V.linep-1]=dg12[V.linep-1]=dg22[V.linep-1]=0.0;
    if(V.linep!=-1) {
      if(caseN[V.linep-1]>0.01){
        dg11[V.linep-1]+=like11*caseN[V.linep-1]*ff[0]/V.temp;
        dg12[V.linep-1]+=like12*caseN[V.linep-1]*ff[1]/V.temp;
        dg22[V.linep-1]+=like22*caseN[V.linep-1]*ff[2]/V.temp;
      }
      if(contrlN[V.linep-1]>0.01){
        dg11[V.linep-1]+=like11*contrlN[V.linep-1]*(1-ff[0])/temp1;
        dg12[V.linep-1]+=like12*contrlN[V.linep-1]*(1-ff[1])/temp1;
        dg22[V.linep-1]+=like22*contrlN[V.linep-1]*(1-ff[2])/temp1;
      }
    }
  } else {
    probpcell(&like11, &V);
    like += like11;
  }
  locip[V.totalloci-1]++;
  for(i=V.totalloci;i>=2;i--)
    if(locip[i-1]>locigeno[i-1]){
      locip[i-1]=1;
      locip[i-2]++;
    }
  locik[V.totalloci-1]++;
  for(i=V.totalloci;i>=2;i--){
    if(locik[i-1]>lociq[i-1]){
      lociq[i-1]++;
      locik[i-1]=1;
      if(lociq[i-1]>loci[i-1]){
        lociq[i-1]=1;
        locik[i-2]++;
        if(i==2){
          if(locik[i-2]>lociq[i-2]){
            lociq[i-2]++;
            locik[i-2]=1;
          }
        }
      }
    }
  }
  if(V.dma) time+=3; else time++;
} while(time<=genonum);
Result=like;
if(!dg) return Result;
if(!ehplus)
{
  strcpy(tmp_NAME, "temp.dat");
  if(tmp!=NULL) tmp=freopen(tmp_NAME,"w",tmp);else tmp=fopen(tmp_NAME, "w");
  fprintf(tmp, "2 ");
  for(i=1;i<V.totalloci;i++)fprintf(tmp," %2d",loci[i]);putc('\n', tmp);
  outind(&tmp,locigeno,dg11,2);
  outind(&tmp,locigeno,dg12,2);
  outind(&tmp,locigeno,dg22,2);
  if(tmp!=NULL) fclose(tmp);
}
else
{
  time=1;
  for(i=1;i<totalloci;i++) time*=locigeno[i];
  for(j=0;j<obscom;j++)
  {
    indivN[j]=dg11[j];
    indivN[obscom+j]=dg12[j];
    indivN[obscom*2+j]=dg22[j];
    id[j]=idsave[j];
    id[obscom+j]=time+idsave[j];
    id[obscom*2+j]=2*time+idsave[j];
  }
}
return Result;
}

struct LOC_gethapN {
  int first, last, linek, lineq, linep;
  linetype vhaplop;
} ;

int multiadd(int i, int *line, struct LOC_gethapN *LINK)
{
int j;
if(locihn[i-1]!=1){
  multiadd(i+1,line,LINK);return;
}
if(i==LINK->first) {
  if(i!=LINK->last) {
    multiadd(i+1,line,LINK);return;
  }
  LINK->linek=linenum(loci,locik,1);
  LINK->lineq=linenum(loci,lociq,1);
  genotype=true;
  LINK->linep=linenum(locigeno,locip,1);
  genotype=false;
  if(LINK->linep==-1) return;
  newhaplo[LINK->linek-1]+=LINK->vhaplop[*line]*indivN[LINK->linep-1];
  newhaplo[LINK->lineq-1]+=LINK->vhaplop[*line]*indivN[LINK->linep-1];
  (*line)++;
  return;
}
for(j=1;j<=2;j++){
  if(j==2) swap__(&locik[i-1],&lociq[i-1]);
  if(i==LINK->last){
    LINK->linek=linenum(loci,locik,1);
    LINK->lineq=linenum(loci,lociq,1);
    genotype=true;
    LINK->linep=linenum(locigeno,locip,1);
    genotype=false;
    if(LINK->linep==-1) continue;
    newhaplo[LINK->linek-1]+=LINK->vhaplop[*line]*indivN[LINK->linep-1];
    newhaplo[LINK->lineq-1]+=LINK->vhaplop[*line]*indivN[LINK->linep-1];
    (*line)++;
  } else multiadd(i+1,line,LINK);
}
swap__(&locik[i-1],&lociq[i-1]);
}

int multigen(int i, int *line, double *totalg, struct LOC_gethapN *LINK)
{
int j;
if(locihn[i-1]!=1) {
  multigen(i+1,line,totalg,LINK);return;
}
if(i==LINK->first) {
  if(i!=LINK->last) {
    multigen(i+1,line,totalg,LINK);return;
  }
  LINK->linek=linenum(loci,locik,1);
  LINK->lineq=linenum(loci,lociq,1);
  LINK->vhaplop[*line]=2*oldhaplo[LINK->linek-1]*oldhaplo[LINK->lineq-1];
  *totalg += LINK->vhaplop[*line];
  (*line)++;
  return;
}
for(j=1;j<=2;j++){
  if(j==2) swap__(&locik[i-1],&lociq[i-1]);
  if(i==LINK->last){
    LINK->linek=linenum(loci,locik,1);
    LINK->lineq=linenum(loci,lociq,1);
    LINK->vhaplop[*line]=2*oldhaplo[LINK->linek-1]*oldhaplo[LINK->lineq-1];
    *totalg+=LINK->vhaplop[*line];
    (*line)++;
  } else multigen(i+1,line,totalg,LINK);
}
swap__(&locik[i-1], &lociq[i-1]);
}

void gethapN()
{
struct LOC_gethapN V;
int i,j,time,h,subn,lines;
iniindex();
newtoold();
time=1;
do {
  for(i=0;i<totalloci;i++) locihn[i]=0;
  lines=linenum(locigeno, locip, 1);
  for(i=0;i<totalloci;i++) if(locik[i]!=lociq[i]) locihn[i]=1;
  h=0;
  for(i=0;i<totalloci;i++) h+=locihn[i];
  if(h<2){
    V.linek=linenum(loci,locik,1);
    V.lineq=linenum(loci,lociq,1);
    genotype=true;
    V.linep=linenum(locigeno,locip,1);
    genotype=false;
    if(V.linep!=-1)
    {
      newhaplo[V.linek-1]+=indivN[V.linep-1];
      newhaplo[V.lineq-1]+=indivN[V.linep-1];
    }
  } else if(h>1) {
    getfirst(&V.first,&V.last);
    lines=0;
    totalg=0.0;
    multigen(V.first,&lines,&totalg,&V);
    subn=exp2(h-1);
    if(totalg<verysmall) for(i=0;i<subn;i++) V.vhaplop[i]=0.0;
    else for(i=0;i<subn;i++) V.vhaplop[i]/=totalg;
    lines=0;
    multiadd(V.first, &lines, &V);
  }
  locip[totalloci-1]++;
  for(i=totalloci;i>=2;i--)
    if(locip[i-1]>locigeno[i-1]){
    locip[i-1]=1;
    locip[i-2]++;
    }
  locik[totalloci-1]++;
  if(totalloci>1) {
    for(i=totalloci;i>=2;i--){
      if(locik[i-1]>lociq[i-1]){
        lociq[i-1]++;
        locik[i-1]=1;
        if(lociq[i-1]>loci[i-1]){
          lociq[i-1]=1;
          locik[i-2]++;
          if(i==2){
            if(locik[i-2]>lociq[i-2]){
              lociq[i-2]++;
              locik[i-2]=1;
            }
          }
        }
      }
    }
  } else {
    if(locik[i-1]>lociq[i-1]){
      lociq[i-1]++;
      locik[i-1]=1;
    }
  }
  time++;
} while(time<=genonum);
}

double getdiff()
{
  int i;
  double temp=0.0;
  for(i=0;i<haplnum;i++) temp+=fabs(oldhaplo[i]-newhaplo[i]);
  return temp;
}

double getp(int *l,int k,int totalloci,double (*p)[maxalle],int *line,int *ll)
{
int i,j;
double n,tempn,temp;
tempn=0;
for(i=0;i<l[k-1];i++) for(j=0;j<i+1;j++)
  if(k==totalloci)
  {
    (*line)++;
    if(*line==id[*ll]) {n=indivN[*ll];(*ll)++;} else n=0;
    p[k-1][i]+=n;
    p[k-1][j]+=n;
    tempn+=n;
  }
  else
  {
    temp=getp(l,k+1,totalloci,p,line,ll);
    p[k-1][i]+=temp;
    p[k-1][j]+=temp;
    tempn+=temp;
  }
return tempn;
}

int checklimit()
/*Check default array upper bounds*/
{
int i,lh,lg,obsbin;
if(totalloci>maxloci) fprintf(stderr,"set maxloci=%d\n",totalloci);
lh=lg=1;
for(i=0;i<totalloci;i++) {
  lh*=loci[i];lg*=loci[i]*(loci[i]+1)/2;
}
if(casestudy) lh*=2;
if(maxall>maxalle) fprintf(stderr,"set maxalle=%d\n",maxall);
if(lh>=maxhap) fprintf(stderr,"set maxhap=%d\n",lh);
if(ehplus) {
  if(obscom>maxobscom) fprintf(stderr,"set maxobscom=%d\n",obscom);
  obsbin=obscom/PRIME;
  if(casestudy) obsbin*=3;
  if(obsbin>maxbinsize) fprintf(stderr,"set maxbinsize=%d\n",obsbin);
} else {
  if(lg>maxposcom) fprintf(stderr,"set maxposcom=%d\n",lg);;
}
return 0;
}

int main(int argc, char **argv)
{
FILE *fi,*fo;
char s[FNL],line[100],rest[100];
int j,j1,j2,n,nloci;
tmp=wf=casef=controlf=NULL;
logo();
ch='n';
fi=stdin;
fo=stdout;
printf("Do you wish to use the case-control sampling option?  [N]\n");
if(seekeoln(&fi)) {
  fscanf(fi,"%*[^\n]");fgetc(fi);
} else {
  fscanf(fi,"%c%*[^\n]", &ch);fgetc(fi);
}
casestudy=(ch=='Y'||ch=='y');
if(ehplus) goto plus;
if(casestudy) {
  do {
    err=false;strcpy(fname, "control.dat");
    printf("Enter name of control data file  [CONTROL.DAT]\n");
    if(seekeoln(&fi)) fgetc(fi);
    else {fgets(s,FNL,fi);sscanf(s,"%s",fname);}
    strcpy(controlf_NAME, fname);
    if(controlf!=NULL) controlf=freopen(controlf_NAME,"r",controlf);
    else controlf=fopen(controlf_NAME, "r");
    if(controlf==NULL) err=true;
  } while(err);
  do {
    err=false;strcpy(fname, "case.dat");
    printf("Enter name of case data file   [CASE.DAT]\n");
    if(seekeoln(&fi))fgetc(fi);
    else {fgets(s,FNL,fi);sscanf(s,"%s",fname);}
    strcpy(casef_NAME, fname);
    if(casef!=NULL) casef=freopen(casef_NAME,"r",casef);
    else casef=fopen(casef_NAME,"r");
    if(casef==NULL) err=true;
  } while(err);
} else
  do {
    err=false;strcpy(fname, "eh.dat");
    printf("Enter name of data file  [EH.DAT]\n");
    if(seekeoln(&fi)) fgetc(fi);
    else {fgets(s,FNL,fi);sscanf(s,"%s",fname);}
    printf("you entered: %s\n", fname);
    strcpy(controlf_NAME, fname);
    if(controlf!=NULL) controlf=freopen(controlf_NAME,"r",controlf);
    else controlf=fopen(controlf_NAME, "r");
    if(controlf==NULL) err=true;
  } while(err);
printf("Enter name of output file.  [EH.OUT]\n");
strcpy(fname, "eh.out");
if(seekeoln(&fi)) fgetc(fi);else {fgets(s,FNL,fi);sscanf(s,"%s",fname);}
strcpy(wf_NAME, fname);
wf=fopen(wf_NAME,"r");fclose(wf);
if(wf) {
  printf("File %s exist. Over write it? [y/n]\007\n", fname);
  if(seekeoln(&fi)) {
    fscanf(fi,"%*[^\n]");fgetc(fi);
    strcpy(wf_NAME, fname);
    if(wf!=NULL) wf=freopen(wf_NAME,"w",wf);else wf=fopen(wf_NAME,"w");
  } else {
    fscanf(fi,"%c%*[^\n]", &ch);fgetc(fi);
    if(ch!='Y'&&ch!='y') {
      programaborted();
      goto _L999;
    }
    strcpy(wf_NAME, fname);
    if(wf!=NULL) wf=freopen(wf_NAME,"w",wf);else wf=fopen(wf_NAME,"w");
  }
} else {
  strcpy(wf_NAME, fname);
  if(wf!=NULL) wf=freopen(wf_NAME,"w",wf);else wf=fopen(wf_NAME,"w");
}
err=false;
logo();
printf("        Program is running\n\n");
indivN=(double *)calloc(maxposcom,sizeof(double));
if(casestudy) {
  caseN=(double*)calloc(maxposcom,sizeof(double));
  contrlN=(double*)calloc(maxposcom,sizeof(double));
  dg11=(double*)calloc(maxposcom,sizeof(double));
  dg12=(double*)calloc(maxposcom,sizeof(double));
  dg22=(double*)calloc(maxposcom,sizeof(double));
  dg33=(double*)calloc(maxposcom,sizeof(double));
  if(!inputdataok(&controlf,&fo,&totalloci,loci)) {
    programaborted();
    printf("        Reason: Input individual genotype should be %d\n",genonum);
    goto _L999;
  }
  if(controlf!=NULL) controlf=freopen(controlf_NAME,"r",controlf);
  else controlf=fopen(controlf_NAME, "r");
  fscanf(controlf, "%*[^\n]");getc(controlf);
  if(!inputdataok(&casef,&fo,&cnloci,caseloci)) {
    programaborted();
    printf("        Reason: Input individual genotype should be %d\n",genonum);
    goto _L999;
  }
  if(casef!=NULL) casef=freopen(casef_NAME,"r",casef);
  else casef=fopen(casef_NAME, "r");
  fscanf(casef, "%*[^\n]");getc(casef);
  if(!dataconsist(loci,caseloci,totalloci,cnloci)) {
    programaborted();
    printf("        Reason: Input data files are not consistant\n");
    goto _L999;
  }
  getlocih(totalloci);
  getfreq(fi);
  i=0;
  casecontrol(&controlf,&casef,loci,1,totalloci,&i);
  strcpy(fname,"tempp.dat");
  strcpy(tmp_NAME,fname);
  if(tmp!=NULL) tmp=freopen(tmp_NAME,"w",tmp);else tmp=fopen(tmp_NAME,"w");
  for(i=0;i<totalloci;i++) fprintf(tmp," %2d",loci[i]);putc('\n',tmp);
  outind(&tmp, locigeno, dg33, 1);
  if(tmp!=NULL) fclose(tmp);
  strcpy(controlf_NAME, "tempp.dat");
  if(controlf!=NULL) controlf=freopen(controlf_NAME, "r", controlf);
  else controlf = fopen(controlf_NAME,"r");
  if(!inputdataok(&controlf,&fo,&totalloci,loci)) {
    programaborted();
    printf("        Reason: Input individual genotype should be %d\n",genonum);
    goto _L999;
  }
  strcpy(controlf_NAME,fname);
  if(controlf!=NULL) controlf=freopen(controlf_NAME,"r",controlf);
  else controlf=fopen(controlf_NAME,"r");
  fscanf(controlf,"%*[^\n]");getc(controlf);
  initdg();
  getlocih(totalloci);
  initial(p);
  i=j=0;
  totalind=getN(&controlf,loci,1,totalloci,p,&i);
  if(controlf!=NULL) fclose(controlf);
  totalall=s2=2*totalind;
  calcula(totalloci);
  gethapP(totalloci);
  onemark=(totalloci==1);
  if(!onemark) {
    loop=0;
    do {
      loop++;
      gethapN();
      for(i=0;i<haplnum;i++) newhaplo[i]/=totalall;
      loglike=likelihood(newhaplo,totalloci,false,false);
      diff=getdiff();
      printf("Iteration=%d  log likelihood=%8.5f   diffence=%8.5f",loop, loglike, diff);
      for(i=1;i<=60;i++) putchar(bs);
      err=(loop>100 ||diff<0.0001);
    } while(!err);
    if(loop>100) printf("Warning: Exceeds iteration limit\007\n");
  }
  for(i=0;i<haplnum;i++) {
    indephaplo[i]=newhaplo[i]*(1-pp);
    indephaplo[i+haplnum]=newhaplo[i]*pp;
  }
  loop=0;
  do {
    strcpy(controlf_NAME,"temp.dat");
    if(controlf!=NULL) controlf=freopen(controlf_NAME,"r",controlf);
    else controlf=fopen(controlf_NAME,"r");
    printf("     Case-control data\n");
    if(!inputdataok(&controlf,&fo,&totalloci,loci)) {
      programaborted();
      printf("        Reason: Input individual genotype should be %d\n",genonum);
      goto _L999;
    }
    if(controlf!=NULL) controlf=freopen(controlf_NAME,"r",controlf);
    else controlf=fopen(controlf_NAME,"r");
    fscanf(controlf, "%*[^\n]");getc(controlf);
    i=0;
    getlocih(totalloci);
    initial(p);
    totalind=getN(&controlf,loci,1,totalloci,p,&i);
    if(controlf!=NULL) fclose(controlf);
    totalall=s2=2*totalind;
    calcula(totalloci);
    if(loop==0) {
      p[0][0]=1-pp;
      p[0][1]=pp;
      outindP(totalloci);
      gethapP(totalloci);
      iniloglike=likelihood(newhaplo,totalloci,true,false);
      indloglike=likelihood(indephaplo,totalloci,true,false);
      loglike=iniloglike;
    }
    oldlog=loglike;
    loop1=0;
    do {
      gethapN();
      for(i=0;i<haplnum;i++) newhaplo[i]/=totalall;
      normalizeH(loci, newhaplo);
      diff=getdiff();
      printf("Iteration=%d  log likelihood=%8.5f   diffence=%8.5f",loop1, loglike, diff);
      for(i=0;i<60;i++) putchar(bs);
      loop1++;
    } while(diff>=0.0001&&loop1<=100);
    loop++;
    loglike=likelihood(newhaplo,totalloci,true,true);
    diff=loglike-oldlog;
  } while(diff>=0.0001&&loop<=15);
} else {
  printf("     Regular data\n");
  if(!inputdataok(&controlf,&fo,&totalloci,loci)) {
    programaborted();
    printf("        Reason: Input individual genotype should be %d\n",genonum);
    goto _L999;
  }
  if(controlf!=NULL) controlf=freopen(controlf_NAME,"r",controlf);
  else controlf=fopen(controlf_NAME,"r");
  fscanf(controlf, "%*[^\n]");getc(controlf);
  getlocih(totalloci);
  initial(p);
  i=0;
  totalind=getN(&controlf,loci,1,totalloci,p,&i);
  if(controlf!=NULL) fclose(controlf);
  totalall=s2=2*totalind;
  calcula(totalloci);
  loop=0;
  outindP(totalloci);
  gethapP(totalloci);
  iniloglike=likelihood(newhaplo,totalloci,false,false);
  do {
    loop++;
    gethapN();
    for(i=0;i<haplnum;i++) newhaplo[i]/=totalall;
    loglike=likelihood(newhaplo,totalloci,false,false);
    diff=getdiff();
    printf("Iteration=%d  log likelihood=%8.5f   diffence=%8.5f",loop,loglike,diff);
    for(i=0;i<60;i++) putchar(bs);
    err=(loop>100||diff<0.0001);
  } while(!err);
}
goto end;

plus:
do {
  err=false;strcpy(fname, "ehplus.dat");
  printf("Enter name of data file  [EHPLUS.DAT]\n");
  if(seekeoln(&fi)) fgetc(fi);
  else {fgets(s,FNL,fi);sscanf(s,"%s",fname);}
  printf("you entered: %s\n", fname);
  strcpy(controlf_NAME, fname);
  if(controlf!=NULL) controlf=freopen(controlf_NAME,"r",controlf);
  else controlf=fopen(controlf_NAME, "r");
  if(controlf==NULL) err=true;
} while(err);
printf("Enter name of output file.  [EHPLUS.OUT]\n");
strcpy(fname, "ehplus.out");
if(seekeoln(&fi)) fgetc(fi);else {fgets(s,FNL,fi);sscanf(s,"%s",fname);}
strcpy(wf_NAME, fname);
wf=fopen(wf_NAME,"r");fclose(wf);
if(wf) {
  printf("File %s exist. Over write it? [y/n]\007\n", fname);
  if(seekeoln(&fi)) {
    fscanf(fi,"%*[^\n]");fgetc(fi);
    strcpy(wf_NAME, fname);
    if(wf!=NULL) wf=freopen(wf_NAME,"w",wf);else wf=fopen(wf_NAME,"w");
  } else {
    fscanf(fi,"%c%*[^\n]", &ch);fgetc(fi);
    if(ch!='Y'&&ch!='y') {
      programaborted();
      goto _L999;
    }
    strcpy(wf_NAME, fname);
    if(wf!=NULL) wf=freopen(wf_NAME,"w",wf);else wf=fopen(wf_NAME,"w");
  }
} else {
  strcpy(wf_NAME, fname);
  if(wf!=NULL) wf=freopen(wf_NAME,"w",wf);else wf=fopen(wf_NAME,"w");
}
logo();
printf("        Program is running\n\n");
indivN=(double *)malloc(3*maxobscom*sizeof(double));
if(casestudy)
{
caseN=(double*)malloc(maxobscom*sizeof(double));
contrlN=(double*)malloc(maxobscom*sizeof(double));
}
controlf=fopen(controlf_NAME,"r");
fgets(line,100,controlf);
maxall=nloci=0;
while(sscanf(line,"%d %[^\n]",&loci[nloci++],rest)>1) strcpy(line,rest);
for(i=0;i<nloci;i++) if(loci[i]>maxall) maxall=loci[i];
totalloci=nloci;
inivar(nloci,loci);
for(i=0;i<PRIME;i++) {
  bin[i].n=0;
  for(j=0;j<maxbinsize;j++) bin[i].id[j]=bin[i].obs[j]=0;
}
totalind=i=0;
if(casestudy) while(fgets(line,100,controlf)&&
sscanf(line,"%ld %lf %lf %lf",&idsave[i],&indivN[i],&caseN[i],&contrlN[i])==4)
{
j=hash(idsave[i]);
j1=bin[j].n;
bin[j].obs[j1]=i;
bin[j].id[j1]=idsave[i];
bin[j].n++;
totalind+=indivN[i];id[i]=idsave[i];
++i;
}
else while(fgets(line,100,controlf)&&
sscanf(line,"%ld %lf ",&idsave[i],&indivN[i])==2)
{
j=hash(idsave[i]);
j1=bin[j].n;
bin[j].obs[j1]=i;
bin[j].id[j1]=idsave[i];
bin[j].n++;
totalind+=indivN[i];id[i]=idsave[i];
++i;
}
fclose(controlf);
sample_size=obscom=i;
for(j=0;j<PRIME;j++) {
/*
  printf("%5d %3d ",j,bin[j].n);
  for(j1=0;j1<bin[j].n;j1++) printf("%ld ",bin[j].id[j1]);printf("\n");
*/
}
checklimit();
if(casestudy)
{
  getfreq(fi);
  dg11=(double *)malloc(obscom*sizeof(double));
  dg12=(double *)malloc(obscom*sizeof(double));
  dg22=(double *)malloc(obscom*sizeof(double));
  for(i=0;i<obscom;i++)
  {
    dg11[i]=casep[0]*caseN[i]+controlp[0]*contrlN[i];
    dg12[i]=casep[1]*caseN[i]+controlp[1]*contrlN[i];
    dg22[i]=casep[2]*caseN[i]+controlp[2]*contrlN[i];
  }
  getlocih(totalloci);
  initial(p);
  i=j=0;
  totalind=getp(loci,1,totalloci,p,&i,&j);
  totalall=s2=2*totalind;
  calcula(totalloci);
  gethapP(totalloci);
  onemark=(boolean)(totalloci==1);
  if(!onemark) {
    loop=0;
    do {
      loop++;
      gethapN();
      for(i=0;i<haplnum;i++) newhaplo[i]/=totalall;
      loglike=likelihood(newhaplo,totalloci,false,false);
      diff=getdiff();
      printf("Iteration=%d  log likelihood=%8.5f   diffence=%8.5f",loop, loglike, diff);
      for(i=0;i<60;i++) putchar(bs);
      err=(loop>100||diff<0.0001);
    } while(!err);
    if(loop>100) printf("Warning: Exceeds iteration limit\007\n");
  }
  for(i=0;i<haplnum;i++) {
    indephaplo[i]=newhaplo[i]*(1-pp);
    indephaplo[i+haplnum]=newhaplo[i]*pp;
  }
  n=1;
  for(i=0;i<totalloci;i++) n*=locigeno[i];
  for(i=totalloci+1;i>0;i--) loci[i]=loci[i-1];loci[0]=2;
  for(j=0;j<obscom;j++)
  {
      indivN[j]=dg11[j];
      indivN[obscom+j]=dg12[j];
      indivN[obscom*2+j]=dg22[j];
      id[j]=idsave[j];
      id[obscom+j]=n+idsave[j];
      id[obscom*2+j]=2*n+idsave[j];
      j1=hash(id[j]);
      j2=bin[j1].n;
      bin[j1].obs[j2]=j;
      bin[j1].id[j2]=id[j];
      bin[j1].n++;
      j1=hash(id[obscom+j]);
      j2=bin[j1].n;
      bin[j1].obs[j2]=obscom+j;
      bin[j1].id[j2]=id[obscom+j];
      bin[j1].n++;
      j1=hash(id[obscom*2+j]);
      j2=bin[j1].n;
      bin[j1].obs[j2]=obscom*2+j;
      bin[j1].id[j2]=id[obscom*2+j];
      bin[j1].n++;
  }
  loop=0;
  totalloci++;nloci++;haplnum*=2;
  inivar(nloci,loci);
  printf("\n     Case-control data\n");
  do {
    getlocih(totalloci);
    initial(p);
    i=j=0;
    totalind=getp(loci,1,totalloci,p,&i,&j);
    totalall=s2=2*totalind;
    calcula(totalloci);
    if(loop==0) {
      p[0][0]=1-pp;
      p[0][1]=pp;
      outindP(totalloci);
      gethapP(totalloci);
      iniloglike=likelihood(newhaplo,totalloci,true,false);
      indloglike=likelihood(indephaplo,totalloci,true,false);
      loglike=iniloglike;
    }
    oldlog=loglike;
    loop1=0;
    do {
      gethapN();
      for(i=0;i<haplnum;i++) newhaplo[i]/=totalall;
      normalizeH(loci,newhaplo);
      diff=getdiff();
      printf("Iteration=%d  log likelihood=%8.5f   diffence=%8.5f",loop1, loglike, diff);
      for(i=0;i<60;i++) putchar(bs);
      loop1++;
    } while(diff>=0.0001&&loop1<=100);
    loop++;
    loglike=likelihood(newhaplo,totalloci,true,true);
    diff=loglike-oldlog;
  } while(diff>=0.0001&&loop<=15);
} else {
  printf("     Regular data\n");
  getlocih(totalloci);
  initial(p);
  i=j=0;
  totalind=getp(loci,1,totalloci,p,&i,&j);
  totalall=s2=2*totalind;
  calcula(totalloci);
  outindP(totalloci);
  gethapP(totalloci);
  loop=0;
  iniloglike=likelihood(newhaplo,totalloci,false,false);
  do {
    loop++;
    gethapN();
    for(i=0;i<haplnum;i++) newhaplo[i]/=totalall;
    loglike=likelihood(newhaplo,totalloci,false,false);
    diff=getdiff();
    printf("Iteration=%d  log likelihood=%8.5f   diffence=%8.5f",loop,loglike,diff);
    for(i=0;i<60;i++) putchar(bs);
    err=(loop>100||diff<0.0001);
  } while(!err);
}
end:
free(indivN);
if(casestudy){
free(caseN);free(contrlN);free(dg11);free(dg12);free(dg22);
if(!ehplus) free(dg33);
}
outvec(loci,newhaplo);
if(wf!=NULL) fclose(wf);
programsuccess();
if(loop>100) printf("Warning: bad data set\007\n");
_L999:
if(controlf!=NULL) fclose(controlf);
if(casef!=NULL) fclose(casef);
if(tmp!=NULL) fclose(tmp);
return 0;
}
