/*                                                                   */
/* GENECOUNTING/PERMUTE - haplotype analysis with missing genotypes  */
/*                        and (global/local) permutation tests       */
/*                                                                   */
/* (C) Copyright JH Zhao, University College London, 2003            */
/*                                                                   */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <time.h>
#include "pgc.h"

#define version 1.2
#ifndef HAPWISE
int     hapwise=0;
#else
int     hapwise=1;
#endif

#ifdef  LINUX
char *gcascmd="gc.linux temp.gc temp1.gco>gc.tmp";
char *gconcmd="gc.linux temp.gc temp2.gco>gc.tmp";
char *gcomcmd="gc.linux temp.gc temp.gco>gc.tmp";
#elif   SUNOS
char *gcascmd="gc.sun temp.gc temp1.gco>gc.tmp";
char *gconcmd="gc.sun temp.gc temp2.gco>gc.tmp";
char *gcomcmd="gc.sun temp.gc temp.gco>gc.tmp";
#else
char *gcascmd="gc temp.gc temp1.gco>gc.tmp";
char *gconcmd="gc temp.gc temp2.gco>gc.tmp";
char *gcomcmd="gc temp.gc temp.gco>gc.tmp";
#endif

int shuffle(int *,int *,int,int);
float getx2ln(float*),getln(char*),getx2het(),getx2block(int,float*);

typedef struct
{
        int id,cc,locus[maxloci][2];
} ind;
ind ind_t,*individuals=0;
enum {ALL,BLOCK1,BLOCK2};
int seed=3000,ip;
double pl,epsh=0.000001;

void *xmalloc(int);

void getnhaps(FILE*,int*,int*);
int geth(FILE*,int,int);
int gethappy(int,int);
int revhapid(double,int*);
double probnorm(double);
double probchi(double, int);

enum {MARKERMAR,MARKERDIS};
enum {CONTROLS,CASES,COMBINED};

int fobs,pobs,fobs0,pobs0;
int fobsca,fobsco,pobsca,pobsco,fobsca0,fobsco0,pobsca0,pobsco0;
int nhaps,nhaps0,nhapcas,nhapcon,nhapcas0,nhapcon0,*hap;
double *h0,*h1,*h1cas,*h1con,*h1caobs,*h1coobs;
double *hid,*hidobs,*hidca,*hidco,*hidcaobs,*hidcoobs;
double *ft,*ftobs,*ftp,*ftcc,*ftccobs,*ftpcc,*fticc,*fticcobs,*ftss,*ftssobs;
long int *nalleles;

int main(int argc, char *argv[])
{
int i,j,k,l,m,ll,irun;
int *a,*b,*a1,*a2;
double p1,p2;
float x2obs,x2,x2obs1,x21,x2obs2,x22,x2obs12,x212,lla,ll1,ll2,p,se;
int dfa,df1,df2,df12,j1,j2;
FILE *gout;
time_t tod;
list t;

printf("\nGENECOUNTING/PERMUTE, %.2f JH Zhao 05/03\n",version);
time(&tod);
printf("%s\n",ctime(&tod));
printf("Maximum number of loci = %d\n",maxloci);
printf("Maximum number of alleles = %d\n\n",maxalleles);

if(argc<3)
{
  fprintf(stderr,"Usage: %s parfile datfile [<outfile> [seed]]\n\n",argv[0]);
  fprintf(stderr,"where:\n\n");
  fprintf(stderr,"parfile and datfile are EHPLUS parameter and data files\n");
  fprintf(stderr,"outfile is output file, omitted when using computer screen\n");
  fprintf(stderr,"seed is an integer for pseudorandom number generator\n");
  exit(1);
}
else
if (argc>4) seed=atoi(argv[4]);
srand(seed);

j=j1=j2=0;
if (!getloci(argv[1]))
{
   printf("\n");
   noid_new(argv[2]);
   individuals=(ind*)malloc(sample_size*sizeof(ind));
   if (!individuals)
   {
      perror("error allocating memory for individuals");
      exit(1);
   }
   t=r;
   i=0;
   do
   {
     individuals[i].id=t->id;
     individuals[i].cc=t->cc;
     for(k=0;k<nloci;k++)
     {
        individuals[i].locus[k][0]=t->locus[k][0];
        individuals[i].locus[k][1]=t->locus[k][1];
     }
     i++;
     t=t->next;
   } while (t);
   a=(int*)malloc(sample_size*sizeof(int));
   b=(int*)malloc(sample_size*sizeof(int));
   if(!a||!b)
   {
      perror("error allocating memory for shuffle arrays");
      exit(1);
   }
   a1=(int*)malloc(sample_size*sizeof(int));
   a2=(int*)malloc(sample_size*sizeof(int));
   nalleles=(long int*)malloc((selected+1)*sizeof(long int));
   hap=(int*)malloc((selected+1)*sizeof(int));
   if (!a1||!a2||!nalleles||!hap)
   {
      perror("error allocating memory for alleles");
      exit(1);
   }
   dfa=df1=df2=df12=1;
   m=0;
   for(k=0;k<nloci;k++) if(sel[k])
   {
      dfa*=alleles[k];
      m++;
   }
   m=selected;
   nalleles[selected]=1;
   for(k=nloci-1;k>0;k--) if(sel[k])
   {
      m--;
      nalleles[m]=nalleles[m+1]*alleles[k];
   }
   dfa--;
   df12=dfa;
   for(k=0;k<nloci;k++) if(sel[k]) dfa-=(alleles[k]-1);
   ip=0;
   if (cc)
   {
      j=0;
      do
      {
        x2=getx2het();
        k=sample_size;
        while (--k>0)
        {
            i=rand()%k;
            ind_t=individuals[i];
            individuals[i]=individuals[k];
            individuals[k]=ind_t;
        }
        for(l=0;l<sample_size;l++) individuals[l].cc=(l<cases)?1:0;
        if(ip==0) x2obs=x2;
        else
        {
          if(x2>=x2obs) j++;
          printf("Permutation: %d/%d, Chi-square=%f\r",ip,npermute,x2);
        }
      } while (ip++<npermute);
   }
   else
   {
      x2=x21=x22=x212=0;
      x2obs=x2obs1=x2obs2=x2obs12=0;
      for (k=0;k<nloci;k++) if(sel[k])
      {
          if(selp[k]) df1*=alleles[k];
          else df2*=alleles[k];
      }
      df1--;
      df2--;
      df12-=df1+df2;
      for (k=0;k<nloci;k++) if(sel[k])
      {
          if(selp[k]) df1-=(alleles[k]-1);
          else df2-=(alleles[k]-1);
      }
      do
      {
        x2=getx2block(ALL,&lla);
        if(selectp>=0) x21=getx2block(BLOCK1,&ll1);
        if (ip==0&&selectn>=0) x22=getx2block(BLOCK2,&ll2);
        if (selectp>=0&&selectn>=0) x212=2*(lla-ll1-ll2);
        if (ip==0)
        {
           x2obs=x2;
           if (selectp>=0) x2obs1=x21;
           if (selectn>=0) x2obs2=x22;
           if (selectp>=0&&selectn>=0) x2obs12=x212;
        }
        else
        {
          if (x2>=x2obs) j++;
          if ((selectp>=0)&&(x21>=x2obs1)) j1++;
          if ((selectp>=0&&selectn>=0)&&(x212>=x2obs12)) j2++;
          for (k=0;k<nloci;k++)
          {
              if(sel[k]) shuffle(a,b,sample_size,k);
              else continue;
              for (l=0;l<sample_size;l++)
              {
                  a1[l]=individuals[l].locus[k][0];
                  a2[l]=individuals[l].locus[k][1];
              }
              for (l=0;l<sample_size;l++)
              {
                  individuals[l].locus[k][0]=a1[a[l]];
                  individuals[l].locus[k][1]=a2[a[l]];
              }
          }
        }
        printf("permutation: %d/%d, Chi-square=%f\r",ip,npermute,x2);
      } while(ip++<npermute);
   }
   if (npermute>0)
   {
      free(a);
      free(b);
      free(a1);
      free(a2);
   }
   free(individuals);
}
remove("gc.tmp");
remove("temp.gc");
remove("temp.gco");
if (cc)
{
   remove("temp1.gco");
   remove("temp2.gco");
}
if(argc<4) gout=stdout;
else
{
  gout=fopen(argv[3],"w");
  if(!gout)
  {
    fprintf(stderr,"\nI can't open file %s for output",argv[3]);
    fprintf(stderr,"\nSo I will write to the screen\n");
    gout=stdout;
  }
}
fprintf(gout,"\nGENECOUNTING/PERMUTE, %.2f JH Zhao 03/03\n",version);
fprintf(gout,"%s\n",ctime(&tod));
fprintf(gout,"Parameter file=%s, data file=%s\n\n\n",argv[1],argv[2]);
fprintf(gout,"*** Global test of association ***\n\n");
fprintf(gout,"Observed Chi-square=%.2f, df=%d\n\n",x2obs,dfa);
if (npermute>1)
{
   p=(float)j/npermute;
   se=sqrt(p*(1-p)/npermute);
   fprintf(gout,"p value=%f, standard error = %f\n\n",p,se);
   fprintf(gout,"Random number seed   = %d\n",seed);
   fprintf(gout,"Number of replicates = %d\n",npermute);
}
if (!cc&&selectp>=0)
{
   fprintf(gout,"Observed Chi-square for blocks: \n");
   if (selectp>=0) fprintf(gout,"block   1=%.2f, df=%d\n",x2obs1,df1);
   if (selectn>=0) fprintf(gout,"block   2=%.2f, df=%d\n",x2obs2,df2);
   if (selectp>=0&&selectn>=0)
      fprintf(gout,"block 1/2=%.2f, df=%d\n\n",x2obs12,df12);
   if (npermute>1)
   {
      if (selectp>=0)
      {
         fprintf(gout,"p value and standard error: \n");
         p=(float)j1/npermute;
         se=sqrt(p*(1-p)/npermute);
         fprintf(gout,"Block   1=%f, se=%f\n",p,se);
         fprintf(gout,"Block 1 contains markers that are permuted\n\n");
      }
      if (selectn>=0)
      {
         p=(float)j2/npermute;
         se=sqrt(p*(1-p)/npermute);
         fprintf(gout,"Block 1/2=%f, se=%f\n",p,se);
         fprintf(gout,"Block 2 contains markers that are not permuted\n");
         fprintf(gout,"Block 1/2 examines association between Blocks 1 & 2\n");
      }
   }
}
/*obtain empirical p value*/
if (hapwise&&npermute>1)
{
   fprintf(gout,"\n\n*** Significance of individual haplotypes ***\n\n");
   switch (cc)
   {
   case 0:
     fprintf(gout,"%d,%d individuals with Full/Partial genotype\n",fobs,pobs);
     fprintf(gout,"%d nonzero haplotypes\n\n",nhaps0);
     fprintf(gout,"      h0       h1     FT       p   Emp. p     haplotypes ID\n\n");
     for (j=0;j<nhaps0;j++)
     {
         ftp[j]/=npermute;
         fprintf(gout,"%f %f %6.2f",h0[j],h1[j],ftobs[j]);
         fprintf(gout,"%8.5f %8.5f ",probnorm(fabs(ftobs[j])),ftp[j]);
         revhapid(hidobs[j],hap);
         for(m=0;m<=selected;m++) fprintf(gout,"%3d",hap[m]);
         fprintf(gout," %.0f\n",hidobs[j]);
     }
     fprintf(gout,"\nh0 = haplotype frequency assuming linkage equilibrium\n");
     fprintf(gout,"h1 = haplotype frequency assuming linkage disequilibrium\n");
     fprintf(gout,"ID = haplotype identifier\n");
     break;
   case 1:
     fprintf(gout,"%d,%d cases with Full/Partial genotype\n",fobsca0,pobsca0);
     fprintf(gout,"%d nonzero haplotypes\n\n",nhapcas0);
     fprintf(gout,"%d,%d controls with Full/Partial genotype\n",fobsco0,pobsco0);
     fprintf(gout,"%d nonzero haplotypes\n\n",nhapcon0);
     fprintf(gout,"  Case h Contr. h       z        p   Emp. p    haplotypes  ID\n\n");
     i=j=ll=0;
     for (k=0;k<nhaps0;k++)
     {
         if (i<nhapcas0||j<nhapcon0)
         if (hidcaobs[i]<hidcoobs[j])
         {
            l=hidcaobs[i];
            p1=h1caobs[i++];
            p2=0;
         }
         else if (hidcaobs[i]>hidcoobs[j])
         {
            l=hidcoobs[j];
            p1=0;
            p2=h1coobs[j++];
         }
         else
         {
            l=hidcaobs[i];
            p1=h1caobs[i++];
            p2=h1coobs[j++];
         }
         if ((p1>epsh||p2>epsh)&&(p1<10||p2<10))
         {
            if(ll==l) continue;
            fprintf(gout,"%f %f %7.2lf",p1,p2,ftccobs[k]);
            fprintf(gout," %8.5lf",probchi(ftccobs[k]*ftccobs[k],1),l);
            fprintf(gout," %8.5lf ",ftpcc[k]/npermute);
            revhapid(l,hap);
            for(m=0;m<=selected;m++) fprintf(gout,"%3d",hap[m]);
            fprintf(gout," %d\n",l);
         }
         ll=l;
     }
     fprintf(gout,"\n  case h = haplotype frequency estimates from cases\n");
     fprintf(gout,"contr. h = haplotype frequency estimates from controls\n");
     fprintf(gout,"      ID = haplotype identifier\n");
     break;
   default:;
   }
}
fclose(gout);

return 0;
}

#define NA 0

int shuffle(int *a,int *b,int array_size,int loci)
/*
 generate index for nonmissing genotypes
 20/02/03 JH Zhao UCL
*/
{
int i,j,k,l,a1,a2;

k=0;
for (i=0;i<array_size;i++) a[i]=b[i]=i;
for (i=0;i<array_size;i++)
{
    a1=individuals[i].locus[loci][0];
    a2=individuals[i].locus[loci][1];
    if (a1!=NA&&a2!=NA) b[k++]=i;
}
while (--k>0)
{
/*
  i=(int)floor(rand()/(float)RAND_MAX*n+1);
*/
    i=rand()%k;
    j=b[i];
    b[i]=b[k];
    b[k]=j;
}
k=0;
for (i=0;i<array_size;i++)
{
    a1=individuals[i].locus[loci][0];
    a2=individuals[i].locus[loci][1];
    if (a1!=NA&&a2!=NA) a[i]=b[k++];
}
return 0;
}

#undef NA

float getx2ln(float *ll)
{
FILE *fp;
char line[301],*ptr;
float x2;
fp=fopen("temp.gco","r");
if(fp) while(fgets(line,300,fp))
{
  ptr=strstr(line,"assuming linkage disequilibrium =");
  if(ptr) sscanf(ptr+strlen("assuming linkage disequilibrium ="),"%f",&x2);
  *ll=x2;
  ptr=strstr(line,"chi-square=");
  if (ptr)
  {
     sscanf(ptr+strlen("chi-square="),"%f",&x2);
     goto done;
  }
}
done:
fclose(fp);
return x2;
}

float getln(char *fname)
{
FILE *fp;
char line[301],*ptr;
float ll;
fp=fopen(fname,"r");
if(fp) while(fgets(line,300,fp))
{
  ptr=strstr(line,"assuming linkage disequilibrium =");
  if(ptr)
  {
     sscanf(ptr+strlen("assuming linkage disequilibrium ="),"%f",&ll);
     goto done;
  }
}
done:
fclose(fp);
return ll;
}

float getx2het()
{
FILE *gout;
int k,l,irun;
float llcas,llcon,llcom;

/*combined*/
gout=fopen("temp.gc","w");
for(k=0;k<nloci;++k) if(sel[k]) fprintf(gout,"%2d ",alleles[k]);
fprintf(gout,"\n");
#ifdef TESTGETSIZE
getsize(gout);
#else
for (l=0;l<sample_size;l++)
{
    fprintf(gout,"%d 1",individuals[l].id);
    for (k=0;k<nloci;k++)
        if(sel[k])
        fprintf(gout," %d %d",individuals[l].locus[k][0],
        individuals[l].locus[k][1]);
    fprintf(gout,"\n");
}
#endif
fclose(gout);
irun=system(gcomcmd);
llcom=getln("temp.gco");

/*cases*/
gout=fopen("temp.gc","w");
for(k=0;k<nloci;++k) if(sel[k]) fprintf(gout,"%2d ",alleles[k]);
fprintf(gout,"\n");
for (l=0;l<sample_size;l++)
{
    if(individuals[l].cc!=1) continue;
    fprintf(gout,"%d 1",individuals[l].id);
    for (k=0;k<nloci;k++)
        if(sel[k])
        fprintf(gout," %d %d",individuals[l].locus[k][0],
        individuals[l].locus[k][1]);
    fprintf(gout,"\n");
}
fclose(gout);
irun=system(gcascmd);
llcas=getln("temp1.gco");

/*controls*/
gout=fopen("temp.gc","w");
for(k=0;k<nloci;++k) if(sel[k]) fprintf(gout,"%2d ",alleles[k]);
fprintf(gout,"\n");
for (l=0;l<sample_size;l++)
{
    if(individuals[l].cc==1) continue;
    fprintf(gout,"%d 1",individuals[l].id);
    for (k=0;k<nloci;k++)
        if(sel[k])
        fprintf(gout," %d %d",individuals[l].locus[k][0],
        individuals[l].locus[k][1]);
    fprintf(gout,"\n");
}
fclose(gout);
irun=system(gconcmd);
llcon=getln("temp2.gco");
if (hapwise) gethappy(MARKERDIS,ip);

return 2*(llcas+llcon-llcom);
}

float getx2block(int op,float *ll)
{
FILE *gout;
int k,l,irun;
float x2;
int a1,a2;

gout=fopen("temp.gc","w");
if(!gout)
{
  fprintf(stderr,"I can't open file temp.pgc for output...exit\n");
  return 1;
}
for(k=0;k<nloci;++k)
if(sel[k]) switch(op)
{
   case ALL:
        fprintf(gout,"%2d ",alleles[k]);
        break;
   case BLOCK1:
        if(selp[k]) fprintf(gout,"%2d ",alleles[k]);
        break;
   case BLOCK2:
        if(!selp[k]) fprintf(gout,"%2d ",alleles[k]);
        break;
   default:
        perror("wrong option in getx2block()");
}
fprintf(gout,"\n");
for (l=0;l<sample_size;l++)
{
    fprintf(gout,"%d 1",individuals[l].id);
    for (k=0;k<nloci;k++)
    {
        if (sel[k])
        {
           a1=individuals[l].locus[k][0];
           a2=individuals[l].locus[k][1];
           switch(op)
           {
              case ALL:
                   fprintf(gout," %d %d",a1,a2);
                   break;
              case BLOCK1:
                   if(selp[k]) fprintf(gout," %d %d",a1,a2);
                   break;
              case BLOCK2:
                   if(!selp[k]) fprintf(gout," %d %d",a1,a2);
                   break;
              default:
                   perror("wrong option in getx2block()");
           }
        }
    }
    fprintf(gout,"\n");
}
fclose(gout);
irun=system(gcomcmd);
if (hapwise) gethappy(MARKERMAR,ip);
x2=getx2ln(ll);
return x2;

}

int getlocim(char *locfile)
/*
** retrieves locus information, incoporating missing genotypes
*/
{
FILE *fp;
int i,j,l,l1,l2,m,n[MAX_LOC],ngtype[MAX_LOC],pg[MAX_LOC],npg[MAX_LOC];
double ll,k1,k2;
char line[501],rest[501];
float kp;

fp=fopen(locfile,"r");
if(!fp)
{
  fprintf(stderr,"Error opening %s",locfile);
  exit(1);
}
fgets(line,500,fp);
sscanf(line,"%d %d %d %d",&nloci,&cc,&permute,&npermute);
if(nloci>=MAX_LOC)
{
  perror("Error: maximum number of loci exceeded");
  exit(1);
}
fgets(line,500,fp);
for(i=0;i<nloci;++i,strcpy(line,rest),*rest='\0')
   if(sscanf(line,"%d %[^\n]",&alleles[i],rest)<1) return 0;
fgets(line,500,fp);
sscanf(line,"%d %d",&isgenotype,&iogenotype);
fgets(line,500,fp);
for(i=0;i<nloci;++i,strcpy(line,rest),*rest='\0')
   if(sscanf(line,"%d %[^\n]",&sel[i],rest)<1) return 0;
fgets(line,500,fp);
for(i=0;i<nloci;++i,strcpy(line,rest),*rest='\0')
   if(sscanf(line,"%d %[^\n]",&selp[i],rest)<1) return 0;
if(fgets(line,500,fp)&&sscanf(line,"%f %f %f %f",&freq,&pen0,&pen1,&pen2)==4)
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
    n[l]=alleles[i]+1;
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

int getloci(char *locfile)
/*
** retrieves locus information, ignoring missing genotypes
*/
{
FILE *fp;
int i,j,l,l1,l2,m,n[MAX_LOC],ngtype[MAX_LOC],pg[MAX_LOC],npg[MAX_LOC];
double ll,k1,k2;
char line[241],rest[241];
float kp;

fp=fopen(locfile,"r");
if(!fp)
{
  fprintf(stderr,"Error opening %s",locfile);
  exit(1);
}
for(i=0;i<MAX_LOC;i++) selidx[i]=selndx[i]=selpdx[i]=0;
fgets(line,240,fp);
sscanf(line,"%d %d %d %d",&nloci,&cc,&permute,&npermute);
if(nloci>=MAX_LOC)
{
  perror("Error: maximum number of loci exceeded");
  exit(1);
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
    selidx[l]=i;
    if(selp[i])
    {
      selpdx[l1]=i;
      pg[l1]=ngtype[l];
      ++l1;
    }
    else
    {
      selndx[l2]=i;
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

int noid_new(char *datfile)
/*9APR2001*/
{
FILE *fp;
int i,j,k,l,n,a1,a2,gid;
char line[1000],rest[1000];
list listi;

if((fp=fopen(datfile,"r"))==NULL) fprintf(stderr,"Error opening %s",datfile);
r = NULL;
i=0;
n=0;
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
        p_t.locus[j][0]=a1;
        p_t.locus[j][1]=a2;
      }
      else
      {
        sscanf(line,"%d %d %[^\n]",&a1,&a2,rest);
        if(a1>a2) _swap_(&a1,&a2);
        p_t.locus[j][0]=a1;
        p_t.locus[j][1]=a2;
        if((a1>alleles[j])||(a2>alleles[j]))
        {
          fprintf(stderr,"Error in record %d,",i+1);
          fprintf(stderr,"reset allele number (or <=0 for missing alleles)\n");
          exit(1);
        }
        p_t.gtype[j]=a2g(a1,a2);
      }
      if(sel[j]&&a1<=0) ++k;
   }
   if (iogenotype)
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
   if(k==selected+1)
   {
     ++n;continue;
   }
   if(!cc) p_t.affection=false;
   else if(p_t.affection==true) ++cases;
        else p_t.affection=false;
   ++i;
   listi = (list)malloc(sizeof(struct newrec));
   if(!listi) exit(1);
   listi->id=i;
   listi->cc=p_t.affection;
   for(j=0;j<nloci;j++)
   {
     listi->k[j]=p_t.gtype[j];
     listi->locus[j][0]=p_t.locus[j][0];
     listi->locus[j][1]=p_t.locus[j][1];
   }
   listi->next = r;
   r = listi;
}
fclose(fp);
sample_size=i;
/*
printf("There are %d cases out of %d individuals\n",cases,sample_size);
if(n>0) printf("%d records with no information have been left out \n",n);
*/
l=0;
for(j=0;j<nloci;j++) if(sel[j]) ++l;
digits=l;
r = sort( r, 0 );
listi = r;
while (listi!=NULL) {
#ifdef DEBUG
      fprintf("%5d",listi->id);
      for(j=0;j<nloci;j++) if(sel[j])
         fprintf(" %2d %2d [%2d]", listi->locus[j][0],listi->locus[j][1],
         listi->k[j]);
      fprintf("\n");
#endif
      listi = listi->next;
}
#ifdef DEBUG
r = sort1(r);
#endif
return 0;
}

int a2g(int l1,int l2)
/*
** converts alleles to genotype
*/
{
int lo,up;

lo=l1;
up=l2;
if (l1>l2)
{
   lo=l2;
   up=l1;
}
if (lo==0) return 0;
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

int revhapid(double hapid, int *hap)
{
long int l;
int j,k;

l=hapid-1;
for (j=0;j<=selected;j++)
{
    hap[j]=(int)(l/nalleles[j]);
    if (j==selected) hap[j]=l;
    else l%=nalleles[j];
    hap[j]++;
}

return 0;
}

list sort( list s, int j )
/*4,5,9APR2001*/
{
int i;
list head[M], t;
struct newrec aux;
extern list Last;
if (s==NULL) return(s);
if (s->next == NULL ) {
   Last = s;
   return(s);
}
if (j>=digits) {
   for (Last=s; Last->next!=NULL; Last = Last->next);
   return( s );
   }
for (i=0; i<M; i++) head[i] = NULL;
while (s != NULL) {
      i = s->k[selidx[j]];
      t = s;
      s = s->next;
      t->next = head[i];
      head[i] = t;
      }
t = &aux;
for (i=0; i<M; i++) if (head[i]!=NULL) {
    t->next = sort( head[i], j+1 );
    t = Last;
    }
return(aux.next);
}

int getsize(FILE *gdat)
/*5,9,10APR2001*/
{
typedef struct {int fid,n,nca,nco,locus[maxloci][2];} ids;
int i=0,j,l,k,k1,k2;
int id=0,s,l1[maxloci],l2[maxloci];
ids *ff=0;
list t;

ff=(ids*)malloc(sample_size*sizeof(ids));
if (!ff) {
   perror("error allocating memory in getsize()");
   exit(1);
}
t=r;
k1=k2=0;
k=t->cc;
if(k==CASE) k1=1;
else if(k==CONTROL) k2=1;
for (l=0;l<nloci;l++) l1[l]=l2[l]=t->k[l];
do {
   id++;
   ff[i].fid=t->id;
   for(l=0;l<nloci;l++)
   {
      ff[i].locus[l][0]=t->locus[l][0];
      ff[i].locus[l][1]=t->locus[l][1];
   }
   if(t->next) {
     for(l=0;l<nloci;l++) l2[l]=t->next->k[l];
     k=t->next->cc;
   } else {
     for(l=0;l<nloci;l++) l2[l]=-999;
     k=-999;
   }
   s=0;
   for(l=0;l<nloci;l++) if(sel[l]&&(l1[l]!=l2[l])) s=1;
   if(s==0)
   {
     if(k==CASE) k1++;
     else if(k==CONTROL) k2++;
   }
   if(s>0)
   {
     ff[i].n=k1+k2;
     ff[i].nca=k1;
     ff[i].nco=k2;
     k1=k2=0;
     if(k==CASE) k1=1;
     else if(k==CONTROL) k2=1;
     for (l=0;l<nloci;l++) l1[l]=l2[l];
     i++;
   }
   t=t->next;
} while(t);
id=0;
for (j=0;j<i;j++)
{
    s=0;
    for (l=0;l<nloci;l++)
        if(sel[l]&&(ff[j].locus[l][0]==0||ff[j].locus[l][1]==0)) ++s;
    if(s>selected) continue;
    ++id;
    fprintf(gdat,"%5d %5d",id,ff[j].n);
    if(cc) fprintf(gdat,"%5d %5d",ff[j].nca,ff[j].nco);
    for (l=0;l<nloci;l++) if(sel[l])
        fprintf(gdat," %2d %2d",ff[j].locus[l][0],ff[j].locus[l][1]);
    fprintf(gdat,"\n");
}
free(ff);
return i;
}

void *xmalloc(int len)
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

int gethappy(int op,int instance)
{
FILE *fp;
int i,j,k,l,ll;
double happ,ss,n1,n2,p1,p2;

switch (cc)
{
case MARKERMAR:
     fp=fopen("temp.gco","r");
     if (!fp)
     {
        perror("no file temp.gco");
        exit(1);
     }
     getnhaps(fp,&fobs,&pobs);
     if (instance==0)
     {
        nhaps0=nhaps;
        fobs0=fobs;
        pobs0=pobs;
     }
     geth(fp,COMBINED,instance);
     for (k=0;k<nhaps0;k++) if (instance>0)
     {
         for(j=0;j<nhaps;j++) if(hid[j]==hidobs[k])
         {
            if(fabs(ft[j])>=fabs(ftobs[k])) ftp[k]++;
            goto carry_mm;
         }
     carry_mm:;
     }
     free(ft);
     free(hid);
     break;
case MARKERDIS:
     fp=fopen("temp1.gco","r");
     if (!fp)
     {
        perror("no file temp1.gco");
        exit(1);
     }
     getnhaps(fp,&fobsca,&pobsca);
     nhapcas=nhaps;
     geth(fp,CASES,instance);
     fp=fopen("temp2.gco","r");
     if (!fp)
     {
        perror("no file temp2.gco");
        exit(1);
     }
     getnhaps(fp,&fobsco,&pobsco);
     geth(fp,CONTROLS,instance);
     nhapcon=nhaps;
     /*calculate chi-square statistic*/
     nhaps=nhapcas+nhapcon;
     if (instance==0)
     {
        nhaps0=nhaps;
        nhapcas0=nhapcas;
        nhapcon0=nhapcon;
        fobsca0=fobsca;
        pobsca0=pobsca;
        fobsco0=fobsco;
        pobsco0=pobsco;
        ftccobs=(double*)xmalloc(nhaps*sizeof(double));
        fticcobs=(double*)xmalloc(nhaps*sizeof(double));
        ftpcc=(double*)xmalloc(nhaps*sizeof(double));
        ftssobs=(double*)xmalloc(nhaps*sizeof(double));
        for (j=0;j<nhaps;j++) ftpcc[j]=0;
     }
     ftss=(double*)xmalloc(nhaps*sizeof(double));
     ftcc=(double*)xmalloc(nhaps*sizeof(double));
     fticc=(double*)xmalloc(nhaps*sizeof(double));
     /*x2 statistic*/
     n1=fobsca+pobsca;
     n2=fobsco+pobsco;
     i=j=ll=0;
     for (k=0;k<nhaps;k++)
     {
         if (i<nhapcas||j<nhapcon)
         if (hidca[i]<hidco[j])
         {
            ss=CASES;
            l=hidca[i];
            p1=h1cas[i++];
            p2=0;
         }
         else if (hidca[i]>hidco[j])
         {
            ss=CONTROLS;
            l=hidco[j];
            p1=0;
            p2=h1con[j++];
         }
         else
         {
            ss=COMBINED;
            l=hidca[i];
            p1=h1cas[i++];
            p2=h1con[j++];
         }
         if ((p1>epsh||p2>epsh)&&(p1<10||p2<10))
         {
            if(ll==l) continue;
            happ=(n1*p1+n2*p2)/(n1+n2);
            happ=(p1-p2)/sqrt(happ*(1-happ)*(0.5/n1+0.5/n2));
            ftcc[k]=happ;
            fticc[k]=i;
            ftss[k]=ss;
            if (instance==0)
            {
               ftssobs[k]=ftss[k];
               ftccobs[k]=ftcc[k];
               fticcobs[k]=fticc[k];
            }
         }
         ll=l;
     }
     for (k=0;k<nhaps0;k++) if (instance>0)
     {
         for(j=0;j<nhaps;j++)
         {
            if (fticcobs[k]!=fticc[j]||ftssobs[k]!=ftss[j]) continue;
            if (fabs(ftcc[j])>=fabs(ftccobs[k])) ftpcc[k]++;
            goto carry_md;
         }
     carry_md:;
     }
     free(hidco);
     free(hidca);
     free(h1con);
     free(h1cas);
     free(ftcc);
     free(fticc);
     break;
default:;
}

return 0;
}

void getnhaps(FILE *fp,int *f, int *p)
/*18-03-2003 in shape*/
{
int i;
char line[256], rest[256], *pstr;

if (fp) while (fgets(line,255,fp))
{
   pstr=strstr(line,"with F/P");
   if (pstr) sscanf(line,"%d/%d",&*f,&*p);
   pstr=strstr(line,"nonzero haplotypes =");
   if (pstr)
   {
      sscanf(pstr+strlen("nonzero haplotypes ="),"%d",&nhaps);
      goto restart;
   }
}
restart:
rewind(fp);
while (fgets(line,255,fp))
{
      pstr=strstr(line,"h1 = frequencies under linkage disequilibrium (");
      if(pstr)
      sscanf(pstr+strlen("h1 = frequencies under linkage disequilibrium ("),"&lf",&pl);
      pstr=strstr(line,"FT = sqrt");
      if (pstr) goto ok;
}
ok:
for (i=0; i<3; i++) fgets(line,255,fp);
}

int geth(FILE *fp, int op, int instance)
{
int i,j,k;
char line[256], rest[256], *pstr;
double x,y,z;

if (instance==0) switch(op)
{
   case CONTROLS:
        h1coobs=(double*)xmalloc((nhaps+1)*sizeof(double));
        hidcoobs=(double*)xmalloc((nhaps+1)*sizeof(double));
        h1coobs[nhaps]=10;
        hidcoobs[nhaps]=ULONG_MAX;
        break;
   case CASES:
        h1caobs=(double*)xmalloc((nhaps+1)*sizeof(double));
        hidcaobs=(double*)xmalloc((nhaps+1)*sizeof(double));
        h1caobs[nhaps]=10;
        hidcaobs[nhaps]=ULONG_MAX;
        break;
   case COMBINED:
        h0=(double*)xmalloc(nhaps*sizeof(double));
        h1=(double*)xmalloc(nhaps*sizeof(double));
        hidobs=(double*)xmalloc(nhaps*sizeof(double));
        ftobs=(double*)xmalloc(nhaps*sizeof(double));
        ftp=(double*)xmalloc(nhaps*sizeof(double));
        for (j=0;j<nhaps;j++) ftp[j]=0;
        break;
   default:;
}
switch (op)
{
       case CONTROLS:
            h1con=(double*)xmalloc((nhaps+1)*sizeof(double));
            hidco=(double*)xmalloc((nhaps+1)*sizeof(double));
            h1con[nhaps]=10;
            hidco[nhaps]=ULONG_MAX;
            break;
       case CASES:
            h1cas=(double*)xmalloc((nhaps+1)*sizeof(double));
            hidca=(double*)xmalloc((nhaps+1)*sizeof(double));
            h1cas[nhaps]=10;
            hidca[nhaps]=ULONG_MAX;
            break;
       case COMBINED:
            ft=(double*)xmalloc(nhaps*sizeof(double));
            hid=(double*)xmalloc(nhaps*sizeof(double));
            break;
}
i=0;
while (fgets(line,255,fp)&&strlen(line)>1)
{
      sscanf(line,"%*lf %lf %lf %lf %[^\n]",&x,&y,&z,rest);
      strcpy(line,rest);
      for (j=0;j<nloci;j++)
      {
          sscanf(line,"%d %[^\n]",&k,rest);
          strcpy(line,rest);
      }
      sscanf(line,"%d",&k);
      if (instance==0) switch(op)
      {
         case CONTROLS:
              h1coobs[i]=x;
              hidcoobs[i]=k;
              break;
         case CASES:
              h1caobs[i]=x;
              hidcaobs[i]=k;
              break;
         case COMBINED:
              h0[i]=y;
              h1[i]=x;
              ftobs[i]=z;
              hidobs[i]=k;
              break;
         default:;
      }
      switch (op)
      {
           case CONTROLS:
                h1con[i]=x;
                hidco[i]=k;
                break;
           case CASES:
                h1cas[i]=x;
                hidca[i]=k;
                break;
           case COMBINED:
                ft[i]=z;
                hid[i]=k;
                break;
           default:;
      }
      ++i;
}
fclose(fp);
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

/*
  History:

  21-2-2003 start implementation (UCL)
  07-3-2003 incorporate pgc code
  11-3-2003 in shape
  12-3-2003 works for case-control and whole block analysis
  13-3-2003 improve I/O and blockwise statistics
  14-3-2003 add random number seed
  20-3-2003 add gethappy
  21-3-2003 intensive check (interface,%d/%lf,merge)
  23-3-2003 pass check with HLA data
  29-3-2003 fix chi-square in case-control option and use normal deviate
  24-5-2003 add haplotype based reverse mixed-radix number function
  31-5-2003 revhap as of 24-5-2003 functioning

*/

