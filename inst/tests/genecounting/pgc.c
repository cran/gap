/*
  Program for preparing haplotype frequencer files
  JH Zhao 13/06/2001, 19/04/2002
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "pgc.h"

#define version 1.0

/*Binary Search Trees*/

node *itree(node *r,double genid)
/*
** insert and sort
*/
{
int i,j;
if (!r)
{
  r=malloc(sizeof(node));
  if (!r)
  {
    printf("out of memory\n");exit(0);
  }
  r->left=r->right=NULL;
  r->genid=genid;
  r->nca=r->nco=0;
  if (p_t.affection) r->nca++;
  else r->nco++;
  j=0;
  for (i=0;i<nloci;i++)
  {
      if (!sel[i]) r->l[i]=r->u[i]=0;
      else {
           r->l[j]=p_t.locus[i][0];
           r->u[j]=p_t.locus[i][1];
           j++;
      }
  }
}
else
if (genid<r->genid) r->left=itree(r->left,genid);
else if (genid>r->genid) r->right=itree(r->right,genid);
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
  if (!t) break;
}
return t;
}

node *dtree(node *t,double key)
/*
** delete
*/
{
node *p,*p2;

if(!t) return t;
if (t->genid==key)
{
  if(t->left==t->right)
  {
    free(t);
    return 0;
  }
  else if (!t->left)
  {
    p=t->right;
    free(t);
    return p;
  }
  else if (!t->right)
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
/*
** left subtree->t->right subtree
*/
{
if (!t) return;
inorder(t->left);
printf("%d ",t->genid);
inorder(t->right);
}

void preorder(node *t)
/*
** t->left subtree->right subtree
*/
{
if (!t) return;
printf("%d ",t->genid);
preorder(t->left);
preorder(t->right);
}

void postorder(node *t)
/*
** left subtree->right subtree->t
*/
{
if (!t) return;
postorder(t->left);
postorder(t->right);
printf("%d ",t->genid);
}

void ptree(node *r,int l,FILE *gdat)
/*
** print tree inorder
*/
{
int i,j;
if (!r) return;
ptree(r->left,l+1,gdat);
fprintf(gdat,"%20.0lf %4d",r->genid,r->nca+r->nco);
if (cc) fprintf(gdat," %4d %4d",r->nca,r->nco);
j=0;
for (i=0;i<nloci;i++)
{
    if (!sel[i]) continue;
    fprintf(gdat,"%3d%3d",r->l[j],r->u[j]);
    j++;
}
fprintf(gdat,"\n");
ptree(r->right,l+1,gdat);
}

void rtree(node *t)
/*
** left subtree->right subtree->t
*/
{
if (!t) return;
rtree(t->left);
rtree(t->right);
free(t);
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

int getdat(char *datfile)
/*
** retrieves data file, ignoring missing genotypes
*/
{
FILE *fp;
int i,j,k,l,n,a1,a2,genotype[MAX_LOC],gid;
char line[1000],rest[1000];
double p;

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
   l=0;
   for(j=0;j<nloci;++j)
   {
     if(!sel[j]) continue;
     genotype[l]=p_t.gtype[j];++l;
   }
   p=position(selected,genotype,0);
   if (!rt) rt=itree(rt,p);
   else itree(rt,p);
   ++i;
}
fclose(fp);
sample_size=i;
printf("There are %d cases out of %d individuals\n",cases,sample_size);
if(n>0) printf("%d records with partial information have been left out \n",n);
return 0;
}

int a2g(int l1,int l2)
/*
** converts alleles to genotype
*/
{
int lo,up;
lo=l1; up=l2;
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

double position(int m,int *genotype,int op)
/*
** (loc1-1)*PROD n[2,...,m]+(loc2-1)*PROD n[3,...,m]+...+locm
*/
{
int l;
double pos,sum;
pos=sum=0;
for(l=0;l<m+1;++l) if(genotype[l]==0) return(pos);
switch(op)
{
case 0:
  for(l=0;l<m;++l) sum+=(genotype[l]-1)*nall[l+1];
  break;
case 1:
  for(l=0;l<m;++l) sum+=(genotype[l]-1)*np[l+1];
  break;
case 2:
  for(l=0;l<m;++l) sum+=(genotype[l]-1)*nnp[l+1];
  break;
default:
  break;
}
pos=sum+genotype[m];
return(pos);
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

int getdatm(char *datfile)
/*
** retrieves data file, incorporating missing genotypes
*/
{
FILE *fp;
int i,j,k,l,n,a1,a2,genotype[MAX_LOC],gid;
char line[1000],rest[1000];
double p;

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
        p_t.locus[j][0]=a1;
        p_t.locus[j][1]=a2;
      }
      else
      {
        sscanf(line,"%d %d %[^\n]",&a1,&a2,rest);
        if(a1>a2) _swap_(&a1,&a2);
        if((a1>alleles[j])||(a2>alleles[j]))
        {
          fprintf(stderr,"Error in record %d,",i+1);
          fprintf(stderr,"reset allele number (or <=0 for missing alleles)\n");
          exit(1);
        }
        p_t.locus[j][0]=a1;
        p_t.locus[j][1]=a2;
      }
      if(sel[j]&&a1<=0) ++k;
      if(a1<=0) a1=alleles[j]+1;
      if(a2<=0) a2=alleles[j]+1;
      p_t.gtype[j]=a2g(a1,a2);
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
   l=0;
   for(j=0;j<nloci;++j)
   {
     if(!sel[j]) continue;
     genotype[l]=p_t.gtype[j];++l;
   }
   p=positionm(selected,genotype,0);
   if (!rt) rt=itree(rt,p);
   else itree(rt,p);
   ++i;
}
fclose(fp);
sample_size=i;
printf("There are %d cases out of %d individuals\n",cases,sample_size);
if(n>0) printf("%d records with no information have been left out \n",n);
return 0;
}

double positionm(int m,int *genotype,int op)
/*
** (loc1-1)*PROD n[2,...,m]+(loc2-1)*PROD n[3,...,m]+...+locm
*/
{
int l;
double pos,sum;
pos=sum=0;
switch(op)
{
case 0:
  for(l=0;l<m;++l) sum+=(genotype[l]-1)*nall[l+1];
  break;
case 1:
  for(l=0;l<m;++l) sum+=(genotype[l]-1)*np[l+1];
  break;
case 2:
  for(l=0;l<m;++l) sum+=(genotype[l]-1)*nnp[l+1];
  break;
default:
  break;
}
pos=sum+genotype[m];
return(pos);
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

list sort1( list r )
{
list head[M], tail[M];
int i, j, h;
for (i=digits-1; i>0; i--) {
    for (j=0; j<M; j++) head[j] = NULL;
    while ( r != NULL ) {
        h = r->k[selidx[i]];
        if ( head[h]==NULL ) head[h] = r;
        else tail[h]->next = r;
        tail[h] = r;
        r = r->next;
        };
    r = NULL;
    for (j=M-1; j>=0; j--)
        if ( head[j] != NULL ) {
           tail[j]->next = r;
           r = head[j];
           }
    };
return( r );
};

int getsize(FILE *gdat)
/*5,9,10APR2001*/
{
int i=0,j,l,k,k1,k2;
int id=0,s,l1[maxloci],l2[maxloci];
typedef struct {int fid,n,nca,nco,locus[maxloci][2];} ids;
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

int noid(char *datfile, FILE *gdat)
/*9APR2001*/
{
FILE *fp;
int i,j,k,l,n,a1,a2,gid;
char line[1000],rest[1000];
double p;
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
printf("There are %d cases out of %d individuals\n",cases,sample_size);
if(n>0) printf("%d records with no information have been left out \n",n);

l=0;
for(j=0;j<nloci;j++) if(sel[j]) ++l;
digits=l;
r = sort( r, 0 );
listi = r;
while (listi!=NULL) {
#ifdef DEBUG
      fprintf(gdat,"%5d",listi->id);
      for(j=0;j<nloci;j++) if(sel[j])
         fprintf(gdat," %2d %2d [%2d]", listi->locus[j][0],listi->locus[j][1],
         listi->k[j]);
      fprintf(gdat,"\n");
#endif
      listi = listi->next;
}
p=getsize(gdat);
printf("There are %.0lf observed multilocus genotypes\n",p);
#ifdef DEBUG
r = sort1(r);
#endif
return 0;
}

int main(int argc,char **argv)
{
FILE *gout;
int j;
time_t t;

printf("Data preparation for GENECOUNTING %.1f JH Zhao 9-7-2002\n",version);
#ifdef NOHM
  printf("(Omit individuals with missing genotypes)\n\n");
#elif ID
  printf("(Collect all information via genotype identifier)\n\n");
#else
  printf("(Collapse over genotypes)\n\n");
#endif
time(&t);
printf("%s\n",ctime(&t));
printf("Maximum number of loci = %d\n",maxloci);
printf("Maximum number of alleles = %d\n",maxalleles);

if(argc<4)
{
  fprintf(stderr,"\nUsage: %s parfile datfile outfile\n",argv[0]);
  fprintf(stderr,"\nwhere parfile/datfile are EHPLUS parameter/data files\n");
  fprintf(stderr,"\noutfile is the target file with condensed information\n");
  exit(1);
}
else
{
  outfile=argv[3];
  gout=fopen(outfile,"w");
  if(!gout)
  {
    fprintf(stderr,"I can't open file %s for output...exit\n",outfile);
    return 1;
  }
#ifdef NOHM /*use only individuals with complete information*/
  if (!getloci(argv[1]))
  {
     getdat(argv[2]);
     for(j=0;j<nloci;++j) if(sel[j]) fprintf(gout,"%2d ",alleles[j]);
     fprintf(gout,"\n");
     ptree(rt,0,gout);
  }
#elif ID /*missing data but with genotype ID*/
  if (!getlocim(argv[1]))
  {
     getdatm(argv[2]);
     for(j=0;j<nloci;++j) if(sel[j]) fprintf(gout,"%2d ",alleles[j]);
     fprintf(gout,"\n");
     ptree(rt,0,gout);
  }
#else /*default, no genotype ID*/
  if (!getloci(argv[1]))
  {
     for(j=0;j<nloci;++j) if(sel[j]) fprintf(gout,"%2d ",alleles[j]);
     fprintf(gout,"\n");
     noid(argv[2],gout);
  }
#endif
  fclose(gout);
  return 0;
}
}
/*
  Compiling options:

    -DNOHM ordinary EM with multilocus genotype identifier
    -DID EM with missing data by multilocus genotype identifier
    default (EM with missing data. no identifier is necessary but slower)
    -mn for Symantec C++

  Possible extensions:
    permutation on oberved genotypes
    housekeeping for various statistics
*/
