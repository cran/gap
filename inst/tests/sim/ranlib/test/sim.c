/******************************************************/
/***J.H. Zhao & Sham P.C. 03MAR97 I.O.P.            ***/
/***see companion documentation for more information***/
/******************************************************/
#include <stdio.h>
#include <math.h>
#include <ranlib.h>

#define MAXLOCI 20
#define MAXALLELES 10
#define MAXFAM 500
#define MAXSIBS 9
#define MAXIND 30
#define MAXALL MAXFAM*MAXIND

#define MALE 1
#define FEMALE 2
#define MU 4
#define NFAM 150
#define NLOCI 5
#define NTRAIT 1
#define NMG 2
#define NPG 2
#define NCE 2
#define NUE 2
#define NDISEQ 1

typedef struct{
  int pid, father, mother, nsibs, sib[MAXSIBS];
} NUC_FAMILY;
typedef struct{
  int pid, id, fid, mid, sex, gtype[MAXLOCI][2];
  float trait[NTRAIT],pg[NTRAIT][NPG],ce[NTRAIT][NCE];
} PERSON;
typedef struct{
  char name;
  int type;
  int nall;
  float freq[MAXALLELES];
  float cumf[MAXALLELES];
} LOCUS;

float r[NLOCI-1]={0.02,0.5,0.0,0.001};
int mg[NMG]={0,1};
float dom[NMG]={1,0.5};
float beta[NTRAIT][NMG+NPG+NCE+NUE]={0.2,0.4,0,0,0,0,0,0};
float kce[NTRAIT][NCE]={0.2,0.2};
int cel[NCE]={1,0};
float kp, pen[3]={0, 0, 1};

/***In this implementation, linkage disequilibrium is allowed.

  diseq[] is indicators for the two biallelic loci in disequilibrium;
  diseqt[] is their marginal/allele frequencies. One of them might be
  disease locus itself. 

  In this case quantitative traits may or may not be unused.

haplotype frequencies  	marginal	conditional
	.2	.1	.3		.2/.3, 1.0
  	.3 	.4	.7		.3/.7, 1.0
  	.5 	.5  	1.0

  For locus with more than 2 loci, haplotype frequencies can be obtained 
  from separate SAS program, as diseqhap below. A comparison with those
  from Weir, B.S. (1990) is pending.

  In the example above, cell(1,2) has linkage disequilibrium coefficient 
  1.0 (complete).

***/
int diseq[NDISEQ]={2};
float diseqt[NDISEQ+1]={0.3,0.5};
/*get haplotype frequencies*/
float diseqhap[2][2]={0.2,0.1,0.3,0.4};
float diseqcum[2][2]={0.2/0.3,1.0,0.3/0.7,1.0};
/*programming accumulation for more allleles pending*/
/*read more than one set of parameters from file*/
char *diseqfile="abc";

main(int argc, char *argv[])
{
FILE *fi,*fo;
float cump[MAXSIBS],tp[MAXSIBS];
static NUC_FAMILY nuc_fam[MAXFAM];
static LOCUS locus[MAXLOCI];
static PERSON person[MAXALL];
static float diseqhap[MAXALLELES][MAXALLELES],diseqcum[MAXALLELES][MAXALLELES];
int pid, id, fid, mid, nid, pbnd, sex, aff, fsize[MAXFAM];
int i,j,k,k1,k2,l,l1,l2,m,m1,m2,type;
float p, s, x, y, mu, pden;
int tl;
float f, g, x1, x2;
float q, psave,t, d, z, e, s1, s2, s3[MAXSIBS];
char line[241];

/*set haplotype with disequilibrium parameters*/
/*specify the natural set number*/
l=0;
if(argc>1) l=atoi(argv[1]);
fi=fopen(diseqfile,"r");
for(i=1;i<l;++i) for (j=0;j<3;++j) fgets(line,240,fi);
printf("allele frequency, penetrance and haplotypes\n");
if(fi==NULL) {printf("error reading disequilibrium parameter file\n");exit(0);}
fgets(line,240,fi);
sscanf(line,"%f",&diseqt[0]);diseqt[1]=0.13;
fgets(line,240,fi);
sscanf(line,"%f %f %f",&pen[0],&pen[1],&pen[2]);
fgets(line,240,fi);
sscanf(line,"%f %f %f %f",&diseqhap[0][0],&diseqhap[0][1],&diseqhap[1][0],&diseqhap[1][1]);
diseqcum[0][0]=diseqhap[0][0]/(diseqhap[0][0]+diseqhap[0][1]);
diseqcum[0][1]=1.0;
diseqcum[1][0]=diseqhap[1][0]/(diseqhap[1][0]+diseqhap[1][1]);
diseqcum[1][1]=1.0;
fclose(fi);
printf("%f\n",diseqt[0]);
printf("%f %f %f\n",pen[0],pen[1],pen[2]);
printf("%f %f %f %f\n",diseqhap[0][0],diseqhap[0][1],diseqhap[1][0],diseqhap[1][1]);

printf("\nkp control\n");
q=diseqt[0];p=1.0-q;
kp=p*p*pen[0]+2.0*p*q*pen[1]+q*q*pen[2];
printf("The calculated prevalence Kp=%f\n",kp);
/*suppose Kp and penetrances are defined instead of allele frequencies*/
s1=pen[0]-2.0*pen[1]+pen[2];
s2=pen[1]*pen[1]-pen[0]*pen[2]+s1*kp;
if(s1==0.0)
  psave=(kp-pen[2])/(2.0*(pen[1]-pen[2]));
else {
  if(s2<0) printf("Error, allele frequencies are imaginary numbers\n");
  x1=(-pen[1]+pen[2]+sqrt(s2))/s1;
  x2=(-pen[1]+pen[2]-sqrt(s2))/s1;
  printf("solution from Kp %f %f\n",x1,x2);
  if((x1<0&&x2<0)||(x1>1&&x2>1)) printf("Error in calculating frequency");
  if(x1>0&&x1<1.0) psave=x1; else psave=x2;
}
p=psave;q=1.0-p;
printf("The penetrance: %.4f %.4f %.4f\n",pen[0],pen[1],pen[2]);
printf("The allele frequency is %f\n",q);
printf("\n");

printf("10 random numbers (uniform, normal):\n");
for(i=0;i<10;i++) printf("%8.5f %8.5f\n",ranf(),snorm());
printf("\n");

mu=MU;p=exp(-mu);pden=1-p;
s=0.0;
for(i=1;i<MAXSIBS;i++) {
   p*=(mu/i);
   tp[i]=p/pden;
   s+=tp[i];
} pden=s;

s=0.0;
cump[0]=0.0;
for(i=1;i<MAXSIBS;i++) {
   tp[i]/=pden;
   s+=tp[i];
   cump[i]=s;
}

printf("Truncated Poisson probabilities:\n");
for(i=1;i<MAXSIBS;i++) {printf("%f %f\n",tp[i],cump[i]);}
printf("\n");

/*Set allele frequencies, rather ad hoc*/
l=0;k=0;
for (i=0;i<NLOCI;++i) {
    locus[i].type=3;locus[i].nall=2;
    if(i==mg[k++])
      for(j=0;j<locus[i].nall;++j) locus[i].freq[j]=1.0/locus[i].nall;
    else if(i==diseq[l]){
      if(l==0) locus[i].type=1;
      locus[i].freq[0]=diseqt[l];
      locus[i].freq[1]=1.0-locus[i].freq[0];
      ++l;
      locus[i+1].freq[0]=diseqt[l];
      locus[i+1].freq[1]=1.0-diseqt[l];
      } else {
         if(i==4) {
           locus[i].nall=2;
           for(j=0;j<locus[i].nall;++j) locus[i].freq[j]=1.0/locus[i].nall;
         }
      }
}

printf("Cumulative allele fractions/locus type:\n");
for(i=0;i<NLOCI;i++){
   s=0;locus[i].cumf[0]=0;
   for(j=0;j<locus[i].nall;j++){
      s+=locus[i].freq[j];
      locus[i].cumf[j]=s;
      printf("%f ",s);
   } printf(" %d\n",locus[i].type);
}
printf("\n");

printf("Recombination fractions:\n");
for(i=0;i<NLOCI-1;i++) printf("%8.5f",r[i]);printf("\n\n");

printf("Nuclear families\n");
id=0;
for(i=0;i<NFAM;i++){ /*get families first, 29/5/97*/
   pid=i+1;
   id++;fid=id;
   nuc_fam[i].pid=pid;
   person[id].pid=pid;
   person[id].id=id;
   person[id].fid=0;
   person[id].mid=0;
   person[id].sex=MALE;
   nuc_fam[i].father=fid;
   id++;mid=id;
   person[id].pid=pid;
   person[id].id=id;
   person[id].fid=0;
   person[id].mid=0;
   person[id].sex=FEMALE;
   nuc_fam[i].mother=mid;
   for(j=0;j<MAXSIBS;j++) if (ranf()>=cump[j]) nuc_fam[i].nsibs=j;
   for(j=0;j<nuc_fam[i].nsibs;++j){
      id++;
      sex=ranf()<0.5?MALE:FEMALE;
      nuc_fam[i].sib[j]=id;
      person[id].pid=pid;
      person[id].id=id;
      person[id].fid=fid;
      person[id].mid=mid;
      person[id].sex=sex;
   }
}

printf("The family size\nFamily ID BOID Size\n");
id=0;k=1;i=0;
for(j=0;j<MAXALL;++j){
   id++;
   if(id==1) 
     l1=person[id].pid;
   else {
     l2=person[id].pid;
     if(l1!=l2&&l1!=0) {/*delete extra printout, 5/8/97*/
       fsize[++i]=k;
       printf("%4d %4d %4d %4d\n",l1,l2,j,k);
       l1=l2;k=1;}
     else ++k;
   }
}

id=0;
for(i=1;i<=NFAM;++i){
   for(l=0;l<fsize[i];++l){
      id++;
      pid=person[id].pid;
      fid=person[id].fid;
      mid=person[id].mid;
      sex=person[id].sex;
      printf("%3d %3d %3d %3d %3d ",pid,id,fid,mid,sex);
      if(fid==0&&mid==0){
        for(j=0;j<NLOCI;++j){
           l1=1;l2=1;
           f=ranf();g=ranf();
           for(k=0;k<locus[j].nall;++k) if (f>=locus[j].cumf[k]) l1++;
           for(k=0;k<locus[j].nall;++k) if (g>=locus[j].cumf[k]) l2++;
           person[id].gtype[j][0]=l1;
           person[id].gtype[j][1]=l2;
        }
/*incorporate linkage disequilibrium*/
        k1=0;k2=0;
        for(j=0;j<NLOCI;++j){
           if(j==diseq[k1]) {if(k1==0) k2=j;++k1;}
           if(j==k2+1) {/*This is the second locus in disequilibrium*/
             l1=1;l2=1;
             m1=person[id].gtype[k2][0]-1;
             f=ranf();
             for(k=0;k<locus[k1].nall;++k) if(f>=diseqcum[m1][k]) l1++;
             m2=person[id].gtype[k2][1]-1;
             g=ranf();
             for(k=0;k<locus[k1].nall;++k) if(g>=diseqcum[m2][k]) l2++;
             person[id].gtype[j][0]=l1;
             person[id].gtype[j][1]=l2;
           }
        }
      }
      else {
        f=ranf();g=ranf();
        l1=person[fid].gtype[0][0];
        if(f>0.5) l1=person[fid].gtype[0][1];
        l2=person[mid].gtype[0][0];
        if(g>0.5) l2=person[mid].gtype[0][1];
        person[id].gtype[0][0]=l1;
        person[id].gtype[0][1]=l2;
        for(k=1;k<NLOCI;++k){
           if (k==1) l1=(f>0.5)?1:0;
           l1=(ranf()<=r[k-1])?1-l1:l1;
           if (k==1) l2=(g>0.5)?1:0;
           l2=(ranf()<=r[k-1])?1-l2:l2;
           person[id].gtype[k][0]=person[fid].gtype[k][l1];
           person[id].gtype[k][1]=person[mid].gtype[k][l2];
        }
      }
      for(k=0;k<NLOCI;++k)printf("  %3d%3d",person[id].gtype[k][0],person[id].gtype[k][1]);
      printf("\n");
   }
}

/*get phenotypes*/
id=0;
for(j=1;j<=NFAM;++j){
   for(l=0;l<fsize[j];++l){
      id++;
      pid=person[id].pid;
      fid=person[id].fid;
      mid=person[id].mid;
      sex=person[id].sex;
      for(k=0;k<NTRAIT;k++) {
         s1=0;
/*Major genes*/
         for(i=0;i<NMG;++i){
            tl=mg[i];
            q=locus[tl].freq[0];p=1-q;
            d=dom[i];
            s=pow(p*q*(q+2*p*d),2)+2*p*q*pow(d-q*q-2*p*q*d,2)+pow(q*(1-q*q-2*p*q*d),2);
            t=sqrt(1/s);
            z=-(pow(q,2)*t+2.0*q*(1-q)*t*d);
/*
            d=0 recessive, d=0.5 additive, d=1 dominant; 
            mean=0, variance=1.
            s=z*pow(p,2)+2*p*q*(z+d*t)+pow(q,2)*(z+t);
            e=pow(p*z,2)+2*p*q*pow(z+d*t,2)+pow(q*(z+t),2);
            printf("Mean/Variance %f %f\n",s,e);
*/
            l1=person[id].gtype[tl][0];
            l2=person[id].gtype[tl][1];
            type=1;
            if((l1==1&&l2==2)||(l1==2&&l2==1)) type=2;
            if(l1==2&&l2==2) type=3;
            d=(type==1)?0:(type==2)?d:1;
            g=z+d*t;
            s1+=g*beta[k][i];
         }
         person[id].trait[k]=s1;
/*polygene*/
         for(i=0;i<NPG;++i){
            if(fid==0&&mid==0) person[id].pg[k][i]=snorm();
            else person[id].pg[k][i]=snorm()/sqrt(2)+(person[fid].pg[k][i]+person[mid].pg[k][i]);
         }
         s1=0;
         for(i=0;i<NPG;++i) s1+=person[id].pg[k][i]*beta[k][NMG+i];
         person[id].trait[k]+=s1;
/*common environment*/
         if(l==0) for(i=0;i<NCE;++i) 
           if(fid==0&&mid==0) {person[id].ce[k][i]=snorm();x=snorm();}
           else person[id].ce[k][i]=kce[k][i]*(person[fid].ce[k][i]+person[mid].ce[k][i])
                                   +sqrt(1-2*pow(kce[k][i],2))*x;
         s1=0;
         for(i=0;i<NCE;++i)
            s1+=person[id].ce[k][cel[i]]*beta[k][NMG+NPG+cel[i]];
         person[id].trait[k]+=s1;
/*unique environment*/
         person[id].trait[k]+=snorm()*beta[k][NMG+NPG+NCE+k];
      }
   }
}

fo=fopen("new1","w");
if(fo==NULL) printf("error openning output\n");
l=0;
for(j=0;j<NLOCI;++j) if(j==mg[l]) {locus[j].type=4;l++;}
id=0;
for(i=1;i<=NFAM;++i)
   for(l=0;l<fsize[i];++l){
      id++;
      pid=person[id].pid;
      fid=person[id].fid;
      mid=person[id].mid;
      sex=person[id].sex;
      fprintf(fo,"%4d %4d %4d %4d %4d",pid,id,fid,mid,sex);
/*    for(k=0;k<NTRAIT;k++) fprintf(fo," %12.5f",person[id].trait[k]);
*/    k=0;
      for(j=0;j<NLOCI;++j) {
         l1=person[id].gtype[j][0];l2=person[id].gtype[j][1];
         m=l1+l2;
         switch(locus[j].type) {
         case 1: /*invariant with disease models*/
           f=ranf();
           switch(m) {
           case 2:aff=1;if(f<pen[0]) aff=2;break;
           case 3:aff=1;if(f<pen[1]) aff=2;break;
           case 4:aff=1;if(f<pen[2]) aff=2;break;
           default: break;
           }
           printf("family %d %d (%d %d %d) is %f aff=%d\n",pid,id,l1,l2,m,f,aff);
           fprintf(fo,"%3d",aff);
           break;
         case 3: 
           fprintf(fo," %3d%3d",l1,l2);
           break;
         case 4: 
/*print MG separately but now do nothing */
           ++k;break;
         default: 
           break;
         }
      } fprintf(fo,"\n");
} fclose(fo);
return 0;
}
