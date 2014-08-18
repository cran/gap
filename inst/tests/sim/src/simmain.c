/***J.H. Zhao & Sham P.C. 03MAR97 I.O.P.            ***/
/******************************************************
 SIM  MAIN PROGRAM
*******************************************************/
#include "sim.h"

/*globals*/
char locfile[50]="samples/loc.tst",
     pedfile[50]="samples/ped.tst",
     confile[50]="samples/problem.dat",
     rawped[50]="sim.out",
     selped[50]="sel.out";

/*internals*/
int ld_sim(int, char**);

main(int argc, char *argv[])
{
ld_sim(argc, argv);
}

int ld_sim(int argc, char **argv)
{
FILE *fi,*fo;
struct truncated_poisson {float mu,p[MAXSIBS],cump[MAXSIBS];} tp;
float kp,c[4],f0,f1,f2;
float p,q,s,s2,pden,ax,delta,f;
int nreps;
int pid,id,fid,mid,offs,nfs,nms,pbnd,sex,aff,liab,nsibs;
int i,j,k,k1,k2,l,l1,l2,m,m1,m2,skiplines=0;
char line[2000],rest[2000];
long seed=3000,seed1,seed2;

/*command line parameters, expanded ?*/
if (argc>1) strcpy(locfile,argv[1]);
if (argc>2) strcpy(pedfile,argv[2]);
if (argc>3) strcpy(confile,argv[3]);
if (argc>4) strcpy(rawped,argv[4]); 

/*allocate working arrays*/
person=(PERSON *)malloc(MAXIND*sizeof(PERSON));
locus=(LOCUS *)malloc(MAXLOCI*sizeof(LOCUS));
typed=imatrix(0,MAXIND,0,MAXLOCI);
locus_order=ivector(0,MAXLOCI);
mg=ivector(0,NMG);
r=vector(0,MAXLOCI-1);
beta=matrix(0,NTRAIT,0,NMG+NPG+NCE+NUE);
kce=matrix(0,NTRAIT,0,NCE);
fsize=ivector(0,MAXFAM);

/*locus file*/
fi=fopen(locfile,"r");
if(fi==NULL) perror("error opening locus file");
locprep(fi);
fclose(fi);

/*control file*/
fi=fopen(confile,"r");
fgets(line,2000,fi);
sscanf(line,"%d",&nreps);
fgets(line,2000,fi);
sscanf(line,"%d %d",&seed1,&seed2);
setall(123,321);
setsd(seed1,seed2);
fclose(fi);

getsd(&seed1,&seed2);
#ifdef DEBUG
/*check random number seeds*/
printf("The current seeds: %d %d\n",seed1,seed2);
#endif

/*get family sizes into fsize*/
nfam=getsize(pedfile);

/*prepare for individual family simulation*/
fo=fopen(rawped,"w");

nobody.pid=0;
nobody.id=0;
nobody.fid=0;
nobody.mid=0;
nobody.sex=UNKNOWN;
nobody.simmed=false;
selected=0;
for (l=0;l<nreps;l++) {
/*doesn't need rewind*/
  fi=fopen(pedfile,"r");
  if (fi==NULL) printf("error opening pedigree file\n");
  for (i=0;i<nfam;++i) {
    ++selected;
    pedsize=fsize[i];
    for (l1=0;l1<MAXIND;l1++) person[l1]=nobody;
    for (k=0;k<5;k++) doubled[k][0]=doubled[k][1]=0;
    for (j=1;j<=pedsize;++j) {
      for (l2=0;l2<MAXLOCI;l2++) typed[l1][l2]=1;
      fgets(line,1000,fi);
#ifdef USEMAKEPED
      sscanf(line,"%d %d %d %d %d %d %d %d %d %[^\n]",
            &pid,&id,&fid,&mid,&offs,&nfs,&nms,&sex,&pbnd,rest);
      person[j].offs=offs;
      person[j].nfs=nfs;
      person[j].nms=nms;
      person[j].pbnd=pbnd;
      if (person[j].pbnd>1) if (fid==0) doubled[person[j].pbnd-1][0]=j;
      else doubled[person[j].pbnd-1][1]=j;
#else
      sscanf(line,"%d %d %d %d %d %[^\n]",&pid,&id,&fid,&mid,&sex,rest);
#endif
      person[j].pid=selected;
      person[j].id=id;
      person[j].fid=fid;
      person[j].mid=mid;
      person[j].sex=sex;
      strcpy(line,rest);
      for (k=0;k<nloci;k++,strcpy(line,rest),*rest='\0') {
        if (locus[k].type==AFFECTION) {
          liab=0;
          if (locus[k].vary.aff.nliability==1) sscanf(line,"%d %[^\n]",&aff,rest);
          else sscanf(line,"%d %d %[^\n]",&aff,&liab,rest);
          person[j].ptype[0].aff.aff=aff;
          person[j].ptype[0].aff.liab=liab;      
        } else sscanf(line,"%d %[^\n]",&typed[j][k],rest);
      }
    }
#ifdef USEMAKEPED
    if (doubled[1][0] != 0 && doubled[1][1] == 0) {
      for (k=1; k<=fsize[j]; k++) if (person[k].pbnd==1) {
        person[k].pbnd=2;
        doubled[1][1]=k;
      }
    } else if (doubled[1][1] != 0 && doubled[1][0] == 0) {
      for (k=1; k<=fsize[j]; k++) if (person[k].pbnd==1) {
        person[k].pbnd=2;
        doubled[1][0]=k;
      }
    }
#endif
    simped();
    for (j=0;j<pedsize;++j) {
      for (k=0;k<nloci;++k) switch(locus[k].type) {
        case AFFECTION:
          if(locus[k].vary.aff.nliability==1) skiplines=0;
          else skiplines=person[id].ptype[k].aff.liab-1;
          /*any other order ?*/
          f2=locus[k].vary.aff.pen[skiplines][0];
          f1=locus[k].vary.aff.pen[skiplines][1];
          f0=locus[k].vary.aff.pen[skiplines][2];
          f=ranf();
          m1=person[id].gtype[k][0];m2=person[id].gtype[k][1];
          m=m1+m2;
          switch(m) {
          case 2:aff=1;if(f<f2) aff=2;break;
          case 3:aff=1;if(f<f1) aff=2;break;
          case 4:aff=1;if(f<f0) aff=2;break;
          default: break;
          }
          person[id].ptype[k].aff.aff=aff;break;
        default:
          break;
      }
    }
    outped(fo);
  }
  fclose(fi);
  if(l%100==0) printf("%d done.\n",l);
}
fclose(fo);

free(person);
free(locus);
free_imatrix(typed,0,MAXIND,0,MAXLOCI);
free_ivector(locus_order,0,MAXLOCI);
free_ivector(mg,0,NMG);
free_vector(r,0,MAXLOCI-1);
free_matrix(beta,0,NTRAIT,0,NMG+NPG+NCE+NUE);
free_matrix(kce,0,NTRAIT,0,NCE);
free_ivector(fsize,0,MAXFAM);

return 0;
}

/*end of simmain.c*/
