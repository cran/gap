/*******************************************************
 SIM/SIMPED: simulate phenotype/genotype for one family
*******************************************************/
#include "sim.h"

int simped()
{
int i,j,k,k1,k2,l,l1,l2,m,m1,m2,tl,type;
int pid,id,fid,mid,sex,unsimmed,strand,ori;
float f,g,s1,q,t,d,p,s,z,x,cumf;

#ifdef USEMAKEPED
int a1,a2,cc,t1,t2;
float mean,std1,std2,zt;
FILE *fo;

/*founder's genotype*/
for (j=0;j<pedsize;j++) {
  id=j+1;
  fid=person[id].fid;
  mid=person[id].mid;
  if (!fid) {
    person[id].simmed=true;
    for (m=0;m<2;m++) for (k=0;k<nloci;k++) {
      x=ranf();
      for (l=1;l<=locus[k].nall;l++) 
        if (x<locus[k].cumf[l-1]) {
          person[id].gtype[k][m]=l;
          break;
        }
    }
    if (sexlink&&person[id].sex==MALE) for (k=0;k<nloci;k++)
      person[id].gtype[k][0]=person[id].gtype[k][1];
  }
}

/*loop around non-founders*/
unsimmed=MAXIND;
while (unsimmed>0) {
unsimmed=0;
for (j=0;j<pedsize;j++) {
  id=j+1;
  fid=person[id].fid;
  mid=person[id].mid;
  if (fid>0) {
    if (!person[id].simmed) if(!person[fid].simmed||!person[mid].simmed) unsimmed++;
    else {
      for (m=0;m<2;m++) {
        ori=(m==0)?fid:mid;
        strand=(ranf()<0.5)?0:1;
        for (k=0;k<nloci;k++) {
          person[id].gtype[k][m]=person[ori].gtype[k][strand];
          if (k<nloci) if (ranf()<r[k]) if (strand==0) strand=1; else strand=0;
        }
      }
      if (sexlink&&person[id].sex==MALE) for (k=0;k<nloci;k++)
        person[id].gtype[k][0]=person[id].gtype[k][1];
      person[id].simmed=true;
    }
  } 
  else if (person[id].pbnd>1) {
    tl=doubled[person[id].pbnd-1][1];
    if (person[tl].simmed) {
      person[id].simmed=true;
      for (k=0;k<nloci;k++) for(m=0;m<2;m++)
        person[id].gtype[k][m]=person[tl].gtype[k][m];
    } else if (!person[id].simmed) unsimmed++;
  }
  }
}

/*temporary checking*/
fo=fopen("simulate.tmp","w");
if (!fo) perror("fail to open simulate.tmp");

for (j=0;j<pedsize;j++) {
  id=j+1;
  fid=person[id].fid;
  mid=person[id].mid;
  sex=person[id].sex;
  fprintf(fo,"%5d%4d%4d%4d%3d",selected,person[id].id,fid,mid,sex);
  if (locus[0].type==AFFECTION) {
    fprintf(fo,"%6d",person[id].ptype[0].aff.aff);
    if (locus[0].vary.aff.nliability>1) fprintf(fo,"%3d",person[id].ptype[0].aff.liab);
  }
  for (k=0;k<nloci;k++) if (typed[j][k]==1) switch(locus[k].type) {
    case QTL:
      t1=person[id].gtype[k][0];
      t2=person[id].gtype[k][1];
      if(t1<=t2) {
        a1=t1;
        a2=t2;
      } else {
        a1=t2;
        a2=t1;
      }
      cc=0;
      for (l1=1;l1<=locus[k].nall;l1++)
      for (l2=l1;l2<=locus[k].nall;l2++) {
        mean=locus[k].vary.quan.mean[cc];
        std1=sqrt(locus[k].vary.quan.variance);
        std2=sqrt(locus[k].vary.quan.multivar)*std1;
        if (a1==l1&&a2==l2) {
          if (a1==a2) zt=gennor(mean,std1);
          else zt=gennor(mean,std2);
          tl=person[id].pbnd-1;
          t1=doubled[tl][0];
          t2=doubled[tl][1];
          if (person[id].pbnd>1) {
            if(t1<j+1) zt=person[t1-1].ptype[k].quan.q;
            else if(t2<j+1) zt=person[t2-1].ptype[k].quan.q;
          }
          person[id].ptype[k].quan.q=zt;
          fprintf(fo,"%12.6f",zt);
        }
        cc++;
      }
      break;
    case AFFECTION:
      t1=person[id].gtype[k][0];
      t2=person[id].gtype[k][1];
      if(t1<=t2) {
        a1=t1;
        a2=t2;
      } else {
        a1=t2;
        a2=t1;
      }
      cc=0;
      for (l1=1;l1<=locus[k].nall;l1++) {
        for (l2=l1;l2<=locus[k].nall;l2++) {
          if (a1==l1&&a2==l2) {
            if (ranf()<locus[k].vary.aff.pen[k][cc]) l=2; else l=1;
            tl=person[id].pbnd-1;
            t1=doubled[tl][0];
            t2=doubled[tl][1];
            if (person[id].pbnd>1) {
              if (t1<j+1) l=person[t1-1].ptype[k].aff.aff;
              else if (t2<j+1) l=person[t2-1].ptype[k].aff.aff;
            }
            person[id].ptype[k].aff.aff=l;
            fprintf(fo,"%8d",l);
          }
          cc++;
        }
      }
      break;
    case BINARY:
      a1=person[id].gtype[k][0];
      a2=person[id].gtype[k][1];
      for (l1=0;l1<locus[k].vary.bin.nfactors;l1++)
        person[id].ptype[k].bin.factor[k][l1]=locus[k].vary.bin.factor[a1-1][l1];
      for (l1=0;l1<locus[k].vary.bin.nfactors;l1++) {
        if (locus[k].vary.bin.factor[a2-1][l1]==1)
          person[id].ptype[k].bin.factor[k][l1]=locus[k].vary.bin.factor[a2-1][l1];
      }
      fprintf(fo,"   ");
      for (l1=0;l1<locus[k].vary.bin.nfactors;l1++)
        fprintf(fo,"%2d",person[id].ptype[k].bin.factor[k][l1]);
      break;
    case NUMBERED:
      fprintf(fo,"%5d%3d",person[id].gtype[k][0],person[id].gtype[k][1]);
      break;
  } else switch(locus[k].type) {
    case QTL:
      fprintf(fo,"%12.6f",0.0);break;
    case AFFECTION:
      fprintf(fo,"%8d",0);break;
    case BINARY:
      fprintf(fo,"   ");
      for (l1=1;l1<=locus[k].vary.bin.nfactors;l1++) fprintf(fo,"%2d",0);
      break;
    case NUMBERED:
      fprintf(fo,"%5d%3d",0,0);break;
    }
  putc('\n',fo);
  person[id].simmed=false;
}
fclose(fo);

#else /*old code*/

id=0;
for (l=0;l<pedsize;++l){
  id++;
  pid=person[id].pid;
  fid=person[id].fid;
  mid=person[id].mid;
  sex=person[id].sex;
#ifdef DEBUG
  printf("%3d %3d %3d %3d %3d ",pid,id,fid,mid,sex);
#endif
  if (!fid&&!mid){
    for (j=0;j<nloci;++j){
      l1=l2=1;
      f=ranf();g=ranf();
      for (k=0;k<locus[j].nall;++k) if (f>=locus[j].cumf[k]) l1++;
      for (k=0;k<locus[j].nall;++k) if (g>=locus[j].cumf[k]) l2++;
      person[id].gtype[j][0]=l1;
      person[id].gtype[j][1]=l2;
    }
/*incorporate linkage disequilibrium*/
    if (ndiseq>0) {
      k1=0;k2=0;
      for (j=0;j<nloci;++j){
         m=diseqloci[k1].first_locus;
         if (j==m) k2=m+1; else
         if (j==k2) {/*This is the second locus in disequilibrium*/
           l1=1;l2=1;
           m1=person[id].gtype[k1][0]-1;
           f=ranf();
           for (k=0;k<locus[k2].nall;++k) if (f>=diseqloci[k1].cumf[m1][k]) l1++;
           m2=person[id].gtype[k1][1]-1;
           g=ranf();
           for (k=0;k<locus[k2].nall;++k) if (g>=diseqloci[k1].cumf[m2][k]) l2++;
           person[id].gtype[j][0]=l1;
           person[id].gtype[j][1]=l2;
           k1++;
         }
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
    for (k=1;k<nloci;++k){
      if (k==1) l1=(f>0.5)?1:0;
      l1=(ranf()<=r[k-1])?1-l1:l1;
      if (k==1) l2=(g>0.5)?1:0;
      l2=(ranf()<=r[k-1])?1-l2:l2;
      person[id].gtype[k][0]=person[fid].gtype[k][l1];
      person[id].gtype[k][1]=person[mid].gtype[k][l2];
    }
  }
#ifdef DEBUG
  for(k=0;k<nloci;++k) printf("  %3d%3d",person[id].gtype[k][0],person[id].gtype[k][1]);
  printf("\n");
#endif
}
#endif
/*get phenotypes*/
id=0;
for (l=0;l<pedsize;++l){
  id++;
  pid=person[id].pid;
  fid=person[id].fid;
  mid=person[id].mid;
  sex=person[id].sex;
  for (k=0;k<ntrait;++k) {
    s1=0;
/*Major genes*/
    for (i=0;i<nmg;++i){
      tl=mg[i];
      q=locus[tl].freq[0];p=1-q;
      d=locus[tl].vary.maj.dominance;
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
    for (i=0;i<npg;++i){
      if(fid==0&&mid==0) person[id].pg[k][i]=snorm();
      else person[id].pg[k][i]=snorm()/sqrt(2)+(person[fid].pg[k][i]+person[mid].pg[k][i]);
    }
    s1=0;
    for (i=0;i<npg;++i) s1+=person[id].pg[k][i]*beta[k][nmg+i];
    person[id].trait[k]+=s1;
/*common environment*/
    if (l==0) x=snorm();
    for (i=0;i<nce;++i) {
      if(fid==0&&mid==0) person[id].ce[k][i]=snorm();
      else person[id].ce[k][i]=kce[k][i]*(person[fid].ce[k][i]+person[mid].ce[k][i])
                               +sqrt(1-2*pow(kce[k][i],2))*x;
    }
    s1=0;
    for (i=0;i<nce;++i)
      s1+=person[id].ce[k][i]*beta[k][nmg+npg+i];
    person[id].trait[k]+=s1;
/*unique environment*/
    for (i=0;i<nue;++i)
      person[id].trait[k]+=snorm()*beta[k][nmg+npg+nce+i];
  }
}

return 0;
}
/*End of simped.c*/
