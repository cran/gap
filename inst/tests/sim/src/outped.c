/*******************************************************
 SIM/OUTPED: output simulated family information
*******************************************************/
#include "sim.h"

int outped(FILE *fo)
{
int i,j,k,l,m,l1,l2,id,aff;
int pid,fid,mid,offs,nfs,nms,sex,pbnd;
float f;

id=0;
for (l=0;l<pedsize;++l) {
  id++;
  pid=selected;
  fid=person[id].fid;
  mid=person[id].mid;
  sex=person[id].sex;
#ifdef USEMAKEPED
  offs=person[id].offs;
  nfs=person[id].nfs;
  nms=person[id].nms;
  pbnd=person[id].pbnd;
  fprintf(fo,"%4d %4d %4d %4d %4d %4d %4d%2d%2d",pid,id,fid,mid,offs,nfs,nms,sex,pbnd);
#else
  fprintf(fo,"%4d %4d %4d %4d %4d",pid,id,fid,mid,sex);
#endif
  for(k=0;k<ntrait;k++) fprintf(fo," %12.5f",person[id].trait[k]);
  k=0;
  for (j=0;j<nloci;++j) {
    l1=person[id].gtype[j][0];l2=person[id].gtype[j][1];
    m=l1+l2;
    switch(locus[j].type) {
    case QTL:
      fprintf(fo," %12.5f",person[id].ptype[j].quan.q);
      break;
    case AFFECTION:
      aff=person[id].ptype[j].aff.aff;
#ifdef DEBUG
      printf("family %d %d (%d %d %d) is aff=%d\n",pid,id,l1,l2,m,aff);
#endif
      fprintf(fo,"%3d",aff);
      if(locus[j].vary.aff.nliability>1) fprintf(fo,"%3d",person[id].ptype[j].aff.liab);
      break;
    case BINARY:
      break;
    case NUMBERED:
      fprintf(fo," %3d%3d",l1,l2);
      break;
    case MAJOR:
      ++k;break;
    default: 
      break;
    }
  } fprintf(fo,"\n");
}
return 0;
}
/*End of outped.c*/
