/*******************************************************
 L O C P R E P
*******************************************************/
#include "sim.h"

/*internals*/

int write_locus();

int locprep(FILE *fi)
{
char line[500],rest[500];
int i,j,k,l;
float s,t;

fgets(line,500,fi);
sscanf(line,"%d %*d %d %*d %[^\n]",&nloci,&sexlink,rest);
sscanf(rest,"%d %d %d %d %d %d",&ntrait,&nmg,&npg,&nce,&nue,&ndiseq);

fgets(line,500,fi);
fgets(line,500,fi);
for (i=0;i<nloci;++i,strcpy(line,rest),*rest='\0') {
  if (sscanf(line,"%d %[^\n]",&locus_order[i],rest)<1)
    printf("check locus order -)\n");
}

l=0;
for (i=0;i<nloci;++i) {
  if (!fgets(line,500,fi)) printf("error in reading locus\n");
  sscanf(line,"%d %d",&locus[i].type,&locus[i].nall);
  fgets(line,500,fi);
  s=0;
  for (j=0;j<locus[i].nall;j++,strcpy(line,rest),*rest='\0')
  {
    if (sscanf(line,"%f %[^\n]",&locus[i].freq[j],rest)<1) printf("check gene frequency -)\n");
    s+=locus[i].freq[j];
    locus[i].cumf[j]=s;
  }
  if (fabs(s-1.0)>0.01) printf("Allele frequencies at locus %d do not sum to 1",i+1);
  switch(locus[i].type){
  case QTL:
    fgets(line,500,fi);
    for(j=0;j<locus[i].nall*(locus[i].nall+1)/2;j++,strcpy(line,rest),*rest='\0')
      if(sscanf(line,"%f %[^\n]",&locus[i].vary.quan.mean[j],rest)<1)
       printf("check means of quatitative trait");
    fgets(line,500,fi);
    sscanf(line,"%f",&locus[i].vary.quan.variance);
    fgets(line,500,fi);
    sscanf(line,"%f",&locus[i].vary.quan.multivar);
    break;
  case AFFECTION:
    fgets(line,500,fi);
    sscanf(line,"%d",&locus[i].vary.aff.nliability);/*add locus[i]*/
    for (j=0;j<locus[i].vary.aff.nliability;++j) {
       fgets(line,500,fi);
       for (k=0;k<3;++k,strcpy(line,rest),*rest='\0')
          if (sscanf(line,"%f %[^\n]",&locus[i].vary.aff.pen[j][k],rest)<1)
            printf("check your penetrances -)\n");
    }
    break;
  case BINARY:
    fgets(line,500,fi);
    sscanf(line,"%d",&locus[i].vary.bin.nfactors);
    for (j=0;j<locus[i].nall;j++)
    {
      fgets(line,500,fi);
      for (k=0;k<locus[i].vary.bin.nfactors;k++,strcpy(line,rest),*rest='\0')
         if (sscanf(line,"%d %[^\n]",&locus[i].vary.bin.factor[j][k],rest)<1)
           printf("check binary factors -)\n");
    }
    break;
  case NUMBERED:
    break;
  case MAJOR:
    fgets(line,500,fi);
    sscanf(line,"%f ",&locus[i].vary.maj.dominance);
    mg[l++]=i;
    break;
  default:
    break;
  }
}
fgets(line,500,fi);
fgets(line,500,fi);
for (i=0;i<nloci-1;++i,strcpy(line,rest),*rest='\0')
  if (sscanf(line,"%f %[^\n]",&r[i],rest)<1)
    printf("check recombination fracs -)\n");
fgets(line,500,fi);
for (i=0;i<ntrait;++i) {
  s=0;
  fgets(line,500,fi);
  for (j=0;j<nmg+npg+nce+nue;++j,strcpy(line,rest),*rest='\0')
  {
   if (sscanf(line,"%f %[^\n]",&beta[i][j],rest)<1)
     printf("check betas -)\n");
   s+=pow(beta[i][j],2);
  }
  for (j=0;j<nmg+npg+nce+nue;j++) beta[i][j]/=s;
}
for (i=0;i<ntrait;++i) {
  fgets(line,500,fi);
  for (j=0;j<nce;++j,strcpy(line,rest),*rest='\0')
    if (sscanf(line,"%f %[^\n]",&kce[i][j],rest)<1)
      printf("check environment transmisson for each trait -)\n");
}
for (k=0;k<ndiseq;k++)
{
  fgets(line,500,fi);
  if (sscanf(line,"%d %[^\n]",&l,rest)<1) printf("check disequilibrium -)\n");
  diseqloci[k].first_locus=l;
  for(i=0;i<locus[l].nall;i++)
  {
    fgets(line,500,fi);
    s=0;
    for (j=0;j<locus[l+1].nall;j++,strcpy(line,rest),*rest='\0')
    {
      if (sscanf(line,"%f %[^\n]",&t,rest)<1) printf("check disequilibrium -)\n");
      diseqloci[k].freq[i][j]=t;
      s+=diseqloci[k].freq[i][j];
      diseqloci[k].cumf[i][j]=s;
    }
    for (j=0;j<locus[l+1].nall;j++) diseqloci[k].cumf[i][j]/=s;
  }
}
#ifdef DEBUG
write_locus();
#endif

return 0;
}

int write_locus()
{
int i,j,k,l;

printf("%d %d %d %d %d %d %d %d\n",nloci,ntrait,nmg,npg,nce,nue,ndiseq,sexlink);
for (i=0;i<nloci;++i) printf("%d ",locus_order[i]);printf("\n");
for (i=0;i<nloci;++i) {
  printf("%d  %d\n",locus[i].type,locus[i].nall);
  for (j=0;j<locus[i].nall;++j) printf("%f ",locus[i].freq[j]);printf("\n");
  switch(locus[i].type){
  case QTL:
    for (j=0;j<locus[i].nall*(locus[i].nall+1)/2;j++)
      printf("%.2f ",locus[i].vary.quan.mean[j]);
    printf("\n");
    printf("%.2f\n",locus[i].vary.quan.variance);
    printf("%.2f\n",locus[i].vary.quan.multivar);
    break;
  case AFFECTION:
    printf("%d\n",locus[i].vary.aff.nliability);
    for (j=0;j<locus[i].vary.aff.nliability;++j) {
       for (k=0;k<3;++k) printf("%f ",locus[i].vary.aff.pen[j][k]); printf("\n");
    }
    break;
  case BINARY:
    printf("%d\n",locus[i].vary.bin.nfactors);
    for(j=0;j<locus[i].nall;j++)
    {
      for(k=0;k<locus[i].vary.bin.nfactors;k++)
        printf("%2d ",locus[i].vary.bin.factor[j][k]);
      printf("\n");
    }
    break;
  case NUMBERED:
    break;
  case MAJOR:
    printf("%f\n",locus[i].vary.maj.dominance);
    break;
  default:
    break;
  }
}
for (i=0;i<nloci-1;++i) printf("%f ",r[i]); printf("\n");
for (i=0;i<ntrait;++i) {
  for (j=0;j<nmg+npg+nce+nue;++j) printf("%f ",beta[i][j]);printf("\n");
}
for (i=0;i<ntrait;++i) {
  for(j=0;j<nce;++j) printf("%f ",kce[i][j]); printf("\n");
}
for (k=0;k<ndiseq;k++)
{
  l=diseqloci[k].first_locus;
  printf("%d\n",l);
  for (i=0;i<locus[l].nall;i++)
  {
    for (j=0;j<locus[l+1].nall;j++)
      printf("%f ",diseqloci[k].cumf[i][j]);printf("\n");
  }
}

return 0;
}

/*End of locprep.c*/
