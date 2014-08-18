/******************************************************
 SIM/GETSIZE -- obtain family sizes
*******************************************************/
#include "sim.h"

#define nl 20
int getsize(char pedinn[nl])
/*decide family size*/
{
FILE *pedin;
int i=0,j=0,k=0,eof=0;
char ll[101],name[nl],id[nl],l1[nl],l2[nl];
typedef struct {char name[nl],id[MAXIND][nl];int n;} ids;
ids *ff;

ff=(ids*)malloc(MAXFAM*sizeof(ids));
if (!ff) perror("error allocating memory in getsize()");
pedin=fopen(pedinn,"r");
while (1) {
  strcpy(ll,"\0");
  if (fgets(ll,100,pedin)==NULL) eof=1;
  if (sscanf(ll,"%s%s",&name,&id)<2) eof=1;
  strcpy(l2,name);
  if (eof) strcpy(l2,"\0");
  if (j==0) {
    strcpy(l1,l2);
    goto next;
  }
  if (k==0) strcpy(ff[i].name,name);
  if (!strcmp(l1,l2)) k++;
  else {
    strcpy(l1,l2);
    ff[i].n=k+1;
    ++i;
    k=0;
  }
next:
  strcpy(ff[i].id[k],id);
  j++;
  if (eof) break;
}
fclose(pedin);
#ifdef DEBUG
printf("\nThere are %d families with sizes:\n",i);
for (j=0;j<i;j++) {
  fsize[j]=ff[j].n;
  printf("%s %d",ff[j].name,ff[j].n);
  for (k=0; k<ff[j].n; k++) printf(" %s",ff[j].id[k]);printf("\n");
}
#endif
free(ff);
return i;
}
#undef nl

#ifdef _useoldcode /*obsolete*/
id=1;k=1;i=0;
l1=person[id].pid;
for(j=0;j<MAXALL;++j){
   id++;
   l2=person[id].pid;
   if(l1!=l2&&l1!=0)
   {
     fsize[++i]=k;
     l1=l2;k=1;
   }
   else ++k;
}
nfam=i;
#endif

/*end of getsize.c*/
