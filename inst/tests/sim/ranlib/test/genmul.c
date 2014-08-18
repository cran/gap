#include <stdio.h>
#include "ranlib.h"
main()
{
int replicates=1000,ncat=4;
int i;
long ix[4];
float p[4]={0.25,0.25,0.25,0.25};
genmul(replicates,p,ncat,ix);
for(i=0;i<4;++i) printf("%ld ",ix[i]);
printf("\n");
return;
}
