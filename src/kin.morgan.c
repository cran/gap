#include "nghds.h"

#define max_size 400
#define STRICT_CHECK

static void nullify(Ind *);
double kinship(Ind *, Ind *),inbreeding(Ind *);

Ind nullnode;

void kin_morgan(int *data, int *pedsize, double *kin)
{
Ind *ped,*t1,*t2;
int i,j,k,id,pa,ma;

nullify(&nullnode);
ped=(Ind *)malloc((max_size+1)*sizeof(Ind));
if (!ped)
{
  printf("\nError to allocate memory for pedigree\n");
  return;
}
for (i=0;i<=max_size;i++) nullify(&ped[i]);

printf("\nThe original family (ID PA MA): \n\n");
i=0;
for(i=1;i<=*pedsize;i++)
{
      id=data[(i-1)*3];
      pa=data[(i-1)*3+1];
      ma=data[(i-1)*3+2];
      ped[i].self=id;
      ped[i].index=i;
      printf("%5d %5d %5d\n",id,pa,ma);
#ifdef STRICT_CHECK
      t1=&ped[pa];
      t2=&ped[ma];
      if ((pa && t1->self==UNKNOWN)||(ma && t2->self==UNKNOWN))
      {
         printf("\nParents not in datafile, quit\n");
         return;
      }
#endif
}
for (i=1;i<=*pedsize;i++)
{
    t1=t2=&nullnode;
    pa=data[(i-1)*3+1];
    ma=data[(i-1)*3+2];
    if (pa) t1=&ped[pa];
    if (ma) t2=&ped[ma];
    ped[i].pa=t1;
    ped[i].ma=t2;
}
k=0;
for (i=1;i<=*pedsize;i++)
{
    printf("%5d ",i);
    for (j=1;j<=i;j++)
    {
        kin[k]=kinship(&(ped[i]),&(ped[j]));
        printf(" %f",kin[k]);
        k++;
    }
    printf("\n");
}

for (i=0;i<=*pedsize;i++) nullify(&ped[i]);
free(ped);

return;

}

static void nullify(Ind *nul)
/*
 * Function to make the nullnode for the neighbourhood routine.
 * Pointers to pa, ma and marriages are set to NULL. It is assumed that a
 * structure  of type Ind has already been set up as nullnode and defined
 * to be a global variable.
 */
{
        nul->self = UNKNOWN;
        nul->index = INVALID_INDEX;
        nul->pa = NULL;
        nul->ma = NULL;
        nul->marriages = NULL;
}

double kinship(Ind * a, Ind * b)
/*
 * Recursive program which returns the kinship coefficient between
 * individuals a and b.
 */
{
    if (a == &nullnode || b == &nullnode)
         return 0.0;

    if (a == b)
        return 0.5 + 0.5 * inbreeding(a);
    else
    if (a->pa->self == 0)
    {
        if (b->index < a->index)
            return 0.0;
        else
        if (b->pa->self == 0)
            return 0.0;
        else
            return (kinship(a, b->pa) + kinship(a, b->ma)) * 0.5;
    } else
    if (b->pa->self == 0)
    {
        if (a->index < b->index)
            return 0.0;
        else
            return (kinship(b, a->pa) + kinship(b, a->ma)) * 0.5;
    } else
    if (a->index < b->index)
        return (kinship(a, b->pa) + kinship(a, b->ma)) * 0.5;
    else
        return (kinship(b, a->pa) + kinship(b, a->ma)) * 0.5;
}

double inbreeding(Ind * a)
/*
 * A recursive program that returns the inbreeding coefficient
 * for individual a.
 */
{
    if (a == &nullnode)
         return 0.0;
    else
        return kinship(a->pa, a->ma);
}
