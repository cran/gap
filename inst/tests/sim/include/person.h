#ifndef FAMILY_H
#define FAMILY_H

typedef struct
{
  int aff,liab;
} affection_p;
typedef struct
{
  int factor[MAXFACTOR][MAXFACTOR];
} binary_p;
typedef struct
{
  float q;
} quantitative_p;

typedef struct
{
  int pid, id, fid, mid, sex, gtype[MAXLOCI][2];
#ifdef USEMAKEPED
  int offs,nfs,nms,pbnd; 
#endif
  int nchildren,nspouse,isproband,simmed;
  float trait[NTRAIT],pg[NTRAIT][NPG],ce[NTRAIT][NCE];
  union {
    affection_p aff;
    binary_p bin;
    quantitative_p quan;
  } ptype[MAXLOCI];
} PERSON;

typedef struct
{
  int pid, father, mother, nsibs, sib[MAXSIBS];
} NUC_FAMILY;

#endif
