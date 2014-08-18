#ifndef LOCUS_H
#define LOCUS_H

enum {UNKNOWN,MALE,FEMALE};
enum {QTL,AFFECTION,BINARY,NUMBERED,MAJOR};

typedef struct
{
  int nliability;
  float pen[MAXLIABILITY][3];
} affection;
typedef struct
{
  int nfactors;
  int factor[MAXFACTOR][MAXFACTOR];
} binary;
typedef struct
{
  float mean[MAXALLELES*(MAXALLELES+1)/2],variance,multivar;
} quantitative;
typedef struct
{
  float dominance;
} major;
typedef struct
{
  char name[10];
  int type,nall;
  float freq[MAXALLELES],cumf[MAXALLELES];
  union {
    affection aff;
    binary bin;
    quantitative quan;
    major maj;
  } vary;
} LOCUS;

typedef struct 
{
int first_locus;
float freq[MAXALLELES][MAXALLELES];
float cumf[MAXALLELES][MAXALLELES];
} DISEQ_LOCI;

#endif
