/***J.H. Zhao & Sham P.C. 03MAR97 I.O.P.           ****
***see companion documentation for more information****
*******************************************************
  S I M . H
******************************************************/
#ifndef SIM_H
#define SIM_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ranlib.h>
#include <string.h>
#include <nrutil.h>

#define DEBUG
#define USEMAKEPED
#undef USEMAKEPED

#define MAXFAM 10
#define MAXSIBS 20
#define MAXIND  50
#define MAXALL MAXFAM*MAXIND
#define SEXRATIO 1.1
#define SELECTED_FAMILY 5

#define MAXLOCI 20
#define MAXLIABILITY 3
#define MAXFACTOR 3
#define MAXALLELES 10
#define MAXDISEQ 5
#define NTRAIT 5
#define NMG 5
#define NPG 5
#define NCE 5
#define NUE 2
#define NDISEQ 2

#include "func.h"
#include "locus.h"
#include "person.h"

int locprep(FILE*);
int simped();
int outped(FILE*);

int selected;
int nfam,nloci,ntrait,nmg,npg,nce,nue,ndiseq,sexlink;
int **typed;

LOCUS *locus;
PERSON *person,nobody;
DISEQ_LOCI diseqloci[MAXDISEQ];
int *locus_order,*mg,*fsize,pedsize,doubled[5][2];
float *r,**beta,**kce;

#endif
