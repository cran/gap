/*SIMULATE (C version) JH Zhao 5/4/1999 IoP*/
/*way to handle maxpeds needs to be revised*/

#include <stdio.h>
#include <string.h>
#include <math.h>

#define version    "2.21" /*22 Apr 1996*/
#define maxpeds    100    /*Max. no. of pedigrees per replicate*/
#define maxloci    32     /*Max. number of loci*/
#define maxind     20     /*Max. no. individuals in one pedigree*/
#define maxall     12     /*Maximum number of alleles at any one locus*/
#define maxgen     (maxall*(maxall+1)/2)
#define maxfact    maxall

typedef enum boolean{false,true} boolean;

FILE *pedfile;
int numgen, numreps, numpeds, numloci;
int liabnum[maxloci], locustype[maxloci], numall[maxloci];
int numind[maxpeds], indnum[maxpeds][maxind], pa[maxpeds][maxind],
    ma[maxpeds][maxind], fo[maxpeds][maxind], nps[maxpeds][maxind],
    nms[maxpeds][maxind], pro[maxpeds][maxind],
    disease[maxpeds][maxind], liabcl[maxpeds][maxind];
int sex[maxpeds][maxind], typed[maxpeds][maxind][maxloci];
int binall[maxloci][maxall][maxfact], facttot[maxloci][maxfact];
int nfact[maxloci], ntrait[maxloci], aphen[maxpeds][maxind][maxloci];
int ss, seed1, seed2, seed3, seedorig1, seedorig2, seedorig3;
int seqrun=-1, numrepdone, howoften=100;
int low, locus[maxpeds][maxind][maxloci][2];
int gasdeviset, sexlink, sexdiff, doubled[maxpeds][5][2];
char pedfilen[30];
double sexratio, gasdevgset, theta[maxloci+1][2], genefreq[maxloci][maxall];
double multvar[maxloci], vari[maxloci];
double mean[maxloci][maxgen], pen[maxloci][maxgen],qphen[maxpeds][maxind][maxloci];
boolean simmed[maxpeds][maxind+1], bad=false;

char peek(FILE *fp)
{
char s;
s=getc(fp);ungetc(s,fp);return(s);
}

boolean eoln(FILE *fp)
{
char s;
s=getc(fp);ungetc(s,fp);
if(feof(fp)) return(feof(fp));return((boolean)(s=='\n'));
}

boolean seekeoln(FILE **f)
{
while((peek(*f)==' '||peek(*f)=='\t')&&!eoln(*f)) if(!feof(*f)) getc(*f);
return (eoln(*f)||feof(*f));
}

double ranuni()
{
  double r;

  seed1 = seed1 % 177 * 171 - seed1 / 177 * 2;
  seed2 = seed2 % 176 * 172 - seed2 / 176 * 35;
  seed3 = seed3 % 178 * 170 - seed3 / 178 * 63;
  if (seed1 < 0) seed1 += 30269;
  if (seed2 < 0) seed2 += 30307;
  if (seed3 < 0) seed3 += 30323;
  r = seed1 / 30269.0 + seed2 / 30307.0 + seed3 / 30323.0;
  return (r - (long)r);
}

double gasdev()
{
  double result, f, r, v1, v2;

  if (gasdeviset != 0) {
    gasdeviset = 0;
    return gasdevgset;
  }
  do {
    v1 = 2.0 * ranuni() - 1.0;
    v2 = 2.0 * ranuni() - 1.0;
    r = v1 * v1 + v2 * v2;
  } while (r >= 1.0 || r <= 0.0);
  f = sqrt(-2.0 * log(r) / r);
  gasdevgset = v1 * f;
  result = v2 * f;
  gasdeviset = 1;
  return result;
}

double rannor(double mu, double sig)
{
  double tmp;

  tmp = sig * gasdev() + mu;
  if (fabs(tmp) < 0.00001) tmp = 0.000001;
  return tmp;
}

void initialize()
{
  FILE *problem, *simdata, *simped;
  char problemn[30], simdatan[30], simpedn[30];
  int i1, i2, i3, i, j, k, l, loopend, pedn, high;
  double totfreq;

  printf("   **********************************************\n");
  printf("   *                                            *\n");
  printf("   *       Program  SIMULATE  version %s      *\n", version);
  printf("   *                                            *\n");
  printf("   *          Written by J. Terwilliger         *\n");
  printf("   *                                            *\n");
  printf("   **********************************************\n\n");
  printf("Input files required:\n");
  printf("  PROBLEM.DAT  3 random seeds and number of replicates\n");
  printf("  SIMDATA.DAT  datafile in LINKAGE format\n");
  printf("  SIMPED.DAT   pedfile in LINKAGE format\n\n");
  printf("Maximum values of parameters:\n");
  printf("%8d pedigrees\n", maxpeds);
  printf("%8d loci\n", maxloci);
  printf("%8d individuals in any one pedigree\n", maxind);
  printf("%8d alleles at any one locus\n", maxall);
  printf("%8d binary factors\n\n", maxfact);
  for (i1 = 0; i1 < maxpeds; i1++) {
    numind[i1]=0;
    for (i2 = 0; i2 < 5; i2++)
      doubled[i1][i2][0] = doubled[i1][i2][1]=0;
    for (i2 = 0; i2 < maxind; i2++) {
      indnum[i1][i2]=pa[i1][i2]=ma[i1][i2]=fo[i1][i2]=nps[i1][i2]=0;
      nms[i1][i2]=pro[i1][i2]=disease[i1][i2]=liabcl[i1][i2]=0;
    }
    for (i2 = 0; i2 <= maxind; i2++) simmed[i1][i2] = false;
    for (i2 = 0; i2 < maxind; i2++) {
      sex[i1][i2] = 1;
      for (i3 = 0; i3 < maxloci; i3++) {
        locus[i1][i2][i3][0]=locus[i1][i2][i3][1]=1;
        typed[i1][i2][i3]=qphen[i1][i2][i3]=aphen[i1][i2][i3]=0;
      }
    }
  }
  for (i1 = 0; i1 <= maxloci; i1++) theta[i1][0]=theta[i1][1]=0.0;
  for (i1 = 0; i1 < maxloci; i1++) {
    liabnum[i1]=locustype[i1]=numall[i1] = 0;
    nfact[i1]=ntrait[i1]=multvar[i1]=vari[i1] = 0.0;
    for (i2 = 0; i2 < maxgen; i2++) mean[i1][i2]=pen[i1][i2] = 0.0;
    for (i2 = 0; i2 < maxfact; i2++) facttot[i1][i2] = 0;
    for (i2 = 0; i2 < maxall; i2++) {
      genefreq[i1][i2] = 0.0;
      for (i3 = 0; i3 < maxfact; i3++) binall[i1][i2][i3] = 0;
    }
  }

/*read PROBLEM.DAT */

  printf("Opening \"problem.dat\" file for input\n");
  strcpy(problemn, "problem.dat");
  problem=fopen(problemn,"r");
  if (problem == NULL) perror("File not found");
  fscanf(problem, "%d%d%d", &seed1, &seed2, &seed3);
  if (seed1<1||seed1>30323||seed2<1||seed2>30323||seed3<1||seed3>30323) {
    printf(" ERROR:  Random seed(s) out of range 0--30324\n");
    printf(" Seeds are  %d %d %d\n", seed1, seed2, seed3);
  }
  seedorig1 = seed1;
  seedorig2 = seed2;
  seedorig3 = seed3;
  fscanf(problem, "%d", &numreps);
  if (!seekeoln(&problem)) fscanf(problem, "%d", &seqrun);
  fclose(problem);

/*read SIMDATA.DAT */

  printf("Opening \"simdata.dat\" file for input\n");
  strcpy(simdatan, "simdata.dat");
  simdata=fopen(simdatan,"r");
  if (simdata==NULL) perror("File not found");
  fscanf(simdata, "%d%d%d%*[^\n]", &numloci, &i, &sexlink);
  getc(simdata);
  fscanf(simdata, "%*[^\n]");
  getc(simdata);
  fscanf(simdata, "%*[^\n]");
  getc(simdata);
  for (i = 0; i < numloci; i++) {
    fscanf(simdata, "%d%d%*[^\n]", &locustype[i], &numall[i]);
    getc(simdata);
    for (j = 0; j < numall[i]; j++) fscanf(simdata, "%lg", &genefreq[i][j]);
    fscanf(simdata, "%*[^\n]");
    getc(simdata);
    totfreq = 0.0;
    for (j = 0; j < numall[i]; j++) totfreq += genefreq[i][j];
    if (fabs(totfreq - 1.0) > 0.01)
      printf("ERROR:  Allele frequencies at locus%2d do not sum to 1.\n",i + 1);
    switch (locustype[i]) {
    case 0:
      fscanf(simdata, "%d%*[^\n]", &ntrait[i]);
      getc(simdata);
      if (ntrait[i] > 1)
        printf("ONLY ONE QUANTITATIVE TRAIT AT A LOCUS ALLOWED!\n");
      numgen = (long)floor(numall[i] * (numall[i] + 1.0) / 2 + 0.5);
      for (j = 0; j < numgen; j++) fscanf(simdata, "%lg", &mean[i][j]);
      fscanf(simdata, "%*[^\n]");
      getc(simdata);
      fscanf(simdata, "%lg%*[^\n]", &vari[i]);
      getc(simdata);
      fscanf(simdata, "%lg%*[^\n]", &multvar[i]);
      getc(simdata);
      break;
    case 1:
      fscanf(simdata, "%d%*[^\n]", &liabnum[i]);
      getc(simdata);
      if (sexlink) loopend = liabnum[i] * 2;
      else loopend = liabnum[i];
      for (j = 1; j <= loopend; j++) {
        if (i > 0) {
          numgen = numall[i] * (numall[i] + 1) / 2;
          if (sexlink) numgen += numall[i];
          for (l=0;l<numgen;l++) fscanf(simdata, "%lg", &pen[i][l]);
          fscanf(simdata, "%*[^\n]");
          getc(simdata);
        } else {
          fscanf(simdata, "%*[^\n]");
          getc(simdata);
        }
      }
      break;
    case 2:
      fscanf(simdata, "%d%*[^\n]", &nfact[i]);
      getc(simdata);
      for (j = 0; j < numall[i]; j++)
        for (k = 0; k < nfact[i]; k++) fscanf(simdata, "%d", &binall[i][j]);
      fscanf(simdata, "%*[^\n]");
      getc(simdata);
      break;
    }
  }
  fscanf(simdata, "%d%*[^\n]", &sexdiff);
  getc(simdata);
  for (i = 1; i < numloci; i++) {
    fscanf(simdata, "%lg", &theta[i][0]);
    theta[i][1] = theta[i][0];
  }
  fscanf(simdata, "%*[^\n]");
  getc(simdata);
  if (sexdiff == 1) {
    fscanf(simdata, "%lg%*[^\n]", &sexratio);
    getc(simdata);
    for (i = 1; i < numloci; i++)
      theta[i][1] =0.5*(1-exp(sexratio*log(1-2*theta[i][0])));
  } else if (sexdiff == 2)
    for (i = 1; i < numloci; i++) fscanf(simdata, "%lg", &theta[i][1]);
  fscanf(simdata, "%*[^\n]");
  getc(simdata);
  fclose(simdata);

/*read SIMPED.DAT */

  printf("Opening \"simped.dat\" file for input\n");
  strcpy(simpedn, "simped.dat");
  if (simped != NULL) simped = freopen(simpedn, "r", simped);
  else simped = fopen(simpedn, "r");
  if (simped == NULL) perror("File not found");
  fscanf(simped, "%d", &numpeds);
  if (numpeds > maxpeds) {
    printf(" ERROR:  number of pedigrees specified too large (>%ld)\n",(long)maxpeds);
    bad = true;
    goto _L90;
  }
  for (i = 0; i < numpeds; i++) fscanf(simped, "%d", &numind[i]);
  fscanf(simped, "%*[^\n]");
  getc(simped);
  if (numind[0] == 1) {
    printf("WARNING:  Pedigree 1 contains only 1 individual.  Perhaps header line\n");
    printf("          in SIMPED.DAT missing.  Press <ENTER> to continue or Ctrl-C\n");
    printf("          to stop.\n");
    scanf("%*[^\n]");
    getchar();
  }
  low = 1;
  for (j = 0; j < maxpeds; j++)
    for (i = 0; i < 5; i++) doubled[j][i][0]=doubled[j][i][1]=0;
  for (j = 0; j < numpeds; j++) {
    for (i = 0; i < numind[j]; i++) {
      if (numind[j] > maxind) {
        printf(" ERROR: number of individuals larger than max of %ld\n",(long)maxind);
        printf("%d-th pedigree\n", j + 1);
        bad = true;
        goto _L90;
      }
      fscanf(simped, "%d%d%d%d%d%d%d%d%d", &pedn, &indnum[j][i],
             &pa[j][i], &ma[j][i], &fo[j][i], &nps[j][i], &nms[j][i],
             &sex[j][i], &pro[j][i]);
      if (indnum[j][i] != i + 1) {
        printf(" ERROR:  INDIVIDUALS MUST BE NUMBERED SEQUENTIALLY IN EACH PEDIGREE,\n");
        printf(" STARTING FROM 1!  Pedigree %d  %d  (top line incorrect?)\n",j + 1, pedn);
        bad = true;
        goto _L90;
      }
      if (pro[j][i] > 1) {
        if (pa[j][i] == 0) doubled[j][pro[j][i] - 1][0] = i + 1;
        else doubled[j][pro[j][i]-1][1] = i + 1;
      }
      if (locustype[0] == 1) {
        fscanf(simped, "%d", &disease[j][i]);
        if (liabnum[0] > 1) fscanf(simped, "%d", &liabcl[j][i]);
        low = 2;
      }
      for (k=low-1;k<numloci;k++) fscanf(simped,"%d",&typed[j][i][k]);
      fscanf(simped, "%*[^\n]");
      getc(simped);
    }
  }
  for (i = 0; i < numpeds; i++) {
    if (doubled[i][1][0] != 0 && doubled[i][1][1] == 0) {
      for (k = 1; k <= numind[i]; k++)
        if (pro[i][k - 1] == 1) {
          pro[i][k - 1] = 2;
          doubled[i][1][1] = k;
        }
    } else if (doubled[i][1][1] != 0 && doubled[i][1][0] == 0) {
      for (k = 1; k <= numind[i]; k++)
        if (pro[i][k - 1] == 1) {
          pro[i][k - 1] = 2;
          doubled[i][1][0] = k;
        }
    }
  }
_L90:
  fclose(simped);
}

void founders()
{
  int i, j, k, l, m;
  double xx, cumf;

  if (locustype[0] == 1) low = 2; else low = 1;
  for (i = 0; i < numpeds; i++) {
    for (j = 0; j < numind[i]; j++) {
      if (pa[i][j] == 0 && pro[i][j] <= 1) {
        simmed[i][j + 1] = true;
        for (m = 0; m <= 1; m++) {
          for (k = low - 1; k < numloci; k++) {
            xx = ranuni();
            cumf = 0.0;
            for (l = 1; l <= numall[k]; l++) {
              if (l == numall[k]) cumf = 2.0; else cumf+=genefreq[k][l-1];
              if (xx < cumf) {
                locus[i][j][k][m] = l;
                goto _L10;
              }
            }
_L10: ;
          }
        }
        if (sexlink==1&&sex[i][j]==1)
          for (k=low-1;k<numloci;k++) locus[i][j][k][0]=locus[i][j][k][1];
      }
    }
  }
}

void finishsimulation()
{
  int i, j, k, m, unsimmed, chro, ori;

  unsimmed = 999;
  for (i = 0; i < numpeds; i++) simmed[i][0] = true;
  while (unsimmed > 0) {
    unsimmed = 0;
    for (i = 0; i < numpeds; i++) for (j = 0; j < numind[i]; j++) {
      if (pa[i][j] > 0) {
        if (!simmed[i][j + 1]) {
          if (!simmed[i][pa[i][j]] || !simmed[i][ma[i][j]]) unsimmed++;
          else {
            for (m = 0; m <= 1; m++) {
              switch (m + 1) {
              case 1: ori = pa[i][j];break;
              case 2: ori = ma[i][j];break;
              }
              chro=(ranuni()<0.5)?1:2;
              for (k=low;k<=numloci;k++) {
                locus[i][j][k-1][m]=locus[i][ori-1][k-1][chro-1];
                if (k<numloci && ranuni()<theta[k][m]) chro=(chro==1)?2:1;
              }
            }
            if (sexlink && sex[i][j] == 1) {
              for (k = low - 1; k < numloci; k++)
                locus[i][j][k][0] = locus[i][j][k][1];
            }
            simmed[i][j + 1] = true;
          }
        }
      } else if (pro[i][j] > 1) {
        if (simmed[i][doubled[i][pro[i][j]-1][1]]) {
          simmed[i][j+1]=true;
          for (k=low-1;k<numloci;k++) {
            for (m=0;m<=1;m++)
              locus[i][j][k][m]=locus[i][doubled[i][pro[i][j]-1][1]-1][k][m];
          }
        } else if (!simmed[i][j + 1]) unsimmed++;
      }
    }
  }
}

void replicate()
{
  int cp, a1, a2, number, i, j, k, l, m;

  founders();
  finishsimulation();
  for (i = 0; i < numpeds; i++) {
    for (j = 0; j < numind[i]; j++) {
      number = (ss-1)*numpeds+i+1;
      fprintf(pedfile, "%5d%4d%4d%4d%4d%4d%4d%3d%3d",
              number, indnum[i][j], pa[i][j], ma[i][j], fo[i][j], nps[i][j],
              nms[i][j], sex[i][j], pro[i][j]);
      if (locustype[0] == 1) {
        fprintf(pedfile, "%6d", disease[i][j]);
        if (liabnum[0] > 1) fprintf(pedfile, "%3d", liabcl[i][j]);
      }
      for (k = low - 1; k < numloci; k++) {
        if (typed[i][j][k] == 1) {
          switch (locustype[k]) {
          case 0:
            if (locus[i][j][k][0] <= locus[i][j][k][1]) {
              a1 = locus[i][j][k][0];
              a2 = locus[i][j][k][1];
            } else {
              a1 = locus[i][j][k][1];
              a2 = locus[i][j][k][0];
            }
            cp = 0;
            for (l=1;l<=numall[k];l++) {
              for (m = l; m <= numall[k]; m++) {
                cp++;
                if (a1 ==l && a2 == m) {
                  if (a1 == a2) qphen[i][j][k]=rannor(mean[k][cp-1],vari[k]);
                  else qphen[i][j][k]=rannor(mean[k][cp-1],vari[k]*multvar[k]);
                  if (pro[i][j] > 1) {
                    if (doubled[i][pro[i][j] - 1][0] < j + 1)
                      qphen[i][j][k] = qphen[i][doubled[i][pro[i][j] - 1][0] - 1][k];
                    else {
                      if (doubled[i][pro[i][j] - 1][1] < j + 1)
                        qphen[i][j][k] = qphen[i][doubled[i][pro[i][j] - 1][1] - 1][k];
                    }
                  }
                  fprintf(pedfile, "%12.6f", qphen[i][j][k]);
                }
              }
            }
            break;
          case 1:
            if (locus[i][j][k][0] <= locus[i][j][k][1]) {
              a1 = locus[i][j][k][0];
              a2 = locus[i][j][k][1];
            } else {
              a1 = locus[i][j][k][1];
              a2 = locus[i][j][k][0];
            }
            cp = 0;
            for (l=1;l<=numall[k];l++) {
              for (m=l;m<=numall[k];m++) {
                cp++;
                if (a1==l&&a2==m) {
                  if (ranuni()<pen[k][cp-1]) aphen[i][j][k] = 2;
                  else aphen[i][j][k] = 1;
                  if (pro[i][j] > 1) {
                    if (doubled[i][pro[i][j]-1][0]<j+1)
                      aphen[i][j][k]=aphen[i][doubled[i][pro[i][j]-1][0]-1][k];
                    else {
                      if (doubled[i][pro[i][j]-1][1]<j+1)
                        aphen[i][j][k]=aphen[i][doubled[i][pro[i][j]-1][1]-1][k];
                    }
                  }
                  fprintf(pedfile, "%8d", aphen[i][j][k]);
                }
              }
            }
            break;
          case 2:
            a1 = locus[i][j][k][0];
            a2 = locus[i][j][k][1];
            for (l=0;l<nfact[k];l++) facttot[k][l]=binall[k][a1-1][l];
            for (l=0;l<nfact[k];l++)
              if (binall[k][a2-1][l]==1) facttot[k][l]=binall[k][a2-1][l];
            fprintf(pedfile, "   ");
            for (l=0;l<nfact[k];l++) fprintf(pedfile, "%2d", facttot[k][l]);
            break;
          case 3:
            fprintf(pedfile, "%5d%3d", locus[i][j][k][0], locus[i][j][k][1]);
            break;
          }
        } else {
          switch (locustype[k]) {
          case 0: fprintf(pedfile, "     0.0    ");break;
          case 1: fprintf(pedfile, "       0");break;
          case 2: fprintf(pedfile, "   ");
            for (l=1;l<=nfact[k];l++) fprintf(pedfile, " 0");
            break;
          case 3:fprintf(pedfile, "    0  0");break;
          }
        }
      }
      putc('\n', pedfile);
      simmed[i][j+1] = false;
    }
  }
}

int main(int argc, char *argv[])
{
  FILE *simout, *problem;
  int i, j=0;
  char simoutn[30], problemn[30];

  initialize();
  if (bad) return;
  printf("\nLoci and interlocus male(/female) recombination fractions:\n");
  for (i = 1; i <= numloci; i++) {
    printf("%3d", i);
    if (i < numloci) printf("%7.3f", theta[i][0]);
    if (sexdiff > 0 && i < numloci) printf("/%7.3f", theta[i][1]);
  }
  printf("\nOpening \"pedfile.dat\" file for output\n");
  strcpy(pedfilen, "pedfile.dat");
  pedfile=fopen(pedfilen,"w");
  if (pedfile == NULL) perror("Open file error");
  printf("Generating replicates.  Please wait...\n");
  numrepdone = 0;
  for (ss = 1; ss <= numreps; ss++)
  {
    if (ss % howoften == 0) {
      numrepdone += howoften;
      printf("%d replicates done\n", numrepdone);
    }
    replicate();
  }
  fclose(pedfile);

/*write SIMOUT.DAT*/

  printf("Opening \"simout.dat\" file for output\n");
  strcpy(simoutn, "simout.dat");
  simout = fopen(simoutn, "w");
  if (simout == NULL) perror("Open file error");
  fprintf(simout," The random number seeds are: %d %d %d\n",seedorig1, seedorig2, seedorig3);
  fprintf(simout," The number of replicates is:%13d\n", numreps);
  fprintf(simout," The requested proportion of unlinked families is:  0.000\n");
  fprintf(simout," The trait locus is number:   1\n");
  fprintf(simout,"    Summary Statistics about simped.dat\n");
  fprintf(simout," Number of Pedigrees%10d\n", numpeds);
  fprintf(simout," Number of People");
  for (i=0;i<numpeds;i++) j+=numind[i];
  fprintf(simout, "%13d\n", j);
  fclose(simout);

/*write PROBLEM.DAT*/

  printf("Opening \"problem.dat\" file for output\n");
  strcpy(problemn, "problem.dat");
  problem=fopen(problemn,"w");
  if (problem==NULL) perror("Open file error");
  fprintf(problem, "%d %d %d", seed1, seed2, seed3);
  fprintf(problem, "  %d", numreps);
  if (seqrun>-1) fprintf(problem, "  %d", ++seqrun);
  putc('\n', problem);
  printf("\nOutput files created (analogous to SLINK files):\n");
  printf("  PEDFILE.DAT   simulated pedigree data\n");
  printf("  SIMOUT.DAT    summary on simulation\n");
  printf("  PROBLEM.DAT   containing updated random seeds and number of replicates\n");
  fclose(problem);
}
