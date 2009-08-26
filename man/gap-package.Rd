\name{gap-package}
\alias{gap-package}
\alias{gap}
\docType{package}
\title{
Genetic analysis package
}
\description{
It is designed as an integrated package for genetic data analysis
of both population and family data. Currently, it contains
functions for sample size calculations of both population-based
and family-based designs, classic twin models, 
probability of familial disease aggregation, kinship calculation, 
some statistics in linkage analysis, and association analysis 
involving one or more genetic markers including haplotype 
analysis with or without environmental covariates.
}
\details{
\tabular{ll}{
Package: \tab gap\cr
Version: \tab 1.0-20\cr
Depends: \tab R(>= 2.0.0)\cr
Suggests: \tab Design, Hmisc, haplo.stats, magic, pedigree, MASS\cr
License: \tab GPL (>=2)\cr
URL: \tab http://www.mrc-epid.cam.ac.uk/~jinghua.zhao\cr
}

Index:
\preformatted{
BFDP                    Bayesian false-discovery probability
FPRP                    False-positive report probability
SNP                     Functions for single nucleotide polymorphisms (SNPs)
aldh2                   ALDH2 markers and Alcoholism
apoeapoc                APOE/APOC1 markers and Schizophrenia
asplot                  Regional association plot
bt                      Bradley-Terry model for contingency table
b2r                     Obtain correlation coefficients and their variance-covariances
ccsize                  Power and sample size for case-cohort design
chow.test               Chow's test for heterogeneity in two regressions
cf                      Cystic Fibrosis data
comp.score              score statistics for testing genetic linkage of quantitative trait
crohn                   Crohn disease data
ESplot                  Effect-size plot
fa                      Friedreich Ataxia data
fbsize                  Sample size for family-based linkage and association design
fsnps                   A case-control data involving four SNPs with missing genotype
gc.em                   Gene counting for haplotype analysis
gcontrol                genomic control
gcontrol2               genomic control based on p values
gcp                     Permutation tests using GENECOUNTING
genecounting            Gene counting for haplotype analysis
gif                     Kinship coefficient and genetic index of familiality
hap                     Haplotype reconstruction
hap.em                  Gene counting for haplotype analysis
hap.score               Score Statistics for Association of Traits with Haplotypes
hla                     HLA markers and Schizophrenia
htr                     Haplotype trend regression
hwe                     Hardy-Weinberg equlibrium test for multiallelic marker
hwe.hardy               Hardy-Weinberg equlibrium test using MCMC
kin.morgan              kinship matrix for simple pedigree
LD22                    LD statistics for two diallelic markers
LDkl                    LD statistics for two multiallelic markers
makeped                 A function to prepare pedigrees in post-MAKEPED format
mao                     A study of Parkinson's disease and MAO gene
metap                   Meta-analysis of p values
metareg                 Fixed and random effects model for meta-analysis
mhtplot                 Manhattan plot of p values
mia                     multiple imputation analysis for hap
mtdt                    Transmission/disequilibrium test of a multiallelic marker
muvar                   Means and variances under 1- and 2- locus (diallelic) QTL model
mvmeta                  Multivariate meta-analysis based on generalized least squares
nep499                  A study of Alzheimer's disease with eight SNPs and APOE
pbsize                  Power for population-based association design
pbsize2                 Power for case-control association design
pedtodot                Converting pedigree(s) to dot file(s)
pfc                     Probability of familial clustering of disease
pfc.sim                 Probability of familial clustering of disease
pgc                     Preparing weight for GENECOUNTING
plot.hap.score          Plot Haplotype Frequencies versus Haplotype Score Statistics
print.hap.score         Print a hap.score object
qqfun                   Quantile-comparison plots
qqunif                  Q-Q plot for uniformly distributed random variable
s2k                     Statistics for 2 by K table
snca                    A study of Parkinson's disease and SNCA makers
tscc                    Power calculation for two-stage case-control design
twinan90                Classic twin models
whscore                 Whittemore-Halpern scores for allele-sharing
}

We have incorporated functions for a wide range of problems. Nevertheless,
this largely remains as a preliminary work to be consolidated in the
near future.

}
\author{
Author: Jing Hua Zhao in collaboration with other colleagues, and with
        help from Kurt Hornik and Brian Ripley of the R core development team

Maintainer: Jing Hua Zhao <jinghua.zhao@mrc-epid.cam.ac.uk>
}
\references{
Zhao JH, gap: genetic analysis package. Journal of Statistical Software 2007, 23(8):1-18
}
\keyword{package}
