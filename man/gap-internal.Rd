\name{gap-internal}
\alias{PARn}
\alias{HapDesign}
\alias{HapFreqSE}
\alias{a2g}
\alias{ci}
\alias{Cox.est}
\alias{Cox.T}
\alias{cov.invlogit}
\alias{Dev}
\alias{Dev.est}
\alias{DevH0dominant}
\alias{DevH0dominant.est}
\alias{DevH0recessive}
\alias{DevH0recessive.est}
\alias{DevHaGdominant}
\alias{DevHaGdominant.est}
\alias{DevHaGrecessive}
\alias{DevHaGrecessive.est}
\alias{DevHa.recessive}
\alias{DevHa.recessive.est}
\alias{ExpN.recessiveH0}
\alias{ExpN.recessiveHa}
\alias{g2a}
\alias{g2a.c}
\alias{gc.control}
\alias{gcode}
\alias{geno.p}
\alias{getb1star}
\alias{getPTE}
\alias{grec2g}
\alias{hap.control}
\alias{hap.score.glm}
\alias{hap.score.podds}
\alias{hmht.control}
\alias{hw}
\alias{invlogit}
\alias{is.miss}
\alias{KCC}
\alias{k}
\alias{logit}
\alias{m2plem}
\alias{mht.control}
\alias{micombine}
\alias{plem2m}
\alias{revhap}
\alias{revhap.i}
\alias{se.exp}
\alias{se.invlogit}
\alias{solve_skol}
\alias{toETDT}
\alias{ungcode}
\alias{x2}
\alias{z}
\title{Internal functions for gap}
\usage{
a2g(a1,a2)
g2a(g)
g2a.c(g)
KCC(model,GRR,p1,K)
revhap(loci,hapid)
}
\description{
These are internal functions.

PARn calculates population attributable risk (PAR) for a list of frequencies and relative risks (RRs).

HapDesign and HapFreqSE both accept a \code{\link[haplo.stats]{haplo.em}} object to derieve a design/dosage
matrix and standard error of haplotype frequency estimates. The former is appropriate for haplotype trend 
regression (HTR), e.g., within the generalized linear model (GLM) framework to be equivllant to a formal 
approach as implemented in the package haplo.stats and hap.score. However, they are expected to be compatible 
with objects from gc.em() \code{\link[gap]{gc.em}} and \code{\link[gap]{hap.em}}. The two functions are 
included as courtesy of Prof Andrea Foulkes from the useR!2008 tutorial.

a2g gives allele-to-genotype conversion.

g2a returns two alleles according to a genotype identifier.

g2a.c is the C version of g2c.

gc.control() is used by gc.em().

gcode is as a2g.

grec2g is undocumented.

hap.control() is used by hap.em().

hap.score.glm, hap.score.podds are used by hap.score().

is.miss is undocumented.

KCC calculates disease prevalences in cases and controls for a given genotype relative risk, 
allele frequency and prevalencen of the disease in the whole population. It is used by tscc
and pbsize2.

k obtains 1st and 2nd order culumants for correlation coefficient.

m2plem is an experimental function for PLEM format.

micombine is used to combine imputation results.

plem2m is also experimental for PLEM format.

revhap recovers the allele indices for a given haplotype ID in a multiallelic system.

revhap.i is similar to revhap.

solve.skol is a function used by tscc.

toETDT a function used to experiment with ETDT.

ungcode recovers alleles from genotype(s).

x2 is a simple chi-squared test of two proportions.

z is a normal z-test of two proportions used by tscc.

}
\arguments{
  \item{a1}{Allele 1}
  \item{a2}{Allele 2}
  \item{g}{A genotype identifier}
  \item{model}{One of "multiplicative", "additive", "recessive", "dominant", "overdominant"}
  \item{GRR}{Genotype relative risk}
  \item{p1}{Frequency of the risk allele}
  \item{K}{Prevalence of disease in the population}
  \item{loci}{A vector of number of alleles at all loci}
  \item{hapid}{Haplotype identifier}
}
\details{These functions are not so frequently called by users}
\keyword{internal}
