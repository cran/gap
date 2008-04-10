\name{gap-internal}
\alias{PARn}
\alias{a2g}
\alias{g2a}
\alias{g2a.c}
\alias{gc.control}
\alias{gcode}
\alias{grec2g}
\alias{hap.control}
\alias{hap.score.glm}
\alias{hap.score.podds}
\alias{is.miss}
\alias{KCC}
\alias{m2plem}
\alias{micombine}
\alias{plem2m}
\alias{revhap}
\alias{revhap.i}
\alias{toETDT}
\alias{ungcode}
\title{internal functions for gap}
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

m2plem is an experimental function for PLEM format.

micombine is used to combine imputation results.

plem2m is also experimental for PLEM format.

revhap recovers the allele indices for a given haplotype ID in a multiallelic system.

revhap.i is similar to revhap.

toETDT a function used to experiment with ETDT.

ungcode recovers alleles from genotype(s).

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
