\name{asplot}
\alias{asplot}
\title{Regional association plot}
\usage{asplot(snp, locusname, chr, locus, gmap, glist, best.pval=NULL, sf=c(3,8), logpmax=9, pch=23)}
\description{
This function obtains regional association plot for a particular locus, based on
the information about recombinatino rates, linkage disequilibria between the
SNP of interest and neighbouring ones, and single-point association tests p values.

Note that the best p value is not necessarily within locus in the original design.
}

\arguments{
 \item{snp}{The SNP name, e.g., "rs10811661"}
 \item{locusname}{The classic term for locus where a gene locates}
 \item{chr}{The name of the chromosome, e.g., "1", \ldots, "X"}
 \item{locus}{The data frame containing association results with columns c("POS", "PVAL", "TYPE", "RSQR"), row.names(locus) contains 
the SNP name, e.g., snp}
 \item{gmap}{The genetic map, e.g., c("position","COMBINED\_rate.cM.Mb.","Genetic\_Map.cM.")}
 \item{glist}{The gene annotation with columns c("START", "STOP", "SIZE", "STRAND", "GENE")}
 \item{best.pval}{The best p value for the locus of interest}
 \item{sf}{scale factors for p values and recombination rates, smaller values are necessary for gene dense regions}
 \item{logpmax}{The maximum value for -log10(p)}
 \item{pch}{Plotting character for the SNPs to be highlighted, e.g., 21 and 23 refer to circle and diamond}
}

\section{References}{
DGI. Whole-genome association analysis identifies novel loci for type 2 diabetes and triglyceride levels.
Science 2007;316(5829):1331-6
}

\examples{
\dontrun{
asplot("rs10811661", "CDKN2A/CDKN2B region", "9", CDKNlocus, CDKNmap, CDKNgenes)
asplot("rs10811661", "CDKN2A/CDKN2B region", "9", CDKNlocus, CDKNmap, CDKNgenes, 5.4e-8, c(3,6))
}
}
\author{Paul de Bakker, Jing Hua Zhao, Shengxu Li}
\keyword{hplot}
