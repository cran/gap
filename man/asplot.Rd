\name{asplot}
\alias{asplot}
\title{Regional association plot}
\usage{asplot(snp, locusname, chr, locus, range, best.pval, pch, recomb, genelist)}
\description{
This function obtains regional association plot for a particular locus, based on
the information about recombinatino rates, linkage disequilibrium between the
SNP of interest and neighbouring ones, and single-point association tests p values.
}

\arguments{
 \item{snp}{The SNP name, e.g., "rs10811661"}
 \item{locusname}{The classic term for locus where a gene locates}
 \item{chr}{The name of the chromosome, e.g., "1", \ldots, "X"}
 \item{locus}{The data frame containing association results}
 \item{range}{The range for -log10(p), e.g., 9}
 \item{best.pval}{The best p value for the locus of interest}
 \item{pch}{Plotting character for the SNPs to be highlighted, e.g., 21 and 23 refer to circle and diamond}
 \item{recomb}{The genetic map}
 \item{genelist}{The gene annotation}
}

\examples{
\dontrun{
data(CDKN)
asplot("rs10811661", "CDKN2A/CDKN2B region", "9", CDKNlocus, 9, 5.4e-8, 23, CDKNmap, CDKNgenes)
}
}
\author{Paul de Bakker, Jing Hua Zhao}
\keyword{hplot}
