\name{mhtplot}
\alias{mhtplot}
\title{Manhattan plot of p values}
\usage{mhtplot(data, usepos=FALSE, logscale=TRUE, base=10, cutoffs=c(4,6,8), colors=NULL, labels=NULL, gap=NULL, ...)}
\description{
To generate Manhattan plot of genomewide significance (p values). It could also be
used for any random variable that is uniformly distributed. By default, a log10-transformation is applied so 
that a different set of cutoffs can be more appropriate. Note that with real chromosomal positions, it is also
appropriate to plot and some but not all chromosomes.

It is possible to specify options such as xlab and ylim when the plot is requested for data in other context.

}
\arguments{
\item{data}{a data frame with three columns representing chromosome, position and p values}
\item{usepos}{a flag to use real chromosomal positions as composed to ordinal positions}
\item{logscale}{a flag to indicate if p value is to be log-transformed}
\item{base}{the base of the logarithm}
\item{cutoffs}{the cut-offs where horizontal line(s) are drawn}
\item{colors}{the color for different chromosome(s), and random if unspecified}
\item{labels}{labels for the x-axis}
\item{gap}{gap between chromosomes}
\item{...}{other options in compatible with the R plot function}
}

\value{
The plot is shown on or saved to the appropriate device.
}

\seealso{
\code{\link[gap]{qqunif}}
}

\examples{
\dontrun{
# foo example
test <- matrix(c(1,1,4,1,1,6,1,10,3,2,1,5,2,2,6,2,4,8),byrow=TRUE,6)
mhtplot(test)
mhtplot(test,logscale=F)

# fake example with Affy500k data
affy <-c(40220, 41400, 33801, 32334, 32056, 31470, 25835, 27457, 22864, 28501, 26273, 
         24954, 19188, 15721, 14356, 15309, 11281, 14881, 6399, 12400, 7125, 6207)
CM <- cumsum(affy)
n.markers <- sum(affy)
n.chr <- length(affy)
test <- data.frame(chr=rep(1:n.chr,affy),pos=1:n.markers,p=runif(n.markers))

# to reduce size of the plot
# bitmap("mhtplot.bmp",res=72*5)
oldpar <- par()
par(las="2",cex=0.6)
colors <- rep(c("blue","green"),11)
# other colors, e.g.
# colors <- c("red","blue","green","cyan","yellow","gray","magenta","red","blue","green",
#             "cyan","yellow","gray","magenta","red","blue","green","cyan","yellow","gray","magenta","red")
mhtplot(test,colors=colors,pch=19,bg=colors)
title("A simulated example according to EPIC-Norfolk QCed SNPs")
par(cex.axis=1.3)
mhtplot(test,usepos=TRUE,colors=colors,gap=10000,pch=19,bg=colors)
title("Real positions with a gap of 10000 bp between chromosomes")
box()
par(oldpar)
# dev.off()
}
}
\author{Jing Hua Zhao}
\keyword{hplot}
