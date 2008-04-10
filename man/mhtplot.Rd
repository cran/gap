\name{mhtplot}
\alias{mhtplot}
\title{Manhattan plot of p values}
\usage{mhtplot(data,logscale=TRUE, base=10, cutoffs=c(3,5,7,9), color=NULL,
              labels=paste(1:22,sep=""),...)}
\description{
To generate Manhattan plot of genomewide significance (p values). It could also be
used for any random variable that is uniformly distributed.
}
\arguments{
\item{data}{a data frame with three columns representing chromosome, position and p values}
\item{logscale}{a flag to indicate if p value is to be log-transformed}
\item{base}{the base of the logarithm}
\item{cutoffs}{the cutt-offs where horizontal line(s) are drawn}
\item{color}{the color for different chromosome(s), and random if unspecified}
\item{labels}{labels for the x-axis}
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
# fake example with Affy500k data
affy <-c(40220, 41400, 33801, 32334, 32056, 31470, 25835, 27457, 22864, 28501, 26273, 
         24954, 19188, 15721, 14356, 15309, 11281, 14881, 6399, 12400, 7125, 6207)
CM <- cumsum(affy)
n.markers <- sum(affy)
n.chr <- length(affy)
test <- data.frame(chr=rep(1:n.chr,affy),pos=1:n.markers,p=runif(n.markers))

#to reduce the massive size of the plot
#bitmap("mhtplot.bmp",res=72*5)
par(las="2",cex=0.6)
colors <- rep(c("blue","green"),11)
mhtplot(test,color=colors,labels=paste(1:n.chr,sep=""),pch=21,bg=colors)
title("A simulated example according to EPIC-Norfolk QCed SNPs")
#dev.off()
}
}
\author{Jing Hua Zhao}
\keyword{hplot}
