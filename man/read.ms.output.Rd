\name{read.ms.output}
\alias{read.ms.output}
\title{A utility function to read ms output}
\usage{
read.ms.output(msout,is.file=TRUE,xpose=TRUE)
}
\description{
This function reads in the output of the program ms, a program to generate
samples under a variety of neutral models.

The argument indicates either a file name or a vector of character strings,
one string for each line of the output of ms. As with the second case, it
is appropriate with system(,intern=TRUE), see example below.
}
\arguments{
\item{msout}{an ms output}
\item{is.file}{a flag indicating ms output as a system file or an R object}
\item{xpose}{a flag to obtain the tranposed format as it is (when TRUE)}
}
\value{
The returned value is a list storing the results.
\item{call}{system call to ms}
\item{seed}{random number seed to ms}
\item{nsam}{number of copies of the locus in each sample}
\item{nreps}{the number of independent samples to generate}
\item{segsites}{a vector of the numbers of segregating sites}
\item{times}{vectors of time to most recent ancester (TMRCA) and total tree lengths}
\item{positions}{positions of polymorphic sites on a scale of (0,1)}
\item{gametes}{a list of haplotype arrays}
\item{probs}{the probability of the specified number of segregating sites 
given the genealogical history of the sample and the value to -t option}
}
\section{References}{
Hudson RR (2002) Generating samples under a Wright-Fisher neutral model. Bioinformatics 18:337-8,

Press WH, SA Teukolsky, WT Vetterling, BP Flannery (1992). Numerical Recipes in C. Cambridge University Press, Cambridge.
}
\examples{
\dontrun{

# Assuming ms is on the path

system("ms 5 4 -s 5 > ms.out")
msout <- read.ms.output("ms.out")

msout <- system("ms 5 4 -s 5 -L", intern=TRUE)
msout <- read.ms.output(msout,FALSE)
}
}
\author{D Davison, RR Hudson, JH Zhao}
\keyword{utilities}
