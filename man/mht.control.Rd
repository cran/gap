% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mhtplot.R
\name{mht.control}
\alias{mht.control}
\title{Controls for mhtplot}
\usage{
mht.control(
  type = "p",
  usepos = FALSE,
  logscale = TRUE,
  base = 10,
  cutoffs = NULL,
  colors = NULL,
  labels = NULL,
  srt = 45,
  gap = NULL,
  cex = 0.4,
  yline = 3,
  xline = 3
)
}
\arguments{
\item{type}{Type of plot.}

\item{usepos}{A flag.}

\item{logscale}{A flag for log-scale.}

\item{base}{Base of log.}

\item{cutoffs}{Cutoffs of P-value, etc.}

\item{colors}{Colours for chromosomes.}

\item{labels}{Labels for chromosomes.}

\item{srt}{Rotation degrees.}

\item{gap}{Gap between data points.}

\item{cex}{Scaling factor of data points.}

\item{yline}{Vertical adjustment.}

\item{xline}{Horiztonal adjustment.}
}
\value{
A list as above.
}
\description{
Parameter specification through function
}
