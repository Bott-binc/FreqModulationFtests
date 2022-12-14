% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hiddenFunctions.R
\name{TimeDomSlepians}
\alias{TimeDomSlepians}
\title{Time domain slepians **Untested**}
\usage{
TimeDomSlepians(N, w, k, passInFullDPSS = NULL, returnDPSS = FALSE)
}
\arguments{
\item{N}{Total number of observations}

\item{w}{bandwidth parameter}

\item{k}{number of tapers}

\item{passInFullDPSS}{this is the object returned by dpss.  both $v and $eigen}

\item{returnDPSS}{returnDPSS for limiting calculations of dpss}
}
\value{
matrix of der for each t = 0, ..., N-1 by k = 0, ..., k-1
}
\description{
Calculates the time domain version of the slepians, it is not used in any of
 the functions but could be usefull to know in the future.
}
