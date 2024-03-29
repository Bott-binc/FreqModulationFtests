% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/userFunctions.R
\name{F3Testpar}
\alias{F3Testpar}
\title{F3TestParallel}
\usage{
F3Testpar(
  xt,
  k,
  p,
  N = length(xt),
  deltat = 1,
  dpss = FALSE,
  undersampleNumber = 100,
  penalty = 1,
  penaltyType = "ScaledExp",
  R = 1,
  cores = 1,
  confLevel = (1 - (1/length(xt))),
  returnFTestVars = FALSE,
  penaltyOnTapersStdInv = FALSE
)
}
\arguments{
\item{xt}{time series}

\item{k}{vector of tapers used in the f test}

\item{p}{Highest degree polynomial you want to test for}

\item{N}{= length(xt)Total number of observations}

\item{deltat}{Time interval between each observation}

\item{dpss}{= FALSE unless you want to use dpss, it will do sine tapers by default}

\item{undersampleNumber}{A numeric of the number the user wants to undersample, usually 100 is a good start}

\item{penalty}{1 is no penalty , 0.2  would give seq(from =  1, to =  1/(0.2*k, length.out = k)
penalty to each respective taper}

\item{penaltyType}{What type of penalty you want to use, "ScaledExp" is the most harsh and the best right now,
"mtm" is for adaptive multitaper weighting, "Cos" is for a cosine weighting scheme, "Clip" is 1 for the number passed into penalty
k's then is 0 for the rest.  The percentage is specified by a fraction in the penalty variable}

\item{R}{= 1 by default: Number of times the specific frequency must be detected to be considered significant}

\item{cores}{must be 1 if on windows, number of cores used for parallelization}

\item{confLevel}{default is 1-1/N, level of confidence used in the Ftest}

\item{returnFTestVars}{= FALSE.  Used for when more information about inner Ftests is needed.
NOTE: this will be very memory intensive and is only implimented for altSig = FALSE}

\item{penaltyOnTapersStdInv}{applys the penalty to all tapers in the calculation instead of just weighting the
eigenCoef's}
}
\value{
$F3testStat, $Freq, $sigFreq, prop, aggrTestResult zero if fail to reject, 1 if rejected at the specified R
}
\description{
This function still works, and I kept it so the older code would still work, but should be using the independent functions in practice for a specific given fTest.
w is chosen by shannons number based on k
}
