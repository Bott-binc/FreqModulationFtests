% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hiddenFunctions.R
\name{reductionSingleKFullComputation}
\alias{reductionSingleKFullComputation}
\title{Full computation of the reduction including the original F3/F4 if not in fast mode.  Note that this only works for linear modulation with the current derivation.}
\usage{
reductionSingleKFullComputation(
  Xt,
  K,
  N,
  penalty = 1,
  penaltyType = "ScaledExp",
  undersample = TRUE,
  undersampleNumber = 100,
  deltat = 1,
  pad = TRUE,
  confLevel = (1 - (1/length(Xt))),
  fast = FALSE
)
}
\arguments{
\item{Xt}{time series that you want to conduct test on}

\item{K}{length of sine tapers}

\item{N}{length of time series}

\item{penalty}{= 1 by default.}

\item{penaltyType}{= "ScaledExp" by default, "Clip" also available}

\item{undersample}{= TRUE, if FALSE and using fast, it will automatically be set to TRUE and use the undersampleNumber imputted}

\item{undersampleNumber}{= 100 by default}

\item{deltat}{= 1 by default}

\item{pad}{= TRUE by default, dont need to pad manually}

\item{confLevel}{= 1-1/N by default}

\item{fast}{= FALSE by default.  Allows for a faster running of the reduction alg.  no F3 is computed and no differences or intermediate terms are returned to the user.}
}
\value{
if fast = FALSE, returns reduction test = FPrime, F3 = FTest, F3 - FPrime =  FtestDiff,
Freq , ratios = list(ratio = ratio,
             ratioPrime = ratioPrime,
             ratioDiff = ratioDiff),
significantFrequencies = list(significantFreqFull = significantFreq,
                             significantFreqPrime = significantFreqRed,
                             FCutOffFull = FCutOff,
                             FCutOffPrime = FCutOffPrime)))
}
\description{
Full computation of the reduction including the original F3/F4 if not in fast mode.  Note that this only works for linear modulation with the current derivation.
}
