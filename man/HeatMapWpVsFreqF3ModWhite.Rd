% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/userFunctions.R
\name{HeatMapWpVsFreqF3ModWhite}
\alias{HeatMapWpVsFreqF3ModWhite}
\title{Simulation of modulation width vs frequency for f3Mod}
\usage{
HeatMapWpVsFreqF3ModWhite(
  K,
  N,
  numSim = 500,
  lengthWp = 200,
  date,
  Amplitude = 1.4,
  FileDiscripter = "",
  DirForSave = "~/",
  cores = 1,
  saveData = TRUE,
  savePlot = TRUE
)
}
\arguments{
\item{K}{Number of tapers}

\item{N}{Length of time series}

\item{numSim}{number of simulations to be conducted}

\item{lengthWp}{number of wp's to check between the range of wp's tested from 4*w to 0.0001}

\item{date}{eg 221206 would be yymmdd}

\item{Amplitude}{1.4 will give SNR of 1}

\item{FileDiscripter}{What makes this sim different than the others done, like SNRHigh}

\item{DirForSave}{Where you would like the data to be saved if you want to use
you're current directory use "."}

\item{cores}{number of cores used, Windows can only use 1}

\item{savePlot}{if you want to save the .png plot to the dir as well}

\item{saveRDS}{if you would like to save the simulation results}
}
\value{
returns the simulation results
}
\description{
This was used to make the heat maps in thesis
}
