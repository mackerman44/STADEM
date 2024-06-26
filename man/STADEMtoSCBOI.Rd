% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/STADEMtoSCOBI.R
\name{STADEMtoSCBOI}
\alias{STADEMtoSCBOI}
\title{Summarise STADEM output for SCOBI}
\usage{
STADEMtoSCBOI(
  stadem_mod = NULL,
  lgr_weekly = NULL,
  saveCSV = F,
  fileName = "SCOBI_input.csv"
)
}
\arguments{
\item{stadem_mod}{results of STADEM model, of class \code{mcmc.list} or `jagsUI`.}

\item{lgr_weekly}{weekly data from Lower Granite dam, compiled by \code{compileGRAdata}.}

\item{saveCSV}{should results be saved as a .csv file? Default is \code{FALSE}.}

\item{fileName}{if \code{saveCSV} is \code{TRUE}, what file name, including the path and extension, should the results be written to?}
}
\description{
Format STADEM output into a form that can be fed into the SCOBI package for further run decomposition.
}
\author{
Kevin See
}
