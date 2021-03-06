% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MODIS_data.R
\name{MODIS_data}
\alias{MODIS_data}
\title{For MODIS EVI data, construct the data object for input into MCMC}
\usage{
MODIS_data(
  siteName,
  lat,
  long,
  startDate = "",
  endDate = "",
  metric,
  startDay = FALSE,
  endDay = FALSE,
  seasonOrder = "FS",
  lastYear = 2018,
  dataDirectory = ""
)
}
\arguments{
\item{siteName}{Site Name}

\item{lat}{Latitude}

\item{long}{Longitude}

\item{startDate}{The start day counted as the day number after 2016-12-31}

\item{endDate}{The end day counted as the day number after 2016-12-31}

\item{metric}{"NDVI" or "EVI"}

\item{dataDirectory}{The data directory}
}
\description{
For MODIS EVI data, construct the data object for input into MCMC
}
