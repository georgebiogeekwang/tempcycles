% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/incoming_data.R
\name{close_stations}
\alias{close_stations}
\title{Find weather stations that are close to a specified location.}
\usage{
close_stations(sample_coords = stop("Sample coords c(lon,lat) required."),
  station_data = stop("Station data with lon, lat cols rqrd."),
  radius = NULL, n = 10)
}
\arguments{
\item{sample_coords}{The \code{cbind(lon, lat)} of the experiment, in decimal
degrees, negative south and west.}

\item{station_data}{A data frame of the weather station metadata, with
\code{lon, lat} columns, in decimal degrees, negative south and west.}

\item{radius}{The search limit in km. Default: no limit.}

\item{n}{The number of stations to return. Default: no limit.}
}
\value{
A subset of the original station_data data frame, with an additional
  \code{distance} column (km), sorted by distance.
}
\description{
This is the underlying distance function for finding close weather records.
If neither \code{radius} or \code{n} is specified, the entire station_data
data frame will be returned.
}
\examples{
require(tempcyclesdata)
close_stations(sample_coords = c(9.0556, 48.52),
               station_data = tempcyclesdata,
               radius = 100,
               n      = 20)
}

