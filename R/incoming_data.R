#' Find weather stations that are close to a specified location.
#'
#' This is the underlying distance function for finding close weather records.
#' If neither \code{radius} or \code{n}, is specified, the entire station_data
#' dataframe will be returned.
#'
#' @param sample_coords The \code{cbind(lon, lat)} of the experiment, in decimal
#'   degrees, negative south and west.
#' @param station_data A data frame of the weather station metadata, with
#'   \code{lon, lat} columns, in decimal degrees, negative south and west.
#' @param radius The search limit in km. Default: no limit.
#' @param n The number of stations to return. Default: no limit.
#' @return The original station_data dataframe, with an additional
#'   \code{distance} column (km), sorted by distance.
#' @examples
#' require(tempcyclesdata)
#' close_stations(cbind(9.0556, 48.52),
#'                tempcyclesdata,
#'                radius = 100,
#'                n      = 20)
close_stations <- function(sample_coords,
                           station_data,
                           radius = NULL,
                           n      = 10) {
  #calculate distances in km
  station_data$distance <- geosphere::distVincentyEllipsoid(cbind(station_data$lon,
                                                                  station_data$lat),
                                                            sample_coords) / 1000
  #sort by distance
  station_data <- station_data[order(station_data$distance),]
  #find close points
  if (!is.null(radius)) {
    within_radius <- which(station_data$distance <= radius)
    if (length(within_radius) == 0) {
      stop(paste("No stations within", radius, "km."))
    }
    station_data <- station_data[within_radius,]
  }
  #limit output number
  if (!is.null(n)) {
    if (dim(station_data)[1] > n) {
      station_data <- station_data[1:n,]
    }
  }
  return(station_data)
}

#' Find closest data from Wang & Dillon 2014 to a specified location and time.
#'
#' Search through the Wang & Dillon cycling dataset, and find the closest
#' stations to a specified sample site.
#' @param samp_coords The \code{cbind(lon, lat)} of the experiment, in decimal
#'   degrees, negative south and west.
#' @param samp_interval period over which to search for data \code{c(start,end),
#'   "mm/dd/yyyy"}. Default returns all.
#' @param ... Additional parameters to limit search (see \code{close_stations}).
#' @return a subset of the cycling data, with and additional column of distance
#'   from sample.
#' @examples
#' close_cycling_data(c(9.0556, 48.52), c("1/1/1998", "4/7/2007"), radius = 50, n = 2)
close_cycling_data <- function(samp_coords,
                               samp_interval = NULL,
                               ...) {
  if (!is.null(samp_interval)) {
    samp_interval <- lubridate::mdy(samp_interval)
    samp_interval <- lubridate::interval(samp_interval[1],
                                         samp_interval[2])
    station_start <- lubridate::mdy(tempcyclesdata::tempcyclesdata$start_date)
    station_end   <- lubridate::mdy(tempcyclesdata::tempcyclesdata$end_date)
    station_ints  <- lubridate::interval(station_start,
                                         station_end)
    overlap_records <- which(lubridate::int_overlaps(samp_interval,
                                                     station_ints))
    records_subset  <- tempcyclesdata::tempcyclesdata[overlap_records,]
  } else {
    records_subset  <- tempcyclesdata::tempcyclesdata
  }
  unique_stations <- records_subset[match(unique(records_subset$id), records_subset$id),]
  station_matches <- close_stations(samp_coords,
                                    unique_stations,
                                    ...)
  record_matches <- records_subset[which(records_subset$id %in% station_matches$id),]
  record_matches$distance <- station_matches$distance[match(record_matches$id,
                                                            station_matches$id)]
  record_matches <- record_matches[order(record_matches$distance),]
  return(record_matches)
}

#' Find closest data from the NOAA ISD dataset to a specified location and time.
#'
#' Search through the Wang & Dillon cycling dataset, and find the closest
#' stations to a specified sample site.
#' @param samp_coords The \code{cbind(lon, lat)} of the experiment, in decimal
#'   degrees, negative south and west.
#' @param samp_interval period over which to search for data \code{c(start,end),
#'   "mm/dd/yyyy"}. Default returns all.
#' @param ... Additional parameters to limit search (see \code{close_stations}).
#' @return a data frame with \code{usaf} and \code{wban} identifiers for use
#'   with \code{get_noaaisd_data}.
#' @examples
#' close_noaa_data(c(9.0556, 48.52), c("1/1/1998", "4/7/2007"), radius = 50, n = 2)
close_noaaisd_data <- function(samp_coords,
                               samp_interval = NULL,
                               ...) {
  station_list <- rnoaa::isd_stations()
  station_list <- station_list[complete.cases(station_list$lat, station_list$lon),}
  if (!is.null(samp_interval)) {
    samp_interval <- lubridate::mdy(samp_interval)
    samp_interval <- lubridate::interval(samp_interval[1],
                                         samp_interval[2])
    station_start <- lubridate::ymd(station_list$begin)
    station_end   <- lubridate::ymd(station_list$end)
    station_ints  <- lubridate::interval(station_start,
                                         station_end)
    overlap_records <- which(lubridate::int_overlaps(samp_interval,
                                                     station_ints))
    station_subset  <- station_list[overlap_records,]
  } else {
    station_subset  <- station_list
  }
  station_matches <- close_stations(samp_coords,
                                    station_subset,
                                    ...)
  station_matches <- station_matches[order(station_matches$distance),]
  return(station_matches)
}

