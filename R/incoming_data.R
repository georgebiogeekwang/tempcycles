#' mm/dd/yyyy to julian date helper function
#'
#' The \pkg{rnoaa} package returns dates in the station list in "mm/dd/yyyy"
#' format, and we (Wang & Dillon) formatted our dates in a similar way. All date
#' and spectral calculations work on multi-year Julian dates. This is a small
#' helper function to convert from "mm/dd/yyyy" to julian day. Because it uses
#' \code{mdy.date()} from the \pkg{date} package, it returns the number of days
#' since January 1, 1960. Inputs are coerced to character if necessary.
#'
#' @param mmddyyyy_dates a vector of dates in "mm/dd/yyyy" format.
#' @return a numeric vector of the number of days since 1/1/1960, dates before
#'   this are negative.
#' @export
#' @examples
#' mmddyyyy_to_numeric(c("3/5/1954", "9/28/2012"))
mmddyyyy_to_numeric <- function(mmddyyyy_dates = stop("Char dates mm/dd/yyyy rqrd.")) {
  if (class(mmddyyyy_dates) != "character") {
    mmddyyyy_dates <- as.character(mmddyyyy_dates)
  }
  splitdates <- strsplit(mmddyyyy_dates,
                         "/")
  l <- lapply(splitdates,
              function(X) c(X, rep(NA, 3 - length(X))))
  mdy_df <- data.frame(t(do.call(cbind, l)))
  names(mdy_df) <- c("m", "d", "y")
  mdy_df$m <- as.numeric(levels(mdy_df$m))[mdy_df$m] #values were converted to factors
  mdy_df$d <- as.numeric(levels(mdy_df$d))[mdy_df$d]
  mdy_df$y <- as.numeric(levels(mdy_df$y))[mdy_df$y]
  jday <- date::mdy.date(mdy_df$m,
                         mdy_df$d,
                         mdy_df$y)
  jday <- as.numeric(jday)
  return(jday)
}

#' yyyy-mm-dd to julian date helper function
#'
#' The \pkg{rnoaa} package returns dates for weather station records in
#' "yyyy-mm-dd" format. All date and spectral calculations work on multi-year
#' Julian dates. This is a small helper function to convert from "yyyy-mm-dd" to
#' julian day. Because it uses \code{mdy.date()} from the \pkg{date} package, it
#' returns the number of days since January 1, 1960. Inputs are coerced to
#' character if necessary.
#'
#' @param yyyymmdd_dates a vector of dates in "yyyy-mm-dd" format.
#' @return a numeric vector of the number of days since 1/1/1960, dates before
#'   this are negative.
#' @export
#' @examples
#' yyyymmdd_to_numeric(c("1954-3-5", "2012-9-28"))
yyyymmdd_to_numeric <- function(yyyymmdd_dates = stop("Char dates yyyy-mm-dd rqrd.")) {
  if (class(yyyymmdd_dates) != "character") {
    yyyymmdd_dates <- as.character(yyyymmdd_dates)
  }
  splitdates <- strsplit(yyyymmdd_dates,
                         "-")
  l <- lapply(splitdates,
              function(X) c(X, rep(NA, 3 - length(X))))
  mdy_df <- data.frame(t(do.call(cbind, l)))
  names(mdy_df) <- c("y", "m", "d")
  mdy_df$m <- as.numeric(levels(mdy_df$m))[mdy_df$m] #values were converted to factors
  mdy_df$d <- as.numeric(levels(mdy_df$d))[mdy_df$d]
  mdy_df$y <- as.numeric(levels(mdy_df$y))[mdy_df$y]
  jday <- date::mdy.date(mdy_df$m,
                         mdy_df$d,
                         mdy_df$y)
  jday <- as.numeric(jday)
  return(jday)
}

#' Find weather stations that are close to a specified location.
#'
#' This is the underlying distance function for finding close weather records.
#' If neither \code{radius} or \code{n} is specified, the entire station_data
#' data frame will be returned.
#'
#' @param sample_coords The \code{cbind(lon, lat)} of the experiment, in decimal
#'   degrees, negative south and west.
#' @param station_data A data frame of the weather station metadata, with
#'   \code{lon, lat} columns, in decimal degrees, negative south and west.
#' @param radius The search limit in km. Default: no limit.
#' @param n The number of stations to return. Default: no limit.
#' @return A subset of the original station_data data frame, with an additional
#'   \code{distance} column (km), sorted by distance.
#' @export
#' @examples
#' require(tempcyclesdata)
#' close_stations(sample_coords = c(9.0556, 48.52),
#'                station_data = tempcyclesdata,
#'                radius = 100,
#'                n      = 20)
close_stations <- function(sample_coords = stop("Sample coords c(lon,lat) required."),
                           station_data = stop("Station data with lon, lat cols rqrd."),
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
      stop(paste("No stations within", radius, "km. Try increasing search radius."))
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

#' Find closest data from Wang and Dillon 2014 to a specified location within a
#' time interval.
#'
#' Search through the Wang and Dillon cycling dataset, and find the closest
#' stations to a specified sample site.
#' @param samp_coords The \code{cbind(lon, lat)} of the search origin, in
#'   decimal degrees, negative south and west.
#' @param samp_interval Period over which to search for data \code{c(start,end),
#'   "mm/dd/yyyy"}. Default returns all.
#' @param ... additional parameters to limit search (see \code{close_stations}).
#' @return a subset of the cycling data, with and additional column of distance
#'   from sample, sorted by distance to the specified location.
#' @seealso \code{\link{close_stations}} for constraining search.
#' @export
#' @examples
#' close_cycling_data(samp_coords = c(9.0556, 48.52),
#'                    samp_interval = c("1/1/1998", "4/7/2007"),
#'                    radius = 50,
#'                    n = 2)
close_cycling_data <- function(samp_coords = stop("Sample coords c(lon,lat) required."),
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
#'
#' @param samp_coords The \code{cbind(lon, lat)} of the search origin, in decimal
#'   degrees, negative south and west.
#' @param samp_interval Period over which to search for data \code{c(start,end),
#'   "mm/dd/yyyy"}. Default returns all.
#' @param ... additional parameters to limit search (see \code{close_stations}).
#' @return a data frame with \code{usaf} and \code{wban} identifiers for use
#'   with \code{get_noaaisd_data}.
#' @seealso \code{\link{close_stations}} for constraining search.
#' @export
#' @examples
#' close_noaaisd_data(samp_coords = c(9.0556, 48.52),
#'                    samp_interval = c("1/1/1998", "4/7/2007"),
#'                    radius = 50,
#'                    n = 2)
close_noaaisd_data <- function(samp_coords = stop("Sample coords c(lon,lat) required."),
                               samp_interval = NULL,
                               ...) {
  station_list <- rnoaa::isd_stations()
  station_list <- station_list[-(grep("BOGUS",
                                      station_list$station_name)),]
  station_list <- station_list[complete.cases(station_list$lat, station_list$lon),]
  station_list <- station_list[-(which(station_list$lon < -180)),]
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

#' Download and temperature data NOAA ISD and format for spectral analysis.
#'
#' A helper script for the \code{rnoaa::isd} function. Download the specified
#' weather station for the specified years, remove errors and low quality data,
#' convert to julian day, and concatenate into a single data frame. Temperature
#' values > 900 are error marks and are removed. Only values which pass all
#' quality control tests, (quality code "1"), are retained.
#'
#' @param usaf A string of the USAF weather station identifier.
#' @param wban A string of the WBAN weather station identifier.
#' @param years A vector of years to download. Default returns all, missing data
#'   results in a warning.
#' @return A data frame with \code{jday} and \code{temperature}.
#' @export
#' @examples
#' get_noaaisd_data(usaf = "702700",
#'                  wban = "00489",
#'                  years = c(2009,2013:2015))
get_noaaisd_data <- function(usaf = stop("6-digit char USAF id required."),
                             wban = stop("5-digit char WBAN id required."),
                             years = NULL) {
  if (is.null(years)) {
    station_list <- rnoaa::isd_stations()
    #Note: rnoaa currently returns station ids as integers!
    if (is.integer(station_list$usaf)) {
      station <- station_list[which(station_list$usaf == as.integer(usaf) &
                                    station_list$wban == as.integer(wban)),]
    } else {
      station <- station_list[which(station_list$usaf == usaf &
                                    station_list$wban == wban),]
    }
    record_start <- floor(station$begin / 10000)
    record_end   <- floor(station$end / 10000)
    years        <- record_start:record_end
  }
  for (year in years) {
    station_data_year <- tryCatch(
      rnoaa::isd(usaf = usaf,
                 wban = wban,
                 year = year),
      error = function(e) {
        warning(paste(year, " not found for ", usaf, "-", wban, sep = ""))
      } )
    if(is.character(station_data_year)) {next}
    year_data <- data.frame(jday = yyyymmdd_to_numeric(station_data_year[[1]]$date) +
                                  ((as.numeric(station_data_year[[1]]$time) / 100) / 24),
                            temperature = station_data_year[[1]]$temperature,
                            quality     = station_data_year[[1]]$temperature_quality)
    year_data <- subset(year_data,
                        year_data$temperature < 900 & year_data$quality == "1")
    year_data$quality <- NULL
    if (exists("station_data")) {
      station_data <- rbind(station_data,
                            year_data)
    } else {
      station_data <- year_data
    }
  }
  return(station_data)
}

#' Preliminary checks for temperature records.
#'
#' Check if a weather record is suitable for spectral analysis. Checks are based
#' on sampling: records must sample at least every six hours on average, be at
#' least 1.5 years in length, and not have a single gap more than ten percent of
#' the record length. As suggested by the function name, these checks are
#' preliminary. The significance of the peaks should also be measured, either by
#' comparison to a white noise model, as in \pkg{nlts}, or against Monte Carlo
#' or Chi-square AR1 levels as in redfit, available in \pkg{dplR}. Most weather
#' records that pass these checks and are not seasonally sampled will have
#' significant DTC and ATC values.
#'
#' @param record_time A time vector, in fractional julian date.
#' @return A logical data frame with check values.
#' @export
#' @examples
#' weather_record <- get_noaaisd_data(usaf = "702700",
#'                                    wban = "00489",
#'                                    years = 2014:2015)
#' prelim_record_check(weather_record$jday)
prelim_record_check <- function(record_time) {
  data_checks <- data.frame(c_mean_diff_lt_6h             = FALSE,
                            c_at_least_one_and_half_years = FALSE,
                            c_no_big_gap                  = FALSE)
  diffvec       <- diff(record_time)
  record_length <- max(record_time) - min(record_time)
  if (mean(diffvec < 0.25)) {
    data_checks$c_mean_diff_lt_6h <- TRUE
  }
  if (record_length > (365 * 1.5)) {
    data_checks$c_at_least_one_and_half_years <- TRUE
  }
  if (max(diffvec) < (0.1 * record_length)) {
    data_checks$c_no_big_gap <- TRUE
  }
  return(data_checks)
}
