#' Generate simulated weather record.
#'
#' Generate time series of temperature measurements
#'
#' @param timebase A vector of dates / times in decimal julian day. Either this
#'   or \code{years} must be supplied.
#' @param years A scalar, the number of years to model. Either this or
#'   \code{timebase} must be supplied.
#' @param spectrum A data frame of \code{cbind(frequency, cyc_range, phase,
#'   tau)}. Day frequency = 1, year frequency = 1/325.25. Cycling range = 2 *
#'   amplitude. Phase is -pi to pi, as output by nlts_plus_phase. tau is lag,
#'   out output by nlts_plus_phase.
#' @param mean The mean temperature. Use \code{t_intercept} for trended data.
#'   Default 0.
#' @param t_int T intercept for use with linear trend.
#' @param t_slope Slope for use with linear trend.
#' @param mean_resid Mean residual. Currently rnorm with sd = 1, without
#'   autocorrelation. Next version will have ARIMA autocorrelation.
#' @return a numeric vector of the number of days since 1/1/1960, dates before
#'   this are negative.
#' @export
#' @examples
#' gen_cycling_rec(years = 3,
#'                 spectrum = data.frame(frequency = c(1, 365.25)
#'                                       cyc_range = c(0.25, )))
gen_cycling_rec <- function(timebase   = NULL,
                            years      = NULL,
  spectrum = stop("data frame of data.frame(frequency, cyc_range, phase, tau) required"),
                            mean       = 0,
                            t_int      = NULL,
                            t_slope    = NULL,
                            mean_resid = 0) {
  if (is.null(timebase)) {
    if (is.null(years)) {stop("Either timebase or years must be supplied")}
    timebase = (1:(years * 365 * 24)) / 24
  }
  #add mean
  temperatures   <- rep(mean, length(timebase))
  #add trend
  if (!is.null(t_int)) {
    if (is.null(t_slope)) {stop("If t_int is supplied, t_slope must also be suppled")}
    temperatures <- temperatures + ((t_slope * timebase) + t_intercept)
  }
  #loop through adding frequencies.
  for (i in 1:dim(spectrum)[1]) {#Note: looped to save memory, a vectorized version could be big.
    temperatures <- temperatures + ((spectrum$cyc_range[i] / 2) *
                                      cos(2 * pi * spectrum$frequency[i]) *
                                      (timebase - spectrum$tau[i]) +
                                      spectrum$phase[i])
  }
  #add residuals
  temperatures <- temperatures + rnorm(length(temperatures),
                                       mean = mean_resid,
                                       sd   = 1)
  return(temperatures)
}
