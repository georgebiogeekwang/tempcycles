#' Generate simulated weather record.
#'
#' Generate time series of temperature measurements
#'
#' @param mmddyyyy_dates a vector of dates in "mm/dd/yyyy" format.
#' @return a numeric vector of the number of days since 1/1/1960, dates before
#'   this are negative.
#' @export
#' @examples
#' gen_cycling_red(c("3/5/1954", "9/28/2012"))
gen_cycling_rec <- function() {
  day.cycle.red       <- workdata$dtc.red * cos(2 * pi * 1 *
                                                  (station.data$jday -
                                                     workdata$day.tau) +
                                                  workdata$day.phase)
  year.cycle.red      <- workdata$atc.red * cos(2 * pi * 1 / 365.25 *
                                                  (station.data$jday -
                                                     workdata$year.tau) +
                                                  workdata$year.phase)
  year.model.red      <- linear.trend + year.cycle.red
  resid.y.red         <- round(sum((station.data$Ta - year.model.red)^2),
                               digits = 2)
  resid.y.red.persamp <- sqrt(resid.y.red) / n.samples
  total.model.red     <- year.model.red + day.cycle.red
  resid.t.red         <- round(sum((station.data$Ta - total.model.red)^2),
                               digits = 2)
  resid.t.red.persamp <- sqrt(resid.t.red) / n.samples
  day.cycle.red       <- workdata$dtc.red * cos(2 * pi * 1 *
                                                  (station.data$jday -
                                                     workdata$day.tau) +
                                                  workdata$day.phase)
  year.cycle.red      <- workdata$atc.red * cos(2 * pi * 1 / 365.25 *
                                                  (station.data$jday -
                                                     workdata$year.tau) +
                                                  workdata$year.phase)
  year.model.red      <- linear.trend + year.cycle.red
  resid.y.red         <- round(sum((station.data$Ta - year.model.red)^2),
                               digits = 2)
  resid.y.red.persamp <- sqrt(resid.y.red) / n.samples
  total.model.red     <- year.model.red + day.cycle.red
  resid.t.red         <- round(sum((station.data$Ta - total.model.red)^2),
                               digits = 2)
  resid.t.red.persamp <- sqrt(resid.t.red) / n.samples
}
