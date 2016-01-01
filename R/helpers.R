#' mm/dd/yyyy to julian date helper function
#'
#' The \code{rnoaa} returns dates in mm/dd/yyyy format, and we (Wang & Dillon)
#' formatted our dates in a similar way. All date and spectral calculations work
#' on multi-year Julian dates. This is a small helper function to convert from
#' mm/dd/yyyy to jdate. Because it uses \code{mdy.date} from the \code{date}
#' package, it returns the number of days since January 1, 1960. Inputs are
#' coerced to character if necessary.
#' @param mdy a vector of dates in "mm/dd/yyyy" format.
#' @return a numeric vector of the number of days since 1/1/1960, dates before
#'   this are negative.
#' @examples
#' mmddy7yy_to_numeric(c("3/5/1954", "9/28/2012"))
mmddyyyy_to_numeric <- function(stringdates){
  if (class(stringdates) != "character") {
    stringdates <- as.character(stringdates)
  }
  splitdates <- strsplit(stringdates,
                         "/")
  l <- lapply(splitdates,
              function(X) c(X, rep(NA, 3 - length(X))))
  mdy_df <- data.frame(t(do.call(cbind, l)))
  names(mdy_df) <- c("m", "d", "y")
  mdy_df$m <- as.numeric(levels(mdy_df$m))[mdy_df$m]    #values were converted to factors
  mdy_df$d <- as.numeric(levels(mdy_df$d))[mdy_df$d]
  mdy_df$y <- as.numeric(levels(mdy_df$y))[mdy_df$y]
  jday <- date::mdy.date(mdy_df$m,
                         mdy_df$d,
                         mdy_df$y)
  jday <- as.numeric(jday)
  return(jday)
}
