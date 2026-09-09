# Output transforms: functions applied to a clock's raw linear prediction to put
# it on the reported scale. A clock references one through its `transform` field.

#' Horvath age transformation (inverse)
#'
#' Several clocks are trained on an age variable that is compressed on a log
#' scale below a reference adult age. This maps a raw prediction back to years.
#'
#' @param x Numeric vector of raw predictions.
#' @param adult_age Reference adult age used by the transformation. Default
#'   \code{20}.
#' @return Numeric vector of ages in years.
#' @examples
#' anti_trafo(c(-0.5, 0, 0.5, 1.5))
#' @export
anti_trafo <- function(x, adult_age = 20) {
    ifelse(x < 0,
           (1 + adult_age) * exp(x) - 1,
           (1 + adult_age) * x + adult_age)
}
