#' methylclock: DNAm age estimation with a declarative clock registry
#'
#' methylclock estimates chronological, gestational and biological DNA
#' methylation (DNAm) age from a panel of methylation-based clocks. Clocks are
#' described declaratively in a registry (one entry per clock, see
#' \code{\link{clock_register}}), and their coefficients are resolved through a
#' single cache-backed resolver (\code{\link{mcd_resource}}). Extending the
#' package with a new clock is therefore a self-contained addition rather than a
#' change spread across the estimation code.
#'
#' Several predictor families are supported: linear, principal-component,
#' surrogate, counter, iterative, neural-network and classifier clocks (see
#' \code{\link{clock_predictors}}).
#'
#' @keywords internal
#' @useDynLib methylclock, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom rlang .data
"_PACKAGE"
