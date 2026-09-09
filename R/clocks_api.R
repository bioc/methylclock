# User-facing estimation API: one engine, three entry points. `methylclock` is
# the unified entry; `DNAmAge` and `DNAmGA` are thin wrappers kept for
# continuity (they are not deprecated).
# INTERNAL: unify inside, don't break outside = D8.

#' Estimate DNA methylation clocks
#'
#' Unified entry point: estimates any registered clocks from a methylation
#' matrix. \code{\link{DNAmAge}} and \code{\link{DNAmGA}} are convenience
#' wrappers that restrict the selection to non-gestational and gestational
#' clocks respectively.
#'
#' @param x A methylation matrix with CpG identifiers as row names and samples
#'   as columns (a data frame with a CpG-identifier column is also accepted), or
#'   a large matrix held on disk as an HDF5-backed matrix.
#' @param clocks Either \code{"all"} (default) or a character vector of clock
#'   names; see \code{\link{clock_list}}.
#' @param generation Optional generation to restrict to (e.g. \code{1}).
#' @param platform Optional native platform to restrict to (e.g. \code{"450K"}),
#'   handy when the data come from a known array and only the clocks built for
#'   it should run.
#' @param storage How the estimates are held: \code{"auto"} (default) backs them
#'   on disk (HDF5) with an in-memory cache when they fit, \code{"hdf5"} on disk
#'   only, \code{"memory"} in memory only (not persisted). See
#'   \code{\link{compute_clocks}} for the \code{file} and \code{compression}
#'   arguments, and \code{\link{persist}} to keep a session result.
#' @param ... Passed to \code{\link{compute_clocks}} (e.g. \code{file},
#'   \code{compression}).
#' @return A \code{"methylclock"} object.
#' @seealso \code{\link{DNAmAge}}, \code{\link{DNAmGA}}, \code{\link{persist}},
#'   \code{\link{clock_list}}
#' @examples
#' data(methylclock_betas)
#' res <- methylclock(methylclock_betas, clocks = c("Horvath", "Levine"))
#' head(as.data.frame(res))
#' @param min.perc Minimum fraction of a clock's CpGs that must be present for
#'   it to be computed (otherwise \code{NA}). Default \code{0.8}.
#' @param impute How missing values among the present CpGs are filled:
#'   \code{"mean"} (default, each CpG's across-sample mean), \code{"reference"}
#'   (the CpG's reference level from the clock's own coefficients, mean
#'   fallback), \code{"none"} (a sample missing any used CpG is returned
#'   \code{NA}), or \code{"knn"} (each from the CpG's nearest co-varying CpGs;
#'   see \code{\link{imputeKNN}}). See \code{\link{compute_clocks}}.
#' @param min.perc.sample Per-sample coverage floor: a sample with fewer than
#'   this fraction of a clock's present CpGs observed is returned \code{NA} for
#'   that clock. Default \code{0} (off).
#' @export
methylclock <- function(x, clocks = "all", generation = NULL, platform = NULL,
                        min.perc = 0.8,
                        impute = c("mean", "reference", "none", "knn"),
                        min.perc.sample = 0, storage = "auto", ...) {
    .alloc_guard(compute_clocks(x, clocks = clocks, target = NULL,
                                generation = generation, platform = platform,
                                min.perc = min.perc, impute = match.arg(impute),
                                min.perc.sample = min.perc.sample,
                                storage = storage, ...))
}

#' Estimate chronological and biological DNAm age
#'
#' Convenience wrapper around \code{\link{methylclock}} restricted to
#' non-gestational clocks.
#'
#' @inheritParams methylclock
#' @return A \code{"methylclock"} object.
#' @seealso \code{\link{methylclock}}, \code{\link{DNAmGA}}
#' @examples
#' data(methylclock_betas)
#' head(as.data.frame(DNAmAge(methylclock_betas, clocks = "Horvath")))
#' @export
DNAmAge <- function(x, clocks = "all", generation = NULL, platform = NULL,
                    min.perc = 0.8, storage = "auto", ...) {
    target <- if (identical(clocks, "all")) "!gestational" else NULL
    .alloc_guard(compute_clocks(x, clocks = clocks, target = target,
                                generation = generation, platform = platform,
                                min.perc = min.perc, storage = storage, ...))
}

#' Estimate gestational DNAm age
#'
#' Convenience wrapper around \code{\link{methylclock}} restricted to
#' gestational clocks.
#'
#' @inheritParams methylclock
#' @return A \code{"methylclock"} object.
#' @seealso \code{\link{methylclock}}, \code{\link{DNAmAge}}
#' @export
DNAmGA <- function(x, clocks = "all", generation = NULL, platform = NULL,
                   min.perc = 0.8, storage = "auto", ...) {
    target <- if (identical(clocks, "all")) "gestational" else NULL
    .alloc_guard(compute_clocks(x, clocks = clocks, target = target,
                                generation = generation, platform = platform,
                                min.perc = min.perc, storage = storage, ...))
}
