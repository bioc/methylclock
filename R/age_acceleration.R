# Epigenetic age acceleration measures, including the cell-composition-adjusted
# (intrinsic) one. Follows the definitions the original package used.

# Residual of y on the columns of df (NA where unfittable), full length.
.mc_resid_on <- function(y, df) {
    r <- rep(NA_real_, length(y))
    ok <- is.finite(y) & stats::complete.cases(df)
    if (sum(ok) >= 3L)
        r[ok] <- stats::residuals(
            stats::lm(y[ok] ~ ., data = df[ok, , drop = FALSE]))
    r
}

#' Epigenetic age acceleration, intrinsic and extrinsic
#'
#' For each age clock, how much older (or younger) a sample looks than its
#' chronological age. Three measures per clock:
#' \describe{
#'   \item{\code{ageAcc}}{The raw difference, predicted minus chronological age.}
#'   \item{\code{residual}}{The classic age-acceleration residual: predicted
#'     regressed on chronological age. (The published EEAA of Chen and
#'     colleagues is a different, cell-weighted measure; \code{\link{EEAA}}
#'     returns it under its usual name.)}
#'   \item{\code{residualCells}}{The residual further adjusted for blood cell
#'     composition (whatever panel you supply), so it reflects ageing
#'     independent of shifts in cell proportions. Only returned when
#'     \code{cell_counts} is supplied. This is the same adjustment idea as
#'     the canonical IEAA of the literature, but that measure is defined
#'     with a specific set of seven immune covariates ---
#'     \code{\link{IEAA}} computes it as defined.}
#' }
#' Estimate the cell proportions with \code{\link{cellCounts}}.
#'
#' @param x A \code{methylclock} result, or a data frame with an \code{id}
#'   column and one numeric column per clock.
#' @param age Numeric vector of chronological ages, one per sample in \code{x}'s
#'   order (or named by sample id).
#' @param cell_counts Optional samples-by-cell-types matrix (from
#'   \code{\link{cellCounts}}); adds the cell-adjusted residual
#'   \code{residualCells}.
#' @param clocks Optional clock names. Default: the age clocks in \code{x}.
#' @return A data frame with columns \code{id}, \code{clock}, \code{ageAcc},
#'   \code{residual} and (with \code{cell_counts}) \code{residualCells}.
#' @seealso \code{\link{cellCounts}}, \code{\link{plotAgeAcceleration}}
#' @examples
#' data(methylclock_demo)
#' head(ageAcceleration(methylclock_demo, age = methylclock_demo$age))
#' @export
ageAcceleration <- function(x, age, cell_counts = NULL, clocks = NULL) {
    v <- .mc_plot_values(x, clocks)
    if (is.null(clocks))
        v <- .mc_plot_values(x, .mc_age_clocks(names(v$values)))
    age <- .mc_align_age(age, v$ids)

    cc <- NULL
    if (!is.null(cell_counts)) {
        cc <- as.matrix(cell_counts)
        if (!is.null(rownames(cc)) && all(v$ids %in% rownames(cc)))
            cc <- cc[v$ids, , drop = FALSE]
        else if (nrow(cc) != length(v$ids))
            stop("`cell_counts` must have one row per sample (or be named by ",
                 "sample id).", call. = FALSE)
        # drop invariant cell types, as the original did (IQR > 1e-6)
        keep <- apply(cc, 2, function(z) stats::IQR(z, na.rm = TRUE) > 1e-6)
        cc <- cc[, keep, drop = FALSE]
    }

    rows <- lapply(names(v$values), function(nm) {
        y <- v$values[[nm]]
        out <- data.frame(id = v$ids, clock = nm, ageAcc = y - age,
                          residual = .mc_resid_on(y, data.frame(age = age)),
                          stringsAsFactors = FALSE)
        if (!is.null(cc))
            out$residualCells <- .mc_resid_on(y, data.frame(age = age, cc,
                                                            check.names = FALSE))
        out
    })
    do.call(rbind, rows)
}
