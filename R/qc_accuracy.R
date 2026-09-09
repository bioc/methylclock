#' How well each clock tracks chronological age
#'
#' A quick quality-control summary for a dataset with known ages: for every age
#' clock, how closely its estimate follows chronological age. It reports the
#' number of samples, the Pearson correlation, the R-squared and the median
#' absolute error (in years). Pass \code{by} to compute the same summary within
#' groups --- sex, case/control, cohort, ethnicity --- which is the usual way to
#' check that a clock behaves across subgroups rather than only on average.
#'
#' A high correlation with a large error means a clock *follows* age but is
#' *miscalibrated* for your data (a constant offset), which the two columns
#' surface separately.
#'
#' @param x A \code{methylclock} result, or a data frame with an \code{id}
#'   column and one numeric column per clock.
#' @param age Numeric vector of chronological ages, one per sample in \code{x}'s
#'   order (or named by sample id).
#' @param by Optional categorical vector (one value per sample, or named by id):
#'   the summary is computed separately within each group.
#' @param clocks Optional clock names to restrict to. Default: the clocks that
#'   estimate an age in years (a comparison with chronological age is only
#'   meaningful for those).
#' @return A data frame with one row per clock (per group), and columns
#'   \code{clock}, optionally \code{group}, \code{n}, \code{r}, \code{R2} and
#'   \code{MAE}.
#' @examples
#' data(methylclock_demo)
#' clockAccuracy(methylclock_demo, age = methylclock_demo$age,
#'               by = methylclock_demo$sex)
#' @export
clockAccuracy <- function(x, age, by = NULL, clocks = NULL) {
    v <- .mc_plot_values(x, clocks)
    if (is.null(clocks))
        v <- .mc_plot_values(x, .mc_age_clocks(names(v$values)))
    age <- .mc_align_age(age, v$ids)
    grp <- if (is.null(by)) rep("all", length(v$ids))
           else as.character(.mc_align_group(by, v$ids))

    summarise <- function(pred, a) {
        ok <- is.finite(pred) & is.finite(a)
        if (sum(ok) < 3L)
            return(list(n = sum(ok), r = NA_real_, R2 = NA_real_,
                        MAE = NA_real_))
        r <- stats::cor(pred[ok], a[ok])
        list(n = sum(ok), r = r, R2 = r^2,
             MAE = stats::median(abs(pred[ok] - a[ok])))
    }

    out <- list()
    for (g in unique(grp)) {
        sel <- grp == g
        for (nm in names(v$values)) {
            s <- summarise(v$values[[nm]][sel], age[sel])
            row <- data.frame(clock = nm, n = s$n, r = round(s$r, 3),
                              R2 = round(s$R2, 3), MAE = round(s$MAE, 2),
                              stringsAsFactors = FALSE)
            if (!is.null(by)) row <- cbind(group = g, row,
                                           stringsAsFactors = FALSE)
            out[[length(out) + 1L]] <- row
        }
    }
    do.call(rbind, out)
}
