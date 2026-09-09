# Per-sample quality control: which samples stand out, and how much of their
# data had to be imputed.

#' Per-sample quality control
#'
#' Flags samples worth a second look: those whose clock estimates are jointly
#' extreme in clock space (a robust sign of a technical outlier), and --- if you
#' pass the methylation matrix --- those missing a large fraction of the clock
#' CpGs (values that had to be imputed). A sample high on either count carries
#' estimates that are closer to an approximation than a measurement.
#'
#' @param x A \code{methylclock} result, or a data frame with an \code{id}
#'   column and one numeric column per clock.
#' @param betas Optional methylation matrix (or accepted container): adds the
#'   per-sample percentage of clock CpGs that were missing (and so imputed).
#' @param clocks Optional clock names to restrict to.
#' @param outlier_sd A sample is flagged an outlier when its mean absolute
#'   z-score across clocks exceeds this. Default 3.
#' @return A data frame with one row per sample: \code{id}, \code{outlier_score}
#'   (mean absolute z-score across clocks), \code{is_outlier},
#'   \code{discordance} (standard deviation of the z-scores across clocks ---
#'   how much the clocks contradict each other on that sample; drawn by
#'   \code{\link{plotSampleDiscordance}}), and --- with \code{betas} ---
#'   \code{pct_missing}.
#' @examples
#' data(methylclock_demo)
#' head(sampleQC(methylclock_demo))
#' @export
sampleQC <- function(x, betas = NULL, clocks = NULL, outlier_sd = 3) {
    v <- .mc_plot_values(x, clocks)
    if (is.null(clocks)) {
        keep <- names(v$values)[.mc_clock_family(names(v$values)) != "other"]
        if (length(keep)) v <- .mc_plot_values(x, keep)
    }
    z <- scale(as.matrix(v$values))
    score <- rowMeans(abs(z), na.rm = TRUE)
    disc <- apply(z, 1, stats::sd, na.rm = TRUE)
    out <- data.frame(id = v$ids, outlier_score = round(score, 2),
                      is_outlier = is.finite(score) & score > outlier_sd,
                      discordance = round(disc, 2),
                      stringsAsFactors = FALSE)
    if (!is.null(betas)) {
        b <- .extract_betas(betas)
        rn <- rownames(b)
        entries <- if (is.null(clocks)) clock_registry()
                   else lapply(clocks, clock_info)
        mk <- unique(unlist(lapply(entries, .clock_cpgs)))
        present <- intersect(mk, rn)
        if (length(present)) {
            sub <- as.matrix(b[present, , drop = FALSE])
            pct <- 100 * colMeans(is.na(sub))
            out$pct_missing <- round(unname(pct[out$id]), 1)
        }
    }
    out
}
