#' Heatmap of the correlation between clocks
#'
#' Correlation of the per-sample estimates across every pair of clocks, as a
#' heatmap. Clocks that measure similar things correlate highly; this is the
#' clearest view of how much the clocks agree (and where they are redundant).
#'
#' @param x A \code{methylclock} result, or a data frame with an \code{id}
#'   column and one numeric column per clock.
#' @param clocks Optional character vector of clock names to include. Default:
#'   all numeric clocks in \code{x}.
#' @param method Correlation method passed to \code{\link[stats]{cor}}
#'   (\code{"pearson"}, \code{"spearman"} or \code{"kendall"}). Default Pearson.
#' @param cluster Reorder clocks by hierarchical clustering so correlated ones
#'   sit together. Default \code{TRUE}.
#' @return A \code{ggplot} object.
#' @examples
#' data(methylclock_demo)
#' plotClockCorrelation(methylclock_demo,
#'                      clocks = c("Horvath", "Hannum", "Levine", "TL"))
#' @export
plotClockCorrelation <- function(x, clocks = NULL, method = "pearson",
                                 cluster = TRUE) {
    v <- .mc_plot_values(x, clocks)
    if (ncol(v$values) < 2L)
        stop("Need at least two clocks to show a correlation.")
    cmat <- stats::cor(as.matrix(v$values), use = "pairwise.complete.obs",
                       method = method)
    if (isTRUE(cluster) && nrow(cmat) > 2L) {
        ord <- tryCatch(
            stats::hclust(stats::as.dist(1 - cmat))$order,
            error = function(e) seq_len(nrow(cmat)))
        cmat <- cmat[ord, ord]
    }
    lv <- rownames(cmat)
    long <- data.frame(
        clock_x = factor(rep(lv, times = length(lv)), levels = lv),
        clock_y = factor(rep(lv, each = length(lv)), levels = rev(lv)),
        r = as.numeric(cmat), stringsAsFactors = FALSE)

    ggplot2::ggplot(long, ggplot2::aes(x = .data$clock_x, y = .data$clock_y,
                                       fill = .data$r)) +
        ggplot2::geom_tile(colour = "white") +
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.2f", .data$r)),
                           size = 2.4) +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "grey97",
            high = "#B2182B", midpoint = 0, limits = c(-1, 1), name = "r") +
        ggplot2::labs(x = NULL, y = NULL) +
        ggplot2::coord_fixed() +
        theme_methylclock() +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
}
