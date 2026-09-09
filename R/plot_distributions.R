#' Plot the distribution of each clock's estimates
#'
#' A violin plus box of every clock's per-sample values, to see the spread and
#' any skew or outliers at a glance. Because clocks live on different scales
#' (years, counts, scores), \code{scale = TRUE} z-scores each clock to
#' one axis; \code{scale = FALSE} facets them on their own scales.
#'
#' @param x A \code{methylclock} result, or a data frame with an \code{id}
#'   column and one numeric column per clock.
#' @param clocks Optional character vector of clock names. Default: all numeric
#'   clocks in \code{x}.
#' @param scale z-score each clock onto one axis (\code{TRUE}) or facet each
#'   on its own scale (\code{FALSE}, the default).
#' @param ncol Facet columns when \code{scale = FALSE}. Default 4.
#' @param color_by Colour the violins by clock family (its registry
#'   \code{target}): \code{"auto"} (the default) colours only when more than one
#'   family is shown, \code{"target"} always, \code{"none"} never.
#' @return A \code{ggplot} object.
#' @examples
#' data(methylclock_demo)
#' plotClockDistributions(methylclock_demo,
#'                        clocks = c("Horvath", "Hannum", "Levine"))
#' @export
plotClockDistributions <- function(x, clocks = NULL, scale = FALSE, ncol = 4,
                                   color_by = c("auto", "target", "none")) {
    color_by <- match.arg(color_by)
    long <- .mc_plot_long(x, clocks)
    if (isTRUE(scale)) {
        long$value <- stats::ave(long$value, long$clock, FUN = function(z) {
            s <- stats::sd(z, na.rm = TRUE)
            if (is.na(s) || s == 0) z - mean(z, na.rm = TRUE)
            else (z - mean(z, na.rm = TRUE)) / s
        })
    }
    fam <- .mc_clock_family(levels(long$clock))
    do_col <- color_by == "target" ||
        (color_by == "auto" && length(unique(fam)) > 1L)

    if (do_col) {
        long$family <- factor(unname(fam[as.character(long$clock)]))
        p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$clock,
                y = .data$value, fill = .data$family)) +
            ggplot2::geom_violin(alpha = 0.4, colour = "grey70",
                                 na.rm = TRUE) +
            ggplot2::geom_boxplot(width = 0.16, outlier.size = 0.6,
                                  outlier.colour = "grey35", alpha = 0.9,
                                  na.rm = TRUE) +
            .mc_fill_qual("Clock type")
    } else {
        p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$clock,
                y = .data$value)) +
            ggplot2::geom_violin(fill = .mc_accent_soft, colour = "grey60",
                                 na.rm = TRUE) +
            ggplot2::geom_boxplot(width = 0.15, outlier.size = 0.6,
                                  fill = .mc_accent_mid, outlier.colour =
                                      "grey35", na.rm = TRUE)
    }
    if (isTRUE(scale)) {
        p <- p + ggplot2::labs(x = NULL, y = "Estimate (z-scored per clock)") +
            theme_methylclock() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                hjust = 1))
    } else {
        # the facet strip already names each clock, so the x label is redundant
        p <- p +
            ggplot2::facet_wrap(~ .data$clock, ncol = ncol, scales = "free") +
            ggplot2::labs(x = NULL, y = "Estimate") +
            theme_methylclock() +
            ggplot2::theme(axis.text.x = ggplot2::element_blank(),
                           axis.ticks.x = ggplot2::element_blank())
    }
    p
}
