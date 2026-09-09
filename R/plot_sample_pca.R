#' PCA of samples in clock space
#'
#' Runs a principal component analysis on the per-sample clock estimates (a
#' samples-by-clocks matrix) and plots the first two components, optionally
#' coloured by a per-sample variable such as age or sex. It shows how samples
#' arrange themselves across the whole panel of clocks: for a healthy age series
#' the leading component typically tracks chronological age.
#'
#' @param x A \code{methylclock} result, or a data frame with an \code{id}
#'   column and one numeric column per clock.
#' @param color Optional per-sample variable to colour points by (numeric for a
#'   continuous scale, factor/character for a discrete one), in \code{x}'s order
#'   or named by sample id.
#' @param color_label Legend title for \code{color}. Default \code{"group"}.
#' @param clocks Optional character vector of clock names to include in the PCA.
#'   Default: all numeric clock columns.
#' @return A \code{ggplot} object.
#' @examples
#' data(methylclock_demo)
#' d <- methylclock_demo
#' plotSamplePCA(d, color = d$age, color_label = "age",
#'               clocks = c("Horvath", "Hannum", "Levine"))
#' @export
plotSamplePCA <- function(x, color = NULL, color_label = "group",
                          clocks = NULL) {
    v <- .mc_plot_values(x, clocks)
    m <- as.matrix(v$values)
    rownames(m) <- v$ids

    # mean-impute stray NAs, drop clocks with no spread (uninformative for PCA)
    for (j in seq_len(ncol(m))) {
        col <- m[, j]
        col[!is.finite(col)] <- mean(col[is.finite(col)])
        m[, j] <- col
    }
    keep <- apply(m, 2, function(col) is.finite(stats::sd(col)) &&
                      stats::sd(col) > 0)
    m <- m[, keep, drop = FALSE]
    if (ncol(m) < 2L)
        stop("Need at least two informative clocks for a PCA.")

    pc <- stats::prcomp(m, center = TRUE, scale. = TRUE)
    pct <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
    df <- data.frame(id = v$ids, PC1 = pc$x[, 1], PC2 = pc$x[, 2],
                     stringsAsFactors = FALSE)

    aes_pt <- ggplot2::aes(x = .data$PC1, y = .data$PC2)
    if (!is.null(color)) {
        df$color <- .mc_align_group(color, v$ids)
        aes_pt <- ggplot2::aes(x = .data$PC1, y = .data$PC2,
                               colour = .data$color)
    }
    p <- ggplot2::ggplot(df, aes_pt) +
        ggplot2::geom_point(size = 2, alpha = 0.85) +
        ggplot2::labs(x = paste0("PC1 (", pct[1], "%)"),
                      y = paste0("PC2 (", pct[2], "%)")) +
        theme_methylclock()
    if (!is.null(color)) {
        p <- p + if (is.numeric(df$color))
            ggplot2::scale_colour_viridis_c(name = color_label)
        else .mc_colour_group(name = color_label)
    }
    p
}
