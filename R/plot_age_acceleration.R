#' Plot epigenetic age acceleration
#'
#' Age acceleration is the residual of the linear fit of predicted age on
#' chronological age: a positive value is a sample that looks older than its
#' years for it. Two views: \code{"box"} shows the spread per clock;
#' \code{"heatmap"} shows every sample and clock, so a
#' sample accelerated across all clocks (a consistent outlier) stands out.
#'
#' @param x A \code{methylclock} result, or a data frame with an \code{id}
#'   column and one numeric column per clock.
#' @param age Numeric vector of chronological ages, one per sample in \code{x}'s
#'   order (or named by sample id).
#' @param clocks Optional character vector of clock names. Default: the age
#'   clocks present in \code{x}.
#' @param type \code{"box"} (per-clock distribution) or \code{"heatmap"}
#'   (samples x clocks). Default \code{"box"}.
#' @param color_by For \code{type = "box"}, colour the boxes by clock family
#'   (its registry \code{target}): \code{"auto"} (default) only when more than
#'   one family is shown, \code{"target"} always, \code{"none"} never.
#' @return A \code{ggplot} object.
#' @examples
#' data(methylclock_demo)
#' d <- methylclock_demo
#' plotAgeAcceleration(d, age = d$age,
#'                     clocks = c("Horvath", "Hannum", "Levine"))
#' @export
plotAgeAcceleration <- function(x, age, clocks = NULL,
                                type = c("box", "heatmap"),
                                color_by = c("auto", "target", "none")) {
    type <- match.arg(type)
    color_by <- match.arg(color_by)
    v <- .mc_plot_values(x, clocks)
    if (is.null(clocks))
        v <- .mc_plot_values(x, .mc_age_clocks(names(v$values)))
    age <- .mc_align_age(age, v$ids)
    rownames(v$values) <- v$ids
    accel <- .mc_age_acceleration(v$values, age)

    if (type == "box") {
        long <- data.frame(
            clock = factor(rep(colnames(accel), each = nrow(accel)),
                           levels = colnames(accel)),
            accel = as.numeric(accel), stringsAsFactors = FALSE)
        fam <- .mc_clock_family(levels(long$clock))
        do_col <- color_by == "target" ||
            (color_by == "auto" && length(unique(fam)) > 1L)
        if (do_col) {
            long$family <- factor(unname(fam[as.character(long$clock)]))
            # colour the box too (not just the thin violin), so a family reads
            # clearly even when the acceleration spread is narrow
            base <- ggplot2::ggplot(long, ggplot2::aes(x = .data$clock,
                    y = .data$accel, fill = .data$family)) +
                ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                                    colour = "grey55") +
                ggplot2::geom_violin(alpha = 0.35, colour = "grey70",
                                     na.rm = TRUE) +
                ggplot2::geom_boxplot(width = 0.16, outlier.size = 0.6,
                    outlier.colour = "grey35", alpha = 0.9, na.rm = TRUE) +
                .mc_fill_qual("Clock type")
        } else {
            base <- ggplot2::ggplot(long, ggplot2::aes(x = .data$clock,
                    y = .data$accel)) +
                ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                                    colour = "grey55") +
                ggplot2::geom_violin(fill = .mc_accent_soft, colour = "grey60",
                                     na.rm = TRUE) +
                ggplot2::geom_boxplot(width = 0.15, outlier.size = 0.6,
                    fill = .mc_accent_mid, outlier.colour = "grey35",
                    na.rm = TRUE)
        }
        return(base +
            ggplot2::labs(x = NULL, y = "Age acceleration (years)") +
            theme_methylclock() +
            ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
                hjust = 1)))
    }

    # heatmap: z-score each clock so clocks on different scales are comparable
    z <- scale(accel)
    long <- data.frame(
        sample = factor(rep(rownames(accel), times = ncol(accel)),
                        levels = rownames(accel)),
        clock = factor(rep(colnames(accel), each = nrow(accel)),
                       levels = colnames(accel)),
        z = as.numeric(z), stringsAsFactors = FALSE)
    ggplot2::ggplot(long, ggplot2::aes(x = .data$clock, y = .data$sample,
                                       fill = .data$z)) +
        ggplot2::geom_tile() +
        ggplot2::scale_fill_gradient2(low = "#2166AC", mid = "grey97",
            high = "#B2182B", midpoint = 0, name = "accel (z)") +
        ggplot2::labs(x = NULL, y = "Sample") +
        theme_methylclock() +
        ggplot2::theme(panel.grid = ggplot2::element_blank(),
            axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
            axis.text.y = ggplot2::element_blank())
}
