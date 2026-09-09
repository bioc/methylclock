#' Plot each subject's clock trajectory over age
#'
#' The longitudinal view: for each clock, one thin line per subject joining
#' its visits (chronological age on the x axis, the clock's estimate on the
#' y axis), with the identity line dashed --- a subject moving parallel to it
#' ages one epigenetic year per calendar year, and a steeper path runs
#' faster. With \code{group}, lines are coloured by group and a thicker line
#' follows each group's median trajectory across the age range it covers.
#'
#' The lines join the observed visits without smoothing: the wobble of a
#' trajectory is a faithful display of how noisy repeated clock estimates
#' are, which is part of what this plot is for. Per-subject rates with their
#' caveats are computed by \code{\link{trajectoryRates}}.
#'
#' @inheritParams clockTrajectories
#' @param group Optional categorical vector (factor or character), one value
#'   per sample in \code{x}'s order, or named by sample id. Samples with a
#'   missing group are dropped.
#' @param ncol Facet columns. Default 3.
#' @return A \code{ggplot} object.
#' @examples
#' data(methylclock_longitudinal)
#' d <- methylclock_longitudinal
#' plotTrajectories(d, subject = d$subject, age = d$age,
#'                  clocks = c("Wu", "PedBE"))
#' @export
plotTrajectories <- function(x, subject, age, group = NULL, clocks = NULL,
                             ncol = 3) {
    tr <- clockTrajectories(x, subject, age, clocks)
    if (!is.null(group)) {
        v <- .mc_plot_values(x, NULL)
        g <- .mc_align_group(group, v$ids)
        subj_all <- as.character(.mc_align_group(subject, v$ids))
        map <- tapply(as.character(g), subj_all, function(z) z[1])
        tr$group <- .mc_group_factor(unname(map[tr$subject]), group)
        tr <- tr[!is.na(tr$group), , drop = FALSE]
    }

    # the thick guide line: medians within age bins, per clock (and group) —
    # with continuous ages a per-age median would zigzag between single
    # subjects instead of tracing the group's course
    f <- if (is.null(group)) list(tr$clock) else list(tr$clock, tr$group)
    med <- do.call(rbind, lapply(split(tr, f, drop = TRUE), function(d) {
        br <- seq(min(d$age), max(d$age), length.out = 9)
        bin <- cut(d$age, unique(br), include.lowest = TRUE)
        out <- data.frame(
            clock = d$clock[1],
            age = as.numeric(tapply(d$age, bin, stats::median)),
            estimate = as.numeric(tapply(d$estimate, bin, stats::median)))
        if (!is.null(d$group)) out$group <- d$group[1]
        out[is.finite(out$age), , drop = FALSE]
    }))

    p <- ggplot2::ggplot(tr, ggplot2::aes(x = .data$age,
                                          y = .data$estimate)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                             colour = "grey60")
    if (is.null(group)) {
        p <- p +
            ggplot2::geom_line(ggplot2::aes(group = .data$subject),
                               colour = .mc_accent, alpha = 0.35,
                               linewidth = 0.4) +
            ggplot2::geom_point(colour = .mc_accent, alpha = 0.5,
                                size = 0.9) +
            ggplot2::geom_line(data = med, colour = .mc_accent,
                               linewidth = 1.1)
    } else {
        p <- p +
            ggplot2::geom_line(
                ggplot2::aes(group = .data$subject, colour = .data$group),
                alpha = 0.35, linewidth = 0.4) +
            ggplot2::geom_point(ggplot2::aes(colour = .data$group),
                                alpha = 0.5, size = 0.9) +
            ggplot2::geom_line(data = med,
                ggplot2::aes(colour = .data$group, group = .data$group),
                linewidth = 1.1) +
            .mc_colour_group()
    }
    p <- p +
        ggplot2::facet_wrap(~ .data$clock, ncol = ncol, scales = "free") +
        ggplot2::labs(x = "Chronological age (years)",
                      y = "Predicted epigenetic age (years)") +
        theme_methylclock() +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "grey80"))
    if (!is.null(group))
        p <- p + ggplot2::theme(legend.position = "bottom")
    p
}
