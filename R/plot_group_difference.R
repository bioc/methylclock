#' Plot the per-clock difference between two groups as a forest plot
#'
#' One row per clock, in one of two layouts. The default forest draws the
#' difference between the two groups' means with its confidence interval, all
#' on one shared axis --- a one-look answer to which clocks separate the
#' groups and by how much. The dashed line at zero means "no difference": an
#' interval that reaches it is compatible with none, and the points of
#' intervals that leave zero out are drawn filled. \code{style = "dumbbell"}
#' draws instead the two group means joined by a segment, so both levels are
#' visible rather than only their difference --- one group being high tells a
#' different story than the other being low, even when the gap is the same.
#'
#' By default the comparison is on the age acceleration (the residual of
#' predicted age regressed on chronological age), which requires \code{age}.
#' When the groups differ in chronological age, the raw estimates differ with
#' it, so the acceleration is the scale on which a group difference means
#' something beyond age. When no ages are known, \code{what = "estimate"}
#' compares the raw estimates instead; any difference in chronological age
#' between the groups is then part of the number.
#'
#' @param x A \code{methylclock} result, or a data frame with an \code{id}
#'   column and one numeric column per clock.
#' @param group Categorical vector (factor or character) with exactly two
#'   levels, one value per sample in \code{x}'s order, or named by sample id.
#'   Samples with a missing group are dropped. The difference is second level
#'   minus first, as the axis label states.
#' @param clocks Optional character vector of clock names. Default: the age
#'   clocks present in \code{x}.
#' @param what Compare the \code{"acceleration"} (the default, requires
#'   \code{age}) or the raw \code{"estimate"}.
#' @param style \code{"forest"} (the default: the difference with its
#'   confidence interval) or \code{"dumbbell"} (the two group means joined by
#'   a segment).
#' @param age Numeric vector of chronological ages, one per sample in
#'   \code{x}'s order (or named by sample id). Only used --- and required ---
#'   when \code{what = "acceleration"}.
#' @param conf Confidence level of the intervals (Welch t). Default 0.95.
#' @return A \code{ggplot} object.
#' @examples
#' data(methylclock_demo)
#' d <- methylclock_demo
#' cl <- c("Horvath", "Hannum", "Levine")
#' plotGroupDifference(d, group = d$sex, age = d$age, clocks = cl)
#' plotGroupDifference(d, group = d$sex, what = "estimate", clocks = cl,
#'                     style = "dumbbell")
#' @export
plotGroupDifference <- function(x, group, clocks = NULL,
                                what = c("acceleration", "estimate"),
                                style = c("forest", "dumbbell"),
                                age = NULL, conf = 0.95) {
    what <- match.arg(what)
    style <- match.arg(style)
    v <- .mc_plot_values(x, clocks)
    if (is.null(clocks))
        v <- .mc_plot_values(x, .mc_age_clocks(names(v$values)))
    group <- .mc_align_group(group, v$ids)
    g <- factor(as.character(group))
    keep <- !is.na(g)
    g <- droplevels(g[keep])
    if (nlevels(g) != 2L)
        stop("`group` must have exactly two levels; got ", nlevels(g), ".")
    vals <- v$values[keep, , drop = FALSE]
    if (what == "acceleration") {
        if (is.null(age))
            stop("`what = \"acceleration\"` needs `age`; without ages, ",
                 "`what = \"estimate\"` compares the raw estimates.")
        age <- .mc_align_age(age, v$ids)[keep]
        rownames(vals) <- v$ids[keep]
        vals <- as.data.frame(.mc_age_acceleration(vals, age))
    }

    lev <- levels(g)
    rows <- lapply(names(vals), function(nm) {
        a <- vals[[nm]][g == lev[1]]
        b <- vals[[nm]][g == lev[2]]
        a <- a[is.finite(a)]
        b <- b[is.finite(b)]
        if (length(a) < 2L || length(b) < 2L)
            return(NULL)
        tt <- stats::t.test(b, a, conf.level = conf)
        data.frame(clock = nm, diff = mean(b) - mean(a),
                   mean_a = mean(a), mean_b = mean(b),
                   lo = tt$conf.int[1], hi = tt$conf.int[2],
                   stringsAsFactors = FALSE)
    })
    skipped <- names(vals)[vapply(rows, is.null, logical(1))]
    fb <- do.call(rbind, rows)
    if (is.null(fb))
        stop("No clock has at least two finite values in each group.")
    if (length(skipped))
        message("Skipped (fewer than two values in a group): ",
                paste(skipped, collapse = ", "))
    fb$clock <- factor(fb$clock, levels = rev(fb$clock))
    fb$excludes_zero <- fb$lo > 0 | fb$hi < 0

    if (style == "dumbbell") {
        pal <- stats::setNames(.mc_group_pal[seq_len(2L)], lev)
        xlab <- if (what == "acceleration")
            "Mean age acceleration by group (years)"
        else "Mean estimate by group"
        return(ggplot2::ggplot(fb) +
            ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                                colour = "grey55") +
            ggplot2::geom_segment(
                ggplot2::aes(x = .data$mean_a, xend = .data$mean_b,
                             y = .data$clock, yend = .data$clock),
                colour = "grey70", linewidth = 0.8) +
            ggplot2::geom_point(
                ggplot2::aes(x = .data$mean_a, y = .data$clock,
                             colour = lev[1]), size = 2.6) +
            ggplot2::geom_point(
                ggplot2::aes(x = .data$mean_b, y = .data$clock,
                             colour = lev[2]), size = 2.6) +
            ggplot2::scale_colour_manual(values = pal, name = NULL) +
            ggplot2::labs(x = xlab, y = NULL) +
            theme_methylclock() +
            ggplot2::theme(
                panel.grid = ggplot2::element_blank(),
                axis.line.x = ggplot2::element_line(colour = "grey80"),
                legend.position = "bottom"))
    }

    pct <- round(conf * 100)
    xlab <- if (what == "acceleration")
        sprintf("Difference in age acceleration, %s - %s (years, %d%% CI)",
                lev[2], lev[1], pct)
    else
        sprintf("Difference in estimate, %s - %s (%d%% CI)", lev[2], lev[1],
                pct)
    ggplot2::ggplot(fb, ggplot2::aes(x = .data$diff, y = .data$clock)) +
        ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                            colour = "grey55") +
        ggplot2::geom_errorbarh(ggplot2::aes(xmin = .data$lo,
                                             xmax = .data$hi),
                                height = 0.18, colour = .mc_accent,
                                linewidth = 0.6) +
        ggplot2::geom_point(ggplot2::aes(shape = .data$excludes_zero),
                            size = 2.4, colour = .mc_accent, fill = "white") +
        ggplot2::scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 21),
                                    guide = "none") +
        ggplot2::labs(x = xlab, y = NULL) +
        theme_methylclock() +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_line(colour = "grey80"))
}
