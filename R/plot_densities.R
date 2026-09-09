#' Plot the density of each clock's estimates, optionally split by a group
#'
#' A smoothed density curve of the per-sample values of every clock. With
#' \code{group}, the curves of the groups are drawn together in each panel, so
#' a shift, a widening or an extra mode between groups shows up as a visible
#' change of shape.
#'
#' Three layouts are available through \code{style}. \code{"overlay"} (the
#' default) superimposes the group curves in one panel per clock.
#' \code{"mirror"} needs exactly two groups and draws one upwards and the
#' other downwards from a shared baseline, so any asymmetry between the two
#' shapes is the figure itself. \code{"ridge"} stacks all the clocks in a
#' single column, one ridge per clock, giving a one-look overview across
#' clocks; because that puts every clock on the same axis, the estimates are
#' z-scored per clock there (accelerations are already in shared units).
#'
#' With \code{what = "acceleration"} the curves are drawn over the age
#' acceleration (the residual of predicted age regressed on chronological age)
#' instead of the raw estimates. The distinction matters when the groups differ
#' in chronological age: the raw estimates then differ with it, so their
#' densities separate for that reason alone, while the acceleration removes
#' that trend and leaves the scales on which the group shapes are comparable.
#'
#' A density smoothed from a handful of points takes its shape mostly from the
#' smoothing, not from the data, so any clock/group with fewer than
#' \code{min_n} values is drawn as ticks along its baseline rather than as a
#' curve.
#'
#' @param x A \code{methylclock} result, or a data frame with an \code{id}
#'   column and one numeric column per clock.
#' @param group Optional categorical vector (factor or character), one value
#'   per sample in \code{x}'s order, or named by sample id. Samples with a
#'   missing group are dropped.
#' @param clocks Optional character vector of clock names. Default: all numeric
#'   clocks in \code{x} (for \code{what = "acceleration"}, the age clocks).
#' @param what Plot the raw \code{"estimate"} (the default) or the
#'   \code{"acceleration"}, which requires \code{age}.
#' @param style \code{"overlay"} (the default), \code{"mirror"} (two groups,
#'   one drawn downwards) or \code{"ridge"} (all clocks in one column).
#' @param age Numeric vector of chronological ages, one per sample in
#'   \code{x}'s order (or named by sample id). Only used --- and required ---
#'   when \code{what = "acceleration"}.
#' @param ncol Facet columns for \code{"overlay"} and \code{"mirror"}.
#'   Default 4.
#' @param min_n Minimum number of values a clock/group needs for a curve;
#'   below it, the values are drawn as ticks. Default 5.
#' @return A \code{ggplot} object.
#' @examples
#' data(methylclock_demo)
#' d <- methylclock_demo
#' cl <- c("Horvath", "Hannum", "Levine")
#' plotClockDensities(d, group = d$sex, clocks = cl)
#' plotClockDensities(d, group = d$sex, clocks = cl, what = "acceleration",
#'                    age = d$age, style = "mirror")
#' plotClockDensities(d, group = d$sex, clocks = cl, style = "ridge")
#' @export
plotClockDensities <- function(x, group = NULL, clocks = NULL,
                               what = c("estimate", "acceleration"),
                               style = c("overlay", "mirror", "ridge"),
                               age = NULL, ncol = 4, min_n = 5L) {
    what <- match.arg(what)
    style <- match.arg(style)
    min_n <- max(2L, as.integer(min_n))
    v <- .mc_plot_values(x, clocks)

    if (what == "acceleration") {
        if (is.null(age))
            stop("`what = \"acceleration\"` needs `age`.")
        if (is.null(clocks))
            v <- .mc_plot_values(x, .mc_age_clocks(names(v$values)))
        age <- .mc_align_age(age, v$ids)
        rownames(v$values) <- v$ids
        vals <- as.data.frame(.mc_age_acceleration(v$values, age))
    } else {
        vals <- v$values
    }

    long <- data.frame(
        id = rep(v$ids, times = ncol(vals)),
        clock = factor(rep(names(vals), each = length(v$ids)),
                       levels = names(vals)),
        value = as.numeric(unlist(vals, use.names = FALSE)),
        stringsAsFactors = FALSE)
    if (!is.null(group)) {
        group <- .mc_align_group(group, v$ids)
        long$group <- .mc_group_factor(
            rep(as.character(group), times = ncol(vals)), group)
        long <- long[!is.na(long$group), , drop = FALSE]
    }
    long <- long[is.finite(long$value), , drop = FALSE]
    if (!nrow(long))
        stop("No finite values to plot.")
    if (style == "mirror" &&
        (is.null(long$group) || nlevels(droplevels(long$group)) != 2L))
        stop("`style = \"mirror\"` needs a `group` with exactly two levels.")

    if (style == "ridge" && what == "estimate") {
        # one shared axis, so each clock is put on it as z-scores
        long$value <- stats::ave(long$value, long$clock, FUN = function(z) {
            s <- stats::sd(z)
            if (is.na(s) || s == 0) z - mean(z) else (z - mean(z)) / s
        })
    }
    xlab <- if (what == "acceleration") "Age acceleration (years)"
        else if (style == "ridge") "Estimate (z-scored per clock)"
        else "Estimate"

    p <- switch(style,
                overlay = .mc_pcd_overlay(long, min_n, ncol),
                mirror = .mc_pcd_mirror(long, min_n, ncol),
                ridge = .mc_pcd_ridge(long, min_n))
    if (what == "acceleration")
        p <- p + ggplot2::geom_vline(xintercept = 0, linetype = "dashed",
                                     colour = "grey55")
    p <- p + ggplot2::labs(x = xlab) + theme_methylclock() +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.line.x = ggplot2::element_line(colour = "grey80"))
    if (!is.null(long$group))
        p <- p + ggplot2::theme(legend.position = "bottom")
    p
}

# Density curves per clock (and group) as ready-to-draw x/y rows; combos with
# fewer than min_n values are returned apart, to show as ticks instead.
.mc_pcd_curves <- function(long, min_n) {
    grouped <- !is.null(long$group)
    f <- if (grouped) interaction(long$clock, long$group, drop = TRUE)
        else factor(long$clock)
    pieces <- split(long, f)
    small <- vapply(pieces, nrow, integer(1)) < min_n
    curves <- lapply(pieces[!small], function(d) {
        den <- stats::density(d$value)
        out <- data.frame(clock = d$clock[1], x = den$x, y = den$y)
        if (grouped) out$group <- d$group[1]
        out
    })
    curves <- if (length(curves)) do.call(rbind, curves) else NULL
    ticks <- do.call(rbind, c(list(long[0, , drop = FALSE]), pieces[small]))
    if (nrow(ticks))
        message("Clock/group combinations with fewer than ", min_n,
                " values are drawn as ticks, not density curves.")
    list(curves = curves, ticks = ticks)
}

# One panel per clock, group curves superimposed.
.mc_pcd_overlay <- function(long, min_n, ncol) {
    cd <- .mc_pcd_curves(long, min_n)
    grouped <- !is.null(long$group)
    p <- ggplot2::ggplot()
    if (!is.null(cd$curves)) {
        if (grouped) {
            p <- p +
                ggplot2::geom_ribbon(data = cd$curves,
                    ggplot2::aes(x = .data$x, ymin = 0, ymax = .data$y,
                                 fill = .data$group),
                    alpha = 0.35, colour = NA) +
                ggplot2::geom_line(data = cd$curves,
                    ggplot2::aes(x = .data$x, y = .data$y,
                                 colour = .data$group),
                    linewidth = 0.5) +
                .mc_fill_group() + .mc_colour_group()
        } else {
            p <- p +
                ggplot2::geom_ribbon(data = cd$curves,
                    ggplot2::aes(x = .data$x, ymin = 0, ymax = .data$y),
                    fill = .mc_accent_soft, colour = NA) +
                ggplot2::geom_line(data = cd$curves,
                    ggplot2::aes(x = .data$x, y = .data$y),
                    colour = .mc_accent, linewidth = 0.5)
        }
    }
    if (nrow(cd$ticks)) {
        p <- if (grouped)
            p + ggplot2::geom_rug(data = cd$ticks,
                ggplot2::aes(x = .data$value, colour = .data$group),
                sides = "b", na.rm = TRUE)
        else
            p + ggplot2::geom_rug(data = cd$ticks,
                ggplot2::aes(x = .data$value), colour = .mc_accent,
                sides = "b", na.rm = TRUE)
    }
    p + ggplot2::facet_wrap(~ .data$clock, ncol = ncol, scales = "free") +
        ggplot2::labs(y = "Density")
}

# One panel per clock, the first group upwards and the second downwards.
.mc_pcd_mirror <- function(long, min_n, ncol) {
    cd <- .mc_pcd_curves(long, min_n)
    lev <- levels(droplevels(long$group))
    cur <- cd$curves
    p <- ggplot2::ggplot()
    if (!is.null(cur)) {
        cur$edge <- ifelse(cur$group == lev[1], cur$y, -cur$y)
        p <- p +
            ggplot2::geom_ribbon(data = cur,
                ggplot2::aes(x = .data$x, ymin = pmin(0, .data$edge),
                             ymax = pmax(0, .data$edge), fill = .data$group),
                alpha = 0.55, colour = NA) +
            ggplot2::geom_line(data = cur,
                ggplot2::aes(x = .data$x, y = .data$edge,
                             colour = .data$group), linewidth = 0.5) +
            .mc_fill_group() + .mc_colour_group()
    }
    if (nrow(cd$ticks))
        p <- p + ggplot2::geom_rug(data = cd$ticks,
            ggplot2::aes(x = .data$value, colour = .data$group),
            sides = "b", na.rm = TRUE)
    p + ggplot2::geom_hline(yintercept = 0, colour = "grey40",
                            linewidth = 0.3) +
        ggplot2::facet_wrap(~ .data$clock, ncol = ncol, scales = "free") +
        ggplot2::scale_y_continuous(
            labels = function(b) format(abs(b), trim = TRUE)) +
        ggplot2::labs(y = paste0("Density (", lev[1], " up, ",
                                 lev[2], " down)"))
}

# All clocks in one column, one ridge per clock (groups overlaid within it).
.mc_pcd_ridge <- function(long, min_n) {
    cd <- .mc_pcd_curves(long, min_n)
    grouped <- !is.null(long$group)
    lev <- levels(droplevels(long$clock))
    base <- stats::setNames(rev(seq_along(lev)), lev)   # first clock on top
    cur <- cd$curves
    if (!is.null(cur)) {
        m <- stats::ave(cur$y, cur$clock, FUN = max)
        cur$ymin <- base[as.character(cur$clock)]
        cur$ymax <- cur$ymin + cur$y / m * 0.85
    }
    p <- ggplot2::ggplot() +
        ggplot2::geom_hline(yintercept = unname(base), colour = "grey88",
                            linewidth = 0.3)
    if (!is.null(cur)) {
        if (grouped) {
            p <- p +
                ggplot2::geom_ribbon(data = cur,
                    ggplot2::aes(x = .data$x, ymin = .data$ymin,
                                 ymax = .data$ymax, fill = .data$group,
                                 group = interaction(.data$clock,
                                                     .data$group)),
                    alpha = 0.35, colour = NA) +
                ggplot2::geom_line(data = cur,
                    ggplot2::aes(x = .data$x, y = .data$ymax,
                                 colour = .data$group,
                                 group = interaction(.data$clock,
                                                     .data$group)),
                    linewidth = 0.45) +
                .mc_fill_group() + .mc_colour_group()
        } else {
            p <- p +
                ggplot2::geom_ribbon(data = cur,
                    ggplot2::aes(x = .data$x, ymin = .data$ymin,
                                 ymax = .data$ymax, group = .data$clock),
                    fill = .mc_accent_soft, alpha = 0.8, colour = NA) +
                ggplot2::geom_line(data = cur,
                    ggplot2::aes(x = .data$x, y = .data$ymax,
                                 group = .data$clock),
                    colour = .mc_accent, linewidth = 0.45)
        }
    }
    if (nrow(cd$ticks)) {
        cd$ticks$y0 <- base[as.character(cd$ticks$clock)]
        p <- if (grouped)
            p + ggplot2::geom_segment(data = cd$ticks,
                ggplot2::aes(x = .data$value, xend = .data$value,
                             y = .data$y0, yend = .data$y0 + 0.15,
                             colour = .data$group))
        else
            p + ggplot2::geom_segment(data = cd$ticks,
                ggplot2::aes(x = .data$value, xend = .data$value,
                             y = .data$y0, yend = .data$y0 + 0.15),
                colour = .mc_accent)
    }
    p + ggplot2::scale_y_continuous(breaks = unname(base),
                                    labels = names(base)) +
        ggplot2::labs(y = NULL)
}
