#' Plot predicted epigenetic age against chronological age
#'
#' For each age clock, a scatter of the predicted age against the samples'
#' chronological age, with the identity line (predicted equal to
#' chronological, drawn dashed), a linear fit and its confidence band, and
#' the fit's \eqn{R^2} in the title. Clocks that do not output an age
#' (mitotic counts, trait scores) are dropped unless named in \code{clocks}.
#'
#' With \code{group}, the points and the fitted lines are drawn per group: a
#' group whose line runs above the other's is predicted older at the same
#' chronological age, and lines with different slopes mean that gap changes
#' along the age range. The outlier highlight is turned off in that case, so
#' its colour cannot be mistaken for a group.
#'
#' @param x A \code{methylclock} result, or a data frame with an \code{id}
#'   column and one numeric column per clock.
#' @param age Numeric vector of chronological ages, one per sample in the same
#'   order as \code{x}'s rows (or named by sample id).
#' @param clocks Optional character vector of clock names to plot. Default: the
#'   age clocks present in \code{x}.
#' @param group Optional categorical vector (factor or character), one value
#'   per sample in \code{x}'s order, or named by sample id. Samples with a
#'   missing group are dropped.
#' @param ncol Number of facet columns. Default 3.
#' @param outlier_sd Points whose residual from the clock's linear fit exceeds
#'   this many standard deviations are highlighted, to flag predictions that fall
#'   well off the trend. \code{NULL} or \code{Inf} turns the highlight off, and
#'   it is not drawn when \code{group} is given. Default 3.
#' @return A \code{ggplot} object.
#' @examples
#' data(methylclock_demo)
#' d <- methylclock_demo
#' plotDNAmAge(d, age = d$age, clocks = c("Horvath", "Hannum", "Levine"))
#' plotDNAmAge(d, age = d$age, group = d$sex, clocks = "Horvath")
#' @export
plotDNAmAge <- function(x, age, clocks = NULL, group = NULL, ncol = 3,
                        outlier_sd = 3) {
    v <- .mc_plot_values(x, clocks)
    if (is.null(clocks))
        v <- .mc_plot_values(x, .mc_age_clocks(names(v$values)))
    age <- .mc_align_age(age, v$ids)
    if (!is.null(group))
        group <- .mc_align_group(group, v$ids)
    if (is.null(outlier_sd)) outlier_sd <- Inf

    # per-clock R^2 and, for the same fit, which points sit far off the trend
    fits <- numeric(ncol(v$values))
    outlier <- logical(0)
    for (j in seq_along(v$values)) {
        y <- v$values[[j]]
        ok <- is.finite(y) & is.finite(age)
        flag <- rep(FALSE, length(y))
        if (sum(ok) >= 3L) {
            fit <- stats::lm(y[ok] ~ age[ok])
            fits[j] <- suppressWarnings(stats::cor(age[ok], y[ok]))^2
            r <- rep(NA_real_, length(y))
            r[ok] <- stats::residuals(fit)
            s <- stats::sd(r, na.rm = TRUE)
            if (is.finite(s) && s > 0)
                flag <- is.finite(r) & abs(r) > outlier_sd * s
        } else fits[j] <- NA_real_
        outlier <- c(outlier, flag)
    }
    labs <- sprintf("%s (R\u00b2 = %.2f)", names(v$values), fits)
    names(labs) <- names(v$values)

    long <- data.frame(
        age = rep(age, times = ncol(v$values)),
        predicted = as.numeric(unlist(v$values, use.names = FALSE)),
        clock = factor(rep(labs[names(v$values)], each = length(v$ids)),
                       levels = labs[names(v$values)]),
        outlier = outlier, stringsAsFactors = FALSE)
    if (!is.null(group)) {
        long$group <- .mc_group_factor(rep(as.character(group),
                                           times = ncol(v$values)), group)
        long <- long[!is.na(long$group), , drop = FALSE]
    }

    p <- ggplot2::ggplot(long, ggplot2::aes(x = .data$age,
                                            y = .data$predicted)) +
        ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed",
                             colour = "grey55")
    if (is.null(group)) {
        p <- p +
            ggplot2::geom_point(alpha = 0.45, size = 1.3, colour = "grey40",
                                na.rm = TRUE) +
            ggplot2::geom_point(data = function(d) d[which(d$outlier), ],
                                colour = .mc_outlier, alpha = 0.8,
                                size = 1.5, na.rm = TRUE) +
            ggplot2::geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                                 colour = .mc_accent, fill = .mc_accent_soft,
                                 linewidth = 0.8, na.rm = TRUE)
    } else {
        p <- p +
            ggplot2::geom_point(ggplot2::aes(colour = .data$group),
                                alpha = 0.35, size = 1.1, na.rm = TRUE) +
            ggplot2::geom_smooth(ggplot2::aes(colour = .data$group,
                                              fill = .data$group),
                                 method = "lm", formula = y ~ x, se = TRUE,
                                 linewidth = 0.8, na.rm = TRUE) +
            .mc_colour_group() + .mc_fill_group()
    }
    p <- p + ggplot2::facet_wrap(~ clock, ncol = ncol, scales = "free_y") +
        ggplot2::labs(x = "Chronological age",
                      y = "Predicted epigenetic age") +
        theme_methylclock()
    if (!is.null(group))
        p <- p + ggplot2::theme(legend.position = "bottom")
    p
}
