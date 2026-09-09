#' Bland-Altman plot of predicted vs. chronological age
#'
#' A method-agreement view that complements the scatter in \code{\link{plotDNAmAge}}:
#' for each age clock, the difference (predicted minus real age) against the mean
#' of the two. The solid line is the mean difference (the clock's bias) and the
#' dashed lines are the 95\% limits of agreement (bias +/- 1.96 SD). A bias away
#' from zero is a miscalibrated clock; wide limits mean noisy predictions.
#'
#' @param x A \code{methylclock} result, or a data frame with an \code{id}
#'   column and one numeric column per clock.
#' @param age Numeric vector of chronological ages, one per sample in \code{x}'s
#'   order (or named by sample id).
#' @param clocks Optional clock names. Default: the age clocks in \code{x}.
#' @param ncol Number of facet columns. Default 3.
#' @return A \code{ggplot} object.
#' @examples
#' data(methylclock_demo)
#' d <- methylclock_demo
#' plotBlandAltman(d, age = d$age, clocks = c("Horvath", "Hannum", "Levine"))
#' @export
plotBlandAltman <- function(x, age, clocks = NULL, ncol = 3) {
    v <- .mc_plot_values(x, clocks)
    if (is.null(clocks))
        v <- .mc_plot_values(x, .mc_age_clocks(names(v$values)))
    age <- .mc_align_age(age, v$ids)

    parts <- lapply(names(v$values), function(nm) {
        y <- v$values[[nm]]
        ok <- is.finite(y) & is.finite(age)
        data.frame(clock = nm, mean = (y[ok] + age[ok]) / 2,
                   diff = y[ok] - age[ok], stringsAsFactors = FALSE)
    })
    long <- do.call(rbind, parts)
    long$clock <- factor(long$clock, levels = names(v$values))

    stats <- do.call(rbind, lapply(split(long, long$clock), function(d) {
        b <- mean(d$diff); s <- stats::sd(d$diff)
        data.frame(clock = d$clock[1], bias = b, lo = b - 1.96 * s,
                   hi = b + 1.96 * s, stringsAsFactors = FALSE)
    }))

    ggplot2::ggplot(long, ggplot2::aes(x = .data$mean, y = .data$diff)) +
        ggplot2::geom_hline(yintercept = 0, colour = "grey80") +
        ggplot2::geom_point(alpha = 0.4, size = 1.2, colour = "grey40",
                            na.rm = TRUE) +
        ggplot2::geom_hline(data = stats, ggplot2::aes(yintercept = .data$bias),
                            colour = .mc_accent, linewidth = 0.7) +
        ggplot2::geom_hline(data = stats, ggplot2::aes(yintercept = .data$lo),
                            colour = .mc_accent, linetype = "dashed") +
        ggplot2::geom_hline(data = stats, ggplot2::aes(yintercept = .data$hi),
                            colour = .mc_accent, linetype = "dashed") +
        ggplot2::facet_wrap(~ clock, ncol = ncol, scales = "free") +
        ggplot2::labs(x = "Mean of predicted and real age",
                      y = "Predicted - real age (years)") +
        theme_methylclock()
}
