#' Plot epigenetic age acceleration by a sample group
#'
#' Splits the per-clock age acceleration (the residual of predicted age on
#' chronological age) by a categorical sample variable, such as sex or a
#' case/control label. A group whose boxes sit consistently above zero looks
#' epigenetically older than its chronological age across clocks, which is the
#' association these plots are meant to surface.
#'
#' @param x A \code{methylclock} result, or a data frame with an \code{id}
#'   column and one numeric column per clock.
#' @param age Numeric vector of chronological ages, one per sample in \code{x}'s
#'   order (or named by sample id).
#' @param group Categorical vector (factor or character), one value per sample
#'   in \code{x}'s order, or named by sample id. Samples with a missing group
#'   are dropped.
#' @param clocks Optional character vector of clock names. Default: the age
#'   clocks present in \code{x}.
#' @return A \code{ggplot} object.
#' @examples
#' data(methylclock_demo)
#' d <- methylclock_demo
#' plotAccelerationByGroup(d, age = d$age, group = d$sex,
#'                         clocks = c("Horvath", "Hannum", "Levine"))
#' @export
plotAccelerationByGroup <- function(x, age, group, clocks = NULL) {
    v <- .mc_plot_values(x, clocks)
    if (is.null(clocks))
        v <- .mc_plot_values(x, .mc_age_clocks(names(v$values)))
    age <- .mc_align_age(age, v$ids)
    group <- .mc_align_group(group, v$ids)
    rownames(v$values) <- v$ids
    accel <- .mc_age_acceleration(v$values, age)

    long <- data.frame(
        clock = factor(rep(colnames(accel), each = nrow(accel)),
                       levels = colnames(accel)),
        group = factor(rep(as.character(group), times = ncol(accel))),
        accel = as.numeric(accel), stringsAsFactors = FALSE)
    long <- long[!is.na(long$group) & is.finite(long$accel), , drop = FALSE]

    ggplot2::ggplot(long, ggplot2::aes(x = .data$clock, y = .data$accel,
                                       fill = .data$group)) +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed",
                            colour = "grey60") +
        ggplot2::geom_boxplot(outlier.size = 0.5, position =
                                  ggplot2::position_dodge(width = 0.8)) +
        .mc_fill_group() +
        ggplot2::labs(x = NULL, y = "Age acceleration (years)") +
        theme_methylclock() +
        ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45,
            hjust = 1), legend.position = "top")
}
