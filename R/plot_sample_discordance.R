#' Plot how much the clocks agree on each sample
#'
#' One point per sample, placed by two summaries of its clock estimates once
#' each clock is z-scored across samples (so clocks on different scales become
#' comparable): the consensus --- the median z-score, where the clocks jointly
#' place the sample relative to the rest --- and the discordance --- the
#' standard deviation of those z-scores, how much the clocks contradict each
#' other on that sample. Most samples form a low band; a sample far above it
#' is one the clocks genuinely disagree about, which can be a technical
#' problem with that sample or a biology worth a second look. The
#' \code{label_n} most discordant samples are highlighted and labelled by id.
#'
#' This complements \code{\link{sampleQC}}, which scores how jointly
#' \emph{extreme} a sample is: a sample can be extreme with every clock in
#' agreement, or unremarkable on average while the clocks disagree.
#'
#' With \code{what = "acceleration"} (requires \code{age}) the same two
#' summaries are computed on the age acceleration instead, so the consensus
#' reads as how jointly accelerated the sample is rather than how old.
#'
#' @param x A \code{methylclock} result, or a data frame with an \code{id}
#'   column and one numeric column per clock.
#' @param clocks Optional character vector of clock names. Default: the clocks
#'   in \code{x} known to the registry (for \code{what = "acceleration"}, the
#'   age clocks). At least three are needed for a disagreement to mean
#'   anything.
#' @param what Summarise the raw \code{"estimate"} (the default) or the
#'   \code{"acceleration"}, which requires \code{age}.
#' @param age Numeric vector of chronological ages, one per sample in
#'   \code{x}'s order (or named by sample id). Only used --- and required ---
#'   when \code{what = "acceleration"}.
#' @param label_n How many of the most discordant samples to highlight and
#'   label. Default 5; 0 turns the highlight off.
#' @return A \code{ggplot} object.
#' @examples
#' data(methylclock_demo)
#' plotSampleDiscordance(methylclock_demo)
#' @export
plotSampleDiscordance <- function(x, clocks = NULL,
                                  what = c("estimate", "acceleration"),
                                  age = NULL, label_n = 5L) {
    what <- match.arg(what)
    label_n <- max(0L, as.integer(label_n))
    v <- .mc_plot_values(x, clocks)
    if (is.null(clocks)) {
        keep <- if (what == "acceleration")
            .mc_age_clocks(names(v$values))
        else
            names(v$values)[.mc_clock_family(names(v$values)) != "other"]
        if (length(keep)) v <- .mc_plot_values(x, keep)
    }
    if (ncol(v$values) < 3L)
        stop("At least three clocks are needed to measure how much they ",
             "disagree; got ", ncol(v$values), ".")
    vals <- v$values
    if (what == "acceleration") {
        if (is.null(age))
            stop("`what = \"acceleration\"` needs `age`.")
        age <- .mc_align_age(age, v$ids)
        rownames(vals) <- v$ids
        vals <- as.data.frame(.mc_age_acceleration(vals, age))
    }

    z <- scale(as.matrix(vals))
    df <- data.frame(
        id = v$ids,
        consensus = apply(z, 1, stats::median, na.rm = TRUE),
        discordance = apply(z, 1, stats::sd, na.rm = TRUE),
        stringsAsFactors = FALSE)
    df <- df[is.finite(df$consensus) & is.finite(df$discordance), ,
             drop = FALSE]
    if (!nrow(df))
        stop("No sample has enough finite values across the clocks.")
    top <- df[utils::head(order(-df$discordance), label_n), , drop = FALSE]

    xlab <- if (what == "acceleration")
        "Consensus of age acceleration (median z-score)"
    else
        "Consensus across clocks (median z-score)"
    p <- ggplot2::ggplot(df, ggplot2::aes(x = .data$consensus,
                                          y = .data$discordance)) +
        ggplot2::geom_point(alpha = 0.5, size = 1.6, colour = "grey40",
                            na.rm = TRUE)
    if (nrow(top))
        p <- p +
            ggplot2::geom_point(data = top, colour = .mc_outlier,
                                size = 2.2) +
            ggplot2::geom_text(data = top,
                ggplot2::aes(label = .data$id), vjust = -0.8, size = 2.8,
                colour = .mc_outlier, check_overlap = TRUE)
    p + ggplot2::labs(x = xlab,
            y = "Discordance across clocks (SD of z-scores)") +
        theme_methylclock() +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "grey80"))
}
