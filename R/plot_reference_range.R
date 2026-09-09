#' Plot samples against the percentile bands of a reference cohort
#'
#' The growth-chart idea applied to epigenetic age: for each clock, the
#' percentile bands of a reference cohort along chronological age (5th--95th
#' light, 25th--75th darker, median line), with the user's samples drawn on
#' top. A sample inside the central band is where most of the reference sits
#' at that age; one outside the light band is where fewer than one in ten of
#' the reference samples fall. Where the reference has too few samples near
#' an age, no band is drawn --- there is no rule to measure against there.
#'
#' The bands describe \emph{the named reference cohort}, nothing more: a
#' reference measures well what resembles it (population, tissue, array,
#' processing), so the most meaningful reference is your own control samples,
#' passed as a data frame. When none is given, the function falls back to a
#' cohort shipped with the package, says so in a message, and the subtitle
#' always states which reference was used and its size. The shipped options
#' (see \code{Details}) are offered so that something can be seen with a
#' stated caveat rather than nothing without data; they are single, specific
#' cohorts and none of them stands for people in general.
#'
#' \itemize{
#'   \item \code{"GSE40279"}: 656 adults aged 19--101 (Hannum 2013, blood,
#'     450K). This is the training cohort of the Hannum, BLUP and EN clocks,
#'     so its bands for those clocks read tighter than a general cohort's
#'     would. Its reported ethnicity splits it into two further references,
#'     \code{"GSE40279-Caucasian"} (426) and \code{"GSE40279-Hispanic"}
#'     (230): a reference measures best the population it comes from.
#'   \item \code{"GSE132203"}: 795 adults aged 18--77 (Grady Trauma Project,
#'     blood, EPIC), a cohort selected for trauma exposure.
#'   \item \code{"GSE36054"}: 134 children aged 1--17 (blood, 450K).
#'   \item \code{"GSE87571"}: 729 Swedish adults aged 14--94 (blood, 450K), a
#'     population-based cohort recorded as \code{disease state: normal}.
#'   \item \code{"GSE210255"}: 1394 African American adults (GENOA, blood,
#'     EPIC), recruited through hypertensive sibships --- enriched for
#'     hypertension, with related samples, and spanning roughly ages 25--90
#'     with most samples between 40 and 75.
#'   \item \code{"GSE116379"}: 79 Han Chinese adults from Changchun (blood,
#'     450K) --- the non-schizophrenia controls of a schizophrenia and
#'     prenatal-famine study, all aged 45--51, so its bands exist only in
#'     that narrow window.
#'   \item \code{"GSE224363"}: 97 mothers from Goma, DR Congo (venous blood,
#'     EPIC), sampled around delivery in a cohort with substantial trauma
#'     exposure, aged 14--42.
#' }
#'
#' @param x A \code{methylclock} result, or a data frame with an \code{id}
#'   column and one numeric column per clock.
#' @param age Numeric vector of chronological ages, one per sample in
#'   \code{x}'s order (or named by sample id). Required: age is the axis.
#' @param reference \code{"auto"} (the default: a shipped cohort picked by the
#'   samples' age range), the name of a shipped cohort (see \code{Details}),
#'   or a data frame with an \code{age} column and one column per clock ---
#'   ideally your own control samples.
#' @param ref_label Label for the reference shown in the subtitle. Defaults to
#'   the shipped cohort's name, or \code{"user reference"} for a data frame.
#' @param clocks Optional character vector of clock names. Default: the age
#'   clocks present in both \code{x} and the reference.
#' @param window Half-width, in years, of the sliding age window the
#'   percentiles are computed in. Default 5.
#' @param min_n Minimum reference samples a window needs for its band to be
#'   drawn. Default 20.
#' @param ncol Facet columns. Default 2.
#' @return A \code{ggplot} object.
#' @examples
#' data(methylclock_demo)
#' pred <- methylclock_demo[seq_len(10), ]
#' plotReferenceRange(pred, age = pred$age)
#' @export
plotReferenceRange <- function(x, age, reference = "auto", ref_label = NULL,
                               clocks = NULL, window = 5, min_n = 20L,
                               ncol = 2) {
    v <- .mc_plot_values(x, clocks)
    if (is.null(clocks))
        v <- .mc_plot_values(x, .mc_age_clocks(names(v$values)))
    age <- .mc_align_age(age, v$ids)
    if (!any(is.finite(age)))
        stop("`age` has no finite values; age is the axis of this plot.")

    ref <- .mc_pick_reference(reference, age)
    if (is.null(ref_label)) ref_label <- ref$label
    rd <- ref$data
    if (!"age" %in% names(rd) || !is.numeric(rd$age))
        stop("The reference needs a numeric `age` column.")

    keep <- intersect(names(v$values),
                      names(rd)[vapply(rd, is.numeric, logical(1))])
    keep <- setdiff(keep, "age")
    if (!length(keep))
        stop("The reference has no clock column in common with `x`.")

    # Percentiles of the reference in a sliding window along age; a window
    # short of samples contributes no band.
    bands <- do.call(rbind, lapply(keep, function(k) {
        ok <- is.finite(rd$age) & is.finite(rd[[k]])
        if (sum(ok) < min_n) return(NULL)   # the reference lacks this clock
        ra <- rd$age[ok]
        rv <- rd[[k]][ok]
        grid <- seq(floor(min(ra)), ceiling(max(ra)), by = 1)
        do.call(rbind, lapply(grid, function(a) {
            w <- rv[abs(ra - a) <= window]
            if (length(w) < min_n) return(NULL)
            q <- stats::quantile(w, c(0.05, 0.25, 0.5, 0.75, 0.95),
                                 names = FALSE)
            data.frame(clock = k, age = a, p5 = q[1], p25 = q[2],
                       p50 = q[3], p75 = q[4], p95 = q[5],
                       stringsAsFactors = FALSE)
        }))
    }))
    if (is.null(bands) || !nrow(bands))
        stop("No age window holds at least `min_n` reference samples; ",
             "no band can be drawn.")
    bands$clock <- factor(bands$clock, levels = keep)

    pts <- do.call(rbind, lapply(keep, function(k) data.frame(
        clock = factor(k, levels = keep), age = age,
        predicted = v$values[[k]], stringsAsFactors = FALSE)))
    pts <- pts[is.finite(pts$age) & is.finite(pts$predicted), , drop = FALSE]

    sub <- sprintf("Reference: %s (n = %d, ages %d-%d)", ref_label,
                   nrow(rd), floor(min(rd$age, na.rm = TRUE)),
                   ceiling(max(rd$age, na.rm = TRUE)))

    ggplot2::ggplot(bands, ggplot2::aes(x = .data$age)) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$p5, ymax = .data$p95),
                             fill = .mc_accent_soft, alpha = 0.55) +
        ggplot2::geom_ribbon(ggplot2::aes(ymin = .data$p25,
                                          ymax = .data$p75),
                             fill = .mc_accent_mid, alpha = 0.45) +
        ggplot2::geom_line(ggplot2::aes(y = .data$p50), colour = .mc_accent,
                           linewidth = 0.6) +
        ggplot2::geom_point(data = pts,
            ggplot2::aes(x = .data$age, y = .data$predicted),
            colour = .mc_outlier, size = 1.8, na.rm = TRUE) +
        ggplot2::facet_wrap(~ .data$clock, ncol = ncol, scales = "free") +
        ggplot2::labs(x = "Chronological age (years)",
                      y = "Predicted epigenetic age (years)",
                      subtitle = sub) +
        theme_methylclock() +
        ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.line = ggplot2::element_line(colour = "grey80"))
}

# Resolve the reference argument: a user data frame is taken as is; "auto"
# picks a shipped cohort by the samples' age range and says so; a shipped
# cohort can also be named directly.
.mc_pick_reference <- function(reference, age) {
    if (is.data.frame(reference))
        return(list(data = reference, label = "user reference"))
    if (!is.character(reference) || length(reference) != 1L)
        stop("`reference` must be \"auto\", a shipped cohort name, ",
             "or a data frame.")
    demo <- get(utils::data("methylclock_demo",
        envir = environment(), package = utils::packageName())[1])
    val <- get(utils::data("methylclock_validation",
        envir = environment(), package = utils::packageName())[1])
    refs <- get(utils::data("methylclock_references",
        envir = environment(), package = utils::packageName())[1])
    shipped <- list(
        GSE40279 = demo,
        `GSE40279-Caucasian` =
            demo[demo$ethnicity == "Caucasian", , drop = FALSE],
        `GSE40279-Hispanic` =
            demo[demo$ethnicity == "Hispanic", , drop = FALSE],
        GSE132203 = val[grepl("adult", val$cohort), , drop = FALSE],
        GSE36054 = val[grepl("pediatric", val$cohort), , drop = FALSE],
        GSE87571 = refs[grepl("GSE87571", refs$cohort), , drop = FALSE],
        GSE210255 = refs[grepl("GSE210255", refs$cohort), , drop = FALSE],
        GSE116379 = refs[grepl("GSE116379", refs$cohort), , drop = FALSE],
        GSE224363 = refs[grepl("GSE224363", refs$cohort), , drop = FALSE])
    if (reference == "auto") {
        pick <- if (stats::median(age, na.rm = TRUE) < 18) "GSE36054"
            else "GSE40279"
        message("No reference given: using the shipped cohort ", pick,
                ". Bands describe that specific cohort; your own control ",
                "samples make the most meaningful reference.")
        reference <- pick
    }
    if (!reference %in% names(shipped))
        stop("Unknown shipped reference \"", reference, "\"; use one of ",
             paste(names(shipped), collapse = ", "),
             ", or pass a data frame.")
    list(data = shipped[[reference]], label = reference)
}
