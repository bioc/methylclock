# Longitudinal support: shape repeated samples of the same subjects into
# trajectories, and summarise each subject's rate of epigenetic aging.

#' Arrange repeated samples into per-subject clock trajectories
#'
#' When the same subjects were sampled at several ages, this pairs their
#' samples into trajectories: one row per subject, visit and clock, ordered
#' in time. It is the table the longitudinal plots and rates are built from,
#' and the one to hand to a mixed-model package for formal analysis --- for
#' example \code{lme4::lmer(estimate ~ age + (age | subject), data = tr)}
#' fits a population slope with per-subject deviations. This package shapes
#' and displays longitudinal data; the inference belongs to those tools.
#'
#' Samples with a missing subject or age are dropped; repeated measurements
#' of a subject at the same recorded age are averaged (whether they are
#' technical replicates or distinct samples whose ages rounded to the same
#' value is for the caller to know --- pass more precise ages to keep them
#' apart); and any subject-and-clock combination seen at fewer than two
#' distinct ages is excluded --- one point has no trajectory. Each exclusion
#' is reported in a message.
#'
#' @param x A \code{methylclock} result, or a data frame with an \code{id}
#'   column and one numeric column per clock.
#' @param subject Subject identifier (character or factor), one per sample in
#'   \code{x}'s order, or named by sample id.
#' @param age Numeric chronological age AT EACH SAMPLE, one per sample in
#'   \code{x}'s order (or named by sample id). This is the time axis.
#' @param clocks Optional character vector of clock names. Default: the age
#'   clocks present in \code{x}.
#' @return A data frame with columns \code{subject}, \code{visit} (1, 2, ...
#'   in age order), \code{age}, \code{clock} and \code{estimate}.
#' @examples
#' data(methylclock_longitudinal)
#' d <- methylclock_longitudinal
#' tr <- clockTrajectories(d, subject = d$subject, age = d$age,
#'                         clocks = "Wu")
#' head(tr)
#' @export
clockTrajectories <- function(x, subject, age, clocks = NULL) {
    v <- .mc_plot_values(x, clocks)
    if (is.null(clocks))
        v <- .mc_plot_values(x, .mc_age_clocks(names(v$values)))
    subject <- as.character(.mc_align_group(subject, v$ids))
    age <- .mc_align_age(age, v$ids)

    ok <- !is.na(subject) & is.finite(age)
    if (any(!ok))
        message(sum(!ok), " sample(s) without subject or age were dropped.")
    vals <- v$values[ok, , drop = FALSE]
    subject <- subject[ok]
    age <- age[ok]

    long <- data.frame(
        subject = rep(subject, times = ncol(vals)),
        age = rep(age, times = ncol(vals)),
        clock = factor(rep(names(vals), each = length(subject)),
                       levels = names(vals)),
        estimate = as.numeric(unlist(vals, use.names = FALSE)),
        stringsAsFactors = FALSE)
    long <- long[is.finite(long$estimate), , drop = FALSE]

    # the same subject measured more than once at the same recorded age
    key <- paste(long$subject, long$age, long$clock, sep = "\r")
    if (anyDuplicated(key)) {
        message("Repeated measurements at the same recorded age were ",
                "averaged.")
        long <- stats::aggregate(estimate ~ subject + age + clock, long,
                                 mean)
    }

    # one point has no trajectory: a trajectory lives per subject AND clock
    key2 <- paste(long$subject, long$clock, sep = "\r")
    n_ages <- tapply(long$age, key2, function(a) length(unique(a)))
    single <- names(n_ages)[n_ages < 2L]
    if (length(single)) {
        message(length(single), " subject-clock trajectorie(s) with a ",
                "single age were excluded.")
        long <- long[!key2 %in% single, , drop = FALSE]
    }
    if (!nrow(long))
        stop("No subject was seen at two or more ages.")

    long <- long[order(long$subject, long$clock, long$age), , drop = FALSE]
    long$visit <- stats::ave(long$age,
                             long$subject, long$clock,
                             FUN = function(a) rank(a, ties.method = "first"))
    rownames(long) <- NULL
    long[, c("subject", "visit", "age", "clock", "estimate")]
}

#' Each subject's rate of epigenetic aging, clock by clock
#'
#' For every subject and clock, the within-subject slope of the estimate:
#' how many epigenetic years it advances per chronological year across that
#' subject's visits. A slope of 1 tracks the calendar. This is a property of
#' the estimate's trajectory, not a validated measure of biological ageing
#' pace --- of the shipped clocks, only DunedinPACE was designed and
#' validated to estimate such a pace, from a single sample.
#' With two visits the slope is the difference quotient and carries no
#' per-subject uncertainty (\code{se} is \code{NA}: read it as orientative);
#' with three or more it is the slope of a least-squares line with its
#' standard error.
#'
#' Read individual slopes with care: the technical replication error of many
#' clocks spans several years (see Higgins-Chen et al. 2022, Nature Aging,
#' on the reliability of clock estimates), and the difference of two
#' measurements carries the noise of both, so over a short follow-up a
#' subject's slope can be mostly measurement. The steadier signal is in the
#' aggregate --- the distribution of slopes in a group, or their group
#' difference (\code{\link{plotGroupDifference}} accepts these rates: one
#' row per subject).
#'
#' @inheritParams clockTrajectories
#' @return A data frame with one row per subject and clock: \code{subject},
#'   \code{clock}, \code{n_visits}, \code{span_years} (age range covered),
#'   \code{slope} and \code{se} (\code{NA} with two visits).
#' @examples
#' data(methylclock_longitudinal)
#' d <- methylclock_longitudinal
#' rates <- trajectoryRates(d, subject = d$subject, age = d$age,
#'                          clocks = "Wu")
#' head(rates)
#' @export
trajectoryRates <- function(x, subject, age, clocks = NULL) {
    tr <- clockTrajectories(x, subject, age, clocks)
    pieces <- split(tr, list(tr$subject, tr$clock), drop = TRUE)
    rows <- lapply(pieces, function(d) {
        n <- nrow(d)
        span <- max(d$age) - min(d$age)
        if (n < 2L || span <= 0)
            return(NULL)
        if (n == 2L) {
            slope <- diff(d$estimate) / diff(d$age)
            se <- NA_real_
        } else {
            fit <- stats::lm(estimate ~ age, data = d)
            slope <- unname(stats::coef(fit)[2])
            se <- summary(fit)$coefficients[2, 2]
        }
        data.frame(subject = d$subject[1], clock = d$clock[1],
                   n_visits = n, span_years = span, slope = slope, se = se,
                   stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows)
    if (is.null(out))
        stop("No subject has two visits at distinct ages.")
    out <- out[order(out$subject, out$clock), , drop = FALSE]
    rownames(out) <- NULL
    out
}
