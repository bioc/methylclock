# Plain-language notes that go with each section of the quality-control report:
# what the numbers in that section are, what they read as on this dataset, and
# what mechanism produces such a pattern. Descriptive throughout --- the report
# states what is there and leaves the judgement to the reader.

# Summary statistics the commentary needs and the tables do not carry: how
# closely the clocks track one another, where each one centres its estimates,
# and the average signed gap to chronological age.
.qc_stats <- function(v, age = NULL, age_clk = NULL) {
    m <- as.matrix(as.data.frame(v$values))
    list(agreement = .qc_agreement(m),
         spread = .qc_spread(m),
         bias = if (!is.null(age) && length(age_clk))
             .qc_bias(m[, intersect(age_clk, colnames(m)), drop = FALSE], age))
}

# Correlation between clocks: the typical pair, the closest and the furthest.
.qc_agreement <- function(m) {
    if (ncol(m) < 2L) return(NULL)
    cm <- suppressWarnings(stats::cor(m, use = "pairwise.complete.obs"))
    diag(cm) <- NA
    if (!any(is.finite(cm))) return(NULL)
    pick <- function(target) {
        i <- which(cm == target, arr.ind = TRUE)[1, ]
        list(r = target,
             clocks = sort(c(rownames(cm)[i[[1]]], colnames(cm)[i[[2]]])))
    }
    list(median_r = stats::median(cm[upper.tri(cm)], na.rm = TRUE),
         closest = pick(max(cm, na.rm = TRUE)),
         furthest = pick(min(cm, na.rm = TRUE)))
}

# Where each clock centres its estimates (lowest and highest median).
.qc_spread <- function(m) {
    med <- apply(m, 2, stats::median, na.rm = TRUE)
    med <- med[is.finite(med)]
    if (!length(med)) return(NULL)
    list(low = names(med)[which.min(med)], low_v = min(med),
         high = names(med)[which.max(med)], high_v = max(med))
}

# Average signed gap to chronological age, and whether that gap drifts across
# the age range (the tilt a Bland-Altman plot shows).
.qc_bias <- function(m, age) {
    if (!ncol(m)) return(NULL)
    rows <- lapply(colnames(m), function(nm) {
        d <- m[, nm] - age
        ok <- is.finite(d)
        if (sum(ok) < 3L) return(NULL)
        avg <- (m[ok, nm] + age[ok]) / 2
        data.frame(clock = nm, bias = mean(d[ok]),
                   trend = suppressWarnings(stats::cor(avg, d[ok])),
                   stringsAsFactors = FALSE)
    })
    rows <- do.call(rbind, rows)
    if (is.null(rows) || !nrow(rows)) NULL else rows
}

# One entry per report section, each a short vector of sentences.
.qc_notes <- function(qc) {
    st <- qc$stats
    notes <- list(overview = .qc_note_overview(qc))
    notes$variables <- .qc_note_variables(qc)
    notes$coverage <- .qc_note_coverage(qc)
    notes$missing <- .qc_note_missing(qc)
    notes$agreement <- c(.qc_note_agreement(st$agreement), .qc_note_shown(qc))
    notes$distributions <- c(.qc_note_spread(st$spread), .qc_note_shown(qc))
    notes$accuracy <- .qc_note_accuracy(qc)
    notes$bland_altman <- .qc_note_bias(st$bias)
    notes$samples <- .qc_note_samples(qc)
    notes[!vapply(notes, is.null, logical(1))]
}

# State which sample variables the report was given, so a report without a
# breakdown by group is not read as a dataset without groups.
.qc_note_variables <- function(qc) {
    have <- c(if (isTRUE(qc$has_age)) "chronological age",
              if (!is.null(qc$groups)) "a grouping variable")
    if (length(have) == 2L)
        return(sprintf(paste("Sample variables supplied: age and a grouping",
            "variable (%s), so the report breaks the accuracy down by group."),
            paste(names(qc$groups), collapse = ", ")))
    if (!length(have))
        return(paste("No sample variables were supplied beyond the estimates",
            "themselves. Passing chronological age adds the error and bias",
            "sections, and a grouping variable (sex, cohort, case/control)",
            "breaks the results down by group."))
    if (is.null(qc$groups))
        return(paste("No grouping variable was supplied, so the results are",
            "not broken down by subgroup. Passing one (sex, cohort,",
            "case/control) adds that breakdown."))
    NULL
}

# Say which clocks the figures leave out, so a shorter plot is not mistaken
# for a shorter analysis.
.qc_note_shown <- function(qc) {
    n <- length(qc$omitted)
    if (!n) return(NULL)
    sprintf(paste("This figure shows the %d clocks that returned values; %d",
        "more (%s) are left out because they returned none for these samples.",
        "The tables above still cover every clock."),
        length(qc$shown), n, .qc_list(qc$omitted))
}

# Name a handful of items and count the rest, so a long list stays readable.
.qc_list <- function(x, n = 6L) {
    if (length(x) <= n) return(paste(x, collapse = ", "))
    sprintf("%s and %d more", paste(utils::head(x, n), collapse = ", "),
            length(x) - n)
}

.qc_note_overview <- function(qc) {
    if (isTRUE(qc$has_age)) {
        a <- qc$samples$age
        msg <- sprintf(paste("This report covers %d samples with chronological",
            "ages from %.0f to %.0f years (median %.0f)."),
            qc$samples$n, a[["Min."]], a[["Max."]], a[["Median"]])
    } else {
        msg <- sprintf(paste("This report covers %d samples. No chronological",
            "age was supplied, so the sections that compare estimates against",
            "age --- error, bias and acceleration --- are not part of it."),
            qc$samples$n)
    }
    ntot <- length(qc$computed)
    part <- names(qc$computed)[qc$computed < qc$samples$n]
    c(msg, if (length(part))
        sprintf(paste("%d of %d clocks returned a value for every sample. The",
            "rest returned NA for some or all of them (%s), which happens when",
            "the CpGs a clock needs are not in the data."),
            ntot - length(part), ntot, .qc_list(part))
      else sprintf("All %d clocks returned a value for every sample.", ntot))
}

.qc_note_coverage <- function(qc) {
    cov <- qc$coverage
    if (is.null(cov) || !nrow(cov)) return(NULL)
    lo <- cov[which.min(cov$pct_present), ]
    full <- sum(cov$pct_present >= 100)
    c(paste("Coverage is the share of a clock's own CpGs that your data",
            "carry."),
      sprintf(paste("It runs from %.0f%% (%s, %d of %d CpGs) to %.0f%%; %d of",
            "%d clocks have all of their CpGs present."),
            min(cov$pct_present), lo$clock, lo$present, lo$n_cpgs,
            max(cov$pct_present), full, nrow(cov)),
      paste("A clock still returns a number when sites are absent --- they are",
            "filled in before the estimate is computed --- so the lower this",
            "share, the more of the estimate rests on filled-in values rather",
            "than on measurements."),
      if (any(cov$pct_present == 0))
          paste("A clock with none of its CpGs present returns NA instead of a",
                "number: there is nothing to fill in from."))
}

.qc_note_missing <- function(qc) {
    mbc <- qc$missing_by_clock
    how <- if (!is.null(qc$impute))
        sprintf("filled by the '%s' method", qc$impute) else "filled in"
    if (is.null(mbc) || !nrow(mbc))
        return(if (!is.null(qc$impute))
            sprintf("Missing values were %s.", how))
    if (!any(mbc$n_missing > 0))
        return(c(sprintf(paste("Among the CpGs that are present, no values",
            "were missing, so nothing was %s."), how)))
    top <- mbc[which.max(mbc$pct_missing), ]
    c(sprintf(paste("Among the CpGs that are present, some values were missing",
            "in the matrix and were %s before the clocks ran."), how),
      sprintf(paste("This peaks at %.1f%% of the values for %s; %d of %d",
            "clocks had none missing."), top$pct_missing, top$clock,
            sum(mbc$n_missing == 0), nrow(mbc)),
      paste("A filled value carries no measurement of its own: it is a",
            "stand-in derived from the rest of the data, so it moves the",
            "estimate towards the dataset rather than away from it."))
}

.qc_note_agreement <- function(ag) {
    if (is.null(ag)) return(NULL)
    neg <- if (ag$furthest$r < 0)
        paste("A negative value means the two move in opposite directions as",
              "samples change.") else NULL
    c(paste("Clocks correlate with one another because they read the same",
            "ageing signal from different sets of CpGs."),
      sprintf(paste("Here the median correlation between pairs of clocks is",
            "%.2f. The pair that tracks each other most closely is %s and %s",
            "(%.2f); the pair furthest apart is %s and %s (%.2f)."),
            ag$median_r,
            ag$closest$clocks[1], ag$closest$clocks[2], ag$closest$r,
            ag$furthest$clocks[1], ag$furthest$clocks[2], ag$furthest$r),
      neg,
      paste("Clocks built for different targets --- calendar age, cell",
            "divisions, a trait --- measure different quantities, so a low",
            "correlation between two of them follows from their design and",
            "not from the data."))
}

.qc_note_spread <- function(sp) {
    if (is.null(sp)) return(NULL)
    c(paste("This shows where each clock places your samples and how widely it",
            "spreads them."),
      sprintf("Median estimates run from %.1f (%s) to %.1f (%s).",
              sp$low_v, sp$low, sp$high_v, sp$high),
      paste("Clocks reported in different units do not share a scale --- years",
            "of age, weeks of gestation and cumulative cell divisions are not",
            "comparable numbers --- so a clock sitting apart from the rest is",
            "often on its own scale rather than disagreeing about the",
            "samples."))
}

.qc_note_accuracy <- function(qc) {
    fin <- qc$accuracy
    if (is.null(fin)) return(NULL)
    fin <- fin[is.finite(fin$R2), , drop = FALSE]
    if (!nrow(fin)) return(NULL)
    lab <- if ("group" %in% names(fin))
        paste0(fin$clock, " [", fin$group, "]") else fin$clock
    hi <- which.max(fin$r); lo <- which.min(fin$r)
    mhi <- which.max(fin$MAE); mlo <- which.min(fin$MAE)
    c(paste("`r` is how tightly a clock's estimate follows chronological age,",
            "`R2` the share of the age variation it accounts for, and `MAE`",
            "the typical gap between estimate and age, in years."),
      sprintf(paste("Across these clocks r runs from %.2f (%s) to %.2f (%s),",
            "and the median absolute error from %.1f years (%s) to %.1f (%s)."),
            fin$r[lo], lab[lo], fin$r[hi], lab[hi],
            fin$MAE[mlo], lab[mlo], fin$MAE[mhi], lab[mhi]),
      paste("The two columns say different things, and a clock can score high",
            "on one and low on the other: a high r with a large MAE is the",
            "signature of a clock that orders the samples by age correctly",
            "while sitting on a different scale for this dataset, which is",
            "what happens when a predictor is applied to a population or",
            "tissue other than the one it was fitted on."))
}

.qc_note_bias <- function(bias) {
    if (is.null(bias)) return(NULL)
    hi <- which.max(bias$bias); lo <- which.min(bias$bias)
    tilt <- is.finite(bias$trend) & abs(bias$trend) > 0.3
    drift <- bias$clock[tilt][order(-abs(bias$trend[tilt]))]
    c(paste("The bias is the average signed gap between estimate and",
            "chronological age: above zero the clock reads older than the",
            "calendar, below zero younger. The dashed lines hold 95% of the",
            "samples."),
      sprintf("It runs from %+.1f years (%s) to %+.1f years (%s).",
              bias$bias[lo], bias$clock[lo], bias$bias[hi], bias$clock[hi]),
      if (length(drift))
          sprintf(paste("For %d of the %d clocks the gap drifts across the age",
              "range rather than staying constant (the cloud tilts); it does",
              "so most in %s. A drifting gap means the offset differs between",
              "the youngest and the oldest samples, which is what a clock",
              "fitted on a narrower age range than yours produces."),
              length(drift), nrow(bias), paste(utils::head(drift, 3),
                                               collapse = ", "))
      else paste("No clock's gap drifts markedly across the age range: the",
                 "offsets are roughly constant from the youngest samples to",
                 "the oldest."))
}

.qc_note_samples <- function(qc) {
    sq <- qc$sample_qc
    if (is.null(sq) || !nrow(sq)) return(NULL)
    n_out <- sum(sq$is_outlier, na.rm = TRUE)
    msg <- c(paste("A sample is flagged when its estimates sit far from the",
            "rest of the dataset --- a mean absolute z-score above 3 across",
            "the clocks."),
      sprintf("%d of %d samples are flagged.", n_out, nrow(sq)),
      paste("The flag is relative to this dataset: it marks a sample as",
            "unusual among these samples, which can come from the assay or",
            "from the sample being genuinely atypical."))
    if (!is.null(sq$pct_missing) && any(is.finite(sq$pct_missing)))
        msg <- c(msg, sprintf(paste("The share of clock CpGs missing per",
            "sample has a median of %.1f%% and a maximum of %.1f%%."),
            stats::median(sq$pct_missing, na.rm = TRUE),
            max(sq$pct_missing, na.rm = TRUE)))
    msg
}

# The short bullet list at the top of the report: one line per section.
.qc_commentary <- function(qc) {
    notes <- qc$notes
    if (is.null(notes)) notes <- .qc_notes(qc)
    key <- list(overview = seq_len(2L), coverage = 2L, missing = 2L,
                agreement = 2L,
                accuracy = 2L, samples = 2L)
    key <- key[names(key) %in% names(notes)]
    out <- unlist(lapply(names(key), function(nm) {
        s <- notes[[nm]]
        s[intersect(key[[nm]], seq_along(s))]
    }), use.names = FALSE)
    if (!is.null(qc$groups))
        out <- c(out, sprintf("Groups compared (sizes): %s.",
            paste(names(qc$groups), as.integer(qc$groups),
                  sep = " n=", collapse = ", ")))
    out
}
