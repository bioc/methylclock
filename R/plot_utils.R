# Shared helpers for the clock plots. The plots consume a methylclock result (or
# any data frame of per-sample clock estimates) and return a ggplot object, so a
# caller can theme, save or combine them.

# Pull a result into a samples x clocks numeric matrix plus the sample ids.
# `x` is a methylclock result or a data frame with id and clock columns;
# `clocks` optionally restricts to a subset (kept in that order).
.mc_plot_values <- function(x, clocks = NULL) {
    df <- as.data.frame(x)
    ids <- if ("id" %in% names(df)) as.character(df$id) else rownames(df)
    df <- df[, setdiff(names(df), "id"), drop = FALSE]
    num <- vapply(df, is.numeric, logical(1))
    df <- df[, num, drop = FALSE]
    if (!is.null(clocks)) {
        keep <- intersect(clocks, names(df))
        if (!length(keep))
            stop("None of the requested clocks are in the result.")
        df <- df[, keep, drop = FALSE]
    }
    if (!ncol(df))
        stop("No numeric clock columns found to plot.")
    list(ids = ids, values = df)
}

# Registry targets whose output is an age in years (a predicted-vs-real-age plot
# only makes sense for these; a mitotic count or a trait score is not an age).
.MC_AGE_TARGETS <- c("chronological", "phenotypic", "gestational", "pace")

# The age-in-years clocks among a set of column names, from the registry. Falls
# back to all of them if the registry is not reachable (e.g. a bare data frame).
.mc_age_clocks <- function(cols) {
    ok <- vapply(cols, function(nm) {
        info <- tryCatch(clock_info(nm), error = function(e) NULL)
        !is.null(info) && info$target %in% .MC_AGE_TARGETS
    }, logical(1))
    if (any(ok)) cols[ok] else cols
}

# Align an age (or any per-sample covariate) to the plot's sample order. A named
# vector is reordered by id; an unnamed one must already match the sample order.
.mc_align_age <- function(age, ids) {
    nm <- names(age)
    age <- as.numeric(age)
    names(age) <- nm
    if (!is.null(nm) && all(ids %in% nm))
        return(unname(age[ids]))
    if (length(age) != length(ids))
        stop("`age` must have one value per sample (", length(ids),
             "), or be named by sample id.")
    age
}

# Align a categorical (or any non-numeric) per-sample covariate to the plot's
# sample order, without coercing its type. A named vector is reordered by id; an
# unnamed one must already match the sample order.
.mc_align_group <- function(group, ids) {
    nm <- names(group)
    if (!is.null(nm) && all(ids %in% nm))
        return(unname(group[ids]))
    if (length(group) != length(ids))
        stop("`group` must have one value per sample (", length(ids),
             "), or be named by sample id.")
    group
}

# Age acceleration per clock: the residual of predicted age regressed on
# chronological age, so a positive value is a sample older than its years.
# `values` is a samples x clocks data frame; `age` is aligned to its rows.
# Returns a samples x clocks matrix (NA where a clock cannot be fitted).
.mc_age_acceleration <- function(values, age) {
    ids <- rownames(values)
    if (is.null(ids)) ids <- as.character(seq_len(nrow(values)))
    accel <- vapply(names(values), function(nm) {
        y <- values[[nm]]
        r <- rep(NA_real_, length(y))
        ok <- is.finite(y) & is.finite(age)
        if (sum(ok) >= 3L)
            r[ok] <- stats::residuals(stats::lm(y[ok] ~ age[ok]))
        r
    }, numeric(nrow(values)))
    matrix(accel, nrow = nrow(values), dimnames = list(ids, names(values)))
}

# Long form (id, clock, value) for the requested clocks, clocks as an ordered
# factor so facets/axes keep the given order.
.mc_plot_long <- function(x, clocks = NULL) {
    v <- .mc_plot_values(x, clocks)
    long <- data.frame(
        id = rep(v$ids, times = ncol(v$values)),
        clock = factor(rep(names(v$values), each = length(v$ids)),
                       levels = names(v$values)),
        value = as.numeric(unlist(v$values, use.names = FALSE)),
        stringsAsFactors = FALSE)
    long
}
