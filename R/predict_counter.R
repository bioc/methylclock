# Counter predictor engine: summarise a CpG set per sample (no weights, no
# intercept), then the clock's transform. Only the set's rows are read; sets are
# small (hundreds of CpGs), so peak memory stays bounded by the set. Missing
# values are dropped per aggregation (na.rm), not imputed, matching how these
# summaries are defined.

# Counter predictor: the mean of a CpG set, one minus that mean (via a
# transform), an upper quantile, or a parametric per-CpG map (epiTOC2).
.predict_counter <- function(entry, betas, min.perc = 0.8) {
    res <- mcd_resource(entry$resource)
    markers <- as.character(res[[.coef_cpg_col(res)]])
    markers <- markers[!markers %in% c("(Intercept)", "Intercept", "intercept")]

    samples <- .src_colnames(betas)
    present <- markers %in% .src_rownames(betas)
    coverage <- mean(present)
    if (coverage <= min.perc) {
        warning(sprintf(
            "clock '%s': %.0f%% of CpGs present (need > %.0f%%); returning NA.",
            entry$name, 100 * coverage, 100 * min.perc))
        out <- rep(NA_real_, length(samples))
        names(out) <- samples
        return(out)
    }
    if (!all(present))
        warning(sprintf("clock '%s': %d of %d CpGs absent; computed on the rest.",
                        entry$name, sum(!present), length(markers)))

    B <- .src_subset(betas, markers[present])
    if (all(c("delta", "beta0") %in% names(res)))
        B <- .counter_parametric(B, res, markers[present])
    pred <- .aggregate_counter(B, entry$agg)
    names(pred) <- samples
    entry$transform(pred)
}

# Parametric mitotic counter (epiTOC2): before aggregating, replace each CpG's
# beta by (beta - beta0) / (delta * (1 - beta0)) using the per-CpG rate `delta`
# and ground-state `beta0` carried in the resource. The clock's `transform`
# supplies the remaining scalar factor (e.g. x * 2). A per-CpG vector is recycled
# down the columns, so the map is applied row-wise (one CpG per row).
.counter_parametric <- function(B, res, cpgs) {
    idx <- match(cpgs, res$CpGmarker)
    delta <- res$delta[idx]
    beta0 <- res$beta0[idx]
    (B - beta0) / (delta * (1 - beta0))
}

# Summarise a CpGs x samples block to one value per sample, dropping NAs.
.aggregate_counter <- function(B, agg) {
    if (is.null(agg) || identical(agg, "mean"))
        return(colMeans(B, na.rm = TRUE))
    if (is.list(agg) && identical(agg$fun, "quantile"))
        return(apply(B, 2L, stats::quantile, probs = agg$p, na.rm = TRUE,
                     names = FALSE))
    stop("Unknown counter aggregation.")
}
