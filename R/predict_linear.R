# Linear predictor engine: intercept + weighted sum of CpG betas, then the
# clock's transform. Reads only the clock's CpG rows from the beta source, in
# row blocks, so a clock that spans most of a disk-backed array still holds only
# one block at a time.

# How missing values among the present CpGs are handled before the weighted sum:
#   "mean"      each missing value -> that CpG's mean across the samples.
#   "reference" each missing value -> that CpG's reference level from the
#               clock's own coefficient table (the level it was trained
#               against), falling back to the across-sample mean when absent.
#   "none"      nothing is filled; any sample missing a used CpG is returned NA
#               rather than given an imputed number.
#   "knn"       each missing value -> the mean of that CpG's nearest CpGs (those
#               that co-vary with it across the observed samples). See
#               impute_knn.R.
# Independently, `min.perc.sample` sets a per-sample coverage floor: a sample
# with fewer than that fraction of the clock's present CpGs observed is returned
# NA for the clock, so a heavily imputed sample is not passed off as a real one.
.MC_IMPUTE <- c("mean", "reference", "none", "knn")

# Linear predictor: intercept + weighted sum of CpG betas, then the transform.
# A clock whose overall CpG coverage does not exceed `min.perc` returns NA for
# every sample; see the file header for `impute` and `min.perc.sample`.
.predict_linear <- function(entry, betas, min.perc = 0.8, impute = "mean",
                            min.perc.sample = 0, knn = .MC_KNN) {
    impute <- match.arg(impute, .MC_IMPUTE)
    coef <- mcd_resource(entry$resource)
    cpg_col <- .coef_cpg_col(coef)
    w_col <- .coef_weight_col(coef, entry$coef_col)
    markers <- as.character(coef[[cpg_col]])
    weights <- as.numeric(coef[[w_col]])
    # some models weight the SQUARE of a probe's beta as well as the beta
    # itself; a `power` column in the coefficient table marks those rows
    pw <- if ("power" %in% names(coef)) as.numeric(coef$power) else
        rep(1, length(markers))

    b0 <- 0
    if (isTRUE(entry$intercept)) {
        ic <- which(markers %in% c("(Intercept)", "Intercept", "intercept"))
        if (!length(ic)) ic <- 1L                 # legacy convention: first row
        b0 <- weights[ic[1]]
        markers <- markers[-ic[1]]
        weights <- weights[-ic[1]]
        pw <- pw[-ic[1]]
    }
    quad <- pw == 2
    mq <- markers[quad]
    wq <- weights[quad]
    markers <- markers[!quad]
    weights <- weights[!quad]

    samples <- .src_colnames(betas)
    present <- markers %in% .src_rownames(betas)
    presq <- mq %in% .src_rownames(betas)
    coverage <- mean(c(present, presq))
    if (coverage <= min.perc) {
        warning(sprintf(
            "clock '%s': %.0f%% of CpGs present (need > %.0f%%); returning NA.",
            entry$name, 100 * coverage, 100 * min.perc))
        out <- rep(NA_real_, length(samples))
        names(out) <- samples
        return(out)
    }
    if (!all(c(present, presq)))
        warning(sprintf("clock '%s': %d of %d CpGs absent; used the rest.",
                        entry$name, sum(!present) + sum(!presq),
                        length(markers) + length(mq)))

    m <- markers[present]
    w <- weights[present]
    mq2 <- mq[presq]
    wq2 <- wq[presq]
    ref <- if (identical(impute, "reference"))
        .coef_ref_means(coef, cpg_col, m) else NULL
    # A clock that uses much of a disk-backed array is read in file order --
    # scattered reads of many rows amplify on a chunked HDF5, a contiguous sweep
    # does not; a small row set reads just its rows. KNN imputation reads the raw
    # rows and fills them from co-varying CpGs before scaling: a small row set is
    # imputed in one pass (.predict_knn), a whole-array one is imputed block by
    # block so the array is never held whole (.predict_knn_seq).
    whole_array <- .src_is_hdf5(betas) &&
        length(m) > .MC_SEQ_COVERAGE * .src_nrow(betas)
    if (identical(impute, "knn"))
        acc <- if (whole_array)
            .predict_knn_seq(betas, m, w, length(samples), knn)
        else
            .predict_knn(betas, m, w, length(samples), knn)
    else if (whole_array)
        acc <- .predict_seq(betas, m, w, length(samples), impute, ref)
    else
        acc <- .predict_scattered(betas, m, w, length(samples), impute, ref)

    # the squared-term rows, when the model has them, add their own weighted
    # sum of beta^2 over the same samples (always a small scattered block)
    if (length(mq2)) {
        accq <- if (identical(impute, "knn"))
            .predict_knn(betas, mq2, wq2, length(samples), knn,
                         square = TRUE)
        else
            .predict_scattered(betas, mq2, wq2, length(samples), impute,
                               ref = NULL, square = TRUE)
        acc$pred <- acc$pred + accq$pred
        acc$na <- acc$na + accq$na
    }

    pred <- acc$pred + b0                           # length = n samples
    # per-sample coverage: honour the floor, and never invent a value under
    # "none" (any missing used CpG makes the sample NA).
    covered <- 1 - acc$na / (length(m) + length(mq2))
    drop <- covered < min.perc.sample
    if (identical(impute, "none")) drop <- drop | acc$na > 0
    pred[drop] <- NA_real_
    if (any(drop) && !identical(impute, "none"))
        warning(sprintf(
            "clock '%s': %d sample(s) below %.0f%% per-sample coverage; NA.",
            entry$name, sum(drop), 100 * min.perc.sample))
    names(pred) <- samples
    entry$transform(pred)
}

# Coverage above which a disk-backed clock is read in file order rather than by
# scattered rows. Only whole-array clocks (e.g. BLUP) cross it.
.MC_SEQ_COVERAGE <- 0.2

# Reference level per CpG from a clock's coefficient table: the training/gold
# standard mean or median the clock was calibrated against, aligned to `m`.
# Returns NULL when the table has no such column (caller then uses the mean).
.coef_ref_means <- function(coef, cpg_col, m) {
    cand <- c("goldstandard_mean", "MeansGS", "Means", "meanByCpG",
              "medianByCpG")
    col <- cand[cand %in% names(coef)]
    if (!length(col)) return(NULL)
    v <- as.numeric(coef[[col[1]]])
    names(v) <- as.character(coef[[cpg_col]])
    v[m]
}

# Mean-impute a block: each missing value -> its CpG's mean across samples (a
# CpG with no observed value contributes nothing). Used by the paths that only
# ever mean-impute (the compiled nn forward and the batched algebra).
.impute_rows <- function(B) .impute_block(B, "mean", NULL)

# Fill a block's missing values per the strategy, in place. "mean" uses each
# CpG's across-sample mean; "reference" uses its reference level (mean
# fallback); "none" leaves a placeholder (those samples are dropped to NA
# afterwards). A CpG with no observed value and no reference is filled with 0.
.impute_block <- function(B, impute, ref) {
    if (!anyNA(B)) return(B)
    fill <- rowMeans(B, na.rm = TRUE)
    if (identical(impute, "reference") && !is.null(ref)) {
        rn <- rownames(B)
        rv <- ref[rn]
        fill <- ifelse(is.finite(rv), rv, fill)
    }
    fill[!is.finite(fill)] <- 0
    na_idx <- which(is.na(B), arr.ind = TRUE)
    B[na_idx] <- fill[na_idx[, 1]]
    B
}

# Weighted sum reading only the clock's CpG rows, a block of markers at a time.
# Best when the clock uses few rows (reads just those). Returns the weighted sum
# and, per sample, how many used CpGs were missing (for the coverage checks).
.predict_scattered <- function(betas, m, w, nsamp, impute = "mean",
                               ref = NULL, square = FALSE) {
    pred <- numeric(nsamp)
    na <- numeric(nsamp)
    for (start in seq(1L, length(m), by = .MC_ROW_BLOCK)) {
        idx <- start:min(start + .MC_ROW_BLOCK - 1L, length(m))
        B <- .src_subset(betas, m[idx])
        na <- na + colSums(is.na(B))
        B <- .impute_block(B, impute, ref)
        # a squared-term block enters the weighted sum as beta^2
        if (isTRUE(square)) B <- B * B
        pred <- pred + as.numeric(crossprod(w[idx], B))
    }
    list(pred = pred, na = na)
}

# Weighted sum reading the array in file order, contiguous block by block, and
# accumulating the markers found in each block. Best when the clock uses much of
# the array on a disk backend: contiguous reads avoid chunk amplification.
.predict_seq <- function(betas, m, w, nsamp, impute = "mean", ref = NULL) {
    wmap <- tapply(w, m, sum)                       # weight per unique CpG
    base <- .src_base(betas)
    scl <- .src_scaling(betas)
    rn <- .src_rownames(betas)
    n <- length(rn)
    want <- rn %in% names(wmap)
    pred <- numeric(nsamp)
    na <- numeric(nsamp)
    for (start in seq(1L, n, by = .MC_ROW_BLOCK)) {
        end <- min(start + .MC_ROW_BLOCK - 1L, n)
        hit <- which(want[start:end])
        if (!length(hit)) next
        B <- .src_range(base, start, end)[hit, , drop = FALSE]
        if (!is.null(scl)) B <- .apply_scale(B, scl$center, scl$scale)
        na <- na + colSums(is.na(B))
        B <- .impute_block(B, impute, ref)
        pred <- pred + as.numeric(crossprod(wmap[rownames(B)], B))
    }
    list(pred = pred, na = na)
}
