# Batched linear predictor over an on-disk methylation matrix, with the algebra
# run through BigDataStatMeth's HDF5 cross-product. When the input is disk-backed,
# every linear clock that shares the same preprocessing is a column of one
# coefficient matrix, and all of them are estimated in a SINGLE on-disk
# cross-product t(coef) %*% betas, reading the shared CpGs once instead of once
# per clock. Clock-specific steps (coverage rule, intercept, output transform)
# stay per clock, applied to the small clocks-by-samples result.
#
# INTERNAL: D20. The algebra follows the input backend: a disk-backed input runs
# INTERNAL: the linear algebra through BigDataStatMeth (crossprod on HDF5); the
# INTERNAL: in-memory path is unchanged. Whole-array clocks (very high coverage,
# INTERNAL: e.g. BLUP) keep the streaming base path where it is lighter; only the
# INTERNAL: sparse groups batch here. See docs/D20_ALGEBRA_BDSM.md.

# Threads for the on-disk algebra (overridable by option).
.mc_bdsm_threads <- function() {
    as.integer(getOption("methylclock.bdsm_threads", 4L))
}

# Cross-product block size in rows per block (overridable by option). A moderate
# value suits the wide (few CpGs x many samples) operands the batches produce.
.mc_bdsm_block <- function() {
    as.integer(getOption("methylclock.bdsm_block", 2000L))
}

# Block size for a tall operand (a whole-array clock's ~10^5 CpGs x samples). A
# larger block cuts the number of reads along the long axis.
.mc_bdsm_block_tall <- function() {
    as.integer(getOption("methylclock.bdsm_block_tall", 10000L))
}

# Coefficient CpG count above which a clock is treated as whole-array and left to
# the streaming path rather than batched (it would pad the shared operand with
# mostly-zero rows). Matches the coverage split used when reading.
.mc_batch_max_markers <- function(nrow_betas) {
    .MC_SEQ_COVERAGE * nrow_betas
}

# Configure BigDataStatMeth's HDF5 algebra for this call. crossprod reads the
# block size from the global options only (not per-call), and its C++ asks the
# OpenMP runtime for the thread count on every call, so both are set here.
.mc_bdsm_set_algebra <- function(block = .mc_bdsm_block(),
                                 threads = .mc_bdsm_threads()) {
    .mc_omp_set_threads(threads)
    BigDataStatMeth::hdf5matrix_options(paral = TRUE, threads = threads,
        block_size = block, compression = 0L)
}

# Per-clock pieces for the linear engine, reusing the shared coefficient helpers:
# CpG markers and weights (intercept split off), the output transform, and how
# many of the clock's CpGs are present in the source.
.mc_linear_pieces <- function(entry, rn) {
    coef <- mcd_resource(entry$resource)
    cpg_col <- .coef_cpg_col(coef)
    w_col <- .coef_weight_col(coef, entry$coef_col)
    markers <- as.character(coef[[cpg_col]])
    weights <- as.numeric(coef[[w_col]])
    ic <- which(markers %in% c("(Intercept)", "Intercept", "intercept"))
    b0 <- 0
    if (isTRUE(entry$intercept)) {
        if (!length(ic)) ic <- 1L
        b0 <- weights[ic[1]]
        markers <- markers[-ic[1]]
        weights <- weights[-ic[1]]
    } else if (length(ic)) {
        markers <- markers[-ic[1]]
        weights <- weights[-ic[1]]
    }
    present <- markers %in% rn
    list(name = entry$name, markers = markers[present],
        weights = weights[present], b0 = b0, transform = entry$transform,
        coverage = mean(present), n_absent = sum(!present),
        n_total = length(markers))
}

# Write an in-memory matrix to an HDF5 dataset without the extra copy that
# creating with data= makes: create empty, then write column blocks.
.mc_stage_operand <- function(m, file, ds, block = 4096L) {
    BigDataStatMeth::hdf5_create_matrix(file, ds, nrow = nrow(m),
        ncol = ncol(m), overwrite = TRUE, compression = 0L)
    h <- BigDataStatMeth::hdf5_matrix(file, ds)
    for (s in seq(1L, ncol(m), by = block)) {
        e <- min(s + block - 1L, ncol(m))
        h[, s:e] <- m[, s:e, drop = FALSE]
    }
    .mc_close(h)
    invisible(file)
}

# Estimate a set of linear clocks that share one beta source in a single on-disk
# cross-product. `src` is a beta source (already scaled if the group
# standardizes); `entries` are its linear clock registry rows. Returns a named
# list of per-sample vectors (an NA vector for a clock below the coverage
# threshold). Coverage and missing-CpG warnings match the per-clock engine.
.predict_linear_batch_bdsm <- function(entries, src, min.perc = 0.8,
                                       block = .mc_bdsm_block(),
                                       threads = .mc_bdsm_threads()) {
    rn <- .src_rownames(src)
    samples <- .src_colnames(src)
    pieces <- lapply(entries, .mc_linear_pieces, rn = rn)
    names(pieces) <- vapply(entries, function(e) e$name, "")

    out <- vector("list", length(pieces))
    names(out) <- names(pieces)

    # coverage rule per clock: at or below the threshold returns NA (with a
    # warning), exactly as the per-clock engine; those clocks skip the matmul.
    keep <- logical(length(pieces))
    for (i in seq_along(pieces)) {
        p <- pieces[[i]]
        if (p$coverage <= min.perc) {
            warning(sprintf(
                "clock '%s': %.0f%% of CpGs present (need > %.0f%%); ",
                p$name, 100 * p$coverage, 100 * min.perc),
                "returning NA.", call. = FALSE)
            na <- rep(NA_real_, length(samples))
            names(na) <- samples
            out[[i]] <- na
        } else {
            if (p$n_absent > 0L) {
                warning(sprintf(
                    "clock '%s': %d of %d CpGs absent; computed on the rest.",
                    p$name, p$n_absent, p$n_total), call. = FALSE)
            }
            keep[i] <- TRUE
        }
    }
    kept <- pieces[keep]
    if (!length(kept)) {
        return(out)
    }

    # shared operand: the union of the kept clocks' present CpGs, read once,
    # missing values imputed to each CpG's mean across samples.
    union_cpg <- unique(unlist(lapply(kept, `[[`, "markers")))
    B <- .impute_rows(.src_subset(src, union_cpg))
    rownames(B) <- union_cpg

    # coefficient matrix: one column per kept clock, its weights on the union
    # rows, zero elsewhere.
    coef <- matrix(0, nrow = length(union_cpg), ncol = length(kept),
        dimnames = list(union_cpg, names(kept)))
    for (j in seq_along(kept)) {
        coef[kept[[j]]$markers, j] <- kept[[j]]$weights
    }

    # single on-disk cross-product: t(coef) %*% B = kept clocks x samples.
    file <- tempfile(fileext = ".h5")
    on.exit(unlink(file), add = TRUE)
    .mc_stage_operand(B, file, "batch/betas")
    .mc_stage_operand(coef, file, "batch/coef")
    rm(B, coef)
    .mc_bdsm_set_algebra(block = block, threads = threads)
    Bh <- BigDataStatMeth::hdf5_matrix(file, "batch/betas")
    Ch <- BigDataStatMeth::hdf5_matrix(file, "batch/coef")
    res <- BigDataStatMeth::crossprod(Ch, Bh, outgroup = "batch",
        outdataset = "out")
    est <- as.matrix(res)
    .mc_close(res)
    .mc_close(Bh)
    .mc_close(Ch)

    # per clock: intercept, then the output transform, on its row of the result.
    for (j in seq_along(kept)) {
        p <- kept[[j]]
        v <- p$transform(est[j, ] + p$b0)
        names(v) <- samples
        out[[p$name]] <- v
    }
    out
}

# A whole-array standardized linear clock (e.g. BLUP) estimated on disk WITHOUT
# materialising a scaled+imputed operand. The per-sample z-score is affine, so
# pred_j = (w . x_j - center_j * sum_w) / scale_j: the cross-product runs on the
# RAW input directly (the coefficient vector, padded over the input's CpGs, is
# staged in a temporary file -- the input is never written to), and the per-sample
# correction is applied to the small 1 x samples result. Exact when there are no
# missing values (imputation would then be a no-op); otherwise this returns NULL
# so the caller falls back to the streaming path, which imputes.
.predict_linear_folded_bdsm <- function(entry, betas, betas_z, min.perc = 0.8,
                                        block = .mc_bdsm_block_tall(),
                                        threads = .mc_bdsm_threads()) {
    # the fold needs a cross-file cross-product (coef in a temp file, betas in the
    # input file); older BigDataStatMeth cannot do it -> decline to the streaming
    # path silently.
    if (utils::packageVersion("BigDataStatMeth") < "2.0.4")
        return(NULL)
    if (!isTRUE(entry$standardize))
        return(NULL)                                 # fold needs per-sample stats
    scl <- .src_scaling(betas_z)
    if (is.null(scl) || !.src_complete(betas_z))
        return(NULL)                                 # missing values -> fall back

    rn <- .src_rownames(betas)
    samples <- .src_colnames(betas)
    p <- .mc_linear_pieces(entry, rn)
    if (p$coverage <= min.perc) {
        warning(sprintf(
            "clock '%s': %.0f%% of CpGs present (need > %.0f%%); ",
            p$name, 100 * p$coverage, 100 * min.perc),
            "returning NA.", call. = FALSE)
        out <- rep(NA_real_, length(samples))
        names(out) <- samples
        return(out)
    }
    if (p$n_absent > 0L)
        warning(sprintf(
            "clock '%s': %d of %d CpGs absent; computed on the rest.",
            p$name, p$n_absent, p$n_total), call. = FALSE)

    sum_w <- sum(p$weights)
    wvec <- numeric(length(rn))
    names(wvec) <- rn
    wvec[p$markers] <- p$weights

    file <- tempfile(fileext = ".h5")
    on.exit(unlink(file), add = TRUE)
    BigDataStatMeth::hdf5_create_matrix(file, "coef/w",
        data = matrix(wvec, ncol = 1), overwrite = TRUE)
    .mc_bdsm_set_algebra(block = block, threads = threads)
    Ch <- BigDataStatMeth::hdf5_matrix(file, "coef/w")
    Bh <- .src_base(betas)$hm                        # raw input, not written to
    res <- BigDataStatMeth::crossprod(Ch, Bh)        # cross-file; result in coef file
    raw <- as.numeric(as.matrix(res))
    .mc_close(res)
    .mc_close(Ch)

    pred <- (raw - scl$center * sum_w) / scl$scale
    out <- p$transform(pred + p$b0)
    names(out) <- samples
    out
}

# Split a linear selection into the groups that batch on disk (sparse clocks that
# share one preprocessing), the whole-array clocks (their own folded cross-product)
# and the names that stay on the per-clock path (quantile, non-linear). Returns a
# list with `groups` (each: entries + their source), `wholearray` (entries) and
# `passthrough` (clock names).
.mc_batch_plan <- function(sel, betas, betas_z) {
    max_markers <- .mc_batch_max_markers(.src_nrow(betas))
    raw <- list()
    zsc <- list()
    whole <- list()
    pass <- character(0)
    for (nm in names(sel)) {
        e <- sel[[nm]]
        quant <- identical(e$normalize, "quantile_gs")
        if (!identical(e$predictor, "linear") || quant) {
            pass <- c(pass, nm)
            next
        }
        # a table with squared-beta terms (power column) is not a plain
        # cross-product: those clocks take the per-clock path, which squares
        # the marked rows
        if ("power" %in% names(mcd_resource(e$resource))) {
            pass <- c(pass, nm)
            next
        }
        n_markers <- nrow(mcd_resource(e$resource)) - isTRUE(e$intercept)
        if (n_markers > max_markers) {           # whole-array: folded path
            whole[[nm]] <- e
            next
        }
        if (isTRUE(e$standardize)) {
            zsc[[nm]] <- e
        } else {
            raw[[nm]] <- e
        }
    }
    groups <- list()
    if (length(raw)) {
        groups <- c(groups, list(list(entries = raw, src = betas)))
    }
    if (length(zsc)) {
        groups <- c(groups, list(list(entries = zsc, src = betas_z)))
    }
    list(groups = groups, wholearray = whole, passthrough = pass)
}
