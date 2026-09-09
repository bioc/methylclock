# A backend-agnostic view over a methylation matrix. The engine asks a source
# for CpG/sample names and for the rows of a given set of CpGs; the source
# serves them from memory or reads just those rows from an HDF5 file, so a clock
# over a few hundred CpGs never loads a full array. A "scaled" source also
# applies a per-sample z-score on the fly, so standardised clocks stay
# out-of-core too.
#
# INTERNAL: the out-of-core input path (D1b/D17/D18). A clock's rows are read by
# INTERNAL: name through .src_subset, in row blocks. An HDF5 input is NEVER
# INTERNAL: loaded whole: z-score stats stream in blocks (.hdf5_col_stats) and
# INTERNAL: each clock reads only its CpG rows, a block at a time. A dense input
# INTERNAL: is worked in RAM by subsets (its rows are the caller's object -- we
# INTERNAL: cannot free it, so we never spill it while it is still resident). We
# INTERNAL: never hold a disk copy and a RAM copy of the same data at once.

# Wrap an already-dense matrix as an in-memory source.
.as_source_from_matrix <- function(m) {
    list(kind = "memory", data = m, rn = rownames(m), cn = colnames(m))
}

# Above this dense size, note that the matrix stays in RAM and point to the
# on-disk route (we cannot free an argument matrix -- the caller holds it).
.MC_DENSE_NOTE_BYTES <- 1024^3

# Turn user input into a beta source. An HDF5-backed matrix stays on disk; any
# other input is coerced to a dense matrix.
.as_beta_source <- function(x) {
    x <- .extract_betas(x)
    if (inherits(x, "HDF5Matrix")) {
        rn <- rownames(x)
        if (is.null(rn))
            stop("The HDF5 methylation matrix needs CpG identifiers as row ",
                 "names.")
        return(list(kind = "hdf5", hm = x, rn = rn, cn = colnames(x)))
    }
    m <- .betas_from_mvalues(.as_beta_matrix(x))
    bytes <- as.numeric(nrow(m)) * ncol(m) * 8
    if (bytes > .MC_DENSE_NOTE_BYTES)
        message("methylclock: large in-memory matrix (",
                round(bytes / 1024^3, 1), " GB). Clocks run block by block ",
                "but it stays in RAM. For lower memory, store it on disk once ",
                "with mc_to_hdf5() and pass that instead.")
    .as_source_from_matrix(m)
}

# A source that applies a per-column (per-sample) center/scale to whatever its
# base source returns.
.as_scaled_source <- function(base, center, scale, complete = TRUE) {
    list(kind = "scaled", base = base, center = center, scale = scale,
         complete = complete,
         rn = .src_rownames(base), cn = .src_colnames(base))
}

# Did the scaled view's stats see a complete (no missing value) input? The folded
# whole-array path is exact only then; otherwise it declines to the streaming
# path, which imputes.
.src_complete <- function(s) !identical(s$kind, "scaled") || isTRUE(s$complete)

.src_rownames <- function(s) s$rn
.src_colnames <- function(s) s$cn
.src_ncol <- function(s) length(s$cn)
.src_nrow <- function(s) length(s$rn)

# Is the source backed on disk (directly or through a scaled view)? On an HDF5
# backend a scattered row read can hit many chunks; a whole-array clock then
# reads it in file order instead (see .predict_linear).
.src_is_hdf5 <- function(s) {
    identical(s$kind, "hdf5") ||
        (identical(s$kind, "scaled") && identical(s$base$kind, "hdf5"))
}

# The unscaled base source, and the (center, scale) of a scaled view (or NULL).
.src_base <- function(s) if (identical(s$kind, "scaled")) s$base else s
.src_scaling <- function(s) {
    if (identical(s$kind, "scaled"))
        list(center = s$center, scale = s$scale)
    else
        NULL
}

# A contiguous range of rows (by position) as a dense matrix. Reads in file
# order, so on an HDF5 backend the read is contiguous (no chunk amplification).
.src_range <- function(s, from, to) {
    switch(s$kind,
        memory = s$data[from:to, , drop = FALSE],
        hdf5 = {
            m <- as.matrix(s$hm[from:to, , drop = FALSE])
            dimnames(m) <- list(s$rn[from:to], s$cn)
            m
        },
        scaled = .apply_scale(.src_range(s$base, from, to), s$center, s$scale),
        stop("Unknown beta source."))
}

# Apply a stored per-column center/scale to a CpGs x samples block, in place
# column by column so no full-size temporary is made (sweep() would copy the
# whole block 2-3 times, which matters for clocks that use most of the array).
.apply_scale <- function(B, center, scale) {
    for (j in seq_len(ncol(B)))
        B[, j] <- (B[, j] - center[j]) / scale[j]
    B
}

# Rows (CpGs) of a source as a dense matrix, loading only those rows. All `cpgs`
# must be present in the source.
.src_subset <- function(s, cpgs) {
    switch(s$kind,
        memory = as.matrix(s$data[cpgs, , drop = FALSE]),
        hdf5 = {
            m <- as.matrix(s$hm[match(cpgs, s$rn), , drop = FALSE])
            dimnames(m) <- list(cpgs, s$cn)
            m
        },
        scaled = .apply_scale(.src_subset(s$base, cpgs), s$center, s$scale),
        stop("Unknown beta source."))
}

# The whole source as a dense matrix; materialises an HDF5 source. Used only
# where a full-matrix pass is unavoidable.
.src_dense <- function(s) {
    switch(s$kind,
        memory = s$data,
        hdf5 = {
            m <- as.matrix(s$hm, force = TRUE)
            if (is.null(rownames(m))) rownames(m) <- s$rn
            if (is.null(colnames(m))) colnames(m) <- s$cn
            m
        },
        scaled = .apply_scale(.src_dense(s$base), s$center, s$scale),
        stop("Unknown beta source."))
}

# Rows (CpGs) processed at a time when a step must sweep many/all of them (disk
# statistics, block-wise prediction), so peak RAM stays bounded to one block.
# Only whole-array clocks span many blocks; a smaller block lowers their peak at
# the cost of more reads (typical clocks fit in a single block either way).
.MC_ROW_BLOCK <- 10000L

# Per-sample (per-column) mean and sd of an HDF5 matrix, computed in row blocks
# so peak RAM stays bounded (one block, not the whole array). Matches scale()'s
# sample sd (divisor n - 1).
# INTERNAL: D18 -- BDSM colMeans/colSds read the full matrix into RAM, which
# INTERNAL: defeats the memory bound; this streams instead.
.hdf5_col_stats <- function(hm, block = .MC_ROW_BLOCK) {
    n <- nrow(hm)
    p <- ncol(hm)
    s1 <- numeric(p)
    s2 <- numeric(p)
    cnt <- numeric(p)
    i <- 1L
    while (i <= n) {
        j <- min(i + block - 1L, n)
        B <- as.matrix(hm[i:j, , drop = FALSE])
        ok <- is.finite(B)
        B[!ok] <- 0                                # NAs drop out of the sums
        s1 <- s1 + colSums(B)
        s2 <- s2 + colSums(B * B)
        cnt <- cnt + colSums(ok)
        i <- j + 1L
    }
    center <- ifelse(cnt > 0, s1 / cnt, 0)
    variance <- ifelse(cnt > 1, (s2 - cnt * center^2) / (cnt - 1), 1)
    variance[!is.finite(variance) | variance <= 0] <- 1
    list(center = center, scale = sqrt(variance), complete = all(cnt == n))
}

# Per-sample mean and sd of a dense matrix, one column at a time so no full-size
# copy is made (x^2 or scale() would copy the whole matrix). Sample sd (n - 1).
.dense_col_stats <- function(x) {
    p <- ncol(x)
    n <- nrow(x)
    center <- numeric(p)
    scale <- numeric(p)
    complete <- TRUE
    for (j in seq_len(p)) {
        cj <- x[, j]
        ok <- is.finite(cj)
        nj <- sum(ok)
        if (nj < n) complete <- FALSE
        m <- if (nj) mean(cj[ok]) else 0
        s <- if (nj > 1) sqrt(sum((cj[ok] - m)^2) / (nj - 1)) else 1
        center[j] <- m
        scale[j] <- if (is.finite(s) && s > 0) s else 1
    }
    list(center = center, scale = scale, complete = complete)
}

# Per-sample z-scored view of a source: compute the per-sample mean/sd once
# (cheaply, without a whole-array copy) and apply them lazily on each subset, so
# a full standardised matrix is never materialised -- standardising clocks add
# no whole-array copy, on any backend.
# INTERNAL: D18. Replaces the earlier spill-to-disk approach, which a memory
# INTERNAL: stress test showed was counterproductive on dense input (it kept the
# INTERNAL: resident matrix AND added disk I/O). Lazy scaling protects memory by
# INTERNAL: construction, with no threshold to tune.
.zscore_source <- function(s) {
    st <- if (identical(s$kind, "hdf5")) .hdf5_col_stats(s$hm)
          else .dense_col_stats(.src_dense(s))
    .as_scaled_source(s, st$center, st$scale, complete = st$complete)
}
