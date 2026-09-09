# K-nearest-neighbour imputation of a methylation matrix.
#
# The matrix is CpGs (rows) by samples (columns). A missing value in row i is
# filled from the k CpGs most similar to i -- by averaged Euclidean distance
# over their co-observed samples -- as the unweighted mean of those neighbours
# observed at the missing column, so co-methylation is respected rather than
# the CpG's marginal mean. This matches the classic gene-expression KNN imputer.
#
# The numerical work runs entirely in the C++ layer (src/impute_knn.cpp): the
# exact kernel, and -- for a row set larger than `maxp` -- the block split with
# row clustering from a BigDataStatMeth SVD embedding. Small row sets are
# imputed exactly; a large one is split into blocks of similar CpGs and each
# block is imputed exactly, so a big matrix never forms a full CpG-by-CpG
# distance matrix. Only the clock's CpGs are read from a disk-backed input, so
# the array is never loaded whole.

# Defaults match the reference KNN imputer so results line up with it.
.MC_KNN <- list(k = 10L, maxp = 1500L, rowmax = 0.5, colmax = 0.8)

# Normalise user-supplied KNN options against the defaults.
.knn_opts <- function(k = NULL, maxp = NULL, rowmax = NULL, colmax = NULL) {
    o <- .MC_KNN
    if (!is.null(k)) o$k <- as.integer(k)
    if (!is.null(maxp)) o$maxp <- as.integer(maxp)
    if (!is.null(rowmax)) o$rowmax <- as.numeric(rowmax)
    if (!is.null(colmax)) o$colmax <- as.numeric(colmax)
    o
}

#' KNN-impute a methylation matrix
#'
#' Fills the missing values of a methylation matrix from each CpG's nearest
#' CpGs -- those that co-vary with it across the observed samples -- rather than
#' from its marginal mean, so co-methylation is respected. A row set within
#' \code{maxp} is imputed exactly; a larger one is first split into blocks of
#' similar CpGs (through a BigDataStatMeth SVD embedding) and each block is
#' imputed exactly. This is the same imputation \code{\link{compute_clocks}}
#' applies with \code{impute = "knn"}; use it to impute a matrix once and reuse
#' it across analyses.
#'
#' @param x A methylation matrix with CpG identifiers as row names and samples
#'   as columns (a data frame with a CpG-identifier column is also accepted).
#' @param k Number of neighbouring CpGs averaged for each imputed value.
#' @param maxp Largest block of CpGs imputed exactly; larger row sets are split
#'   into similar blocks first.
#' @param rowmax If a CpG has a greater fraction of samples missing than this,
#'   it is filled with the per-sample mean instead of by KNN.
#' @param colmax If a sample has a greater fraction of CpGs missing than this, a
#'   warning is raised (the sample is very sparse); imputation still proceeds.
#' @return The matrix with missing values filled, same shape and dimnames.
#' @seealso \code{\link{compute_clocks}}
#' @examples
#' data(methylclock_betas)
#' m <- methylclock_betas[seq_len(200), ]
#' m[5, 2] <- NA
#' m[40, c(1, 3)] <- NA
#' filled <- imputeKNN(m)
#' filled[c(5, 40), seq_len(3)]
#' @export
imputeKNN <- function(x, k = 10, maxp = 1500, rowmax = 0.5, colmax = 0.8) {
    m <- .betas_from_mvalues(.as_beta_matrix(.extract_betas(x)))
    if (!anyNA(m)) return(m)
    storage.mode(m) <- "double"
    mc_knn_impute_mem(m, k = as.integer(k), maxp = as.integer(maxp),
                      rowmax = rowmax, colmax = colmax)
}

# One out-of-core KNN block is held in RAM at a time; its row count is derived
# from the sample count against this byte budget, so peak memory stays flat
# whatever the array's shape. Override with options(methylclock.knn_block_mb=).
.MC_KNN_BLOCK_MB <- 512
# Bytes held per block beyond the raw block itself (its in-kernel copy and the
# embedding's zero-filled copy), as a multiple of the block size.
.MC_KNN_BLOCK_SAFETY <- 3L

# Rows per out-of-core KNN block for `nsamp` samples: the byte budget divided by
# one block-row's cost, floored at `maxp` (so a block still clusters) and capped
# at the array height. Flat memory by construction: a block always costs about
# the budget, whatever the total number of CpGs.
.knn_block_rows <- function(nsamp, maxp, nrow_total) {
    mb <- getOption("methylclock.knn_block_mb", .MC_KNN_BLOCK_MB)
    budget <- as.numeric(mb) * 1024^2
    b <- floor(budget / (.MC_KNN_BLOCK_SAFETY * nsamp * 8))
    as.integer(max(maxp, min(b, nrow_total)))
}

# Weighted sum of a whole-array clock after KNN imputation, out-of-core. Sweeps
# the array in file order in RAM-bounded blocks (their row count derived from the
# sample count), reading each contiguous block once; keeps the clock's CpGs it
# holds, KNN-imputes those raw rows in the C++ engine, applies the source's
# per-sample scaling, and accumulates the weighted sum. Peak RAM is one block, so
# a clock spanning the array is never loaded whole. Neighbours are drawn from
# within a block rather than the entire array, matching how the in-memory engine
# splits a large row set into blocks. Returns the weighted sum and per-sample
# missing counts, like the other .predict_* helpers.
.predict_knn_seq <- function(betas, m, w, nsamp, opts = .MC_KNN) {
    wmap <- tapply(w, m, sum)                       # weight per unique CpG
    base <- .src_base(betas)
    scl <- .src_scaling(betas)
    rn <- .src_rownames(betas)
    n <- length(rn)
    want <- rn %in% names(wmap)
    step <- .knn_block_rows(nsamp, opts$maxp, n)
    pred <- numeric(nsamp)
    na <- numeric(nsamp)
    for (start in seq(1L, n, by = step)) {
        end <- min(start + step - 1L, n)
        hit <- which(want[start:end])
        if (!length(hit)) next
        B <- .src_range(base, start, end)[hit, , drop = FALSE]
        na <- na + colSums(is.na(B))
        if (anyNA(B)) {
            storage.mode(B) <- "double"
            B <- mc_knn_impute_mem(B, k = opts$k, maxp = opts$maxp,
                                   rowmax = opts$rowmax, colmax = opts$colmax)
        }
        if (!is.null(scl)) B <- .apply_scale(B, scl$center, scl$scale)
        pred <- pred + as.numeric(crossprod(wmap[rownames(B)], B))
    }
    list(pred = pred, na = na)
}

# Weighted sum of a clock's CpGs after KNN imputation (used by .predict_linear
# when impute == "knn"). Reads the clock's raw rows from the beta source -- only
# those rows, so a disk-backed array is not loaded whole -- imputes them in the
# C++ engine, applies the source's per-sample scaling, then forms the weighted
# sum. Imputing the raw betas before standardising matches how the input is
# imputed once, up front. Returns the weighted sum and per-sample missing
# counts, like the other .predict_* helpers.
.predict_knn <- function(betas, m, w, nsamp, opts = .MC_KNN,
                         square = FALSE) {
    base <- .src_base(betas)
    scl <- .src_scaling(betas)
    B <- .src_subset(base, m)
    na <- colSums(is.na(B))
    storage.mode(B) <- "double"
    if (anyNA(B))
        B <- mc_knn_impute_mem(B, k = opts$k, maxp = opts$maxp,
                               rowmax = opts$rowmax, colmax = opts$colmax)
    if (!is.null(scl)) B <- .apply_scale(B, scl$center, scl$scale)
    # a squared-term block enters the weighted sum as beta^2
    if (isTRUE(square)) B <- B * B
    list(pred = as.numeric(crossprod(w, B)), na = na)
}
