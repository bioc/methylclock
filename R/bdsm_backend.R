# On-disk (HDF5) backend for a methylclock result, via BigDataStatMeth. HDF5 is
# the source of truth for every result; a copy is kept in memory as a read cache
# only when it fits comfortably. The result table (samples x clocks) is small,
# so this is about durability and consistency with the rest of the stack.
#
# INTERNAL: D7/D12-D16. BigDataStatMeth is reached ONLY through this file, so
# INTERNAL: later phases (out-of-core input) have one call site. Calls are
# INTERNAL: qualified `BigDataStatMeth::` to avoid masking (cor/var/scale...).
# INTERNAL: Never rhdf5. Durability (D15): open/close in the same scope, never
# INTERNAL: store a live handle, verify every write (read-after-write), never
# INTERNAL: silently swallow a close failure.

# Where a result's values live inside its HDF5 file.
.MC_H5_GROUP <- "methylclock"
.MC_H5_DATASET <- "values"

# Keep a result in the in-memory read cache unless its values would occupy more
# than this in RAM. Only trips on very large sample x clock counts.
# INTERNAL: D12 (memory cache threshold).
.MC_CACHE_BYTES <- 256 * 1024^2

# Compression by size: small results uncompressed (faster, no overhead); larger
# ones use level 3 (a good compute/retrieval vs. disk trade-off in BDSM).
# INTERNAL: D13.
.MC_COMPRESS_BYTES <- 16 * 1024^2

# Is BigDataStatMeth usable? It is an Import, so normally yes; guarded so an
# on-disk request fails with a clear message rather than an obscure one.
.bdsm_available <- function() {
    requireNamespace("BigDataStatMeth", quietly = TRUE)
}

.bdsm_require <- function() {
    if (!.bdsm_available())
        stop("On-disk (hdf5) storage needs the 'BigDataStatMeth' package.")
}

# Rough in-memory size (bytes) of a result's numeric values.
.mc_values_bytes <- function(values) {
    k <- length(setdiff(names(values), "id"))
    as.numeric(nrow(values)) * k * 8
}

# Split a values data frame (leading `id` column + numeric clock columns) into
# the numeric matrix (samples x clocks) that goes to HDF5, with sample ids as
# row names.
.mc_values_to_matrix <- function(values) {
    num <- values[, setdiff(names(values), "id"), drop = FALSE]
    m <- as.matrix(num)
    storage.mode(m) <- "double"
    rownames(m) <- as.character(values$id)
    m
}

# Close an HDF5 handle, warning (never silently) if it does not close cleanly.
# The global net BigDataStatMeth::hdf5_close_all() is a last resort.
# INTERNAL: D15 (never swallow a close failure).
.mc_close <- function(h) {
    err <- tryCatch({ close(h); NULL }, error = function(e) conditionMessage(e))
    if (!is.null(err))
        warning("methylclock: an HDF5 handle did not close cleanly: ", err)
    invisible(is.null(err))
}

# Read a result's values back into a matrix (samples x clocks), closing the read
# handle in the same scope.
.mc_hdf5_read_matrix <- function(file, path) {
    hm <- BigDataStatMeth::hdf5_matrix(file, path)
    on.exit(.mc_close(hm), add = TRUE)
    as.matrix(hm, force = TRUE)
}

# Compression level for a matrix, by size; NULL means "choose by size".
.mc_compression <- function(m, compression = NULL) {
    if (!is.null(compression)) return(compression)
    if (as.numeric(length(m)) * 8 > .MC_COMPRESS_BYTES) 3L else 0L
}

# Write a result's values to an HDF5 file and return a handle descriptor (text
# only -- never a live handle). Verifies the write by reading it back before
# returning: once this returns, the data are on disk, flushed and checked.
# INTERNAL: D15 (read-after-write verification).
.mc_hdf5_write <- function(values, file, compression = NULL,
                           group = .MC_H5_GROUP, dataset = .MC_H5_DATASET) {
    .bdsm_require()
    m <- .mc_values_to_matrix(values)
    comp <- .mc_compression(m, compression)
    path <- paste(group, dataset, sep = "/")
    h <- BigDataStatMeth::hdf5_create_matrix(
        filename = file, dataset = path, data = m,
        overwrite = TRUE, compression = comp)
    .mc_close(h)                                    # flush + release the writer

    back <- .mc_hdf5_read_matrix(file, path)        # read-after-write check
    ok <- identical(dim(back), dim(m)) &&
        identical(dimnames(back), dimnames(m)) &&
        isTRUE(all.equal(back, m, tolerance = 0))
    if (!ok)
        stop("HDF5 write verification failed for '", file, "'.")

    list(file = file, group = group, dataset = dataset,
         clocks = colnames(m), nrow = nrow(m), ncol = ncol(m),
         compression = comp)
}

# Read a result's values back from its HDF5 handle into a data frame, restoring
# the leading `id` column from the stored row names.
.mc_hdf5_read <- function(handle) {
    .bdsm_require()
    path <- paste(handle$group, handle$dataset, sep = "/")
    m <- .mc_hdf5_read_matrix(handle$file, path)
    df <- data.frame(id = rownames(m), stringsAsFactors = FALSE)
    for (cl in colnames(m)) df[[cl]] <- as.numeric(m[, cl])
    rownames(df) <- NULL
    df
}

# A unique HDF5 file name (no timestamp/random needed: tempfile gives a fresh
# token even under a user-chosen directory).
.mc_h5_name <- function() {
    paste0("methylclock_", basename(tempfile("")), ".h5")
}

# Resolve where a result's HDF5 file goes and whether it persists:
#   file given            -> that path, persistent
#   option hdf5_dir set   -> a file in that directory, persistent (user consent)
#   otherwise             -> a session temp file, NOT persistent (announced)
# Writing outside the session tempdir requires explicit user input, per the
# Bioconductor policy against writing to the working directory by default.
# INTERNAL: D16.
.mc_resolve_file <- function(file) {
    if (!is.null(file))
        return(list(path = file, persistent = TRUE))
    dir <- getOption("methylclock.hdf5_dir", default = NULL)
    if (!is.null(dir)) {
        if (!dir.exists(dir))
            dir.create(dir, recursive = TRUE, showWarnings = FALSE)
        return(list(path = file.path(dir, .mc_h5_name()), persistent = TRUE))
    }
    list(path = tempfile(fileext = ".h5"), persistent = FALSE)
}

# INTERNAL: the storage policy (D12 source-of-truth + cache / D14 fallback / D16
# INTERNAL: location) is applied inline in compute_clocks, where results are
# INTERNAL: streamed to disk clock by clock (D18).
