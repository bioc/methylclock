# Promote a result to a persistent on-disk file.
#
# INTERNAL: D16. A result already lives in an HDF5 file (D12); when that file is
# INTERNAL: a session temporary, `persist()` moves it to a user path -- a file
# INTERNAL: copy, not a recomputation. A memory-only result (disk fallback, D14)
# INTERNAL: is written out fresh. Named `persist` (not `save`/`export`): `save`
# INTERNAL: is not a generic and `export` clashes with rtracklayer/rio.

#' Persist a result to a file
#'
#' Save an object to a file so it outlives the current session.
#'
#' @param x The object to persist.
#' @param file Destination file path.
#' @param ... Further arguments passed to methods.
#' @return The object, updated to reference the persistent file (invisibly).
#' @export
persist <- function(x, file, ...) UseMethod("persist")

#' @rdname persist
#' @param overwrite Overwrite \code{file} if it already exists. Default
#'   \code{FALSE}.
#' @details For a \code{"methylclock"} result the estimates already live in an
#'   HDF5 file; when that file is a session temporary, persisting copies it to
#'   \code{file} (no recomputation) and verifies the copy. A result held only in
#'   memory is written out to \code{file}. The returned object references the
#'   persistent file.
#' @examples
#' data(methylclock_betas)
#' res <- methylclock(methylclock_betas, clocks = "Horvath")
#' f <- tempfile(fileext = ".h5")
#' res <- persist(res, f)
#' file.exists(f)
#' unlink(f)
#' @export
persist.methylclock <- function(x, file, overwrite = FALSE, ...) {
    if (!is.character(file) || length(file) != 1L)
        stop("`file` must be a single file path.")
    if (file.exists(file) && !overwrite)
        stop("`", file, "` already exists; pass overwrite = TRUE to replace.")

    if (identical(x$storage, "hdf5") && !is.null(x$handle)) {
        if (identical(normalizePath(x$handle$file, mustWork = FALSE),
                      normalizePath(file, mustWork = FALSE))) {
            x$persistent <- TRUE
            return(invisible(x))
        }
        ok <- file.copy(x$handle$file, file, overwrite = overwrite)
        if (!isTRUE(ok))
            stop("Could not copy the result file to '", file, "'.")
        h <- x$handle
        path <- paste(h$group, h$dataset, sep = "/")
        back <- .mc_hdf5_read_matrix(file, path)   # verify the copy
        # INTERNAL: D15 (verify after copy).
        if (!identical(dim(back), c(h$nrow, h$ncol)) ||
            !identical(colnames(back), h$clocks))
            stop("The persisted copy '", file, "' does not match the result.")
        x$handle$file <- file
    } else {
        x$handle <- .mc_hdf5_write(as.data.frame(x), file = file)
        x$storage <- "hdf5"
    }
    x$persistent <- TRUE
    message("methylclock: result persisted to '", file, "'.")
    invisible(x)
}
