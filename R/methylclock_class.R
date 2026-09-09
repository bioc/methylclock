# The `methylclock` result object (S3).
#
# INTERNAL: D7/D12. Backed on disk (HDF5, via BigDataStatMeth -- never rhdf5)
# INTERNAL: whenever possible; `values` holds an in-memory read cache only when
# INTERNAL: it fits. The `handle` slot is text only (file/group/dataset), never
# INTERNAL: a live connection, so the object is safe to keep and pass on.

#' Construct a methylclock result
#'
#' Low-level constructor for the object returned by the clock estimation
#' functions. Most users obtain one from \code{\link{methylclock}},
#' \code{\link{DNAmAge}} or \code{\link{DNAmGA}} rather than calling this
#' directly.
#'
#' @param values A data frame of estimates (one row per sample, one column per
#'   clock) when it is held in memory; \code{NULL} when the estimates live
#'   only on disk and are reached through \code{handle}.
#' @param clocks Character vector of the clocks that were computed.
#' @param storage Where the authoritative estimates live: \code{"memory"} or
#'   \code{"hdf5"}.
#' @param handle For \code{storage = "hdf5"}, a descriptor of the on-disk
#'   dataset (file, group, dataset); \code{NULL} otherwise.
#' @param persistent \code{TRUE} if the on-disk file is meant to outlive the
#'   session (a user-supplied path); \code{FALSE} for a session temporary file.
#' @param meta Named list of metadata (e.g. sample count, call, warnings).
#' @return An object of class \code{"methylclock"}.
#' @seealso \code{\link{methylclock}}
#' @examples
#' est <- data.frame(id = c("s1", "s2"), Horvath = c(41.2, 55.8))
#' res <- new_methylclock(est, clocks = "Horvath")
#' as.data.frame(res)
#' @export
new_methylclock <- function(values, clocks, storage = "memory",
                            handle = NULL, persistent = FALSE, meta = list()) {
    storage <- match.arg(storage, c("memory", "hdf5"))
    structure(
        list(values = values, clocks = clocks, storage = storage,
             handle = handle, persistent = persistent, meta = meta),
        class = "methylclock"
    )
}

#' Coerce a methylclock result to a data frame
#'
#' Returns the estimates as a data frame, materialising them from disk first if
#' they are not already held in memory. The backend is invisible to the caller.
#'
#' @param x A \code{"methylclock"} object.
#' @param ... Unused.
#' @return A data frame of estimates.
#' @export
as.data.frame.methylclock <- function(x, ...) {
    if (!is.null(x$values))
        return(as.data.frame(x$values))            # in-memory (or cache)
    if (identical(x$storage, "hdf5"))
        return(.mc_hdf5_read(x$handle))            # materialise from disk
    stop("This methylclock result holds no estimates.")
}

#' @export
print.methylclock <- function(x, ...) {
    df <- if (!is.null(x$values)) x$values else NULL
    n <- if (!is.null(df)) nrow(df) else x$meta$n %||% NA
    where <- if (identical(x$storage, "hdf5")) {
        if (isTRUE(x$persistent)) paste0("hdf5: ", x$handle$file)
        else "hdf5, session file"
    } else {
        "memory"
    }
    cat(sprintf("<methylclock> %s sample(s) x %d clock(s) [%s]\n",
                as.character(n), length(x$clocks), where))
    cat("  clocks:", paste(x$clocks, collapse = ", "), "\n")
    if (identical(x$storage, "hdf5") && !isTRUE(x$persistent))
        cat("  session file (removed when R exits) -- to keep it: ",
            "persist(result, \"file.h5\")\n", sep = "")
    if (!is.null(df))
        print(utils::head(df))
    invisible(x)
}

#' @export
summary.methylclock <- function(object, ...) {
    df <- as.data.frame(object)
    num <- vapply(df, is.numeric, logical(1))
    summary(df[, num, drop = FALSE])
}
