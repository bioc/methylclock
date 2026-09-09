# The estimation orchestrator: select clocks, run each predictor, assemble the
# result. Shared by the user-facing entry points.

#' Estimate clocks from a methylation matrix
#'
#' Computes the selected clocks over a methylation beta matrix and returns a
#' \code{\link{new_methylclock}} result. This is the engine shared by
#' \code{\link{methylclock}}, \code{\link{DNAmAge}} and \code{\link{DNAmGA}}.
#'
#' @param x A methylation matrix with CpG identifiers as row names and samples
#'   as columns (a data frame with a CpG-identifier column is also accepted). A
#'   large matrix may be supplied on disk as an HDF5-backed matrix, in which
#'   case each clock reads only the CpGs it needs rather than the whole array.
#' @param clocks Either \code{"all"} or a character vector of clock names.
#' @param target Optional biological target to restrict to (e.g.
#'   \code{"gestational"}); \code{"!gestational"} selects everything else.
#' @param generation Optional generation to restrict to (e.g. \code{1}).
#' @param platform Optional native platform to restrict to (e.g. \code{"450K"}):
#'   keeps clocks whose native platform includes any of the given arrays.
#' @param min.perc Minimum fraction of a clock's CpGs that must be present. A
#'   clock is only computed when its coverage is strictly greater than this;
#'   otherwise it is returned as \code{NA}. Default \code{0.8}.
#' @param impute How missing values among the present CpGs are filled before
#'   the weighted sum. \code{"mean"} (default) fills each with that CpG's mean
#'   across samples; \code{"reference"} fills it with the CpG's reference level
#'   from the clock's own coefficients (mean fallback), staying closer to how
#'   the clock was trained; \code{"none"} fills nothing, so a sample missing any
#'   used CpG is returned \code{NA} instead of an imputed value; \code{"knn"}
#'   fills each from the CpG's nearest co-varying CpGs (see \code{\link{imputeKNN}}).
#' @param k Neighbours averaged when \code{impute = "knn"}. Default 10.
#' @param maxp Largest block imputed by exact KNN when \code{impute = "knn"};
#'   larger row sets are split by 2-means first. Default 1500.
#' @param min.perc.sample Per-sample coverage floor. For each clock, a sample
#'   with fewer than this fraction of the clock's present CpGs observed is
#'   returned \code{NA}, so a heavily imputed sample is not reported as if it
#'   were a real estimate. Default \code{0} (off).
#' @param storage How the estimates are held: \code{"auto"} (default) backs them
#'   on disk (HDF5, through BigDataStatMeth) and keeps an in-memory cache when
#'   they fit; \code{"hdf5"} backs them on disk without an in-memory cache;
#'   \code{"memory"} keeps them in memory only (they will not persist). When a
#'   disk write is not possible the estimates fall back to memory, with a
#'   warning.
#' @param file HDF5 file to write the estimates to. When \code{NULL} (default) a
#'   session temporary file is used (or the directory in
#'   \code{options(methylclock.hdf5_dir=)} if set); pass a path to keep the
#'   result, or promote a session result later with \code{\link{persist}}.
#' @param compression HDF5 compression level (0-9); \code{NULL} (default) picks
#'   it by size (none for small results, moderate for large ones).
#' @param ... Reserved for future inputs (e.g. age, sex, cell composition).
#' @return A \code{"methylclock"} object.
#' @seealso \code{\link{methylclock}}, \code{\link{persist}},
#'   \code{\link{clock_list}}
#' @examples
#' data(methylclock_betas)
#' res <- compute_clocks(methylclock_betas, clocks = c("Horvath", "Levine"),
#'                       storage = "memory")
#' head(as.data.frame(res))
#' @export
compute_clocks <- function(x, clocks = "all", target = NULL,
                           generation = NULL, platform = NULL,
                           min.perc = 0.8, impute = c("mean", "reference",
                           "none", "knn"), min.perc.sample = 0, k = 10,
                           maxp = 1500,
                           storage = c("auto", "hdf5", "memory"),
                           file = NULL, compression = NULL, ...) {
    storage <- match.arg(storage)
    impute <- match.arg(impute)
    knn <- .knn_opts(k = k, maxp = maxp)
    betas <- .as_beta_source(x)
    sel <- .select_clocks(clocks, target, generation, platform)
    if (length(sel) == 0L)
        stop("No clocks selected.")
    ids <- .src_colnames(betas)

    # Some clocks need per-sample standardized input; build it once as a lazily
    # scaled source (per-sample mean/sd computed without a whole-array copy,
    # applied per subset), so standardising clocks add no full matrix.
    betas_z <- NULL
    if (any(vapply(sel, function(e) isTRUE(e$standardize), logical(1))))
        betas_z <- .zscore_source(betas)

    # On-disk streaming target: as each clock finishes, the running result is
    # written to its HDF5 file, so a mid-run failure keeps what completed.
    loc <- if (identical(storage, "memory")) NULL else .mc_resolve_file(file)
    stream_ok <- !is.null(loc)
    handle <- NULL

    # A clock named explicitly must succeed; one pulled in by a filter is
    # skipped (with a note) if its engine is not ready or it cannot be computed.
    explicit <- !(length(clocks) == 1L && identical(clocks, "all"))
    warns <- character(0)
    failed <- character(0)
    cols <- list()

    # On a disk-backed input, estimate the batchable linear clocks (those sharing
    # a preprocessing, not whole-array) through BigDataStatMeth: one on-disk
    # cross-product per group. The rest stay on the per-clock path below. A group
    # that fails leaves its clocks to that path, so the result is unaffected.
    # The batched/folded algebra implements only the default fill (each missing
    # value -> the CpG's across-sample mean) with no per-sample floor. When
    # either is changed, every clock takes the per-clock path, which honours
    # both.
    batch_vals <- list()
    if (.src_is_hdf5(betas) && identical(impute, "mean") &&
        min.perc.sample <= 0) {
        plan <- .mc_batch_plan(sel, betas, betas_z)
        capture_warn <- function(w) {
            warns[[length(warns) + 1L]] <<- conditionMessage(w)
            invokeRestart("muffleWarning")
        }
        for (g in plan$groups) {
            got <- tryCatch(
                withCallingHandlers(
                    .predict_linear_batch_bdsm(g$entries, g$src,
                                               min.perc = min.perc),
                    warning = capture_warn),
                error = function(err) {
                    warning("Batched algebra failed (", conditionMessage(err),
                            "); computing these clocks individually.")
                    NULL
                })
            if (!is.null(got)) batch_vals <- c(batch_vals, got)
        }
        # whole-array clocks: one folded cross-product each. A NULL means the
        # folded path declined (missing values); the clock then falls to the
        # per-clock streaming path below.
        for (nm in names(plan$wholearray)) {
            got <- tryCatch(
                withCallingHandlers(
                    .predict_linear_folded_bdsm(plan$wholearray[[nm]], betas,
                                                betas_z, min.perc = min.perc),
                    warning = capture_warn),
                error = function(err) {
                    warning("Whole-array algebra failed (",
                            conditionMessage(err),
                            "); computing this clock individually.")
                    NULL
                })
            if (!is.null(got)) batch_vals[[nm]] <- got
        }
    }

    for (nm in names(sel)) {
        mat <- if (isTRUE(sel[[nm]]$standardize)) betas_z else betas
        val <- if (nm %in% names(batch_vals)) {
            batch_vals[[nm]]
        } else tryCatch(
            withCallingHandlers(
                .predict_clock(sel[[nm]], mat, min.perc = min.perc,
                               impute = impute,
                               min.perc.sample = min.perc.sample, knn = knn),
                warning = function(w) {
                    warns[[length(warns) + 1L]] <<- conditionMessage(w)
                    invokeRestart("muffleWarning")
                }),
            error = function(err) {
                if (explicit) stop(err)
                failed[[length(failed) + 1L]] <<-
                    sprintf("%s (%s)", nm, conditionMessage(err))
                NULL
            })
        if (is.null(val)) next
        cols[[nm]] <- val
        if (stream_ok) {
            partial <- data.frame(id = ids, stringsAsFactors = FALSE)
            for (c2 in names(cols)) partial[[c2]] <- cols[[c2]]
            handle <- tryCatch(
                .mc_hdf5_write(partial, file = loc$path,
                               compression = compression),
                error = function(e) {
                    warning("Could not write results to disk (",
                            conditionMessage(e),
                            "); keeping them in memory only.")
                    stream_ok <<- FALSE
                    NULL
                })
        }
    }
    if (length(failed))
        warning("Skipped ", length(failed), " clock(s): ",
                paste(failed, collapse = "; "))
    if (length(cols) == 0L)
        stop("None of the selected clocks could be computed.")

    values <- data.frame(id = ids, stringsAsFactors = FALSE)
    for (nm in names(cols)) values[[nm]] <- cols[[nm]]
    meta <- list(n = length(ids), call = match.call(),
                 warnings = warns, skipped = failed)

    if (!is.null(loc) && stream_ok && !is.null(handle)) {
        cache <- if (identical(storage, "auto") &&
                     .mc_values_bytes(values) <= .MC_CACHE_BYTES)
            values else NULL
        return(new_methylclock(values = cache, clocks = names(cols),
                               storage = "hdf5", handle = handle,
                               persistent = loc$persistent, meta = meta))
    }
    new_methylclock(values = values, clocks = names(cols), storage = "memory",
                    handle = NULL, persistent = FALSE, meta = meta)
}
