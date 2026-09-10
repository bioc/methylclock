# Manifest-backed resource resolver: the single access point for clock
# coefficient objects. A resource identifier is resolved against a manifest
# through a configurable chain of backends and cached in a private environment.
#
# INTERNAL: replaces loadClocks.R / loadClocksGA.R /
# INTERNAL: meffilListCellTypeReferences.R, which assigned to .GlobalEnv 24
# INTERNAL: times (B-06, B-07). See DESIGN.md sec 7 and
# INTERNAL: DATA_LAYER_AND_NEW_CLOCKS.md. All three backends are wired: `local`
# INTERNAL: reads the data mirror; `eh` reads ExperimentHub by the manifest
# INTERNAL: `eh_id`; `zenodo` downloads by the manifest `zenodo_doi`. A resource
# INTERNAL: whose manifest row has no id for a given source returns NULL there,
# INTERNAL: so the chain falls through. Publishing (filling eh_id/zenodo_doi,
# INTERNAL: creating methylclockData) is the remaining step (D6) and needs no
# INTERNAL: change here.

# Private cache.
.mcd_cache <- new.env(parent = emptyenv())

#' Configure the resource backend chain
#'
#' Sets or returns the ordered list of sources consulted to resolve a resource.
#' The chain is stored in \code{options(methylclockData.backend=)}.
#'
#' @param backends Optional character vector to set, e.g. \code{c("local")}.
#'   With no argument, returns the current chain.
#' @return The current backend chain; returned invisibly when setting.
#' @details Order matters: the first backend that resolves an identifier wins.
#'   Supported backends:
#'   \describe{
#'     \item{\code{"local"}}{a directory of \code{.rda} files with an
#'       accompanying \code{MANIFEST.csv} (set its location with
#'       \code{options(methylclockData.local_root=)} or the
#'       \code{METHYLCLOCKDATA_MIRROR} environment variable).}
#'     \item{\code{"eh"}}{ExperimentHub; resolves the resource by the manifest's
#'       \code{eh_id} column. Requires the \pkg{ExperimentHub} package.}
#'     \item{\code{"zenodo"}}{downloads and caches the resource for the record in
#'       the manifest's \code{zenodo_doi} column.}
#'   }
#'   A backend returns nothing for a resource whose manifest row has no id for
#'   that source, so the chain falls through to the next backend.
#' @seealso \code{\link{mcd_resource}}
#' @examples
#' mcd_backends()                 # the current chain
#' old <- mcd_backends()
#' mcd_backends(c("local", "eh")) # set a chain
#' mcd_backends(old)              # and restore it
#' @export
mcd_backends <- function(backends) {
    if (missing(backends))
        return(getOption("methylclockData.backend",
                         default = c("local", "eh", "zenodo")))
    stopifnot(is.character(backends), length(backends) >= 1L)
    options(methylclockData.backend = backends)
    invisible(backends)
}

# Resolve the local root (a directory holding rda/ + MANIFEST.csv).
# Precedence: option -> env var -> nearby data directory. With required = FALSE
# an unset root returns NULL so a backend can decline and let the chain move on
# to the next source instead of failing the whole resolution.
.mcd_local_root <- function(required = TRUE) {
    root <- getOption("methylclockData.local_root",
                      default = Sys.getenv("METHYLCLOCKDATA_MIRROR", unset = NA))
    if (is.na(root) || !nzchar(root)) {
        # INTERNAL: development fallback to the verified data-mirror tree.
        for (cand in c("../data-mirror", "data-mirror",
                       file.path("..", "..", "data-mirror")))
            if (dir.exists(cand)) { root <- cand; break }
    }
    if (is.na(root) || !dir.exists(root)) {
        if (!required) return(NULL)
        stop("Local data root not found. Set options(methylclockData.local_root=) ",
             "or the METHYLCLOCKDATA_MIRROR env var to a directory with rda/ + ",
             "MANIFEST.csv.")
    }
    root
}

#' The resource manifest
#'
#' Reads and caches the manifest that describes the available resources
#' (identifier, stored object name, class, dimensions, columns, size, checksum
#' and source). A copy ships with the package; a configured local data root
#' (see \code{\link{mcd_backends}}) takes precedence when present.
#' @return A data frame, one row per resource identifier.
#' @seealso \code{\link{mcd_resource}}
#' @examples
#' man <- mcd_manifest()
#' head(man$id)
#' @export
mcd_manifest <- function() {
    if (!is.null(.mcd_cache$.manifest))
        return(.mcd_cache$.manifest)
    # A configured local root takes precedence (it may carry newer rows);
    # otherwise fall back to the copy shipped with the package.
    root <- .mcd_local_root(required = FALSE)
    path <- if (!is.null(root)) file.path(root, "MANIFEST.csv") else ""
    if (!file.exists(path))
        path <- system.file("extdata", "MANIFEST.csv", package = "methylclock")
    if (!nzchar(path) || !file.exists(path))
        stop("MANIFEST.csv not found (neither a local data root nor the ",
             "installed package provides one).")
    man <- utils::read.csv(path, stringsAsFactors = FALSE)
    .mcd_cache$.manifest <- man
    man
}

# local backend: load `<root>/rda/<id>.rda` and return the named object.
.mcd_load_local <- function(id, row) {
    root <- .mcd_local_root(required = FALSE)
    if (is.null(root))
        return(NULL)                         # no mirror -> try next backend
    path <- file.path(root, "rda", paste0(id, ".rda"))
    if (!file.exists(path))
        return(NULL)                         # not resolvable here -> try next backend
    e <- new.env(parent = emptyenv())
    objs <- load(path, envir = e)
    obj_name <- if (!is.na(row$object) && nzchar(row$object)) row$object else objs[1]
    if (!obj_name %in% objs)
        obj_name <- objs[1]
    get(obj_name, envir = e)
}

# A manifest field may be absent, NA or empty; treat all three as "not set".
.mcd_field <- function(row, name) {
    if (is.null(row) || !name %in% names(row)) return(NA_character_)
    v <- row[[name]]
    if (length(v) != 1L || is.na(v) || !nzchar(as.character(v))) return(NA_character_)
    as.character(v)
}

# Some data hubs return the object directly; others return a path to an .rda
# that must be loaded. Normalize both to the stored object (named by the
# manifest's `object` column when the archive holds several).
.mcd_materialize <- function(res, row) {
    if (is.character(res) && length(res) == 1L && file.exists(res)) {
        e <- new.env(parent = emptyenv())
        objs <- load(res, envir = e)
        want <- .mcd_field(row, "object")
        obj_name <- if (!is.na(want) && want %in% objs) want else objs[1]
        return(get(obj_name, envir = e))
    }
    res
}

# ExperimentHub backend: fetch the resource named by the manifest's `eh_id`.
# Returns NULL when the row carries no ExperimentHub identifier, so the backend
# chain falls through to the next entry.
.mcd_load_eh <- function(id, row) {
    eh_id <- .mcd_field(row, "eh_id")
    if (is.na(eh_id)) return(NULL)
    if (!requireNamespace("ExperimentHub", quietly = TRUE)) {
        warning("ExperimentHub is not installed; cannot use the 'eh' backend.")
        return(NULL)
    }
    hub <- ExperimentHub::ExperimentHub()
    res <- tryCatch(hub[[eh_id]], error = function(e) NULL)
    if (is.null(res)) return(NULL)
    .mcd_materialize(res, row)
}

# Zenodo backend: download and cache the resource for a record addressed by the
# manifest's `zenodo_doi`. Dependency-light (no zen4R): the file is fetched by
# convention as `<id>.rda` from the record. Returns NULL when no DOI is set or
# the download fails, so the chain can fall through.
.mcd_zenodo_fetch <- function(doi, id) {
    cache <- tools::R_user_dir("methylclock", which = "cache")
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    dest <- file.path(cache, paste0(id, ".rda"))
    if (file.exists(dest) && file.size(dest) > 0L) return(dest)
    rec <- sub(".*zenodo\\.", "", doi)       # 10.5281/zenodo.<rec> -> <rec>
    if (!grepl("^[0-9]+$", rec)) return(NULL)
    url <- sprintf("https://zenodo.org/records/%s/files/%s.rda?download=1",
                   rec, id)
    ok <- tryCatch({
        utils::download.file(url, dest, mode = "wb", quiet = TRUE)
        file.exists(dest) && file.size(dest) > 0L
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!isTRUE(ok)) {
        if (file.exists(dest)) unlink(dest)
        return(NULL)
    }
    dest
}

.mcd_load_zenodo <- function(id, row) {
    doi <- .mcd_field(row, "zenodo_doi")
    if (is.na(doi)) return(NULL)
    path <- .mcd_zenodo_fetch(doi, id)
    if (is.null(path)) return(NULL)
    .mcd_materialize(path, row)
}

# Zenodo backend for FILE resources (e.g. an .hdf5 of network weights): same
# record convention as .mcd_zenodo_fetch, with the file's own extension.
.mcd_file_zenodo <- function(id, row) {
    doi <- .mcd_field(row, "zenodo_doi")
    if (is.na(doi)) return(NULL)
    ext <- if (identical(.mcd_field(row, "class"), "hdf5")) ".hdf5" else
        ".rda"
    cache <- tools::R_user_dir("methylclock", which = "cache")
    dir.create(cache, recursive = TRUE, showWarnings = FALSE)
    dest <- file.path(cache, paste0(id, ext))
    if (file.exists(dest) && file.size(dest) > 0L) return(dest)
    rec <- sub(".*zenodo\\.", "", doi)
    if (!grepl("^[0-9]+$", rec)) return(NULL)
    url <- sprintf("https://zenodo.org/records/%s/files/%s%s?download=1",
                   rec, id, ext)
    ok <- tryCatch({
        utils::download.file(url, dest, mode = "wb", quiet = TRUE)
        file.exists(dest) && file.size(dest) > 0L
    }, error = function(e) FALSE, warning = function(w) FALSE)
    if (!isTRUE(ok)) {
        if (file.exists(dest)) unlink(dest)
        return(NULL)
    }
    dest
}

# Optional checksum verification (requires the digest package).
.mcd_verify_sha256 <- function(obj, row, id) {
    if (is.null(row$sha256) || is.na(row$sha256) || !nzchar(row$sha256))
        return(invisible(TRUE))
    if (!requireNamespace("digest", quietly = TRUE)) {
        warning("digest not installed; skipping checksum for '", id, "'.")
        return(invisible(NA))
    }
    got <- digest::digest(obj, algo = "sha256", serialize = TRUE)
    # INTERNAL: MANIFEST sha256 was computed over the on-disk artifact; object-
    # INTERNAL: level hashing is a placeholder until the manifest records a hash
    # INTERNAL: we can reproduce here. Hook only (verify defaults to FALSE).
    invisible(identical(got, row$sha256))
}

#' Resolve a clock resource by identifier
#'
#' The single access point for coefficient objects. Resolves \code{id} through
#' the backend chain (see \code{\link{mcd_backends}}), caches the result in a
#' private environment and returns it.
#'
#' @param id Resource identifier, e.g. \code{"coefHorvath"} (a row of
#'   \code{\link{mcd_manifest}}).
#' @param verify Logical; verify the object's checksum when possible. Default
#'   \code{FALSE}.
#' @param cache Logical; use and populate the private cache. Default
#'   \code{TRUE}.
#' @return The resource object in canonical form (see \code{\link{normalize_coef}}).
#' @seealso \code{\link{mcd_manifest}}, \code{\link{mcd_backends}},
#'   \code{\link{clock_register}}
#' @examples
#' coef <- mcd_resource("coefHorvath")
#' head(coef)
#' @export
mcd_resource <- function(id, verify = FALSE, cache = TRUE) {
    stopifnot(is.character(id), length(id) == 1L)
    key <- paste0("res:", id)
    if (cache && !is.null(.mcd_cache[[key]]))
        return(.mcd_cache[[key]])

    man <- mcd_manifest()
    idx <- match(id, man$id)
    if (is.na(idx))
        stop(sprintf("Unknown resource '%s' (not in manifest).", id))
    row <- man[idx, , drop = FALSE]

    obj <- NULL
    for (backend in mcd_backends()) {
        obj <- switch(backend,
            local  = .mcd_load_local(id, row),
            eh     = .mcd_load_eh(id, row),
            zenodo = .mcd_load_zenodo(id, row),
            stop("Unknown backend '", backend, "'.")
        )
        if (!is.null(obj)) break
    }
    if (is.null(obj))
        stop(sprintf("Resource '%s' could not be resolved by backends: %s",
                     id, paste(mcd_backends(), collapse = ", ")))

    if (isTRUE(verify)) {
        ok <- .mcd_verify_sha256(obj, row, id)
        if (isFALSE(ok))
            warning("Checksum mismatch for resource '", id, "'.")
    }

    obj <- normalize_coef(obj, id = id, row = row)
    if (cache) .mcd_cache[[key]] <- obj
    obj
}

# local backend for a file resource: an HDF5 bundle under <root>/altumage/.
.mcd_file_local <- function(id, row) {
    root <- .mcd_local_root(required = FALSE)
    if (is.null(root))
        return(NULL)                         # no mirror -> try next backend
    path <- file.path(root, "altumage", paste0(id, ".hdf5"))
    if (file.exists(path)) path else NULL
}

# ExperimentHub backend for a file resource: download (or reuse) the hub's
# cached copy of the file named by the manifest's `eh_id` and return its path.
# ExperimentHub::cache() gives the on-disk file without loading it into R.
.mcd_file_eh <- function(id, row) {
    eh_id <- .mcd_field(row, "eh_id")
    if (is.na(eh_id)) return(NULL)
    if (!requireNamespace("ExperimentHub", quietly = TRUE))
        return(NULL)
    path <- tryCatch({
        hub <- ExperimentHub::ExperimentHub()
        unname(ExperimentHub::cache(hub[eh_id]))
    }, error = function(e) NULL)
    if (is.null(path) || !length(path) || !file.exists(path[1]))
        return(NULL)
    path[1]
}

#' Resolve a file resource to a local path
#'
#' Some resources are files read directly by the C++ layer (e.g. an HDF5 bundle
#' of neural-network weights) rather than objects loaded into R. This resolves
#' \code{id} through the backend chain and returns a local path: the file on
#' disk for the \code{local} backend, or a downloaded-and-cached copy for a
#' remote backend.
#'
#' @param id Resource identifier (a row of \code{\link{mcd_manifest}} whose
#'   class marks it as a file resource, e.g. \code{"coefAltumAge"}).
#' @return A path to the resource file.
#' @seealso \code{\link{mcd_resource}}
#' @export
mcd_resource_file <- function(id) {
    stopifnot(is.character(id), length(id) == 1L)
    key <- paste0("file:", id)
    if (!is.null(.mcd_cache[[key]]))
        return(.mcd_cache[[key]])
    man <- mcd_manifest()
    idx <- match(id, man$id)
    if (is.na(idx))
        stop(sprintf("Unknown resource '%s' (not in manifest).", id))
    row <- man[idx, , drop = FALSE]
    path <- NULL
    for (backend in mcd_backends()) {
        path <- switch(backend,
            local  = .mcd_file_local(id, row),
            eh     = .mcd_file_eh(id, row),
            zenodo = .mcd_file_zenodo(id, row),
            stop("Unknown backend '", backend, "'.")
        )
        if (!is.null(path)) break
    }
    if (is.null(path))
        stop(sprintf("File resource '%s' could not be resolved by backends: %s",
                     id, paste(mcd_backends(), collapse = ", ")))
    .mcd_cache[[key]] <- path
    path
}

#' Clear the resource cache
#'
#' Empties the private cache of resolved resources. Resources are re-resolved
#' through the backend chain the next time they are needed.
#' @return Invisibly \code{NULL}.
#' @seealso \code{\link{mcd_resource}}, \code{\link{mcd_backends}}
#' @examples
#' mcd_cache_clear()
#' @export
mcd_cache_clear <- function() {
    rm(list = ls(.mcd_cache), envir = .mcd_cache)
    invisible(NULL)
}
