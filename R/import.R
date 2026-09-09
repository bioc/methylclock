# Convert a methylation matrix (in memory or in a text file) to an HDF5-backed
# matrix, so it can be estimated on disk in blocks rather than held in RAM.
#
# INTERNAL: D18. The genuine low-memory path needs the data to start on disk: a
# INTERNAL: matrix passed as an argument stays resident (caller holds it), so
# INTERNAL: converting once and dropping the dense copy is what lowers peak RAM.
# INTERNAL: TODO(bdsm): bdImportTextFile_hdf5 with rownames=TRUE errors
# INTERNAL: off-by-one ("index N/N in SET_STRING_ELT") in BDSM 2.0.3, so a text
# INTERNAL: file is read with base R for now (fits-in-RAM only). Fixing the BDSM
# INTERNAL: importer would let files too big for RAM stream in.

#' Store a methylation matrix on disk as HDF5
#'
#' Writes a methylation matrix to an HDF5 file and returns it as an on-disk
#' matrix that \code{\link{methylclock}} can estimate block by block, without
#' holding the whole array in memory. Useful for arrays too large to work with
#' comfortably in RAM: convert once, then pass the result (and let the in-memory
#' copy go).
#'
#' @param x A methylation matrix or data frame (CpG identifiers as row names / a
#'   CpG-identifier column, samples as columns), or the path to a delimited text
#'   file laid out the same way. A text file is imported without being loaded
#'   into memory first.
#' @param file Destination HDF5 file. Pass a real path to keep it; a session
#'   temporary file is used when \code{NULL} (default).
#' @param overwrite Overwrite \code{file} if it exists. Default \code{FALSE}.
#' @param sep For a text file, the field separator (default: auto by extension,
#'   comma for \code{.csv} else tab).
#' @param header,row.names For a text file, whether it has a header row / a
#'   leading column of CpG identifiers. Both default to \code{TRUE}.
#' @param compression HDF5 compression level (0-9); \code{NULL} (default) picks
#'   it by size.
#' @return An on-disk (HDF5-backed) matrix, open for reading; pass it to
#'   \code{\link{methylclock}} and \code{close()} it when done.
#' @seealso \code{\link{methylclock}}
#' @examples
#' data(methylclock_betas)
#' f <- tempfile(fileext = ".h5")
#' hm <- mc_to_hdf5(methylclock_betas, f)
#' res <- methylclock(hm, clocks = "Horvath")
#' head(as.data.frame(res))
#' close(hm)
#' unlink(f)
#' @export
mc_to_hdf5 <- function(x, file = NULL, overwrite = FALSE, sep = NULL,
                       header = TRUE, row.names = TRUE, compression = NULL) {
    .bdsm_require()
    if (is.null(file))
        file <- tempfile(fileext = ".h5")
    else if (file.exists(file) && !overwrite)
        stop("`", file, "` already exists; pass overwrite = TRUE to replace.")
    path <- paste("betas", "values", sep = "/")

    if (is.character(x) && length(x) == 1L && file.exists(x)) {
        if (is.null(sep))
            sep <- if (grepl("\\.csv$", x, ignore.case = TRUE)) "," else "\t"
        df <- utils::read.delim(x, sep = sep, header = header,
                                row.names = if (isTRUE(row.names)) 1L else NULL,
                                check.names = FALSE)
        m <- as.matrix(df)
    } else {
        m <- .betas_from_mvalues(.as_beta_matrix(.extract_betas(x)))
    }
    comp <- .mc_compression(m, compression)
    h <- BigDataStatMeth::hdf5_create_matrix(
        filename = file, dataset = path, data = m,
        overwrite = TRUE, compression = comp)
    .mc_close(h)
    BigDataStatMeth::hdf5_matrix(file, path)
}

# Run the estimation, turning an out-of-memory failure into an actionable error
# that points at the on-disk route. `expr` is a promise, forced inside tryCatch.
# INTERNAL: D18. We do NOT auto-spill a resident matrix: writing it to disk
# INTERNAL: needs it in RAM, so the spill peak is >= computing in place
# INTERNAL: (measured). The real low-peak route is data that never becomes a
# INTERNAL: dense matrix (import a text file with mc_to_hdf5, or make HDF5
# INTERNAL: upstream).
.alloc_guard <- function(expr) {
    tryCatch(expr, error = function(e) {
        if (grepl("cannot allocate", conditionMessage(e), fixed = TRUE))
            stop("methylclock: not enough memory to compute this in RAM. Give ",
                 "the input on disk instead -- import it with mc_to_hdf5() (a ",
                 "text file is not loaded whole) and pass that.", call. = FALSE)
        stop(e)
    })
}
