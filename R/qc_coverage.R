# How many of each clock's CpGs are present in a dataset --- the first thing to
# check, since a clock computed from few of its CpGs is on thin ice.

# A clock's CpG markers, from its coefficient resource (or, for the compiled
# neural-network clock, its auxiliary CpG set). Intercept and non-cg rows drop out.
.clock_cpgs <- function(entry) {
    from <- function(res) {
        if (is.null(res) || (length(res) == 1 && is.na(res)))
            return(character(0))
        co <- tryCatch(mcd_resource(res), error = function(e) NULL)
        if (is.null(co)) return(character(0))
        cc <- tryCatch(.coef_cpg_col(co), error = function(e) NULL)
        if (is.null(cc)) return(character(0))
        grep("^cg", as.character(co[[cc]]), value = TRUE)
    }
    m <- from(entry$resource)
    if (!length(m)) m <- from(entry$aux)
    m
}

#' CpG coverage of each clock in a dataset
#'
#' For every clock, how many of the CpGs it needs are present in your data. A
#' clock built from only a fraction of its CpGs should be read with caution (and
#' below the coverage threshold it is not computed at all). This needs the
#' methylation matrix, not the estimates, so it can look up which CpGs you have.
#'
#' @param betas A methylation matrix (or data frame, HDF5-backed matrix, or one
#'   of the accepted Bioconductor containers): only its CpG identifiers (row
#'   names) are used.
#' @param clocks Optional clock names to restrict to. Default: all clocks.
#' @return A data frame with one row per clock and columns \code{clock},
#'   \code{n_cpgs} (CpGs the clock uses), \code{present} (how many are in your
#'   data) and \code{pct_present}.
#' @examples
#' data(methylclock_betas)
#' clockCoverage(methylclock_betas, clocks = c("Horvath", "Hannum", "Levine"))
#' @export
clockCoverage <- function(betas, clocks = NULL) {
    x <- .extract_betas(betas)
    rn <- rownames(x)
    if (is.null(rn) && is.data.frame(x))
        rn <- as.character(x[[1]])
    if (is.null(rn))
        stop("`betas` needs CpG identifiers as row names.", call. = FALSE)
    entries <- if (is.null(clocks)) clock_registry()
               else lapply(clocks, clock_info)
    rows <- lapply(entries, function(e) {
        mk <- .clock_cpgs(e)
        if (!length(mk)) return(NULL)
        present <- sum(mk %in% rn)
        data.frame(clock = e$name, n_cpgs = length(mk), present = present,
                   pct_present = round(100 * present / length(mk), 1),
                   stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
    out[order(-out$pct_present), ]
}

# Per-clock missingness: among the CpGs a clock uses that are present in the
# data, the fraction of values (over all samples) that were missing --- and so
# imputed before the clock ran. Complements clockCoverage(), which counts the
# CpGs absent from the array altogether. Reads only each clock's CpG rows.
.clock_missing <- function(betas, clocks = NULL) {
    b <- .extract_betas(betas)
    rn <- rownames(b)
    if (is.null(rn) && is.data.frame(b)) {
        rn <- as.character(b[[1]])
        b <- b[, -1, drop = FALSE]
        rownames(b) <- rn
    }
    if (is.null(rn)) return(NULL)
    entries <- if (is.null(clocks)) clock_registry()
               else lapply(clocks, clock_info)
    rows <- lapply(entries, function(e) {
        mk <- intersect(.clock_cpgs(e), rn)
        if (!length(mk)) return(NULL)
        sub <- as.matrix(b[mk, , drop = FALSE])
        data.frame(clock = e$name, present = length(mk),
                   n_missing = sum(is.na(sub)),
                   pct_missing = round(100 * mean(is.na(sub)), 2),
                   stringsAsFactors = FALSE)
    })
    out <- do.call(rbind, rows[!vapply(rows, is.null, logical(1))])
    if (is.null(out)) return(NULL)
    out[order(-out$pct_missing), ]
}
