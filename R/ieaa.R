# Canonical intrinsic epigenetic age acceleration (IEAA): the Horvath clock
# regressed on chronological age and the seven blood cell covariates of Chen
# et al. (2016) -- naive CD8+ T cells, exhausted CD8+ T cells, plasmablasts,
# CD4+ T cells, natural killer cells, monocytes and granulocytes. The three
# rare subsets come from the same published CpG estimators the canonical
# software uses (shipped as internal data, see `EEAA()`); the four common
# ones come from reference-based deconvolution.

#' Intrinsic epigenetic age acceleration (canonical IEAA)
#'
#' The IEAA of the literature (Chen et al. 2016): the Horvath clock's
#' estimate regressed on chronological age and seven blood immune cell
#' covariates --- naive CD8+ T cells, exhausted CD8+ T cells
#' (CD28-CD45RA-), plasmablasts, CD4+ T cells, natural killer cells,
#' monocytes and granulocytes --- whose residual reflects cell-intrinsic
#' ageing independent of shifts in blood cell composition. This is the
#' specific measure researchers mean by "IEAA" (its extrinsic sibling is
#' \code{\link{EEAA}}). The same adjustment idea for any clock and any cell
#' panel is available as the \code{residualCells} column of
#' \code{\link{ageAcceleration}}, under that plainer name.
#'
#' The three rare subsets are computed from the methylation data itself, as
#' the published linear combinations of CpGs the original software uses
#' (parameters from European patent EP 3 494 210 B1; their additive
#' intercepts are not published, but constants are absorbed by the
#' regression, so the residual is unaffected). The four common cell types
#' are taken from \code{cell_counts} if you provide them, or estimated with
#' \code{\link{cellCounts}} otherwise; granulocytes are used directly when
#' a \code{Gran} column exists, and otherwise formed as neutrophils plus
#' eosinophils. Reference-based deconvolution estimates are not numerically
#' identical to the original calculator's own estimates of these four, so
#' the measure is canonical in structure and covariates while small
#' numerical differences with that software remain possible. If some of the
#' estimator CpGs are absent from your array a warning reports it.
#'
#' @param x A \code{methylclock} result, or a data frame with an \code{id}
#'   column and a numeric \code{Horvath} column.
#' @param age Numeric vector of chronological ages, one per sample in
#'   \code{x}'s order (or named by sample id).
#' @param betas The methylation data the estimate came from: a matrix (or
#'   data frame, HDF5-backed matrix, or accepted Bioconductor container)
#'   with CpG identifiers as row names and one column per sample; blood
#'   450K/EPIC data. The rare-subset estimators read only their CpG rows.
#' @param cell_counts Optional samples-by-cell-types matrix (for example
#'   from \code{\link{cellCounts}}) supplying the common cell types:
#'   \code{CD4T}, \code{NK}, \code{Mono} and \code{Gran} (or \code{Neu},
#'   plus \code{Eos} when available). When omitted, \code{cellCounts(betas)}
#'   is run internally, which loads the full matrix into memory.
#' @return A data frame with columns \code{id} and \code{IEAA}.
#' @references Chen BH, Marioni RE, Colicino E, et al. (2016). DNA
#'   methylation-based measures of biological age. \emph{Aging}
#'   8(9):1844--1865.
#'   European patent EP 3 494 210 B1 (2024). DNA methylation based
#'   biological age prediction.
#' @seealso \code{\link{EEAA}}, \code{\link{ageAcceleration}},
#'   \code{\link{cellCounts}}
#' @examples
#' \dontrun{
#' res <- DNAmAge(betas, clocks = "Horvath")
#' head(IEAA(res, age = pheno$age, betas = betas))
#' }
#' @export
IEAA <- function(x, age, betas, cell_counts = NULL) {
    v <- .mc_plot_values(x, "Horvath")
    horvath <- v$values[["Horvath"]]
    age <- .mc_align_age(age, v$ids)

    p <- .mc_eeaa_params
    src <- .as_beta_source(betas)
    cn <- .src_colnames(src)
    if (!is.null(cn) && all(v$ids %in% cn)) {
        sel <- match(v$ids, cn)
    } else if (.src_ncol(src) == length(v$ids)) {
        sel <- seq_along(v$ids)
    } else {
        stop("`betas` does not cover the samples in `x` (match them by ",
             "column name, or pass the same samples in the same order).",
             call. = FALSE)
    }

    cpgs <- unique(c(names(p$plasmablast), names(p$cd8exhausted),
                     names(p$cd8naive)))
    present <- intersect(cpgs, .src_rownames(src))
    if (!length(present))
        stop("None of the immune-estimator CpGs are present in `betas`.",
             call. = FALSE)
    B <- .src_subset(src, present)[, sel, drop = FALSE]

    cc <- if (is.null(cell_counts)) cellCounts(betas) else
        as.matrix(cell_counts)
    if (!is.null(rownames(cc)) && all(v$ids %in% rownames(cc)))
        cc <- cc[v$ids, , drop = FALSE]
    else if (nrow(cc) != length(v$ids))
        stop("`cell_counts` must have one row per sample (or be named by ",
             "sample id).", call. = FALSE)
    need <- c("CD4T", "NK", "Mono")
    miss <- setdiff(need, colnames(cc))
    if (length(miss) || !any(c("Gran", "Neu") %in% colnames(cc)))
        stop("The common cell covariates need columns ",
             "'CD4T', 'NK', 'Mono' and 'Gran' (or 'Neu'); missing: ",
             paste0("'", c(miss,
                           if (!any(c("Gran", "Neu") %in% colnames(cc)))
                               "Gran/Neu"), "'", collapse = ", "), ".",
             call. = FALSE)
    gran <- if ("Gran" %in% colnames(cc)) cc[, "Gran"]
            else cc[, "Neu"] +
                (if ("Eos" %in% colnames(cc)) cc[, "Eos"] else 0)

    covars <- data.frame(
        age = age,
        CD8.naive = .mc_eeaa_cell_score(B, p$cd8naive,
                                        "naive CD8+ T cell"),
        CD8pCD28nCD45RAn = .mc_eeaa_cell_score(B, p$cd8exhausted,
                                               "exhausted CD8+ T cell"),
        PlasmaBlast = .mc_eeaa_cell_score(B, p$plasmablast, "plasmablast"),
        CD4T = cc[, "CD4T"], NK = cc[, "NK"], Mono = cc[, "Mono"],
        Gran = gran)

    data.frame(id = v$ids,
               IEAA = .mc_resid_on(horvath, covars),
               stringsAsFactors = FALSE)
}
