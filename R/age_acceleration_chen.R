# Chen-style extrinsic age acceleration: an age clock blended with immune
# blood cell measures, with the blend weights re-estimated from the data at
# hand. Kept separate from the simple residual measures in `ageAcceleration()`
# and from the fixed-weight canonical measure in `EEAA()`.

# Klemera--Doubal biological age from several biomarkers that are each linear in
# chronological age. Each marker x_j is regressed on age to get slope k_j,
# intercept q_j and residual scale s_j; the combined estimate weights markers by
# k_j / s_j^2 (sharper, more age-informative markers count more). Returns NA for
# samples missing any marker.
.mc_kd_bioage <- function(markers, age) {
    markers <- as.matrix(markers)
    fits <- lapply(seq_len(ncol(markers)), function(j) {
        x <- markers[, j]
        ok <- is.finite(x) & is.finite(age)
        if (sum(ok) < 3L) return(NULL)
        f <- stats::lm(x[ok] ~ age[ok])
        k <- unname(stats::coef(f)[2L]); q <- unname(stats::coef(f)[1L])
        s <- stats::sigma(f)
        if (!is.finite(k) || !is.finite(s) || s <= 0) return(NULL)
        list(k = k, q = q, s = s)
    })
    keep <- !vapply(fits, is.null, logical(1))
    if (sum(keep) < 1L) return(rep(NA_real_, nrow(markers)))
    fits <- fits[keep]; M <- markers[, keep, drop = FALSE]
    k <- vapply(fits, `[[`, numeric(1), "k")
    q <- vapply(fits, `[[`, numeric(1), "q")
    s <- vapply(fits, `[[`, numeric(1), "s")
    num <- sweep(sweep(M, 2, q, "-"), 2, k / s^2, "*")
    ba <- rowSums(num) / sum(k^2 / s^2)
    ba[!stats::complete.cases(M)] <- NA_real_
    ba
}

#' Chen-style extrinsic age acceleration with data-driven weights
#'
#' An immune-aware measure of epigenetic age acceleration in the spirit of
#' Chen et al. (2016): a Klemera--Doubal weighted average of the Hannum clock
#' with immune cell types that track age (by default naive CD8+ T cells,
#' exhausted CD8+ T cells and plasmablasts) is formed, then regressed on
#' chronological age, and the returned \code{chenAcc} is that residual.
#'
#' The defining trait of this function is that the Klemera--Doubal weights
#' are re-estimated from YOUR data, and the immune inputs are whatever cell
#' measures you supply. The published EEAA instead fixes both: its weights
#' were fixed once in the Women's Health Initiative cohorts and its cell
#' inputs are three specific CpG-based estimators (the parameters of that
#' fixed recipe are published in European patent EP 3 494 210 B1, and
#' \code{\link{EEAA}} implements it). The patent documents the re-estimated
#' flavour as well (its "dynamic" variant) and reports the two to be highly
#' correlated; still, the two are numerically different measures, which is
#' why this function's output column is named \code{chenAcc} and not EEAA.
#'
#' @section Where the immune proportions can come from:
#' The subsets this measure was defined on --- naive CD8+ T cells, exhausted
#' CD8+ T cells (CD28-CD45RA-) and plasmablasts --- are not resolved by the
#' reference panels bundled here (see \code{\link{listCellReferences}}), which
#' report total CD8+ T cells. Two external routes are commonly used, and both
#' produce values you can pass straight to \code{cells}:
#' \itemize{
#'   \item the \emph{advanced analysis} of the online epigenetic-age calculator
#'     of Horvath (2013), which reports these subsets directly;
#'   \item an extended immune deconvolution such as the EPIC IDOL-Ext library of
#'     Salas et al. (2022), distributed as the
#'     \pkg{FlowSorted.BloodExtended.EPIC} package under its own licence from
#'     the authors' institution --- with the caveat that it resolves naive and
#'     memory CD8+ T cells and B cells but NOT exhausted CD8+ T cells or
#'     plasmablasts, so proportions from it are stand-ins: the result is a
#'     measure in the spirit of Chen's, not the one defined on the original
#'     subsets.
#' }
#' Users who obtain proportions this way should cite the source they used.
#' Note that the rarer subsets (plasmablasts, exhausted CD8+ T cells) are
#' reported by some algorithms as scores linearly related to abundance rather
#' than as cell-type proportions, which is worth keeping in mind when comparing
#' values across sources.
#'
#' @param hannum Numeric vector of Hannum-clock estimated ages, one per
#'   sample. May be named by sample id.
#' @param age Numeric vector of chronological ages, one per sample.
#' @param cells Samples-by-cell-types matrix with the immune cell measures;
#'   it must contain every column named in \code{immune}.
#' @param immune Column names in \code{cells} for the immune subsets that
#'   enter the weighted average. All of them must be present: the blend is
#'   only the intended measure with its full set of inputs, so a missing
#'   column is an error rather than a silent change of definition.
#' @return A data frame with columns \code{id}, \code{age} and
#'   \code{chenAcc}.
#' @references Chen BH, Marioni RE, Colicino E, et al. (2016). DNA
#'   methylation-based measures of biological age. \emph{Aging}
#'   8(9):1844--1865.
#'   Klemera P, Doubal S (2006). A new approach to the concept and computation
#'   of biological age. \emph{Mech Ageing Dev} 127(3):240--248.
#'   European patent EP 3 494 210 B1 (2024). DNA methylation based biological
#'   age prediction.
#' @seealso \code{\link{EEAA}} for the canonical fixed-weight measure,
#'   \code{\link{IEAA}}, \code{\link{ageAcceleration}},
#'   \code{\link{cellCounts}}
#' @examples
#' \dontrun{
#' # `subsets` holding naive CD8+, exhausted CD8+ and plasmablast measures
#' # per sample, e.g. from the Horvath calculator's advanced analysis:
#' chen <- ageAccelerationChen(res$Hannum, age = pheno$age, cells = subsets)
#' }
#' @export
ageAccelerationChen <- function(hannum, age, cells,
                                immune = c("CD8nv", "CD8exhausted",
                                           "PlasmaBlast")) {
    n <- length(age)
    if (length(hannum) != n)
        stop("`hannum` and `age` must have the same length.", call. = FALSE)
    cells <- as.matrix(cells)
    if (nrow(cells) != n)
        stop("`cells` must have one row per sample.", call. = FALSE)
    miss <- setdiff(immune, colnames(cells))
    if (length(miss))
        stop("`cells` is missing the immune column(s) ",
             paste0("'", miss, "'", collapse = ", "),
             ". The blend needs every column named in `immune`; see the ",
             "help page for where these subsets can come from.",
             call. = FALSE)
    id <- if (!is.null(names(hannum))) names(hannum)
          else if (!is.null(rownames(cells))) rownames(cells)
          else as.character(seq_len(n))

    ba <- .mc_kd_bioage(cbind(Hannum = hannum,
                              cells[, immune, drop = FALSE]), age)
    data.frame(id = id, age = age,
               chenAcc = .mc_resid_on(ba, data.frame(age = age)),
               stringsAsFactors = FALSE)
}
