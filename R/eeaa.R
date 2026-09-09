# Canonical extrinsic epigenetic age acceleration (EEAA): the fixed-weight
# blend of the Hannum clock with three CpG-based immune cell estimators, as
# defined by Chen et al. (2016) with the weights fixed in the Women's Health
# Initiative. The numeric parameters (cell-estimator CpG coefficients,
# standardisation centers/scales and blend weights) are the ones published in
# European patent EP 3 494 210 B1; they ship as internal data.

# One CpG-based cell score: the published linear combination of beta values.
# The patent's estimators carry an intercept that it does not print; a missing
# additive constant shifts every sample equally and is absorbed by the final
# regression on age, so the residual (EEAA) is unaffected. Scattered missing
# values are filled with the CpG's own mean across samples.
.mc_eeaa_cell_score <- function(B, coefs, what) {
    present <- intersect(names(coefs), rownames(B))
    if (length(present) < length(coefs))
        warning(sprintf("%s: only %d of %d CpGs present; ", what,
                        length(present), length(coefs)),
                "the score deviates from the published estimator.",
                call. = FALSE)
    if (length(present) < 0.5 * length(coefs))
        stop("Too few CpGs for the ", what, " estimator (",
             length(present), " of ", length(coefs), ").", call. = FALSE)
    M <- B[present, , drop = FALSE]
    if (anyNA(M)) {
        mu <- rowMeans(M, na.rm = TRUE)
        idx <- which(is.na(M), arr.ind = TRUE)
        M[idx] <- mu[idx[, 1L]]
    }
    colSums(M * coefs[present])
}

#' Extrinsic epigenetic age acceleration (canonical EEAA)
#'
#' The EEAA of the literature (Chen et al. 2016): the Hannum clock blended
#' with three immune cell estimators --- plasmablasts, exhausted CD8+ T cells
#' (CD28-CD45RA-) and naive CD8+ T cells --- with weights fixed once in the
#' Women's Health Initiative, then regressed on chronological age; EEAA is
#' the residual. By construction it rises with the immune shifts of ageing,
#' which is what makes it "extrinsic", and it is the specific measure
#' researchers mean by "EEAA" (its intrinsic sibling is \code{\link{IEAA}}).
#'
#' The three cell estimators are computed here from the methylation data
#' itself, as linear combinations of published CpG coefficients --- the same
#' route the original software takes --- so no external deconvolution is
#' needed; \code{betas} should be blood 450K/EPIC data. All numeric
#' parameters (the estimators' CpG coefficients and the fixed
#' centers, scales and weights of the blend) are the ones published in
#' European patent EP 3 494 210 B1. The estimators' additive intercepts are
#' not published; they shift every sample equally and cancel in the residual,
#' which is why this function returns EEAA and not the intermediate
#' biological-age value. If some of the estimator CpGs are absent from your
#' array a warning reports it: the scores then deviate from the published
#' estimators and comparability with other cohorts' EEAA weakens.
#'
#' @param x A \code{methylclock} result, or a data frame with an \code{id}
#'   column and a numeric \code{Hannum} column.
#' @param age Numeric vector of chronological ages, one per sample in
#'   \code{x}'s order (or named by sample id).
#' @param betas The methylation data the estimate came from: a matrix (or
#'   data frame, HDF5-backed matrix, or accepted Bioconductor container) with
#'   CpG identifiers as row names and one column per sample. Only the rows of
#'   the estimator CpGs are read, so an on-disk matrix is not loaded whole.
#' @return A data frame with columns \code{id} and \code{EEAA}.
#' @references Chen BH, Marioni RE, Colicino E, et al. (2016). DNA
#'   methylation-based measures of biological age. \emph{Aging}
#'   8(9):1844--1865.
#'   European patent EP 3 494 210 B1 (2024). DNA methylation based biological
#'   age prediction.
#' @seealso \code{\link{IEAA}}, \code{\link{ageAcceleration}},
#'   \code{\link{ageAccelerationChen}}
#' @examples
#' \dontrun{
#' res <- DNAmAge(betas, clocks = "Hannum")
#' head(EEAA(res, age = pheno$age, betas = betas))
#' }
#' @export
EEAA <- function(x, age, betas) {
    v <- .mc_plot_values(x, "Hannum")
    hannum <- v$values[["Hannum"]]
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

    scores <- cbind(
        DNAmAge = hannum,
        PlasmaBlast = .mc_eeaa_cell_score(B, p$plasmablast, "plasmablast"),
        CD8pCD28nCD45RAn = .mc_eeaa_cell_score(B, p$cd8exhausted,
                                               "exhausted CD8+ T cell"),
        CD8.naive = .mc_eeaa_cell_score(B, p$cd8naive,
                                        "naive CD8+ T cell"))
    std <- sweep(sweep(scores, 2, p$centers[colnames(scores)], "-"),
                 2, p$scales[colnames(scores)], "/")
    bioage <- as.numeric(std %*% p$weights[colnames(std)])

    data.frame(id = v$ids,
               EEAA = .mc_resid_on(bioage, data.frame(age = age)),
               stringsAsFactors = FALSE)
}
