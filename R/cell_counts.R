# Reference-based blood cell-type deconvolution (Houseman constrained least
# squares), the prerequisite for intrinsic age acceleration. Ported from the
# reference-panel machinery the original package vendored from meffil. The
# purified-cell reference panels live in the data layer, not the package.

# Available cell-type reference panels (by name).
#' Cell-type reference panels available for deconvolution
#'
#' The reference panels shipped with the package (via the data layer), each a
#' name accepted by \code{\link{cellCounts}}'s \code{reference} argument.
#' \code{cellCounts()} also accepts the higher-resolution
#' \code{"FlowSorted.Blood.EPIC"} and \code{"FlowSorted.BloodExtended.EPIC"}
#' panels, which are not shipped and are used only from your own installation of
#' those packages (see \code{\link{cellCounts}}).
#'
#' @return A character vector of the shipped reference-panel names.
#' @examples
#' listCellReferences()
#' @export
listCellReferences <- function() {
    c("andrews and bakulski cord blood", "blood gse35069",
      "blood gse35069 chen", "blood gse35069 complete", "combined cord blood",
      "cord blood gse68456", "gervin and lyle cord blood", "guintivano dlpfc",
      "saliva gse48472")
}

# Externally licensed panels supported by run-time detection only: we never ship
# their data, the user installs the package (under its own licence) and we use it
# if present. FlowSorted.Blood.EPIC is GPL-3; FlowSorted.BloodExtended.EPIC (the
# 12-cell high-resolution panel, Salas 2022) carries a restrictive non-commercial
# licence, so it can only be used from the user's own installation.
.MC_FLOWSORTED <- c("FlowSorted.Blood.EPIC", "FlowSorted.BloodExtended.EPIC")

# The 12 immune populations of the extended panel.
.MC_EXTENDED_CELLS <- c("Bas", "Bmem", "Bnv", "CD4mem", "CD4nv", "CD8mem",
                        "CD8nv", "Eos", "Mono", "Neu", "NK", "Treg")

# Build a CpGs x cell-types reference matrix from an installed FlowSorted package,
# so the same Houseman projection can run on it. Not shipped; loaded at run time.
.cell_reference_flowsorted <- function(reference, array = "EPIC") {
    if (!requireNamespace(reference, quietly = TRUE)) {
        hint <- if (grepl("Extended", reference))
            paste0("install it from github.com/immunomethylomics/",
                   reference, " under its (Dartmouth, non-commercial academic) ",
                   "licence")
        else paste0("install it with BiocManager::install(\"", reference, "\")")
        stop("Reference '", reference, "' needs the ", reference,
             " package, which methylclock does not ship; ", hint, ".",
             call. = FALSE)
    }
    # the objects may be exported, lazy-data, or ExperimentHub-backed; try each.
    got <- function(obj) {
        v <- tryCatch(getExportedValue(reference, obj), error = function(x) NULL)
        if (!is.null(v)) return(v)
        e <- new.env()
        loaded <- tryCatch(utils::data(list = obj, package = reference,
                                       envir = e), error = function(x) NULL)
        if (!is.null(loaded) && exists(obj, envir = e))
            return(get(obj, envir = e))
        ns <- tryCatch(getNamespace(reference), error = function(x) NULL)
        if (!is.null(ns) && exists(obj, envir = ns, inherits = FALSE))
            return(get(obj, envir = ns, inherits = FALSE))
        NULL
    }
    if (grepl("Extended", reference)) {
        is450 <- array == "450k"
        markers <- got(if (is450) "IDOLOptimizedCpGsBloodExtended450k"
                       else "IDOLOptimizedCpGsBloodExtended")
        comp <- got(if (is450) "FlowSorted.BloodExtended.450klegacy.compTable"
                    else "FlowSorted.BloodExtended.EPIC.compTable")
        cells <- .MC_EXTENDED_CELLS
    } else {
        is450 <- array == "450k"
        markers <- got(if (is450) "IDOLOptimizedCpGs450klegacy"
                       else "IDOLOptimizedCpGs")
        comp <- got(if (is450) "IDOLOptimizedCpGs450klegacy.compTable"
                    else "IDOLOptimizedCpGs.compTable")
        cells <- c("CD8T", "CD4T", "NK", "Bcell", "Mono", "Neu")
    }
    if (is.null(comp))
        stop("The reference table from ", reference, " could not be loaded ",
             "(its data objects were not found via getExportedValue(), data() ",
             "or ExperimentHub). Check that the package's data are installed ",
             "and loadable in your session.", call. = FALSE)
    comp <- as.matrix(comp)
    cells <- intersect(cells, colnames(comp))
    if (length(cells) < 2L)
        stop("The reference table from ", reference,
             " does not have the expected cell-type columns.", call. = FALSE)
    B <- comp[, cells, drop = FALSE]
    if (!is.null(markers))
        B <- B[intersect(as.character(markers), rownames(B)), , drop = FALSE]
    list(beta = B, quantiles = NULL, subsets = NULL, external = TRUE)
}

# Load a reference panel (a list with beta / quantiles / subsets) from the data
# mirror. Kept out of mcd_resource, which normalises to the coefficient contract.
.cell_reference <- function(name) {
    if (!name %in% listCellReferences())
        stop("Unknown cell reference '", name, "'. See listCellReferences().",
             call. = FALSE)
    id <- gsub("[^a-z0-9]+", "_", tolower(name))
    path <- file.path(.mcd_local_root(), "cellref", paste0(id, ".rds"))
    if (!file.exists(path))
        stop("Cell reference '", name, "' not found in the data mirror.",
             call. = FALSE)
    readRDS(path)
}

# Quantile-normalise the CpG subsets of a beta matrix to the reference's stored
# quantiles (as the reference deconvolution expects).
.qnorm_to_reference <- function(beta, subsets, quantiles) {
    for (sn in names(subsets)) {
        subset <- intersect(subsets[[sn]], rownames(beta))
        if (!length(subset)) next
        target <- quantiles[[sn]]$beta
        target <- stats::approx(seq_along(target), target,
                                seq_along(subset))$y
        tryCatch(
            beta[subset, ] <- preprocessCore::normalize.quantiles.use.target(
                beta[subset, , drop = FALSE], target),
            error = function(e) message("cellCounts: quantile normalisation ",
                "failed (", conditionMessage(e), ")."))
    }
    beta
}

# One sample: the non-negative cell proportions minimising ||b - B r||^2, via a
# quadratic program (Houseman et al., PMID 22568884).
.deconv_sample <- function(b, B) {
    ok <- which(!is.na(b))
    D <- crossprod(B[ok, , drop = FALSE])
    d <- crossprod(B[ok, , drop = FALSE], b[ok])
    r <- quadprog::solve.QP(D, d, diag(ncol(B)), rep(0, ncol(B)))$solution
    names(r) <- colnames(B)
    r
}

#' Estimate blood cell-type proportions from methylation
#'
#' Reference-based deconvolution: given a methylation matrix, it estimates the
#' proportion of each blood cell type per sample, by fitting the sample's profile
#' as a non-negative mix of purified-cell reference profiles (Houseman et al.
#' 2012). These proportions are what the intrinsic age acceleration adjusts for.
#'
#' @param betas A methylation matrix (or data frame, HDF5-backed matrix, or an
#'   accepted Bioconductor container) with CpG identifiers as row names.
#' @param reference Reference panel name (see \code{\link{listCellReferences}}).
#'   Default \code{"blood gse35069 complete"} (whole-blood, 7 cell types). The
#'   high-resolution FlowSorted panels are supported only from the user's own
#'   installation of the corresponding package (see Details).
#' @param normalize Quantile-normalise the CpGs to a shipped reference before
#'   fitting (ignored for the FlowSorted panels). Default \code{TRUE}.
#' @param array \code{"EPIC"} or \code{"450k"}; selects the marker set for the
#'   extended FlowSorted panel. Default \code{"EPIC"}.
#' @details The panels named in \code{\link{listCellReferences}} are shipped via
#'   the data layer. Two higher-resolution panels are also supported but not
#'   shipped, because of their licence: \code{"FlowSorted.Blood.EPIC"} (6 cells,
#'   GPL-3) and \code{"FlowSorted.BloodExtended.EPIC"} (12 immune cells, Salas et
#'   al. 2022, a restrictive non-commercial academic licence). Install the
#'   package yourself and it is used at run time; nothing licensed is
#'   redistributed by methylclock.
#' @return A matrix of samples (rows) by cell types (columns), each row summing
#'   to approximately one.
#' @references Houseman EA et al. (2012) \emph{BMC Bioinformatics} 13:86.
#'   Salas LA et al. (2022) \emph{Nature Communications} 13:761.
#' @examples
#' \dontrun{
#' cellCounts(betas)
#' cellCounts(betas, reference = "FlowSorted.BloodExtended.EPIC")
#' }
#' @export
cellCounts <- function(betas, reference = "blood gse35069 complete",
                       normalize = TRUE, array = c("EPIC", "450k")) {
    array <- match.arg(array)
    beta <- as.matrix(.extract_betas(betas))
    if (reference %in% .MC_FLOWSORTED) {
        ref <- .cell_reference_flowsorted(reference, array)
        normalize <- FALSE
    } else {
        ref <- .cell_reference(reference)
    }
    if (isTRUE(normalize))
        beta <- .qnorm_to_reference(beta, ref$subsets, ref$quantiles)
    cpgs <- intersect(rownames(beta), rownames(ref$beta))
    if (!length(cpgs))
        stop("None of the reference CpGs are present in the data.",
             call. = FALSE)
    if (length(cpgs) < 0.5 * nrow(ref$beta))
        warning(sprintf("only %d of %d reference CpGs present; ",
                        length(cpgs), nrow(ref$beta)),
                "cell-count estimates may be unreliable.", call. = FALSE)
    B <- ref$beta[cpgs, , drop = FALSE]
    out <- t(apply(beta[cpgs, , drop = FALSE], 2, .deconv_sample, B))
    colnames(out) <- colnames(ref$beta)
    out
}
