# Input handling: accept the several shapes methylation data arrive in, and
# accept M-values as well as betas. Kept separate from the source machinery so
# the "what can I pass in" rules live in one place.

# Pull a CpGs x samples matrix out of the containers people hold after processing
# IDATs, so methylclock() takes more than a bare matrix. A matrix, data frame or
# HDF5-backed matrix passes through untouched (handled downstream); the Bioconductor
# containers are read through their own accessors, loaded only if present.
.extract_betas <- function(x) {
    if (is.matrix(x) || is.data.frame(x) || inherits(x, "HDF5Matrix"))
        return(x)
    need <- function(pkg, what) {
        if (!requireNamespace(pkg, quietly = TRUE))
            stop("Input looks like ", what, "; install '", pkg,
                 "' to use it, or pass a beta matrix.", call. = FALSE)
    }
    if (inherits(x, c("MethylSet", "GenomicMethylSet", "RatioSet",
                      "GenomicRatioSet"))) {
        need("minfi", "a minfi object")
        return(minfi::getBeta(x))
    }
    if (inherits(x, "ExpressionSet")) {
        need("Biobase", "an ExpressionSet")
        return(Biobase::exprs(x))
    }
    if (inherits(x, "SummarizedExperiment")) {
        need("SummarizedExperiment", "a SummarizedExperiment")
        return(as.matrix(SummarizedExperiment::assay(x)))
    }
    stop("Unsupported input. Provide a beta (or M-value) matrix or data frame ",
         "with CpG identifiers as row names, an HDF5-backed matrix, or a minfi, ",
         "SummarizedExperiment or ExpressionSet object.", call. = FALSE)
}

# Accept M-values as well as betas. Beta values lie in [0, 1]; M-values are an
# unbounded log-ratio. If the data fall outside [0, 1] we treat them as M-values
# and map them back to betas (beta = 2^M / (2^M + 1)), as the original package did.
.betas_from_mvalues <- function(m) {
    rng <- suppressWarnings(range(m, na.rm = TRUE))
    if (all(is.finite(rng)) && (rng[1] < -0.1 || rng[2] > 1.1)) {
        message("methylclock: values fall outside [0, 1]; treating the input as ",
                "M-values and converting to betas.")
        m <- 2^m / (2^m + 1)
    }
    m
}
