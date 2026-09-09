# Neural-network predictor engine. Two kinds of network:
#   - weights in an HDF5 bundle (`resource`) with a CpG-order reference (`aux`):
#     the forward pass runs in C++/Eigen reading the weights from HDF5 (AltumAge);
#   - weights compiled into the package C++ (`source == "compiled"`, no
#     resource): a baked-in forward over a fixed feature set (BNN).
.predict_nn <- function(entry, betas, min.perc = 0.8) {
    if (identical(entry$source, "compiled"))
        return(.predict_nn_compiled(entry, betas, min.perc))
    if (is.null(entry$resource) || is.null(entry$aux))
        stop(sprintf("Predictor 'nn' for clock '%s' is not implemented in this ",
                     entry$name), "version.")
    ref <- as.character(mcd_resource(entry$aux)$CpGmarker)
    h5  <- mcd_resource_file(entry$resource)

    samples <- .src_colnames(betas)
    present <- ref %in% .src_rownames(betas)
    coverage <- mean(present)
    if (coverage <= min.perc) {
        warning(sprintf(
            "clock '%s': %.0f%% of CpGs present (need > %.0f%%); returning NA.",
            entry$name, 100 * coverage, 100 * min.perc))
        out <- rep(NA_real_, length(samples))
        names(out) <- samples
        return(out)
    }
    if (!all(present))
        warning(sprintf("clock '%s': %d of %d CpGs absent; filled with 0.",
                        entry$name, sum(!present), length(ref)))

    X <- matrix(0, nrow = length(ref), ncol = length(samples),
                dimnames = list(ref, samples))
    m <- ref[present]
    X[m, ] <- .src_subset(betas, m)
    X[is.na(X)] <- 0
    pred <- altumage_forward(h5, X)                 # C++/Eigen forward, HDF5 weights
    names(pred) <- samples
    entry$transform(pred)
}

# Compiled neural network (BNN): the weights are baked into the package C++, not a
# data resource. Its fixed feature set is the reference CpGs named by `aux`, in
# that order; the compiled forward needs all of them present, so a clock missing
# any feature CpG returns NA. Missing values among the present CpGs are imputed to
# each CpG's mean across samples.
.predict_nn_compiled <- function(entry, betas, min.perc = 0.8) {
    ref <- as.character(mcd_resource(entry$aux)$CpGmarker)
    ref <- ref[!ref %in% c("(Intercept)", "Intercept", "intercept")]
    samples <- .src_colnames(betas)
    present <- ref %in% .src_rownames(betas)
    if (!all(present)) {
        warning(sprintf(
            "clock '%s': %d of %d CpGs absent; the compiled model needs all of ",
            entry$name, sum(!present), length(ref)), "them; returning NA.",
            call. = FALSE)
        out <- rep(NA_real_, length(samples))
        names(out) <- samples
        return(out)
    }
    X <- .impute_rows(.src_subset(betas, ref))      # CpGs x samples, ref order
    pred <- bnn_predict(X)                          # baked-in C++ forward
    names(pred) <- samples
    entry$transform(pred)
}
