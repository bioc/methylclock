# Predictor dispatch and the coefficient-table helpers shared by the engines.
# Each engine lives in its own file: predict_linear.R, predict_counter.R,
# predict_nn.R. `.predict_clock` routes a clock to its engine by `predictor`.

# Locate the CpG-identifier column in a coefficient table.
.coef_cpg_col <- function(coef) {
    nm <- names(coef)
    hit <- nm[tolower(nm) %in% c("cpgmarker", "cpg", "probeid", "id")]
    if (length(hit)) return(hit[1])
    chr <- nm[vapply(coef, is.character, logical(1))]
    if (length(chr)) return(chr[1])
    stop("No CpG-identifier column found in coefficient table.")
}

# Locate the weight column: the one named by the registry, else a default.
.coef_weight_col <- function(coef, coef_col) {
    if (!is.null(coef_col) && coef_col %in% names(coef))
        return(coef_col)
    cand <- names(coef)[tolower(names(coef)) %in%
                        c("coefficienttraining", "coefficient", "weight", "beta")]
    if (length(cand)) return(cand[1])
    num <- names(coef)[vapply(coef, is.numeric, logical(1))]
    if (length(num)) return(num[1])
    stop("No coefficient column found in coefficient table.")
}

# Dispatch one clock to its predictor family. `betas` is a beta source.
.predict_clock <- function(entry, betas, covars = NULL, min.perc = 0.8,
                           impute = "mean", min.perc.sample = 0,
                           knn = .MC_KNN) {
    switch(entry$predictor,
        linear = {
            src <- betas
            if (identical(entry$normalize, "quantile_gs"))
                src <- .as_source_from_matrix(
                    .quantile_normalize_gs(betas, entry$aux))
            .predict_linear(entry, src, min.perc = min.perc, impute = impute,
                            min.perc.sample = min.perc.sample, knn = knn)
        },
        counter = .predict_counter(entry, betas, min.perc = min.perc),
        nn = .predict_nn(entry, betas, min.perc = min.perc),
        stop(sprintf("Predictor '%s' (clock '%s') is not implemented in this ",
                     entry$predictor, entry$name),
             "version.")
    )
}
