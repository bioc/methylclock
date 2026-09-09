# Input handling: turn a user request into the clocks to run and a clean
# methylation matrix to run them on.

# Targets that are not aging clocks (trait/exposure and disease-risk predictors).
# They are opt-in: a plain clocks = "all" (with no target filter) does not include
# them, keeping the default estimation age-focused. They are still reachable by
# name or with an explicit target = "trait" / "risk".
.MC_OPTIN_TARGETS <- c("trait", "risk")

# Resolve a clock selection to registry entries. `clocks` is "all" or a vector
# of names; the target/generation/platform filters narrow an "all" selection.
# Platform matching uses each clock's native platform only (conversion to
# non-native platforms is a separate, later feature).
.select_clocks <- function(clocks = "all", target = NULL,
                           generation = NULL, platform = NULL) {
    reg <- clock_registry()
    keep <- function(e) {
        t_ok <- is.null(target) ||
            (identical(target, "!gestational") && e$target != "gestational") ||
            (!identical(target, "!gestational") && e$target == target)
        g_ok <- is.null(generation) || (e$generation %in% generation)
        p_ok <- is.null(platform) || any(platform %in% e$platform)
        t_ok && g_ok && p_ok
    }
    reg <- reg[vapply(reg, keep, logical(1))]
    if (length(clocks) == 1L && identical(clocks, "all")) {
        # a plain "all" is age-focused: drop opt-in targets unless the caller
        # asked for one explicitly through `target`.
        if (is.null(target))
            reg <- reg[vapply(reg,
                function(e) !e$target %in% .MC_OPTIN_TARGETS, logical(1))]
        return(reg)
    }
    miss <- setdiff(clocks, names(reg))
    if (length(miss))
        stop("Unknown or filtered-out clock(s): ", paste(miss, collapse = ", "),
             ". See clock_list().")
    reg[clocks]
}

# Quantile-normalize each sample to a gold-standard reference distribution.
# `aux_id` resolves to a table with a CpG column and a reference-means column;
# only those CpGs are read from the source and each sample (column) is mapped to
# the reference distribution. Returns a dense CpGs x samples matrix.
.quantile_normalize_gs <- function(betas, aux_id) {
    gs <- mcd_resource(aux_id)
    gs_cpg <- as.character(gs[[.coef_cpg_col(gs)]])
    means_col <- names(gs)[tolower(names(gs)) %in%
                           c("meansgs", "means", "target")]
    if (!length(means_col))
        stop("Reference means column not found in '", aux_id, "'.")
    target <- as.numeric(gs[[means_col[1]]])
    obs <- intersect(gs_cpg, .src_rownames(betas))
    X <- .src_subset(betas, obs)                   # reads only the aux CpGs
    Xn <- preprocessCore::normalize.quantiles.use.target(X, target = target)
    dimnames(Xn) <- dimnames(X)
    Xn
}

# Coerce input to a CpG x sample matrix with CpG row names.
.as_beta_matrix <- function(x) {
    if (is.data.frame(x)) {
        chr <- vapply(x, is.character, logical(1))
        if (any(chr)) {
            rn <- as.character(x[[which(chr)[1]]])
            x <- as.matrix(x[, !chr, drop = FALSE])
            rownames(x) <- rn
        } else {
            x <- as.matrix(x)
        }
    }
    if (is.null(rownames(x)))
        stop("The methylation matrix needs CpG identifiers as row names.")
    x
}
