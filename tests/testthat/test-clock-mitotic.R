# Mitotic counter clocks (epiTOC1, HypoClock, stemTOC). Their reference is the
# EpiMitClocks package: epiTOC1 = mean of its CpG set, HypoClock = 1 - that mean,
# stemTOC = the 0.95 quantile. Here we check the registration and the arithmetic
# against a synthetic input; the exact match against the reference implementation
# on its own example data is verified in the build/validation script (not shipped).

has_mirror <- function() {
    root <- tryCatch(methylclock:::.mcd_local_root(), error = function(e) NA)
    !is.na(root) && dir.exists(root)
}

test_that("the three mitotic clocks are registered as counters", {
    reg <- clock_registry()
    for (nm in c("epiTOC1", "HypoClock", "stemTOC")) {
        expect_true(nm %in% names(reg))
        expect_identical(reg[[nm]]$predictor, "counter")
        expect_identical(reg[[nm]]$target, "mitotic")
    }
    expect_identical(reg[["stemTOC"]]$agg, list(fun = "quantile", p = 0.95))
})

test_that("their CpG-set resources load with the expected sizes", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    mcd_cache_clear()
    expect_equal(nrow(mcd_resource("coefEpiTOC1")), 385L)
    expect_equal(nrow(mcd_resource("coefHypoClock")), 678L)
    expect_equal(nrow(mcd_resource("coefStemTOC")), 371L)
})

test_that("engine matches the reference summaries on a synthetic input", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    mcd_cache_clear()
    e1 <- mcd_resource("coefEpiTOC1")$CpGmarker
    eh <- mcd_resource("coefHypoClock")$CpGmarker
    es <- mcd_resource("coefStemTOC")$CpGmarker
    cpgs <- unique(c(e1, eh, es))
    set.seed(3)
    B <- matrix(runif(length(cpgs) * 5), nrow = length(cpgs),
                dimnames = list(cpgs, paste0("S", 1:5)))
    input <- data.frame(ProbeID = cpgs, B, check.names = FALSE)
    res <- as.data.frame(methylclock(
        input, clocks = c("epiTOC1", "HypoClock", "stemTOC"), min.perc = 0))

    expect_equal(res$epiTOC1, unname(colMeans(B[e1, ])), tolerance = 1e-12)
    expect_equal(res$HypoClock, unname(1 - colMeans(B[eh, ])), tolerance = 1e-12)
    expect_equal(res$stemTOC,
                 unname(apply(B[es, ], 2, stats::quantile, 0.95, names = FALSE)),
                 tolerance = 1e-12)
})

test_that("epiTOC2 is a parametric counter matching its closed-form score", {
    reg <- clock_registry()
    expect_true("epiTOC2" %in% names(reg))
    expect_identical(reg[["epiTOC2"]]$predictor, "counter")
    expect_identical(reg[["epiTOC2"]]$target, "mitotic")

    skip_if_not(has_mirror(), "data-mirror not reachable")
    mcd_cache_clear()
    res <- mcd_resource("coefEpiTOC2")
    expect_equal(nrow(res), 163L)
    expect_true(all(c("delta", "beta0") %in% names(res)))
    cpgs <- res$CpGmarker
    set.seed(6)
    B <- matrix(runif(length(cpgs) * 4), nrow = length(cpgs),
                dimnames = list(cpgs, paste0("S", 1:4)))
    # tnsc = 2 * mean_i[(beta - beta0) / (delta * (1 - beta0))]
    manual <- 2 * colMeans((B - res$beta0) / (res$delta * (1 - res$beta0)),
                           na.rm = TRUE)
    input <- data.frame(ProbeID = cpgs, B, check.names = FALSE)
    mine <- as.data.frame(
        methylclock(input, clocks = "epiTOC2", min.perc = 0))$epiTOC2
    expect_equal(mine, unname(manual), tolerance = 1e-12)
})

test_that("RepliTali is a linear mitotic clock matching its weighted sum", {
    reg <- clock_registry()
    expect_true("RepliTali" %in% names(reg))
    expect_identical(reg[["RepliTali"]]$predictor, "linear")
    expect_identical(reg[["RepliTali"]]$target, "mitotic")

    skip_if_not(has_mirror(), "data-mirror not reachable")
    mcd_cache_clear()
    coef <- mcd_resource("coefRepliTali")
    ic <- grepl("intercept", coef$CpGmarker, ignore.case = TRUE)
    expect_equal(sum(ic), 1L)
    expect_equal(sum(!ic), 87L)
    b0 <- coef$CoefficientTraining[ic]
    cpgs <- coef$CpGmarker[!ic]
    w <- coef$CoefficientTraining[!ic]
    set.seed(4)
    B <- matrix(runif(length(cpgs) * 3), nrow = length(cpgs),
                dimnames = list(cpgs, paste0("S", 1:3)))
    manual <- as.vector(w %*% B) + b0
    input <- data.frame(ProbeID = cpgs, B, check.names = FALSE)
    mine <- as.data.frame(
        methylclock(input, clocks = "RepliTali", min.perc = 0))$RepliTali
    expect_equal(mine, manual, tolerance = 1e-9)
})
