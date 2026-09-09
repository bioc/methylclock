# Vidal-Bralo (2016): an 8-CpG whole-blood linear clock. As a direct linear
# regression its reference is exact arithmetic, so we check the engine against a
# hand-computed weighted sum. Needs the data mirror; skips cleanly without it.

has_mirror <- function() {
    root <- tryCatch(methylclock:::.mcd_local_root(), error = function(e) NA)
    !is.na(root) && dir.exists(root)
}

vb_cpgs <- c("cg16386080", "cg24768561", "cg19761273", "cg25809905",
             "cg09809672", "cg02228185", "cg17471102", "cg10917602")
vb_coef <- c(59.5, 33.9, -44.0, -19.7, -22.8, -16.8, -17.7, -11.4)
vb_intercept <- 84.7

test_that("VidalBralo is registered as a 1st-generation linear clock", {
    reg <- clock_registry()
    expect_true("VidalBralo" %in% names(reg))
    e <- reg[["VidalBralo"]]
    expect_identical(e$predictor, "linear")
    expect_identical(e$resource, "coefVidalBralo")
    expect_equal(e$generation, 1L)
})

test_that("its coefficient resource has the intercept plus 8 CpGs", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    mcd_cache_clear()
    coef <- mcd_resource("coefVidalBralo")
    expect_true(all(vb_cpgs %in% coef$CpGmarker))
    expect_true(any(grepl("[Ii]ntercept", coef$CpGmarker)))
    expect_equal(nrow(coef), length(vb_cpgs) + 1L)
})

test_that("the engine reproduces the hand-computed weighted sum", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    set.seed(42)
    B <- matrix(runif(length(vb_cpgs) * 4, 0.1, 0.9),
                nrow = length(vb_cpgs),
                dimnames = list(vb_cpgs, paste0("S", 1:4)))
    input <- data.frame(ProbeID = vb_cpgs, B, check.names = FALSE)
    manual <- vb_intercept + as.vector(vb_coef %*% B)
    res <- as.data.frame(methylclock(input, clocks = "VidalBralo"))
    expect_equal(res$VidalBralo, manual, tolerance = 1e-10)
})
