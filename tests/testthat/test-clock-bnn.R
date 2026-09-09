# BNN: a neural-network clock whose weights are compiled into the package C++
# (not a data resource). Its features are the Horvath CpGs (aux), fed as a
# CpGs x samples matrix to the baked-in forward. Validated against oracle B, the
# original package's output on the shared fixture. Needs the data mirror.

bnn_root <- function() {
    root <- tryCatch(methylclock:::.mcd_local_root(), error = function(e) NA)
    if (is.na(root) || !dir.exists(root)) return(NA_character_)
    root
}

test_that("BNN is registered as a compiled nn clock over the Horvath CpGs", {
    reg <- clock_registry()
    expect_true("BNN" %in% names(reg))
    expect_identical(reg[["BNN"]]$predictor, "nn")
    expect_identical(reg[["BNN"]]$source, "compiled")
    expect_null(reg[["BNN"]]$resource)
    expect_identical(reg[["BNN"]]$aux, "coefHorvath")
})

test_that("BNN reproduces oracle B (the original package) on the fixture", {
    root <- bnn_root()
    skip_if_not(!is.na(root), "data-mirror not reachable")
    orc <- file.path(root, "..", "oracles", "github", "oracle_synthetic.rds")
    fix <- file.path(root, "..", "oracles", "fixtures", "synthetic.rds")
    skip_if_not(file.exists(orc) && file.exists(fix), "oracle/fixture not found")
    mcd_cache_clear()
    res <- as.data.frame(suppressWarnings(
        methylclock(readRDS(fix), clocks = "BNN")))
    oracle <- readRDS(orc)$DNAmAge$BNN
    expect_equal(res$BNN, oracle, tolerance = 1e-9)
})

test_that("BNN returns NA when a feature CpG is absent (compiled needs all)", {
    root <- bnn_root()
    skip_if_not(!is.na(root), "data-mirror not reachable")
    fix <- file.path(root, "..", "oracles", "fixtures", "synthetic.rds")
    skip_if_not(file.exists(fix), "fixture not found")
    mcd_cache_clear()
    m <- methylclock:::.as_beta_matrix(readRDS(fix))
    hv <- mcd_resource("coefHorvath")$CpGmarker
    drop <- hv[hv != "(Intercept)"][1]                # remove one feature CpG
    m2 <- m[rownames(m) != drop, , drop = FALSE]
    res <- as.data.frame(suppressWarnings(
        methylclock(m2, clocks = "BNN")))
    expect_true(all(is.na(res$BNN)))
})
