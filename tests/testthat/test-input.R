# Input handling: M-values are detected and converted, and the common
# Bioconductor containers are accepted as well as a bare matrix.

betas_small <- function() {
    set.seed(1)
    m <- matrix(runif(12, 0.05, 0.95), nrow = 4,
                dimnames = list(paste0("cg", 1:4), paste0("s", 1:3)))
    m
}

test_that("values outside [0,1] are treated as M-values and converted", {
    b <- betas_small()
    M <- log2(b / (1 - b))
    conv <- suppressMessages(methylclock:::.betas_from_mvalues(M))
    expect_equal(conv, b, tolerance = 1e-12)
    # a genuine beta matrix is left untouched
    expect_identical(methylclock:::.betas_from_mvalues(b), b)
})

test_that(".extract_betas pulls a matrix from Bioconductor containers", {
    b <- betas_small()
    expect_identical(methylclock:::.extract_betas(b), b)

    skip_if_not_installed("SummarizedExperiment")
    se <- SummarizedExperiment::SummarizedExperiment(assays = list(b))
    expect_equal(methylclock:::.extract_betas(se), b)

    skip_if_not_installed("Biobase")
    es <- Biobase::ExpressionSet(assayData = b)
    expect_equal(methylclock:::.extract_betas(es), b)
})

test_that("an unsupported input is rejected with a clear message", {
    expect_error(methylclock:::.extract_betas(list(1, 2)),
                 "Unsupported input")
})
