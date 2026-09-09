# Canonical IEAA: Horvath residual on age plus the seven Chen covariates
# (three rare subsets from the published CpG estimators, four common cell
# types supplied or deconvolved). Synthetic betas carry the estimator CpGs.

.ieaa_fixture <- function(n = 40) {
    set.seed(11)
    p <- methylclock:::.mc_eeaa_params
    cpgs <- unique(c(names(p$plasmablast), names(p$cd8exhausted),
                     names(p$cd8naive)))
    betas <- matrix(runif(length(cpgs) * n, 0.1, 0.9),
                    nrow = length(cpgs),
                    dimnames = list(cpgs, sprintf("s%02d", seq_len(n))))
    age <- runif(n, 30, 80)
    x <- data.frame(id = colnames(betas),
                    Horvath = age + rnorm(n, 0, 4),
                    stringsAsFactors = FALSE)
    cc <- matrix(runif(n * 5, 0, 0.4), n, 5,
                 dimnames = list(colnames(betas),
                                 c("CD4T", "NK", "Mono", "Neu", "Eos")))
    list(betas = betas, age = age, x = x, cc = cc)
}

test_that("IEAA returns a residual orthogonal to age and the covariates", {
    f <- .ieaa_fixture()
    res <- IEAA(f$x, age = f$age, betas = f$betas, cell_counts = f$cc)
    expect_named(res, c("id", "IEAA"))
    expect_equal(res$id, f$x$id)
    expect_lt(abs(mean(res$IEAA)), 1e-6)
    expect_gt(stats::sd(res$IEAA), 0)
    expect_lt(abs(stats::cor(res$IEAA, f$age)), 1e-6)
    expect_lt(abs(stats::cor(res$IEAA, f$cc[, "CD4T"])), 1e-6)
})

test_that("a Gran column and Neu+Eos are interchangeable", {
    f <- .ieaa_fixture()
    r1 <- IEAA(f$x, age = f$age, betas = f$betas, cell_counts = f$cc)
    cc2 <- cbind(f$cc[, c("CD4T", "NK", "Mono"), drop = FALSE],
                 Gran = f$cc[, "Neu"] + f$cc[, "Eos"])
    r2 <- IEAA(f$x, age = f$age, betas = f$betas, cell_counts = cc2)
    expect_equal(r1$IEAA, r2$IEAA, tolerance = 1e-10)
})

test_that("missing common covariates are an error", {
    f <- .ieaa_fixture()
    expect_error(
        IEAA(f$x, age = f$age, betas = f$betas,
             cell_counts = f$cc[, c("CD4T", "NK"), drop = FALSE]),
        "'Mono'")
    expect_error(
        IEAA(f$x, age = f$age, betas = f$betas,
             cell_counts = f$cc[, c("CD4T", "NK", "Mono"), drop = FALSE]),
        "Gran")
})

test_that("IEAA differs from the arbitrary-panel residualCells", {
    # the seven canonical covariates are not the generic panel adjustment
    f <- .ieaa_fixture()
    i <- IEAA(f$x, age = f$age, betas = f$betas, cell_counts = f$cc)
    aa <- ageAcceleration(f$x, age = f$age, cell_counts = f$cc,
                          clocks = "Horvath")
    expect_gt(max(abs(i$IEAA - aa$residualCells)), 1e-8)
})

test_that("samples are matched by column name when betas are shuffled", {
    f <- .ieaa_fixture()
    r1 <- IEAA(f$x, age = f$age, betas = f$betas, cell_counts = f$cc)
    r2 <- IEAA(f$x, age = f$age,
               betas = f$betas[, rev(colnames(f$betas))],
               cell_counts = f$cc)
    expect_equal(r1$IEAA, r2$IEAA)
})
