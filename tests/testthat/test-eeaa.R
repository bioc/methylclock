# Canonical EEAA: the fixed-parameter blend. Checked on synthetic betas that
# carry the estimator CpGs, plus integrity checks of the shipped parameters.

.eeaa_fixture <- function(n = 40) {
    set.seed(7)
    p <- methylclock:::.mc_eeaa_params
    cpgs <- unique(c(names(p$plasmablast), names(p$cd8exhausted),
                     names(p$cd8naive)))
    betas <- matrix(runif(length(cpgs) * n, 0.1, 0.9),
                    nrow = length(cpgs),
                    dimnames = list(cpgs, sprintf("s%02d", seq_len(n))))
    age <- runif(n, 30, 80)
    x <- data.frame(id = colnames(betas),
                    Hannum = age + rnorm(n, 0, 4),
                    stringsAsFactors = FALSE)
    list(betas = betas, age = age, x = x)
}

test_that("the shipped EEAA parameters are the published ones", {
    p <- methylclock:::.mc_eeaa_params
    expect_length(p$plasmablast, 10L)
    expect_length(p$cd8exhausted, 35L)
    expect_length(p$cd8naive, 39L)
    expect_equal(unname(p$weights),
                 c(0.931238800, 0.009650902, 0.041428925, 0.017681373))
    expect_equal(unname(p$centers[1]), 10.398623)
    expect_equal(unname(p$scales[4]), -1.137203262)
    # spot checks against the source tables
    expect_equal(unname(p$plasmablast["cg24735235"]), -1.22037)
    expect_equal(unname(p$cd8exhausted["cg13608166"]), 51.38004579)
    expect_equal(unname(p$cd8naive["cg10493055"]), 410.2807)
})

test_that("EEAA returns a per-sample residual on age", {
    f <- .eeaa_fixture()
    res <- EEAA(f$x, age = f$age, betas = f$betas)
    expect_named(res, c("id", "EEAA"))
    expect_equal(res$id, f$x$id)
    expect_lt(abs(mean(res$EEAA)), 1e-6)
    expect_gt(stats::sd(res$EEAA), 0)
    # the residual is orthogonal to age
    expect_lt(abs(stats::cor(res$EEAA, f$age)), 1e-6)
})

test_that("a constant shift of an estimator's CpGs leaves EEAA unchanged", {
    # the estimators' unpublished intercepts shift every sample equally, and
    # any such constant must be absorbed by the regression on age
    f <- .eeaa_fixture()
    r1 <- EEAA(f$x, age = f$age, betas = f$betas)
    p <- methylclock:::.mc_eeaa_params
    shifted <- f$betas
    shifted[names(p$plasmablast), ] <- shifted[names(p$plasmablast), ] + 0.05
    r2 <- EEAA(f$x, age = f$age, betas = shifted)
    expect_equal(r1$EEAA, r2$EEAA, tolerance = 1e-8)
})

test_that("samples are matched by column name when betas are shuffled", {
    f <- .eeaa_fixture()
    r1 <- EEAA(f$x, age = f$age, betas = f$betas)
    r2 <- EEAA(f$x, age = f$age,
               betas = f$betas[, rev(colnames(f$betas))])
    expect_equal(r1$EEAA, r2$EEAA)
})

test_that("absent estimator CpGs warn, and too few stop", {
    f <- .eeaa_fixture()
    p <- methylclock:::.mc_eeaa_params
    drop_one <- f$betas[setdiff(rownames(f$betas),
                                names(p$plasmablast)[1]), ]
    expect_warning(EEAA(f$x, age = f$age, betas = drop_one),
                   "plasmablast")
    drop_most <- f$betas[setdiff(rownames(f$betas),
                                 names(p$plasmablast)[1:6]), ]
    expect_error(suppressWarnings(
        EEAA(f$x, age = f$age, betas = drop_most)),
        "Too few CpGs")
})

test_that("scattered NAs are tolerated via per-CpG means", {
    f <- .eeaa_fixture()
    holed <- f$betas
    holed[1, 3] <- NA
    holed[20, 10] <- NA
    res <- EEAA(f$x, age = f$age, betas = holed)
    expect_true(all(is.finite(res$EEAA)))
})
