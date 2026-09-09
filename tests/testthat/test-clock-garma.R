# The Garma 2024 cross-platform clock: a linear model whose table also
# weights the SQUARE of some probes' betas (power column). The reference is
# exact arithmetic on the stored coefficients, squared terms included.
# Needs the data mirror; skips cleanly without it.

has_mirror <- function() {
    root <- tryCatch(methylclock:::.mcd_local_root(), error = function(e) NA)
    !is.na(root) && dir.exists(root)
}

test_that("Garma is registered as a native cross-platform linear clock", {
    e <- clock_info("Garma")
    expect_identical(e$predictor, "linear")
    expect_identical(e$target, "chronological")
    expect_true(all(c("450K", "EPIC", "EPICv2") %in% e$platform))
})

test_that("its resource carries the intercept and the squared-term rows", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    mcd_cache_clear()
    coef <- mcd_resource("coefGarma")
    expect_identical(nrow(coef), 4963L)
    expect_identical(coef$CpGmarker[1], "(Intercept)")
    expect_identical(sum(coef$power == 2), 168L)
    expect_false(anyNA(coef$CoefficientTraining))
})

test_that("the engine equals exact arithmetic, squared terms included", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    coef <- mcd_resource("coefGarma")
    b0 <- coef$CoefficientTraining[1]
    tab <- coef[-1, ]
    cpgs <- unique(tab$CpGmarker)
    set.seed(11)
    betas <- matrix(runif(length(cpgs) * 8, 0.05, 0.95),
                    nrow = length(cpgs),
                    dimnames = list(cpgs, paste0("s", 1:8)))
    res <- suppressWarnings(methylclock(betas, clocks = "Garma"))
    got <- as.data.frame(res)$Garma

    lin <- tab[tab$power == 1, ]
    quad <- tab[tab$power == 2, ]
    manual <- b0 +
        as.numeric(crossprod(lin$CoefficientTraining,
                             betas[lin$CpGmarker, ])) +
        as.numeric(crossprod(quad$CoefficientTraining,
                             betas[quad$CpGmarker, ]^2))
    expect_equal(got, unname(manual), tolerance = 1e-12)

    # the squared block genuinely contributes: dropping it must change the
    # result (guards against the power column being silently ignored)
    wrong <- b0 + as.numeric(crossprod(lin$CoefficientTraining,
                                       betas[lin$CpGmarker, ]))
    expect_gt(max(abs(got - wrong)), 0.1)
})

test_that("an HDF5 input equals RAM, partial coverage included", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    skip_if_not_installed("BigDataStatMeth")
    coef <- mcd_resource("coefGarma")
    cpgs <- unique(coef$CpGmarker[-1])
    set.seed(13)
    # partial coverage plus filler rows: the shape that routes an HDF5 input
    # through the batched algebra, which must hand this clock back to the
    # per-clock path (its squared terms are not a plain cross-product)
    keep <- sample(cpgs, floor(length(cpgs) * 0.9))
    rows <- c(keep, paste0("cg9999", seq_len(2000)))
    betas <- matrix(runif(length(rows) * 5, 0.05, 0.95),
                    nrow = length(rows),
                    dimnames = list(rows, paste0("s", 1:5)))
    ram <- suppressWarnings(methylclock(betas, clocks = "Garma"))
    f <- tempfile(fileext = ".h5")
    h <- BigDataStatMeth::hdf5_create_matrix(
        filename = f, dataset = "in/betas", data = betas, overwrite = TRUE)
    close(h)
    hm <- BigDataStatMeth::hdf5_matrix(f, "in/betas")
    on.exit({ close(hm); unlink(f) }, add = TRUE)
    dsk <- suppressWarnings(methylclock(hm, clocks = "Garma"))
    expect_equal(as.data.frame(dsk)$Garma, as.data.frame(ram)$Garma,
                 tolerance = 1e-12)
})

test_that("impute='none' with full coverage matches, and NA path behaves", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    coef <- mcd_resource("coefGarma")
    cpgs <- unique(coef$CpGmarker[-1])
    set.seed(12)
    betas <- matrix(runif(length(cpgs) * 4, 0.05, 0.95),
                    nrow = length(cpgs),
                    dimnames = list(cpgs, paste0("s", 1:4)))
    full <- suppressWarnings(methylclock(betas, clocks = "Garma"))
    none <- suppressWarnings(methylclock(betas, clocks = "Garma",
                                         impute = "none"))
    expect_equal(as.data.frame(full)$Garma, as.data.frame(none)$Garma,
                 tolerance = 1e-12)
    betas[cpgs[1], 2] <- NA          # a missing used CpG under "none" -> NA
    none2 <- suppressWarnings(methylclock(betas, clocks = "Garma",
                                          impute = "none"))
    expect_true(is.na(as.data.frame(none2)$Garma[2]))
    expect_false(anyNA(as.data.frame(none2)$Garma[-2]))
})
