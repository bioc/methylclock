# T-07: imputation strategies and per-sample coverage floor. Uses a small
# in-memory beta matrix so it needs no data mirror for the matrix itself, but
# the clock coefficients come from the mirror; skip when it is unreachable.

mirror_ok <- function() {
    root <- tryCatch(methylclock:::.mcd_local_root(), error = function(e) NA)
    !is.na(root) && dir.exists(root)
}

# A beta matrix covering a clock's CpGs, with reproducible values in [0, 1].
clock_betas <- function(clock = "Horvath", nsamp = 8) {
    co <- mcd_resource(methylclock:::clock_registry()[[clock]]$resource)
    cpgs <- as.character(co$CpGmarker)
    cpgs <- setdiff(cpgs, c("(Intercept)", "Intercept", "intercept"))
    m <- matrix(seq(0.1, 0.9, length.out = length(cpgs) * nsamp),
                nrow = length(cpgs),
                dimnames = list(cpgs, paste0("s", seq_len(nsamp))))
    m
}

test_that("with no missing values the three strategies agree", {
    skip_if_not(mirror_ok(), "data-mirror not reachable")
    m <- clock_betas()
    base <- as.data.frame(suppressWarnings(
        methylclock(m, clocks = "Horvath")))$Horvath
    ref <- as.data.frame(suppressWarnings(
        methylclock(m, clocks = "Horvath", impute = "reference")))$Horvath
    none <- as.data.frame(suppressWarnings(
        methylclock(m, clocks = "Horvath", impute = "none")))$Horvath
    expect_equal(base, ref)
    expect_equal(base, none)
})

test_that("impute='none' returns NA only for samples missing a used CpG", {
    skip_if_not(mirror_ok(), "data-mirror not reachable")
    m <- clock_betas()
    bad <- c(2L, 5L)
    m[2, bad] <- NA                                  # one used CpG, two samples
    mean_v <- as.data.frame(suppressWarnings(
        methylclock(m, clocks = "Horvath", impute = "mean")))$Horvath
    none_v <- as.data.frame(suppressWarnings(
        methylclock(m, clocks = "Horvath", impute = "none")))$Horvath
    expect_equal(which(is.na(none_v)), bad)
    expect_false(any(is.na(mean_v)))                 # mean imputes, no NA
    expect_equal(mean_v[-bad], none_v[-bad])         # untouched samples agree
})

test_that("impute='reference' differs from 'mean' on the imputed cells", {
    skip_if_not(mirror_ok(), "data-mirror not reachable")
    m <- clock_betas()
    bad <- c(1L, 3L)
    m[2, bad] <- NA
    mean_v <- as.data.frame(suppressWarnings(
        methylclock(m, clocks = "Horvath", impute = "mean")))$Horvath
    ref_v <- as.data.frame(suppressWarnings(
        methylclock(m, clocks = "Horvath", impute = "reference")))$Horvath
    expect_false(any(is.na(ref_v)))
    expect_false(isTRUE(all.equal(mean_v[bad], ref_v[bad])))
    expect_equal(mean_v[-bad], ref_v[-bad])          # unaffected samples agree
})

test_that("min.perc.sample drops under-covered samples to NA", {
    skip_if_not(mirror_ok(), "data-mirror not reachable")
    m <- clock_betas()
    half <- seq_len(floor(nrow(m) / 2))
    m[half, 4] <- NA                                 # sample 4: ~50% missing
    v <- as.data.frame(suppressWarnings(
        methylclock(m, clocks = "Horvath", min.perc.sample = 0.8)))$Horvath
    expect_true(is.na(v[4]))
    expect_equal(sum(is.na(v)), 1L)                  # only that sample
})
