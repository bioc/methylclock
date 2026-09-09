# Trait/exposure EpiScores (McCartney 2018, via biolearn). Each is a pure linear
# weighted sum with no intercept, so the reference is exact arithmetic on the
# stored coefficients. They are opt-in: a plain clocks = "all" is age-focused and
# must not include them. Needs the data mirror; skips cleanly without it.

trait_has_mirror <- function() {
    root <- tryCatch(methylclock:::.mcd_local_root(), error = function(e) NA)
    !is.na(root) && dir.exists(root)
}

mccartney <- c("McCartney.Smoking", "McCartney.Alcohol", "McCartney.BMI",
               "McCartney.BodyFat", "McCartney.Education", "McCartney.HDL",
               "McCartney.LDL", "McCartney.TotalChol")

test_that("the McCartney EpiScores are registered as trait linear clocks", {
    reg <- clock_registry()
    for (nm in mccartney) {
        expect_true(nm %in% names(reg), info = nm)
        expect_identical(reg[[nm]]$predictor, "linear", info = nm)
        expect_identical(reg[[nm]]$target, "trait", info = nm)
        expect_false(isTRUE(reg[[nm]]$intercept), info = nm)
    }
})

test_that("trait clocks are opt-in: not in a default 'all', reachable otherwise", {
    all_default <- names(methylclock:::.select_clocks("all"))
    expect_false(any(mccartney %in% all_default))
    by_target <- names(methylclock:::.select_clocks("all", target = "trait"))
    expect_true(all(mccartney %in% by_target))
    expect_true("mCigarette" %in% by_target)
    by_name <- names(methylclock:::.select_clocks(c("McCartney.BMI", "Horvath")))
    expect_setequal(by_name, c("McCartney.BMI", "Horvath"))
})

test_that("mCigarette is a trait smoking clock, opt-in, matching its sum", {
    reg <- clock_registry()
    expect_true("mCigarette" %in% names(reg))
    expect_identical(reg[["mCigarette"]]$target, "trait")
    expect_false("mCigarette" %in% names(methylclock:::.select_clocks("all")))

    skip_if_not(trait_has_mirror(), "data-mirror not reachable")
    mcd_cache_clear()
    coef <- mcd_resource("coefMcCigarette")
    expect_false(any(grepl("intercept", coef$CpGmarker, ignore.case = TRUE)))
    cpgs <- coef$CpGmarker
    w <- coef$CoefficientTraining
    set.seed(9)
    B <- matrix(runif(length(cpgs) * 3), nrow = length(cpgs),
                dimnames = list(cpgs, paste0("S", 1:3)))
    manual <- as.vector(w %*% B)
    input <- data.frame(ProbeID = cpgs, B, check.names = FALSE)
    mine <- as.data.frame(
        methylclock(input, clocks = "mCigarette", min.perc = 0))$mCigarette
    expect_equal(mine, manual, tolerance = 1e-9)
})

test_that("the engine reproduces the hand-computed weighted sum (no intercept)", {
    skip_if_not(trait_has_mirror(), "data-mirror not reachable")
    mcd_cache_clear()
    for (nm in mccartney) {
        coef <- mcd_resource(clock_info(nm)$resource)
        expect_false(any(grepl("intercept", coef$CpGmarker, ignore.case = TRUE)),
                     info = nm)
        cpgs <- coef$CpGmarker
        w <- coef$CoefficientTraining
        set.seed(7)
        B <- matrix(runif(length(cpgs) * 3), nrow = length(cpgs),
                    dimnames = list(cpgs, paste0("S", 1:3)))
        manual <- as.vector(w %*% B)
        input <- data.frame(ProbeID = cpgs, B, check.names = FALSE)
        mine <- as.data.frame(
            methylclock(input, clocks = nm, min.perc = 0))[[nm]]
        expect_equal(mine, manual, tolerance = 1e-9, info = nm)
    }
})
