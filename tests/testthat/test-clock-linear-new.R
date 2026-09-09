# Linear clocks added from external weight sources: CausAge/DamAge/AdaptAge
# (biolearn, Ying 2024) and Lin (2016 supplement). Each is a direct linear model
# (intercept + weighted sum), so the reference is exact arithmetic on the stored
# coefficients. Needs the data mirror; skips cleanly without it.

has_mirror <- function() {
    root <- tryCatch(methylclock:::.mcd_local_root(), error = function(e) NA)
    !is.na(root) && dir.exists(root)
}

new_linear <- c("CausAge", "DamAge", "AdaptAge", "Lin", "Weidner")

test_that("the new linear clocks are registered", {
    reg <- clock_registry()
    for (nm in new_linear) {
        expect_true(nm %in% names(reg))
        expect_identical(reg[[nm]]$predictor, "linear")
    }
    expect_identical(clock_info("CausAge")$target, "causal")
})

test_that("their resources carry an intercept and the expected CpG count", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    mcd_cache_clear()
    counts <- c(CausAge = 585L, DamAge = 1089L, AdaptAge = 999L, Lin = 99L,
                Weidner = 3L)
    for (nm in new_linear) {
        coef <- mcd_resource(paste0("coef", nm))
        ic <- grepl("intercept", coef$CpGmarker, ignore.case = TRUE)
        expect_equal(sum(ic), 1L, info = nm)
        expect_equal(sum(!ic), counts[[nm]], info = nm)
    }
})

test_that("the engine reproduces the hand-computed weighted sum", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    mcd_cache_clear()
    for (nm in new_linear) {
        coef <- mcd_resource(paste0("coef", nm))
        ic <- grepl("intercept", coef$CpGmarker, ignore.case = TRUE)
        b0 <- coef$CoefficientTraining[ic]
        cpgs <- coef$CpGmarker[!ic]
        w <- coef$CoefficientTraining[!ic]
        set.seed(11)
        B <- matrix(runif(length(cpgs) * 3), nrow = length(cpgs),
                    dimnames = list(cpgs, paste0("S", 1:3)))
        manual <- as.vector(w %*% B) + b0
        input <- data.frame(ProbeID = cpgs, B, check.names = FALSE)
        mine <- as.data.frame(methylclock(input, clocks = nm, min.perc = 0))[[nm]]
        expect_equal(mine, manual, tolerance = 1e-9, info = nm)
    }
})
