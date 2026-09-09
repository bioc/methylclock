# Engine + API. These need the data mirror; they skip if it is unreachable.

has_mirror <- function() {
    root <- tryCatch(methylclock:::.mcd_local_root(), error = function(e) NA)
    !is.na(root) && dir.exists(root)
}

fixture <- function() {
    p <- file.path(methylclock:::.mcd_local_root(), "..", "oracles",
                   "fixtures", "synthetic.rds")
    skip_if_not(file.exists(p), "synthetic fixture not found")
    readRDS(p)
}

test_that("methylclock() computes a linear clock and returns an S3 result", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    res <- methylclock(fixture(), clocks = c("Levine", "Hannum"))
    expect_s3_class(res, "methylclock")
    df <- as.data.frame(res)
    expect_true(all(c("id", "Levine", "Hannum") %in% names(df)))
    expect_equal(nrow(df), 20L)
    expect_true(is.numeric(df$Levine))
})

test_that("DNAmGA restricts to gestational clocks", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    res <- suppressWarnings(DNAmGA(fixture()))
    expect_true(all(vapply(res$clocks, function(c)
        clock_info(c)$target == "gestational", logical(1))))
})

test_that("DNAmAge excludes gestational clocks", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    res <- suppressWarnings(DNAmAge(fixture()))
    expect_false(any(vapply(res$clocks, function(c)
        clock_info(c)$target == "gestational", logical(1))))
})

test_that("selection by native platform keeps only matching clocks", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    res <- suppressWarnings(methylclock(fixture(), platform = "EPIC"))
    expect_true(all(vapply(res$clocks, function(c)
        "EPIC" %in% clock_info(c)$platform, logical(1))))
})

test_that("selection by generation works", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    res <- suppressWarnings(DNAmAge(fixture(), generation = 1))
    expect_true(all(vapply(res$clocks, function(c)
        clock_info(c)$generation == 1L, logical(1))))
})

test_that("a filtered-in clock whose engine is unavailable is skipped, not fatal", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    # a throwaway clock on an as-yet unimplemented engine: a plain filter that
    # pulls it in must skip it (with a note), not abort the whole run.
    clock_register("FakeUnimplPC", target = "chronological", generation = 1L,
                   predictor = "pc")
    on.exit(suppressWarnings(rm("FakeUnimplPC",
                                envir = methylclock:::.clock_env)), add = TRUE)
    res <- suppressWarnings(DNAmAge(fixture(), generation = 1))
    expect_false("FakeUnimplPC" %in% res$clocks)
    expect_true(length(res$meta$skipped) >= 1L)
    expect_true("BNN" %in% res$clocks)              # the real nn clock computes
})

test_that("an explicitly named clock with an unavailable engine errors clearly", {
    expect_error(
        methylclock:::.predict_clock(
            list(name = "FakePC", predictor = "pc"), NULL), "not implemented")
})

test_that("DunedinPACE (quantile-normalized) reproduces the oracle", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    skip_if_not_installed("preprocessCore")
    op <- file.path(methylclock:::.mcd_local_root(), "..", "oracles",
                    "github", "oracle_synthetic.rds")
    skip_if_not(file.exists(op), "oracle not found")
    oracle <- readRDS(op)$DNAmAge$DUNEDIN
    got <- as.data.frame(methylclock(fixture(), clocks = "DunedinPACE"))$DunedinPACE
    expect_equal(got, oracle, tolerance = 1e-8)
})

test_that("hdf5 storage produces the same estimates as memory", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    skip_if_not_installed("BigDataStatMeth")
    f <- tempfile(fileext = ".h5")
    on.exit(unlink(f), add = TRUE)
    mem <- methylclock(fixture(), clocks = "Levine", storage = "memory")
    h5 <- methylclock(fixture(), clocks = "Levine", storage = "hdf5", file = f)
    expect_identical(h5$storage, "hdf5")
    expect_equal(as.data.frame(h5), as.data.frame(mem), tolerance = 1e-12)
})
