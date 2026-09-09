# AltumAge: a neural-network clock (predictor "nn"). Its forward pass runs in
# C++/Eigen with the layer weights read from HDF5 via BigDataStatMeth. The
# reference is biolearn's published output on GSE41169 samples; the fixture holds
# the input subset + expected ages. Needs the data mirror (HDF5 weights + ref).

has_altumage <- function() {
    root <- tryCatch(methylclock:::.mcd_local_root(), error = function(e) NA)
    if (is.na(root) || !dir.exists(root)) return(FALSE)
    file.exists(file.path(root, "altumage", "coefAltumAge.hdf5")) &&
        file.exists(file.path(root, "rda", "coefAltumAgeRef.rda"))
}
fixture_path <- function()
    testthat::test_path("..", "..", "..", "oracles", "fixtures",
                        "altumage_fixture.rds")

test_that("AltumAge is registered as an nn clock with a file + ref resource", {
    e <- clock_info("AltumAge")
    expect_identical(e$predictor, "nn")
    expect_identical(e$resource, "coefAltumAge")
    expect_identical(e$aux, "coefAltumAgeRef")
})

test_that("the weight bundle resolves to an HDF5 path and the reference loads", {
    skip_if_not(has_altumage(), "AltumAge data not in the mirror")
    mcd_cache_clear()
    p <- mcd_resource_file("coefAltumAge")
    expect_true(file.exists(p))
    expect_match(p, "\\.hdf5$")
    expect_equal(nrow(mcd_resource("coefAltumAgeRef")), 20318L)
})

test_that("methylclock reproduces biolearn's AltumAge output end to end", {
    skip_if_not(has_altumage(), "AltumAge data not in the mirror")
    skip_if_not(file.exists(fixture_path()), "AltumAge fixture not present")
    fx <- readRDS(fixture_path())
    input <- data.frame(ProbeID = rownames(fx$betas), fx$betas,
                        check.names = FALSE)
    res <- as.data.frame(methylclock(input, clocks = "AltumAge"))
    expect_equal(res$AltumAge, unname(fx$expected), tolerance = 1e-3)
})
