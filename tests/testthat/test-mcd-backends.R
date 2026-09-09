# Backend chain: fall-through behaviour (offline, unit) and an end-to-end check
# that ExperimentHub returns the same object as the local mirror (skips offline).

test_that("backend chain is stored under methylclockData.backend", {
    old <- getOption("methylclockData.backend")
    on.exit(options(methylclockData.backend = old), add = TRUE)
    mcd_backends(c("local", "eh"))
    expect_identical(getOption("methylclockData.backend"), c("local", "eh"))
    expect_identical(mcd_backends(), c("local", "eh"))
})

test_that("a manifest field that is absent, NA or empty reads as not-set", {
    f <- methylclock:::.mcd_field
    expect_true(is.na(f(data.frame(id = "x"), "eh_id")))          # absent
    expect_true(is.na(f(data.frame(id = "x", eh_id = NA), "eh_id")))
    expect_true(is.na(f(data.frame(id = "x", eh_id = ""), "eh_id")))
    expect_identical(f(data.frame(id = "x", eh_id = "EH1"), "eh_id"), "EH1")
})

test_that("eh backend falls through when the row has no eh_id (no network)", {
    expect_null(methylclock:::.mcd_load_eh("x", data.frame(id = "x", eh_id = NA)))
    expect_null(methylclock:::.mcd_load_eh("x", data.frame(id = "x")))
})

test_that("zenodo backend falls through when the row has no doi (no network)", {
    expect_null(methylclock:::.mcd_load_zenodo("x",
        data.frame(id = "x", zenodo_doi = NA)))
    expect_null(methylclock:::.mcd_load_zenodo("x", data.frame(id = "x")))
})

test_that(".mcd_materialize loads a path but passes an object through", {
    obj <- data.frame(CpGmarker = "cg1", CoefficientTraining = 0.5)
    expect_identical(methylclock:::.mcd_materialize(obj, NULL), obj)
    tmp <- tempfile(fileext = ".rda")
    coefX <- obj; save(coefX, file = tmp)
    got <- methylclock:::.mcd_materialize(tmp, data.frame(object = "coefX"))
    expect_identical(got, obj)
    unlink(tmp)
})

# End-to-end: the eh backend must return exactly what local returns. Needs
# network + ExperimentHub; skips cleanly without them.
test_that("eh backend returns the same object as local", {
    skip_if_not(requireNamespace("ExperimentHub", quietly = TRUE),
                "ExperimentHub not installed")
    skip_if_offline()
    root <- tryCatch(methylclock:::.mcd_local_root(), error = function(e) NA)
    skip_if(is.na(root) || !dir.exists(root), "data-mirror not reachable")

    old <- getOption("methylclockData.backend")
    on.exit({ options(methylclockData.backend = old); mcd_cache_clear() }, add = TRUE)

    mcd_cache_clear(); mcd_backends("local"); a <- mcd_resource("coefHorvath")
    mcd_cache_clear(); mcd_backends("eh")
    b <- tryCatch(mcd_resource("coefHorvath"), error = function(e) NULL)
    skip_if(is.null(b), "ExperimentHub resource not retrievable")
    expect_identical(a, b)
})

# End-to-end: the zenodo backend serves users with no mirror and any R
# version (the hub's snapshot gate does not apply). Needs network.
test_that("zenodo backend returns the same objects as local", {
    skip_if_offline()
    root <- tryCatch(methylclock:::.mcd_local_root(), error = function(e) NA)
    skip_if(is.na(root) || !dir.exists(root), "data-mirror not reachable")

    old <- getOption("methylclockData.backend")
    on.exit({ options(methylclockData.backend = old); mcd_cache_clear() },
            add = TRUE)

    mcd_cache_clear(); mcd_backends("local"); a <- mcd_resource("coefLin")
    mcd_cache_clear(); mcd_backends("zenodo")
    b <- tryCatch(mcd_resource("coefLin"), error = function(e) NULL)
    skip_if(is.null(b), "Zenodo not retrievable")
    expect_identical(a, b)

    # file resources (the AltumAge weights) resolve through zenodo too
    p <- tryCatch(mcd_resource_file("coefAltumAge"), error = function(e) NULL)
    skip_if(is.null(p), "Zenodo file not retrievable")
    expect_true(file.exists(p) && file.size(p) > 1e6)
})
