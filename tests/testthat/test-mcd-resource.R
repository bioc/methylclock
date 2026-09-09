# Data resolver against the local backend (data-mirror). These tests need the
# data mirror; they skip cleanly if it is not reachable from the test wd.

has_mirror <- function() {
    root <- tryCatch(methylclock:::.mcd_local_root(), error = function(e) NA)
    !is.na(root) && dir.exists(root)
}

test_that("manifest loads and lists known resources", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    man <- mcd_manifest()
    expect_true(all(c("id", "object", "sha256") %in% names(man)))
    expect_true("coefHorvath" %in% man$id)
})

test_that("mcd_resource resolves a coefficient object via local backend", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    mcd_cache_clear()
    coef <- mcd_resource("coefHorvath")
    expect_s3_class(coef, "data.frame")
    expect_true(nrow(coef) > 300)
})

test_that("mcd_resource caches (second call is the same object)", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    a <- mcd_resource("coefHannum")
    b <- mcd_resource("coefHannum")
    expect_identical(a, b)
})

test_that("unknown resource errors clearly", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    expect_error(mcd_resource("coefNope"), "Unknown resource")
})

test_that("every built-in clock resource exists in the manifest", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    man_ids <- mcd_manifest()$id
    ids <- unlist(lapply(clock_registry(), function(e) c(e$resource, e$aux)))
    ids <- unique(ids[!is.na(ids)])
    expect_true(all(ids %in% man_ids),
                info = paste("missing:", paste(setdiff(ids, man_ids), collapse = ", ")))
})
