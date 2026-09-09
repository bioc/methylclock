# Coefficient contract normalization.

test_that("a bare CpG vector becomes a one-column CpGmarker table", {
    out <- normalize_coef(c("cg001", "cg002"))
    expect_s3_class(out, "data.frame")
    expect_identical(names(out), "CpGmarker")
    expect_type(out$CpGmarker, "character")
})

test_that("a tibble-like table is coerced to a plain data.frame", {
    tb <- structure(
        list(CpGmarker = c("cg1", "cg2"), CoefficientTraining = c(1.5, 2.5)),
        class = c("tbl_df", "tbl", "data.frame"), row.names = 1:2)
    out <- normalize_coef(tb)
    expect_identical(class(out), "data.frame")
    expect_type(out$CpGmarker, "character")
})

test_that("the CpG-identifier column is renamed to CpGmarker", {
    df <- data.frame(ProbeID = c("cg1", "cg2"), w = c(0.1, 0.2),
                     stringsAsFactors = FALSE)
    out <- normalize_coef(df)
    expect_true("CpGmarker" %in% names(out))
    expect_false("ProbeID" %in% names(out))
})

test_that("non-table resources pass through unchanged", {
    ref <- list(a = 1, b = 2)
    expect_identical(normalize_coef(ref), ref)
})

has_mirror <- function() {
    root <- tryCatch(methylclock:::.mcd_local_root(), error = function(e) NA)
    !is.na(root) && dir.exists(root)
}

test_that("resolved resources are in canonical form", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    mcd_cache_clear()
    epic <- mcd_resource("coefEPIC")          # was a tibble
    expect_identical(class(epic), "data.frame")
    expect_type(epic$CpGmarker, "character")
    bn <- mcd_resource("cpgs.bn")             # was a character vector
    expect_s3_class(bn, "data.frame")
    expect_true("CpGmarker" %in% names(bn))
})
