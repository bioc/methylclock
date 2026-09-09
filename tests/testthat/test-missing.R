# Missing-CpG handling: coverage threshold (min.perc) and imputation.

has_mirror <- function() {
    root <- tryCatch(methylclock:::.mcd_local_root(), error = function(e) NA)
    !is.na(root) && dir.exists(root)
}

beta_matrix <- function() {
    p <- file.path(methylclock:::.mcd_local_root(), "..", "oracles",
                   "fixtures", "synthetic.rds")
    skip_if_not(file.exists(p), "fixture not found")
    df <- readRDS(p)
    m <- as.matrix(df[, -1])
    rownames(m) <- df$ProbeID
    m
}

clock_cpgs <- function(name) {
    e <- clock_info(name)
    as.character(mcd_resource(e$resource)$CpGmarker[-1])
}

test_that("a clock below the coverage threshold returns NA with a clear note", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    bm <- beta_matrix()
    cpgs <- clock_cpgs("Horvath")
    half <- bm[rownames(bm) %in% cpgs[seq_len(floor(length(cpgs) * 0.5))], ,
               drop = FALSE]
    res <- suppressWarnings(methylclock(half, clocks = "Horvath"))
    expect_true(all(is.na(as.data.frame(res)$Horvath)))
    expect_match(res$meta$warnings[[1]], "returning NA")
})

test_that("min.perc is configurable", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    bm <- beta_matrix()
    cpgs <- clock_cpgs("Horvath")
    most <- bm[rownames(bm) %in% cpgs[seq_len(floor(length(cpgs) * 0.6))], ,
               drop = FALSE]
    # 60% coverage: rejected at the 0.8 default, accepted at 0.5.
    r_default <- suppressWarnings(methylclock(most, clocks = "Horvath"))
    r_low <- suppressWarnings(methylclock(most, clocks = "Horvath",
                                          min.perc = 0.5))
    expect_true(all(is.na(as.data.frame(r_default)$Horvath)))
    expect_false(any(is.na(as.data.frame(r_low)$Horvath)))
})

test_that("NA values among present CpGs are imputed, not propagated", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    bm <- beta_matrix()
    cpgs <- clock_cpgs("Horvath")
    bm[cpgs[1:5], 1] <- NA
    res <- suppressWarnings(methylclock(bm, clocks = "Horvath"))
    expect_false(is.na(as.data.frame(res)$Horvath[1]))
})

test_that("clockCoverage reports how many of each clock's CpGs are present", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    bm <- beta_matrix()
    cov <- clockCoverage(bm, clocks = c("Horvath", "Hannum"))
    expect_setequal(names(cov), c("clock", "n_cpgs", "present", "pct_present"))
    expect_true(all(cov$present <= cov$n_cpgs))
    expect_true(all(cov$pct_present >= 0 & cov$pct_present <= 100))
    # dropping half of Horvath's CpGs halves its coverage
    cpgs <- clock_cpgs("Horvath")
    half <- bm[!(rownames(bm) %in% cpgs[seq_len(floor(length(cpgs) / 2))]), ]
    cov2 <- clockCoverage(half, clocks = "Horvath")
    expect_lt(cov2$pct_present, cov$pct_present[cov$clock == "Horvath"])
})

test_that("a present CpG that is all-NA is treated as absent, not NaN", {
    # A real cohort can carry a probe present in the array but blank for every
    # sample (QC). Its (undefined) mean must not poison the whole clock.
    skip_if_not(has_mirror(), "data-mirror not reachable")
    bm <- beta_matrix()
    cpgs <- clock_cpgs("Horvath")
    bm[cpgs[1], ] <- NA                              # one present CpG, all-NA
    res <- suppressWarnings(methylclock(bm, clocks = "Horvath"))
    expect_false(any(is.na(as.data.frame(res)$Horvath)))
})

test_that("a standardized clock survives an all-NA CpG (per-sample stats)", {
    # BLUP z-scores each sample; an all-NA CpG must not turn its mean/sd into NaN
    # and blank every sample.
    skip_if_not(has_mirror(), "data-mirror not reachable")
    bm <- beta_matrix()
    cpgs <- clock_cpgs("BLUP")
    present <- intersect(cpgs, rownames(bm))
    skip_if(length(present) < 2, "BLUP CpGs not in fixture")
    clean <- as.data.frame(methylclock(bm, clocks = "BLUP"))$BLUP
    bm[present[1], ] <- NA
    dirty <- suppressWarnings(
        as.data.frame(methylclock(bm, clocks = "BLUP"))$BLUP)
    expect_false(any(is.na(dirty)))
    expect_equal(dirty, clean, tolerance = 0.5)      # one dropped CpG barely moves it
})
