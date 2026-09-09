# INTERNAL: out-of-core input path (D1b/D17). An HDF5-backed methylation matrix
# INTERNAL: must give byte-identical estimates to the same matrix in memory,
# INTERNAL: reading only each clock's CpGs rather than the whole array.

mc_in_betas <- function() {
    root <- tryCatch(methylclock:::.mcd_local_root(), error = function(e) NA)
    skip_if_not(!is.na(root) && dir.exists(root), "data-mirror not reachable")
    p <- file.path(root, "..", "oracles", "fixtures", "synthetic.rds")
    skip_if_not(file.exists(p), "synthetic fixture not found")
    df <- readRDS(p)
    m <- as.matrix(df[, -1])
    rownames(m) <- df[[1]]
    m
}

# Write a matrix to a fresh HDF5 file and open it as a BigDataStatMeth matrix.
mc_in_hdf5 <- function(m) {
    f <- tempfile(fileext = ".h5")
    h <- BigDataStatMeth::hdf5_create_matrix(
        filename = f, dataset = "in/betas", data = m, overwrite = TRUE)
    close(h)
    list(hm = BigDataStatMeth::hdf5_matrix(f, "in/betas"), file = f)
}

test_that("an HDF5-backed input reproduces the in-memory result exactly", {
    skip_if_not_installed("BigDataStatMeth")
    m <- mc_in_betas()
    ram <- as.data.frame(suppressWarnings(methylclock(m, storage = "memory")))

    ih <- mc_in_hdf5(m)
    on.exit({ close(ih$hm); unlink(ih$file) }, add = TRUE)
    disk <- as.data.frame(
        suppressWarnings(methylclock(ih$hm, storage = "memory")))

    expect_equal(disk[order(disk$id), ], ram[order(ram$id), ],
                 tolerance = 1e-12)
})

test_that("a single linear clock matches from an HDF5 input", {
    skip_if_not_installed("BigDataStatMeth")
    m <- mc_in_betas()
    ram <- methylclock(m, clocks = "Horvath", storage = "memory")
    ih <- mc_in_hdf5(m)
    on.exit({ close(ih$hm); unlink(ih$file) }, add = TRUE)
    disk <- methylclock(ih$hm, clocks = "Horvath", storage = "memory")
    expect_equal(as.data.frame(disk), as.data.frame(ram), tolerance = 1e-12)
})

test_that("a quantile-normalised clock matches from an HDF5 input", {
    skip_if_not_installed("BigDataStatMeth")
    m <- mc_in_betas()
    ram <- methylclock(m, clocks = "DunedinPACE", storage = "memory")
    ih <- mc_in_hdf5(m)
    on.exit({ close(ih$hm); unlink(ih$file) }, add = TRUE)
    disk <- methylclock(ih$hm, clocks = "DunedinPACE", storage = "memory")
    expect_equal(as.data.frame(disk), as.data.frame(ram), tolerance = 1e-12)
})

test_that("a standardized clock is z-scored out-of-core from an HDF5 input", {
    skip_if_not_installed("BigDataStatMeth")
    m <- mc_in_betas()
    ram <- methylclock(m, clocks = c("EN", "BLUP"), storage = "memory")
    ih <- mc_in_hdf5(m)
    on.exit({ close(ih$hm); unlink(ih$file) }, add = TRUE)
    disk <- methylclock(ih$hm, clocks = c("EN", "BLUP"), storage = "memory")
    expect_equal(as.data.frame(disk), as.data.frame(ram), tolerance = 1e-9)
})

test_that("standardization never builds a full z-scored matrix", {
    # The scaled source applies the z-score lazily per subset; its stats match
    # base scale() computed the eager way.
    m <- mc_in_betas()
    st <- methylclock:::.dense_col_stats(m)
    z <- scale(m)
    expect_equal(st$center, attr(z, "scaled:center"),
                 tolerance = 1e-12, ignore_attr = TRUE)
    expect_equal(st$scale, attr(z, "scaled:scale"),
                 tolerance = 1e-12, ignore_attr = TRUE)
    # A subset scaled with the stored stats equals the eager z-score subset.
    src <- methylclock:::.as_scaled_source(
        methylclock:::.as_source_from_matrix(m), st$center, st$scale)
    cpgs <- rownames(m)[c(1, 100, 5000)]
    got <- methylclock:::.src_subset(src, cpgs)
    expect_equal(got, z[cpgs, ], tolerance = 1e-12)
})

test_that("mc_to_hdf5 stores a matrix on disk and reproduces the result", {
    skip_if_not_installed("BigDataStatMeth")
    clk <- c("Horvath", "Levine", "EN")
    m <- mc_in_betas()
    ram <- methylclock(m, clocks = clk, storage = "memory")
    f <- tempfile(fileext = ".h5")
    on.exit(unlink(f), add = TRUE)
    hm <- mc_to_hdf5(m, f)
    on.exit(close(hm), add = TRUE)
    disk <- methylclock(hm, clocks = clk, storage = "memory")
    expect_equal(as.data.frame(disk), as.data.frame(ram), tolerance = 1e-9)
})

test_that("mc_to_hdf5 imports a delimited text file", {
    skip_if_not_installed("BigDataStatMeth")
    m <- mc_in_betas()[seq_len(2000), , drop = FALSE]   # a manageable slice
    csv <- tempfile(fileext = ".csv")
    utils::write.csv(data.frame(cpg = rownames(m), m, check.names = FALSE),
                     csv, row.names = FALSE)
    f <- tempfile(fileext = ".h5")
    on.exit(unlink(c(csv, f)), add = TRUE)
    hm <- mc_to_hdf5(csv, f)
    on.exit(close(hm), add = TRUE)
    back <- as.matrix(hm, force = TRUE)
    expect_equal(dim(back), dim(m))
    expect_equal(back, m, tolerance = 1e-9, ignore_attr = TRUE)
})

test_that("a group of linear clocks batches into one on-disk cross-product", {
    # INTERNAL: D20. On a disk-backed input the batchable linear clocks are
    # INTERNAL: estimated together through BigDataStatMeth; result must equal the
    # INTERNAL: per-clock base path, and a below-coverage clock stays NA.
    skip_if_not_installed("BigDataStatMeth")
    clk <- c("Horvath", "Levine", "Hannum", "PedBE", "TL")   # raw linear group
    m <- mc_in_betas()
    ram <- as.data.frame(methylclock(m, clocks = clk, storage = "memory"))
    ih <- mc_in_hdf5(m)
    on.exit({ close(ih$hm); unlink(ih$file) }, add = TRUE)
    disk <- as.data.frame(suppressWarnings(
        methylclock(ih$hm, clocks = clk, storage = "memory")))
    expect_equal(disk[order(disk$id), ], ram[order(ram$id), ],
                 tolerance = 1e-12)

    # the batch engine is what produced these on the disk input
    plan <- methylclock:::.mc_batch_plan(
        methylclock:::.select_clocks(clk),
        methylclock:::.as_beta_source(ih$hm), NULL)
    expect_true(length(plan$groups) >= 1L)
    expect_length(plan$passthrough, 0L)
})

test_that("a whole-array clock is folded on disk and falls back on NAs", {
    # INTERNAL: D20. BLUP (standardized, ~10^5 CpGs) is estimated with a cross-file
    # INTERNAL: cross-product on the RAW input plus a per-sample z-score fold, with
    # INTERNAL: no scaled operand materialised. Requires BigDataStatMeth's cross-file
    # INTERNAL: crossprod (>= 2.0.4). With missing values the fold declines (NULL).
    skip_if_not_installed("BigDataStatMeth")
    m <- mc_in_betas()
    ram <- as.data.frame(methylclock(m, clocks = "BLUP", storage = "memory"))
    ih <- mc_in_hdf5(m)
    on.exit({ close(ih$hm); unlink(ih$file) }, add = TRUE)

    betas <- methylclock:::.as_beta_source(ih$hm)
    betas_z <- methylclock:::.zscore_source(betas)
    entry <- clock_info("BLUP")
    folded <- suppressWarnings(
        methylclock:::.predict_linear_folded_bdsm(entry, betas, betas_z))
    expect_false(is.null(folded))                      # the fold ran (not fallback)
    expect_equal(unname(folded), ram$BLUP, tolerance = 1e-9)

    # with a missing value the per-sample stats carry NaN -> fold declines (NULL)
    mna <- m; mna[1, 1] <- NA
    ihna <- mc_in_hdf5(mna)
    on.exit({ close(ihna$hm); unlink(ihna$file) }, add = TRUE)
    bna <- methylclock:::.as_beta_source(ihna$hm)
    bzna <- methylclock:::.zscore_source(bna)
    expect_null(methylclock:::.predict_linear_folded_bdsm(entry, bna, bzna))
})

test_that("results are streamed to disk clock by clock (D18)", {
    skip_if_not_installed("BigDataStatMeth")
    m <- mc_in_betas()
    f <- tempfile(fileext = ".h5")
    on.exit(unlink(f), add = TRUE)
    res <- suppressWarnings(
        methylclock(m, clocks = c("Horvath", "Levine", "EN"), file = f))
    expect_identical(res$storage, "hdf5")
    expect_true(file.exists(f))
    # The on-disk copy holds every computed clock.
    expect_setequal(res$handle$clocks, res$clocks)
})
