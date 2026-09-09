# INTERNAL: D12-D16 -- HDF5 is the source of truth for every result; memory is a
# INTERNAL: read cache or a fallback. Persisting is a one-call promotion.

# A small methylation matrix that covers a few linear clocks. Skips (rather than
# errors) when the data mirror is unreachable, e.g. under R CMD check.
mc_h5_betas <- function() {
    root <- tryCatch(methylclock:::.mcd_local_root(), error = function(e) NA)
    skip_if_not(!is.na(root) && dir.exists(root), "data-mirror not reachable")
    p <- file.path(root, "..", "oracles", "fixtures", "synthetic.rds")
    skip_if_not(file.exists(p), "synthetic fixture not found")
    readRDS(p)
}

mc_h5_clocks <- c("Horvath", "Hannum", "Levine", "Knight", "DunedinPACE")

test_that("auto backs a result on disk and keeps an in-memory cache (D12)", {
    skip_if_not_installed("BigDataStatMeth")
    betas <- mc_h5_betas()
    res <- suppressWarnings(methylclock(betas, clocks = mc_h5_clocks))
    on.exit(unlink(res$handle$file), add = TRUE)

    expect_identical(res$storage, "hdf5")
    expect_false(res$persistent)                   # session file
    expect_false(is.null(res$values))              # cache present
    expect_true(file.exists(res$handle$file))
    # Cache and the on-disk copy must agree.
    disk <- res
    disk$values <- NULL
    expect_equal(as.data.frame(res), as.data.frame(disk), tolerance = 1e-12)
})

test_that("hdf5 estimates match a pure in-memory computation", {
    skip_if_not_installed("BigDataStatMeth")
    betas <- mc_h5_betas()
    mem <- methylclock(betas, clocks = mc_h5_clocks, storage = "memory")
    f <- tempfile(fileext = ".h5")
    on.exit(unlink(f), add = TRUE)
    h5 <- methylclock(betas, clocks = mc_h5_clocks, storage = "hdf5", file = f)

    expect_identical(mem$storage, "memory")
    expect_identical(h5$storage, "hdf5")
    expect_true(h5$persistent)                     # user-supplied path
    expect_null(h5$values)                         # hdf5 -> no cache
    expect_equal(as.data.frame(h5), as.data.frame(mem), tolerance = 1e-12)
})

test_that("persist promotes a session result without recomputing (D16)", {
    skip_if_not_installed("BigDataStatMeth")
    betas <- mc_h5_betas()
    res <- suppressWarnings(methylclock(betas, clocks = mc_h5_clocks))
    on.exit(unlink(res$handle$file), add = TRUE)
    dst <- tempfile(fileext = ".h5")
    on.exit(unlink(dst), add = TRUE)

    before <- as.data.frame(res)
    kept <- persist(res, dst)
    expect_true(kept$persistent)
    expect_identical(kept$handle$file, dst)
    expect_true(file.exists(dst))
    expect_equal(as.data.frame(kept), before, tolerance = 1e-12)
})

test_that("persist refuses to clobber an existing file unless told to", {
    skip_if_not_installed("BigDataStatMeth")
    betas <- mc_h5_betas()
    res <- suppressWarnings(methylclock(betas, clocks = mc_h5_clocks))
    on.exit(unlink(res$handle$file), add = TRUE)
    dst <- tempfile(fileext = ".h5")
    file.create(dst)
    on.exit(unlink(dst), add = TRUE)
    expect_error(persist(res, dst), "already exists")
    expect_message(kept <- persist(res, dst, overwrite = TRUE), "persisted")
})

test_that("memory storage keeps a result off disk (opt-out)", {
    skip_if_not_installed("BigDataStatMeth")
    betas <- mc_h5_betas()
    res <- methylclock(betas, clocks = mc_h5_clocks, storage = "memory")
    expect_identical(res$storage, "memory")
    expect_null(res$handle)
    expect_false(is.null(res$values))
})

test_that("a disk write failure falls back to memory with a warning (D14)", {
    skip_if_not_installed("BigDataStatMeth")
    betas <- mc_h5_betas()
    bad <- file.path("/no", "such", "dir", "x.h5")
    expect_warning(
        res <- methylclock(betas, clocks = mc_h5_clocks, file = bad),
        "memory only")
    expect_identical(res$storage, "memory")
    expect_null(res$handle)
})

test_that("summary() works over an hdf5-backed result", {
    skip_if_not_installed("BigDataStatMeth")
    betas <- mc_h5_betas()
    f <- tempfile(fileext = ".h5")
    on.exit(unlink(f), add = TRUE)
    h5 <- methylclock(betas, clocks = mc_h5_clocks, storage = "hdf5", file = f)
    s <- summary(h5)
    expect_true(is.table(s) || is.matrix(s))
})
