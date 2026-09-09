# T-10: out-of-core KNN imputation. Parity with the reference KNN imputer,
# equality between the in-memory and HDF5 (out-of-core) paths, and the streamed
# whole-array regime. The clock coefficients come from the mirror; skip when it
# is unreachable, and skip the parity/HDF5 pieces without their optional packages.

mirror_ok <- function() {
    root <- tryCatch(methylclock:::.mcd_local_root(), error = function(e) NA)
    !is.na(root) && dir.exists(root)
}

# A CpGs x samples beta matrix over a clock's markers, with a latent structure so
# nearby CpGs genuinely co-vary (KNN neighbours then carry signal), plus a
# scattering of missing values.
knn_betas <- function(clock = "Horvath", nsamp = 16, nmiss = 60, seed = 1) {
    co <- mcd_resource(methylclock:::clock_registry()[[clock]]$resource)
    cpgs <- setdiff(as.character(co$CpGmarker),
                    c("(Intercept)", "Intercept", "intercept"))
    set.seed(seed)
    lat <- matrix(runif(6 * nsamp), 6, nsamp)
    load <- matrix(rnorm(length(cpgs) * 6), length(cpgs), 6)
    m <- stats::plogis(load %*% lat +
                       matrix(rnorm(length(cpgs) * nsamp, 0, 0.3),
                              length(cpgs), nsamp))
    dimnames(m) <- list(cpgs, paste0("s", seq_len(nsamp)))
    m[cbind(sample(length(cpgs), nmiss), sample(nsamp, nmiss, TRUE))] <- NA
    m
}

# Open a matrix as a BigDataStatMeth HDF5 matrix.
knn_hdf5 <- function(m) {
    f <- tempfile(fileext = ".h5")
    h <- BigDataStatMeth::hdf5_create_matrix(
        filename = f, dataset = "in/betas", data = m, overwrite = TRUE)
    close(h)
    list(hm = BigDataStatMeth::hdf5_matrix(f, "in/betas"), file = f)
}

test_that("imputeKNN matches the reference KNN imputer", {
    skip_if_not_installed("impute")
    set.seed(3)
    m <- matrix(runif(300 * 14), 300, 14,
                dimnames = list(paste0("cg", 1:300), paste0("s", 1:14)))
    m[cbind(sample(300, 30), sample(14, 30, TRUE))] <- NA
    ours <- imputeKNN(m, k = 10)
    ref <- suppressMessages(impute::impute.knn(m, k = 10)$data)
    expect_equal(unname(ours), unname(ref), tolerance = 1e-10)
    expect_false(anyNA(ours))
})

test_that("imputeKNN leaves a complete matrix untouched", {
    m <- matrix(runif(50 * 6), 50, 6,
                dimnames = list(paste0("cg", 1:50), paste0("s", 1:6)))
    expect_equal(imputeKNN(m), m)
})

test_that("the C++ kernel matches impute.knn cell for cell (no split)", {
    skip_if_not_installed("impute")
    cfgs <- list(c(80, 12, 30, 10, 0.5), c(200, 20, 120, 10, 0.5),
                 c(150, 16, 60, 5, 0.5), c(300, 25, 200, 15, 0.5),
                 c(120, 14, 90, 10, 0.3))
    for (cf in cfgs) {
        set.seed(cf[1])
        m <- matrix(runif(cf[1] * cf[2], .05, .95), cf[1], cf[2])
        m[cbind(sample(cf[1], cf[3], TRUE),
                sample(cf[2], cf[3], TRUE))] <- NA
        ref <- suppressMessages(impute::impute.knn(
            m, k = cf[4], rowmax = cf[5], colmax = .9, maxp = 1e6)$data)
        ours <- suppressWarnings(methylclock:::mc_knn_impute_mem(
            m, k = cf[4], maxp = 5000L, rowmax = cf[5], colmax = .9))
        expect_equal(unname(ours), unname(ref), tolerance = 1e-10)
    }
})

test_that("the block split stays close to the exact result", {
    set.seed(5)
    lat <- matrix(runif(8 * 24), 8, 24)
    load <- matrix(rnorm(800 * 8), 800, 8)
    m <- stats::plogis(load %*% lat + matrix(rnorm(800 * 24, 0, .3), 800, 24))
    m[cbind(sample(800, 400, TRUE), sample(24, 400, TRUE))] <- NA
    exact <- suppressWarnings(
        methylclock:::mc_knn_impute_mem(m, k = 10, maxp = 5000L))
    split <- suppressWarnings(
        methylclock:::mc_knn_impute_mem(m, k = 10, maxp = 150L))
    na <- is.na(m)
    expect_false(anyNA(split))
    expect_gt(stats::cor(split[na], exact[na]), 0.95)   # clustering keeps them
})

test_that("the HDF5 C++ read path equals the in-memory result", {
    skip_if_not_installed("BigDataStatMeth")
    set.seed(7)
    m <- matrix(runif(300 * 20, .05, .95), 300, 20,
                dimnames = list(paste0("cg", 1:300), paste0("s", 1:20)))
    m[cbind(sample(300, 150, TRUE), sample(20, 150, TRUE))] <- NA
    f <- tempfile(fileext = ".h5")
    h <- BigDataStatMeth::hdf5_create_matrix(
        filename = f, dataset = "in/betas", data = m, overwrite = TRUE)
    close(h)
    on.exit(unlink(f), add = TRUE)
    mem <- suppressWarnings(
        methylclock:::mc_knn_impute_mem(m, k = 10, maxp = 5000L))
    disk <- suppressWarnings(methylclock:::mc_knn_impute_hdf5(
        f, "in", "betas", k = 10, maxp = 5000L))
    expect_equal(unname(mem), unname(disk), tolerance = 1e-12)
    expect_equal(m[!is.na(m)], unname(disk)[!is.na(m)])  # observed preserved
})

test_that("impute='knn' fills every gap and differs from mean imputation", {
    skip_if_not(mirror_ok(), "data-mirror not reachable")
    skip_if_not_installed("impute")
    m <- knn_betas()
    knn <- as.data.frame(suppressWarnings(
        methylclock(m, clocks = "Horvath", impute = "knn")))$Horvath
    mean_v <- as.data.frame(suppressWarnings(
        methylclock(m, clocks = "Horvath", impute = "mean")))$Horvath
    expect_false(anyNA(knn))
    expect_gt(max(abs(knn - mean_v)), 1e-6)          # genuinely different fill
})

test_that("with no missing values knn agrees with mean imputation", {
    skip_if_not(mirror_ok(), "data-mirror not reachable")
    skip_if_not_installed("impute")
    m <- knn_betas(nmiss = 0)
    knn <- as.data.frame(suppressWarnings(
        methylclock(m, clocks = "Horvath", impute = "knn")))$Horvath
    mean_v <- as.data.frame(suppressWarnings(
        methylclock(m, clocks = "Horvath", impute = "mean")))$Horvath
    expect_equal(knn, mean_v)
})

test_that("the HDF5 (out-of-core) knn path equals the in-memory result", {
    skip_if_not(mirror_ok(), "data-mirror not reachable")
    skip_if_not_installed("impute")
    skip_if_not_installed("BigDataStatMeth")
    m <- knn_betas()
    ram <- as.data.frame(suppressWarnings(
        methylclock(m, clocks = "Horvath", impute = "knn",
                    storage = "memory")))$Horvath
    ih <- knn_hdf5(m)
    on.exit({ close(ih$hm); unlink(ih$file) }, add = TRUE)
    disk <- as.data.frame(suppressWarnings(
        methylclock(ih$hm, clocks = "Horvath", impute = "knn",
                    storage = "memory")))$Horvath
    expect_equal(disk, ram, tolerance = 1e-10)
})

test_that("min.perc.sample still drops under-covered samples under knn", {
    skip_if_not(mirror_ok(), "data-mirror not reachable")
    skip_if_not_installed("impute")
    m <- knn_betas(nmiss = 0)
    half <- seq_len(floor(nrow(m) / 2))
    m[half, 4] <- NA                                 # sample 4: ~50% missing
    v <- as.data.frame(suppressWarnings(
        methylclock(m, clocks = "Horvath", impute = "knn",
                    min.perc.sample = 0.8)))$Horvath
    expect_true(is.na(v[4]))
    expect_equal(sum(is.na(v)), 1L)
})

# The streamed whole-array regime fires when a clock spans a large fraction of a
# disk-backed array. An array holding only the clock's CpGs has coverage 1, so it
# routes there: a single block (maxp above the row count) must reproduce the exact
# in-memory result, and a forced multi-block split must stay well correlated with
# it (block-boundary effects match the reference imputer's own divide-and-conquer).
test_that("the streamed whole-array regime reproduces the exact result", {
    skip_if_not(mirror_ok(), "data-mirror not reachable")
    skip_if_not_installed("impute")
    skip_if_not_installed("BigDataStatMeth")
    m <- knn_betas(nsamp = 20, nmiss = 60, seed = 5)
    mem <- as.data.frame(suppressWarnings(
        methylclock(m, clocks = "Horvath", impute = "knn",
                    storage = "memory")))$Horvath
    ih <- knn_hdf5(m)                                # array == clock rows only
    on.exit({ close(ih$hm); unlink(ih$file) }, add = TRUE)
    single <- as.data.frame(suppressWarnings(
        methylclock(ih$hm, clocks = "Horvath", impute = "knn", maxp = 5000,
                    storage = "memory")))$Horvath
    blocked <- as.data.frame(suppressWarnings(
        methylclock(ih$hm, clocks = "Horvath", impute = "knn", maxp = 100,
                    storage = "memory")))$Horvath
    expect_equal(single, mem, tolerance = 1e-10)     # single block == exact
    expect_false(anyNA(blocked))
    expect_gt(stats::cor(blocked, mem), 0.95)        # clustering keeps neighbours
})

# The whole-array sweep bounds RAM by reading the array in blocks whose row count
# is derived from a byte budget; shrinking the budget only changes how many blocks
# are read, not the answer. At a fixed maxp, one block and many blocks must stay
# closely correlated: the only difference is which CpGs share a block as neighbours.
test_that("the whole-array knn sweep is stable across block counts", {
    skip_if_not(mirror_ok(), "data-mirror not reachable")
    skip_if_not_installed("impute")
    skip_if_not_installed("BigDataStatMeth")
    m <- knn_betas(nsamp = 20, nmiss = 60, seed = 5)
    ih <- knn_hdf5(m)                                # array == clock rows only
    on.exit({ close(ih$hm); unlink(ih$file) }, add = TRUE)
    one <- as.data.frame(suppressWarnings(           # default budget: one block
        methylclock(ih$hm, clocks = "Horvath", impute = "knn", maxp = 50,
                    storage = "memory")))$Horvath
    old <- options(methylclock.knn_block_mb = 0.03)  # tiny budget: several blocks
    on.exit(options(old), add = TRUE)
    many <- as.data.frame(suppressWarnings(
        methylclock(ih$hm, clocks = "Horvath", impute = "knn", maxp = 50,
                    storage = "memory")))$Horvath
    expect_false(anyNA(many))
    expect_gt(stats::cor(many, one), 0.95)
})
