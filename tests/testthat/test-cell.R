# Cell-type deconvolution and intrinsic age acceleration.

has_mirror <- function() {
    root <- tryCatch(methylclock:::.mcd_local_root(), error = function(e) NA)
    !is.na(root) && dir.exists(root)
}

test_that("cellCounts recovers pure cell profiles as an identity", {
    skip_if_not(has_mirror(), "data-mirror not reachable")
    skip_if_not_installed("quadprog")
    ref <- methylclock:::.cell_reference("blood gse35069 complete")
    cc <- cellCounts(ref$beta, normalize = FALSE)
    # deconvolving each pure profile returns ~1 for its own type, ~0 elsewhere
    expect_equal(unname(diag(cc)), rep(1, ncol(cc)), tolerance = 1e-3)
    expect_equal(unname(rowSums(cc)), rep(1, nrow(cc)), tolerance = 1e-3)
})

test_that("listCellReferences names include the default blood panel", {
    expect_true("blood gse35069 complete" %in% listCellReferences())
})

test_that("FlowSorted panel not installed gives a clear, actionable error", {
    skip_if(requireNamespace("FlowSorted.BloodExtended.EPIC", quietly = TRUE),
            "FlowSorted.BloodExtended.EPIC is installed")
    m <- matrix(0.5, 3, 2, dimnames = list(c("cg1", "cg2", "cg3"), c("a", "b")))
    expect_error(cellCounts(m, reference = "FlowSorted.BloodExtended.EPIC"),
                 "does not ship")
})

test_that("the FlowSorted.Blood.EPIC (6-cell) panel deconvolves correctly", {
    skip_if_not_installed("FlowSorted.Blood.EPIC")
    ct <- getExportedValue("FlowSorted.Blood.EPIC", "IDOLOptimizedCpGs.compTable")
    cc <- cellCounts(ct, reference = "FlowSorted.Blood.EPIC", normalize = FALSE)
    # deconvolving the reference's own profiles returns the identity
    expect_equal(unname(diag(cc)), rep(1, ncol(cc)), tolerance = 1e-3)
    expect_equal(unname(rowSums(cc)), rep(1, nrow(cc)), tolerance = 1e-3)
})

test_that("ageAcceleration returns the expected measures", {
    set.seed(1)
    n <- 50
    age <- runif(n, 20, 80)
    df <- data.frame(id = paste0("S", seq_len(n)),
                     Horvath = age + rnorm(n, 0, 4),
                     Hannum = age + rnorm(n, 0, 6))
    aa <- ageAcceleration(df, age)
    expect_setequal(names(aa), c("id", "clock", "ageAcc", "residual"))
    expect_false("residualCells" %in% names(aa))
    expect_equal(aa$ageAcc[1], df$Horvath[1] - age[1])
    # the residual column is centred at zero by construction
    expect_lt(abs(mean(aa$residual[aa$clock == "Horvath"], na.rm = TRUE)), 1e-6)

    cc <- matrix(runif(n * 4), n, 4,
                 dimnames = list(df$id, c("A", "B", "C", "D")))
    aa2 <- ageAcceleration(df, age, cell_counts = cc)
    expect_true("residualCells" %in% names(aa2))
})
