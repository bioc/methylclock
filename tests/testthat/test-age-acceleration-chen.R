test_that("ageAccelerationChen returns a per-sample chenAcc residual", {
    data(methylclock_cells)
    cells <- as.matrix(methylclock_cells[, c("CD8T", "CD4T", "NK", "Mono")])
    res <- ageAccelerationChen(methylclock_cells$Hannum,
                               age = methylclock_cells$age,
                               cells = cells,
                               immune = c("CD8T", "NK", "Mono"))
    expect_equal(nrow(res), nrow(methylclock_cells))
    expect_named(res, c("id", "age", "chenAcc"))
    # a (near) mean-zero residual that varies across samples
    expect_lt(abs(mean(res$chenAcc, na.rm = TRUE)), 1e-6)
    expect_gt(stats::sd(res$chenAcc, na.rm = TRUE), 0)
})

test_that("a missing immune column is an error, not a silent change", {
    data(methylclock_cells)
    cells <- as.matrix(methylclock_cells[, c("CD4T", "Mono"), drop = FALSE])
    expect_error(
        ageAccelerationChen(methylclock_cells$Hannum,
                            age = methylclock_cells$age, cells = cells),
        "missing the immune column")
    # even one absent column out of an explicit request must stop
    expect_error(
        ageAccelerationChen(methylclock_cells$Hannum,
                            age = methylclock_cells$age, cells = cells,
                            immune = c("CD4T", "NK")),
        "'NK'")
})

test_that("mismatched lengths are caught", {
    data(methylclock_cells)
    cells <- as.matrix(methylclock_cells[, c("CD8T", "NK")])
    expect_error(
        ageAccelerationChen(methylclock_cells$Hannum[-1],
                            age = methylclock_cells$age, cells = cells,
                            immune = c("CD8T", "NK")),
        "same length")
})
