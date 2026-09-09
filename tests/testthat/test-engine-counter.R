# Counter engine: aggregates a CpG set per sample (mitotic-like clocks). No data
# mirror needed -- the resource is mocked and inputs are synthetic.

test_that("aggregation computes mean and quantile per sample", {
    B <- matrix(c(1, 2, 3, 10, 20, 30), nrow = 3,
                dimnames = list(paste0("cg", 1:3), c("A", "B")))
    expect_equal(methylclock:::.aggregate_counter(B, "mean"), c(A = 2, B = 20))
    q <- methylclock:::.aggregate_counter(B, list(fun = "quantile", p = 0.5))
    expect_equal(unname(q), c(2, 20))
    expect_error(methylclock:::.aggregate_counter(B, "median"), "aggregation")
})

test_that("engine reads a CpG set and aggregates it (mean, 1-mean, quantile)", {
    cpgs <- paste0("cg", 1:20)
    set.seed(7)
    B <- matrix(runif(20 * 4), 20, 4, dimnames = list(cpgs, paste0("S", 1:4)))
    src <- methylclock:::.as_source_from_matrix(B)
    testthat::local_mocked_bindings(
        mcd_resource = function(id, ...)
            data.frame(CpGmarker = cpgs, stringsAsFactors = FALSE),
        .package = "methylclock")

    mk <- function(agg, tr) structure(
        list(name = "T", predictor = "counter", resource = "x",
             agg = agg, transform = tr), class = "clock_entry")

    p_mean <- methylclock:::.predict_counter(mk("mean", identity), src)
    expect_equal(p_mean, colMeans(B), tolerance = 1e-12)

    p_hypo <- methylclock:::.predict_counter(mk("mean", function(x) 1 - x), src)
    expect_equal(p_hypo, 1 - colMeans(B), tolerance = 1e-12)

    p_q <- methylclock:::.predict_counter(
        mk(list(fun = "quantile", p = 0.95), identity), src)
    manual <- apply(B, 2L, stats::quantile, 0.95, names = FALSE)
    expect_equal(p_q, setNames(manual, colnames(B)), tolerance = 1e-12)
})

test_that("counter clock below coverage returns NA with a warning", {
    cpgs <- paste0("cg", 1:20)
    B <- matrix(0.5, nrow = 2, ncol = 3,
                dimnames = list(cpgs[1:2], paste0("S", 1:3)))   # only 2 of 20
    src <- methylclock:::.as_source_from_matrix(B)
    testthat::local_mocked_bindings(
        mcd_resource = function(id, ...)
            data.frame(CpGmarker = cpgs, stringsAsFactors = FALSE),
        .package = "methylclock")
    e <- structure(list(name = "T", predictor = "counter", resource = "x",
                        agg = "mean", transform = identity), class = "clock_entry")
    expect_warning(p <- methylclock:::.predict_counter(e, src), "returning NA")
    expect_true(all(is.na(p)))
})

test_that("clock_register validates the agg spec", {
    on.exit({ clock_reset(); register_builtin_clocks() }, add = TRUE)
    expect_error(
        clock_register("BadAgg", target = "mitotic", predictor = "counter",
                       resource = "x", agg = "median"), "must be")
    e <- clock_register("GoodQ", target = "mitotic", predictor = "counter",
                        resource = "x", agg = list(fun = "quantile", p = 0.95))
    expect_equal(e$agg$p, 0.95)
})
