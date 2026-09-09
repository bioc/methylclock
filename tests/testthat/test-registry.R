# Registry mechanics + built-in registrations.

test_that("built-in registry is populated with the 22 unified clocks", {
    reg <- clock_registry()
    expect_gte(length(reg), 22L)
    expect_true(all(c("Horvath", "Hannum", "Levine", "DunedinPACE",
                      "Knight", "Lee.RPC", "BNN") %in% names(reg)))
})

test_that("Levine is registered as the phenotypic (PhenoAge) clock", {
    e <- clock_info("Levine")
    expect_identical(e$target, "phenotypic")
    expect_identical(e$generation, 2L)
})

test_that("DunedinPACE is kept but flagged as license-risk (policy D9)", {
    e <- clock_info("DunedinPACE")
    expect_identical(e$status, "license-risk")
    expect_identical(e$aux, "coefDunedinPACEGS")
})

test_that("BNN uses an internal (compiled) predictor, no public resource", {
    e <- clock_info("BNN")
    expect_identical(e$predictor, "nn")
    expect_null(e$resource)
    expect_identical(e$license, "internal")
})

test_that("clock_register validates predictor and target", {
    on.exit(register_builtin_clocks(), add = TRUE)  # restore
    expect_error(clock_register("X", target = "chronological", predictor = "bogus"))
    expect_error(clock_register("X", target = "bogus"))
})

test_that("clock_list filters by target and predictor", {
    gest <- clock_list(target = "gestational")
    expect_true(all(gest$target == "gestational"))
    expect_true("Knight" %in% gest$name)
})

test_that("every clock declares the tissue it was trained on", {
    reg <- clock_registry()
    tis <- vapply(reg, function(e) e$tissue, character(1))
    expect_false(anyNA(tis))
    expect_identical(unname(tis["PedBE"]), "buccal epithelium (children)")
    expect_identical(unname(tis["Mayne"]), "placenta")
    expect_identical(unname(tis["Horvath"]), "multi-tissue")
})
