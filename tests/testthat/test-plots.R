# Clock plots. They consume a result (or a plain data frame of estimates) and
# return a ggplot; checked here on a small synthetic frame so they run without
# the data mirror.

fake_result <- function(n = 40) {
    set.seed(1)
    age <- runif(n, 20, 80)
    data.frame(
        id = paste0("S", seq_len(n)),
        Horvath = age + rnorm(n, 0, 4),
        Hannum = age + rnorm(n, 0, 6),
        Levine = age + rnorm(n, 0, 8),
        stringsAsFactors = FALSE)
}

test_that("the plot functions return ggplot objects", {
    skip_if_not_installed("ggplot2")
    df <- fake_result()
    age <- df$Horvath - rnorm(nrow(df), 0, 3)
    expect_s3_class(plotDNAmAge(df, age), "ggplot")
    expect_s3_class(plotClockCorrelation(df), "ggplot")
    expect_s3_class(plotAgeAcceleration(df, age, type = "box"), "ggplot")
    expect_s3_class(plotAgeAcceleration(df, age, type = "heatmap"), "ggplot")
    expect_s3_class(plotClockDistributions(df), "ggplot")
    expect_s3_class(plotClockDistributions(df, scale = TRUE), "ggplot")
})

test_that("the Block-2 plots (group acceleration, sample PCA) return ggplots", {
    skip_if_not_installed("ggplot2")
    df <- fake_result()
    age <- df$Horvath - rnorm(nrow(df), 0, 3)
    grp <- ifelse(age >= median(age), "old", "young")
    expect_s3_class(plotAccelerationByGroup(df, age, grp), "ggplot")
    expect_s3_class(plotSamplePCA(df), "ggplot")
    expect_s3_class(plotSamplePCA(df, color = age, color_label = "age"),
                    "ggplot")
    expect_s3_class(plotSamplePCA(df, color = grp), "ggplot")
})

test_that("plotClockDensities draws estimates and acceleration, by group", {
    skip_if_not_installed("ggplot2")
    df <- fake_result()
    age <- df$Horvath - rnorm(nrow(df), 0, 3)
    grp <- ifelse(age >= median(age), "old", "young")
    expect_s3_class(plotClockDensities(df), "ggplot")
    expect_s3_class(plotClockDensities(df, group = grp), "ggplot")
    expect_s3_class(plotClockDensities(df, group = grp,
        what = "acceleration", age = age), "ggplot")
    # acceleration is a residual, so it needs the chronological age
    expect_error(plotClockDensities(df, what = "acceleration"), "age")
    # a named group is aligned by sample id, like the other plots
    named <- setNames(grp, df$id)
    expect_s3_class(plotClockDensities(df, group = rev(named)), "ggplot")
})

test_that("plotDNAmAge draws per-group points and fits", {
    skip_if_not_installed("ggplot2")
    df <- fake_result()
    age <- df$Horvath - rnorm(nrow(df), 0, 3)
    grp <- ifelse(age >= median(age), "old", "young")
    expect_s3_class(plotDNAmAge(df, age, group = grp), "ggplot")
    grp[1:4] <- NA                      # missing groups are dropped
    expect_s3_class(plotDNAmAge(df, age, group = grp), "ggplot")
})

test_that("plotGroupDifference builds the two-group forest", {
    skip_if_not_installed("ggplot2")
    df <- fake_result()
    age <- df$Horvath - rnorm(nrow(df), 0, 3)
    grp <- ifelse(age >= median(age), "old", "young")
    expect_s3_class(plotGroupDifference(df, grp, age = age), "ggplot")
    expect_s3_class(plotGroupDifference(df, grp, age = age,
                                        style = "dumbbell"), "ggplot")
    # without ages the raw estimates can still be compared, stated as such
    expect_s3_class(plotGroupDifference(df, grp, what = "estimate"), "ggplot")
    expect_error(plotGroupDifference(df, grp), "age")
    grp3 <- sample(c("a", "b", "c"), nrow(df), replace = TRUE)
    expect_error(plotGroupDifference(df, grp3, age = age), "two levels")
    # a clock with too few values in a group is skipped, with a message
    df$Levine[grp == "old"] <- NA
    expect_message(plotGroupDifference(df, grp, age = age), "Levine")
})

test_that("plotClockDensities offers mirror and ridge layouts", {
    skip_if_not_installed("ggplot2")
    df <- fake_result()
    age <- df$Horvath - rnorm(nrow(df), 0, 3)
    grp <- ifelse(age >= median(age), "old", "young")
    expect_s3_class(plotClockDensities(df, group = grp, style = "mirror"),
                    "ggplot")
    expect_s3_class(plotClockDensities(df, group = grp, style = "mirror",
        what = "acceleration", age = age), "ggplot")
    expect_s3_class(plotClockDensities(df, style = "ridge"), "ggplot")
    expect_s3_class(plotClockDensities(df, group = grp, style = "ridge"),
                    "ggplot")
    # a mirror has an up side and a down side: exactly two groups
    expect_error(plotClockDensities(df, style = "mirror"), "two levels")
    grp3 <- sample(c("a", "b", "c"), nrow(df), replace = TRUE)
    expect_error(plotClockDensities(df, group = grp3, style = "mirror"),
                 "two levels")
})

test_that("plotClockDensities falls back to ticks for tiny groups", {
    skip_if_not_installed("ggplot2")
    df <- fake_result()
    grp <- rep(c("big", "tiny"), c(nrow(df) - 3L, 3L))
    expect_message(p <- plotClockDensities(df, group = grp), "ticks")
    expect_s3_class(p, "ggplot")
    # samples with a missing group are dropped, not plotted as a group
    grp[1:5] <- NA
    expect_s3_class(suppressMessages(plotClockDensities(df, group = grp)),
                    "ggplot")
})

test_that("clockAccuracy summarises fit, overall and by group", {
    df <- fake_result()
    age <- df$Horvath - rnorm(nrow(df), 0, 3)
    a <- clockAccuracy(df, age)
    expect_true(all(c("clock", "n", "r", "R2", "MAE") %in% names(a)))
    expect_true(all(a$R2 >= 0 & a$R2 <= 1, na.rm = TRUE))
    # Horvath tracks age closest here, so its error is the smallest.
    expect_lt(a$MAE[a$clock == "Horvath"], a$MAE[a$clock == "Levine"])

    grp <- ifelse(age >= median(age), "old", "young")
    ag <- clockAccuracy(df, age, by = grp)
    expect_true("group" %in% names(ag))
    expect_setequal(unique(ag$group), c("old", "young"))
    expect_equal(nrow(ag), 2L * nrow(a))
})

test_that("qcReport bundles descriptives, coverage and accuracy", {
    df <- fake_result()
    age <- df$Horvath - rnorm(nrow(df), 0, 3)
    grp <- ifelse(age >= median(age), "old", "young")
    q <- qcReport(df, age, by = grp)
    expect_s3_class(q, "methylclock_qc")
    expect_equal(q$samples$n, nrow(df))
    expect_setequal(names(q$groups), c("old", "young"))
    expect_true(all(q$computed == nrow(df)))     # all clocks computed here
    expect_true("group" %in% names(q$accuracy))
    expect_output(print(q), "QC report")
    # observational commentary, and no judgement words
    expect_true(length(q$commentary) >= 3)
    expect_false(any(grepl("good|bad|excellent|poor|accurate",
                           q$commentary, ignore.case = TRUE)))
})

test_that("qcReport writes a reading of each section", {
    df <- fake_result()
    age <- df$Horvath - rnorm(nrow(df), 0, 3)
    q <- qcReport(df, age)
    # one entry per section that the report actually has
    expect_true(all(c("overview", "agreement", "distributions", "accuracy",
                      "bland_altman", "samples") %in% names(q$notes)))
    expect_null(q$notes$coverage)                # no betas were supplied
    txt <- unlist(q$notes)
    expect_true(all(nzchar(txt)))
    # the notes describe; they do not rate a clock or tell the reader what to do
    expect_false(any(grepl(
        "good|bad|excellent|poor|accurate|reliable|should|recommend|caution",
        txt, ignore.case = TRUE)))
    # and they quote this dataset's own numbers, not a generic template
    expect_match(paste(q$notes$overview, collapse = " "),
                 paste(nrow(df), "samples"), fixed = TRUE)
    expect_output(print(q), "Reading of each section")
    expect_output(print(q), "Agreement")
})

test_that("the figures drop clocks that returned nothing", {
    df <- fake_result(30)
    df$Levine <- NA_real_                        # a clock with no values at all
    q <- qcReport(df)
    expect_false("Levine" %in% q$shown)
    expect_equal(q$omitted, "Levine")
    expect_match(paste(q$notes$distributions, collapse = " "), "Levine")
    expect_match(paste(q$notes$agreement, collapse = " "), "left out")
    # the tables still account for every clock
    expect_true("Levine" %in% names(q$computed))
})

test_that("the report says which sample variables it was given", {
    df <- fake_result()
    expect_match(paste(qcReport(df)$notes$variables, collapse = " "),
                 "No sample variables")
    age <- df$Horvath - rnorm(nrow(df), 0, 3)
    expect_match(paste(qcReport(df, age)$notes$variables, collapse = " "),
                 "No grouping variable")
    grp <- ifelse(age >= median(age), "old", "young")
    expect_match(paste(qcReport(df, age, by = grp)$notes$variables,
                       collapse = " "), "age and a grouping variable")
})

test_that("the section reading follows the numbers it is given", {
    df <- fake_result(40)
    # two clocks made identical: they must come out as the closest pair
    df$Hannum <- df$Horvath
    q <- qcReport(df)
    expect_equal(sort(q$stats$agreement$closest$clocks),
                 c("Hannum", "Horvath"))
    expect_match(paste(q$notes$agreement, collapse = " "), "Hannum and Horvath")
})

test_that("plotBlandAltman returns a ggplot", {
    skip_if_not_installed("ggplot2")
    df <- fake_result()
    age <- df$Horvath - rnorm(nrow(df), 0, 3)
    expect_s3_class(plotBlandAltman(df, age), "ggplot")
})

test_that("sampleQC flags per-sample outliers", {
    df <- fake_result(60)
    # a sample extreme across every clock (a technical outlier), not just one
    df[1, c("Horvath", "Hannum", "Levine")] <-
        df[1, c("Horvath", "Hannum", "Levine")] + 90
    sq <- sampleQC(df)
    expect_true(all(c("id", "outlier_score", "is_outlier", "discordance")
                    %in% names(sq)))
    expect_equal(nrow(sq), nrow(df))
    expect_true(sq$is_outlier[1])
})

test_that("discordance separates disagreement from joint extremeness", {
    df <- fake_result(60)
    df$Levine[2] <- df$Levine[2] + 90     # one clock wildly off: disagreement
    sq <- sampleQC(df)
    expect_equal(which.max(sq$discordance), 2L)
    expect_s3_class(plotSampleDiscordance(df), "ggplot")
    expect_s3_class(plotSampleDiscordance(df, label_n = 0), "ggplot")
    age <- df$Horvath - rnorm(nrow(df), 0, 3)
    expect_s3_class(plotSampleDiscordance(df, what = "acceleration",
                                          age = age), "ggplot")
    expect_error(plotSampleDiscordance(df, what = "acceleration"), "age")
    # a disagreement needs at least three opinions
    expect_error(plotSampleDiscordance(df[, c("id", "Horvath", "Hannum")]),
                 "three clocks")
})

test_that("plotReferenceRange draws bands from a declared reference", {
    skip_if_not_installed("ggplot2")
    data(methylclock_demo, envir = environment())
    pred <- methylclock_demo[seq_len(12), ]
    # no reference: a shipped cohort is used, and the choice is announced
    expect_message(p <- plotReferenceRange(pred, age = pred$age), "GSE40279")
    expect_s3_class(p, "ggplot")
    # a named shipped cohort, and a user's own data frame, are quiet
    expect_silent(plotReferenceRange(pred, age = pred$age,
                                     reference = "GSE132203"))
    # the population-specific splits of GSE40279 are shipped references too
    expect_silent(plotReferenceRange(pred, age = pred$age,
                                     reference = "GSE40279-Hispanic"))
    expect_silent(plotReferenceRange(pred, age = pred$age,
                                     reference = "GSE40279-Caucasian"))
    # so are the reference cohorts shipped in methylclock_references
    expect_silent(plotReferenceRange(pred, age = pred$age,
                                     reference = "GSE87571"))
    expect_silent(plotReferenceRange(pred, age = pred$age,
                                     reference = "GSE210255"))
    expect_silent(plotReferenceRange(pred, age = pred$age,
                                     reference = "GSE224363"))
    # the Chinese controls live in a ~6-year window: bands need samples
    # of a matching age, and the reference declares itself in the subtitle
    p6 <- plotReferenceRange(pred, age = pred$age, reference = "GSE116379",
                             min_n = 10)
    expect_s3_class(p6, "ggplot")
    own <- methylclock_demo[100:300, c("age", "Horvath", "Hannum")]
    expect_silent(p2 <- plotReferenceRange(pred, age = pred$age,
                                           reference = own))
    expect_s3_class(p2, "ggplot")
    expect_error(plotReferenceRange(pred, age = pred$age,
                                    reference = "nope"), "shipped")
    noage <- own[, c("Horvath", "Hannum")]
    expect_error(plotReferenceRange(pred, age = pred$age, reference = noage),
                 "age")
    expect_error(plotReferenceRange(pred, age = rep(NA_real_, nrow(pred))),
                 "finite")
})

test_that("qcReport works without ages (no accuracy)", {
    df <- fake_result()
    q <- qcReport(df)
    expect_null(q$accuracy)
    expect_null(q$samples$age)
    expect_false(isTRUE(q$has_age))
    expect_match(paste(q$commentary, collapse = " "),
                 "[Nn]o chronological age")
    expect_output(print(q), "not computed")
})

test_that("qcReport renders an HTML report", {
    skip_if_not_installed("rmarkdown")
    skip_if(!nzchar(Sys.which("pandoc")), "pandoc not available")
    df <- fake_result()
    age <- df$Horvath - rnorm(nrow(df), 0, 3)
    f <- tempfile(fileext = ".html")
    p <- qcReport(df, age, output = "html", file = f)
    expect_true(file.exists(p))
    expect_gt(file.info(p)$size, 1000)
})

test_that("a PCA needs at least two informative clocks", {
    df <- data.frame(id = c("a", "b", "c"), Horvath = c(50, 60, 70),
                     Flat = c(1, 1, 1), stringsAsFactors = FALSE)
    expect_error(plotSamplePCA(df), "two informative")
})

test_that("a categorical covariate is aligned without losing its type", {
    ids <- c("a", "b", "c")
    expect_identical(
        methylclock:::.mc_align_group(c(c = "M", a = "F", b = "M"), ids),
        c("F", "M", "M"))
    expect_error(methylclock:::.mc_align_group(c("F", "M"), ids), "one value")
})

test_that("age is aligned by id when named, else by position", {
    ids <- c("a", "b", "c")
    expect_equal(methylclock:::.mc_align_age(c(c = 3, a = 1, b = 2), ids),
                 c(1, 2, 3))
    expect_equal(methylclock:::.mc_align_age(c(10, 20, 30), ids),
                 c(10, 20, 30))
    expect_error(methylclock:::.mc_align_age(c(1, 2), ids), "one value")
})

test_that("only numeric clock columns are used, id is dropped", {
    df <- data.frame(id = c("a", "b"), Horvath = c(50, 60),
                     note = c("x", "y"), stringsAsFactors = FALSE)
    v <- methylclock:::.mc_plot_values(df)
    expect_identical(names(v$values), "Horvath")
    expect_identical(v$ids, c("a", "b"))
})

test_that("a correlation needs at least two clocks", {
    df <- data.frame(id = c("a", "b"), Horvath = c(1, 2))
    expect_error(plotClockCorrelation(df), "two clocks")
})
