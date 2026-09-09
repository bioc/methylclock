# Longitudinal shaping and rates, checked on synthetic subjects whose true
# slopes are known, so the machinery must recover them.

fake_long <- function(n_subj = 20, visits = 3, slope = 1.5, noise = 0.5) {
    set.seed(42)
    subject <- rep(sprintf("P%02d", seq_len(n_subj)), each = visits)
    start <- rep(runif(n_subj, 30, 60), each = visits)
    age <- start + rep(seq_len(visits) - 1, times = n_subj) * 2
    est <- 5 + slope * age + rnorm(length(age), 0, noise)
    list(df = data.frame(id = paste0("s", seq_along(age)),
                         Horvath = est, stringsAsFactors = FALSE),
         subject = subject, age = age)
}

test_that("clockTrajectories orders visits and excludes singletons", {
    f <- fake_long()
    tr <- clockTrajectories(f$df, f$subject, f$age)
    expect_setequal(names(tr), c("subject", "visit", "age", "clock",
                                 "estimate"))
    expect_equal(nrow(tr), 60L)
    expect_true(all(tr$visit %in% 1:3))
    # visits are numbered along age within each subject
    p1 <- tr[tr$subject == "P01", ]
    expect_equal(p1$visit[order(p1$age)], c(1, 2, 3))

    # a subject with a single age is excluded, with a message
    subj <- f$subject
    subj[subj == "P01"] <- c("P01", "S1", "S2")
    expect_message(tr2 <- clockTrajectories(f$df, subj, f$age), "single age")
    expect_false(any(c("P01", "S1", "S2") %in% tr2$subject))
})

test_that("replicates at the same age are averaged", {
    f <- fake_long(n_subj = 2, visits = 2)
    df <- rbind(f$df, f$df[1, ])            # duplicate the first sample
    df$id <- paste0("s", seq_len(nrow(df)))
    df$Horvath[5] <- df$Horvath[1] + 2      # a differing technical replicate
    expect_message(
        tr <- clockTrajectories(df, c(f$subject, f$subject[1]),
                                c(f$age, f$age[1])),
        "averaged")
    v <- tr$estimate[tr$subject == "P01" & tr$visit == 1]
    expect_equal(v, f$df$Horvath[1] + 1)    # the mean of the two runs
})

test_that("trajectoryRates recovers a known slope", {
    f <- fake_long(n_subj = 30, visits = 4, slope = 1.5, noise = 0.3)
    r <- trajectoryRates(f$df, f$subject, f$age)
    expect_equal(nrow(r), 30L)
    expect_true(all(r$n_visits == 4L))
    # the median recovered slope sits near the simulated truth
    expect_lt(abs(median(r$slope) - 1.5), 0.1)
    expect_true(all(is.finite(r$se)))       # >= 3 visits carry their SE
})

test_that("two-visit subjects are kept, marked without an SE", {
    f <- fake_long(n_subj = 10, visits = 2)
    r <- trajectoryRates(f$df, f$subject, f$age)
    expect_equal(nrow(r), 10L)
    expect_true(all(r$n_visits == 2L))
    expect_true(all(is.na(r$se)))
})

test_that("plotTrajectories draws spaghetti, plain and by group", {
    skip_if_not_installed("ggplot2")
    f <- fake_long(n_subj = 12, visits = 3)
    expect_s3_class(plotTrajectories(f$df, f$subject, f$age), "ggplot")
    grp <- rep(rep(c("a", "b"), each = 3), length.out = length(f$subject))
    expect_s3_class(plotTrajectories(f$df, f$subject, f$age, group = grp),
                    "ggplot")
})

test_that("the shipped longitudinal example flows through the tools", {
    data(methylclock_longitudinal, envir = environment())
    d <- methylclock_longitudinal
    # replicates at the same age are averaged, and every child has >= 2 ages
    tr <- suppressMessages(
        clockTrajectories(d, d$subject, d$age, clocks = "Wu"))
    expect_true(all(tapply(tr$age, tr$subject, length) >= 2))
    r <- suppressMessages(trajectoryRates(d, d$subject, d$age,
                                          clocks = "Wu"))
    # in fast-changing early childhood the pediatric clock advances with age
    expect_gt(median(r$slope), 0.5)
    skip_if_not_installed("ggplot2")
    expect_s3_class(suppressMessages(
        plotTrajectories(d, d$subject, d$age, clocks = c("Wu", "PedBE"))),
        "ggplot")
})

test_that("the rates feed the group forest directly", {
    skip_if_not_installed("ggplot2")
    f <- fake_long(n_subj = 24, visits = 3)
    r <- trajectoryRates(f$df, f$subject, f$age)
    wide <- data.frame(id = r$subject, Horvath = r$slope,
                       stringsAsFactors = FALSE)
    grp <- rep(c("a", "b"), length.out = nrow(wide))
    expect_s3_class(plotGroupDifference(wide, grp, what = "estimate",
                                        clocks = "Horvath"), "ggplot")
})
