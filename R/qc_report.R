# A one-call quality-control report over a set of clock estimates: sample
# descriptives, which clocks were computable, and how well each tracks age. Built
# on clockAccuracy(); returns the pieces and prints a short report.

#' Quality-control report for clock estimates
#'
#' Bundles the checks worth running before trusting clock estimates into one call:
#' sample descriptives, how many clocks were computable, how well each age clock
#' tracks chronological age (overall and within groups), and --- if you pass the
#' methylation matrix --- how many of each clock's CpGs your data carry. The
#' chronological ages are optional: without them the accuracy and acceleration
#' cannot be computed, and the report falls back to coverage, agreement and
#' distributions.
#'
#' @param x A \code{methylclock} result, or a data frame with an \code{id} column
#'   and one numeric column per clock.
#' @param age Optional numeric vector of chronological ages (one per sample in
#'   \code{x}'s order, or named by sample id). Needed only for accuracy
#'   (\eqn{R^2}, MAE) and age acceleration.
#' @param by Optional categorical vector (one value per sample, or named by id):
#'   the accuracy is broken down within each group.
#' @param clocks Optional clock names to restrict to. Default: all clocks in
#'   \code{x}.
#' @param betas Optional methylation matrix (or accepted container): if given,
#'   the report adds each clock's CpG coverage (see \code{\link{clockCoverage}})
#'   and how many of each clock's CpG values had to be imputed.
#' @param impute Optional label of the imputation method used to compute
#'   \code{x} (\code{"mean"}, \code{"knn"}, \code{"reference"} or \code{"none"}),
#'   recorded in the report so the reader knows how missing values were filled.
#' @param output \code{"object"} (the default) returns the report object;
#'   \code{"html"} also renders a small self-contained HTML report to \code{file}
#'   and returns its path invisibly.
#' @param file Where to write the HTML report when \code{output = "html"}.
#'   Default: a temporary file.
#' @return An object of class \code{methylclock_qc}: a list with \code{samples},
#'   \code{groups}, \code{computed}, \code{accuracy} (or \code{NULL} without
#'   \code{age}), \code{coverage} and \code{missing_by_clock} (or \code{NULL}
#'   without \code{betas}), \code{impute} (the recorded method, if any),
#'   \code{notes} (a written reading of each section, in plain language) and
#'   \code{commentary} (the one-line-per-section summary). With
#'   \code{output = "html"} the path to the rendered report is returned
#'   invisibly.
#' @seealso \code{\link{clockAccuracy}}, \code{\link{clockCoverage}}
#' @examples
#' data(methylclock_demo)
#' qcReport(methylclock_demo, age = methylclock_demo$age,
#'          by = methylclock_demo$sex)
#' @export
qcReport <- function(x, age = NULL, by = NULL, clocks = NULL, betas = NULL,
                     impute = NULL, output = c("object", "html"), file = NULL) {
    output <- match.arg(output)
    has_age <- !is.null(age)
    v0 <- .mc_plot_values(x, clocks)
    cl <- if (is.null(clocks))                       # drop stray numeric columns
        names(v0$values)[.mc_clock_family(names(v0$values)) != "other"]
        else names(v0$values)
    if (!length(cl)) cl <- names(v0$values)
    v <- .mc_plot_values(x, cl)
    computed <- vapply(v$values, function(col) sum(is.finite(col)), integer(1))
    groups <- if (!is.null(by))
        table(group = as.character(.mc_align_group(by, v$ids))) else NULL
    age_clk <- .mc_age_clocks(cl)
    acc <- if (has_age && length(age_clk))
        clockAccuracy(x, age, by = by, clocks = age_clk) else NULL
    coverage <- if (!is.null(betas)) clockCoverage(betas, clocks = cl) else NULL
    missing_by_clock <- if (!is.null(betas))
        tryCatch(.clock_missing(betas, clocks = cl), error = function(e) NULL)
    sample_qc <- tryCatch(sampleQC(x, betas = betas, clocks = cl),
                          error = function(e) NULL)
    # The figures show only the clocks that returned enough values to draw:
    # a clock that is all NA would otherwise take up a panel showing nothing.
    # The tables above keep every clock, so nothing is hidden.
    keep <- computed >= min(3L, length(v$ids))
    shown <- names(computed)[keep]
    if (!length(shown)) shown <- cl
    v_shown <- if (identical(shown, cl)) v else .mc_plot_values(x, shown)
    qc <- structure(list(
        samples = list(n = length(v$ids),
                       age = if (has_age) summary(.mc_align_age(age, v$ids))),
        groups = groups, computed = computed,
        shown = shown, omitted = setdiff(cl, shown),
        accuracy = acc, coverage = coverage, sample_qc = sample_qc,
        missing_by_clock = missing_by_clock,
        impute = if (!is.null(impute)) as.character(impute)[1],
        has_age = has_age),
        class = "methylclock_qc")
    qc$stats <- .qc_stats(v_shown, if (has_age) .mc_align_age(age, v$ids),
                          age_clk)
    qc$notes <- .qc_notes(qc)
    qc$commentary <- .qc_commentary(qc)

    if (output == "html") {
        plots <- .qc_plots(x, if (has_age) age, shown,
                           intersect(age_clk, shown))
        path <- .qc_render_html(qc, plots, file)
        message("methylclock QC report written to ", path)
        return(invisible(path))
    }
    qc
}

# Build the report's figures (as ggplots): agreement and distributions always,
# predicted-vs-real age only when ages are available.
.qc_plots <- function(x, age, cl, age_clk) {
    list(
        correlation = if (length(cl) >= 2L)
            tryCatch(plotClockCorrelation(x, clocks = cl),
                     error = function(e) NULL),
        distributions = tryCatch(plotClockDistributions(x, clocks = cl),
                                 error = function(e) NULL),
        age = if (!is.null(age) && length(age_clk))
            tryCatch(plotDNAmAge(x, age, clocks = age_clk),
                     error = function(e) NULL),
        bland_altman = if (!is.null(age) && length(age_clk))
            tryCatch(plotBlandAltman(x, age, clocks = age_clk),
                     error = function(e) NULL))
}

# Render the report object to a small self-contained HTML file.
.qc_render_html <- function(qc, plots = NULL, file = NULL) {
    if (!requireNamespace("rmarkdown", quietly = TRUE))
        stop("Rendering the HTML report needs the 'rmarkdown' package.",
             call. = FALSE)
    tmpl <- system.file("rmd", "qc_report.Rmd", package = "methylclock")
    if (!nzchar(tmpl))
        tmpl <- file.path("inst", "rmd", "qc_report.Rmd")   # dev (load_all)
    if (is.null(file)) file <- tempfile(fileext = ".html")
    rmarkdown::render(tmpl, output_file = normalizePath(file, mustWork = FALSE),
                      params = list(qc = qc, plots = plots), quiet = TRUE,
                      envir = new.env(parent = globalenv()))
    normalizePath(file)
}

#' @export
print.methylclock_qc <- function(x, ...) {
    cat("methylclock QC report\n")
    cat("- samples: ", x$samples$n, "\n", sep = "")
    if (!is.null(x$samples$age))
        cat("- age (years): ",
            paste(names(x$samples$age), round(x$samples$age, 1),
                  sep = " ", collapse = "  "), "\n", sep = "")
    if (!is.null(x$groups)) {
        cat("- groups: ",
            paste(names(x$groups), as.integer(x$groups), sep = "=",
                  collapse = "  "), "\n", sep = "")
    }
    miss <- names(x$computed)[x$computed < x$samples$n]
    cat("- clocks computed on all samples: ",
        sum(x$computed == x$samples$n), " of ", length(x$computed),
        if (length(miss))
            paste0(" (partial/NA: ", paste(miss, collapse = ", "), ")") else "",
        "\n", sep = "")
    if (!is.null(x$coverage))
        cat("- CpG coverage: ", min(x$coverage$pct_present), "%-",
            max(x$coverage$pct_present), "% across clocks\n", sep = "")
    if (!is.null(x$accuracy)) {
        cat("- accuracy vs chronological age:\n")
        print(x$accuracy, row.names = FALSE)
    } else cat("- accuracy: not computed (no chronological age supplied)\n")
    .qc_print_notes(x$notes)
    invisible(x)
}

# Section titles used when the report is printed at the console.
.QC_SECTIONS <- c(overview = "Overview", variables = "Sample variables",
                  coverage = "CpG coverage",
                  missing = "Missing values", agreement = "Agreement",
                  distributions = "Distributions",
                  accuracy = "Accuracy against age",
                  bland_altman = "Bias", samples = "Per-sample quality")

# Print the written reading of each section, wrapped to the console width.
.qc_print_notes <- function(notes) {
    if (!length(notes)) return(invisible(NULL))
    cat("\nReading of each section\n")
    for (nm in names(notes)) {
        cat("\n", .QC_SECTIONS[[nm]], "\n", sep = "")
        txt <- gsub("`", "", paste(notes[[nm]], collapse = " "))
        cat(strwrap(txt, width = 0.9 * getOption("width", 80), prefix = "  "),
            sep = "\n")
    }
    invisible(NULL)
}
