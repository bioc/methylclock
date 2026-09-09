# Shared look for the package plots: one clean theme and a colour-blind-safe
# palette, so the default output is modern and consistent and the user can reuse
# the same theme on their own figures.

# Two distinct colour-blind-safe palettes for two distinct jobs, so a reader
# never confuses a clock family with a sample group: Okabe-Ito is kept for clock
# *type* (its registry target); Paul Tol's "bright" palette is used for sample
# *variables* (sex, cohort, smoking status...) — its first colours contrast
# strongly, so two overlaid groups stay apart at a glance.
.mc_qual <- c("#0072B2", "#E69F00", "#009E73", "#D55E00", "#56B4E9",
              "#CC79A7", "#F0E442", "#999999")
.mc_group_pal <- c("#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE",
                   "#AA3377", "#BBBBBB")

# Single accent used where one series is drawn (points, a lone violin/box).
.mc_accent <- "#0072B2"
.mc_accent_soft <- "#c6dbef"   # light blue, for a violin or CI band
.mc_accent_mid <- "#6baed6"    # medium blue, for a filled box
.mc_outlier <- "#D55E00"       # vermilion, to flag points off the fit

# A discrete fill/colour scale from the palette (recycled if there are more
# levels than colours).
.mc_fill_qual <- function(name = NULL, n = length(.mc_qual)) {
    ggplot2::scale_fill_manual(values = rep_len(.mc_qual, max(n, 1)),
                               name = name)
}
.mc_colour_qual <- function(name = NULL, n = length(.mc_qual)) {
    ggplot2::scale_colour_manual(values = rep_len(.mc_qual, max(n, 1)),
                                 name = name)
}

# Scales for a sample *variable* (a separate palette from the clock-type one).
.mc_fill_group <- function(name = NULL, n = length(.mc_group_pal)) {
    ggplot2::scale_fill_manual(values = rep_len(.mc_group_pal, max(n, 1)),
                               name = name)
}
.mc_colour_group <- function(name = NULL, n = length(.mc_group_pal)) {
    ggplot2::scale_colour_manual(values = rep_len(.mc_group_pal, max(n, 1)),
                                 name = name)
}

# Turn a per-sample grouping into a factor, keeping the caller's level order
# when one was given: the order decides which palette colour each group gets,
# so a user's factor keeps its colours stable across plots.
.mc_group_factor <- function(values, original) {
    if (is.factor(original))
        factor(as.character(values), levels = levels(original))
    else factor(values)
}

#' A clean theme for methylclock plots
#'
#' The theme the package's plotting functions use by default: a light background,
#' restrained gridlines and clear facet labels. It returns a standard
#' \code{ggplot2} theme, so you can add it to your own plots or override any part
#' of it --- every plotting function returns a \code{ggplot} object, so the
#' defaults are a starting point, not a constraint.
#'
#' @param base_size Base font size in points. Default 12.
#' @param base_family Base font family. Default "" (the device default).
#' @return A \code{ggplot2} theme object.
#' @examples
#' library(ggplot2)
#' ggplot(mtcars, aes(wt, mpg)) + geom_point() + theme_methylclock()
#' @export
theme_methylclock <- function(base_size = 12, base_family = "") {
    ggplot2::theme_minimal(base_size = base_size, base_family = base_family) +
        ggplot2::theme(
            panel.grid.minor = ggplot2::element_blank(),
            panel.grid.major = ggplot2::element_line(colour = "grey92",
                                                     linewidth = 0.3),
            panel.background = ggplot2::element_rect(fill = "white",
                                                     colour = NA),
            plot.background = ggplot2::element_rect(fill = "white",
                                                    colour = NA),
            axis.title = ggplot2::element_text(colour = "grey25"),
            axis.text = ggplot2::element_text(colour = "grey35"),
            strip.background = ggplot2::element_rect(fill = "grey95",
                                                     colour = NA),
            strip.text = ggplot2::element_text(face = "bold", colour = "grey20",
                margin = ggplot2::margin(4, 4, 4, 4)),
            legend.key = ggplot2::element_blank(),
            plot.title = ggplot2::element_text(face = "bold"),
            plot.margin = ggplot2::margin(8, 10, 8, 8))
}

# The registry target ("chronological", "mitotic", ...) of each clock, for
# colouring boxes by clock family when several families are shown together.
# Falls back to "other" for names not in the registry.
.mc_clock_family <- function(clocks) {
    fam <- vapply(clocks, function(nm) {
        info <- tryCatch(clock_info(nm), error = function(e) NULL)
        if (is.null(info) || is.null(info$target)) "other" else info$target
    }, character(1))
    fam
}
