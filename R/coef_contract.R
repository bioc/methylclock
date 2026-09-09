# Coefficient contract: the single place where a resource is coerced to the
# canonical layout the estimation engine expects. Downstream code sees only the
# canonical shape.
#
# Canonical form:
#   * a data.frame (never a tibble/matrix)
#   * a character column named `CpGmarker` with the CpG identifiers
#   * one or more numeric coefficient columns (a clock's registry entry names
#     the one it uses via `coef_col`); the intercept, when present, is kept as a
#     row so the predictor can strip it consistently
#
# INTERNAL: the legacy contract was implicit and false -- it assumed every coef*
# INTERNAL: object is a data.frame with CpGmarker / CoefficientTraining and row
# INTERNAL: 1 = intercept. 10 of 20 objects break that (DESIGN.md sec 7.1):
# INTERNAL: coefEPIC is a tibble, cpgs.bn is character, coefLeeGA has several
# INTERNAL: coefficient columns and none named CoefficientTraining, etc. Legacy
# INTERNAL: DNAmGA.R rebuilt a fake data.frame 3x and stripped the intercept by
# INTERNAL: hand in ~20 places inconsistently; all of that collapses here.

# Find the CpG-identifier column of a table, or NULL if there is none.
.find_cpg_col <- function(df) {
    nm <- names(df)
    hit <- nm[tolower(nm) %in% c("cpgmarker", "cpg", "probeid", "id", "name")]
    if (length(hit)) return(hit[1])
    chr <- nm[vapply(df, function(col) is.character(col) || is.factor(col),
                     logical(1))]
    if (length(chr)) return(chr[1])
    NULL
}

#' Normalize a coefficient resource to canonical form
#'
#' Coerces a coefficient object to the single canonical layout used by the
#' estimation engine (a data frame with a character \code{CpGmarker} column and
#' numeric coefficient columns), so that downstream code does not handle the
#' shape differences between clocks. Called by \code{\link{mcd_resource}} when a
#' resource is loaded. Objects that are not coefficient tables (e.g. cell-type
#' references) are returned unchanged.
#'
#' @param obj The raw resource object as loaded by a backend.
#' @param id The resource identifier (unused; kept for dispatch).
#' @param row The manifest row for \code{id} (unused; kept for dispatch).
#' @return The object in canonical form.
#' @seealso \code{\link{mcd_resource}}
#' @examples
#' # a table whose CpG column is not called CpGmarker
#' raw <- data.frame(probeID = c("(Intercept)", "cg00075967"),
#'                   weight = c(0.5, -0.02))
#' normalize_coef(raw)
#'
#' # a bare CpG vector becomes a one-column table
#' normalize_coef(c("cg00075967", "cg00374717"))
#' @export
normalize_coef <- function(obj, id = NULL, row = NULL) {
    # A bare CpG vector (a clock's probe list) becomes a one-column table.
    if (is.character(obj) || is.factor(obj))
        return(data.frame(CpGmarker = as.character(obj),
                          stringsAsFactors = FALSE))

    # Table-like input becomes a plain data.frame with a `CpGmarker` column.
    if (is.data.frame(obj) || is.matrix(obj)) {
        df <- as.data.frame(obj, stringsAsFactors = FALSE)
        cpg <- .find_cpg_col(df)
        if (!is.null(cpg)) {
            if (cpg != "CpGmarker")
                names(df)[names(df) == cpg] <- "CpGmarker"
            df$CpGmarker <- as.character(df$CpGmarker)
        }
        return(df)
    }

    # Other resources (lists, references) are not coefficient tables.
    obj
}
