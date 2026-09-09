# Mitotic clocks. Most use the counter engine (a per-sample summary of a CpG
# set, no weights); RepliTali is an ordinary linear model (intercept + weights).

.register_mitotic <- function() {
    clock_register("epiTOC1", target = "mitotic", predictor = "counter",
        tissue = "multi-tissue (proliferative)",
        platform = c("450K", "EPIC"), resource = "coefEpiTOC1",
        intercept = FALSE, agg = "mean",
        license = "CC-BY", source = "EpiMitClocks (Yang 2016)",
        citation = "yang2016")
    clock_register("HypoClock", target = "mitotic", predictor = "counter",
        tissue = "multi-tissue (proliferative)",
        platform = c("450K", "EPIC"), resource = "coefHypoClock",
        intercept = FALSE, agg = "mean", transform = function(x) 1 - x,
        license = "CC-BY", source = "EpiMitClocks (Teschendorff 2020)",
        citation = "teschendorff2020")
    clock_register("stemTOC", target = "mitotic", predictor = "counter",
        tissue = "multi-tissue (proliferative)",
        platform = c("450K", "EPIC"), resource = "coefStemTOC",
        intercept = FALSE, agg = list(fun = "quantile", p = 0.95),
        license = "CC-BY", source = "EpiMitClocks (Zhu 2024)",
        citation = "zhu2024")
    clock_register("epiTOC2", target = "mitotic", predictor = "counter",
        tissue = "multi-tissue (proliferative)",
        platform = c("450K", "EPIC"), resource = "coefEpiTOC2",
        intercept = FALSE, agg = "mean", transform = function(x) 2 * x,
        license = "CC-BY", source = "EpiMitClocks (Teschendorff 2020)",
        citation = "teschendorff2020")
    # INTERNAL: epiTOC2 tnsc = 2 * mean_i[(b_i-beta0_i)/(delta_i*(1-beta0_i))];
    # INTERNAL: the per-CpG delta/beta0 map is applied by .counter_parametric.
    clock_register("RepliTali", target = "mitotic", generation = NA_integer_,
        tissue = "cultured cells",
        platform = c("450K", "EPIC"), resource = "coefRepliTali",
        license = "CC-BY", source = "EpiMitClocks (Endicott 2022)",
        citation = "endicott2022")
    invisible(NULL)
}
