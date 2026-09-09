# NEO chronological-age predictors (platform-specific, standardized input).

.register_neo <- function() {
    clock_register("NEOaPMA450K", target = "chronological",
        tissue = "buccal (preterm infants)",
        resource = "coefNEOaPMA450K", platform = "450K",
        standardize = TRUE, citation = "graw2021neoage")
    clock_register("NEOaPNA450K", target = "chronological",
        tissue = "buccal (preterm infants)",
        resource = "coefNEOaPNA450K", platform = "450K",
        standardize = TRUE, citation = "graw2021neoage")
    clock_register("NEOaPMAEPIC", target = "chronological",
        tissue = "buccal (preterm infants)",
        resource = "coefNEOaPMAEPIC", platform = "EPIC",
        standardize = TRUE, citation = "graw2021neoage")
    clock_register("NEOaPNAEPIC", target = "chronological",
        tissue = "buccal (preterm infants)",
        resource = "coefNEOaPNAEPIC", platform = "EPIC",
        standardize = TRUE, citation = "graw2021neoage")
    invisible(NULL)
}
