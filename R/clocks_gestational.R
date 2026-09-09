# Gestational-age clocks. The three Lee predictors share one coefficient table,
# picking their weight column with `coef_col`.

.register_gestational <- function() {
    clock_register("Knight", target = "gestational", resource = "coefKnightGA",
        tissue = "cord blood",
        platform = "450K", citation = "knight2016")
    clock_register("Bohlin", target = "gestational", resource = "coefBohlin",
        tissue = "cord blood",
        platform = "450K", transform = function(x) x / 7, citation = "bohlin2016")
    # INTERNAL: Bohlin prediction is in days -> /7 to weeks (legacy DNAmGA.R:134).
    # INTERNAL: Bohlin -> D-02 (percentage fixed at 100% in legacy checks).
    clock_register("Mayne", target = "gestational", resource = "coefMayneGA",
        tissue = "placenta",
        platform = "450K", citation = "mayne2017")
    clock_register("EPIC", target = "gestational", resource = "coefEPIC",
        tissue = "cord blood",
        platform = "EPIC", citation = "haftorn2021")
    # INTERNAL: coefEPIC resource is a tibble (contract sec 7.1).
    clock_register("Lee.RPC", target = "gestational", resource = "coefLeeGA",
        tissue = "placenta",
        coef_col = "Coefficient_RPC", platform = c("450K", "EPIC"),
        citation = "lee2019")
    clock_register("Lee.CPC", target = "gestational", resource = "coefLeeGA",
        tissue = "placenta",
        coef_col = "Coefficient_CPC", platform = c("450K", "EPIC"),
        citation = "lee2019")
    clock_register("Lee.refRPC", target = "gestational", resource = "coefLeeGA",
        tissue = "placenta",
        coef_col = "Coefficient_refined_RPC", platform = c("450K", "EPIC"),
        citation = "lee2019")
    invisible(NULL)
}
