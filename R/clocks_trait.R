# Trait / exposure EpiScores (McCartney 2018): blood DNAm predictors of a trait
# rather than of age. They are opt-in -- a default clocks = "all" selection is
# age-focused and does not include them; request them by name or with
# target = "trait". Pure weighted sums (no intercept), linear engine.

.register_trait <- function() {
    mcc <- function(name, resource) {
        clock_register(name, target = "trait", generation = NA_integer_,
            tissue = "whole blood",
            platform = c("450K", "EPIC"), resource = resource,
            intercept = FALSE,
            license = "CC0/CC-BY (via biolearn, BSD-3)",
            source = "McCartney 2018 Genome Biol 19:136 (biolearn)",
            citation = "mccartney2018")
    }
    mcc("McCartney.Smoking",   "coefMcCartneySmoking")
    mcc("McCartney.Alcohol",   "coefMcCartneyAlcohol")
    mcc("McCartney.BMI",       "coefMcCartneyBMI")
    mcc("McCartney.BodyFat",   "coefMcCartneyBodyFat")
    mcc("McCartney.Education", "coefMcCartneyEducation")
    mcc("McCartney.HDL",       "coefMcCartneyHDL")
    mcc("McCartney.LDL",       "coefMcCartneyLDL")
    mcc("McCartney.TotalChol", "coefMcCartneyTotalChol")

    # mCigarette: elastic-net smoking (pack-years) EpiScore, no intercept.
    clock_register("mCigarette", target = "trait", generation = NA_integer_,
        tissue = "whole blood",
        platform = "EPIC", resource = "coefMcCigarette", intercept = FALSE,
        license = "MIT", source = "Chybowska 2025 Nat Commun 16:3210 (github)",
        citation = "chybowska2025")
    invisible(NULL)
}
