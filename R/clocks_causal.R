# Causal-aging clocks (Ying 2024): damage, adaptation and net causal age.

.register_causal <- function() {
    clock_register("CausAge", target = "causal", generation = NA_integer_,
        tissue = "whole blood",
        platform = c("450K", "EPIC"), resource = "coefCausAge",
        license = "BSD-3 (biolearn)", source = "biolearn (Ying 2024)",
        citation = "ying2024")
    clock_register("DamAge", target = "causal", generation = NA_integer_,
        tissue = "whole blood",
        platform = c("450K", "EPIC"), resource = "coefDamAge",
        license = "BSD-3 (biolearn)", source = "biolearn (Ying 2024)",
        citation = "ying2024")
    clock_register("AdaptAge", target = "causal", generation = NA_integer_,
        tissue = "whole blood",
        platform = c("450K", "EPIC"), resource = "coefAdaptAge",
        license = "BSD-3 (biolearn)", source = "biolearn (Ying 2024)",
        citation = "ying2024")
    invisible(NULL)
}
