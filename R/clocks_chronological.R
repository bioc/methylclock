# Chronological-age clocks (1st generation). Metadata is verified; `transform`
# defaults to identity unless a clock needs one, and `predictor` names the engine
# family.

.register_chronological <- function() {
    clock_register("Horvath", target = "chronological", generation = 1L,
        tissue = "multi-tissue",
        platform = c("27K", "450K"), resource = "coefHorvath",
        transform = anti_trafo, citation = "horvath2013")
    clock_register("Hannum", target = "chronological", generation = 1L,
        tissue = "whole blood",
        platform = "450K", resource = "coefHannum", intercept = FALSE,
        citation = "hannum2013")
    clock_register("skinHorvath", target = "chronological", generation = 1L,
        tissue = "skin and blood",
        platform = c("450K", "EPIC"), resource = "coefSkin",
        transform = anti_trafo, citation = "horvath2018")
    # INTERNAL: skinHorvath -> D-01 (column label under cell.count).
    clock_register("PedBE", target = "chronological", generation = 1L,
        tissue = "buccal epithelium (children)",
        platform = "450K", resource = "coefPedBE",
        transform = anti_trafo, citation = "mcewen2020")
    clock_register("Wu", target = "chronological", generation = 1L,
        tissue = "whole blood (children)",
        platform = c("27K", "450K"), resource = "coefWu",
        transform = function(x) anti_trafo(x, adult_age = 48) / 12,
        citation = "wu2019")
    # INTERNAL: Wu -> D-08. Using adult.age = 48 (oracle B, the likely-correct
    # INTERNAL: value); deviates from oracle A (default 20). Confirm vs paper.
    clock_register("EN", target = "chronological", generation = 1L,
        tissue = "whole blood",
        platform = "450K", resource = "coefEN", standardize = TRUE,
        citation = "zhang2019")
    clock_register("BLUP", target = "chronological", generation = 1L,
        tissue = "whole blood",
        platform = "450K", resource = "coefBlup", standardize = TRUE,
        citation = "zhang2019")
    clock_register("BNN", target = "chronological", generation = 1L,
        tissue = "whole blood",
        predictor = "nn", platform = "450K", resource = NULL,
        aux = "coefHorvath", source = "compiled", license = "internal",
        citation = "alfonso2020bayesian")
    # INTERNAL: BNN weights compiled in src/bnn*.cpp (codegen, VLAs fixed);
    # INTERNAL: features = the Horvath CpGs (aux), fed CpGs x samples to the
    # INTERNAL: baked-in bnn_predict().
    clock_register("AltumAge", target = "chronological", generation = 1L,
        tissue = "multi-tissue",
        predictor = "nn", platform = c("27K", "450K", "EPIC"),
        resource = "coefAltumAge", aux = "coefAltumAgeRef",
        license = "MIT", source = "biolearn / AltumAge (de Lima Camillo 2022)",
        citation = "delimacamillo2022")
    # INTERNAL: nn forward in src/altumage.cpp (Eigen), weights from HDF5.
    clock_register("VidalBralo", target = "chronological", generation = 1L,
        tissue = "whole blood",
        platform = c("27K", "450K"), resource = "coefVidalBralo",
        license = "CC-BY", source = "Vidal-Bralo 2016 Front Genet 7:126",
        citation = "vidalbralo2016")
    clock_register("Lin", target = "chronological", generation = 1L,
        tissue = "whole blood",
        platform = "450K", resource = "coefLin",
        license = "CC-BY", source = "Lin 2016 Aging 8:394 (Suppl Table 1)",
        citation = "lin2016")
    clock_register("Weidner", target = "chronological", generation = 1L,
        tissue = "whole blood",
        platform = c("27K", "450K"), resource = "coefWeidner",
        license = "CC-BY", source = "Weidner 2014 Genome Biol 15:R24 (3-CpG)",
        citation = "weidner2014")
    # Trained only on probes shared by the 450K, EPICv1 and EPICv2 arrays, so
    # it runs natively on all three; its table weights the square of some
    # probes' betas as well as the betas (`power` column).
    clock_register("Garma", target = "chronological", generation = 1L,
        tissue = "whole blood",
        platform = c("450K", "EPIC", "EPICv2"), resource = "coefGarma",
        license = "CC-BY-NC-ND (non-commercial, no derivatives)",
        source = "Garma & Quintela-Fandino 2024 Genome Med 16:116 (Table S5)",
        citation = "garma2024")
    invisible(NULL)
}
