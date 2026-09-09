# methylclock

`methylclock` estimates chronological, gestational and biological DNA
methylation (DNAm) age, trait scores and epigenetic age acceleration from
methylation beta values.

Version 2.0 is a rewrite of the package around a clock registry: every clock
is a declarative entry describing what it estimates, which CpGs and weights
it uses, and how its output is transformed, so all 44 clocks run through the
same engines and new clocks can be added without touching the computation
code.

## Overview

- **44 published clocks** across eight families (see the tables below),
  selectable by name, family or array platform, including the
  EPICv2-native cross-platform Garma clock.
- **Flexible input**: a matrix or data frame of betas, M-values (detected
  and converted), `SummarizedExperiment`/`ExpressionSet`/minfi objects, or
  an HDF5 file. Arrays larger than memory are computed **out-of-core**,
  block by block, through
  [BigDataStatMeth](https://cran.r-project.org/package=BigDataStatMeth).
- **Missing-value imputation**: per-CpG mean, reference values from the
  clock's own training data, none, or KNN (compiled code, and block-wise on
  disk-backed input).
- **Cell-type deconvolution** (`cellCounts()`, Houseman's reference-based
  method) and **epigenetic age acceleration**: `ageAcceleration()` for
  residual and cell-adjusted measures, and the canonical `EEAA()` and
  `IEAA()` of Chen et al. (2016), plus `ageAccelerationChen()`
  (Klemera-Doubal weighting).
- **Quality control**: `clockCoverage()` (which clocks your data can
  support), `clockAccuracy()` (agreement with chronological age) and
  `qcReport()` (an HTML report with a written reading of each section).
- **Plots**: predicted-vs-age panels, clock correlations, distributions and
  densities by group, acceleration heatmaps, Bland-Altman, forest/dumbbell
  group differences, per-sample discordance, reference percentile bands and
  longitudinal trajectories. All return `ggplot` objects.
- **Longitudinal helpers**: `clockTrajectories()`, `trajectoryRates()` and
  `plotTrajectories()` for repeated measures per subject.

Clock coefficients are not stored in this package: they are distributed
through [methylclockData](https://bioconductor.org/packages/methylclockData/)
(ExperimentHub, with a Zenodo fallback), so no manual data setup is needed.

## Installation

From Bioconductor devel (R devel required):

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("methylclock", version = "devel")
```

From GitHub (the `devel` branch holds the current source):

```r
BiocManager::install("methylclockData")
BiocManager::install("isglobal-brge/methylclock@devel")
```

## Quick start

```r
library(methylclock)

# A small bundled beta matrix (919 CpGs x 16 samples)
data(methylclock_betas)

# Every age clock the data cover; returns a "methylclock" object
res <- methylclock(methylclock_betas)
head(as.data.frame(res))

# Specific clocks, with KNN imputation of missing values
res <- methylclock(methylclock_betas, clocks = c("Horvath", "Levine"),
                   impute = "knn")

# Which clocks can your data support?
clockCoverage(methylclock_betas)
```

Age acceleration and plotting, on a bundled cohort of precomputed
estimates (656 blood samples with age and sex):

```r
data(methylclock_demo)

# Predicted vs chronological age, per clock
plotDNAmAge(methylclock_demo, age = methylclock_demo$age,
            clocks = c("Horvath", "Hannum", "Levine"))

# Age acceleration (residual-based; EEAA() and IEAA() compute the
# canonical cell-weighted measures from the beta values)
acc <- ageAcceleration(methylclock_demo, age = methylclock_demo$age)
```

Gestational age clocks have their own wrapper, `DNAmGA()`, and trait
EpiScores are opt-in by name, e.g.
`methylclock(betas, clocks = "McCartney.Smoking")`.

## Implemented clocks

The bundled `clock_catalog` data frame holds this catalogue (CpG counts,
native arrays, training tissue, licence and citation per clock);
`clock_list()` queries the live registry.

**Chronological age**

| Clock | Estimates | Tissue | Arrays | CpGs | Reference |
|---|---|---|---|---|---|
| AltumAge | Age (years) | multi-tissue | 27K/450K/EPIC | 20318 | de Lima Camillo et al. (2022) |
| BLUP | Age (years) | whole blood | 450K | 318395 | Zhang et al. (2019) |
| BNN | Age (years) | whole blood | 450K | 353 | Alfonso et al. (2020) |
| EN | Age (years) | whole blood | 450K | 514 | Zhang et al. (2019) |
| Garma | Age (years) | whole blood | 450K/EPIC/EPICv2 | 4962 | Garma et al. (2024) |
| Hannum | Age (years) | whole blood | 450K | 69 | Hannum et al. (2013) |
| Horvath | Age (years) | multi-tissue | 27K/450K | 353 | Horvath (2013) |
| Lin | Age (years) | whole blood | 450K | 99 | Lin et al. (2016) |
| NEOaPMA450K | Age (years) | buccal (preterm infants) | 450K | 408 | Graw et al. (2021) |
| NEOaPMAEPIC | Age (years) | buccal (preterm infants) | EPIC | 520 | Graw et al. (2021) |
| NEOaPNA450K | Age (years) | buccal (preterm infants) | 450K | 302 | Graw et al. (2021) |
| NEOaPNAEPIC | Age (years) | buccal (preterm infants) | EPIC | 508 | Graw et al. (2021) |
| PedBE | Age (years) | buccal epithelium (children) | 450K | 94 | McEwen et al. (2020) |
| skinHorvath | Age (years) | skin and blood | 450K/EPIC | 391 | Horvath et al. (2018) |
| VidalBralo | Age (years) | whole blood | 27K/450K | 8 | Vidal-Bralo et al. (2016) |
| Weidner | Age (years) | whole blood | 27K/450K | 3 | Weidner et al. (2014) |
| Wu | Age (years) | whole blood (children) | 27K/450K | 111 | Wu et al. (2019) |

**Gestational age**

| Clock | Estimates | Tissue | Arrays | CpGs | Reference |
|---|---|---|---|---|---|
| Bohlin | Gestational age (weeks) | cord blood | 450K | 96 | Bohlin et al. (2016) |
| EPIC | Gestational age (weeks) | cord blood | EPIC | 176 | Haftorn et al. (2021) |
| Knight | Gestational age (weeks) | cord blood | 450K | 148 | Knight et al. (2016) |
| Lee.CPC | Gestational age (weeks) | placenta | 450K/EPIC | 1125 | Lee et al. (2019) |
| Lee.refRPC | Gestational age (weeks) | placenta | 450K/EPIC | 1125 | Lee et al. (2019) |
| Lee.RPC | Gestational age (weeks) | placenta | 450K/EPIC | 1125 | Lee et al. (2019) |
| Mayne | Gestational age (weeks) | placenta | 450K | 62 | Mayne et al. (2017) |

**Phenotypic age**

| Clock | Estimates | Tissue | Arrays | CpGs | Reference |
|---|---|---|---|---|---|
| Levine | Biological age (years) | whole blood | 450K | 513 | Levine et al. (2018) |

**Pace of aging**

| Clock | Estimates | Tissue | Arrays | CpGs | Reference |
|---|---|---|---|---|---|
| DunedinPACE | Pace of aging (rate) | whole blood | 450K/EPIC | 173 | Belsky et al. (2022) |

**Telomere length**

| Clock | Estimates | Tissue | Arrays | CpGs | Reference |
|---|---|---|---|---|---|
| TL | Telomere length | whole blood | 450K | 140 | Lu et al. (2019) |

**Mitotic history**

| Clock | Estimates | Tissue | Arrays | CpGs | Reference |
|---|---|---|---|---|---|
| epiTOC1 | Cell divisions (count) | multi-tissue (proliferative) | 450K/EPIC | 385 | Yang et al. (2016) |
| epiTOC2 | Cell divisions (count) | multi-tissue (proliferative) | 450K/EPIC | 163 | Teschendorff (2020) |
| HypoClock | Cell divisions (count) | multi-tissue (proliferative) | 450K/EPIC | 678 | Teschendorff (2020) |
| RepliTali | Cell divisions (count) | cultured cells | 450K/EPIC | 87 | Endicott et al. (2022) |
| stemTOC | Cell divisions (count) | multi-tissue (proliferative) | 450K/EPIC | 371 | Zhu et al. (2024) |

**Causal age components**

| Clock | Estimates | Tissue | Arrays | CpGs | Reference |
|---|---|---|---|---|---|
| AdaptAge | Age (years, causal) | whole blood | 450K/EPIC | 998 | Ying et al. (2024) |
| CausAge | Age (years, causal) | whole blood | 450K/EPIC | 581 | Ying et al. (2024) |
| DamAge | Age (years, causal) | whole blood | 450K/EPIC | 1089 | Ying et al. (2024) |

**Trait EpiScores** (opt-in: request by name)

| Clock | Estimates | Tissue | Arrays | CpGs | Reference |
|---|---|---|---|---|---|
| McCartney.Alcohol | Trait score | whole blood | 450K/EPIC | 450 | McCartney et al. (2018) |
| McCartney.BMI | Trait score | whole blood | 450K/EPIC | 1109 | McCartney et al. (2018) |
| McCartney.BodyFat | Trait score | whole blood | 450K/EPIC | 968 | McCartney et al. (2018) |
| McCartney.Education | Trait score | whole blood | 450K/EPIC | 373 | McCartney et al. (2018) |
| McCartney.HDL | Trait score | whole blood | 450K/EPIC | 737 | McCartney et al. (2018) |
| McCartney.LDL | Trait score | whole blood | 450K/EPIC | 233 | McCartney et al. (2018) |
| McCartney.Smoking | Trait score | whole blood | 450K/EPIC | 233 | McCartney et al. (2018) |
| McCartney.TotalChol | Trait score | whole blood | 450K/EPIC | 204 | McCartney et al. (2018) |
| mCigarette | Trait score | whole blood | EPIC | 1255 | Chybowska et al. (2025) |

## Working with large arrays

An HDF5-backed matrix can be passed directly, and `mc_to_hdf5()` converts a
matrix or a delimited text file into one. Clocks then read only the CpGs
they need, in file order, and whole-array clocks (BLUP, EN) run block-wise
under a configurable memory budget, so memory use follows the budget rather
than the size of the array.

```r
hm <- mc_to_hdf5("betas_big.csv", "betas.h5")
res <- methylclock(hm)
close(hm)
```

## Documentation

The vignette walks through the clock catalogue, input formats, imputation,
age acceleration, quality control, group comparisons, reference ranges and
longitudinal analysis:

```r
vignette("methylclock")
```

## Citation

Dolors Pelegri-Siso, Paula de Prado, Justiina Ronkainen, Mariona Bustamante,
Juan R Gonzalez. methylclock: a Bioconductor package to estimate DNA
methylation age. *Bioinformatics* 37(12):1759-1760, 2021.
<https://doi.org/10.1093/bioinformatics/btaa825>

## License

MIT. Individual clock coefficient sets keep their original licenses, stated
per clock in `clock_catalog` and in the methylclockData documentation.
