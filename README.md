# methylclock

Estimation of chronological, gestational and biological DNA methylation
(DNAm) age, trait scores and epigenetic age acceleration from methylation
beta values.

Version 2.0 is a rewrite of the package: 44 published clocks selectable by
name, family or array platform; input as matrices, data frames, Bioconductor
containers or HDF5 files, with block-wise out-of-core computation through
[BigDataStatMeth](https://cran.r-project.org/package=BigDataStatMeth) when
the data exceed memory; configurable missing-value imputation (mean,
reference or KNN), cell-type deconvolution, intrinsic and extrinsic
epigenetic age acceleration (IEAA, EEAA), quality-control reports and
plotting functions. See `NEWS.md` for the full changelog.

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

Clock coefficients are distributed through the
[methylclockData](https://bioconductor.org/packages/methylclockData/)
package (ExperimentHub), with a Zenodo fallback, so no manual data setup is
needed.

## Usage

```r
library(methylclock)
data(methylclock_betas)
res <- methylclock(methylclock_betas, clocks = c("Horvath", "Levine"))
head(as.data.frame(res))
```

The vignette (`vignette("methylclock")`) documents the clock catalogue, the
input formats, imputation, age acceleration, quality control and the
plotting suite.

## Citation

Dolors Pelegri-Siso, Paula de Prado, Justiina Ronkainen, Mariona Bustamante,
Juan R Gonzalez. methylclock: a Bioconductor package to estimate DNA
methylation age. *Bioinformatics* 37(12):1759-1760, 2021.
<https://doi.org/10.1093/bioinformatics/btaa825>

## License

MIT. Individual clock coefficient sets keep their original licenses, stated
per clock in the catalogue (`clock_catalog`) and in the methylclockData
documentation.
