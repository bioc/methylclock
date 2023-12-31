% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DNAmAge.R
\name{DNAmAge}
\alias{DNAmAge}
\title{DNAm age estimation using different DNA methylation clocks.}
\usage{
DNAmAge(
  x,
  clocks = "all",
  toBetas = FALSE,
  fastImp = FALSE,
  normalize = FALSE,
  age,
  cell.count = TRUE,
  cell.count.reference = "blood gse35069 complete",
  min.perc = 0.8,
  ...
)
}
\arguments{
\item{x}{data.frame (Individual in columns, CpGs in rows, CpG names in 
first colum - i.e. Horvath's format), matrix (individuals in columns and 
Cpgs in rows having CpG names in the rownames), ExpressionSet or 
GenomicRatioSet.}

\item{clocks}{the methods used for estimating DNAmAge. Currrently 
"Horvath", "Hannum", "Levine", "BNN", "skinHorvath", "PedBE", "Wu", "TL",
"BLUP", "EN" and "all" 
are available. Default is "all" and all clocks are estimated.}

\item{toBetas}{Should data be transformed to beta values? Default is FALSE. 
If TRUE, it implies data are M values.}

\item{fastImp}{Is fast imputation performed if necessary? (see details). 
Default is FALSE}

\item{normalize}{Is Horvath's normalization performed? By default is FALSE}

\item{age}{individual's chronological age.}

\item{cell.count}{Are cell counts estimated? Default is TRUE.}

\item{cell.count.reference}{Used when 'cell.count' is TRUE. Default is 
"blood gse35069 complete".See 'meffil::meffil.list.cell.count.references()' 
for possible values.}

\item{min.perc}{Indicates the minimum conicidence percentage required 
between CpGs in or dataframee x and CpGs in clock coefficients to perform 
the calculation. If min.prec is too low, the estimated gestational DNAm 
age can be poor}

\item{...}{Other arguments to be passed through impute package}
}
\value{
The estimated chronological and biological mDNA age
}
\description{
DNAm age estimation using different DNA methylation clocks.
}
\details{
Imputation is performed when having missing data.
Fast imputation is performed by ...
what about imputing only when CpGs for the clock are missing?
}
\examples{
MethylationData <- get_MethylationDataExample()
age.example55 <- DNAmAge(MethylationData)
}
