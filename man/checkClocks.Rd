% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkClocks.R
\name{checkClocks}
\alias{checkClocks}
\title{Check wheter input data contains the required 
CpGs for the implemented clocks.}
\usage{
checkClocks(x, ...)
}
\arguments{
\item{x}{data.frame or tibble (Individual in columns, 
CpGs in rows, CpG names in first colum - i.e. Horvath's format), 
ExpressionSet or GenomicRatioSet. 
A matrix is also possible having the CpG names in the rownames.}

\item{...}{other parameters}
}
\value{
a list with the different clocks when there are more than 80% of 
the required CpGs
}
\description{
Check wheter input data contains the required 
CpGs for the implemented clocks.
}
\details{
To be supplied
}
\examples{
TestDataset <- get_TestDataset()
checkClocks(TestDataset)
}
