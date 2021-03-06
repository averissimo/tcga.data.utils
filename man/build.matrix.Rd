% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{build.matrix}
\alias{build.matrix}
\title{Build list of matrices with different types of tissues}
\usage{
build.matrix(source.name, assay.data)
}
\arguments{
\item{source.name}{data type}

\item{assay.data}{data from MultiAssayExperiment}
}
\value{
data divided by tissue type, with clinical and source.name assay
}
\description{
Build list of matrices with different types of tissues
}
\examples{
data('multiAssay')
fpkm.data <- build.matrix('RNASeqFPKM', multiAssay)
fpkm.per.tissue <- fpkm.data$data
fpkm.clinical   <- fpkm.data$clinical
names(fpkm.per.tissue)
}
