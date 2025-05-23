% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotSignifBarcodeProportion.R
\name{plotSignifBarcodeProportion}
\alias{plotSignifBarcodeProportion}
\title{Plot Barcode contributions as mean Barcode proportion across samples}
\usage{
plotSignifBarcodeProportion(barbieQ)
}
\arguments{
\item{barbieQ}{A \code{SummarizedExperiment} object created by the \link{createBarbieQ} function,
updated with Barcode test results by calling the \link{testBarcodeSignif} function.}
}
\value{
A \code{ggplot} S3 class object displaying Barcode contributions in a
bar plot.
}
\description{
\code{plotBarcodeProportion()} visualizes the overall proportion of Barcodes
within groups of significance, as determined by \link{testBarcodeSignif}, for each sample.
}
\examples{
Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
Treat <- factor(rep(c('ctrl', 'drug'), each = 6))
Time <- rep(rep(seq_len(2), each = 3), 2)
nbarcodes <- 50
nsamples <- 12
count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
rownames(count) <- paste0('Barcode', seq_len(nbarcodes))
barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))
testBB <- testBarcodeSignif(barbieQ, sampleGroup = 'Treat')
plotSignifBarcodeProportion(testBB)
}
