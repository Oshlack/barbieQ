% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotBarcodePareto.R
\name{plotBarcodePareto}
\alias{plotBarcodePareto}
\title{Plot contributions of Barcodes split by \emph{top} v.s. \emph{bottom}
in a circular bar plot.}
\usage{
plotBarcodePareto(barbieQ, absoluteProportion = FALSE)
}
\arguments{
\item{barbieQ}{A \code{SummarizedExperiment} object created by the \link{createBarbieQ} function.}

\item{absoluteProportion}{A logical value indicating whether to present
absolute Barcode mean proportion (across samples) or relative values across Barcodes,
Defaults to FALSE,
which means it will present the percentage of (Barcode mean proportion) across Barcodes.}
}
\value{
A \code{ggplot} S3 class object displaying a circular bar plot,
highlighting the relative total proportion of each Barcode across samples.
}
\description{
After the \link{tagTopBarcodes} function tags Barcodes as either
\emph{top} or \emph{bottom}, \code{plotBarcodePareto()} visualizes the
proportion of each Barcode and separates them by these tags
in a circular bar plot, also known as a Pareto plot.
}
\note{
To save the plot with its original aspect ratio, use the
\code{ggplot2::ggsave} function and set \code{width = 8} and \code{height = 6},
or an equivalent size.
}
\examples{
## sample conditions and color palettes
sampleConditions <- data.frame(
  Treat = factor(rep(c("ctrl", "drug"), each = 6)),
  Time = rep(rep(seq_len(2), each = 3), 2)
)
conditionColor <- list(
  Treat = c(ctrl = "#999999", drug = "#112233"),
  Time = c("1" = "#778899", "2" = "#998877")
)
## Barcode count data
nbarcodes <- 50
nsamples <- 12
barcodeCount <- abs(matrix(10, nbarcodes, nsamples))
barcodeCount[seq(21, 50), ] <- 0.0001
rownames(barcodeCount) <- paste0("Barcode", seq_len(nbarcodes))
## create a `barbieQ` object
myBarbieQ <- createBarbieQ(barcodeCount, sampleConditions, conditionColor)
myBarbieQ <- tagTopBarcodes(myBarbieQ)
plotBarcodePareto(myBarbieQ)
}
