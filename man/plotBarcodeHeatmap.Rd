% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotBarcodeHeatmap.R
\name{plotBarcodeHeatmap}
\alias{plotBarcodeHeatmap}
\title{Plot Barcode CPM (or occurrence) across samples in a Heatmap}
\usage{
plotBarcodeHeatmap(
  barbieQ,
  value = "CPM",
  splitSamples = FALSE,
  targets = NULL,
  sampleGroups = NULL,
  barcodeAnnotation = NULL,
  sampleAnnotation = NULL
)
}
\arguments{
\item{barbieQ}{A \code{barbieQ} object created by the \link{createBarbieQ} function.}

\item{value}{A string indicating what to visualize.
Defaults to 'CPM'. Options include: 'CPM' and 'occurrence'.}

\item{splitSamples}{A logical value deciding whether to split samples
into slices. Defaults to FALSE.}

\item{targets}{A \code{matrix} or \code{data.frame} of sample conditions,
where each factor is represented in a separate column. Defaults to NULL,
in which case sample conditions are inherited from \code{barbieQ$metadata}.}

\item{sampleGroups}{A string representing the name of a factor from the
sample conditions passed by \code{barbieQ} or \code{targets}, or a vector of
sample conditions, indicating the primary factor to split sample slices.}

\item{barcodeAnnotation}{A row Annotation object created by the
\link[ComplexHeatmap:rowAnnotation]{ComplexHeatmap::rowAnnotation} function. Defaults to NULL, which means
no Barcode annotation will be displayed in the Heatmap.}

\item{sampleAnnotation}{A column Annotation object created by the
\link[ComplexHeatmap:HeatmapAnnotation]{ComplexHeatmap::HeatmapAnnotation} function. Defaults to NULL, which means
the sample annotations are generated from the sample conditions provided by
\code{barbieQ} and \code{targets}.}
}
\value{
A \code{Heatmap} S4 object displaying the heatmap of Barcode
data across samples, optionally annotated with sample and Barcode
information.
}
\description{
\code{plotBarcodeHeatmap()} visualizes Barcode output across samples with
the option to include sample annotations. The Heatmap can display either:
\itemize{
\item CPM: log2(CPM+1)
\item occurrence: 0 or 1
}
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
plotBarcodeHeatmap(myBarbieQ)
}
