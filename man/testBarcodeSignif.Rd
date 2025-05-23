% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/testBarcodeSignif.R
\name{testBarcodeSignif}
\alias{testBarcodeSignif}
\title{Test significant change in Barcode proportion or occurrence
between sample groups}
\usage{
testBarcodeSignif(
  barbieQ,
  method = "diffProp",
  sampleMetadata = NULL,
  sampleGroup = NULL,
  contrastFormula = NULL,
  designFormula = NULL,
  designMatrix = NULL,
  block = NULL,
  transformation = "asin-sqrt",
  regularization = "firth"
)
}
\arguments{
\item{barbieQ}{A \code{barbieQ} object created by the \link{createBarbieQ} function.}

\item{method}{A string specifying what is to be tested.
Options: 'diffProp' and 'diffOcc'. Defaults to 'diffProp'.}

\item{sampleMetadata}{A \code{matrix} or \code{data.frame} of sample conditions,
where each factor is represented in a separate column. Defaults to NULL,
in which case sample conditions are inherited from \code{barbieQ$metadata}.}

\item{sampleGroup}{A string representing the name of a factor from the
sample conditions passed by \code{barbieQ} or \code{sampleMetadata}, or a vector of
sample conditions, indicating the primary factor to be tested.
Defaults to the first factor in the sample conditions.}

\item{contrastFormula}{A string indicating the contrast between sample conditions.
Defaults to contrast bewteen the original levels of conditions in the specified factor.}

\item{designFormula}{A formula to compute the \code{designMatrix}, generated by
the \link[stats:formula]{stats::as.formula} function. Defaults to include all factors
in the sample conditions.}

\item{designMatrix}{A numeric matrix standardizing \code{sampleMetadata}, generated by
the \link[stats:model.matrix]{stats::model.matrix} function. Defaults to be generated by
\code{designFormula}.}

\item{block}{A vector (array) indicating sample duplicates. Defaults to
no duplicates among the samples.}

\item{transformation}{A string specifying the transformation method when
testing 'diffProp'. Options include: 'asin-sqrt', 'logit', and 'none'.
Defaults to 'asin-sqrt'.}

\item{regularization}{A string specifying the regularization method when
testing 'diffOcc'. Options: 'firth' and 'none'. Defaults to 'firth'.}
}
\value{
A \code{barbieQ} object updated with a \code{testBarcodes} component, adding
a \code{list} named by the test name, containing:
\itemize{
\item \code{results}: a \code{data.frame} of statistical test results,
including p-values, etc.
\item \code{methods}: a \code{list} indicating the specific statistical test method used.
\item \code{targets}: a numeric \code{matrix} of the standardized design matrix
used in the test.
}
}
\description{
\code{testBarcodeSignif()} tests differential proportion (diffProp) or
occurrence (diffOcc) for each Barcode between sample groups,
with the option account for multiple factors using regression models.
The results can be visualized by \link{plotBarcodePValue},
\link{plotSignifBarcodeHeatmap}, \link{plotBarcodeMA}, and \link{plotSignifBarcodeProportion}.
}
\examples{
Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
Treat <- factor(rep(seq_len(2), each = 6))
Time <- rep(rep(seq_len(2), each = 3), 2)
nbarcodes <- 50
nsamples <- 12
count <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
rownames(count) <- paste0('Barcode', seq_len(nbarcodes))
barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))
testBarcodeSignif(barbieQ, sampleGroup = 'Treat')
}
