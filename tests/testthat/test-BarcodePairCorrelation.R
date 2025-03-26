test_that("extracting Barcode pairs works", {
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) %>% abs()
  rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
  myBarbieQ <- createBarbieQ(count)
  corBarbieQ <- extractBarcodePairs(myBarbieQ)
  corDF <- SummarizedExperiment::rowData(corBarbieQ)$barcodeCorrelation
  ## confirm extracting a square correlation mat
  expect_equal(dim(corDF), c(50, 50))
  barcodeArray <- c(seq_len(40), rep(41, 10))
  ## confirm can take group array
  corBarbieQ <- extractBarcodePairs(myBarbieQ, preDefinedCluster = barcodeArray)
  corDF <- SummarizedExperiment::rowData(corBarbieQ)$barcodeCorrelation
  knownPairDf <- S4Vectors::metadata(corDF)$preDefinedBarcodePair
  expect_equal(nrow(knownPairDf), choose(10, 2))
  ## confirm can take list of groups
  barcodeList <- list(group1 = c("Barcode1", "Barcode2", "Barcode3"))
  corBarbieQ <- extractBarcodePairs(myBarbieQ, preDefinedCluster = barcodeList)
  corDF <- SummarizedExperiment::rowData(corBarbieQ)$barcodeCorrelation
  knownPairDf <- S4Vectors::metadata(corDF)$preDefinedBarcodePair
  expect_equal(nrow(knownPairDf), choose(3, 2))
})

test_that("plotting Barcode pairs works", {
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) %>% abs()
  rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
  myBarbieQ <- createBarbieQ(count)
  p <- plotBarcodePairCorrelation(myBarbieQ)
  ## checking method
  expect_error(plotBarcodePairCorrelation(myBarbieQ, method = "wrongName"))
  ## confirm object class
  expect_s3_class(p, "ggplot")
  expect_true(all(p$data$correlationGroup == "non-Corr"))
  ## confirm correct group info
  barcodeArray <- c(seq_len(40), rep(41, 10))
  p <- plotBarcodePairCorrelation(myBarbieQ, preDefinedCluster = barcodeArray)
  expect_equal(sum(p$data$correlationGroup == "pre-Defined"), choose(10, 2))
})

test_that("clustering Barcodes based on correlation works", {
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) %>% abs()
  rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
  myBarbieQ <- createBarbieQ(count)
  ## confirm correct cluster dimension
  myBarbieQ <- clusterCorrelatingBarcodes(myBarbieQ)
  expect_equal(
    SummarizedExperiment::rowData(myBarbieQ)$barcodeCorrelatedCluster$cluster, 
    -seq_len(50), ignore_attr = TRUE
    )
  ## confirm known Barcode groups are taken
  barcodeArray <- c(seq_len(40), rep(41, 10))
  myBarbieQ <- clusterCorrelatingBarcodes(myBarbieQ, preDefinedCluster = barcodeArray)
  expect_equal(
    SummarizedExperiment::rowData(myBarbieQ)$barcodeCorrelatedCluster$cluster, 
    c(-seq_len(40), rep(1, 10)), ignore_attr = TRUE
    )
  expect_message(
    clusterCorrelatingBarcodes(myBarbieQ, preDefinedCluster = barcodeArray),
    "identified 1 clusters, including 10 Barcodes."
    )
})
