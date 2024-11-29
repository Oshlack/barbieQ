test_that("extracting Barcode pairs works", {
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) %>% abs()
  rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
  myBarbieQ <- createBarbieQ(count)
  pairInfo <- extractBarcodePairs(myBarbieQ)
  ## confirm extracting all pairs
  expect_equal(
    nrow(pairInfo$corTestResults),
    choose(nrow(myBarbieQ$assay), 2)
  )
  barcodeArray <- c(seq_len(40), rep(41, 10))
  ## confirm can take group array
  pairInfo <- extractBarcodePairs(myBarbieQ, BarcodeClusters = barcodeArray)
  expect_equal(nrow(pairInfo$knownPairDf), choose(10, 2))
  ## confirm can take list of groups
  barcodeList <- list(group1 = c("Barcode1", "Barcode2", "Barcode3"))
  pairInfo <- extractBarcodePairs(myBarbieQ, BarcodeClusters = barcodeList)
  expect_equal(nrow(pairInfo$knownPairDf), choose(3, 2))
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
  expect_true(all(p$data$knownCorrelating == FALSE))
  ## confirm correct group info
  barcodeArray <- c(seq_len(40), rep(41, 10))
  p <- plotBarcodePairCorrelation(myBarbieQ, BarcodeClusters = barcodeArray)
  expect_equal(sum(p$data$knownCorrelating), choose(10, 2))
})

test_that("clustering Barcodes based on correlation works", {
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) %>% abs()
  rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
  myBarbieQ <- createBarbieQ(count)
  ## confirm correct cluster dimension
  myBarbieQ <- clusterCorrelatingBarcodes(myBarbieQ)
  expect_equal(myBarbieQ$BarcodeCluster |> nrow(), 50)
  ## confirm known Barcode groups are taken
  barcodeArray <- c(seq_len(40), rep(41, 10))
  myBarbieQ <- clusterCorrelatingBarcodes(myBarbieQ, BarcodeClusters = barcodeArray)
  expect_equal(myBarbieQ$BarcodeCluster$corCluster, c(-seq_len(40), rep(1, 10)),
    ignore_attr = TRUE
  )
  expect_message(
    clusterCorrelatingBarcodes(myBarbieQ, BarcodeClusters = barcodeArray),
    "predicting 1 clusters, including 10 Barcodes."
  )
})
