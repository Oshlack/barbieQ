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
  barcodeList <- list(group1 = c("Barcode1", "Barcode2", "Barcode3", "BarcodeOdd"))
  corBarbieQ <- extractBarcodePairs(myBarbieQ, preDefinedCluster = barcodeList)
  corDF <- SummarizedExperiment::rowData(corBarbieQ)$barcodeCorrelation
  knownPairDf <- S4Vectors::metadata(corDF)$preDefinedBarcodePair
  expect_equal(nrow(knownPairDf), choose(4, 2))
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
  p <-plotBarcodePairCorrelation(myBarbieQ, preDefinedCluster = barcodeArray, showRawProportion = T)
  ## confirm barcodes not-exsiting in the barbieQ will not trigger errors
  barcodeList <- list(group1 = c("Barcode1", "Barcode2", "Barcode3", "BarcodeOdd"))
  p <- plotBarcodePairCorrelation(myBarbieQ, preDefinedCluster = barcodeList)
  expect_equal(sum(p$data$correlationGroup == "pre-Defined"), choose(3, 2))
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
  ## confirm barcodes not-exsiting in the barbieQ will not trigger errors
  barcodeList <- list(group1 = c("Barcode1", "Barcode2", "Barcode3", "BarcodeOdd"))
  expect_message(
    clusterCorrelatingBarcodes(myBarbieQ, preDefinedCluster = barcodeList),
    "identified 1 clusters, including 3 Barcodes."
  )
})

test_that("plotting barcode cluster structure works", {
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) %>% abs()
  rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
  myBarbieQ <- createBarbieQ(count)
  ## add known Barcode groups
  barcodeArray <- c(seq_len(40), rep(41, 10))
  myBarbieQ <- clusterCorrelatingBarcodes(myBarbieQ, preDefinedCluster = barcodeArray)
  ## expect all elements to be ggplot objects
  p_list <- inspectCorrelatingBarcodes(myBarbieQ)
  sapply(p_list, function(p) expect_s3_class(p, "ggplot"))
  ## confirm barcodes not-exsiting in the barbieQ will not trigger errors
  barcodeList <- list(group1 = c("Barcode1", "Barcode2", "Barcode3", "BarcodeOdd"))
  myBarbieQ <- clusterCorrelatingBarcodes(myBarbieQ, preDefinedCluster = barcodeList)
  p_list <- inspectCorrelatingBarcodes(myBarbieQ)
  sapply(p_list, function(p) expect_s3_class(p, "ggplot"))
})

test_that("detailed plotting of barcode cluster proportions works", {
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) %>% abs()
  rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
  myBarbieQ <- createBarbieQ(count)
  ## add known barcode groups
  barcodeArray <- c(seq_len(46), rep(47, 2), rep(49, 2))
  myBarbieQ <- clusterCorrelatingBarcodes(myBarbieQ, preDefinedCluster = barcodeArray)
  
  ## expect a ggplot object with default arguments
  p <- inspectCorrelatingBarcodesDetailed(myBarbieQ)
  expect_s3_class(p, "ggplot")
  
  ## expect ggplot with each transformation
  p_asin <- inspectCorrelatingBarcodesDetailed(myBarbieQ, transformation = "asin-sqrt")
  expect_s3_class(p_asin, "ggplot")
  
  p_logit <- inspectCorrelatingBarcodesDetailed(myBarbieQ, transformation = "logit")
  expect_s3_class(p_logit, "ggplot")
  
  ## expect ggplot with manual ncol
  p_ncol <- inspectCorrelatingBarcodesDetailed(myBarbieQ, ncol = 1)
  expect_s3_class(p_ncol, "ggplot")
  
  ## expect error when cluster info is missing
  myBarbieQ_empty <- createBarbieQ(count)
  expect_error(
    inspectCorrelatingBarcodesDetailed(myBarbieQ_empty),
    "Cluster information does not exist"
  )
  
  ## expect error on invalid transformation
  expect_error(
    inspectCorrelatingBarcodesDetailed(myBarbieQ, transformation = "invalid"),
    "transformation must be one of"
  )
  
  ## confirm barcodes not-exsiting in the barbieQ will not trigger errors
  barcodeList <- list(group1 = c("Barcode1", "Barcode2", "Barcode3", "BarcodeOdd"))
  myBarbieQ <- clusterCorrelatingBarcodes(myBarbieQ, preDefinedCluster = barcodeList)
  p <- inspectCorrelatingBarcodesDetailed(myBarbieQ)
  expect_s3_class(p, "ggplot")
})

test_that("merging correlating barcodes works", {
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) %>% abs()
  rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
  myBarbieQ <- createBarbieQ(count)
  ## add known Barcode groups
  barcodeArray <- c(seq_len(40), rep(41, 5), rep(42,5))
  myBarbieQ <- clusterCorrelatingBarcodes(myBarbieQ, preDefinedCluster = barcodeArray)
  ## merge correlating barcodes
  myBarbieQ_merged <- mergeCorrelatingBarcodes(myBarbieQ)
  ## expect the unchaned barcodes
  expect_equal(rownames(myBarbieQ_merged)[seq_len(40)], paste0("Barcode", seq_len(40)))
  ## expect the representative barcode with max mean CPM
  expect_equal(
    rownames(myBarbieQ_merged)[41], 
    rowMeans(myBarbieQ@assays@data$CPM)[seq(41,45)] %>% which.max() %>% names())
  expect_equal(
    rownames(myBarbieQ_merged)[42], 
    rowMeans(myBarbieQ@assays@data$CPM)[seq(46,50)] %>% which.max() %>% names())
  ## test merge a different barbieQ object
  count2 <- matrix(rnorm(2 * nsamples), 2, nsamples) %>% abs()
  rownames(count2) <- paste0("Barcode", seq(51,52))
  myBarbieQ2 <- createBarbieQ(rbind(count, count2))
  myBarbieQ2_merged <- mergeCorrelatingBarcodes(
    barbieQ_clustered = myBarbieQ, barbieQ_toMerge = myBarbieQ2)
  ## expect the unchaned barcodes
  expect_equal(rownames(myBarbieQ2_merged)[seq_len(40)], paste0("Barcode", seq_len(40)))
  expect_equal(rownames(myBarbieQ2_merged)[seq(43,44)], paste0("Barcode", seq(51,52)))
  ## expect the representative barcode with max mean CPM
  expect_equal(
    rownames(myBarbieQ2_merged)[41], 
    rowMeans(myBarbieQ@assays@data$CPM)[seq(41,45)] %>% which.max() %>% names())
  expect_equal(
    rownames(myBarbieQ2_merged)[42], 
    rowMeans(myBarbieQ@assays@data$CPM)[seq(46,50)] %>% which.max() %>% names())
  
  ## confirm barcodes not-exsiting in the barbieQ will not trigger errors
  barcodeList <- list(group1 = c("Barcode1", "Barcode2", "Barcode3", "BarcodeOdd"))
  myBarbieQ <- clusterCorrelatingBarcodes(myBarbieQ, preDefinedCluster = barcodeList)
  ## merge correlating barcodes
  myBarbieQ_merged <- mergeCorrelatingBarcodes(myBarbieQ)
  ## expect the unchaned barcodes
  expect_equal(rownames(myBarbieQ_merged), paste0("Barcode", seq(3,50)))
})
