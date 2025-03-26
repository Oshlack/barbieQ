test_that("extracting count matrix, sampleMetadata file and factor colors works", {
  ## Sample conditions and color palettes
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
  barcodeCount <- abs(matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples))
  rownames(barcodeCount) <- paste0("Barcode", seq_len(nbarcodes))
  colnames(barcodeCount) <- paste0("Sample", seq_len(nsamples))

  ## when `object` is a numeric matrix
  ## Passing `object` is necessary, but `sampleMetadata` and `factorColors` are optional
  object1 <- createBarbieQ(object = barcodeCount)
  ## object is S4 SE class
  expect_s4_class(object1, "SummarizedExperiment")
  ## correct Barcode count passed
  expect_equal(as.matrix(SummarizedExperiment::assay(object1)), barcodeCount, ignore_attr = TRUE)
  ## homogeneous sample conditions set up
  expect_equal(as.matrix(SummarizedExperiment::colData(object1)$sampleMetadata), rep(1, 12), ignore_attr = TRUE)
  object2 <- createBarbieQ(object = barcodeCount, sampleMetadata = sampleConditions)
  ## correct sample conditions passed
  expect_equal(SummarizedExperiment::colData(object2)$sampleMetadata, 
               S4Vectors::DataFrame(sampleConditions), ignore_attr = TRUE)
  object3 <- createBarbieQ(
    object = barcodeCount, sampleMetadata = sampleConditions,
    factorColors = conditionColor
  )
  ## correct color palettes passed
  expect_equal(S4Vectors::metadata(object3)$factorColors, conditionColor, ignore_attr = TRUE)

  ## when `object` is a `barbieQ` SE object
  object4 <- createBarbieQ(barcodeCount, sampleConditions, conditionColor)
  SummarizedExperiment::colData(object4)$foo <- list(tester = seq_len(12))
  object5 <- createBarbieQ(object = object4)
  ## correct Barcode count, sample condition and color palettes inherited
  expect_equal(SummarizedExperiment::assay(object4), SummarizedExperiment::assay(object5))
  ## cannot add problematic `barbieQ` component
  object0 <- object4
  expect_error(SummarizedExperiment::assay(object0) <- SummarizedExperiment::assay(object0)[-1, ])
  ## row subsetting works
  object00 <- createBarbieQ(object = object0[-1,])
  expect_equal(colSums(SummarizedExperiment::assays(object00)$proportion), rep(1, 12), ignore_attr = TRUE)
  expect_equal(colSums(SummarizedExperiment::assays(object00)$CPM), rep(1e6, 12), ignore_attr = TRUE)
  expect_equal(dim(SummarizedExperiment::assays(object00)$occurrence), c(49, 12))
  expect_equal(max(SummarizedExperiment::assays(object00)$rank), 49)
  ## correct sample condition and color palettes inherited
  expect_equal(SummarizedExperiment::colData(object4)$sampleMetadata, 
               SummarizedExperiment::colData(object5)$sampleMetadata)
  expect_equal(S4Vectors::metadata(object4)$factorColors, 
               S4Vectors::metadata(object5)$factorColors)
  ## correct other components inherited
  expect_equal(SummarizedExperiment::colData(object5)$foo[[1]], seq_len(12), ignore_attr = TRUE)
  ## Updating a `barbieQ` object by passing new `sampleMetadata` and `factorColors`
  object6 <- createBarbieQ(
    object = object4, sampleMetadata = data.frame(Mouse = rep(seq_len(4), each = 3))
  )
  ## success update sample conditions
  expect_equal(SummarizedExperiment::colData(object6)$sampleMetadata, 
               S4Vectors::DataFrame(Mouse = rep(seq_len(4), each = 3)), ignore_attr = TRUE)
  object7 <- createBarbieQ(
    object = object4, sampleMetadata = data.frame(Mouse = rep(seq_len(4), each = 3)),
    factorColors = list(
      Mouse = c("1" = "#111199", "2" = "#112200", "3" = "#441111", "4" = "#000000")
    )
  )
  ## success update color palettes
  expect_equal(
    S4Vectors::metadata(object7)$factorColors[[1]],
    c("1" = "#111199", "2" = "#112200", "3" = "#441111", "4" = "#000000")
  )
})
