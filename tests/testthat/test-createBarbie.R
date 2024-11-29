test_that("extracting count matrix, target file and factor colors works", {
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

  ## when `object` is a numeric matrix
  ## Passing `object` is necessary, but `target` and `factorColors` are optional
  object1 <- createBarbieQ(object = barcodeCount)
  ## correct Barcode count passed
  expect_equal(as.matrix(object1$assay), barcodeCount, ignore_attr = TRUE)
  ## homogeneous sample conditions set up
  expect_equal(as.matrix(object1$metadata), rep(1, 12), ignore_attr = TRUE)
  object2 <- createBarbieQ(object = barcodeCount, target = sampleConditions)
  ## correct sample conditions passed
  expect_equal(object2$metadata, sampleConditions, ignore_attr = TRUE)
  object3 <- createBarbieQ(
    object = barcodeCount, target = sampleConditions,
    factorColors = conditionColor
  )
  ## correct color palettes passed
  expect_equal(object3$factorColors, conditionColor, ignore_attr = TRUE)

  ## when `object` is a `barbieQ` object
  object4 <- createBarbieQ(barcodeCount, sampleConditions, conditionColor)
  object4$foo <- list(tester = seq_len(12))
  object5 <- createBarbieQ(object = object4)
  ## correct Barcode count, sample condition and color palettes inherited
  expect_equal(object4$assay, object5$assay)
  ## cannot pass problematic `barbieQ` object
  object0 <- object4
  object0$assay <- object0$assay[-1, ]
  expect_error(createBarbieQ(object = object0))
  ## processed information updated
  object0$proportion <- object0$proportion[-1, ]
  object0$CPM <- object0$CPM[-1, ]
  object0$occurrence <- object0$occurrence[-1, ]
  object0$rank <- object0$rank[-1, ]
  object00 <- createBarbieQ(object = object0)
  expect_equal(colSums(object00$proportion), rep(1, 12), ignore_attr = TRUE)
  expect_equal(colSums(object00$CPM), rep(1e6, 12), ignore_attr = TRUE)
  expect_equal(dim(object00$occurrence), c(49, 12))
  expect_equal(max(object00$rank), 49)
  ## correct sample condition and color palettes inherited
  expect_equal(object4$metadata, object5$metadata)
  expect_equal(object4$factorColors, object5$factorColors)
  ## correct other components inherited
  expect_equal(object5$foo[[1]], seq_len(12), ignore_attr = TRUE)
  ## Updating a `barbieQ` object by passing new `target` and `factorColors`
  object6 <- createBarbieQ(
    object = object4, target = data.frame(Mouse = rep(seq_len(4), each = 3))
  )
  ## success update sample conditions
  expect_equal(object6$metadata, data.frame(Mouse = rep(seq_len(4), each = 3)))
  object7 <- createBarbieQ(
    object = object4, target = data.frame(Mouse = rep(seq_len(4), each = 3)),
    factorColors = list(
      Mouse = c("1" = "#111199", "2" = "#112200", "3" = "#441111", "4" = "#000000")
    )
  )
  ## success update color palettes
  expect_equal(
    object7$factorColors[[1]],
    c("1" = "#111199", "2" = "#112200", "3" = "#441111", "4" = "#000000")
  )
})
