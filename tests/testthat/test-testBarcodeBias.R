test_that("testing differential proportions works", {
  Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
  Treat <- factor(rep(c("ctrl", "drug"), each = 6))
  Time <- rep(rep(1:2, each = 3), 2)
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) |> abs()
  rownames(count) <- paste0("Barcode", 1:nbarcodes)
  barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))

  resultStat <- testDiffProp(
    barbieQ = barbieQ,
    mycontrasts = c(-1, 1, 0),
    contrastLevels = c("ctrl", "drug"),
    designMatrix = model.matrix(~ 0 + Treat + Time),
    transformation = "asin-sqrt",
    block = Block
  )
  expect_equal(rownames(resultStat), rownames(barbieQ$proportion))

  resultStat2 <- testDiffProp(
    barbieQ = barbieQ,
    mycontrasts = c(0, 0, 1),
    designMatrix = model.matrix(~ 0 + Treat + Time),
    transformation = "asin-sqrt",
    block = Block
  )
  expect_equal(rownames(resultStat2), rownames(barbieQ$proportion))
})

test_that("testing differential occurrence works", {
  Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
  Treat <- factor(rep(c("ctrl", "drug"), each = 6))
  Time <- rep(rep(1:2, each = 3), 2)
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) |> abs()
  rownames(count) <- paste0("Barcode", 1:nbarcodes)
  barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))

  resultStat <- testDiffOcc(
    barbieQ = barbieQ,
    regularization = "firth",
    mycontrasts = c(-1, 1, 0),
    contrastLevels = c("ctrl", "drug"),
    designMatrix = model.matrix(~ 0 + Treat + Time)
  )
  expect_equal(rownames(resultStat), rownames(barbieQ$occurrence))

  ## check one factor
  resultStat <- testDiffOcc(
    barbieQ = barbieQ,
    regularization = "firth",
    mycontrasts = c(-1, 1),
    contrastLevels = c("ctrl", "drug"),
    designMatrix = model.matrix(~ 0 + Treat)
  )

  resultStat <- testDiffOcc(
    barbieQ = barbieQ,
    regularization = "firth",
    mycontrasts = c(1),
    designMatrix = model.matrix(~ 0 + Time)
  )

  ## this will cause warings, because of
  resultStat <- testDiffOcc(
    barbieQ = barbieQ,
    regularization = "firth",
    mycontrasts = c(1, 0, 0),
    designMatrix = model.matrix(~ 0 + Time + Treat)
  )
})

test_that("barcode test extracting correct arguments, dispatching right function", {
  Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
  Treat <- factor(rep(c("ctrl", "drug"), each = 6))
  Time <- rep(rep(seq_len(2), each = 3), 2)
  nbarcodes <- 50
  nsamples <- 12
  ## Generate the base counts
  base_counts <- rnorm(nbarcodes, mean = 10, sd = 5) |> abs()
  count <- matrix(0, nbarcodes, nsamples)
  ## Decrease counts based on time
  for (i in seq_len(nbarcodes)) {
    count[i, ] <- base_counts[i] * (1 - (Time - 1) * 0.1)
  }
  rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
  barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))

  testBB1 <- testBarcodeBias(barbieQ, sampleGroups = "Treat")
  testBB11 <- testBarcodeBias(
    barbieQ,
    sampleGroups = "Treat",
    contrastLevels = c("ctrl", "drug")
  )
  testBB111 <- testBarcodeBias(
    barbieQ,
    sampleGroups = "Treat",
    contrastLevels = c("drug", "ctrl")
  )
  expect_equal(
    testBB11$testBarcodes$diffProp_Treat$results$direction,
    testBB111$testBarcodes$diffProp_Treat$results$direction
  )

  testBB1111 <- testBarcodeBias(barbieQ,
    sampleGroups = "Treat", contrastLevels = c("drug", "ctrl"),
    designFormula = formula("~0 + Treat + Time")
  )
  testBB11111 <- testBarcodeBias(barbieQ,
    sampleGroups = "Treat", contrastLevels = c("drug", "ctrl"),
    designMatrix = model.matrix(~ 0 + Treat + Time)
  )
  expect_equal(
    testBB1111$testBarcodes$diffProp_Treat$result,
    testBB11111$testBarcodes$diffProp_Treat$result
  )

  ## confirm `targets` is updated; design matrix and formula are updated accordingly
  testBB12 <- testBarcodeBias(barbieQ,
    sampleGroups = "Treat", contrastLevels = c("drug", "ctrl"),
    targets = data.frame(Treat = Treat)
  )
  expect_equal(colnames(testBB12$testBarcodes$diffProp_Treat$targets), "Treat")
  expect_equal(
    as.character(testBB12$testBarcodes$diffProp_Treat$methods$formula),
    "~0 + Treat"
  )

  ## with specified `designMatrix` used for test, should expect not using default formula
  testBB13 <- testBarcodeBias(barbieQ,
    sampleGroups = "Treat", contrastLevels = c("drug", "ctrl"),
    designMatrix = model.matrix(~ 0 + Treat)
  )
  expect_equal(testBB13$testBarcodes$diffProp_Treat$methods$formula, "NA")
  expect_equal(
    testBB13$testBarcodes$diffProp_Treat$methods$design,
    model.matrix(~ 0 + Treat),
    ignore_attr = TRUE
  )

  ## with specified `designFormula`, should expect updated formula and `designMatrix`
  testBB14 <- testBarcodeBias(barbieQ,
    sampleGroups = "Treat", contrastLevels = c("drug", "ctrl"),
    designFormula = formula("~ 0 + Treat")
  )
  expect_equal(
    testBB14$testBarcodes$diffProp_Treat$methods$design,
    model.matrix(~ 0 + Treat),
    ignore_attr = TRUE
  )

  testBB2 <- testBarcodeBias(barbieQ, sampleGroups = "Time")

  testBB3 <- testBarcodeBias(
    barbieQ,
    sampleGroups = "Time", method = "diffOcc",
    designFormula = formula("~ 0 + Time + Treat")
  )

  testBB4 <- testBarcodeBias(
    barbieQ,
    sampleGroups = "Time", method = "diffOcc",
    designFormula = formula("~ 0 + Time")
  )

  testBB <- testBarcodeBias(barbieQ, sampleGroups = rep(seq_len(4), each = 3))
})
