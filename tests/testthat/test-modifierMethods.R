test_that("subset barbieQ by metadata works", {
  myBarbieQ <- createBarbieQ(
    object = matrix(seq_len(20), nrow = 5, ncol = 4),
    target = data.frame(fac = c(1, 1, 2, 2))
  )
  subBarbieQ <- subsetSamplesByMetadata(
    barbieQ = myBarbieQ, factor = "fac", specifiedConditions = 1, keep = TRUE
  )
  expect_equal(subBarbieQ$metadata$fac, c(1, 1))
})

test_that("subset barbieQ samples works", {
  myBarbieQ <- createBarbieQ(
    object = matrix(seq_len(20), nrow = 5, ncol = 4),
    target = data.frame(fac = c(1, 1, 2, 2))
  )
  subBarbieQ <- subsetSamples(barbieQ = myBarbieQ, retainedColumns = c(TRUE, TRUE, FALSE, FALSE))
  expect_equal(subBarbieQ$metadata$fac, c(1, 1))
})

test_that("combining barbieQ samples works", {
  barbieQ1 <- createBarbieQ(
    object = matrix(seq_len(20), nrow = 5, ncol = 4),
    target = data.frame(
      a = c("a", "a", "b", "b"),
      b = c("x", "x", "x", "x")
    ),
    factorColors = list(
      a = c("a" = "#FFFFFF", "b" = "#CCCCFF"),
      b = c("x" = "#000099")
    )
  )
  barbieQ2 <- createBarbieQ(
    object = matrix(seq(21, 30), nrow = 5, ncol = 2),
    target = data.frame(
      a = c("a", "c"),
      c = c("u", "u")
    ),
    factorColors = list(
      a = c("a" = "#FFFFFF", "c" = "#000000"),
      c = c("u" = "#CCCCCC")
    )
  )
  barbieQ3 <- createBarbieQ(
    object = matrix(seq(31, 35), nrow = 5, ncol = 1),
    target = data.frame(b = c("x")),
    factorColors = list(b = c("x" = "#999999"))
  )

  combinedBarbieQ <- combineBarbieQs(barbieQ1, barbieQ2, barbieQ3)

  manualBarbieQ <- createBarbieQ(
    object = matrix(seq_len(35), nrow = 5, ncol = 7),
    target = data.frame(
      a = c("a", "a", "b", "b", "a", "c", ""),
      b = c("x", "x", "x", "x", "", "", "x"),
      c = c("", "", "", "", "u", "u", "")
    ),
    factorColors = list(
      a = c("a" = "#FFFFFF", "b" = "#CCCCFF", "c" = "#000000"),
      b = c("x" = "#000099", "x" = "#999999"),
      c = c("u" = "#CCCCCC")
    )
  )
  rownames(manualBarbieQ$metadata) <- rownames(combinedBarbieQ$metadata)
  expect_equal(combinedBarbieQ$metadata, manualBarbieQ$metadata)

  barbieQ4 <- createBarbieQ(
    object = matrix(seq(31, 34), nrow = 4, ncol = 1)
  )
  expect_error(combineBarbieQs(barbieQ1, barbieQ4))
})

test_that("merging lists works", {
  list1 <- list(
    a = c("a" = "#FFFFFF", "b" = "#CCCCFF"),
    b = c("x" = "#000099")
  )
  list2 <- list(
    a = c("a" = "#FFFFFF", "c" = "#000000"),
    c = c("u" = "#CCCCCC")
  )
  list3 <- list(b = c("x" = "#999999"))
  mergedList <- mergeLists(list1, list2, list3)
  expectdcList <- list(
    a = c("a" = "#FFFFFF", "b" = "#CCCCFF", "c" = "#000000"),
    b = c("x" = "#000099", "x" = "#999999"),
    c = c("u" = "#CCCCCC")
  )
  expect_equal(mergedList, expectdcList)
})

test_that("subset barbieQ Barcodes works", {
  myBarbieQ <- createBarbieQ(
    object = matrix(seq_len(20), nrow = 5, ncol = 4),
    target = data.frame(fac = c(1, 1, 2, 2))
  )
  subBarbieQ <- subsetBarcodes(
    barbieQ = myBarbieQ, retainedRows = c(TRUE, TRUE, FALSE, FALSE, FALSE)
  )
  expectedAssay <- as.data.frame(
    matrix(c(1, 2, 6, 7, 11, 12, 16, 17), nrow = 2, ncol = 4)
  )
  expect_equal(subBarbieQ$assay, expectedAssay)
})

test_that("collapsing barbieQ samples works", {
  myBarbieQ <- createBarbieQ(
    object = matrix(seq_len(30), nrow = 5, ncol = 6),
    target = data.frame(fac = letters[seq_len(6)])
  )
  collapsedBarbieQ <- collapseSamples(
    barbieQ = myBarbieQ, groupArray = c(8, 8, 8, 9, 9, 9)
  )
  expect_equal(collapsedBarbieQ$metadata[1, 1], "a.b.c")
})

test_that("collapsing columns works", {
  mymatrix <- matrix(seq_len(30), nrow = 5, ncol = 6)
  collapsedMat <- collapseColumnsByArray(
    mymatrix,
    groupArray = c(8, 8, 8, 9, 9, 9), method = mean
  )
  expectedMat <- matrix(c(seq(6, 10), seq(21, 25)), nrow = 5, ncol = 2)
  colnames(expectedMat) <- c(8, 9)
  expect_equal(collapsedMat, expectedMat)
})

test_that("collapsing barbieQ Barcodes works", {
  myBarbieQ <- createBarbieQ(
    object = matrix(seq_len(30), nrow = 6, ncol = 5),
    target = data.frame(fac = letters[seq_len(5)])
  )
  collapsedBarbieQ <- collapseBarcodes(
    barbieQ = myBarbieQ, groupArray = c(8, 8, 8, 9, 9, 9)
  )
  expect_equal(collapsedBarbieQ$metadata, myBarbieQ$metadata)
  expectMat <- matrix(
    c(3, 9, 15, 21, 27, 6, 12, 18, 24, 30),
    ncol = 2, nrow = 5
  ) |>
    t() |>
    as.data.frame()
  expect_equal(collapsedBarbieQ$assay, expectMat, ignore_attr = TRUE)
  myBarbieQ <- tagTopBarcodes(myBarbieQ)
  collapsedBarbieQ <- collapseBarcodes(
    barbieQ = myBarbieQ, groupArray = c(7, 8, 8, 9, 9, 9)
  )
  expect_equal(collapsedBarbieQ$isTop$vec, rep(TRUE, 3), ignore_attr = TRUE)
  expect_equal(
    collapsedBarbieQ$isTop$mat,
    matrix(TRUE, nrow = 3, ncol = 5),
    ignore_attr = TRUE
  )
})

test_that("collapsing rows works", {
  mymatrix <- matrix(seq_len(30), nrow = 6, ncol = 5)
  collapsedMat <- collapseRowsByArray(
    mymatrix,
    groupArray = c(8, 8, 8, 9, 9, 9), method = mean
  )
  expectedMat <- matrix(
    c(2, 8, 14, 20, 26, 5, 11, 17, 23, 29),
    nrow = 5, ncol = 2
  ) |> t()
  expect_equal(collapsedMat, expectedMat, ignore_attr = TRUE)
})
