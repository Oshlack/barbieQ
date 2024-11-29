test_that("converting first non-numeric column into rownames works", {
  mat <- data.frame(id = letters[seq_len(5)], matrix(seq_len(25), 5, 5))
  expectMat <- as.matrix(
    data.frame(matrix(seq_len(25), 5, 5), row.names = letters[seq_len(5)])
  )
  expect_equal(returnNumMat(mat), expectMat)
})

test_that("check barbieQ dimensions works", {
  Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
  Treat <- factor(rep(seq_len(2), each = 6))
  Time <- rep(rep(seq_len(2), each = 3), 2)
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) |> abs()
  rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
  barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))
  expect_true(checkBarbieQDimensions(barbieQ))
  barbieQ2 <- barbieQ
  barbieQ2$metadata <- barbieQ2$metadata[, -1, drop = FALSE]
  expect_true(checkBarbieQDimensions(barbieQ2))
  barbieQ3 <- barbieQ
  barbieQ3$metadata <- barbieQ3$metadata[-1, , drop = FALSE]
  expect_error(checkBarbieQDimensions(barbieQ3))
})

test_that("extracting targets and sampleGroups works", {
  Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
  Treat <- factor(rep(seq_len(2), each = 6))
  Time <- rep(rep(seq_len(2), each = 3), 2)
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) |> abs()
  rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
  barbieQ <- createBarbieQ(count, data.frame(Treat = Treat, Time = Time))

  expect_message(
    extractTargetsAndPrimaryFactor(barbieQ = barbieQ),
    "setting Treat as the primary effector of sample conditions."
  )
  expect_message(
    extractTargetsAndPrimaryFactor(
      barbieQ = barbieQ,
      targets = data.frame(Treat = Treat, Time = Time, Block = Block)
    ),
    "setting Treat as the primary effector of sample conditions."
  )
  expect_message(
    extractTargetsAndPrimaryFactor(barbieQ = barbieQ, sampleGroups = "Treat"),
    "setting Treat as the primary effector of sample conditions."
  )
  expect_equal(
    extractTargetsAndPrimaryFactor(
      barbieQ = barbieQ,
      targets = data.frame(Treat = Treat, Time = Time, Block = Block),
      sampleGroups = "Treat"
    ),
    list(
      mytargets = data.frame(Treat = Treat, Time = Time, Block = Block),
      pointer = 1
    )
  )
  expect_error(
    extractTargetsAndPrimaryFactor(barbieQ = barbieQ, sampleGroups = "treat")
  )
  expect_error(
    extractTargetsAndPrimaryFactor(barbieQ = barbieQ, sampleGroups = seq_len(6))
  )
})
