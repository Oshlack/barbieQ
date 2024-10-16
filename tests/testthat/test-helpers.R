test_that("converting first non-numeric column into rownames works", {
  mat <- data.frame(id = letters[seq_len(5)], matrix(seq_len(25), 5, 5))
  expectMat <- as.matrix(
    data.frame(matrix(seq_len(25), 5, 5), row.names = letters[seq_len(5)])
  )
  expect_equal(returnNumMat(mat), expectMat)
})

test_that("check Barbie dimensions works", {
  Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
  Treat <- factor(rep(seq_len(2), each = 6))
  Time <- rep(rep(seq_len(2), each = 3), 2)
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) |> abs()
  rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
  Barbie <- createBarbie(count, data.frame(Treat = Treat, Time = Time))
  expect_true(checkBarbieDimensions(Barbie))
  Barbie2 <- Barbie
  Barbie2$metadata <- Barbie2$metadata[, -1, drop = FALSE]
  expect_true(checkBarbieDimensions(Barbie2))
  Barbie3 <- Barbie
  Barbie3$metadata <- Barbie3$metadata[-1, , drop = FALSE]
  expect_error(checkBarbieDimensions(Barbie3))
})

test_that("extracting targets and sampleGroups works", {
  Block <- c(1, 1, 2, 3, 3, 4, 1, 1, 2, 3, 3, 4)
  Treat <- factor(rep(seq_len(2), each = 6))
  Time <- rep(rep(seq_len(2), each = 3), 2)
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes * nsamples), nbarcodes, nsamples) |> abs()
  rownames(count) <- paste0("Barcode", seq_len(nbarcodes))
  Barbie <- createBarbie(count, data.frame(Treat = Treat, Time = Time))

  expect_message(
    extractTargetsAndPrimaryFactor(Barbie = Barbie),
    "setting Treat as the primary effector of sample conditions."
  )
  expect_message(
    extractTargetsAndPrimaryFactor(
      Barbie = Barbie,
      targets = data.frame(Treat = Treat, Time = Time, Block = Block)
    ),
    "setting Treat as the primary effector of sample conditions."
  )
  expect_message(
    extractTargetsAndPrimaryFactor(Barbie = Barbie, sampleGroups = "Treat"),
    "setting Treat as the primary effector of sample conditions."
  )
  expect_equal(
    extractTargetsAndPrimaryFactor(
      Barbie = Barbie,
      targets = data.frame(Treat = Treat, Time = Time, Block = Block),
      sampleGroups = "Treat"
    ),
    list(
      mytargets = data.frame(Treat = Treat, Time = Time, Block = Block),
      pointer = 1
    )
  )
  expect_error(
    extractTargetsAndPrimaryFactor(Barbie = Barbie, sampleGroups = "treat")
  )
  expect_error(
    extractTargetsAndPrimaryFactor(Barbie = Barbie, sampleGroups = seq_len(6))
  )
})
