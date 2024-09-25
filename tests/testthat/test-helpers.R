test_that("converting first non-numeric column into rownames works", {
  mat <- data.frame(id = letters[1:5], matrix(1:25,5,5))
  expectMat <- as.matrix(data.frame(matrix(1:25,5,5), row.names=letters[1:5]))
  expect_equal(returnNumMat(mat), expectMat)
})

test_that("check Barbie dimensions works", {
  Block <- c(1,1,2,3,3,4,1,1,2,3,3,4)
  Treat <- factor(rep(1:2, each=6))
  Time <- rep(rep(1:2, each=3), 2)
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes*nsamples), nbarcodes, nsamples) %>% abs()
  rownames(count) <- paste0("Barcode", 1:nbarcodes)
  Barbie <- Barbie::createBarbie(count, data.frame(Treat=Treat, Time=Time))
  expect_true(checkBarbieDimensions(Barbie))
  Barbie2 <- Barbie
  Barbie2$metadata <- Barbie2$metadata[,-1,drop=FALSE]
  expect_true(checkBarbieDimensions(Barbie2))
  Barbie3 <- Barbie
  Barbie3$metadata <- Barbie3$metadata[-1,,drop=FALSE]
  expect_error(checkBarbieDimensions(Barbie3))
})

test_that("extarcting targets and groupBy works", {
  Block <- c(1,1,2,3,3,4,1,1,2,3,3,4)
  Treat <- factor(rep(1:2, each=6))
  Time <- rep(rep(1:2, each=3), 2)
  nbarcodes <- 50
  nsamples <- 12
  count <- matrix(rnorm(nbarcodes*nsamples), nbarcodes, nsamples) %>% abs()
  rownames(count) <- paste0("Barcode", 1:nbarcodes)
  Barbie <- Barbie::createBarbie(count, data.frame(Treat=Treat, Time=Time))

  expect_message(
    extarctTargetsAndPrimaryFactor(Barbie = Barbie),
    "no properly specified 'groupBy'. setting samples by homogenenous group.")
  expect_message(
    extarctTargetsAndPrimaryFactor(Barbie = Barbie,
                                   targets = data.frame(Treat=Treat, Time=Time, Block = Block)),
    "no properly specified 'groupBy'. setting samples by homogenenous group."
  )
  expect_message(
    extarctTargetsAndPrimaryFactor(Barbie = Barbie, groupBy = "Treat"),
    "setting Treat as the primary effector of sample conditions."
  )
  expect_equal(
    extarctTargetsAndPrimaryFactor(Barbie = Barbie,
                                   targets = data.frame(Treat=Treat, Time=Time, Block = Block),
                                   groupBy = "Treat"),
    list(mytargets = data.frame(Treat=Treat, Time=Time, Block = Block),
         pointer = 1)
    )
  expect_error(extarctTargetsAndPrimaryFactor(Barbie = Barbie, groupBy = "treat"))
  expect_error(extarctTargetsAndPrimaryFactor(Barbie = Barbie, groupBy = c(1:6)))
})
