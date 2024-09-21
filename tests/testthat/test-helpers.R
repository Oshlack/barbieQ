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
