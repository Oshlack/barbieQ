test_that("converting first non-numeric column into rownames works", {
  mat <- data.frame(id = letters[1:5], matrix(1:25,5,5))
  expectMat <- as.matrix(data.frame(matrix(1:25,5,5), row.names=letters[1:5]))
  expect_equal(returnNumMat(mat), expectMat)
})
