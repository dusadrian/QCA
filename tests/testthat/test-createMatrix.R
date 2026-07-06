local_edition(3)
noflevels <- c(2, 2, 2)

test_that("tests for createMatrix() have the same output", {
  expect_snapshot(createMatrix(noflevels))
  expect_snapshot(createMatrix(noflevels + 1))
  expect_snapshot(createMatrix(c(2, 3, 2)))
})
