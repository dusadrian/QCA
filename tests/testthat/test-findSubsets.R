local_edition(3)
noflevels <- c(2, 2, 2)

test_that("tests for findSubsets() have the same output", {
  expect_snapshot(findSubsets(input = 2, noflevels = noflevels + 1))
  expect_snapshot(findSubsets(input = 2, noflevels = noflevels + 1, stop = 20))
  expect_snapshot(findSubsets(input = c(8, 79), noflevels = rep(3, 4)))
})
