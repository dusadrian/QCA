local_edition(3)
noflevels <- c(2, 2, 2)

test_that("tests for allExpressions() have the same output", {
  expect_snapshot(allExpressions(noflevels))
  expect_snapshot(allExpressions(noflevels, arrange = TRUE))
  expect_snapshot(allExpressions(noflevels, raw = TRUE))
})
