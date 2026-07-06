local_edition(3)
PI1 <- c("A", "B", "~C", "D")
CO1 <- c("A~BCD", "A~BC~D", "~AB~C~D", "~ABCD")
chart1 <- makeChart(PI1, CO1, snames = c(A, B, C, D))
PI2 <- "A, B, ~C"
CO2 <- "ABC, A~B~C, A~BC, ~AB~C"
chart2 <- makeChart(PI2, CO2, snames = "A,B,C")
PI3 <- c("AB", "BC", "A~C", "~AC", "~A~B~D", "~B~C~D")
CO3 <- c("ABCD", "ABC~D", "AB~CD", "AB~C~D", "A~B~CD", "A~B~C~D", "~ABCD", "~ABC~D", "~A~BCD", "~A~BC~D", "~A~B~C~D")
chart3 <- makeChart(PI3, CO3, snames = c(A, B, C, D))
PI4 <- c("EF", "~GH", "IJ")
CO4 <- c("~EF*~GH*IJ", "EF*GH*~IJ", "~EF*GH*IJ", "EF*~GH*~IJ")
chart4 <- makeChart(PI4, CO4)


test_that("tests for solveChart() have the same output", {
  expect_snapshot(chart1)
  expect_snapshot(solveChart(chart1))
  expect_snapshot(solveChart(chart1, type = "lagrangian"))
  expect_snapshot(solveChart(chart1, all.sol = TRUE))
  expect_snapshot(chart2)
  expect_snapshot(solveChart(chart2))
  expect_snapshot(solveChart(chart3))
  expect_snapshot(solveChart(chart4))
})
