local_edition(3)
data(d.educate)
cc <- causalChain(d.educate)
data(d.women)
data(d.pban)
data(d.autonomy)
dat2 <- d.autonomy[15:30, c("AU","RE", "CN", "DE")]

test_that("tests for causalChain() have the same output", {
  expect_snapshot(cc$E$IC)
  expect_snapshot(causalChain(d.women))
  expect_snapshot(causalChain(d.pban, ordering = "C, F, T, V < PB", sol.cov = 0.95))
  expect_snapshot(causalChain(d.pban, ordering = "C, F, T, V < PB", pi.cons = 0.93, sol.cons = 0.95))
  expect_snapshot(causalChain(dat2, ordering = "AU", sol.cons = 0.9, pi.cons = 0.85, sol.cov = 0.85))
})
