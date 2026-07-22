local_edition(3)


make_adversarial_Ft <- function(t) {
  block <- 2L * t

  on <- rbind(
    c(rep(0L, block), 0L, 0L, 0L, 0L),
    c(rep(0L, block), 1L, 1L, 1L, 1L)
  )

  off <- t(apply(combn(block, t), 2L, function(T) {
    row <- rep(1L, block)
    row[T] <- 0L
    c(row, 0L, 1L, 0L, 1L)
  }))

  data <- as.data.frame(rbind(on, off))
  names(data) <- paste0("X", seq_len(ncol(data)))
  data$Y <- c(1L, 1L, rep(0L, nrow(off)))
  data
}


test_that("adaptive stopping reaches the deeper F_t minimum", {
  for (t in 2:4) {
    tt <- truthTable(
      make_adversarial_Ft(t),
      outcome = "Y",
      incl.cut = 1
    )

    first <- minimize(
      tt,
      include = "?",
      first.min = TRUE
    )

    expect_equal(
      length(first$solution[[1L]]),
      1L,
      info = paste0("F", t)
    )
  }
})


test_that("default minimization enumerates every F_t minimum model", {
  for (t in 2:4) {
    tt <- truthTable(
      make_adversarial_Ft(t),
      outcome = "Y",
      incl.cut = 1
    )

    result <- minimize(tt, include = "?")

    expect_true(
      all(lengths(result$solution) == 1L),
      info = paste0("F", t)
    )
    expect_equal(
      length(result$solution),
      as.integer(choose(2L * t, t + 1L)),
      info = paste0("F", t)
    )
  }
})


test_that("a warned F_t search cannot end uncertified at pi.depth", {
  for (t in 2:4) {
    tt <- truthTable(
      make_adversarial_Ft(t),
      outcome = "Y",
      incl.cut = 1
    )

    expect_error(
      do.call(
        minimize,
        list(
          input = tt,
          include = "?",
          first.min = TRUE,
          pi.depth = t
        )
      ),
      "before the adaptive stopping certificate",
      info = paste0("F", t)
    )
  }
})
