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


test_that("serial and parallel PI generation retain the same chart and models", {
  previous <- Sys.getenv("QCA_NUM_THREADS", unset = NA_character_)
  on.exit({
    if (is.na(previous)) {
      Sys.unsetenv("QCA_NUM_THREADS")
    } else {
      Sys.setenv(QCA_NUM_THREADS = previous)
    }
  }, add = TRUE)

  data(d.represent, package = "QCA")

  Sys.setenv(QCA_NUM_THREADS = "1")
  serial <- minimize(
    d.represent,
    outcome = "WNP",
    include = "?",
    details = TRUE
  )

  Sys.setenv(QCA_NUM_THREADS = "4")
  parallel <- minimize(
    d.represent,
    outcome = "WNP",
    include = "?",
    details = TRUE
  )

  expect_equal(parallel$PIchart, serial$PIchart)
  expect_equal(parallel$solution, serial$solution)
})


test_that("CCubes transfers a validated cover between complexity levels", {
  data(d.represent, package = "QCA")
  .Call("C_resetScpProfile", PACKAGE = "QCA")

  minimize(
    d.represent,
    outcome = "WNP",
    include = "?",
    details = TRUE,
    first.min = TRUE
  )
  profile <- .Call("C_getScpProfile", PACKAGE = "QCA")

  expect_true(profile$incumbent_requested)
  expect_true(profile$incumbent_accepted)
})
