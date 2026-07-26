test_that("subject.neighborhood.stats computes expected diagnostics for conductance weights", {
  fit <- list(
    fitted.values = c(0.1, 0.2, 0.3, 0.4),
    graph = list(
      adj.list = list(
        c(2L, 3L),
        c(1L, 3L),
        c(1L, 2L, 4L),
        c(3L)
      ),
      edge.list = rbind(
        c(1L, 2L),
        c(1L, 3L),
        c(2L, 3L),
        c(3L, 4L)
      ),
      edge.densities = c(2, 4, 1, 0.5),
      vertex.densities = c(1, 4, 1, 2)
    )
  )

  subj.id <- c("A", "A", "B", "C")

  stats <- subject.neighborhood.stats(
    fit,
    subj.id,
    weight.type = "conductance",
    include.self.in.R = TRUE
  )

  expect_equal(stats$R, c(1.5, 1.5, 4 / 3, 1))
  expect_equal(stats$p.max, c(2 / 3, 2 / 3, 0.6153846154, 1), tolerance = 1e-10)
  expect_equal(stats$s.eff, c(1.8, 1.8, 1.8988764045, 1), tolerance = 1e-10)
})

test_that("subject.neighborhood.stats supports mass.sym weights and optimal.fit wrapper", {
  fit <- list(
    fitted.values = c(0.1, 0.2, 0.3, 0.4),
    graph = list(
      adj.list = list(
        c(2L, 3L),
        c(1L, 3L),
        c(1L, 2L, 4L),
        c(3L)
      ),
      edge.list = rbind(
        c(1L, 2L),
        c(1L, 3L),
        c(2L, 3L),
        c(3L, 4L)
      ),
      edge.densities = c(2, 4, 1, 0.5),
      vertex.densities = c(1, 4, 1, 2)
    )
  )

  wrapped <- list(optimal.fit = fit)
  subj.id <- c("A", "A", "B", "C")

  stats <- subject.neighborhood.stats(
    wrapped,
    subj.id,
    weight.type = "mass.sym",
    include.self.in.R = FALSE
  )

  expect_equal(stats$R, c(1, 1, 1.5, 1))
  expect_equal(stats$p.max, c(0.5, 2 / 3, 0.6534537935, 1), tolerance = 1e-8)
  expect_equal(stats$s.eff, c(2, 1.8, 1.8278312164, 1), tolerance = 1e-6)
})

test_that("subject.neighborhood.stats validates subject.id length", {
  fit <- list(
    fitted.values = c(0.1, 0.2),
    graph = list(
      adj.list = list(c(2L), c(1L)),
      edge.list = rbind(c(1L, 2L)),
      edge.densities = 1,
      vertex.densities = c(1, 1)
    )
  )

  expect_error(
    subject.neighborhood.stats(fit, subject.id = "s1"),
    "Length of subject.id"
  )
})
