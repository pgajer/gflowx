is_truthy_env <- function(value) {
  tolower(trimws(value)) %in% c("1", "true", "t", "yes", "y", "on")
}

skip_if_slow_runtime_tests_disabled <- function() {
  if (!is_truthy_env(Sys.getenv("GFLOW_RUN_SLOW_TESTS", unset = ""))) {
    skip("Slow runtime tests disabled. Set GFLOW_RUN_SLOW_TESTS=true to enable.")
  }
}

run_rdgraph_runtime_probe <- function(n,
                                      d,
                                      k,
                                      seed,
                                      max.ratio.threshold,
                                      threshold.percentile) {
  set.seed(seed)
  X <- matrix(rnorm(n * d), nrow = n, ncol = d)
  y <- rnorm(n)

  iknn.elapsed <- as.numeric(system.time({
    iknn <- dgraphs::create.single.iknn.graph(
      X,
      k = as.integer(k),
      max.path.edge.ratio.deviation.thld = max.ratio.threshold,
      threshold.percentile = threshold.percentile,
      compute.full = FALSE,
      verbose = FALSE
    )
  })["elapsed"])

  fit.elapsed <- as.numeric(system.time({
    fit <- fit.rdgraph.regression(
      X,
      y,
      k = as.integer(k),
      max.iterations = 1L,
      n.eigenpairs = min(90L, as.integer(max(40L, floor(n * 0.25)))),
      pca.dim = NULL,
      apply.geometric.pruning = max.ratio.threshold > 0,
      max.ratio.threshold = max.ratio.threshold,
      threshold.percentile = threshold.percentile,
      dense.fallback = "never",
      verbose.level = 0L
    )
  })["elapsed"])

  list(
    iknn.edges = as.integer(iknn$n_edges_in_pruned_graph),
    fit.edges = as.integer(fit$graph$n.edges),
    iknn.elapsed = iknn.elapsed,
    fit.elapsed = fit.elapsed
  )
}

test_that("small fixed dataset runtime and edge-count guardrails hold", {
  skip_on_cran()
  skip_if_slow_runtime_tests_disabled()

  res <- run_rdgraph_runtime_probe(
    n = 500L,
    d = 8L,
    k = 12L,
    seed = 7101L,
    max.ratio.threshold = 0.1,
    threshold.percentile = 0.2
  )

  expect_identical(res$iknn.edges, res$fit.edges)
  expect_true(res$iknn.edges >= 9000L)
  expect_true(res$iknn.edges <= 18000L)

  expect_true(res$iknn.elapsed <= 60, info = sprintf("iknn elapsed=%.3fs", res$iknn.elapsed))
  expect_true(res$fit.elapsed <= 120, info = sprintf("fit elapsed=%.3fs", res$fit.elapsed))
})

test_that("medium fixed dataset runtime and edge-count guardrails hold", {
  skip_on_cran()
  skip_if_slow_runtime_tests_disabled()

  res <- run_rdgraph_runtime_probe(
    n = 1000L,
    d = 8L,
    k = 16L,
    seed = 7102L,
    max.ratio.threshold = 0.1,
    threshold.percentile = 0.2
  )

  expect_identical(res$iknn.edges, res$fit.edges)
  expect_true(res$iknn.edges >= 25000L)
  expect_true(res$iknn.edges <= 65000L)

  expect_true(res$iknn.elapsed <= 120, info = sprintf("iknn elapsed=%.3fs", res$iknn.elapsed))
  expect_true(res$fit.elapsed <= 180, info = sprintf("fit elapsed=%.3fs", res$fit.elapsed))
})
