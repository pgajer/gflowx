as_edge_table <- function(adj.list, weight.list) {
  n <- length(adj.list)
  from <- rep.int(seq_len(n), lengths(adj.list))
  to <- as.integer(unlist(adj.list, use.names = FALSE))
  wt <- as.numeric(unlist(weight.list, use.names = FALSE))

  keep <- is.finite(to) & to >= 1L & to <= n & is.finite(wt)
  from <- from[keep]
  to <- to[keep]
  wt <- wt[keep]

  i <- pmin(from, to)
  j <- pmax(from, to)
  keep2 <- i != j
  i <- i[keep2]
  j <- j[keep2]
  wt <- wt[keep2]

  key <- paste(i, j, sep = "_")
  ord <- order(key, wt)
  key <- key[ord]
  i <- i[ord]
  j <- j[ord]
  wt <- wt[ord]

  first <- !duplicated(key)
  data.frame(
    i = i[first],
    j = j[first],
    w = wt[first],
    stringsAsFactors = FALSE
  )
}

is_truthy_env <- function(value) {
  tolower(trimws(value)) %in% c("1", "true", "t", "yes", "y", "on")
}

skip_if_slow_triangle_policy_tests_disabled <- function() {
  if (!is_truthy_env(Sys.getenv("GFLOW_RUN_SLOW_TESTS", unset = ""))) {
    skip("Slow triangle-policy tests disabled. Set GFLOW_RUN_SLOW_TESTS=true to enable.")
  }
}

collect_warnings <- function(expr) {
  out <- character(0)
  val <- withCallingHandlers(
    expr,
    warning = function(w) {
      out <<- c(out, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(value = val, warnings = out)
}

test_that("triangle.policy argument is validated", {
  set.seed(2101)
  X <- matrix(rnorm(100 * 5), nrow = 100, ncol = 5)
  y <- rnorm(100)

  expect_error(
    fit.rdgraph.regression(
      X = X,
      y = y,
      k = 8L,
      max.iterations = 1L,
      n.eigenpairs = 60L,
      pca.dim = NULL,
      triangle.policy = "invalid-mode",
      verbose.level = 0L
    ),
    "should be one of"
  )
})

test_that("triangle.policy auto and never agree when response.penalty.exp is zero", {
  set.seed(2102)
  X <- matrix(rnorm(140 * 6), nrow = 140, ncol = 6)
  y <- rnorm(140)

  fit.auto <- fit.rdgraph.regression(
    X = X,
    y = y,
    k = 9L,
    max.iterations = 1L,
    n.eigenpairs = 100L,
    pca.dim = NULL,
    apply.geometric.pruning = FALSE,
    response.penalty.exp = 0,
    triangle.policy = "auto",
    dense.fallback = "never",
    verbose.level = 0L
  )

  fit.never <- fit.rdgraph.regression(
    X = X,
    y = y,
    k = 9L,
    max.iterations = 1L,
    n.eigenpairs = 100L,
    pca.dim = NULL,
    apply.geometric.pruning = FALSE,
    response.penalty.exp = 0,
    triangle.policy = "never",
    dense.fallback = "never",
    verbose.level = 0L
  )

  edges.auto <- as_edge_table(fit.auto$graph$adj.list, fit.auto$graph$edge.length.list)
  edges.never <- as_edge_table(fit.never$graph$adj.list, fit.never$graph$edge.length.list)

  expect_equal(edges.auto, edges.never, tolerance = 1e-12)
  expect_equal(fit.auto$fitted.values, fit.never$fitted.values, tolerance = 1e-12)

  # Cleanup check: these modes are now included in fit$parameters.
  expect_identical(fit.never$parameters$dense.fallback, "never")
  expect_identical(fit.never$parameters$triangle.policy, "never")
})

test_that("triangle.policy never disables triangles needed for response modulation", {
  set.seed(2103)
  X <- matrix(rnorm(130 * 6), nrow = 130, ncol = 6)
  y <- rnorm(130)

  never.res <- collect_warnings(
    fit.rdgraph.regression(
      X = X,
      y = y,
      k = 10L,
      max.iterations = 1L,
      n.eigenpairs = 100L,
      pca.dim = NULL,
      apply.geometric.pruning = FALSE,
      response.penalty.exp = 0.8,
      triangle.policy = "never",
      dense.fallback = "never",
      verbose.level = 0L
    )
  )

  always.res <- collect_warnings(
    fit.rdgraph.regression(
      X = X,
      y = y,
      k = 10L,
      max.iterations = 1L,
      n.eigenpairs = 100L,
      pca.dim = NULL,
      apply.geometric.pruning = FALSE,
      response.penalty.exp = 0.8,
      triangle.policy = "always",
      dense.fallback = "never",
      verbose.level = 0L
    )
  )

  expect_true(any(grepl("requires triangles for off-diagonal modulation", never.res$warnings, fixed = TRUE)))
  expect_false(any(grepl("requires triangles for off-diagonal modulation", always.res$warnings, fixed = TRUE)))
})

test_that("triangle.policy auto skips triangles on large-edge graph (slow)", {
  skip_on_cran()
  skip_if_slow_triangle_policy_tests_disabled()

  set.seed(2104)
  n <- 1900L
  X <- matrix(rnorm(n * 6), nrow = n, ncol = 6)
  y <- rnorm(n)

  auto.res <- collect_warnings(
    fit.rdgraph.regression(
      X = X,
      y = y,
      k = 20L,
      max.iterations = 1L,
      n.eigenpairs = 30L,
      pca.dim = NULL,
      apply.geometric.pruning = FALSE,
      response.penalty.exp = 0.8,
      triangle.policy = "auto",
      dense.fallback = "never",
      verbose.level = 0L
    )
  )

  expect_gt(auto.res$value$graph$n.edges, 200000L)
  expect_true(any(grepl("n_edges=.*exceeds 200000", auto.res$warnings)))
  expect_true(any(grepl("requires triangles for off-diagonal modulation", auto.res$warnings, fixed = TRUE)))
})
