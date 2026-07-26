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

test_that("fit.rdgraph.regression honors explicit geometric pruning off switch", {
  set.seed(1101)
  X <- matrix(rnorm(90 * 6), nrow = 90, ncol = 6)
  y <- rnorm(90)

  fit.off.a <- fit.rdgraph.regression(
    X = X,
    y = y,
    k = 8L,
    max.iterations = 1L,
    n.eigenpairs = 20L,
    pca.dim = NULL,
    apply.geometric.pruning = FALSE,
    max.ratio.threshold = 0.10,
    verbose.level = 0L
  )

  fit.off.b <- fit.rdgraph.regression(
    X = X,
    y = y,
    k = 8L,
    max.iterations = 1L,
    n.eigenpairs = 20L,
    pca.dim = NULL,
    apply.geometric.pruning = FALSE,
    max.ratio.threshold = 0.00,
    verbose.level = 0L
  )

  edges.a <- as_edge_table(fit.off.a$graph$adj.list, fit.off.a$graph$edge.length.list)
  edges.b <- as_edge_table(fit.off.b$graph$adj.list, fit.off.b$graph$edge.length.list)

  expect_equal(edges.a, edges.b, tolerance = 1e-12)
  expect_equal(fit.off.a$fitted.values, fit.off.b$fitted.values, tolerance = 1e-12)
})

test_that("geometric pruning switch can only reduce (or preserve) edge count", {
  set.seed(1102)
  theta <- seq(0, 2 * pi, length.out = 96L)
  X <- cbind(cos(theta), sin(theta)) + matrix(rnorm(96 * 2, sd = 0.04), ncol = 2)
  y <- sin(theta) + rnorm(96, sd = 0.03)

  fit.on <- fit.rdgraph.regression(
    X = X,
    y = y,
    k = 10L,
    max.iterations = 1L,
    n.eigenpairs = 20L,
    pca.dim = NULL,
    apply.geometric.pruning = TRUE,
    max.ratio.threshold = 0.10,
    verbose.level = 0L
  )

  fit.off <- fit.rdgraph.regression(
    X = X,
    y = y,
    k = 10L,
    max.iterations = 1L,
    n.eigenpairs = 20L,
    pca.dim = NULL,
    apply.geometric.pruning = FALSE,
    max.ratio.threshold = 0.10,
    verbose.level = 0L
  )

  n.edges.on <- nrow(as_edge_table(fit.on$graph$adj.list, fit.on$graph$edge.length.list))
  n.edges.off <- nrow(as_edge_table(fit.off$graph$adj.list, fit.off$graph$edge.length.list))

  expect_lte(n.edges.on, n.edges.off)
})

test_that("deviation threshold zero behaves like geometric pruning off", {
  set.seed(1104)
  X <- matrix(rnorm(96 * 5), nrow = 96, ncol = 5)
  y <- rnorm(96)

  fit.zero <- fit.rdgraph.regression(
    X = X,
    y = y,
    k = 10L,
    max.iterations = 1L,
    n.eigenpairs = 20L,
    pca.dim = NULL,
    apply.geometric.pruning = TRUE,
    max.ratio.threshold = 0.00,
    verbose.level = 0L
  )

  fit.off <- fit.rdgraph.regression(
    X = X,
    y = y,
    k = 10L,
    max.iterations = 1L,
    n.eigenpairs = 20L,
    pca.dim = NULL,
    apply.geometric.pruning = FALSE,
    max.ratio.threshold = 0.10,
    verbose.level = 0L
  )

  edges.zero <- as_edge_table(fit.zero$graph$adj.list, fit.zero$graph$edge.length.list)
  edges.off <- as_edge_table(fit.off$graph$adj.list, fit.off$graph$edge.length.list)

  expect_equal(edges.zero, edges.off, tolerance = 1e-12)
  expect_equal(fit.zero$fitted.values, fit.off$fitted.values, tolerance = 1e-12)
})

test_that("apply.geometric.pruning argument is validated", {
  set.seed(1103)
  X <- matrix(rnorm(80 * 5), nrow = 80, ncol = 5)
  y <- rnorm(80)

  expect_error(
    fit.rdgraph.regression(
      X = X,
      y = y,
      k = 8L,
      max.iterations = 1L,
      n.eigenpairs = 20L,
      pca.dim = NULL,
      apply.geometric.pruning = NA,
      verbose.level = 0L
    ),
    "apply.geometric.pruning must be TRUE or FALSE"
  )
})
