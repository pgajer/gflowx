canonicalize_undirected_edges <- function(adj.list, weight.list) {
  expect_true(is.list(adj.list))
  expect_true(is.list(weight.list))
  expect_identical(length(adj.list), length(weight.list))

  u <- integer(0)
  v <- integer(0)
  w <- numeric(0)

  for (i in seq_along(adj.list)) {
    nbrs <- as.integer(adj.list[[i]])
    wts <- as.double(weight.list[[i]])

    expect_identical(length(nbrs), length(wts))

    keep <- nbrs > i
    if (!any(keep)) {
      next
    }

    n_keep <- sum(keep)
    u <- c(u, rep.int(i, n_keep))
    v <- c(v, nbrs[keep])
    w <- c(w, wts[keep])
  }

  if (!length(u)) {
    return(data.frame(u = integer(0), v = integer(0), w = numeric(0)))
  }

  ord <- order(u, v)
  data.frame(u = u[ord], v = v[ord], w = w[ord])
}

assert_cross_api_graph_invariance <- function(X,
                                              y,
                                              k,
                                              apply.geometric.pruning,
                                              max.ratio.threshold,
                                              threshold.percentile,
                                              path.edge.ratio.percentile = 0.5,
                                              case.label = "case") {
  fit <- fit.rdgraph.regression(
    X,
    y,
    k = as.integer(k),
    max.iterations = 1L,
    n.eigenpairs = 100L,
    pca.dim = NULL,
    apply.geometric.pruning = isTRUE(apply.geometric.pruning),
    max.ratio.threshold = max.ratio.threshold,
    path.edge.ratio.percentile = path.edge.ratio.percentile,
    threshold.percentile = threshold.percentile,
    dense.fallback = "never",
    verbose.level = 0L
  )

  iknn_geom_threshold <- if (isTRUE(apply.geometric.pruning)) max.ratio.threshold else 0
  iknn <- dgraphs::create.single.iknn.graph(
    X,
    k = as.integer(k),
    max.path.edge.ratio.deviation.thld = iknn_geom_threshold,
    path.edge.ratio.percentile = path.edge.ratio.percentile,
    threshold.percentile = threshold.percentile,
    compute.full = FALSE,
    verbose = FALSE
  )

  fit_edges <- canonicalize_undirected_edges(
    fit$graph$adj.list,
    fit$graph$edge.length.list
  )
  iknn_edges <- canonicalize_undirected_edges(
    iknn$pruned_adj_list,
    iknn$pruned_weight_list
  )

  expect_identical(nrow(fit_edges), fit$graph$n.edges, info = case.label)
  expect_equal(nrow(iknn_edges), iknn$n_edges_in_pruned_graph, info = case.label)
  expect_identical(nrow(fit_edges), nrow(iknn_edges), info = case.label)
  expect_equal(fit_edges[, c("u", "v")], iknn_edges[, c("u", "v")], info = case.label)
  expect_equal(fit_edges$w, iknn_edges$w, tolerance = 1e-12, info = case.label)
}

test_that("fit path and ikNN API produce identical unpruned winner graph", {
  set.seed(11)
  X <- matrix(rnorm(120 * 7), nrow = 120, ncol = 7)
  y <- rnorm(nrow(X))

  assert_cross_api_graph_invariance(
    X = X,
    y = y,
    k = 9L,
    apply.geometric.pruning = FALSE,
    max.ratio.threshold = 0,
    threshold.percentile = 0,
    case.label = "unpruned"
  )
})

test_that("fit path and ikNN API produce identical graph with geometric pruning only", {
  set.seed(12)
  X <- matrix(rnorm(120 * 7), nrow = 120, ncol = 7)
  y <- rnorm(nrow(X))

  assert_cross_api_graph_invariance(
    X = X,
    y = y,
    k = 9L,
    apply.geometric.pruning = TRUE,
    max.ratio.threshold = 0.1,
    threshold.percentile = 0,
    case.label = "geometric-only"
  )
})

test_that("fit path and ikNN API produce identical graph with quantile pruning only", {
  set.seed(13)
  X <- matrix(rnorm(120 * 7), nrow = 120, ncol = 7)
  y <- rnorm(nrow(X))

  assert_cross_api_graph_invariance(
    X = X,
    y = y,
    k = 9L,
    apply.geometric.pruning = FALSE,
    max.ratio.threshold = 0,
    threshold.percentile = 0.2,
    case.label = "quantile-only"
  )
})

test_that("fit path and ikNN API produce identical graph with geometric and quantile pruning", {
  set.seed(14)
  X <- matrix(rnorm(120 * 7), nrow = 120, ncol = 7)
  y <- rnorm(nrow(X))

  assert_cross_api_graph_invariance(
    X = X,
    y = y,
    k = 9L,
    apply.geometric.pruning = TRUE,
    max.ratio.threshold = 0.1,
    threshold.percentile = 0.2,
    case.label = "geometric+quantile"
  )
})

test_that("fit-path geometric threshold boundary (deviation=0) matches ikNN no-prune", {
  set.seed(15)
  X <- matrix(rnorm(120 * 7), nrow = 120, ncol = 7)
  y <- rnorm(nrow(X))

  assert_cross_api_graph_invariance(
    X = X,
    y = y,
    k = 9L,
    apply.geometric.pruning = TRUE,
    max.ratio.threshold = 0.0,
    threshold.percentile = 0.0,
    case.label = "geometric-boundary-zero"
  )
})
