as_edge_table_undirected <- function(adj.list, weight.list) {
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
    if (!any(keep)) next

    u <- c(u, rep.int(i, sum(keep)))
    v <- c(v, nbrs[keep])
    w <- c(w, wts[keep])
  }

  if (!length(u)) {
    return(data.frame(u = integer(0), v = integer(0), w = numeric(0)))
  }

  ord <- order(u, v)
  data.frame(u = u[ord], v = v[ord], w = w[ord])
}

make_complete_precomputed_graph <- function(n) {
  adj.list <- vector("list", n)
  weight.list <- vector("list", n)

  for (i in seq_len(n)) {
    nbrs <- setdiff(seq_len(n), i)
    adj.list[[i]] <- as.integer(nbrs)
    weight.list[[i]] <- as.numeric(abs(nbrs - i) + 1)
  }

  list(adj.list = adj.list, weight.list = weight.list)
}

test_that("precomputed graph mode rejects density reference and accepts counting measure", {
  set.seed(2614)
  n <- 24L
  graph <- make_complete_precomputed_graph(n)
  X <- matrix(rnorm(n * 3), nrow = n, ncol = 3)
  y <- sin(seq_len(n) / 3)

  expect_error(
    fit.rdgraph.regression(
      X,
      y,
      k = 5L,
      adj.list = graph$adj.list,
      weight.list = graph$weight.list,
      use.counting.measure = FALSE,
      max.iterations = 1L,
      n.eigenpairs = 10L,
      pca.dim = NULL,
      apply.geometric.pruning = FALSE,
      max.ratio.threshold = 0,
      threshold.percentile = 0,
      response.penalty.exp = 0,
      dense.fallback = "auto",
      verbose.level = 0L
    ),
    "precomputed adj.list/weight.list graphs currently require use.counting.measure = TRUE"
  )

  fit <- suppressWarnings(
    fit.rdgraph.regression(
      X,
      y,
      k = 5L,
      adj.list = graph$adj.list,
      weight.list = graph$weight.list,
      use.counting.measure = TRUE,
      max.iterations = 1L,
      n.eigenpairs = 10L,
      pca.dim = NULL,
      apply.geometric.pruning = FALSE,
      max.ratio.threshold = 0,
      threshold.percentile = 0,
      response.penalty.exp = 0,
      dense.fallback = "auto",
      verbose.level = 0L
    )
  )

  expect_identical(fit$parameters$graph.source, "precomputed")
  expect_true(all(is.finite(fit$fitted.values)))
})

test_that("dense.fallback='never' disables forced dense spectral path", {
  set.seed(2615)
  n <- 80L
  X <- matrix(rnorm(n * 4), nrow = n, ncol = 4)
  y <- sin(seq_len(n) / 5)

  out <- capture.output({
    fit <- suppressWarnings(
      fit.rdgraph.regression(
        X,
        y,
        k = 9L,
        max.iterations = 1L,
        n.eigenpairs = 20L,
        pca.dim = NULL,
        apply.geometric.pruning = FALSE,
        max.ratio.threshold = 0,
        threshold.percentile = 0,
        response.penalty.exp = 0,
        dense.fallback = "never",
        verbose.level = 2L
      )
    )
  })

  expect_false(any(grepl("Using dense solver", out, fixed = TRUE)))
  expect_true(any(grepl("Using sparse iterative solver", out, fixed = TRUE)))
  expect_identical(fit$parameters$dense.fallback, "never")
  expect_true(all(is.finite(fit$fitted.values)))
})

test_that("fit accepts precomputed graph and matches standard fit-path result", {
  set.seed(2611)
  X <- matrix(rnorm(120 * 7), nrow = 120, ncol = 7)
  y <- rnorm(nrow(X))
  k <- 9L

  iknn <- dgraphs::create.single.iknn.graph(
    X,
    k = k,
    max.path.edge.ratio.deviation.thld = 0,
    threshold.percentile = 0,
    compute.full = FALSE,
    verbose = FALSE
  )

  fit.standard <- fit.rdgraph.regression(
    X,
    y,
    k = k,
    max.iterations = 1L,
    n.eigenpairs = 100L,
    pca.dim = NULL,
    apply.geometric.pruning = FALSE,
    max.ratio.threshold = 0,
    threshold.percentile = 0,
    dense.fallback = "never",
    verbose.level = 0L
  )

  fit.injected <- fit.rdgraph.regression(
    X,
    y,
    k = k,
    adj.list = iknn$pruned_adj_list,
    weight.list = iknn$pruned_weight_list,
    max.iterations = 1L,
    n.eigenpairs = 100L,
    pca.dim = NULL,
    apply.geometric.pruning = FALSE,
    max.ratio.threshold = 0,
    threshold.percentile = 0,
    dense.fallback = "never",
    verbose.level = 0L
  )

  expect_identical(fit.injected$parameters$graph.source, "precomputed")
  expect_identical(fit.standard$graph$n.edges, fit.injected$graph$n.edges)

  edges.standard <- as_edge_table_undirected(
    fit.standard$graph$adj.list,
    fit.standard$graph$edge.length.list
  )
  edges.injected <- as_edge_table_undirected(
    fit.injected$graph$adj.list,
    fit.injected$graph$edge.length.list
  )

  expect_equal(edges.standard[, c("u", "v")], edges.injected[, c("u", "v")])
  expect_equal(edges.standard$w, edges.injected$w, tolerance = 1e-12)
  expect_identical(length(fit.standard$fitted.values), length(fit.injected$fitted.values))
  expect_true(all(is.finite(fit.injected$fitted.values)))
  expect_true(all(is.finite(fit.injected$residuals)))
})

test_that("fit with same injected graph reproduces internal fitted estimates", {
  set.seed(2613)
  X <- matrix(rnorm(140 * 8), nrow = 140, ncol = 8)
  y <- rnorm(nrow(X))
  k <- 11L

  fit.internal <- fit.rdgraph.regression(
    X,
    y,
    k = k,
    max.iterations = 1L,
    n.eigenpairs = 100L,
    pca.dim = NULL,
    apply.geometric.pruning = FALSE,
    max.ratio.threshold = 0,
    threshold.percentile = 0,
    dense.fallback = "never",
    verbose.level = 0L
  )

  fit.precomputed <- fit.rdgraph.regression(
    X,
    y,
    k = k,
    adj.list = fit.internal$graph$adj.list,
    weight.list = fit.internal$graph$edge.length.list,
    max.iterations = 1L,
    n.eigenpairs = 100L,
    pca.dim = NULL,
    apply.geometric.pruning = FALSE,
    max.ratio.threshold = 0,
    threshold.percentile = 0,
    dense.fallback = "never",
    verbose.level = 0L
  )

  expect_identical(fit.precomputed$parameters$graph.source, "precomputed")
  expect_equal(fit.internal$fitted.values, fit.precomputed$fitted.values, tolerance = 1e-12)
  expect_equal(fit.internal$residuals, fit.precomputed$residuals, tolerance = 1e-12)
  expect_equal(
    fit.internal$graph$edge.densities,
    fit.precomputed$graph$edge.densities,
    tolerance = 1e-12
  )
})

test_that("fit rejects malformed precomputed graph inputs", {
  set.seed(2612)
  X <- matrix(rnorm(80 * 6), nrow = 80, ncol = 6)
  y <- rnorm(nrow(X))
  k <- 8L

  iknn <- dgraphs::create.single.iknn.graph(
    X,
    k = k,
    max.path.edge.ratio.deviation.thld = 0,
    threshold.percentile = 0,
    compute.full = FALSE,
    verbose = FALSE
  )

  expect_error(
    fit.rdgraph.regression(
      X,
      y,
      k = k,
      adj.list = iknn$pruned_adj_list,
      weight.list = NULL,
      max.iterations = 1L,
      n.eigenpairs = 40L,
      pca.dim = NULL,
      verbose.level = 0L
    ),
    "provided together"
  )

  adj.bad <- iknn$pruned_adj_list
  w.bad <- iknn$pruned_weight_list
  i <- which(lengths(adj.bad) > 0L)[1L]
  j <- adj.bad[[i]][1L]
  rev.idx <- match(i, adj.bad[[j]])
  adj.bad[[j]] <- adj.bad[[j]][-rev.idx]
  w.bad[[j]] <- w.bad[[j]][-rev.idx]

  expect_error(
    fit.rdgraph.regression(
      X,
      y,
      k = k,
      adj.list = adj.bad,
      weight.list = w.bad,
      max.iterations = 1L,
      n.eigenpairs = 40L,
      pca.dim = NULL,
      apply.geometric.pruning = FALSE,
      max.ratio.threshold = 0,
      threshold.percentile = 0,
      dense.fallback = "never",
      verbose.level = 0L
    ),
    "undirected|reciprocal"
  )
})
