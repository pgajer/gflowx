make_bridge_graph_with_empty_local_overlap <- function() {
  n <- 12L
  adj.list <- vector("list", n)
  weight.list <- vector("list", n)

  add.edge <- function(i, j, w) {
    adj.list[[i]] <<- c(adj.list[[i]], j)
    weight.list[[i]] <<- c(weight.list[[i]], w)
    adj.list[[j]] <<- c(adj.list[[j]], i)
    weight.list[[j]] <<- c(weight.list[[j]], w)
  }

  for (block in list(1:6, 7:12)) {
    edge.pairs <- utils::combn(block, 2L)
    for (col.idx in seq_len(ncol(edge.pairs))) {
      i <- edge.pairs[1L, col.idx]
      j <- edge.pairs[2L, col.idx]
      add.edge(i, j, 1 + abs(i - j) / 100)
    }
  }

  # The bridge keeps the graph connected, but it is too long to enter either
  # endpoint's local neighborhood when k = 2 (C++ neighborhood size = k + 1).
  # Its neighborhood intersection is therefore empty, exercising the iterative
  # edge-density clamp added for density-reference fits.
  add.edge(6L, 7L, 100)

  list(adj.list = adj.list, weight.list = weight.list)
}

test_that("density-reference fits clamp empty-overlap iterative edge densities", {
  graph <- make_bridge_graph_with_empty_local_overlap()
  n <- length(graph$adj.list)
  X <- cbind(seq_len(n), rep(c(0, 10), each = n / 2))
  y <- sin(seq_len(n))
  adj.list.0based <- lapply(graph$adj.list, function(v) as.integer(v - 1L))
  weight.list.cpp <- lapply(graph$weight.list, as.double)

  expect_warning(
    fit <- .Call(
      get("S_fit_rdgraph_regression", envir = asNamespace("gflowx")),
      X,
      as.double(y),
      NULL,
      as.integer(3L),
      adj.list.0based,
      weight.list.cpp,
      FALSE,
      FALSE,
      0.95,
      1000L,
      1L,
      FALSE,
      0,
      0,
      0,
      0,
      0.5,
      0.1,
      0L,
      1.25,
      10L,
      "heat_kernel",
      1e-4,
      1e-4,
      1L,
      0,
      0.5,
      0,
      1.5,
      1e-10,
      FALSE,
      10,
      1000,
      1e12,
      FALSE,
      0.95,
      30L,
      NULL,
      0L,
      0L,
      0L,
      0L,
      PACKAGE = "gflowx"
    ),
    regexp = "update_edge_densities_from_vertices: clamped"
  )

  expect_type(fit, "list")
  expect_true(all(is.finite(fit$graph$edge.densities)))
  expect_gt(min(fit$graph$edge.densities), 0)
})
