rdgraph.weighted.circle.graph <- function(theta) {
  theta <- as.double(theta) %% (2 * pi)
  n <- length(theta)
  if (n < 3L) stop("Need at least three vertices for a circle graph.")

  ord <- order(theta)
  adj.list <- vector("list", n)
  weight.list <- vector("list", n)

  add.edge <- function(i, j, w) {
    adj.list[[i]] <<- c(adj.list[[i]], as.integer(j))
    weight.list[[i]] <<- c(weight.list[[i]], as.double(w))
    adj.list[[j]] <<- c(adj.list[[j]], as.integer(i))
    weight.list[[j]] <<- c(weight.list[[j]], as.double(w))
  }

  for (pos in seq_len(n)) {
    i <- ord[pos]
    j <- ord[if (pos == n) 1L else pos + 1L]
    raw.gap <- if (pos == n) {
      theta[ord[1L]] + 2 * pi - theta[i]
    } else {
      theta[j] - theta[i]
    }
    add.edge(i, j, raw.gap)
  }

  list(adj.list = adj.list, weight.list = weight.list)
}

rdgraph.weighted.edge.graph <- function(X, edges) {
  X <- as.matrix(X)
  n <- nrow(X)
  adj.list <- vector("list", n)
  weight.list <- vector("list", n)
  if (is.null(edges) || !length(edges)) {
    return(list(adj.list = adj.list, weight.list = weight.list))
  }
  edges <- as.matrix(edges)
  if (ncol(edges) != 2L) stop("edges must have two columns.")

  add.edge <- function(i, j) {
    i <- as.integer(i)
    j <- as.integer(j)
    w <- sqrt(sum((X[i, ] - X[j, ])^2))
    adj.list[[i]] <<- c(adj.list[[i]], j)
    weight.list[[i]] <<- c(weight.list[[i]], w)
    adj.list[[j]] <<- c(adj.list[[j]], i)
    weight.list[[j]] <<- c(weight.list[[j]], w)
  }

  for (e in seq_len(nrow(edges))) {
    add.edge(edges[e, 1L], edges[e, 2L])
  }
  list(adj.list = adj.list, weight.list = weight.list)
}

rdgraph.case.oracle.graph <- function(case, data) {
  if (!is.null(data$oracle.graph)) {
    return(list(
      model = data$oracle.graph$model %||% "planted_graph",
      adj.list = data$oracle.graph$adj.list,
      weight.list = data$oracle.graph$weight.list,
      command = data$oracle.graph$command %||% paste(
        "oracle.graph <- planted weighted graph returned by the case generator",
        "oracle.fit <- gflowx::fit.rdgraph.regression(",
        "  X = data$X, y = data$responses[[1]], k = 2L,",
        "  adj.list = oracle.graph$adj.list,",
        "  weight.list = oracle.graph$weight.list,",
        "  use.counting.measure = TRUE",
        ")",
        sep = "\n"
      )
    ))
  }

  if (identical(data$type, "line")) {
    g <- create.chain.graph(x = data$coord)
    return(list(
      model = "weighted_path_graph",
      adj.list = g$adj.list,
      weight.list = g$edge.lengths,
      command = paste(
        "path.graph <- gflow::create.chain.graph(x = data$coord)",
        "path.fit <- gflowx::fit.rdgraph.regression(",
        "  X = data$X, y = data$responses[[1]], k = 2L,",
        "  adj.list = path.graph$adj.list,",
        "  weight.list = path.graph$edge.lengths,",
        "  use.counting.measure = TRUE",
        ")",
        sep = "\n"
      )
    ))
  }

  if (identical(data$type, "circle")) {
    g <- rdgraph.weighted.circle.graph(data$theta)
    return(list(
      model = "weighted_circle_graph",
      adj.list = g$adj.list,
      weight.list = g$weight.list,
      command = paste(
        "circle.graph <- weighted circle graph from observed sorted angles",
        "circle.fit <- gflowx::fit.rdgraph.regression(",
        "  X = data$X, y = data$responses[[1]], k = 2L,",
        "  adj.list = circle.graph$adj.list,",
        "  weight.list = circle.graph$weight.list,",
        "  use.counting.measure = TRUE",
        ")",
        sep = "\n"
      )
    ))
  }

  NULL
}

rdgraph.edge.table <- function(adj.list, weight.list = NULL) {
  n <- length(adj.list)
  out <- vector("list", n)
  ptr <- 0L
  for (i in seq_len(n)) {
    nbrs <- as.integer(adj.list[[i]])
    if (!length(nbrs)) next
    wts <- if (is.null(weight.list)) rep(1, length(nbrs)) else as.double(weight.list[[i]])
    keep <- nbrs > i
    if (!any(keep)) next
    ptr <- ptr + 1L
    out[[ptr]] <- data.frame(
      from = rep.int(i, sum(keep)),
      to = nbrs[keep],
      weight = wts[keep]
    )
  }
  if (!ptr) {
    return(data.frame(from = integer(), to = integer(), weight = numeric()))
  }
  do.call(rbind, out[seq_len(ptr)])
}

rdgraph.graph.summary <- function(adj.list, weight.list = NULL, k = NA_integer_) {
  deg <- lengths(adj.list)
  cc <- graph.connected.components(adj.list)
  comp.sizes <- as.integer(tabulate(cc))
  comp.sizes <- comp.sizes[comp.sizes > 0L]
  comp.sizes.sorted <- sort(comp.sizes, decreasing = TRUE)
  n.edges <- sum(deg) / 2
  data.frame(
    k = as.integer(k),
    edge_count = as.integer(n.edges),
    component_count = length(comp.sizes),
    largest_component_size = max(comp.sizes),
    largest_component_fraction = max(comp.sizes) / length(adj.list),
    second_largest_component_size = if (length(comp.sizes.sorted) >= 2L) comp.sizes.sorted[[2L]] else 0L,
    mean_degree = mean(deg),
    min_degree = min(deg),
    max_degree = max(deg),
    stringsAsFactors = FALSE
  )
}

rdgraph.largest.component.indices <- function(adj.list) {
  cc <- graph.connected.components(adj.list)
  comp.sizes <- tabulate(cc)
  largest <- which.max(comp.sizes)
  which(cc == largest)
}

rdgraph.induced.graph <- function(adj.list, weight.list, indices) {
  indices <- as.integer(indices)
  map <- integer(length(adj.list))
  map[indices] <- seq_along(indices)

  adj.out <- vector("list", length(indices))
  weight.out <- vector("list", length(indices))
  for (pos in seq_along(indices)) {
    old.i <- indices[[pos]]
    nbrs <- as.integer(adj.list[[old.i]])
    wts <- as.double(weight.list[[old.i]])
    keep <- nbrs %in% indices
    adj.out[[pos]] <- as.integer(map[nbrs[keep]])
    weight.out[[pos]] <- wts[keep]
  }

  list(adj.list = adj.out, weight.list = weight.out)
}

rdgraph.degree.js <- function(adj1, adj2) {
  d1 <- lengths(adj1)
  d2 <- lengths(adj2)
  max.degree <- max(c(d1, d2))
  p <- tabulate(d1 + 1L, nbins = max.degree + 1L)
  q <- tabulate(d2 + 1L, nbins = max.degree + 1L)
  jensen.shannon.divergence(p, q)
}

rdgraph.gcv.value <- function(fit) {
  if (is.null(fit$gcv) || is.null(fit$gcv$gcv.optimal)) return(NA_real_)
  idx <- fit$optimal.iteration %||% which.min(fit$gcv$gcv.optimal)
  as.double(fit$gcv$gcv.optimal[[idx]])
}
