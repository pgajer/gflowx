rdgraph.default.fit.args <- function(k, n.eigenpairs) {
  list(
    k = as.integer(k),
    max.iterations = 5L,
    n.eigenpairs = as.integer(n.eigenpairs),
    pca.dim = NULL,
    apply.geometric.pruning = FALSE,
    max.ratio.threshold = 0,
    threshold.percentile = 0,
    response.penalty.exp = 0,
    use.counting.measure = TRUE,
    filter.type = "heat_kernel",
    dense.fallback = "never",
    compute.extremality = FALSE,
    verbose.level = 0L
  )
}

rdgraph.gaussian.signal <- function(x, centers, heights, sds) {
  out <- numeric(length(x))
  for (i in seq_along(centers)) {
    out <- out + heights[[i]] * exp(-0.5 * ((x - centers[[i]]) / sds[[i]])^2)
  }
  out
}

rdgraph.circle.xy <- function(theta, noise = 0, radius = 1) {
  X <- cbind(x = radius * cos(theta), y = radius * sin(theta))
  if (noise > 0) {
    X <- X + matrix(stats::rnorm(length(theta) * 2L, sd = noise), ncol = 2L)
  }
  X
}

rdgraph.path.edges <- function(indices) {
  indices <- as.integer(indices)
  if (length(indices) < 2L) return(matrix(integer(0), ncol = 2L))
  cbind(indices[-length(indices)], indices[-1L])
}

rdgraph.knn.edges <- function(X, indices, k = 4L) {
  indices <- as.integer(indices)
  if (length(indices) < 2L) return(matrix(integer(0), ncol = 2L))
  X.sub <- as.matrix(X[indices, , drop = FALSE])
  d <- as.matrix(stats::dist(X.sub))
  edges <- list()
  ptr <- 0L
  for (i in seq_along(indices)) {
    ord <- order(d[i, ])
    ord <- ord[ord != i]
    nbrs <- head(ord, min(k, length(ord)))
    for (j in nbrs) {
      ptr <- ptr + 1L
      edges[[ptr]] <- sort(c(indices[[i]], indices[[j]]))
    }
  }
  edges <- unique(do.call(rbind, edges))
  matrix(as.integer(edges), ncol = 2L)
}

rdgraph.y.geometry <- function(noise = 0.025) {
  stem.t <- seq(0, 1, length.out = 50L)
  left.t <- seq(0.02, 1, length.out = 45L)
  right.t <- seq(0.02, 1, length.out = 45L)

  stem <- cbind(x = 0, y = -1 + stem.t)
  left <- cbind(x = -0.85 * left.t, y = 0.85 * left.t)
  right <- cbind(x = 0.85 * right.t, y = 0.85 * right.t)
  X <- rbind(stem, left, right)
  if (noise > 0) {
    X <- X + matrix(stats::rnorm(nrow(X) * 2L, sd = noise), ncol = 2L)
  }

  stem.idx <- seq_len(nrow(stem))
  left.idx <- length(stem.idx) + seq_len(nrow(left))
  right.idx <- length(stem.idx) + length(left.idx) + seq_len(nrow(right))
  edges <- rbind(
    rdgraph.path.edges(stem.idx),
    cbind(tail(stem.idx, 1L), head(left.idx, 1L)),
    cbind(tail(stem.idx, 1L), head(right.idx, 1L)),
    rdgraph.path.edges(left.idx),
    rdgraph.path.edges(right.idx)
  )

  list(
    X = X,
    edges = edges,
    branch = c(rep("stem", length(stem.idx)), rep("left", length(left.idx)), rep("right", length(right.idx))),
    branch.t = c(stem.t, left.t, right.t)
  )
}

rdgraph.case.1d.gaussian.mixture <- function() {
  list(
    id = "1d_gaussian_mixture_recovers_smooth_signal",
    title = "1D Gaussian Mixture: Smooth Signal Recovery",
    group = "signal_recovery",
    seed = 1001L,
    k.grid = 3:10,
    fit.args = rdgraph.default.fit.args(k = 6L, n.eigenpairs = 100L),
    data.command = paste(
      "set.seed(1001L)",
      "gm <- gflow::generate.1d.gaussian.mixture(",
      "  n.points = 120L,",
      "  x.knots = c(-4, 0.5, 5),",
      "  y.knots = c(1.4, 3.2, 1.7),",
      "  sd.knot = 0.9, x.offset = 1.5, add.noise = FALSE",
      ")",
      "X <- cbind(x = gm$x)",
      "y.true <- gm$y.true",
      "y <- y.true + rnorm(length(gm$x), sd = 0.35)",
      sep = "\n"
    ),
    fit.command = paste(
      "fit <- gflowx::fit.rdgraph.regression(",
      "  X = X, y = y, k = 6L, max.iterations = 5L,",
      "  n.eigenpairs = 100L, pca.dim = NULL,",
      "  apply.geometric.pruning = FALSE, max.ratio.threshold = 0,",
      "  threshold.percentile = 0, response.penalty.exp = 0,",
      "  use.counting.measure = TRUE, filter.type = \"heat_kernel\",",
      "  dense.fallback = \"never\", compute.extremality = FALSE,",
      "  verbose.level = 0L",
      ")",
      sep = "\n"
    ),
    generate = function() {
      set.seed(1001L)
      gm <- generate.1d.gaussian.mixture(
        n.points = 120L,
        x.knots = c(-4, 0.5, 5),
        y.knots = c(1.4, 3.2, 1.7),
        sd.knot = 0.9,
        x.offset = 1.5,
        add.noise = FALSE
      )
      y <- as.double(gm$y.true + stats::rnorm(length(gm$x), sd = 0.35))

      list(
        type = "line",
        X = cbind(x = as.double(gm$x)),
        coord = as.double(gm$x),
        y.true = as.double(gm$y.true),
        responses = list(noisy = y),
        metadata = list(
          noise_sd = 0.35,
          signal = "three-component 1D Gaussian mixture"
        )
      )
    }
  )
}

rdgraph.case.nonuniform.circle <- function() {
  list(
    id = "nonuniform_circle_recovers_periodic_signal",
    title = "Nonuniform Circle: Periodic Signal Under Density Imbalance",
    group = "geometry_density_stress",
    seed = 2001L,
    k.grid = 3:12,
    fit.args = modifyList(
      rdgraph.default.fit.args(k = 6L, n.eigenpairs = 120L),
      list(
        apply.geometric.pruning = TRUE,
        max.ratio.threshold = 0.1
      )
    ),
    data.command = paste(
      "set.seed(2001L)",
      "n <- 180L",
      "component <- sample(1:3, n, replace = TRUE, prob = c(0.55, 0.30, 0.15))",
      "theta <- numeric(n)",
      "theta[component == 1] <- rnorm(sum(component == 1), mean = 0.65, sd = 0.28)",
      "theta[component == 2] <- rnorm(sum(component == 2), mean = 3.05, sd = 0.45)",
      "theta[component == 3] <- runif(sum(component == 3), min = 0, max = 2 * pi)",
      "theta <- theta %% (2 * pi)",
      "X <- cbind(x = cos(theta), y = sin(theta)) + matrix(rnorm(2 * n, sd = 0.055), ncol = 2)",
      "y.true <- gflow::circular.synthetic.mixture.of.gaussians(",
      "  x = theta, x.knot = c(0.45, 2.35, 5.35),",
      "  y.knot = c(1.25, 0.85, 1.7), sd.knot = 0.32",
      ")",
      "y <- y.true + rnorm(n, sd = 0.12)",
      sep = "\n"
    ),
    fit.command = paste(
      "fit <- gflowx::fit.rdgraph.regression(",
      "  X = X, y = y, k = 6L, max.iterations = 5L,",
      "  n.eigenpairs = 120L, pca.dim = NULL,",
      "  apply.geometric.pruning = TRUE, max.ratio.threshold = 0.1,",
      "  threshold.percentile = 0, response.penalty.exp = 0,",
      "  use.counting.measure = TRUE, filter.type = \"heat_kernel\",",
      "  dense.fallback = \"never\", compute.extremality = FALSE,",
      "  verbose.level = 0L",
      ")",
      sep = "\n"
    ),
    generate = function() {
      set.seed(2001L)
      n <- 180L
      component <- sample(1:3, n, replace = TRUE, prob = c(0.55, 0.30, 0.15))
      theta <- numeric(n)
      theta[component == 1L] <- stats::rnorm(sum(component == 1L), mean = 0.65, sd = 0.28)
      theta[component == 2L] <- stats::rnorm(sum(component == 2L), mean = 3.05, sd = 0.45)
      theta[component == 3L] <- stats::runif(sum(component == 3L), min = 0, max = 2 * pi)
      theta <- as.double(theta %% (2 * pi))
      X <- rdgraph.circle.xy(theta, noise = 0.055)
      y.true <- circular.synthetic.mixture.of.gaussians(
        x = theta,
        x.knot = c(0.45, 2.35, 5.35),
        y.knot = c(1.25, 0.85, 1.7),
        sd.knot = 0.32
      )
      y <- as.double(y.true + stats::rnorm(n, sd = 0.12))

      list(
        type = "circle",
        X = X,
        coord = theta,
        theta = theta,
        y.true = as.double(y.true),
        responses = list(noisy = y),
        metadata = list(
          noise_sd = 0.12,
          geometry_noise = 0.055,
          signal = "periodic Gaussian mixture on a nonuniformly sampled circle"
        )
      )
    }
  )
}

rdgraph.case.circle.with.gap <- function() {
  list(
    id = "circle_with_gap_reports_connectivity_and_wraparound_limits",
    title = "Circle With Gap: Connectivity And Wraparound Limits",
    group = "geometry_topology_stress",
    seed = 2002L,
    k.grid = 3:12,
    fit.args = modifyList(
      rdgraph.default.fit.args(k = 9L, n.eigenpairs = 110L),
      list(
        apply.geometric.pruning = TRUE,
        max.ratio.threshold = 0.1
      )
    ),
    data.command = paste(
      "set.seed(2002L)",
      "n <- 150L",
      "theta <- runif(n, min = 0.70, max = 2 * pi - 0.70)",
      "theta <- sort(theta)",
      "X <- cbind(x = cos(theta), y = sin(theta)) + matrix(rnorm(2 * n, sd = 0.05), ncol = 2)",
      "y.true <- gflow::circular.synthetic.mixture.of.gaussians(",
      "  x = theta, x.knot = c(0.85, 2.70, 5.15),",
      "  y.knot = c(1.45, 0.8, 1.55), sd.knot = 0.28",
      ")",
      "y <- y.true + rnorm(n, sd = 0.12)",
      sep = "\n"
    ),
    fit.command = paste(
      "fit <- gflowx::fit.rdgraph.regression(",
      "  X = X, y = y, k = 9L, max.iterations = 5L,",
      "  n.eigenpairs = 110L, pca.dim = NULL,",
      "  apply.geometric.pruning = TRUE, max.ratio.threshold = 0.1,",
      "  threshold.percentile = 0, response.penalty.exp = 0,",
      "  use.counting.measure = TRUE, filter.type = \"heat_kernel\",",
      "  dense.fallback = \"never\", compute.extremality = FALSE,",
      "  verbose.level = 0L",
      ")",
      sep = "\n"
    ),
    generate = function() {
      set.seed(2002L)
      n <- 150L
      theta <- sort(stats::runif(n, min = 0.70, max = 2 * pi - 0.70))
      X <- rdgraph.circle.xy(theta, noise = 0.05)
      y.true <- circular.synthetic.mixture.of.gaussians(
        x = theta,
        x.knot = c(0.85, 2.70, 5.15),
        y.knot = c(1.45, 0.8, 1.55),
        sd.knot = 0.28
      )
      y <- as.double(y.true + stats::rnorm(n, sd = 0.12))

      list(
        type = "circle",
        X = X,
        coord = theta,
        theta = theta,
        y.true = as.double(y.true),
        responses = list(noisy = y),
        metadata = list(
          noise_sd = 0.12,
          geometry_noise = 0.05,
          angular_gap = "no samples in [0, 0.70) union (2*pi - 0.70, 2*pi]",
          signal = "periodic Gaussian mixture with a sampling gap around the wrap point"
        )
      )
    }
  )
}

rdgraph.case.1d.boundary.peaks <- function() {
  list(
    id = "1d_boundary_peaks_exposes_endpoint_bias",
    title = "1D Boundary Peaks: Endpoint Bias Stress",
    group = "boundary_stress",
    seed = 2003L,
    k.grid = 3:10,
    fit.args = rdgraph.default.fit.args(k = 6L, n.eigenpairs = 110L),
    data.command = paste(
      "set.seed(2003L)",
      "n <- 140L",
      "x <- sort(runif(n, min = 0, max = 1))",
      "y.true <- 0.25 + 1.4 * exp(-0.5 * ((x - 0.06) / 0.045)^2) +",
      "  0.75 * exp(-0.5 * ((x - 0.52) / 0.16)^2) +",
      "  1.2 * exp(-0.5 * ((x - 0.94) / 0.055)^2)",
      "X <- cbind(x = x)",
      "y <- y.true + rnorm(n, sd = 0.12)",
      sep = "\n"
    ),
    fit.command = paste(
      "fit <- gflowx::fit.rdgraph.regression(",
      "  X = X, y = y, k = 6L, max.iterations = 5L,",
      "  n.eigenpairs = 110L, pca.dim = NULL,",
      "  apply.geometric.pruning = FALSE, max.ratio.threshold = 0,",
      "  threshold.percentile = 0, response.penalty.exp = 0,",
      "  use.counting.measure = TRUE, filter.type = \"heat_kernel\",",
      "  dense.fallback = \"never\", compute.extremality = FALSE,",
      "  verbose.level = 0L",
      ")",
      sep = "\n"
    ),
    generate = function() {
      set.seed(2003L)
      n <- 140L
      x <- sort(stats::runif(n, min = 0, max = 1))
      y.true <- 0.25 + rdgraph.gaussian.signal(
        x,
        centers = c(0.06, 0.52, 0.94),
        heights = c(1.4, 0.75, 1.2),
        sds = c(0.045, 0.16, 0.055)
      )
      y <- as.double(y.true + stats::rnorm(n, sd = 0.12))

      list(
        type = "line",
        X = cbind(x = as.double(x)),
        coord = as.double(x),
        y.true = as.double(y.true),
        responses = list(noisy = y),
        metadata = list(
          noise_sd = 0.12,
          boundary_fraction = 0.10,
          signal = "1D smooth signal with narrow peaks near both endpoints"
        )
      )
    }
  )
}

rdgraph.case.pure.noise.response <- function() {
  list(
    id = "pure_noise_response_does_not_look_structured",
    title = "Pure Noise Response: Negative Control",
    group = "negative_control",
    seed = 2004L,
    k.grid = 3:10,
    fit.args = rdgraph.default.fit.args(k = 6L, n.eigenpairs = 100L),
    data.command = paste(
      "set.seed(2004L)",
      "n <- 120L",
      "x <- sort(runif(n, min = -3, max = 3))",
      "X <- cbind(x = x)",
      "y.true <- rep(0, n)",
      "y <- rnorm(n, sd = 1)",
      sep = "\n"
    ),
    fit.command = paste(
      "fit <- gflowx::fit.rdgraph.regression(",
      "  X = X, y = y, k = 6L, max.iterations = 5L,",
      "  n.eigenpairs = 100L, pca.dim = NULL,",
      "  apply.geometric.pruning = FALSE, max.ratio.threshold = 0,",
      "  threshold.percentile = 0, response.penalty.exp = 0,",
      "  use.counting.measure = TRUE, filter.type = \"heat_kernel\",",
      "  dense.fallback = \"never\", compute.extremality = FALSE,",
      "  verbose.level = 0L",
      ")",
      sep = "\n"
    ),
    generate = function() {
      set.seed(2004L)
      n <- 120L
      x <- sort(stats::runif(n, min = -3, max = 3))
      y.true <- rep(0, n)
      y <- as.double(stats::rnorm(n, sd = 1))

      list(
        type = "line",
        X = cbind(x = as.double(x)),
        coord = as.double(x),
        y.true = as.double(y.true),
        responses = list(noise = y),
        metadata = list(
          noise_sd = 1,
          signal = "zero conditional expectation with independent Gaussian response noise"
        )
      )
    }
  )
}

rdgraph.case.two.components.independent.signals <- function() {
  list(
    id = "two_components_independent_smooth_signals",
    title = "Two Components: Independent Smooth Signals",
    group = "topology_stress",
    seed = 3001L,
    k.grid = 3:14,
    fit.args = modifyList(
      rdgraph.default.fit.args(k = 10L, n.eigenpairs = 120L),
      list(
        apply.geometric.pruning = TRUE,
        max.ratio.threshold = 0.1
      )
    ),
    data.command = paste(
      "set.seed(3001L)",
      "n.left <- 70L; n.right <- 70L",
      "x.left <- sort(runif(n.left, -1.45, -0.55))",
      "x.right <- sort(runif(n.right, 0.55, 1.45))",
      "X <- rbind(cbind(x = x.left, y = rnorm(n.left, sd = 0.055)),",
      "           cbind(x = x.right, y = rnorm(n.right, sd = 0.055)))",
      "left.t <- (x.left - min(x.left)) / diff(range(x.left))",
      "right.t <- (x.right - min(x.right)) / diff(range(x.right))",
      "y.true <- c(0.65 + 0.55 * sin(pi * left.t),",
      "            1.45 - 0.45 * cos(pi * right.t))",
      "y <- y.true + rnorm(length(y.true), sd = 0.10)",
      sep = "\n"
    ),
    fit.command = paste(
      "fit <- gflowx::fit.rdgraph.regression(",
      "  X = X, y = y, k = 10L, max.iterations = 5L,",
      "  n.eigenpairs = 120L, pca.dim = NULL,",
      "  apply.geometric.pruning = TRUE, max.ratio.threshold = 0.1,",
      "  threshold.percentile = 0, response.penalty.exp = 0,",
      "  use.counting.measure = TRUE, filter.type = \"heat_kernel\",",
      "  dense.fallback = \"never\", compute.extremality = FALSE,",
      "  verbose.level = 0L",
      ")",
      sep = "\n"
    ),
    generate = function() {
      set.seed(3001L)
      n.left <- 70L
      n.right <- 70L
      x.left <- sort(stats::runif(n.left, -1.45, -0.55))
      x.right <- sort(stats::runif(n.right, 0.55, 1.45))
      X <- rbind(
        cbind(x = x.left, y = stats::rnorm(n.left, sd = 0.055)),
        cbind(x = x.right, y = stats::rnorm(n.right, sd = 0.055))
      )
      left.t <- (x.left - min(x.left)) / diff(range(x.left))
      right.t <- (x.right - min(x.right)) / diff(range(x.right))
      y.true <- c(
        0.65 + 0.55 * sin(pi * left.t),
        1.45 - 0.45 * cos(pi * right.t)
      )
      y <- as.double(y.true + stats::rnorm(length(y.true), sd = 0.10))

      list(
        type = "plane",
        X = X,
        coord = c(left.t, right.t),
        y.true = as.double(y.true),
        responses = list(noisy = y),
        metadata = list(
          noise_sd = 0.10,
          topology = "two well-separated components with independent smooth signals"
        )
      )
    }
  )
}

rdgraph.case.weak.bridge.between.components <- function() {
  list(
    id = "weak_bridge_between_components_tests_bottleneck_leakage",
    title = "Weak Bridge: Bottleneck Leakage Stress",
    group = "topology_stress",
    seed = 3002L,
    k.grid = 3:12,
    fit.args = modifyList(
      rdgraph.default.fit.args(k = 7L, n.eigenpairs = 120L),
      list(
        apply.geometric.pruning = TRUE,
        max.ratio.threshold = 0.1
      )
    ),
    data.command = paste(
      "set.seed(3002L)",
      "n.blob <- 55L; n.bridge <- 28L",
      "left <- cbind(x = rnorm(n.blob, -1.15, 0.16), y = rnorm(n.blob, 0, 0.12))",
      "right <- cbind(x = rnorm(n.blob, 1.15, 0.16), y = rnorm(n.blob, 0, 0.12))",
      "bridge <- cbind(x = seq(-0.85, 0.85, length.out = n.bridge),",
      "                 y = rnorm(n.bridge, 0, 0.025))",
      "X <- rbind(left, bridge, right)",
      "y.true <- 0.45 + 0.95 / (1 + exp(-4 * X[, 'x'])) +",
      "  0.18 * exp(-0.5 * ((X[, 'x'] + 1.15) / 0.22)^2) -",
      "  0.14 * exp(-0.5 * ((X[, 'x'] - 1.15) / 0.25)^2)",
      "y <- y.true + rnorm(length(y.true), sd = 0.10)",
      sep = "\n"
    ),
    fit.command = paste(
      "fit <- gflowx::fit.rdgraph.regression(",
      "  X = X, y = y, k = 7L, max.iterations = 5L,",
      "  n.eigenpairs = 120L, pca.dim = NULL,",
      "  apply.geometric.pruning = TRUE, max.ratio.threshold = 0.1,",
      "  threshold.percentile = 0, response.penalty.exp = 0,",
      "  use.counting.measure = TRUE, filter.type = \"heat_kernel\",",
      "  dense.fallback = \"never\", compute.extremality = FALSE,",
      "  verbose.level = 0L",
      ")",
      sep = "\n"
    ),
    generate = function() {
      set.seed(3002L)
      n.blob <- 55L
      n.bridge <- 28L
      left <- cbind(x = stats::rnorm(n.blob, -1.15, 0.16), y = stats::rnorm(n.blob, 0, 0.12))
      bridge <- cbind(x = seq(-0.85, 0.85, length.out = n.bridge), y = stats::rnorm(n.bridge, 0, 0.025))
      right <- cbind(x = stats::rnorm(n.blob, 1.15, 0.16), y = stats::rnorm(n.blob, 0, 0.12))
      X <- rbind(left, bridge, right)
      y.true <- 0.45 + 0.95 / (1 + exp(-4 * X[, "x"])) +
        0.18 * exp(-0.5 * ((X[, "x"] + 1.15) / 0.22)^2) -
        0.14 * exp(-0.5 * ((X[, "x"] - 1.15) / 0.25)^2)
      y <- as.double(y.true + stats::rnorm(length(y.true), sd = 0.10))

      left.idx <- seq_len(n.blob)
      bridge.idx <- n.blob + seq_len(n.bridge)
      right.idx <- n.blob + n.bridge + seq_len(n.blob)
      left.anchor <- left.idx[which.min(abs(X[left.idx, "x"] - X[head(bridge.idx, 1L), "x"]))]
      right.anchor <- right.idx[which.min(abs(X[right.idx, "x"] - X[tail(bridge.idx, 1L), "x"]))]
      edges <- rbind(
        rdgraph.knn.edges(X, left.idx, k = 4L),
        rdgraph.path.edges(bridge.idx),
        rdgraph.knn.edges(X, right.idx, k = 4L),
        cbind(left.anchor, head(bridge.idx, 1L)),
        cbind(tail(bridge.idx, 1L), right.anchor)
      )
      oracle <- rdgraph.weighted.edge.graph(X, edges)

      list(
        type = "plane",
        X = X,
        coord = X[, "x"],
        y.true = as.double(y.true),
        responses = list(noisy = y),
        oracle.graph = list(
          model = "planted_weak_bridge_graph",
          adj.list = oracle$adj.list,
          weight.list = oracle$weight.list,
          command = paste(
            "oracle.graph <- planted graph with within-blob kNN edges,",
            "  a sparse bridge path, and one connector at each bridge end",
            "oracle.fit <- gflowx::fit.rdgraph.regression(",
            "  X = X, y = y, k = 2L, adj.list = oracle.graph$adj.list,",
            "  weight.list = oracle.graph$weight.list, use.counting.measure = TRUE",
            ")",
            sep = "\n"
          )
        ),
        metadata = list(
          noise_sd = 0.10,
          topology = "two dense blobs joined by a sparse bridge"
        )
      )
    }
  )
}

rdgraph.case.y.branching.geometry <- function() {
  list(
    id = "y_shaped_branching_geometry_branch_specific_signal",
    title = "Y-Shaped Branching Geometry: Branch-Specific Signal",
    group = "topology_stress",
    seed = 3003L,
    k.grid = 3:10,
    fit.args = modifyList(
      rdgraph.default.fit.args(k = 6L, n.eigenpairs = 120L),
      list(
        apply.geometric.pruning = TRUE,
        max.ratio.threshold = 0.1
      )
    ),
    data.command = paste(
      "set.seed(3003L)",
      "geom <- rdgraph.y.geometry(noise = 0.025)",
      "X <- geom$X",
      "y.true <- ifelse(geom$branch == 'stem',",
      "  0.45 + 0.35 * geom$branch.t,",
      "  ifelse(geom$branch == 'left',",
      "    0.82 + 0.42 * geom$branch.t,",
      "    1.18 - 0.36 * geom$branch.t))",
      "y <- y.true + rnorm(length(y.true), sd = 0.09)",
      sep = "\n"
    ),
    fit.command = paste(
      "fit <- gflowx::fit.rdgraph.regression(",
      "  X = X, y = y, k = 6L, max.iterations = 5L,",
      "  n.eigenpairs = 120L, pca.dim = NULL,",
      "  apply.geometric.pruning = TRUE, max.ratio.threshold = 0.1,",
      "  threshold.percentile = 0, response.penalty.exp = 0,",
      "  use.counting.measure = TRUE, filter.type = \"heat_kernel\",",
      "  dense.fallback = \"never\", compute.extremality = FALSE,",
      "  verbose.level = 0L",
      ")",
      sep = "\n"
    ),
    generate = function() {
      set.seed(3003L)
      geom <- rdgraph.y.geometry(noise = 0.025)
      X <- geom$X
      y.true <- ifelse(
        geom$branch == "stem",
        0.45 + 0.35 * geom$branch.t,
        ifelse(
          geom$branch == "left",
          0.82 + 0.42 * geom$branch.t,
          1.18 - 0.36 * geom$branch.t
        )
      )
      y <- as.double(y.true + stats::rnorm(length(y.true), sd = 0.09))
      oracle <- rdgraph.weighted.edge.graph(X, geom$edges)

      list(
        type = "plane",
        X = X,
        coord = geom$branch.t,
        y.true = as.double(y.true),
        responses = list(noisy = y),
        oracle.graph = list(
          model = "planted_y_tree_graph",
          adj.list = oracle$adj.list,
          weight.list = oracle$weight.list,
          command = paste(
            "oracle.graph <- planted Y-shaped tree graph",
            "oracle.fit <- gflowx::fit.rdgraph.regression(",
            "  X = X, y = y, k = 2L, adj.list = oracle.graph$adj.list,",
            "  weight.list = oracle.graph$weight.list, use.counting.measure = TRUE",
            ")",
            sep = "\n"
          )
        ),
        metadata = list(
          noise_sd = 0.09,
          topology = "Y-shaped tree with branch-specific smooth responses"
        )
      )
    }
  )
}

rdgraph.case.y.branch.localized.peak <- function() {
  list(
    id = "y_branch_localized_peak_tests_branch_leakage",
    title = "Y-Shaped Branching Geometry: Localized Branch Feature",
    group = "topology_stress",
    seed = 3004L,
    k.grid = 3:10,
    fit.args = modifyList(
      rdgraph.default.fit.args(k = 6L, n.eigenpairs = 120L),
      list(
        apply.geometric.pruning = TRUE,
        max.ratio.threshold = 0.1
      )
    ),
    data.command = paste(
      "set.seed(3004L)",
      "geom <- rdgraph.y.geometry(noise = 0.025)",
      "X <- geom$X",
      "right.peak <- geom$branch == 'right'",
      "y.true <- 0.35 + 0.18 * pmax(geom$branch.t, 0) +",
      "  ifelse(right.peak, 1.25 * exp(-0.5 * ((geom$branch.t - 0.68) / 0.08)^2), 0)",
      "y <- y.true + rnorm(length(y.true), sd = 0.08)",
      sep = "\n"
    ),
    fit.command = paste(
      "fit <- gflowx::fit.rdgraph.regression(",
      "  X = X, y = y, k = 6L, max.iterations = 5L,",
      "  n.eigenpairs = 120L, pca.dim = NULL,",
      "  apply.geometric.pruning = TRUE, max.ratio.threshold = 0.1,",
      "  threshold.percentile = 0, response.penalty.exp = 0,",
      "  use.counting.measure = TRUE, filter.type = \"heat_kernel\",",
      "  dense.fallback = \"never\", compute.extremality = FALSE,",
      "  verbose.level = 0L",
      ")",
      sep = "\n"
    ),
    generate = function() {
      set.seed(3004L)
      geom <- rdgraph.y.geometry(noise = 0.025)
      X <- geom$X
      right.peak <- geom$branch == "right"
      y.true <- 0.35 + 0.18 * pmax(geom$branch.t, 0) +
        ifelse(right.peak, 1.25 * exp(-0.5 * ((geom$branch.t - 0.68) / 0.08)^2), 0)
      y <- as.double(y.true + stats::rnorm(length(y.true), sd = 0.08))
      oracle <- rdgraph.weighted.edge.graph(X, geom$edges)

      list(
        type = "plane",
        X = X,
        coord = geom$branch.t,
        y.true = as.double(y.true),
        responses = list(noisy = y),
        oracle.graph = list(
          model = "planted_y_tree_graph",
          adj.list = oracle$adj.list,
          weight.list = oracle$weight.list,
          command = paste(
            "oracle.graph <- planted Y-shaped tree graph",
            "oracle.fit <- gflowx::fit.rdgraph.regression(",
            "  X = X, y = y, k = 2L, adj.list = oracle.graph$adj.list,",
            "  weight.list = oracle.graph$weight.list, use.counting.measure = TRUE",
            ")",
            sep = "\n"
          )
        ),
        metadata = list(
          noise_sd = 0.08,
          topology = "Y-shaped tree with a narrow peak on the right arm"
        )
      )
    }
  )
}

rdgraph.case.two.components.componentwise.signals <- function() {
  list(
    id = "two_components_componentwise_smooth_signals",
    title = "Two Components 3B: Component-Wise Smooth Signals",
    group = "topology_stress",
    seed = 3501L,
    k.grid = 3:14,
    fit.args = modifyList(
      rdgraph.default.fit.args(k = 10L, n.eigenpairs = 120L),
      list(apply.geometric.pruning = TRUE, max.ratio.threshold = 0.1)
    ),
    data.command = paste(
      "set.seed(3501L)",
      "n.left <- 80L; n.right <- 65L",
      "x.left <- sort(runif(n.left, -1.55, -0.55))",
      "x.right <- sort(runif(n.right, 0.55, 1.55))",
      "X <- rbind(cbind(x = x.left, y = rnorm(n.left, sd = 0.045)),",
      "           cbind(x = x.right, y = rnorm(n.right, sd = 0.045)))",
      "left.t <- (x.left - min(x.left)) / diff(range(x.left))",
      "right.t <- (x.right - min(x.right)) / diff(range(x.right))",
      "y.true <- c(0.55 + 0.50 * sin(pi * left.t),",
      "            1.25 + 0.35 * sin(2 * pi * right.t))",
      "y <- y.true + rnorm(length(y.true), sd = 0.09)",
      sep = "\n"
    ),
    fit.command = paste(
      "fit <- gflowx::fit.rdgraph.regression(",
      "  X = X, y = y, k = 10L, max.iterations = 5L,",
      "  n.eigenpairs = 120L, pca.dim = NULL,",
      "  apply.geometric.pruning = TRUE, max.ratio.threshold = 0.1,",
      "  threshold.percentile = 0, response.penalty.exp = 0,",
      "  use.counting.measure = TRUE, filter.type = \"heat_kernel\",",
      "  dense.fallback = \"never\", compute.extremality = FALSE,",
      "  verbose.level = 0L",
      ")",
      sep = "\n"
    ),
    generate = function() {
      set.seed(3501L)
      n.left <- 80L
      n.right <- 65L
      x.left <- sort(stats::runif(n.left, -1.55, -0.55))
      x.right <- sort(stats::runif(n.right, 0.55, 1.55))
      X <- rbind(
        cbind(x = x.left, y = stats::rnorm(n.left, sd = 0.045)),
        cbind(x = x.right, y = stats::rnorm(n.right, sd = 0.045))
      )
      left.t <- (x.left - min(x.left)) / diff(range(x.left))
      right.t <- (x.right - min(x.right)) / diff(range(x.right))
      y.true <- c(0.55 + 0.50 * sin(pi * left.t), 1.25 + 0.35 * sin(2 * pi * right.t))
      y <- as.double(y.true + stats::rnorm(length(y.true), sd = 0.09))
      list(
        type = "plane",
        X = X,
        coord = c(left.t, right.t),
        y.true = as.double(y.true),
        responses = list(noisy = y),
        metadata = list(
          noise_sd = 0.09,
          topology = "two intentionally disconnected components requiring component-wise interpretation"
        )
      )
    }
  )
}

rdgraph.case.long.sparse.bridge.sharp.transition <- function() {
  list(
    id = "long_sparse_bridge_sharp_transition",
    title = "Long Sparse Bridge 3B: Sharp Transition",
    group = "topology_stress",
    seed = 3502L,
    k.grid = 2:12,
    fit.args = modifyList(
      rdgraph.default.fit.args(k = 5L, n.eigenpairs = 120L),
      list(apply.geometric.pruning = TRUE, max.ratio.threshold = 0.08)
    ),
    data.command = paste(
      "set.seed(3502L)",
      "n.blob <- 55L; n.bridge <- 14L",
      "left <- cbind(x = rnorm(n.blob, -1.35, 0.13), y = rnorm(n.blob, 0, 0.11))",
      "right <- cbind(x = rnorm(n.blob, 1.35, 0.13), y = rnorm(n.blob, 0, 0.11))",
      "bridge <- cbind(x = seq(-1.05, 1.05, length.out = n.bridge),",
      "                 y = rnorm(n.bridge, 0, 0.018))",
      "X <- rbind(left, bridge, right)",
      "y.true <- 0.55 + 0.95 / (1 + exp(-9 * X[, 'x']))",
      "y <- y.true + rnorm(length(y.true), sd = 0.08)",
      sep = "\n"
    ),
    fit.command = paste(
      "fit <- gflowx::fit.rdgraph.regression(",
      "  X = X, y = y, k = 5L, max.iterations = 5L,",
      "  n.eigenpairs = 120L, pca.dim = NULL,",
      "  apply.geometric.pruning = TRUE, max.ratio.threshold = 0.08,",
      "  threshold.percentile = 0, response.penalty.exp = 0,",
      "  use.counting.measure = TRUE, filter.type = \"heat_kernel\",",
      "  dense.fallback = \"never\", compute.extremality = FALSE,",
      "  verbose.level = 0L",
      ")",
      sep = "\n"
    ),
    generate = function() {
      set.seed(3502L)
      n.blob <- 55L
      n.bridge <- 14L
      left <- cbind(x = stats::rnorm(n.blob, -1.35, 0.13), y = stats::rnorm(n.blob, 0, 0.11))
      bridge <- cbind(x = seq(-1.05, 1.05, length.out = n.bridge), y = stats::rnorm(n.bridge, 0, 0.018))
      right <- cbind(x = stats::rnorm(n.blob, 1.35, 0.13), y = stats::rnorm(n.blob, 0, 0.11))
      X <- rbind(left, bridge, right)
      y.true <- 0.55 + 0.95 / (1 + exp(-9 * X[, "x"]))
      y <- as.double(y.true + stats::rnorm(length(y.true), sd = 0.08))

      left.idx <- seq_len(n.blob)
      bridge.idx <- n.blob + seq_len(n.bridge)
      right.idx <- n.blob + n.bridge + seq_len(n.blob)
      left.anchor <- left.idx[which.min(abs(X[left.idx, "x"] - X[head(bridge.idx, 1L), "x"]))]
      right.anchor <- right.idx[which.min(abs(X[right.idx, "x"] - X[tail(bridge.idx, 1L), "x"]))]
      edges <- rbind(
        rdgraph.knn.edges(X, left.idx, k = 4L),
        rdgraph.path.edges(bridge.idx),
        rdgraph.knn.edges(X, right.idx, k = 4L),
        cbind(left.anchor, head(bridge.idx, 1L)),
        cbind(tail(bridge.idx, 1L), right.anchor)
      )
      oracle <- rdgraph.weighted.edge.graph(X, edges)
      list(
        type = "plane",
        X = X,
        coord = X[, "x"],
        y.true = as.double(y.true),
        responses = list(noisy = y),
        oracle.graph = list(
          model = "planted_long_sparse_bridge_graph",
          adj.list = oracle$adj.list,
          weight.list = oracle$weight.list,
          command = "oracle.graph <- planted long sparse bridge graph"
        ),
        metadata = list(noise_sd = 0.08, topology = "two blobs connected by a long sparse bridge with a sharp response transition")
      )
    }
  )
}

rdgraph.case.y.branch.contrast.metrics <- function() {
  list(
    id = "y_branch_contrast_topology_baseline",
    title = "Y Branch Contrast 3B: Planted Topology Baseline",
    group = "topology_stress",
    seed = 3503L,
    k.grid = 3:10,
    fit.args = modifyList(
      rdgraph.default.fit.args(k = 5L, n.eigenpairs = 120L),
      list(apply.geometric.pruning = TRUE, max.ratio.threshold = 0.08)
    ),
    data.command = paste(
      "set.seed(3503L)",
      "geom <- rdgraph.y.geometry(noise = 0.022)",
      "X <- geom$X",
      "y.true <- ifelse(geom$branch == 'stem', 0.55 + 0.28 * geom$branch.t,",
      "  ifelse(geom$branch == 'left', 0.95 + 0.25 * geom$branch.t,",
      "         1.45 - 0.22 * geom$branch.t))",
      "y <- y.true + rnorm(length(y.true), sd = 0.08)",
      sep = "\n"
    ),
    fit.command = paste(
      "fit <- gflowx::fit.rdgraph.regression(",
      "  X = X, y = y, k = 5L, max.iterations = 5L,",
      "  n.eigenpairs = 120L, pca.dim = NULL,",
      "  apply.geometric.pruning = TRUE, max.ratio.threshold = 0.08,",
      "  threshold.percentile = 0, response.penalty.exp = 0,",
      "  use.counting.measure = TRUE, filter.type = \"heat_kernel\",",
      "  dense.fallback = \"never\", compute.extremality = FALSE,",
      "  verbose.level = 0L",
      ")",
      sep = "\n"
    ),
    generate = function() {
      set.seed(3503L)
      geom <- rdgraph.y.geometry(noise = 0.022)
      X <- geom$X
      y.true <- ifelse(
        geom$branch == "stem",
        0.55 + 0.28 * geom$branch.t,
        ifelse(geom$branch == "left", 0.95 + 0.25 * geom$branch.t, 1.45 - 0.22 * geom$branch.t)
      )
      y <- as.double(y.true + stats::rnorm(length(y.true), sd = 0.08))
      oracle <- rdgraph.weighted.edge.graph(X, geom$edges)
      list(
        type = "plane",
        X = X,
        coord = geom$branch.t,
        y.true = as.double(y.true),
        responses = list(noisy = y),
        oracle.graph = list(model = "planted_y_tree_graph", adj.list = oracle$adj.list, weight.list = oracle$weight.list),
        metadata = list(noise_sd = 0.08, topology = "Y-shaped tree with stronger branch contrast")
      )
    }
  )
}

rdgraph.case.y.branch.narrow.peak.leakage <- function() {
  list(
    id = "y_branch_narrow_peak_leakage_stress",
    title = "Y Branch Narrow Peak 3B: Leakage Stress",
    group = "topology_stress",
    seed = 3504L,
    k.grid = 3:10,
    fit.args = modifyList(
      rdgraph.default.fit.args(k = 5L, n.eigenpairs = 120L),
      list(apply.geometric.pruning = TRUE, max.ratio.threshold = 0.08)
    ),
    data.command = paste(
      "set.seed(3504L)",
      "geom <- rdgraph.y.geometry(noise = 0.022)",
      "X <- geom$X",
      "right.peak <- geom$branch == 'right'",
      "y.true <- 0.35 + 0.12 * geom$branch.t +",
      "  ifelse(right.peak, 1.45 * exp(-0.5 * ((geom$branch.t - 0.70) / 0.055)^2), 0)",
      "y <- y.true + rnorm(length(y.true), sd = 0.07)",
      sep = "\n"
    ),
    fit.command = paste(
      "fit <- gflowx::fit.rdgraph.regression(",
      "  X = X, y = y, k = 5L, max.iterations = 5L,",
      "  n.eigenpairs = 120L, pca.dim = NULL,",
      "  apply.geometric.pruning = TRUE, max.ratio.threshold = 0.08,",
      "  threshold.percentile = 0, response.penalty.exp = 0,",
      "  use.counting.measure = TRUE, filter.type = \"heat_kernel\",",
      "  dense.fallback = \"never\", compute.extremality = FALSE,",
      "  verbose.level = 0L",
      ")",
      sep = "\n"
    ),
    generate = function() {
      set.seed(3504L)
      geom <- rdgraph.y.geometry(noise = 0.022)
      X <- geom$X
      right.peak <- geom$branch == "right"
      y.true <- 0.35 + 0.12 * geom$branch.t +
        ifelse(right.peak, 1.45 * exp(-0.5 * ((geom$branch.t - 0.70) / 0.055)^2), 0)
      y <- as.double(y.true + stats::rnorm(length(y.true), sd = 0.07))
      oracle <- rdgraph.weighted.edge.graph(X, geom$edges)
      list(
        type = "plane",
        X = X,
        coord = geom$branch.t,
        y.true = as.double(y.true),
        responses = list(noisy = y),
        oracle.graph = list(model = "planted_y_tree_graph", adj.list = oracle$adj.list, weight.list = oracle$weight.list),
        metadata = list(noise_sd = 0.07, topology = "Y-shaped tree with a very narrow right-arm peak")
      )
    }
  )
}

rdgraph.case.circular.gaussian.mixture <- function() {
  list(
    id = "circular_gaussian_mixture_recovers_periodic_signal",
    title = "Circular Gaussian Mixture: Periodic Signal Recovery",
    group = "signal_recovery",
    seed = 1002L,
    k.grid = 3:10,
    fit.args = modifyList(
      rdgraph.default.fit.args(k = 5L, n.eigenpairs = 100L),
      list(
        apply.geometric.pruning = TRUE,
        max.ratio.threshold = 0.1
      )
    ),
    data.command = paste(
      "set.seed(1002L)",
      "X.df <- gflow::generate.circle.data(",
      "  n = 160L, radius = 1, noise = 0.08,",
      "  type = \"random\", noise.type = \"normal\", seed = 1002L",
      ")",
      "X <- as.matrix(X.df[, c(\"x\", \"y\")])",
      "theta <- X.df$angles",
      "y.true <- gflow::circular.synthetic.mixture.of.gaussians(",
      "  x = theta,",
      "  x.knot = c(0.45, 2.35, 5.35),",
      "  y.knot = c(1.25, 0.85, 1.7),",
      "  sd.knot = 0.32",
      ")",
      "y <- y.true + rnorm(length(theta), sd = 0.12)",
      sep = "\n"
    ),
    fit.command = paste(
      "fit <- gflowx::fit.rdgraph.regression(",
      "  X = X, y = y, k = 5L, max.iterations = 5L,",
      "  n.eigenpairs = 100L, pca.dim = NULL,",
      "  apply.geometric.pruning = TRUE, max.ratio.threshold = 0.1,",
      "  threshold.percentile = 0, response.penalty.exp = 0,",
      "  use.counting.measure = TRUE, filter.type = \"heat_kernel\",",
      "  dense.fallback = \"never\", compute.extremality = FALSE,",
      "  verbose.level = 0L",
      ")",
      sep = "\n"
    ),
    generate = function() {
      set.seed(1002L)
      n <- 160L
      X.df <- generate.circle.data(
        n = n,
        radius = 1,
        noise = 0.08,
        type = "random",
        noise.type = "normal",
        seed = 1002L
      )
      X <- as.matrix(X.df[, c("x", "y")])
      theta <- as.double(X.df$angles)

      y.true <- circular.synthetic.mixture.of.gaussians(
        x = theta,
        x.knot = c(0.45, 2.35, 5.35),
        y.knot = c(1.25, 0.85, 1.7),
        sd.knot = 0.32
      )
      y <- as.double(y.true + stats::rnorm(n, sd = 0.12))

      list(
        type = "circle",
        X = X,
        coord = theta,
        theta = theta,
        y.true = as.double(y.true),
        responses = list(noisy = y),
        metadata = list(
          noise_sd = 0.12,
          geometry_noise = 0.08,
          signal = "three-component circular Gaussian mixture"
        )
      )
    }
  )
}

rdgraph.case.shuffled.response <- function() {
  list(
    id = "shuffled_response_does_not_look_oracular",
    title = "Shuffled Response: Negative Control",
    group = "negative_control",
    seed = 1003L,
    k.grid = 3:10,
    fit.args = rdgraph.default.fit.args(k = 6L, n.eigenpairs = 100L),
    data.command = paste(
      "set.seed(1003L)",
      "gm <- gflow::generate.1d.gaussian.mixture(",
      "  n.points = 120L,",
      "  x.knots = c(-4, 0.5, 5),",
      "  y.knots = c(1.4, 3.2, 1.7),",
      "  sd.knot = 0.9, x.offset = 1.5, add.noise = FALSE",
      ")",
      "X <- cbind(x = gm$x)",
      "y.true <- gm$y.true",
      "y <- y.true + rnorm(length(gm$x), sd = 0.35)",
      "y.shuffled <- sample(y, length(y), replace = FALSE)",
      sep = "\n"
    ),
    fit.command = paste(
      "real.fit <- gflowx::fit.rdgraph.regression(X, y, k = 6L, ...)",
      "shuffled.fit <- gflowx::fit.rdgraph.regression(X, y.shuffled, k = 6L, ...)",
      sep = "\n"
    ),
    generate = function() {
      set.seed(1003L)
      gm <- generate.1d.gaussian.mixture(
        n.points = 120L,
        x.knots = c(-4, 0.5, 5),
        y.knots = c(1.4, 3.2, 1.7),
        sd.knot = 0.9,
        x.offset = 1.5,
        add.noise = FALSE
      )
      y <- as.double(gm$y.true + stats::rnorm(length(gm$x), sd = 0.35))
      y.shuffled <- sample(y, length(y), replace = FALSE)

      list(
        type = "line",
        X = cbind(x = as.double(gm$x)),
        coord = as.double(gm$x),
        y.true = as.double(gm$y.true),
        responses = list(real = y, shuffled = y.shuffled),
        metadata = list(
          noise_sd = 0.35,
          signal = "1D Gaussian mixture with response permutation negative control"
        )
      )
    }
  )
}

rdgraph.case.2d.gaussian.mixture <- function() {
  list(
    id = "2d_gaussian_mixture_recovers_smooth_signal",
    title = "2D Gaussian Mixture: Smooth Surface Recovery",
    group = "signal_recovery",
    seed = 1004L,
    k.grid = 4:11,
    fit.args = modifyList(
      rdgraph.default.fit.args(k = 7L, n.eigenpairs = 100L),
      list(
        apply.geometric.pruning = TRUE,
        max.ratio.threshold = 0.1
      )
    ),
    data.command = paste(
      "set.seed(1004L)",
      "mix <- gflow::create.gaussian.mixture(",
      "  x1 = 0.25, y1 = 0.25, x2 = 0.75, y2 = 0.72,",
      "  A1 = 0.75, A2 = 1.4, sigma = 0.18",
      ")",
      "grid <- expand.grid(x = seq(0, 1, length.out = 34),",
      "                    y = seq(0, 1, length.out = 34))",
      "idx <- sample(seq_len(nrow(grid)), size = 180L)",
      "X <- as.matrix(grid[idx, ])",
      "y.true <- apply(X, 1L, function(pt) mix$f(pt[1], pt[2]))",
      "y <- y.true + rnorm(length(y.true), sd = 0.10)",
      sep = "\n"
    ),
    fit.command = paste(
      "fit <- gflowx::fit.rdgraph.regression(",
      "  X = X, y = y, k = 7L, max.iterations = 5L,",
      "  n.eigenpairs = 100L, pca.dim = NULL,",
      "  apply.geometric.pruning = TRUE, max.ratio.threshold = 0.1,",
      "  threshold.percentile = 0, response.penalty.exp = 0,",
      "  use.counting.measure = TRUE, filter.type = \"heat_kernel\",",
      "  dense.fallback = \"never\", compute.extremality = FALSE,",
      "  verbose.level = 0L",
      ")",
      sep = "\n"
    ),
    generate = function() {
      set.seed(1004L)
      mix <- create.gaussian.mixture(
        x1 = 0.25, y1 = 0.25,
        x2 = 0.75, y2 = 0.72,
        A1 = 0.75, A2 = 1.4,
        sigma = 0.18
      )
      grid <- expand.grid(
        x = seq(0, 1, length.out = 34L),
        y = seq(0, 1, length.out = 34L)
      )
      idx <- sample(seq_len(nrow(grid)), size = 180L)
      X <- as.matrix(grid[idx, ])
      y.true <- apply(X, 1L, function(pt) mix$f(pt[1], pt[2]))
      y <- as.double(y.true + stats::rnorm(length(y.true), sd = 0.10))

      list(
        type = "plane",
        X = X,
        coord = seq_along(y.true),
        y.true = as.double(y.true),
        responses = list(noisy = y),
        metadata = list(
          noise_sd = 0.10,
          signal = "two-component 2D Gaussian mixture"
        )
      )
    }
  )
}

rdgraph.regression.correctness.cases <- function() {
  cases <- list(
    rdgraph.case.1d.gaussian.mixture(),
    rdgraph.case.circular.gaussian.mixture(),
    rdgraph.case.shuffled.response(),
    rdgraph.case.2d.gaussian.mixture(),
    rdgraph.case.nonuniform.circle(),
    rdgraph.case.circle.with.gap(),
    rdgraph.case.1d.boundary.peaks(),
    rdgraph.case.pure.noise.response(),
    rdgraph.case.two.components.independent.signals(),
    rdgraph.case.weak.bridge.between.components(),
    rdgraph.case.y.branching.geometry(),
    rdgraph.case.y.branch.localized.peak(),
    rdgraph.case.two.components.componentwise.signals(),
    rdgraph.case.long.sparse.bridge.sharp.transition(),
    rdgraph.case.y.branch.contrast.metrics(),
    rdgraph.case.y.branch.narrow.peak.leakage()
  )
  rdgraph.assign.case.rounds(cases)
}

rdgraph.correctness.rounds <- function() {
  list(
    round_1 = list(
      id = "round_1",
      title = "Round 1: Initial Synthetic Examples",
      file = "round_1_initial_synthetic.html",
      description = "Baseline synthetic recovery, periodic recovery, shuffled-response negative control, and 2D smooth surface recovery."
    ),
    round_2 = list(
      id = "round_2",
      title = "Round 2: Geometry, Boundary, And Negative-Control Stress",
      file = "round_2_geometry_boundary_negative_controls.html",
      description = "Focused stress cases for density imbalance, sampling gaps, endpoint bias, and pure-noise response behavior."
    ),
    round_3 = list(
      id = "round_3",
      title = "Round 3: Topology And Bridge Stress",
      file = "round_3_topology_bridge_stress.html",
      description = "Topology stress cases for disconnected components, weak bridges, branching geometry, and localized branch features."
    ),
    round_3b = list(
      id = "round_3b",
      title = "Round 3B: Revised Topology Diagnostics",
      file = "round_3b_revised_topology_diagnostics.html",
      description = "Revised topology cases using component-wise disconnected-graph diagnostics, sharper bridge stress, and stronger branch-leakage challenges."
    )
  )
}

rdgraph.assign.case.rounds <- function(cases) {
  round.1.ids <- c(
    "1d_gaussian_mixture_recovers_smooth_signal",
    "circular_gaussian_mixture_recovers_periodic_signal",
    "shuffled_response_does_not_look_oracular",
    "2d_gaussian_mixture_recovers_smooth_signal"
  )
  round.2.ids <- c(
    "nonuniform_circle_recovers_periodic_signal",
    "circle_with_gap_reports_connectivity_and_wraparound_limits",
    "1d_boundary_peaks_exposes_endpoint_bias",
    "pure_noise_response_does_not_look_structured"
  )
  round.3.ids <- c(
    "two_components_independent_smooth_signals",
    "weak_bridge_between_components_tests_bottleneck_leakage",
    "y_shaped_branching_geometry_branch_specific_signal",
    "y_branch_localized_peak_tests_branch_leakage"
  )
  round.3b.ids <- c(
    "two_components_componentwise_smooth_signals",
    "long_sparse_bridge_sharp_transition",
    "y_branch_contrast_topology_baseline",
    "y_branch_narrow_peak_leakage_stress"
  )
  rounds <- rdgraph.correctness.rounds()

  lapply(cases, function(case) {
    if (case$id %in% round.1.ids) {
      case$round <- rounds$round_1
    } else if (case$id %in% round.2.ids) {
      case$round <- rounds$round_2
    } else if (case$id %in% round.3.ids) {
      case$round <- rounds$round_3
    } else if (case$id %in% round.3b.ids) {
      case$round <- rounds$round_3b
    } else {
      case$round <- list(
        id = "unassigned",
        title = "Unassigned Cases",
        file = "unassigned_cases.html",
        description = "Cases without explicit round metadata."
      )
    }
    case
  })
}
