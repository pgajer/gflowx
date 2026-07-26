#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x

find.repo.root <- function(start = getwd()) {
  cur <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(cur, "DESCRIPTION"))) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) stop("Could not locate repository root (missing DESCRIPTION).", call. = FALSE)
    cur <- parent
  }
}

parse.cli.args <- function(args) {
  out <- list()
  if (length(args) == 0L) return(out)
  for (a in args) {
    if (!grepl("^--[^=]+=", a)) next
    key <- sub("^--([^=]+)=.*$", "\\1", a)
    val <- sub("^--[^=]+=(.*)$", "\\1", a)
    out[[key]] <- val
  }
  out
}

parse.char.vec <- function(x, default) {
  if (is.null(x) || !nzchar(x)) return(default)
  v <- trimws(unlist(strsplit(x, ",", fixed = TRUE)))
  v <- v[nzchar(v)]
  if (length(v) == 0L) return(default)
  v
}

as.int.arg <- function(x, default) {
  if (is.null(x)) return(as.integer(default))
  y <- suppressWarnings(as.integer(x))
  if (is.na(y) || y < 1L) return(as.integer(default))
  y
}

as.num.arg <- function(x, default) {
  if (is.null(x)) return(as.numeric(default))
  y <- suppressWarnings(as.numeric(x))
  if (!is.finite(y)) return(as.numeric(default))
  y
}

format.seconds <- function(sec) {
  if (!is.finite(sec) || sec < 0) return("NA")
  h <- floor(sec / 3600)
  m <- floor((sec %% 3600) / 60)
  s <- floor(sec %% 60)
  sprintf("%02dh:%02dm:%02ds", h, m, s)
}

compute.filter.weights <- function(eigenvalues, eta, filter.type = "heat_kernel") {
  switch(
    filter.type,
    "heat_kernel" = exp(-eta * eigenvalues),
    "tikhonov" = 1.0 / (1.0 + eta * eigenvalues),
    "cubic_spline" = 1.0 / (1.0 + eta * eigenvalues^2),
    "gaussian" = exp(-eta * eigenvalues^2),
    "exponential" = exp(-eta * sqrt(pmax(eigenvalues, 0))),
    "butterworth" = {
      ratio <- eigenvalues / max(eta, 1e-15)
      1.0 / (1.0 + ratio^4)
    },
    exp(-eta * eigenvalues)
  )
}

set.model.eta <- function(fit, eta) {
  fit2 <- fit
  filter.type <- fit2$spectral$filter.type %||% "heat_kernel"
  eigenvalues <- as.double(fit2$spectral$eigenvalues)
  fit2$spectral$eta.optimal <- as.double(eta)
  fit2$spectral$filtered.eigenvalues <- compute.filter.weights(
    eigenvalues = eigenvalues,
    eta = as.double(eta),
    filter.type = filter.type
  )
  fit2
}

safe.col.scales <- function(X.train) {
  med <- apply(X.train, 2L, stats::mad, constant = 1.4826)
  sdv <- apply(X.train, 2L, stats::sd)
  out <- med
  bad <- !is.finite(out) | out <= sqrt(.Machine$double.eps)
  out[bad] <- sdv[bad]
  bad <- !is.finite(out) | out <= sqrt(.Machine$double.eps)
  out[bad] <- 1.0
  as.double(out)
}

weighted.mean.safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  sum(x[ok] * w[ok]) / sum(w[ok])
}

weighted.se.safe <- function(x, w) {
  ok <- is.finite(x) & is.finite(w) & w > 0
  if (!any(ok)) return(NA_real_)
  x <- x[ok]
  w <- w[ok]
  sw <- sum(w)
  if (sw <= 0) return(NA_real_)
  mu <- sum(w * x) / sw
  var.w <- sum(w * (x - mu)^2) / sw
  ess <- (sw^2) / sum(w^2)
  if (!is.finite(ess) || ess <= 1) return(0.0)
  sqrt(var.w / ess)
}

sqdist2 <- function(A, B) {
  A <- as.matrix(A)
  B <- as.matrix(B)
  d2 <- outer(rowSums(A^2), rowSums(B^2), "+") - 2 * (A %*% t(B))
  pmax(d2, 0)
}

predict.knn.mean <- function(X.train, ord.idx, k) {
  X.train <- as.matrix(X.train)
  k.use <- min(max(1L, as.integer(round(k))), nrow(X.train))
  if (is.null(dim(ord.idx))) ord.idx <- matrix(ord.idx, nrow = 1L)
  n.test <- nrow(ord.idx)
  pred <- matrix(NA_real_, nrow = n.test, ncol = ncol(X.train))
  for (i in seq_len(n.test)) {
    idx <- ord.idx[i, seq_len(k.use)]
    pred[i, ] <- colMeans(X.train[idx, , drop = FALSE])
  }
  pred
}

predict.gaussian.kernel <- function(X.train, d2, bandwidth) {
  X.train <- as.matrix(X.train)
  h <- max(as.numeric(bandwidth), 1e-6)
  W <- exp(-d2 / (2 * h * h))
  sw <- rowSums(W)
  bad <- !is.finite(sw) | sw <= 1e-12
  sw[bad] <- 1.0
  pred <- sweep(W %*% X.train, 1L, sw, "/")
  if (any(bad)) {
    nn.idx <- max.col(-d2, ties.method = "first")
    pred[bad, ] <- X.train[nn.idx[bad], , drop = FALSE]
  }
  pred
}

predict.classical.full <- function(X, method.id, param.value) {
  X <- as.matrix(X)
  n <- nrow(X)
  if (n < 2L) return(X)

  if (identical(method.id, "classical_knn_mean_cv")) {
    d2 <- sqdist2(X, X)
    diag(d2) <- Inf
    ord.idx <- t(apply(d2, 1L, order))
    if (is.null(dim(ord.idx))) ord.idx <- matrix(ord.idx, nrow = 1L)
    return(predict.knn.mean(X.train = X, ord.idx = ord.idx, k = param.value))
  }

  if (identical(method.id, "classical_gaussian_kernel_cv")) {
    d2 <- sqdist2(X, X)
    diag(d2) <- Inf
    nn.dist <- sqrt(apply(d2, 1L, min))
    base.nn <- stats::median(nn.dist[is.finite(nn.dist)], na.rm = TRUE)
    if (!is.finite(base.nn) || base.nn <= 0) base.nn <- 0.25
    h <- max(as.numeric(param.value) * base.nn, 1e-6)
    W <- exp(-d2 / (2 * h * h))
    diag(W) <- 0
    sw <- rowSums(W)
    bad <- !is.finite(sw) | sw <= 1e-12
    sw[bad] <- 1.0
    pred <- sweep(W %*% X, 1L, sw, "/")
    if (any(bad)) {
      nn.idx <- max.col(-d2, ties.method = "first")
      pred[bad, ] <- X[nn.idx[bad], , drop = FALSE]
    }
    return(pred)
  }

  stop("Unsupported method id: ", method.id, call. = FALSE)
}

evaluate.classical.cv <- function(X, splits, classical.grid) {
  X <- as.matrix(X)
  method.spec <- list(
    classical_knn_mean_cv = list(param_name = "k", param_grid = as.double(classical.grid$knn_k_grid)),
    classical_gaussian_kernel_cv = list(param_name = "bw_mult", param_grid = as.double(classical.grid$gaussian_bw_mult_grid))
  )

  candidate.out <- list()
  selection.out <- list()
  cand.ptr <- 0L
  sel.ptr <- 0L

  for (method.id in names(method.spec)) {
    grid <- method.spec[[method.id]]$param_grid
    param.name <- method.spec[[method.id]]$param_name
    if (length(grid) == 0L) next

    scaled.mat <- matrix(NA_real_, nrow = length(splits), ncol = length(grid))
    raw.mat <- matrix(NA_real_, nrow = length(splits), ncol = length(grid))
    w.vec <- numeric(length(splits))

    for (sidx in seq_along(splits)) {
      split <- splits[[sidx]]
      X.train <- X[split$train.idx, , drop = FALSE]
      X.test <- X[split$test.idx, , drop = FALSE]
      w.vec[sidx] <- nrow(X.test)

      d2 <- sqdist2(X.test, X.train)
      ord.idx <- t(apply(d2, 1L, order))
      if (is.null(dim(ord.idx))) ord.idx <- matrix(ord.idx, nrow = 1L)

      d2.train <- sqdist2(X.train, X.train)
      diag(d2.train) <- Inf
      nn.dist <- sqrt(apply(d2.train, 1L, min))
      base.nn <- stats::median(nn.dist[is.finite(nn.dist)], na.rm = TRUE)
      if (!is.finite(base.nn) || base.nn <= 0) base.nn <- 0.25

      for (pidx in seq_along(grid)) {
        param.val <- grid[pidx]
        pred <- switch(
          method.id,
          classical_knn_mean_cv = predict.knn.mean(X.train = X.train, ord.idx = ord.idx, k = param.val),
          classical_gaussian_kernel_cv = predict.gaussian.kernel(X.train = X.train, d2 = d2, bandwidth = param.val * base.nn),
          stop("Unsupported method: ", method.id, call. = FALSE)
        )

        resid <- X.test - pred
        resid.scaled <- sweep(resid, 2L, split$col.scale, "/")
        scaled.mat[sidx, pidx] <- mean(resid.scaled^2)
        raw.mat[sidx, pidx] <- mean(resid^2)
      }
    }

    cand <- data.frame(
      method = method.id,
      param_name = param.name,
      param_value = grid,
      scaled_mse_mean = vapply(seq_along(grid), function(j) weighted.mean.safe(scaled.mat[, j], w.vec), numeric(1)),
      scaled_mse_se = vapply(seq_along(grid), function(j) weighted.se.safe(scaled.mat[, j], w.vec), numeric(1)),
      raw_mse_mean = vapply(seq_along(grid), function(j) weighted.mean.safe(raw.mat[, j], w.vec), numeric(1)),
      n_eval_rows = vapply(seq_along(grid), function(j) sum(is.finite(scaled.mat[, j])), integer(1)),
      stringsAsFactors = FALSE
    )

    cand <- cand[is.finite(cand$scaled_mse_mean), , drop = FALSE]
    cand <- cand[order(cand$scaled_mse_mean, cand$raw_mse_mean), , drop = FALSE]
    if (nrow(cand) == 0L) next

    cand.ptr <- cand.ptr + 1L
    candidate.out[[cand.ptr]] <- cand

    best <- cand[1L, , drop = FALSE]
    sel.ptr <- sel.ptr + 1L
    selection.out[[sel.ptr]] <- data.frame(
      method = method.id,
      selected = "min",
      param_name = best$param_name,
      param_value = best$param_value,
      scaled_mse_mean = best$scaled_mse_mean,
      scaled_mse_se = best$scaled_mse_se,
      raw_mse_mean = best$raw_mse_mean,
      stringsAsFactors = FALSE
    )
  }

  list(
    candidates = if (cand.ptr > 0L) do.call(rbind, candidate.out[seq_len(cand.ptr)]) else data.frame(),
    selection = if (sel.ptr > 0L) do.call(rbind, selection.out[seq_len(sel.ptr)]) else data.frame()
  )
}

evaluate.graph.split <- function(model, X, split, model.idx, block.size = NULL) {
  model.idx <- as.integer(model.idx)
  X.model <- X[model.idx, , drop = FALSE]

  train.local <- match(split$train.idx, model.idx, nomatch = 0L)
  train.local <- as.integer(train.local[train.local > 0L])
  if (length(train.local) < 2L) stop("Too few training rows inside modeled component.")

  test.local.all <- match(split$test.idx, model.idx, nomatch = 0L)
  modeled.mask <- as.logical(test.local.all > 0L)
  test.local <- as.integer(test.local.all[modeled.mask])
  if (length(test.local) < 1L) stop("No test rows inside modeled component.")

  refit <- gflowx::refit.rdgraph.regression(
    fitted.model = model,
    y.new = X.model,
    y.vertices = train.local,
    per.column.gcv = FALSE,
    block.size = block.size,
    verbose = FALSE
  )

  X.hat.model <- as.matrix(refit$fitted.values)
  X.pred.test <- X[split$test.idx, , drop = FALSE]
  X.pred.test[modeled.mask, ] <- X.hat.model[test.local, , drop = FALSE]

  resid <- X[split$test.idx, , drop = FALSE] - X.pred.test
  resid.scaled <- sweep(resid, 2L, split$col.scale, "/")

  list(
    scaled.mse = mean(resid.scaled^2),
    raw.mse = mean(resid^2),
    n_test_total = length(split$test.idx),
    n_test_modeled = length(test.local),
    n_train_total = length(split$train.idx),
    n_train_modeled = length(train.local),
    modeled_fraction = length(test.local) / max(length(split$test.idx), 1L)
  )
}

draw.radial.noise <- function(n, noise, noise.type = "normal") {
  if (!is.finite(noise) || noise <= 0) return(rep(0.0, n))
  if (identical(noise.type, "laplace")) {
    u <- stats::runif(n, min = -0.5, max = 0.5)
    return(-(noise / sqrt(2)) * sign(u) * log1p(-2 * abs(u)))
  }
  stats::rnorm(n, mean = 0, sd = noise)
}

orthonormal.two.cols <- function(p, seed) {
  set.seed(as.integer(seed))
  M <- matrix(stats::rnorm(p * 2L), nrow = p, ncol = 2L)
  Q <- qr.Q(qr(M))
  Q[, 1:2, drop = FALSE]
}

generate.highdim.circle <- function(n, p, radius, sigma, embedding.type,
                                    ambient.noise.factor = 1.0,
                                    noise.type = "normal",
                                    seed = NULL) {
  if (!is.null(seed)) set.seed(as.integer(seed))

  angles <- sort(stats::runif(n, min = 0, max = 2 * pi))
  eps <- draw.radial.noise(n, noise = sigma, noise.type = noise.type)

  s.true <- cbind(radius * cos(angles), radius * sin(angles))
  s.obs <- cbind((radius + eps) * cos(angles), (radius + eps) * sin(angles))

  if (identical(embedding.type, "sparse")) {
    basis <- matrix(0, nrow = p, ncol = 2)
    basis[1, 1] <- 1
    basis[2, 2] <- 1
  } else if (identical(embedding.type, "mixed")) {
    basis <- orthonormal.two.cols(p = p, seed = as.integer(seed + 17L))
  } else {
    stop("Unsupported embedding.type: ", embedding.type, call. = FALSE)
  }

  X.true <- s.true %*% t(basis)
  X.signal.obs <- s.obs %*% t(basis)

  ambient.sd <- sigma * ambient.noise.factor
  if (identical(embedding.type, "sparse")) {
    ambient <- matrix(0, nrow = n, ncol = p)
    if (p > 2L && ambient.sd > 0) {
      ambient[, 3:p] <- matrix(stats::rnorm(n * (p - 2L), mean = 0, sd = ambient.sd), nrow = n, ncol = p - 2L)
    }
  } else {
    ambient <- if (ambient.sd > 0) {
      matrix(stats::rnorm(n * p, mean = 0, sd = ambient.sd), nrow = n, ncol = p)
    } else {
      matrix(0, nrow = n, ncol = p)
    }
  }

  X.obs <- X.signal.obs + ambient

  list(
    X.obs = X.obs,
    X.true = X.true,
    basis = basis,
    angles = angles,
    s.true = s.true,
    s.obs = s.obs
  )
}

compute.oracle.metrics <- function(X.obs, X.hat, X.true, basis,
                                   method, embedding.type, p, sigma, replicate,
                                   notes = NA_character_, modeled_fraction = NA_real_) {
  X.obs <- as.matrix(X.obs)
  X.hat <- as.matrix(X.hat)
  X.true <- as.matrix(X.true)
  basis <- as.matrix(basis)

  latent.obs <- X.obs %*% basis
  latent.hat <- X.hat %*% basis
  latent.true <- X.true %*% basis

  data.frame(
    embedding_type = embedding.type,
    p = as.integer(p),
    sigma = as.double(sigma),
    replicate = as.integer(replicate),
    method = method,
    n_rows = nrow(X.hat),
    n_cols = ncol(X.hat),
    modeled_fraction = modeled_fraction,
    rmse_all_obs = sqrt(mean((X.obs - X.true)^2)),
    rmse_all_hat = sqrt(mean((X.hat - X.true)^2)),
    rmse_latent_obs = sqrt(mean((latent.obs - latent.true)^2)),
    rmse_latent_hat = sqrt(mean((latent.hat - latent.true)^2)),
    var_ratio_mean = mean(apply(X.hat, 2L, stats::var) / pmax(apply(X.obs, 2L, stats::var), 1e-12)),
    notes = notes,
    stringsAsFactors = FALSE
  )
}

append.csv <- function(df, path) {
  if (is.null(df) || nrow(df) == 0L) return(invisible(NULL))
  write.table(
    df,
    file = path,
    append = file.exists(path),
    sep = ",",
    row.names = FALSE,
    col.names = !file.exists(path),
    quote = TRUE
  )
}

mode.param <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  tab <- sort(table(x), decreasing = TRUE)
  as.numeric(names(tab)[1L])
}

build.scenario.grid <- function(config) {
  sg <- expand.grid(
    embedding_type = config$embedding_grid,
    p = as.integer(config$p_grid),
    sigma = as.double(config$sigma_grid),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )

  sg$n_replicates <- ifelse(
    sg$p >= config$ultra_p_threshold,
    config$n_replicates_ultra,
    ifelse(sg$p >= config$high_p_threshold, config$n_replicates_high, config$n_replicates_main)
  )

  sg <- sg[order(match(sg$embedding_type, config$embedding_grid), sg$p, sg$sigma), ]
  rownames(sg) <- NULL
  sg
}

repo.root <- find.repo.root()
setwd(repo.root)

suppressPackageStartupMessages({
  loaded <- FALSE
  if (requireNamespace("pkgload", quietly = TRUE)) {
    loaded <- isTRUE(tryCatch({
      pkgload::load_all(repo.root, quiet = TRUE, export_all = FALSE)
      TRUE
    }, error = function(e) FALSE))
  }
  if (!loaded) {
    library(gflow, quietly = TRUE, warn.conflicts = FALSE)
  }
})

args <- parse.cli.args(commandArgs(trailingOnly = TRUE))
mode <- args$mode %||% "full"
base.seed <- as.int.arg(args$base_seed, 20260218L)

config <- list(
  experiment_name = "noisy_circle_cv_phase3_highdim",
  mode = mode,
  base_seed = base.seed,
  n_points = 220L,
  radius = 1.0,
  noise_type = "normal",
  embedding_grid = c("sparse", "mixed"),
  p_grid = c(100L, 500L, 1000L, 10000L),
  sigma_grid = c(0.05, 0.10, 0.20),
  n_replicates_main = 8L,
  n_replicates_high = 6L,
  n_replicates_ultra = 4L,
  high_p_threshold = 1000L,
  ultra_p_threshold = 10000L,
  n_folds = 5L,
  n_repeats = 3L,
  k_grid = c(5L, 9L, 13L, 17L, 21L),
  n_eigenpairs_grid = c(25L, 50L),
  eta_grid = exp(seq(log(1e-3), log(1.5), length.out = 10L)),
  classical_knn_k_grid = c(5L, 9L, 15L, 25L, 35L),
  classical_gaussian_bw_mult_grid = c(0.40, 0.70, 1.00, 1.50, 2.20),
  filter_type = "heat_kernel",
  use_counting_measure = TRUE,
  response_penalty_exp = 0.0,
  t_scale_factor = 0.5,
  beta_coef_factor = 0.1,
  max_iterations = 8L,
  pca_dim = 100L,
  variance_explained = 0.99,
  max_path_edge_ratio_deviation_thld = 0.1,
  path_edge_ratio_percentile = 0.5,
  threshold_percentile = 0.0,
  ambient_noise_factor = 1.0,
  block_size = NULL,
  log_every_tasks = 100L
)

if (mode == "pilot") {
  config$p_grid <- c(100L, 500L)
  config$n_replicates_main <- 2L
  config$n_replicates_high <- 1L
  config$n_replicates_ultra <- 1L
  config$n_repeats <- 2L
  config$n_folds <- 3L
  config$eta_grid <- exp(seq(log(1e-3), log(1.5), length.out = 6L))
}

if (!is.null(args$embedding_grid)) config$embedding_grid <- parse.char.vec(args$embedding_grid, config$embedding_grid)
if (!is.null(args$p_grid)) config$p_grid <- as.integer(parse.char.vec(args$p_grid, config$p_grid))
if (!is.null(args$sigma_grid)) config$sigma_grid <- as.numeric(parse.char.vec(args$sigma_grid, config$sigma_grid))
if (!is.null(args$n_replicates_main)) config$n_replicates_main <- as.int.arg(args$n_replicates_main, config$n_replicates_main)
if (!is.null(args$n_replicates_high)) config$n_replicates_high <- as.int.arg(args$n_replicates_high, config$n_replicates_high)
if (!is.null(args$n_replicates_ultra)) config$n_replicates_ultra <- as.int.arg(args$n_replicates_ultra, config$n_replicates_ultra)
if (!is.null(args$n_repeats)) config$n_repeats <- as.int.arg(args$n_repeats, config$n_repeats)
if (!is.null(args$n_folds)) config$n_folds <- as.int.arg(args$n_folds, config$n_folds)
if (!is.null(args$ambient_noise_factor)) config$ambient_noise_factor <- as.num.arg(args$ambient_noise_factor, config$ambient_noise_factor)
if (!is.null(args$log_every_tasks)) config$log_every_tasks <- as.int.arg(args$log_every_tasks, config$log_every_tasks)

scenario.grid <- build.scenario.grid(config)
config$scenario_grid <- scenario.grid

run.id <- args$run_id %||% format(Sys.time(), "%Y%m%d_%H%M%S")
run.dir <- file.path(repo.root, "tests", "manual", "results", "noisy_circle_cv_phase3_highdim", run.id)
dir.create(run.dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run.dir, "logs"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run.dir, "tables"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run.dir, "checkpoints"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run.dir, "artifacts"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(run.dir, "config"), recursive = TRUE, showWarnings = FALSE)

log.path <- file.path(run.dir, "logs", "progress.log")
log.con <- file(log.path, open = "a")
on.exit(close(log.con), add = TRUE)

log.msg <- function(...) {
  txt <- paste(...)
  line <- sprintf("[%s] %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), txt)
  cat(line, "\n")
  writeLines(line, con = log.con, sep = "\n")
  flush(log.con)
  flush.console()
}

saveRDS(config, file.path(run.dir, "config", "experiment_config.rds"))
writeLines(capture.output(str(config, max.level = 2L)), file.path(run.dir, "config", "experiment_config.txt"))
write.csv(scenario.grid, file.path(run.dir, "config", "scenario_grid.csv"), row.names = FALSE)

candidate.grid <- expand.grid(
  k = as.integer(config$k_grid),
  n_eigenpairs = as.integer(config$n_eigenpairs_grid),
  eta = as.double(config$eta_grid),
  KEEP.OUT.ATTRS = FALSE,
  stringsAsFactors = FALSE
)
candidate.grid <- candidate.grid[order(candidate.grid$k, candidate.grid$n_eigenpairs, candidate.grid$eta), ]

pair.grid <- unique(candidate.grid[, c("k", "n_eigenpairs")])
split.count <- as.integer(config$n_repeats * config$n_folds)
tasks.per.replicate.graph <- nrow(candidate.grid) * split.count
total.replicates <- sum(as.integer(scenario.grid$n_replicates))
total.tasks <- total.replicates * tasks.per.replicate.graph

log.msg("Starting phase-3 high-dimensional experiment run.")
log.msg("Run directory:", run.dir)
log.msg("Graph CV tasks:", total.tasks,
        "| scenarios:", nrow(scenario.grid),
        "| total replicates:", total.replicates,
        "| candidates:", nrow(candidate.grid),
        "| splits/replicate:", split.count)

cv.fold.path <- file.path(run.dir, "tables", "graph_cv_fold_scores.csv")
cv.summary.path <- file.path(run.dir, "tables", "graph_cv_candidate_summary.csv")
selection.graph.path <- file.path(run.dir, "tables", "graph_selection_by_replicate.csv")
selection.classical.path <- file.path(run.dir, "tables", "classical_selection_by_replicate.csv")
cv.classical.path <- file.path(run.dir, "tables", "classical_cv_candidate_summary.csv")
oracle.path <- file.path(run.dir, "tables", "oracle_metrics_by_replicate.csv")
runtime.path <- file.path(run.dir, "tables", "runtime_by_replicate.csv")

global.task.idx <- 0L
run.start.time <- Sys.time()

for (scenario.idx in seq_len(nrow(scenario.grid))) {
  embedding.type <- as.character(scenario.grid$embedding_type[scenario.idx])
  p.dim <- as.integer(scenario.grid$p[scenario.idx])
  sigma <- as.double(scenario.grid$sigma[scenario.idx])
  n.rep <- as.integer(scenario.grid$n_replicates[scenario.idx])

  for (rep.idx in seq_len(n.rep)) {
    rep.seed <- as.integer(base.seed + scenario.idx * 100000L + rep.idx)
    set.seed(rep.seed)

    log.msg(sprintf("Embedding %s | p=%d | sigma %.4f | replicate %d/%d | seed %d",
                    embedding.type, p.dim, sigma, rep.idx, n.rep, rep.seed))

    gen <- generate.highdim.circle(
      n = config$n_points,
      p = p.dim,
      radius = config$radius,
      sigma = sigma,
      embedding.type = embedding.type,
      ambient.noise.factor = config$ambient_noise_factor,
      noise.type = config$noise_type,
      seed = rep.seed
    )

    X <- as.matrix(gen$X.obs)
    X.true <- as.matrix(gen$X.true)
    basis <- as.matrix(gen$basis)

    proxy.y <- stats::prcomp(X, center = TRUE, scale. = FALSE, rank. = 1L)$x[, 1L]

    splits <- vector("list", split.count)
    split.ptr <- 0L
    for (repeat.id in seq_len(config$n_repeats)) {
      split.seed <- as.integer(rep.seed + repeat.id * 1000L)
      set.seed(split.seed)
      fold.id <- sample(rep(seq_len(config$n_folds), length.out = nrow(X)))
      for (fold.idx in seq_len(config$n_folds)) {
        split.ptr <- split.ptr + 1L
        test.idx <- which(fold.id == fold.idx)
        train.idx <- which(fold.id != fold.idx)
        splits[[split.ptr]] <- list(
          repeat_id = repeat.id,
          fold_id = fold.idx,
          train.idx = as.integer(train.idx),
          test.idx = as.integer(test.idx),
          col.scale = safe.col.scales(X[train.idx, , drop = FALSE])
        )
      }
    }

    graph.start <- Sys.time()
    rep.rows <- vector("list", tasks.per.replicate.graph)
    rep.row.ptr <- 0L
    fit.cache <- list()

    for (pair.idx in seq_len(nrow(pair.grid))) {
      k.val <- as.integer(pair.grid$k[pair.idx])
      eig.val <- as.integer(pair.grid$n_eigenpairs[pair.idx])
      pair.key <- sprintf("%d::%d", k.val, eig.val)

      log.msg(sprintf("  Graph base fit | emb=%s p=%d sigma=%.4f | k=%d | n.eigenpairs=%d",
                      embedding.type, p.dim, sigma, k.val, eig.val))

      fit.fail.status <- "fit_error"
      fit.fail.message <- NA_character_
      model.idx <- seq_len(nrow(X))

      pca.dim.use <- as.integer(min(config$pca_dim, ncol(X), nrow(X) - 1L))
      if (!is.finite(pca.dim.use) || pca.dim.use < 2L) pca.dim.use <- NULL

      fit.base <- tryCatch(
        gflowx::fit.rdgraph.regression(
          X = X,
          y = proxy.y,
          k = k.val,
          pca.dim = pca.dim.use,
          variance.explained = config$variance_explained,
          max.iterations = config$max_iterations,
          n.eigenpairs = eig.val,
          filter.type = config$filter_type,
          t.scale.factor = config$t_scale_factor,
          beta.coef.factor = config$beta_coef_factor,
          response.penalty.exp = config$response_penalty_exp,
          use.counting.measure = config$use_counting_measure,
          max.ratio.threshold = config$max_path_edge_ratio_deviation_thld,
          path.edge.ratio.percentile = config$path_edge_ratio_percentile,
          threshold.percentile = config$threshold_percentile,
          compute.extremality = FALSE,
          verbose.level = 0
        ),
        error = function(e) e
      )

      fit.ok <- !inherits(fit.base, "error")
      if (!fit.ok) {
        fit.fail.message <- conditionMessage(fit.base)
        if (grepl("connected components", fit.fail.message, ignore.case = TRUE)) {
          fit.fail.status <- "graph_disconnected"
        }
        log.msg(sprintf("    Graph base fit failed | k=%d | n.eigenpairs=%d | %s",
                        k.val, eig.val, fit.fail.message))
      } else {
        fit.cache[[pair.key]] <- list(fit.base = fit.base, model.idx = model.idx)
      }

      eta.values <- as.double(config$eta_grid)
      for (eta.val in eta.values) {
        fit.eta <- if (fit.ok) set.model.eta(fit.base, eta.val) else NULL

        for (split in splits) {
          global.task.idx <- global.task.idx + 1L
          rep.row.ptr <- rep.row.ptr + 1L

          status <- "ok"
          err.msg <- NA_character_
          scaled.mse <- NA_real_
          raw.mse <- NA_real_
          n.test.total <- length(split$test.idx)
          n.test.modeled <- if (fit.ok) NA_integer_ else 0L
          n.train.total <- length(split$train.idx)
          n.train.modeled <- if (fit.ok) NA_integer_ else 0L
          modeled.frac <- if (fit.ok) NA_real_ else 0.0

          if (!fit.ok) {
            status <- fit.fail.status
            err.msg <- fit.fail.message
          } else {
            eval.res <- tryCatch(
              evaluate.graph.split(
                model = fit.eta,
                X = X,
                split = split,
                model.idx = model.idx,
                block.size = config$block_size
              ),
              error = function(e) e
            )

            if (inherits(eval.res, "error")) {
              status <- "refit_error"
              err.msg <- conditionMessage(eval.res)
              n.test.modeled <- 0L
              n.train.modeled <- 0L
              modeled.frac <- 0.0
            } else {
              scaled.mse <- eval.res$scaled.mse
              raw.mse <- eval.res$raw.mse
              n.test.total <- as.integer(eval.res$n_test_total)
              n.test.modeled <- as.integer(eval.res$n_test_modeled)
              n.train.total <- as.integer(eval.res$n_train_total)
              n.train.modeled <- as.integer(eval.res$n_train_modeled)
              modeled.frac <- as.double(eval.res$modeled_fraction)
            }
          }

          rep.rows[[rep.row.ptr]] <- data.frame(
            embedding_type = embedding.type,
            p = p.dim,
            sigma = sigma,
            replicate = rep.idx,
            repeat_id = split$repeat_id,
            fold_id = split$fold_id,
            k = k.val,
            n_eigenpairs = eig.val,
            eta = eta.val,
            scaled_mse = scaled.mse,
            raw_mse = raw.mse,
            status = status,
            n_test_total = as.integer(n.test.total),
            n_test_modeled = as.integer(n.test.modeled),
            n_train_total = as.integer(n.train.total),
            n_train_modeled = as.integer(n.train.modeled),
            modeled_fraction = as.double(modeled.frac),
            error = err.msg,
            stringsAsFactors = FALSE
          )

          if ((global.task.idx %% config$log_every_tasks) == 0L || global.task.idx == total.tasks) {
            elapsed <- as.numeric(difftime(Sys.time(), run.start.time, units = "secs"))
            rate <- global.task.idx / max(elapsed, 1e-6)
            eta.sec <- (total.tasks - global.task.idx) / max(rate, 1e-6)
            pct <- 100 * global.task.idx / total.tasks
            log.msg(sprintf("Progress %.2f%% | task %d/%d | rate %.2f tasks/s | ETA %s",
                            pct, global.task.idx, total.tasks, rate, format.seconds(eta.sec)))
          }
        }
      }
    }

    graph.elapsed <- as.numeric(difftime(Sys.time(), graph.start, units = "secs"))

    rep.df <- do.call(rbind, rep.rows[seq_len(rep.row.ptr)])
    append.csv(rep.df, cv.fold.path)

    valid.df <- rep.df[is.finite(rep.df$scaled_mse) & rep.df$status == "ok", , drop = FALSE]
    sel.rows <- NULL

    if (nrow(valid.df) > 0L) {
      key.vec <- sprintf("%d::%d::%.16g", valid.df$k, valid.df$n_eigenpairs, valid.df$eta)
      key.uniq <- unique(key.vec)
      cand.rows <- vector("list", length(key.uniq))

      for (ii in seq_along(key.uniq)) {
        kk <- valid.df[key.vec == key.uniq[ii], , drop = FALSE]
        w <- pmax(kk$n_test_total, 1)
        cand.rows[[ii]] <- data.frame(
          embedding_type = embedding.type,
          p = p.dim,
          sigma = sigma,
          replicate = rep.idx,
          k = as.integer(kk$k[1L]),
          n_eigenpairs = as.integer(kk$n_eigenpairs[1L]),
          eta = as.double(kk$eta[1L]),
          scaled_mse_mean = weighted.mean.safe(kk$scaled_mse, w),
          scaled_mse_se = weighted.se.safe(kk$scaled_mse, w),
          raw_mse_mean = weighted.mean.safe(kk$raw_mse, w),
          modeled_fraction_mean = weighted.mean.safe(kk$modeled_fraction, w),
          n_eval_rows = nrow(kk),
          stringsAsFactors = FALSE
        )
      }

      cand.df <- do.call(rbind, cand.rows)
      cand.df <- cand.df[order(cand.df$scaled_mse_mean, cand.df$raw_mse_mean), ]
      append.csv(cand.df, cv.summary.path)

      best.min <- cand.df[1L, ]
      threshold.1se <- best.min$scaled_mse_mean + best.min$scaled_mse_se
      cand.1se <- cand.df[cand.df$scaled_mse_mean <= threshold.1se, , drop = FALSE]
      cand.1se <- cand.1se[order(cand.1se$n_eigenpairs, -cand.1se$k, -cand.1se$eta), , drop = FALSE]
      best.1se <- cand.1se[1L, ]

      sel.rows <- rbind(
        data.frame(
          embedding_type = embedding.type,
          p = p.dim,
          sigma = sigma,
          replicate = rep.idx,
          selected = "min",
          k = best.min$k,
          n_eigenpairs = best.min$n_eigenpairs,
          eta = best.min$eta,
          scaled_mse_mean = best.min$scaled_mse_mean,
          scaled_mse_se = best.min$scaled_mse_se,
          raw_mse_mean = best.min$raw_mse_mean,
          modeled_fraction_mean = best.min$modeled_fraction_mean,
          stringsAsFactors = FALSE
        ),
        data.frame(
          embedding_type = embedding.type,
          p = p.dim,
          sigma = sigma,
          replicate = rep.idx,
          selected = "one_se",
          k = best.1se$k,
          n_eigenpairs = best.1se$n_eigenpairs,
          eta = best.1se$eta,
          scaled_mse_mean = best.1se$scaled_mse_mean,
          scaled_mse_se = best.1se$scaled_mse_se,
          raw_mse_mean = best.1se$raw_mse_mean,
          modeled_fraction_mean = best.1se$modeled_fraction_mean,
          stringsAsFactors = FALSE
        )
      )
      append.csv(sel.rows, selection.graph.path)
    } else {
      log.msg(sprintf("  No valid graph CV rows for emb=%s p=%d sigma=%.4f rep=%d",
                      embedding.type, p.dim, sigma, rep.idx))
    }

    if (!is.null(sel.rows) && nrow(sel.rows) > 0L) {
      for (sel.idx in seq_len(nrow(sel.rows))) {
        k.sel <- as.integer(sel.rows$k[sel.idx])
        eig.sel <- as.integer(sel.rows$n_eigenpairs[sel.idx])
        eta.sel <- as.double(sel.rows$eta[sel.idx])
        pair.key <- sprintf("%d::%d", k.sel, eig.sel)
        fit.rec <- fit.cache[[pair.key]]

        if (is.null(fit.rec) || is.null(fit.rec$fit.base)) {
          log.msg(sprintf("  Selected graph candidate missing in cache | emb=%s p=%d sigma=%.4f rep=%d sel=%s",
                          embedding.type, p.dim, sigma, rep.idx, sel.rows$selected[sel.idx]))
          next
        }

        fit.sel <- set.model.eta(fit.rec$fit.base, eta.sel)
        refit.full <- tryCatch(
          gflowx::refit.rdgraph.regression(
            fitted.model = fit.sel,
            y.new = X,
            per.column.gcv = FALSE,
            verbose = FALSE
          ),
          error = function(e) e
        )

        if (inherits(refit.full, "error")) {
          log.msg(sprintf("  Graph oracle refit failed | emb=%s p=%d sigma=%.4f rep=%d sel=%s | %s",
                          embedding.type, p.dim, sigma, rep.idx, sel.rows$selected[sel.idx], conditionMessage(refit.full)))
          next
        }

        X.hat.full <- as.matrix(refit.full$fitted.values)
        om <- compute.oracle.metrics(
          X.obs = X,
          X.hat = X.hat.full,
          X.true = X.true,
          basis = basis,
          method = paste0("cv_", sel.rows$selected[sel.idx]),
          embedding.type = embedding.type,
          p = p.dim,
          sigma = sigma,
          replicate = rep.idx,
          notes = sprintf("k=%d,n.eigenpairs=%d,eta=%.8f", k.sel, eig.sel, eta.sel),
          modeled_fraction = 1.0
        )
        append.csv(om, oracle.path)
      }
    }

    classical.start <- Sys.time()
    class.res <- evaluate.classical.cv(
      X = X,
      splits = splits,
      classical.grid = list(
        knn_k_grid = config$classical_knn_k_grid,
        gaussian_bw_mult_grid = config$classical_gaussian_bw_mult_grid
      )
    )
    class.elapsed <- as.numeric(difftime(Sys.time(), classical.start, units = "secs"))

    if (nrow(class.res$candidates) > 0L) {
      class.cand <- class.res$candidates
      class.cand$embedding_type <- embedding.type
      class.cand$p <- p.dim
      class.cand$sigma <- sigma
      class.cand$replicate <- rep.idx
      class.cand <- class.cand[, c("embedding_type", "p", "sigma", "replicate",
                                   "method", "param_name", "param_value",
                                   "scaled_mse_mean", "scaled_mse_se", "raw_mse_mean", "n_eval_rows")]
      append.csv(class.cand, cv.classical.path)
    }

    if (nrow(class.res$selection) > 0L) {
      class.sel <- class.res$selection
      class.sel$embedding_type <- embedding.type
      class.sel$p <- p.dim
      class.sel$sigma <- sigma
      class.sel$replicate <- rep.idx
      class.sel <- class.sel[, c("embedding_type", "p", "sigma", "replicate",
                                 "method", "selected", "param_name", "param_value",
                                 "scaled_mse_mean", "scaled_mse_se", "raw_mse_mean")]
      append.csv(class.sel, selection.classical.path)

      for (jj in seq_len(nrow(class.sel))) {
        method.id <- as.character(class.sel$method[jj])
        param.value <- as.double(class.sel$param_value[jj])
        X.hat <- tryCatch(
          predict.classical.full(X = X, method.id = method.id, param.value = param.value),
          error = function(e) e
        )

        if (inherits(X.hat, "error")) {
          log.msg(sprintf("  Classical oracle prediction failed | emb=%s p=%d sigma=%.4f rep=%d method=%s | %s",
                          embedding.type, p.dim, sigma, rep.idx, method.id, conditionMessage(X.hat)))
          next
        }

        om <- compute.oracle.metrics(
          X.obs = X,
          X.hat = X.hat,
          X.true = X.true,
          basis = basis,
          method = method.id,
          embedding.type = embedding.type,
          p = p.dim,
          sigma = sigma,
          replicate = rep.idx,
          notes = sprintf("%s=%.6f", class.sel$param_name[jj], param.value),
          modeled_fraction = 1.0
        )
        append.csv(om, oracle.path)
      }
    }

    baseline.obs <- compute.oracle.metrics(
      X.obs = X,
      X.hat = X,
      X.true = X.true,
      basis = basis,
      method = "observed_no_smoothing",
      embedding.type = embedding.type,
      p = p.dim,
      sigma = sigma,
      replicate = rep.idx,
      notes = "identity baseline",
      modeled_fraction = 0.0
    )
    append.csv(baseline.obs, oracle.path)

    runtime.rows <- rbind(
      data.frame(
        embedding_type = embedding.type,
        p = p.dim,
        sigma = sigma,
        replicate = rep.idx,
        method_group = "graph_cv_total",
        elapsed_sec = graph.elapsed,
        stringsAsFactors = FALSE
      ),
      data.frame(
        embedding_type = embedding.type,
        p = p.dim,
        sigma = sigma,
        replicate = rep.idx,
        method_group = "classical_cv_total",
        elapsed_sec = class.elapsed,
        stringsAsFactors = FALSE
      )
    )
    append.csv(runtime.rows, runtime.path)

    ckpt.path <- file.path(
      run.dir, "checkpoints",
      sprintf("checkpoint_%s_p_%d_sigma_%s_rep_%03d.rds",
              embedding.type,
              p.dim,
              formatC(sigma, digits = 3, format = "f"),
              rep.idx)
    )

    saveRDS(
      list(
        embedding_type = embedding.type,
        p = p.dim,
        sigma = sigma,
        replicate = rep.idx,
        seed = rep.seed,
        config = config,
        graph_selection = sel.rows,
        classical_selection = class.res$selection
      ),
      ckpt.path
    )

    elapsed.rep <- as.numeric(difftime(Sys.time(), run.start.time, units = "secs"))
    log.msg(sprintf("Completed emb=%s p=%d sigma=%.4f rep=%d | elapsed %s",
                    embedding.type, p.dim, sigma, rep.idx, format.seconds(elapsed.rep)))
  }
}

selection.graph.df <- if (file.exists(selection.graph.path)) utils::read.csv(selection.graph.path, stringsAsFactors = FALSE) else data.frame()
if (nrow(selection.graph.df) > 0L) {
  out.rows <- list()
  out.ptr <- 0L
  for (emb in unique(selection.graph.df$embedding_type)) {
    for (pp in sort(unique(selection.graph.df$p))) {
      for (s in sort(unique(selection.graph.df$sigma))) {
        for (sel in c("min", "one_se")) {
          sub <- selection.graph.df[
            selection.graph.df$embedding_type == emb &
              selection.graph.df$p == pp &
              selection.graph.df$sigma == s &
              selection.graph.df$selected == sel,
            , drop = FALSE
          ]
          if (nrow(sub) == 0L) next
          out.ptr <- out.ptr + 1L
          out.rows[[out.ptr]] <- data.frame(
            embedding_type = emb,
            p = pp,
            sigma = s,
            selected = sel,
            k_mode = mode.param(sub$k),
            n_eigenpairs_mode = mode.param(sub$n_eigenpairs),
            eta_mode = mode.param(sub$eta),
            scaled_mse_mean_avg = mean(sub$scaled_mse_mean, na.rm = TRUE),
            scaled_mse_mean_sd = stats::sd(sub$scaled_mse_mean, na.rm = TRUE),
            n_replicates = nrow(sub),
            stringsAsFactors = FALSE
          )
        }
      }
    }
  }
  selected.by.group <- do.call(rbind, out.rows)
  write.csv(selected.by.group,
            file = file.path(run.dir, "tables", "graph_selected_params_by_group.csv"),
            row.names = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(run.dir, "artifacts", "session_info.txt"))

elapsed.total <- as.numeric(difftime(Sys.time(), run.start.time, units = "secs"))
log.msg("Run complete.")
log.msg("Total elapsed:", format.seconds(elapsed.total))
log.msg("Outputs written under:", run.dir)
