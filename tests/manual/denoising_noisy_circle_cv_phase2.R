#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x

find.repo.root <- function(start = getwd()) {
  cur <- normalizePath(start, winslash = "/", mustWork = TRUE)
  repeat {
    if (file.exists(file.path(cur, "DESCRIPTION"))) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) {
      stop("Could not locate repository root (missing DESCRIPTION).", call. = FALSE)
    }
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

draw.radial.noise <- function(n, noise, noise.type) {
  if (!is.finite(noise) || noise <= 0) return(rep(0.0, n))
  if (identical(noise.type, "laplace")) {
    rlaplace(n, location = 0, scale = noise / sqrt(2))
  } else {
    stats::rnorm(n, mean = 0, sd = noise)
  }
}

generate.circle.nonuniform <- function(n, radius, noise, noise.type = "normal", seed = NULL,
                                       arc1.weight = 0.45, arc2.weight = 0.35,
                                       arc.halfwidth = pi / 6) {
  if (!is.null(seed)) set.seed(seed)
  u <- stats::runif(n)
  angles <- numeric(n)

  idx.arc1 <- u <= arc1.weight
  idx.arc2 <- u > arc1.weight & u <= (arc1.weight + arc2.weight)
  idx.bg <- !(idx.arc1 | idx.arc2)

  angles[idx.arc1] <- stats::runif(sum(idx.arc1), -arc.halfwidth, arc.halfwidth)
  angles[idx.arc2] <- pi + stats::runif(sum(idx.arc2), -arc.halfwidth, arc.halfwidth)
  angles[idx.bg] <- stats::runif(sum(idx.bg), 0, 2 * pi)
  angles <- (angles %% (2 * pi))

  eps <- draw.radial.noise(n, noise = noise, noise.type = noise.type)
  radial <- pmax(radius + eps, 1e-6)

  x.noise.free <- radius * cos(angles)
  y.noise.free <- radius * sin(angles)

  data.frame(
    x = radial * cos(angles),
    y = radial * sin(angles),
    x.noise.free = x.noise.free,
    y.noise.free = y.noise.free,
    angles = angles,
    is_outlier = FALSE,
    stringsAsFactors = FALSE
  )
}

generate.circle.outliers <- function(n, radius, noise, noise.type = "normal", seed = NULL,
                                     outlier.frac = 0.05,
                                     outlier.radial.shift = 0.8,
                                     outlier.radial.sd = 0.25) {
  if (!is.null(seed)) set.seed(seed)
  angles <- stats::runif(n, 0, 2 * pi)
  eps <- draw.radial.noise(n, noise = noise, noise.type = noise.type)
  radial <- pmax(radius + eps, 1e-6)

  n.out <- max(1L, as.integer(round(outlier.frac * n)))
  out.idx <- sort(sample.int(n, size = n.out, replace = FALSE))
  bump <- outlier.radial.shift + abs(stats::rnorm(n.out, mean = 0, sd = outlier.radial.sd))
  signs <- sample(c(-1, 1), size = n.out, replace = TRUE)
  radial[out.idx] <- pmax(radial[out.idx] + signs * bump, 1e-6)

  x.noise.free <- radius * cos(angles)
  y.noise.free <- radius * sin(angles)
  is.outlier <- rep(FALSE, n)
  is.outlier[out.idx] <- TRUE

  data.frame(
    x = radial * cos(angles),
    y = radial * sin(angles),
    x.noise.free = x.noise.free,
    y.noise.free = y.noise.free,
    angles = angles,
    is_outlier = is.outlier,
    stringsAsFactors = FALSE
  )
}

generate.dgp.data <- function(dgp.id, config, sigma, seed) {
  if (identical(dgp.id, "circle_random")) {
    out <- generate.circle.data(
      n = config$n_points,
      radius = config$radius,
      noise = sigma,
      type = "random",
      noise.type = config$noise_type,
      seed = seed
    )
    out$is_outlier <- FALSE
    return(out)
  }

  if (identical(dgp.id, "circle_nonuniform")) {
    return(generate.circle.nonuniform(
      n = config$n_points,
      radius = config$radius,
      noise = sigma,
      noise.type = config$noise_type,
      seed = seed,
      arc1.weight = config$nonuniform_arc1_weight,
      arc2.weight = config$nonuniform_arc2_weight,
      arc.halfwidth = config$nonuniform_arc_halfwidth
    ))
  }

  if (identical(dgp.id, "circle_outliers")) {
    return(generate.circle.outliers(
      n = config$n_points,
      radius = config$radius,
      noise = sigma,
      noise.type = config$noise_type,
      seed = seed,
      outlier.frac = config$outlier_fraction,
      outlier.radial.shift = config$outlier_radial_shift,
      outlier.radial.sd = config$outlier_radial_sd
    ))
  }

  stop("Unsupported DGP id: ", dgp.id, call. = FALSE)
}

build.connectivity.profile <- function(X, config) {
  k.min <- min(config$k_grid)
  k.max <- max(config$k_grid)

  gtmp <- create.iknn.graphs(
    X = X,
    kmin = k.min,
    kmax = k.max,
    max.path.edge.ratio.deviation.thld = config$max_path_edge_ratio_deviation_thld,
    path.edge.ratio.percentile = config$path_edge_ratio_percentile,
    threshold.percentile = config$threshold_percentile,
    compute.full = TRUE,
    pca.dim = config$pca_dim,
    variance.explained = config$variance_explained,
    n.cores = 1L,
    verbose = FALSE
  )

  st <- NULL
  invisible(capture.output({
    st <- summary(gtmp)
  }))

  if (is.null(st) || !all(c("k", "n_ccomp") %in% names(st))) {
    stop("Connectivity pre-screen failed: missing k/n_ccomp in summary(create.iknn.graphs()).")
  }

  k.vec <- as.integer(st$k)
  n.cc <- as.integer(st$n_ccomp)
  g.list <- gtmp$geom_pruned_graphs

  if (is.null(g.list) || length(g.list) != length(k.vec)) {
    stop("Connectivity pre-screen failed: geom_pruned_graphs missing or length mismatch.")
  }

  profile <- vector("list", length(k.vec))
  names(profile) <- as.character(k.vec)

  for (i in seq_along(k.vec)) {
    g.obj <- g.list[[i]]
    adj <- g.obj$adj_list
    if (is.null(adj)) {
      keep.idx <- seq_len(nrow(X))
      n.comp <- 1L
    } else {
      comp <- graph.connected.components(adj)
      tab <- table(comp)
      keep.label <- as.integer(names(tab)[which.max(tab)])
      keep.idx <- which(comp == keep.label)
      n.comp <- length(tab)
    }

    profile[[i]] <- list(
      k = as.integer(k.vec[i]),
      n_components = as.integer(n.comp),
      lcc_size = as.integer(length(keep.idx)),
      lcc_frac = as.double(length(keep.idx) / nrow(X)),
      keep_idx = as.integer(keep.idx)
    )
  }

  names(profile) <- as.character(k.vec)

  list(
    profile = profile,
    connected_k = as.integer(k.vec[n.cc == 1L]),
    n_ccomp_by_k = stats::setNames(n.cc, as.character(k.vec)),
    lcc_frac_by_k = stats::setNames(vapply(profile, `[[`, numeric(1), "lcc_frac"), as.character(k.vec)),
    lcc_size_by_k = stats::setNames(vapply(profile, `[[`, integer(1), "lcc_size"), as.character(k.vec))
  )
}

evaluate.split <- function(model, X, split, model.idx, block.size = NULL) {
  model.idx <- as.integer(model.idx)
  X.model <- X[model.idx, , drop = FALSE]

  train.local <- match(split$train.idx, model.idx, nomatch = 0L)
  train.local <- as.integer(train.local[train.local > 0L])
  if (length(train.local) < 2L) {
    stop("Too few training rows inside modeled component.")
  }

  test.local.all <- match(split$test.idx, model.idx, nomatch = 0L)
  modeled.mask <- as.logical(test.local.all > 0L)
  test.local <- as.integer(test.local.all[modeled.mask])
  if (length(test.local) < 1L) {
    stop("No test rows inside modeled component.")
  }

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

compute.oracle.metrics <- function(X.obs, X.hat, X.true, method, dgp.id, sigma, replicate,
                                   notes = NA_character_, modeled_fraction = NA_real_) {
  X.obs <- as.matrix(X.obs)
  X.hat <- as.matrix(X.hat)
  X.true <- as.matrix(X.true)
  rad.obs <- sqrt(rowSums(X.obs^2))
  rad.hat <- sqrt(rowSums(X.hat^2))
  rad.true <- sqrt(rowSums(X.true^2))
  data.frame(
    dgp_id = dgp.id,
    sigma = as.double(sigma),
    replicate = as.integer(replicate),
    method = method,
    n_rows = nrow(X.hat),
    modeled_fraction = modeled_fraction,
    rmse_xy_obs = sqrt(mean((X.obs - X.true)^2)),
    rmse_xy_hat = sqrt(mean((X.hat - X.true)^2)),
    rmse_rad_obs = sqrt(mean((rad.obs - rad.true)^2)),
    rmse_rad_hat = sqrt(mean((rad.hat - rad.true)^2)),
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
    dgp_id = config$dgp_grid,
    sigma = as.double(config$sigma_grid),
    KEEP.OUT.ATTRS = FALSE,
    stringsAsFactors = FALSE
  )
  sg$n_replicates <- ifelse(
    abs(sg$sigma - config$edge_sigma) <= 1e-12,
    config$n_replicates_edge,
    config$n_replicates_main
  )
  sg <- sg[order(match(sg$dgp_id, config$dgp_grid), sg$sigma), ]
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
  experiment_name = "noisy_circle_cv_phase2",
  mode = mode,
  base_seed = base.seed,
  n_points = 250L,
  radius = 1.0,
  noise_type = "normal",
  dgp_grid = c("circle_random", "circle_nonuniform", "circle_outliers"),
  sigma_grid = c(0.02, 0.05, 0.10, 0.20, 0.30),
  edge_sigma = 0.02,
  n_replicates_edge = 8L,
  n_replicates_main = 24L,
  n_folds = 5L,
  n_repeats = 5L,
  k_grid = seq.int(5L, 35L, by = 2L),
  n_eigenpairs_grid = c(25L, 50L, 100L),
  eta_grid = exp(seq(log(1e-3), log(2.0), length.out = 12L)),
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
  block_size = NULL,
  use_connectivity_prescreen = TRUE,
  connectivity_policy = "largest_component",
  lcc_min_fraction = 0.95,
  lcc_min_train_vertices = 20L,
  lcc_min_total_vertices = 80L,
  outlier_fraction = 0.05,
  outlier_radial_shift = 0.8,
  outlier_radial_sd = 0.25,
  nonuniform_arc1_weight = 0.45,
  nonuniform_arc2_weight = 0.35,
  nonuniform_arc_halfwidth = pi / 6,
  log_every_tasks = 100L,
  run_data_smoother_baseline = TRUE
)

if (mode == "pilot") {
  config$n_replicates_edge <- 1L
  config$n_replicates_main <- 2L
  config$n_repeats <- 2L
  config$n_folds <- 3L
}

if (!is.null(args$dgp_grid)) config$dgp_grid <- parse.char.vec(args$dgp_grid, config$dgp_grid)
if (!is.null(args$sigma_grid)) config$sigma_grid <- as.numeric(parse.char.vec(args$sigma_grid, config$sigma_grid))
if (!is.null(args$n_replicates_edge)) config$n_replicates_edge <- as.int.arg(args$n_replicates_edge, config$n_replicates_edge)
if (!is.null(args$n_replicates_main)) config$n_replicates_main <- as.int.arg(args$n_replicates_main, config$n_replicates_main)
if (!is.null(args$n_repeats)) config$n_repeats <- as.int.arg(args$n_repeats, config$n_repeats)
if (!is.null(args$n_folds)) config$n_folds <- as.int.arg(args$n_folds, config$n_folds)
if (!is.null(args$lcc_min_fraction)) config$lcc_min_fraction <- as.num.arg(args$lcc_min_fraction, config$lcc_min_fraction)
if (!is.null(args$lcc_min_train_vertices)) config$lcc_min_train_vertices <- as.int.arg(args$lcc_min_train_vertices, config$lcc_min_train_vertices)
if (!is.null(args$lcc_min_total_vertices)) config$lcc_min_total_vertices <- as.int.arg(args$lcc_min_total_vertices, config$lcc_min_total_vertices)
if (!is.null(args$log_every_tasks)) config$log_every_tasks <- as.int.arg(args$log_every_tasks, config$log_every_tasks)

scenario.grid <- build.scenario.grid(config)
config$scenario_grid <- scenario.grid

run.id <- args$run_id %||% format(Sys.time(), "%Y%m%d_%H%M%S")
run.dir <- file.path(repo.root, "tests", "manual", "results", "noisy_circle_cv_phase2", run.id)
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
tasks.per.replicate <- nrow(candidate.grid) * split.count
total.replicates <- sum(as.integer(scenario.grid$n_replicates))
total.tasks <- total.replicates * tasks.per.replicate

log.msg("Starting experiment run.")
log.msg("Run directory:", run.dir)
log.msg("Total tasks:", total.tasks,
        "| DGP-sigma scenarios:", nrow(scenario.grid),
        "| total replicates:", total.replicates,
        "| candidates:", nrow(candidate.grid),
        "| splits/replicate:", split.count)

cv.fold.path <- file.path(run.dir, "tables", "cv_fold_scores.csv")
cv.summary.path <- file.path(run.dir, "tables", "cv_candidate_summary.csv")
selection.path <- file.path(run.dir, "tables", "selection_by_replicate.csv")
oracle.path <- file.path(run.dir, "tables", "oracle_metrics_by_replicate.csv")

global.task.idx <- 0L
run.start.time <- Sys.time()

for (scenario.idx in seq_len(nrow(scenario.grid))) {
  dgp.id <- as.character(scenario.grid$dgp_id[scenario.idx])
  sigma <- as.double(scenario.grid$sigma[scenario.idx])
  n.rep <- as.integer(scenario.grid$n_replicates[scenario.idx])

  for (rep.idx in seq_len(n.rep)) {
    rep.seed <- as.integer(base.seed + scenario.idx * 100000L + rep.idx)
    set.seed(rep.seed)

    log.msg(sprintf("DGP %s | sigma %.4f | replicate %d/%d | seed %d",
                    dgp.id, sigma, rep.idx, n.rep, rep.seed))

    X.df <- generate.dgp.data(
      dgp.id = dgp.id,
      config = config,
      sigma = sigma,
      seed = rep.seed
    )
    X <- as.matrix(X.df[, c("x", "y"), drop = FALSE])
    X.true <- as.matrix(X.df[, c("x.noise.free", "y.noise.free"), drop = FALSE])

    connected.k <- NULL
    n.ccomp.by.k <- NULL
    lcc.profile <- NULL

    if (isTRUE(config$use_connectivity_prescreen)) {
      pre.res <- tryCatch(
        build.connectivity.profile(X = X, config = config),
        error = function(e) e
      )

      if (inherits(pre.res, "error")) {
        log.msg(sprintf("  Connectivity pre-screen failed; falling back to direct fit checks | error=%s",
                        conditionMessage(pre.res)))
      } else {
        connected.k <- pre.res$connected_k
        n.ccomp.by.k <- pre.res$n_ccomp_by_k
        lcc.profile <- pre.res$profile
        disconnected.k <- setdiff(as.integer(config$k_grid), connected.k)
        if (length(disconnected.k) > 0L) {
          lcc.msg <- vapply(disconnected.k, function(kv) {
            p <- lcc.profile[[as.character(kv)]]
            sprintf("k=%d(n_ccomp=%d,lcc=%.3f)", kv, p$n_components, p$lcc_frac)
          }, FUN.VALUE = character(1))
          log.msg(sprintf("  Connectivity pre-screen: disconnected k values detected: %s",
                          paste(lcc.msg, collapse = ", ")))
        } else {
          log.msg("  Connectivity pre-screen: all k values connected.")
        }
      }
    }

    proxy.y <- stats::prcomp(X, center = TRUE, scale. = FALSE)$x[, 1L]

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

    rep.rows <- vector("list", tasks.per.replicate)
    rep.row.ptr <- 0L
    fit.cache <- list()

    for (pair.idx in seq_len(nrow(pair.grid))) {
      k.val <- as.integer(pair.grid$k[pair.idx])
      eig.val <- as.integer(pair.grid$n_eigenpairs[pair.idx])
      pair.key <- sprintf("%d::%d", k.val, eig.val)

      log.msg(sprintf("  Fitting base model | dgp=%s | sigma=%.4f | k=%d | n.eigenpairs=%d",
                      dgp.id, sigma, k.val, eig.val))

      fit.fail.status <- "fit_error"
      fit.fail.message <- NA_character_
      graph.policy <- "full_graph"
      model.idx <- seq_len(nrow(X))
      n.components <- 1L
      lcc.size <- nrow(X)
      lcc.frac <- 1.0

      if (!is.null(n.ccomp.by.k)) {
        n.components <- as.integer(n.ccomp.by.k[[as.character(k.val)]] %||% NA_integer_)
      }
      if (!is.null(lcc.profile) && !is.null(lcc.profile[[as.character(k.val)]])) {
        p <- lcc.profile[[as.character(k.val)]]
        lcc.size <- as.integer(p$lcc_size)
        lcc.frac <- as.double(p$lcc_frac)
      }

      if (!is.null(connected.k) && !(k.val %in% connected.k)) {
        graph.policy <- "graph_disconnected"
        if (identical(config$connectivity_policy, "largest_component") &&
            is.finite(lcc.frac) &&
            lcc.frac >= config$lcc_min_fraction &&
            lcc.size >= config$lcc_min_total_vertices) {
          model.idx <- as.integer(lcc.profile[[as.character(k.val)]]$keep_idx)
          graph.policy <- "largest_component"
          log.msg(sprintf("    Disconnected graph accepted by LCC policy | k=%d | n.ccomp=%d | lcc.frac=%.3f | n.model=%d",
                          k.val, n.components, lcc.frac, length(model.idx)))
        } else {
          fit.base <- simpleError(sprintf(
            "Pre-screen: graph at k=%d disconnected (n_ccomp=%d, lcc_frac=%.3f). Base fit skipped by policy.",
            k.val,
            as.integer(n.components),
            as.double(lcc.frac)
          ))
          fit.ok <- FALSE
          fit.fail.status <- "graph_disconnected"
          fit.fail.message <- conditionMessage(fit.base)
          log.msg(sprintf("    Base fit skipped | k=%d | n.eigenpairs=%d | error=%s",
                          k.val, eig.val, fit.fail.message))
        }
      }

      if (!exists("fit.ok", inherits = FALSE)) {
        fit.X <- X[model.idx, , drop = FALSE]
        fit.y <- proxy.y[model.idx]

        if (length(model.idx) < config$lcc_min_train_vertices) {
          fit.base <- simpleError(sprintf(
            "Modeled component too small for stable fit: n.model=%d < lcc_min_train_vertices=%d",
            length(model.idx), config$lcc_min_train_vertices
          ))
          fit.ok <- FALSE
          fit.fail.status <- "graph_disconnected_small_lcc"
          fit.fail.message <- conditionMessage(fit.base)
          log.msg(sprintf("    Base fit skipped | k=%d | n.eigenpairs=%d | error=%s",
                          k.val, eig.val, fit.fail.message))
        } else {
          fit.base <- tryCatch(
            gflowx::fit.rdgraph.regression(
              X = fit.X,
              y = fit.y,
              k = k.val,
              pca.dim = config$pca_dim,
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
            log.msg(sprintf("    Base fit failed | k=%d | n.eigenpairs=%d | error=%s",
                            k.val, eig.val, fit.fail.message))
          } else {
            fit.cache[[pair.key]] <- list(
              fit.base = fit.base,
              model.idx = as.integer(model.idx),
              graph.policy = graph.policy,
              n_components = as.integer(n.components),
              lcc_size = as.integer(lcc.size),
              lcc_frac = as.double(lcc.frac)
            )
          }
        }
      }

      eta.values <- as.double(config$eta_grid)
      for (eta.val in eta.values) {
        if (fit.ok) {
          fit.eta <- set.model.eta(fit.base, eta.val)
        } else {
          fit.eta <- NULL
        }

        for (split in splits) {
          global.task.idx <- global.task.idx + 1L
          rep.row.ptr <- rep.row.ptr + 1L

          status <- "ok"
          err.msg <- NA_character_
          scaled.mse <- NA_real_
          raw.mse <- NA_real_
          eval.time <- Sys.time()
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
              evaluate.split(
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
            dgp_id = dgp.id,
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
            graph_policy = graph.policy,
            n_components = as.integer(n.components),
            lcc_size = as.integer(lcc.size),
            lcc_frac = as.double(lcc.frac),
            n_rows_model = as.integer(length(model.idx)),
            n_test_total = as.integer(n.test.total),
            n_test_modeled = as.integer(n.test.modeled),
            n_train_total = as.integer(n.train.total),
            n_train_modeled = as.integer(n.train.modeled),
            modeled_fraction = as.double(modeled.frac),
            error = err.msg,
            eval_started = format(eval.time, "%Y-%m-%d %H:%M:%S"),
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

      if (exists("fit.ok", inherits = FALSE)) rm(fit.ok)
    }

    rep.df <- do.call(rbind, rep.rows[seq_len(rep.row.ptr)])
    append.csv(rep.df, cv.fold.path)

    valid.df <- rep.df[is.finite(rep.df$scaled_mse) & rep.df$status == "ok", , drop = FALSE]

    if (nrow(valid.df) == 0L) {
      log.msg(sprintf("  No valid CV rows for dgp=%s sigma=%.4f replicate=%d; skipping selections.", dgp.id, sigma, rep.idx))
      sel.row <- data.frame(
        dgp_id = dgp.id,
        sigma = sigma,
        replicate = rep.idx,
        selected = "none",
        k = NA_integer_,
        n_eigenpairs = NA_integer_,
        eta = NA_real_,
        scaled_mse_mean = NA_real_,
        scaled_mse_se = NA_real_,
        raw_mse_mean = NA_real_,
        modeled_fraction_mean = NA_real_,
        stringsAsFactors = FALSE
      )
      append.csv(sel.row, selection.path)
      next
    }

    key.vec <- sprintf("%d::%d::%.16g", valid.df$k, valid.df$n_eigenpairs, valid.df$eta)
    key.uniq <- unique(key.vec)
    cand.rows <- vector("list", length(key.uniq))

    for (ii in seq_along(key.uniq)) {
      kk <- valid.df[key.vec == key.uniq[ii], , drop = FALSE]
      w <- pmax(kk$n_test_total, 1)
      cand.rows[[ii]] <- data.frame(
        dgp_id = dgp.id,
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

    append.csv(cand.df[, c("dgp_id", "sigma", "replicate", "k", "n_eigenpairs", "eta",
                           "scaled_mse_mean", "scaled_mse_se", "raw_mse_mean",
                           "modeled_fraction_mean", "n_eval_rows")], cv.summary.path)

    best.min <- cand.df[1L, ]
    threshold.1se <- best.min$scaled_mse_mean + best.min$scaled_mse_se
    cand.1se <- cand.df[cand.df$scaled_mse_mean <= threshold.1se, , drop = FALSE]
    cand.1se <- cand.1se[order(cand.1se$n_eigenpairs, -cand.1se$k, -cand.1se$eta), ]
    best.1se <- cand.1se[1L, ]

    sel.rows <- rbind(
      data.frame(
        dgp_id = dgp.id,
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
        dgp_id = dgp.id,
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
    append.csv(sel.rows, selection.path)

    for (sel.idx in seq_len(nrow(sel.rows))) {
      k.sel <- as.integer(sel.rows$k[sel.idx])
      eig.sel <- as.integer(sel.rows$n_eigenpairs[sel.idx])
      eta.sel <- as.double(sel.rows$eta[sel.idx])
      pair.key <- sprintf("%d::%d", k.sel, eig.sel)
      fit.rec <- fit.cache[[pair.key]]

      if (is.null(fit.rec) || is.null(fit.rec$fit.base)) {
        log.msg(sprintf("  Selected candidate not available in fit cache | dgp=%s sigma=%.4f rep=%d sel=%s",
                        dgp.id, sigma, rep.idx, sel.rows$selected[sel.idx]))
        next
      }

      fit.sel <- set.model.eta(fit.rec$fit.base, eta.sel)
      X.model <- X[fit.rec$model.idx, , drop = FALSE]
      refit.full <- tryCatch(
        gflowx::refit.rdgraph.regression(
          fitted.model = fit.sel,
          y.new = X.model,
          per.column.gcv = FALSE,
          verbose = FALSE
        ),
        error = function(e) e
      )

      if (inherits(refit.full, "error")) {
        log.msg(sprintf("  Oracle refit failed | dgp=%s sigma=%.4f rep=%d sel=%s | %s",
                        dgp.id, sigma, rep.idx, sel.rows$selected[sel.idx], conditionMessage(refit.full)))
        next
      }

      X.hat.full <- X
      X.hat.full[fit.rec$model.idx, ] <- as.matrix(refit.full$fitted.values)
      om <- compute.oracle.metrics(
        X.obs = X,
        X.hat = X.hat.full,
        X.true = X.true,
        method = paste0("cv_", sel.rows$selected[sel.idx]),
        dgp.id = dgp.id,
        sigma = sigma,
        replicate = rep.idx,
        notes = sprintf("k=%d,n.eigenpairs=%d,eta=%.8f,graph.policy=%s,lcc.frac=%.4f",
                        k.sel, eig.sel, eta.sel, fit.rec$graph.policy, fit.rec$lcc_frac),
        modeled_fraction = length(fit.rec$model.idx) / nrow(X)
      )
      append.csv(om, oracle.path)
    }

    baseline.obs <- compute.oracle.metrics(
      X.obs = X,
      X.hat = X,
      X.true = X.true,
      method = "observed_no_smoothing",
      dgp.id = dgp.id,
      sigma = sigma,
      replicate = rep.idx,
      notes = "identity baseline",
      modeled_fraction = 0.0
    )
    append.csv(baseline.obs, oracle.path)

    if (isTRUE(config$run_data_smoother_baseline)) {
      ds <- tryCatch(
        {
          ds.out <- NULL
          invisible(capture.output({
            ds.out <- data.smoother(
              X = X,
              kmin = min(config$k_grid),
              kmax = max(config$k_grid),
              proxy.response = "pc1",
              max.iterations = config$max_iterations,
              n.eigenpairs = 50L,
              filter.type = config$filter_type,
              t.scale.factor = config$t_scale_factor,
              beta.coef.factor = config$beta_coef_factor,
              use.counting.measure = config$use_counting_measure,
              max.path.edge.ratio.deviation.thld = config$max_path_edge_ratio_deviation_thld,
              path.edge.ratio.percentile = config$path_edge_ratio_percentile,
              threshold.percentile = config$threshold_percentile,
              n.cores = 1L,
              verbose = FALSE
            )
          }))
          ds.out
        },
        error = function(e) e
      )

      if (inherits(ds, "error")) {
        log.msg(sprintf("  data.smoother baseline failed | dgp=%s sigma=%.4f rep=%d | %s",
                        dgp.id, sigma, rep.idx, conditionMessage(ds)))
      } else {
        X.hat.full <- X
        keep.idx <- if (isTRUE(ds$trimmed)) as.integer(ds$kept.rows) else seq_len(nrow(X))
        X.hat.full[keep.idx, ] <- as.matrix(ds$X.smoothed)
        ds.metrics <- compute.oracle.metrics(
          X.obs = X,
          X.hat = X.hat.full,
          X.true = X.true,
          method = "data_smoother_default",
          dgp.id = dgp.id,
          sigma = sigma,
          replicate = rep.idx,
          notes = sprintf("k.best=%d,trimmed=%s,n_model=%d",
                          ds$k.best, isTRUE(ds$trimmed), length(keep.idx)),
          modeled_fraction = length(keep.idx) / nrow(X)
        )
        append.csv(ds.metrics, oracle.path)
      }
    }

    ckpt.path <- file.path(
      run.dir, "checkpoints",
      sprintf("checkpoint_%s_sigma_%s_rep_%03d.rds",
              dgp.id,
              formatC(sigma, digits = 3, format = "f"),
              rep.idx)
    )
    saveRDS(
      list(
        dgp_id = dgp.id,
        sigma = sigma,
        replicate = rep.idx,
        seed = rep.seed,
        config = config,
        summary = cand.df,
        selection = sel.rows
      ),
      ckpt.path
    )

    elapsed.rep <- as.numeric(difftime(Sys.time(), run.start.time, units = "secs"))
    log.msg(sprintf("Completed dgp %s sigma %.4f replicate %d | elapsed %s",
                    dgp.id, sigma, rep.idx, format.seconds(elapsed.rep)))
  }
}

selection.df <- if (file.exists(selection.path)) utils::read.csv(selection.path, stringsAsFactors = FALSE) else data.frame()

if (nrow(selection.df) > 0L) {
  out.rows <- list()
  out.ptr <- 0L
  for (dgp.id in unique(selection.df$dgp_id)) {
    for (s in sort(unique(selection.df$sigma))) {
      for (sel in c("min", "one_se")) {
        sub <- selection.df[
          selection.df$dgp_id == dgp.id &
            selection.df$sigma == s &
            selection.df$selected == sel,
          , drop = FALSE
        ]
        if (nrow(sub) == 0L) next
        out.ptr <- out.ptr + 1L
        out.rows[[out.ptr]] <- data.frame(
          dgp_id = dgp.id,
          sigma = s,
          selected = sel,
          k_mode = mode.param(sub$k),
          n_eigenpairs_mode = mode.param(sub$n_eigenpairs),
          eta_mode = mode.param(sub$eta),
          scaled_mse_mean_avg = mean(sub$scaled_mse_mean, na.rm = TRUE),
          scaled_mse_mean_sd = stats::sd(sub$scaled_mse_mean, na.rm = TRUE),
          raw_mse_mean_avg = mean(sub$raw_mse_mean, na.rm = TRUE),
          modeled_fraction_mean_avg = mean(sub$modeled_fraction_mean, na.rm = TRUE),
          n_replicates = nrow(sub),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  selected.by.group <- do.call(rbind, out.rows)
  write.csv(selected.by.group,
            file = file.path(run.dir, "tables", "selected_params_by_dgp_sigma.csv"),
            row.names = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(run.dir, "artifacts", "session_info.txt"))

elapsed.total <- as.numeric(difftime(Sys.time(), run.start.time, units = "secs"))
log.msg("Run complete.")
log.msg("Total elapsed:", format.seconds(elapsed.total))
log.msg("Outputs written under:", run.dir)
