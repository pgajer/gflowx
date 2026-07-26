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

evaluate.split <- function(model, X, split, block.size = NULL) {
  refit <- gflowx::refit.rdgraph.regression(
    fitted.model = model,
    y.new = X,
    y.vertices = split$train.idx,
    per.column.gcv = FALSE,
    block.size = block.size,
    verbose = FALSE
  )
  X.hat <- as.matrix(refit$fitted.values)
  resid <- X[split$test.idx, , drop = FALSE] - X.hat[split$test.idx, , drop = FALSE]
  resid.scaled <- sweep(resid, 2L, split$col.scale, "/")
  list(
    scaled.mse = mean(resid.scaled^2),
    raw.mse = mean(resid^2)
  )
}

compute.oracle.metrics <- function(X.obs, X.hat, X.true, method, sigma, replicate, notes = NA_character_) {
  X.obs <- as.matrix(X.obs)
  X.hat <- as.matrix(X.hat)
  X.true <- as.matrix(X.true)
  rad.obs <- sqrt(rowSums(X.obs^2))
  rad.hat <- sqrt(rowSums(X.hat^2))
  rad.true <- sqrt(rowSums(X.true^2))
  data.frame(
    sigma = as.double(sigma),
    replicate = as.integer(replicate),
    method = method,
    n_rows = nrow(X.hat),
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

se.fun <- function(x) {
  x <- x[is.finite(x)]
  n <- length(x)
  if (n <= 1L) return(0.0)
  stats::sd(x) / sqrt(n)
}

mode.param <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  tab <- sort(table(x), decreasing = TRUE)
  as.numeric(names(tab)[1L])
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
  experiment_name = "noisy_circle_cv_phase1",
  mode = mode,
  base_seed = base.seed,
  n_points = 250L,
  radius = 1.0,
  noise_type = "normal",
  sampling_type = "random",
  sigma_grid = c(0.02, 0.05, 0.10, 0.20),
  n_replicates = 30L,
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
  log_every_tasks = 100L,
  run_data_smoother_baseline = TRUE
)

if (mode == "pilot") {
  config$n_replicates <- 2L
  config$n_repeats <- 2L
  config$n_folds <- 3L
}

if (!is.null(args$n_replicates)) config$n_replicates <- as.int.arg(args$n_replicates, config$n_replicates)
if (!is.null(args$n_repeats)) config$n_repeats <- as.int.arg(args$n_repeats, config$n_repeats)
if (!is.null(args$n_folds)) config$n_folds <- as.int.arg(args$n_folds, config$n_folds)
if (!is.null(args$log_every_tasks)) config$log_every_tasks <- as.int.arg(args$log_every_tasks, config$log_every_tasks)

run.id <- args$run_id %||% format(Sys.time(), "%Y%m%d_%H%M%S")
run.dir <- file.path(repo.root, "tests", "manual", "results", "noisy_circle_cv_phase1", run.id)
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
total.tasks <- length(config$sigma_grid) * config$n_replicates * tasks.per.replicate

log.msg("Starting experiment run.")
log.msg("Run directory:", run.dir)
log.msg("Total tasks:", total.tasks,
        "| sigma:", length(config$sigma_grid),
        "| replicates:", config$n_replicates,
        "| candidates:", nrow(candidate.grid),
        "| splits/replicate:", split.count)

cv.fold.path <- file.path(run.dir, "tables", "cv_fold_scores.csv")
cv.summary.path <- file.path(run.dir, "tables", "cv_candidate_summary.csv")
selection.path <- file.path(run.dir, "tables", "selection_by_replicate.csv")
oracle.path <- file.path(run.dir, "tables", "oracle_metrics_by_replicate.csv")

global.task.idx <- 0L
run.start.time <- Sys.time()

for (sigma.idx in seq_along(config$sigma_grid)) {
  sigma <- as.double(config$sigma_grid[sigma.idx])
  for (rep.idx in seq_len(config$n_replicates)) {
    rep.seed <- as.integer(base.seed + sigma.idx * 100000L + rep.idx)
    set.seed(rep.seed)

    log.msg(sprintf("Sigma %.4f | replicate %d/%d | seed %d",
                    sigma, rep.idx, config$n_replicates, rep.seed))

    X.df <- generate.circle.data(
      n = config$n_points,
      radius = config$radius,
      noise = sigma,
      type = config$sampling_type,
      noise.type = config$noise_type,
      seed = rep.seed
    )
    X <- as.matrix(X.df[, c("x", "y"), drop = FALSE])
    X.true <- as.matrix(X.df[, c("x.noise.free", "y.noise.free"), drop = FALSE])

    connected.k <- NULL
    n.ccomp.by.k <- NULL
    if (isTRUE(config$use_connectivity_prescreen)) {
      k.min <- min(config$k_grid)
      k.max <- max(config$k_grid)
      pre.res <- tryCatch(
        {
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
          if (is.null(st) || !("k" %in% names(st)) || !("n_ccomp" %in% names(st))) {
            stop("Connectivity pre-screen failed: summary(create.iknn.graphs()) missing k or n_ccomp.")
          }
          list(
            connected.k = as.integer(st$k[st$n_ccomp == 1L]),
            n.ccomp.by.k = stats::setNames(as.integer(st$n_ccomp), as.character(as.integer(st$k)))
          )
        },
        error = function(e) e
      )

      if (inherits(pre.res, "error")) {
        log.msg(sprintf("  Connectivity pre-screen failed; falling back to direct fit checks | error=%s",
                        conditionMessage(pre.res)))
      } else {
        connected.k <- pre.res$connected.k
        n.ccomp.by.k <- pre.res$n.ccomp.by.k
        disconnected.k <- setdiff(as.integer(config$k_grid), connected.k)
        if (length(disconnected.k) > 0L) {
          log.msg(sprintf("  Connectivity pre-screen: disconnected k values skipped: %s",
                          paste(disconnected.k, collapse = ",")))
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
    fit.cache <- vector("list", nrow(pair.grid))

    for (pair.idx in seq_len(nrow(pair.grid))) {
      k.val <- as.integer(pair.grid$k[pair.idx])
      eig.val <- as.integer(pair.grid$n_eigenpairs[pair.idx])

      log.msg(sprintf("  Fitting base model | k=%d | n.eigenpairs=%d", k.val, eig.val))
      fit.fail.status <- "fit_error"
      fit.fail.message <- NA_character_

      if (!is.null(connected.k) && !(k.val %in% connected.k)) {
        n.ccomp <- n.ccomp.by.k[[as.character(k.val)]] %||% NA_integer_
        fit.base <- simpleError(sprintf(
          "Pre-screen: graph at k=%d disconnected (n_ccomp=%d). Base fit skipped.",
          k.val, as.integer(n.ccomp)
        ))
        fit.ok <- FALSE
        fit.fail.status <- "graph_disconnected"
        fit.fail.message <- conditionMessage(fit.base)
        log.msg(sprintf("    Base fit skipped | k=%d | n.eigenpairs=%d | error=%s",
                        k.val, eig.val, conditionMessage(fit.base)))
      } else {
        fit.base <- tryCatch(
          gflowx::fit.rdgraph.regression(
            X = X,
            y = proxy.y,
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
          log.msg(sprintf("    Base fit failed | k=%d | n.eigenpairs=%d | error=%s",
                          k.val, eig.val, fit.fail.message))
        } else {
          fit.cache[[pair.idx]] <- fit.base
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

          if (!fit.ok) {
            status <- fit.fail.status
            err.msg <- fit.fail.message
          } else {
            eval.res <- tryCatch(
              evaluate.split(
                model = fit.eta,
                X = X,
                split = split,
                block.size = config$block_size
              ),
              error = function(e) e
            )

            if (inherits(eval.res, "error")) {
              status <- "refit_error"
              err.msg <- conditionMessage(eval.res)
            } else {
              scaled.mse <- eval.res$scaled.mse
              raw.mse <- eval.res$raw.mse
            }
          }

          rep.rows[[rep.row.ptr]] <- data.frame(
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
    }

    rep.df <- do.call(rbind, rep.rows[seq_len(rep.row.ptr)])
    append.csv(rep.df, cv.fold.path)

    valid.df <- rep.df[is.finite(rep.df$scaled_mse), , drop = FALSE]

    if (nrow(valid.df) == 0L) {
      log.msg(sprintf("  No valid CV rows for sigma=%.4f replicate=%d; skipping selections.", sigma, rep.idx))
      sel.row <- data.frame(
        sigma = sigma,
        replicate = rep.idx,
        selected = "none",
        k = NA_integer_,
        n_eigenpairs = NA_integer_,
        eta = NA_real_,
        scaled_mse_mean = NA_real_,
        scaled_mse_se = NA_real_,
        raw_mse_mean = NA_real_,
        stringsAsFactors = FALSE
      )
      append.csv(sel.row, selection.path)
      next
    }

    mean.df <- aggregate(
      cbind(scaled_mse, raw_mse) ~ k + n_eigenpairs + eta,
      data = valid.df,
      FUN = mean
    )
    se.df <- aggregate(
      scaled_mse ~ k + n_eigenpairs + eta,
      data = valid.df,
      FUN = se.fun
    )
    names(mean.df)[names(mean.df) == "scaled_mse"] <- "scaled_mse_mean"
    names(mean.df)[names(mean.df) == "raw_mse"] <- "raw_mse_mean"
    names(se.df)[names(se.df) == "scaled_mse"] <- "scaled_mse_se"
    cand.df <- merge(mean.df, se.df, by = c("k", "n_eigenpairs", "eta"), all = TRUE)
    cand.df$sigma <- sigma
    cand.df$replicate <- rep.idx
    cand.df <- cand.df[order(cand.df$scaled_mse_mean, cand.df$raw_mse_mean), ]

    append.csv(cand.df[, c("sigma", "replicate", "k", "n_eigenpairs", "eta",
                           "scaled_mse_mean", "scaled_mse_se", "raw_mse_mean")], cv.summary.path)

    best.min <- cand.df[1L, ]
    threshold.1se <- best.min$scaled_mse_mean + best.min$scaled_mse_se
    cand.1se <- cand.df[cand.df$scaled_mse_mean <= threshold.1se, , drop = FALSE]
    cand.1se <- cand.1se[order(cand.1se$n_eigenpairs, -cand.1se$k, -cand.1se$eta), ]
    best.1se <- cand.1se[1L, ]

    sel.rows <- rbind(
      data.frame(
        sigma = sigma, replicate = rep.idx, selected = "min",
        k = best.min$k, n_eigenpairs = best.min$n_eigenpairs, eta = best.min$eta,
        scaled_mse_mean = best.min$scaled_mse_mean, scaled_mse_se = best.min$scaled_mse_se,
        raw_mse_mean = best.min$raw_mse_mean, stringsAsFactors = FALSE
      ),
      data.frame(
        sigma = sigma, replicate = rep.idx, selected = "one_se",
        k = best.1se$k, n_eigenpairs = best.1se$n_eigenpairs, eta = best.1se$eta,
        scaled_mse_mean = best.1se$scaled_mse_mean, scaled_mse_se = best.1se$scaled_mse_se,
        raw_mse_mean = best.1se$raw_mse_mean, stringsAsFactors = FALSE
      )
    )
    append.csv(sel.rows, selection.path)

    for (sel.idx in seq_len(nrow(sel.rows))) {
      k.sel <- as.integer(sel.rows$k[sel.idx])
      eig.sel <- as.integer(sel.rows$n_eigenpairs[sel.idx])
      eta.sel <- as.double(sel.rows$eta[sel.idx])
      key.idx <- which(pair.grid$k == k.sel & pair.grid$n_eigenpairs == eig.sel)[1L]
      fit.base <- fit.cache[[key.idx]]
      fit.sel <- set.model.eta(fit.base, eta.sel)
      refit.full <- gflowx::refit.rdgraph.regression(
        fitted.model = fit.sel,
        y.new = X,
        per.column.gcv = FALSE,
        verbose = FALSE
      )
      X.hat <- as.matrix(refit.full$fitted.values)
      om <- compute.oracle.metrics(
        X.obs = X,
        X.hat = X.hat,
        X.true = X.true,
        method = paste0("cv_", sel.rows$selected[sel.idx]),
        sigma = sigma,
        replicate = rep.idx,
        notes = sprintf("k=%d,n.eigenpairs=%d,eta=%.8f", k.sel, eig.sel, eta.sel)
      )
      append.csv(om, oracle.path)
    }

    baseline.obs <- compute.oracle.metrics(
      X.obs = X,
      X.hat = X,
      X.true = X.true,
      method = "observed_no_smoothing",
      sigma = sigma,
      replicate = rep.idx,
      notes = "identity baseline"
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
        log.msg(sprintf("  data.smoother baseline failed | sigma=%.4f rep=%d | %s",
                        sigma, rep.idx, conditionMessage(ds)))
      } else {
        keep.idx <- if (isTRUE(ds$trimmed)) as.integer(ds$kept.rows) else seq_len(nrow(X))
        X.obs.sub <- X[keep.idx, , drop = FALSE]
        X.true.sub <- X.true[keep.idx, , drop = FALSE]
        X.hat.sub <- as.matrix(ds$X.smoothed)
        ds.metrics <- compute.oracle.metrics(
          X.obs = X.obs.sub,
          X.hat = X.hat.sub,
          X.true = X.true.sub,
          method = "data_smoother_default",
          sigma = sigma,
          replicate = rep.idx,
          notes = sprintf("k.best=%d,trimmed=%s,n_rows=%d",
                          ds$k.best, isTRUE(ds$trimmed), nrow(X.hat.sub))
        )
        append.csv(ds.metrics, oracle.path)
      }
    }

    ckpt.path <- file.path(
      run.dir, "checkpoints",
      sprintf("checkpoint_sigma_%s_rep_%03d.rds", formatC(sigma, digits = 3, format = "f"), rep.idx)
    )
    saveRDS(
      list(
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
    log.msg(sprintf("Completed sigma %.4f replicate %d | elapsed %s",
                    sigma, rep.idx, format.seconds(elapsed.rep)))
  }
}

selection.df <- if (file.exists(selection.path)) utils::read.csv(selection.path, stringsAsFactors = FALSE) else data.frame()

if (nrow(selection.df) > 0L) {
  out.rows <- list()
  out.ptr <- 0L
  for (s in unique(selection.df$sigma)) {
    for (sel in c("min", "one_se")) {
      sub <- selection.df[selection.df$sigma == s & selection.df$selected == sel, , drop = FALSE]
      if (nrow(sub) == 0L) next
      out.ptr <- out.ptr + 1L
      out.rows[[out.ptr]] <- data.frame(
        sigma = s,
        selected = sel,
        k_mode = mode.param(sub$k),
        n_eigenpairs_mode = mode.param(sub$n_eigenpairs),
        eta_mode = mode.param(sub$eta),
        scaled_mse_mean_avg = mean(sub$scaled_mse_mean, na.rm = TRUE),
        scaled_mse_mean_sd = stats::sd(sub$scaled_mse_mean, na.rm = TRUE),
        raw_mse_mean_avg = mean(sub$raw_mse_mean, na.rm = TRUE),
        n_replicates = nrow(sub),
        stringsAsFactors = FALSE
      )
    }
  }
  selected.by.sigma <- do.call(rbind, out.rows)
  write.csv(selected.by.sigma,
            file = file.path(run.dir, "tables", "selected_params_by_sigma.csv"),
            row.names = FALSE)
}

writeLines(capture.output(sessionInfo()), file.path(run.dir, "artifacts", "session_info.txt"))

elapsed.total <- as.numeric(difftime(Sys.time(), run.start.time, units = "secs"))
log.msg("Run complete.")
log.msg("Total elapsed:", format.seconds(elapsed.total))
log.msg("Outputs written under:", run.dir)
