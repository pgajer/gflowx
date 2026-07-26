rdgraph.run.k.sweep <- function(case, data) {
  k.grid <- case$k.grid %||% integer(0)
  if (!length(k.grid)) return(NULL)

  response.name <- names(data$responses)[1L]
  y <- as.double(data$responses[[response.name]])
  rows <- vector("list", length(k.grid))
  fits <- vector("list", length(k.grid))
  graphs <- vector("list", length(k.grid))
  lcc.fits <- vector("list", length(k.grid))
  component.fits <- vector("list", length(k.grid))
  names(fits) <- paste0("k", k.grid)
  names(graphs) <- paste0("k", k.grid)
  names(lcc.fits) <- paste0("k", k.grid)
  names(component.fits) <- paste0("k", k.grid)

  for (i in seq_along(k.grid)) {
    k <- as.integer(k.grid[[i]])
    fit.args <- modifyList(case$fit.args, list(k = k))
    elapsed <- system.time({
      fit <- tryCatch(
        do.call(gflowx::fit.rdgraph.regression, c(list(X = data$X, y = y), fit.args)),
        error = function(e) e
      )
    })[["elapsed"]]

    if (inherits(fit, "error")) {
      constructed.graph <- rdgraph.k.sweep.construct.graph(data$X, fit.args)
      graph.row <- if (is.null(constructed.graph$error)) {
        graphs[[i]] <- constructed.graph$graph
        rdgraph.graph.summary(
          constructed.graph$graph$adj.list,
          constructed.graph$graph$edge.length.list,
          k = k
        )
      } else {
        data.frame(
          k = k,
          edge_count = NA_integer_,
          component_count = NA_integer_,
          largest_component_size = NA_integer_,
          largest_component_fraction = NA_real_,
          second_largest_component_size = NA_integer_,
          mean_degree = NA_real_,
          min_degree = NA_integer_,
          max_degree = NA_integer_,
          stringsAsFactors = FALSE
        )
      }
      lcc.row <- if (is.null(constructed.graph$error)) {
        rdgraph.k.sweep.lcc.row(
          graph = constructed.graph$graph,
          fit.args = fit.args,
          data = data,
          y = y,
          k = k
        )
      } else {
        list(
          summary = rdgraph.k.sweep.lcc.empty(
            NA_character_,
            message = conditionMessage(fit)
          ),
          lcc.fit = NULL
        )
      }
      lcc.fits[i] <- list(lcc.row$lcc.fit)
      component.row <- if (is.null(constructed.graph$error)) {
        rdgraph.k.sweep.componentwise.row(
          graph = constructed.graph$graph,
          fit.args = fit.args,
          data = data,
          y = y,
          k = k
        )
      } else {
        list(summary = rdgraph.k.sweep.componentwise.empty(NA_character_), component.fit = NULL)
      }
      component.fits[i] <- list(component.row$component.fit)

      rows[[i]] <- cbind(data.frame(
        k = k,
        fit_status = "error",
        error_message = conditionMessage(fit),
        rmse_fit = NA_real_,
        cor_fit_truth = NA_real_,
        gcv = NA_real_,
        elapsed_sec = elapsed,
        stringsAsFactors = FALSE
      ),
      graph.row[, setdiff(names(graph.row), "k"), drop = FALSE],
      lcc.row$summary,
      component.row$summary)
      next
    }

    y.hat <- as.double(fit$fitted.values)
    graph.row <- rdgraph.graph.summary(
      fit$graph$adj.list,
      fit$graph$edge.length.list,
      k = k
    )
    lcc.row <- rdgraph.k.sweep.lcc.row(
      graph = list(
        adj.list = fit$graph$adj.list,
        edge.length.list = fit$graph$edge.length.list
      ),
      fit.args = fit.args,
      data = data,
      y = y,
      k = k
    )
    lcc.fits[i] <- list(lcc.row$lcc.fit)
    component.row <- rdgraph.k.sweep.componentwise.row(
      graph = list(
        adj.list = fit$graph$adj.list,
        edge.length.list = fit$graph$edge.length.list
      ),
      fit.args = fit.args,
      data = data,
      y = y,
      k = k
    )
    component.fits[i] <- list(component.row$component.fit)
    graphs[[i]] <- list(
      adj.list = fit$graph$adj.list,
      edge.length.list = fit$graph$edge.length.list
    )

    rows[[i]] <- cbind(
      data.frame(
        k = k,
        fit_status = "ok",
        error_message = NA_character_,
        rmse_fit = rdgraph.rmse(y.hat, data$y.true),
        cor_fit_truth = rdgraph.safe.cor(y.hat, data$y.true),
        gcv = rdgraph.gcv.value(fit),
        elapsed_sec = elapsed,
        stringsAsFactors = FALSE
      ),
      graph.row[, setdiff(names(graph.row), "k"), drop = FALSE],
      lcc.row$summary,
      component.row$summary
    )
    fits[[i]] <- fit
  }

  summary <- do.call(rbind, rows)
  ok <- summary$fit_status == "ok"
  if (any(ok)) {
    summary$min_norm_gcv <- NA_real_
    best.gcv <- min(summary$gcv[ok], na.rm = TRUE)
    summary$min_norm_gcv[ok] <- summary$gcv[ok] / best.gcv
    med.gcv <- stats::median(summary$gcv[ok], na.rm = TRUE)
    mad.gcv <- stats::mad(summary$gcv[ok], constant = 1.4826, na.rm = TRUE)
    if (!is.finite(mad.gcv) || mad.gcv <= sqrt(.Machine$double.eps)) {
      mad.gcv <- stats::IQR(summary$gcv[ok], na.rm = TRUE)
    }
    if (!is.finite(mad.gcv) || mad.gcv <= sqrt(.Machine$double.eps)) {
      mad.gcv <- stats::sd(summary$gcv[ok], na.rm = TRUE)
    }
    if (!is.finite(mad.gcv) || mad.gcv <= sqrt(.Machine$double.eps)) {
      mad.gcv <- 1
    }
    summary$robust_norm_gcv <- NA_real_
    summary$robust_norm_gcv[ok] <- (summary$gcv[ok] - med.gcv) / mad.gcv
  }

  stability <- rdgraph.k.sweep.stability(graphs, k.grid)
  selection <- rdgraph.k.sweep.selection(summary)

  list(
    response_variant = response.name,
    summary = summary,
    stability = stability,
    selection = selection,
    fits = fits,
    graphs = graphs,
    lcc.fits = lcc.fits,
    component.fits = component.fits
  )
}

rdgraph.k.sweep.construct.graph <- function(X, fit.args) {
  graph <- tryCatch({
    graph.args <- list(
      X = X,
      kmin = fit.args$k,
      kmax = fit.args$k,
      max.path.edge.ratio.deviation.thld = fit.args$max.ratio.threshold %||% 0.1,
      path.edge.ratio.percentile = fit.args$path.edge.ratio.percentile %||% 0.5,
      threshold.percentile = fit.args$threshold.percentile %||% 0,
      compute.full = TRUE,
      pca.dim = fit.args$pca.dim,
      variance.explained = fit.args$variance.explained %||% 0.95,
      n.cores = 1L,
      verbose = FALSE,
      knn.cache.path = fit.args$knn.cache.path %||% NULL,
      knn.cache.mode = fit.args$knn.cache.mode %||% "none"
    )
    graphs <- do.call(create.iknn.graphs, graph.args)
    g <- graphs$geom_pruned_graphs[[1L]]
    list(adj.list = g$adj_list, edge.length.list = g$weight_list)
  }, error = function(e) e)

  if (inherits(graph, "error")) {
    list(graph = NULL, error = graph)
  } else {
    list(graph = graph, error = NULL)
  }
}

rdgraph.k.sweep.lcc.empty <- function(status,
                                      message = NA_character_,
                                      n = NA_integer_,
                                      fraction = NA_real_,
                                      elapsed = NA_real_) {
  data.frame(
    lcc_fit_status = status,
    lcc_error_message = message,
    lcc_n = as.integer(n),
    lcc_fraction = as.double(fraction),
    lcc_rmse_fit = NA_real_,
    lcc_cor_fit_truth = NA_real_,
    lcc_gcv = NA_real_,
    lcc_elapsed_sec = as.double(elapsed),
    lcc_edge_count = NA_integer_,
    stringsAsFactors = FALSE
  )
}

rdgraph.k.sweep.lcc.row <- function(graph, fit.args, data, y, k) {
  cc <- graph.connected.components(graph$adj.list)
  comp.sizes <- tabulate(cc)
  comp.sizes <- comp.sizes[comp.sizes > 0L]
  lcc.indices <- rdgraph.largest.component.indices(graph$adj.list)
  lcc.n <- length(lcc.indices)
  lcc.frac <- lcc.n / length(graph$adj.list)

  if (length(comp.sizes) <= 1L) {
    return(list(
      summary = rdgraph.k.sweep.lcc.empty("connected_graph", n = lcc.n, fraction = lcc.frac),
      lcc.fit = NULL
    ))
  }

  if (lcc.n < 10L) {
    return(list(
      summary = rdgraph.k.sweep.lcc.empty(
        "skipped_small_component",
        message = sprintf("Largest connected component has n=%d; gflowx::fit.rdgraph.regression() requires at least 10 observations.", lcc.n),
        n = lcc.n,
        fraction = lcc.frac
      ),
      lcc.fit = NULL
    ))
  }

  lcc.graph <- rdgraph.induced.graph(
    graph$adj.list,
    graph$edge.length.list,
    lcc.indices
  )
  lcc.args <- modifyList(fit.args, list(
    k = min(max(2L, as.integer(k)), lcc.n - 1L),
    adj.list = lcc.graph$adj.list,
    weight.list = lcc.graph$weight.list,
    use.counting.measure = TRUE,
    apply.geometric.pruning = FALSE,
    max.ratio.threshold = 0,
    threshold.percentile = 0
  ))
  if (!is.null(lcc.args$n.eigenpairs)) {
    max.eigenpairs <- if (lcc.n > 12L) lcc.n - 2L else lcc.n
    lcc.args$n.eigenpairs <- min(as.integer(lcc.args$n.eigenpairs), max(10L, max.eigenpairs))
  }

  elapsed <- system.time({
    lcc.fit <- tryCatch(
      do.call(
        gflowx::fit.rdgraph.regression,
        c(
          list(
            X = data$X[lcc.indices, , drop = FALSE],
            y = y[lcc.indices]
          ),
          lcc.args
        )
      ),
      error = function(e) e
    )
  })[["elapsed"]]

  if (inherits(lcc.fit, "error")) {
    return(list(
      summary = rdgraph.k.sweep.lcc.empty(
        "error",
        message = conditionMessage(lcc.fit),
        n = lcc.n,
        fraction = lcc.frac,
        elapsed = elapsed
      ),
      lcc.fit = NULL
    ))
  }

  y.hat <- as.double(lcc.fit$fitted.values)
  graph.row <- rdgraph.graph.summary(
    lcc.fit$graph$adj.list,
    lcc.fit$graph$edge.length.list,
    k = k
  )
  summary <- data.frame(
    lcc_fit_status = "ok",
    lcc_error_message = NA_character_,
    lcc_n = as.integer(lcc.n),
    lcc_fraction = as.double(lcc.frac),
    lcc_rmse_fit = rdgraph.rmse(y.hat, data$y.true[lcc.indices]),
    lcc_cor_fit_truth = rdgraph.safe.cor(y.hat, data$y.true[lcc.indices]),
    lcc_gcv = rdgraph.gcv.value(lcc.fit),
    lcc_elapsed_sec = as.double(elapsed),
    lcc_edge_count = graph.row$edge_count,
    stringsAsFactors = FALSE
  )

  list(
    summary = summary,
    lcc.fit = list(
      fit = lcc.fit,
      indices = lcc.indices,
      k = as.integer(k),
      fraction = as.double(lcc.frac)
    )
  )
}

rdgraph.k.sweep.componentwise.empty <- function(status,
                                                message = NA_character_,
                                                n = NA_integer_,
                                                fraction = NA_real_,
                                                n.components = NA_integer_,
                                                elapsed = NA_real_) {
  data.frame(
    component_fit_status = status,
    component_fit_error_message = message,
    component_fit_components = as.integer(n.components),
    component_fit_n = as.integer(n),
    component_fit_fraction = as.double(fraction),
    component_fit_rmse = NA_real_,
    component_fit_cor_truth = NA_real_,
    component_fit_elapsed_sec = as.double(elapsed),
    stringsAsFactors = FALSE
  )
}

rdgraph.k.sweep.fit.component <- function(indices, graph, fit.args, data, y, k) {
  component.graph <- rdgraph.induced.graph(
    graph$adj.list,
    graph$edge.length.list,
    indices
  )
  n.component <- length(indices)
  component.args <- modifyList(fit.args, list(
    k = min(max(2L, as.integer(k)), n.component - 1L),
    adj.list = component.graph$adj.list,
    weight.list = component.graph$weight.list,
    use.counting.measure = TRUE,
    apply.geometric.pruning = FALSE,
    max.ratio.threshold = 0,
    threshold.percentile = 0
  ))
  if (!is.null(component.args$n.eigenpairs)) {
    max.eigenpairs <- if (n.component > 12L) n.component - 2L else n.component
    component.args$n.eigenpairs <- min(as.integer(component.args$n.eigenpairs), max(10L, max.eigenpairs))
  }

  fit <- tryCatch(
    do.call(
      gflowx::fit.rdgraph.regression,
      c(
        list(
          X = data$X[indices, , drop = FALSE],
          y = y[indices]
        ),
        component.args
      )
    ),
    error = function(e) e
  )
  fit
}

rdgraph.k.sweep.componentwise.row <- function(graph, fit.args, data, y, k, min.component.n = 10L) {
  cc <- graph.connected.components(graph$adj.list)
  comp.ids <- sort(unique(cc))
  comp.indices <- lapply(comp.ids, function(id) which(cc == id))
  comp.indices <- comp.indices[lengths(comp.indices) >= min.component.n]

  if (length(comp.ids) <= 1L) {
    return(list(
      summary = rdgraph.k.sweep.componentwise.empty(
        "connected_graph",
        n = length(graph$adj.list),
        fraction = 1,
        n.components = 1L
      ),
      component.fit = NULL
    ))
  }
  if (!length(comp.indices)) {
    return(list(
      summary = rdgraph.k.sweep.componentwise.empty(
        "skipped_small_components",
        message = sprintf("No connected component has at least %d observations.", min.component.n),
        n = 0L,
        fraction = 0,
        n.components = 0L
      ),
      component.fit = NULL
    ))
  }

  elapsed <- system.time({
    component.results <- lapply(comp.indices, function(indices) {
      fit <- rdgraph.k.sweep.fit.component(indices, graph, fit.args, data, y, k)
      list(indices = indices, fit = fit)
    })
  })[["elapsed"]]

  ok <- !vapply(component.results, function(x) inherits(x$fit, "error"), logical(1L))
  if (!any(ok)) {
    messages <- vapply(component.results, function(x) {
      if (inherits(x$fit, "error")) conditionMessage(x$fit) else NA_character_
    }, character(1L))
    return(list(
      summary = rdgraph.k.sweep.componentwise.empty(
        "error",
        message = paste(stats::na.omit(messages), collapse = " | "),
        n = 0L,
        fraction = 0,
        n.components = 0L,
        elapsed = elapsed
      ),
      component.fit = NULL
    ))
  }

  component.results <- component.results[ok]
  y.hat.all <- rep(NA_real_, nrow(data$X))
  component.table <- do.call(rbind, lapply(seq_along(component.results), function(i) {
    res <- component.results[[i]]
    indices <- res$indices
    y.hat <- as.double(res$fit$fitted.values)
    y.hat.all[indices] <<- y.hat
    data.frame(
      component = i,
      n = length(indices),
      fraction = length(indices) / nrow(data$X),
      rmse_fit = rdgraph.rmse(y.hat, data$y.true[indices]),
      cor_fit_truth = rdgraph.safe.cor(y.hat, data$y.true[indices]),
      stringsAsFactors = FALSE
    )
  }))
  fitted <- is.finite(y.hat.all)
  summary <- data.frame(
    component_fit_status = "ok",
    component_fit_error_message = NA_character_,
    component_fit_components = as.integer(length(component.results)),
    component_fit_n = as.integer(sum(fitted)),
    component_fit_fraction = mean(fitted),
    component_fit_rmse = rdgraph.rmse(y.hat.all[fitted], data$y.true[fitted]),
    component_fit_cor_truth = rdgraph.safe.cor(y.hat.all[fitted], data$y.true[fitted]),
    component_fit_elapsed_sec = as.double(elapsed),
    stringsAsFactors = FALSE
  )

  list(
    summary = summary,
    component.fit = list(
      k = as.integer(k),
      y.hat = y.hat.all,
      components = component.results,
      component.table = component.table,
      fraction = mean(fitted)
    )
  )
}

rdgraph.k.sweep.selection <- function(summary) {
  ok <- summary$fit_status == "ok"
  if (!any(ok)) return(data.frame())
  s <- summary[ok, , drop = FALSE]
  criteria <- list(
    min.rmse_fit = c("rmse_fit", "lower"),
    max.cor_fit_truth = c("cor_fit_truth", "higher"),
    min.gcv = c("gcv", "lower"),
    min.min_norm_gcv = c("min_norm_gcv", "lower"),
    min.robust_norm_gcv = c("robust_norm_gcv", "lower")
  )
  rows <- lapply(names(criteria), function(nm) {
    field <- criteria[[nm]][[1L]]
    direction <- criteria[[nm]][[2L]]
    vals <- s[[field]]
    finite <- is.finite(vals)
    if (!any(finite)) {
      return(data.frame(
        criterion = nm,
        k_selected = NA_integer_,
        value = NA_real_,
        stringsAsFactors = FALSE
      ))
    }
    idx.local <- if (identical(direction, "higher")) {
      which.max(vals[finite])
    } else {
      which.min(vals[finite])
    }
    idx <- which(finite)[[idx.local]]
    data.frame(
      criterion = nm,
      k_selected = s$k[[idx]],
      value = vals[[idx]],
      stringsAsFactors = FALSE
    )
  })
  do.call(rbind, rows)
}

rdgraph.k.sweep.stability <- function(fits, k.grid) {
  rows <- list()
  ptr <- 0L
  for (i in seq_len(length(k.grid) - 1L)) {
    graph1 <- fits[[i]]
    graph2 <- fits[[i + 1L]]
    if (is.null(graph1) || is.null(graph2)) next
    ptr <- ptr + 1L
    rows[[ptr]] <- data.frame(
      transition_from_k = k.grid[[i]],
      transition_to_k = k.grid[[i + 1L]],
      edit_distance = graph.edit.distance(
        graph1$adj.list,
        graph1$edge.length.list,
        graph2$adj.list,
        graph2$edge.length.list,
        edge.cost = 1,
        weight.cost.factor = 0
      ),
      js_divergence_degree = rdgraph.degree.js(graph1$adj.list, graph2$adj.list),
      stringsAsFactors = FALSE
    )
  }
  if (!ptr) return(data.frame())
  do.call(rbind, rows)
}
