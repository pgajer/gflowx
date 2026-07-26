rdgraph.rmse <- function(x, truth) {
  sqrt(mean((as.double(x) - as.double(truth))^2, na.rm = TRUE))
}

rdgraph.safe.cor <- function(x, truth) {
  x <- as.double(x)
  truth <- as.double(truth)
  if (length(x) != length(truth) || length(x) < 2L) return(NA_real_)
  if (!all(is.finite(x)) || !all(is.finite(truth))) return(NA_real_)
  if (stats::sd(x) <= sqrt(.Machine$double.eps)) return(NA_real_)
  if (stats::sd(truth) <= sqrt(.Machine$double.eps)) return(NA_real_)
  stats::cor(x, truth)
}

rdgraph.safe.var.ratio <- function(x, truth) {
  truth.var <- stats::var(as.double(truth))
  if (!is.finite(truth.var) || truth.var <= sqrt(.Machine$double.eps)) {
    return(NA_real_)
  }
  stats::var(as.double(x)) / truth.var
}

rdgraph.boundary.rmse <- function(y.hat, y.true, data, width = NULL) {
  if (identical(data$type, "circle") && !is.null(data$theta)) {
    width <- width %||% 0.25
    theta <- as.double(data$theta) %% (2 * pi)
    boundary <- theta <= width | theta >= (2 * pi - width)
  } else if (identical(data$type, "line") && !is.null(data$coord)) {
    width <- width %||% 0.10
    coord <- as.double(data$coord)
    rng <- range(coord, finite = TRUE)
    if (!all(is.finite(rng)) || diff(rng) <= sqrt(.Machine$double.eps)) {
      return(NA_real_)
    }
    margin <- width * diff(rng)
    boundary <- coord <= rng[[1L]] + margin | coord >= rng[[2L]] - margin
  } else {
    return(NA_real_)
  }

  if (!any(boundary)) return(NA_real_)
  rdgraph.rmse(y.hat[boundary], y.true[boundary])
}

rdgraph.variant.metrics <- function(case.result, variant.name, variant.result) {
  truth <- as.double(case.result$data$y.true)
  response <- as.double(variant.result$y)

  if (!is.null(variant.result$error)) {
    return(data.frame(
      case_id = case.result$case$id,
      case_title = case.result$case$title,
      round_id = case.result$case$round$id %||% NA_character_,
      round_title = case.result$case$round$title %||% NA_character_,
      group = case.result$case$group,
      response_variant = variant.name,
      graph_model = variant.result$graph.model %||% NA_character_,
      fit_status = "error",
      error_message = conditionMessage(variant.result$error),
      n = length(truth),
      rmse_raw = rdgraph.rmse(response, truth),
      rmse_fit = NA_real_,
      rmse_improvement = NA_real_,
      rmse_improvement_fraction = NA_real_,
      cor_raw_truth = rdgraph.safe.cor(response, truth),
      cor_fit_truth = NA_real_,
      var_ratio_fit_truth = NA_real_,
      residual_mean = NA_real_,
      residual_sd = NA_real_,
      residual_max_abs = NA_real_,
      boundary_rmse_fit = NA_real_,
      fit_elapsed_sec = as.double(variant.result$elapsed),
      k = case.result$case$fit.args$k %||% NA_integer_,
      n_eigenpairs = case.result$case$fit.args$n.eigenpairs %||% NA_integer_,
      n_edges = NA_integer_,
      stringsAsFactors = FALSE
    ))
  }

  y.hat <- as.double(variant.result$fit$fitted.values)
  residual <- y.hat - truth
  rmse.raw <- rdgraph.rmse(response, truth)
  rmse.fit <- rdgraph.rmse(y.hat, truth)
  n.edges <- variant.result$fit$graph$n.edges %||% NA_integer_

  data.frame(
    case_id = case.result$case$id,
    case_title = case.result$case$title,
    round_id = case.result$case$round$id %||% NA_character_,
    round_title = case.result$case$round$title %||% NA_character_,
    group = case.result$case$group,
    response_variant = variant.name,
    graph_model = variant.result$graph.model %||% NA_character_,
    fit_status = "ok",
    error_message = NA_character_,
    n = length(truth),
    rmse_raw = rmse.raw,
    rmse_fit = rmse.fit,
    rmse_improvement = rmse.raw - rmse.fit,
    rmse_improvement_fraction = (rmse.raw - rmse.fit) / max(rmse.raw, .Machine$double.eps),
    cor_raw_truth = rdgraph.safe.cor(response, truth),
    cor_fit_truth = rdgraph.safe.cor(y.hat, truth),
    var_ratio_fit_truth = rdgraph.safe.var.ratio(y.hat, truth),
    residual_mean = mean(residual),
    residual_sd = stats::sd(residual),
    residual_max_abs = max(abs(residual)),
    boundary_rmse_fit = rdgraph.boundary.rmse(y.hat, truth, case.result$data),
    fit_elapsed_sec = as.double(variant.result$elapsed),
    k = case.result$case$fit.args$k %||% NA_integer_,
    n_eigenpairs = case.result$case$fit.args$n.eigenpairs %||% NA_integer_,
    n_edges = as.integer(n.edges),
    stringsAsFactors = FALSE
  )
}

rdgraph.case.metrics <- function(case.result) {
  rows <- Map(
    f = function(name, variant) rdgraph.variant.metrics(case.result, name, variant),
    name = names(case.result$variants),
    variant = case.result$variants
  )
  do.call(rbind, rows)
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}
