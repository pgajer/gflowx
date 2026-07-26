#' Local Correlation over Posterior rdgraph Refits
#'
#' Archive adapter that refits one or more vertex responses with an existing
#' `knn.riem.fit` object, requests posterior draws, and delegates the
#' model-independent local-correlation summary to [gflow::lcor.with.posterior()].
#'
#' @param fitted.model A `knn.riem.fit` object.
#' @param Z.abundances Numeric response vector or matrix to refit.
#' @param y.hat Fixed response field. Defaults to
#'   `fitted.model$fitted.values`.
#' @param per.column.gcv,n.candidates,n.cores Arguments passed to
#'   [refit.rdgraph.regression()].
#' @param n.posterior.samples,credible.level,seed Posterior sampling controls.
#' @param lcor.type Local-correlation weighting passed to
#'   [gflow::lcor.with.posterior()].
#' @param return.samples Logical; retain the local-correlation draws.
#' @param verbose Logical; report progress.
#'
#' @return A `lcor.posterior` object returned by
#'   [gflow::lcor.with.posterior()].
#'
#' @export
lcor.with.rdgraph.posterior <- function(
    fitted.model,
    Z.abundances,
    y.hat = fitted.model$fitted.values,
    per.column.gcv = TRUE,
    n.candidates = 40L,
    n.cores = 1L,
    n.posterior.samples = 500L,
    credible.level = 0.95,
    seed = 12345L,
    lcor.type = c("derivative", "unit", "sign"),
    return.samples = FALSE,
    verbose = TRUE
) {
    if (!inherits(fitted.model, "knn.riem.fit")) {
        stop("fitted.model must be a 'knn.riem.fit' object.", call. = FALSE)
    }
    graph <- fitted.model$graph
    if (is.null(graph$adj.list) || is.null(graph$edge.length.list)) {
        stop(
            "fitted.model$graph must contain adj.list and edge.length.list.",
            call. = FALSE
        )
    }
    if (is.null(y.hat) || length(y.hat) != length(graph$adj.list)) {
        stop("y.hat must contain one value per graph vertex.", call. = FALSE)
    }

    lcor.posterior <- .gflowx.get.namespace.export(
        "gflow", "lcor.with.posterior",
        api = "lcor.with.rdgraph.posterior()",
        install_hint = "Install gflow to summarize local associations."
    )
    lcor.type <- match.arg(lcor.type)

    refit <- refit.rdgraph.regression(
        fitted.model = fitted.model,
        y.new = Z.abundances,
        per.column.gcv = per.column.gcv,
        n.candidates = n.candidates,
        n.cores = n.cores,
        verbose = verbose,
        with.posterior = TRUE,
        return.posterior.samples = TRUE,
        credible.level = credible.level,
        n.posterior.samples = n.posterior.samples,
        posterior.seed = seed
    )

    lcor.posterior(
        adj.list = graph$adj.list,
        weight.list = graph$edge.length.list,
        y.hat = y.hat,
        Z.hat.samples = refit$posterior$samples,
        lcor.type = lcor.type,
        credible.level = credible.level,
        return.samples = return.samples,
        verbose = verbose
    )
}
