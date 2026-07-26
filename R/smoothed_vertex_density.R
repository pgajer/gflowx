#' Compute Smoothed Vertex Density with an Archived Rdgraph Fit
#'
#' Computes the inverse nearest-edge-length density and refits an archived
#' rdgraph model to estimate its graph-conditional expectation.
#'
#' @param fit An object returned by [fit.rdgraph.regression()].
#' @param adj.list Optional 1-based adjacency list. If omitted, it is read from
#'   `fit$adj.list` or `fit$graph$adj.list`.
#' @param weight.list Optional edge-length list parallel to `adj.list`.
#'
#' @return Numeric vector normalized to `[0, 1]`.
#' @export
compute.smoothed.density <- function(fit, adj.list = NULL,
                                     weight.list = NULL) {
    if (is.null(adj.list)) {
        adj.list <- fit$adj.list
        if (is.null(adj.list) && !is.null(fit$graph)) {
            adj.list <- fit$graph$adj.list
        }
    }
    if (is.null(weight.list)) {
        weight.list <- fit$weight.list
        if (is.null(weight.list) && !is.null(fit$graph)) {
            weight.list <- fit$graph$weight.list
        }
    }
    if (!is.list(adj.list) || !is.list(weight.list) ||
        length(adj.list) != length(weight.list)) {
        stop("A parallel adj.list and weight.list must be supplied or stored in fit.")
    }

    raw.rho <- vapply(seq_along(adj.list), function(i) {
        weights <- as.double(weight.list[[i]])
        if (!length(weights)) {
            return(NA_real_)
        }
        valid <- weights[is.finite(weights) & weights > 0]
        if (!length(valid)) {
            return(NA_real_)
        }
        1 / min(valid)
    }, numeric(1L))
    if (anyNA(raw.rho)) {
        stop("Cannot estimate density for isolated vertices or invalid weights.")
    }

    density.fit <- refit.rdgraph.regression(fit, raw.rho)
    smoothed.rho <- as.numeric(density.fit$fitted.values)
    rho.range <- range(smoothed.rho)
    if (diff(rho.range) > 0) {
        (smoothed.rho - rho.range[1L]) / diff(rho.range)
    } else {
        rep(1, length(smoothed.rho))
    }
}
