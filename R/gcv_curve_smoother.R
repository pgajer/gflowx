.gflowx.fit.gcv.smoother <- function(x, y, params = list()) {
    x <- as.numeric(x)
    y <- as.numeric(y)
    keep <- is.finite(x) & is.finite(y)
    x <- x[keep]
    y <- y[keep]
    if (length(x) < 4L || length(unique(x)) < 4L) {
        stop("At least four finite, distinct p values are required.")
    }

    grid.size <- params$grid.size
    if (is.null(grid.size)) {
        grid.size <- max(200L, 5L * length(unique(x)))
    }
    xgrid <- seq(min(x), max(x), length.out = as.integer(grid.size))
    fit.args <- list(x = x, y = y)
    if (!is.null(params$spar)) {
        fit.args$spar <- params$spar
    } else if (!is.null(params$df)) {
        fit.args$df <- params$df
    } else {
        fit.args$cv <- isTRUE(params$use.gcv)
    }
    fit <- do.call(stats::smooth.spline, fit.args)
    prediction <- stats::predict(fit, xgrid)$y

    residual.scale <- sqrt(mean((y - stats::predict(fit, x)$y)^2))
    interval <- 1.96 * residual.scale
    out <- list(
        xgrid = xgrid,
        gpredictions = prediction,
        gpredictions.CrI = rbind(prediction - interval, prediction + interval),
        xg = xgrid,
        Eyg = prediction,
        fit = fit,
        fit.method = "stats::smooth.spline",
        selected.spar = fit$spar,
        selected.df = fit$df
    )
    class(out) <- c("gflowx_spline_fit", "list")
    out
}
