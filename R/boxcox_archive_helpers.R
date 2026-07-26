.gflowx.boxcox.transform <- function(y, lambda) {
    if (any(y <= 0 | !is.finite(y))) {
        stop("Box-Cox requires y > 0 and finite.")
    }
    if (abs(lambda) < 1e-8) {
        return(log(y))
    }
    (y^lambda - 1) / lambda
}

.gflowx.boxcox.loglik <- function(lambda, y, Xmat) {
    transformed <- .gflowx.boxcox.transform(y, lambda)
    fit <- stats::lm.fit(x = Xmat, y = transformed)
    rss <- sum(fit$residuals^2)
    n <- length(y)
    -(n / 2) * log(rss / n) + (lambda - 1) * sum(log(y))
}

.gflowx.boxcox.mle <- function(
    formula,
    data = environment(formula),
    lambdas = seq(-2, 2, by = 0.1),
    refine = TRUE
) {
    frame <- stats::model.frame(formula, data = data)
    y <- stats::model.response(frame)
    if (any(y <= 0 | !is.finite(y))) {
        stop("Box-Cox requires a positive, finite response.")
    }
    Xmat <- stats::model.matrix(attr(frame, "terms"), data = frame)
    loglik <- vapply(
        lambdas, .gflowx.boxcox.loglik, numeric(1), y = y, Xmat = Xmat
    )
    lambda <- lambdas[which.max(loglik)]
    if (isTRUE(refine)) {
        step <- if (length(lambdas) > 1L) max(diff(lambdas)) else 0.25
        interval <- c(
            max(min(lambdas), lambda - step),
            min(max(lambdas), lambda + step)
        )
        lambda <- stats::optimize(
            function(value) {
                -.gflowx.boxcox.loglik(value, y = y, Xmat = Xmat)
            },
            interval = interval
        )$minimum
    }
    list(lambda = lambda)
}
