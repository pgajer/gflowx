pca.optimal.components <- function(X, variance.threshold = 0.99,
                                   max.components = NULL, center = TRUE,
                                   scale = FALSE) {
    if (is.null(max.components) || max.components > ncol(X)) {
        max.components <- ncol(X)
    }
    pca.result <- stats::prcomp(X, center = center, scale. = scale)
    cumulative.variance <-
        cumsum(pca.result$sdev^2) / sum(pca.result$sdev^2)
    n.components <- which(cumulative.variance >= variance.threshold)[1L]
    if (is.na(n.components)) {
        n.components <- length(cumulative.variance)
    }
    n.components <- min(n.components, max.components)
    list(
        n.components = n.components,
        variance.explained = cumulative.variance[n.components],
        cumulative.variance = cumulative.variance,
        pca.result = pca.result
    )
}

pca.project <- function(X, pca.result, n.components) {
    out <- pca.result$x[, seq_len(n.components), drop = FALSE]
    rownames(out) <- rownames(X)
    out
}
