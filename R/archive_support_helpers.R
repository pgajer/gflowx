.maxp_no_check <- function(i, nbrs, y, rho) {
    if (!length(nbrs)) {
        return(NA_real_)
    }
    total <- sum(rho[nbrs])
    if (total < .Machine$double.eps) {
        return(NA_real_)
    }
    sum(rho[nbrs[y[i] > y[nbrs]]]) / total
}

.minp_no_check <- function(i, nbrs, y, rho) {
    if (!length(nbrs)) {
        return(NA_real_)
    }
    total <- sum(rho[nbrs])
    if (total < .Machine$double.eps) {
        return(NA_real_)
    }
    sum(rho[nbrs[y[i] < y[nbrs]]]) / total
}

generate.eta.grid <- function(eigenvalues, filter.type,
                              n.candidates = 40L) {
    n.candidates <- max(5L, as.integer(n.candidates))
    eigenvalues <- sort(as.numeric(eigenvalues))
    eigenvalues <- eigenvalues[
        is.finite(eigenvalues) & eigenvalues > .Machine$double.eps
    ]
    if (!length(eigenvalues)) {
        return(seq(1e-3, 1, length.out = n.candidates))
    }
    lower <- 1 / max(eigenvalues)
    upper <- 1 / min(eigenvalues)
    if (!is.finite(lower) || !is.finite(upper) ||
        lower <= 0 || upper <= lower) {
        return(seq(1e-3, 1, length.out = n.candidates))
    }
    exp(seq(log(max(lower, 1e-6)), log(upper),
            length.out = n.candidates))
}

right.winsorize <- function(y, p = 0.01, verbose = FALSE) {
    if (!is.numeric(p) || length(p) != 1L || p < 0 || p >= 0.25) {
        stop("p must be in [0, 0.25).")
    }
    if (p == 0 || !length(y) || all(is.na(y))) {
        return(y)
    }
    threshold <- stats::quantile(y, 1 - p, na.rm = TRUE, names = FALSE)
    if (isTRUE(verbose)) {
        message("Right winsorization threshold: ", threshold)
    }
    y[!is.na(y) & y > threshold] <- threshold
    y
}

winsorize <- function(x, p = 0.01, verbose = FALSE) {
    if (!is.numeric(x)) {
        stop("x must be numeric.")
    }
    if (!is.numeric(p) || length(p) != 1L || p <= 0 || p >= 0.25) {
        stop("p must be in (0, 0.25).")
    }
    if (!length(x) || all(is.na(x))) {
        return(x)
    }
    thresholds <- stats::quantile(
        x, c(p, 1 - p), na.rm = TRUE, names = FALSE
    )
    if (isTRUE(verbose)) {
        message(
            "Winsorization interval: [",
            thresholds[1L], ", ", thresholds[2L], "]"
        )
    }
    x[!is.na(x) & x < thresholds[1L]] <- thresholds[1L]
    x[!is.na(x) & x > thresholds[2L]] <- thresholds[2L]
    x
}

# Archived diagnostic used only to count classical graph extrema and annotate
# probabilistic extrema. Basin/extrema construction remains in gflow.
compute.extrema.hop.nbhds <- function(adj.list, weight.list, y) {
    if (!is.list(adj.list) || !is.list(weight.list) ||
        length(adj.list) != length(weight.list) ||
        length(y) != length(adj.list)) {
        stop("adj.list, weight.list, and y must describe the same graph.")
    }
    y <- as.numeric(y)
    maxima <- vapply(seq_along(adj.list), function(i) {
        neighbors <- as.integer(adj.list[[i]])
        !length(neighbors) || all(y[i] >= y[neighbors])
    }, logical(1L))
    minima <- vapply(seq_along(adj.list), function(i) {
        neighbors <- as.integer(adj.list[[i]])
        !length(neighbors) || all(y[i] <= y[neighbors])
    }, logical(1L))

    max.vertices <- which(maxima)
    min.vertices <- which(minima)
    data <- rbind(
        data.frame(
            label = paste0("M", seq_along(max.vertices)),
            vertex = max.vertices,
            value = y[max.vertices],
            type = "max",
            hop_idx = ifelse(
                y[max.vertices] == max(y), Inf, 1
            )
        ),
        data.frame(
            label = paste0("m", seq_along(min.vertices)),
            vertex = min.vertices,
            value = y[min.vertices],
            type = "min",
            hop_idx = ifelse(
                y[min.vertices] == min(y), Inf, 1
            )
        )
    )
    rownames(data) <- NULL
    list(extrema_df = data)
}
