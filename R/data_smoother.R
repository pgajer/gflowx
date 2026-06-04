utils::globalVariables("k.val")

#' Data smoothing via RD-graph regression with k selection by GCV
#'
#' @description
#' Smooths a multivariate dataset \code{X} over a sequence of geometrically pruned
#' intersection kNN (ikNN) graphs by (i) selecting a graph scale parameter \code{k}
#' using GCV on a proxy scalar response \code{y}, then (ii) refitting the learned
#' smoothing operator to all columns of \code{X}.
#'
#' The algorithm:
#' \enumerate{
#'   \item Build ikNN graphs for \code{k = kmin:kmax} using \code{create.iknn.graphs()}
#'         with geometric pruning (and optional quantile pruning).
#'   \item If all graphs in the k-range are disconnected, trim rows of \code{X} to the
#'         largest connected component using the graph at \code{k = kmin}, then rebuild
#'         the graph sequence.
#'   \item Compute the terminal connected tail threshold \code{k.cc}, defined as the
#'         smallest k such that all graphs for k' >= k are connected. If \code{k.cc > kmin},
#'         restrict k-candidates to \code{k.cc:kmax}.
#'   \item Choose a proxy response \code{y} from \code{X} according to \code{proxy.response}:
#'         either the column with maximum variance (\code{"max.variance"}) or the first
#'         principal component score (\code{"pc1"}).
#'   \item For each candidate k, fit \code{fit.rdgraph.regression(X, y, k)} and record
#'         the model's GCV at the optimal iteration. Select \code{k.best} minimizing GCV.
#'   \item Smooth the full matrix via \code{refit.rdgraph.regression(fit.best, X)}.
#' }
#'
#' @param X Numeric matrix (n x p). Rows are vertices/observations, columns are features.
#' @param kmin Integer >= 1. Minimum k for the ikNN graph sequence.
#' @param kmax Integer >= kmin. Maximum k for the ikNN graph sequence. Requires \code{nrow(X) > kmax}.
#'
#' @param max.path.edge.ratio.deviation.thld Numeric in \eqn{[0, 0.2)}. Geometric pruning parameter
#'   passed to \code{create.iknn.graphs()}.
#' @param path.edge.ratio.percentile Numeric in \eqn{[0, 1]}. Percentile threshold used in geometric
#'   pruning (passed to \code{create.iknn.graphs()}).
#' @param threshold.percentile Numeric in \eqn{[0, 0.5]}. Optional quantile-based edge-length pruning
#'   (passed to \code{create.iknn.graphs()}).
#' @param knn.cache.path Optional character scalar path to a binary kNN cache file
#'   passed to \code{create.iknn.graphs()}.
#' @param knn.cache.mode Character cache mode passed to \code{create.iknn.graphs()}:
#'   \itemize{
#'     \item \code{"none"}: always compute kNN; do not read/write cache.
#'     \item \code{"read"}: read kNN from cache only; if cache is missing/invalid, error.
#'     \item \code{"write"}: always compute kNN and atomically write/overwrite cache.
#'     \item \code{"readwrite"}: read valid cache when available; if cache file is missing,
#'       compute and write cache; if cache exists but is invalid, error.
#'   }
#'
#' @param n.cores Integer >= 1. Number of CPU cores to use for fitting across k values.
#'   If \code{n.cores > 1}, requires the \pkg{foreach} and \pkg{doParallel} packages.
#'
#' @param block.size Optional positive integer. Passed to
#'   \code{refit.rdgraph.regression()} to process response columns in blocks
#'   during the final refit step. Useful for large multi-response matrices to
#'   reduce peak memory use. Default \code{NULL} (no explicit block size).
#'
#' @param proxy.response Character. Proxy response used to select \code{k.best} by GCV.
#'   Options:
#'   \itemize{
#'     \item \code{"max.variance"}: use the column of \code{X} with the largest sample variance.
#'     \item \code{"pc1"}: use the first PCA score (projection of samples onto the first principal axis)
#'           computed by \code{stats::prcomp()}.
#'   }
#' @param pca.center Logical. If \code{proxy.response = "pc1"}, whether to center columns of \code{X}
#'   before PCA (passed to \code{stats::prcomp()}). Default \code{TRUE}.
#' @param pca.scale Logical. If \code{proxy.response = "pc1"}, whether to scale columns of \code{X}
#'   before PCA (passed to \code{stats::prcomp()}). Default \code{FALSE}.
#'
#' @param pca.dim Integer. Passed to \code{fit.rdgraph.regression()} (default 100).
#' @param variance.explained Numeric in (0,1]. Passed to \code{fit.rdgraph.regression()} (default 0.99).
#' @param max.iterations Integer. Passed to \code{fit.rdgraph.regression()} (default 10).
#' @param n.eigenpairs Integer. Passed to \code{fit.rdgraph.regression()} (default 10).
#' @param filter.type Character. Passed to \code{fit.rdgraph.regression()} (default \code{"heat_kernel"}).
#' @param t.scale.factor Numeric. Passed to \code{fit.rdgraph.regression()} (default 0.5).
#' @param beta.coef.factor Numeric. Passed to \code{fit.rdgraph.regression()} (default 0.1).
#' @param dense.fallback Character. Passed to \code{fit.rdgraph.regression()}
#'   to control dense direct-solver fallback in diffusion solves. Options are
#'   \code{"auto"} (default), \code{"never"}, and \code{"always"}.
#' @param use.counting.measure Logical. If TRUE, uses uniform vertex weights
#'     (counting measure). If FALSE, uses distance-based weights inversely
#'     proportional to local k-NN density: \eqn{w(x) = (\epsilon +
#'     d_k(x))^{-\alpha}}. Distance-based weights are useful when sampling
#'     density varies across the feature space. Default: TRUE.
#' @param verbose Logical. If TRUE, prints progress messages.
#'
#' @return A list with elements:
#' \itemize{
#'   \item \code{X.smoothed}: numeric matrix, the smoothed version of \code{X} (possibly trimmed if all k graphs were disconnected).
#'   \item \code{k.best}: integer, selected k minimizing GCV.
#'   \item \code{k.values.used}: integer vector of k values evaluated.
#'   \item \code{gcv}: numeric vector of GCV values aligned with \code{k.values.used}.
#'   \item \code{selected.column}: integer index of the proxy column if \code{proxy.response="max.variance"}, otherwise \code{NA}.
#'   \item \code{proxy.y}: numeric vector used as the proxy response.
#'   \item \code{fit.best}: fitted object returned by \code{fit.rdgraph.regression()} at \code{k.best}.
#'   \item \code{trimmed}: logical, TRUE if rows were trimmed to the largest connected component.
#'   \item \code{kept.rows}: integer vector of original row indices kept if \code{trimmed=TRUE}, otherwise \code{NULL}.
#' }
#'
#' @export
data.smoother <- function(
    X,
    kmin,
    kmax,
    ## selected create.iknn.graphs() parameters
    max.path.edge.ratio.deviation.thld = 0.1,
    path.edge.ratio.percentile = 0.5,
    threshold.percentile = 0,
    knn.cache.path = NULL,
    knn.cache.mode = c("none", "read", "write", "readwrite"),
    n.cores = 1L,
    block.size = NULL,
    proxy.response = c("pc1","max.variance"),
    pca.center = TRUE,
    pca.scale = FALSE,
    ## selected fit.rdgraph.regression() parameters (posterior + other excluded params kept at defaults)
    pca.dim = 100,
    variance.explained = 0.99,
    max.iterations = 50,
    n.eigenpairs = 50,
    filter.type = "heat_kernel",
    t.scale.factor = 0.5,
    beta.coef.factor = 0.1,
    dense.fallback = c("auto", "never", "always"),
    use.counting.measure = TRUE,
    verbose = FALSE
) {

    ## -------------------------------------------------------------------------
    ## Helpers
    ## -------------------------------------------------------------------------

    if (!requireNamespace("gflow", quietly = TRUE)) {
        stop("data.smoother() currently requires gflow for its legacy graph ",
             "construction and rdgraph-regression backend.", call. = FALSE)
    }

    .is.int.scalar <- function(x) is.numeric(x) && length(x) == 1 && is.finite(x) && x == floor(x)

    .normalize.knn.cache.path <- function(knn.cache.path, knn.cache.mode) {
        utils::getFromNamespace(".normalize.knn.cache.path", "gflow")(
            knn.cache.path,
            knn.cache.mode
        )
    }

    .create.iknn.graphs <- function(...) {
        gflow::create.iknn.graphs(...)
    }

    .largest.cc.indices <- function(adj.list) {
        ## Returns 1-based vertex indices in the largest connected component
        n <- length(adj.list)
        if (n == 0) return(integer(0))

        visited <- rep(FALSE, n)
        comp.id <- integer(n)
        comp.sizes <- integer(0)
        current.comp <- 0L

        for (v in seq_len(n)) {
            if (visited[v]) next
            current.comp <- current.comp + 1L

            ## BFS
            q <- v
            visited[v] <- TRUE
            comp.id[v] <- current.comp
            size <- 0L

            while (length(q) > 0) {
                u <- q[1]
                q <- q[-1]
                size <- size + 1L
                nbrs <- adj.list[[u]]
                if (length(nbrs) == 0) next
                for (w in nbrs) {
                    if (!visited[w]) {
                        visited[w] <- TRUE
                        comp.id[w] <- current.comp
                        q <- c(q, w)
                    }
                }
            }

            comp.sizes[current.comp] <- size
        }

        if (length(comp.sizes) == 0) return(integer(0))
        main.comp <- which.max(comp.sizes)
        which(comp.id == main.comp)
    }

    .extract.gcv <- function(fit.res) {
        ## Your preferred extractor:
        ## gcv.optimal[optimal.iteration]
        if (!is.null(fit.res$gcv) &&
            !is.null(fit.res$gcv$gcv.optimal) &&
            !is.null(fit.res$optimal.iteration)) {
            ii <- fit.res$optimal.iteration
            g <- fit.res$gcv$gcv.optimal
            if (length(g) >= ii && ii >= 1) return(as.numeric(g[ii]))
        }

        ## Common fallbacks
        if (!is.null(fit.res$gcv.optimal) && is.numeric(fit.res$gcv.optimal)) {
            return(as.numeric(fit.res$gcv.optimal[1]))
        }
        if (!is.null(fit.res$gcv) && is.numeric(fit.res$gcv)) {
            return(as.numeric(fit.res$gcv[1]))
        }
        NA_real_
    }

    ## -------------------------------------------------------------------------
    ## Validate inputs
    ## -------------------------------------------------------------------------

    if (!is.matrix(X)) {
        X <- try(as.matrix(X), silent = TRUE)
        if (inherits(X, "try-error")) stop("X must be a matrix or coercible to a numeric matrix.")
    }
    if (!is.numeric(X)) stop("X must be numeric.")
    if (any(!is.finite(X))) stop("X cannot contain NA/NaN/Inf.")

    if (!.is.int.scalar(kmin) || kmin < 1) stop("kmin must be a positive integer.")
    if (!.is.int.scalar(kmax) || kmax < kmin) stop("kmax must be an integer not smaller than kmin.")
    kmin <- as.integer(kmin)
    kmax <- as.integer(kmax)

    if (nrow(X) <= kmax) stop("nrow(X) must be greater than kmax.")

    if (!.is.int.scalar(n.cores) || n.cores < 1) stop("n.cores must be a positive integer.")
    n.cores <- as.integer(n.cores)
    if (!is.null(block.size)) {
        if (!.is.int.scalar(block.size) || block.size < 1) {
            stop("block.size must be NULL or a positive integer.")
        }
        block.size <- as.integer(block.size)
    }
    knn.cache.mode <- match.arg(knn.cache.mode)
    knn.cache.path <- .normalize.knn.cache.path(knn.cache.path, knn.cache.mode)

    ## Keep X as double
    if (!is.double(X)) storage.mode(X) <- "double"
    dense.fallback <- match.arg(dense.fallback)

    ## -------------------------------------------------------------------------
    ## Build graph sequence (geom pruned)
    ## -------------------------------------------------------------------------

    if (verbose) cat("Building ikNN graphs (geom pruned) for k in [", kmin, ", ", kmax, "]\n", sep = "")

    ## Keep PCA/pruning params aligned with fit.rdgraph.regression() to avoid connectivity mismatches.
    X.graphs <- .create.iknn.graphs(
        X,
        kmin = kmin,
        kmax = kmax,
        max.path.edge.ratio.deviation.thld = max.path.edge.ratio.deviation.thld,
        path.edge.ratio.percentile = path.edge.ratio.percentile,
        threshold.percentile = threshold.percentile,
        compute.full = TRUE,
        pca.dim = pca.dim,
        variance.explained = variance.explained,
        knn.cache.path = knn.cache.path,
        knn.cache.mode = knn.cache.mode,
        n.cores = n.cores,
        verbose = isTRUE(verbose)
    )

    graphs <- X.graphs$geom_pruned_graphs
    if (is.null(graphs)) stop("create.iknn.graphs() did not return geom_pruned_graphs (compute.full=TRUE required).")

    X.graphs.stats <- summary(X.graphs)
    if (is.null(X.graphs.stats) || !("n_ccomp" %in% names(X.graphs.stats))) {
        stop("summary(create.iknn.graphs()) did not return an n_ccomp column.")
    }

    ## -------------------------------------------------------------------------
    ## Outlier trimming rule (only if all graphs are disconnected)
    ## -------------------------------------------------------------------------

    trimmed <- FALSE
    kept.rows <- seq_len(nrow(X))

    any.connected <- any(X.graphs.stats$n_ccomp == 1)

    if (!any.connected) {
        if (verbose) {
            cat("All graphs in k range have >1 connected component.\n")
            cat("Applying outlier trimming (largest CC) using graph at k = kmin.\n")
        }

        g0 <- graphs[[1]]
        if (is.null(g0$adj_list)) stop("geom_pruned_graphs[[1]] is missing adj_list.")

        keep.idx <- .largest.cc.indices(g0$adj_list)
        if (length(keep.idx) < 2) stop("Largest connected component after trimming is too small.")

        ## Trim X
        X <- X[keep.idx, , drop = FALSE]
        kept.rows <- kept.rows[keep.idx]
        trimmed <- TRUE

        ## Rebuild graphs after trimming
        if (verbose) cat("Rebuilding ikNN graph sequence (after trimming)\n")
        rebuild.knn.cache.path <- knn.cache.path
        rebuild.knn.cache.mode <- knn.cache.mode
        if (!is.null(rebuild.knn.cache.path) && rebuild.knn.cache.mode %in% c("read", "readwrite")) {
            if (verbose) {
                cat("Trimming changed matrix dimensions; disabling kNN cache for rebuilt graph sequence.\n")
            }
            rebuild.knn.cache.path <- NULL
            rebuild.knn.cache.mode <- "none"
        }
        X.graphs <- .create.iknn.graphs(
            X,
            kmin = kmin,
            kmax = kmax,
            max.path.edge.ratio.deviation.thld = max.path.edge.ratio.deviation.thld,
            path.edge.ratio.percentile = path.edge.ratio.percentile,
            threshold.percentile = threshold.percentile,
            compute.full = TRUE,
            pca.dim = pca.dim,
            variance.explained = variance.explained,
            knn.cache.path = rebuild.knn.cache.path,
            knn.cache.mode = rebuild.knn.cache.mode,
            n.cores = n.cores,
            verbose = isTRUE(verbose)
        )
        graphs <- X.graphs$geom_pruned_graphs
        if (is.null(graphs)) stop("After trimming, create.iknn.graphs() did not return geom_pruned_graphs.")
        X.graphs.stats <- summary(X.graphs)
    }

    ## -------------------------------------------------------------------------
    ## Compute terminal connected tail k.cc and restrict k range if needed
    ## -------------------------------------------------------------------------

    k.values.full <- kmin:kmax
    n.ccomp.vec <- X.graphs.stats$n_ccomp

    if (length(n.ccomp.vec) != length(k.values.full)) {
        stop("Length mismatch: X.graphs.stats$n_ccomp must match kmin:kmax")
    }

    disc.idx <- which(n.ccomp.vec > 1)

    if (length(disc.idx) == 0) {
        k.cc <- as.integer(kmin)
    } else if (max(disc.idx) == length(k.values.full)) {
        k.cc <- NA_integer_
    } else {
        k.cc <- as.integer(k.values.full[max(disc.idx) + 1L])
    }

    if (!is.na(k.cc)) {
        tail.idx <- which(k.values.full >= k.cc)
        if (any(n.ccomp.vec[tail.idx] > 1)) {
            stop("Internal error: computed k.cc does not define a connected tail")
        }
    }

    if (is.na(k.cc)) {
        ## No terminal connected tail; keep full range (but warn)
        if (verbose) cat("No terminal connected tail found in [kmin, kmax]; using full k range.\n")
        k.values.use <- k.values.full
    } else if (k.cc > kmin) {
        if (verbose) {
            cat("Found connected graphs starting at k = ", k.cc, "\n", sep = "")
            cat("Restricting stability analysis to k in [", k.cc, ", ", kmax, "]\n", sep = "")
        }
        k.values.use <- k.cc:kmax
    } else {
        k.values.use <- k.values.full
    }

    ## -------------------------------------------------------------------------
    ## Select proxy response
    ## -------------------------------------------------------------------------

    proxy.response <- match.arg(proxy.response)

    if (proxy.response == "max.variance") {
        var.vec <- apply(X, 2, stats::var)
        if (all(!is.finite(var.vec))) stop("All column variances are non-finite; cannot select proxy response.")
        col.idx <- which.max(var.vec)
        y <- as.double(X[, col.idx])

        if (verbose) {
            cat("Proxy response: max-variance column ", col.idx,
                " (var = ", signif(var.vec[col.idx], 4), ")\n", sep = "")
        }
    } else if (proxy.response == "pc1") {
        ## First PC score (projection of samples onto first principal axis)
        ## Note: prcomp() returns scores in $x (n x r); we use the first column.
        pca.res <- stats::prcomp(X, center = isTRUE(pca.center), scale. = isTRUE(pca.scale))

        if (is.null(pca.res$x) || ncol(pca.res$x) < 1L) {
            stop("PCA did not return PC scores; cannot use proxy.response='pc1'.")
        }

        y <- as.double(pca.res$x[, 1L])
        col.idx <- NA_integer_

        if (verbose) {
            ## Optional diagnostic: fraction variance explained by PC1
            sdev2 <- pca.res$sdev^2
            frac <- if (length(sdev2) > 0 && sum(sdev2) > 0) sdev2[1L] / sum(sdev2) else NA_real_
            cat("Proxy response: PC1 score (fraction variance explained = ",
                if (is.finite(frac)) signif(frac, 4) else "NA", ")\n", sep = "")
        }
    } else {
        stop("Unknown proxy.response value.")
    }

    ## -------------------------------------------------------------------------
    ## Fit rdgraph for each k and collect GCV
    ## -------------------------------------------------------------------------

    do.parallel <- (n.cores > 1L)

    fit.worker <- function(k.val) {
        fit.res <- fit.rdgraph.regression(
            X = X,
            y = y,
            k = as.integer(k.val),
            ## keep defaults for excluded params; enforce verbose.level = 0
            pca.dim = pca.dim,
            variance.explained = variance.explained,
            max.iterations = max.iterations,
            n.eigenpairs = n.eigenpairs,
            filter.type = filter.type,
            t.scale.factor = t.scale.factor,
            beta.coef.factor = beta.coef.factor,
            dense.fallback = dense.fallback,
            use.counting.measure = use.counting.measure,
            max.ratio.threshold = max.path.edge.ratio.deviation.thld,
            path.edge.ratio.percentile = path.edge.ratio.percentile,
            threshold.percentile = threshold.percentile,
            verbose.level = 0
        )
        list(
            k = as.integer(k.val),
            fit = fit.res,
            gcv = .extract.gcv(fit.res)
        )
    }

    if (verbose) cat("Fitting fit.rdgraph.regression() across k values and extracting GCV\n")

    res.list <- NULL

    if (do.parallel) {
        if (!requireNamespace("foreach", quietly = TRUE) ||
            !requireNamespace("doParallel", quietly = TRUE)) {
            warning("foreach/doParallel not available; running sequentially.")
            do.parallel <- FALSE
        }
    }

    if (do.parallel) {
        cl <- parallel::makeCluster(n.cores)
        on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)
        doParallel::registerDoParallel(cl)

        res.list <- foreach::`%dopar%`(foreach::foreach(
            k.val = k.values.use,
            .packages = character(0),
            .export = character(0)
        ), {
            fit.worker(k.val)
        })
    } else {
        res.list <- lapply(k.values.use, fit.worker)
    }

    gcv.vec <- vapply(res.list, function(z) z$gcv, numeric(1))
    k.vec <- vapply(res.list, function(z) z$k, integer(1))

    if (all(!is.finite(gcv.vec))) {
        warning("All GCV values are NA/non-finite; selecting the smallest k in k.values.use.")
        best.idx <- 1L
    } else {
        best.idx <- which.min(gcv.vec)
    }

    k.best <- k.vec[best.idx]
    fit.best <- res.list[[best.idx]]$fit

    if (verbose) {
        cat("Selected k = ", k.best, " (min GCV = ", signif(gcv.vec[best.idx], 6), ")\n", sep = "")
    }

    ## -------------------------------------------------------------------------
    ## Refit smoother to full matrix
    ## -------------------------------------------------------------------------

    refit.res <- refit.rdgraph.regression(
        fit.best,
        X,
        block.size = block.size
    )
    X.smoothed <- refit.res$fitted.values

    ## -------------------------------------------------------------------------
    ## Return richer object
    ## -------------------------------------------------------------------------

    out <- list(
        X.smoothed = X.smoothed,
        k.best = as.integer(k.best),
        k.values.used = as.integer(k.vec),
        gcv = as.numeric(gcv.vec),
        selected.column = as.integer(col.idx),
        proxy.y = y,
        fit.best = fit.best,
        trimmed = isTRUE(trimmed),
        kept.rows = if (isTRUE(trimmed)) as.integer(kept.rows) else NULL
    )

    out
}
