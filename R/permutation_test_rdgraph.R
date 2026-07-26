#' Vertex-wise permutation test for rdgraph regression fitted values
#'
#' @description
#' Performs a vertex-wise permutation test for a fitted rdgraph regression model.
#' The response vector \code{y} is permuted across vertices \code{n.perms} times,
#' the model is refit using \code{refit.rdgraph.regression()}, and a p-value is
#' computed at each vertex by comparing the observed statistic to the permutation
#' null distribution.
#'
#' The default statistic is \eqn{T_v = |\hat y(v) - \bar y|}, where \eqn{\bar y}
#' is the global mean of the original \code{y}. This is appropriate for binary
#' outcomes (incidence) and yields a two-sided test against departures from the
#' global mean incidence at each vertex.
#'
#' @param fitted.model A fitted model object returned by
#'   \code{\link{fit.rdgraph.regression}}.
#'   Must contain \code{fitted.values}.
#'
#' @param y Numeric vector of length n (original response used to produce
#'   \code{fitted.model}). For binary incidence, use 0/1 (numeric or integer).
#'
#' @param n.perms Integer \eqn{\ge 1}. Number of permutations.
#'
#' @param y.vertices Optional integer vector of labeled vertex indices (1-based).
#'   If provided, only these vertices contribute to the refit objective; residuals
#'   at unlabeled vertices are returned as \code{NA}. Values of \code{y.new} at
#'   unlabeled vertices are ignored and may be \code{NA}.
#' 
#' @param subj.id Optional vector of length n giving subject IDs (e.g., woman IDs)
#'   for repeated-measures data. If not \code{NULL}, permutations are performed
#'   at the subject level (block permutation): the response label is permuted
#'   across unique subjects and then assigned to all vertices belonging to that
#'   subject. This preserves within-subject dependence and is recommended when
#'   multiple samples per subject share the same outcome.
#'
#'   The current implementation requires \code{y} to be constant within subject.
#'
#' @param adj.list Adjacency list
#'
#' @param with.cluster Set to TRUE to compute cluster p-values
#'
#' @param cluster.threshold.quantile Cluster quantile threshold.
#' 
#' @param per.column.gcv Logical. Passed to \code{\link{refit.rdgraph.regression}}.
#'   If TRUE, reselect eta by GCV for each permutation column.
#'
#' @param n.cores Integer. Passed to \code{\link{refit.rdgraph.regression}}.
#'
#' @param seed Integer or NULL. If not NULL, sets RNG seed for reproducibility.
#'
#' @param stat Character scalar specifying the global test statistic to compute
#'     in addition to the vertex-wise p-values. Must be one of \code{"l2"},
#'     \code{"l1"}, \code{"lp"}, \code{"log"}, \code{"logit"}, or \code{"max"}.
#'
#'   The global statistic is computed from the fitted field \eqn{\hat y} and
#'   the global mean incidence \eqn{\bar y}:
#'   \itemize{
#'     \item \code{"l1"}: \eqn{T_1(\hat y) = \sum_{v\in V} |\hat y(v) - \bar y|}
#'     \item \code{"l2"}: \eqn{T_2(\hat y) = \sum_{v\in V} (\hat y(v) - \bar y)^2}
#'     \item \code{"lp"}: \eqn{T_p(\hat y) = \sum_{v\in V} |\hat y(v) - \bar y|^p}
#'     \item \code{"log"}: \eqn{T_{2,\log}(\hat y) = \sum_{v\in V} (\log(\hat y_\epsilon(v)) - \log(\bar y_\epsilon))^2}
#'     \item \code{"logit"}: \eqn{T_{2,\mathrm{logit}}(\hat y) = \sum_{v\in V} (\mathrm{logit}(\hat y_\epsilon(v)) - \mathrm{logit}(\bar y_\epsilon))^2}
#'     \item \code{"max"}: \eqn{T_{\mathrm{max}}(\hat y) = \max_{v\in V} |\hat y(v) - \bar y|}
#'   }
#'
#'   For \code{"log"} and \code{"logit"}, fitted values are clipped to
#'   \eqn{[\epsilon, 1-\epsilon]} with \eqn{\epsilon = \max(10^{-6}, 0.5/n)}
#'   to avoid infinities when \eqn{\hat y(v)} is numerically 0 or 1.
#'
#'   The global permutation p-value is computed as
#'   \eqn{p = (1 + \#\{b: T^{(b)} \ge T^{\mathrm{obs}}\})/(B+1)}.
#'
#' @param lp.test.exponent Exponent in the Lp test. Default: 0.5.
#' 
#' @param two.sided Logical. If TRUE (default), use absolute-value statistic.
#'   If FALSE, uses one-sided statistic \eqn{\hat y(v) - \bar y} (upper-tail).
#'
#' @param return.perm.fits Logical. If TRUE, return the matrix of permuted
#'   fitted values (can be large). Default FALSE.
#'
#' @param return.perm.refit Logical. If TRUE, return the output list of the refit function. Default FALSE.
#'
#' @param verbose Logical. If TRUE, print progress messages. Default TRUE.
#'
#' @return A list with components:
#'   \describe{
#'     \item{\code{p.value}}{Numeric vector of length n, vertex-wise permutation p-values.}
#'     \item{\code{q.value}}{BH-adjusted p-values (FDR), same length as \code{p.value}.}
#'     \item{\code{stat.obs}}{Observed test statistic at each vertex.}
#'     \item{\code{stat.perm}}{Optional matrix of permutation statistics (n x n.perms) if requested.}
#'     \item{\code{y.mean}}{Global mean \eqn{\bar y} used for centering.}
#'     \item{\code{fitted.obs}}{Observed fitted values \eqn{\hat y}.}
#'     \item{\code{perm.refit}}{The refit object returned by \code{refit.rdgraph.regression}.}
#'     \item{\code{perm.fitted.values}}{Optional matrix of permuted fitted values if requested.}
#'     \item{\code{call}}{Matched call.}
#'   }
#'
#' @examples
#' \dontrun{
#' y.fit <- fit.rdgraph.regression(X0.smoothed, y = as.double(sptb.int.im), k = k.opt)
#' perm.res <- permutation.test.rdgraph(
#'   fitted.model = y.fit,
#'   y = as.double(sptb.int.im),
#'   n.perms = 200,
#'   per.column.gcv = TRUE,
#'   n.cores = 14,
#'   seed = 1
#' )
#' summary(perm.res$p.value)
#' }
#'
#' @export
permutation.test.rdgraph <- function(fitted.model,
                                     y,
                                     n.perms,
                                     y.vertices = NULL,
                                     subj.id = NULL,
                                     adj.list = NULL,
                                     with.cluster = FALSE,
                                     cluster.threshold.quantile = 0.95,
                                     per.column.gcv = TRUE,
                                     n.cores = 1L,
                                     seed = 12345L,
                                     stat = c("l2","l1","lp","log","logit","max"),
                                     lp.test.exponent = 0.5,
                                     two.sided = TRUE,
                                     return.perm.fits = FALSE,
                                     return.perm.refit = FALSE,
                                     verbose = TRUE) {
    ## ---------------------------------------------------------------------
    ## Input validation
    ## ---------------------------------------------------------------------
    if (is.null(fitted.model)) stop("fitted.model cannot be NULL.")
    if (is.null(fitted.model$fitted.values)) {
        stop("fitted.model must contain fitted.values.")
    }

    if (!is.numeric(y)) stop("y must be numeric.")
    y <- as.double(y)

    n <- length(y)
    if (n < 2) stop("y must have length >= 2.")

    fitted.obs <- fitted.model$fitted.values
    if (!is.numeric(fitted.obs)) stop("fitted.model$fitted.values must be numeric.")
    if (length(fitted.obs) != n) {
        stop("Length mismatch: length(y) != length(fitted.model$fitted.values).")
    }

    if (!is.numeric(n.perms) || length(n.perms) != 1 || n.perms < 1 || n.perms != floor(n.perms)) {
        stop("n.perms must be a positive integer.")
    }
    n.perms <- as.integer(n.perms)

    if (!is.null(subj.id)) {
        if (length(subj.id) != n) {
            stop("subj.id must have length equal to length(y).")
        }
        if (anyNA(subj.id)) {
            stop("subj.id contains NA; please impute/remove or treat missing IDs explicitly.")
        }

        n.subj.tmp <- length(unique(as.character(subj.id)))
        if (n.subj.tmp < 2L) stop("Need at least 2 unique subjects for subj.id permutation.")
    }
    
    if (!is.logical(per.column.gcv) || length(per.column.gcv) != 1) {
        stop("per.column.gcv must be TRUE/FALSE.")
    }

    if (!is.numeric(n.cores) || length(n.cores) != 1 || n.cores < 1 || n.cores != floor(n.cores)) {
        stop("n.cores must be a positive integer.")
    }
    n.cores <- as.integer(n.cores)

    if (!is.logical(two.sided) || length(two.sided) != 1) {
        stop("two.sided must be TRUE/FALSE.")
    }
    if (!is.logical(return.perm.fits) || length(return.perm.fits) != 1) {
        stop("return.perm.fits must be TRUE/FALSE.")
    }
    if (!is.logical(return.perm.refit) || length(return.perm.refit) != 1) {
        stop("return.perm.refit must be TRUE/FALSE.")
    }
    if (!is.logical(verbose) || length(verbose) != 1) {
        stop("verbose must be TRUE/FALSE.")
    }

    if (!is.null(seed)) {
        if (!is.numeric(seed) || length(seed) != 1 || seed != floor(seed)) {
            stop("seed must be an integer or NULL.")
        }
        seed <- as.integer(seed)
        set.seed(seed)
    }

    stat <- match.arg(stat)

    ## ---------------------------------------------------------------------
    ## Build permutation matrix (n x n.perms)
    ## ---------------------------------------------------------------------
    if (verbose) {
        message(sprintf("Generating %d permutations for n = %d vertices...", n.perms, n))
    }

    n.subj <- NULL
    if (is.null(subj.id)) {

        ## ------------------------------------------------------------
        ## Vertex-level permutations (original behavior)
        ## ------------------------------------------------------------
        perm.idx <- replicate(n.perms, sample.int(n, size = n, replace = FALSE))
        perm.y <- matrix(y[perm.idx], nrow = n, ncol = n.perms)

    } else {

        ## ------------------------------------------------------------
        ## Subject-level (block) permutations
        ## ------------------------------------------------------------
        subj <- as.character(subj.id)
        subj.levels <- unique(subj)
        n.subj <- length(subj.levels)

        ## Map each vertex to a subject index in 1..n.subj
        subj.idx <- match(subj, subj.levels)
        if (anyNA(subj.idx)) stop("Internal error: NA in subj.idx mapping.")

        ## Require y be constant within subject (typical for woman-level outcomes)
        ## If you want to allow within-subject variation, block-permutation needs a different design.

        y.by.subj <- tapply(y, subj, function(z) z[1])
        y.by.subj <- as.double(y.by.subj[subj.levels])

        y.min <- tapply(y, subj, min)
        y.max <- tapply(y, subj, max)
        y.min <- as.double(y.min[subj.levels])
        y.max <- as.double(y.max[subj.levels])
        
        if (any(y.min != y.max)) {
            stop("y is not constant within subject; block permutation is ambiguous. ",
                 "Either aggregate to subject-level first or revise the permutation scheme.")
        }

        ## Subject permutation indices: n.subj x n.perms
        perm.subj.idx <- replicate(n.perms, sample.int(n.subj, size = n.subj, replace = FALSE))

        ## Expand to vertex-level permuted y: n x n.perms
        idx.mat <- perm.subj.idx[subj.idx, , drop = FALSE]  ## n x n.perms integer matrix
        perm.y <- matrix(as.double(y.by.subj)[idx.mat], nrow = n, ncol = n.perms)
    }
    
    ## ---------------------------------------------------------------------
    ## Refit using cached spectral structure
    ## ---------------------------------------------------------------------
    if (verbose) {
        message("Refitting permuted responses with refit.rdgraph.regression()...")
    }

    perm.refit <- refit.rdgraph.regression(
        fitted.model = fitted.model,
        y.new = perm.y,
        y.vertices = y.vertices,
        per.column.gcv = per.column.gcv,
        n.cores = n.cores,
        verbose = verbose
    )

    perm.fitted <- perm.refit$fitted.values
    if (is.null(perm.fitted)) stop("refit.rdgraph.regression returned NULL fitted.values.")
    if (!is.matrix(perm.fitted) || nrow(perm.fitted) != n || ncol(perm.fitted) != n.perms) {
        stop("Unexpected shape of permuted fitted.values; expected n x n.perms matrix.")
    }

    ## ---------------------------------------------------------------------
    ## Compute vertex-wise permutation p-values
    ## ---------------------------------------------------------------------

    y.mean.obs <- mean(y)
    y.mean.perm <- colMeans(perm.y)
    
    ## ---------------------------------------------------------------------
    ## Global permutation test(s): T2 on probability / log / logit scale
    ## ---------------------------------------------------------------------
    ## Note: We compute one global statistic chosen by `stat`, using the same
    ##       permuted fitted fields already available in `perm.fitted`.
    ##       This is independent of the vertex-wise statistic controlled by `two.sided`.
    eps <- max(1e-6, 0.5 / n)

    ## helper: clip to [eps, 1-eps] while preserving dimensions
    clip01 <- function(x, eps) {
        if (is.matrix(x) || length(dim(x)) == 2L || !is.null(dim(x))) {
            x2 <- x
            x2[x2 < eps] <- eps
            x2[x2 > 1 - eps] <- 1 - eps
            return(x2)
        }

        ## vector
        pmin(1 - eps, pmax(eps, x))
    }

    y.mean.obs.eps <- clip01(y.mean.obs, eps)

    if (stat == "l2") {
        ## T2 = sum (yhat - ybar)^2
        global.stat.obs <- sum((fitted.obs - y.mean.obs)^2)
        global.stat.perm <- colSums(sweep(perm.fitted, 2, y.mean.perm, FUN = "-")^2)
    } else if (stat == "l1") {
        global.stat.obs <- sum(abs(fitted.obs - y.mean.obs))
        global.stat.perm <- colSums(abs(sweep(perm.fitted, 2, y.mean.perm, FUN = "-")))
    } else if (stat == "lp") {
        if (!is.numeric(lp.test.exponent) || length(lp.test.exponent) != 1L ||
            !is.finite(lp.test.exponent) || lp.test.exponent <= 0) {
            stop("lp.test.exponent must be a finite numeric scalar > 0.")
        }
        global.stat.obs <- sum(abs(fitted.obs - y.mean.obs)^lp.test.exponent)
        global.stat.perm <- colSums(abs(sweep(perm.fitted, 2, y.mean.perm, FUN = "-"))^lp.test.exponent)

    } else if (stat == "log") {
        fitted.obs.eps <- clip01(fitted.obs, eps)
        y.mean.obs.eps <- clip01(y.mean.obs, eps)

        global.stat.obs <- sum((log(fitted.obs.eps) - log(y.mean.obs.eps))^2)

        y.mean.perm.eps <- clip01(y.mean.perm, eps)
        lp <- log(clip01(perm.fitted, eps))
        global.stat.perm <- colSums(sweep(lp, 2, log(y.mean.perm.eps), FUN = "-")^2)

    } else if (stat == "logit") {
        logit <- function(p) log(p / (1 - p))

        fitted.obs.eps <- clip01(fitted.obs, eps)
        y.mean.obs.eps <- clip01(y.mean.obs, eps)
        global.stat.obs <- sum((logit(fitted.obs.eps) - logit(y.mean.obs.eps))^2)

        y.mean.perm.eps <- clip01(y.mean.perm, eps)
        lp <- logit(clip01(perm.fitted, eps))
        global.stat.perm <- colSums(sweep(lp, 2, logit(y.mean.perm.eps), FUN = "-")^2)
        
    } else if (stat == "max") {
        ## T2 = max |yhat - ybar|
        global.stat.obs <- max(abs(fitted.obs - y.mean.obs))
        global.stat.perm <- apply(abs(sweep(perm.fitted, 2, y.mean.perm, FUN = "-")), 2, max)
    }

    ## Upper-tail p-value with +1 correction
    global.p.value <- (sum(global.stat.perm >= global.stat.obs) + 1) / (n.perms + 1)

    if (two.sided) {
        stat.obs <- abs(fitted.obs - y.mean.obs)
        ## subtract permutation-specific means columnwise
        stat.perm <- abs(sweep(perm.fitted, 2, y.mean.perm, FUN = "-"))
        ge.count <- rowSums(stat.perm >= stat.obs)
    } else {
        stat.obs <- fitted.obs - y.mean.obs
        stat.perm <- sweep(perm.fitted, 2, y.mean.perm, FUN = "-")
        ge.count <- rowSums(stat.perm >= stat.obs)
    }
    
    ## Phipson-Smyth style +1 correction for exactness
    p.value <- (ge.count + 1) / (n.perms + 1)

    ## Multiple testing correction (BH/FDR)
    q.value <- stats::p.adjust(p.value, method = "BH")

    rng.kind <- RNGkind()
    out <- list(
        n.subj = n.subj,
        y.mean.obs = y.mean.obs,
        y.mean.perm = y.mean.perm,
        
        p.value = p.value,
        q.value = q.value,
        stat.obs = stat.obs,

        stat = stat,
        global.stat.perm = global.stat.perm,
        global.stat.obs = global.stat.obs,
        global.p.value = global.p.value,
        fitted.obs = fitted.obs,

        eps = eps,
        n.perms = n.perms,
        two.sided = two.sided,
        seed = seed,
        lp.test.exponent = lp.test.exponent,
        per.column.gcv = per.column.gcv,
        with.cluster = with.cluster,
        cluster.threshold.quantile = cluster.threshold.quantile,
        rng.kind = rng.kind,
        call = match.call()
    )

    if (with.cluster) {
        if (is.null(adj.list)) stop("adj.list must be provided when with.cluster = TRUE.")
        if (!is.list(adj.list) || length(adj.list) != n) stop("adj.list must be a list of length n.")

        stat.obs.v <- stat.obs
        stat.perm.mat <- stat.perm

        threshold <- stats::quantile(as.double(stat.perm.mat),
                                     probs = cluster.threshold.quantile,
                                     na.rm = TRUE,
                                     names = FALSE)

        max.size.obs <- cluster.max.size(adj.list, stat.obs.v, threshold)
        max.size.perm <- apply(stat.perm.mat, 2, function(x) cluster.max.size(adj.list, x, threshold))
        cluster.size.p.value <- (sum(max.size.perm >= max.size.obs) + 1) / (n.perms + 1)

        max.mass.obs <- cluster.max.mass(adj.list, stat.obs.v, threshold)
        max.mass.perm <- apply(stat.perm.mat, 2, function(x) cluster.max.mass(adj.list, x, threshold))
        cluster.mass.p.value <- (sum(max.mass.perm >= max.mass.obs) + 1) / (n.perms + 1)

        out$threshold <- threshold
        
        out$max.mass.obs <- max.mass.obs
        out$max.mass.perm <- max.mass.perm
        out$cluster.mass.p.value <- cluster.mass.p.value

        out$max.size.obs <- max.size.obs
        out$max.size.perm <- max.size.perm
        out$cluster.size.p.value <- cluster.size.p.value
    }
    
    if (return.perm.fits) {
        out$perm.fitted.values <- perm.fitted
        out$stat.perm <- stat.perm
    }

    if (return.perm.refit) {
        out$perm.refit <- perm.refit
    }

    class(out) <- c("rdgraph_permutation_test", "list")
    return(out)
}

##' Maximum connected component size above a threshold on a graph
##'
##' @description
##' Given a graph adjacency list and a per-vertex statistic \code{stat.v}, this
##' function thresholds vertices by \code{stat.v >= threshold}, computes
##' connected components in the induced subgraph, and returns the size of the
##' largest component (cluster).
##'
##' This is intended for cluster-based permutation inference, where the observed
##' field and each permuted field produce a set of suprathreshold vertices and
##' one compares the observed maximum cluster size to the permutation null
##' distribution of maximum cluster sizes.
##'
##' @param adj.list List of length n. Each element is an integer vector of
##'   neighbor vertex indices in \code{1..n}. The graph is treated as undirected
##'   for connectivity (i.e., connectivity is evaluated using the provided
##'   neighbor links; if the underlying graph is directed, supply a symmetrized
##'   adjacency list).
##'
##' @param stat.v Numeric vector of length n giving a per-vertex statistic
##'   (e.g., \code{fitted.obs - y.mean} for one-sided tests, or
##'   \code{abs(fitted.obs - y.mean)} for two-sided tests).
##'
##' @param threshold Numeric scalar. Vertices with \code{stat.v >= threshold}
##'   are included in the suprathreshold induced subgraph.
##'
##' @return Integer scalar. Size (number of vertices) of the largest connected
##'   component among suprathreshold vertices. Returns 0L if no vertex exceeds
##'   \code{threshold}.
##'
##' @examples
##' \dontrun{
##' ## observed max cluster size above a threshold c
##' c <- stats::quantile(as.double(stat.perm.mat), probs = 0.95, na.rm = TRUE)
##' max.size.obs <- cluster.max.size(adj.list, stat.obs.v, threshold = c)
##'
##' ## permutation null of max cluster sizes
##' max.size.perm <- apply(stat.perm.mat, 2, function(x) {
##'     cluster.max.size(adj.list, x, threshold = c)
##' })
##'
##' ## p-value (upper tail) with +1 correction
##' p.cluster <- (sum(max.size.perm >= max.size.obs) + 1) / (length(max.size.perm) + 1)
##' }
cluster.max.size <- function(adj.list, stat.v, threshold) {
    ## ------------------------------------------------------------
    ## Input validation
    ## ------------------------------------------------------------
    if (!is.list(adj.list)) stop("adj.list must be a list.")
    n <- length(adj.list)

    if (!is.numeric(stat.v)) stop("stat.v must be numeric.")
    stat.v <- as.double(stat.v)
    if (length(stat.v) != n) {
        stop("Length mismatch: length(stat.v) must equal length(adj.list).")
    }

    if (!is.numeric(threshold) || length(threshold) != 1L || !is.finite(threshold)) {
        stop("threshold must be a finite numeric scalar.")
    }

    ## ------------------------------------------------------------
    ## Identify suprathreshold vertices
    ## ------------------------------------------------------------
    keep <- which(stat.v >= threshold)
    if (length(keep) == 0L) return(0L)

    ## visited: TRUE means "do not visit / already handled"
    ## mark all non-keep vertices as visited so BFS stays inside keep-set
    visited <- rep(FALSE, n)
    visited[-keep] <- TRUE

    max.size <- 0L

    ## ------------------------------------------------------------
    ## BFS/DFS over induced subgraph to find component sizes
    ## ------------------------------------------------------------
    for (s in keep) {
        if (visited[s]) next

        ## start a new component
        visited[s] <- TRUE
        queue <- s
        comp.size <- 0L

        while (length(queue) > 0L) {
            v <- queue[[1]]
            queue <- queue[-1]

            comp.size <- comp.size + 1L

            nb <- adj.list[[v]]
            if (length(nb) == 0L) next

            ## keep only neighbors not yet visited
            nb <- nb[!visited[nb]]

            if (length(nb) > 0L) {
                visited[nb] <- TRUE
                queue <- c(queue, nb)
            }
        }

        if (comp.size > max.size) max.size <- comp.size
    }

    max.size
}

cluster.max.mass <- function(adj.list, stat.v, threshold) {
    ## stat.v: numeric vector length n
    ## threshold: numeric scalar; keep vertices with stat.v >= threshold
    n <- length(stat.v)
    keep <- which(stat.v >= threshold)
    if (length(keep) == 0L) return(0)

    ## Build induced subgraph components on keep-set
    ## Use a fast BFS/DFS on adj.list to avoid igraph dependency
    visited <- rep(FALSE, n)
    visited[-keep] <- TRUE

    max.mass <- 0
    for (s in keep) {
        if (visited[s]) next
        ## BFS
        queue <- s
        visited[s] <- TRUE
        comp.mass <- 0
        while (length(queue) > 0L) {
            v <- queue[[1]]
            queue <- queue[-1]
            comp.mass <- comp.mass + stat.v[v]
            nb <- adj.list[[v]]
            nb <- nb[!visited[nb]]
            if (length(nb) > 0L) {
                visited[nb] <- TRUE
                queue <- c(queue, nb)
            }
        }
        if (comp.mass > max.mass) max.mass <- comp.mass
    }
    max.mass
}

#' Create a compact audit record from a permutation test result
#'
#' Reduce the output of \code{permutation.test.rdgraph()} (or a compatible list)
#' to a compact "audit stub" suitable for saving to disk while preserving
#' reproducibility-critical information (e.g., observed statistic, permutation
#' null, test options, and RNG metadata).
#'
#' @param res A list, typically the return value of \code{permutation.test.rdgraph()}.
#'   The function is tolerant to missing components and will keep only the fields
#'   that are present and relevant.
#' @param keep.vertex Logical scalar. If \code{TRUE}, retain vertex-wise outputs
#'   when available (e.g., \code{p.value}, \code{q.value}, \code{stat.obs},
#'   \code{fitted.obs}). Default is \code{FALSE}.
#' @param keep.global.null Logical scalar. If \code{TRUE}, retain the permutation
#'   null for the global statistic (typically \code{global.stat.perm}) when present.
#'   This enables post-hoc diagnostics (null histograms, Monte Carlo uncertainty)
#'   and verification of the reported p-value without rerunning refits. Default is
#'   \code{TRUE}.
#' @param keep.cluster Logical scalar. If \code{TRUE}, retain cluster-level
#'   summaries when available (e.g., max cluster size/mass distributions and
#'   corresponding p-values). Default is \code{TRUE}.
#'
#' @return
#' A named list with class \code{c("rdgraph_permutation_test_audit", "list")}.
#' The audit record typically includes:
#' \itemize{
#'   \item \code{global.p.value}: global permutation p-value (if present in \code{res})
#'   \item \code{global.stat.obs}: observed global test statistic (if present)
#'   \item \code{global.stat.perm}: permutation null vector (optional; controlled by \code{keep.global.null})
#'   \item \code{n.perms}: number of permutations inferred from \code{global.stat.perm} (if present)
#'   \item \code{stat}, \code{y.mean.obs}, \code{eps}: common test descriptors (if present)
#'   \item \code{seed}, \code{two.sided}: extracted from \code{res$call} when available
#'   \item \code{rng.kind}: current \code{RNGkind()} at time of auditing
#'   \item \code{call}: call stored in \code{res} (if present)
#'   \item optional sublists \code{vertex} and/or \code{cluster} (controlled by flags)
#' }
#'
#' @details
#' This helper is intended to be called immediately after
#' \code{permutation.test.rdgraph()} to avoid saving large intermediate objects
#' (e.g., permutation refits or fitted-value matrices) when only summary inference
#' is needed. For strict reproducibility, ensure the permutation test itself is
#' run with an explicit seed, and save the audited object alongside the inputs
#' that define the test (e.g., graph, response vector, and relevant model settings).
#'
#' @examples
#' ## Minimal mock object with the most common fields
#' set.seed(1)
#' res <- list(
#'   global.p.value = 0.05,
#'   global.stat.obs = 1.23,
#'   global.stat.perm = rnorm(1000),
#'   stat = "L2",
#'   y.mean.obs = 0.4,
#'   eps = 0,
#'   call = quote(permutation.test.rdgraph(seed = 1, two.sided = TRUE))
#' )
#'
#' audit <- perm.test.audit(res)
#' audit$global.p.value
#' length(audit$global.stat.perm)
#'
#' audit.small <- perm.test.audit(res, keep.global.null = FALSE)
#' names(audit.small)
#'
#' @seealso \code{\link{RNGkind}}
#' @export
perm.test.audit <- function(res,
                            keep.vertex = FALSE,
                            keep.global.null = TRUE,
                            keep.cluster = TRUE) {
    stopifnot(is.list(res))

    out <- list(
        global.p.value = res$global.p.value,
        global.stat.obs = res$global.stat.obs,
        global.stat.perm = if (keep.global.null) res$global.stat.perm else NULL,
        n.perms = if (!is.null(res$global.stat.perm)) length(res$global.stat.perm) else NA_integer_,
        stat = res$stat,
        y.mean.obs = res$y.mean.obs,
        eps = res$eps,
        call = res$call
    )

    ## Pull reproducibility knobs out of the call (since your current output
    ## does not store them explicitly)
    if (!is.null(res$call$seed)) out$seed <- as.integer(res$call$seed)
    if (!is.null(res$call$two.sided)) out$two.sided <- isTRUE(res$call$two.sided)
    out$rng.kind <- RNGkind()

    if (keep.vertex) {
        out$vertex <- list(
            p.value = res$p.value,
            q.value = res$q.value,
            stat.obs = res$stat.obs,
            fitted.obs = res$fitted.obs
        )
    }

    if (keep.cluster && !is.null(res$cluster.size.p.value)) {
        out$cluster <- list(
            threshold = res$threshold,
            max.size.obs = res$max.size.obs,
            max.size.perm = res$max.size.perm,
            cluster.size.p.value = res$cluster.size.p.value,
            max.mass.obs = res$max.mass.obs,
            max.mass.perm = res$max.mass.perm,
            cluster.mass.p.value = res$cluster.mass.p.value
        )
    }

    class(out) <- c("rdgraph_permutation_test_audit", "list")
    out
}
