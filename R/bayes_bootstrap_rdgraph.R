#' Bayesian bootstrap uncertainty for rdgraph fitted values
#'
#' @description
#' Computes Bayesian bootstrap (Rubin-style) uncertainty summaries for the
#' rdgraph regression fitted values \eqn{\hat y}. The Bayesian bootstrap draws
#' random weights on vertices, forms weighted pseudo-responses, and refits
#' using \code{\link{refit.rdgraph.regression}} with cached spectral structure.
#'
#' This provides posterior-like credible intervals and tail probabilities for
#' \eqn{\hat y(v)} (or derived contrasts) under uncertainty in the empirical
#' distribution of observed vertices. It is distinct from a permutation test:
#' Bayesian bootstrap quantifies estimator uncertainty rather than testing an
#' exchangeability null.
#'
#' @param fitted.model A fitted model object returned by
#'   \code{\link{fit.rdgraph.regression}}. Must contain \code{fitted.values}.
#'
#' @param y Numeric vector of length n (original response).
#'
#' @param n.draws Integer \eqn{\ge 1}. Number of Bayesian bootstrap draws.
#'
#' @param per.column.gcv Logical. Passed to \code{\link{refit.rdgraph.regression}}.
#'   If FALSE (default), uses the original filter for all draws (linear map).
#'   If TRUE, reselects smoothing per draw (more variable; slower).
#'
#' @param n.cores Integer. Passed to \code{\link{refit.rdgraph.regression}}.
#'
#' @param seed Integer or NULL. If not NULL, sets RNG seed.
#'
#' @param stat Character scalar specifying the global Bayesian bootstrap
#'   statistic to compute from fitted fields. Must be one of \code{"l2"},
#'   \code{"log"}, or \code{"logit"}.
#'
#'   Let \eqn{\hat y} be the fitted field and \eqn{\bar y} the global mean of
#'   the original \code{y}. The global statistic is computed for the observed
#'   fit and for each Bayesian bootstrap draw:
#'   \itemize{
#'     \item \code{"l2"}: \eqn{T_2(\hat y) = \sum_{v\in V} (\hat y(v) - \bar y)^2}
#'     \item \code{"log"}: \eqn{T_{2,\log}(\hat y) = \sum_{v\in V} (\log(\hat y_\epsilon(v)) - \log(\bar y_\epsilon))^2}
#'     \item \code{"logit"}: \eqn{T_{2,\mathrm{logit}}(\hat y) = \sum_{v\in V} (\mathrm{logit}(\hat y_\epsilon(v)) - \mathrm{logit}(\bar y_\epsilon))^2}
#'   }
#'
#'   For \code{"log"} and \code{"logit"}, fitted values are clipped to
#'   \eqn{[\epsilon, 1-\epsilon]} with \eqn{\epsilon = \max(10^{-6}, 0.5/n)}
#'   to avoid infinities when fitted values are numerically 0 or 1.
#'
#' @param credible.level Numeric in (0, 1). Default 0.95.
#'
#' @param baseline Character specifying baseline for enrichment probabilities.
#'   One of \code{"mean"} (global mean of y) or \code{"none"}.
#'
#' @param return.draws Logical. If TRUE, returns the matrix of draws
#'   (n x n.draws). Default FALSE to save memory.
#'
#' @param verbose Logical. If TRUE, prints progress messages.
#'
#' @return A list of class \code{"rdgraph_bayes_bootstrap"} with:
#'   \describe{
#'     \item{\code{fitted.obs}}{Observed fitted values \eqn{\hat y}.}
#'     \item{\code{lower}}{Lower credible bound per vertex.}
#'     \item{\code{upper}}{Upper credible bound per vertex.}
#'     \item{\code{median}}{Posterior median per vertex.}
#'     \item{\code{sd}}{Posterior SD per vertex (Monte Carlo).}
#'     \item{\code{prob.gt.baseline}}{If \code{baseline="mean"}, posterior
#'       probability \eqn{\Pr(\hat y(v) > \bar y)} per vertex.}
#'     \item{\code{bb.pseudo.p}}{A two-sided BB tail-area summary for
#'       \eqn{|\hat y(v)-\bar y|} if \code{baseline="mean"}.}
#'     \item{\code{refit}}{The refit object from \code{\link{refit.rdgraph.regression}}.}
#'     \item{\code{draws}}{Optional matrix of fitted draws (n x n.draws).}
#'     \item{\code{call}}{Matched call.}
#'   }
#'
#' @export
bayes.bootstrap.rdgraph <- function(fitted.model,
                                   y,
                                   n.draws,
                                   per.column.gcv = FALSE,
                                   n.cores = 1L,
                                   seed = 12345L,
                                   stat = c("l2", "log", "logit"),
                                   credible.level = 0.95,
                                   baseline = c("mean", "none"),
                                   return.draws = FALSE,
                                   verbose = TRUE) {
    ## -----------------------------
    ## Input validation
    ## -----------------------------
    if (is.null(fitted.model)) stop("fitted.model cannot be NULL.")
    if (is.null(fitted.model$fitted.values)) stop("fitted.model must contain fitted.values.")

    if (!is.numeric(y)) stop("y must be numeric.")
    y <- as.double(y)

    n <- length(y)
    if (n < 2) stop("y must have length >= 2.")

    fitted.obs <- fitted.model$fitted.values
    if (!is.numeric(fitted.obs)) stop("fitted.model$fitted.values must be numeric.")
    if (length(fitted.obs) != n) stop("Length mismatch: length(y) != length(fitted.model$fitted.values).")

    if (!is.numeric(n.draws) || length(n.draws) != 1 || n.draws < 1 || n.draws != floor(n.draws)) {
        stop("n.draws must be a positive integer.")
    }
    n.draws <- as.integer(n.draws)

    if (!is.numeric(n.cores) || length(n.cores) != 1 || n.cores < 1 || n.cores != floor(n.cores)) {
        stop("n.cores must be a positive integer.")
    }
    n.cores <- as.integer(n.cores)

    if (!is.numeric(credible.level) || length(credible.level) != 1 ||
        credible.level <= 0 || credible.level >= 1) {
        stop("credible.level must be in (0, 1).")
    }

    baseline <- match.arg(baseline)

    if (!is.logical(return.draws) || length(return.draws) != 1) stop("return.draws must be TRUE/FALSE.")
    if (!is.logical(verbose) || length(verbose) != 1) stop("verbose must be TRUE/FALSE.")

    if (!is.null(seed)) {
        if (!is.numeric(seed) || length(seed) != 1 || seed != floor(seed)) stop("seed must be an integer or NULL.")
        set.seed(as.integer(seed))
    }

    stat <- match.arg(stat)

    ## -----------------------------
    ## Bayesian bootstrap weights
    ## Dirichlet(1,...,1) via Exp(1)/sum
    ## Rescale to mean 1 (sum to n): w.n = n*w
    ## Pseudo-response: y.bb = w.n * y
    ## -----------------------------
    if (verbose) message(sprintf("Generating %d Bayesian bootstrap draws (n = %d)...", n.draws, n))

    w.raw <- matrix(stats::rexp(n * n.draws, rate = 1), nrow = n, ncol = n.draws)
    w <- w.raw / rep(colSums(w.raw), each = n)       ## columns sum to 1
    w.n <- n * w                                     ## columns sum to n

    y.bb <- w.n * y                                  ## elementwise by row recycling

    ## -----------------------------
    ## Refit (cached spectral)
    ## -----------------------------
    if (verbose) message("Refitting Bayesian bootstrap pseudo-responses with refit.rdgraph.regression()...")

    refit <- refit.rdgraph.regression(
        fitted.model = fitted.model,
        y.new = y.bb,
        per.column.gcv = per.column.gcv,
        n.cores = n.cores,
        verbose = verbose
    )

    draws <- refit$fitted.values
    if (is.null(draws) || !is.matrix(draws) || nrow(draws) != n || ncol(draws) != n.draws) {
        stop("Unexpected shape of refit$fitted.values; expected n x n.draws matrix.")
    }

    ## -----------------------------
    ## Global Bayesian bootstrap statistic draws
    ## -----------------------------
    y.mean.global <- mean(y)

    eps <- max(1e-6, 0.5 / n)

    ## helper: clip to [eps, 1-eps] while preserving dimensions
    clip01 <- function(x, eps) {
        x2 <- x
        x2[x2 < eps] <- eps
        x2[x2 > 1 - eps] <- 1 - eps
        x2
    }

    if (stat == "l2") {
        ## Observed
        global.stat.obs <- sum((fitted.obs - y.mean.global)^2)
        ## Draws: column-wise sum
        global.stat.draws <- colSums((draws - y.mean.global)^2)
    } else if (stat == "log") {
        fitted.obs.eps <- clip01(fitted.obs, eps)
        draws.eps <- clip01(draws, eps)
        y.mean.eps <- clip01(y.mean.global, eps)

        global.stat.obs <- sum((log(fitted.obs.eps) - log(y.mean.eps))^2)

        ld <- log(draws.eps)  ## n x n.draws
        global.stat.draws <- colSums((ld - log(y.mean.eps))^2)
    } else if (stat == "logit") {
        logit <- function(p) log(p / (1 - p))

        fitted.obs.eps <- clip01(fitted.obs, eps)
        draws.eps <- clip01(draws, eps)
        y.mean.eps <- clip01(y.mean.global, eps)

        global.stat.obs <- sum((logit(fitted.obs.eps) - logit(y.mean.eps))^2)

        ld <- logit(draws.eps)
        global.stat.draws <- colSums((ld - logit(y.mean.eps))^2)
    }

    ## posterior summaries for the global statistic
    alpha <- (1 - credible.level) / 2
    global.stat.mean <- mean(global.stat.draws)
    global.stat.ci <- as.double(stats::quantile(
                                           global.stat.draws,
                                           probs = c(alpha, 1 - alpha),
                                           na.rm = TRUE,
                                           names = FALSE
                                       ))

    ## optional: BB tail-area analog vs observed global statistic
    global.prob.ge.obs <- mean(global.stat.draws >= global.stat.obs)

    ## -----------------------------
    ## Summaries
    ## -----------------------------
    probs <- c(alpha, 0.5, 1 - alpha)

    lower <- apply(draws, 1, stats::quantile, probs = probs[1], na.rm = TRUE, names = FALSE)
    median <- apply(draws, 1, stats::quantile, probs = probs[2], na.rm = TRUE, names = FALSE)
    upper <- apply(draws, 1, stats::quantile, probs = probs[3], na.rm = TRUE, names = FALSE)
    sd <- apply(draws, 1, stats::sd, na.rm = TRUE)

    out <- list(
        fitted.obs = as.double(fitted.obs),
        stat = stat,
        y.mean.global = y.mean.global,
        global.stat.obs = as.double(global.stat.obs),
        global.stat.draws.mean = as.double(global.stat.mean),
        global.stat.draws.ci = as.double(global.stat.ci),
        global.prob.ge.obs = as.double(global.prob.ge.obs),
        lower = as.double(lower),
        median = as.double(median),
        upper = as.double(upper),
        sd = as.double(sd),
        refit = refit,
        call = match.call()
    )

    if (baseline == "mean") {
        y.mean <- mean(y)
        prob.gt <- rowMeans(draws > y.mean)
        ## BB tail-area for |yhat - y.mean|
        stat.obs <- abs(fitted.obs - y.mean)
        stat.bb <- abs(draws - y.mean)
        tail.ge <- rowMeans(stat.bb >= stat.obs)
        tail.le <- rowMeans(stat.bb <= stat.obs)
        bb.pseudo.p <- 2 * pmin(tail.ge, tail.le)
        bb.pseudo.p <- pmin(bb.pseudo.p, 1)

        out$baseline <- "mean"
        out$y.mean <- y.mean
        out$prob.gt.baseline <- as.double(prob.gt)
        out$bb.pseudo.p <- as.double(bb.pseudo.p)
    } else {
        out$baseline <- "none"
    }

    if (return.draws) {
        out$draws <- draws
        out$global.stat.draws <- as.double(global.stat.draws)
    }

    class(out) <- c("rdgraph_bayes_bootstrap", "list")
    return(out)
}
