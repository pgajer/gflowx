## ============================================================================
## Local Slope (lslope) Association Test with Paired Bayesian Bootstrap
## ============================================================================
##
## This implementation uses PAIRED Bayesian bootstrap weights for signal and
## null samples, ensuring that any difference is purely due to the permutation
## destroying the y-z association, not due to different random weights.
##
## Key features:
## - Paired BB weights: same lambda used for signal and null
## - Three test options: paired.t, weighted.pvalue, wilcoxon
## - Automatic Box-Cox transformation with normality testing
## - Vertex-wise testing with FDR correction
##
## @author Pawel Gajer
## @date 2025
## ============================================================================


#' Local Slope Association Test with Paired Bayesian Bootstrap
#'
#' Perform vertex-wise hypothesis testing for local slope (lslope) associations
#' using paired Bayesian bootstrap. At each vertex, the test compares the
#' distribution of lslope values under signal (original y) versus null
#' (permuted y), using the same BB weights for both to ensure proper pairing.
#'
#' @param fitted.model A fitted model object from \code{gflowx::fit.rdgraph.regression()}
#'   containing the spectral decomposition and graph structure.
#' @param y Numeric vector: original (unsmoothed) response variable.
#' @param z Numeric vector: feature to test for association. Can be raw
#'   abundances or pre-smoothed values.
#' @param z.is.smoothed Logical. If \code{TRUE}, \code{z} is treated as already
#'   smoothed. If \code{FALSE} (default), \code{z} is smoothed using the fitted
#'   model's spectral decomposition before testing.
#' @param test.type Character string specifying the test for inference. Options:
#'   \describe{
#'     \item{"paired.t"}{Paired t-test on differences (signal - null). Most
#'       powerful when differences are approximately normal (default).}
#'     \item{"weighted.pvalue"}{Weighted p-value integrating over signal
#'       uncertainty. Assumes null distribution is approximately normal.}
#'     \item{"wilcoxon"}{Wilcoxon signed-rank test on differences.
#'       Distribution-free, most robust to non-normality.}
#'   }
#' @param boxcox Character string controlling Box-Cox transformation:
#'   \describe{
#'     \item{"auto"}{Apply Box-Cox if Shapiro-Wilk rejects normality (default)}
#'     \item{"always"}{Always apply Box-Cox transformation}
#'     \item{"never"}{Never apply Box-Cox transformation}
#'   }
#' @param boxcox.alpha Significance level for Shapiro-Wilk normality test when
#'   boxcox = "auto". Default is 0.05.
#' @param n.BB Integer: number of paired Bayesian bootstrap samples.
#'   Default is 500. Each BB sample generates one signal and one null lslope
#'   value using the same weights.
#' @param lslope.type Character: type of lslope measure. One of:
#'   \describe{
#'     \item{"normalized"}{Sigmoid-normalized slope (default, bounded \eqn{[-1,1]})}
#'     \item{"slope"}{Raw gradient slope Delta z / Delta y (unbounded)}
#'     \item{"sign"}{Sign of Delta z along gradient: -1, 0, or +1}
#'   }
#' @param ascending Logical. If \code{TRUE} (default), use ascending gradient
#'   direction (toward higher y values).
#' @param fdr.method Character: method for FDR correction across vertices.
#'   Default is "BH" (Benjamini-Hochberg). See \code{\link[stats]{p.adjust}}.
#' @param seed Integer: random seed for reproducibility. Default is NULL.
#' @param n.cores Integer: number of parallel cores. Default is 1.
#' @param verbose Logical: print progress messages. Default is TRUE.
#'
#' @return An object of class \code{"lslope.test"} containing:
#'   \describe{
#'     \item{p.value}{Numeric vector: raw p-value at each vertex}
#'     \item{p.adjusted}{Numeric vector: FDR-adjusted p-value at each vertex}
#'     \item{signal.mean}{Numeric vector: mean of signal lslope at each vertex}
#'     \item{signal.sd}{Numeric vector: SD of signal lslope at each vertex}
#'     \item{null.mean}{Numeric vector: mean of null lslope at each vertex}
#'     \item{null.sd}{Numeric vector: SD of null lslope at each vertex}
#'     \item{null.median}{Numeric vector: median of null lslope at each vertex}
#'     \item{null.mad}{Numeric vector: MAD of null lslope at each vertex (scaled)}
#'     \item{diff.mean}{Numeric vector: mean of (signal - null) at each vertex}
#'     \item{diff.sd}{Numeric vector: SD of (signal - null) at each vertex}
#'     \item{lslope.point}{Numeric vector: point estimate of lslope at each vertex}
#'     \item{lslope.z}{Numeric vector: normalized effect size (lslope - null.mean) / null.sd}
#'     \item{lslope.robust.z}{Numeric vector: robust normalized effect size using median/MAD}
#'     \item{shapiro.pvalue}{Numeric vector: Shapiro-Wilk p-value for differences}
#'     \item{boxcox.applied}{Logical vector: whether Box-Cox was applied}
#'     \item{boxcox.lambda}{Numeric vector: Box-Cox lambda (if applied)}
#'     \item{signal.samples}{Matrix (n.BB x n.vertices): signal lslope samples}
#'     \item{null.samples}{Matrix (n.BB x n.vertices): null lslope samples}
#'     \item{n.vertices}{Number of vertices}
#'     \item{n.BB}{Number of BB samples}
#'     \item{test.type}{Test type used}
#'     \item{fdr.method}{FDR method used}
#'   }
#'
#' @details
#' The test evaluates whether lslope(y.hat, z.hat) at each vertex reflects a
#' genuine association or could arise by chance. The paired BB approach:
#'
#' \enumerate{
#'   \item Generates shared Dirichlet weights lambda_1, ..., lambda_B
#'   \item For each b = 1, ..., B:
#'     \itemize{
#'       \item Signal: smooth y with weights lambda_b, compute lslope
#'       \item Null: permute y, smooth with SAME lambda_b, compute lslope
#'       \item Difference: d_b = signal_b - null_b
#'     }
#'   \item Test H0: \eqn{E[d] = 0} at each vertex using paired test
#' }
#'
#' The pairing ensures that any difference is due to the permutation breaking
#' the y-z association, not due to different random BB weights.
#'
#' @examples
#' \dontrun{
#' ## Fit a model
#' fit <- gflowx::fit.rdgraph.regression(X, y, k = 15)
#'
#' ## Test association with a single phylotype
#' result <- lslope.test(
#'     fitted.model = fit,
#'     y = y,
#'     z = Z[, "Lactobacillus_crispatus"],
#'     n.BB = 500,
#'     test.type = "paired.t",
#'     seed = 42
#' )
#'
#' ## Significant vertices (FDR < 0.05)
#' sig.vertices <- which(result$p.adjusted < 0.05)
#'
#' ## Plot results
#' plot(result)
#' }
#'
#' @seealso \code{\link{lslope}} for local slope computation,
#'   \code{\link{fassoc1.test}} for functional association testing
#'
#' @importFrom stats shapiro.test t.test wilcox.test pnorm sd p.adjust
#' @export
lslope.test <- function(fitted.model,
                        y,
                        z,
                        z.is.smoothed = FALSE,
                        test.type = c("paired.t", "weighted.pvalue", "wilcoxon"),
                        boxcox = c("auto", "always", "never"),
                        boxcox.alpha = 0.05,
                        n.BB = 500L,
                        lslope.type = c("normalized", "slope", "sign"),
                        ascending = TRUE,
                        fdr.method = "BH",
                        seed = NULL,
                        n.cores = 1L,
                        verbose = TRUE) {

    lslope.fun <- .gflowx.get.namespace.export(
        "gflow", "lslope", api = "lslope.test()",
        install_hint = "Install gflow to evaluate local-slope statistics."
    )
    weighted.p.fun <- .gflowx.get.namespace.export(
        "gflow", "weighted.p.value", api = "lslope.test()"
    )

    ## ========================================================================
    ## Input validation
    ## ========================================================================

    test.type <- match.arg(test.type)
    boxcox <- match.arg(boxcox)
    lslope.type <- match.arg(lslope.type)

    if (!inherits(fitted.model, "knn.riem.fit")) {
        stop("fitted.model must be of class 'knn.riem.fit' from gflowx::fit.rdgraph.regression()")
    }

    if (!is.numeric(y) || !is.numeric(z)) {
        stop("y and z must be numeric vectors")
    }

    n <- length(y)
    if (length(z) != n) {
        stop("y and z must have the same length")
    }

    if (!is.numeric(n.BB) || n.BB < 50) {
        stop("n.BB must be at least 50")
    }

    if (!is.numeric(boxcox.alpha) || boxcox.alpha <= 0 || boxcox.alpha >= 1) {
        stop("boxcox.alpha must be between 0 and 1")
    }

    n.BB <- as.integer(n.BB)
    n.cores <- as.integer(max(1, n.cores))

    if (!is.null(seed)) {
        set.seed(seed)
    }

    ## ========================================================================
    ## Extract model components
    ## ========================================================================

    if (verbose) {
        routine.ptm <- proc.time()
        message("Local Slope Association Test (Paired BB)")
        message(sprintf("  n = %d, n.BB = %d, test.type = %s", n, n.BB, test.type))
    }

    ## Extract spectral decomposition
    V <- fitted.model$V
    eigenvalues <- fitted.model$eigenvalues
    filter.weights <- fitted.model$filter.weights

    ## Extract graph structure
    adj.list <- fitted.model$adj.list
    weight.list <- fitted.model$weight.list

    ## ========================================================================
    ## Smooth z if needed
    ## ========================================================================

    if (!z.is.smoothed) {
        if (verbose) message("Smoothing z...")

        ## Apply spectral smoothing to z
        z.hat <- as.vector(V %*% (filter.weights * (t(V) %*% z)))
    } else {
        z.hat <- z
    }

    ## ========================================================================
    ## Compute point estimate of lslope
    ## ========================================================================

    if (verbose) message("Computing point estimate of lslope...")

    ## Smooth y for point estimate
    y.hat <- as.vector(V %*% (filter.weights * (t(V) %*% y)))

    ## Compute lslope point estimate
    lslope.point <- lslope.fun(
        adj.list, weight.list, y.hat, z.hat,
        type = lslope.type, ascending = ascending,
        instrumented = FALSE
    )

    ## ========================================================================
    ## Generate shared BB weights
    ## ========================================================================

    if (verbose) {
        message("Generating Dirichlet weights...")
        ptm <- proc.time()
    }

    lambda <- matrix(stats::rexp(n * n.BB), nrow = n, ncol = n.BB)
    lambda <- sweep(lambda, 2L, colSums(lambda), "/")

    if (verbose) {
        message(sprintf("  Done (%.2f sec)", (proc.time() - ptm)[3]))
    }

    ## ========================================================================
    ## Generate paired signal and null samples
    ## ========================================================================

    if (verbose) {
        message("Computing paired signal and null samples...")
        ptm <- proc.time()
    }

    ## Storage: rows = BB samples, columns = vertices
    signal.samples <- matrix(NA_real_, nrow = n.BB, ncol = n)
    null.samples <- matrix(NA_real_, nrow = n.BB, ncol = n)

    ## Generate one permutation to use for all null samples
    ## (Alternative: different permutation per BB sample)
    y.perm <- sample(y)

    for (b in seq_len(n.BB)) {
        ## Get weights for this BB sample
        w.b <- lambda[, b]

        ## ---- Signal: original y with BB weights ----

        ## Weighted spectral smoothing: y.hat^(w) = V * F * V^T * (w * y)
        ## where w is applied as observation weights
        y.weighted <- w.b * y
        y.hat.signal <- as.vector(V %*% (filter.weights * (t(V) %*% y.weighted)))

        ## Compute lslope for signal
        signal.samples[b, ] <- lslope.fun(
            adj.list, weight.list, y.hat.signal, z.hat,
            type = lslope.type, ascending = ascending,
            instrumented = FALSE
        )

        ## ---- Null: permuted y with SAME BB weights ----

        y.perm.weighted <- w.b * y.perm
        y.hat.null <- as.vector(V %*% (filter.weights * (t(V) %*% y.perm.weighted)))

        null.samples[b, ] <- lslope.fun(
            adj.list, weight.list, y.hat.null, z.hat,
            type = lslope.type, ascending = ascending,
            instrumented = FALSE
        )
    }

    if (verbose) {
        message(sprintf("  Done (%.2f sec)", (proc.time() - ptm)[3]))
    }

    ## ========================================================================
    ## Compute paired differences and perform vertex-wise tests
    ## ========================================================================

    if (verbose) {
        message("Performing vertex-wise hypothesis tests...")
        ptm <- proc.time()
    }

    ## Storage for results
    p.value <- numeric(n)
    shapiro.pvalue <- numeric(n)
    boxcox.applied <- logical(n)
    boxcox.lambda.vec <- rep(NA_real_, n)

    ## Summary statistics
    signal.mean <- colMeans(signal.samples, na.rm = TRUE)
    signal.sd <- apply(signal.samples, 2, sd, na.rm = TRUE)
    null.mean <- colMeans(null.samples, na.rm = TRUE)
    null.sd <- apply(null.samples, 2, sd, na.rm = TRUE)
    null.median <- apply(null.samples, 2, median, na.rm = TRUE)
    null.mad <- apply(null.samples, 2, function(x) mad(x, constant = 1.4826, na.rm = TRUE))

    ## Paired differences
    diff.samples <- signal.samples - null.samples
    diff.mean <- colMeans(diff.samples, na.rm = TRUE)
    diff.sd <- apply(diff.samples, 2, sd, na.rm = TRUE)

    ## Normalized effect sizes (per vertex)
    ## z-score: how many SDs is the point estimate from the null mean?
    lslope.z <- ifelse(null.sd > .Machine$double.eps,
                       (lslope.point - null.mean) / null.sd,
                       NA_real_)

    ## Robust z-score: using median and MAD
    lslope.robust.z <- ifelse(null.mad > .Machine$double.eps,
                              (lslope.point - null.median) / null.mad,
                              NA_real_)

    ## Test at each vertex
    for (v in seq_len(n)) {
        diff.v <- diff.samples[, v]

        ## Skip if all NA or constant
        if (all(is.na(diff.v)) || sd(diff.v, na.rm = TRUE) < .Machine$double.eps) {
            p.value[v] <- NA
            shapiro.pvalue[v] <- NA
            next
        }

        ## Test normality
        shapiro.result <- tryCatch(
            shapiro.test(diff.v),
            error = function(e) list(p.value = NA)
        )
        shapiro.pvalue[v] <- shapiro.result$p.value

        ## Decide on Box-Cox
        apply.bc <- FALSE
        if (boxcox == "always") {
            apply.bc <- TRUE
        } else if (boxcox == "auto") {
            if (!is.na(shapiro.pvalue[v]) && shapiro.pvalue[v] < boxcox.alpha) {
                apply.bc <- TRUE
            }
        }

        diff.for.test <- diff.v

        if (apply.bc) {
            ## Shift to positive if needed
            diff.for.bc <- diff.v
            shift <- 0
            if (min(diff.for.bc, na.rm = TRUE) <= 0) {
                shift <- abs(min(diff.for.bc, na.rm = TRUE)) +
                         0.01 * sd(diff.for.bc, na.rm = TRUE)
                diff.for.bc <- diff.for.bc + shift
            }

            ## Box-Cox transformation
            bc.fit <- tryCatch(
                .gflowx.boxcox.mle(diff.for.bc ~ 1),
                error = function(e) NULL
            )

            if (!is.null(bc.fit)) {
                boxcox.lambda.vec[v] <- bc.fit$lambda
                boxcox.applied[v] <- TRUE

                if (abs(bc.fit$lambda) > .Machine$double.eps) {
                    diff.for.test <- (diff.for.bc^bc.fit$lambda - 1) / bc.fit$lambda
                } else {
                    diff.for.test <- log(diff.for.bc)
                }
            }
        }

        ## Perform test
        if (test.type == "paired.t") {
            ## Paired t-test: H0: mean(diff) = 0, H1: mean(diff) != 0
            ## Use two-sided for lslope since association can be positive or negative
            t.result <- tryCatch(
                t.test(diff.for.test, mu = 0, alternative = "two.sided"),
                error = function(e) list(p.value = NA)
            )
            p.value[v] <- t.result$p.value

        } else if (test.type == "weighted.pvalue") {
            ## Weighted p-value approach
            ## Fit normal to null, compute weighted p-value over signal

            if (boxcox.applied[v]) {
                ## Transform null and signal
                null.v <- null.samples[, v]
                signal.v <- signal.samples[, v]

                shift <- 0
                if (min(null.v, na.rm = TRUE) <= 0) {
                    shift <- abs(min(null.v, na.rm = TRUE)) +
                             0.01 * sd(null.v, na.rm = TRUE)
                }

                null.for.test <- null.v + shift
                signal.for.test <- signal.v + shift

                if (abs(boxcox.lambda.vec[v]) > .Machine$double.eps) {
                    null.for.test <- (null.for.test^boxcox.lambda.vec[v] - 1) /
                                     boxcox.lambda.vec[v]
                    signal.for.test <- (signal.for.test^boxcox.lambda.vec[v] - 1) /
                                       boxcox.lambda.vec[v]
                } else {
                    null.for.test <- log(null.for.test)
                    signal.for.test <- log(signal.for.test)
                }
            } else {
                null.for.test <- null.samples[, v]
                signal.for.test <- signal.samples[, v]
            }

            mu <- mean(null.for.test, na.rm = TRUE)
            sigma <- sd(null.for.test, na.rm = TRUE)

            p.value[v] <- weighted.p.fun(
                signal.for.test, mu, sigma, alternative = "two.sided"
            )

        } else if (test.type == "wilcoxon") {
            ## Wilcoxon signed-rank test: H0: median(diff) = 0
            wilcox.result <- tryCatch(
                wilcox.test(diff.v, mu = 0, alternative = "two.sided"),
                error = function(e) list(p.value = NA)
            )
            p.value[v] <- wilcox.result$p.value
        }
    }

    ## FDR correction
    p.adjusted <- p.adjust(p.value, method = fdr.method)

    if (verbose) {
        message(sprintf("  Done (%.2f sec)", (proc.time() - ptm)[3]))

        n.sig.raw <- sum(p.value < 0.05, na.rm = TRUE)
        n.sig.fdr <- sum(p.adjusted < 0.05, na.rm = TRUE)
        message(sprintf("  Significant vertices: %d (raw), %d (FDR < 0.05)",
                        n.sig.raw, n.sig.fdr))

        message(sprintf("Total time: %.2f sec",
                        (proc.time() - routine.ptm)[3]))
    }

    ## ========================================================================
    ## Prepare output
    ## ========================================================================

    output <- list(
        p.value = p.value,
        p.adjusted = p.adjusted,
        signal.mean = signal.mean,
        signal.sd = signal.sd,
        null.mean = null.mean,
        null.sd = null.sd,
        null.median = null.median,
        null.mad = null.mad,
        diff.mean = diff.mean,
        diff.sd = diff.sd,
        lslope.point = lslope.point,
        lslope.z = lslope.z,
        lslope.robust.z = lslope.robust.z,
        shapiro.pvalue = shapiro.pvalue,
        boxcox.applied = boxcox.applied,
        boxcox.lambda = boxcox.lambda.vec,
        signal.samples = signal.samples,
        null.samples = null.samples,
        n.vertices = n,
        n.BB = n.BB,
        test.type = test.type,
        boxcox = boxcox,
        lslope.type = lslope.type,
        fdr.method = fdr.method,
        call = match.call()
    )

    class(output) <- "lslope.test"

    return(output)
}


## ============================================================================
## S3 Methods
## ============================================================================

#' Print Method for lslope.test Objects
#'
#' @param x An object of class "lslope.test".
#' @param digits Number of significant digits to display.
#' @param ... Additional arguments (currently ignored).
#'
#' @return Invisible x.
#'
#' @examples
#' result <- list(
#'   call = quote(lslope.test(fit, y, z)),
#'   n.vertices = 4L,
#'   n.BB = 8L,
#'   test.type = "paired",
#'   lslope.type = "forward",
#'   fdr.method = "BH",
#'   p.value = c(0.01, 0.20, 0.03, 0.50),
#'   p.adjusted = c(0.04, 0.27, 0.06, 0.50),
#'   lslope.point = c(0.8, -0.1, 0.5, 0.0),
#'   lslope.z = c(2.5, -0.4, 1.9, 0.0)
#' )
#' class(result) <- "lslope.test"
#'
#' print(result)
#' @method print lslope.test
#' @export
print.lslope.test <- function(x, digits = 4, ...) {
    cat("\nLocal Slope Association Test (Paired BB)\n")
    cat("=========================================\n\n")

    cat("Call:\n")
    print(x$call)
    cat("\n")

    cat("Configuration:\n")
    cat("  Vertices: ", x$n.vertices, "\n", sep = "")
    cat("  BB samples: ", x$n.BB, "\n", sep = "")
    cat("  Test type: ", x$test.type, "\n", sep = "")
    cat("  lslope type: ", x$lslope.type, "\n", sep = "")
    cat("  FDR method: ", x$fdr.method, "\n\n", sep = "")

    cat("Results Summary:\n")
    n.tested <- sum(!is.na(x$p.value))
    n.sig.05 <- sum(x$p.adjusted < 0.05, na.rm = TRUE)
    n.sig.01 <- sum(x$p.adjusted < 0.01, na.rm = TRUE)
    n.sig.001 <- sum(x$p.adjusted < 0.001, na.rm = TRUE)

    cat("  Vertices tested: ", n.tested, "\n", sep = "")
    cat("  Significant (FDR < 0.05): ", n.sig.05, "\n", sep = "")
    cat("  Significant (FDR < 0.01): ", n.sig.01, "\n", sep = "")
    cat("  Significant (FDR < 0.001): ", n.sig.001, "\n", sep = "")

    if (n.sig.05 > 0) {
        cat("\n  Top significant vertices:\n")
        top.idx <- head(order(x$p.adjusted), 10)
        for (i in top.idx) {
            if (!is.na(x$p.adjusted[i]) && x$p.adjusted[i] < 0.05) {
                cat(sprintf("    Vertex %d: p.adj = %.4g, lslope = %.4f, z = %.2f\n",
                            i, x$p.adjusted[i], x$lslope.point[i], x$lslope.z[i]))
            }
        }
    }

    invisible(x)
}


#' Summary Method for lslope.test Objects
#'
#' @param object An object of class "lslope.test".
#' @param ... Additional arguments (currently ignored).
#'
#' @return A list of class "summary.lslope.test" with summary statistics.
#'
#' @examples
#' result <- list(
#'   call = quote(lslope.test(fit, y, z)),
#'   n.vertices = 4L,
#'   n.BB = 8L,
#'   test.type = "paired",
#'   lslope.type = "forward",
#'   fdr.method = "BH",
#'   p.value = c(0.01, 0.20, 0.03, 0.50),
#'   p.adjusted = c(0.04, 0.27, 0.06, 0.50),
#'   lslope.point = c(0.8, -0.1, 0.5, 0.0),
#'   diff.mean = c(0.4, -0.1, 0.3, 0.0)
#' )
#' class(result) <- "lslope.test"
#'
#' summary(result)
#' @method summary lslope.test
#' @export
summary.lslope.test <- function(object, ...) {
    out <- list(
        call = object$call,
        n.vertices = object$n.vertices,
        n.BB = object$n.BB,
        test.type = object$test.type,
        lslope.type = object$lslope.type,
        fdr.method = object$fdr.method,
        n.tested = sum(!is.na(object$p.value)),
        n.sig.05 = sum(object$p.adjusted < 0.05, na.rm = TRUE),
        n.sig.01 = sum(object$p.adjusted < 0.01, na.rm = TRUE),
        n.sig.001 = sum(object$p.adjusted < 0.001, na.rm = TRUE),
        p.value.summary = summary(object$p.value),
        p.adjusted.summary = summary(object$p.adjusted),
        lslope.summary = summary(object$lslope.point),
        diff.mean.summary = summary(object$diff.mean)
    )

    class(out) <- "summary.lslope.test"
    return(out)
}


#' @method print summary.lslope.test
#' @examples
#' result <- list(
#'   call = quote(lslope.test(fit, y, z)),
#'   n.vertices = 4L,
#'   n.BB = 8L,
#'   test.type = "paired",
#'   lslope.type = "forward",
#'   fdr.method = "BH",
#'   p.value = c(0.01, 0.20, 0.03, 0.50),
#'   p.adjusted = c(0.04, 0.27, 0.06, 0.50),
#'   lslope.point = c(0.8, -0.1, 0.5, 0.0),
#'   diff.mean = c(0.4, -0.1, 0.3, 0.0)
#' )
#' class(result) <- "lslope.test"
#'
#' result.summary <- summary(result)
#' result.summary
#' @export
print.summary.lslope.test <- function(x, digits = 4, ...) {
    cat("\nLocal Slope Association Test Summary\n")
    cat("====================================\n\n")

    cat("Call:\n")
    print(x$call)
    cat("\n")

    cat("Configuration:\n")
    cat("  Vertices: ", x$n.vertices, "\n", sep = "")
    cat("  BB samples: ", x$n.BB, "\n", sep = "")
    cat("  Test type: ", x$test.type, "\n\n", sep = "")

    cat("Significance:\n")
    cat("  Tested: ", x$n.tested, "\n", sep = "")
    cat("  FDR < 0.05: ", x$n.sig.05,
        sprintf(" (%.1f%%)", 100 * x$n.sig.05 / x$n.tested), "\n", sep = "")
    cat("  FDR < 0.01: ", x$n.sig.01, "\n", sep = "")
    cat("  FDR < 0.001: ", x$n.sig.001, "\n\n", sep = "")

    cat("P-value Distribution:\n")
    print(signif(x$p.value.summary, digits))
    cat("\n")

    cat("lslope Point Estimate Distribution:\n")
    print(signif(x$lslope.summary, digits))

    invisible(x)
}


#' Coef Method for lslope.test Objects
#'
#' Extracts lslope point estimates from a lslope.test object.
#'
#' @param object An object of class "lslope.test".
#' @param ... Additional arguments (currently ignored).
#'
#' @return A named numeric vector of lslope point estimates.
#' @method coef lslope.test
#' @export
coef.lslope.test <- function(object, ...) {
    setNames(object$lslope.point, seq_along(object$lslope.point))
}


#' Plot Method for lslope.test Objects
#'
#' Creates diagnostic plots for lslope.test results.
#'
#' @param x An object of class "lslope.test".
#' @param type Character string specifying plot type. Options are:
#'   \describe{
#'     \item{"volcano"}{Volcano plot: lslope vs -log10(p.adjusted) (default)}
#'     \item{"manhattan"}{Manhattan plot: -log10(p.value) by vertex}
#'     \item{"qq"}{QQ plot of p-values}
#'     \item{"histogram"}{Histogram of p-values}
#'     \item{"diff"}{Histogram of mean differences at each vertex}
#'   }
#' @param fdr.threshold Numeric: FDR threshold for highlighting (default 0.05).
#' @param ... Additional arguments passed to plotting functions.
#'
#' @examples
#' result <- list(
#'   n.vertices = 4L,
#'   p.value = c(0.01, 0.20, 0.03, 0.50),
#'   p.adjusted = c(0.04, 0.27, 0.06, 0.50),
#'   lslope.point = c(0.8, -0.1, 0.5, 0.0),
#'   diff.mean = c(0.4, -0.1, 0.3, 0.0)
#' )
#' class(result) <- "lslope.test"
#'
#' plot(result, type = "volcano")
#' @method plot lslope.test
#' @export
plot.lslope.test <- function(x,
                             type = c("volcano", "manhattan", "qq",
                                      "histogram", "diff"),
                             fdr.threshold = 0.05,
                             ...) {
    type <- match.arg(type)

    if (type == "volcano") {
        ## Volcano plot
        neg.log.p <- -log10(pmax(x$p.adjusted, 1e-300))
        sig <- x$p.adjusted < fdr.threshold

        plot(x$lslope.point, neg.log.p,
             xlab = "lslope point estimate",
             ylab = expression(-log[10](p.adjusted)),
             main = "Volcano Plot",
             col = ifelse(sig, "red", "gray50"),
             pch = ifelse(sig, 19, 1),
             cex = ifelse(sig, 0.8, 0.5),
             las = 1, ...)

        abline(h = -log10(fdr.threshold), col = "blue", lty = 2)
        abline(v = 0, col = "gray", lty = 3)

        legend("topright",
               legend = c(sprintf("FDR < %.2f (n=%d)", fdr.threshold, sum(sig)),
                          "Not significant"),
               col = c("red", "gray50"),
               pch = c(19, 1), cex = 0.8, bty = "n")

    } else if (type == "manhattan") {
        ## Manhattan plot
        neg.log.p <- -log10(pmax(x$p.value, 1e-300))
        sig <- x$p.adjusted < fdr.threshold

        plot(seq_len(x$n.vertices), neg.log.p,
             xlab = "Vertex",
             ylab = expression(-log[10](p.value)),
             main = "Manhattan Plot",
             col = ifelse(sig, "red", "gray50"),
             pch = 19, cex = 0.5,
             las = 1, ...)

        abline(h = -log10(0.05), col = "blue", lty = 2)

    } else if (type == "qq") {
        ## QQ plot of p-values
        p.obs <- sort(x$p.value[!is.na(x$p.value)])
        n.p <- length(p.obs)
        p.exp <- (seq_len(n.p) - 0.5) / n.p

        plot(-log10(p.exp), -log10(p.obs),
             xlab = expression(Expected ~ -log[10](p)),
             ylab = expression(Observed ~ -log[10](p)),
             main = "QQ Plot of P-values",
             pch = 19, cex = 0.5, col = "darkblue",
             las = 1, ...)

        abline(0, 1, col = "red", lty = 2)

    } else if (type == "histogram") {
        ## Histogram of p-values
        hist(x$p.value, breaks = 20,
             main = "Distribution of P-values",
             xlab = "P-value",
             col = "lightblue", border = "white",
             las = 1, ...)

        abline(v = 0.05, col = "red", lty = 2)

    } else if (type == "diff") {
        ## Histogram of mean differences
        hist(x$diff.mean, breaks = 30,
             main = "Distribution of Mean Differences",
             xlab = "Mean(signal - null)",
             col = "lightgreen", border = "white",
             las = 1, ...)

        abline(v = 0, col = "red", lty = 2)
        abline(v = mean(x$diff.mean, na.rm = TRUE), col = "blue", lwd = 2)
    }

    invisible(NULL)
}


#' Extract Significant Vertices from lslope.test
#'
#' Convenience function to extract vertices with significant lslope associations.
#'
#' @param object An object of class "lslope.test".
#' @param fdr.threshold Numeric: FDR threshold (default 0.05).
#' @param return.details Logical: return full details or just indices (default FALSE).
#'
#' @return If return.details is FALSE, a vector of significant vertex indices.
#'   If TRUE, a data.frame with columns: vertex, lslope, lslope.z, lslope.robust.z,
#'   p.value, p.adjusted, diff.mean, diff.sd.
#'
#' @examples
#' \dontrun{
#' result <- lslope.test(fit, y, z)
#'
#' # Get significant vertex indices
#' sig.idx <- significant.vertices(result)
#'
#' # Get detailed results
#' sig.details <- significant.vertices(result, return.details = TRUE)
#' }
#'
#' @export
significant.vertices <- function(object,
                                  fdr.threshold = 0.05,
                                  return.details = FALSE) {

    if (!inherits(object, "lslope.test")) {
        stop("object must be of class 'lslope.test'")
    }

    sig.idx <- which(object$p.adjusted < fdr.threshold)

    if (!return.details) {
        return(sig.idx)
    }

    if (length(sig.idx) == 0) {
        return(data.frame(
            vertex = integer(0),
            lslope = numeric(0),
            lslope.z = numeric(0),
            lslope.robust.z = numeric(0),
            p.value = numeric(0),
            p.adjusted = numeric(0),
            diff.mean = numeric(0),
            diff.sd = numeric(0)
        ))
    }

    data.frame(
        vertex = sig.idx,
        lslope = object$lslope.point[sig.idx],
        lslope.z = object$lslope.z[sig.idx],
        lslope.robust.z = object$lslope.robust.z[sig.idx],
        p.value = object$p.value[sig.idx],
        p.adjusted = object$p.adjusted[sig.idx],
        diff.mean = object$diff.mean[sig.idx],
        diff.sd = object$diff.sd[sig.idx]
    )
}
