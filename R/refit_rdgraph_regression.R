#' Refit Riemannian Graph Regression with New Response(s)
#'
#' Apply spectral filtering from a fitted model to new response variable(s).
#' Optionally, perform semi-supervised refitting when responses are observed
#' only on a subset of vertices (labeled set).
#'
#' This function uses the cached eigendecomposition from an existing fit,
#' avoiding the computational cost of rebuilding the geometric structure.
#'
#' When \code{y.vertices} is provided, the refit uses a masked (semi-supervised)
#' spectral Tikhonov formulation in the eigenbasis. Only labeled vertices
#' contribute to the data-fit term; fitted values are still returned for all
#' vertices.
#'
#' @param fitted.model A fitted model object of class \code{"knn.riem.fit"}
#'   returned by \code{\link{fit.rdgraph.regression}}. Must contain the
#'   \code{spectral} component with cached eigendecomposition, including
#'   \code{eigenvectors}, \code{eigenvalues}, and \code{eta.optimal}.
#'
#' @param y.new New response variable(s). Either a numeric vector of length n,
#'   or a numeric matrix-like object (base \code{matrix} or \code{Matrix}
#'   class) with n rows where each column is a separate response.
#'
#' @param y.vertices Optional integer vector of labeled vertex indices (1-based).
#'   If provided, only these vertices contribute to the refit objective; residuals
#'   at unlabeled vertices are returned as \code{NA}. Values of \code{y.new} at
#'   unlabeled vertices are ignored and may be \code{NA}.
#'
#' @param per.column.gcv Logical. If \code{TRUE}, select optimal smoothing
#'   parameter eta independently for each column via GCV. If \code{FALSE}
#'   (default), use the filter/eta from the original fit.
#'
#'   \strong{Restriction:} \code{per.column.gcv = TRUE} is not supported when
#'   \code{y.vertices} is provided.
#'
#' @param n.candidates Integer. Number of candidate eta values for GCV grid
#'   search. Only used when \code{per.column.gcv = TRUE}. Default is 40.
#'
#' @param n.cores Integer. Number of CPU cores for parallel processing.
#'   Default is 1 (sequential). Used only for per-column GCV.
#'
#' @param block.size Optional positive integer. When provided and
#'   \code{per.column.gcv = FALSE}, responses are processed in column blocks of
#'   this size to reduce peak memory pressure for large multi-response inputs.
#'   Ignored for single-response input.
#'
#' @param verbose Logical. If \code{TRUE}, print progress information. Default
#'   is \code{FALSE}.
#'
#' @param with.posterior Logical. (Kept for backward compatibility.) If \code{TRUE},
#'   compute posterior credible intervals (not supported when \code{y.vertices}
#'   is provided).
#'
#' @param return.posterior.samples Logical. Posterior option (see original
#'   implementation).
#'
#' @param credible.level Numeric in (0, 1). Posterior option.
#'
#' @param n.posterior.samples Integer. Posterior option.
#'
#' @param posterior.seed Integer. Posterior option.
#'
#' @return A list of class \code{"knn.riem.refit"} containing:
#'   \describe{
#'     \item{\code{fitted.values}}{Smoothed response(s). Vector if single
#'       response, matrix if multiple.}
#'     \item{\code{residuals}}{Residuals (y.new - fitted.values). In semi-supervised
#'       mode, residuals at unlabeled vertices are \code{NA}.}
#'     \item{\code{n.responses}}{Number of response columns processed.}
#'     \item{\code{per.column.gcv}}{Logical indicating whether per-column GCV
#'       was used.}
#'     \item{\code{y.vertices}}{Labeled vertex indices used (or \code{NULL}).}
#'     \item{\code{eta.used}}{Eta used for the refit (scalar) when \code{per.column.gcv = FALSE}.}
#'   }
#' 
#' \strong{Semi-supervised refit (y.vertices):}
#' Let \eqn{L} be the labeled vertex set and let \eqn{V} be the cached eigenvectors
#' (with eigenvalues \eqn{\lambda_1,\dots,\lambda_m}). The refit solves for spectral
#' coefficients \eqn{a}:
#' \deqn{(A + \eta \Lambda)a = b,}
#' where \eqn{A = V_L^T V_L}, \eqn{b = V_L^T y_L}, and \eqn{\Lambda = \mathrm{diag}(\lambda_1,\dots,\lambda_m)}.
#' The fitted field is \eqn{\hat y = V a} over all vertices.
#'
#' @export
refit.rdgraph.regression <- function(fitted.model,
                                     y.new,
                                     y.vertices = NULL,
                                     per.column.gcv = FALSE,
                                     n.candidates = 40L,
                                     n.cores = 1L,
                                     verbose = FALSE,
                                     with.posterior = FALSE,
                                     return.posterior.samples = FALSE,
                                     credible.level = 0.95,
                                     n.posterior.samples = 1000L,
                                     posterior.seed = 12345L,
                                     block.size = NULL) {

    ## ================================================================
    ## INPUT VALIDATION
    ## ================================================================

    if (!inherits(fitted.model, "knn.riem.fit")) {
        stop("fitted.model must be a 'knn.riem.fit' object from fit.rdgraph.regression()")
    }
    if (is.null(fitted.model$spectral)) {
        stop("Fitted model does not contain spectral component. ",
             "Cannot refit without cached eigendecomposition.")
    }

    ## Validate and extract spectral components
    V <- fitted.model$spectral$eigenvectors
    eigenvalues <- fitted.model$spectral$eigenvalues
    f.lambda <- fitted.model$spectral$filtered.eigenvalues
    eta.used <- fitted.model$spectral$eta.optimal

    if (is.null(V) || !is.matrix(V)) stop("fitted.model$spectral$eigenvectors must be a matrix.")
    n <- nrow(V)
    m <- ncol(V)

    if (is.null(eigenvalues) || !is.numeric(eigenvalues)) {
        stop("fitted.model$spectral$eigenvalues must be present and numeric.")
    }
    eigenvalues <- as.double(eigenvalues)
    if (length(eigenvalues) < m) {
        stop("Length of spectral$eigenvalues is smaller than ncol(eigenvectors).")
    }
    eigenvalues <- eigenvalues[seq_len(m)]

    if (is.null(eta.used) || !is.numeric(eta.used) || length(eta.used) != 1L || !is.finite(eta.used)) {
        stop("fitted.model$spectral$eta.optimal must be a finite numeric scalar.")
    }
    eta.used <- as.double(eta.used)

    ## Validate y.vertices
    if (!is.null(y.vertices)) {
        if (per.column.gcv) {
            stop("per.column.gcv = TRUE is not supported when y.vertices is provided.")
        }
        if (with.posterior) {
            stop("with.posterior = TRUE is not supported when y.vertices is provided.")
        }
        if (!is.numeric(y.vertices) && !is.integer(y.vertices)) {
            stop("y.vertices must be an integer (or numeric coercible to integer) vector, or NULL.")
        }
        y.vertices <- as.integer(y.vertices)
        if (length(y.vertices) < 1L) stop("y.vertices must have positive length when provided.")
        if (anyNA(y.vertices)) stop("y.vertices cannot contain NA.")
        if (any(y.vertices < 1L | y.vertices > n)) stop("y.vertices must be in 1..n.")
        if (anyDuplicated(y.vertices)) stop("y.vertices cannot contain duplicates.")
        y.vertices <- sort(y.vertices)
    }

    ## Validate n.cores
    n.cores <- as.integer(n.cores)
    if (n.cores < 1L) {
        warning("n.cores must be >= 1; setting to 1", call. = FALSE)
        n.cores <- 1L
    }

    ## Validate block.size
    if (is.null(block.size)) {
        block.size <- NA_integer_
    } else {
        if (!is.numeric(block.size) || length(block.size) != 1L ||
            is.na(block.size) || !is.finite(block.size) ||
            block.size < 1 || block.size != floor(block.size)) {
            stop("block.size must be NULL or a positive integer.")
        }
        block.size <- as.integer(block.size)
    }

    .as_dense_matrix <- function(x) {
        if (inherits(x, "Matrix")) as.matrix(x) else x
    }

    .make_block_index <- function(p, bsize) {
        if (is.na(bsize) || p <= 1L) {
            return(list(seq_len(p)))
        }
        split(seq_len(p), ceiling(seq_len(p) / bsize))
    }

    ## Prepare y.new
    is.matrix.input <- is.matrix(y.new) || inherits(y.new, "Matrix")
    is.sparse.input <- inherits(y.new, "sparseMatrix")

    if (is.matrix.input) {
        if (nrow(y.new) != n) {
            stop(sprintf("Number of rows in y.new (%d) must match fitted model (%d)",
                         nrow(y.new), n))
        }
        n.responses <- ncol(y.new)
        col.names <- colnames(y.new)
    } else {
        if (length(y.new) != n) {
            stop(sprintf("Length of y.new (%d) must match fitted model (%d)",
                         length(y.new), n))
        }
        y.new <- matrix(as.double(y.new), ncol = 1)
        n.responses <- 1L
        col.names <- NULL
    }

    block.index <- .make_block_index(n.responses, block.size)
    blockwise.enabled <- (n.responses > 1L && length(block.index) > 1L)

    ## Finite-value checks
    if (is.null(y.vertices)) {
        ## Supervised: all entries must be finite
        if (is.sparse.input) {
            if (any(!is.finite(y.new@x))) {
                bad.cols <- which(vapply(
                    seq_len(n.responses),
                    function(j) {
                        xj <- y.new[, j, drop = FALSE]@x
                        any(!is.finite(xj))
                    },
                    logical(1L)
                ))
                stop(sprintf("y.new contains non-finite sparse entries in column(s): %s",
                             paste(bad.cols, collapse = ", ")))
            }
        } else {
            y.new.dense <- .as_dense_matrix(y.new)
            if (any(!is.finite(y.new.dense))) {
                bad.cols <- which(vapply(
                    seq_len(n.responses),
                    function(j) any(!is.finite(y.new.dense[, j])),
                    logical(1L)
                ))
                stop(sprintf("y.new contains non-finite values in column(s): %s",
                             paste(bad.cols, collapse = ", ")))
            }
        }
    } else {
        ## Semi-supervised: require finiteness only on labeled vertices
        y.lab <- y.new[y.vertices, , drop = FALSE]
        if (inherits(y.lab, "sparseMatrix")) {
            if (any(!is.finite(y.lab@x))) {
                bad.cols <- which(vapply(
                    seq_len(n.responses),
                    function(j) {
                        xj <- y.lab[, j, drop = FALSE]@x
                        any(!is.finite(xj))
                    },
                    logical(1L)
                ))
                stop(sprintf("y.new contains non-finite sparse entries on labeled vertices in column(s): %s",
                             paste(bad.cols, collapse = ", ")))
            }
        } else {
            y.lab.dense <- .as_dense_matrix(y.lab)
            if (any(!is.finite(y.lab.dense))) {
                bad.cols <- which(vapply(
                    seq_len(n.responses),
                    function(j) any(!is.finite(y.lab.dense[, j])),
                    logical(1L)
                ))
                stop(sprintf("y.new contains non-finite values on labeled vertices in column(s): %s",
                             paste(bad.cols, collapse = ", ")))
            }
        }
    }

    ## ================================================================
    ## BRANCH: Per-column GCV vs. fixed filter
    ## ================================================================

    if (per.column.gcv) {

        if (blockwise.enabled && isTRUE(verbose)) {
            message("block.size is ignored when per.column.gcv = TRUE.")
        }

        start.time <- Sys.time()

        ## Extract raw eigenvalues and filter type from fitted model
        eigenvalues <- fitted.model$spectral$eigenvalues
        filter.type <- fitted.model$spectral$filter.type

        ## Default to heat_kernel if filter.type not stored (older models)
        if (is.null(filter.type)) {
            filter.type <- "heat_kernel"
            if (verbose) {
                message("Filter type not stored in model; using 'heat_kernel'")
            }
        }

        if (verbose) {
            message(sprintf("Using '%s' filter for per-column GCV", filter.type))
        }

        ## Project all responses onto eigenbasis once: V^T Y (m x p)
        Vt.Y <- crossprod(V, y.new)

        ## Generate eta grid
        eta.grid <- generate.eta.grid(eigenvalues, filter.type, n.candidates)

        ## Pre-compute filter weights for all eta values: m x n.candidates
        filter.weights.matrix <- compute.filter.weights.matrix(
            eigenvalues, eta.grid, filter.type
        )

        ## Pre-compute trace(S) for all eta (used in all columns)
        trace.S.all <- colSums(filter.weights.matrix)

        ## ============================================================
        ## PARALLEL vs SEQUENTIAL PROCESSING
        ## ============================================================

        use.parallel <- (n.cores > 1L) && (n.responses > 1L)

        parallel.completed <- FALSE

        if (use.parallel) {
            ## Check for required packages
            if (!requireNamespace("foreach", quietly = TRUE)) {
                warning("Package 'foreach' not available; falling back to sequential processing")
                use.parallel <- FALSE
            }
            if (use.parallel && !requireNamespace("doParallel", quietly = TRUE)) {
                warning("Package 'doParallel' not available; falling back to sequential processing")
                use.parallel <- FALSE
            }
        }

        if (use.parallel) {

            ## --------------------------------------------------------
            ## PARALLEL EXECUTION
            ## --------------------------------------------------------

            ## Determine actual number of cores to use.
            ## On some systems detectCores() can return NA; guard and fall back.
            max.cores <- suppressWarnings(parallel::detectCores(logical = FALSE))
            if (!is.finite(max.cores) || max.cores < 1L) {
                max.cores <- suppressWarnings(parallel::detectCores(logical = TRUE))
            }
            if (!is.finite(max.cores) || max.cores < 1L) {
                max.cores <- n.cores
                if (verbose) {
                    message(sprintf("parallel::detectCores() unavailable; using requested n.cores=%d", n.cores))
                }
            }
            max.cores <- as.integer(max.cores)
            n.cores.use <- as.integer(min(n.cores, max.cores, n.responses))
            if (!is.finite(n.cores.use) || n.cores.use < 1L) {
                n.cores.use <- 1L
            }

            if (verbose) {
                message(sprintf("Processing %d columns in parallel using %d cores...",
                                n.responses, n.cores.use))
            }

            ## Register parallel backend
            cl <- tryCatch(
                parallel::makeCluster(n.cores.use),
                error = function(e) {
                    warning(
                        sprintf("Could not start parallel PSOCK cluster (%s); falling back to sequential processing.",
                                conditionMessage(e)),
                        call. = FALSE
                    )
                    NULL
                }
            )

            if (!is.null(cl)) {
                doParallel::registerDoParallel(cl)
                on.exit(parallel::stopCluster(cl), add = TRUE)
                `%dopar%` <- foreach::`%dopar%`

                ## Define iterator variable for foreach
                j <- NULL  # Avoid R CMD check NOTE about undefined global

                ## Parallel loop over columns
                results.list <- foreach::foreach(
                                             j = seq_len(n.responses),
                                             .combine = "rbind",
                                             .packages = character(0),
                                             .export = c("select.eta.gcv.single.fast")
                                         ) %dopar% {

                                             y.obs.j <- y.new[, j]
                                             y.spectral.j <- Vt.Y[, j]

                                             ## Fast GCV selection (no intermediate y.hat storage)
                                             gcv.result <- select.eta.gcv.single.fast(
                                                 y.obs.j, y.spectral.j, V,
                                                 filter.weights.matrix, eta.grid, trace.S.all
                                             )

                                             ## Return as single-row data frame for rbind combining
                                             data.frame(
                                                 col.idx = j,
                                                 eta.optimal = gcv.result$eta.optimal,
                                                 gcv.min = gcv.result$gcv.min,
                                                 effective.df = gcv.result$effective.df,
                                                 best.idx = gcv.result$best.idx
                                             )
                                         }

                ## Reconstruct Y.hat using optimal indices
                Y.hat <- matrix(0, nrow = n, ncol = n.responses)
                for (i in seq_len(nrow(results.list))) {
                    j <- results.list$col.idx[i]
                    best.idx <- results.list$best.idx[i]
                    y.spectral.j <- Vt.Y[, j]
                    filtered.spectral <- y.spectral.j * filter.weights.matrix[, best.idx]
                    Y.hat[, j] <- V %*% filtered.spectral
                }

                eta.optimal <- results.list$eta.optimal
                gcv.scores <- results.list$gcv.min
                effective.df <- results.list$effective.df
                best.idx.vec <- results.list$best.idx
                n.cores.used <- n.cores.use
                parallel.completed <- TRUE
            }
        }

        if (!parallel.completed) {

            ## --------------------------------------------------------
            ## SEQUENTIAL EXECUTION
            ## --------------------------------------------------------

            n.cores.used <- 1L

            ## Initialize output storage
            Y.hat <- matrix(0, nrow = n, ncol = n.responses)
            eta.optimal <- numeric(n.responses)
            gcv.scores <- numeric(n.responses)
            effective.df <- numeric(n.responses)
            best.idx.vec <- integer(n.responses)

            if (verbose && n.responses > 1L) {
                message(sprintf("Selecting eta via GCV for %d response(s)...",
                                n.responses))
                pb <- utils::txtProgressBar(min = 0, max = n.responses, style = 3)
            }

            ## Process each column
            for (j in seq_len(n.responses)) {

                y.spectral <- Vt.Y[, j]

                ## Compute GCV scores for all eta candidates (vectorized)
                gcv.result <- select.eta.gcv.single(
                    y.new[, j], y.spectral, V, filter.weights.matrix, eta.grid
                )

                Y.hat[, j] <- gcv.result$y.hat
                eta.optimal[j] <- gcv.result$eta.optimal
                gcv.scores[j] <- gcv.result$gcv.min
                effective.df[j] <- gcv.result$effective.df
                best.idx.vec[j] <- gcv.result$best.idx

                if (verbose && n.responses > 1L) utils::setTxtProgressBar(pb, j)
            }

            if (verbose && n.responses > 1L) {
                close(pb)
            }
        }

        elapsed.time <- as.numeric(difftime(Sys.time(), start.time, units = "secs"))

        if (verbose) {
            message(sprintf("Done. Processed %d columns in %.1f seconds (%.1f cols/sec)",
                            n.responses, elapsed.time, n.responses / elapsed.time))
            message(sprintf("Eta range: [%.2e, %.2e]",
                            min(eta.optimal), max(eta.optimal)))
        }

        ## Preserve column names if present
        if (!is.null(col.names)) {
            colnames(Y.hat) <- col.names
            names(eta.optimal) <- col.names
            names(gcv.scores) <- col.names
            names(effective.df) <- col.names
        }

        ## Compute residuals
        residuals <- y.new - Y.hat

        ## ============================================================
        ## POSTERIOR INFERENCE (PER-COLUMN GCV CASE)
        ## ============================================================

        posterior <- NULL

        if (with.posterior) {

            if (verbose) {
                message("Computing posterior credible intervals...")
            }

            ## Initialize storage for posterior summaries
            posterior.lower <- matrix(0, nrow = n, ncol = n.responses)
            posterior.upper <- matrix(0, nrow = n, ncol = n.responses)
            posterior.sd <- matrix(0, nrow = n, ncol = n.responses)
            posterior.sigma <- numeric(n.responses)

            ## Optionally store samples
            if (return.posterior.samples) {
                posterior.samples <- vector("list", n.responses)
            }

            if (verbose && n.responses > 1L) {
                pb.post <- utils::txtProgressBar(min = 0, max = n.responses, style = 3)
            }

            for (j in seq_len(n.responses)) {

                ## Get the filtered eigenvalues for this column's optimal eta
                filtered.eigenvalues.j <- filter.weights.matrix[, best.idx.vec[j]]

                ## Call C++ posterior computation
                post.j <- .Call(
                    "S_compute_posterior_summary",
                    V,
                    eigenvalues,
                    filtered.eigenvalues.j,
                    y.new[, j],
                    Y.hat[, j],
                    eta.optimal[j],
                    credible.level,
                    n.posterior.samples,
                    posterior.seed + j - 1L,  # Different seed per column
                    return.posterior.samples,
                    PACKAGE = "gflowx"
                )

                posterior.lower[, j] <- post.j$lower
                posterior.upper[, j] <- post.j$upper
                posterior.sd[, j] <- post.j$sd
                posterior.sigma[j] <- post.j$sigma

                if (return.posterior.samples) {
                    posterior.samples[[j]] <- post.j$samples
                }

                if (verbose && n.responses > 1L) {
                    utils::setTxtProgressBar(pb.post, j)
                }
            }

            if (verbose && n.responses > 1L) {
                close(pb.post)
            }

            ## Preserve column names
            if (!is.null(col.names)) {
                colnames(posterior.lower) <- col.names
                colnames(posterior.upper) <- col.names
                colnames(posterior.sd) <- col.names
                names(posterior.sigma) <- col.names
                if (return.posterior.samples) {
                    names(posterior.samples) <- col.names
                }
            }

            ## Assemble posterior component
            posterior <- list(
                lower = if (n.responses == 1L) as.vector(posterior.lower) else posterior.lower,
                upper = if (n.responses == 1L) as.vector(posterior.upper) else posterior.upper,
                sd = if (n.responses == 1L) as.vector(posterior.sd) else posterior.sd,
                sigma = if (n.responses == 1L) posterior.sigma[1] else posterior.sigma,
                credible.level = credible.level,
                n.samples = n.posterior.samples
            )

            if (return.posterior.samples) {
                posterior$samples <- if (n.responses == 1L) {
                                         posterior.samples[[1]]
                                     } else {
                                         posterior.samples
                                     }
            }

            if (verbose) {
                message("Posterior computation complete.")
            }
        }

        ## Build result
        result <- list(
            fitted.values = if (n.responses == 1L) as.vector(Y.hat) else Y.hat,
            residuals = if (n.responses == 1L) as.vector(residuals) else residuals,
            n.responses = n.responses,
            per.column.gcv = TRUE,
            filter.type = filter.type,
            eta.optimal = if (n.responses == 1L) eta.optimal[1] else eta.optimal,
            gcv.scores = if (n.responses == 1L) gcv.scores[1] else gcv.scores,
            effective.df = if (n.responses == 1L) effective.df[1] else effective.df,
            eta.grid = eta.grid,
            y.vertices = NULL,
            eta.used = NULL,
            block.size.used = NA_integer_,
            n.cores.used = n.cores.used,
            elapsed.time = elapsed.time
        )

        if (with.posterior) {
            result$posterior <- posterior
        }

        class(result) <- c("knn.riem.refit", "list")
        return(result)
    }

    ## ================================================================
    ## FIXED-ETA REFIT
    ## ================================================================

    if (is.null(y.vertices)) {

        ## ------------------------------------------------------------
        ## ORIGINAL BEHAVIOR: Apply cached spectral filter weights
        ## ------------------------------------------------------------
        if (is.null(f.lambda)) {
            stop("Fitted model does not contain filtered.eigenvalues. ",
                 "Cannot refit without spectral filter information.")
        }
        f.lambda <- as.double(f.lambda)
        if (length(f.lambda) != m) {
            stop("Length mismatch: length(filtered.eigenvalues) must equal ncol(eigenvectors).")
        }

        if (blockwise.enabled && isTRUE(verbose)) {
            message(sprintf(
                "Fixed-eta refit in %d column blocks (block.size=%d).",
                length(block.index), block.size
            ))
        }

        Y.hat <- matrix(0.0, nrow = n, ncol = n.responses)
        residuals <- matrix(0.0, nrow = n, ncol = n.responses)
        if (!is.null(col.names) && n.responses > 1L) {
            colnames(Y.hat) <- col.names
            colnames(residuals) <- col.names
        }

        for (bi in seq_along(block.index)) {
            cols <- block.index[[bi]]

            y.block <- y.new[, cols, drop = FALSE]

            ## Apply spectral filter blockwise:
            ## Y_hat_block = V diag(f.lambda) V^T Y_block
            Vt.Y.block <- crossprod(V, y.block)
            filtered.block <- f.lambda * Vt.Y.block
            Y.hat.block <- .as_dense_matrix(V %*% filtered.block)

            Y.hat[, cols] <- Y.hat.block
            residuals[, cols] <- .as_dense_matrix(y.block) - Y.hat.block

            if (blockwise.enabled && isTRUE(verbose)) {
                message(sprintf("  processed block %d/%d", bi, length(block.index)))
            }
        }

        result <- list(
            fitted.values = if (n.responses == 1L) as.vector(Y.hat) else Y.hat,
            residuals = if (n.responses == 1L) as.vector(residuals) else residuals,
            n.responses = n.responses,
            per.column.gcv = FALSE,
            y.vertices = NULL,
            eta.used = eta.used,
            block.size.used = if (blockwise.enabled) block.size else NA_integer_
        )

        class(result) <- c("knn.riem.refit", "list")
        return(result)

    } else {

        ## ------------------------------------------------------------
        ## SEMI-SUPERVISED: Masked spectral Tikhonov refit
        ## ------------------------------------------------------------
        ## Solve (A + eta*Lambda)a = b where A = V_L^T V_L, b = V_L^T y_L
        if (blockwise.enabled && isTRUE(verbose)) {
            message(sprintf(
                "Semi-supervised fixed-eta refit in %d column blocks (block.size=%d).",
                length(block.index), block.size
            ))
        }

        V.lab <- V[y.vertices, , drop = FALSE]          ## |L| x m

        A.mat <- crossprod(V.lab)                      ## m x m

        ## M = A + eta * diag(eigenvalues)
        M.mat <- A.mat
        diag(M.mat) <- diag(M.mat) + eta.used * eigenvalues

        ## Solve efficiently via Cholesky; fallback to solve() if needed
        chol.M <- tryCatch(chol(M.mat), error = function(e) NULL)

        if (is.null(chol.M)) {
            ## Rare numerical fallback: add tiny ridge then retry
            ridge <- sqrt(.Machine$double.eps)
            if (verbose) {
                message(sprintf("Cholesky failed; adding ridge %.3e to diagonal and retrying.", ridge))
            }
            M2 <- M.mat
            diag(M2) <- diag(M2) + ridge
            chol.M <- tryCatch(chol(M2), error = function(e) NULL)
        }

        Y.hat <- matrix(0.0, nrow = n, ncol = n.responses)
        residuals <- matrix(NA_real_, nrow = n, ncol = n.responses)

        for (bi in seq_along(block.index)) {
            cols <- block.index[[bi]]
            y.lab.block <- y.new[y.vertices, cols, drop = FALSE]
            b.block <- .as_dense_matrix(crossprod(V.lab, y.lab.block))

            if (!is.null(chol.M)) {
                ## a = M^{-1} b using triangular solves
                a.block <- backsolve(chol.M, forwardsolve(t(chol.M), b.block))
            } else {
                ## Final fallback
                if (verbose) message("Cholesky retry failed; using solve(M, b).")
                a.block <- solve(M.mat, b.block)
            }

            Y.hat.block <- .as_dense_matrix(V %*% a.block)
            Y.hat[, cols] <- Y.hat.block
            residuals[y.vertices, cols] <- .as_dense_matrix(y.lab.block) - Y.hat.block[y.vertices, , drop = FALSE]

            if (blockwise.enabled && isTRUE(verbose)) {
                message(sprintf("  processed block %d/%d", bi, length(block.index)))
            }
        }

        if (!is.null(col.names) && n.responses > 1L) {
            colnames(Y.hat) <- col.names
            colnames(residuals) <- col.names
        }

        result <- list(
            fitted.values = if (n.responses == 1L) as.vector(Y.hat) else Y.hat,
            residuals = if (n.responses == 1L) as.vector(residuals) else residuals,
            n.responses = n.responses,
            per.column.gcv = FALSE,
            y.vertices = as.integer(y.vertices),
            eta.used = eta.used,
            block.size.used = if (blockwise.enabled) block.size else NA_integer_
        )

        class(result) <- c("knn.riem.refit", "list")
        return(result)
    }
}

#' Blockwise Refit Wrapper for Riemannian Graph Regression
#'
#' Convenience wrapper around \code{\link{refit.rdgraph.regression}} that enables
#' blockwise processing by default.
#'
#' @inheritParams refit.rdgraph.regression
#' @param block.size Positive integer block size. Default is 250.
#' @param ... Additional arguments passed to \code{\link{refit.rdgraph.regression}}.
#'
#' @return A \code{"knn.riem.refit"} object returned by
#'   \code{\link{refit.rdgraph.regression}}.
refit.blocked <- function(fitted.model,
                          y.new,
                          block.size = 250L,
                          ...) {
    refit.rdgraph.regression(
        fitted.model = fitted.model,
        y.new = y.new,
        block.size = block.size,
        ...
    )
}

#' Compute filter weights for all eigenvalues and all eta candidates
#'
#' @param eigenvalues Vector of eigenvalues
#' @param eta.grid Vector of eta values
#' @param filter.type Filter type string
#' @return Matrix of filter weights (m x length(eta.grid))
#' @keywords internal
compute.filter.weights.matrix <- function(eigenvalues, eta.grid, filter.type) {

    m <- length(eigenvalues)
    n.eta <- length(eta.grid)

    ## Pre-allocate matrix (eigenvalues in rows, eta in columns)
    weights <- matrix(0, nrow = m, ncol = n.eta)

    for (k in seq_len(n.eta)) {
        eta <- eta.grid[k]

        weights[, k] <- switch(filter.type,
            "heat_kernel" = exp(-eta * eigenvalues),
            "tikhonov" = 1.0 / (1.0 + eta * eigenvalues),
            "cubic_spline" = 1.0 / (1.0 + eta * eigenvalues^2),
            "gaussian" = exp(-eta * eigenvalues^2),
            "exponential" = exp(-eta * sqrt(pmax(eigenvalues, 0))),
            "butterworth" = {
                ## Order-2 Butterworth: 1/(1 + (lambda/eta)^4)
                x <- eigenvalues / eta
                1.0 / (1.0 + x^4)
            },
            ## Default to heat kernel
            exp(-eta * eigenvalues)
        )
    }

    return(weights)
}

#' Single-column GCV selection with full fitted values
#'
#' @param y.obs Observed response vector
#' @param y.spectral Spectral coefficients V^T y
#' @param V Eigenvector matrix
#' @param filter.weights.matrix Pre-computed filter weights (m x n.eta)
#' @param eta.grid Vector of eta candidates
#' @return List with y.hat, eta.optimal, gcv.min, effective.df
#' @keywords internal
select.eta.gcv.single <- function(y.obs, y.spectral, V,
                                   filter.weights.matrix, eta.grid) {

    n <- length(y.obs)
    m <- length(y.spectral)
    n.eta <- length(eta.grid)

    ## Compute filtered spectral coefficients for all eta: m x n.eta
    filtered.spectral <- y.spectral * filter.weights.matrix

    ## Compute predictions for all eta: n x n.eta
    Y.hat.all <- V %*% filtered.spectral

    ## Compute RSS and effective df for all eta
    residuals.all <- y.obs - Y.hat.all
    rss <- colSums(residuals.all^2)
    trace.S <- colSums(filter.weights.matrix)

    ## GCV score = RSS / (n - trace(S))^2
    denom <- pmax(n - trace.S, 1e-10)
    gcv.scores <- rss / (denom^2)

    ## Find optimal
    best.idx <- which.min(gcv.scores)

    return(list(
        y.hat = Y.hat.all[, best.idx],
        eta.optimal = eta.grid[best.idx],
        gcv.min = gcv.scores[best.idx],
        effective.df = trace.S[best.idx],
        best.idx = best.idx
    ))
}

#' Fast single-column GCV selection (returns only scalars)
#'
#' Optimized version that doesn't return full y.hat vector.
#' Used in parallel processing where we reconstruct y.hat afterward.
#'
#' @param y.obs Observed response vector
#' @param y.spectral Spectral coefficients V^T y
#' @param V Eigenvector matrix
#' @param filter.weights.matrix Pre-computed filter weights (m x n.eta)
#' @param eta.grid Vector of eta candidates
#' @param trace.S.all Pre-computed trace(S) for all eta candidates
#' @return List with eta.optimal, gcv.min, effective.df, best.idx
#' @keywords internal
select.eta.gcv.single.fast <- function(y.obs, y.spectral, V,
                                        filter.weights.matrix, eta.grid,
                                        trace.S.all) {

    n <- length(y.obs)
    n.eta <- length(eta.grid)

    ## Compute filtered spectral coefficients for all eta: m x n.eta
    filtered.spectral <- y.spectral * filter.weights.matrix

    ## Compute predictions for all eta: n x n.eta
    Y.hat.all <- V %*% filtered.spectral

    ## Compute RSS for all eta (vectorized)
    residuals.all <- y.obs - Y.hat.all
    rss <- colSums(residuals.all^2)

    ## GCV score = RSS / (n - trace(S))^2
    denom <- pmax(n - trace.S.all, 1e-10)
    gcv.scores <- rss / (denom^2)

    ## Find optimal
    best.idx <- which.min(gcv.scores)

    return(list(
        eta.optimal = eta.grid[best.idx],
        gcv.min = gcv.scores[best.idx],
        effective.df = trace.S.all[best.idx],
        best.idx = best.idx
    ))
}

## ============================================================
## PRINT AND SUMMARY METHODS
## ============================================================

#' @export
print.knn.riem.refit <- function(x, ...) {
    cat("\nRefitted Riemannian Graph Regression\n")
    cat("====================================\n\n")

    if (x$n.responses == 1) {
        cat(sprintf("Single response variable (n = %d)\n", length(x$fitted.values)))
        cat(sprintf("RMSE: %.4f\n", sqrt(mean(x$residuals^2, na.rm=TRUE))))
        cat(sprintf("MAE:  %.4f\n", mean(abs(x$residuals), na.rm=TRUE)))

        if (isTRUE(x$per.column.gcv)) {
            cat(sprintf("\nPer-column GCV (%s filter):\n", x$filter.type))
            cat(sprintf("  Optimal eta:    %.4e\n", x$eta.optimal))
            cat(sprintf("  Effective df:   %.1f\n", x$effective.df))
            cat(sprintf("  GCV score:      %.4e\n", x$gcv.scores))
        }

        if (!is.null(x$posterior)) {
            cat(sprintf("\nPosterior inference (%.0f%% credible intervals):\n",
                        x$posterior$credible.level * 100))
            cat(sprintf("  Residual sigma: %.4f\n", x$posterior$sigma))
            ci.width <- x$posterior$upper - x$posterior$lower
            cat(sprintf("  Mean CI width:  %.4f\n", mean(ci.width, na.rm=TRUE)))
        }
    } else {
        cat(sprintf("Multiple responses (n = %d, p = %d)\n",
                    nrow(x$fitted.values), x$n.responses))
        rmse.by.col <- apply(x$residuals, 2, function(r) sqrt(mean(r^2, na.rm=TRUE)))
        cat(sprintf("Mean RMSE across responses: %.4f\n", mean(rmse.by.col, na.rm=TRUE)))
        cat(sprintf("RMSE range: [%.4f, %.4f]\n", min(rmse.by.col), max(rmse.by.col)))

        if (isTRUE(x$per.column.gcv)) {
            cat(sprintf("\nPer-column GCV selection (%s filter):\n", x$filter.type))
            cat(sprintf("  Eta range:          [%.2e, %.2e]\n",
                        min(x$eta.optimal), max(x$eta.optimal)))
            cat(sprintf("  Effective df range: [%.1f, %.1f]\n",
                        min(x$effective.df), max(x$effective.df)))

            if (!is.null(x$n.cores.used) && !is.null(x$elapsed.time)) {
                cat(sprintf("\nComputation: %.1f sec using %d core%s (%.1f cols/sec)\n",
                            x$elapsed.time, x$n.cores.used,
                            if (x$n.cores.used > 1) "s" else "",
                            x$n.responses / x$elapsed.time))
            }
        }

        if (!is.null(x$posterior)) {
            cat(sprintf("\nPosterior inference (%.0f%% credible intervals):\n",
                        x$posterior$credible.level * 100))
            cat(sprintf("  Sigma range:    [%.4f, %.4f]\n",
                        min(x$posterior$sigma), max(x$posterior$sigma)))
            ci.width <- x$posterior$upper - x$posterior$lower
            mean.ci.width <- colMeans(ci.width)
            cat(sprintf("  Mean CI width:  [%.4f, %.4f]\n",
                        min(mean.ci.width), max(mean.ci.width)))
        }
    }

    invisible(x)
}

#' @export
summary.knn.riem.refit <- function(object, ...) {
    n.resp <- object$n.responses
    n.obs <- if (n.resp == 1) length(object$fitted.values) else nrow(object$fitted.values)

    if (n.resp == 1) {
        ## Single response summary
        y.hat <- object$fitted.values
        resid <- object$residuals
        y.obs <- y.hat + resid

        rss <- sum(resid^2)
        tss <- sum((y.obs - mean(y.obs, na.rm=TRUE))^2)
        r.squared <- 1 - rss / tss
        rmse <- sqrt(mean(resid^2, na.rm=TRUE))
        mae <- mean(abs(resid), na.rm=TRUE)

        result <- list(
            n.obs = n.obs,
            n.responses = 1L,
            r.squared = r.squared,
            rmse = rmse,
            mae = mae,
            per.column.gcv = isTRUE(object$per.column.gcv),
            eta.optimal = object$eta.optimal,
            effective.df = object$effective.df,
            filter.type = object$filter.type,
            n.cores.used = object$n.cores.used,
            elapsed.time = object$elapsed.time
        )

        if (!is.null(object$posterior)) {
            result$posterior.summary <- list(
                credible.level = object$posterior$credible.level,
                sigma = object$posterior$sigma,
                mean.ci.width = mean(object$posterior$upper - object$posterior$lower, na.rm=TRUE)
            )
        }
    } else {
        ## Multiple response summary
        Y.hat <- object$fitted.values
        residuals <- object$residuals
        Y.obs <- Y.hat + residuals

        ## Per-column metrics
        rmse.by.col <- apply(residuals, 2, function(r) sqrt(mean(r^2, na.rm=TRUE)))
        mae.by.col <- apply(residuals, 2, function(r) mean(abs(r), na.rm=TRUE))

        r.squared.by.col <- sapply(seq_len(n.resp), function(j) {
            rss <- sum(residuals[, j]^2)
            tss <- sum((Y.obs[, j] - mean(Y.obs[, j]))^2, na.rm=TRUE)
            if (tss > 1e-15) 1 - rss / tss else NA_real_
        })

        result <- list(
            n.obs = n.obs,
            n.responses = n.resp,
            rmse.by.column = rmse.by.col,
            mae.by.column = mae.by.col,
            r.squared.by.column = r.squared.by.col,
            per.column.gcv = isTRUE(object$per.column.gcv),
            eta.optimal = object$eta.optimal,
            effective.df = object$effective.df,
            filter.type = object$filter.type,
            n.cores.used = object$n.cores.used,
            elapsed.time = object$elapsed.time
        )

        if (!is.null(object$posterior)) {
            ci.width <- object$posterior$upper - object$posterior$lower
            result$posterior.summary <- list(
                credible.level = object$posterior$credible.level,
                sigma.range = range(object$posterior$sigma),
                mean.ci.width.by.col = colMeans(ci.width)
            )
        }
    }

    class(result) <- "summary.knn.riem.refit"
    return(result)
}

#' @export
print.summary.knn.riem.refit <- function(x, digits = 4, ...) {
    cat("\nSummary: Refitted Riemannian Graph Regression\n")
    cat("==============================================\n\n")

    cat(sprintf("Observations: %d\n", x$n.obs))
    cat(sprintf("Responses: %d\n", x$n.responses))

    if (x$n.responses == 1) {
        cat(sprintf("\nFit quality:\n"))
        cat(sprintf("  R-squared: %.*f\n", digits, x$r.squared))
        cat(sprintf("  RMSE:      %.*f\n", digits, x$rmse))
        cat(sprintf("  MAE:       %.*f\n", digits, x$mae))

        if (x$per.column.gcv) {
            cat(sprintf("\nGCV-selected parameters (%s filter):\n", x$filter.type))
            cat(sprintf("  Optimal eta:    %.*e\n", digits, x$eta.optimal))
            cat(sprintf("  Effective df:   %.1f\n", x$effective.df))
        }

        if (!is.null(x$posterior.summary)) {
            cat(sprintf("\nPosterior inference (%.0f%% CI):\n",
                        x$posterior.summary$credible.level * 100))
            cat(sprintf("  Residual sigma: %.*f\n", digits, x$posterior.summary$sigma))
            cat(sprintf("  Mean CI width:  %.*f\n", digits, x$posterior.summary$mean.ci.width))
        }
    } else {
        cat(sprintf("\nFit quality summary across %d responses:\n", x$n.responses))

        print.quantile.summary <- function(vals, name) {
            q <- quantile(vals, probs = c(0, 0.25, 0.5, 0.75, 1), na.rm = TRUE)
            cat(sprintf("  %s: min=%.4f, Q1=%.4f, med=%.4f, Q3=%.4f, max=%.4f\n",
                        name, q[1], q[2], q[3], q[4], q[5]))
        }

        print.quantile.summary(x$r.squared.by.column, "R-sq")
        print.quantile.summary(x$rmse.by.column, "RMSE")

        if (x$per.column.gcv) {
            cat(sprintf("\nGCV-selected parameters (%s filter):\n", x$filter.type))
            print.quantile.summary(x$eta.optimal, "Eta ")
            print.quantile.summary(x$effective.df, "Eff.df")

            if (!is.null(x$n.cores.used) && !is.null(x$elapsed.time)) {
                cat(sprintf("\nComputation: %.1f sec using %d core%s\n",
                            x$elapsed.time, x$n.cores.used,
                            if (x$n.cores.used > 1) "s" else ""))
            }
        }

        if (!is.null(x$posterior.summary)) {
            cat(sprintf("\nPosterior inference (%.0f%% CI):\n",
                        x$posterior.summary$credible.level * 100))
            cat(sprintf("  Sigma range: [%.4f, %.4f]\n",
                        x$posterior.summary$sigma.range[1],
                        x$posterior.summary$sigma.range[2]))
            print.quantile.summary(x$posterior.summary$mean.ci.width.by.col, "CI width")
        }
    }

    invisible(x)
}
