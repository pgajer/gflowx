#' Compute ikNN Graphs and Riemannian Regression Under Power Transformation
#'
#' For each power transformation exponent p, transforms the abundance matrix
#' \eqn{X \mapsto X^p}, constructs ikNN graphs, fits Riemannian spectral
#' regression for the response y, and computes k-selection diagnostics. Results
#' are saved to disk and optionally computed in parallel.
#'
#' @param X Numeric matrix of relative abundances (samples x phylotypes).
#'   Should be non-negative with no duplicate rows.
#' @param y Numeric vector of response values (e.g., binary sPTB indicator).
#' @param p Numeric vector of power transformation exponents in (0, 1].
#' @param k.min Integer; minimum number of neighbors for ikNN graph.
#' @param k.max Integer; maximum number of neighbors for ikNN graph.
#' @param L1.normalize Logical; if TRUE, apply L1 normalization so that
#'   the sum of each row is 1.
#' @param out.dir Directory for saving results. If NULL, results are returned
#'   but not saved.
#' @param plot.dir Directory for saving diagnostic plots. If NULL, no plots
#'   are generated.
#' @param data.tag Character string identifying the dataset (used in filenames).
#' @param skip.existing Logical; if TRUE, skip p values with existing results.
#' @param n.cores Integer; number of cores for parallel computation.
#'   If 1, computation is sequential.
#' @param iknn.params Named list of additional parameters passed to
#'   \code{dgraphs::create.iknn.graphs()}.
#' @param regression.params Named list of additional parameters passed to
#'   \code{\link{fit.rdgraph.regression}}.
#' @param n.embeddings Integer; number of 3D graph embeddings to compute.
#' @param verbose Logical; print progress messages.
#'
#' @return A named list (keyed by p.tag) of results for each p value. Each
#'   element contains:
#'   \describe{
#'     \item{p}{Power transformation exponent.}
#'     \item{p.tag}{Tag string for p (e.g., "p065").}
#'     \item{status}{One of "success", "failed", "skipped", or "error".}
#'     \item{opt.k}{Optimal k selected by GCV (if status is "success").}
#'     \item{gcv}{GCV score at optimal k (if status is "success").}
#'     \item{n.extrema}{Number of local extrema at optimal k (if status is "success").}
#'     \item{n.connected.k}{Number of k values yielding connected graphs (if status is "success").}
#'     \item{file}{Path to saved results file (if out.dir provided and status is "success").}
#'     \item{reason}{Description of failure (if status is "failed", "skipped", or "error").}
#'   }
#'
#' @details
#' The function applies the power transformation \eqn{X \mapsto X^p} to the
#' abundance matrix without re-normalization (unless \code{L1.normalize = TRUE}).
#' This boosts the contribution of low-abundance phylotypes to inter-sample
#' distances. For each p value, the function:
#' \enumerate{
#'   \item Constructs ikNN graphs for k in \eqn{[k.min, k.max]}
#'   \item Identifies k values yielding connected graphs
#'   \item Fits Riemannian spectral regression for each valid k
#'   \item Selects optimal k via GCV
#'   \item Computes 3D graph embeddings at optimal k
#' }
#'
#' Full results are saved to disk (if out.dir is provided) to manage memory
#' for large datasets. The returned list contains lightweight summaries.
#'
#' @note Parallel execution requires the \pkg{doParallel} and \pkg{foreach}
#'   packages. If these are not installed and \code{n.cores > 1}, an error
#'   is raised.
#'
#' @seealso \code{dgraphs::create.iknn.graphs()} for graph construction,
#'   \code{\link{fit.rdgraph.regression}} for the regression method,
#'   \code{\link{k.diagnostics.plots}} for bandwidth selection diagnostics.
#'
#' @examples
#' \dontrun{
#' ## Load data (project-specific)
#' sptb.zmb.data <- load.sptb.data("Zambia")
#'
#' ## Run analysis for a grid of p values
#' p.grid <- seq(0.5, 1, by = 0.05)
#' results <- compute.Lp.graphs(
#'   X = sptb.zmb.data$X,
#'   y = sptb.zmb.data$y,
#'   p = p.grid,
#'   k.min = 6,
#'   k.max = 12,
#'   L1.normalize = FALSE,
#'   out.dir = "~/analysis/output",
#'   plot.dir = "~/analysis/figures",
#'   data.tag = sptb.zmb.data$data.tag,
#'   n.cores = 4
#' )
#' }
#'
#' @importFrom grDevices dev.off pdf
#' @importFrom utils modifyList
#' @export
compute.Lp.graphs <- function(X,
                              y,
                              p,
                              k.min = 6,
                              k.max = 12,
                              L1.normalize = FALSE,
                              out.dir = NULL,
                              plot.dir = NULL,
                              data.tag = "data",
                              skip.existing = TRUE,
                              n.cores = 1,
                              iknn.params = list(),
                              regression.params = list(),
                              n.embeddings = 5,
                              verbose = TRUE) {

    ## ==========================================================================
    ## Input validation
    ## ==========================================================================
    if (!is.matrix(X) && !is.data.frame(X)) {
        stop("X must be a matrix or data frame")
    }
    X <- as.matrix(X)

    if (!is.numeric(X) || any(X < 0)) {
        stop("X must be a non-negative numeric matrix")
    }

    if (any(duplicated(X))) {
        stop("X contains duplicate rows; please dereplicate before calling")
    }

    y <- as.double(y)
    if (length(y) != nrow(X)) {
        stop("length(y) must equal nrow(X)")
    }

    if (!is.numeric(p) || any(p <= 0) || any(p > 1)) {
        stop("p must be numeric with all values in (0, 1]")
    }

    if (!is.null(out.dir) && !dir.exists(out.dir)) {
        dir.create(out.dir, recursive = TRUE)
    }

    if (!is.null(plot.dir) && !dir.exists(plot.dir)) {
        dir.create(plot.dir, recursive = TRUE)
    }

    ## ==========================================================================
    ## Helper functions
    ## ==========================================================================
    make.p.tag <- function(p) {
        gsub("\\.", "", sprintf("p%g", p))
    }

    make.rds.filename <- function(p.tag, data.tag, L1.normalize) {
        if (L1.normalize) {
            sprintf("Lp_experiment_L1_normalized_%s_%s.rds", p.tag, data.tag)
        } else {
            sprintf("Lp_experiment_%s_%s.rds", p.tag, data.tag)
        }
    }

    ## Generate p.tags
    p.tags <- make.p.tag(p)
    names(p.tags) <- as.character(p)

    ## ==========================================================================
    ## Identify p values needing computation
    ## ==========================================================================
    if (skip.existing && !is.null(out.dir)) {
        needs.computation <- sapply(seq_along(p), function(i) {
            p.tag <- p.tags[i]
            p.rds.file <- file.path(out.dir,
                                    make.rds.filename(p.tag, data.tag, L1.normalize))
            !file.exists(p.rds.file)
        })
        p.to.compute <- p[needs.computation]
        p.tags.to.compute <- p.tags[needs.computation]

        if (verbose && any(!needs.computation)) {
            message("Skipping existing results for p = ",
                    paste(p[!needs.computation], collapse = ", "))
        }
    } else {
        p.to.compute <- p
        p.tags.to.compute <- p.tags
    }

    if (length(p.to.compute) == 0) {
        if (verbose) message("All results already exist. Nothing to compute.")
        return(invisible(list()))
    }

    if (verbose) {
        message("Computing for p = ", paste(round(p.to.compute, 3), collapse = ", "))
        message("Using ", n.cores, ifelse(n.cores == 1, " core", " cores"))
    }

    ## ==========================================================================
    ## Default parameters for ikNN and regression
    ## ==========================================================================
    iknn.defaults <- list(
        max.path.edge.ratio.deviation.thld = 0.1,
        path.edge.ratio.percentile = 0.5,
        compute.full = TRUE,
        pca.dim = 100,
        variance.explained = 0.99
    )
    iknn.params <- modifyList(iknn.defaults, iknn.params)

    regression.defaults <- list(
        pca.dim = 100,
        max.iterations = 10,
        n.eigenpairs = 50,
        compute.extremality = FALSE,
        t.scale.factor = 0.05,
        beta.coef.factor = 0.075,
        use.counting.measure = FALSE
    )
    regression.params <- modifyList(regression.defaults, regression.params)

    ## ==========================================================================
    ## Single p computation function
    ## ==========================================================================
    single.p.fn <- function(single.p,
                            p.tag,
                            X,
                            y,
                            k.min,
                            k.max,
                            L1.normalize,
                            out.dir,
                            plot.dir,
                            data.tag,
                            skip.existing,
                            iknn.params,
                            regression.params,
                            n.embeddings) {

        ## Check for existing results
        if (!is.null(out.dir)) {
            p.rds.file <- file.path(out.dir,
                                    make.rds.filename(p.tag, data.tag, L1.normalize))

            if (skip.existing && file.exists(p.rds.file)) {
                return(list(
                    p = single.p,
                    p.tag = p.tag,
                    status = "skipped",
                    reason = "results exist"
                ))
            }
        } else {
            p.rds.file <- NULL
        }

        ## Apply power transformation
        X.transformed <- X^single.p

        if (L1.normalize) {
            X.transformed <- X.transformed / rowSums(X.transformed)
        }

        ## Build ikNN graphs
        iknn.args <- c(
            list(
                X = X.transformed,
                kmin = k.min,
                kmax = k.max,
                n.cores = 1,
                verbose = FALSE
            ),
            iknn.params
        )
        graphs <- do.call(dgraphs::create.iknn.graphs, iknn.args)

        graphs.stats <- summary(graphs)

        ## Identify k values with connected graphs
        idx <- graphs.stats$n_ccomp == 1
        good.indices <- graphs.stats$idx[idx]
        good.k.values <- graphs.stats$k[idx]
        good.k <- list(indices = good.indices, k.values = good.k.values)

        if (length(good.k.values) == 0) {
            return(list(
                p = single.p,
                p.tag = p.tag,
                status = "failed",
                reason = "no connected k values"
            ))
        }

        ## Fit regression for each good k
        rcx.res <- list()
        for (i in seq_along(good.indices)) {
            reg.args <- c(
                list(
                    X = X.transformed,
                    y = y,
                    k = good.k.values[i],
                    verbose.level = 0L
                ),
                regression.params
            )
            rcx.res[[i]] <- do.call(fit.rdgraph.regression, reg.args)
        }
        names(rcx.res) <- paste0("k", good.k.values)

        ## k-selection diagnostics
        k.diag <- k.diagnostics.plots(
            rcx.res,
            good.k.values,
            y,
            main.title = sprintf("Bandwidth Selection (p = %g): %s", single.p, data.tag),
            mark.optimal = TRUE
        )

        opt.k.info <- get.rcx.optimal.k(k.diag)
        opt.k <- opt.k.info$k

        ## Save diagnostic plot if plot.dir provided
        if (!is.null(plot.dir)) {
            k.diag.pdf <- file.path(plot.dir,
                                    sprintf("k_diagnostic_plot_%s_%s.pdf", p.tag, data.tag))
            grDevices::pdf(k.diag.pdf, width = 12, height = 12)
            k.diagnostics.plots(
                rcx.res,
                good.k.values,
                y,
                main.title = sprintf("Bandwidth Selection (p = %g): %s", single.p, data.tag)
            )
            grDevices::dev.off()
        }

        ## 3D embeddings at optimal k
        opt.k.idx <- which(good.k.values == opt.k)
        graph.at.opt <- graphs$geom_pruned_graphs[[good.indices[opt.k.idx]]]
        adj.list <- graph.at.opt$adj_list
        edgelen.list <- graph.at.opt$weight_list

        embed.3d <- list()
        if (n.embeddings > 0) {
            for (i in seq_len(n.embeddings)) {
                embed.3d[[i]] <- dgraphs::graph.embedding(
                    adj.list = adj.list,
                    weights.list = edgelen.list,
                    invert.weights = TRUE,
                    dim = 3,
                    method = "fr"
                )
            }
        }

        ## Assemble full results
        p.experiment <- list(
            p = single.p,
            p.tag = p.tag,
            status = "success",
            graphs = graphs,
            graphs.stats = graphs.stats,
            good.k = good.k,
            rcx.res = rcx.res,
            k.diag = k.diag,
            opt.k = opt.k,
            adj.list = adj.list,
            edgelen.list = edgelen.list,
            embed.3d = embed.3d
        )

        ## Save full results
        if (!is.null(p.rds.file)) {
            saveRDS(p.experiment, p.rds.file)
        }

        ## Return lightweight summary (with diagnostics across all k)
        list(
            p = single.p,
            p.tag = p.tag,
            status = "success",
            opt.k = opt.k,
            good.k = good.k,
            k.diag = k.diag,
            gcv = k.diag$gcv[k.diag$k == opt.k],
            n.extrema = k.diag$n.extrema[k.diag$k == opt.k],
            n.connected.k = length(good.k.values),
            file = p.rds.file
        )
    }

    ## ==========================================================================
    ## Execute computation (parallel or sequential)
    ## ==========================================================================
    if (n.cores > 1 && length(p.to.compute) > 1) {
        if (!requireNamespace("doParallel", quietly = TRUE)) {
            stop("Package 'doParallel' is required for parallel execution. ",
                 "Install it or set n.cores = 1.")
        }
        if (!requireNamespace("foreach", quietly = TRUE)) {
            stop("Package 'foreach' is required for parallel execution. ",
                 "Install it or set n.cores = 1.")
        }

        n.cores.use <- min(n.cores, length(p.to.compute))

        cl <- parallel::makeCluster(n.cores.use)
        doParallel::registerDoParallel(cl)
        on.exit(parallel::stopCluster(cl), add = TRUE)

        `%dopar%` <- foreach::`%dopar%`

        results <- foreach::foreach(
                                i = seq_along(p.to.compute),
                                .errorhandling = "pass",
                                .packages = c("gflowx", "dgraphs"),
                                .verbose = FALSE
                            ) %dopar% {
                                single.p <- p.to.compute[i]
                                p.tag <- p.tags.to.compute[i]

                                tryCatch(
                                    single.p.fn(
                                        single.p = single.p,
                                        p.tag = p.tag,
                                        X = X,
                                        y = y,
                                        k.min = k.min,
                                        k.max = k.max,
                                        L1.normalize = L1.normalize,
                                        out.dir = out.dir,
                                        plot.dir = plot.dir,
                                        data.tag = data.tag,
                                        skip.existing = skip.existing,
                                        iknn.params = iknn.params,
                                        regression.params = regression.params,
                                        n.embeddings = n.embeddings
                                    ),
                                    error = function(e) {
                                        list(
                                            p = single.p,
                                            p.tag = p.tag,
                                            status = "error",
                                            reason = conditionMessage(e)
                                        )
                                    }
                                )
                            }
        names(results) <- p.tags.to.compute
    } else {
        ## Sequential execution
        results <- list()
        for (i in seq_along(p.to.compute)) {
            single.p <- p.to.compute[i]
            p.tag <- p.tags.to.compute[i]

            if (verbose) {
                message("\n=== Processing p = ", single.p, " (", p.tag, ") ===")
            }

            results[[p.tag]] <- tryCatch(
                single.p.fn(
                    single.p = single.p,
                    p.tag = p.tag,
                    X = X,
                    y = y,
                    k.min = k.min,
                    k.max = k.max,
                    L1.normalize = L1.normalize,
                    out.dir = out.dir,
                    plot.dir = plot.dir,
                    data.tag = data.tag,
                    skip.existing = skip.existing,
                    iknn.params = iknn.params,
                    regression.params = regression.params,
                    n.embeddings = n.embeddings
                ),
                error = function(e) {
                    list(
                        p = single.p,
                        p.tag = p.tag,
                        status = "error",
                        reason = conditionMessage(e)
                    )
                }
            )

            if (verbose && results[[p.tag]]$status == "success") {
                message("  Optimal k: ", results[[p.tag]]$opt.k)
                message("  GCV: ", results[[p.tag]]$gcv)
                message("  Saved to: ", results[[p.tag]]$file)
            }
        }
    }

    ## Summary of results
    if (verbose) {
        message("\n=== Summary ===")
        status.counts <- table(sapply(results, function(x) x$status))
        print(status.counts)

        ## Report any failures
        failures <- Filter(function(x) x$status != "success", results)
        if (length(failures) > 0) {
            message("\nFailed computations:")
            for (f in failures) {
                message("  p = ", f$p, ": ", f$reason)
            }
        }
    }

    class(results) <- c("Lp.results", "list")
    invisible(results)
}

#' Load Results from Power Transformation Experiment
#'
#' Loads and combines results saved by \code{\link{compute.Lp.graphs}}
#' for multiple power transformation exponents.
#'
#' @param p.values Numeric vector of power transformation exponents whose
#'   results should be loaded.
#' @param out.dir Directory containing the saved RDS files.
#' @param data.tag Character string identifying the dataset (must match the
#'   tag used when results were saved).
#' @param L1.normalize Logical; if TRUE, load results generated with L1
#'   normalization. Must match the setting used in
#'   \code{\link{compute.Lp.graphs}}.
#'
#' @return A named list (keyed by p.tag, e.g., "p05", "p075") where each
#'   element contains the full results object saved by
#'   \code{\link{compute.Lp.graphs}}, including graphs, regression fits,
#'   and diagnostics. Returns an empty list if no results are found.
#'
#' @seealso \code{\link{compute.Lp.graphs}} for generating results,
#'   \code{\link{summary.Lp.results}} for extracting summary statistics.
#'
#' @examples
#' \dontrun{
#' p.grid <- seq(0.5, 1, by = 0.05)
#'
#' ## Load results without L1 normalization
#' results <- load.Lp.results(p.grid,
#'                            out.dir = "~/analysis/output",
#'                            data.tag = "Zambia",
#'                            L1.normalize = FALSE)
#'
#' ## Load results with L1 normalization
#' results.L1 <- load.Lp.results(p.grid,
#'                               out.dir = "~/analysis/output",
#'                               data.tag = "Zambia",
#'                               L1.normalize = TRUE)
#' }
#'
#' @export
load.Lp.results <- function(p.values,
                            out.dir,
                            data.tag,
                            L1.normalize = FALSE) {

    make.p.tag <- function(p) gsub("\\.", "", sprintf("p%g", p))

    make.rds.filename <- function(p.tag, data.tag, L1.normalize) {
        if (L1.normalize) {
            sprintf("Lp_experiment_L1_normalized_%s_%s.rds", p.tag, data.tag)
        } else {
            sprintf("Lp_experiment_%s_%s.rds", p.tag, data.tag)
        }
    }

    results <- list()
    for (p in p.values) {
        p.tag <- make.p.tag(p)
        rds.file <- file.path(out.dir, make.rds.filename(p.tag, data.tag, L1.normalize))
        if (file.exists(rds.file)) {
            results[[p.tag]] <- readRDS(rds.file)
        } else {
            warning(sprintf("Results for p = %g not found: %s", p, rds.file))
        }
    }

    class(results) <- c("Lp.results", "list")
    results
}

#' Summarize Power Transformation Experiment Results
#'
#' Extracts key diagnostic metrics from results loaded by
#' \code{\link{load.Lp.results}} into a data frame suitable for comparison
#' and visualization.
#'
#' @param object A list of results as returned by \code{\link{load.Lp.results}}.
#' @param ... Additional arguments (currently unused).
#'
#' @return A data frame with one row per successful experiment, containing:
#'   \describe{
#'     \item{p}{Power transformation exponent.}
#'     \item{p.tag}{Tag string for p (e.g., "p065").}
#'     \item{opt.k}{Optimal k selected by GCV.}
#'     \item{gcv}{GCV score at optimal k.}
#'     \item{n.extrema}{Number of local extrema at optimal k.}
#'     \item{n.connected.k}{Number of k values yielding connected graphs.}
#'   }
#'   Rows are ordered by increasing p. Returns NULL if no successful
#'   results are found.
#'
#' @seealso \code{\link{load.Lp.results}} for loading results,
#'   \code{\link{compute.Lp.graphs}} for generating results.
#'
#' @examples
#' \dontrun{
#' all.results <- load.Lp.results(seq(0.5, 1, by = 0.05))
#' summary.df <- summary(all.results)
#' print(summary.df)
#'
#' ## Visualize GCV across p values
#' plot(summary.df$p, summary.df$gcv, type = "b",
#'      xlab = "Power exponent p", ylab = "GCV")
#' }
#'
#' @export
summary.Lp.results <- function(object, ...) {

    if (length(object) == 0) {
        warning("No results to summarize")
        return(NULL)
    }

    if (!inherits(object, "Lp.results")) {
        stop("'object' must be of class 'Lp.results', as returned by load.Lp.results()")
    }

    summary.list <- lapply(object, function(res) {
        if (is.null(res) || is.null(res$status) || res$status != "success") {
            return(NULL)
        }
        data.frame(
            p = res$p,
            p.tag = res$p.tag,
            opt.k = res$opt.k,
            gcv = res$k.diag$gcv[res$k.diag$k == res$opt.k],
            n.extrema = res$k.diag$n.extrema[res$k.diag$k == res$opt.k],
            n.connected.k = length(res$good.k$k.values),
            stringsAsFactors = FALSE
        )
    })

    summary.df <- do.call(rbind, summary.list)

    if (is.null(summary.df) || nrow(summary.df) == 0) {
        warning("No successful results to summarize")
        return(NULL)
    }

    rownames(summary.df) <- NULL
    summary.df <- summary.df[order(summary.df$p), ]

    summary.df
}
