#' Extract Optimal Power Transformation Graph
#'
#' Given GCV values across a grid of power transformation exponents,
#' finds the optimal p and extracts the corresponding graph, fitted values,
#' and diagnostics. The optimal p can be determined in three ways: (1) by
#' fitting an internal spline smoother to the GCV curve and finding its minimum,
#' (2) by providing a pre-fitted smoother object (including legacy "magelo"),
#' or (3) by directly specifying a target p value.
#'
#' When a target p is specified directly, the function bypasses smoother fitting
#' entirely and extracts results for the grid point closest to the specified
#' value. This is useful when the user has already identified the optimal p
#' through other means or wishes to examine results at a specific transformation
#' exponent.
#'
#' When a pre-fitted smoother object is provided, the function uses its smoothed
#' predictions to determine the optimal p, avoiding redundant computation when
#' the smoother has already been fit for exploratory purposes.
#'
#' @param summary.df Data frame with columns 'p' and 'gcv' (and optionally
#'     others). Required unless \code{p} is specified directly.
#' @param out.dir Directory containing Lp_experiment_*.rds files.
#' @param data.tag Tag identifying the dataset (e.g., "5k_Zambia").
#' @param p Numeric scalar specifying the target power transformation exponent
#'     directly. If provided, smoothing is bypassed and the function
#'     extracts results for the grid point closest to this value. Default is
#'     NULL.
#' @param magelo.fit A pre-fitted smoother object (e.g., a legacy object of
#'     class "magelo"). If provided (and \code{p} is NULL), this pre-fitted
#'     smoother is used instead of fitting a new one. The object must contain
#'     at minimum the fields \code{xgrid}, \code{gpredictions}, and
#'     \code{gpredictions.CrI}. Default is NULL.
#' @param L1.normalize If TRUE, extracts optimal graph using L1 normalized
#'     (after the power transformation) data.
#' @param plot.dir Directory for saving plots (NULL to skip plotting).
#' @param magelo.params List of optional smoothing parameters used when fitting
#'     a new spline smoother (\code{df}, \code{spar}, \code{grid.size}).
#'     Ignored if \code{p} or \code{magelo.fit} is provided.
#' @param save.results Logical; if TRUE, save extracted results to out.dir.
#' @param verbose Logical; print progress messages.
#'
#' @return A list containing:
#'   \item{p.opt}{Selected optimal p (closest grid point to target).}
#'   \item{p.opt.tag}{Tag string for optimal p.}
#'   \item{p.smoothed.min}{Smoothed GCV minimum location (NULL if p was
#'       specified directly).}
#'   \item{p.equivalent.range}{Range of p values with overlapping CrI (NULL if
#'       p was specified directly).}
#'   \item{k.opt}{Optimal k at selected p.}
#'   \item{gcv}{GCV at optimal (p, k).}
#'   \item{n.extrema}{Number of local extrema at optimal.}
#'   \item{fitted.values}{sPTB conditional expectation estimates.}
#'   \item{adj.list}{Adjacency list of optimal graph.}
#'   \item{edgelen.list}{Edge lengths of optimal graph.}
#'   \item{rcx.fit}{Full regression fit object at optimal.}
#'   \item{embed.3d}{List of 3D embeddings.}
#'   \item{magelo.fit}{The smoother object (NULL if p was specified directly).}
#'   \item{summary.df}{Input summary data frame (for reference).}
#'
#' @examples
#' \dontrun{
#' ## Example 1: Automatic optimal p selection via spline smoothing
#' opt.analysis <- extract.optimal.Lp.graph(
#'   summary.df = summary.df,
#'   out.dir = out.dir,
#'   data.tag = "5k_Zambia"
#' )
#'
#' ## Example 2: Re-use a pre-fitted smoother object
#' gcv.fit <- gflowx:::.gflowx.fit.gcv.smoother(summary.df$p, summary.df$gcv)
#' opt.analysis <- extract.optimal.Lp.graph(
#'   summary.df = summary.df,
#'   out.dir = out.dir,
#'   data.tag = "5k_Zambia",
#'   magelo.fit = gcv.fit
#' )
#'
#' ## Example 3: Directly specify target p
#' opt.analysis <- extract.optimal.Lp.graph(
#'   summary.df = summary.df,
#'   out.dir = out.dir,
#'   data.tag = "5k_Zambia",
#'   p = 0.25
#' )
#' }
#'
#' @export
extract.optimal.Lp.graph <- function(summary.df,
                                     out.dir,
                                     data.tag,
                                     p = NULL,
                                     magelo.fit = NULL,
                                     L1.normalize = FALSE,
                                     plot.dir = NULL,
                                     magelo.params = list(),
                                     save.results = TRUE,
                                     verbose = TRUE) {

    ## ================================================================
    ## Input validation
    ## ================================================================
    if (!dir.exists(out.dir)) {
        stop("Output directory does not exist: ", out.dir)
    }

    ## Validate p parameter if provided
    if (!is.null(p)) {
        if (!is.numeric(p) || length(p) != 1 || !is.finite(p)) {
            stop("'p' must be a single finite numeric value")
        }
        if (p <= 0) {
            stop("'p' must be positive")
        }
    }

    ## Validate magelo.fit if provided
    if (!is.null(magelo.fit)) {
        has.required.fields <- is.list(magelo.fit) &&
            all(c("xgrid", "gpredictions") %in% names(magelo.fit))
        if (!(inherits(magelo.fit, "magelo") || inherits(magelo.fit, "gflow_spline_fit") || has.required.fields)) {
            stop("'magelo.fit' must be a smoother object with fields 'xgrid' and 'gpredictions'")
        }
        if (inherits(magelo.fit, "magelo")) {
            warning("Legacy 'magelo' fit detected in 'magelo.fit'; support is maintained for compatibility but will be deprecated.",
                    call. = FALSE)
        }
        required.magelo.fields <- c("xgrid", "gpredictions")
        missing.fields <- setdiff(required.magelo.fields, names(magelo.fit))
        if (length(missing.fields) > 0) {
            stop("'magelo.fit' is missing required fields: ",
                 paste(missing.fields, collapse = ", "))
        }
        ## Check for CrI field (warn if missing but don't fail)
        if (!"gpredictions.CrI" %in% names(magelo.fit)) {
            warning("'magelo.fit' does not contain 'gpredictions.CrI'; ",
                    "equivalent p range will not be computed")
        }
    }

    ## Validate summary.df only when needed (not using direct p specification)
    if (is.null(p)) {
        if (missing(summary.df) || is.null(summary.df)) {
            stop("'summary.df' is required when 'p' is not specified directly")
        }
        required.cols <- c("p", "gcv")
        missing.cols <- setdiff(required.cols, names(summary.df))
        if (length(missing.cols) > 0) {
            stop("summary.df missing required columns: ",
                 paste(missing.cols, collapse = ", "))
        }
    }

    ## Determine mode of operation
    use.direct.p <- !is.null(p)
    use.prefitted.smoother <- !is.null(magelo.fit) && !use.direct.p

    ## ================================================================
    ## Determine optimal p based on mode
    ## ================================================================
    gcv.smooth <- NULL
    p.smoothed.min <- NULL
    p.equivalent.range <- NULL

    if (use.direct.p) {
        ## Mode 1: Direct p specification - bypass smoother fitting entirely
        if (verbose) {
            cat("Using directly specified p =", p, "\n")
        }

        ## Determine p grid from available files if summary.df not provided
        if (missing(summary.df) || is.null(summary.df)) {
            ## Extract available p values from filenames
            pattern <- if (L1.normalize) {
                           sprintf("Lp_experiment_L1_normalized_p.*_%s\\.rds$", data.tag)
                       } else {
                           sprintf("Lp_experiment_p.*_%s\\.rds$", data.tag)
                       }
            available.files <- list.files(out.dir, pattern = pattern)
            if (length(available.files) == 0) {
                stop("No Lp experiment files found in: ", out.dir)
            }
            ## Extract p values from filenames
            p.pattern <- "p([0-9]+)"
            p.matches <- regmatches(available.files, regexpr(p.pattern, available.files))
            p.grid <- as.numeric(gsub("p", "0.", p.matches))
            p.grid <- sort(unique(p.grid))
            if (verbose) {
                cat("  Found", length(p.grid), "available p values:",
                    paste(head(p.grid, 5), collapse = ", "),
                    if (length(p.grid) > 5) "...", "\n")
            }
        } else {
            p.grid <- summary.df$p
        }

        ## Find closest grid point
        p.opt.idx <- which.min(abs(p.grid - p))
        p.opt <- p.grid[p.opt.idx]

        if (verbose) {
            cat("  Closest grid point: p =", p.opt, "\n")
            cat("  Distance from target:", round(abs(p.opt - p), 5), "\n")
        }

    } else if (use.prefitted.smoother) {
        ## Mode 2: Use pre-fitted smoother object
        if (verbose) cat("Using pre-fitted smoother object...\n")

        gcv.smooth <- magelo.fit

        ## Extract smoothed minimum
        min.idx <- which.min(gcv.smooth$gpredictions)
        p.smoothed.min <- gcv.smooth$xgrid[min.idx]

        if (verbose) {
            cat("  Smoothed GCV minimum at p =", round(p.smoothed.min, 4), "\n")
        }

        ## Find closest grid point to smoothed minimum
        p.grid <- summary.df$p
        p.opt.idx <- which.min(abs(p.grid - p.smoothed.min))
        p.opt <- p.grid[p.opt.idx]

        if (verbose) {
            cat("  Closest grid point: p =", p.opt, "\n")
            cat("  Distance from smoothed min:",
                round(abs(p.opt - p.smoothed.min), 5), "\n")
        }

        ## Compute equivalent p range using credible intervals if available
        if ("gpredictions.CrI" %in% names(gcv.smooth) &&
            !is.null(gcv.smooth$gpredictions.CrI)) {
            gcv.lower <- gcv.smooth$gpredictions.CrI[1, ]
            gcv.upper <- gcv.smooth$gpredictions.CrI[2, ]
            gcv.upper.at.min <- gcv.upper[min.idx]

            ## p values whose lower CI is below upper CI at minimum (overlapping CrI)
            equivalent.to.min <- gcv.lower <= gcv.upper.at.min
            p.equivalent.range <- range(gcv.smooth$xgrid[equivalent.to.min])

            if (verbose) {
                cat("  Equivalent p range (CrI overlap): [",
                    round(p.equivalent.range[1], 3), ",",
                    round(p.equivalent.range[2], 3), "]\n")
            }
        }

    } else {
        ## Mode 3: Fit new spline smoother
        if (verbose) cat("Fitting spline smoother to GCV curve...\n")

        gcv.smooth <- .gflowx.fit.gcv.smoother(
            x = summary.df$p,
            y = summary.df$gcv,
            params = magelo.params
        )

        ## Extract smoothed minimum
        min.idx <- which.min(gcv.smooth$gpredictions)
        p.smoothed.min <- gcv.smooth$xgrid[min.idx]

        if (verbose) {
            cat("  Smoothed GCV minimum at p =", round(p.smoothed.min, 4), "\n")
        }

        ## Find closest grid point to smoothed minimum
        p.grid <- summary.df$p
        p.opt.idx <- which.min(abs(p.grid - p.smoothed.min))
        p.opt <- p.grid[p.opt.idx]

        if (verbose) {
            cat("  Closest grid point: p =", p.opt, "\n")
            cat("  Distance from smoothed min:",
                round(abs(p.opt - p.smoothed.min), 5), "\n")
        }

        ## Compute equivalent p range using interval overlap if available
        if ("gpredictions.CrI" %in% names(gcv.smooth) &&
            !is.null(gcv.smooth$gpredictions.CrI)) {
            gcv.lower <- gcv.smooth$gpredictions.CrI[1, ]
            gcv.upper <- gcv.smooth$gpredictions.CrI[2, ]
            gcv.upper.at.min <- gcv.upper[min.idx]

            ## p values whose lower CI is below upper CI at minimum (overlapping CrI)
            equivalent.to.min <- gcv.lower <= gcv.upper.at.min
            p.equivalent.range <- range(gcv.smooth$xgrid[equivalent.to.min])

            if (verbose) {
                cat("  Equivalent p range (CrI overlap): [",
                    round(p.equivalent.range[1], 3), ",",
                    round(p.equivalent.range[2], 3), "]\n")
            }
        }
    }

    ## ================================================================
    ## Generate p tag and load results file
    ## ================================================================
    make.p.tag <- function(p.val) gsub("\\.", "", sprintf("p%g", p.val))
    p.opt.tag <- make.p.tag(p.opt)

    if (verbose) {
        cat("  Tag:", p.opt.tag, "\n")
    }

    ## ================================================================
    ## Load full results for optimal p
    ## ================================================================
    make.rds.filename <- function(p.tag, d.tag, L1.norm) {
        if (L1.norm) {
            sprintf("Lp_experiment_L1_normalized_%s_%s.rds", p.tag, d.tag)
        } else {
            sprintf("Lp_experiment_%s_%s.rds", p.tag, d.tag)
        }
    }

    opt.results.file <- file.path(out.dir,
                                  make.rds.filename(p.opt.tag, data.tag,
                                                    L1.normalize))

    if (!file.exists(opt.results.file)) {
        stop("Results file not found: ", opt.results.file,
             "\nAvailable p values may not include p = ", p.opt)
    }

    if (verbose) cat("Loading results from:", basename(opt.results.file), "\n")

    opt.experiment <- readRDS(opt.results.file)

    if (opt.experiment$status != "success") {
        stop("Experiment at p = ", p.opt, " has status: ", opt.experiment$status)
    }

    ## ================================================================
    ## Extract key components
    ## ================================================================
    ##     p.experiment <- list(
    ##     p = single.p,
    ##     p.tag = p.tag,
    ##     status = "success",
    ##     graphs = graphs,
    ##     graphs.stats = graphs.stats,
    ##     good.k = good.k,
    ##     rcx.res = rcx.res,
    ##     k.diag = k.diag,
    ##     opt.k = opt.k,
    ##     adj.list = adj.list,
    ##     edgelen.list = edgelen.list,
    ##     embed.3d = embed.3d
    ## )

    k.opt <- opt.experiment$opt.k
    k.opt.idx <- which(opt.experiment$good.k$k.values == k.opt)

    ## Graph structure
    opt.graph <- opt.experiment$graphs$geom_pruned_graphs[[
                                                              opt.experiment$good.k$indices[k.opt.idx]
                                                              ]]
    adj.list <- opt.graph$adj_list
    edgelen.list <- opt.graph$weight_list

    ## Regression fit and fitted values
    rcx.fit <- opt.experiment$rcx.res[[k.opt.idx]]
    fitted.values <- rcx.fit$fitted.values
    rel.fitted.values <- fitted.values / mean(rcx.fit$y)
    r.winsorized.rel.fitted.values <- right.winsorize(rel.fitted.values, p = 0.025)
    winsorized.rel.fitted.values <- winsorize(rel.fitted.values, p = 0.025)

    ## Diagnostics
    k.diag <- opt.experiment$k.diag
    gcv.opt <- k.diag$gcv[k.diag$k == k.opt]
    n.extrema.opt <- k.diag$n.extrema[k.diag$k == k.opt]

    ## Embeddings
    embed.3d <- opt.experiment$embed.3d

    if (verbose) {
        cat("\nOptimal configuration:\n")
        cat("  p =", p.opt, ", k =", k.opt, "\n")
        cat("  GCV =", format(gcv.opt, scientific = TRUE, digits = 4), "\n")
        cat("  Local extrema:", n.extrema.opt, "\n")
        cat("  Fitted values range: [",
            round(min(fitted.values), 4), ",",
            round(max(fitted.values), 4), "]\n")
        cat("  Relative fitted values range: [",
            round(min(rel.fitted.values), 4), ",",
            round(max(rel.fitted.values), 4), "]\n")
        cat("  Winsorized relative fitted values range: [",
            round(min(winsorized.rel.fitted.values), 4), ",",
            round(max(winsorized.rel.fitted.values), 4), "]\n")
        cat("  Graph: n =", length(adj.list), "vertices,",
            sum(sapply(adj.list, length)) / 2, "edges\n")
    }

    ## ================================================================
    ## Assemble output object
    ## ================================================================
    result <- list(
        ## Optimal parameters
        p.opt = p.opt,
        p.opt.tag = p.opt.tag,
        p.smoothed.min = p.smoothed.min,
        p.equivalent.range = p.equivalent.range,
        k.opt = k.opt,

        ## Diagnostics at optimal
        gcv = gcv.opt,
        n.extrema = n.extrema.opt,

        ## Graph structure
        adj.list = adj.list,
        edgelen.list = edgelen.list,
        n.vertices = length(adj.list),
        n.edges = sum(sapply(adj.list, length)) / 2,

        ## Regression results
        fitted.values = fitted.values,
        rel.fitted.values = rel.fitted.values,
        r.winsorized.rel.fitted.values = r.winsorized.rel.fitted.values,
        winsorized.rel.fitted.values = winsorized.rel.fitted.values,
        rcx.fit = rcx.fit,

        ## Embeddings
        embed.3d = embed.3d,

        ## Smoother results (for reproducibility and further analysis)
        magelo.fit = gcv.smooth,

        ## Reference data
        summary.df = if (!missing(summary.df)) summary.df else NULL,
        data.tag = data.tag,
        results.file = opt.results.file,

        ## Mode information
        mode = if (use.direct.p) "direct_p" else if (use.prefitted.smoother)
                                                "prefitted_smoother" else "fitted_spline"
    )

    ## ================================================================
    ## Generate plots if requested
    ## ================================================================
    if (!is.null(plot.dir) && !is.null(gcv.smooth)) {
        if (!dir.exists(plot.dir)) {
            dir.create(plot.dir, recursive = TRUE)
        }

        ## Plot 1: GCV curve with optimal region
        pdf.gcv <- file.path(plot.dir,
                             sprintf("GCV_optimal_p_%s.pdf", data.tag))
        pdf(pdf.gcv, width = 8, height = 6)

        op <- par(mar = c(4, 5.5, 2, 1), mgp = c(2.5, 0.5, 0), tcl = -0.3)

        if (inherits(gcv.smooth, "magelo")) {
            plot(gcv.smooth, with.pts = TRUE, pts.col = "red",
                 xlab = "p", ylab = "", legend.position = "topleft")
        } else {
            plot(summary.df$p, summary.df$gcv,
                 xlab = "p", ylab = "", pch = 16, cex = 0.8,
                 col = "red")
            lines(gcv.smooth$xgrid, gcv.smooth$gpredictions, col = "black", lwd = 2)
        }
        mtext("GCV", 2, line = 4.5)

        ## Shade equivalent region if available
        if (!is.null(p.equivalent.range)) {
            usr <- par("usr")
            rect(p.equivalent.range[1], usr[3],
                 p.equivalent.range[2], usr[4],
                 col = rgb(0, 0.5, 0, 0.1), border = NA)
        }

        ## Mark smoothed minimum and selected grid point
        if (!is.null(p.smoothed.min)) {
            abline(v = p.smoothed.min, col = "darkgreen", lty = 2, lwd = 1.5)
            mtext(round(p.smoothed.min, 3), 1, at = p.smoothed.min,
                  line = 0.5, cex = 0.8)
        }
        abline(v = p.opt, col = "blue", lty = 3, lwd = 1.5)

        legend.items <- c("Selected grid point")
        legend.cols <- c("blue")
        legend.ltys <- c(3)
        legend.pchs <- c(NA)

        if (!is.null(p.smoothed.min)) {
            legend.items <- c("Smoothed minimum", legend.items)
            legend.cols <- c("darkgreen", legend.cols)
            legend.ltys <- c(2, legend.ltys)
            legend.pchs <- c(NA, legend.pchs)
        }

        if (!is.null(p.equivalent.range)) {
            legend.items <- c(legend.items, "Equivalent region")
            legend.cols <- c(legend.cols, rgb(0, 0.5, 0, 0.5))
            legend.ltys <- c(legend.ltys, NA)
            legend.pchs <- c(legend.pchs, 15)
        }

        legend("topright",
               legend.items,
               col = legend.cols,
               lty = legend.ltys, lwd = 1.5,
               pch = legend.pchs, pt.cex = 2,
               cex = 0.8, bg = "white")

        par(op)
        dev.off()

        if (verbose) cat("\nGCV plot saved to:", pdf.gcv, "\n")

        ## Plot 2: Fitted values diagnostics
        pdf.fitted <- file.path(plot.dir,
                                sprintf("sptb_fitted_optimal_p_%s.pdf", data.tag))
        pdf(pdf.fitted, width = 10, height = 4)

        par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))

        ## Need y for this plot - extract from rcx.fit if available
        if (!is.null(rcx.fit$y)) {
            y <- rcx.fit$y
        } else {
            ## Attempt to infer from fitted values (won't work, just placeholder)
            y <- NULL
        }

        ## Histogram of fitted values
        hist(fitted.values, breaks = 30, col = "lightblue", border = "white",
             main = sprintf("Fitted sPTB Risk (p = %g, k = %d)", p.opt, k.opt),
             xlab = "Fitted probability", ylab = "Frequency")
        if (!is.null(y)) {
            abline(v = mean(y), col = "red", lty = 2, lwd = 2)
            legend("topright", "Observed prevalence", col = "red", lty = 2, lwd = 2)
        }

        ## Fitted vs observed (if y available)
        if (!is.null(y)) {
            plot(jitter(y, amount = 0.03), fitted.values,
                 pch = 16, col = rgb(0, 0, 0, 0.3), cex = 0.8,
                 xlab = "Observed sPTB (jittered)", ylab = "Fitted probability",
                 main = "Fitted vs Observed")
            abline(h = mean(fitted.values[y == 0]), col = "blue", lty = 2)
            abline(h = mean(fitted.values[y == 1]), col = "red", lty = 2)
            legend("topleft", c("Mean (controls)", "Mean (cases)"),
                   col = c("blue", "red"), lty = 2, cex = 0.8)
        } else {
            plot(fitted.values, pch = 16, col = rgb(0, 0, 0, 0.3),
                 xlab = "Sample index", ylab = "Fitted probability",
                 main = "Fitted Values")
        }

        dev.off()

        if (verbose) cat("Fitted values plot saved to:", pdf.fitted, "\n")
    } else if (!is.null(plot.dir) && is.null(gcv.smooth)) {
        if (verbose) {
            cat("\nNote: GCV plot not generated (no smoother fit available).\n")
        }
    }

    ## ================================================================
    ## Save results if requested
    ## ================================================================
    if (save.results) {
        save.file <- file.path(out.dir,
                               sprintf("sptb_optimal_Lp_analysis_%s_%s.rds", p.opt.tag, data.tag))
        saveRDS(result, save.file)
        if (verbose) cat("\nResults saved to:", save.file, "\n")
        result$save.file <- save.file
    }

    invisible(result)
}
