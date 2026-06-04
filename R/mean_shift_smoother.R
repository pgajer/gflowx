#' Mean Shift Data Smoother
#'
#' @description
#' Performs mean shift smoothing on a dataset using various algorithms. Mean shift is a
#' non-parametric feature-space analysis technique for locating the maxima of a density
#' function. This implementation offers multiple variants including adaptive step sizes,
#' gradient field averaging, and momentum-based updates.
#'
#' @param X A numeric matrix or data frame where rows represent points and columns
#'   represent features. Will be coerced to a matrix if necessary.
#'
#' @param k An integer specifying the number of nearest neighbors to consider for
#'   gradient estimation. Must be positive and less than the number of points.
#'
#' @param density.k An integer specifying the number of nearest neighbors to consider
#'   for density estimation. Must be positive and less than the number of points.
#'   Default is 1.
#'
#' @param n.steps An integer specifying the number of smoothing steps to perform.
#'   Must be positive. Default is 10.
#'
#' @param step.size A numeric value specifying the step size for updating point positions.
#'   Should be between 0 and 1. Default is 0.1.
#'
#' @param ikernel An integer specifying the kernel function to use:
#'   \itemize{
#'     \item 1: Epanechnikov kernel
#'     \item 2: Triangular kernel
#'     \item 3: Truncated Exponential kernel
#'     \item 4: Normal (Gaussian) kernel
#'   }
#'   Default is 1.
#'
#' @param dist.normalization.factor A numeric value specifying the scaling factor for
#'   distance normalization. Should be greater than 1. Default is 1.01.
#'
#' @param method A character string or integer specifying the smoothing method:
#'   \describe{
#'     \item{"basic" (0)}{Basic mean shift algorithm}
#'     \item{"precomputed" (1)}{Mean shift with precomputed nearest neighbors (default)}
#'     \item{"grad_field" (2)}{Mean shift with gradient field averaging}
#'     \item{"grad_field_dir" (3)}{Mean shift with gradient field averaging (direction only)}
#'     \item{"adaptive" (4)}{Adaptive mean shift}
#'     \item{"adaptive_dir" (5)}{Adaptive mean shift (direction only)}
#'     \item{"knn_adaptive" (6)}{KNN adaptive mean shift}
#'     \item{"knn_adaptive_grad" (7)}{KNN adaptive with gradient field averaging}
#'     \item{"knn_adaptive_grad_dir" (8)}{KNN adaptive with gradient field averaging (direction only)}
#'     \item{"adaptive_momentum" (9)}{Adaptive with gradient field averaging and momentum}
#'     \item{"adaptive_momentum_dir" (10)}{Adaptive with gradient field averaging and momentum (direction only)}
#'   }
#'   Default is "precomputed".
#'
#' @param average.direction.only A logical value indicating whether to average only the
#'   directions of gradients (TRUE) or full gradient vectors (FALSE). If NULL (default),
#'   this is automatically determined based on the method:
#'   \itemize{
#'     \item Methods ending in "_dir" use TRUE
#'     \item Other methods use FALSE
#'   }
#'
#' @param momentum A numeric value specifying the momentum factor for methods
#'   "adaptive_momentum" and "adaptive_momentum_dir". Should be between 0 and 1.
#'   Default is 0.9. Ignored for other methods.
#'
#' @param increase.factor A numeric value specifying the factor by which to increase
#'   the step size in adaptive methods. Should be greater than 1. Default is 1.2.
#'   Only used for adaptive methods (4-10).
#'
#' @param decrease.factor A numeric value specifying the factor by which to decrease
#'   the step size in adaptive methods. Should be between 0 and 1. Default is 0.5.
#'   Only used for adaptive methods (4-10).
#'
#' @return A list of class "MSD" containing:
#'   \describe{
#'     \item{X}{The original input dataset}
#'     \item{X.traj}{A list of matrices, each representing the smoothed dataset at each step}
#'     \item{median.kdistances}{A numeric vector of median k-distances for each step}
#'     \item{opt.step}{The step number with minimum median k-distance}
#'     \item{dX}{The smoothed dataset at the optimal step}
#'     \item{method}{The method used (as a string)}
#'     \item{params}{A list of all parameters used}
#'   }
#'
#' @details
#' The mean shift algorithm iteratively moves each data point towards the mode of the
#' density estimated in its neighborhood. Different variants offer various improvements:
#'
#' \itemize{
#'   \item \strong{Basic/Precomputed}: Standard mean shift with optional nearest neighbor precomputation
#'   \item \strong{Gradient Field}: Averages gradients across the dataset for smoother updates
#'   \item \strong{Adaptive}: Dynamically adjusts step sizes based on convergence behavior
#'   \item \strong{KNN Adaptive}: Uses k-nearest neighbors for adaptive bandwidth selection
#'   \item \strong{Momentum}: Incorporates previous update directions for faster convergence
#' }
#'
#' The "direction only" variants normalize gradient vectors before averaging, which can
#' be more robust to outliers but may converge more slowly.
#'
#' @examples
#' \dontrun{
#' # Generate sample data
#' set.seed(123)
#' X <- matrix(rnorm(200), ncol = 2)
#'
#' # Basic mean shift smoothing
#' result1 <- meanshift.data.smoother(X, k = 5)
#'
#' # Adaptive mean shift with gradient field averaging
#' result2 <- meanshift.data.smoother(X, k = 5, method = "adaptive")
#'
#' # With momentum
#' result3 <- meanshift.data.smoother(X, k = 5, method = "adaptive_momentum",
#'                                   momentum = 0.95)
#'
#' # Plot the trajectory
#' plot(result1$median.kdistances, type = 'l',
#'      xlab = 'Step', ylab = 'Median k-distance')
#' abline(v = result1$opt.step, col = 'red', lty = 2)
#'
#' # Compare original and smoothed data
#' par(mfrow = c(1, 2))
#' plot(X, main = "Original", pch = 19, cex = 0.5)
#' plot(result1$dX, main = "Smoothed", pch = 19, cex = 0.5)
#' }
#'
#' @references
#' Comaniciu, D., & Meer, P. (2002). Mean shift: A robust approach toward feature
#' space analysis. IEEE Transactions on Pattern Analysis and Machine Intelligence,
#' 24(5), 603-619.
#'
#' @seealso
#' \code{\link{plot.MSD}} for plotting the results
#'
#' @export
meanshift.data.smoother <- function(X,
                                    k,
                                    density.k = 1,
                                    n.steps = 10,
                                    step.size = 0.1,
                                    ikernel = 1,
                                    dist.normalization.factor = 1.01,
                                    method = "precomputed",
                                    momentum = 0.9,
                                    increase.factor = 1.2,
                                    decrease.factor = 0.5,
                                    average.direction.only = NULL
                                    ) {

    # Input validation for X
    if (!is.matrix(X)) {
        X <- try(as.matrix(X), silent = TRUE)
        if (inherits(X, "try-error")) {
            stop("X must be a matrix or coercible to a matrix")
        }
    }

    if (!is.numeric(X)) {
        stop("X must contain numeric values")
    }

    if (any(is.na(X)) || any(is.infinite(X))) {
        stop("X cannot contain NA, NaN, or Inf values")
    }

    if (!is.double(X)) {
        storage.mode(X) <- "double"
    }

    if (!requireNamespace("gflow", quietly = TRUE)) {
        stop("meanshift.data.smoother() currently requires gflow for its ",
             "native mean-shift backend.", call. = FALSE)
    }

    n.points <- nrow(X)
    n.features <- ncol(X)

    if (n.points < 2) {
        stop("X must contain at least 2 points")
    }
    if (n.features < 1) {
        stop("X must have at least 1 feature")
    }

    # Parameter validation
    if (!is.numeric(k) || length(k) != 1 || k != round(k) || k < 1 || k >= n.points) {
        stop("k must be a positive integer less than the number of points in X")
    }

    if (!is.numeric(density.k) || length(density.k) != 1 || density.k != round(density.k) ||
        density.k < 1 || density.k >= n.points) {
        stop("density.k must be a positive integer less than the number of points in X")
    }

    if (!is.numeric(n.steps) || length(n.steps) != 1 || n.steps != round(n.steps) ||
        n.steps < 1) {
        stop("n.steps must be a positive integer")
    }

    if (!is.numeric(step.size) || length(step.size) != 1 || step.size <= 0 || step.size > 1) {
        stop("step.size must be a numeric value between 0 and 1")
    }

    if (!is.numeric(ikernel) || length(ikernel) != 1 || ikernel != round(ikernel) ||
        ikernel < 1 || ikernel > 4) {
        stop("ikernel must be an integer between 1 and 4")
    }

    if (!is.numeric(dist.normalization.factor) || length(dist.normalization.factor) != 1 ||
        dist.normalization.factor <= 1) {
        stop("dist.normalization.factor must be a numeric value greater than 1")
    }

    # Method validation and conversion
    method.map <- c(
        "basic" = 0,
        "precomputed" = 1,
        "grad_field" = 2,
        "grad_field_dir" = 3,
        "adaptive" = 4,
        "adaptive_dir" = 5,
        "knn_adaptive" = 6,
        "knn_adaptive_grad" = 7,
        "knn_adaptive_grad_dir" = 8,
        "adaptive_momentum" = 9,
        "adaptive_momentum_dir" = 10
    )

    method.string <- method  # Store original method for output

    if (is.character(method)) {
        if (!(method %in% names(method.map))) {
            stop("Invalid method. Choose from: ", paste(names(method.map), collapse = ", "))
        }
        method.num <- method.map[method]
    } else if (is.numeric(method)) {
        if (!(method %in% 0:10)) {
            stop("method must be an integer between 0 and 10")
        }
        method.num <- method
        # Get string name for output
        method.string <- names(method.map)[which(method.map == method)]
    } else {
        stop("method must be a character string or integer")
    }

    # Determine average.direction.only if not specified
    if (is.null(average.direction.only)) {
        # Methods 3, 5, 8, 10 use direction only
        average.direction.only <- method.num %in% c(3, 5, 8, 10)
    }

    if (!is.logical(average.direction.only)) {
        stop("average.direction.only must be a logical value")
    }

    # Validate method-specific parameters
    adaptive.methods <- 4:10
    momentum.methods <- 9:10

    if (!(method.num %in% momentum.methods)) {
        if (!is.numeric(momentum) || length(momentum) != 1 || momentum <= 0 || momentum >= 1) {
            # Only warn if momentum was explicitly provided
            if (!missing(momentum)) {
                warning("momentum parameter is only used for methods 'adaptive_momentum' and 'adaptive_momentum_dir'")
            }
        }
    } else {
        if (!is.numeric(momentum) || length(momentum) != 1 || momentum <= 0 || momentum >= 1) {
            stop("momentum must be a numeric value between 0 and 1")
        }
    }

    if (!(method.num %in% adaptive.methods)) {
        if (!missing(increase.factor) || !missing(decrease.factor)) {
            warning("increase.factor and decrease.factor are only used for adaptive methods")
        }
    } else {
        if (!is.numeric(increase.factor) || length(increase.factor) != 1 || increase.factor <= 1) {
            stop("increase.factor must be a numeric value greater than 1")
        }

        if (!is.numeric(decrease.factor) || length(decrease.factor) != 1 ||
            decrease.factor <= 0 || decrease.factor >= 1) {
            stop("decrease.factor must be a numeric value between 0 and 1")
        }
    }

    # Call the appropriate C function based on method
    if (method.num == 0) {
        ## Basic method
        result <- .Call("S_mean_shift_data_smoother",
                        X,
                        as.integer(k),
                        as.integer(density.k),
                        as.integer(n.steps),
                        as.double(step.size),
                        as.integer(ikernel),
                        as.double(dist.normalization.factor),
                        as.integer(method.num),
                        as.double(momentum),
                        as.double(increase.factor),
                        as.double(decrease.factor),
                        PACKAGE = "gflow")

    } else if (method.num == 1) {
        ## Precomputed method
        stop("Precomputed method not implemented yet")
        ## result <- .Call("S_mean_shift_data_smoother_precomputed",
        ##                X,
        ##                as.integer(k),
        ##                as.integer(n.steps),
        ##                as.double(step.size),
        ##                as.integer(ikernel),
        ##                as.double(dist.normalization.factor))

    } else if (method.num %in% c(2, 3)) {
        # Gradient field averaging methods
        result <- .Call("S_mean_shift_data_smoother_with_grad_field_averaging",
                       X,
                       as.integer(k),
                       as.integer(density.k),
                       as.integer(n.steps),
                       as.double(step.size),
                       as.integer(ikernel),
                       as.double(dist.normalization.factor),
                       as.logical(average.direction.only),
                       PACKAGE = "gflow")
    } else if (method.num %in% c(4, 5)) {
        # Adaptive methods
        result <- .Call("S_mean_shift_data_smoother_adaptive",
                       X,
                       as.integer(k),
                       as.integer(density.k),
                       as.integer(n.steps),
                       as.double(step.size),
                       as.integer(ikernel),
                       as.double(dist.normalization.factor),
                       as.logical(average.direction.only),
                       PACKAGE = "gflow")
    } else if (method.num == 6) {
        ## KNN adaptive method
        stop("knn_adaptive_mean_shift_smoother not implemented yet")
        ## result <- .Call("S_knn_adaptive_mean_shift_smoother",
        ##                X,
        ##                as.integer(k),
        ##                as.integer(density.k),
        ##                as.integer(n.steps),
        ##                as.double(step.size),
        ##                as.integer(ikernel),
        ##                as.double(dist.normalization.factor),
        ##                as.double(increase.factor),
        ##                as.double(decrease.factor))

    } else if (method.num %in% c(7, 8)) {
        ## KNN adaptive with gradient field averaging
        result <- knn_adaptive_mean_shift_gfa(X,
                                              as.integer(k),
                                              as.integer(density.k),
                                              as.integer(n.steps),
                                              as.double(step.size),
                                              as.integer(ikernel),
                                              as.double(dist.normalization.factor),
                                              as.logical(average.direction.only))

    } else if (method.num %in% c(9, 10)) {
        ## Adaptive with gradient field averaging and momentum
        result <- adaptive_mean_shift_gfa(X,
                                          as.integer(k),
                                          as.integer(density.k),
                                          as.integer(n.steps),
                                          as.double(step.size),
                                          as.integer(ikernel),
                                          as.double(dist.normalization.factor),
                                          as.logical(average.direction.only),
                                          as.double(momentum),
                                          as.double(increase.factor),
                                          as.double(decrease.factor))
    }

    # Process results
    names(result) <- c("X.traj", "median.kdistances")

    # Find optimal step (minimum median k-distance, excluding first step)
    if (n.steps > 1) {
        steps <- 2:n.steps
        opt.step <- steps[which.min(result[["median.kdistances"]][steps])]
    } else {
        opt.step <- 1
    }

    # Store results
    result[["opt.step"]] <- opt.step
    result[["dX"]] <- result[["X.traj"]][[opt.step]]
    result[["X"]] <- X
    result[["method"]] <- method.string

    # Store all parameters for reproducibility
    result[["params"]] <- list(
        k = k,
        density.k = density.k,
        n.steps = n.steps,
        step.size = step.size,
        ikernel = ikernel,
        dist.normalization.factor = dist.normalization.factor,
        average.direction.only = average.direction.only,
        momentum = if (method.num %in% momentum.methods) momentum else NULL,
        increase.factor = if (method.num %in% adaptive.methods) increase.factor else NULL,
        decrease.factor = if (method.num %in% adaptive.methods) decrease.factor else NULL
    )

    class(result) <- "MSD"

    return(result)
}

#' Plot Results of Mean Shift Data Denoising
#'
#' This function plots the results of the Mean Shift Data Denoising function
#' \code{\link{meanshift.data.smoother}}. It can produce three types of plots based on the
#' 'type' parameter.
#'
#' @param x A list of class "MSD" containing the results from meanshift.data.smoother function.
#' @param type A character string specifying the type of plot. Must be one of
#'   "dX" (default), "dXi", or "kdists".
#' @param i An integer specifying which trajectory step to plot when type = "dXi".
#'   Default is 2.
#' @param rg A numeric value specifying the range for x and y axes in "dX" and
#'   "dXi" plots. If NULL (default), ranges are computed from the data.
#' @param mar Numeric vector of length 4 specifying plot margins.
#' @param mgp Numeric vector of length 3 for axis label positions.
#' @param tcl Numeric value for tick mark length.
#' @param ... Additional arguments to be passed to the plot function.
#'
#' @details
#' The function produces different plots based on the 'type' parameter:
#' \itemize{
#'   \item "dX": Plots original data (X) and denoised data (dX) side by side.
#'   \item "dXi": Plots original data (X) and a specific trajectory step (X.traj\\[\\[i\\]\\])
#'     side by side.
#'   \item "kdists": Plots median k-distances with a vertical line at the optimal step.
#' }
#'
#' @return This function does not return a value; it produces a plot as a side effect.
#'
#' @examples
#' \dontrun{
#' # Generate sample data
#' X <- matrix(rnorm(200), ncol = 2)
#'
#' # Run mean shift smoothing
#' res <- meanshift.data.smoother(X, k = 5)
#'
#' # Plot original vs denoised data
#' plot(res, type = "dX")
#'
#' # Plot trajectory at step 3
#' plot(res, type = "dXi", i = 3)
#'
#' # Plot median k-distances
#' plot(res, type = "kdists")
#' }
#'
#' @seealso \code{\link{meanshift.data.smoother}}
#'
#' @export
plot.MSD <- function(x,
                     type = "dX",
                     i = 2,
                     rg = NULL,
                     mar = c(2.5, 2.5, 2, 0.5),
                     mgp = c(2.5, 0.5, 0),
                     tcl = -0.3,
                     ...) {

    # Validate input
    if (!inherits(x, "MSD")) {
        stop("x must be an object of class 'MSD' from meanshift.data.smoother()")
    }

    # Use x directly (no need for renaming)
    res <- x

    # Validate type
    type <- match.arg(type, c("dX", "dXi", "kdists"))

    # Save and restore graphical parameters
    old.par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old.par), add = TRUE)

    switch(type,
        "dX" = {
            # Validate that dX exists
            if (is.null(res$dX)) {
                stop("Denoised data (dX) not found in results")
            }

            # Calculate range if not provided
            if (is.null(rg)) {
                rg <- max(abs(range(c(res$X, res$dX))))
                rg <- rg * 1.1  # Add 10% padding
            }

            graphics::par(mfrow = c(1, 2), mar = mar, mgp = mgp, tcl = tcl)

            # Plot original data
            graphics::plot(res$X,
                 main = "Original Data",
                 xlab = "X1",
                 ylab = "X2",
                 las = 1,
                 xlim = c(-rg, rg),
                 ylim = c(-rg, rg),
                 ...)

            # Plot denoised data
            graphics::plot(res$dX,
                 main = paste("Denoised Data (step", res$opt.step, ")"),
                 xlab = "X1",
                 ylab = "X2",
                 las = 1,
                 xlim = c(-rg, rg),
                 ylim = c(-rg, rg),
                 ...)
        },

        "dXi" = {
            # Validate trajectory index
            if (!is.numeric(i) || i < 1 || i > length(res$X.traj)) {
                stop("i must be between 1 and ", length(res$X.traj))
            }

            # Calculate range if not provided
            if (is.null(rg)) {
                rg <- max(abs(range(c(res$X, res$X.traj[[i]]))))
                rg <- rg * 1.1  # Add 10% padding
            }

            graphics::par(mfrow = c(1, 2), mar = mar, mgp = mgp, tcl = tcl)

            # Plot original data
            graphics::plot(res$X,
                 main = "Original Data",
                 xlab = "X1",
                 ylab = "X2",
                 las = 1,
                 xlim = c(-rg, rg),
                 ylim = c(-rg, rg),
                 ...)

            # Plot trajectory at step i
            graphics::plot(res$X.traj[[i]],
                 main = paste("Trajectory at step", i),
                 xlab = "X1",
                 ylab = "X2",
                 las = 1,
                 xlim = c(-rg, rg),
                 ylim = c(-rg, rg),
                 ...)
        },

        "kdists" = {
            # Validate median.kdistances exists
            if (is.null(res$median.kdistances)) {
                stop("median.kdistances not found in results")
            }

            # Single plot
            graphics::par(mar = c(4, 4, 2, 1), mgp = mgp, tcl = tcl)

            graphics::plot(res$median.kdistances,
                 main = "Median k-distances by Step",
                 xlab = "Step",
                 ylab = "Median k-distance",
                 las = 1,
                 type = 'b',
                 pch = 19,
                 ...)

            # Add vertical line at optimal step
            if (!is.null(res$opt.step)) {
                graphics::abline(v = res$opt.step, col = "red", lty = 2, lwd = 2)
                graphics::mtext(paste("Optimal step:", res$opt.step),
                                side = 3, line = 0.5, at = res$opt.step,
                                col = "red", cex = 0.8)
            }
        }
    )
}


#' Fully Adaptive Mean Shift with Gradient Field Averaging
#'
#' Implements a fully adaptive mean-shift smoother with gradient-field averaging.
#' This version adapts both neighborhood selection (kNN recomputed each step)
#' and per-point step sizes via momentum and step increase/decrease factors.
#'
#' @param X Numeric matrix (n x d), one row per point.
#' @param k Integer > 0, k-NN for gradient estimation (excl. self).
#' @param density_k Integer > 0, k-NN for density estimation (excl. self).
#' @param n_steps Integer > 0, number of iterations.
#' @param initial_step_size Positive numeric, initial per-point step size.
#' @param ikernel Integer kernel id; 1=Gaussian, 2=Epanechnikov, 3=Quartic.
#' @param dist_normalization_factor Positive numeric distance scaling (default 1.01).
#' @param average_direction_only Logical; if TRUE, average directions only.
#' @param momentum Numeric in (0,1], influence of previous gradients (default 0.9).
#' @param increase_factor Numeric > 1, step growth factor (default 1.2).
#' @param decrease_factor Numeric in (0,1), step shrink factor (default 0.5).
#'
#' @return A list with:
#' \itemize{
#'   \item \code{X_traj}: list of length \code{n_steps+1} (typically),
#'         each an n x d matrix of point positions per iteration.
#'   \item \code{median_kdistances}: numeric vector of median k-distances per step.
#' }
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(200), ncol = 2)
#' out <- adaptive_mean_shift_gfa(
#'   X, k = 10, density_k = 10, n_steps = 5, initial_step_size = 0.1
#' )
#' length(out$X_traj)
#' @export
adaptive_mean_shift_gfa <- function(
  X,
  k,
  density_k,
  n_steps,
  initial_step_size,
  ikernel = 1L,
  dist_normalization_factor = 1.01,
  average_direction_only = FALSE,
  momentum = 0.9,
  increase_factor = 1.2,
  decrease_factor = 0.5
) {
  if (!requireNamespace("gflow", quietly = TRUE)) {
    stop("adaptive_mean_shift_gfa() currently requires gflow for its native ",
         "mean-shift backend.", call. = FALSE)
  }
  X <- as.matrix(X)
  utils::getFromNamespace("rcpp_adaptive_mean_shift_gfa", "gflow")(
    X, as.integer(k), as.integer(density_k), as.integer(n_steps),
    as.numeric(initial_step_size),
    as.integer(ikernel), as.numeric(dist_normalization_factor),
    as.logical(average_direction_only),
    as.numeric(momentum), as.numeric(increase_factor), as.numeric(decrease_factor)
  )
}

#' kNN-Adaptive Mean Shift with Gradient Field Averaging
#'
#' Adaptive in the sense that kNN is recomputed as points move, but with a
#' constant step size over iterations.
#'
#' @inheritParams adaptive_mean_shift_gfa
#' @param step_size Positive numeric constant step size per iteration.
#'
#' @return A list with \code{X_traj} and \code{median_kdistances}; see
#'   \code{\link{adaptive_mean_shift_gfa}}.
#'
#' @examples
#' set.seed(1)
#' X <- matrix(rnorm(200), ncol = 2)
#' out <- knn_adaptive_mean_shift_gfa(
#'   X, k = 10, density_k = 10, n_steps = 5, step_size = 0.1
#' )
#' sapply(out$X_traj, dim)
#' @export
knn_adaptive_mean_shift_gfa <- function(
  X,
  k,
  density_k,
  n_steps,
  step_size,
  ikernel = 1L,
  dist_normalization_factor = 1.01,
  average_direction_only = FALSE
) {
  if (!requireNamespace("gflow", quietly = TRUE)) {
    stop("knn_adaptive_mean_shift_gfa() currently requires gflow for its ",
         "native mean-shift backend.", call. = FALSE)
  }
  X <- as.matrix(X)
  utils::getFromNamespace("rcpp_knn_adaptive_mean_shift_gfa", "gflow")(
    X, as.integer(k), as.integer(density_k), as.integer(n_steps),
    as.numeric(step_size),
    as.integer(ikernel), as.numeric(dist_normalization_factor),
    as.logical(average_direction_only)
  )
}
