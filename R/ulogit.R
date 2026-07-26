#' Fit Univariate Logistic Regression Model
#'
#' @description
#' Fits a univariate logistic regression model to binary response data using
#' Newton-Raphson optimization with optional ridge regularization. The function
#' implements numerical safeguards to ensure stable convergence and provides
#' leave-one-out cross-validation errors for model assessment.
#'
#' @param x Numeric vector of predictor values. Must contain only finite values
#'   (no NA, NaN, or Inf).
#' @param y Binary response vector containing only 0 or 1 values. Must be the
#'   same length as \code{x}.
#' @param w Optional numeric vector of non-negative observation weights. If
#'   \code{NULL} (default), all observations are given equal weight. Must be
#'   the same length as \code{x} if provided.
#' @param max.iterations Maximum number of iterations for Newton-Raphson
#'   optimization. Must be a positive integer. Default is 100.
#' @param ridge.lambda Ridge regularization parameter to improve numerical
#'   stability. Must be non-negative. Default is 0.002. Higher values provide
#'   more regularization but may increase bias.
#' @param max.beta Maximum allowed absolute value for coefficient estimates.
#'   Used to prevent numerical overflow. Must be positive. Default is 100.0.
#' @param tolerance Convergence tolerance for optimization. The algorithm stops
#'   when the relative change in log-likelihood is less than this value. Must
#'   be positive. Default is 1e-8.
#' @param verbose Logical flag to enable detailed output during optimization.
#'   If \code{TRUE}, prints iteration progress. Default is \code{FALSE}.
#'
#' @return A list of class \code{"ulogit"} containing:
#' \describe{
#'   \item{predictions}{Numeric vector of fitted probabilities for each observation}
#'   \item{errors}{Numeric vector of leave-one-out cross-validation errors}
#'   \item{weights}{Numeric vector of weights used in model fitting}
#'   \item{converged}{Logical indicating whether the algorithm converged}
#'   \item{iterations}{Integer giving the number of iterations used}
#'   \item{call}{The matched call}
#' }
#'
#' @details
#' The function fits a logistic regression model of the form:
#' \deqn{logit(p_i) = \beta_0 + \beta_1 x_i}
#' where \eqn{p_i} is the probability of success for observation \eqn{i}.
#'
#' The Newton-Raphson algorithm is used for maximum likelihood estimation with
#' optional ridge regularization. The ridge penalty adds \eqn{\lambda \beta^2}
#' to the negative log-likelihood, which helps stabilize the estimation when
#' the data are nearly separable.
#'
#' Leave-one-out cross-validation (LOOCV) errors are computed efficiently using
#' the hat matrix diagonal elements, providing a measure of predictive accuracy
#' without requiring repeated model fitting.
#'
#' @references
#' Hastie, T., Tibshirani, R., & Friedman, J. (2009). The Elements of
#' Statistical Learning (2nd ed.). Springer.
#'
#' @examples
#' # Basic usage with simulated data
#' set.seed(123)
#' x <- seq(0, 1, length.out = 100)
#' true_prob <- 1/(1 + exp(-(2*x - 1)))
#' y <- rbinom(100, 1, prob = true_prob)
#' fit <- ulogit(x, y)
#'
#' # Plot results
#' plot(x, y, pch = 16, col = ifelse(y == 1, "blue", "red"),
#'      main = "Univariate Logistic Regression")
#' lines(x, fit$predictions, lwd = 2)
#' legend("topleft", c("y = 1", "y = 0", "Fitted"),
#'        col = c("blue", "red", "black"),
#'        pch = c(16, 16, NA), lty = c(NA, NA, 1))
#'
#' # Example with weights
#' w <- runif(100, 0.5, 1.5)
#' fit_weighted <- ulogit(x, y, w = w)
#'
#' # Example with increased regularization
#' fit_regularized <- ulogit(x, y, ridge.lambda = 0.1)
#'
#' # Compare LOOCV errors
#' cat("Standard model LOOCV error:", mean(fit$errors), "\n")
#' cat("Weighted model LOOCV error:", mean(fit_weighted$errors), "\n")
#' cat("Regularized model LOOCV error:", mean(fit_regularized$errors), "\n")
#'
#' @seealso
#' \code{\link{eigen.ulogit}} for an alternative implementation using Eigen,
#' \code{\link[stats]{glm}} for multivariate logistic regression
#'
#' @export
ulogit <- function(x,
                   y,
                   w = NULL,
                   max.iterations = 100L,
                   ridge.lambda = 0.002,
                   max.beta = 100.0,
                   tolerance = 1e-8,
                   verbose = FALSE) {

    # Store the call
    cl <- match.call()

    # Input validation with informative error messages
    if (!is.numeric(x)) {
        stop("'x' must be a numeric vector", call. = FALSE)
    }
    if (!is.numeric(y)) {
        stop("'y' must be a numeric vector", call. = FALSE)
    }
    if (length(x) != length(y)) {
        stop(sprintf("'x' and 'y' must have the same length (x: %d, y: %d)",
                     length(x), length(y)), call. = FALSE)
    }
    if (length(x) < 3L) {
        stop("At least 3 observations are required for model fitting", call. = FALSE)
    }
    if (!all(y %in% c(0, 1))) {
        stop("'y' must contain only 0 or 1 values", call. = FALSE)
    }
    if (!all(is.finite(x))) {
        stop("'x' contains non-finite values (NA, NaN, or Inf)", call. = FALSE)
    }

    # Parameter validation
    max.iterations <- as.integer(max.iterations)
    if (length(max.iterations) != 1L || is.na(max.iterations) || max.iterations <= 0L) {
        stop("'max.iterations' must be a single positive integer", call. = FALSE)
    }

    if (!is.numeric(ridge.lambda) || length(ridge.lambda) != 1L ||
        is.na(ridge.lambda) || ridge.lambda < 0) {
        stop("'ridge.lambda' must be a single non-negative number", call. = FALSE)
    }

    if (!is.numeric(max.beta) || length(max.beta) != 1L ||
        is.na(max.beta) || max.beta <= 0) {
        stop("'max.beta' must be a single positive number", call. = FALSE)
    }

    if (!is.numeric(tolerance) || length(tolerance) != 1L ||
        is.na(tolerance) || tolerance <= 0) {
        stop("'tolerance' must be a single positive number", call. = FALSE)
    }

    if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
        stop("'verbose' must be a single logical value (TRUE or FALSE)", call. = FALSE)
    }

    # Weight validation and processing
    if (is.null(w)) {
        w <- rep(1.0, length(x))
    } else {
        if (!is.numeric(w)) {
            stop("'w' must be a numeric vector", call. = FALSE)
        }
        if (length(w) != length(x)) {
            stop(sprintf("'w' must have the same length as 'x' (w: %d, x: %d)",
                         length(w), length(x)), call. = FALSE)
        }
        if (!all(is.finite(w))) {
            stop("'w' contains non-finite values (NA, NaN, or Inf)", call. = FALSE)
        }
        if (any(w < 0)) {
            stop("All weights must be non-negative", call. = FALSE)
        }
        if (all(w == 0)) {
            stop("At least one weight must be positive", call. = FALSE)
        }
        # Normalize weights for numerical stability
        w <- w / sum(w) * length(w)
    }

    # Call the C implementation
    result <- .Call("S_ulogit",
                    as.double(x),
                    as.double(y),
                    as.double(w),
                    as.integer(max.iterations),
                    as.double(ridge.lambda),
                    as.double(max.beta),
                    as.double(tolerance),
                    isTRUE(verbose),
                    PACKAGE = "gflowx")
    
    # Add call to result
    result$call <- cl

    # Set class for S3 methods
    class(result) <- "ulogit"

    return(result)
}


#' Fit Univariate Logistic Regression Using Eigen
#'
#' @description
#' Fits a univariate logistic regression model using the Eigen linear algebra
#' library for enhanced numerical stability. This function supports both linear
#' and quadratic models with optional observation weights and provides efficient
#' leave-one-out cross-validation.
#'
#' @param x Numeric vector of predictor values. Must contain only finite values
#'   (no NA, NaN, or Inf).
#' @param y Binary response vector containing only 0 or 1 values. Must be the
#'   same length as \code{x}.
#' @param w Optional numeric vector of non-negative observation weights. If
#'   \code{NULL} (default), all observations are given equal weight. Must be
#'   the same length as \code{x} if provided.
#' @param fit.quadratic Logical; if \code{TRUE}, includes a quadratic term
#'   (\eqn{x^2}) in the model. Default is \code{FALSE}.
#' @param with.errors Logical indicating whether to compute leave-one-out
#'   cross-validation errors. Default is \code{TRUE}. Setting to \code{FALSE}
#'   can improve performance for large datasets.
#' @param max.iterations Maximum number of Newton-Raphson iterations. Must be
#'   a positive integer. Default is 100.
#' @param ridge.lambda Ridge regularization parameter for numerical stability.
#'   Must be non-negative. Default is 0.002. Applied to all non-intercept
#'   coefficients.
#' @param tolerance Convergence tolerance for Newton-Raphson algorithm. The
#'   algorithm stops when the relative change in log-likelihood is less than
#'   this value. Must be positive. Default is 1e-8.
#'
#' @return A list of class \code{"eigen.ulogit"} containing:
#' \describe{
#'   \item{predictions}{Numeric vector of fitted probabilities}
#'   \item{converged}{Logical indicating whether the algorithm converged}
#'   \item{iterations}{Integer giving the number of iterations used}
#'   \item{beta}{Numeric vector of coefficients: intercept, linear coefficient,
#'     and (if \code{fit.quadratic = TRUE}) quadratic coefficient}
#'   \item{errors}{If \code{with.errors = TRUE}, numeric vector of leave-one-out
#'     cross-validation errors; otherwise \code{NULL}}
#'   \item{loglik}{Final log-likelihood value}
#'   \item{aic}{Akaike Information Criterion}
#'   \item{bic}{Bayesian Information Criterion}
#'   \item{call}{The matched call}
#'   \item{model}{Character string indicating model type ("linear" or "quadratic")}
#' }
#'
#' @details
#' For the linear model (\code{fit.quadratic = FALSE}), the function fits:
#' \deqn{logit(p_i) = \beta_0 + \beta_1 x_i}
#'
#' For the quadratic model (\code{fit.quadratic = TRUE}), it fits:
#' \deqn{logit(p_i) = \beta_0 + \beta_1 x_i + \beta_2 x_i^2}
#'
#' The implementation uses the Eigen C++ library for linear algebra operations,
#' providing enhanced numerical stability compared to standard BLAS/LAPACK
#' routines, particularly for ill-conditioned problems.
#'
#' Ridge regularization is applied to non-intercept coefficients to prevent
#' overfitting and improve numerical stability. The penalty term added to the
#' negative log-likelihood is \eqn{\lambda (\beta_1^2 + \beta_2^2)}.
#'
#' @references
#' McCullagh, P., & Nelder, J. A. (1989). Generalized Linear Models (2nd ed.).
#' Chapman and Hall/CRC.
#'
#' @examples
#' # Generate example data with non-linear relationship
#' set.seed(456)
#' x <- runif(150, -3, 3)
#' true_prob <- 1/(1 + exp(-(1 + 2*x - 0.5*x^2)))
#' y <- rbinom(150, 1, prob = true_prob)
#'
#' # Fit linear model
#' fit_linear <- eigen.ulogit(x, y)
#'
#' # Fit quadratic model
#' fit_quad <- eigen.ulogit(x, y, fit.quadratic = TRUE)
#'
#' # Compare models using AIC
#' cat("Linear model AIC:", fit_linear$aic, "\n")
#' cat("Quadratic model AIC:", fit_quad$aic, "\n")
#'
#' # Plot both fits
#' ord <- order(x)
#' plot(x, y, pch = 16, col = ifelse(y == 1, "blue", "red"),
#'      main = "Linear vs Quadratic Logistic Regression")
#' lines(x[ord], fit_linear$predictions[ord], lwd = 2, col = "green")
#' lines(x[ord], fit_quad$predictions[ord], lwd = 2, col = "purple")
#' legend("topleft", c("y = 1", "y = 0", "Linear", "Quadratic"),
#'        col = c("blue", "red", "green", "purple"),
#'        pch = c(16, 16, NA, NA), lty = c(NA, NA, 1, 1))
#'
#' # Example with weights and without errors for efficiency
#' w <- runif(150, 0.5, 2)
#' fit_fast <- eigen.ulogit(x, y, w = w, with.errors = FALSE)
#'
#' # Access coefficients
#' cat("Quadratic model coefficients:\n")
#' cat("  Intercept:", fit_quad$beta[1], "\n")
#' cat("  Linear:", fit_quad$beta[2], "\n")
#' cat("  Quadratic:", fit_quad$beta[3], "\n")
#'
#' @seealso
#' \code{\link{ulogit}} for the standard implementation,
#' \code{\link[stats]{glm}} for general linear models,
#' \code{\link{predict.eigen.ulogit}} for predictions on new data
#'
#' @export
eigen.ulogit <- function(x,
                         y,
                         w = NULL,
                         fit.quadratic = FALSE,
                         with.errors = TRUE,
                         max.iterations = 100L,
                         ridge.lambda = 0.002,
                         tolerance = 1e-8) {

    # Store the call
    cl <- match.call()

    # Input validation
    if (!is.numeric(x)) {
        stop("'x' must be a numeric vector", call. = FALSE)
    }
    if (!is.numeric(y)) {
        stop("'y' must be a numeric vector", call. = FALSE)
    }
    if (length(x) != length(y)) {
        stop(sprintf("'x' and 'y' must have the same length (x: %d, y: %d)",
                     length(x), length(y)), call. = FALSE)
    }

    # Check minimum sample size
    min_n <- if (fit.quadratic) 4L else 3L
    if (length(x) < min_n) {
        stop(sprintf("At least %d observations are required for %s model",
                     min_n, if (fit.quadratic) "quadratic" else "linear"),
             call. = FALSE)
    }

    if (!all(y %in% c(0, 1))) {
        stop("'y' must contain only 0 or 1 values", call. = FALSE)
    }
    if (!all(is.finite(x))) {
        stop("'x' contains non-finite values (NA, NaN, or Inf)", call. = FALSE)
    }

    # Parameter validation
    if (!is.logical(fit.quadratic) || length(fit.quadratic) != 1L ||
        is.na(fit.quadratic)) {
        stop("'fit.quadratic' must be a single logical value (TRUE or FALSE)",
             call. = FALSE)
    }

    if (!is.logical(with.errors) || length(with.errors) != 1L ||
        is.na(with.errors)) {
        stop("'with.errors' must be a single logical value (TRUE or FALSE)",
             call. = FALSE)
    }

    max.iterations <- as.integer(max.iterations)
    if (length(max.iterations) != 1L || is.na(max.iterations) ||
        max.iterations <= 0L) {
        stop("'max.iterations' must be a single positive integer", call. = FALSE)
    }

    if (!is.numeric(ridge.lambda) || length(ridge.lambda) != 1L ||
        is.na(ridge.lambda) || ridge.lambda < 0) {
        stop("'ridge.lambda' must be a single non-negative number", call. = FALSE)
    }

    if (!is.numeric(tolerance) || length(tolerance) != 1L ||
        is.na(tolerance) || tolerance <= 0) {
        stop("'tolerance' must be a single positive number", call. = FALSE)
    }

    # Weight validation and processing
    if (is.null(w)) {
        w <- rep(1.0, length(x))
    } else {
        if (!is.numeric(w)) {
            stop("'w' must be a numeric vector", call. = FALSE)
        }
        if (length(w) != length(x)) {
            stop(sprintf("'w' must have the same length as 'x' (w: %d, x: %d)",
                         length(w), length(x)), call. = FALSE)
        }
        if (!all(is.finite(w))) {
            stop("'w' contains non-finite values (NA, NaN, or Inf)", call. = FALSE)
        }
        if (any(w < 0)) {
            stop("All weights must be non-negative", call. = FALSE)
        }
        if (all(w == 0)) {
            stop("At least one weight must be positive", call. = FALSE)
        }
        # Normalize weights
        w <- w / sum(w) * length(w)
    }

    # Call the C implementation
    result <- .Call("S_eigen_ulogit",
                    as.double(x),
                    as.double(y),
                    as.double(w),
                    as.logical(fit.quadratic),
                    as.logical(with.errors),
                    as.integer(max.iterations),
                    as.double(ridge.lambda),
                    as.double(tolerance),
                    PACKAGE = "gflowx")

    # Calculate AIC and BIC
    n <- length(y)
    k <- length(result$beta)  # Number of parameters
    result$aic <- -2 * result$loglik + 2 * k
    result$bic <- -2 * result$loglik + log(n) * k

    # Add additional information
    result$call <- cl
    result$model <- if (fit.quadratic) "quadratic" else "linear"

    # Set class for S3 methods
    class(result) <- "eigen.ulogit"

    return(result)
}


#' Print Method for ulogit Objects
#'
#' @param x An object of class \code{"ulogit"}
#' @param ... Additional arguments (currently ignored)
#'
#' @return Invisibly returns the input object
#'
#' @examples
#' fit <- list(
#'   call = quote(ulogit(x, y)),
#'   converged = TRUE,
#'   iterations = 6L,
#'   errors = c(0.10, 0.08, 0.12),
#'   beta = c(0.5, 1.2)
#' )
#' class(fit) <- "ulogit"
#'
#' fit
#' @export
#' @method print ulogit
print.ulogit <- function(x, ...) {
    cat("\nUnivariate Logistic Regression Model\n")
    cat("Call:\n")
    print(x$call)
    cat("\nConvergence:", ifelse(x$converged, "Yes", "No"))
    cat("\nIterations:", x$iterations)
    cat("\nMean LOOCV error:", round(mean(x$errors), 4))
    cat("\n")
    invisible(x)
}


#' Print Method for eigen.ulogit Objects
#'
#' @param x An object of class \code{"eigen.ulogit"}
#' @param ... Additional arguments (currently ignored)
#'
#' @return Invisibly returns the input object
#'
#' @examples
#' fit <- list(
#'   call = quote(eigen.ulogit(x, y, fit.quadratic = TRUE)),
#'   model = "quadratic",
#'   beta = c(0.5, 1.2, -0.3),
#'   converged = TRUE,
#'   iterations = 5L,
#'   loglik = -12.34,
#'   aic = 30.68,
#'   bic = 34.12,
#'   errors = c(0.09, 0.11, 0.10)
#' )
#' class(fit) <- "eigen.ulogit"
#'
#' fit
#' @export
#' @method print eigen.ulogit
print.eigen.ulogit <- function(x, ...) {
    cat("\nUnivariate Logistic Regression Model (Eigen Implementation)\n")
    cat("Call:\n")
    print(x$call)
    cat("\nModel type:", x$model)
    cat("\nCoefficients:\n")
    coef_names <- c("(Intercept)", "x")
    if (x$model == "quadratic") coef_names <- c(coef_names, "x^2")
    names(x$beta) <- coef_names
    print(round(x$beta, 4))
    cat("\nConvergence:", ifelse(x$converged, "Yes", "No"))
    cat("\nIterations:", x$iterations)
    cat("\nLog-likelihood:", round(x$loglik, 4))
    cat("\nAIC:", round(x$aic, 2))
    cat("\nBIC:", round(x$bic, 2))
    if (!is.null(x$errors)) {
        cat("\nMean LOOCV error:", round(mean(x$errors), 4))
    }
    cat("\n")
    invisible(x)
}


#' Predict Method for eigen.ulogit Objects
#'
#' @description
#' Obtain predictions from a fitted eigen.ulogit model
#'
#' @param object An object of class \code{"eigen.ulogit"}
#' @param newdata Optional numeric vector of new x values for prediction.
#'   If omitted, fitted values are returned.
#' @param type Character string specifying the type of prediction. Either
#'   \code{"response"} for probabilities or \code{"link"} for linear predictors.
#' @param ... Additional arguments (currently ignored)
#'
#' @return Numeric vector of predictions
#'
#' @export
#' @method predict eigen.ulogit
predict.eigen.ulogit <- function(object, newdata = NULL,
                                 type = c("response", "link"), ...) {
    type <- match.arg(type)

    if (is.null(newdata)) {
        if (type == "response") {
            return(object$predictions)
        } else {
            # Calculate linear predictor from probabilities
            p <- object$predictions
            p <- pmax(pmin(p, 1 - 1e-15), 1e-15)  # Avoid log(0) or log(1)
            return(log(p / (1 - p)))
        }
    }

    # Validate newdata
    if (!is.numeric(newdata)) {
        stop("'newdata' must be numeric", call. = FALSE)
    }
    if (!all(is.finite(newdata))) {
        stop("'newdata' contains non-finite values", call. = FALSE)
    }

    # Calculate linear predictor
    eta <- object$beta[1] + object$beta[2] * newdata
    if (object$model == "quadratic") {
        eta <- eta + object$beta[3] * newdata^2
    }

    if (type == "response") {
        return(1 / (1 + exp(-eta)))
    } else {
        return(eta)
    }
}
