#' gflowx: Archived Geometric Smoothing Methods
#'
#' `gflowx` preserves retired graph-regression, response-smoothing, and
#' conditional-expectation estimators formerly shipped by `gflow`. Graph
#' construction is provided by `dgraphs`; the archived numerical estimators are
#' implemented in this package.
#'
#' @keywords internal
#' @useDynLib gflowx, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom graphics abline grid hist legend lines mtext par polygon rect text
#' @importFrom grDevices adjustcolor rgb
#' @importFrom stats cor prcomp quantile sd
#' @importFrom utils head
"_PACKAGE"
