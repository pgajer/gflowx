.gflowx_gflow_function <- function(name) {
    if (!requireNamespace("gflow", quietly = TRUE)) {
        stop("gflow is required for the legacy rdgraph backend.", call. = FALSE)
    }

    ns <- asNamespace("gflow")
    if (exists(name, envir = ns, inherits = FALSE)) {
        return(get(name, envir = ns, inherits = FALSE))
    }

    getExportedValue("gflow", name)
}

#' Bridge to `gflow::fit.rdgraph.regression()`
#'
#' These functions are the initial `gflowx` landing surface for legacy rdgraph
#' regression workflows. They currently delegate to the legacy backend in
#' `gflow`; this keeps existing behavior fixed while the native rdgraph engine
#' and its graph utilities are extracted in smaller audited slices.
#'
#' @param ... Arguments forwarded unchanged to the corresponding `gflow`
#'   function.
#'
#' @return The object returned by the delegated `gflow` implementation.
#'
#' @name rdgraph-bridges
NULL

#' @rdname rdgraph-bridges
#' @export
fit.rdgraph.regression <- function(...) {
    .gflowx_gflow_function("fit.rdgraph.regression")(...)
}

#' @rdname rdgraph-bridges
#' @export
refit.rdgraph.regression <- function(...) {
    .gflowx_gflow_function("refit.rdgraph.regression")(...)
}

#' @rdname rdgraph-bridges
#' @export
permutation.test.rdgraph <- function(...) {
    .gflowx_gflow_function("permutation.test.rdgraph")(...)
}

#' @rdname rdgraph-bridges
#' @export
perm.test.audit <- function(...) {
    .gflowx_gflow_function("perm.test.audit")(...)
}

#' @rdname rdgraph-bridges
#' @export
bayes.bootstrap.rdgraph <- function(...) {
    .gflowx_gflow_function("bayes.bootstrap.rdgraph")(...)
}

#' @rdname rdgraph-bridges
#' @export
lcor.with.posterior <- function(...) {
    .gflowx_gflow_function("lcor.with.posterior")(...)
}

#' @rdname rdgraph-bridges
#' @export
extremality.summary <- function(...) {
    .gflowx_gflow_function("extremality.summary")(...)
}

#' @rdname rdgraph-bridges
#' @export
label.extremality.3d <- function(...) {
    .gflowx_gflow_function("label.extremality.3d")(...)
}

#' @rdname rdgraph-bridges
#' @export
compute.pextrema.nbhds <- function(...) {
    .gflowx_gflow_function("compute.pextrema.nbhds")(...)
}

#' @rdname rdgraph-bridges
#' @export
compute.cluster.summary <- function(...) {
    .gflowx_gflow_function("compute.cluster.summary")(...)
}

#' @rdname rdgraph-bridges
#' @export
extract.cluster.representatives <- function(...) {
    .gflowx_gflow_function("extract.cluster.representatives")(...)
}

#' @rdname rdgraph-bridges
#' @export
compute.dbscan.cluster.summary <- function(...) {
    .gflowx_gflow_function("compute.dbscan.cluster.summary")(...)
}

#' @rdname rdgraph-bridges
#' @export
extract.dbscan.cluster.representatives <- function(...) {
    .gflowx_gflow_function("extract.dbscan.cluster.representatives")(...)
}

#' @rdname rdgraph-bridges
#' @export
analyze.dbscan.noise <- function(...) {
    .gflowx_gflow_function("analyze.dbscan.noise")(...)
}

#' @rdname rdgraph-bridges
#' @export
subject.neighborhood.stats <- function(...) {
    .gflowx_gflow_function("subject.neighborhood.stats")(...)
}

#' @rdname rdgraph-bridges
#' @export
rdgraph.neighbor.weights <- function(...) {
    .gflowx_gflow_function("rdgraph.neighbor.weights")(...)
}
