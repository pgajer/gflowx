.gflowx_gflow_export <- function(name) {
    getExportedValue("gflow", name)
}

#' Bridge to `gflow::fit.rdgraph.regression()`
#'
#' These functions are the initial `gflowx` landing surface for legacy rdgraph
#' regression workflows. They currently delegate to the implementations in
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
    .gflowx_gflow_export("fit.rdgraph.regression")(...)
}

#' @rdname rdgraph-bridges
#' @export
refit.rdgraph.regression <- function(...) {
    .gflowx_gflow_export("refit.rdgraph.regression")(...)
}

#' @rdname rdgraph-bridges
#' @export
permutation.test.rdgraph <- function(...) {
    .gflowx_gflow_export("permutation.test.rdgraph")(...)
}

#' @rdname rdgraph-bridges
#' @export
perm.test.audit <- function(...) {
    .gflowx_gflow_export("perm.test.audit")(...)
}

#' @rdname rdgraph-bridges
#' @export
bayes.bootstrap.rdgraph <- function(...) {
    .gflowx_gflow_export("bayes.bootstrap.rdgraph")(...)
}

#' @rdname rdgraph-bridges
#' @export
lcor.with.posterior <- function(...) {
    .gflowx_gflow_export("lcor.with.posterior")(...)
}

#' @rdname rdgraph-bridges
#' @export
extremality.summary <- function(...) {
    .gflowx_gflow_export("extremality.summary")(...)
}

#' @rdname rdgraph-bridges
#' @export
label.extremality.3d <- function(...) {
    .gflowx_gflow_export("label.extremality.3d")(...)
}

#' @rdname rdgraph-bridges
#' @export
compute.pextrema.nbhds <- function(...) {
    .gflowx_gflow_export("compute.pextrema.nbhds")(...)
}

#' @rdname rdgraph-bridges
#' @export
compute.cluster.summary <- function(...) {
    .gflowx_gflow_export("compute.cluster.summary")(...)
}

#' @rdname rdgraph-bridges
#' @export
extract.cluster.representatives <- function(...) {
    .gflowx_gflow_export("extract.cluster.representatives")(...)
}

#' @rdname rdgraph-bridges
#' @export
compute.dbscan.cluster.summary <- function(...) {
    .gflowx_gflow_export("compute.dbscan.cluster.summary")(...)
}

#' @rdname rdgraph-bridges
#' @export
extract.dbscan.cluster.representatives <- function(...) {
    .gflowx_gflow_export("extract.dbscan.cluster.representatives")(...)
}

#' @rdname rdgraph-bridges
#' @export
analyze.dbscan.noise <- function(...) {
    .gflowx_gflow_export("analyze.dbscan.noise")(...)
}

#' @rdname rdgraph-bridges
#' @export
subject.neighborhood.stats <- function(...) {
    .gflowx_gflow_export("subject.neighborhood.stats")(...)
}

#' @rdname rdgraph-bridges
#' @export
rdgraph.neighbor.weights <- function(...) {
    .gflowx_gflow_export("rdgraph.neighbor.weights")(...)
}
