.gflowx.get.namespace.export <- function(pkg, name, api = NULL,
                                         install_hint = NULL) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        label <- if (is.null(api)) name else api
        hint <- if (is.null(install_hint)) "" else paste0(" ", install_hint)
        stop(label, " requires optional package '", pkg, "'.", hint,
             call. = FALSE)
    }
    ns <- asNamespace(pkg)
    if (!exists(name, envir = ns, inherits = FALSE)) {
        stop("Optional package '", pkg, "' does not provide ", name, "().",
             call. = FALSE)
    }
    get(name, envir = ns, inherits = FALSE)
}
