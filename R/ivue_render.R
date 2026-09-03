.gflowx.render.ivue <- function(kind, args) {
    if (!requireNamespace("ivue", quietly = TRUE) ||
        utils::packageVersion("ivue") < "0.0.0.9001") {
        stop("Install ivue >= 0.0.0.9001 for pipeline plotting.", call. = FALSE)
    }
    file <- args$output.file
    selfcontained <- if (is.null(args$selfcontained)) TRUE else args$selfcontained
    open.browser <- isTRUE(args$open.browser)
    args[c("output.file", "selfcontained", "open.browser")] <- NULL
    fun <- switch(kind, plain = ivue::plot3D.plain,
                  cont = ivue::plot3D.cont, cltrs = ivue::plot3D.cltrs,
                  stop("Unknown pipeline plot kind.", call. = FALSE))
    if (identical(kind, "cont") && is.null(args$scale)) {
        args$scale <- ivue::color.scale.cont(args$values, mode = "binned",
            winsor.p = 0.01, palette = function(n) grDevices::rainbow(n, start = 1/6, end = 0))
    }
    widget <- do.call(fun, args)
    if (!is.null(file)) {
        if (!is.character(file) || length(file) != 1L || is.na(file) || !nzchar(file))
            stop("output.file must be one nonempty path.", call. = FALSE)
        file <- path.expand(file)
        dir.create(dirname(file), recursive = TRUE, showWarnings = FALSE)
        htmlwidgets::saveWidget(widget, file, selfcontained = selfcontained)
        if (open.browser) utils::browseURL(file)
    }
    widget
}
