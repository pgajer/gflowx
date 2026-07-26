#' Build iKNN Graphs, Lay Out Components, and Fit rdgraph Regression
#'
#' @description
#' Convenience pipeline that combines:
#' \itemize{
#'   \item \code{dgraphs::build.iknn.graphs.and.selectk()}
#'   \item \code{grip::grip.layout}
#'   \item \code{\link{fit.rdgraph.regression}}
#' }
#' for a single data matrix \code{X} and response \code{y}.
#'
#' The function first builds/selects an iKNN graph sequence, chooses one graph
#' (by \code{selected.k} or internal selection fields), splits it into connected
#' components, computes 3D layouts for each component, and then:
#' \itemize{
#'   \item if the selected graph has multiple components: generates a whole-graph
#'         3D HTML view colored by connected component (default)
#'   \item if the selected graph is connected: fits rdgraph regression on the selected
#'         component and generates continuous 3D HTML for \eqn{\hat y - \bar y}
#' }
#'
#' @param X Numeric matrix-like object (samples in rows, features in columns).
#' @param y Numeric response vector of length \code{nrow(X)} (binary allowed).
#' @param k.min Integer minimum k passed to
#'   \code{dgraphs::build.iknn.graphs.and.selectk()}.
#' @param k.max Integer maximum k passed to
#'   \code{dgraphs::build.iknn.graphs.and.selectk()}.
#' @param build.args Named list passed to
#'   \code{dgraphs::build.iknn.graphs.and.selectk()}.
#'   \code{X}, \code{kmin}, and \code{kmax} are always supplied by this wrapper.
#'   If \code{build.args$kmin}/\code{build.args$kmax} are provided, top-level
#'   \code{k.min}/\code{k.max} take precedence.
#' @param grip.args Named list passed to \code{grip::grip.layout()}.
#'   Defaults are:
#'   \code{dim=3L, rounds=200, final_rounds=200, num_init=10, num_nbrs=30,
#'   r=0.1, s=1.0, tinit_factor=6, seed=6L}.
#' @param fit.args Named list passed to \code{fit.rdgraph.regression()} in the
#'   connected-graph case. \code{X}, \code{y}, and \code{k} are injected by this
#'   wrapper and take precedence.
#' @param plain.plot.args Named list passed to
#'   \code{gflow::plot3D.plain.widget()}.
#' @param cltr.plot.args Named list passed to
#'   \code{gflow::plot3D.cltrs.widget()}
#'   for the multi-component whole-graph plot.
#' @param cont.plot.args Named list passed to
#'   \code{gflow::plot3D.cont.widget()}.
#' @param multi.comp.plot Plot mode used when selected graph has multiple
#'   connected components: \code{"cltrs"} (default) for component-colored whole-graph
#'   plot, \code{"plain"} for whole-graph plain plot, or \code{"both"}.
#' @param selected.k Optional integer. If supplied, selects that k (or nearest
#'   available k). If \code{NULL}, chooses in this order:
#'   \code{k.opt.edit}, \code{k.opt.mixing}, \code{k.cc.edit}, \code{k.cc.mixing},
#'   then first connected k, then largest-LCC k.
#' @param out.dir Optional output directory. If non-NULL, directory is created
#'   and generated objects are saved as timestamped \code{.rds} files. HTML files
#'   from plotting functions are also saved there unless explicit \code{output.file}
#'   is already provided in plot args.
#' @param file.prefix Prefix used in saved filenames.
#' @param timestamp Optional timestamp string used in saved filenames. If NULL,
#'   current time in \code{"\%Y\%m\%d_\%H\%M\%S"} format is used.
#' @param save.rds Logical. If TRUE and \code{out.dir} is set, saves created
#'   objects as timestamped \code{.rds}.
#' @param verbose Logical.
#' @param suppress.rplots Logical. If TRUE, temporarily routes the default plot
#'   device to a null PDF device while this function runs and removes a newly
#'   created \code{Rplots.pdf} in the working directory.
#'
#' @return A list with:
#' \itemize{
#'   \item \code{graphs}: full output of
#'         \code{dgraphs::build.iknn.graphs.and.selectk()}
#'   \item \code{connectivity}: connectivity table
#'   \item \code{k.values}: k sequence
#'   \item \code{selected.k}, \code{selected.k.index}
#'   \item \code{selected.graph}: chosen graph object
#'   \item \code{n.ccomp}: number of connected components in selected graph
#'   \item \code{components}: per-component objects (layout, optional fit/HTML)
#'   \item \code{html.objects}: list of htmlwidget objects produced
#'   \item \code{saved.files}: named vector of saved \code{.rds} files
#'   \item \code{notes}: character vector with run notes
#' }
#'
#' @details
#' In the connected case, this wrapper passes the selected component graph
#' directly to \code{fit.rdgraph.regression()} via \code{adj.list} and
#' \code{weight.list}, avoiding redundant graph reconstruction from \code{X}.
#'
#' @export
iknn.graph.response.pipeline <- function(
    X,
    y,
    k.min = NULL,
    k.max = NULL,
    build.args = list(),
    grip.args = list(
        dim = 3L,
        rounds = 200,
        final_rounds = 200,
        num_init = 10,
        num_nbrs = 30,
        r = 0.1,
        s = 1.0,
        tinit_factor = 6,
        seed = 6L
    ),
    fit.args = list(
        max.iterations = 20L,
        n.eigenpairs = 50L,
        t.scale.factor = 1.0,
        beta.coef.factor = 0.5,
        verbose.level = 1L
    ),
    plain.plot.args = list(),
    cltr.plot.args = list(),
    cont.plot.args = list(),
    multi.comp.plot = c("cltrs", "plain", "both"),
    selected.k = NULL,
    out.dir = NULL,
    file.prefix = "iknn_graph_response_pipeline",
    timestamp = NULL,
    save.rds = TRUE,
    verbose = TRUE,
    suppress.rplots = TRUE
) {
    grip.layout <- .gflowx.get.namespace.export(
        pkg = "grip",
        name = "grip.layout",
        api = "grip.layout()",
        install_hint = "Install it first (e.g. remotes::install_github('pgajer/grip'))."
    )

    get.gflow.fn <- function(name) {
        if (exists(name, mode = "function", inherits = TRUE)) {
            return(get(name, mode = "function", inherits = TRUE))
        }
        if (requireNamespace("gflow", quietly = TRUE) &&
            exists(name, envir = asNamespace("gflow"), inherits = FALSE)) {
            return(get(name, envir = asNamespace("gflow"), inherits = FALSE))
        }
        stop("Function not found: ", name)
    }

    f.build <- dgraphs::build.iknn.graphs.and.selectk
    f.cc <- dgraphs::graph.connected.components
    f.fit <- fit.rdgraph.regression
    f.plot.plain <- get.gflow.fn("plot3D.plain.widget")
    f.plot.cltrs <- get.gflow.fn("plot3D.cltrs.widget")
    f.plot.cont <- get.gflow.fn("plot3D.cont.widget")
    multi.comp.plot <- match.arg(multi.comp.plot)

    if (!is.matrix(X)) {
        X <- tryCatch(as.matrix(X), error = function(e) NULL)
    }
    if (is.null(X) || !is.numeric(X)) stop("X must be numeric or coercible to numeric matrix.")
    if (nrow(X) < 3L || ncol(X) < 1L) stop("X must have at least 3 rows and 1 column.")

    y <- as.numeric(y)
    if (length(y) != nrow(X)) stop("length(y) must equal nrow(X).")
    if (anyNA(y) || any(!is.finite(y))) stop("y cannot contain NA/Inf.")

    if (is.null(k.min) && !is.null(build.args$kmin)) {
        k.min <- build.args$kmin
        if (isTRUE(verbose)) warning("Using deprecated build.args$kmin; prefer top-level k.min.")
    }
    if (is.null(k.max) && !is.null(build.args$kmax)) {
        k.max <- build.args$kmax
        if (isTRUE(verbose)) warning("Using deprecated build.args$kmax; prefer top-level k.max.")
    }

    if (is.null(k.min) || is.null(k.max)) {
        stop("Both k.min and k.max must be provided.")
    }

    k.min <- as.integer(k.min)
    k.max <- as.integer(k.max)
    if (!is.finite(k.min) || !is.finite(k.max) || k.min < 1L || k.max < k.min) {
        stop("Invalid k range: require integers with 1 <= k.min <= k.max.")
    }

    if (is.null(timestamp) || !nzchar(timestamp)) {
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    }
    timestamp <- gsub("[^0-9A-Za-z_\\-]", "_", as.character(timestamp))
    file.prefix <- gsub("[^0-9A-Za-z_\\-]", "_", as.character(file.prefix))

    if (!is.null(out.dir)) {
        out.dir <- path.expand(out.dir)
        dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
    }

    old.wd <- getwd()
    tmp.wd <- NULL
    if (isTRUE(suppress.rplots)) {
        tmp.wd <- file.path(
            tempdir(),
            sprintf("iknn_graph_response_pipeline_%d_%d", Sys.getpid(), as.integer(stats::runif(1) * 1e9))
        )
        dir.create(tmp.wd, recursive = TRUE, showWarnings = FALSE)
        setwd(tmp.wd)
    }
    on.exit({
        if (!is.null(tmp.wd) && dir.exists(tmp.wd)) {
            setwd(old.wd)
            unlink(file.path(tmp.wd, "Rplots.pdf"), force = TRUE)
            unlink(tmp.wd, recursive = TRUE, force = TRUE)
        } else if (getwd() != old.wd) {
            setwd(old.wd)
        }
    }, add = TRUE)

    old.device.opt <- getOption("device")
    rplots.path <- file.path(old.wd, "Rplots.pdf")
    had.rplots.before <- file.exists(rplots.path)
    if (isTRUE(suppress.rplots)) {
        options(device = function(...) grDevices::pdf(file = NULL))
    }
    on.exit({
        if (isTRUE(suppress.rplots)) {
            options(device = old.device.opt)
            if (!had.rplots.before && file.exists(rplots.path)) {
                unlink(rplots.path, force = TRUE)
            }
        }
    }, add = TRUE)

    save.object <- function(object, stem) {
        if (is.null(out.dir) || !isTRUE(save.rds)) return(NA_character_)
        file <- file.path(out.dir, sprintf("%s_%s_%s.rds", file.prefix, timestamp, stem))
        saveRDS(object, file = file)
        file
    }

    `%or%` <- function(a, b) if (!is.null(a)) a else b

    if (is.null(build.args$verbose)) build.args$verbose <- isTRUE(verbose)
    build.args$kmin <- NULL
    build.args$kmax <- NULL
    build.args$k.min <- NULL
    build.args$k.max <- NULL

    graphs <- do.call(f.build, c(list(X = X, kmin = k.min, kmax = k.max), build.args))

    if (is.null(graphs$X.graphs) || is.null(graphs$X.graphs$geom_pruned_graphs)) {
        stop("build.iknn.graphs.and.selectk() did not return X.graphs$geom_pruned_graphs.")
    }
    g.list <- graphs$X.graphs$geom_pruned_graphs
    if (length(g.list) == 0L) stop("No graphs in X.graphs$geom_pruned_graphs.")

    k.values <- as.integer(graphs$k.values %or% seq_along(g.list))
    conn.df <- graphs$connectivity
    if (is.null(conn.df) || !is.data.frame(conn.df) || nrow(conn.df) == 0L) {
        conn.df <- data.frame(
            k = k.values,
            n.components = rep(NA_integer_, length(k.values)),
            stringsAsFactors = FALSE
        )
    }
    if (!("k" %in% names(conn.df))) conn.df$k <- k.values

    choose.k.index <- function() {
        if (!is.null(selected.k) && is.finite(as.numeric(selected.k))) {
            sk <- as.integer(selected.k)
            hit <- which(k.values == sk)
            if (length(hit) > 0L) return(hit[[1L]])
            if (isTRUE(verbose)) {
                warning("selected.k=", sk, " not available. Using nearest k in sequence.")
            }
            return(which.min(abs(k.values - sk)))
        }

        cand <- c(graphs$k.opt.edit, graphs$k.opt.mixing, graphs$k.cc.edit, graphs$k.cc.mixing)
        cand <- as.integer(cand[is.finite(as.numeric(cand))])
        cand <- unique(cand)
        for (kk in cand) {
            hit <- which(k.values == kk)
            if (length(hit) > 0L) return(hit[[1L]])
        }

        ncomp.col <- NULL
        if ("n.components" %in% names(conn.df)) ncomp.col <- "n.components"
        if (is.null(ncomp.col) && "n_ccomp" %in% names(conn.df)) ncomp.col <- "n_ccomp"
        if (!is.null(ncomp.col)) {
            hit <- which(as.integer(conn.df[[ncomp.col]]) == 1L)
            if (length(hit) > 0L) return(hit[[1L]])
        }

        if ("lcc.frac" %in% names(conn.df) && any(is.finite(conn.df$lcc.frac))) {
            return(which.max(conn.df$lcc.frac))
        }

        1L
    }

    selected.idx <- as.integer(choose.k.index())
    selected.k.value <- as.integer(k.values[[selected.idx]])
    selected.graph <- g.list[[selected.idx]]

    if (is.null(selected.graph$adj_list) || !is.list(selected.graph$adj_list)) {
        stop("Selected graph does not contain a valid adj_list.")
    }
    sanitize.weights <- function(w) {
        w <- as.double(w)
        bad <- !is.finite(w) | (w <= 0)
        if (any(bad)) {
            good <- w[is.finite(w) & (w > 0)]
            fill <- if (length(good) > 0L) min(good) else 1.0
            w[bad] <- fill
        }
        w
    }
    sanitize.weight.list <- function(weight.list, adj.list) {
        n <- length(adj.list)
        if (is.null(weight.list) || !is.list(weight.list) || length(weight.list) != n) {
            return(lapply(adj.list, function(nb) rep(1.0, length(nb))))
        }
        out <- vector("list", n)
        for (i in seq_len(n)) {
            nb <- as.integer(adj.list[[i]])
            wi <- weight.list[[i]]
            if (is.null(wi) || length(wi) != length(nb)) {
                out[[i]] <- rep(1.0, length(nb))
            } else {
                out[[i]] <- sanitize.weights(wi)
            }
        }
        out
    }
    selected.graph$weight_list <- sanitize.weight.list(
        weight.list = selected.graph$weight_list,
        adj.list = selected.graph$adj_list
    )

    cc <- f.cc(selected.graph$adj_list)
    component.ids <- sort(unique(as.integer(cc)))
    n.ccomp <- length(component.ids)
    component.sizes <- as.integer(vapply(component.ids, function(cid) sum(cc == cid), integer(1L)))
    names(component.sizes) <- as.character(component.ids)

    extract.component.graph <- function(adj.list, weight.list, vertices) {
        vertices <- sort(unique(as.integer(vertices)))
        map <- seq_along(vertices)
        names(map) <- as.character(vertices)

        adj.sub <- vector("list", length(vertices))
        w.sub <- vector("list", length(vertices))

        for (ii in seq_along(vertices)) {
            v <- vertices[[ii]]
            nb <- as.integer(adj.list[[v]])
            if (length(nb) == 0L) {
                adj.sub[[ii]] <- integer(0)
                w.sub[[ii]] <- numeric(0)
                next
            }

            keep <- nb %in% vertices
            nb.keep <- nb[keep]
            if (length(nb.keep) == 0L) {
                adj.sub[[ii]] <- integer(0)
                w.sub[[ii]] <- numeric(0)
                next
            }

            adj.sub[[ii]] <- as.integer(unname(map[as.character(nb.keep)]))

            if (is.null(weight.list) || length(weight.list) < v || is.null(weight.list[[v]])) {
                w.sub[[ii]] <- rep(1, length(adj.sub[[ii]]))
            } else {
                wv <- as.double(weight.list[[v]])
                if (length(wv) != length(nb)) {
                    w.sub[[ii]] <- rep(1, length(adj.sub[[ii]]))
                } else {
                    w.sub[[ii]] <- sanitize.weights(wv[keep])
                }
            }
        }

        list(
            adj_list = adj.sub,
            weight_list = w.sub,
            vertex.idx = vertices
        )
    }

    html.objects <- list()
    components <- vector("list", n.ccomp)
    saved.files <- c(
        graphs = save.object(graphs, "graphs"),
        connectivity = save.object(conn.df, "connectivity"),
        selected_graph = save.object(selected.graph, sprintf("selected_graph_k%02d", selected.k.value))
    )
    notes <- character(0)
    layout.multi <- NULL

    if (n.ccomp == 1L) {
        notes <- c(
            notes,
            "Connected graph fit used precomputed adj.list/weight.list from selected iKNN graph."
        )
    } else {
        notes <- c(
            notes,
            sprintf(
                "Multi-component graph (n.ccomp=%d; sizes=%s) rendered as whole-graph layout with mode='%s'.",
                n.ccomp,
                paste(component.sizes, collapse = ","),
                multi.comp.plot
            )
        )
    }

    layout.args.base <- utils::modifyList(
        list(
            dim = 3L,
            rounds = 200,
            final_rounds = 200,
            num_init = 10,
            num_nbrs = 30,
            r = 0.1,
            s = 1.0,
            tinit_factor = 6,
            seed = 6L
        ),
        grip.args
    )

    if (n.ccomp > 1L) {
        layout.args.multi <- layout.args.base
        layout.args.multi$adj_list <- selected.graph$adj_list
        layout.args.multi$weight_list <- selected.graph$weight_list
        if (is.null(layout.args.multi$disconnected)) {
            layout.args.multi$disconnected <- "components"
        }

        layout.multi <- do.call(grip.layout, layout.args.multi)
        layout.multi <- as.matrix(layout.multi)
        if (ncol(layout.multi) != 3L || nrow(layout.multi) != length(selected.graph$adj_list)) {
            stop(
                "Whole-graph grip.layout() output must be n x 3 with n=",
                length(selected.graph$adj_list),
                "; got ",
                nrow(layout.multi),
                " x ",
                ncol(layout.multi),
                "."
            )
        }
        saved.files[paste0("layout_multi_k", sprintf("%02d", selected.k.value))] <-
            save.object(layout.multi, paste0("layout_multi_k", sprintf("%02d", selected.k.value)))
    }

    for (ii in seq_along(component.ids)) {
        cid <- component.ids[[ii]]
        cseq <- as.integer(ii)
        vertex.idx <- which(cc == cid)
        g.comp <- extract.component.graph(
            adj.list = selected.graph$adj_list,
            weight.list = selected.graph$weight_list,
            vertices = vertex.idx
        )

        comp.tag <- sprintf("k%02d_comp%02d", selected.k.value, cseq)
        fit.obj <- NULL
        html.obj <- NULL
        html.file <- NA_character_

        if (n.ccomp > 1L) {
            layout.3d <- layout.multi[g.comp$vertex.idx, , drop = FALSE]
        } else {
            layout.args <- layout.args.base
            layout.args$adj_list <- g.comp$adj_list
            layout.args$weight_list <- g.comp$weight_list

            layout.3d <- do.call(grip.layout, layout.args)
            layout.3d <- as.matrix(layout.3d)
            if (ncol(layout.3d) != 3L) {
                stop("grip.layout() must return a 3-column embedding. Got ncol=", ncol(layout.3d))
            }

            x.comp <- X[g.comp$vertex.idx, , drop = FALSE]
            y.comp <- as.double(y[g.comp$vertex.idx])

            fit.args.use <- fit.args
            fit.args.use$X <- x.comp
            fit.args.use$y <- y.comp
            fit.args.use$adj.list <- g.comp$adj_list
            fit.args.use$weight.list <- g.comp$weight_list

            if (!is.null(fit.args.use$verbose) && is.null(fit.args.use$verbose.level)) {
                fit.args.use$verbose.level <- if (isTRUE(fit.args.use$verbose)) 1L else 0L
            }
            fit.args.use$verbose <- NULL

            k.fit <- as.integer(fit.args.use$k %or% selected.k.value)
            if (!is.finite(k.fit)) k.fit <- as.integer(max(2L, min(10L, nrow(x.comp) - 1L)))
            if (k.fit >= nrow(x.comp)) {
                k.fit <- nrow(x.comp) - 1L
                if (isTRUE(verbose)) {
                    warning("Adjusted k in fit.rdgraph.regression to ", k.fit,
                            " for component size n=", nrow(x.comp), ".")
                }
            }
            if (k.fit < 2L) {
                stop("Cannot fit rdgraph on component with < 3 vertices.")
            }
            fit.args.use$k <- k.fit

            fit.obj <- do.call(f.fit, fit.args.use)
            saved.files[paste0("fit_", comp.tag)] <- save.object(fit.obj, paste0("fit_", comp.tag))

            y.plot <- as.double(fit.obj$fitted.values) - mean(y.comp)
            cont.args <- utils::modifyList(
                list(X = layout.3d, y = y.plot, legend.title = "fitted(y)-mean(y)"),
                cont.plot.args
            )
            if (is.null(cont.args$output.file) && !is.null(out.dir)) {
                html.file <- file.path(
                    out.dir,
                    sprintf("%s_%s_%s_cont.html", file.prefix, timestamp, comp.tag)
                )
                cont.args$output.file <- html.file
            } else {
                html.file <- as.character(cont.args$output.file %or% NA_character_)
            }

            html.obj <- do.call(f.plot.cont, cont.args)
        }

        saved.files[paste0("layout_", comp.tag)] <- save.object(layout.3d, paste0("layout_", comp.tag))
        if (!is.null(html.obj)) {
            html.objects[[length(html.objects) + 1L]] <- html.obj
            names(html.objects)[[length(html.objects)]] <- comp.tag
            saved.files[paste0("html_widget_", comp.tag)] <- save.object(html.obj, paste0("html_widget_", comp.tag))
        }

        components[[ii]] <- list(
            component.seq = cseq,
            component.id = as.integer(cid),
            n.vertices = length(g.comp$vertex.idx),
            vertex.idx = g.comp$vertex.idx,
            graph = g.comp,
            layout.3d = layout.3d,
            fit = fit.obj,
            html.object = html.obj,
            html.file = html.file
        )
    }

    if (n.ccomp > 1L) {
        comp.labels <- sprintf(
            "comp%02d (n=%d)",
            seq_along(component.ids),
            component.sizes
        )
        cltr.vec <- factor(
            cc,
            levels = component.ids,
            labels = comp.labels
        )

        if (multi.comp.plot %in% c("cltrs", "both")) {
            cltr.args <- utils::modifyList(
                list(
                    X = layout.multi,
                    cltr = cltr.vec,
                    legend.title = "Connected Component"
                ),
                cltr.plot.args
            )
            cltr.tag <- sprintf("k%02d_components_cltrs", selected.k.value)
            cltr.file <- as.character(cltr.args$output.file %or% NA_character_)
            if (is.null(cltr.args$output.file) && !is.null(out.dir)) {
                cltr.file <- file.path(
                    out.dir,
                    sprintf("%s_%s_%s.html", file.prefix, timestamp, cltr.tag)
                )
                cltr.args$output.file <- cltr.file
            }
            html.cltr <- do.call(f.plot.cltrs, cltr.args)
            html.objects[[length(html.objects) + 1L]] <- html.cltr
            names(html.objects)[[length(html.objects)]] <- cltr.tag
            saved.files[paste0("html_widget_", cltr.tag)] <- save.object(html.cltr, paste0("html_widget_", cltr.tag))
            saved.files[paste0("html_file_", cltr.tag)] <- cltr.file
        }

        if (multi.comp.plot %in% c("plain", "both")) {
            plain.args <- utils::modifyList(list(X = layout.multi), plain.plot.args)
            plain.tag <- sprintf("k%02d_components_plain", selected.k.value)
            plain.file <- as.character(plain.args$output.file %or% NA_character_)
            if (is.null(plain.args$output.file) && !is.null(out.dir)) {
                plain.file <- file.path(
                    out.dir,
                    sprintf("%s_%s_%s.html", file.prefix, timestamp, plain.tag)
                )
                plain.args$output.file <- plain.file
            }
            html.plain <- do.call(f.plot.plain, plain.args)
            html.objects[[length(html.objects) + 1L]] <- html.plain
            names(html.objects)[[length(html.objects)]] <- plain.tag
            saved.files[paste0("html_widget_", plain.tag)] <- save.object(html.plain, paste0("html_widget_", plain.tag))
            saved.files[paste0("html_file_", plain.tag)] <- plain.file
        }
    }

    result <- list(
        graphs = graphs,
        connectivity = conn.df,
        k.values = as.integer(k.values),
        selected.k = as.integer(selected.k.value),
        selected.k.index = as.integer(selected.idx),
        selected.graph = selected.graph,
        n.ccomp = as.integer(n.ccomp),
        component.sizes = component.sizes,
        layout.multi = layout.multi,
        multi.comp.plot = multi.comp.plot,
        components = components,
        html.objects = html.objects,
        saved.files = saved.files,
        notes = notes
    )

    saved.files["result"] <- save.object(result, "result")
    result$saved.files <- saved.files
    result
}
