#' Compute Subject Concentration in Local rdgraph Neighborhoods
#'
#' @description
#' Computes per-vertex diagnostics for repeated-measures concentration in
#' neighborhoods of a graph returned by \code{gflowx::fit.rdgraph.regression()}.
#'
#' \itemize{
#'   \item \eqn{R(v) = |\hat N(v)| / \#\{\text{subjects in }\hat N(v)\}}
#'   with \eqn{\hat N(v) = N(v) \cup \{v\}} (optional self-inclusion).
#'   \item \code{p.max.conductance(v)}, \code{s.eff.conductance(v)}: subject
#'   concentration in \eqn{N(v)} using conductance weights
#'   \eqn{w_{ij} = 1/\rho_1(i,j)}.
#'   \item \code{p.max.operator(v)}, \code{s.eff.operator(v)}: subject
#'   concentration in \eqn{N(v)} using operator weights derived from
#'   \code{L0.mass.sym} (off-diagonals of
#'   \eqn{L_0^{\text{mass.sym}} = M^{-1/2} L_0 M^{-1/2}}) when available;
#'   otherwise from the mass-sym reconstruction
#'   \eqn{\{1/\rho_1(i,j)\}/\sqrt{\rho_0(i)\rho_0(j)}}.
#' }
#'
#' For backward compatibility, the legacy columns \code{p.max} and \code{s.eff}
#' are set according to \code{weight.type} (either conductance or operator).
#'
#' @param fitted.model Fitted object from \code{gflowx::fit.rdgraph.regression()} (class
#'   \code{"knn.riem.fit"}) or a list containing \code{$optimal.fit}.
#' @param subject.id Vector of subject IDs of length \eqn{n} (one ID per vertex).
#' @param weight.type Character scalar controlling the legacy columns
#'   \code{n.subjects.weighted}, \code{weight.sum}, \code{p.max}, \code{s.eff}.
#'   Choices: \code{"conductance"} (default) or \code{"mass.sym"}.
#' @param include.self.in.R Logical scalar; if \code{TRUE} (default), include
#'   vertex \eqn{v} in \eqn{\hat N(v)} for \code{R(v)}.
#' @param edge.mass.floor Positive numeric floor applied to edge density before
#'   inversion. Default \code{1e-10}.
#' @param vertex.mass.floor Positive numeric floor for vertex densities used in
#'   the mass-sym reconstruction. Default \code{1e-15}.
#' @param L0.mass.sym Optional \eqn{n\times n} matrix (dense or sparse) giving
#'   \code{L0_mass_sym} used by the spectral smoother. If \code{NULL} (default),
#'   operator weights are approximated via the mass-sym reconstruction from
#'   \code{fitted.model$graph} densities.
#' @param operator.offdiag.tol Non-negative tolerance for flagging positive
#'   off-diagonals in \code{L0.mass.sym}. Used only for diagnostics.
#'
#' @return Data frame with one row per vertex and columns:
#' \itemize{
#'   \item \code{vertex}: vertex index (1-based)
#'   \item \code{n.neighbors}: \eqn{|N(v)|} based on \code{fitted.model$graph$adj.list}
#'   \item \code{n.hat}: \eqn{|\hat N(v)|} used in \code{R(v)}
#'   \item \code{n.subjects.hat}: number of distinct subjects in \eqn{\hat N(v)}
#'   \item \code{R}: multiplicity ratio
#'   \item \code{n.subjects.conductance}, \code{weight.sum.conductance},
#'   \code{p.max.conductance}, \code{s.eff.conductance}: conductance-based stats
#'   \item \code{n.subjects.operator}, \code{weight.sum.operator},
#'   \code{p.max.operator}, \code{s.eff.operator}: operator-based stats
#'   \item \code{n.subjects.weighted}, \code{weight.sum}, \code{p.max}, \code{s.eff}:
#'   legacy columns controlled by \code{weight.type}
#' }
#'
#' @export
subject.neighborhood.stats <- function(
    fitted.model,
    subject.id,
    weight.type = c("conductance", "mass.sym"),
    include.self.in.R = TRUE,
    edge.mass.floor = 1e-10,
    vertex.mass.floor = 1e-15,
    L0.mass.sym = NULL,
    operator.offdiag.tol = 1e-12
) {
    weight.type <- match.arg(weight.type)

    if (!is.logical(include.self.in.R) || length(include.self.in.R) != 1L || is.na(include.self.in.R)) {
        stop("include.self.in.R must be a single non-NA logical value.")
    }
    if (!is.numeric(edge.mass.floor) || length(edge.mass.floor) != 1L ||
        is.na(edge.mass.floor) || !is.finite(edge.mass.floor) || edge.mass.floor <= 0) {
        stop("edge.mass.floor must be a single finite positive number.")
    }
    if (!is.numeric(vertex.mass.floor) || length(vertex.mass.floor) != 1L ||
        is.na(vertex.mass.floor) || !is.finite(vertex.mass.floor) || vertex.mass.floor <= 0) {
        stop("vertex.mass.floor must be a single finite positive number.")
    }
    if (!is.numeric(operator.offdiag.tol) || length(operator.offdiag.tol) != 1L ||
        is.na(operator.offdiag.tol) || !is.finite(operator.offdiag.tol) || operator.offdiag.tol < 0) {
        stop("operator.offdiag.tol must be a single finite non-negative number.")
    }

    ## Always compute conductance-based neighbor weights
    wobj.conductance <- rdgraph.neighbor.weights(
        fitted.model = fitted.model,
        weight.type = "conductance",
        edge.mass.floor = edge.mass.floor,
        vertex.mass.floor = vertex.mass.floor,
        return.weights.df = FALSE
    )

    adj.list <- wobj.conductance$adj.list
    weight.list.conductance <- wobj.conductance$weight.list
    n <- length(adj.list)

    if (length(subject.id) != n) {
        stop(sprintf(
            "Length of subject.id (%d) must match number of vertices (%d).",
            length(subject.id), n
        ))
    }
    if (anyNA(subject.id)) {
        stop("subject.id cannot contain NA values.")
    }

    ## Map subjects to 1..S indices (stable)
    subject.index <- as.integer(factor(subject.id))

    ## Resolve L0_mass_sym if available (exact operator weights)
    fit <- .resolve.rdgraph_fit_for_subject_stats(fitted.model)
    L0.mass.sym.use <- L0.mass.sym
    if (is.null(L0.mass.sym.use)) {
        if (!is.null(fit$spectral) && is.list(fit$spectral)) {
            if (!is.null(fit$spectral$L0_mass_sym)) {
                L0.mass.sym.use <- fit$spectral$L0_mass_sym
            } else if (!is.null(fit$spectral$L0.mass.sym)) {
                L0.mass.sym.use <- fit$spectral$L0.mass.sym
            }
        }
    }

    operator.weight.source <- "mass.sym.reconstructed"
    if (!is.null(L0.mass.sym.use)) {
        d <- dim(L0.mass.sym.use)
        if (length(d) != 2L || d[1L] != n || d[2L] != n) {
            stop(sprintf(
                "L0.mass.sym must be an %dx%d matrix (got %s).",
                n, n, paste(d, collapse = "x")
            ))
        }
        operator.weight.source <- "L0.mass.sym"
    }

    ## Fallback operator weights: mass.sym reconstruction from densities
    ## (used when L0.mass.sym is not provided / not stored in fit)
    weight.list.mass.sym <- NULL
    if (is.null(L0.mass.sym.use)) {
        weight.list.mass.sym <- tryCatch(
            rdgraph.neighbor.weights(
                fitted.model = fitted.model,
                weight.type = "mass.sym",
                edge.mass.floor = edge.mass.floor,
                vertex.mass.floor = vertex.mass.floor,
                return.weights.df = FALSE
            )$weight.list,
            error = function(e) {
                warning(
                    "Could not reconstruct mass.sym operator weights from fitted.model; ",
                    "p.max.operator and s.eff.operator will be NA. Reason: ", conditionMessage(e),
                    call. = FALSE
                )
                NULL
            }
        )
    }

    ## Small helper: subject concentration stats for a given (neighbors, weights)
    .subject.weight.stats <- function(neighbors, weights) {
        if (length(neighbors) == 0L) {
            return(list(n.subjects = 0L, weight.sum = 0.0, p.max = NA_real_, s.eff = NA_real_))
        }
        if (length(weights) != length(neighbors)) {
            stop("Internal error: neighbors and weights length mismatch.")
        }

        subj.weights <- tapply(weights, subject.index[neighbors], sum)
        total.w <- sum(subj.weights)

        n.subjects <- length(subj.weights)

        if (!is.finite(total.w) || total.w <= 0.0) {
            return(list(n.subjects = n.subjects, weight.sum = total.w, p.max = NA_real_, s.eff = NA_real_))
        }

        p <- as.numeric(subj.weights / total.w)

        list(
            n.subjects = n.subjects,
            weight.sum = total.w,
            p.max = max(p),
            s.eff = 1.0 / sum(p^2)
        )
    }

    n.neighbors <- lengths(adj.list)
    n.hat <- n.neighbors + if (include.self.in.R) 1L else 0L
    n.subjects.hat <- integer(n)
    ratio.R <- rep(NA_real_, n)

    n.subjects.conductance <- integer(n)
    weight.sum.conductance <- numeric(n)
    p.max.conductance <- rep(NA_real_, n)
    s.eff.conductance <- rep(NA_real_, n)

    n.subjects.operator <- integer(n)
    weight.sum.operator <- numeric(n)
    p.max.operator <- rep(NA_real_, n)
    s.eff.operator <- rep(NA_real_, n)

    max.pos.offdiag.operator <- NA_real_
    if (!is.null(L0.mass.sym.use)) {
        max.pos.offdiag.operator <- -Inf
    }

    for (i in seq_len(n)) {
        nhat <- if (include.self.in.R) c(i, adj.list[[i]]) else adj.list[[i]]
        n.subj.hat.i <- length(unique(subject.index[nhat]))
        n.subjects.hat[i] <- n.subj.hat.i
        if (n.subj.hat.i > 0L) {
            ratio.R[i] <- length(nhat) / n.subj.hat.i
        }

        if (n.neighbors[i] == 0L) {
            next
        }

        neighbors <- as.integer(adj.list[[i]])

        ## Conductance weights: w_ij = 1 / rho_1(i,j)
        weights.c <- as.numeric(weight.list.conductance[[i]])
        stats.c <- .subject.weight.stats(neighbors, weights.c)
        n.subjects.conductance[i] <- stats.c$n.subjects
        weight.sum.conductance[i] <- stats.c$weight.sum
        p.max.conductance[i] <- stats.c$p.max
        s.eff.conductance[i] <- stats.c$s.eff

        ## Operator weights: derived from L0_mass_sym off-diagonals when available,
        ## otherwise from mass.sym reconstruction.
        weights.o <- NULL
        if (!is.null(L0.mass.sym.use)) {
            x <- as.numeric(L0.mass.sym.use[i, neighbors])
            if (any(is.finite(x) & x > operator.offdiag.tol)) {
                max.pos.offdiag.operator <- max(max.pos.offdiag.operator, max(x[is.finite(x)], na.rm = TRUE))
            }
            weights.o <- pmax(0.0, -x)
        } else if (!is.null(weight.list.mass.sym)) {
            weights.o <- as.numeric(weight.list.mass.sym[[i]])
        }

        if (!is.null(weights.o)) {
            stats.o <- .subject.weight.stats(neighbors, weights.o)
            n.subjects.operator[i] <- stats.o$n.subjects
            weight.sum.operator[i] <- stats.o$weight.sum
            p.max.operator[i] <- stats.o$p.max
            s.eff.operator[i] <- stats.o$s.eff
        }
    }

    ## Legacy columns preserved for backward compatibility, controlled by weight.type
    if (weight.type == "conductance") {
        n.subjects.weighted <- n.subjects.conductance
        weight.sum <- weight.sum.conductance
        p.max <- p.max.conductance
        s.eff <- s.eff.conductance
    } else {
        n.subjects.weighted <- n.subjects.operator
        weight.sum <- weight.sum.operator
        p.max <- p.max.operator
        s.eff <- s.eff.operator
    }

    out <- data.frame(
        vertex = seq_len(n),
        n.neighbors = n.neighbors,
        n.hat = n.hat,
        n.subjects.hat = n.subjects.hat,
        R = ratio.R,
        n.subjects.weighted = n.subjects.weighted,
        weight.sum = weight.sum,
        p.max = p.max,
        s.eff = s.eff,
        n.subjects.conductance = n.subjects.conductance,
        weight.sum.conductance = weight.sum.conductance,
        p.max.conductance = p.max.conductance,
        s.eff.conductance = s.eff.conductance,
        n.subjects.operator = n.subjects.operator,
        weight.sum.operator = weight.sum.operator,
        p.max.operator = p.max.operator,
        s.eff.operator = s.eff.operator,
        stringsAsFactors = FALSE
    )

    attr(out, "weight.type") <- weight.type
    attr(out, "include.self.in.R") <- include.self.in.R
    attr(out, "operator.weight.source") <- operator.weight.source
    attr(out, "operator.offdiag.tol") <- operator.offdiag.tol
    if (!is.null(L0.mass.sym.use)) {
        attr(out, "operator.max.pos.offdiag") <- if (is.finite(max.pos.offdiag.operator)) max.pos.offdiag.operator else NA_real_
    }
    out
}

#' Reconstruct Per-Vertex Laplacian Neighbor Weights
#'
#' @description
#' Returns explicit per-vertex neighbor weights \eqn{w_{vj}} parallel to
#' \code{fitted.model$graph$adj.list}, reconstructed from the same graph quantities
#' used by \code{gflowx::fit.rdgraph.regression()}.
#'
#' This helper is intended for direct inspection/debugging of local mixing weights.
#'
#' @param fitted.model Fitted object from \code{gflowx::fit.rdgraph.regression()} (class
#'   \code{"knn.riem.fit"}) or a list containing \code{$optimal.fit}.
#' @param weight.type Character scalar: \code{"conductance"} (default) or
#'   \code{"mass.sym"}.
#' @param edge.mass.floor Positive numeric floor applied to edge density before
#'   inversion. Default \code{1e-10}.
#' @param vertex.mass.floor Positive numeric floor for vertex densities when
#'   \code{weight.type = "mass.sym"}. Default \code{1e-15}.
#' @param return.weights.df Logical scalar; if \code{TRUE} (default), also return
#'   a long-format edge-direction table with columns \code{vertex}, \code{neighbor},
#'   and \code{weight}.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{adj.list}: neighbor indices by vertex
#'   \item \code{weight.list}: numeric weights parallel to \code{adj.list}
#'   \item \code{edge.list}: undirected edge endpoints from fit object
#'   \item \code{edge.weights}: one weight per row of \code{edge.list}
#'   \item \code{weights.df}: optional long-format directed table (if requested)
#' }
#'
rdgraph.neighbor.weights <- function(
    fitted.model,
    weight.type = c("conductance", "mass.sym"),
    edge.mass.floor = 1e-10,
    vertex.mass.floor = 1e-15,
    return.weights.df = TRUE
) {
    weight.type <- match.arg(weight.type)

    if (!is.numeric(edge.mass.floor) || length(edge.mass.floor) != 1L ||
        is.na(edge.mass.floor) || !is.finite(edge.mass.floor) || edge.mass.floor <= 0) {
        stop("edge.mass.floor must be a single finite positive number.")
    }
    if (!is.numeric(vertex.mass.floor) || length(vertex.mass.floor) != 1L ||
        is.na(vertex.mass.floor) || !is.finite(vertex.mass.floor) || vertex.mass.floor <= 0) {
        stop("vertex.mass.floor must be a single finite positive number.")
    }
    if (!is.logical(return.weights.df) || length(return.weights.df) != 1L || is.na(return.weights.df)) {
        stop("return.weights.df must be a single non-NA logical value.")
    }

    wdata <- .build.rdgraph_neighbor_weights(
        fitted.model = fitted.model,
        weight.type = weight.type,
        edge.mass.floor = edge.mass.floor,
        vertex.mass.floor = vertex.mass.floor
    )

    out <- list(
        adj.list = wdata$adj.list,
        weight.list = wdata$weight.list,
        edge.list = wdata$edge.list,
        edge.weights = wdata$edge.weights
    )

    if (return.weights.df) {
        n.rows <- sum(lengths(wdata$adj.list))
        v.col <- integer(n.rows)
        n.col <- integer(n.rows)
        w.col <- numeric(n.rows)
        pos <- 1L

        for (i in seq_along(wdata$adj.list)) {
            neighbors <- wdata$adj.list[[i]]
            weights <- wdata$weight.list[[i]]
            k <- length(neighbors)
            if (k == 0L) next
            idx <- pos:(pos + k - 1L)
            v.col[idx] <- i
            n.col[idx] <- neighbors
            w.col[idx] <- weights
            pos <- pos + k
        }

        out$weights.df <- data.frame(
            vertex = v.col,
            neighbor = n.col,
            weight = w.col,
            stringsAsFactors = FALSE
        )
    }

    attr(out, "weight.type") <- weight.type
    attr(out, "edge.mass.floor") <- edge.mass.floor
    attr(out, "vertex.mass.floor") <- vertex.mass.floor

    out
}

# gflowx::fit.rdgraph.regression() results can appear directly or nested in $optimal.fit
# in downstream wrappers.
.resolve.rdgraph_fit_for_subject_stats <- function(x) {
    if (is.list(x) && !is.null(x$graph)) {
        return(x)
    }
    if (is.list(x) && !is.null(x$optimal.fit) && is.list(x$optimal.fit)) {
        if (!is.null(x$optimal.fit$graph)) {
            return(x$optimal.fit)
        }
    }
    stop(
        "fitted.model must be a gflowx::fit.rdgraph.regression() result ",
        "or a list containing $optimal.fit with a $graph component."
    )
}

.build.rdgraph_neighbor_weights <- function(
    fitted.model,
    weight.type,
    edge.mass.floor,
    vertex.mass.floor
) {
    fit <- .resolve.rdgraph_fit_for_subject_stats(fitted.model)
    graph <- fit$graph

    required.graph.fields <- c("adj.list", "edge.list", "edge.densities")
    missing.fields <- setdiff(required.graph.fields, names(graph))
    if (length(missing.fields) > 0L) {
        stop(sprintf(
            "fitted.model$graph is missing required fields: %s",
            paste(missing.fields, collapse = ", ")
        ))
    }
    if (weight.type == "mass.sym" && !("vertex.densities" %in% names(graph))) {
        stop("weight.type = 'mass.sym' requires fitted.model$graph$vertex.densities.")
    }

    adj.list <- graph$adj.list
    if (!is.list(adj.list)) {
        stop("fitted.model$graph$adj.list must be a list of integer neighbor vectors.")
    }
    adj.list <- lapply(adj.list, as.integer)
    n <- length(adj.list)
    if (n < 1L) {
        stop("fitted.model$graph$adj.list is empty.")
    }

    edge.densities <- as.numeric(graph$edge.densities)
    if (any(!is.finite(edge.densities))) {
        stop("fitted.model$graph$edge.densities contains non-finite values.")
    }

    vertex.densities <- NULL
    if (weight.type == "mass.sym") {
        vertex.densities <- as.numeric(graph$vertex.densities)
        if (length(vertex.densities) != n) {
            stop(sprintf(
                "Length of vertex.densities (%d) must match number of vertices (%d).",
                length(vertex.densities), n
            ))
        }
        if (any(!is.finite(vertex.densities))) {
            stop("fitted.model$graph$vertex.densities contains non-finite values.")
        }
    }

    edge.list <- graph$edge.list
    if (!is.matrix(edge.list)) {
        stop("fitted.model$graph$edge.list must be an integer matrix.")
    }
    if (ncol(edge.list) != 2L) {
        if (nrow(edge.list) == 2L) {
            edge.list <- t(edge.list)
        } else {
            stop("fitted.model$graph$edge.list must have 2 columns (edge endpoints).")
        }
    }
    storage.mode(edge.list) <- "integer"

    # Backward-compatibility fallback:
    # Some fit objects contain edge.densities but an empty edge.list.
    # Reconstruct in canonical order by traversing adj.list rows and keeping j > i.
    if (nrow(edge.list) == 0L && length(edge.densities) > 0L) {
        edge.i <- integer(0)
        edge.j <- integer(0)
        for (i in seq_len(n)) {
            neighbors <- adj.list[[i]]
            if (length(neighbors) == 0L) next
            keep <- neighbors > i
            if (any(keep)) {
                edge.i <- c(edge.i, rep.int(i, sum(keep)))
                edge.j <- c(edge.j, neighbors[keep])
            }
        }
        edge.list <- cbind(edge.i, edge.j)
        storage.mode(edge.list) <- "integer"
    }

    if (length(edge.densities) != nrow(edge.list)) {
        stop(sprintf(
            "Length mismatch: edge.densities has length %d, edge.list has %d rows.",
            length(edge.densities), nrow(edge.list)
        ))
    }

    neighbor.position <- vector("list", n)
    for (i in seq_len(n)) {
        neighbors <- adj.list[[i]]
        if (length(neighbors) == 0L) {
            neighbor.position[[i]] <- integer(0)
            next
        }
        if (any(neighbors < 1L | neighbors > n)) {
            stop(sprintf("adj.list[[%d]] contains out-of-range vertex indices.", i))
        }
        if (any(neighbors == i)) {
            stop(sprintf("adj.list[[%d]] should not contain self-loop index %d.", i, i))
        }
        if (anyDuplicated(neighbors)) {
            stop(sprintf("adj.list[[%d]] contains duplicated neighbor indices.", i))
        }
        pos <- seq_along(neighbors)
        names(pos) <- as.character(neighbors)
        neighbor.position[[i]] <- pos
    }

    weight.list <- lapply(adj.list, function(neighbors) rep(NA_real_, length(neighbors)))
    edge.weights <- numeric(nrow(edge.list))

    for (e in seq_len(nrow(edge.list))) {
        i <- as.integer(edge.list[e, 1L])
        j <- as.integer(edge.list[e, 2L])

        if (is.na(i) || is.na(j) || i < 1L || j < 1L || i > n || j > n) {
            stop(sprintf("edge.list row %d contains invalid vertex indices.", e))
        }
        if (i == j) {
            stop(sprintf("edge.list row %d contains a self-edge (%d, %d).", e, i, j))
        }

        edge.mass <- max(edge.densities[e], edge.mass.floor)
        w <- 1.0 / edge.mass

        if (weight.type == "mass.sym") {
            mass.i <- max(vertex.densities[i], vertex.mass.floor)
            mass.j <- max(vertex.densities[j], vertex.mass.floor)
            w <- w / sqrt(mass.i * mass.j)
        }

        idx.ij <- neighbor.position[[i]][as.character(j)]
        idx.ji <- neighbor.position[[j]][as.character(i)]

        if (length(idx.ij) == 0L || is.na(idx.ij)) {
            stop(sprintf(
                "edge.list row %d (%d, %d) is not present in adj.list[[%d]].",
                e, i, j, i
            ))
        }
        if (length(idx.ji) == 0L || is.na(idx.ji)) {
            stop(sprintf(
                "edge.list row %d (%d, %d) is not present in adj.list[[%d]].",
                e, i, j, j
            ))
        }

        weight.list[[i]][idx.ij] <- w
        weight.list[[j]][idx.ji] <- w
        edge.weights[e] <- w
    }

    for (i in seq_len(n)) {
        if (length(weight.list[[i]]) > 0L && anyNA(weight.list[[i]])) {
            stop(sprintf(
                "Failed to assign Laplacian weights to all neighbors in adj.list[[%d]].",
                i
            ))
        }
    }

    list(
        adj.list = adj.list,
        weight.list = weight.list,
        edge.list = edge.list,
        edge.weights = edge.weights
    )
}
