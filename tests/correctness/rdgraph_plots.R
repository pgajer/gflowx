rdgraph.color.values <- function(z, palette = grDevices::hcl.colors(100L, "Viridis")) {
  z <- as.double(z)
  out <- rep("gray80", length(z))
  finite <- is.finite(z)
  if (!any(finite)) return(out)

  rng <- range(z[finite])
  if (diff(rng) <= sqrt(.Machine$double.eps)) {
    out[finite] <- palette[ceiling(length(palette) / 2)]
    return(out)
  }

  idx <- 1L + floor((z[finite] - rng[1]) / diff(rng) * (length(palette) - 1L))
  idx <- pmin(pmax(idx, 1L), length(palette))
  out[finite] <- palette[idx]
  out
}

rdgraph.color.values.limits <- function(z,
                                         limits,
                                         palette = grDevices::hcl.colors(100L, "Viridis")) {
  z <- as.double(z)
  out <- rep("gray80", length(z))
  finite <- is.finite(z)
  if (!any(finite)) return(out)

  rng <- as.double(limits)
  if (length(rng) != 2L || !all(is.finite(rng)) || diff(rng) <= sqrt(.Machine$double.eps)) {
    return(rdgraph.color.values(z, palette = palette))
  }

  idx <- 1L + floor((z[finite] - rng[1]) / diff(rng) * (length(palette) - 1L))
  idx <- pmin(pmax(idx, 1L), length(palette))
  out[finite] <- palette[idx]
  out
}

rdgraph.fit.k.value <- function(name) {
  as.integer(sub("^k", "", name))
}

rdgraph.ok.k.fits <- function(ks) {
  fits <- Filter(Negate(is.null), ks$fits)
  fits[order(vapply(names(fits), rdgraph.fit.k.value, integer(1L)))]
}

rdgraph.ok.lcc.fits <- function(ks) {
  fits <- Filter(Negate(is.null), ks$lcc.fits %||% list())
  if (!length(fits)) return(fits)
  fits[order(vapply(fits, function(x) as.integer(x$k), integer(1L)))]
}

rdgraph.ok.component.fits <- function(ks) {
  fits <- Filter(Negate(is.null), ks$component.fits %||% list())
  if (!length(fits)) return(fits)
  fits[order(vapply(fits, function(x) as.integer(x$k), integer(1L)))]
}

rdgraph.selected.k.labels <- function(ks) {
  if (is.null(ks$selection) || !nrow(ks$selection)) {
    return(data.frame(k = integer(), label = character()))
  }
  split.criteria <- split(ks$selection$criterion, ks$selection$k_selected)
  data.frame(
    k = as.integer(names(split.criteria)),
    label = vapply(split.criteria, paste, character(1L), collapse = ", "),
    stringsAsFactors = FALSE
  )
}

rdgraph.write.png <- function(path, width = 1200L, height = 800L, expr) {
  grDevices::png(path, width = width, height = height, res = 130)
  on.exit(grDevices::dev.off(), add = TRUE)
  force(expr)
}

rdgraph.plot.metric.series <- function(x, y, xlab, ylab, main) {
  finite <- is.finite(x) & is.finite(y)
  if (!any(finite)) {
    graphics::plot.new()
    graphics::title(main = main)
    graphics::text(0.5, 0.5, "No finite values", col = "gray40")
    return(invisible(NULL))
  }
  graphics::plot(x[finite], y[finite], type = "b", pch = 16,
                 xlab = xlab, ylab = ylab, main = main)
}

rdgraph.variant.yhat <- function(variant.result) {
  if (!is.null(variant.result$error)) return(NULL)
  as.double(variant.result$fit$fitted.values)
}

rdgraph.plot.signal.overlay <- function(case.result, out.dir) {
  path <- file.path(out.dir, paste0(case.result$case$id, "_signal_overlay.png"))
  data <- case.result$data
  variants <- case.result$variants
  n.variants <- length(variants)
  coord <- as.double(data$coord)
  ord <- order(coord)
  xlab <- if (identical(data$type, "circle")) "angle (radians)" else "coordinate"

  rdgraph.write.png(path, width = 1300L, height = if (n.variants > 1L) 650L else 800L, {
    old.par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old.par), add = TRUE)
    graphics::par(
      mfrow = c(1L, n.variants),
      mar = c(4.2, 4.2, 3.0, 1.0),
      mgp = c(2.4, 0.7, 0)
    )

    for (variant.name in names(variants)) {
      variant <- variants[[variant.name]]
      response <- as.double(variant$y)
      y.hat <- rdgraph.variant.yhat(variant)
      ylim <- range(c(data$y.true, response, y.hat), finite = TRUE)

      graphics::plot(
        coord[ord],
        data$y.true[ord],
        type = "l",
        lwd = 2,
        col = "black",
        xlab = xlab,
        ylab = "response",
        ylim = ylim,
        main = variant.name
      )
      graphics::points(coord, response, pch = 16, cex = 0.45, col = grDevices::adjustcolor("gray45", 0.55))
      if (!is.null(y.hat)) {
        graphics::lines(coord[ord], y.hat[ord], lwd = 2, col = "#D55E00")
      }
      graphics::legend(
        "topright",
        legend = c("truth", "observed", "fit"),
        col = c("black", "gray45", "#D55E00"),
        lty = c(1, NA, 1),
        pch = c(NA, 16, NA),
        bty = "n",
        cex = 0.85
      )
    }
  })

  list(path = path, title = "Signal overlay")
}

rdgraph.plot.residuals <- function(case.result, out.dir) {
  path <- file.path(out.dir, paste0(case.result$case$id, "_residuals.png"))
  data <- case.result$data
  variants <- case.result$variants
  n.variants <- length(variants)
  coord <- as.double(data$coord)
  xlab <- if (identical(data$type, "circle")) "angle (radians)" else "coordinate"

  rdgraph.write.png(path, width = 1300L, height = if (n.variants > 1L) 650L else 800L, {
    old.par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old.par), add = TRUE)
    graphics::par(
      mfrow = c(1L, n.variants),
      mar = c(4.2, 4.2, 3.0, 1.0),
      mgp = c(2.4, 0.7, 0)
    )

    for (variant.name in names(variants)) {
      variant <- variants[[variant.name]]
      y.hat <- rdgraph.variant.yhat(variant)
      if (is.null(y.hat)) {
        graphics::plot.new()
        graphics::title(main = paste(variant.name, "fit failed"))
        next
      }
      residual <- y.hat - data$y.true
      graphics::plot(
        coord,
        residual,
        pch = 16,
        cex = 0.55,
        col = grDevices::adjustcolor("#0072B2", 0.7),
        xlab = xlab,
        ylab = "fit - truth",
        main = variant.name
      )
      graphics::abline(h = 0, lty = 2, col = "gray40")
    }
  })

  list(path = path, title = "Fit residuals")
}

rdgraph.plot.circle.geometry <- function(case.result, out.dir) {
  if (!identical(case.result$data$type, "circle") && !identical(case.result$data$type, "plane")) return(NULL)

  variant.name <- names(case.result$variants)[1L]
  variant <- case.result$variants[[variant.name]]
  y.hat <- rdgraph.variant.yhat(variant)
  if (is.null(y.hat)) return(NULL)

  path <- file.path(out.dir, paste0(case.result$case$id, "_geometry_panels.png"))
  data <- case.result$data
  response <- as.double(variant$y)
  residual <- y.hat - data$y.true
  z.list <- list(
    truth = data$y.true,
    observed = response,
    fit = y.hat,
    residual = residual
  )

  rdgraph.write.png(path, width = 1100L, height = 1000L, {
    old.par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old.par), add = TRUE)
    graphics::par(mfrow = c(2L, 2L), mar = c(2.5, 2.5, 3.0, 1.0), mgp = c(1.8, 0.5, 0))

    for (nm in names(z.list)) {
      z <- z.list[[nm]]
      graphics::plot(
        data$X,
        asp = 1,
        pch = 16,
        cex = 0.85,
        col = rdgraph.color.values(z),
        xlab = "x",
        ylab = "y",
        main = nm
      )
    }
  })

  list(path = path, title = "Geometry panels")
}

rdgraph.plot.k.metrics <- function(case.result, out.dir) {
  ks <- case.result$k.sweep
  if (is.null(ks) || is.null(ks$summary) || !nrow(ks$summary)) return(NULL)
  path <- file.path(out.dir, paste0(case.result$case$id, "_k_metrics.png"))
  s <- ks$summary
  ok <- s$fit_status == "ok"

  rdgraph.write.png(path, width = 1300L, height = 900L, {
    old.par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old.par), add = TRUE)
    graphics::par(mfrow = c(2L, 2L), mar = c(4.2, 4.2, 3.0, 1.0), mgp = c(2.4, 0.7, 0))

    rdgraph.plot.metric.series(s$k[ok], s$rmse_fit[ok],
                               xlab = "k", ylab = "RMSE(fit, truth)", main = "Oracle RMSE")
    rdgraph.plot.metric.series(s$k[ok], s$gcv[ok],
                               xlab = "k", ylab = "GCV", main = "GCV")
    rdgraph.plot.metric.series(s$k[ok], s$min_norm_gcv[ok],
                               xlab = "k", ylab = "min-normalized GCV", main = "Min-normalized GCV")
    rdgraph.plot.metric.series(s$k[ok], s$cor_fit_truth[ok],
                               xlab = "k", ylab = "cor(fit, truth)", main = "Correlation")
  })

  list(path = path, title = "k selection metrics")
}

rdgraph.plot.k.predictions.all.panels <- function(case.result, out.dir) {
  ks <- case.result$k.sweep
  if (is.null(ks) || !length(ks$fits)) return(NULL)
  path <- file.path(out.dir, paste0(case.result$case$id, "_k_predictions_all_panels.png"))
  data <- case.result$data
  fits <- rdgraph.ok.k.fits(ks)
  if (!length(fits)) return(NULL)
  selected <- rdgraph.selected.k.labels(ks)
  selected.k <- selected$k

  if (identical(data$type, "line") || identical(data$type, "circle")) {
    coord <- as.double(data$coord)
    ord <- order(coord)
    xlab <- if (identical(data$type, "circle")) "angle (radians)" else "coordinate"
    n.panels <- length(fits)
    nr <- ceiling(sqrt(n.panels))
    nc <- ceiling(n.panels / nr)
    ylim <- range(c(data$y.true, unlist(lapply(fits, `[[`, "fitted.values"))), finite = TRUE)
    rdgraph.write.png(path, width = 1500L, height = 1200L, {
      old.par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old.par), add = TRUE)
      graphics::par(mfrow = c(nr, nc), mar = c(3.4, 3.4, 2.7, 0.8), mgp = c(2.0, 0.6, 0))
      for (nm in names(fits)) {
        k <- rdgraph.fit.k.value(nm)
        y.hat <- as.double(fits[[nm]]$fitted.values)
        row <- ks$summary[ks$summary$k == k, , drop = FALSE]
        title <- sprintf("k = %d | RMSE %.3g | GCV %.3g", k, row$rmse_fit, row$gcv)
        if (k %in% selected.k) title <- paste0("* ", title)
        graphics::plot(coord[ord], data$y.true[ord], type = "l", lwd = 2,
                       col = "black", xlab = xlab, ylab = "response",
                       ylim = ylim, main = title)
        graphics::lines(coord[ord], y.hat[ord], lwd = if (k %in% selected.k) 2.4 else 1.6,
                        col = if (k %in% selected.k) "#D55E00" else "#0072B2")
      }
    })
  } else {
    n.panels <- length(fits)
    nr <- ceiling(sqrt(n.panels))
    nc <- ceiling(n.panels / nr)
    z.limits <- range(unlist(lapply(fits, `[[`, "fitted.values")), finite = TRUE)
    rdgraph.write.png(path, width = 1400L, height = 1200L, {
      old.par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old.par), add = TRUE)
      graphics::par(mfrow = c(nr, nc), mar = c(2.5, 2.5, 2.8, 0.8), mgp = c(1.7, 0.5, 0))
      for (nm in names(fits)) {
        k <- rdgraph.fit.k.value(nm)
        row <- ks$summary[ks$summary$k == k, , drop = FALSE]
        z <- as.double(fits[[nm]]$fitted.values)
        title <- sprintf("k = %d | RMSE %.3g | GCV %.3g", k, row$rmse_fit, row$gcv)
        if (k %in% selected.k) title <- paste0("* ", title)
        graphics::plot(data$X, asp = 1, pch = 16, cex = 0.75,
                       col = rdgraph.color.values.limits(z, z.limits),
                       xlab = "x", ylab = "y", main = title)
      }
    })
  }

  list(path = path, title = "Predictions across k: all panels")
}

rdgraph.plot.k.predictions.selected <- function(case.result, out.dir) {
  ks <- case.result$k.sweep
  if (is.null(ks) || !length(ks$fits)) return(NULL)
  selected <- rdgraph.selected.k.labels(ks)
  if (!nrow(selected)) return(NULL)
  fits <- rdgraph.ok.k.fits(ks)
  selected <- selected[selected$k %in% vapply(names(fits), rdgraph.fit.k.value, integer(1L)), , drop = FALSE]
  if (!nrow(selected)) return(NULL)

  path <- file.path(out.dir, paste0(case.result$case$id, "_k_predictions_selected.png"))
  data <- case.result$data

  if (identical(data$type, "line") || identical(data$type, "circle")) {
    coord <- as.double(data$coord)
    ord <- order(coord)
    xlab <- if (identical(data$type, "circle")) "angle (radians)" else "coordinate"
    selected.names <- paste0("k", selected$k)
    ylim <- range(c(data$y.true, unlist(lapply(fits, `[[`, "fitted.values"))), finite = TRUE)

    rdgraph.write.png(path, width = 1400L, height = 850L, {
      old.par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old.par), add = TRUE)
      graphics::par(mar = c(4.2, 4.2, 3.0, 1.0), mgp = c(2.4, 0.7, 0))
      graphics::plot(coord[ord], data$y.true[ord], type = "l", lwd = 3,
                     col = "black", xlab = xlab, ylab = "response",
                     ylim = ylim, main = "Selected k fits")
      for (nm in names(fits)) {
        y.hat <- as.double(fits[[nm]]$fitted.values)
        graphics::lines(coord[ord], y.hat[ord], col = grDevices::adjustcolor("gray65", 0.35), lwd = 0.8)
      }
      cols <- grDevices::hcl.colors(nrow(selected), "Dark 3")
      for (i in seq_len(nrow(selected))) {
        nm <- selected.names[[i]]
        y.hat <- as.double(fits[[nm]]$fitted.values)
        graphics::lines(coord[ord], y.hat[ord], col = cols[[i]], lwd = 2.4)
      }
      graphics::legend(
        "topright",
        legend = c("truth", sprintf("k=%d: %s", selected$k, selected$label)),
        col = c("black", cols),
        lty = 1,
        lwd = c(3, rep(2.4, nrow(selected))),
        bty = "n",
        cex = 0.72
      )
    })
  } else {
    selected.names <- paste0("k", selected$k)
    selected.fits <- fits[selected.names]
    z.limits <- range(c(data$y.true, unlist(lapply(selected.fits, `[[`, "fitted.values"))), finite = TRUE)
    n.panels <- length(selected.fits) + 1L
    nr <- ceiling(sqrt(n.panels))
    nc <- ceiling(n.panels / nr)

    rdgraph.write.png(path, width = 1400L, height = 1000L, {
      old.par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old.par), add = TRUE)
      graphics::par(mfrow = c(nr, nc), mar = c(2.5, 2.5, 2.8, 0.8), mgp = c(1.7, 0.5, 0))
      graphics::plot(data$X, asp = 1, pch = 16, cex = 0.75,
                     col = rdgraph.color.values.limits(data$y.true, z.limits),
                     xlab = "x", ylab = "y", main = "truth")
      for (i in seq_along(selected.fits)) {
        z <- as.double(selected.fits[[i]]$fitted.values)
        title <- sprintf("k=%d: %s", selected$k[[i]], selected$label[[i]])
        graphics::plot(data$X, asp = 1, pch = 16, cex = 0.75,
                       col = rdgraph.color.values.limits(z, z.limits),
                       xlab = "x", ylab = "y", main = title)
      }
    })
  }

  list(path = path, title = "Selected k fits")
}

rdgraph.plot.k.residuals.selected <- function(case.result, out.dir) {
  ks <- case.result$k.sweep
  if (is.null(ks) || !length(ks$fits)) return(NULL)
  selected <- rdgraph.selected.k.labels(ks)
  if (!nrow(selected)) return(NULL)
  fits <- rdgraph.ok.k.fits(ks)
  selected <- selected[selected$k %in% vapply(names(fits), rdgraph.fit.k.value, integer(1L)), , drop = FALSE]
  if (!nrow(selected)) return(NULL)

  path <- file.path(out.dir, paste0(case.result$case$id, "_k_residuals_selected.png"))
  data <- case.result$data
  selected.names <- paste0("k", selected$k)
  selected.fits <- fits[selected.names]

  if (identical(data$type, "line") || identical(data$type, "circle")) {
    coord <- as.double(data$coord)
    xlab <- if (identical(data$type, "circle")) "angle (radians)" else "coordinate"
    residuals <- lapply(selected.fits, function(fit) as.double(fit$fitted.values) - data$y.true)
    ylim <- range(unlist(residuals), finite = TRUE)
    n.panels <- length(selected.fits)
    rdgraph.write.png(path, width = 1400L, height = if (n.panels > 2L) 900L else 650L, {
      old.par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old.par), add = TRUE)
      graphics::par(mfrow = c(ceiling(n.panels / 2L), min(2L, n.panels)),
                    mar = c(4.0, 4.0, 3.0, 1.0), mgp = c(2.3, 0.7, 0))
      for (i in seq_along(selected.fits)) {
        graphics::plot(coord, residuals[[i]], pch = 16, cex = 0.55,
                       col = grDevices::adjustcolor("#0072B2", 0.75),
                       xlab = xlab, ylab = "fit - truth", ylim = ylim,
                       main = sprintf("k=%d: %s", selected$k[[i]], selected$label[[i]]))
        graphics::abline(h = 0, lty = 2, col = "gray40")
      }
    })
  } else {
    residuals <- lapply(selected.fits, function(fit) as.double(fit$fitted.values) - data$y.true)
    lim <- max(abs(unlist(residuals)), na.rm = TRUE)
    n.panels <- length(selected.fits)
    nr <- ceiling(sqrt(n.panels))
    nc <- ceiling(n.panels / nr)
    rdgraph.write.png(path, width = 1300L, height = 1000L, {
      old.par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old.par), add = TRUE)
      graphics::par(mfrow = c(nr, nc), mar = c(2.5, 2.5, 2.8, 0.8), mgp = c(1.7, 0.5, 0))
      for (i in seq_along(selected.fits)) {
        graphics::plot(data$X, asp = 1, pch = 16, cex = 0.75,
                       col = rdgraph.color.values.limits(residuals[[i]], c(-lim, lim),
                                                         palette = grDevices::hcl.colors(100L, "Blue-Red 3")),
                       xlab = "x", ylab = "y",
                       main = sprintf("k=%d residual: %s", selected$k[[i]], selected$label[[i]]))
      }
    })
  }

  list(path = path, title = "Selected k residuals")
}

rdgraph.plot.k.lcc.predictions <- function(case.result, out.dir) {
  ks <- case.result$k.sweep
  if (is.null(ks)) return(NULL)
  lcc.fits <- rdgraph.ok.lcc.fits(ks)
  if (!length(lcc.fits)) return(NULL)

  path <- file.path(out.dir, paste0(case.result$case$id, "_k_lcc_predictions.png"))
  data <- case.result$data

  if (identical(data$type, "line") || identical(data$type, "circle")) {
    coord <- as.double(data$coord)
    ord <- order(coord)
    xlab <- if (identical(data$type, "circle")) "angle (radians)" else "coordinate"
    yhat.values <- unlist(lapply(lcc.fits, function(x) x$fit$fitted.values))
    ylim <- range(c(data$y.true, yhat.values), finite = TRUE)
    n.panels <- length(lcc.fits)

    rdgraph.write.png(path, width = 1400L, height = if (n.panels > 3L) 900L else 650L, {
      old.par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old.par), add = TRUE)
      graphics::par(mfrow = c(ceiling(n.panels / 3L), min(3L, n.panels)),
                    mar = c(4.0, 4.0, 3.0, 1.0), mgp = c(2.3, 0.7, 0))

      for (lcc in lcc.fits) {
        row <- ks$summary[ks$summary$k == lcc$k, , drop = FALSE]
        idx <- lcc$indices
        idx.ord <- idx[order(coord[idx])]
        graphics::plot(coord[ord], data$y.true[ord], type = "l", lwd = 2,
                       col = grDevices::adjustcolor("black", 0.5),
                       xlab = xlab, ylab = "response", ylim = ylim,
                       main = sprintf("k=%d | LCC %d/%d (%.1f%%) | RMSE %.3g",
                                      lcc$k, row$lcc_n, nrow(data$X),
                                      100 * row$lcc_fraction, row$lcc_rmse_fit))
        graphics::points(coord[-idx], data$y.true[-idx], pch = 16, cex = 0.45,
                         col = grDevices::adjustcolor("gray60", 0.45))
        graphics::points(coord[idx], data$y.true[idx], pch = 16, cex = 0.45,
                         col = grDevices::adjustcolor("black", 0.55))
        graphics::lines(coord[idx.ord], as.double(lcc$fit$fitted.values)[order(coord[idx])],
                        col = "#D55E00", lwd = 2.3)
      }
    })
  } else {
    yhat.values <- unlist(lapply(lcc.fits, function(x) x$fit$fitted.values))
    z.limits <- range(c(data$y.true, yhat.values), finite = TRUE)
    n.panels <- length(lcc.fits) + 1L
    nr <- ceiling(sqrt(n.panels))
    nc <- ceiling(n.panels / nr)

    rdgraph.write.png(path, width = 1400L, height = 1100L, {
      old.par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old.par), add = TRUE)
      graphics::par(mfrow = c(nr, nc), mar = c(2.5, 2.5, 2.8, 0.8), mgp = c(1.7, 0.5, 0))
      graphics::plot(data$X, asp = 1, pch = 16, cex = 0.75,
                     col = rdgraph.color.values.limits(data$y.true, z.limits),
                     xlab = "x", ylab = "y", main = "truth")

      for (lcc in lcc.fits) {
        row <- ks$summary[ks$summary$k == lcc$k, , drop = FALSE]
        z <- rep(NA_real_, nrow(data$X))
        z[lcc$indices] <- as.double(lcc$fit$fitted.values)
        graphics::plot(data$X, asp = 1, pch = 16, cex = 0.75,
                       col = rdgraph.color.values.limits(z, z.limits),
                       xlab = "x", ylab = "y",
                       main = sprintf("k=%d | LCC %d/%d (%.1f%%) | RMSE %.3g",
                                      lcc$k, row$lcc_n, nrow(data$X),
                                      100 * row$lcc_fraction, row$lcc_rmse_fit))
      }
    })
  }

  list(path = path, title = "Largest-component predictions for disconnected k")
}

rdgraph.plot.k.lcc.residuals <- function(case.result, out.dir) {
  ks <- case.result$k.sweep
  if (is.null(ks)) return(NULL)
  lcc.fits <- rdgraph.ok.lcc.fits(ks)
  if (!length(lcc.fits)) return(NULL)

  path <- file.path(out.dir, paste0(case.result$case$id, "_k_lcc_residuals.png"))
  data <- case.result$data

  if (identical(data$type, "line") || identical(data$type, "circle")) {
    coord <- as.double(data$coord)
    xlab <- if (identical(data$type, "circle")) "angle (radians)" else "coordinate"
    residuals <- lapply(lcc.fits, function(lcc) {
      as.double(lcc$fit$fitted.values) - data$y.true[lcc$indices]
    })
    ylim <- range(unlist(residuals), finite = TRUE)
    n.panels <- length(lcc.fits)

    rdgraph.write.png(path, width = 1400L, height = if (n.panels > 3L) 900L else 650L, {
      old.par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old.par), add = TRUE)
      graphics::par(mfrow = c(ceiling(n.panels / 3L), min(3L, n.panels)),
                    mar = c(4.0, 4.0, 3.0, 1.0), mgp = c(2.3, 0.7, 0))

      for (i in seq_along(lcc.fits)) {
        lcc <- lcc.fits[[i]]
        row <- ks$summary[ks$summary$k == lcc$k, , drop = FALSE]
        graphics::plot(coord[lcc$indices], residuals[[i]], pch = 16, cex = 0.6,
                       col = grDevices::adjustcolor("#0072B2", 0.75),
                       xlab = xlab, ylab = "LCC fit - truth", ylim = ylim,
                       main = sprintf("k=%d | LCC %.1f%% | RMSE %.3g",
                                      lcc$k, 100 * row$lcc_fraction, row$lcc_rmse_fit))
        graphics::abline(h = 0, lty = 2, col = "gray40")
      }
    })
  } else {
    residual.values <- unlist(lapply(lcc.fits, function(lcc) {
      as.double(lcc$fit$fitted.values) - data$y.true[lcc$indices]
    }))
    lim <- max(abs(residual.values), na.rm = TRUE)
    n.panels <- length(lcc.fits)
    nr <- ceiling(sqrt(n.panels))
    nc <- ceiling(n.panels / nr)

    rdgraph.write.png(path, width = 1300L, height = 1000L, {
      old.par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old.par), add = TRUE)
      graphics::par(mfrow = c(nr, nc), mar = c(2.5, 2.5, 2.8, 0.8), mgp = c(1.7, 0.5, 0))
      for (lcc in lcc.fits) {
        row <- ks$summary[ks$summary$k == lcc$k, , drop = FALSE]
        z <- rep(NA_real_, nrow(data$X))
        z[lcc$indices] <- as.double(lcc$fit$fitted.values) - data$y.true[lcc$indices]
        graphics::plot(data$X, asp = 1, pch = 16, cex = 0.75,
                       col = rdgraph.color.values.limits(z, c(-lim, lim),
                                                         palette = grDevices::hcl.colors(100L, "Blue-Red 3")),
                       xlab = "x", ylab = "y",
                       main = sprintf("k=%d | LCC %.1f%% | RMSE %.3g",
                                      lcc$k, 100 * row$lcc_fraction, row$lcc_rmse_fit))
      }
    })
  }

  list(path = path, title = "Largest-component residuals for disconnected k")
}

rdgraph.plot.k.componentwise.predictions <- function(case.result, out.dir) {
  ks <- case.result$k.sweep
  if (is.null(ks)) return(NULL)
  component.fits <- rdgraph.ok.component.fits(ks)
  if (!length(component.fits)) return(NULL)

  path <- file.path(out.dir, paste0(case.result$case$id, "_k_componentwise_predictions.png"))
  data <- case.result$data
  z.limits <- range(c(data$y.true, unlist(lapply(component.fits, `[[`, "y.hat"))), finite = TRUE)
  n.panels <- length(component.fits) + 1L
  nr <- ceiling(sqrt(n.panels))
  nc <- ceiling(n.panels / nr)

  rdgraph.write.png(path, width = 1400L, height = 1100L, {
    old.par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old.par), add = TRUE)
    graphics::par(mfrow = c(nr, nc), mar = c(2.5, 2.5, 2.8, 0.8), mgp = c(1.7, 0.5, 0))
    graphics::plot(data$X, asp = if (ncol(data$X) == 2L) 1 else NA,
                   pch = 16, cex = 0.75,
                   col = rdgraph.color.values.limits(data$y.true, z.limits),
                   xlab = "x", ylab = if (ncol(data$X) == 2L) "y" else "",
                   main = "truth")

    for (component.fit in component.fits) {
      row <- ks$summary[ks$summary$k == component.fit$k, , drop = FALSE]
      title <- sprintf(
        "k=%d | %d comps | %.1f%% | RMSE %.3g",
        component.fit$k,
        row$component_fit_components,
        100 * row$component_fit_fraction,
        row$component_fit_rmse
      )
      graphics::plot(data$X, asp = if (ncol(data$X) == 2L) 1 else NA,
                     pch = 16, cex = 0.75,
                     col = rdgraph.color.values.limits(component.fit$y.hat, z.limits),
                     xlab = "x", ylab = if (ncol(data$X) == 2L) "y" else "",
                     main = title)
    }
  })

  list(path = path, title = "Component-wise predictions for disconnected k")
}

rdgraph.plot.k.componentwise.residuals <- function(case.result, out.dir) {
  ks <- case.result$k.sweep
  if (is.null(ks)) return(NULL)
  component.fits <- rdgraph.ok.component.fits(ks)
  if (!length(component.fits)) return(NULL)

  path <- file.path(out.dir, paste0(case.result$case$id, "_k_componentwise_residuals.png"))
  data <- case.result$data
  residual.values <- unlist(lapply(component.fits, function(component.fit) {
    component.fit$y.hat - data$y.true
  }))
  lim <- max(abs(residual.values), na.rm = TRUE)
  n.panels <- length(component.fits)
  nr <- ceiling(sqrt(n.panels))
  nc <- ceiling(n.panels / nr)

  rdgraph.write.png(path, width = 1300L, height = 1000L, {
    old.par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old.par), add = TRUE)
    graphics::par(mfrow = c(nr, nc), mar = c(2.5, 2.5, 2.8, 0.8), mgp = c(1.7, 0.5, 0))
    for (component.fit in component.fits) {
      row <- ks$summary[ks$summary$k == component.fit$k, , drop = FALSE]
      residual <- component.fit$y.hat - data$y.true
      title <- sprintf(
        "k=%d | %d comps | RMSE %.3g",
        component.fit$k,
        row$component_fit_components,
        row$component_fit_rmse
      )
      graphics::plot(data$X, asp = if (ncol(data$X) == 2L) 1 else NA,
                     pch = 16, cex = 0.75,
                     col = rdgraph.color.values.limits(
                       residual,
                       c(-lim, lim),
                       palette = grDevices::hcl.colors(100L, "Blue-Red 3")
                     ),
                     xlab = "x", ylab = if (ncol(data$X) == 2L) "y" else "",
                     main = title)
    }
  })

  list(path = path, title = "Component-wise residuals for disconnected k")
}

rdgraph.plot.k.predictions.ranked <- function(case.result, out.dir, n.show = 6L) {
  ks <- case.result$k.sweep
  if (is.null(ks) || !length(ks$fits) || is.null(ks$summary)) return(NULL)
  fits <- rdgraph.ok.k.fits(ks)
  s <- ks$summary[ks$summary$fit_status == "ok", , drop = FALSE]
  if (!length(fits) || !nrow(s)) return(NULL)
  s <- s[order(s$rmse_fit), , drop = FALSE]
  s <- head(s, n.show)
  fit.names <- paste0("k", s$k)
  ranked.fits <- fits[fit.names]
  path <- file.path(out.dir, paste0(case.result$case$id, "_k_predictions_ranked_rmse.png"))
  data <- case.result$data

  if (identical(data$type, "line") || identical(data$type, "circle")) {
    coord <- as.double(data$coord)
    ord <- order(coord)
    xlab <- if (identical(data$type, "circle")) "angle (radians)" else "coordinate"
    ylim <- range(c(data$y.true, unlist(lapply(ranked.fits, `[[`, "fitted.values"))), finite = TRUE)
    n.panels <- length(ranked.fits)
    rdgraph.write.png(path, width = 1400L, height = if (n.panels > 3L) 900L else 650L, {
      old.par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old.par), add = TRUE)
      graphics::par(mfrow = c(ceiling(n.panels / 3L), min(3L, n.panels)),
                    mar = c(4.0, 4.0, 3.0, 1.0), mgp = c(2.3, 0.7, 0))
      for (i in seq_along(ranked.fits)) {
        y.hat <- as.double(ranked.fits[[i]]$fitted.values)
        graphics::plot(coord[ord], data$y.true[ord], type = "l", lwd = 2,
                       col = "black", xlab = xlab, ylab = "response",
                       ylim = ylim,
                       main = sprintf("#%d: k=%d | RMSE %.3g", i, s$k[[i]], s$rmse_fit[[i]]))
        graphics::lines(coord[ord], y.hat[ord], col = "#D55E00", lwd = 2)
      }
    })
  } else {
    z.limits <- range(c(data$y.true, unlist(lapply(ranked.fits, `[[`, "fitted.values"))), finite = TRUE)
    n.panels <- length(ranked.fits) + 1L
    nr <- ceiling(sqrt(n.panels))
    nc <- ceiling(n.panels / nr)
    rdgraph.write.png(path, width = 1400L, height = 1100L, {
      old.par <- graphics::par(no.readonly = TRUE)
      on.exit(graphics::par(old.par), add = TRUE)
      graphics::par(mfrow = c(nr, nc), mar = c(2.5, 2.5, 2.8, 0.8), mgp = c(1.7, 0.5, 0))
      graphics::plot(data$X, asp = 1, pch = 16, cex = 0.75,
                     col = rdgraph.color.values.limits(data$y.true, z.limits),
                     xlab = "x", ylab = "y", main = "truth")
      for (i in seq_along(ranked.fits)) {
        z <- as.double(ranked.fits[[i]]$fitted.values)
        graphics::plot(data$X, asp = 1, pch = 16, cex = 0.75,
                       col = rdgraph.color.values.limits(z, z.limits),
                       xlab = "x", ylab = "y",
                       main = sprintf("#%d: k=%d | RMSE %.3g", i, s$k[[i]], s$rmse_fit[[i]]))
      }
    })
  }

  list(path = path, title = "RMSE-ranked prediction gallery")
}

rdgraph.plot.k.graphs <- function(case.result, out.dir) {
  ks <- case.result$k.sweep
  if (is.null(ks)) return(NULL)
  graphs <- Filter(Negate(is.null), ks$graphs %||% list())
  if (!length(graphs)) return(NULL)
  path <- file.path(out.dir, paste0(case.result$case$id, "_k_graphs.png"))
  data <- case.result$data
  n.panels <- length(graphs)
  nr <- ceiling(sqrt(n.panels))
  nc <- ceiling(n.panels / nr)

  rdgraph.write.png(path, width = 1400L, height = 1200L, {
    old.par <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(old.par), add = TRUE)
    graphics::par(mfrow = c(nr, nc), mar = c(2.5, 2.5, 2.8, 0.8), mgp = c(1.7, 0.5, 0))
    for (nm in names(graphs)) {
      graph <- graphs[[nm]]
      row <- ks$summary[ks$summary$k == rdgraph.fit.k.value(nm), , drop = FALSE]
      edges <- rdgraph.edge.table(graph$adj.list, graph$edge.length.list)
      lcc.indices <- if (row$component_count > 1L) {
        rdgraph.largest.component.indices(graph$adj.list)
      } else {
        seq_len(nrow(data$X))
      }
      point.col <- rep(grDevices::adjustcolor("gray60", 0.45), nrow(data$X))
      point.col[lcc.indices] <- "#1f2933"
      title <- if (row$component_count > 1L) {
        sprintf("%s | %d comps | LCC %.1f%%", nm, row$component_count, 100 * row$largest_component_fraction)
      } else {
        paste0(nm, " | connected")
      }
      if (ncol(data$X) == 1L) {
        graphics::plot(data$X[, 1], rep(0, nrow(data$X)), pch = 16, cex = 0.45,
                       col = point.col, xlab = "x", ylab = "", yaxt = "n", main = title)
        for (e in seq_len(nrow(edges))) {
          graphics::segments(data$X[edges$from[e], 1], 0, data$X[edges$to[e], 1], 0,
                             col = grDevices::adjustcolor("#0072B2", 0.35))
        }
        graphics::points(data$X[, 1], rep(0, nrow(data$X)), pch = 16, cex = 0.45, col = point.col)
      } else {
        graphics::plot(data$X, asp = 1, pch = 16, cex = 0.45, col = point.col,
                       xlab = "x1", ylab = "x2", main = title)
        for (e in seq_len(nrow(edges))) {
          graphics::segments(data$X[edges$from[e], 1], data$X[edges$from[e], 2],
                             data$X[edges$to[e], 1], data$X[edges$to[e], 2],
                             col = grDevices::adjustcolor("#0072B2", 0.35))
        }
        graphics::points(data$X, pch = 16, cex = 0.45, col = point.col)
      }
    }
  })

  list(path = path, title = "Graphs across k")
}

rdgraph.write.k.sweep.plots <- function(case.result, out.dir) {
  plots <- list(
    rdgraph.plot.k.metrics(case.result, out.dir),
    rdgraph.plot.k.predictions.all.panels(case.result, out.dir),
    rdgraph.plot.k.predictions.selected(case.result, out.dir),
    rdgraph.plot.k.residuals.selected(case.result, out.dir),
    rdgraph.plot.k.lcc.predictions(case.result, out.dir),
    rdgraph.plot.k.lcc.residuals(case.result, out.dir),
    rdgraph.plot.k.componentwise.predictions(case.result, out.dir),
    rdgraph.plot.k.componentwise.residuals(case.result, out.dir),
    rdgraph.plot.k.predictions.ranked(case.result, out.dir),
    rdgraph.plot.k.graphs(case.result, out.dir)
  )
  Filter(Negate(is.null), plots)
}

rdgraph.write.case.plots <- function(case.result, out.dir) {
  if (!dir.exists(out.dir)) dir.create(out.dir, recursive = TRUE)
  plots <- list(
    rdgraph.plot.signal.overlay(case.result, out.dir),
    rdgraph.plot.residuals(case.result, out.dir),
    rdgraph.plot.circle.geometry(case.result, out.dir)
  )
  Filter(Negate(is.null), plots)
}
