test_that("pipeline exports survive its temporary graphics directory", {
    skip_if_not_installed("ivue")
    skip_if_not_installed("rgl")
    skip_if_not_installed("grip")
    old.wd <- getwd()
    root <- tempfile("pipeline-paths-")
    dir.create(root)
    on.exit({setwd(old.wd); unlink(root, recursive = TRUE)}, add = TRUE)
    setwd(root)
    set.seed(812)
    X <- matrix(rnorm(80 * 3), 80, 3)
    y <- X[, 1]^2 - X[, 2]^2
    args <- list(X = X, y = y, k.min = 6L, k.max = 8L, selected.k = 8L,
      build.args = list(method = "none", trim.disconnected = FALSE),
      grip.args = list(rounds = 5L, final_rounds = 5L, num_init = 3L, num_nbrs = 15L),
      fit.args = list(max.iterations = 1L, n.eigenpairs = 10L, verbose.level = 0L),
      cont.plot.args = list(point.type = "point", selfcontained = FALSE),
      out.dir = "relative", timestamp = "fixture", save.rds = FALSE, verbose = FALSE)
    # The installed grip compatibility entry point emits a deprecation warning.
    result <- suppressWarnings(do.call(iknn.graph.response.pipeline, args))
    expect_true(length(result$html.objects) > 0L)
    expect_true(length(list.files("relative", pattern = "[.]html$")) > 0L)
    expect_identical(getwd(), normalizePath(root))
    args$cont.plot.args$output.file <- "explicit/scene.html"
    result <- suppressWarnings(do.call(iknn.graph.response.pipeline, args))
    expect_true(file.exists("explicit/scene.html"))
    expect_true(all(vapply(result$components, function(x) file.exists(x$html.file), logical(1))))
})
