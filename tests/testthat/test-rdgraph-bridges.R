test_that("rdgraph bridge fits and refits through gflow", {
    skip_if_not_installed("gflow")

    set.seed(1)
    X <- matrix(runif(60), ncol = 3)
    y <- sin(X[, 1]) + X[, 2]

    fit <- fit.rdgraph.regression(
        X,
        y,
        k = 4,
        n.eigenpairs = 10,
        max.iterations = 1,
        verbose.level = 0,
        apply.geometric.pruning = FALSE
    )

    expect_s3_class(fit, "knn.riem.fit")
    expect_length(fit$fitted.values, nrow(X))

    refit <- refit.rdgraph.regression(fit, y)
    expect_s3_class(refit, "knn.riem.refit")
    expect_length(refit$fitted.values, nrow(X))
})
