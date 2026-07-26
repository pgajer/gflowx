test_that("mean-shift smoother uses the local archived native backend", {
    set.seed(1)
    X <- matrix(runif(24), ncol = 2)

    fit <- meanshift.data.smoother(
        X,
        k = 3,
        density.k = 2,
        n.steps = 2,
        method = "basic"
    )

    expect_s3_class(fit, "MSD")
    expect_equal(dim(fit$dX), dim(X))
})

test_that("data smoother is available from gflowx", {
    expect_true(is.function(data.smoother))
})
