source("/Users/pgajer/current_projects/gflow/tests/manual/00-setup.R")
test.data <- generate.test.data.1d(n.pts = 50, noise_sd = 1)
plot.test.data(test.data)

fit <- gflowx::fit.rdgraph.regression(
    test.data$X,
    test.data$y,
    k = 3,
    response.penalty.exp = 0.1,
    density.uniform.pull = 0.1,
    n.eigenpairs = 10,
    max.iterations = 10,
    ##filter.type = "cubic_spline",
    use.counting.measure = TRUE,
    density.normalization = 0,
    t.diffusion = 10,
    epsilon.y = 1e-4,
    epsilon.rho = 1e-4,
    max.ratio.threshold = 0.1,
    threshold.percentile = 0.5,
    test.stage = -1
)

fit <- gflowx::fit.rdgraph.regression(
    test.data$X,
    test.data$y,
    k = 3,
    with.posterior = TRUE,
    return.posterior.samples = TRUE,
    credible.level = 0.95,
    n.posterior.samples = 1000,
    posterior.seed = NULL,
    n.eigenpairs = 10,
    max.iterations = 10)

summary(fit)

plot(test.data$X[,1], test.data$y.smooth, xlab = "", ylab = "", type = "l", ylim = range(c(test.data$y, fit$fitted.values)))
points(test.data$X[,1], test.data$y)
lines(test.data$X[,1], fit$raw.fitted.values, col = "red")

plot(test.data$X[,1], test.data$y.smooth, xlab = "", ylab = "", type = "l", ylim = range(c(test.data$y, fit$fitted.values, fit$posterior$lower, fit$posterior$upper)))
points(test.data$X[,1], test.data$y)
lines(test.data$X[,1], fit$fitted.values, col = "red")
lines(test.data$X[,1], fit$posterior$lower, col = "gray")
lines(test.data$X[,1], fit$posterior$upper, col = "gray")
lines(test.data$X[,1], fit$posterior$samples[,1], col = "blue")

## Get summary statistics only (memory efficient)
fit1 <- gflowx::fit.rdgraph.regression(
    X, y, k = 20,
    with.posterior = TRUE,
    n.posterior.samples = 1000,
    return.posterior.samples = FALSE
)

## Access: fit1$posterior$lower, fit1$posterior$upper, fit1$posterior$sd

## Get individual samples (for downstream analysis)
fit2 <- gflowx::fit.rdgraph.regression(
    X, y, k = 20,
    with.posterior = TRUE,
    n.posterior.samples = 100,
    return.posterior.samples = TRUE
)

## Access: fit2$posterior$samples is an n × 100 matrix
## Each column is one posterior realization
posterior.mean <- rowMeans(fit2$posterior$samples)
posterior.quantiles <- apply(fit2$posterior$samples, 1, quantile, probs = c(0.025, 0.975))



plot(fit$gamma.selection$gamma.grid, fit$gamma.selection$gcv.scores, type = "l")

# Summary GCV information
fit$gcv$eta.optimal    # Vector of optimal η per iteration
fit$gcv$gcv.optimal    # Vector of optimal GCV scores per iteration

# Full GCV curves for detailed analysis
fit$gcv$eta.grid       # List of η grids (one per iteration)
fit$gcv$gcv.scores     # List of GCV score vectors (one per iteration)

# Plot GCV evolution
plot(1:length(fit$gcv$gcv.optimal), fit$gcv$gcv.optimal,
     type = "b", xlab = "Iteration", ylab = "Optimal GCV Score",
     main = "GCV Score Evolution")

# Plot GCV curves for specific iterations
op <- par(mfrow = c(2, 2), mar=c(3.5, 3.5, 1.5, 0.5), mgp=c(2.25,0.5,0),tcl = -0.3)
plot(1:length(fit$gcv$gcv.optimal), fit$gcv$gcv.optimal,
     type = "b", xlab = "Iteration", ylab = "Optimal GCV Score",
     main = "GCV Score Evolution")
for (iter in c(1, 5, 10, 20)) {
    if (iter <= length(fit$gcv$eta.grid)) {
        plot(log10(fit$gcv$eta.grid[[iter]]), fit$gcv$gcv.scores[[iter]],
             type = "l", xlab = "log10(eta)", ylab = "GCV Score",
             main = paste("Iteration", iter))
        abline(v = log10(fit$gcv$eta.optimal[iter]), col = "red", lty = 2)
    }
}
par(op)

timestamp <- generate.timestamp()
file <- "~/current_projects/gflow/tests/manual/pics/gcv_evolultion__"
(file <- paste0(file, timestamp, ".pdf"))
pdf(file, width=6, height=6)
op <- par(mfrow = c(2, 2), mar=c(3.5, 3.5, 1.5, 0.5), mgp=c(2.25,0.5,0),tcl = -0.3)
plot(1:length(fit$gcv$gcv.optimal), fit$gcv$gcv.optimal,
     type = "b", xlab = "Iteration", ylab = "Optimal GCV Score",
     main = "GCV Score Evolution")
for (iter in c(1, 5, 10, 20)) {
    if (iter <= length(fit$gcv$eta.grid)) {
        plot(log10(fit$gcv$eta.grid[[iter]]), fit$gcv$gcv.scores[[iter]],
             type = "l", xlab = "log10(eta)", ylab = "GCV Score",
             main = paste("Iteration", iter))
        abline(v = log10(fit$gcv$eta.optimal[iter]), col = "red", lty = 2)
    }
}
par(op)
dev.off()
system(paste0("open ",file))
## "~/current_projects/gflow/tests/manual/pics/gcv_evolultion__2025_10_07_214301.pdf"



## Compare GCV trajectories for different γ values
X <- test.data$X
y <- test.data$y
gamma_values <- c(0.1, 0.3, 0.5, 0.8)
gcv_trajectories <- lapply(gamma_values, function(g) {
    fit <- gflowx::fit.rdgraph.regression(X, y,
                                         response.penalty.exp = g,
                                         k = 3,
                                         n.eigenpairs = 10,
                                         max.iterations = 10,
                                         filter.type = "cubic_spline",
                                         use.counting.measure = TRUE,
                                         density.normalization = 0,
                                         t.diffusion = 10,
                                         density.uniform.pull = 0.01,
                                         epsilon.y = 1e-4,
                                         epsilon.rho = 1e-4,
                                         max.ratio.threshold = 0.1,
                                         threshold.percentile = 0.5,
                                         test.stage = -1
                                         )
    fit$gcv$gcv.optimal
})

matplot(do.call(cbind, gcv_trajectories), type = "l",
        xlab = "Iteration", ylab = "GCV Score",
        main = "GCV vs gamma")
legend("topright", legend = paste("gamma =", gamma_values),
       col = 1:4, lty = 1:4)

op <- par(mfrow = c(2, 2), mar=c(3.5, 3.5, 1.5, 0.5), mgp=c(2.25,0.5,0),tcl = -0.3)
for (i in 1:4) {
    plot(gcv_trajectories[[i]], type = "b", main = paste0("gamma: ", gamma_values[i]))
}
par(op)


timestamp <- generate.timestamp()
file <- "~/current_projects/gflow/tests/manual/pics/gcv_vs_gamma_"
(file <- paste0(file, timestamp, ".pdf"))
pdf(file, width=6, height=6)
op <- par(mfrow = c(2, 2), mar=c(3.5, 3.5, 1.5, 0.5), mgp=c(2.25,0.5,0),tcl = -0.3)
for (i in 1:4) {
    plot(gcv_trajectories[[i]], type = "b", main = paste0("gamma: ", gamma_values[i]))
}
par(op)
dev.off()
system(paste0("open ",file))
## "~/current_projects/gflow/tests/manual/pics/gcv_vs_gamma_2025_10_07_215358.pdf"
