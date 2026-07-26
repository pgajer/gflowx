# gflowx

`gflowx` is the installable archive for retired geometric smoothing and
conditional-expectation methods formerly shipped by `gflow`.

The package owns its rdgraph regression/refit engine, posterior and resampling
diagnostics, mean-shift/data smoothers, ulogit implementations, and legacy
response-graph selection helpers. It does not call private `gflow` code.
Generic graph construction is delegated to public `dgraphs` APIs.

The methods are preserved for reproducibility and comparison. They are not
under active methodological development.

```r
library(gflowx)

set.seed(1)
X <- matrix(rnorm(120), ncol = 3)
y <- X[, 1] - 0.5 * X[, 2] + rnorm(nrow(X), sd = 0.2)

fit <- fit.rdgraph.regression(
    X, y,
    k = 5,
    n.eigenpairs = 10,
    max.iterations = 2,
    verbose.level = 0
)

refit <- refit.rdgraph.regression(fit, y)
```

The optional `iknn.graph.response.pipeline()` wrapper still uses visualization
and layout functions from `gflow` and `grip`; those packages are only suggested
dependencies and are not required by the estimators.
