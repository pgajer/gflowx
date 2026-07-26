# Rdgraph Regression Correctness Suite

## Purpose

This document sketches an extensible correctness suite for
`gflowx::fit.rdgraph.regression()` conditional expectation estimates. The first
goal is not to create brittle numerical unit tests, but to generate
human-readable reports that make estimator behavior easy to inspect against
known synthetic truth. Stable report metrics can later be promoted into
`testthat` regression tests.

The suite should answer questions such as:

- Does the fitted conditional expectation recover known smooth signal better
  than the noisy response?
- Does it preserve important geometric features such as peaks, troughs,
  periodic wraparound, and local variation?
- Does it fail gracefully on negative controls such as shuffled responses or
  pure-noise outcomes?
- Do future changes to `gflowx::fit.rdgraph.regression()` produce visible or measurable
  behavior drift?

## Guiding Principles

- Start visual-first: every correctness case should produce plots comparing
  truth, noisy observations, fitted estimates, and residuals.
- Emit structured metrics for every visual case, even before those metrics are
  used as pass/fail criteria.
- Keep the report runner deterministic: fixed seeds, fixed case parameters, and
  explicit package versions/session information.
- Separate correctness reports from normal CRAN-facing checks. The suite can be
  run locally or in extended CI without making `make check-fast` fragile.
- Prefer conservative, interpretable metrics over exact fitted-value snapshots.
- Promote only stable, well-understood metrics into `testthat` tests.

## Proposed Directory Layout

```text
tests/correctness/
  rdgraph-regression-correctness-suite.md
  rdgraph_regression_report.R
  rdgraph_cases.R
  rdgraph_metrics.R
  rdgraph_plots.R
  rdgraph_report_utils.R
```

Suggested generated outputs:

```text
tests/manual/reports/rdgraph-regression-correctness/
  rdgraph_regression_correctness.html
  metrics.csv
  plots/
```

Generated reports, images, and metric tables should usually stay out of
functional commits unless the task is explicitly about publishing or comparing
report artifacts.

## Case Registry Design

Each correctness case should be registered as a small object with a common
interface. The report runner should not need case-specific logic beyond calling
standard hooks.

Example shape:

```r
list(
  id = "1d_gaussian_mixture_recovers_smooth_signal",
  title = "1D Gaussian Mixture: Smooth Signal Recovery",
  group = "signal_recovery",
  seed = 1001L,
  generate = function() {
    # Returns X, y, y.true, plotting coordinates, and metadata.
  },
  fit.args = list(
    k = 5L,
    max.iterations = 5L,
    n.eigenpairs = 30L,
    pca.dim = NULL,
    verbose.level = 0L
  ),
  plot = function(result, out.dir) {
    # Writes one or more plot files and returns metadata for the report.
  },
  metrics = function(result) {
    # Returns a one-row data.frame of metrics.
  }
)
```

The runner should:

1. Load the package with `pkgload::load_all(".", quiet = TRUE)`.
2. Iterate over registered cases.
3. Generate data.
4. Fit `gflowx::fit.rdgraph.regression()` or `gflowx::refit.rdgraph.regression()` as specified.
5. Fit any case-specific oracle graph models, such as an ordered path graph for
   1D data or a weighted circle graph for circular data.
6. Run optional k-sweeps over ikNN graph construction parameters.
7. Compute metrics.
8. Write plots.
9. Assemble a self-contained HTML report.
10. Write the combined metrics table.

## Initial Cases

### 1D Gaussian Mixture Recovers Smooth Signal

Use `generate.1d.gaussian.mixture()` or a direct equivalent to create a known
smooth response, then add Gaussian noise.

Visuals:

- Truth, noisy response, and fitted response over the 1D coordinate.
- Residual `y.hat - y.true` over the 1D coordinate.
- Optional graph layout colored by truth and fitted values.

Metrics:

- `rmse_raw = RMSE(y, y.true)`
- `rmse_fit = RMSE(y.hat, y.true)`
- `rmse_improvement = rmse_raw - rmse_fit`
- `cor_fit_truth = cor(y.hat, y.true)`
- `var_ratio_fit_truth = var(y.hat) / var(y.true)`

### Circular Gaussian Mixture Recovers Periodic Signal

Use `generate.circle.data()` plus `circular.synthetic.mixture.of.gaussians()` to
create a periodic response field on noisy circle samples.

Visuals:

- Truth, noisy response, and fitted response over angle.
- Residual over angle, with attention to the `0 / 2pi` boundary.
- Geometry panels colored by `y.true`, `y`, `y.hat`, and residual.

Metrics:

- Same recovery metrics as the 1D case.
- Boundary residual summaries near angles close to `0` and `2pi`.

### Shuffled Response Does Not Look Oracular

Use one of the same signal-generating datasets, then compare a fit to the real
response with a fit to a shuffled response.

Visuals:

- Real-response fit versus truth.
- Shuffled-response fit versus truth.
- Residual comparison panels.

Metrics:

- Real-response recovery metrics.
- Shuffled-response recovery metrics.
- Difference between real and shuffled RMSE/correlation.

This case is a negative control. It should show that the estimator is not
recovering the oracle signal when the response-to-geometry relationship has been
destroyed.

### 2D Gaussian Mixture Recovers Smooth Signal

Use `create.gaussian.mixture()` on a 2D grid/subsample with known truth and
add Gaussian response noise.

Visuals:

- Geometry panels colored by truth, noisy response, fitted response, and
  residuals.
- k-sweep panels showing ikNN graphs and predictions across k.

Metrics:

- Same recovery metrics as the 1D case.
- k-sweep diagnostics over a small graph-construction range.

## k Selection Diagnostics

For synthetic correctness cases, the report can show both oracle diagnostics and
model-internal diagnostics.

- **Oracle RMSE**: \(\sqrt{n^{-1}\sum_i(\hat y_i-y_i^{true})^2}\). Lower is
  better; available only because these cases have known truth.
- **GCV**: generalized cross-validation reported by
  `gflowx::fit.rdgraph.regression()` at the selected iteration. Lower is better.
- **Min-normalized GCV**: `GCV(k) / min_k GCV(k)`. Lower is better; the best
  possible value is 1.
- **Robust-normalized GCV**: `(GCV(k) - median_k GCV(k)) / MAD_k GCV(k)`, with
  IQR, SD, and constant fallbacks for degenerate scales. Lower is better.
- **Graph structural stability**: adjacent-k graph edit distance and
  Jensen-Shannon divergence of degree distributions. Smaller values indicate
  more stable adjacent graph structure.
- **Largest connected component fit**: if a k graph is disconnected, refit
  `gflowx::fit.rdgraph.regression()` on the largest connected component using the
  induced weighted graph. Report the component count, largest component size,
  largest component fraction, second-largest component size, and LCC-only
  recovery metrics. These metrics are computed only on vertices in the largest
  component and should be read as a diagnostic companion to the full-graph fit,
  not as a replacement for the disconnected-graph warning.

## Future Case Groups

Signal recovery:

- 2D Gaussian mixture on uniform samples.
- 2D Gaussian mixture on nonuniform samples.
- Branching or Y-shaped graph with piecewise smooth response.
- Multi-scale signal with one broad and one narrow feature.

Negative controls:

- Pure noise response.
- Shuffled response at multiple noise levels.
- Signal generated from a coordinate unrelated to the fitted geometry.

Geometry stress:

- Nonuniform density on a circle.
- Outliers near a smooth manifold.
- Sparse bridge between components.
- Boundary-heavy samples in a square or interval.

Parameter sensitivity:

- Multiple `k` values.
- Multiple `filter.type` settings.
- Multiple `n.eigenpairs` values.
- Multiple `response.penalty.exp` values.

Regression guards:

- Compare current metrics against saved baseline ranges.
- Flag large changes in fitted RMSE, correlation, variance ratio, selected
  tuning parameters, or graph edge counts.

## Relationship To `testthat`

The visual suite should come first. After report metrics are stable across
local runs and package changes, selected cases can be promoted into small,
conservative `testthat` checks.

Good future `testthat` assertions:

- Fitted RMSE is lower than raw noisy RMSE for a deterministic easy case.
- Fitted values are finite and not collapsed to a constant.
- Real-response fit correlates with truth substantially more than shuffled fit.
- Precomputed-graph and internally built graph paths agree when given the same
  graph.

Avoid early brittle assertions:

- Exact fitted values.
- Overly tight RMSE thresholds from a single seed.
- Thresholds that depend heavily on BLAS, sparse eigensolver behavior, or small
  stochastic differences.

## Suggested Make Target

The report should eventually be runnable with a Makefile target:

```make
check-correctness-report:
	Rscript tests/correctness/rdgraph_regression_report.R
```

This target should remain separate from CRAN-facing package checks unless a
later release decision explicitly folds part of the suite into fast tests.
