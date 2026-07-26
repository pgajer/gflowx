# Rdgraph Regression Correctness Lessons Learned

This document records empirical lessons learned while building and reviewing the
visual correctness suite for `gflowx::fit.rdgraph.regression()`. It is intended
to be updated after each new round of examples.

The companion design document describes the intended suite architecture. This
document is a running lab notebook: what the examples actually showed, what
diagnostics were useful, what was confusing, and how the next round of examples
should be shaped.

## Round 1: Initial Synthetic Examples

### Cases Added

The first round introduced four synthetic visual correctness cases:

- `1d_gaussian_mixture_recovers_smooth_signal`
- `circular_gaussian_mixture_recovers_periodic_signal`
- `shuffled_response_does_not_look_oracular`
- `2d_gaussian_mixture_recovers_smooth_signal`

These cases were intentionally small enough to run quickly but structured
enough to expose several different aspects of conditional expectation recovery:

- recovery of a smooth 1D response with known truth,
- recovery of a periodic response on noisy circular geometry,
- a negative control where the response-to-geometry relationship is destroyed,
- recovery of a smooth 2D surface from a subsampled grid.

### Main Lessons

#### Visual correctness needs multiple complementary views

The initial "Predictions across k" overlay was too dense. When many fitted
curves were drawn on the same axes, it became hard to see which k values were
actually close to the truth and whether the selected k values were visually
reasonable.

The report became more useful after replacing the single crowded overlay with
several complementary views:

- all-k small multiples,
- selected-k comparison,
- selected-k residual panels,
- oracle-RMSE-ranked prediction gallery.

These plots answer different questions:

- all-k panels show the global behavior of the sweep,
- selected-k plots show what the automated criteria would actually choose,
- residual panels show where estimates fail spatially,
- RMSE-ranked galleries show whether visual quality agrees with oracle metrics.

No single plot should be treated as sufficient for future cases. The report
should preserve this multi-view structure.

#### Graph construction is part of estimator correctness

The circular case showed that low k values can produce disconnected graphs, and
in those cases the regression fit may fail before conditional expectation
quality can be evaluated. This is not merely a nuisance. It is a correctness
signal about whether the graph construction supports the estimator.

The report now treats graph connectivity as a first-class diagnostic. For each
k, it records:

- number of connected components,
- largest connected component size,
- largest connected component fraction,
- second-largest component size.

This distinction matters because a poor estimate may be caused by the smoother,
by an unsuitable graph, or by both. Future cases should be designed so that the
report can separate these causes whenever possible.

#### Largest connected component fits are useful but must stay diagnostic

For disconnected k graphs, fitting on the largest connected component gives a
useful companion diagnostic. It answers:

> Given the part of the graph that actually supports connected diffusion, does
> the conditional expectation estimate behave reasonably?

This is especially helpful in cases like the noisy circle, where low-k graphs
may fragment. In the first circular example, disconnected low-k graphs could
still produce reasonable largest-component fits, but those fits used only a
subset of points. Therefore, largest connected component results must not be
silently substituted for full-graph results.

The report should keep the following separation:

- full-graph metrics describe the original requested fit,
- LCC metrics describe only vertices in the largest connected component,
- points outside the LCC are visibly marked or excluded in LCC plots,
- LCC fits are presented as diagnostic evidence, not as a successful full-data
  recovery.

This separation prevents a fragmented graph from looking better than it is.

#### Oracle graph baselines are highly informative

The weighted path graph for 1D examples and the weighted circle graph for the
circular example provide important anchors. They show what happens when the
graph is close to the known geometry rather than inferred only through the ikNN
construction.

Future examples should include an oracle or near-oracle graph whenever the
synthetic geometry allows one. Possible oracle graph baselines include:

- weighted path graph for ordered 1D samples,
- weighted circle graph for periodic samples,
- grid graph for structured 2D grids,
- known branch graph for Y-shaped or tree-like examples,
- planted component graph for disconnected examples,
- known bridge graph for weak-bridge topology examples.

Oracle graph fits are not meant to replace ikNN fits. They are reference points
that help interpret whether failures arise from graph construction or from the
regression method itself.

#### k selection should be evaluated as a disagreement surface

The first examples made it clear that k selection should not be reduced to a
single winning criterion too early. Oracle RMSE, GCV, normalized GCV, robust
normalized GCV, structural stability, visual smoothness, connectivity, and
residual structure can disagree.

The report is most useful when it makes these disagreements visible:

- selection table showing each criterion's selected k,
- metric curves across k,
- selected-k prediction plots,
- selected-k residual plots,
- RMSE-ranked prediction galleries,
- graph panels across k.

As the suite grows, disagreements between criteria should be treated as useful
signals. They can reveal cases where internal criteria select an oversmoothed
fit, a fragmented graph, or a visually plausible but biased estimate.

#### Negative controls are necessary

The shuffled-response case is important because it asks whether the method can
avoid looking oracular when the response-to-geometry relationship has been
destroyed. A correctness suite that only includes recoverable signals can be too
forgiving: a smoother may look impressive simply because all cases reward
smoothness.

Future negative controls should include:

- shuffled responses,
- pure noise responses,
- responses generated from a coordinate unrelated to the fitted geometry,
- response signals assigned to a deliberately wrong graph.

These cases should be reviewed visually and metrically. The desired behavior is
not "excellent recovery," but a report that clearly shows the absence of a
meaningful geometry-response signal.

### Report Design Lessons

The first round established a useful report contract. Each case should be
self-contained and reproducible from the report.

Every case should include:

- a human-readable title and stable case id,
- the synthetic data generation command with parameter values,
- the main `gflowx::fit.rdgraph.regression()` command with parameter values,
- any oracle graph construction and oracle fit commands,
- a metric table for the primary response variants,
- geometry/truth/response/fit/residual plots where applicable,
- k-sweep diagnostics over a small but informative k range,
- graph summaries for every k,
- selected-k diagnostics,
- graph panels across k,
- prediction panels across k,
- selected-k prediction comparison,
- selected-k residual panels,
- oracle-RMSE-ranked prediction gallery,
- LCC table and LCC plots whenever disconnected graphs occur.

The report should continue to prioritize visual review. Numeric metrics are
important, but in this stage they should support interpretation rather than act
as brittle pass/fail gates. The same cases can later be converted into
`testthat` checks once we know which metric ranges are stable and meaningful.

### Metrics Worth Keeping

The following metrics were useful in the first round and should remain in the
core report:

- raw response RMSE against truth,
- fitted response RMSE against truth,
- RMSE improvement,
- RMSE improvement fraction,
- correlation between fit and truth,
- variance ratio of fit to truth,
- residual mean,
- residual standard deviation,
- maximum absolute residual,
- boundary RMSE when a case has meaningful boundary or wraparound structure,
- GCV from `gflowx::fit.rdgraph.regression()`,
- min-normalized GCV,
- robust-normalized GCV,
- adjacent-k graph edit distance,
- adjacent-k degree Jensen-Shannon divergence,
- component count,
- largest component size,
- largest component fraction,
- second-largest component size,
- LCC-only RMSE, correlation, GCV, and edge count for disconnected graphs.

### Plot Types Worth Keeping

The following plot types should be treated as part of the standard visual
contract for future cases:

- signal overlay: truth, observed response, fitted response,
- residual plot: fitted response minus truth,
- geometry panels for 2D or manifold-like cases,
- k metric curves,
- all-k prediction small multiples,
- selected-k prediction comparison,
- selected-k residual panels,
- RMSE-ranked prediction gallery,
- graph panels across k,
- LCC prediction panels for disconnected k values,
- LCC residual panels for disconnected k values.

For 1D and circular cases, line plots against coordinate or angle are usually
most interpretable. For 2D cases, geometry-colored panels are usually more
informative. For future higher-dimensional examples, each case should define a
low-dimensional diagnostic coordinate or embedding for visualization.

## Round 2: Geometry, Boundary, And Negative-Control Stress

### Cases Added

The second round added four focused stress cases:

- `nonuniform_circle_recovers_periodic_signal`
- `circle_with_gap_reports_connectivity_and_wraparound_limits`
- `1d_boundary_peaks_exposes_endpoint_bias`
- `pure_noise_response_does_not_look_structured`

These cases moved beyond baseline smooth recovery and exercised density
imbalance, sampling gaps, endpoint bias, and constant-truth negative-control
behavior.

### Main Lessons

#### Round-specific reports are needed

After Round 2, the combined visual report became long enough that reading it as
a single artifact was no longer the best review experience. A combined report is
still useful for periodic full-suite review, but normal development should use
round-specific reports.

The report workflow should therefore produce:

- an index page linking all reports,
- one focused HTML report per round,
- a combined report for full-suite review,
- one metrics CSV for the full suite,
- one metrics CSV per round.

This prevents earlier examples from creating visual inflation as new case
families are added. Each round can be reviewed as a coherent experiment with a
bounded purpose.

#### Disconnected low-k behavior is common in stress cases

The Round 2 cases made it clear that disconnected graphs are not an edge case.
They occurred in every new k sweep:

- nonuniform circle: disconnected at `k = 3, 4`,
- circle with gap: disconnected at `k = 3, 4, 5, 6, 7`,
- 1D boundary peaks: disconnected at `k = 3, 4`,
- pure-noise response: disconnected at `k = 3, 4, 5`.

This confirms that LCC diagnostics should remain part of the standard report
contract. Connectivity is a graph-construction diagnostic, but it also governs
whether the regression problem is globally meaningful.

#### The circle-with-gap case is a strong topology stress test

The circle-with-gap case separates local recoverability from global graph
adequacy. Low k values produced disconnected graphs and therefore unsuitable
full-data fits, while higher k values eventually produced connected graphs with
interpretable recovery metrics.

This case is useful because it makes the following distinction visible:

- LCC-only fits can look reasonable on connected subsets,
- full-graph fits require a globally connected graph,
- a periodic oracle graph can recover well even when ikNN graph construction is
  sensitive to the sampling gap.

Future topology stress tests should preserve this pattern: show the graph
failure, show what can still be recovered on the LCC, and show an oracle graph
baseline when the intended topology is known.

#### Oracle graph baselines remain essential

Round 2 again showed that oracle graph baselines are not decorative. They are
often the most direct way to distinguish estimator behavior from graph
construction behavior.

In the boundary peaks case, the weighted path graph fit had much lower RMSE and
boundary RMSE than the ikNN fit. In the circle-with-gap and nonuniform-circle
cases, the weighted circle graph provided a useful reference for the intended
periodic geometry.

This supports a standing rule:

> If the synthetic geometry has a known graph, include an oracle or near-oracle
> graph fit.

When ikNN and oracle graph fits diverge, the report should treat that divergence
as a primary result, not as a nuisance.

#### Boundary stress is doing real work

The 1D boundary peaks case exposed a clear difference between the ikNN fit and
the weighted path graph baseline, especially in boundary RMSE. This suggests
that endpoint behavior should become a dedicated subfamily of correctness
examples rather than a single case.

Future boundary examples should vary:

- peak width near endpoints,
- endpoint signal amplitude,
- sample density near boundaries,
- whether the signal has endpoint peaks, endpoint slopes, or both.

The report should keep boundary RMSE for line cases, not only circular
wraparound cases.

#### Negative controls expose reporting assumptions

The pure-noise response case exposed two reporting fragilities:

- correlation with truth can be undefined when the true conditional expectation
  is constant,
- selection criteria can be all `NA` or non-finite.

The report was hardened so unavailable criteria are reported as unavailable and
metric plots can show "No finite values" rather than failing. This is a useful
lesson: negative controls test the report machinery as well as the estimator.

Future negative controls should deliberately include cases where some standard
metrics are undefined. The report should treat undefined metrics as valid
information, not exceptional failure.

#### Primary fit k and k-sweep k have different roles

The circle-with-gap example initially used a primary k that belonged to the
disconnected part of the sweep. That made the top-level case summary less
useful, because the headline ikNN fit errored even though the sweep later showed
connected, interpretable k values.

The better pattern is:

- choose the primary case fit at a plausible connected k,
- keep disconnected low-k values in the sweep,
- let the k sweep and LCC diagnostics document the graph failure region.

This gives the report both a readable headline fit and a useful stress-test
sweep.

### Implications For Round 3

Round 3 should probably focus on one of two directions:

- **2D density and boundary stress**: nonuniform 2D sampling, peaks near edges,
  and sparse boundary coverage.
- **Topology stress**: two components, weak bridge, and Y-shaped branching
  geometry.

The second direction may be more informative next because Round 2 showed how
central connectivity and graph topology are to interpreting regression
correctness. A weak-bridge or planted-component family would extend the LCC
diagnostic from "graph fragmented" to "graph connected but topologically
fragile."

## Proposed Modular Case Families For Next Rounds

Future examples should be added as modular families. Each family should expose
one or two specific stress signatures rather than trying to test everything at
once.

### Geometry Stress

Purpose: test whether graph construction and regression handle imperfect or
challenging geometry.

Candidate cases:

- nonuniform circle with periodic smooth truth,
- circle with a sampling gap,
- noisy circle with outliers,
- ellipse or distorted circle with periodic truth,
- two-dimensional manifold sampled with anisotropic noise.

Expected diagnostics:

- connectivity across k,
- LCC behavior for fragmented low-k graphs,
- wraparound residuals,
- comparison to weighted circle or known manifold graph when available.

### Boundary Stress

Purpose: test whether estimates degrade near geometric boundaries.

Candidate cases:

- 1D interval with peaks near endpoints,
- 1D interval with steep endpoint slope,
- 2D square with Gaussian peaks near edges or corners,
- nonuniform 2D sample with sparse boundary coverage.

Expected diagnostics:

- boundary RMSE,
- residuals near endpoints or edges,
- comparison of selected k against oracle RMSE,
- sensitivity to oversmoothing.

### Density Stress

Purpose: test whether estimates remain reasonable when sample density varies.

Candidate cases:

- same smooth truth under uniform and clustered X samples,
- 1D mixture with dense and sparse regions,
- 2D Gaussian mixture on nonuniform samples,
- circle with angular density imbalance.

Expected diagnostics:

- residuals stratified by dense versus sparse regions,
- variance ratio of fit to truth,
- graph degree summaries,
- selected k stability under density imbalance.

### Topology Stress

Purpose: test whether estimates behave around branches, bridges, and
disconnected or nearly disconnected structures.

Candidate cases:

- two separated components with independent smooth truths,
- two components with a weak bridge,
- Y-shaped graph with branch-specific responses,
- branching manifold with shared stem and diverging arms,
- planted components where some components are deliberately small.

Expected diagnostics:

- component count,
- largest and second-largest component sizes,
- LCC-only fits,
- bridge leakage or cross-component smoothing,
- oracle graph comparison when the planted topology is known.

### Signal Stress

Purpose: test response functions that challenge smoothing assumptions.

Candidate cases:

- broad plus narrow peaks,
- two nearby Gaussian peaks,
- steep but continuous transition,
- periodic signal with narrow wraparound feature,
- multi-scale response on a 2D surface.

Expected diagnostics:

- peak attenuation,
- oversmoothing,
- residual localization,
- RMSE-ranked gallery,
- disagreement between GCV-selected k and oracle-RMSE-selected k.

### Negative Controls

Purpose: test whether the method avoids producing misleadingly oracle-looking
fits when no valid geometry-response relationship exists.

Candidate cases:

- shuffled response,
- pure noise response,
- response generated from an unrelated latent coordinate,
- response generated on one graph but fitted on another,
- random labels on clustered geometry.

Expected diagnostics:

- weak or absent correlation with truth,
- limited RMSE improvement,
- visually non-oracular fits,
- comparison to the corresponding non-shuffled positive case.

### Parameter Sensitivity

Purpose: test robustness of the estimator to important tuning parameters.

Candidate cases:

- sweep over `filter.type`,
- sweep over `n.eigenpairs`,
- sweep over `response.penalty.exp`,
- sweep over geometric pruning parameters,
- sweep over counting measure versus density-related options when applicable.

Expected diagnostics:

- parameter-specific metric curves,
- selected parameter comparisons,
- visual change in residual structure,
- runtime and failure status summaries.

## Standard Report Contract For Future Cases

Every new case should satisfy this contract unless there is a clear reason to
deviate.

### Reproducibility

Each case must include:

- fixed random seed,
- data generation command with all important parameter values,
- response generation command with noise level and truth definition,
- main `gflowx::fit.rdgraph.regression()` command with all important parameter values,
- oracle graph construction command when applicable,
- oracle graph fit command when applicable.

### Data Fields

Each case should return:

- `type`: visualization type such as `line`, `circle`, `plane`, or a future
  type,
- `X`: fitted geometry matrix,
- `coord`: diagnostic coordinate for 1D plotting or ordering,
- `theta`: circular coordinate when applicable,
- `y.true`: known synthetic conditional expectation,
- `responses`: one or more observed response variants,
- `metadata`: noise levels, signal description, and geometry notes.

### Fit Variants

At minimum, each case should include:

- primary ikNN-based fit,
- oracle graph fit when a known graph is available,
- negative-control response variant when appropriate.

### k Sweep

Each case should define a small but informative `k.grid`. The grid should be
wide enough to include:

- under-connected or under-smoothed values when possible,
- plausible selected values,
- oversmoothed or overly dense values.

For each k, the report should record:

- fit status,
- error message if any,
- elapsed time,
- oracle RMSE,
- correlation with truth,
- GCV,
- graph edge count,
- component count,
- largest component size and fraction,
- second-largest component size,
- degree summaries,
- LCC-only fit status and metrics when disconnected.

### Visuals

Each case should produce:

- truth/response/fit overlay when meaningful,
- residual plot,
- geometry panels for 2D or manifold-like cases,
- k metric curves,
- graph panels across k,
- all-k prediction small multiples,
- selected-k prediction comparison,
- selected-k residual panels,
- RMSE-ranked prediction gallery,
- LCC prediction and residual panels for disconnected k values.

### Interpretation

After each round of examples, this lessons-learned document should be updated
with:

- cases added,
- failure modes exposed,
- diagnostics that were useful,
- diagnostics that were confusing or redundant,
- metric disagreements worth investigating,
- report design changes made,
- implications for the next round of cases.

This keeps the correctness suite from becoming a pile of examples. The goal is
to build an increasingly explicit understanding of when
`gflowx::fit.rdgraph.regression()` behaves well, when it fails, and which diagnostics
make those behaviors visible.
