# Round 3 Lessons Learned: Topology And Bridge Stress

Round 3 added topology-focused examples for disconnected components, weak
bridges, Y-shaped branching geometry, and branch-localized signal features.

## Cases Reviewed

- `two_components_independent_smooth_signals`
- `weak_bridge_between_components_tests_bottleneck_leakage`
- `y_shaped_branching_geometry_branch_specific_signal`
- `y_branch_localized_peak_tests_branch_leakage`

## Main Lessons

### LCC diagnostics are not enough for meaningful multi-component data

The two-component example stayed disconnected across the full k sweep. That is
not merely a graph-construction failure. It is the correct topology for that
synthetic data. The previous largest-connected-component diagnostic was useful
for accidental fragmentation, but it discarded the second component even when
both components were meaningful.

Round 3B therefore adds component-wise diagnostics: for disconnected graphs, fit
every connected component with at least 10 vertices and compute aggregate
component-wise metrics over all fitted vertices.

### Disconnected full-graph fits need interpretation

When a graph is intentionally disconnected, a full-graph
`gflowx::fit.rdgraph.regression()` error is informative. The report should distinguish:

- accidental fragmentation,
- intentional multi-component geometry,
- weakly connected geometry,
- connected but topologically fragile geometry.

Round 3B keeps the intentionally disconnected component case but shifts the
diagnostic emphasis from "largest component" to "all sufficiently large
components."

### The weak-bridge case was too easy

The original weak-bridge example fit well across the k sweep. The planted bridge
baseline did not outperform ikNN, which suggests that the bridge was not a hard
enough bottleneck and/or the planted graph was too sparse.

Round 3B replaces this with a longer, sparser bridge and a sharper response
transition. The goal is to make leakage and bottleneck smoothing visible.

### Planted topology baselines are hypothesis baselines

The planted Y-tree graph helped the Y-shaped examples, but the planted
weak-bridge graph did not help the original bridge example. This means planted
graphs should not automatically be interpreted as oracle truth. They are known
topology baselines: useful for diagnosis, but still subject to edge-density and
weighting choices.

### Branch-localized features need stronger stress

The original branch-localized peak case was useful, but Round 3B makes the peak
narrower and increases branch contrast. This should better reveal leakage from
one arm into another and oversmoothing near the branch point.

## Round 3B Design

Round 3B adds four revised topology examples:

- `two_components_componentwise_smooth_signals`
- `long_sparse_bridge_sharp_transition`
- `y_branch_contrast_topology_baseline`
- `y_branch_narrow_peak_leakage_stress`

The report machinery is extended with:

- component-wise fit status,
- number of fitted components,
- total fitted component size and fraction,
- aggregate component-wise RMSE and correlation,
- component-wise prediction panels,
- component-wise residual panels.

These additions let the report answer a new question:

> If the graph is disconnected and multiple components are meaningful, does the
> estimator behave correctly on each component rather than only on the largest
> one?
