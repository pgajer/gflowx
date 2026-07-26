#' Fit k-Nearest Neighbor Riemannian Graph Regression
#'
#' Constructs a 1-skeleton k-nearest neighbor complex and iteratively adapts
#' its Riemannian geometry to reflect response structure, enabling conditional
#' expectation estimation that respects both feature topology and response
#' boundaries.
#'
#' @param X Numeric matrix (n x d) of features, or sparse matrix (dgCMatrix).
#'   Each row represents an observation, each column a feature. Features
#'   should be appropriately scaled if they have different units.
#'
#' @param y Numeric vector of length n containing response values. Cannot
#'   contain NA or infinite values.
#'
#' @param y.vertices Indices of labeled vertices; only these contribute to the
#'     data-fit term; fitted field is returned on all vertices. NULL or an
#'     integer vector of vertex indices over which the conditional expectation
#'     of y will be estimated. If y.vertices = NULL, the conditional expectation
#'     of y will be estimted over all vertices.
#'
#' @param k Integer scalar giving the number of nearest neighbors. Must satisfy
#'   \eqn{2 \le k < n}. Larger k produces smoother fits but may oversmooth
#'   fine-scale features. Typical values are in the range 5 to 30. If NULL,
#'   selects k via cross-validation (NOT YET IMPLEMENTED). When
#'   \code{adj.list}/\code{weight.list} are supplied, graph construction is
#'   skipped. The current precomputed-graph API supports only
#'   \code{use.counting.measure = TRUE}; distance-reference initialization needs
#'   true kNN neighborhoods and distances rather than graph adjacency alone.
#'
#' @param adj.list Optional precomputed graph adjacency list (1-based indexing).
#'   If provided together with \code{weight.list}, the fit uses this graph
#'   directly and skips internal iKNN graph construction/pruning. Currently only
#'   supported with \code{use.counting.measure = TRUE}.
#'
#' @param weight.list Optional precomputed edge-length list parallel to
#'   \code{adj.list}. Must be provided together with \code{adj.list}.
#'
#' @param with.posterior Logical indicating whether to compute Bayesian posterior
#'   credible intervals for fitted values. When \code{TRUE}, the function generates
#'   Monte Carlo samples from the posterior distribution of the response surface
#'   and returns vertex-wise credible interval bounds and posterior standard
#'   deviations. This enables direct probabilistic statements about fitted values
#'   rather than relying on frequentist p-values. Posterior inference adds minimal
#'   computational overhead as it reuses the cached spectral decomposition. Default
#'   is \code{FALSE}.
#'
#' @param return.posterior.samples Set to TRUE to add posterior samples to the output list.
#'
#' @param credible.level Numeric value in (0, 1) specifying the posterior
#'   probability mass contained in the credible intervals. For example,
#'   \code{credible.level = 0.95} produces 95\% credible intervals such that
#'   P(y_lower < y_true < y_upper | data) = 0.95 at each vertex. Unlike
#'   frequentist confidence intervals, these admit direct probability
#'   interpretation: given the observed data and smoothness prior, there is
#'   95\% posterior probability that the true response lies within the interval.
#'   Typical values are 0.90, 0.95, or 0.99. Only used when
#'   \code{with.posterior = TRUE}. Default is 0.95.
#'
#' @param n.posterior.samples Integer specifying the number of Monte Carlo samples
#'   to draw from the posterior distribution for computing credible intervals and
#'   posterior standard deviations. Larger values provide more accurate quantile
#'   estimates but increase computation time linearly. Values between 1000-2000
#'   typically provide good accuracy. For n vertices and m eigenpairs, generating
#'   B samples requires O(Bnm) operations, which is negligible compared to the
#'   initial spectral decomposition for typical B values. Only used when
#'   \code{with.posterior = TRUE}. Must be at least 100. Default is 1000.
#'
#' @param posterior.seed Integer seed for the random number generator used in
#'   posterior sampling, enabling reproducible credible intervals across runs.
#'   If \code{NULL}, a random seed is generated automatically. Specify an integer
#'   value to ensure identical posterior samples and credible bounds in repeated
#'   analyses. Only used when \code{with.posterior = TRUE}. Default is \code{NULL}.
#'
#' @param pca.dim Positive integer or `NULL`. If not `NULL` and `ncol(X) >
#'     pca.dim`, PCA is used to reduce to at most `pca.dim` components.
#'
#' @param variance.explained Numeric in \eqn{(0,1]} or `NULL`. If not `NULL`,
#'     choose the smallest number of PCs whose cumulative variance explained
#'     exceeds this threshold, capped by `pca.dim`.
#'
#' @param max.iterations Integer scalar, positive. Maximum number of iterations.
#'   Default 50. Typical convergence occurs within 5-10 iterations for
#'   well-chosen parameters. If maximum iterations reached without convergence,
#'   a warning is issued.
#'
#' @param n.eigenpairs Integer scalar in the range 10 to n. Number of Laplacian
#'   eigenpairs to compute for spectral filtering. Default 200. Larger values
#'   capture more frequency components but increase computation time. Must be
#'   less than the number of observations n.
#'
#' @param filter.type Character string specifying the spectral filter type for
#'   response smoothing. The filter determines how high-frequency (rapid variation)
#'   components are attenuated relative to low-frequency (smooth) components.
#'   The smoothing parameter is selected automatically via generalized cross-validation.
#'   Available options:
#'   \describe{
#'     \item{\code{"heat_kernel"}}{Exponential decay \eqn{\exp(-\eta \lambda)}. General-purpose
#'       filter with smooth, progressive attenuation. Most commonly used and
#'       recommended as the default choice. Equivalent to heat diffusion for
#'       time \eqn{\eta}.}
#'     \item{\code{"tikhonov"}}{Rational decay \eqn{1/(1+\eta \lambda)}. Corresponds to
#'       Tikhonov regularization. Gentler than heat kernel, preserves more
#'       mid-frequency detail and linear trends. Good when moderate smoothing
#'       is sufficient.}
#'     \item{\code{"cubic_spline"}}{Spline smoothness \eqn{1/(1+\eta \lambda^2)}. Minimizes
#'       second derivatives of the fitted function. Produces very smooth,
#'       spline-like results. Excellent when the response should vary gradually
#'       and continuously.}
#'     \item{\code{"gaussian"}}{Super-exponential decay \eqn{\exp(-\eta \lambda^2)}. The most
#'       aggressive smoothing among exponential filters. Produces extremely
#'       smooth results with minimal oscillations. Use when maximum smoothness
#'       is desired and fine-scale features are noise.}
#'     \item{\code{"exponential"}}{Intermediate decay \eqn{\exp(-\eta \sqrt{\lambda})}. Less aggressive
#'       than heat kernel. Maintains more detail in mid-frequency range while
#'       still providing meaningful smoothing. Useful for preserving moderate-scale
#'       features while removing high-frequency noise.}
#'     \item{\code{"butterworth"}}{Smooth cutoff \eqn{1/(1+(\lambda/\eta)^4)}. Fourth-order
#'       rational filter providing clear frequency separation with reduced
#'       ringing artifacts compared to ideal low-pass filters. Good balance
#'       between sharp cutoff and smooth transition.}
#'   }
#'   Default: \code{"heat_kernel"}
#'
#' @param t.diffusion Diffusion time parameter (non-negative). Controls how much
#'     vertex density evolves in each iteration. If 0 (default), automatically
#'     selected as \eqn{0.5/\lambda_2} where \eqn{\lambda_2} is the spectral
#'     gap. Larger values produce more aggressive density updates; smaller
#'     values are more conservative. Typical range is \eqn{0.1/\lambda_2} to
#'     \eqn{1.0/\lambda_2}. Set to 0 for automatic selection (recommended).
#'
#' @param density.uniform.pull Restoring force strength (non-negative). Controls
#'     how strongly vertex densities are pulled back toward uniform distribution
#'     during diffusion, preventing excessive concentration. If 0 (default),
#'     automatically selected as 0.1/t.diffusion to maintain a 10:1 ratio of
#'     diffusion to restoration. Larger values prevent concentration but may
#'     over-smooth; smaller values allow more geometric structure but risk
#'     collapse. Set to 0 for automatic selection (recommended). Former
#'     beta.damping.
#'
#' @param response.penalty.exp Response penalty exponent that controls the rate
#'     of edge weight reduction as response variation increases, via the penalty
#'     function \eqn{\Gamma(\Delta) = (1 + \Delta^2/\sigma^2)^{-\gamma}}. Larger
#'     values create steeper penalties (sharper adaptation to response changes);
#'     smaller values create gentler penalties. For smooth responses use
#'     0.5-0.75; for discontinuous responses use 1.5-2.0. When set to 0,
#'     disables response-coherence modulation entirely (geometry remains fixed
#'     based on feature space only). When negative (e.g., -1), triggers
#'     automatic selection via first-iteration GCV evaluation over a default
#'     grid c(0.05, 0.1, 0.2, ..., 0.9, 1.0, 1.2, ..., 2.0). Automatic selection
#'     efficiently determines the optimal value by evaluating GCV after a single
#'     iteration for each candidate, then uses the selected value for the full
#'     iterative refinement. Recommended when the response smoothness
#'     characteristics are unknown or when exploring new datasets. When
#'     positive, must be less than 2.0 and Default 0.
#'
#' @param t.scale.factor Diffusion scale factor (positive numeric). Controls the
#'     dimensionless product \eqn{t \cdot \lambda_2} when \code{t.diffusion = 0}
#'     (automatic selection), where \eqn{\lambda_2} is the spectral gap of the
#'     graph Laplacian. This parameter determines how aggressively density
#'     diffuses per iteration relative to the graph's intrinsic geometric scale.
#'     Conservative values \eqn{(0.1-0.2)} produce slow but stable density
#'     evolution, requiring more iterations but reducing risk of oscillation.
#'     Moderate values \eqn{(0.4-0.6)} balance speed and stability for typical
#'     applications. Aggressive values \eqn{(0.8-1.0)} accelerate convergence but may
#'     cause oscillatory behavior when geometric-response feedback is strong.
#'     The scale factor is ignored when \code{t.diffusion > 0} (manual
#'     specification). Typical range: \eqn{[0.1, 1.0]}. Default: 0.5
#'
#' @param beta.coef.factor Damping coefficient factor (positive numeric).
#'     Controls the product \eqn{\beta \cdot t} when \code{density.uniform.pull
#'     = 0} (automatic selection). This parameter determines the strength of the
#'     restoring force toward uniform distribution during each iteration step.
#'     The damping prevents excessive concentration of density in small regions
#'     while allowing genuine accumulation in well-connected areas. Light
#'     damping (0.05) permits strong concentration and is appropriate when the
#'     data distribution is expected to be highly non-uniform. Moderate damping
#'     (0.1) provides balanced geometric adaptation for most applications.
#'     Strong damping (0.2-0.3) prevents concentration aggressively and is
#'     appropriate when maintaining dispersed density is important or when the
#'     response varies smoothly without sharp boundaries. The coefficient factor
#'     is ignored when \code{density.uniform.pull > 0} (manual specification). A
#'     coefficient of 0.1 means the restoring force contributes approximately
#'     10\% as much as the identity term in the system matrix. Typical range:
#'     \eqn{[0.05, 0.3]}. Default: 0.1
#'
#' @param t.update Character scalar. Either \code{"fixed"} (default) or
#'     \code{"per_iteration"}. If \code{"per_iteration"}, the diffusion time
#'     \code{t.diffusion} is updated each iteration to keep the diffusion scale
#'     \eqn{s = t \lambda_2} approximately constant using the current
#'     \eqn{\lambda_2}. Only applied when \code{t.diffusion <= 0} (auto-selected
#'     t).
#'
#' @param t.update.max.mult Numeric scalar \eqn{\ge 1}. Maximum multiplicative
#'     change allowed for \code{t.diffusion} in a single iteration when
#'     \code{t.update="per_iteration"} (default 1.25).
#'
#' @param use.counting.measure Logical scalar. If TRUE, uses uniform vertex
#'     weights (counting measure). If FALSE, uses distance-based weights
#'     inversely proportional to local k-NN density: \eqn{w(x) = (\epsilon +
#'     d_k(x))^{-\alpha}}. Distance-based weights are useful when sampling
#'     density varies across the feature space. Must be TRUE when
#'     \code{adj.list}/\code{weight.list} are supplied. Default: TRUE.
#'
#' @param density.normalization Numeric scalar, non-negative. Specifies target
#'   sum for normalized vertex densities. If 0 (default), densities are
#'   normalized to sum to n (giving average vertex density = 1). If positive,
#'   densities sum to the specified value. Most users should use the default.
#'   Set to 1.0 for probability interpretation, or to a large value (e.g., 1e6)
#'   for very large graphs to avoid numerical underflow.
#'
#' @param density.alpha Exponent alpha in \eqn{[1,2]} for distance-based
#'     reference measure formula \eqn{\mu(x) = (\epsilon + d_k(x))^{-\alpha}}.
#'     Larger values (near 2) create stronger density-adaptive weighting,
#'     concentrating measure on densely sampled regions. Smaller values (near 1)
#'     produce more uniform weighting. Only used when use_counting_measure =
#'     FALSE. Default: 1.5.
#'
#' @param density.epsilon Regularization parameter epsilon > 0 in reference
#'     measure formula \eqn{\mu(x) = (\epsilon + d_k(x))^{-\alpha}}. Prevents
#'     numerical issues when nearest neighbor distances are very small. Should
#'     be small relative to typical k-NN distances but large enough for
#'     stability. Only used when use_counting_measure = FALSE. Default: 1e-10.
#'
#' @param clamp.dk Logical scalar. If \code{TRUE}, clamp the k-th neighbor distances
#'   \eqn{d_k} used to construct the distance-based reference measure before converting
#'   them to weights. Clamping uses median-factor bounds controlled by
#'   \code{dk.clamp.median.factor}. If \code{FALSE}, no clamping is applied (bounds
#'   may still be computed for diagnostics).
#'
#' @param dk.clamp.median.factor Numeric scalar \eqn{> 1}. Median-factor used to define
#'   clamping bounds when \code{clamp.dk = TRUE}:
#'   \eqn{[ \mathrm{median}(d_k) / f,\ \mathrm{median}(d_k) \cdot f ]}.
#'   Larger values reduce clamping; smaller values increase clamping.
#'
#' @param target.weight.ratio Numeric scalar. Upper bound on the dynamic range of the
#'   reference measure weights, enforcing \eqn{\max(w)/\min(w) \le R} via compression
#'   when needed. Larger values reduce compression; smaller values increase compression.
#'
#' @param pathological.ratio.threshold Numeric scalar. Threshold for triggering a more
#'   robust compression path when the initial weight ratio is extremely large or non-finite.
#'   When \eqn{\max(w)/\min(w)} exceeds this value, an additional log-compression step
#'   may be applied before enforcing \code{target.weight.ratio}.
#' 
#' @param compute.extremality Logical; compute Riemannian extremality scores
#'   during fitting? When TRUE, extremality scores are computed for all vertices
#'   at the optimal GCV iteration using the full Riemannian metric via the
#'   coboundary operator. When FALSE (default), extremality analysis is skipped,
#'   reducing computation time and memory usage. Set to FALSE for applications
#'   like k-selection where only fitted values and GCV scores are needed. Set to
#'   TRUE for comprehensive extrema detection and characterization. See Details
#'   for extremality score computation. Default: TRUE.
#'
#' @param p.threshold Numeric in \eqn{[0, 1]}; extremality threshold for hop radius
#'   computation. Only used when compute.extremality = TRUE. When p.threshold = 0,
#'   only extremality scores are computed without hop radii (faster). When
#'   p.threshold > 0, hop-extremality radii are computed for all vertices with
#'   |extremality| >= p.threshold, measuring spatial persistence of extremal
#'   character. Non-candidates are assigned NA. Typical values are 0.85-0.95.
#'   Higher thresholds select fewer, stronger extrema. See compute.gextrema.nbhds()
#'   for interpretation of hop radii. Default: 0.95.
#'
#' @param max.hop Integer; maximum hop distance for extremality radius computation.
#'   Only used when compute.extremality = TRUE and p.threshold > 0. Controls how
#'   far the BFS traversal extends when computing hop-extremality radii. Vertices
#'   maintaining |extremality| >= p.threshold beyond max.hop are assigned radius
#'   = max.hop (not infinity, to distinguish from global extrema). Larger values
#'   capture more persistent extrema but increase computation. Typical range is
#'   10-50 depending on graph diameter and application. For dense graphs or when
#'   computational cost is a concern, use smaller values (10-20). For sparse graphs
#'   or when long-range persistence is important, use larger values (30-50).
#'   Default: 30.
#'
#' @param max.ratio.threshold Numeric in the range 0 to 0.2.
#'     Geometric pruning removes an edge (i,j) when there exists an
#'     alternative path between i and j whose path/edge length ratio
#'     minus 1.0 is less than this threshold. This is a deviation
#'     threshold in the range 0 to 0.2. Internally we compare the
#'     path-to-edge ratio R to \eqn{1 + \delta}. A value of 0 means
#'     no geometric pruning (ratio threshold 1.0). Default 0.1.
#'
#' @param apply.geometric.pruning Logical; if \code{TRUE} (default), apply
#'     geometric pruning using \code{max.ratio.threshold}. If \code{FALSE},
#'     skip geometric pruning entirely in the rdgraph fit path. This provides
#'     an explicit off-switch and avoids relying on threshold conventions.
#'
#' @param path.edge.ratio.percentile Numeric in \eqn{[0,1]}. Only edges with
#'     length above this percentile are considered for geometric pruning.
#'
#' @param threshold.percentile Numeric in \eqn{[0, 0.5]}. Percentile threshold for
#'     quantile-based edge length pruning. Default is 0 (no quantile pruning).
#'     When greater than 0, edges in the top (1 - threshold.percentile) quantile
#'     by length are removed if their removal preserves graph connectivity.
#'     For example, threshold.percentile = 0.9 removes edges in the top 10\% by length.
#'     This pruning is applied after geometric pruning and targets unusually long edges
#'     based on absolute edge lengths rather than path-to-edge ratios.
#'
#' @param knn.cache.path Optional character scalar path to a binary kNN cache file.
#'   Used only when \code{knn.cache.mode != "none"}.
#' @param knn.cache.mode Character cache mode:
#'   \itemize{
#'     \item \code{"none"}: always compute kNN; do not read/write cache.
#'     \item \code{"read"}: read kNN from cache only; if cache is missing/invalid, error.
#'     \item \code{"write"}: always compute kNN and atomically write/overwrite cache.
#'     \item \code{"readwrite"}: read valid cache when available; if cache file is missing,
#'       compute and write cache; if cache exists but is invalid, error.
#'   }
#'
#' @param epsilon.y Numeric scalar, positive. Relative convergence threshold
#'   for response. Iteration stops when the relative change in fitted values
#'   \eqn{||\hat{y}^{(\ell)} - \hat{y}^{(\ell-1)}|| / ||\hat{y}^{(\ell-1)}||}
#'   falls below this value. Default 1e-4. Smaller values require more iterations
#'   but give more precise convergence.
#'
#' @param epsilon.rho Numeric scalar, positive. Relative convergence threshold
#'   for densities. Iteration stops when
#'   \eqn{\max_p ||\rho_p^{(\ell)} - \rho_p^{(\ell-1)}|| / ||\rho_p^{(\ell-1)}||}
#'   falls below this value across all dimensions. Default 1e-4.
#'
#' @param store.projected.X Logical scalar. If \code{TRUE}, store the exact
#'     feature matrix used for fitting as \code{attr(fit, "X.used")}. When PCA
#'     is applied, also stores the PCA-projected matrix as
#'     \code{attr(fit, "pca")$X_projected}. Default is \code{FALSE} to avoid
#'     increasing object size.
#'
#' @param dense.fallback Character scalar controlling dense direct-solver fallback
#'   in the diffusion linear solve and spectral decomposition. \code{"auto"}
#'   (default) preserves current behavior (direct solve for small systems and
#'   dense fallback after iterative-solver failure on larger systems).
#'   \code{"never"} disables dense fallback and raises a detailed iterative-solver
#'   error instead of attempting sparse-to-dense allocation. \code{"always"}
#'   forces dense direct solve (debugging mode; can require \eqn{O(n^2)} memory).
#'
#' @param triangle.policy Character scalar controlling triangle construction in
#'   simplicial initialization. \code{"auto"} (default) builds triangles only
#'   when they can affect response-coherence modulation, \code{"never"} disables
#'   triangle construction, and \code{"always"} forces legacy behavior.
#'
#' @param verbose.level Integer in \code{c(0, 1, 2, 3)} controlling reporting.
#' @param ... Additional advanced options forwarded to lower-level fitting routines.
#'
#' @return An object of class \code{c("knn.riem.fit", "list")}, which is
#'   a list containing:
#'   \describe{
#'     \item{fitted.values}{Numeric vector of fitted values selected from the
#'       iteration with minimum GCV.}
#'     \item{residuals}{Numeric vector of residuals.}
#'     \item{optimal.iteration}{Integer indicating which iteration's fitted
#'       values were selected (1-indexed).}
#'     \item{graph}{List with graph structure information.}
#'     \item{iteration}{List with convergence information.}
#'     \item{parameters}{List of input parameters.}
#'     \item{y}{Original response values.}
#'     \item{gcv}{List with GCV information for each iteration.}
#'     \item{density}{List with density history.}
#'     \item{extremality.scores}{(Only if \code{compute.extremality = TRUE})
#'       List containing extremality diagnostics; includes \code{scores} and,
#'       when \code{p.threshold > 0}, \code{hop_extremality_radii}.}
#'   }
#'   When \code{store.projected.X = TRUE}, additional attributes are attached:
#'   \describe{
#'     \item{attr(fit, "X.used")}{Exact matrix passed to the C++ fitter
#'       (PCA-projected if PCA was applied, otherwise validated input \code{X}).}
#'     \item{attr(fit, "pca")$X_projected}{PCA-projected matrix (only when
#'       PCA is applied).}
#'   }
#'
#' @details
#' \strong{Algorithm Overview}
#'
#' The algorithm proceeds in two main phases:
#'
#' \emph{Phase I: Initialization}
#' \enumerate{
#'   \item Constructs 1-skeleton k-NN complex from features using ANN library
#'   \item Computes initial vertex and edge densities from neighborhood overlaps
#'   \item Builds mass matrices \eqn{M_0} (diagonal) and \eqn{M_1} (full, via
#'         Gram determinants)
#'   \item Assembles Hodge Laplacian \eqn{L_0 = B_1 M_1^{-1} B_1^T M_0}
#'   \item Performs initial spectral filtering: \eqn{\hat{y}^{(0)} = V F_\eta(\Lambda) V^T y}
#' }
#'
#' \emph{Phase II: Iterative Refinement}
#'
#' Repeats until convergence or max.iterations reached:
#' \enumerate{
#'   \item \strong{Density diffusion}: Evolve vertex density via damped heat
#'         equation \eqn{(I + t(L_0 + \beta I))\rho_0 = \rho_{0,prev} + t\beta u}
#'   \item \strong{Edge density update}: Recompute edge densities from evolved
#'         vertex masses via neighborhood intersections
#'   \item \strong{Response-coherence modulation}: Reduce edge masses where
#'         response varies rapidly: \eqn{\rho_1([i,j]) \leftarrow \rho_1([i,j]) \cdot (1 + \Delta_{ij}^2/\sigma^2)^{-\gamma}}
#'   \item \strong{Metric update}: Rebuild \eqn{M_0}, \eqn{M_1} from modified densities
#'   \item \strong{Laplacian reassembly}: Recompute \eqn{L_0} with updated metric
#'   \item \strong{Response smoothing}: Re-filter response using new geometry
#'   \item \strong{Convergence check}: Test if response and densities stabilized
#' }
#'
#' \strong{Extremality Analysis During Fitting}
#'
#' When compute.extremality = TRUE, the function computes Riemannian extremality
#' scores at the optimal GCV iteration, optionally with spatial persistence
#' measures (hop-extremality radii). This integrated approach is more efficient
#' than post-processing because:
#' \itemize{
#'   \item The internal graph structures are already in memory
#'   \item Extremality scores use the same fitted values selected by GCV
#'   \item Results are guaranteed to be consistent with the returned fit
#' }
#'
#' The extremality score at vertex v is computed via the coboundary operator:
#' \deqn{\text{extr}(v; \hat{y}, g) = -\frac{\sum_{e \ni v} \rho_1(e) \cdot [\partial_1^*(\hat{y})](e)}{\sum_{e \ni v} \rho_1(e) \cdot |[\partial_1^*(\hat{y})](e)|}}
#'
#' where \eqn{\partial_1^*(\hat{y})} is the coboundary computed by solving
#' \eqn{M_1 x = B_1^T M_0 \hat{y}} with the full (potentially non-diagonal)
#' edge mass matrix \eqn{M_1}. Values near +1 indicate local maxima, values
#' near -1 indicate local minima, and values near 0 indicate saddle points or
#' plateaus.
#'
#' When p.threshold > 0, hop-extremality radii are computed for vertices with
#' |\code{extremality}| >= p.threshold. For a vertex v, the hop-extremality radius is:
#' \deqn{r_\tau(v) = \max\{h \geq 1 : |\text{extr}^{(h)}(v; \hat{y})| \geq \tau\}}
#'
#' where \eqn{\text{extr}^{(h)}(v; \hat{y})} is the extremality computed over the
#' h-hop neighborhood. Large hop radii indicate persistent, stable extrema while
#' small hop radii suggest local noise fluctuations.
#'
#' \strong{Computational Cost}
#'
#' Extremality score computation: \eqn{O(m^{1.5})} for first call (includes solving
#' \eqn{M_1 x = B_1^T M_0 \hat{y}} via Cholesky factorization or iterative CG),
#' where m is the number of edges. For typical kNN graphs with m = O(kn), this
#' is \eqn{O((kn)^{1.5})}.
#'
#' Hop radius computation: \eqn{O(k d h)} where \eqn{k} is the number of candidates
#' (typically \eqn{O(\sqrt{n})}), \eqn{d} is average degree, and \eqn{h} is average hop radius at
#' termination (typically << max.hop due to early stopping).
#'
#' For large-scale problems (n > 10,000), consider:
#' \itemize{
#'   \item Set compute.extremality = FALSE during k-selection
#'   \item Use p.threshold = 0 to skip hop radius computation
#'   \item Increase p.threshold to reduce candidate count (e.g., 0.95 instead of 0.90)
#'   \item Decrease max.hop if only local persistence matters
#' }
#'
#' \strong{Using the Extremality Results}
#'
#' When compute.extremality = TRUE, the returned fit object contains an
#' additional component \code{extremality.scores}, which is a list with:
#' \itemize{
#'   \item \code{scores}: Numeric vector of extremality scores for all vertices
#'   \item \code{hop_extremality_radii}: Numeric vector of hop radii (only if p.threshold > 0)
#'         Values: hop radius for candidates, NA for non-candidates, -1 for global extrema
#' }
#'
#' Access these via:
#' \preformatted{
#' fit <- fit.rdgraph.regression(X, y, k = 50, compute.extremality = TRUE)
#' extremality_scores <- fit$extremality.scores$scores
#' hop_radii <- fit$extremality.scores$hop_extremality_radii  # if p.threshold > 0
#'
#' # Identify strong maxima
#' strong_maxima <- which(extremality_scores >= 0.90)
#' persistent_maxima <- strong_maxima[hop_radii[strong_maxima] >= 5]
#' }
#'
#' \strong{Geometric Interpretation}
#'
#' The method creates a Riemannian metric on the k-NN graph where:
#' \itemize{
#'   \item Edges within response-coherent regions have high mass (short distance)
#'   \item Edges crossing response boundaries have low mass (long distance)
#'   \item Diffusion is rapid within coherent regions, slow across boundaries
#'   \item Spectral filtering respects this geometry, producing smooth fits
#'         within regions while preserving discontinuities at boundaries
#' }
#'
#' \strong{Parameter Selection Guidelines}
#'
#' For most applications, default parameters work well. Consider tuning:
#' \itemize{
#'   \item \strong{k}: Primary tuning parameter. Start with \eqn{k = \lceil\log(n)\rceil}
#'         and adjust based on cross-validation. Larger k smooths more.
#'   \item \strong{response.penalty.exp}: Increase (up to 2.0) for sharper response
#'         boundaries; decrease (down to 0.5) for gentler transitions.
#'   \item \strong{filter.type}: "heat_kernel" (default) for general use;
#'         "tikhonov" for better linear trend preservation; "cubic_spline"
#'         for smoother results on chain-like graphs.
#' }
#'
#' \strong{Automatic Parameter Selection}
#'
#' Several key parameters are automatically selected if not specified:
#' \itemize{
#'   \item \strong{t.diffusion}: Set to \eqn{0.5/\lambda_2} where \eqn{\lambda_2}
#'         is the spectral gap of the initial Laplacian
#'   \item \strong{density.uniform.pull}: Set to \eqn{0.1/t} to maintain 10:1
#'         diffusion to damping ratio
#'   \item \strong{response.penalty.exp}: When set to negative value (e.g., -1),
#'         automatically selected via first-iteration GCV evaluation across a
#'         grid of candidate values. This efficiently determines the optimal
#'         response-coherence strength by evaluating GCV after one iteration
#'         for each candidate, avoiding the computational cost of full convergence
#'         for each value
#' }
#'
#' The first-iteration GCV approach for gamma selection provides 5-10x speedup
#' compared to full cross-validation while maintaining comparable selection quality.
#'
#' @section Computational Complexity:
#' \itemize{
#'   \item k-NN search: \eqn{O(n \log n)} using kd-trees (ANN library)
#'   \item Edge construction: \eqn{O(n k^2)} for intersection computations
#'   \item Per iteration: \eqn{O(n k^2 + m n_{eig})} where \eqn{m \approx n k}
#'         is edge count
#'   \item Total: \eqn{O(n k^2 + iter \times m n_{eig})}, typically manageable
#'         for \eqn{n \le 10,000}
#' }
#'
#' @section Memory Usage:
#' \itemize{
#'   \item Feature matrix: \eqn{O(n d)}
#'   \item Complex storage: \eqn{O(n k)} for vertices and edges
#'   \item Mass matrices: \eqn{O(n + m^2)} but \eqn{M_1} is sparse with
#'         \eqn{O(m k^2)} nonzeros
#'   \item Eigenpairs: \eqn{O(n \times n_{eigenpairs})}
#'   \item Peak usage: dominated by Laplacian operations, typically < 1GB
#'         for n = 5000
#' }
#'
#' @note \strong{Graph Connectivity Requirement:}
#'   The k-nearest neighbor graph must be fully connected (form a single
#'   connected component). If the graph is disconnected, the function will
#'   terminate with an error. To explore connectivity structure before fitting:
#'   \itemize{
#'     \item Use \code{dgraphs::create.iknn.graphs()} to examine connectivity across
#'           multiple k values
#'     \item Use \code{dgraphs::create.single.iknn.graph()} to inspect the graph
#'           structure for a specific k
#'     \item Examine connected components in the returned graph objects
#'   }
#'
#'   Common solutions for disconnected graphs:
#'   \itemize{
#'     \item Increase k to add more edges between vertices
#'     \item Remove isolated outlier observations
#'     \item Verify that your distance metric is appropriate for the data
#'   }
#'
#' @references
#' Lim, L. H. (2020). Hodge Laplacians on graphs. \emph{SIAM Review}, 62(3), 685-715.
#' (Hodge theory on simplicial complexes)
#'
#' Belkin, M., & Niyogi, P. (2003). Laplacian eigenmaps for dimensionality
#' reduction and data representation. \emph{Neural Computation}, 15(6), 1373-1396.
#' (Graph-based geometric methods)
#'
#' @examples
#' \dontrun{
#' # Example 1: Simple nonlinear regression
#' set.seed(123)
#' n <- 200
#' X <- matrix(rnorm(n * 3), n, 3)
#' y <- sin(2 * pi * X[,1]) + 0.5 * X[,2]^2 + rnorm(n, sd = 0.2)
#'
#' fit <- fit.rdgraph.regression(X, y, k = 15)
#' y.hat <- fitted(fit)
#' plot(y, y.hat, asp = 1, main = "Fitted vs Actual")
#' abline(0, 1, col = "red")
#'
#' # Example 2: Regression with sharp boundaries
#' X <- matrix(runif(n * 2, -1, 1), n, 2)
#' y <- ifelse(X[,1]^2 + X[,2]^2 < 0.5, 1, 0) + rnorm(n, sd = 0.1)
#'
#' # Use larger gamma for sharper boundary detection
#' fit <- fit.rdgraph.regression(X, y, k = 20, response.penalty.exp = 1.5)
#' y.hat <- fitted(fit)
#'
#' # Example 3: Using sparse matrix input
#' library(Matrix)
#' X.sparse <- Matrix(X, sparse = TRUE)
#' fit <- fit.rdgraph.regression(X.sparse, y, k = 15)
#'
#' # Example 4: Automatic response.penalty.exp selection
#' # When optimal gamma is unknown, use automatic selection
#' set.seed(456)
#' X <- matrix(rnorm(n * 2), n, 2)
#' # Mix of smooth and discontinuous regions
#' y <- ifelse(X[,1] > 0, sin(3*X[,2]), 2*X[,2]^2) + rnorm(n, sd = 0.1)
#'
#' # Let algorithm select optimal gamma
#' fit.auto <- fit.rdgraph.regression(X, y, k = 20,
#'                                            response.penalty.exp = -1)
#' cat("Selected gamma:", fit.auto$parameters$response.penalty.exp, "\n")
#' cat("Optimal iteration:", fit.auto$optimal.iteration, "\n")
#'
#' # Example 5: Computing extremality scores during fitting
#' set.seed(789)
#' n <- 300
#' X <- matrix(rnorm(n * 2), n, 2)
#' # Function with clear local maxima and minima
#' y <- sin(2 * pi * X[,1]) * cos(2 * pi * X[,2]) + rnorm(n, sd = 0.1)
#'
#' # Fit with extremality analysis
#' fit_with_extr <- fit.rdgraph.regression(X, y, k = 20,
#'                                         compute.extremality = TRUE,
#'                                         p.threshold = 0.90,
#'                                         max.hop = 20)
#'
#' # Extract extremality scores
#' extr_scores <- fit_with_extr$extremality
#' hop_radii <- fit_with_extr$extremality_extremality_radii
#'
#' # Identify persistent extrema
#' strong_maxima <- which(extr_scores >= 0.90 & !is.na(hop_radii) & hop_radii >= 5)
#' strong_minima <- which(extr_scores <= -0.90 & !is.na(hop_radii) & hop_radii >= 5)
#'
#' cat("Found", length(strong_maxima), "persistent maxima\n")
#' cat("Found", length(strong_minima), "persistent minima\n")
#'
#' # Visualize extremality scores
#' plot(X[,1], X[,2], col = rainbow(100)[cut(extr_scores, 100)],
#'      pch = 19, main = "Extremality Scores")
#' points(X[strong_maxima, 1], X[strong_maxima, 2], pch = 24, cex = 1.5, col = "red")
#' points(X[strong_minima, 1], X[strong_minima, 2], pch = 25, cex = 1.5, col = "blue")
#'
#' # Example 6: K-selection without extremality (faster)
#' # For large-scale k-selection, skip extremality computation
#' k_values <- c(10, 15, 20, 25, 30)
#' gcv_scores <- numeric(length(k_values))
#'
#' for (i in seq_along(k_values)) {
#'     fit <- fit.rdgraph.regression(X, y, k = k_values[i],
#'                                   compute.extremality = FALSE,  # Skip for speed
#'                                   verbose.level = 1)
#'     gcv_scores[i] <- fit$gcv$gcv.optimal[fit$optimal.iteration]
#' }
#'
#' optimal_k <- k_values[which.min(gcv_scores)]
#' cat("Optimal k:", optimal_k, "\n")
#'
#' # Now fit with optimal k and compute extremality
#' final_fit <- fit.rdgraph.regression(X, y, k = optimal_k,
#'                                     compute.extremality = TRUE,
#'                                     p.threshold = 0.90)
#' }
#'
#' @seealso
#' \code{\link{fitted}} for extracting fitted values,
#' \code{\link{residuals}} for extracting residuals,
#' \code{\link{summary}} for model summaries
#'
#' @export
fit.rdgraph.regression <- function(
    X,
    y,
    y.vertices = NULL,
    k,
    adj.list = NULL,
    weight.list = NULL,
    with.posterior = FALSE,
    return.posterior.samples = FALSE,
    credible.level = 0.95,
    n.posterior.samples = 1000,
    posterior.seed = NULL,
    pca.dim = 50,
    variance.explained = 0.95,
    max.iterations = 10,
    n.eigenpairs = 10,
    filter.type = "heat_kernel",
    t.diffusion = 0,
    density.uniform.pull = 0,
    response.penalty.exp = 0,
    t.scale.factor = 0.5,
    beta.coef.factor = 0.1,
    t.update = c("fixed", "per_iteration"),
    t.update.max.mult = 1.25,
    use.counting.measure = TRUE,
    density.normalization = 0,
    density.alpha = 1.5,
    density.epsilon = 1e-10,
    clamp.dk = FALSE,
    dk.clamp.median.factor = 10.0,
    target.weight.ratio = 1000,
    pathological.ratio.threshold = 1e12,
    compute.extremality = FALSE,
    p.threshold = 0.95,
    max.hop = 30L,
    max.ratio.threshold = 0.1,
    apply.geometric.pruning = TRUE,
    path.edge.ratio.percentile = 0.5,
    threshold.percentile = 0,
    knn.cache.path = NULL,
    knn.cache.mode = c("none", "read", "write", "readwrite"),
    epsilon.y = 1e-4,
    epsilon.rho = 1e-4,
    dense.fallback = c("auto", "never", "always"),
    triangle.policy = c("auto", "never", "always"),
    store.projected.X = FALSE,
    verbose.level = 1,
    ...
    ) {

    ## ==================== Dots / Deprecated Arguments ====================

    dots <- list(...)

    if ("verbose" %in% names(dots)) {
        stop("fit.rdgraph.regression() does not have a logical 'verbose' argument. Use integer 'verbose.level'.")
    }

    if (length(dots) > 0) {
        stop(sprintf("fit.rdgraph.regression(): unused argument(s): %s",
                     paste(names(dots), collapse = ", ")))
    }

    ## ==================== Feature Matrix Validation ====================

    ## Check X type
    is.dense <- is.matrix(X)
    is.sparse <- inherits(X, "dgCMatrix")

    if (!is.dense && !is.sparse) {
        stop("X must be a numeric matrix or sparse matrix (dgCMatrix from Matrix package)")
    }

    ## For dense matrix, check numeric
    if (is.dense && !is.numeric(X)) {
        stop("X must be numeric")
    }

    ## Extract dimensions
    n <- nrow(X)
    d <- ncol(X)

    if (length(y) != n) {
        stop("Length mismatch: length(y) != nrow(X).")
    }

    ## Check dimensions
    if (is.null(n) || is.null(d) || n < 1 || d < 1) {
        stop("X must have valid dimensions (positive rows and columns)")
    }

    if (n < 10) {
        stop(sprintf("X must have at least 10 observations (got n=%d)", n))
    }

    ## Check for NA/Inf in X
    if (is.dense) {
        if (anyNA(X)) {
            stop("X cannot contain NA values")
        }
        if (any(!is.finite(X))) {
            stop("X cannot contain infinite values")
        }
        ## Ensure storage mode is double
        storage.mode(X) <- "double"
    } else {
        ## For sparse matrix, check x slot
        if (anyNA(X@x)) {
            stop("X cannot contain NA values")
        }
        if (any(!is.finite(X@x))) {
            stop("X cannot contain infinite values")
        }
    }

    ## ==================== Response Vector Validation ====================

    ## Semi-supervised option:
    ## - if y.vertices is NULL: y must be length n (current behavior)
    ## - if y.vertices is provided: y over unlabeled entries will be set to 0 and ignored by mask

    if (!is.numeric(y)) {
        stop("y must be numeric")
    }

    if (is.null(y.vertices)) {

        if (length(y) != n) {
            stop(sprintf("Length of y (%d) must equal number of rows in X (%d)",
                         length(y), n))
        }

        if (anyNA(y)) {
            stop("y cannot contain NA values")
        }
        if (any(!is.finite(y))) {
            stop("y cannot contain infinite values")
        }

        y <- as.vector(y)

    } else {

        ## Validate y.vertices
        if (!is.numeric(y.vertices) && !is.integer(y.vertices)) {
            stop("y.vertices must be an integer (or numeric coercible to integer) vector")
        }
        y.vertices <- as.integer(y.vertices)

        if (length(y.vertices) < 1L) {
            stop("y.vertices must have positive length")
        }
        if (anyNA(y.vertices)) {
            stop("y.vertices cannot contain NA")
        }
        if (any(y.vertices < 1L | y.vertices > n)) {
            stop("y.vertices must be in 1..n (vertex indices are 1-based)")
        }
        if (anyDuplicated(y.vertices)) {
            stop("y.vertices cannot contain duplicates")
        }

        if (anyNA(y)) {
            stop("y cannot contain NA values (provide y only for labeled vertices)")
        }
        if (any(!is.finite(y))) {
            stop("y cannot contain infinite values")
        }

        y <- as.vector(y)

        ## Embed into full vector for C++ (unlabeled entries are ignored there)
        y.full <- numeric(n)
        y.full[y.vertices] <- y[y.vertices]
        y <- y.full
    }

    ## ==================== Parameter k Validation ====================

    if (!is.numeric(k) || length(k) != 1) {
        stop("k must be a single numeric value")
    }

    ## Convert to integer
    k <- as.integer(k)

    ## Check range
    if (k < 2) {
        stop(sprintf("k must be at least 2 (got k=%d)", k))
    }

    if (k >= n) {
        stop(sprintf("k must be less than n (got k=%d, n=%d)", k, n))
    }

    ## ==================== Optional Precomputed Graph Validation ====================

    has.precomputed.graph <- (!is.null(adj.list) || !is.null(weight.list))
    adj.list.0based <- NULL
    weight.list.cpp <- NULL

    if (has.precomputed.graph &&
        is.logical(use.counting.measure) &&
        length(use.counting.measure) == 1L &&
        !is.na(use.counting.measure) &&
        !use.counting.measure) {
        stop(
            "fit.rdgraph.regression(): precomputed adj.list/weight.list graphs ",
            "currently require use.counting.measure = TRUE. Distance-reference ",
            "initialization requires true kNN neighborhoods and distances; cached ",
            "graph adjacency alone is not a valid substitute."
        )
    }

    if (has.precomputed.graph) {
        if (is.null(adj.list) || is.null(weight.list)) {
            stop("adj.list and weight.list must be provided together (or both NULL).")
        }
        if (!is.list(adj.list)) {
            stop("adj.list must be a list when provided.")
        }
        if (!is.list(weight.list)) {
            stop("weight.list must be a list when provided.")
        }
        if (length(adj.list) != n || length(weight.list) != n) {
            stop(sprintf("adj.list/weight.list length mismatch: expected %d entries (nrow(X)).", n))
        }

        adj.list.norm <- vector("list", n)
        weight.list.norm <- vector("list", n)
        tol <- 1e-10

        for (i in seq_len(n)) {
            nbrs <- adj.list[[i]]
            wts <- weight.list[[i]]

            if (!is.numeric(nbrs) && !is.integer(nbrs)) {
                stop(sprintf("adj.list[[%d]] must be numeric/integer.", i))
            }
            if (!is.numeric(wts)) {
                stop(sprintf("weight.list[[%d]] must be numeric.", i))
            }

            nbrs <- as.integer(nbrs)
            wts <- as.double(wts)

            if (length(nbrs) != length(wts)) {
                stop(sprintf("Length mismatch at vertex %d: adj.list=%d, weight.list=%d.",
                             i, length(nbrs), length(wts)))
            }
            if (anyNA(nbrs)) {
                stop(sprintf("adj.list[[%d]] contains NA.", i))
            }
            if (any(nbrs < 1L | nbrs > n)) {
                stop(sprintf("adj.list[[%d]] has indices outside 1..n.", i))
            }
            if (any(nbrs == i)) {
                stop(sprintf("adj.list[[%d]] contains self-loops; self indices are not allowed.", i))
            }
            if (anyDuplicated(nbrs)) {
                stop(sprintf("adj.list[[%d]] contains duplicate neighbors.", i))
            }
            if (any(!is.finite(wts))) {
                stop(sprintf("weight.list[[%d]] contains non-finite values.", i))
            }
            if (any(wts <= 0)) {
                stop(sprintf("weight.list[[%d]] contains non-positive values.", i))
            }

            adj.list.norm[[i]] <- nbrs
            weight.list.norm[[i]] <- wts
        }

        ## Undirected/symmetric graph validation with reciprocal weight agreement
        for (i in seq_len(n)) {
            nbrs <- adj.list.norm[[i]]
            wts <- weight.list.norm[[i]]
            if (!length(nbrs)) {
                next
            }
            for (j.idx in seq_along(nbrs)) {
                j <- nbrs[[j.idx]]
                rev.idx <- match(i, adj.list.norm[[j]])
                if (is.na(rev.idx)) {
                    stop(sprintf("Graph must be undirected: edge %d -> %d has no reciprocal entry.", i, j))
                }
                w.ij <- wts[[j.idx]]
                w.ji <- weight.list.norm[[j]][[rev.idx]]
                if (abs(w.ij - w.ji) > tol * max(1, abs(w.ij), abs(w.ji))) {
                    stop(sprintf("Reciprocal edge weights mismatch for (%d, %d): %.12g vs %.12g.",
                                 i, j, w.ij, w.ji))
                }
            }
        }

        adj.list.0based <- lapply(adj.list.norm, function(v) as.integer(v - 1L))
        weight.list.cpp <- lapply(weight.list.norm, as.double)
    }

    ## ==================== Posterior Parameters ====================

    if (!is.logical(with.posterior) || length(with.posterior) != 1)
        stop("with.posterior must be TRUE or FALSE")

    if (!is.numeric(credible.level) || length(credible.level) != 1 ||
        credible.level <= 0 || credible.level >= 1)
        stop("credible.level must be in (0, 1)")

    if (!is.numeric(n.posterior.samples) || length(n.posterior.samples) != 1 ||
        n.posterior.samples < 100)
        stop("n.posterior.samples must be >= 100")

    if (is.null(posterior.seed)) {
        posterior.seed <- sample.int(.Machine$integer.max, 1)
    }

    ## ==================== Parameter density.alpha ====================

    if (!is.numeric(density.alpha) || length(density.alpha) != 1) {
        stop("density.alpha must be a single numeric value")
    }
    if (!is.finite(density.alpha) || density.alpha < 0.1 || density.alpha > 2.0) {
        stop("density.alpha must be finite and in [0.1, 2], got ", density.alpha)
    }

    ## ==================== Parameter density.epsilon ====================

    if (!is.numeric(density.epsilon) || length(density.epsilon) != 1) {
        stop("density.epsilon must be a single numeric value")
    }
    if (!is.finite(density.epsilon) || density.epsilon <= 0) {
        stop("density.epsilon must be finite and positive, got ", density.epsilon)
    }

    ## ---- reference measure: d_k clamping + weight compression controls ----

    if (!is.logical(clamp.dk) || length(clamp.dk) != 1L || is.na(clamp.dk)) {
        stop("clamp.dk must be a single, non-NA logical.")
    }

    if (!is.numeric(dk.clamp.median.factor) || length(dk.clamp.median.factor) != 1L ||
        is.na(dk.clamp.median.factor) || !is.finite(dk.clamp.median.factor) ||
        dk.clamp.median.factor <= 1.0) {
        stop(sprintf("dk.clamp.median.factor must be a single finite number > 1 (got %s).",
                     format(dk.clamp.median.factor)))
    }

    if (!is.numeric(target.weight.ratio) || length(target.weight.ratio) != 1L ||
        is.na(target.weight.ratio) || !is.finite(target.weight.ratio) ||
        target.weight.ratio <= 10.0) {
        stop(sprintf("target.weight.ratio must be a single finite number > 10 (got %s).",
                     format(target.weight.ratio)))
    }

    if (!is.numeric(pathological.ratio.threshold) || length(pathological.ratio.threshold) != 1L ||
        is.na(pathological.ratio.threshold) || !is.finite(pathological.ratio.threshold) ||
        pathological.ratio.threshold <= 1e6) {
        stop(sprintf("pathological.ratio.threshold must be a single finite number > 1e6 (got %s).",
                     format(pathological.ratio.threshold)))
    }

    ## Strongly recommended for consistency with intended behavior
    if (pathological.ratio.threshold <= target.weight.ratio) {
        stop(sprintf("pathological.ratio.threshold must be > target.weight.ratio (got %.3e <= %.3e).",
                     pathological.ratio.threshold, target.weight.ratio))
    }
    
    ## ==================== Parameter compute.extremality Validation ====================

    if (!is.logical(compute.extremality) || length(compute.extremality) != 1) {
        stop("compute.extremality must be TRUE or FALSE")
    }

    if (is.na(compute.extremality)) {
        stop("compute.extremality cannot be NA")
    }

    ## ==================== Parameter p.threshold Validation ====================

    if (!is.numeric(p.threshold) || length(p.threshold) != 1) {
        stop("p.threshold must be a single numeric value")
    }

    if (is.na(p.threshold) || !is.finite(p.threshold)) {
        stop("p.threshold must be finite")
    }

    ## p.threshold = 0 is allowed to skip hop radius computation
    ## Otherwise it must be in (0, 1]
    if (p.threshold < 0 || p.threshold > 1) {
        stop(sprintf("p.threshold must be in [0, 1], got %.3f", p.threshold))
    }
    ## ==================== Parameter max.hop Validation ====================

    if (!is.numeric(max.hop) || length(max.hop) != 1) {
        stop("max.hop must be a single numeric value")
    }

    ## Convert to integer
    max.hop <- as.integer(max.hop)

    if (is.na(max.hop)) {
        stop("max.hop cannot be NA")
    }

    if (max.hop < 1) {
        stop(sprintf("max.hop must be at least 1 (got %d)", max.hop))
    }

    ## Practical upper bound check
    if (max.hop > 1000) {
        warning(sprintf("max.hop = %d is unusually large; typical values are 10-50", max.hop))
    }

    ## ==================== Parameter use.counting.measure Validation ====================

    if (!is.logical(use.counting.measure) || length(use.counting.measure) != 1) {
        stop("use.counting.measure must be TRUE or FALSE")
    }

    if (is.na(use.counting.measure)) {
        stop("use.counting.measure cannot be NA")
    }

    ## ==================== Numeric Parameter Validation ====================

    ## density.normalization
    if (!is.numeric(density.normalization) || length(density.normalization) != 1) {
        stop("density.normalization must be a single numeric value")
    }

    if (is.na(density.normalization) || !is.finite(density.normalization)) {
        stop("density.normalization cannot be NA or infinite")
    }

    if (density.normalization < 0) {
        stop(sprintf("density.normalization must be non-negative (got %.3f)",
                     density.normalization))
    }

    ## t.diffusion
    if (!is.numeric(t.diffusion) || length(t.diffusion) != 1) {
        stop("t.diffusion must be a single numeric value")
    }

    if (is.na(t.diffusion) || !is.finite(t.diffusion)) {
        stop("t.diffusion cannot be NA or infinite")
    }

    if (t.diffusion < 0) {
        stop(sprintf("t.diffusion must be non-negative (got %.3f; use 0 for auto)",
                     t.diffusion))
    }

    ## density.uniform.pull
    if (!is.numeric(density.uniform.pull) || length(density.uniform.pull) != 1) {
        stop("density.uniform.pull must be a single numeric value")
    }

    if (is.na(density.uniform.pull) || !is.finite(density.uniform.pull)) {
        stop("density.uniform.pull cannot be NA or infinite")
    }

    if (density.uniform.pull < 0) {
        stop(sprintf("density.uniform.pull must be non-negative (got %.3f; use 0 for auto)",
                     density.uniform.pull))
    }

    ## response.penalty.exp
    if (!is.numeric(response.penalty.exp) || length(response.penalty.exp) != 1) {
        stop("response.penalty.exp must be a single numeric value")
    }

    if (is.na(response.penalty.exp) || !is.finite(response.penalty.exp)) {
        stop("response.penalty.exp cannot be NA or infinite")
    }

    if (response.penalty.exp > 2.0) {
        warning(sprintf("response.penalty.exp=%.3f must be less than 2.0.\n",
                        response.penalty.exp))
    }

    ## t.scale.factor
    if (!is.numeric(t.scale.factor) || length(t.scale.factor) != 1) {
        stop("t.scale.factor must be a single numeric value")
    }

    if (t.scale.factor <= 0) {
        stop("t.scale.factor must be positive (got ", t.scale.factor,
             "). Typical range: [0.1, 1.0]")
    }

    if (t.scale.factor > 2.0) {
        warning("Very large t.scale.factor (", t.scale.factor,
                "): may cause oscillations. Consider reducing to below 1.0.")
    }

    ## beta.coef.factor
    if (!is.numeric(beta.coef.factor) || length(beta.coef.factor) != 1) {
        stop("beta.coef.factor must be a single numeric value")
    }

    if (beta.coef.factor <= 0) {
        stop("beta.coef.factor must be positive (got ", beta.coef.factor,
             "). Typical range: [0.05, 0.3]")
    }

    if (beta.coef.factor < 0.01) {
        warning("Very small beta.coef.factor (", beta.coef.factor,
                "): density may collapse. Consider increasing to at least 0.05.")
    }

    if (beta.coef.factor > 0.5) {
        warning("Very large beta.coef.factor (", beta.coef.factor,
                "): may over-suppress structure. Consider reducing to below 0.3.")
    }

    ## n.eigenpairs
    if (!is.numeric(n.eigenpairs) || length(n.eigenpairs) != 1) {
        stop("n.eigenpairs must be a single numeric value")
    }

    n.eigenpairs <- as.integer(n.eigenpairs)

    if (n.eigenpairs < 10) {
        stop(sprintf("n.eigenpairs must be at least 10 (got %d)", n.eigenpairs))
    }

    if (n.eigenpairs > n) {
        stop(sprintf("n.eigenpairs cannot exceed n (got %d, n=%d)", n.eigenpairs, n))
    }

    ## Warn if n.eigenpairs seems too small
    if (n.eigenpairs < 50 && n > 100) {
        warning(sprintf(paste0("n.eigenpairs=%d may be too small for n=%d observations.\n",
                              "Consider n.eigenpairs >= 100 for better frequency resolution."),
                       n.eigenpairs, n))
    }

    ## epsilon.y
    if (!is.numeric(epsilon.y) || length(epsilon.y) != 1) {
        stop("epsilon.y must be a single numeric value")
    }

    if (is.na(epsilon.y) || !is.finite(epsilon.y)) {
        stop("epsilon.y cannot be NA or infinite")
    }

    if (epsilon.y <= 0) {
        stop(sprintf("epsilon.y must be positive (got %.3e)", epsilon.y))
    }

    if (epsilon.y >= 0.1) {
        warning(sprintf("epsilon.y=%.3f is very large; convergence may be premature", epsilon.y))
    }

    ## epsilon.rho
    if (!is.numeric(epsilon.rho) || length(epsilon.rho) != 1) {
        stop("epsilon.rho must be a single numeric value")
    }

    if (is.na(epsilon.rho) || !is.finite(epsilon.rho)) {
        stop("epsilon.rho cannot be NA or infinite")
    }

    if (epsilon.rho <= 0) {
        stop(sprintf("epsilon.rho must be positive (got %.3e)", epsilon.rho))
    }

    if (epsilon.rho >= 0.1) {
        warning(sprintf("epsilon.rho=%.3f is very large; convergence may be premature", epsilon.rho))
    }

    ## max.iterations
    if (!is.numeric(max.iterations) || length(max.iterations) != 1) {
        stop("max.iterations must be a single numeric value")
    }

    max.iterations <- as.integer(max.iterations)

    if (max.iterations < 1) {
        stop(sprintf("max.iterations must be at least 1 (got %d)", max.iterations))
    }

    if (max.iterations > 1000) {
        warning(sprintf(paste0("max.iterations=%d is very large.\n",
                              "Typical convergence occurs within 5-10 iterations."),
                       max.iterations))
    }

    ## max.ratio.threshold
    if (!is.numeric(max.ratio.threshold) || length(max.ratio.threshold) != 1) {
        stop("max.ratio.threshold must be numeric.")
    }
        
    if (max.ratio.threshold < 0 || max.ratio.threshold >= 0.2) {
        stop("max.ratio.threshold must be in [0, 0.2).")
    }

    if (!is.logical(apply.geometric.pruning) ||
        length(apply.geometric.pruning) != 1L ||
        is.na(apply.geometric.pruning)) {
        stop("apply.geometric.pruning must be TRUE or FALSE.")
    }

    if (!is.numeric(path.edge.ratio.percentile) || length(path.edge.ratio.percentile) != 1 ||
        path.edge.ratio.percentile < 0 || path.edge.ratio.percentile > 1)
        stop("path.edge.ratio.percentile must be in [0, 1].")

    ## threshold.percentile
    if (!is.numeric(threshold.percentile) || length(threshold.percentile) != 1 ||
        threshold.percentile < 0 || threshold.percentile > 0.5) {
        stop("threshold.percentile must be in [0, 0.5].")
    }

    knn.cache.mode <- match.arg(knn.cache.mode)
    if (has.precomputed.graph) {
        if (!identical(knn.cache.mode, "none") || !is.null(knn.cache.path)) {
            warning("adj.list/weight.list provided: knn.cache.path/mode are ignored.")
        }
        knn.cache.path <- NULL
        knn.cache.mode <- "none"
    }
    knn.cache.path <- .normalize.knn.cache.path(knn.cache.path, knn.cache.mode)
    knn.cache.mode.id <- switch(knn.cache.mode,
                                none = 0L,
                                read = 1L,
                                write = 2L,
                                readwrite = 3L)

    ## verbose.level
    if (length(verbose.level) != 1 || !(verbose.level %in% 0:3) )
        stop("verbose.level must be an integer in {0,1,2,3}.")

    ## pca.dim
    if (!is.null(pca.dim)) {
        if (!is.numeric(pca.dim) || length(pca.dim) != 1 || pca.dim < 1 || pca.dim != floor(pca.dim))
            stop("pca.dim must be a positive integer or NULL.")
    }

    ## variance.explained
    if (!is.null(variance.explained)) {
        if (!is.numeric(variance.explained) || length(variance.explained) != 1 ||
            variance.explained <= 0 || variance.explained > 1)
            stop("variance.explained must be in (0, 1], or NULL.")
    }

    ## store.projected.X
    if (!is.logical(store.projected.X) || length(store.projected.X) != 1 || is.na(store.projected.X)) {
        stop("store.projected.X must be TRUE or FALSE.")
    }

    if (has.precomputed.graph && !is.null(pca.dim) && ncol(X) > pca.dim) {
        pca.dim <- NULL
    }

    ## PCA (optional)
    pca_info <- NULL
    if (!is.null(pca.dim) && ncol(X) > pca.dim) {
        if (verbose.level >= 1) message("High-dimensional data detected. Performing PCA.")
        original_dim <- ncol(X)
        if (!is.null(variance.explained)) {
            pca_analysis <- pca.optimal.components(
                X, variance.threshold = variance.explained, max.components = pca.dim
            )
            n_components <- pca_analysis$n.components
            if (verbose.level == 1) {
                message(sprintf("Using %d PCs (explains %.2f%% variance)",
                                n_components, 100 * pca_analysis$variance.explained))
            }
            X <- pca.project(X, pca_analysis$pca.result, n_components)
            pca_info <- list(
                original_dim = original_dim,
                n_components = n_components,
                variance_explained = pca_analysis$variance.explained,
                cumulative_variance = pca_analysis$cumulative.variance
            )
        } else {
            if (verbose.level >= 1) message(sprintf("Projecting to first %d PCs", pca.dim))
            pca_result <- prcomp(X)
            X <- pca.project(X, pca_result, pca.dim)
            variance_explained <- sum(pca_result$sdev[1:pca.dim]^2) / sum(pca_result$sdev^2)
            pca_info <- list(
                original_dim = original_dim,
                n_components = pca.dim,
                variance_explained = variance_explained
            )
        }

        if (isTRUE(store.projected.X)) {
            pca_info$X_projected <- X
        }
    }

    ## gflowx owns only the archived estimator. Generic graph construction is
    ## delegated to dgraphs and the native backend always receives an explicit
    ## graph.
    if (!has.precomputed.graph) {
        if (verbose.level >= 1L) {
            message("Constructing the archived estimator graph with dgraphs.")
        }
        archive.graph <- dgraphs::create.single.iknn.graph(
            X = X,
            k = k,
            max.path.edge.ratio.deviation.thld =
                if (isTRUE(apply.geometric.pruning)) max.ratio.threshold else 0,
            prune.method =
                if (isTRUE(apply.geometric.pruning)) "global.geodesic" else "none",
            path.edge.ratio.percentile = path.edge.ratio.percentile,
            threshold.percentile = threshold.percentile,
            compute.full = FALSE,
            pca.dim = NULL,
            variance.explained = NULL,
            knn.cache.path = knn.cache.path,
            knn.cache.mode = knn.cache.mode,
            verbose = verbose.level >= 1L
        )
        adj.list.0based <- lapply(
            archive.graph$adj_list,
            function(v) as.integer(v - 1L)
        )
        weight.list.cpp <- lapply(archive.graph$weight_list, as.double)
    }

        
    ## ==================== String Parameter Validation ====================

    ## Match filter.type
    filter.type <- match.arg(filter.type)

    ## Match t.update
    t.update <- match.arg(t.update)

    ## Match dense.fallback
    dense.fallback <- match.arg(dense.fallback)
    dense.fallback.mode <- switch(
        dense.fallback,
        auto = 0L,
        never = 1L,
        always = 2L
    )

    ## Match triangle.policy
    triangle.policy <- match.arg(triangle.policy)
    triangle.policy.mode <- switch(
        triangle.policy,
        auto = 0L,
        never = 1L,
        always = 2L
    )

    if (!is.numeric(t.update.max.mult) || length(t.update.max.mult) != 1L ||
        !is.finite(t.update.max.mult) || t.update.max.mult < 1.0) {
        stop("t.update.max.mult must be a finite numeric scalar >= 1.0")
    }

    t.update.mode <- if (t.update == "fixed") 0L else 1L

    ## ==================== Call C++ Function ====================

    fit <- .Call(
        S_fit_rdgraph_regression,
        X,
        as.double(y),
        if (is.null(y.vertices)) NULL else as.integer(y.vertices),
        as.integer(k + 1L), # this is to account for the fact that ANN library is set up to return for query point as the first elements of the list of kNN's
        if (is.null(adj.list.0based)) NULL else adj.list.0based,
        if (is.null(weight.list.cpp)) NULL else weight.list.cpp,
        as.logical(with.posterior),
        as.logical(return.posterior.samples),
        as.double(credible.level),
        as.integer(n.posterior.samples),
        as.integer(posterior.seed),
        as.logical(use.counting.measure),
        as.double(density.normalization),
        as.double(t.diffusion),
        as.double(density.uniform.pull),
        as.double(response.penalty.exp),
        as.double(t.scale.factor),
        as.double(beta.coef.factor),
        as.integer(t.update.mode),
        as.double(t.update.max.mult),
        as.integer(n.eigenpairs),
        as.character(filter.type),
        as.double(epsilon.y),
        as.double(epsilon.rho),
        as.integer(max.iterations),
        as.double(if (isTRUE(apply.geometric.pruning)) max.ratio.threshold + 1.0 else 0.0),
        as.double(path.edge.ratio.percentile),
        as.double(threshold.percentile),
        as.double(density.alpha),
        as.double(density.epsilon),
        as.logical(clamp.dk),
        as.double(dk.clamp.median.factor),
        as.double(target.weight.ratio),
        as.double(pathological.ratio.threshold),
        as.logical(compute.extremality),
        as.double(p.threshold),
        as.integer(max.hop),
        if (is.null(knn.cache.path)) NULL else as.character(knn.cache.path),
        as.integer(knn.cache.mode.id),
        as.integer(dense.fallback.mode),
        as.integer(triangle.policy.mode),
        as.integer(verbose.level),
        PACKAGE = "gflowx"
    )

    ## ==================== Post-Processing ====================

    ## Add class
    class(fit) <- c("knn.riem.fit", "riem.dcx")

    ## Store call for reproducibility
    attr(fit, "call") <- match.call()

    ## Store key parameters for summary()
    attr(fit, "n") <- n
    attr(fit, "d") <- d
    attr(fit, "k") <- k
    if (is.list(fit$parameters)) {
        fit$parameters$graph.source <-
            if (has.precomputed.graph) "precomputed" else "dgraphs"
    }
    if (!is.null(pca_info)) attr(fit, "pca") <- pca_info
    if (isTRUE(store.projected.X)) attr(fit, "X.used") <- X

    if (!is.null(y.vertices)) {
        attr(fit, "y.vertices") <- as.integer(y.vertices)
    }

    return(fit)
}

#' Semi-supervised Riemannian Graph Regression (subset-labeled y)
#'
#' Convenience wrapper around \code{fit.rdgraph.regression()} for the case when
#' \code{y} is provided only on a subset of vertices specified by \code{y.vertices}.
#'
#' @param X Feature matrix (dense numeric matrix or dgCMatrix).
#' @param y Numeric vector of response values for labeled vertices only.
#' @param y.vertices Integer vector of labeled vertex indices (1-based).
#' @param ... Passed to \code{fit.rdgraph.regression()}; must include \code{k} and other tuning parameters.
#'
#' @return An object of class \code{c("knn.riem.fit", "riem.dcx")}.
fit.rdgraph.regression.semiy <- function(X, y, y.vertices, ...) {
    if (missing(y.vertices) || is.null(y.vertices)) {
        stop("y.vertices must be provided (integer 1-based indices of labeled vertices)")
    }
    fit.rdgraph.regression(X = X, y = y, y.vertices = y.vertices, ...)
}

#' Extract Fitted Values from kNN Riemannian Regression Fit
#'
#' @param object Object of class \code{knn.riem.fit}
#' @param ... Additional arguments (currently unused)
#' @return Numeric vector of fitted values (length n)
#' @export
fitted.knn.riem.fit <- function(object, ...) {
    object$fitted.values
}

#' @export
residuals.knn.riem.fit <- function(object, ...) {
    object$residuals
}

#' @export
summary.knn.riem.fit <- function(object, ...) {
    # Extract key information
    resid <- object$residuals
    n <- length(object$y)

    # Compute fit quality metrics
    rss <- sum(resid^2)
    tss <- sum((object$y - mean(object$y))^2)
    r_squared <- 1 - rss / tss
    adj_r_squared <- 1 - (1 - r_squared) * (n - 1) / (n - object$parameters$k - 1)
    rmse <- sqrt(mean(resid^2))
    mae <- mean(abs(resid))

    # Get optimal iteration info
    optimal_iter <- object$optimal.iteration
    total_iters <- object$iteration$n.iterations

    # Get GCV information
    gcv_optimal <- object$gcv$gcv.optimal[optimal_iter]
    gcv_initial <- object$gcv$gcv.optimal[1]
    gcv_improvement <- (gcv_initial - gcv_optimal) / gcv_initial * 100

    # Check if parameters were auto-selected
    params <- object$parameters
    auto_selected <- list(
        t.diffusion = (params$t.diffusion > 0),  # Was auto-selected if positive
        density.uniform.pull = (params$density.uniform.pull > 0),
        response.penalty.exp = attr(object, "gamma.auto.selected")  # Need to store this
    )

    # Graph structure info
    n_vertices <- object$graph$n.vertices
    n_edges <- object$graph$n.edges
    avg_degree <- 2 * n_edges / n_vertices

    # Convergence diagnostics
    final_response_change <- if (length(object$iteration$response.changes) > 0) {
        utils::tail(object$iteration$response.changes, 1)
    } else {
        NA
    }

    structure(
        list(
            call = attr(object, "call"),

            # Data dimensions
            n = n,
            d = attr(object, "d"),

            # Model parameters
            k = params$k,
            response.penalty.exp = params$response.penalty.exp,
            filter.type = params$filter.type,
            n.eigenpairs = params$n.eigenpairs,

            # Auto-selection flags
            auto.selected = auto_selected,

            # Graph structure
            graph = list(
                n.vertices = n_vertices,
                n.edges = n_edges,
                avg.degree = avg_degree,
                density = n_edges / choose(n_vertices, 2)
            ),

            # Convergence
            converged = object$iteration$converged,
            n.iterations = total_iters,
            optimal.iteration = optimal_iter,
            final.response.change = final_response_change,

            # Fit quality
            residuals = resid,
            sigma = sd(resid),
            rmse = rmse,
            mae = mae,
            r.squared = r_squared,
            adj.r.squared = adj_r_squared,

            # GCV diagnostics
            gcv.optimal = gcv_optimal,
            gcv.initial = gcv_initial,
            gcv.improvement.pct = gcv_improvement,

            # Response characteristics
            y.range = range(object$y),
            y.mean = mean(object$y),
            y.sd = sd(object$y)
        ),
        class = "summary.knn.riem.fit"
    )
}

#' @export
print.summary.knn.riem.fit <- function(x, digits = 4, ...) {
    cat("\n")
    cat("========================================================\n")
    cat("k-NN Riemannian Graph Regression Summary\n")
    cat("========================================================\n\n")

    # Call
    cat("Call:\n")
    print(x$call)
    cat("\n")

    # Data and Model
    cat("Data and Model:\n")
    cat(sprintf("  Observations:        n = %d\n", x$n))
    cat(sprintf("  Features:            d = %d\n", x$d))
    cat(sprintf("  Neighbors:           k = %d\n", x$k))
    cat(sprintf("  Filter type:         %s\n", x$filter.type))
    cat(sprintf("  Eigenpairs:          %d\n", x$n.eigenpairs))
    cat("\n")

    # Parameter selection
    cat("Parameter Selection:\n")
    cat(sprintf("  Response penalty:    gamma = %.3f", x$response.penalty.exp))
    if (!is.null(x$auto.selected$response.penalty.exp) && x$auto.selected$response.penalty.exp) {
        cat(" [auto-selected]")
    }
    cat("\n")
    cat(sprintf("  Diffusion time:      t = %.4f", x$k))  # Need actual value
    if (!is.null(x$auto.selected$t.diffusion) && x$auto.selected$t.diffusion) {
        cat(" [auto-selected]")
    }
    cat("\n\n")

    # Graph structure
    cat("Graph Structure:\n")
    cat(sprintf("  Vertices:            %d\n", x$graph$n.vertices))
    cat(sprintf("  Edges:               %d\n", x$graph$n.edges))
    cat(sprintf("  Avg degree:          %.1f\n", x$graph$avg.degree))
    cat(sprintf("  Graph density:       %.4f\n", x$graph$density))
    cat("\n")

    # Convergence
    cat("Convergence:\n")
    cat(sprintf("  Status:              %s\n",
                if(x$converged) "CONVERGED" else "Max iterations reached"))
    cat(sprintf("  Total iterations:    %d\n", x$n.iterations))
    cat(sprintf("  Optimal iteration:   %d (min GCV)\n", x$optimal.iteration))
    if (!is.na(x$final.response.change)) {
        cat(sprintf("  Final change:        %.2e\n", x$final.response.change))
    }
    cat("\n")

    # Fit quality
    cat("Fit Quality:\n")
    cat(sprintf("  R-squared:           %.4f\n", x$r.squared))
    cat(sprintf("  Adj R-squared:       %.4f\n", x$adj.r.squared))
    cat(sprintf("  RMSE:                %.4f\n", x$rmse))
    cat(sprintf("  MAE:                 %.4f\n", x$mae))
    cat(sprintf("  Residual SD:         %.4f\n", x$sigma))
    cat("\n")

    # GCV diagnostics
    cat("GCV Diagnostics:\n")
    cat(sprintf("  Optimal GCV:         %.4e\n", x$gcv.optimal))
    cat(sprintf("  Initial GCV:         %.4e\n", x$gcv.initial))
    cat(sprintf("  Improvement:         %.1f%%\n", x$gcv.improvement.pct))
    cat("\n")

    # Residual summary
    cat("Residuals:\n")
    print(summary(x$residuals), digits = digits)
    cat("\n")

    # Response summary
    cat("Response:\n")
    cat(sprintf("  Range:               [%.3f, %.3f]\n",
                x$y.range[1], x$y.range[2]))
    cat(sprintf("  Mean:                %.3f\n", x$y.mean))
    cat(sprintf("  SD:                  %.3f\n", x$y.sd))

    cat("\n")
    cat("========================================================\n")

    invisible(x)
}

#' Diagnostic Plots for k-NN Riemannian Graph Regression
#'
#' @param results.list List of fitted model objects from fit.rdgraph.regression
#' @param k.values Vector of k values corresponding to results.list
#' @param response Vector of observed response values (for computing relative ranges)
#' @param main.title Optional main title for the entire plot grid
#' @param mark.optimal Logical; if TRUE, mark the k selected by GCV in red
#'
#' @return Invisibly returns a data frame with all diagnostic metrics
k.diagnostics.plots <- function(results.list,
                                k.values,
                                response,
                                main.title = "k-NN Bandwidth Diagnostics",
                                mark.optimal = TRUE) {

  n.models <- length(results.list)
  if (length(k.values) != n.models) {
    stop("Length of k.values must match length of results.list")
  }

    ## Compute all diagnostics
    gcv.scores <- numeric(n.models)
    range.lower <- numeric(n.models)
    range.upper <- numeric(n.models)
    stability.corr <- rep(NA_real_, n.models)
    n.extrema <- numeric(n.models)

    y.mean <- mean(response)

    for (j in seq_len(n.models)) {
        fit <- results.list[[j]]

        ## GCV score
        gcv.scores[j] <- fit$gcv$gcv.optimal[fit$optimal.iteration]

        ## Fitted values range (winsorized)
        fitted.wins <- winsorize(fit$fitted.values, p = 0.02)
        range.lower[j] <- min(fitted.wins) / y.mean
        range.upper[j] <- max(fitted.wins) / y.mean

        ## Stability correlation
        if (j > 1) {
            fit.prev <- results.list[[j - 1]]
            stability.corr[j] <- cor(fit.prev$fitted.values, fit$fitted.values)
        }

        ## Number of local extrema
        extrema.res <- compute.extrema.hop.nbhds(
            fit$graph$adj.list,
            fit$graph$edge.length.list,
            fit$fitted.values
        )
        n.extrema[j] <- nrow(extrema.res$extrema_df)
    }

    ## Find optimal k by GCV
    optimal.idx <- which.min(gcv.scores)
    optimal.k <- k.values[optimal.idx]

    ## Set up 2x2 plotting grid
    old.par <- par(mfrow = c(2, 2),
                   mar = c(4, 4, 2.5, 1),
                   oma = c(0, 0, 2, 0))
    on.exit(par(old.par))

    ## Plot 1: GCV scores (expressed as deviation from minimum)
    gcv.min <- min(gcv.scores)
    gcv.relative <- (gcv.scores - gcv.min) / gcv.min * 100  # percent increase from minimum

    plot(k.values, gcv.relative,
         type = "b",
         pch = 19,
         col = "steelblue",
         xlab = "Neighborhood size k",
         ylab = "GCV increase from minimum (%)",
         main = "Cross-Validation Criterion",
         las = 1,
         panel.first = grid(col = "gray90", lty = 1))

    if (mark.optimal) {
        abline(v = optimal.k, col = "firebrick", lwd = 2, lty = 2)
        abline(h = 0, col = "gray40", lty = 1)
        text(optimal.k, max(gcv.relative),
             labels = paste0("k = ", optimal.k),
             pos = 4, col = "firebrick", cex = 0.9)
    }

    ## Plot 2: Fitted values range
    y.range <- range(c(range.lower, range.upper))
    plot(k.values, range.lower,
         type = "n",
         ylim = y.range,
         xlab = "Neighborhood size k",
         ylab = "Relative range (normalized by mean response)",
         main = "Prediction Range Stability",
         las = 1,
         panel.first = grid(col = "gray90", lty = 1))

    ## Add shaded region
    polygon(c(k.values, rev(k.values)),
            c(range.lower, rev(range.upper)),
            col = adjustcolor("steelblue", alpha.f = 0.3),
            border = NA)

    lines(k.values, range.lower, type = "b", pch = 19, col = "steelblue")
    lines(k.values, range.upper, type = "b", pch = 19, col = "steelblue")
    abline(h = 1, col = "gray40", lwd = 1.5, lty = 1)

    if (mark.optimal) {
        abline(v = optimal.k, col = "firebrick", lwd = 2, lty = 2)
    }

    ## Add legend
    legend("topright",
           legend = c("Lower bound", "Upper bound", "Mean response"),
           col = c("steelblue", "steelblue", "gray40"),
           lty = c(1, 1, 1),
           lwd = c(1, 1, 1.5),
           pch = c(19, 19, NA),
           bty = "n",
           cex = 0.8)

    ## Plot 3: Stability correlation
    plot(k.values, stability.corr,
         type = "b",
         pch = 19,
         col = "steelblue",
         xlab = "Neighborhood size k",
         ylab = "Correlation with k-1 fit",
         main = "Fit Stability Across k",
         las = 1,
         ylim = c(0, 1),
         panel.first = grid(col = "gray90", lty = 1))

    if (mark.optimal) {
        abline(v = optimal.k, col = "firebrick", lwd = 2, lty = 2)
    }

    ## Plot 4: Number of extrema
    plot(k.values, n.extrema,
         type = "b",
         pch = 19,
         col = "steelblue",
         xlab = "Neighborhood size k",
         ylab = "Count",
         main = "Geometric Complexity (Local Extrema)",
         las = 1,
         panel.first = grid(col = "gray90", lty = 1))

    if (mark.optimal) {
        abline(v = optimal.k, col = "firebrick", lwd = 2, lty = 2)
    }

    ## Add overall title
    mtext(main.title, outer = TRUE, cex = 1.2, font = 2)

    ## Return diagnostic data frame
    diagnostics.df <- data.frame(
        k = k.values,
        gcv = gcv.scores,
        range.lower = range.lower,
        range.upper = range.upper,
        stability.corr = stability.corr,
        n.extrema = n.extrema,
        is.optimal = k.values == optimal.k
    )

    invisible(diagnostics.df)
}
#' Compute Effective Degrees for All Vertices
#'
#' @description
#' Helper function that computes effective degrees by summing edge densities
#' incident to each vertex. This measures connectivity strength accounting for
#' edge weights.
#'
#' @param adj.list List of integer vectors giving neighbor indices for each vertex
#' @param edge_densities List of numeric vectors giving edge densities parallel
#'   to adj.list structure
#'
#' @return Numeric vector of effective degrees, one per vertex
#'
#' @keywords internal
compute.effective.degrees <- function(adj.list, edge_densities) {
    n <- length(adj.list)
    eff.deg <- numeric(n)

    for (i in seq_len(n)) {
        if (length(adj.list[[i]]) > 0 && length(edge_densities[[i]]) > 0) {
            eff.deg[i] <- sum(edge_densities[[i]], na.rm = TRUE)
        }
    }

    return(eff.deg)
}

#' Extract Optimal k from Diagnostics Data Frame
#'
#' @param diag.df Data frame returned by k.diagnostics.plots
#'
#' @return A list with optimal k value and its index
get.rcx.optimal.k <- function(diag.df) {
  if (!any(diag.df$is.optimal)) {
    stop("No optimal k found in diagnostics data frame")
  }

  optimal.idx <- which(diag.df$is.optimal)
  optimal.k <- diag.df$k[optimal.idx]

  list(
    k = optimal.k,
    idx = optimal.idx,
    gcv = diag.df$gcv[optimal.idx]
  )
}

#' Compute Probabilistic Extrema Neighborhoods
#'
#' @description
#' Identifies local extrema using probabilistic scoring (maxp/minp) and computes
#' their spatial persistence via hop-extremp radius. This statistical approach
#' distinguishes genuine extrema from spurious peaks caused by discretization
#' and estimation noise.
#'
#' @details
#' Classical extrema detection on discrete graphs uses a binary condition: a
#' vertex either is or is not a local extremum. This rigidity creates problems
#' when the graph is derived from continuous data, where sampling artifacts and
#' estimation errors introduce spurious extrema. The probabilistic framework
#' provides quantitative measures that enable flexible, data-driven filtering.
#'
#' \strong{The Three-Stage Detection Pipeline}
#'
#' The function computes three complementary measures for each extremum candidate,
#' allowing downstream analysis to apply domain-specific filtering criteria.
#'
#' Stage 1 (Probabilistic Scoring): For each vertex, compute the extremum
#' probability score, which measures what fraction of the neighborhood density
#' mass satisfies the extremum condition. For maxima, this is the maxp score;
#' for minima, the minp score. Vertices are selected if their score exceeds
#' a high quantile threshold (default 95th percentile), ensuring we capture
#' vertices where most neighbors satisfy the condition.
#'
#' Stage 2 (Connectivity Filtering): Compute the effective degree for each
#' vertex, which sums the density weights of incident edges. This measures
#' how well-connected the vertex is to the graph structure. Vertices below
#' a low threshold (default 10th percentile) are discarded as poorly connected
#' outliers. This removes obvious outliers while retaining genuine extrema in
#' well-sampled regions.
#'
#' Stage 3 (Spatial Persistence): For vertices passing the first two filters,
#' compute the hop-extremp radius, which measures how far the extremum property
#' extends into the neighborhood. A vertex with large hop-extremp radius maintains
#' its extremum character over many hops, indicating a stable, genuine peak or
#' valley. Vertices with small radius may be local noise fluctuations.
#'
#' The function returns all three measures, allowing users to apply custom
#' filtering logic based on their domain knowledge and analytical needs.
#'
#' \strong{Relationship to Classical Detection}
#'
#' The function can optionally compute classical hop_idx for comparison. The
#' classical measure requires all neighbors to strictly satisfy the inequality,
#' corresponding to extremp = 1.0. The hop-extremp radius generalizes this by
#' allowing a threshold p < 1.0, making it more robust to outliers. For genuine
#' extrema, hop-extremp radius should be similar to or larger than hop_idx.
#' Large discrepancies indicate vertices where density weighting significantly
#' alters the assessment.
#'
#' \strong{Parameter Selection Guidelines}
#'
#' The extremp_quantile parameter controls initial sensitivity. Higher values
#' (0.99) produce fewer candidates but with higher confidence. Lower values
#' (0.90) produce more candidates, useful for exploratory analysis. The default
#' 0.95 balances sensitivity and specificity for most applications.
#'
#' The hop_extremp_threshold parameter determines the strictness of the spatial
#' persistence test. Higher values (0.95) require very strong evidence across
#' the neighborhood. Lower values (0.85) are more permissive. The default 0.90
#' works well when combined with the 95th percentile initial filter.
#'
#' The eff_degree_quantile parameter filters poorly connected vertices. The
#' default 0.10 removes the bottom 10% of vertices by connectivity, targeting
#' clear outliers without being overly aggressive. Increasing this value (0.20)
#' produces more conservative filtering.
#'
#' \strong{Computational Considerations}
#'
#' The function is optimized for efficiency through staged filtering. By
#' computing expensive hop-extremp radii only for vertices passing initial
#' filters, computational cost scales with the number of genuine extrema rather
#' than total vertices. For a graph with n vertices and average degree k:
#' \itemize{
#'   \item Stage 1 (extremp scores): `O(n*k)` - single pass over all vertices
#'   \item Stage 2 (effective degrees): `O(n*k)` - already computed, just filter
#'   \item Stage 3 (hop-extremp radii): `O(m*(n+E))` where m is number of candidates
#' }
#'
#' Typically m << n (perhaps 20-50 extrema in a graph of 2000 vertices), making
#' Stage 3 tractable even though each hop-extremp computation requires BFS.
#'
#' @param dcx Output from \code{fit.rdgraph.regression()}, a list containing
#'   the fitted model components including \code{vertex_cofaces}, \code{y_hat},
#'   \code{vertex_density}, and \code{graph$adj_list}.
#' @param extremp_quantile Numeric scalar in \eqn{[0,1]} giving the quantile threshold
#'   for initial extremp filtering. Vertices with extremp score above this
#'   quantile are considered candidates. Default: 0.95 (top 5%).
#' @param hop_extremp_threshold Numeric scalar in (0,1] giving the extremp
#'   threshold for hop-extremp radius computation. Default: 0.90.
#' @param eff_degree_quantile Numeric scalar in \eqn{[0,1]} giving the quantile
#'   threshold for effective degree filtering. Vertices below this quantile
#'   are discarded as poorly connected. Default: 0.10 (bottom 10%).
#' @param max_hop Integer giving maximum hop distance to explore. Limits
#'   computation time on large graphs. Default: 20.
#' @param compute_hop_idx Logical indicating whether to compute classical hop_idx
#'   for comparison. Default: FALSE (skip for efficiency).
#'
#' @return Data frame with one row per detected extremum (both maxima and minima),
#'   containing the following columns:
#' \describe{
#'   \item{vertex}{Integer vertex index (1-based)}
#'   \item{value}{Numeric function value at the vertex}
#'   \item{type}{Character: "max" for local maximum, "min" for local minimum}
#'   \item{extremp}{Numeric extremum probability score (maxp or minp)}
#'   \item{eff_degree}{Numeric effective degree of the vertex}
#'   \item{hop_extremp_radius}{Numeric or Inf: hop-extremp radius. Inf indicates
#'         global extremum. Values of 0 indicate vertex failed 1-hop condition.}
#'   \item{hop_idx}{Numeric or Inf or NA: classical hop index if
#'         \code{compute_hop_idx = TRUE}, otherwise NA. Inf indicates global
#'         extremum. NA indicates computation was skipped or failed.}
#'   \item{label}{Character: extremum label. Maxima labeled M1, M2, ... in
#'         descending order by value. Minima labeled m1, m2, ... in ascending
#'         order by value.}
#' }
#'
#' The returned data frame is sorted with all maxima first (descending by value),
#' then all minima (ascending by value). Users can apply additional filtering
#' based on the quantitative measures to identify high-confidence extrema.
#'
#' @examples
#' \dontrun{
#' # Fit regression model
#' dcx <- fit.rdgraph.regression(X, y, k = 50)
#'
#' # Detect probabilistic extrema with default parameters
#' pextrema <- compute.pextrema.nbhds(dcx)
#'
#' # Use more stringent filtering
#' pextrema_strict <- compute.pextrema.nbhds(
#'   dcx,
#'   extremp_quantile = 0.99,
#'   hop_extremp_threshold = 0.95
#' )
#'
#' # Apply custom filtering based on returned measures
#' high_confidence <- pextrema[
#'   pextrema$extremp > 0.90 &
#'   pextrema$eff_degree > quantile(pextrema$eff_degree, 0.25) &
#'   pextrema$hop_extremp_radius >= 3,
#' ]
#'
#' # Compare probabilistic vs classical measures
#' pextrema_with_classical <- compute.pextrema.nbhds(
#'   dcx,
#'   compute_hop_idx = TRUE
#' )
#'
#' # Examine discrepancies
#' with(pextrema_with_classical, {
#'   plot(hop_idx, hop_extremp_radius,
#'        main = "Classical vs Probabilistic Persistence")
#'   abline(0, 1, col = "red")
#' })
#' }
#'
#' @seealso [fit.rdgraph.regression()]
#'
#' @keywords internal
compute.pextrema.nbhds <- function(dcx,
                                   extremp_quantile = 0.95,
                                   hop_extremp_threshold = 0.90,
                                   eff_degree_quantile = 0.10,
                                   max_hop = 20L,
                                   compute_hop_idx = FALSE) {

    # ================================================================
    # ARGUMENT VALIDATION
    # ================================================================

    # Validate dcx is a riem.dcx object
    if (!inherits(dcx, "riem.dcx")) {
        stop("dcx must be an object of class 'riem.dcx' (output from fit.rdgraph.regression)")
    }

    required_fields <- c("fitted.values", "graph")
    missing_fields <- setdiff(required_fields, names(dcx$optimal.fit))
    if (length(missing_fields) > 0) {
        stop(sprintf("dcx$optimal.fit is missing required fields: %s",
                     paste(missing_fields, collapse = ", ")))
    }

    required_graph_fields <- c("adj.list", "vertex.densities", "edge.densities")
    missing_graph_fields <- setdiff(required_graph_fields, names(dcx$optimal.fit$graph))
    if (length(missing_graph_fields) > 0) {
        stop(sprintf("dcx$optimal.fit$graph is missing required fields: %s",
                     paste(missing_graph_fields, collapse = ", ")))
    }

    # Extract components from the actual structure
    yhat <- dcx$optimal.fit$fitted.values
    Ey <- mean(dcx$optimal.fit$y)
    rho <- dcx$optimal.fit$graph$vertex.densities
    adj.list <- dcx$optimal.fit$graph$adj.list

    n <- length(adj.list)

    if (length(yhat) != n) {
        stop(sprintf("Length of yhat (%d) must match length of adj.list (%d)",
                     length(yhat), n))
    }

    if (length(rho) != n) {
        stop(sprintf("Length of rho (%d) must match length of adj.list (%d)",
                     length(rho), n))
    }

    if (extremp_quantile < 0 || extremp_quantile > 1) {
        stop(sprintf("extremp_quantile must be in [0,1], got %.3f", extremp_quantile))
    }

    if (hop_extremp_threshold <= 0 || hop_extremp_threshold > 1) {
        stop(sprintf("hop_extremp_threshold must be in (0,1], got %.3f",
                     hop_extremp_threshold))
    }

    if (eff_degree_quantile < 0 || eff_degree_quantile > 1) {
        stop(sprintf("eff_degree_quantile must be in [0,1], got %.3f",
                     eff_degree_quantile))
    }

    if (max_hop <= 0) {
        stop(sprintf("max_hop must be positive, got %d", max_hop))
    }

    # ================================================================
    # STAGE 1: COMPUTE EXTREMP SCORES FOR ALL VERTICES
    # ================================================================

    vertex_maxp <- numeric(n)
    vertex_minp <- numeric(n)

    for (i in 1:n) {
        vertex_maxp[i] <- .maxp_no_check(i, adj.list[[i]], yhat, rho)
        vertex_minp[i] <- .minp_no_check(i, adj.list[[i]], yhat, rho)
    }

    # ================================================================
    # STAGE 2: COMPUTE EFFECTIVE DEGREES
    # ================================================================

    # Compute effective degrees directly from adjacency list and edge densities
    # For each vertex, sum the densities of all incident edges
    eff_degree <- numeric(n)
    edge_idx <- 1  # Track position in edge.densities vector

    for (i in 1:n) {
        total_edge_mass <- 0.0
        n_neighbors <- length(adj.list[[i]])

        if (n_neighbors > 0) {
            # Sum densities of edges incident to vertex i
            # Edge densities are stored in order of traversing adj.list
            edge_densities_i <- dcx$optimal.fit$graph$edge.densities[edge_idx:(edge_idx + n_neighbors - 1)]
            total_edge_mass <- sum(edge_densities_i)
            edge_idx <- edge_idx + n_neighbors
        }

        eff_degree[i] <- total_edge_mass
    }

    # ================================================================
    # SELECT CANDIDATES
    # ================================================================

    # Threshold for extremp scores
    maxp_thresh <- quantile(vertex_maxp, probs = extremp_quantile, na.rm = TRUE)
    minp_thresh <- quantile(vertex_minp, probs = extremp_quantile, na.rm = TRUE)

    # Threshold for effective degree (remove bottom percentile)
    eff_deg_thresh <- quantile(eff_degree, probs = eff_degree_quantile, na.rm = TRUE)

    # Select maxima candidates
    max_candidates <- which(vertex_maxp >= maxp_thresh & eff_degree >= eff_deg_thresh)

    # Select minima candidates
    min_candidates <- which(vertex_minp >= minp_thresh & eff_degree >= eff_deg_thresh)

    # ================================================================
    # STAGE 3: COMPUTE HOP-EXTREMP RADII
    # ================================================================

    # Initialize results data frame
    results <- data.frame(
        vertex = integer(),
        value = numeric(),
        rel_value = numeric(),
        type = character(),
        extremp = numeric(),
        eff_degree = numeric(),
        hop_extremp_radius = numeric(),
        hop_idx = numeric(),
        stringsAsFactors = FALSE
    )

    # Process maxima candidates
    if (length(max_candidates) > 0) {
        # Call C++ function for batch hop-extremp radius computation
        hop_maxp_radii <- .Call("S_compute_hop_extremp_radii_batch",
                                adj.list,
                                dcx$optimal.fit$graph$edge.densities,
                                rho,
                                as.integer(max_candidates),
                                yhat,
                                as.numeric(hop_extremp_threshold),
                                TRUE,  # detect_maxima
                                as.integer(max_hop),
                                PACKAGE = "gflowx")

        # Convert -1 (global extrema) to Inf
        hop_maxp_radii[hop_maxp_radii == -1] <- Inf

        hop_idx_max <- rep(NA_real_, length(max_candidates))

        # Build maxima data frame
        max_df <- data.frame(
            vertex = max_candidates,
            value = yhat[max_candidates],
            rel_value = yhat[max_candidates] / Ey,
            type = "max",
            extremp = vertex_maxp[max_candidates],
            eff_degree = eff_degree[max_candidates],
            hop_extremp_radius = hop_maxp_radii,
            hop_idx = hop_idx_max,
            stringsAsFactors = FALSE
        )

        results <- rbind(results, max_df)
    }

    # Process minima candidates
    if (length(min_candidates) > 0) {
        # Call C++ function for batch hop-extremp radius computation
        hop_minp_radii <- .Call("S_compute_hop_extremp_radii_batch",
                                adj.list,
                                dcx$optimal.fit$graph$edge.densities,
                                rho,
                                as.integer(min_candidates),
                                yhat,
                                as.numeric(hop_extremp_threshold),
                                FALSE,  # detect_maxima
                                as.integer(max_hop),
                                PACKAGE = "gflowx")

        # Convert -1 (global extrema) to Inf
        hop_minp_radii[hop_minp_radii == -1] <- Inf

        if (compute_hop_idx) {
            hop_idx_min <- rep(NA_real_, length(min_candidates))
        } else {
            hop_idx_min <- rep(NA_real_, length(min_candidates))
        }

        # Build minima data frame
        min_df <- data.frame(
            vertex = min_candidates,
            value = yhat[min_candidates],
            rel_value = yhat[min_candidates] / Ey,
            type = "min",
            extremp = vertex_minp[min_candidates],
            eff_degree = eff_degree[min_candidates],
            hop_extremp_radius = hop_minp_radii,
            hop_idx = hop_idx_min,
            stringsAsFactors = FALSE
        )

        results <- rbind(results, min_df)
    }

    # ================================================================
    # SORT AND LABEL
    # ================================================================

    if (nrow(results) == 0) {
        # No extrema found, return empty data frame with correct structure
        results$label <- character()
        return(results)
    }

    # Separate maxima and minima for sorting and labeling
    max_subset <- results[results$type == "max", ]
    min_subset <- results[results$type == "min", ]

    # Sort maxima by value (descending) and label M1, M2, ...
    if (nrow(max_subset) > 0) {
        max_subset <- max_subset[order(-max_subset$value), ]
        max_subset$label <- paste0("M", seq_len(nrow(max_subset)))
    }

    # Sort minima by value (ascending) and label m1, m2, ...
    if (nrow(min_subset) > 0) {
        min_subset <- min_subset[order(min_subset$value), ]
        min_subset$label <- paste0("m", seq_len(nrow(min_subset)))
    }

    # Combine: maxima first, then minima
    results <- rbind(max_subset, min_subset)
    rownames(results) <- NULL

    ## call to classical compute.extrema.hop.nbhds to compute classical hop_idx
    if (compute_hop_idx) {
        extr.res <- compute.extrema.hop.nbhds(adj.list,
                                              dcx$optimal.fit$graph$edge.length.list,
                                              yhat)

        ## Select vertices of results that are in extr.res and populate
        ## results$hop_idx with the corresponding value in extr.res$extrema_df
        classical.extr <- results$vertex[results$vertex %in% extr.res$extrema_df$vertex]
        for (v in classical.extr) {
            results$hop_idx[results$vertex == v] <- extr.res$extrema_df$hop_idx[extr.res$extrema_df$vertex == v]
        }
    }

    return(results)
}
