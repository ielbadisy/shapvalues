# shapvalues 0.1.2

* `sv$phi` is always a matrix (`n_obs x n_feature`) for single-output models,
  also when one row is explained (it used to collapse to a vector).
* Internal table operations use 'basetable' (row binding, the
  observation-feature merge in `sv_table()`, and the pathway counts in
  `sv_pathway()`). Values are unchanged.
* `sv_table()` rows for multi-output models (classification, survival) now
  follow a fixed order: observation, then feature, then output. Before, the
  output order within an observation-feature pair was arbitrary.

# shapvalues 0.1.1

* `theme_shapvalues()` is exported. Every `sv_*` plot returns a plain `ggplot`
  object, so you can re-theme it with any `ggplot2::theme_*()` via `+`.
* New multivariate attribution views:
  - `sv_compass()` places each observation where its attribution mass points;
    distance from the centre is the explanation concentration index `EC`.
  - `sv_concentration()` returns the compass numbers (`cx`, `cy`, `ec`,
    `angle`, `dominant`, `total_abs`), optionally over feature groups.
  - `sv_fingerprint()` draws a radial `|SHAP|` signature per observation or per
    k-means cluster of the SHAP profiles.
  - `sv_pathway()` is an alluvial of the top-`rank` contributors (node =
    feature, flow width = observations, flow colour = dominant SHAP sign);
    needs the `ggalluvial` package.
* One colour system for every plot, with selectable styles. `sv_style()` sets
  the active style; `sv_palette()` shows the colours. The default `"shapvalues"`
  style is the package's own diverging pair; the others are the `viridisLite`
  maps `"viridis"`, `"magma"`, `"inferno"`, `"plasma"`, `"cividis"`. Each plot
  also takes a `style =` argument.
* New: pairwise SHAP interaction values. `shapvalues(..., interactions = TRUE)`
  returns an `n x n_feature x n_feature` array in `$interactions` (the halved,
  symmetrised Grabisch-Roubens index; each row sums to the SHAP value).
  Single-output models only. Plot with `sv_interaction()`.
* Removed `sv_force()`. Use `sv_waterfall()` for a single-instance breakdown.
* `sv_waterfall()` restyled to match the Python `shap` waterfall: `value =
  feature` y-axis labels, value labels inside wide arrows / outside narrow
  ones, `+0` / `-0` for near-zero contributions, value axis at the bottom.

* Reworked the plot layer to match the Python `shap` library: white ground, no
  gridlines or panel border, light-grey axis lines only, the shap blue-to-red
  colour map (`#008bfb` to `#ff0d57`), and a slim unframed colour bar.
* `sv_beeswarm()` uses a deterministic centre-out swarm offset (no random
  jitter) with dashed per-feature row separators and the "SHAP value (impact on
  model output)" axis label.
* `sv_waterfall()` is now the shap waterfall: fat red/blue arrows, `E[f(x)]`
  and `f(x)` annotations, top axis, dashed connectors.
* `sv_force()` is the shap additive force layout (red features push above the
  baseline, blue pull back, blocks meet at the prediction).
* Every `sv_*` plot gains a `title` argument (default `NULL`);
  `sv_dependence()` gains `smooth` for a loess trend line and reports its y
  axis as "SHAP value for <feature>".
* Fixed `sv_dependence(colour_by = ...)`: the `colour` aesthetic was being
  silently overridden, so points rendered monochrome with no colour bar. The
  interaction colouring and colour bar now display, with smaller points to
  reduce overplotting.
* `sv_waterfall()` bars are now crisp triangle-tipped rectangles (the blob
  arrowheads are gone), near-zero contributions label as `0.00` rather than
  `-0.00`, and the connector lines line up.
* `sv_force()` rebuilt with interlocking arrow blocks, a lower default
  `max_features` (5), greedy non-overlapping labels, and `feature = value`
  sublabels only on the dominant block of each side, so labels no longer
  collide.

# shapvalues 0.1.0

* First release.
* `shapvalues()`: model-agnostic sampling (permutation) SHAP via a single
  batched prediction call, with the marginal-contribution reduction in
  parallel C++ (`RcppParallel`).
* Regression, binary and multiclass classification, and right-censored
  survival (time-dependent SHAP) tasks, auto-detected from a probe prediction.
* Plots: `sv_beeswarm()`, `sv_importance()`, `sv_waterfall()`, `sv_force()`,
  `sv_dependence()`, `sv_surv()`, plus `plot()` / native `ggplot2`.
* Tables: `sv_table()`, `sv_importance_table()`, `as.data.frame()`.
* `inst/benchmarks/bench-shapvalues.R`: reproducible comparison against
  `kernelshap`, `fastshap`, `iml` and `treeshap` (runtime, accuracy vs exact,
  convergence, thread scaling), with precomputed results and figures.
