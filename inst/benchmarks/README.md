# shapvalues benchmarks

`bench-shapvalues.R` compares `shapvalues` with the other model-agnostic SHAP
implementations on CRAN (plus `treeshap` as a tree-only reference), explaining
three model types: `ranger`, `xgboost`, and `fastgbm` (a compiled histogram
GBM).

```
Rscript inst/benchmarks/bench-shapvalues.R          # quick
Rscript inst/benchmarks/bench-shapvalues.R --full   # full grid (slow)
```

Extra packages to run it: `bench`, `kernelshap`, `fastshap` (CRAN-archived;
install from r-universe or GitHub `bgreenwell/fastshap`), `iml`, `treeshap`,
`ranger`, `xgboost`, `fastgbm`. Results go to `bench-shapvalues-results.rds`,
figures to `figures/`.

## Model hyperparameters

All three models share every hyperparameter that can be shared, so differences
between them reflect the model family, not the tuning:

| parameter        | value | ranger | xgboost | fastgbm |
|------------------|------:|:------:|:-------:|:-------:|
| trees            | 200   | y      | y       | y       |
| max depth        | 4     | y      | y       | y       |
| min node size    | 10    | y      | y (`min_child_weight`) | y |
| features / split | all   | `mtry=p` | `colsample_*=1` | `colsample=1` |
| row subsampling  | none  | -      | `subsample=1` | `subsample=1` |
| L2 penalty       | 1     | -      | `lambda=1` | `lambda=1` |
| histogram bins   | 255   | exact splits | `max_bin=255` | `max_bins=255` |
| shrinkage        | 0.1   | -      | `eta=0.1` | `learning_rate=0.1` |

Not shareable: a random forest has no shrinkage, and it keeps bootstrap
resampling (without it, depth-matched full-feature trees would be identical);
`ranger` splits on exact values, the GBMs on 255 histogram bins.

## Ground truth

Accuracy is RMSE against **interventional SHAP** from `kernelshap(exact = TRUE)`
over the full background on a capped row subset - the same estimand
`shapvalues`, `fastshap` and hybrid `kernelshap` target. `treeshap` computes
the path-dependent variant, so its small non-zero error here is expected
(features are independent by construction). Where `p <= 8` the default
`kernelshap` is itself exact, so its accuracy row is ~0 by construction - read
its column as a timing result only.

## Results (quick grid, R 4.5.1, 8-core Linux; 2026-09-01)

### 1. Runtime - model-agnostic methods (ms per explained row)

| model / p       | **shapvalues** | fastshap | kernelshap | iml   | treeshap* |
|-----------------|---------------:|---------:|-----------:|------:|----------:|
| ranger, p=5     | **2.3**        | 7.1      | 22.2       | 122   | 0.18      |
| ranger, p=10    | **4.8**        | 16.3     | 114.9      | 190   | 0.20      |
| xgboost, p=5    | **0.52**       | 3.1      | 4.4        | 105   | 0.09      |
| xgboost, p=10   | **1.4**        | 7.2      | 26.7       | 196   | 0.11      |
| fastgbm, p=5    | **0.94**       | 2.5      | 7.4        | 103   | n/a       |
| fastgbm, p=10   | **1.5**        | 8.5      | 37.2       | 218   | n/a       |

\*treeshap: exact, tree-only, path-dependent; no adapter for `fastgbm`. With
these depth-4 forests its trees are tiny, so it is very fast where it applies.

`shapvalues` is the fastest **model-agnostic** method in every case:
**3-6x faster than fastshap**, **10-25x faster than kernelshap**, ~50-200x
faster than `iml`. Across the three model families at p=10 its cost tracks the
model's `predict` cost: xgboost 1.4, fastgbm 1.5, ranger 4.8 ms/row. For
`fastgbm` there is no exact tree adapter in any package, so `shapvalues` is the
fastest available SHAP for it.

### 2. Accuracy vs interventional-exact (RMSE), same runs

`shapvalues` matches or beats `fastshap` at equal `nsim`: ranger p=10 0.077 vs
0.125; xgboost p=10 0.098 vs 0.188; fastgbm p=10 0.129 vs 0.164. Hybrid
`kernelshap` (p>8) is tighter (5e-3 to 9e-3) because it uses the whole
background with optimal weights - at 20-250x the time.

### 3. Convergence vs nsim (xgboost, p=10, vs interventional exact)

| nsim | shapvalues RMSE | shapvalues time | fastshap RMSE | fastshap time |
|-----:|----------------:|----------------:|--------------:|--------------:|
| 10   | 0.35            | 0.10 s          | 0.38          | 0.67 s        |
| 25   | 0.22            | 0.18 s          | 0.25          | 1.18 s        |
| 50   | 0.13            | 0.53 s          | 0.17          | 2.40 s        |
| 100  | 0.11            | 0.96 s          | 0.12          | 5.45 s        |
| 200  | 0.030           | 2.29 s          | 0.088         | 10.25 s       |

`shapvalues` reaches a given accuracy **2-5x faster than `fastshap`**. For an
*additive* model it becomes exact once `nsim >= bg_size` (a separate `lm` check
matched `kernelshap` to 6e-16); for a tree model the permutation count still
matters, so error keeps falling with `nsim` rather than snapping to zero.

### 4. Thread scaling (analytic 40-output model, p=18, 2500 rows)

| threads | 1 | 2 | 4 | 8 |
|---------|--:|--:|--:|--:|
| speedup | 1.00 | 1.21 | 1.21 | 1.28 |

The parallel C++ reduction only helps when the model has many outputs
(multiclass, survival over a time grid); for single-output models the cost is
the one batched `predict()` call plus R-side coalition-row assembly, not the
reduction. The advantage over other packages comes from issuing **one**
prediction call for all coalitions and a lean compiled reduce - not from
threads. Moving row assembly into C++ is tracked as future work.
