# Centered-DARMA

Code and locked-window results for:

> Katz, Harrison (2026). *Centered-Innovation MA for Bayesian Dirichlet ARMA:
> Theoretical Equivalence and an Application to Bank-Asset Shares.*

The repository contains the analysis scripts, frozen input data, and the
exact locked-window results reported in Tables 1--4 and Figures 1--6 of the
manuscript.

## Repository layout

```
.
|-- README.md                          # this file
|-- centered_DARMA_main.R              # main analysis (loans reference)
|-- centered_DARMA_sensitivity.R       # four-reference sensitivity analysis
|-- data/                              # frozen FRED CSV snapshots (see Data section)
|   |-- TLAACBW027SBOG.csv             # total assets
|   |-- CASACBW027SBOG.csv             # cash
|   |-- SBCACBW027SBOG.csv             # securities
|   `-- TOTLL.csv                      # loans
|-- results/                           # locked-window run outputs
|   |-- h8_weekly_seed_symmetric_diagnostic_loans_main/   # main analysis
|   |-- 20260508_142126_h8_sensitivity_ref_cash/          # cash sensitivity
|   |-- 20260508_142710_h8_sensitivity_ref_other/         # other sensitivity
|   `-- 20260511_141952_h8_sensitivity_ref_securities/    # securities sensitivity
`-- scripts/
    |-- build_sensitivity_figures.py   # regenerate cross-reference figures (Python)
    `-- build_sensitivity_figures.R    # same, in R
```

## Analysis configuration

Both scripts use:

- **Locked data window**: October 7, 2015 through October 1, 2025
  (T = 522 weekly observations, T_train = 418, T_test = 104). The script
  fails fast if the realized window does not match.
- **Seed-symmetric protocol**: identical random seed at every rolling
  origin for both Centered-MA and Raw-MA, so HMC differences are
  attributable to posterior geometry rather than to seed variability.
- **No auto-refits**: a divergence-triggered refit policy is disabled.
  Both specifications use a fixed sampler at each origin.
- **Rolling sampler**: 2 chains, 1,200 iterations, 600 warmup,
  `adapt_delta = 0.95`, `max_treedepth = 12`, `init = 0`.

## Data

The `data/` directory contains direct CSV downloads from FRED for the four
H.8 series used in the analysis. To regenerate from FRED, download each at:

- `https://fred.stlouisfed.org/graph/fredgraph.csv?id=TLAACBW027SBOG`  (total assets)
- `https://fred.stlouisfed.org/graph/fredgraph.csv?id=CASACBW027SBOG`  (cash)
- `https://fred.stlouisfed.org/graph/fredgraph.csv?id=SBCACBW027SBOG`  (securities)
- `https://fred.stlouisfed.org/graph/fredgraph.csv?id=TOTLL`           (loans)

Each CSV has two columns (`DATE`/`observation_date`, value). Place the
files in `data/` with the names above.

## Reproducing the manuscript results

### Run the main analysis (loans reference)

```r
source("centered_DARMA_main.R")
```

Produces a timestamped directory under `results/` with `config.txt`,
`tables/`, `figs/`, and the merged input `h8_weekly_composition.csv`.

### Run the four-reference sensitivity analysis

```r
# Edit the ref_label variable in centered_DARMA_sensitivity.R, then:
source("centered_DARMA_sensitivity.R")
```

The script writes its output to a timestamped directory under `results/`
named `<timestamp>_h8_sensitivity_ref_<ref_label>` where `<ref_label>` is
one of `cash`, `securities`, `other`, or `loans`.

### Regenerate the cross-reference figures

The two figures comparing all four references are not produced by the R
scripts above (which run one reference at a time). Instead:

```bash
python3 scripts/build_sensitivity_figures.py
```

or

```r
source("scripts/build_sensitivity_figures.R")
```

Both write `sensitivity_cumelpd_by_ref.png` and
`sensitivity_divergences_by_ref.png` into a new `figs_sensitivity/`
directory at the repo root.

## Headline numbers

Across the four ALR references on the locked 104-origin rolling window:

| Reference     | Centered divs | Raw divs | Raw / Cent | Cum ELPD diff |
|---------------|--------------:|---------:|-----------:|---------------:|
| Loans (main)  | 34            | 446      | 13.1       | +0.37          |
| Cash          | 28            | 29       | 1.04       | +0.12          |
| Securities    | 11            | 33       | 3.00       | +0.10          |
| Other         | 60            | 186      | 3.10       | +0.23          |

All numbers can be read directly from
`results/<run>/tables/rolling_diagnostic_summary.csv` and
`results/<run>/tables/rolling_elpd_summary.csv`.

## Key result files per run

Each `results/<run>/tables/` directory contains:

| File | Contents |
|------|----------|
| `rolling_diagnostic_summary.csv` | Per-spec totals: divergences, R-hat, ESS, treedepth hits |
| `rolling_diagnostics.csv` | Per-origin per-spec diagnostics |
| `rolling_elpd_summary.csv` | Cumulative ELPD difference and win counts |
| `rolling_elpd_cov95.csv` | Per-origin Centered/Raw ELPD and componentwise coverage |
| `rolling_per_origin_totals.csv` | Per-origin per-spec total-share RMSE and MAE |
| `metrics_*_fixed.csv` | Fixed-holdout metrics (retained for completeness; not reported in the manuscript) |
