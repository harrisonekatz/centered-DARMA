# Centered-DARMA: Locked-window results bundle

This bundle contains the locked-window rolling-origin results for the
Centered-MA vs Raw-MA comparison reported in:

> Katz, Harrison (2026). *Centered-Innovation MA for Bayesian Dirichlet ARMA:
> Theoretical Equivalence and an Application to Bank-Asset Shares.*

All runs use:

- **Locked data window**: October 7, 2015 through October 1, 2025 (T = 522
  weekly observations, T_train = 418, T_test = 104).
- **Seed-symmetric protocol**: identical random seed at every rolling origin
  for both the Centered-MA and Raw-MA specifications, so that any difference
  in HMC behavior is attributable to posterior geometry, not seed-induced
  sampling variability.
- **No auto-refits**: a divergence-triggered refit policy is disabled. Both
  specifications use a fixed 2-chain, 1,200-iteration sampler with
  `adapt_delta = 0.90` and `max_treedepth = 12`.

## Directory layout

```
.
|-- README.md                                  # this file
|-- results/
|   |-- h8_weekly_seed_symmetric_diagnostic_loans_main/   # main analysis (loans reference)
|   |-- 20260508_142126_h8_sensitivity_ref_cash/          # cash sensitivity run
|   |-- 20260511_141952_h8_sensitivity_ref_securities/    # securities sensitivity run
|   `-- 20260508_142710_h8_sensitivity_ref_other/         # other sensitivity run
|-- scripts/
|   |-- build_sensitivity_figures.py           # regenerates cross-reference figures (Python)
|   `-- build_sensitivity_figures.R            # same, in R
```

Each `results/<run>/` directory contains:

- `config.txt`: the full configuration used for that run (data window, sampler
  settings, priors, seed protocol, sim name).
- `h8_weekly_composition.csv`: the frozen weekly composition input data
  (cash, securities, loans, other shares) for the locked window. Identical
  across all four runs.
- `tables/`: CSV outputs from the run.
- `figs/`: PNG outputs from the run (per-reference cumulative ELPD curve,
  per-origin divergence plot, total-share RMSE plot).

## Key tables in each `results/<run>/tables/` directory

| File | Contents |
|------|----------|
| `rolling_diagnostic_summary.csv` | Per-spec totals: divergences, R-hat, ESS, treedepth hits |
| `rolling_diagnostics.csv` | Per-origin per-spec diagnostics |
| `rolling_elpd_summary.csv` | Cumulative ELPD difference (Centered minus Raw), per-origin mean and SD, win counts |
| `rolling_elpd_cov95.csv` | Per-origin Centered/Raw ELPD and 95% coverage indicator |
| `rolling_per_origin_totals.csv` | Per-origin per-spec total-share RMSE and MAE |
| `metrics_*_fixed.csv` | Fixed-holdout metrics (retained for completeness; not reported in the manuscript) |

## Headline numbers (from `rolling_diagnostic_summary.csv` and
`rolling_elpd_summary.csv` across the four references)

| Reference     | Centered divs | Raw divs | Raw / Cent | Cum ELPD diff |
|---------------|--------------:|---------:|-----------:|---------------:|
| Loans (main)  | 34            | 446      | 13.1       | +0.37          |
| Cash          | 28            | 29       | 1.04       | +0.12          |
| Securities    | 11            | 33       | 3.00       | +0.10          |
| Other         | 60            | 186      | 3.10       | +0.23          |

## Reproducing the cross-reference figures

### Python
```
pip install pandas matplotlib
python3 scripts/build_sensitivity_figures.py
```

### R
```
Rscript scripts/build_sensitivity_figures.R
```

Both scripts write `sensitivity_cumelpd_by_ref.png` and
`sensitivity_divergences_by_ref.png` into a new `figs_sensitivity/` directory
at the bundle root. The numerical inputs come from the CSVs in
`results/<run>/tables/`; no Stan re-fitting is required.

## License

This bundle accompanies the manuscript and its associated public
repository. See the main repository for license terms.
