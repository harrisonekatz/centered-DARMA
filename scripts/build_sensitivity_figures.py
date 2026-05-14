"""
build_sensitivity_figures.py

Cross-reference comparison figures for the centered-vs-raw MA sensitivity
analysis. Reads the per-run output produced by sensitivity_one_reference.R
(tables/rolling_elpd_cov95.csv, tables/rolling_diagnostic_summary.csv) for
each ALR reference, builds:

  1) sensitivity_cumelpd_by_ref.png
  2) sensitivity_divergences_by_ref.png

Assumes the script is run from the bundle root, with a `results/` directory
containing one subdirectory per reference run.

Usage:
  python3 scripts/build_sensitivity_figures.py
"""
import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates

# Paths are relative to the bundle root.
HERE = os.path.dirname(os.path.abspath(__file__))
BUNDLE_ROOT = os.path.abspath(os.path.join(HERE, ".."))
RESULTS_DIR = os.path.join(BUNDLE_ROOT, "results")

RUN_DIRS = {
    "Loans (main)": os.path.join(RESULTS_DIR, "h8_weekly_seed_symmetric_diagnostic_loans_main"),
    "Cash":         os.path.join(RESULTS_DIR, "20260508_142126_h8_sensitivity_ref_cash"),
    "Securities":   os.path.join(RESULTS_DIR, "20260511_141952_h8_sensitivity_ref_securities"),
    "Other":        os.path.join(RESULTS_DIR, "20260508_142710_h8_sensitivity_ref_other"),
}
OUTPUT_DIR = os.path.join(BUNDLE_ROOT, "figs_sensitivity")
os.makedirs(OUTPUT_DIR, exist_ok=True)

COLORS = {
    "Loans (main)": "#1f77b4",
    "Cash":         "#d62728",
    "Securities":   "#2ca02c",
    "Other":        "#9467bd",
}

# ---- Load per-origin ELPD and per-run divergence totals ------------------
elpd_frames = []
div_rows = []
for ref, run_dir in RUN_DIRS.items():
    elpd_path = os.path.join(run_dir, "tables", "rolling_elpd_cov95.csv")
    df = pd.read_csv(elpd_path, parse_dates=["date"]).sort_values("date").reset_index(drop=True)
    df["reference"] = ref
    df["cum_elpd_diff"] = df["ELPD_diff"].cumsum()
    elpd_frames.append(df)

    diag_path = os.path.join(run_dir, "tables", "rolling_diagnostic_summary.csv")
    diag = pd.read_csv(diag_path)
    for _, row in diag.iterrows():
        spec = "Centered" if "Centered" in str(row["model"]) else "Raw"
        div_rows.append({"reference": ref, "spec": spec,
                         "total_divergences": int(row["total_divergences"])})

elpd_long = pd.concat(elpd_frames, ignore_index=True)
div_long  = pd.DataFrame(div_rows)

ref_order = list(RUN_DIRS.keys())

# ---- Figure 1: cumulative ELPD diff by reference ------------------------
fig, ax = plt.subplots(figsize=(8, 4.5))
ax.axhline(0, color="gray", linestyle="--", linewidth=0.8, alpha=0.7)
for ref in ref_order:
    sub = elpd_long[elpd_long["reference"] == ref]
    ax.plot(sub["date"], sub["cum_elpd_diff"],
            label=ref, color=COLORS[ref], linewidth=1.4)
ax.set_xlabel("")
ax.set_ylabel("Cumulative ELPD difference (Centered \u2212 Raw), nats")
ax.legend(loc="upper left", frameon=False, ncol=len(ref_order))
ax.xaxis.set_major_locator(mdates.MonthLocator(interval=3))
ax.xaxis.set_major_formatter(mdates.DateFormatter("%Y-%m"))
plt.setp(ax.get_xticklabels(), rotation=45, ha="right")
ax.grid(True, linestyle=":", alpha=0.4)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
plt.tight_layout()
fig.patch.set_facecolor("white")
ax.set_facecolor("white")
fig.savefig(os.path.join(OUTPUT_DIR, "sensitivity_cumelpd_by_ref.png"),
            dpi=200, bbox_inches="tight", facecolor="white", transparent=False)
plt.close(fig)

# ---- Figure 2: total rolling divergences, Centered vs Raw ---------------
fig, ax = plt.subplots(figsize=(7, 4.5))
x_positions = range(len(ref_order))
width = 0.38

centered_totals = [
    int(div_long[(div_long["reference"] == r) & (div_long["spec"] == "Centered")]["total_divergences"].values[0])
    for r in ref_order
]
raw_totals = [
    int(div_long[(div_long["reference"] == r) & (div_long["spec"] == "Raw")]["total_divergences"].values[0])
    for r in ref_order
]

bars_c = ax.bar([x - width/2 for x in x_positions], centered_totals,
                width=width, label="Centered", color="#1f77b4", edgecolor="white")
bars_r = ax.bar([x + width/2 for x in x_positions], raw_totals,
                width=width, label="Raw", color="#d62728", edgecolor="white")

for bars, totals in [(bars_c, centered_totals), (bars_r, raw_totals)]:
    for bar, total in zip(bars, totals):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(totals)*0.015,
                str(total), ha="center", va="bottom", fontsize=10)

ax.set_xticks(list(x_positions))
ax.set_xticklabels(ref_order)
ax.set_xlabel("ALR reference")
ax.set_ylabel("Total HMC divergent transitions (rolling fits)")
ax.legend(loc="upper left", frameon=False)
ax.grid(True, axis="y", linestyle=":", alpha=0.4)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.set_ylim(0, max(centered_totals + raw_totals) * 1.18)
plt.tight_layout()
fig.patch.set_facecolor("white")
ax.set_facecolor("white")
fig.savefig(os.path.join(OUTPUT_DIR, "sensitivity_divergences_by_ref.png"),
            dpi=200, bbox_inches="tight", facecolor="white", transparent=False)
plt.close(fig)

print("Wrote:")
print(" ", os.path.join(OUTPUT_DIR, "sensitivity_cumelpd_by_ref.png"))
print(" ", os.path.join(OUTPUT_DIR, "sensitivity_divergences_by_ref.png"))
