# -*- coding: utf-8 -*-
"""
Generic analysis of omegasome diameter over time
Version with configurable y-limits, QC filters, PNG export, per-track duration
summaries, manual CIs (bands) for time courses, and Prism/Excel-friendly exports.

Script version: Analysis_new_V7_2026-02-26
"""

from __future__ import division
import shutil
from datetime import datetime
from pathlib import Path

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import glob

# -------------------------------------------------------------------
# SCRIPT VERSION
# -------------------------------------------------------------------

SCRIPT_VERSION = "Analysis_new_V7_2026-02-26"

# -------------------------------------------------------------------
# USER CONFIG
# -------------------------------------------------------------------

# Root directory that contains <date>/<strain>/Measurements_2/*.csv
DATA_ROOT = r"C:\DataAnalysis\DFCP1\Analysis"

# Subdirectories (e.g. dates / experiments)
DIRNAMES = ["251016", "251024", "251028"]  # ["251028"] #

# Strain IDs as stored in the CSV ("Strain" column)
# and their human-readable names
STRAIN_IDS = [243, 244, 245, 246]
STRAIN_NAME_MAP = {
    243: "WT",
    244: "244",
    245: "245",
    246: "246",
}

# Which strain to use as reference for vertical time lines
REFERENCE_STRAIN_ID = 243  # must be in STRAIN_IDS

# Order in plots
ORDER_STRAIN_NUM = STRAIN_IDS
ORDER_STRAIN_NAME = [STRAIN_NAME_MAP[s] for s in ORDER_STRAIN_NUM]

# Folder names for strains in the filesystem
CELL_LINES = [str(s) for s in STRAIN_IDS]

# Physical pixel size and frame interval
PIXEL_SIZE = 0.08      # µm per pixel
FRAME_INTERVAL = 2.0   # seconds per frame

# Size bins (MaxDiameter in µm)
SIZE_BIN_EDGES = [0, 1, 1.5, 2, np.inf]
SIZE_BIN_LABELS = ["<1", "1-1.5", "1.5-2", ">2"]

# Lifetime bins (MaxTime in s)
TIME_BIN_EDGES = [0, 100, 200, 300, 400, 500, 600, np.inf]
TIME_BIN_LABELS = ["<100", "100-200", "200-300", "300-400",
                   "400-500", "500-600", ">600"]

# Thresholds for constriction analysis (in µm)
CONSTRICTION_THRESHOLDS = [1.0, 0.8, 0.6]

# Where to save plots
PLOT_ROOT = r"C:\DataAnalysis\DFCP1\Analysis\plots"

# Also save PNG copies of all plots
SAVE_PNG = True

# How to draw the vertical reference line
# options: "max", "median", "quantile"
REF_TIME_METHOD = "max"
REF_TIME_QUANTILE = 0.95   # used only if REF_TIME_METHOD == "quantile"

# -------------------------------------------------------------------
# Y-LIMIT CONFIG
# -------------------------------------------------------------------
# mode: "fixed" -> use fixed tuple; "dynamic" -> factor × max(column)

YLIM_MODE_DIAMETER = "dynamic"      # mean diameter curves
YLIM_MODE_GROUPED_SIZE = "dynamic"  # grouped by size
YLIM_MODE_GROUPED_DURATION = "dynamic"
YLIM_MODE_TRACK_DURATION = "fixed"
YLIM_MODE_TRACK_MAX_DIAM = "dynamic"

# Fixed baselines (lower bound + default upper bound)
YLIM_DIAMETER_FIXED = (0, 2.0)
YLIM_GROUPED_SIZE_FIXED = (0, 2.0)
YLIM_GROUPED_DURATION_FIXED = (0, 1.7)
YLIM_TRACK_DURATION_FIXED = (0, 900)
YLIM_TRACK_MAX_DIAM_FIXED = (0, 3.5)

# Factors for dynamic mode (upper = factor × max(column))
YLIM_DIAMETER_FACTOR = 1.1
YLIM_GROUPED_SIZE_FACTOR = 1.1
YLIM_GROUPED_DURATION_FACTOR = 1.1
YLIM_TRACK_MAX_DIAM_FACTOR = 1.1

# -------------------------------------------------------------------
# QC CONFIG FOR PLOTTED CURVES
# -------------------------------------------------------------------
# Minimum number of movies per bin (Strain, size / time_class)
MIN_MOVIES_PER_SIZE_BIN = 2
MIN_MOVIES_PER_TIMECLASS_BIN = 2

# Minimum number of tracks contributing at a given timepoint
MIN_TRACKS_PER_TIME_SIZE = 2
MIN_TRACKS_PER_TIME_TIMECLASS = 2

# Apply QC also to the per-track summary plots?
APPLY_QC_TO_TRACK_SUMMARIES = True

# -------------------------------------------------------------------
# SETUP: PLOTS FOLDER, STYLE, AND VERSION INFO
# -------------------------------------------------------------------

print(f"Running analysis script version: {SCRIPT_VERSION}")

plot_root_path = Path(PLOT_ROOT)
plot_root_path.mkdir(parents=True, exist_ok=True)

# Empty plots folder at each run
for child in plot_root_path.iterdir():
    if child.is_file() or child.is_symlink():
        child.unlink()
    elif child.is_dir():
        shutil.rmtree(child)

# Write a small version / provenance file
version_md = plot_root_path / "analysis_version.md"
with version_md.open("w", encoding="utf-8") as f:
    f.write("# Analysis script version\n\n")
    f.write(f"- Script version: {SCRIPT_VERSION}\n")
    f.write(f"- Generated at: {datetime.now().isoformat()}\n")
    f.write(f"- DATA_ROOT: {DATA_ROOT}\n")
    f.write(f"- DIRNAMES: {', '.join(DIRNAMES)}\n")
    f.write(f"- STRAIN_IDS: {', '.join(str(s) for s in STRAIN_IDS)}\n")

plt.rcParams["figure.figsize"] = (3.5, 2.5)
sns.set_style("white")

PAL_LINE = sns.color_palette("bright", n_colors=len(STRAIN_IDS))
PAL_STRIP = sns.color_palette("bright", n_colors=len(STRAIN_IDS))

ORDER_SIZE = SIZE_BIN_LABELS
ORDER_TIME_CLASS = TIME_BIN_LABELS

# Global CI factor for 95% CI (normal approx)
CI_FACTOR = 1.96

# -------------------------------------------------------------------
# HELPER: COMPUTE YLIM
# -------------------------------------------------------------------

def compute_ylim(data, column, fixed, mode="fixed", factor=1.1):
    """
    Decide y-limits for a plot.

    data:   DataFrame used for the plot
    column: column name containing the y-data
    fixed:  (ymin, ymax) default tuple
    mode:   "fixed" or "dynamic"
    factor: multiplier for dynamic upper limit
    """
    if mode == "fixed":
        return fixed
    elif mode == "dynamic":
        if column not in data.columns:
            return fixed
        y_max = data[column].max()
        if pd.isna(y_max):
            return fixed
        upper = y_max * factor
        # keep lower bound from fixed
        return (fixed[0], upper)
    else:
        # fallback
        return fixed

# -------------------------------------------------------------------
# HELPER FOR FACETGRID OVERLAYS
# -------------------------------------------------------------------

def add_ref_vline_to_facets(grid, ref_time):
    """Add a vertical reference line with the same x position in all facets."""
    if ref_time is None:
        return
    try:
        axes = grid.axes.flat
    except AttributeError:
        axes = [grid.ax]
    for ax in axes:
        ax.axvline(ref_time)

def plot_mean_ci(data, x_col, mean_col, lower_col, upper_col, color=None, **kwargs):
    """
    Helper for FacetGrid: plot mean line + CI band with a visible edge.
    """
    ax = plt.gca()

    data = (
        data.sort_values(x_col)
            .replace([np.inf, -np.inf], np.nan)
            .dropna(subset=[mean_col, lower_col, upper_col])
    )
    if data.empty:
        return

    if color is None:
        color = "black"

    # CI band with a faint edge
    ax.fill_between(
        data[x_col],
        data[lower_col],
        data[upper_col],
        facecolor=color,
        edgecolor=color,      # <- edge back
        linewidth=0.7,
        alpha=0.2,
        zorder=1,
    )

    # Mean line on top as clear boundary
    ax.plot(
        data[x_col],
        data[mean_col],
        color=color,
        linewidth=1.2,
        alpha=1.0,
        zorder=2,
    )
    
    
# -------------------------------------------------------------------
# LOAD AND PREP DATA
# -------------------------------------------------------------------

all_tracks = []

for dir_idx, directory in enumerate(DIRNAMES):
    print("Processing directory:", directory)
    for cellline in CELL_LINES:
        pattern = f"{DATA_ROOT}/{directory}/{cellline}/Measurements_2/*.csv"
        files = glob.glob(pattern)
        if not files:
            print("  No files for", pattern)
            continue

        for csvfile in files:
            df = pd.read_csv(csvfile, decimal=".")
            # Time axis from Timepoint in seconds
            df["Time (s)"] = (df["Timepoint"] - df["Timepoint"].iloc[0]) * FRAME_INTERVAL
            df["Experiment"] = dir_idx + 1

            # Index (row) of max diameter within this track
            max_idx = df["Max Diameter"].idxmax()

            # Proper time of maximum in seconds from the Time (s) column
            time_of_max = df.loc[max_idx, "Time (s)"]
            df["Time of maximum (s)"] = time_of_max

            # Constriction time = duration from max until end of track (in seconds)
            df["Tracklength from max (s)"] = df["Time (s)"].max() - time_of_max

            all_tracks.append(df)

# Concatenate all tracks
df1 = pd.concat(all_tracks).reset_index(drop=True)

# Ensure Strain is numeric so it matches STRAIN_IDS / REFERENCE_STRAIN_ID
df1["Strain"] = pd.to_numeric(df1["Strain"], errors="coerce")

# Convert to physical units
df1["Diameter (µm)"] = df1["Max Diameter"] * PIXEL_SIZE
df1["Area (µm^2)"] = df1["Area"] * (PIXEL_SIZE ** 2)

# Map strains to names
df1["Strainname"] = df1["Strain"].map(STRAIN_NAME_MAP)

print("Strain mapping:")
for k, v in STRAIN_NAME_MAP.items():
    print(f"  {k}: {v}")

# -------------------------------------------------------------------
# MEAN DIAMETER OVER TIME (per strain) – line + CI band
# -------------------------------------------------------------------

plt.rcParams["figure.figsize"] = (5, 4)

df_mean_time = (
    df1
    .groupby(["Strain", "Strainname", "Time (s)"])["Diameter (µm)"]
    .agg(["mean", "std", "count"])
    .reset_index()
)

df_mean_time["sem"] = df_mean_time["std"] / np.sqrt(df_mean_time["count"])
df_mean_time["ci_lower"] = df_mean_time["mean"] - CI_FACTOR * df_mean_time["sem"]
df_mean_time["ci_upper"] = df_mean_time["mean"] + CI_FACTOR * df_mean_time["sem"]

fig, ax = plt.subplots()

for strain, color in zip(ORDER_STRAIN_NUM, PAL_LINE):
    sub = df_mean_time[df_mean_time["Strain"] == strain].sort_values("Time (s)")
    if sub.empty:
        continue
    sub = (
        sub.replace([np.inf, -np.inf], np.nan)
           .dropna(subset=["mean", "ci_lower", "ci_upper"])
    )
    if sub.empty:
        continue

    label = STRAIN_NAME_MAP.get(strain, str(strain))

    ax.fill_between(
        sub["Time (s)"],
        sub["ci_lower"],
        sub["ci_upper"],
        facecolor=color,
        edgecolor=color,
        linewidth=0.7,
        alpha=0.2,
        zorder=1,
    )
    ax.plot(
        sub["Time (s)"],
        sub["mean"],
        color=color,
        linewidth=1.5,
        label=label,
        zorder=2,
    )

sns.despine()
ax.set_xlabel("Time (s)")
ax.set_ylabel("Diameter (µm)")
ax.legend()

ylim_mean = compute_ylim(
    df1, "Diameter (µm)", YLIM_DIAMETER_FIXED,
    mode=YLIM_MODE_DIAMETER, factor=YLIM_DIAMETER_FACTOR
)
ax.set_ylim(*ylim_mean)

plt.savefig(f"{PLOT_ROOT}/mean_diameter_over_time.svg")
if SAVE_PNG:
    plt.savefig(f"{PLOT_ROOT}/mean_diameter_over_time.png", dpi=300)
plt.show()

# -------------------------------------------------------------------
# SPLIT LIFETIME BY MAXIMUM SIZE
# -------------------------------------------------------------------

df2 = df1.groupby(["Strain", "Experiment", "Filename"]).max().reset_index()
df2 = df2.rename(
    columns={"Diameter (µm)": "MaxDiameter (µm)", "Time (s)": "MaxTime (s)"}
)

df2_simple = df2[["Filename", "MaxDiameter (µm)", "MaxTime (s)"]]
df1_1 = df1.merge(df2_simple, on="Filename", how="left")

# Size class
df1_1["size"] = pd.cut(
    df1_1["MaxDiameter (µm)"],
    bins=SIZE_BIN_EDGES,
    labels=SIZE_BIN_LABELS,
    include_lowest=True,
    right=True,
)

# ---------------- QC FILTERS FOR SIZE-GROUPED CURVES ----------------

track_size = df1_1[["Filename", "Strainname", "size"]].drop_duplicates()
bin_counts_size = (
    track_size
    .groupby(["Strainname", "size"])["Filename"]
    .nunique()
    .reset_index(name="n_movies")
)
print("\nMovies per (Strainname, size) bin:")
print(bin_counts_size.sort_values(["Strainname", "size"]))

valid_bins_size = bin_counts_size[
    bin_counts_size["n_movies"] >= MIN_MOVIES_PER_SIZE_BIN
]

df1_size = df1_1.merge(
    valid_bins_size[["Strainname", "size"]],
    on=["Strainname", "size"],
    how="inner"
)

per_time_counts_size = (
    df1_size
    .groupby(["Strainname", "size", "Time (s)"])["Filename"]
    .nunique()
    .reset_index(name="n_tracks")
)

df1_size = df1_size.merge(
    per_time_counts_size,
    on=["Strainname", "size", "Time (s)"],
    how="left"
)

df1_size_plot = df1_size[df1_size["n_tracks"] >= MIN_TRACKS_PER_TIME_SIZE].copy()

df_size_curves = (
    df1_size_plot
    .groupby(["Strainname", "Strain", "size", "Time (s)"])["Diameter (µm)"]
    .agg(["mean", "std", "count"])
    .reset_index()
)
df_size_curves["sem"] = df_size_curves["std"] / np.sqrt(df_size_curves["count"])
df_size_curves["ci_lower"] = df_size_curves["mean"] - CI_FACTOR * df_size_curves["sem"]
df_size_curves["ci_upper"] = df_size_curves["mean"] + CI_FACTOR * df_size_curves["sem"]

# ---------------- REFERENCE TIME AFTER QC (SIZE-GROUPED) -----------

if REFERENCE_STRAIN_ID in df1_size_plot["Strain"].unique():
    df_ref_size = df1_size_plot[df1_size_plot["Strain"] == REFERENCE_STRAIN_ID]
else:
    df_ref_size = pd.DataFrame()

if df_ref_size.empty:
    ref_time_size = None
else:
    track_ends_size = df_ref_size.groupby("Filename")["Time (s)"].max()

    if REF_TIME_METHOD == "max":
        ref_time_size = track_ends_size.max()
    elif REF_TIME_METHOD == "median":
        ref_time_size = track_ends_size.median()
    elif REF_TIME_METHOD == "quantile":
        ref_time_size = track_ends_size.quantile(REF_TIME_QUANTILE)
    else:
        ref_time_size = track_ends_size.max()

    coverage_size = (track_ends_size >= ref_time_size).mean()
    print(
        f"[SIZE] Using ref_time_size = {ref_time_size:.1f} s "
        f"(method={REF_TIME_METHOD}, coverage={coverage_size:.2f} of WT movies AFTER QC)"
    )

ylim_grouped_size = compute_ylim(
    df1_size_plot, "Diameter (µm)", YLIM_GROUPED_SIZE_FIXED,
    mode=YLIM_MODE_GROUPED_SIZE, factor=YLIM_GROUPED_SIZE_FACTOR
)

# Facet by strain, color by size – mean + CI, with horizontal + vertical line
g = sns.FacetGrid(
    df_size_curves,
    col="Strainname",
    hue="size",
    col_order=ORDER_STRAIN_NAME,
    hue_order=SIZE_BIN_LABELS,
    sharex=True,
    sharey=True,
    palette="bright",
)
g.map_dataframe(
    plot_mean_ci,
    x_col="Time (s)",
    mean_col="mean",
    lower_col="ci_lower",
    upper_col="ci_upper",
)
g.map(plt.axhline, y=0.5)
add_ref_vline_to_facets(g, ref_time_size)
for ax in g.axes.flat:
    ax.set_ylim(*ylim_grouped_size)
g.savefig(f"{PLOT_ROOT}/grouped_by_diameter_includinglines.svg")
if SAVE_PNG:
    g.savefig(f"{PLOT_ROOT}/grouped_by_diameter_includinglines.png", dpi=300)
plt.show()

# Same data, no horizontal line
g = sns.FacetGrid(
    df_size_curves,
    col="Strainname",
    hue="size",
    col_order=ORDER_STRAIN_NAME,
    hue_order=SIZE_BIN_LABELS,
    sharex=True,
    sharey=True,
    palette="bright",
)
g.map_dataframe(
    plot_mean_ci,
    x_col="Time (s)",
    mean_col="mean",
    lower_col="ci_lower",
    upper_col="ci_upper",
)
add_ref_vline_to_facets(g, ref_time_size)
for ax in g.axes.flat:
    ax.set_ylim(*ylim_grouped_size)
g.savefig(f"{PLOT_ROOT}/grouped_by_diameter_includingline_xaxis.svg")
if SAVE_PNG:
    g.savefig(f"{PLOT_ROOT}/grouped_by_diameter_includingline_xaxis.png", dpi=300)
plt.show()

# Horizontal line only (no vertical line)
g = sns.FacetGrid(
    df_size_curves,
    col="Strainname",
    hue="size",
    col_order=ORDER_STRAIN_NAME,
    hue_order=SIZE_BIN_LABELS,
    sharex=True,
    sharey=True,
    palette="bright",
)
g.map_dataframe(
    plot_mean_ci,
    x_col="Time (s)",
    mean_col="mean",
    lower_col="ci_lower",
    upper_col="ci_upper",
)
g.map(plt.axhline, y=0.5)
for ax in g.axes.flat:
    ax.set_ylim(*ylim_grouped_size)
g.savefig(f"{PLOT_ROOT}/grouped_by_diameter_includingline_yaxis.svg")
if SAVE_PNG:
    g.savefig(f"{PLOT_ROOT}/grouped_by_diameter_includingline_yaxis.png", dpi=300)
plt.show()

# -------------------------------------------------------------------
# CONSTRICTION DURATION (MEAN CURVE, PRINTED ONLY)
# -------------------------------------------------------------------

print("Calculating the duration of the constriction (mean curves)...")
df12 = df1_1.groupby(["Strain", "Time (s)"]).mean(numeric_only=True)

for threshold in CONSTRICTION_THRESHOLDS:
    print(f"Threshold: {threshold}")
    for strain in STRAIN_IDS:
        if strain not in df12.index.levels[0]:
            print(f"  Strain {strain} missing in data, skipping.")
            continue

        df13 = df12.xs(strain, level="Strain")  # index = Time (s)
        if df13.empty:
            print(f"  Strain {strain}: no data after grouping, skipping.")
            continue

        maximum = df13["Diameter (µm)"].max()
        timepoint = df13["Diameter (µm)"].idxmax()

        print("___________________")
        print(f"Strain: {strain} ({STRAIN_NAME_MAP.get(strain, 'NA')})")
        print(f"Maximum diameter (mean): {maximum}")
        print(f"Time of maximum diameter: {timepoint}")

        df_after = df13[df13.index > timepoint]
        df_below = df_after[df_after["Diameter (µm)"] < threshold]

        if df_below.empty:
            print("No time point below threshold – no completed constriction.")
            print("___________________")
            continue

        constricted = df_below.index[0]
        print(f"Time of completed constriction: {constricted}")
        print(f"Duration of constriction: {constricted - timepoint}")
        print("___________________")

# -------------------------------------------------------------------
# GROUP BY DURATION OF MAXIMUM (TIME CLASS)
# -------------------------------------------------------------------

df1_1["time_class"] = pd.cut(
    df1_1["MaxTime (s)"],
    bins=TIME_BIN_EDGES,
    labels=TIME_BIN_LABELS,
    include_lowest=True,
    right=True,
)

# ------------- QC FILTERS FOR TIME-CLASS-GROUPED CURVES -------------

track_time = df1_1[["Filename", "Strainname", "time_class"]].drop_duplicates()
bin_counts_time = (
    track_time
    .groupby(["Strainname", "time_class"])["Filename"]
    .nunique()
    .reset_index(name="n_movies")
)
print("\nMovies per (Strainname, time_class) bin:")
print(bin_counts_time.sort_values(["Strainname", "time_class"]))

valid_bins_time = bin_counts_time[
    bin_counts_time["n_movies"] >= MIN_MOVIES_PER_TIMECLASS_BIN
]

df1_time = df1_1.merge(
    valid_bins_time[["Strainname", "time_class"]],
    on=["Strainname", "time_class"],
    how="inner"
)

per_time_counts_time = (
    df1_time
    .groupby(["Strainname", "time_class", "Time (s)"])["Filename"]
    .nunique()
    .reset_index(name="n_tracks")
)

df1_time = df1_time.merge(
    per_time_counts_time,
    on=["Strainname", "time_class", "Time (s)"],
    how="left"
)

df1_time_plot = df1_time[df1_time["n_tracks"] >= MIN_TRACKS_PER_TIME_TIMECLASS].copy()

df_time_curves = (
    df1_time_plot
    .groupby(["Strainname", "Strain", "time_class", "Time (s)"])["Diameter (µm)"]
    .agg(["mean", "std", "count"])
    .reset_index()
)
df_time_curves["sem"] = df_time_curves["std"] / np.sqrt(df_time_curves["count"])
df_time_curves["ci_lower"] = df_time_curves["mean"] - CI_FACTOR * df_time_curves["sem"]
df_time_curves["ci_upper"] = df_time_curves["mean"] + CI_FACTOR * df_time_curves["sem"]

# ------------- REFERENCE TIME AFTER QC (TIME-CLASS-GROUPED) --------

if REFERENCE_STRAIN_ID in df1_time_plot["Strain"].unique():
    df_ref_timeclass = df1_time_plot[df1_time_plot["Strain"] == REFERENCE_STRAIN_ID]
else:
    df_ref_timeclass = pd.DataFrame()

if df_ref_timeclass.empty:
    ref_time_timeclass = None
else:
    track_ends_timeclass = df_ref_timeclass.groupby("Filename")["Time (s)"].max()

    if REF_TIME_METHOD == "max":
        ref_time_timeclass = track_ends_timeclass.max()
    elif REF_TIME_METHOD == "median":
        ref_time_timeclass = track_ends_timeclass.median()
    elif REF_TIME_METHOD == "quantile":
        ref_time_timeclass = track_ends_timeclass.quantile(REF_TIME_QUANTILE)
    else:
        ref_time_timeclass = track_ends_timeclass.max()

    coverage_timeclass = (track_ends_timeclass >= ref_time_timeclass).mean()
    print(
        f"[TIMECLASS] Using ref_time_timeclass = {ref_time_timeclass:.1f} s "
        f"(method={REF_TIME_METHOD}, "
        f"coverage={coverage_timeclass:.2f} of WT movies AFTER QC)"
    )

ylim_grouped_duration = compute_ylim(
    df1_time_plot, "Diameter (µm)", YLIM_GROUPED_DURATION_FIXED,
    mode=YLIM_MODE_GROUPED_DURATION, factor=YLIM_GROUPED_DURATION_FACTOR
)

g = sns.FacetGrid(
    df_time_curves,
    col="Strainname",
    hue="time_class",
    col_order=ORDER_STRAIN_NAME,
    hue_order=TIME_BIN_LABELS,
    sharex=True,
    sharey=True,
    palette="bright",
)
g.map_dataframe(
    plot_mean_ci,
    x_col="Time (s)",
    mean_col="mean",
    lower_col="ci_lower",
    upper_col="ci_upper",
)
for ax in g.axes.flat:
    ax.set_ylim(*ylim_grouped_duration)
g.savefig(f"{PLOT_ROOT}/grouped_by_duration.svg")
if SAVE_PNG:
    g.savefig(f"{PLOT_ROOT}/grouped_by_duration.png", dpi=300)
plt.show()

g = sns.FacetGrid(
    df_time_curves,
    col="Strainname",
    hue="time_class",
    col_order=ORDER_STRAIN_NAME,
    hue_order=TIME_BIN_LABELS,
    sharex=True,
    sharey=True,
    palette="bright",
)
g.map_dataframe(
    plot_mean_ci,
    x_col="Time (s)",
    mean_col="mean",
    lower_col="ci_lower",
    upper_col="ci_upper",
)
g.map(plt.axhline, y=0.5)
add_ref_vline_to_facets(g, ref_time_timeclass)
for ax in g.axes.flat:
    ax.set_ylim(*ylim_grouped_duration)
g.savefig(f"{PLOT_ROOT}/grouped_by_duration_lines.svg")
if SAVE_PNG:
    g.savefig(f"{PLOT_ROOT}/grouped_by_duration_lines.png", dpi=300)
plt.show()

# -------------------------------------------------------------------
# DEFINE "GOOD" TRACKS FOR SUMMARY PLOTS (OPTIONAL QC)
# -------------------------------------------------------------------

if APPLY_QC_TO_TRACK_SUMMARIES:
    GOOD_FILES_DURATION = set(df1_time_plot["Filename"].unique())
    GOOD_FILES_DIAMETER = set(df1_size_plot["Filename"].unique())
    print(
        f"\nTracks kept for duration summaries after QC: {len(GOOD_FILES_DURATION)}"
    )
    print(
        f"Tracks kept for diameter summaries after QC: {len(GOOD_FILES_DIAMETER)}"
    )
else:
    GOOD_FILES_DURATION = None
    GOOD_FILES_DIAMETER = None

# -------------------------------------------------------------------
# PER-TRACK SUMMARY: DURATIONS AND MAX DIAMETER
# -------------------------------------------------------------------

df_tracks = df1.groupby(["Strain", "Experiment", "Filename"]).max().reset_index()
df_tracks["Strain name"] = df_tracks["Strain"].map(STRAIN_NAME_MAP)

df_tracks["Constriction time (s)"] = df_tracks["Tracklength from max (s)"]
df_tracks["Track length total (s)"] = df_tracks["Time (s)"]
df_tracks["Time to max (s)"] = df_tracks["Time of maximum (s)"]

df_tracks["Diameter (µm)"] = df_tracks["Max Diameter"] * PIXEL_SIZE

if GOOD_FILES_DURATION is not None:
    df_tracks_dur = df_tracks[df_tracks["Filename"].isin(GOOD_FILES_DURATION)].copy()
    print(
        f"Per-track duration summaries: using {len(df_tracks_dur)} tracks "
        f"out of {len(df_tracks)} total after time-class QC."
    )
else:
    df_tracks_dur = df_tracks

def plot_strip_with_errorbars(df, value_col, filename_base, ylabel):
    """
    Stripplot + mean ± SEM errorbars per strain (these are the *only* error bars).
    """
    summary = (
        df.groupby("Strain name")[value_col]
        .agg(["mean", "std", "count"])
        .reset_index()
    )
    summary["sem"] = summary["std"] / np.sqrt(summary["count"])

    fig, ax = plt.subplots()

    sns.stripplot(
        x="Strain name",
        y=value_col,
        data=df,
        zorder=2,
        order=ORDER_STRAIN_NAME,
        palette=PAL_STRIP,
        ax=ax,
    )

    x_positions = np.arange(len(ORDER_STRAIN_NAME))
    means = []
    sems = []
    for name in ORDER_STRAIN_NAME:
        row = summary[summary["Strain name"] == name]
        if row.empty:
            means.append(np.nan)
            sems.append(np.nan)
        else:
            means.append(row["mean"].iloc[0])
            sems.append(row["sem"].iloc[0])

    ax.errorbar(
        x_positions,
        means,
        yerr=sems,
        fmt="o",
        color="black",
        zorder=3,
        capsize=4,
    )

    ax.set_xticks(x_positions)
    ax.set_xticklabels(ORDER_STRAIN_NAME)
    ax.set_ylabel(ylabel)
    sns.despine()
    plt.savefig(f"{PLOT_ROOT}/{filename_base}.svg")
    if SAVE_PNG:
        plt.savefig(f"{PLOT_ROOT}/{filename_base}.png", dpi=300)
    plt.show()


# 1) Constriction time: from max diameter to end of track
plot_strip_with_errorbars(
    df_tracks_dur,
    value_col="Constriction time (s)",
    filename_base="ConstrictionTime_max_to_end",
    ylabel="Constriction time (s)",
)

# 2) Overall track length: from track start to end
plot_strip_with_errorbars(
    df_tracks_dur,
    value_col="Track length total (s)",
    filename_base="TrackLength_total",
    ylabel="Track length total (s)",
)

# 3) Time from beginning of track to maximum diameter
plot_strip_with_errorbars(
    df_tracks_dur,
    value_col="Time to max (s)",
    filename_base="TimeToMax_from_start",
    ylabel="Time to max (s)",
)

if GOOD_FILES_DIAMETER is not None:
    df_tracks_diam = df_tracks[df_tracks["Filename"].isin(GOOD_FILES_DIAMETER)].copy()
    print(
        f"Per-track diameter summaries: using {len(df_tracks_diam)} tracks "
        f"out of {len(df_tracks)} total after size QC."
    )
else:
    df_tracks_diam = df_tracks

plot_strip_with_errorbars(
    df_tracks_diam,
    value_col="Diameter (µm)",
    filename_base="Diameter_01",
    ylabel="Max diameter (µm)",
)

# -------------------------------------------------------------------
# QC SUMMARY (TRACK COUNTS)
# -------------------------------------------------------------------

print("\n================ QC SUMMARY ================")
total_tracks = len(df_tracks)

if GOOD_FILES_DURATION is not None:
    n_dur_used = len(GOOD_FILES_DURATION)
    n_dur_excl = total_tracks - n_dur_used
    print("Time-class (duration-based QC):")
    print(f"  Total tracks:               {total_tracks}")
    print(f"  Used in duration analyses:  {n_dur_used}")
    print(f"  Excluded by time-class QC:  {n_dur_excl}")
else:
    print("Time-class QC for summaries:  DISABLED")

if GOOD_FILES_DIAMETER is not None:
    n_size_used = len(GOOD_FILES_DIAMETER)
    n_size_excl = total_tracks - n_size_used
    print("\nSize-class (max-diameter QC):")
    print(f"  Total tracks:              {total_tracks}")
    print(f"  Used in size analyses:     {n_size_used}")
    print(f"  Excluded by size QC:       {n_size_excl}")
else:
    print("\nSize QC for summaries:       DISABLED")

print("============================================")

# -------------------------------------------------------------------
# EXPORTS FOR GRAPHPRISM / EXCEL (semicolon CSV, wide format)
# -------------------------------------------------------------------

export_root = plot_root_path / "exports"
export_root.mkdir(parents=True, exist_ok=True)

def make_per_track_wide_df(df, metric_col):
    """
    Build a wide per-track DataFrame for a single metric.

    Wide format:
      - columns: '<Strain>_Exp<Experiment>'
      - rows: track index within that (Experiment, Strain) group
    """
    tmp = df[["Experiment", "Strain name", metric_col]].copy()
    tmp["TrackIndex"] = tmp.groupby(["Experiment", "Strain name"]).cumcount()

    wide = tmp.pivot(
        index="TrackIndex",
        columns=["Experiment", "Strain name"],
        values=metric_col,
    )

    wide.columns = [
        f"{strain}_Exp{exp}" for (exp, strain) in wide.columns
    ]

    wide = wide.reset_index()
    return wide

per_track_constriction = make_per_track_wide_df(
    df_tracks_dur, "Constriction time (s)"
)
per_track_tracklen = make_per_track_wide_df(
    df_tracks_dur, "Track length total (s)"
)
per_track_timemax = make_per_track_wide_df(
    df_tracks_dur, "Time to max (s)"
)
per_track_diam = make_per_track_wide_df(
    df_tracks_diam, "Diameter (µm)"
)

per_track_constriction.to_csv(
    export_root / "PerTrack_ConstrictionTime_max_to_end_wide.csv",
    sep=";",
    decimal=",",
    index=False,
)
per_track_tracklen.to_csv(
    export_root / "PerTrack_TrackLength_total_wide.csv",
    sep=";",
    decimal=",",
    index=False,
)
per_track_timemax.to_csv(
    export_root / "PerTrack_TimeToMax_from_start_wide.csv",
    sep=";",
    decimal=",",
    index=False,
)
per_track_diam.to_csv(
    export_root / "PerTrack_MaxDiameter_wide.csv",
    sep=";",
    decimal=",",
    index=False,
)

per_exp_dur = (
    df_tracks_dur
    .groupby(["Experiment", "Strain name"], as_index=False)[
        ["Constriction time (s)",
         "Track length total (s)",
         "Time to max (s)"]
    ]
    .mean()
)

constriction_wide = per_exp_dur.pivot(
    index="Experiment",
    columns="Strain name",
    values="Constriction time (s)",
).reset_index()

tracklen_wide = per_exp_dur.pivot(
    index="Experiment",
    columns="Strain name",
    values="Track length total (s)",
).reset_index()

timemax_wide = per_exp_dur.pivot(
    index="Experiment",
    columns="Strain name",
    values="Time to max (s)",
).reset_index()

constriction_wide.to_csv(
    export_root / "ConstrictionTime_max_to_end_perExperiment_wide.csv",
    sep=";",
    decimal=",",
    index=False,
)
tracklen_wide.to_csv(
    export_root / "TrackLength_total_perExperiment_wide.csv",
    sep=";",
    decimal=",",
    index=False,
)
timemax_wide.to_csv(
    export_root / "TimeToMax_from_start_perExperiment_wide.csv",
    sep=";",
    decimal=",",
    index=False,
)
per_exp_diam = (
    df_tracks_diam
    .groupby(["Experiment", "Strain name"], as_index=False)[
        ["Diameter (µm)"]
    ]
    .mean()
)

diam_wide = per_exp_diam.pivot(
    index="Experiment",
    columns="Strain name",
    values="Diameter (µm)",
).reset_index()

diam_wide.to_csv(
    export_root / "MaxDiameter_perExperiment_wide.csv",
    sep=";",
    decimal=",",
    index=False,
)

print(f"\nCSV export files written to: {export_root}")