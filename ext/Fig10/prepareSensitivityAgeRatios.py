from __future__ import annotations

"""
Sensitivity (age effect) on public schema:
1) Compute duct/(duct+lobule) in densest square per ROI and size -> ratios.csv
2) For each (NImages, Size): fit spline vs linear and do LRT -> panel_stats.csv + panel_series.csv

Author: Eric Lee

Inputs:
- data/derived/cells.parquet        (ImageID, CellID, CenterX, CenterY)
- data/derived/cellContext.parquet  (ImageID, CellID, TissueStructure)
- data/derived/clinical.parquet     (ImageID, Age)

Outputs:
- scratch/ext/Fig10/ratios.csv
- scratch/ext/Fig10/panel_stats.csv
- scratch/ext/Fig10/panel_series.csv

Dependencies:
- Python >= 3.9
- numpy
- pandas
- scipy  (scipy.stats.gaussian_kde, scipy.optimize.curve_fit, scipy.interpolate.LSQUnivariateSpline, scipy.stats.chi2)
- pyarrow  (for reading parquet via pandas)

This script can be run with the provided sensitivity.yaml conda environment. The output can be plotted with plotSensitivityAgeRatios.R.
"""

from pathlib import Path
import random
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde, chi2
from scipy.optimize import curve_fit
from scipy.interpolate import LSQUnivariateSpline


STRUCT_COL = "TissueStructure"
XCOL, YCOL = "CenterX", "CenterY"

SIZES = [400, 600, 800, 1000, 1200]
IMAGE_NUMS_FOR_FIT = [50, 100, 400, 700, 1400]

RANDOM_SEED = 420
KDE_GRID_N = 100
NOT_BIG_ENOUGH_ROI_IDXS = [75, 93, 134, 515, 620, 1170, 1273]


def repo_root_from_this_file() -> Path:
    return Path(__file__).resolve().parents[3]


def find_densest_point(x: np.ndarray, y: np.ndarray) -> tuple[float, float] | None:
    # gaussian_kde needs n_samples >= n_dims (here 2), and non-singular covariance.
    if x.size < 2:
        return None
    if np.allclose(x, x[0]) and np.allclose(y, y[0]):
        return None  # all points identical -> singular covariance

    xy = np.vstack([x, y])
    try:
        kde = gaussian_kde(xy)
    except Exception:
        return None

    x_grid = np.linspace(x.min(), x.max(), KDE_GRID_N)
    y_grid = np.linspace(y.min(), y.max(), KDE_GRID_N)
    xg, yg = np.meshgrid(x_grid, y_grid)

    grid_coords = np.vstack([xg.ravel(), yg.ravel()])
    z = kde(grid_coords).reshape(xg.shape)

    densest_x, densest_y = grid_coords[:, int(np.argmax(z))]
    return float(densest_x), float(densest_y)

def densest_square_bounds(dense_x: float, dense_y: float, x: np.ndarray, y: np.ndarray, s: float):
    half = s / 2.0
    left, right = dense_x - half, dense_x + half
    down, up = dense_y - half, dense_y + half

    xmin, xmax = float(x.min()), float(x.max())
    ymin, ymax = float(y.min()), float(y.max())

    if right > xmax:
        right = xmax
        left = xmax - s
    if left < xmin:
        left = xmin
        right = xmin + s

    if up > ymax:
        up = ymax
        down = ymax - s
    if down < ymin:
        down = ymin
        up = ymin + s

    return left, right, up, down


def linear_func(x, a, b):
    return a * x + b


def loglik_from_rss(rss: float, n: int) -> float:
    sigma2 = rss / n
    return -n / 2 * np.log(2 * np.pi * sigma2) - n / 2


def main():
    repo_root = repo_root_from_this_file()
    out_dir = repo_root / "scratch" / "ext" / "Fig10"
    out_dir.mkdir(parents=True, exist_ok=True)

    cells_fp = repo_root / "data" / "derived" / "cells.parquet"
    ctx_fp = repo_root / "data" / "derived" / "cellContext.parquet"
    clinical_fp = repo_root / "data" / "derived" / "clinical.parquet"

    cells = pd.read_parquet(cells_fp, columns=["ImageID", "CellID", XCOL, YCOL])
    ctx = pd.read_parquet(ctx_fp, columns=["ImageID", "CellID", STRUCT_COL])
    clinical = pd.read_parquet(clinical_fp, columns=["ImageID", "Age"])

    df = (
        cells.merge(ctx, on=["ImageID", "CellID"], how="inner")
             .merge(clinical, on=["ImageID"], how="inner")
    )

    df[STRUCT_COL] = df[STRUCT_COL].astype(str).str.strip().str.lower()
    df = df[df[STRUCT_COL].isin(["duct", "lobule"])].copy()

    unique_rois = list(pd.unique(df["ImageID"]))
    all_idxs = list(set(range(len(unique_rois))) - set(NOT_BIG_ENOUGH_ROI_IDXS))

    random.seed(RANDOM_SEED)
    random.shuffle(all_idxs)

    roi_order = [unique_rois[i] for i in all_idxs]

    rows = []
    for roi in roi_order:
        roi_df = df[df["ImageID"] == roi]
        x = roi_df[XCOL].to_numpy(dtype=float)
        y = roi_df[YCOL].to_numpy(dtype=float)

        if len(roi_df) == 0:
            continue

        dense = find_densest_point(x, y)
        if dense is None:
            continue
        dense_x, dense_y = dense


        for s in SIZES:
            left, right, up, down = densest_square_bounds(dense_x, dense_y, x, y, float(s))
            sub = roi_df[
                (roi_df[XCOL] > left) & (roi_df[XCOL] < right) &
                (roi_df[YCOL] > down) & (roi_df[YCOL] < up)
            ]

            if len(sub) == 0:
                continue

            counts = sub[STRUCT_COL].value_counts()
            duct = int(counts.get("duct", 0))
            lob = int(counts.get("lobule", 0))
            denom = duct + lob
            ratio = duct / denom if denom > 0 else np.nan
            age = float(sub["Age"].iloc[0])

            rows.append([roi, int(s), age, float(ratio)])

    ratios = pd.DataFrame(rows, columns=["ImageID", "Size", "Age", "Ratio"]).dropna()
    ratios.to_csv(out_dir / "ratios.csv", index=False)

    ratios["ImageID"] = pd.Categorical(ratios["ImageID"], categories=roi_order, ordered=True)

    stats_rows = []
    series_rows = []

    for n_imgs in IMAGE_NUMS_FOR_FIT:
        rois_n = set(roi_order[: min(n_imgs, len(roi_order))])
        sub_n = ratios[ratios["ImageID"].isin(rois_n)].copy()

        for s in SIZES:
            sub = sub_n[sub_n["Size"] == s].copy()
            if sub.empty:
                stats_rows.append([s, n_imgs, 0, np.nan])
                continue

            sub = sub.sort_values("Age")
            grp = sub[["Age", "Ratio"]].groupby("Age").mean()
            x = grp.index.to_numpy(dtype=float)
            y = grp["Ratio"].to_numpy(dtype=float)

            spline = LSQUnivariateSpline(x, y, [50.0], k=2)
            y_spline = spline(x)

            (a, b), _ = curve_fit(linear_func, x, y)
            y_linear = linear_func(x, a, b)

            rss_lin = float(np.sum((y - y_linear) ** 2))
            rss_spl = float(np.sum((y - y_spline) ** 2))

            n = len(y)
            test_stat = 2 * (loglik_from_rss(rss_spl, n) - loglik_from_rss(rss_lin, n))
            pval = float(1 - chi2.cdf(test_stat, 2))

            stats_rows.append([int(s), int(n_imgs), int(len(x)), pval])

            for age_i, sy, ly in zip(x, y_spline, y_linear):
                series_rows.append([int(s), int(n_imgs), float(age_i), float(sy), float(ly)])

    pd.DataFrame(stats_rows, columns=["Size", "NImages", "n_points", "p_value"]).to_csv(
        out_dir / "panel_stats.csv", index=False
    )
    pd.DataFrame(series_rows, columns=["Size", "NImages", "Age", "spline_y", "linear_y"]).to_csv(
        out_dir / "panel_series.csv", index=False
    )


if __name__ == "__main__":
    main()
