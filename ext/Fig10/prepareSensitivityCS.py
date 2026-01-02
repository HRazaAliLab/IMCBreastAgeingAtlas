"""
Sensitivity analysis (CS): Squidpy neighborhood-enrichment stability.
Author: Eric Lee

- Builds spatial neighbor graph from coordinates using squidpy.
- Computes neighborhood enrichment z-score matrices (sq.gr.nhood_enrichment).
- Uses the base window (largest size) to define the label universe + reformat mapping.
- Compares reduced windows to base using cosine similarity (row-wise mean) and
  normalized Frobenius similarity (1 - normalized Frobenius distance).

Inputs (repo-relative):
- data/derived/cells.parquet

Outputs (repo-relative):
- scratch/ext/Fig10/cs_loss.csv
- scratch/ext/Fig10/fn_loss.csv
- scratch/ext/Fig10/cs_master.csv
- scratch/ext/Fig10/cs_roc.csv

Dependencies:
Python >=3.9, numpy, pandas, scipy, anndata, squidpy, pyarrow

Conda environment is also available in the provided yaml file.
"""

from __future__ import annotations

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

from anndata import AnnData
import squidpy as sq


def repo_root_from_this_file() -> Path:
    return Path(__file__).resolve().parents[3]


def find_densest_point(df: pd.DataFrame) -> tuple[float, float]:
    x = df["X"].values
    y = df["Y"].values
    xy = np.vstack([x, y])
    kde = gaussian_kde(xy)

    x_grid = np.linspace(x.min(), x.max(), 100)
    y_grid = np.linspace(y.min(), y.max(), 100)
    x_grid, y_grid = np.meshgrid(x_grid, y_grid)
    grid_coords = np.vstack([x_grid.ravel(), y_grid.ravel()])
    z = kde(grid_coords).reshape(x_grid.shape)

    max_density_idx = np.argmax(z)
    densest_x, densest_y = grid_coords[:, max_density_idx]
    return float(densest_x), float(densest_y)


def find_densest_square(dense_x: float, dense_y: float, df: pd.DataFrame, s: float) -> tuple[float, float, float, float]:
    half_s = s / 2.0
    left = dense_x - half_s
    right = dense_x + half_s
    up = dense_y + half_s
    down = dense_y - half_s

    if df["X"].max() < right:
        right = float(df["X"].max())
        left = right - s
    if df["X"].min() > left:
        left = float(df["X"].min())
        right = left + s
    if df["Y"].max() < up:
        up = float(df["Y"].max())
        down = up - s
    if df["Y"].min() > down:
        down = float(df["Y"].min())
        up = down + s

    return left, right, up, down


def cosine_similarity(p: np.ndarray, q: np.ndarray) -> np.ndarray:
    #cosine similarity per each row
    p_norm = p / np.linalg.norm(p, axis=1, keepdims=True)
    q_norm = q / np.linalg.norm(q, axis=1, keepdims=True)
    return np.einsum("ij,ij->i", p_norm, q_norm)


def normalized_frobenius_norm(p: np.ndarray, q: np.ndarray) -> float:
    #normalize frobenius norm by sum of norms
    f_norm = np.linalg.norm(p - q, "fro")
    norm_const = np.linalg.norm(p, "fro") + np.linalg.norm(q, "fro")
    if norm_const == 0:
        return 0.0
    return float(f_norm / norm_const)


def similarity_ensemble(p: np.ndarray, q: np.ndarray, alpha: float = 0.5) -> tuple[float, float]:
    #returns (cosine_similarity_mean, 1 - normalized_frobenius_norm)
    cs = float(np.nanmean(cosine_similarity(p, q)))
    fn = normalized_frobenius_norm(p, q)
    return cs, float(1.0 - fn)


def get_int_matrix(df_dict: dict[int, pd.DataFrame], sizes: list[int]) -> dict[int, np.ndarray]:
    int_dict: dict[int, np.ndarray] = {}

    # ----- base -----
    print(f"  Processing window size {sizes[0]}")
    df0 = df_dict[sizes[0]]

    adata0 = AnnData(df0.iloc[:, 2:-3].values, var=pd.DataFrame([], index=df0.columns[2:-3]))
    adata0.obsm["spatial"] = df0[["X", "Y"]].values

    cluster_key = "key"
    adata0.obs[cluster_key] = df0["label"].values.astype(int)

    unique_labels = pd.unique(adata0.obs[cluster_key])
    reformat_dict = dict(zip(sorted(unique_labels), range(len(unique_labels))))

    adata0.obs[cluster_key] = [reformat_dict[i] for i in adata0.obs[cluster_key]]
    adata0.obs[cluster_key] = adata0.obs[cluster_key].astype("category")

    sq.gr.spatial_neighbors(adata0)
    sq.gr.nhood_enrichment(adata0, cluster_key=cluster_key)

    z_df0 = pd.DataFrame(
        adata0.uns[f"{cluster_key}_nhood_enrichment"]["zscore"],
        columns=list(adata0.obs[cluster_key].value_counts(sort=False).index),
        index=list(adata0.obs[cluster_key].value_counts(sort=False).index),
    )

    C = len(unique_labels)
    mtx0 = np.zeros((C, C), dtype=float)
    for i in z_df0.index:
        for j in z_df0.columns:
            mtx0[int(i)][int(j)] = float(z_df0.loc[i, j])

    int_dict[sizes[0]] = mtx0

    # ----- rest -----
    for s in sizes[1:]:
        print(f"  Processing window size {s}")
        df = df_dict[s]

        adata = AnnData(df.iloc[:, 2:-3].values, var=pd.DataFrame([], index=df.columns[2:-3]))
        adata.obsm["spatial"] = df[["X", "Y"]].values
        adata.obs[cluster_key] = df["label"].values.astype(int)

        adata.obs[cluster_key] = [reformat_dict[i] for i in adata.obs[cluster_key]]
        adata.obs[cluster_key] = adata.obs[cluster_key].astype("category")

        sq.gr.spatial_neighbors(adata)
        sq.gr.nhood_enrichment(adata, cluster_key=cluster_key)

        z_df = pd.DataFrame(
            adata.uns[f"{cluster_key}_nhood_enrichment"]["zscore"],
            columns=list(adata.obs[cluster_key].value_counts(sort=False).index),
            index=list(adata.obs[cluster_key].value_counts(sort=False).index),
        )

        mtx = np.zeros((C, C), dtype=float)
        for i in z_df.index:
            for j in z_df.columns:
                mtx[int(i)][int(j)] = float(z_df.loc[i, j])

        int_dict[s] = mtx

    return int_dict


def main():
    repo_root = repo_root_from_this_file()
    out_dir = repo_root / "scratch" / "ext" / "Fig10"
    out_dir.mkdir(parents=True, exist_ok=True)

    sizes = [1200, 1100, 1000, 900, 800, 700, 600, 500, 400]
    img_wanted = 200

    cells_fp = repo_root / "data" / "derived" / "cells.parquet"
    cols_needed = ["ImageID", "CellID", "ClusterID", "CenterX", "CenterY"]
    slide_df = pd.read_parquet(cells_fp, columns=cols_needed).copy()
    slide_df = slide_df.rename(columns={"CenterX": "X", "CenterY": "Y"})

    slide_df["label"] = slide_df["ClusterID"].astype(int)

    import random
    random.seed(420)
    unique_rois = pd.unique(slide_df["ImageID"])
    total_imgs = len(unique_rois)

    if img_wanted > total_imgs:
        img_wanted = total_imgs

    ran_smp = random.sample(range(0, total_imgs), img_wanted)

    cs_loss_aggregate = []
    fn_loss_aggregate = []

    for i in ran_smp:
        df = slide_df[slide_df["ImageID"] == unique_rois[i]].reset_index(drop=True)
        print(f"Processing ROI {df.ImageID.iloc[0]}")

        if df["X"].max() < sizes[0] and df["Y"].max() < sizes[0]:
            print(f"  {df.ImageID.iloc[0]} is smaller than the required max dimension; skipping")
            continue

        dense_x, dense_y = find_densest_point(df)
        df_dict: dict[int, pd.DataFrame] = {}
        for s in sizes:
            left, right, up, down = find_densest_square(dense_x, dense_y, df, s)
            df_dict[s] = df[(df["X"] < right) & (df["X"] > left) & (df["Y"] < up) & (df["Y"] > down)]

        int_dict = get_int_matrix(df_dict, sizes)

        cs_losses = [similarity_ensemble(int_dict[s], int_dict[sizes[0]])[0] for s in sizes]
        fn_losses = [similarity_ensemble(int_dict[s], int_dict[sizes[0]])[1] for s in sizes]

        cs_loss_aggregate.append(np.array(cs_losses))
        fn_loss_aggregate.append(np.array(fn_losses))

    cs_loss = pd.DataFrame(np.array(cs_loss_aggregate), columns=sizes)
    fn_loss = pd.DataFrame(np.array(fn_loss_aggregate), columns=sizes)

    cs_loss.to_csv(out_dir / "cs_loss.csv", index=False)
    fn_loss.to_csv(out_dir / "fn_loss.csv", index=False)

    cs_master = cs_loss.melt(var_name="variable", value_name="value")
    cs_master["variable"] = cs_master["variable"].astype(int)
    cs_master.to_csv(out_dir / "cs_master.csv", index=False)

    cs_med = cs_loss.median(axis=0).reindex(sizes)
    cs_roc = pd.DataFrame({"V1": sizes, "cosine_sim": cs_med.diff().fillna(0.0).to_numpy()})
    cs_roc.to_csv(out_dir / "cs_roc.csv", index=False)

    print(f"\nWrote CS outputs to: {out_dir.resolve()}")


if __name__ == "__main__":
    main()
