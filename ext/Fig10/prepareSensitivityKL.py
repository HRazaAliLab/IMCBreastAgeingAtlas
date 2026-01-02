"""
Sensitivity analysis (KL divergence): cell type composition stability.
Author: Eric Lee

This script quantifies how cell type distributions change as IMC images are
progressively restricted to smaller spatial windows. For each image, it
computes KL divergence between the full-size distribution and reduced-size
distributions for epithelial, immune, stromal, and coarse cell-type groupings.

Inputs (repo-relative):
- data/derived/cells.parquet
- data/derived/cellClusterAnnotation.csv

Outputs (repo-relative):
- scratch/ext/Fig2/kl_epi_loss.csv
- scratch/ext/Fig2/kl_imm_loss.csv
- scratch/ext/Fig2/kl_str_loss.csv
- scratch/ext/Fig2/kl_type_loss.csv

Dependencies:
Python >=3.9, numpy, pandas, scipy, pyarrow
(compatible with the `sensitivity` conda environment provided as a yaml)

Data are plotted with plotSensitivityTests.R
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import gaussian_kde

def kl_divergence(p: np.ndarray, q: np.ndarray) -> float:
    # Guard q==0 to avoid inf explosions in sparse categories
    q_safe = np.where(q == 0, 1e-12, q)
    return float(np.sum(np.where(p != 0, p * np.log(p / q_safe), 0.0)))


def get_simplex_distr(df: pd.DataFrame, k: int) -> tuple[np.ndarray, np.ndarray, np.ndarray]:

    count = df["label"].value_counts()
    count_np = np.zeros(k, dtype=float)
    for i in count.index:
        # labels are 1..k
        count_np[int(i) - 1] = float(count[i])

    epi = count_np[:11]
    imm = count_np[11:19]
    strom = count_np[19:]

    epi_distr = epi / epi.sum() if epi.sum() > 0 else np.zeros_like(epi)
    imm_distr = imm / imm.sum() if imm.sum() > 0 else np.zeros_like(imm)
    str_distr = strom / strom.sum() if strom.sum() > 0 else np.zeros_like(strom)
    return epi_distr, imm_distr, str_distr

def get_simplex_type_distr(df: pd.DataFrame, k: int = 3) -> np.ndarray:
    count = df["label_type"].value_counts()
    count_np = np.zeros(k, dtype=float)
    for i in count.index:
        count_np[int(i)] = float(count[i])
    return count_np / count_np.sum() if count_np.sum() > 0 else np.zeros_like(count_np)


def find_densest_point(df: pd.DataFrame) -> tuple[float, float]:
    x = df["X"].to_numpy()
    y = df["Y"].to_numpy()

    xy = np.vstack([x, y])
    kde = gaussian_kde(xy)

    x_grid = np.linspace(x.min(), x.max(), 100)
    y_grid = np.linspace(y.min(), y.max(), 100)
    x_grid, y_grid = np.meshgrid(x_grid, y_grid)
    grid_coords = np.vstack([x_grid.ravel(), y_grid.ravel()])
    z = kde(grid_coords).reshape(x_grid.shape)

    max_density_idx = int(np.argmax(z))
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


def repo_root_from_this_file() -> Path:
    return Path(__file__).resolve().parents[3]


def load_cluster_annotation(fp: Path) -> pd.DataFrame:
    ann = pd.read_csv(fp)
    if "ClusterID" not in ann.columns:
        raise ValueError("cellClusterAnnotation.csv must contain a 'ClusterID' column.")

    ann["ClusterID"] = ann["ClusterID"].astype(int)

    if "isEpithelial" in ann.columns:
        is_epi = ann["isEpithelial"].astype(bool)
    else:
        is_epi = pd.Series(False, index=ann.index)

    typ = ann["Type"].astype(str).str.lower() if "Type" in ann.columns else pd.Series("", index=ann.index)
    is_imm = typ.str.contains("imm")

    ann["label_type"] = np.where(is_epi, 0, np.where(is_imm, 1, 2)).astype(int)

    if "PrintOrder" in ann.columns:
        ann = ann.sort_values(["PrintOrder", "ClusterID"], kind="mergesort")
    else:
        ann = ann.sort_values(["ClusterID"], kind="mergesort")

    #create a 1..k label index for the ordering
    ann["label_index"] = np.arange(1, len(ann) + 1, dtype=int)
    return ann[["ClusterID", "label_type", "label_index"]].copy()


def main():
    repo_root = repo_root_from_this_file()

    cells_fp = repo_root / "data" / "derived" / "cells.parquet"
    ann_fp = repo_root / "data" / "derived" / "cellClusterAnnotation.csv"

    out_dir = repo_root / "scratch" / "ext" / "Fig10"
    out_dir.mkdir(parents=True, exist_ok=True)

    ann = load_cluster_annotation(ann_fp)
    k = int(ann["label_index"].max())

    cols_needed = ["ImageID", "ClusterID", "CenterX", "CenterY"]
    cells = pd.read_parquet(cells_fp, columns=cols_needed).copy()
    cells = cells.rename(columns={"CenterX": "X", "CenterY": "Y"})
    cells["ClusterID"] = cells["ClusterID"].astype(int)

    cells = cells.merge(ann, on="ClusterID", how="inner")
    cells = cells.rename(columns={"label_index": "label"})

    sizes = [1200, 1100, 1000, 900, 800, 700, 600, 500, 400]
    base_size = sizes[0]

    unique_rois = pd.unique(cells["ImageID"])
    total_imgs = len(unique_rois)

    kl_epi_loss_aggregate = []
    kl_imm_loss_aggregate = []
    kl_str_loss_aggregate = []
    kl_type_loss_aggregate = []
    kept_image_ids = []

    for i in range(total_imgs):
        df = cells[cells["ImageID"] == unique_rois[i]].reset_index(drop=True)
        if df.empty:
            continue

        print(f"Processing {df['ImageID'].iloc[0]}")

        if df["X"].max() < base_size and df["Y"].max() < base_size:
            print(f"{df['ImageID'].iloc[0]} is smaller than the required max dimension!")
            continue

        dense_x, dense_y = find_densest_point(df)

        df_dict = {}
        for s in sizes:
            left, right, up, down = find_densest_square(dense_x, dense_y, df, s)
            df_dict[s] = df[(df["X"] < right) & (df["X"] > left) & (df["Y"] < up) & (df["Y"] > down)]

        distr_dict = {}
        distr_type_dict = {}
        for s in sizes:
            distr_dict[s] = get_simplex_distr(df_dict[s], k=k)
            distr_type_dict[s] = get_simplex_type_distr(df_dict[s], k=3)

        epi_kl_losses = [kl_divergence(distr_dict[s][0], distr_dict[base_size][0]) for s in sizes]
        imm_kl_losses = [kl_divergence(distr_dict[s][1], distr_dict[base_size][1]) for s in sizes]
        str_kl_losses = [kl_divergence(distr_dict[s][2], distr_dict[base_size][2]) for s in sizes]
        type_kl_losses = [kl_divergence(distr_type_dict[s], distr_type_dict[base_size]) for s in sizes]

        kl_epi_loss_aggregate.append(np.array(epi_kl_losses, dtype=float))
        kl_imm_loss_aggregate.append(np.array(imm_kl_losses, dtype=float))
        kl_str_loss_aggregate.append(np.array(str_kl_losses, dtype=float))
        kl_type_loss_aggregate.append(np.array(type_kl_losses, dtype=float))
        kept_image_ids.append(df["ImageID"].iloc[0])

    epi_mat = pd.DataFrame(np.array(kl_epi_loss_aggregate), columns=sizes)
    imm_mat = pd.DataFrame(np.array(kl_imm_loss_aggregate), columns=sizes)
    str_mat = pd.DataFrame(np.array(kl_str_loss_aggregate), columns=sizes)
    typ_mat = pd.DataFrame(np.array(kl_type_loss_aggregate), columns=sizes)

    epi_mat.to_csv(out_dir / "kl_epi_loss.csv", index=False)
    imm_mat.to_csv(out_dir / "kl_imm_loss.csv", index=False)
    str_mat.to_csv(out_dir / "kl_str_loss.csv", index=False)
    typ_mat.to_csv(out_dir / "kl_type_loss.csv", index=False)

    pd.Series(kept_image_ids, name="ImageID").to_csv(out_dir / "kl_images_used.csv", index=False)

    def to_master(mat: pd.DataFrame, which: str) -> pd.DataFrame:
        long = mat.melt(var_name="variable", value_name="value")
        long["which"] = which
        long["variable"] = long["variable"].astype(int)
        return long[["which", "variable", "value"]]

    kl_master = pd.concat(
        [
            to_master(epi_mat, "epi"),
            to_master(imm_mat, "imm"),
            to_master(str_mat, "str"),
            to_master(typ_mat, "type"),
        ],
        ignore_index=True,
    )
    kl_master.to_csv(out_dir / "kl_master.csv", index=False)

    def roc_curve(mat: pd.DataFrame) -> pd.Series:
        med = mat.median(axis=0)
        med = med.reindex(sizes)
        d = med.diff().fillna(0.0)
        return d

    kl_roc = pd.DataFrame({"V1": sizes})
    kl_roc["epi"] = roc_curve(epi_mat).to_numpy()
    kl_roc["imm"] = roc_curve(imm_mat).to_numpy()
    kl_roc["str"] = roc_curve(str_mat).to_numpy()
    kl_roc["type"] = roc_curve(typ_mat).to_numpy()
    kl_roc.to_csv(out_dir / "kl_roc.csv", index=False)

    print(f"\nWrote KL outputs to: {out_dir.resolve()}")


if __name__ == "__main__":
    main()
