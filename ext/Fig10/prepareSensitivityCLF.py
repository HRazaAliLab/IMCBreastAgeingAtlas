from __future__ import annotations
"""
Sensitivity analysis (classifier): tissue compartment prediction robustness.
Author: Eric Lee

This script trains a multiclass classifier to predict tissue compartment
(e.g. duct vs lobule) from IMC-derived cell phenotypes and evaluates how
classification accuracy degrades as the number of images used for training
is progressively reduced.

Inputs (repo-relative):
- data/derived/cells.parquet
- data/derived/cellClusterAnnotation.csv
- data/derived/cellAdjacency.parquet
- data/derived/cellContext.parquet

Outputs (repo-relative):
- scratch/ext/Fig10/clf_loss_label.csv
- scratch/ext/Fig10/clf_loss_structure.csv

Dependencies:
Python >=3.9, numpy, pandas, lightgbm, scikit-learn, pyarrow
Conda environment is provided as sensitivity.yaml.

Data are plotted with plotSensitivityTests.R.
"""

from pathlib import Path
import random
import numpy as np
import pandas as pd
import lightgbm as lgb
from sklearn.metrics import accuracy_score


def repo_root_from_this_file() -> Path:
    return Path(__file__).resolve().parents[3]


MARKER_COLS = [
    "CD163", "CD20", "Cyclin D1", "CD56", "CD45", "CD8",
    "GATA3", "CD11c", "CD3", "ER", "SMA", "PR", "HER2",
    "PD-1", "IDO", "AR", "PD-L1", "GZMB", "Ki67", "CD4",
    "CK5/14", "TCF1", "PDGFRbeta", "CD31", "CK7", "PDPN",
    "HLA-ABC", "FOXA1", "panCK", "pH2AX", "CK8/18",
    "Vimentin", "Calponin", "Caveolin-1", "CD15", "MPO",
    "HLA-DR", "CD68", "CD79a", "CA9", "LDHA",
]


def run_one_target(
    slide_df: pd.DataFrame,
    y_col: str,
    image_nums: list[int],
    total_imgs: list
) -> list[float]:

    print(f"  Preparing data for target: {y_col}")

    df = slide_df.reset_index(drop=True)

    X_df = df[MARKER_COLS].apply(pd.to_numeric, errors="coerce")
    keep = ~X_df.isna().any(axis=1)
    df = df.loc[keep].reset_index(drop=True)

    print(f"    Retained {len(df):,} cells after NA filtering")

    y_raw = df[y_col].to_numpy()
    classes = pd.unique(y_raw)
    class_map = {c: i for i, c in enumerate(classes)}

    y_all = np.array([class_map[v] for v in y_raw], dtype=int)
    num_class = int(y_all.max()) + 1

    print(f"    Number of classes: {num_class}")

    params = {
        "objective": "multiclass",
        "num_class": num_class,
        "metric": "multi_logloss",
        "verbose": -1,
        "seed": 420,
        "num_threads": 10,
    }

    X_all = df[MARKER_COLS].to_numpy(dtype=float)

    accs = []
    for i, n in enumerate(image_nums, start=1):
        print(f"    [{i}/{len(image_nums)}] Training with first {n} images")

        keep_imgs = set(total_imgs[: min(int(n), len(total_imgs))])
        sub = df[df["ImageID"].isin(keep_imgs)].reset_index(drop=True)

        if len(sub) == 0:
            print("      No cells after filtering — skipping")
            accs.append(np.nan)
            continue

        X_train = sub[MARKER_COLS].to_numpy(dtype=float)
        y_train = np.array([class_map[v] for v in sub[y_col].to_numpy()], dtype=int)

        if len(np.unique(y_train)) < 2:
            print("      <2 classes in training subset — skipping")
            accs.append(np.nan)
            continue

        print(f"      Training on {len(sub):,} cells")

        train_data = lgb.Dataset(X_train, label=y_train)
        model = lgb.train(params, train_data, num_boost_round=100)

        print("      Predicting on full dataset")
        pred = model.predict(X_all, num_iteration=model.best_iteration)
        pred_lab = np.argmax(pred, axis=1)

        acc = float(accuracy_score(y_all, pred_lab))
        accs.append(acc)

        print(f"      Accuracy = {acc:.4f}")

    return accs


def main():
    repo_root = repo_root_from_this_file()
    out_dir = repo_root / "scratch" / "ext" / "Fig10"
    out_dir.mkdir(parents=True, exist_ok=True)

    print("Loading input data")

    cells_fp = repo_root / "data" / "derived" / "cells.parquet"
    ctx_fp   = repo_root / "data" / "derived" / "cellContext.parquet"

    cols_cells = ["ImageID", "CellID", "CenterX", "CenterY",
                  "ClusterID", "SpatialClusterID"] + MARKER_COLS
    cells = pd.read_parquet(cells_fp, columns=cols_cells)
    ctx = pd.read_parquet(ctx_fp, columns=["ImageID", "CellID", "TissueStructure"])

    print(f"  Loaded {len(cells):,} cells")

    slide_df = pd.merge(cells, ctx, on=["ImageID", "CellID"], how="inner")
    slide_df = slide_df.rename(columns={"CenterX": "X", "CenterY": "Y"})

    needed = ["ImageID", "CellID", "X", "Y",
              "ClusterID", "SpatialClusterID", "TissueStructure"] + MARKER_COLS
    slide_df = slide_df.dropna(subset=needed).reset_index(drop=True)

    print(f"  Retained {len(slide_df):,} cells after merge and filtering")

    random.seed(420)
    np.random.seed(420)

    unique_rois = pd.unique(slide_df["ImageID"])
    print(f"  Found {len(unique_rois)} unique images")

    all_imgs = list(range(len(unique_rois)))
    random.shuffle(all_imgs)
    total_imgs = [unique_rois[i] for i in all_imgs]

    image_nums = [
        1400, 1300, 1200, 1100, 1000, 900, 800, 700, 600, 500,
        450, 400, 350, 300, 250, 200, 150, 100, 50, 40, 30, 20, 10, 5
    ]

    TARGETS = [
        ("label", "ClusterID"),
        ("spatial_label", "SpatialClusterID"),
        ("structure", "TissueStructure"),
    ]

    master_rows = []
    roc_cols = {"V1": image_nums}

    for which_name, y_col in TARGETS:
        print(f"\n=== Running target: {which_name} ({y_col}) ===")

        accs = run_one_target(
            slide_df,
            y_col=y_col,
            image_nums=image_nums,
            total_imgs=total_imgs
        )

        master_rows.extend(
            [{"which": which_name, "variable": n, "value": v}
             for n, v in zip(image_nums, accs)]
        )

        roc_cols[which_name] = pd.Series(accs).diff().fillna(0.0).to_list()

        out_fp = out_dir / f"clf_loss_{y_col}.csv"
        pd.DataFrame(np.array(accs), index=image_nums).to_csv(out_fp)
        print(f"  Wrote {out_fp}")

    pd.DataFrame(master_rows).to_csv(out_dir / "clf_master.csv", index=False)
    pd.DataFrame(roc_cols).to_csv(out_dir / "clf_roc.csv", index=False)

    print("\nAll classifier sensitivity analyses complete.")


if __name__ == "__main__":
    main()
