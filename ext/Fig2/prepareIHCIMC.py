"""
prepareIHCIMC.py
AUTHOR: Eric Lee

Dependencies
------------
This script requires the following Python packages:
    - python >= 3.9
    - numpy
    - pandas
    - scipy
    - pyarrow   (for reading Parquet IMC files)

The conda environment ihc-imc.yaml is also provided in this folder.

Input Data Requirements
-----------------------
This script assumes that the following dataset has been downloaded and placed
in the repository at:

    data/derived/NormalBreast13_TMA_raw_data.csv

This file corresponds to the published NormalBreast13 tissue microarray (TMA)
dataset originally described in:

    Osako et al.
    "Age-correlated protein and transcript expression in breast cancer and
    normal breast tissues is dominated by host endocrine effects."
    Nature Cancer (2020).

The original data are also publicly available at:
    https://github.com/BCCRCMO/BrCa_AgeAssociations/blob/master/BrCa_Age_Associated_TMA/NormalBreast13/Analysis/BuildData/NormalBreast13_TMA_raw_data.csv

Patient Identifiers
-------------------
The PatientID field in this dataset (e.g., "166-03", "87-08") matches the
PatientIDs used in the NormalBreast13 dataset and is preserved during processing.

Overview
--------
This script:
  1) Cleans the published IHC TMA data using logic matching the original analysis.
  2) Aggregates IMC cell-level data to patient-level means.
  3) Joins IHC and IMC measurements using block and within-block identifiers.
  4) Writes marker-specific CSV files for downstream figure generation in R.

Output
------
Processed files are written to:
    scratch/ext/Fig2/

including marker-specific CSVs (e.g., ck5.csv, er.csv, ki67.csv) and an
intermediate joined table for reproducibility. The data are then plotted with plotIHCValidation.R.
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


def build_section_dict():
    sections = [f"R{i}C{j}" for i in range(2, 9) for j in range(1, 13)]
    WITHIN_ID = [
        1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6,
        12, 12, 11, 11, 10, 10, 9, 9, 8, 8, 7, 7,
        13, 13, 14, 14, 15, 15, 16, 16, 17, 17, 18, 18,
        24, 24, 23, 23, 22, 22, 21, 21, 20, 20, 19, 19,
        25, 25, 26, 26, 27, 27, 28, 28, 29, 29, 30, 30,
        36, 36, 35, 35, 34, 34, 33, 33, 32, 32, 31, 31,
        37, 37, 38, 38, 39, 39, 40, 40, 41, 41, 42, 42
    ]
    assert len(sections) == len(WITHIN_ID), "Section mapping length mismatch."
    return dict(zip(sections, WITHIN_ID))


def clean_ihc_raw_eric_exact(raw_fp: Path, out_dir: Path) -> pd.DataFrame:
    ihc_raw = pd.read_csv(raw_fp)

    keep = [
        "nb13_id", "TMA.Block", "within.TMA.Block.ID", "Patient.ID", "Age", "Laterality",
        "PctPosCells.ck5", "Intensity.ck5", "PctPosCells.esr1", "Intensity.esr1",
        "PctPosCells.pgr", "Intensity.pgr", "PctPosCells.her2", "Intensity.her2",
        "PctPosCells.foxa1", "Intensity.foxa1", "PctPosCells.mki67", "Intensity.mki67",
        "PctPosCells.ecad", "Intensity.ecad", "PctPosCells.egfr", "Intensity.egfr",
        "PctPosCells.bcl2", "Intensity.bcl2", "PctPosCells.tp53", "Intensity.tp53",
        "PctPosCells.ezh2", "Intensity.ezh2", "PctPosCells.h3k27me3", "Intensity.h3k27me3"
    ]
    ihc_cleaned = ihc_raw[keep].copy()
    ihc_cleaned.columns = [
        "nb13_id", "block_id", "within_block_id", "patient_id", "age", "laterality",
        "pp_CK5", "CK5", "pp_ER", "ER", "pp_PR", "PR", "pp_HER2", "HER2",
        "pp_FOXA1", "FOXA1", "pp_Ki67", "Ki67",
        "pp_ECAD", "ECAD", "pp_EGFR", "EGFR",
        "pp_BCL2", "BCL2", "pp_p53", "p53",
        "pp_EZH2", "EZH2", "pp_H3K27me3", "H3K27me3"
    ]

    paired_fp = out_dir / "nb13_tma_paired.csv"
    ihc_cleaned.to_csv(paired_fp, index=False)

    ihc = pd.read_csv(paired_fp)

    for c in ihc.columns[6:]:
        for r in range(len(ihc[c])):
            if str(ihc[c][r]).isalpha():
                ihc.at[r, c] = np.nan

    ihc["PR"] = ["3" if i == "0" else np.nan for i in ihc["PR"]]
    ihc["FOXA1"] = ["3" if i == "0" else np.nan for i in ihc["FOXA1"]]
    ihc["Ki67"] = ["3" if i == "0" else np.nan for i in ihc["Ki67"]]
    ihc["pp_EGFR"] = np.nan

    for c in ihc.columns[6:]:
        for r in range(len(ihc[c])):
            if str(ihc[c][r]) == "<1":
                ihc.at[r, c] = "0.5"
            elif str(ihc[c][r]) == "0c":
                ihc.at[r, c] = "0"
            elif not str(ihc[c][r]).isdigit():
                ihc.at[r, c] = np.nan

    cleaned_fp = out_dir / "nb13_tma_paired_cleaned.csv"
    ihc.to_csv(cleaned_fp, index=False)

    print("\n=== IHC cleaning diagnostics ===")
    print(f"Raw rows: {len(ihc_raw):,}")
    print(f"Paired rows: {len(ihc):,}")
    for k in ["pp_CK5", "pp_ER", "pp_PR", "pp_HER2", "pp_FOXA1", "pp_Ki67"]:
        n = pd.to_numeric(ihc[k], errors="coerce").notna().sum()
        print(f"{k:10s} non-NA: {n:,}")

    return ihc


def main():
    # repo root (script is at code/ext/Fig2/prepareIHCIMC.py)
    repo_root = Path(__file__).resolve().parents[3]

    ihc_raw_fp = repo_root / "data" / "derived" / "NormalBreast13_TMA_raw_data.csv"
    
    out_dir = repo_root / "scratch" / "ext" / "Fig2"
    out_dir.mkdir(parents=True, exist_ok=True)

    imc_fp = repo_root / "data" / "derived" / "cells.parquet"

    section_dict = build_section_dict()

    ihc_clean = clean_ihc_raw_eric_exact(ihc_raw_fp, out_dir)

    ihc_patient = ihc_clean.copy()
    ihc_patient["ID"] = ihc_patient["block_id"].astype(str) + "_" + ihc_patient["within_block_id"].astype(str)
    ihc_patient = ihc_patient.drop(columns=["block_id", "within_block_id"])
    ihc_patient = ihc_patient[["ID", "pp_CK5", "pp_ER", "pp_PR", "pp_HER2", "pp_FOXA1", "pp_Ki67"]].copy()

    imc = pd.read_parquet(imc_fp, columns=["ImageID", "CK5/14", "ER", "PR", "HER2", "FOXA1", "Ki67"]).copy()
    marker_cols = ["CK5/14", "ER", "PR", "HER2", "FOXA1", "Ki67"]

    for col in range(len(marker_cols)):
        imc.iloc[:, imc.columns.get_loc(marker_cols[col])] = np.arcsinh(imc[marker_cols[col]].astype(float))

    img = imc["ImageID"].astype(str)
    tokens = img.str.split("_")

    prefix = tokens.str[0]
    section_token_raw = tokens.str[1]  # "R5C11" or "R5C11_split"
    section_token = section_token_raw.str.replace(r"_split.*$", "", regex=True)

    imc["block_id"] = prefix.str.extract(r"^NormalBreastImg([A-Za-z])$", expand=False)
    if imc["block_id"].isna().any():
        bad = imc.loc[imc["block_id"].isna(), "ImageID"].head(10).tolist()
        raise ValueError(f"Couldn't parse block letter from ImageID prefix (up to 10): {bad}")

    imc["within_block_id"] = section_token.map(section_dict)
    if imc["within_block_id"].isna().any():
        bad_tokens = section_token_raw[imc["within_block_id"].isna()].value_counts().head(10)
        bad_ids = imc.loc[imc["within_block_id"].isna(), "ImageID"].head(10).tolist()
        raise ValueError(
            "Unmapped section tokens (after Eric-emulation).\n"
            f"Top bad raw tokens:\n{bad_tokens}\n"
            f"Example ImageIDs (up to 10): {bad_ids}"
        )

    imc_patient = []
    for b in pd.unique(imc["block_id"]):
        block_subset = imc[imc["block_id"] == b]
        for w in pd.unique(block_subset["within_block_id"]):
            w_subset = block_subset[block_subset["within_block_id"] == w]
            w_exp = np.nanmean(w_subset[marker_cols].to_numpy(), axis=0)
            imc_patient.append([f"{b}_{int(w)}"] + list(w_exp))

    imc_patient = pd.DataFrame(imc_patient, columns=["ID"] + marker_cols)

    joint = pd.merge(imc_patient, ihc_patient, on=["ID"], how="inner")

    print("\n=== Join diagnostics ===")
    print(f"IHC patients (with IDs): {len(ihc_patient):,}")
    print(f"IMC patients (with IDs): {len(imc_patient):,}")
    print(f"Inner-joined:            {len(joint):,}")

    # write 6 small CSVs for R (IHCpercentage, IMCintensity) to scratch/ext/Fig2
    pairs = {
        "ck5.csv":   ("pp_CK5",  "CK5/14"),
        "er.csv":    ("pp_ER",   "ER"),
        "pr.csv":    ("pp_PR",   "PR"),
        "her2.csv":  ("pp_HER2", "HER2"),
        "foxa1.csv": ("pp_FOXA1","FOXA1"),
        "ki67.csv":  ("pp_Ki67", "Ki67"),
    }

    print("\n=== Writing scatter inputs ===")
    for fn, (xcol, ycol) in pairs.items():
        df = joint[[xcol, ycol]].dropna().copy()
        df.columns = ["IHCpercentage", "IMCintensity"]
        df.to_csv(out_dir / fn, index=False)

        rho, p = spearmanr(df["IHCpercentage"], df["IMCintensity"])
        print(f"{fn:10s}  n={len(df):4d}  rho={rho: .3f}  p={p: .3e}")

    joint.to_csv(out_dir / "joint_ihc_imc_patient_means.csv", index=False)
    print(f"\nWrote outputs to: {out_dir.resolve()}")


if __name__ == "__main__":
    main()
