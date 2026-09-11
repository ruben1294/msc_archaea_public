#!/usr/bin/env python3
"""
Script: 09_combine_env_metadata.py
Description: Combine Physicochemical_properties.xlsx (depth-aggregated means)
             and Depths.xlsx (per-station x depth metadata) into a single
             unified workbook with two sheets, normalizing depth labels to
             match MAGs_metadata.tsv conventions (1000, BTM, ROMZ).

             Sheet1 "Physicochemical_by_depth"  -> joined by Depth
             Sheet2 "Stations_depths"           -> joined by Station + Depth

             Output overwrites the input Physicochemical_properties.xlsx.
             A backup of the original is created next to it.

Requires pandas and openpyxl.
"""

import os
import shutil
from datetime import datetime

import pandas as pd

from pathlib import Path

# Project root, from the environment or the usual location.
PROJECT_DIR = Path(os.environ.get("MSC_ARCHAEA_DIR", Path.home() / "msc_archaea"))

# The environmental tables live under data/ in the published tree.
PHYSCHEM_XLSX = str(PROJECT_DIR / "data" / "Physicochemical_properties.xlsx")
DEPTHS_XLSX = str(PROJECT_DIR / "data" / "Depths.xlsx")
OUTPUT_XLSX = PHYSCHEM_XLSX  # overwrite

SHEET1_NAME = "Physicochemical_by_depth"
SHEET2_NAME = "Stations_depths"

# Depth label normalization (matches MAGs_metadata.tsv canonical values)
DEPTH_NORMALIZATION = {
    # From Physicochemical_properties.xlsx Sheet1 ("Depth" column)
    "ROMZ": "ROMZ",
    "1000 m": "1000",
    "1000m": "1000",
    "1000": "1000",
    "BTM": "BTM",
    "FONDO": "BTM",
    "Fondo": "BTM",
    # From Depths.xlsx ("Profundidad" column)
    "ZMOR": "ROMZ",
    "ZMOR1": "ROMZ",  # subscript "1" already removed by user; kept for safety
    "800 m": "800",
    "800m": "800",
}

# Spanish display labels for figures (R scripts will use these)
# Mapping canonical -> Spanish
DEPTH_ES_LABELS = {
    "ROMZ": "ZMOR",
    "1000": "1000 m",
    "BTM": "Fondo",
}

# Restrict to the MAG cohort. Anything outside these sets is dropped from
# the combined workbook to keep it tightly scoped to the thesis dataset.
MAG_STATIONS = {"A4", "A6", "C22", "F37", "PO1", "Y6"}
MAG_DEPTHS = {"1000", "BTM", "ROMZ"}

# Column rename map for Depths.xlsx -> tidy schema
DEPTHS_COL_RENAME = {
    "Estación": "Station",
    "Estacion": "Station",
    "Profundidad": "Depth",
    "Campaña oceanográfica": "Cruise",
    "Campania oceanografica": "Cruise",
    "Profundidad del fondo marino (m)": "Seafloor_depth_m",
    "Densidad": "Density_kg_m3",
    "Profundidad de colecta (m)": "Collection_depth_m",
}


# ==========================================
# HELPERS
# ==========================================

def normalize_depth(value):
    if value is None:
        return None
    s = str(value).strip()
    if s == "" or s.lower() == "nan":
        return None
    return DEPTH_NORMALIZATION.get(s, s)


def backup(path):
    if not os.path.exists(path):
        return None
    ts = datetime.now().strftime("%Y%m%d_%H%M%S")
    bak = f"{path}.bak_{ts}"
    shutil.copy2(path, bak)
    print(f"[backup] {path} -> {bak}")
    return bak


# ==========================================
# READERS
# ==========================================

def read_physchem(path):
    """Sheet1 of Physicochemical_properties.xlsx: 1 row per Depth aggregated."""
    df = pd.read_excel(path, sheet_name=0, engine="openpyxl")
    # Normalize Depth column
    if "Depth" not in df.columns:
        raise ValueError(f"Expected 'Depth' column in {path}; got {list(df.columns)}")
    df["Depth"] = df["Depth"].apply(normalize_depth)
    # Restrict to MAG cohort depths
    before = len(df)
    df = df[df["Depth"].isin(MAG_DEPTHS)].reset_index(drop=True)
    dropped = before - len(df)
    if dropped > 0:
        print(f"[info]   dropped {dropped} rows outside MAG depth cohort")
    # Ensure Spanish display column matches the canonical mapping
    df["Depth_es"] = df["Depth"].map(DEPTH_ES_LABELS)
    return df


def read_depths(path):
    """Sheet1 of Depths.xlsx: per-station x depth metadata.
    Skips the units pseudo-row (row 2 in the source) by dropping rows where
    the Station column is empty.
    """
    df = pd.read_excel(path, sheet_name=0, engine="openpyxl")
    # Drop the spacer rows the spreadsheet carries, where Station is blank
    station_col = None
    for cand in ("Estación", "Estacion", "Station"):
        if cand in df.columns:
            station_col = cand
            break
    if station_col is None:
        raise ValueError(f"Could not find Station column in {path}; got {list(df.columns)}")
    df = df[df[station_col].notna() & (df[station_col].astype(str).str.strip() != "")]
    # Rename to tidy schema
    df = df.rename(columns=DEPTHS_COL_RENAME)
    # Normalize Depth
    df["Depth"] = df["Depth"].apply(normalize_depth)
    # Drop rows with empty Depth after normalization
    df = df[df["Depth"].notna()]
    # Coerce numeric columns
    for col in ("Seafloor_depth_m", "Density_kg_m3", "Collection_depth_m"):
        if col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
    # Restrict to MAG cohort: stations and depths
    before = len(df)
    df = df[df["Station"].astype(str).isin(MAG_STATIONS)]
    df = df[df["Depth"].isin(MAG_DEPTHS)].reset_index(drop=True)
    dropped = before - len(df)
    if dropped > 0:
        print(f"[info]   dropped {dropped} rows outside MAG cohort "
              f"(station or depth)")
    # Spanish display label for figures
    df["Depth_es"] = df["Depth"].map(DEPTH_ES_LABELS)
    return df


# ==========================================
# MAIN
# ==========================================

def main():
    print(f"[info] reading {PHYSCHEM_XLSX}")
    physchem_df = read_physchem(PHYSCHEM_XLSX)
    print(f"[info]   rows: {len(physchem_df)}, cols: {list(physchem_df.columns)}")

    print(f"[info] reading {DEPTHS_XLSX}")
    depths_df = read_depths(DEPTHS_XLSX)
    print(f"[info]   rows: {len(depths_df)}, cols: {list(depths_df.columns)}")

    # Sanity check: report station coverage vs MAG metadata
    mag_stations = {"A4", "A6", "C22", "F37", "PO1", "Y6"}
    depths_stations = set(depths_df["Station"].astype(str).unique())
    print(f"[info] MAG stations missing in Depths.xlsx: "
          f"{sorted(mag_stations - depths_stations)}")
    print(f"[info] Depths.xlsx stations not in MAG dataset: "
          f"{sorted(depths_stations - mag_stations)}")

    # Sanity check: depth coverage
    physchem_depths = set(physchem_df["Depth"].dropna().astype(str).unique())
    depths_depths = set(depths_df["Depth"].dropna().astype(str).unique())
    print(f"[info] Depth labels in Physicochemical: {sorted(physchem_depths)}")
    print(f"[info] Depth labels in Depths: {sorted(depths_depths)}")
    expected_mag_depths = {"1000", "BTM", "ROMZ"}
    missing_in_depths = expected_mag_depths - depths_depths
    if missing_in_depths:
        print(f"[warn] MAG Depth labels missing from Depths.xlsx: "
              f"{sorted(missing_in_depths)} -- RDA join will produce NA for these")

    # Backup and write combined workbook
    backup(PHYSCHEM_XLSX)
    print(f"[info] writing combined workbook to {OUTPUT_XLSX}")
    with pd.ExcelWriter(OUTPUT_XLSX, engine="openpyxl") as writer:
        physchem_df.to_excel(writer, sheet_name=SHEET1_NAME, index=False)
        depths_df.to_excel(writer, sheet_name=SHEET2_NAME, index=False)

    print(f"[done] sheets written: {SHEET1_NAME!r}, {SHEET2_NAME!r}")
    print(f"[note] R reader (script 29):")
    print(f"  physchem <- readxl::read_xlsx('{OUTPUT_XLSX}', sheet = '{SHEET1_NAME}')")
    print(f"  stations <- readxl::read_xlsx('{OUTPUT_XLSX}', sheet = '{SHEET2_NAME}')")


if __name__ == "__main__":
    main()
