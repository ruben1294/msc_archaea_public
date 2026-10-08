#!/usr/bin/env python3
"""Build the MAG metadata table from the MAG filenames and the GTDB-Tk taxonomy.

The filenames encode where each MAG came from: <station>.<depth>.megahit.bin.<n>,
so A4.BTM.megahit.bin.25 is bin 25 from station A4 at the bottom depth. This
splits that apart and joins it against the classify_wf summary of step 03.

The curated result is data/MAGs_metadata.tsv, the spine of the whole project.
Every R script keys off its NameWithoutExtension and DisplayName columns, and
pipeline steps 04, 18, 20, 23 and 25 read it too. This script regenerates it
into results/ rather than over the curated copy, so the two can be compared:

    python3 scripts/utils/06_generate_metadata.py
    diff data/MAGs_metadata.tsv results/MAGs_metadata.generated.tsv

The oxygen minimum zone samples carry OMZ in their filenames and ROMZ in the
Depth, DisplayName and PanShortName columns, the label the rest of the analysis
uses. OriginalName keeps each MAG's name from before utils/05_rename_mags.sh
turned its underscores into dots.
"""

import argparse
import os
import re
import sys
from pathlib import Path

import pandas as pd

REPO_ROOT = Path(__file__).resolve().parents[2]
PROJECT_DIR = Path(os.environ.get("PROJECT_DIR")
                   or os.environ.get("MSC_ARCHAEA_DIR")
                   or REPO_ROOT)
RESULTS_DIR = Path(os.environ.get("RESULTS_DIR") or PROJECT_DIR / "results")

EXTENSION = ".fa"

# Filename depth code -> the label used everywhere else.
DEPTH_LABELS = {"OMZ": "ROMZ"}


def default_paths():
    filtered = os.environ.get("FILTERED_SUBDIR", "01_mag_quality/filtered")
    taxonomy = os.environ.get("TAXONOMY_SUBDIR", "03_taxonomy_assigned")
    prefix = os.environ.get("PREFIX", "tax_")
    return {
        "mag_dir": RESULTS_DIR / filtered,
        "originals_dir": Path(os.environ.get("MAGS_DIR") or PROJECT_DIR / "data" / "mags"),
        "taxonomy": RESULTS_DIR / taxonomy / "classify" / f"{prefix}.ar53.summary.tsv",
        "output": RESULTS_DIR / "MAGs_metadata.generated.tsv",
    }


def extract_taxonomy_level(taxonomy_string, level_prefix):
    """Return the taxon at one GTDB level ("g__", "s__"), or None if empty."""
    if pd.isna(taxonomy_string) or taxonomy_string == "":
        return None
    for level in taxonomy_string.split(";"):
        level = level.strip()
        if level.startswith(level_prefix):
            taxon = level[len(level_prefix):].strip()
            return taxon or None
    return None


def create_display_name(row):
    """Species, else genus, then station, depth and bin: the label on figures.

    For example "Thalassarchaeum sp002507125 A4 1000_38" or
    "MGIIb-Q1 A4 1000_28".
    """
    species = row.get("Species")
    genus = row.get("Genus")
    if species and species != "N/A":
        taxon_name = species
    elif genus and genus != "N/A":
        taxon_name = genus
    else:
        taxon_name = "Unclassified"
    return f"{taxon_name} {row['Station']} {row['Depth']}_{row['BinID']}"


def original_names(originals_dir):
    """Map renamed MAG names back to the input filenames in data/mags/.

    utils/05_rename_mags.sh only turns underscores into dots, so the renamed
    form is the key. Returns an empty map when the inputs are not at hand.
    """
    if not originals_dir.is_dir():
        return {}
    mapping = {}
    for path in originals_dir.glob(f"*{EXTENSION}"):
        stem = path.name[: -len(EXTENSION)]
        mapping[stem.replace("_", ".")] = stem
    return mapping


def guess_original_name(name):
    """Undo the rename for names of the form <station>_<depth>_megahit_bin.<n>."""
    return re.sub(r"_bin_(\d+)$", r"_bin.\1", name.replace(".", "_"))


def parse_mags(mag_dir, originals):
    """One row per MAG file, with the fields its name encodes."""
    files = sorted(f for f in os.listdir(mag_dir) if f.endswith(EXTENSION))
    if not files:
        print(f"[WARN] No {EXTENSION} files in {mag_dir}")
        return pd.DataFrame()

    rows = []
    for filename in files:
        name = filename[: -len(EXTENSION)]
        parts = name.split(".")

        if len(parts) >= 5:
            # Station.Depth.Assembler.bin.BinNumber
            station, depth, assembler = parts[0], parts[1], parts[2]
            bin_id_raw = ".".join(parts[4:])
        elif len(parts) == 4:
            # Same, but with no "bin" field
            station, depth, assembler = parts[0], parts[1], parts[2]
            bin_id_raw = parts[3]
        else:
            # Anything else: keep what can be read and flag the rest
            station = parts[0] if len(parts) > 0 else "unknown"
            depth = parts[1] if len(parts) > 1 else "unknown"
            assembler = parts[2] if len(parts) > 2 else "unknown"
            bin_id_raw = ".".join(parts[3:]) if len(parts) > 3 else "unknown"

        depth = DEPTH_LABELS.get(depth, depth)
        bin_id = bin_id_raw.replace("bin.", "")

        rows.append({
            "FullName": filename,
            "NameWithoutExtension": name,
            "Station": station,
            "Depth": depth,
            "Assembler": assembler,
            "BinID": bin_id,
            "OriginalName": originals.get(name, guess_original_name(name)),
            "PanShortName": f"{station}_{depth}_{bin_id}",
        })

    return pd.DataFrame(rows)


def add_taxonomy(df_mags, taxonomy_path):
    """Join the GTDB-Tk classification on the MAG name."""
    if not taxonomy_path.is_file():
        print(f"[WARN] Taxonomy file not found: {taxonomy_path}")
        print("  Continuing without taxonomy...")
        for col in ("FullTaxonomy", "Genus", "Species"):
            df_mags[col] = "N/A"
        return df_mags

    print(f"Loading taxonomy from: {taxonomy_path}")
    df_tax = pd.read_csv(taxonomy_path, sep="\t")
    mag_col, tax_col = df_tax.columns[0], df_tax.columns[1]
    print(f"[OK] Loaded taxonomy for {len(df_tax)} MAGs ('{mag_col}', '{tax_col}')")

    # GTDB-Tk reports the names it was given: dots when it runs after the
    # rename, as the pipeline now has it, underscores when it ran before.
    df_tax = df_tax.rename(columns={tax_col: "FullTaxonomy"})
    df_tax["NameWithoutExtension"] = df_tax[mag_col].astype(str).str.replace("_", ".", regex=False)
    df_tax["Genus"] = df_tax["FullTaxonomy"].apply(lambda x: extract_taxonomy_level(x, "g__"))
    df_tax["Species"] = df_tax["FullTaxonomy"].apply(lambda x: extract_taxonomy_level(x, "s__"))

    df_mags = df_mags.merge(
        df_tax[["NameWithoutExtension", "FullTaxonomy", "Genus", "Species"]],
        on="NameWithoutExtension",
        how="left",
    )
    for col in ("FullTaxonomy", "Genus", "Species"):
        df_mags[col] = df_mags[col].fillna("N/A")

    n = len(df_mags)
    print("[OK] Taxonomy merged:")
    print(f"  MAGs with taxonomy: {(df_mags['FullTaxonomy'] != 'N/A').sum()}/{n}")
    print(f"  MAGs with genus: {(df_mags['Genus'] != 'N/A').sum()}/{n}")
    print(f"  MAGs with species: {(df_mags['Species'] != 'N/A').sum()}/{n}")
    return df_mags


def main():
    defaults = default_paths()
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--mag-dir", type=Path, default=defaults["mag_dir"],
                        help="Filtered, renamed MAGs (default: %(default)s)")
    parser.add_argument("--originals-dir", type=Path, default=defaults["originals_dir"],
                        help="Input MAGs under their original names (default: %(default)s)")
    parser.add_argument("--taxonomy", type=Path, default=defaults["taxonomy"],
                        help="GTDB-Tk ar53 summary (default: %(default)s)")
    parser.add_argument("--output", type=Path, default=defaults["output"],
                        help="Where to write the table (default: %(default)s)")
    args = parser.parse_args()

    print("=" * 70)
    print("MAG metadata from filenames and taxonomy")
    print("=" * 70)

    if not args.mag_dir.is_dir():
        print(f"[FAIL] MAG directory not found: {args.mag_dir}")
        return 1
    if args.output.resolve() == (REPO_ROOT / "data" / "MAGs_metadata.tsv").resolve():
        print("[FAIL] Refusing to overwrite the curated data/MAGs_metadata.tsv; "
              "write elsewhere and compare.")
        return 1

    print(f"MAG directory: {args.mag_dir}")
    df_mags = parse_mags(args.mag_dir, original_names(args.originals_dir))
    if df_mags.empty:
        print("[FAIL] No MAG files processed")
        return 1
    print(f"[OK] Parsed {len(df_mags)} MAG files")

    df_mags = add_taxonomy(df_mags, args.taxonomy)
    df_mags["DisplayName"] = df_mags.apply(create_display_name, axis=1)

    column_order = [
        "FullName", "NameWithoutExtension", "DisplayName", "Station", "Depth",
        "Assembler", "BinID", "Genus", "Species", "FullTaxonomy",
        "OriginalName", "PanShortName",
    ]
    df_mags = df_mags[column_order]

    args.output.parent.mkdir(parents=True, exist_ok=True)
    df_mags.to_csv(args.output, sep="\t", index=False)

    print(f"[OK] {len(df_mags)} MAGs written to {args.output}")
    print(f"  Stations: {', '.join(sorted(df_mags['Station'].unique()))}")
    print(f"  Depths: {', '.join(sorted(df_mags['Depth'].unique()))}")
    print("Example display names:")
    for name in df_mags["DisplayName"].head(5):
        print(f"  {name}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
