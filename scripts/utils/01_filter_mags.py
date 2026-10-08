#!/usr/bin/env python3
"""Copy the MAGs that clear the MIMAG thresholds into results/01_mag_quality/filtered.

Reads the CheckM2 quality_report.tsv from pipeline step 01 and sorts each MAG
into one of three MIMAG bands:

    high    completeness >= 90 and contamination <= 5
    medium  completeness >= 50 and contamination < 10
    low     anything else, which is left behind

Everything from step 02 onward reads the filtered directory, so this is what
fixes the cohort. Of the 64 MAGs originally recovered, 54 clear the medium bar.
Run utils/05_rename_mags.sh right after it, before step 02.

Paths follow the project layout: the input MAGs in data/mags and the results
under results/, both relative to the repository root. PROJECT_DIR, MAGS_DIR and
RESULTS_DIR from configs/server.env move them when the configuration is loaded,
MSC_ARCHAEA_DIR moves the root, and every path can be given on the command line.
The thresholds come from the same configuration, with the MIMAG values as the
fallback.
"""

import argparse
import logging
import os
import shutil
import sys
from datetime import datetime
from pathlib import Path

import pandas as pd

# MIMAG bands (Bowers et al. 2017 Nat Biotechnol). Medium quality takes strictly
# less than 10% contamination, which is why a MAG at 9.9% is kept and one at
# exactly 10% is not.
MIN_COMP_HIGH = float(os.environ.get("HIGH_QUALITY_COMPLETENESS", 90))
MAX_CONT_HIGH = float(os.environ.get("HIGH_QUALITY_CONTAMINATION", 5))
MIN_COMP_MEDIUM = float(os.environ.get("MIN_COMPLETENESS", 50))
MAX_CONT_MEDIUM = float(os.environ.get("MAX_CONTAMINATION", 10))

# scripts/utils/ sits two levels below the repository root.
REPO_ROOT = Path(__file__).resolve().parents[2]

logger = logging.getLogger("filter_mags")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--base", type=Path,
                        help="project root (default: $PROJECT_DIR, $MSC_ARCHAEA_DIR "
                             "or the repository this script lives in)")
    parser.add_argument("--raw-dir", type=Path,
                        help="directory holding the unfiltered MAGs "
                             "(default: <base>/data/mags)")
    parser.add_argument("--out-dir", type=Path,
                        help="destination for the MAGs that pass "
                             "(default: <base>/results/01_mag_quality/filtered)")
    parser.add_argument("--report", type=Path,
                        help="CheckM2 quality_report.tsv "
                             "(default: <base>/results/01_mag_quality/quality_report.tsv)")
    parser.add_argument("--log-dir", type=Path,
                        help="where to write the run log (default: <base>/logs)")

    args = parser.parse_args()
    if args.base is not None:
        # An explicit root wins over the configuration.
        results = args.base / "results"
        mags = args.base / "data" / "mags"
    else:
        args.base = Path(os.environ.get("PROJECT_DIR")
                         or os.environ.get("MSC_ARCHAEA_DIR")
                         or REPO_ROOT)
        results = Path(os.environ.get("RESULTS_DIR") or args.base / "results")
        mags = Path(os.environ.get("MAGS_DIR") or args.base / "data" / "mags")
    quality_dir = results / os.environ.get("QUALITY_SUBDIR", "01_mag_quality")
    args.raw_dir = args.raw_dir or mags
    args.out_dir = args.out_dir or results / os.environ.get("FILTERED_SUBDIR", "01_mag_quality/filtered")
    args.report = args.report or quality_dir / "quality_report.tsv"
    args.log_dir = args.log_dir or args.base / os.environ.get("LOGS_SUBDIR", "logs")
    return args


def setup_logging(log_dir: Path) -> Path:
    log_dir.mkdir(parents=True, exist_ok=True)
    log_file = log_dir / f"01_filter_MAGs_{datetime.now():%Y%m%d_%H%M%S}.log"
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s | %(levelname)s | %(message)s",
        handlers=[logging.FileHandler(log_file), logging.StreamHandler(sys.stdout)],
    )
    return log_file


def classify(completeness: float, contamination: float) -> tuple[bool, str]:
    """Return whether the MAG is worth keeping, and which MIMAG band it lands in."""
    if completeness >= MIN_COMP_HIGH and contamination <= MAX_CONT_HIGH:
        return True, "HIGH"
    if completeness >= MIN_COMP_MEDIUM and contamination < MAX_CONT_MEDIUM:
        return True, "MEDIUM"
    return False, "LOW"


def main() -> int:
    args = parse_args()
    log_file = setup_logging(args.log_dir)

    try:
        logger.info("MAG quality filtering")
        logger.info("  Raw MAGs:       %s", args.raw_dir)
        logger.info("  Output:         %s", args.out_dir)
        logger.info("  Quality report: %s", args.report)

        for path, label in ((args.raw_dir, "raw MAGs directory"),
                            (args.report, "quality report")):
            if not path.exists():
                raise FileNotFoundError(f"{label} not found: {path}")

        args.out_dir.mkdir(parents=True, exist_ok=True)

        df = pd.read_csv(args.report, sep="\t")
        logger.info("Loaded %d MAGs from the quality report", len(df))

        required = {"Name", "Completeness", "Contamination"}
        missing = required - set(df.columns)
        if missing:
            raise RuntimeError(f"Quality report is missing columns: {sorted(missing)}")

        kept = {"HIGH": 0, "MEDIUM": 0}
        rejected = 0
        missing_files = []

        for _, row in df.iterrows():
            mag_name = str(row["Name"])
            completeness = row["Completeness"]
            contamination = row["Contamination"]

            # CheckM2 reports the name without the extension.
            if not mag_name.endswith(".fa"):
                mag_name += ".fa"

            acceptable, band = classify(completeness, contamination)

            if not acceptable:
                rejected += 1
                logger.debug("[LOW] %s (comp %.2f%%, cont %.2f%%)",
                             mag_name, completeness, contamination)
                continue

            src = args.raw_dir / mag_name
            if not src.exists():
                missing_files.append(mag_name)
                logger.warning("Listed in the report but not on disk: %s", mag_name)
                continue

            shutil.copy2(src, args.out_dir / mag_name)
            kept[band] += 1
            logger.info("[%s] %s (comp %.2f%%, cont %.2f%%)",
                        band, mag_name, completeness, contamination)

        logger.info("Summary")
        logger.info("  MAGs analysed:  %d", len(df))
        logger.info("  High quality:   %d  (comp >= %d%%, cont <= %d%%)",
                    kept["HIGH"], MIN_COMP_HIGH, MAX_CONT_HIGH)
        logger.info("  Medium quality: %d  (comp >= %d%%, cont < %d%%)",
                    kept["MEDIUM"], MIN_COMP_MEDIUM, MAX_CONT_MEDIUM)
        logger.info("  Kept in total:  %d", kept["HIGH"] + kept["MEDIUM"])
        logger.info("  Rejected:       %d", rejected)
        logger.info("  Files missing:  %d", len(missing_files))

        for name in missing_files:
            logger.warning("  missing: %s", name)

        logger.info("Filtered MAGs: %s", args.out_dir)
        logger.info("Log file:      %s", log_file)
        return 0

    except Exception as exc:
        logger.error("%s: %s", type(exc).__name__, exc, exc_info=True)
        return 1

if __name__ == "__main__":
    sys.exit(main())
