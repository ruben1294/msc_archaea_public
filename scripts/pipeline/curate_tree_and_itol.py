#!/usr/bin/env python3
"""
curate_tree_and_itol.py

Curates the GTDB-Tk de_novo_wf phylogenomic tree for 54 archaeal MAGs:

  1. Selects GTDB reference genomes based on ANI proximity and
     phylogenetic coverage (see REFERENCE_SELECTION_RATIONALE below).
  2. Extracts alignment columns for selected taxa from the GTDB-Tk MSA.
  3. Re-infers a tree with FastTree (LG+GAMMA) on the curated alignment.
  4. Renames all leaf labels (MAGs + refs) so the tree and annotation files
     share exactly the same identifiers (required by iTOL).
  5. Generates iTOL annotation files matching the format of
     generate_itol_annotations.py (files 01–10) for the full 54-MAG set.

Leaf-name convention (same in tree AND annotation files):
  MAGs  to DisplayName from MAGs_metadata.tsv with spaces become underscores,
           e.g. "MGIIb-Q1_A4_1000_28"  (no Newick quoting needed)
  Refs  to REF_DISPLAY label, e.g. "Thalassarchaeum_betae_B"
           (underscores throughout)

Reference Selection Methodology
Based on Parks et al. (2022) GTDB r207 and Chaumeil et al. (2022) GTDB-Tk v2:
  • Per genus: 2–3 GTDB references with highest ANI to the MAGs of that genus,
    prioritising species cluster representatives and named type strains.
  • Full genus coverage: if a genus has >1 GTDB species, include ≥1 ref per
    species to anchor within-genus diversity.
  • Outgroup: 1 Halobacteriota genome (Halobacterium salinarum NRC-1,
    RS_GCF_000025685.1), a well-characterised phylum used routinely as
    outgroup for Thermoplasmatota/Thermoproteota trees.

Usage, in the gtdbtk conda environment of env/gtdbtk.yml, which provides
FastTree:
  python3 scripts/pipeline/curate_tree_and_itol.py [--project-dir PATH]

Output tree:  results/SO_02/phylogenomics/curated_54mags.fasttree.nwk
iTOL files:   results/SO_02/phylogenomics/itol_annotations/

Author : Rubén Castañeda-Martínez
Created: 2026-03-27
"""

import argparse
import csv
import gzip
import logging
import os
import re
import subprocess
import sys

# SELF-BOOTSTRAP: re-exec inside the conda env that provides FastTree, so the
# script also runs from an editor's own Python interpreter. CURATE_CONDA_ENV
# names the environment; nothing happens when FastTree is already on PATH.
_CONDA_ENV = os.environ.get("CURATE_CONDA_ENV", "gtdbtk")


def _ensure_conda_env() -> None:
    import shutil
    if os.environ.get("CONDA_DEFAULT_ENV") == _CONDA_ENV or shutil.which("FastTree"):
        return
    conda_exe = os.environ.get("CONDA_EXE") or ""
    if not conda_exe:
        for candidate in [
            os.path.expanduser("~/miniconda3/bin/conda"),
            os.path.expanduser("~/anaconda3/bin/conda"),
            os.path.expanduser("~/miniforge3/bin/conda"),
            "/opt/conda/bin/conda",
        ]:
            if os.path.isfile(candidate):
                conda_exe = candidate
                break
    if not conda_exe:
        print("[WARN] conda not found — proceeding without env switch "
              "(FastTree may not be in PATH)", flush=True)
        return
    print(f"[INFO] Re-execing inside conda env '{_CONDA_ENV}' …", flush=True)
    os.execv(conda_exe, [conda_exe, "run", "--no-capture-output",
                         "-n", _CONDA_ENV, "python3", __file__] + sys.argv[1:])

_ensure_conda_env()

# LOGGING
logging.basicConfig(
    format="[%(asctime)s] %(levelname)s  %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
    level=logging.INFO,
    stream=sys.stdout,
)
log = logging.getLogger(__name__)

# PATHS (results resolved relative to --project-dir; the metadata table is part
# of the repository and always read from it)
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
DEFAULT_PROJECT_DIR = (os.environ.get("PROJECT_DIR")
                       or os.environ.get("MSC_ARCHAEA_DIR")
                       or REPO_ROOT)


def build_paths(project_dir: str) -> dict:
    r = project_dir
    return {
        "de_novo_tree":    os.path.join(r, "results", "11_de_novo_tree"),
        "msa_full_gz":     os.path.join(r, "results", "11_de_novo_tree", "align",
                                        "54_archaea.ar53.msa.fasta.gz"),
        "taxonomy_tsv":    os.path.join(r, "results", "03_taxonomy_assigned",
                                        "tax_.ar53.summary.tsv"),
        "metadata_tsv":    os.path.join(REPO_ROOT, "data", "MAGs_metadata.tsv"),
        "quality_tsv":     os.path.join(r, "results", "01_mag_quality",
                                        "quality_report.tsv"),
        "product_tsv":     os.path.join(r, "results", "10_MAG_distill",
                                        "VISUALIZE", "product.tsv"),
        "output_dir":      os.path.join(r, "results", "SO_02", "phylogenomics"),
        "itol_dir":        os.path.join(r, "results", "SO_02", "phylogenomics",
                                        "itol_annotations"),
    }

# REFERENCE SELECTION
# Each tuple: (GTDB_leaf_in_de_novo_tree, display_label, genus, justification)
#
# GTDB-Tk prefixes: GB_GCA_* (GenBank) or RS_GCF_* (RefSeq).
# ANI values are taken from tax_.ar53.summary.tsv column
# "other_related_references"; the closest_genome_reference column is also used.
#
# Literature basis:
#   Parks et al. (2022) Nat Biotechnol, GTDB classification method
#   Chaumeil et al. (2022) Bioinformatics, GTDB-Tk v2
#   Rinke et al. (2019) Nat Microbiol, Marine Group II/III Euryarchaeota
#   Orellana et al. (2022) ISME J, Thalassarchaeum betae type species
#   Vanni et al. (2022) eLife, MAG annotation best practices

SELECTED_REFERENCES = [
    # ---- Thalassarchaeum (Thermoplasmatota / Poseidoniales / Thalassarchaeaceae) ----
    # 16 MAGs; closest by ANI for the majority (~92–93 %).
    # GCF_022448665.1 = Thalassarchaeum sp002507125 (RefSeq, highest completeness)
    ("RS_GCF_022448665.1",
     "Thalassarchaeum_sp002507125",
     "Thalassarchaeum",
     "Closest GTDB reference (ANI ~93 %) for most Thalassarchaeum MAGs; "
     "species cluster representative (Parks et al. 2022)"),

    # GCA_022447095.1 = Thalassarchaeum betae_B (second lineage within genus)
    ("GB_GCA_022447095.1",
     "Thalassarchaeum_betae_B",
     "Thalassarchaeum",
     "Second-closest ANI (~92–93 %); represents the betae_B species cluster, "
     "ensuring within-genus diversity coverage"),

    # ---- MGIIb-Q1 (Thermoplasmatota / Poseidoniales / Thalassarchaeaceae) ----
    # 7 MAGs; no species assignment (novel species cluster).
    # GCA_963580235.1, closest placement reference (cp_ref) for A4_1000_bin.28
    # NOTE: GCA_002505495.1 and GCA_012964335.1 absent from GTDB r226 MSA;
    #       GCA_003193925.1 (Thalassarchaeum_betae) also absent from MSA.
    #       MGIIb-O1 sole ref (GCA_002498525.1) likewise absent from MSA.
    ("GB_GCA_963580235.1",
     "MGIIb-Q1_sp963580235",
     "MGIIb-Q1",
     "Closest placement reference (ANI 93–95 %) for 3 of 7 MGIIb-Q1 MAGs "
     "(A4_1000_bin.28, A6_BTM_bin.49, Y6_1000_bin.23)"),

    # ---- JAKUTA01 (Thermoplasmatota / Poseidoniales / Thalassarchaeaceae) ----
    # 2 MAGs; newly described genus in Thalassarchaeaceae.
    ("GB_GCA_950068635.1",
     "JAKUTA01_sp950068635",
     "JAKUTA01",
     "Closest genome reference by ANI (~93 %) for both JAKUTA01 MAGs; "
     "type representative of this genus"),

    ("GB_GCA_022447885.1",
     "JAKUTA01_sp022447885",
     "JAKUTA01",
     "Second-closest ANI (~93 %) for JAKUTA01 MAGs; represents an "
     "independent JAKUTA01 lineage for within-genus resolution"),

    # ---- UBA36 (Thermoplasmatota / Poseidoniales / Thalassarchaeaceae) ----
    # 2 MAGs; genus within Thalassarchaeaceae (Parks et al. 2022).
    ("GB_GCA_036153975.1",
     "UBA36_sp036153975",
     "UBA36",
     "Closest GTDB reference by ANI (~93 %) for both UBA36 MAGs; "
     "recently deposited UBA36 genome"),

    ("GB_GCA_036336595.1",
     "UBA36_sp036336595",
     "UBA36",
     "Second-closest ANI (~92 %); independent UBA36 species cluster, "
     "ensures within-genus diversity"),

    # ---- UBA8886 (Thermoplasmatota / Poseidoniia / MGIII / CG-Epi1) ----
    # 6 MAGs; related to Marine Group III (MGIII) archaea.
    # NOTE: GCA_003193815.1 (UBA8886_sp003193815) absent from GTDB r226 MSA.
    ("GB_GCA_030754655.1",
     "UBA8886_sp030754655",
     "UBA8886",
     "Highest ANI (~94–95 %) for the majority of UBA8886 MAGs; "
     "best available UBA8886 representative"),

    ("GB_GCA_963578705.1",
     "UBA8886_sp963578705",
     "UBA8886",
     "Second-highest ANI (~92–93 %) across UBA8886 MAGs; "
     "independent species cluster"),

    # ---- UBA57 (Thermoproteota / Nitrososphaeria / Nitrososphaerales / UBA57) ----
    # 4 MAGs; novel family within Nitrososphaerales.
    # NOTE: GCA_013306035.1 (UBA57_sp013306035) and GCA_013378385.1
    #       (UBA57_sp013378385) absent from GTDB r226 MSA.
    ("GB_GCA_963577725.1",
     "UBA57_sp963577725",
     "UBA57",
     "Highest ANI (~94–95 %) for 3 of 4 UBA57 MAGs; "
     "best available UBA57 representative in GTDB r226 MSA"),

    # ---- UBA60 (Thermoproteota / Nitrososphaeria / Nitrososphaerales / UBA60) ----
    # 4 MAGs; sister to UBA57 within Nitrososphaerales.
    # NOTE: GCA_002503395.1 (UBA60_sp002503395) absent from GTDB r226 MSA.
    ("GB_GCA_028283345.1",
     "UBA60_sp028283345",
     "UBA60",
     "Best available UBA60 reference in GTDB r226 MSA (~89–92 % ANI); "
     "GCA_002503395.1 (primary cg_ref) is absent from this MSA"),

    # ---- Nitrosopelagicus / Undinarchaeum ----
    # NOTE: All Nitrosopelagicus refs (GCA_013390375.1, GCA_000402075.1) and
    #       all Undinarchaeum refs (GCA_014189675.1, GCA_002495465.1) are absent
    #       from the GTDB r226 de_novo_wf MSA. Their 3+3 MAGs appear in the
    #       tree without a reference genome from their own genus.

    # ---- Outgroup (Halobacteriota / Halobacteria / Halobacteriales) ----
    # Halobacterium salinarum NRC-1: well-characterised, complete genome,
    # widely used as outgroup for Thermoplasmatota/Thermoproteota trees.
    # In the de_novo_wf this phylum defines the root (--outgroup_taxon
    # p__Euryarchaeota, now p__Halobacteriota in GTDB r220).
    ("RS_GCF_000025685.1",
     "Halobacterium_salinarum_NRC-1",
     "Outgroup_Halobacteriota",
     "Halobacterium salinarum NRC-1 (complete genome, RefSeq); "
     "Halobacteriota is the outgroup phylum used in the de_novo_wf rooting; "
     "standard reference for archaeal phylogenomics (Ng et al. 2000)"),
]

# Set of GTDB leaf IDs (as they appear in the tree) for quick lookup
SELECTED_REF_IDS = {r[0] for r in SELECTED_REFERENCES}

# Display name map: GTDB_leaf_id to display_label
REF_DISPLAY = {r[0]: r[1] for r in SELECTED_REFERENCES}

# Genus map: GTDB_leaf_id to genus
REF_GENUS = {r[0]: r[2] for r in SELECTED_REFERENCES}

# COLORS AND TAXONOMY
# Genus-level color palette (consistent with DRAM2 analysis heatmaps)
GENUS_COLORS = {
    # Thalassarchaeaceae / Poseidoniales
    "Thalassarchaeum":           "#1d6a96",
    "MGIIb-Q1":                  "#e76f51",
    "MGIIb-O1":                  "#f4a261",
    "JAKUTA01":                  "#e9c46a",
    "UBA36":                     "#2a9d8f",
    # MGIII
    "UBA8886":                   "#8ecae6",
    # Thermoproteota / Nitrososphaerales
    "UBA57":                     "#6a4c93",
    "UBA60":                     "#9b72cf",
    "Nitrosopelagicus":          "#264653",
    # Undinarchaeota
    "Undinarchaeum":             "#a8dadc",
    # Outgroup
    "Outgroup_Halobacteriota":   "#adb5bd",
}

# Station colors (same as existing scripts)
STATION_COLORS = {
    "A4":  "#e63946",
    "A6":  "#2a9d8f",
    "PO1": "#e9c46a",
    "Y6":  "#6a4c93",
    "F37": "#457b9d",
    "C22": "#f4a261",
}

# Depth colors: code to (hex, Spanish label)
# OMZ samples carry Depth="ROMZ" in MAGs_metadata.tsv, not "OMZ". ROMZ is the
# Spanish abbreviation for the oxygen minimum zone.
DEPTH_COLORS = {
    "BTM":  ("#1d3557", "Fondo"),
    "1000": ("#a8dadc", "1000 m"),
    "ROMZ": ("#e9c46a", "ZMOR"),
}

# ---- DRAM2 metabolic field definitions (from product.tsv) ----
DRAM2_FIELDS = [
    ("Fijación de carbono (CBB)",
     ["Reductive pentose phosphate cycle (Calvin cycle)"], 0.5),
    ("Fijación de carbono (rTCA)",
     ["Reductive citrate cycle (Arnon-Buchanan cycle)"], 0.5),
    ("Fijación de carbono (WL)",
     ["Reductive acetyl-CoA pathway (Wood-Ljungdahl pathway)"], 0.5),
    ("Metanogénesis",
     ["Methanogenesis, CO2 => methane"], 0.5),
    ("Fijación de nitrógeno",
     ["Nitrogen: nitrogen => ammonia"], None),
    ("Oxidación de amoniaco",
     ["Nitrogen: Bacterial/Archaeal ammonia oxidation",
      "Nitrogen: ammonia => nitrite"], None),
    ("Reducción disimatoria de nitrato",
     ["Nitrogen: nitrate => nitrite",
      "Nitrogen: Dissimilatory nitrite reduction to ammonia (DNRA)"], None),
    ("Reducción asimilatoria de nitrato",
     ["Nitrogen: nitrate => nitrite"], None),
    ("Desnitrificación",
     ["Nitrogen: nitrite => nitric oxide",
      "Nitrogen: nitric oxide => nitrous oxide",
      "Nitrogen: nitrous oxide => nitrogen"], None),
    ("Oxidación de azufre",
     ["Sulfur: Thiosulfate oxidation by SOX complex, thiosulfate => sulfate"], None),
    ("Reducción de sulfato",
     ["Sulfur: dissimilatory sulfate reduction (and oxidation) sulfate => sulfide"],
     None),
    ("Respiración aeróbica",
     ["Complex IV L: Cytochrome c oxidase",
      "Complex IV L: Cytochrome c oxidase, prokaryotes",
      "Complex IV H: Cytochrome c oxidase, cbb3-type",
      "Complex IV L: Cytochrome o ubiquinol oxidase",
      "Complex IV H: Cytochrome bd ubiquinol oxidase"], 0.5),
    ("Fermentación",
     ["SCFA/alcohol conversions: Alcohol production",
      "SCFA/alcohol conversions: Butyrate, pt 1",
      "SCFA/alcohol conversions: Butyrate, pt 2",
      "SCFA/alcohol conversions: lactate D",
      "SCFA/alcohol conversions: lactate L",
      "SCFA/alcohol conversions: acetate, pt 1",
      "SCFA/alcohol conversions: acetate, pt 2"], None),
    ("Fototrofía",
     ["Photo: Photosystem I",
      "Photo: Photosystem II"], None),
]

DRAM2_FIELD_COLORS = [
    "#2a9d8f",  # CBB
    "#2a9d8f",  # rTCA
    "#2a9d8f",  # WL
    "#6a4c93",  # Methanogenesis
    "#e9c46a",  # N fixation
    "#e9c46a",  # Ammonia oxidation
    "#e9c46a",  # Dissim. NO3
    "#e9c46a",  # Assim. NO3
    "#e9c46a",  # Denitrification
    "#f4a261",  # Sulfur oxidation
    "#f4a261",  # Sulfate reduction
    "#264653",  # Aerobic respiration
    "#264653",  # Fermentation
    "#e76f51",  # Phototrophy
]

# UTILITIES


def tsv_name_to_dot(tsv_name: str) -> str:
    """A4_BTM_megahit_bin.28  to  A4.BTM.megahit.bin.28"""
    m = re.match(r'^([^_]+)_([^_]+)_([^_]+)_bin\.(.+)$', tsv_name)
    if m:
        return f"{m.group(1)}.{m.group(2)}.{m.group(3)}.bin.{m.group(4)}"
    return tsv_name


def dot_name_to_display(dot_name: str, display_map: dict) -> str:
    """Returns the display name (spaces become underscores) from metadata, or dot_name."""
    return display_map.get(dot_name, dot_name.replace(".", "_"))


def _open_out(path: str):
    return open(path, "w", newline="\n", encoding="utf-8")

# DATA LOADING


def load_metadata(path: str):
    """
    Returns:
        mag_order  : list[str] , NameWithoutExtension in file order
        stations   : dict[str, str]
        depths     : dict[str, str]
        display    : dict[str, str] , display name with spaces become underscores
        genera     : dict[str, str]
    """
    mag_order, stations, depths, display, genera = [], {}, {}, {}, {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            name = row["NameWithoutExtension"].strip()
            if not name:
                continue
            mag_order.append(name)
            stations[name] = row.get("Station", "").strip() or name.split(".")[0]
            d = row.get("Depth", "").strip()
            depths[name] = d or name.split(".")[1]
            disp = row.get("DisplayName", "").strip()
            display[name] = disp if disp else name.replace(".", "_")
            genera[name] = row.get("Genus", "Unknown").strip()
    return mag_order, stations, depths, display, genera


def load_quality(path: str, target: set):
    """Returns completeness, contamination dicts keyed by dot_name."""
    comp, cont = {}, {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            name = tsv_name_to_dot(row.get("Name", "").strip())
            if name not in target:
                continue
            try:
                comp[name] = float(row["Completeness"])
                cont[name] = float(row["Contamination"])
            except (KeyError, ValueError):
                log.warning("Quality: could not parse %s", name)
    return comp, cont


def load_ani(path: str, target: set):
    """Best ANI from closest_genome_ani to closest_placement_ani, else -1."""
    ani = {}
    with open(path, newline="", encoding="utf-8") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            name = tsv_name_to_dot(row.get("user_genome", "").strip())
            if name not in target:
                continue
            val = -1.0
            for col in ("closest_genome_ani", "closest_placement_ani"):
                raw = row.get(col, "").strip()
                if raw and raw not in ("N/A", "", "None", "nan"):
                    try:
                        val = float(raw)
                        break
                    except ValueError:
                        pass
            ani[name] = val
    return ani


def load_product(path: str, mag_order: list, target: set):
    """Returns dict mag_name to {col: value} for target MAGs."""
    data = {}
    with open(path, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for i, row in enumerate(reader):
            if i >= len(mag_order):
                break
            name = mag_order[i]
            if name in target:
                data[name] = dict(row)
    return data


def eval_field(row: dict, columns: list, threshold) -> int:
    for col in columns:
        raw = row.get(col, "").strip()
        if not raw:
            continue
        if threshold is None:
            if raw.lower() in ("true", "1"):
                return 1
        else:
            try:
                if float(raw) >= threshold:
                    return 1
            except ValueError:
                pass
    return 0

# REFERENCE SELECTION REPORT


def write_reference_selection_table(out_path: str):
    """Writes a TSV documenting the selected references and their justification."""
    with _open_out(out_path) as fh:
        fh.write("gtdb_leaf_id\tdisplay_label\tgenus\tjustification\n")
        for leaf_id, label, genus, reason in SELECTED_REFERENCES:
            fh.write(f"{leaf_id}\t{label}\t{genus}\t{reason}\n")
    log.info("Reference selection table to %s", out_path)

# MSA EXTRACTION


def extract_msa(msa_gz_path: str, keep_ids: set, out_fasta_path: str) -> int:
    """
    Reads the gzipped full MSA (MAGs + GTDB refs) and writes a FASTA
    containing only sequences whose header is in keep_ids.

    GTDB-Tk MSA headers for refs are 'GB_GCA_*' or 'RS_GCF_*';
    for MAGs they are 'A4.1000.megahit.bin.28' etc.

    Returns the number of sequences written.
    """
    written = 0
    current_id = None
    current_seq = []

    def flush(fh):
        nonlocal written
        if current_id and current_id in keep_ids:
            fh.write(f">{current_id}\n")
            fh.write("".join(current_seq) + "\n")
            written += 1

    with gzip.open(msa_gz_path, "rt") as gz, open(out_fasta_path, "w") as out:
        for line in gz:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush(out)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        flush(out)

    log.info("MSA extraction: %d sequences written to %s", written, out_fasta_path)
    return written

# TREE INFERENCE (FastTree)


def run_fasttree(msa_fasta: str, out_tree: str, log_file: str,
                 n_threads: int = 4, prot_model: str = "LG") -> None:
    """
    Runs FastTree on the curated MSA.
    Model: LG + CAT approximation (default for FastTree protein mode).
    """
    cmd = [
        "FastTree",
        "-lg",       # LG substitution model
        "-gamma",    # Gamma rate variation
        "-log", log_file,
        msa_fasta,
    ]
    log.info("Running FastTree: %s", " ".join(cmd))
    with open(out_tree, "w") as tree_out, open(log_file, "w") as log_out:
        proc = subprocess.run(
            cmd,
            stdout=tree_out,
            stderr=subprocess.PIPE,
            text=True,
        )
        log_out.write(proc.stderr)

    if proc.returncode != 0:
        log.error("FastTree failed (exit %d). See %s", proc.returncode, log_file)
        sys.exit(proc.returncode)
    log.info("FastTree completed to %s", out_tree)

# TREE PRUNING (DendroPy), used to create a GTDB-topology pruned tree
# as a complement/validation of the re-inferred FastTree


def prune_decorated_tree(decorated_tree_path: str, keep_ids: set,
                          out_path: str) -> None:
    """
    Uses DendroPy to prune the GTDB-Tk decorated tree.
    Retains only leaves in keep_ids and writes the pruned Newick.
    """
    try:
        import dendropy
    except ImportError:
        log.warning("DendroPy not available; skipping pruned-original-tree output.")
        return

    log.info("Loading decorated tree with DendroPy …")
    tree = dendropy.Tree.get(
        path=decorated_tree_path,
        schema="newick",
        preserve_underscores=True,
    )

    # Build taxon set to retain
    taxa_to_keep = set()
    for taxon in tree.taxon_namespace:
        if taxon.label in keep_ids:
            taxa_to_keep.add(taxon)

    missing = keep_ids - {t.label for t in taxa_to_keep}
    if missing:
        log.warning("Leaves not found in decorated tree: %s", missing)

    log.info("Pruning tree: retaining %d of %d leaves …",
             len(taxa_to_keep), len(tree.taxon_namespace))
    tree.retain_taxa(taxa_to_keep)
    tree.purge_taxon_namespace()

    tree.write(path=out_path, schema="newick")
    log.info("Pruned GTDB-topology tree to %s", out_path)

# DISPLAY NAME MAPPING for iTOL


def build_leaf_display_map(mag_order: list, display: dict) -> dict:
    """
    Returns dot_name to leaf_id for all leaves.
    MAGs: dot_name to DisplayName with spaces replaced by underscores,
          e.g. A4.1000.megahit.bin.28 to MGIIb-Q1_A4_1000_28.
    Refs: GTDB_leaf_id to REF_DISPLAY label (already uses underscores).
    The same IDs are used in both the tree and annotation files (iTOL requirement).
    """
    m = {mag: display[mag].replace(" ", "_") for mag in mag_order}
    m.update(REF_DISPLAY)
    return m


def rename_tree_leaves(tree_path: str, out_path: str,
                       dot_to_display: dict) -> None:
    """
    Renames leaf labels in a Newick file using the provided mapping:
      - MAGs: dot-notation (A4.1000.megahit.bin.28) to underscore ID (MGIIb-Q1_A4_1000_28)
      - Refs: GTDB prefix (GB_GCA_022447095.1) to REF_DISPLAY label (Thalassarchaeum_betae_B)

    All output labels use underscores (no spaces), so no Newick quoting is needed
    and IDs match the annotation files exactly (iTOL requirement).
    """
    with open(tree_path, encoding="utf-8") as fh:
        tree_str = fh.read()
    for dot_name, leaf_id in dot_to_display.items():
        # In Newick, leaf names are immediately followed by ':'
        tree_str = tree_str.replace(f"{dot_name}:", f"{leaf_id}:")
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(tree_str)
    log.info("Leaf labels renamed to %s", out_path)


def translate(name: str, leaf_map: dict) -> str:
    return leaf_map.get(name, name)

# iTOL FILE WRITERS


def write_01_source_colorstrip(out_dir: str, mag_order: list,
                                leaf_map: dict) -> str:
    path = os.path.join(out_dir, "01_source_colorstrip.txt")
    with _open_out(path) as f:
        f.write("DATASET_COLORSTRIP\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tFuente\n")
        f.write("COLOR\t#e63946\n")
        f.write("STRIP_WIDTH\t25\n")
        f.write("MARGIN\t5\n")
        f.write("BORDER_WIDTH\t1\n")
        f.write("BORDER_COLOR\t#ffffff\n")
        f.write("SHOW_STRIP_LABELS\t0\n")
        f.write("COLOR_BRANCHES\t0\n")
        f.write("LEGEND_TITLE\tFuente\n")
        f.write("LEGEND_SHAPES\t1\t1\n")
        f.write("LEGEND_COLORS\t#e63946\t#457b9d\n")
        f.write("LEGEND_LABELS\tMAG\tGenoma de referencia (GTDB)\n")
        f.write("DATA\n")
        for mag in mag_order:
            f.write(f"{translate(mag, leaf_map)}\t#e63946\tMAG\n")
        for ref_id in SELECTED_REF_IDS:
            f.write(f"{translate(ref_id, leaf_map)}\t#457b9d\tReferencia\n")
    return path


def write_02_taxonomy_colorstrip(out_dir: str, mag_order: list,
                                  genera: dict, leaf_map: dict) -> str:
    path = os.path.join(out_dir, "02_taxonomy_colorstrip.txt")

    # Build ordered legend (unique genera in appearance order)
    seen, leg_colors, leg_labels = [], [], []
    for mag in mag_order:
        g = genera.get(mag, "Unknown")
        if g not in seen:
            seen.append(g)
            leg_colors.append(GENUS_COLORS.get(g, "#cccccc"))
            leg_labels.append(g)
    for ref_id in SELECTED_REF_IDS:
        g = REF_GENUS.get(ref_id, "Unknown")
        if g not in seen:
            seen.append(g)
            leg_colors.append(GENUS_COLORS.get(g, "#cccccc"))
            leg_labels.append(g)

    with _open_out(path) as f:
        f.write("DATASET_COLORSTRIP\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tGénero (GTDB)\n")
        f.write("COLOR\t#e76f51\n")
        f.write("STRIP_WIDTH\t30\n")
        f.write("MARGIN\t5\n")
        f.write("BORDER_WIDTH\t1\n")
        f.write("BORDER_COLOR\t#ffffff\n")
        f.write("SHOW_STRIP_LABELS\t0\n")
        f.write("LEGEND_TITLE\tGénero (GTDB)\n")
        f.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(seen)) + "\n")
        f.write("LEGEND_COLORS\t" + "\t".join(leg_colors) + "\n")
        f.write("LEGEND_LABELS\t" + "\t".join(leg_labels) + "\n")
        f.write("DATA\n")
        for mag in mag_order:
            g = genera.get(mag, "Unknown")
            color = GENUS_COLORS.get(g, "#cccccc")
            f.write(f"{translate(mag, leaf_map)}\t{color}\t{g}\n")
        for ref_id in SELECTED_REF_IDS:
            g = REF_GENUS.get(ref_id, "Unknown")
            color = GENUS_COLORS.get(g, "#cccccc")
            f.write(f"{translate(ref_id, leaf_map)}\t{color}\t{g}\n")
    return path


def write_03_completeness_bar(out_dir: str, mag_order: list,
                               comp: dict, leaf_map: dict) -> str:
    path = os.path.join(out_dir, "03_completeness_bar.txt")
    with _open_out(path) as f:
        f.write("DATASET_SIMPLEBAR\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tCompletitud CheckM2 (%)\n")
        f.write("COLOR\t#2a9d8f\n")
        f.write("WIDTH\t150\n")
        f.write("MARGIN\t10\n")
        f.write("DATASET_SCALE\t0,20,40,60,80,100\n")
        f.write("HEIGHT_FACTOR\t0.8\n")
        f.write("BAR_SHIFT\t0\n")
        f.write("LEGEND_TITLE\tCompletitud (MIMAG)\n")
        f.write("LEGEND_SHAPES\t1\t1\t1\n")
        f.write("LEGEND_COLORS\t#2a9d8f\t#a8dadc\t#f4a261\n")
        f.write("LEGEND_LABELS\t≥90 % (alta)\t50–89 % (media)\t<50 % (baja)\n")
        f.write("DATA\n")
        for mag in mag_order:
            val = comp.get(mag)
            if val is None:
                continue
            color = "#2a9d8f" if val >= 90 else ("#a8dadc" if val >= 50 else "#f4a261")
            f.write(f"{translate(mag, leaf_map)}\t{val:.2f}\t{color}\n")
    return path


def write_04_contamination_bar(out_dir: str, mag_order: list,
                                cont: dict, leaf_map: dict) -> str:
    path = os.path.join(out_dir, "04_contamination_bar.txt")
    with _open_out(path) as f:
        f.write("DATASET_SIMPLEBAR\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tContaminación CheckM2 (%)\n")
        f.write("COLOR\t#e63946\n")
        f.write("WIDTH\t80\n")
        f.write("MARGIN\t5\n")
        f.write("DATASET_SCALE\t0,2,4,6,8,10\n")
        f.write("HEIGHT_FACTOR\t0.8\n")
        f.write("LEGEND_TITLE\tContaminación (MIMAG)\n")
        f.write("LEGEND_SHAPES\t1\t1\t1\n")
        f.write("LEGEND_COLORS\t#2a9d8f\t#f4a261\t#e63946\n")
        f.write("LEGEND_LABELS\t≤5 % (alta)\t5–10 % (media)\t>10 % (baja)\n")
        f.write("DATA\n")
        for mag in mag_order:
            val = cont.get(mag)
            if val is None:
                continue
            color = "#2a9d8f" if val <= 5 else ("#f4a261" if val <= 10 else "#e63946")
            f.write(f"{translate(mag, leaf_map)}\t{val:.2f}\t{color}\n")
    return path


def write_05_ani_heatmap(out_dir: str, mag_order: list,
                          ani: dict, leaf_map: dict) -> str:
    path = os.path.join(out_dir, "05_ani_heatmap.txt")
    with _open_out(path) as f:
        f.write("DATASET_HEATMAP\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tANI a referencia más cercana (%)\n")
        f.write("COLOR\t#1d3557\n")
        f.write("FIELD_LABELS\tANI (%)\n")
        f.write("COLOR_MIN\t#f1faee\n")
        f.write("COLOR_MAX\t#1d3557\n")
        f.write("USE_MID_COLOR\t1\n")
        f.write("COLOR_MID\t#a8dadc\n")
        f.write("MIN_VALUE\t80\n")
        f.write("MID_VALUE\t90\n")
        f.write("MAX_VALUE\t100\n")
        f.write("COLOR_NAN\t#dee2e6\n")
        f.write("STRIP_WIDTH\t30\n")
        f.write("MARGIN\t5\n")
        f.write("BORDER_WIDTH\t1\n")
        f.write("BORDER_COLOR\t#ffffff\n")
        f.write("DATA\n")
        for mag in mag_order:
            val = ani.get(mag, -1.0)
            if val < 0:
                f.write(f"{translate(mag, leaf_map)}\t\n")
            else:
                f.write(f"{translate(mag, leaf_map)}\t{val:.2f}\n")
    return path


def write_06_label_styles(out_dir: str, mag_order: list, leaf_map: dict) -> str:
    path = os.path.join(out_dir, "06_label_styles.txt")
    with _open_out(path) as f:
        f.write("DATASET_STYLE\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tLabel styles\n")
        f.write("COLOR\t#000000\n")
        f.write("DATA\n")
        for mag in mag_order:
            f.write(f"{translate(mag, leaf_map)}\tlabel\tnode\t#e63946\t1.1\tbold\t\n")
        for ref_id, _, genus, _ in SELECTED_REFERENCES:
            color = "#adb5bd" if "Outgroup" in genus else "#457b9d"
            f.write(f"{translate(ref_id, leaf_map)}\tlabel\tnode\t{color}\t1.0\titalic\t\n")
    return path


def write_07_dram2_binary(out_dir: str, mag_order: list,
                           product_data: dict, leaf_map: dict) -> str:
    path = os.path.join(out_dir, "07_dram2_binary.txt")

    # Determine active fields (have at least one MAG with value=1)
    n_total = len(DRAM2_FIELDS)
    active_idx = []
    for fi in range(n_total):
        _, cols, thresh = DRAM2_FIELDS[fi]
        for mag in mag_order:
            if mag in product_data and eval_field(product_data[mag], cols, thresh):
                active_idx.append(fi)
                break

    active_fields  = [DRAM2_FIELDS[i] for i in active_idx]
    active_colors  = [DRAM2_FIELD_COLORS[i] for i in active_idx]
    field_labels   = [f[0] for f in active_fields]
    n_active       = len(active_fields)

    results = {}
    for mag in mag_order:
        if mag not in product_data:
            log.warning("DRAM2: %s not found in product.tsv", mag)
            results[mag] = [0] * n_active
            continue
        row = product_data[mag]
        results[mag] = [eval_field(row, cols, thresh)
                        for _, cols, thresh in active_fields]

    with _open_out(path) as f:
        f.write("DATASET_BINARY\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tPotencial metabólico (DRAM2)\n")
        f.write("COLOR\t#264653\n")
        f.write("FIELD_COLORS\t" + "\t".join(active_colors) + "\n")
        f.write("FIELD_LABELS\t" + "\t".join(field_labels) + "\n")
        f.write("FIELD_SHAPES\t" + "\t".join(["1"] * n_active) + "\n")
        f.write("WIDTH\t350\n")
        f.write("MARGIN\t10\n")
        f.write("LEGEND_TITLE\tPotencial metabólico (DRAM2)\n")
        f.write("LEGEND_SHAPES\t1\t1\t1\t1\n")
        f.write("LEGEND_COLORS\t#2a9d8f\t#e9c46a\t#f4a261\t#264653\n")
        f.write("LEGEND_LABELS\t"
                "Fijación de carbono\t"
                "Metabolismo del nitrógeno\t"
                "Metabolismo del azufre\t"
                "Conservación de energía\n")
        f.write("DATA\n")
        for mag in mag_order:
            vals = results.get(mag, [0] * n_active)
            f.write(translate(mag, leaf_map) + "\t"
                    + "\t".join(str(v) for v in vals) + "\n")
    return path


def write_08_clade_colors(out_dir: str, mag_order: list,
                           genera: dict, leaf_map: dict) -> str:
    """
    TREE_COLORS file: one MRCA range per major clade.
    Each range is defined by the first and last MAG (or ref) from that clade.
    """
    path = os.path.join(out_dir, "08_branch_support_colors.txt")

    # Group MAGs by their major taxonomic lineage
    CLADE_ORDER = [
        # (clade_label, set_of_genera, background_color)
        ("Thalassarchaeaceae\n(Poseidoniales)",
         {"Thalassarchaeum", "MGIIb-Q1", "MGIIb-O1", "JAKUTA01", "UBA36"},
         "#fde8e0"),
        ("MGIII\n(CG-Epi1)",
         {"UBA8886"},
         "#d4edda"),
        ("Nitrososphaerales\n(Thermoproteota)",
         {"UBA57", "UBA60", "Nitrosopelagicus"},
         "#d6eaf8"),
        ("Undinarchaeota",
         {"Undinarchaeum"},
         "#fff3cd"),
    ]

    # For each clade, find two representative leaves (first MAG + last MAG)
    clade_ranges = []
    for label, genus_set, color in CLADE_ORDER:
        members = [
            translate(m, leaf_map)
            for m in mag_order
            if genera.get(m, "") in genus_set
        ]
        if len(members) >= 2:
            clade_ranges.append((members[0], members[-1], color, label))
        elif len(members) == 1:
            clade_ranges.append((members[0], members[0], color, label))

    with _open_out(path) as f:
        f.write("TREE_COLORS\n")
        f.write("SEPARATOR TAB\n")
        f.write("# Resaltado de clados principales mediante rangos MRCA.\n")
        f.write("# Par de hojas delimita el MRCA del clado deseado.\n")
        f.write("# Ajustar hoja1|hoja2 si la topología del árbol re-inferido\n")
        f.write("# modifica la posición de algún MAG respecto al esperado.\n")
        f.write("DATA\n")
        for leaf1, leaf2, color, label in clade_ranges:
            safe_label = label.replace("\n", " ")
            f.write(f"{leaf1}|{leaf2}\trange\t{color}\t{safe_label}\n")
    return path


def write_09_station_colorstrip(out_dir: str, mag_order: list,
                                 stations: dict, leaf_map: dict) -> str:
    path = os.path.join(out_dir, "09_station_colorstrip.txt")
    present = list(dict.fromkeys(stations.get(m, m.split(".")[0]) for m in mag_order))
    leg_colors = [STATION_COLORS.get(s, "#999999") for s in present]

    with _open_out(path) as f:
        f.write("DATASET_COLORSTRIP\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tEstación de muestreo\n")
        f.write("COLOR\t#e63946\n")
        f.write("STRIP_WIDTH\t20\n")
        f.write("MARGIN\t5\n")
        f.write("BORDER_WIDTH\t1\n")
        f.write("BORDER_COLOR\t#ffffff\n")
        f.write("SHOW_STRIP_LABELS\t0\n")
        f.write("LEGEND_TITLE\tEstación de muestreo\n")
        f.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(present)) + "\n")
        f.write("LEGEND_COLORS\t" + "\t".join(leg_colors) + "\n")
        f.write("LEGEND_LABELS\t" + "\t".join(present) + "\n")
        f.write("DATA\n")
        for mag in mag_order:
            st = stations.get(mag, mag.split(".")[0])
            color = STATION_COLORS.get(st, "#999999")
            f.write(f"{translate(mag, leaf_map)}\t{color}\t{st}\n")
    return path


def write_10_depth_colorstrip(out_dir: str, mag_order: list,
                               depths: dict, leaf_map: dict) -> str:
    path = os.path.join(out_dir, "10_depth_colorstrip.txt")
    present = list(dict.fromkeys(depths.get(m, "") for m in mag_order if depths.get(m)))
    leg_colors = [DEPTH_COLORS.get(d, ("#999999", d))[0] for d in present]
    leg_labels = [DEPTH_COLORS.get(d, ("#999999", d))[1] for d in present]

    with _open_out(path) as f:
        f.write("DATASET_COLORSTRIP\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tProfundidad\n")
        f.write("COLOR\t#1d3557\n")
        f.write("STRIP_WIDTH\t20\n")
        f.write("MARGIN\t5\n")
        f.write("BORDER_WIDTH\t1\n")
        f.write("BORDER_COLOR\t#ffffff\n")
        f.write("SHOW_STRIP_LABELS\t0\n")
        f.write("LEGEND_TITLE\tProfundidad\n")
        f.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(present)) + "\n")
        f.write("LEGEND_COLORS\t" + "\t".join(leg_colors) + "\n")
        f.write("LEGEND_LABELS\t" + "\t".join(leg_labels) + "\n")
        f.write("DATA\n")
        for mag in mag_order:
            d = depths.get(mag, "")
            color, label = DEPTH_COLORS.get(d, ("#999999", d))
            f.write(f"{translate(mag, leaf_map)}\t{color}\t{label}\n")
    return path

# VALIDATION


def validate_itol_coverage(itol_dir: str,
                            mag_leaves: set,
                            all_leaves: set) -> None:
    """
    Validates iTOL file coverage.

    File-type expectations:
      • 01, 02, 06       to must cover ALL leaves (MAGs + refs)
      • 03–05, 07, 09–10 to must cover MAG-only leaves (refs have no data)
      • 08               to TREE_COLORS range file; only lists clade-boundary
                           leaf pairs — skip completeness check, just report counts
    """
    log.info("iTOL coverage validation")
    for fname in sorted(os.listdir(itol_dir)):
        if not fname.endswith(".txt"):
            continue

        # Select expected leaf set for this file
        prefix = fname.split("_")[0]
        if prefix in ("01", "02"):
            expected = all_leaves
        elif prefix == "08":
            expected = None   # range file, no full-coverage expectation
        else:
            expected = mag_leaves

        referenced = set()
        in_data = False
        with open(os.path.join(itol_dir, fname), encoding="utf-8") as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line == "DATA":
                    in_data = True
                    continue
                if in_data and line and not line.startswith("#"):
                    leaf = line.split("\t")[0]
                    for part in leaf.split("|"):
                        referenced.add(part.strip())

        if expected is None:
            # TREE_COLORS: just report how many boundary leaves are listed
            log.info("  %s: %d range entries (clade-boundary file, no full coverage expected)",
                     fname, len(referenced))
            continue

        covered = referenced & expected
        missing = expected - referenced
        extra   = referenced - expected
        status  = "OK" if not missing and not extra else "WARN"
        log.info("  %s: %d/%d covered [%s]",
                 fname, len(covered), len(expected), status)
        for m in sorted(missing):
            log.warning("    MISSING in file: %s", m)
        for e in sorted(extra):
            log.warning("    EXTRA in file (not in expected leaf set): %s", e)

# MAIN


def main() -> None:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "--project-dir", "-p",
        default=DEFAULT_PROJECT_DIR,
        help="Root of the msc_archaea project (default: auto-detected)",
    )
    parser.add_argument(
        "--threads", "-t", type=int, default=12,
        help="CPU threads for FastTree (default: 12)",
    )
    parser.add_argument(
        "--skip-fasttree", action="store_true",
        help="Skip FastTree re-inference (only prune + iTOL)",
    )
    parser.add_argument(
        "--skip-prune", action="store_true",
        help="Skip DendroPy pruning of the decorated tree",
    )
    args = parser.parse_args()

    P = build_paths(args.project_dir)
    os.makedirs(P["output_dir"], exist_ok=True)
    os.makedirs(P["itol_dir"],   exist_ok=True)

    # ---- Validate inputs ----
    required = ["msa_full_gz", "taxonomy_tsv", "metadata_tsv",
                "quality_tsv", "product_tsv"]
    missing_files = [P[k] for k in required if not os.path.isfile(P[k])]
    if missing_files:
        for f in missing_files:
            log.error("Input file not found: %s", f)
        sys.exit(1)

    # ---- Write reference selection report ----
    write_reference_selection_table(
        os.path.join(P["output_dir"], "reference_selection.tsv"))

    # ---- Load metadata ----
    log.info("Loading metadata …")
    mag_order, stations, depths, display, genera = load_metadata(P["metadata_tsv"])
    target_set = set(mag_order)
    log.info("  %d MAGs loaded from metadata", len(mag_order))

    # ---- Load quality / ANI / DRAM2 ----
    log.info("Loading quality metrics …")
    comp, cont = load_quality(P["quality_tsv"], target_set)
    # Called for the coverage warnings it logs; the table itself goes unused
    # here, since the ANI annotation file is written from product instead.
    load_ani(P["taxonomy_tsv"], target_set)
    product     = load_product(P["product_tsv"], mag_order, target_set)

    # ---- Build leaf display map ----
    leaf_map = build_leaf_display_map(mag_order, display)
    # Set of display names used in iTOL files (for validation)
    itol_mag_leaves = {leaf_map[m] for m in mag_order}
    itol_ref_leaves = {leaf_map[r] for r in SELECTED_REF_IDS}

    # ---- Define keep_ids for MSA extraction & tree pruning ----
    # MAG IDs as they appear in the GTDB-Tk MSA headers (dot-notation)
    # Ref IDs as GB_/RS_ prefixed GTDB leaf IDs
    keep_ids = set(mag_order) | SELECTED_REF_IDS

    # ---- Extract curated MSA ----
    curated_msa = os.path.join(P["output_dir"], "curated_54mags.msa.fasta")
    log.info("Extracting curated MSA (%d taxa) …", len(keep_ids))
    n_written = extract_msa(P["msa_full_gz"], keep_ids, curated_msa)
    if n_written < len(keep_ids):
        log.warning(
            "Only %d/%d requested sequences found in MSA. "
            "Some reference IDs may not be in the alignment.",
            n_written, len(keep_ids),
        )

    # Mapping used to rename tree leaves: covers MAGs (dot-notation to DisplayName)
    # AND refs (GB_/RS_ GTDB prefix to REF_DISPLAY label).  Both tree files and
    # all iTOL annotation files use these final labels, so they must match exactly.
    full_rename_map = {**{mag: display[mag].replace(" ", "_") for mag in mag_order}, **REF_DISPLAY}

    # ---- Re-infer tree with FastTree ----
    fasttree_nwk = os.path.join(P["output_dir"], "curated_54mags.fasttree.nwk")
    fasttree_log = os.path.join(P["output_dir"], "curated_54mags.fasttree.log")
    if not args.skip_fasttree:
        run_fasttree(curated_msa, fasttree_nwk, fasttree_log, args.threads)
        rename_tree_leaves(fasttree_nwk, fasttree_nwk, full_rename_map)
    else:
        log.info("--skip-fasttree: skipping tree re-inference.")

    # ---- Prune GTDB-topology decorated tree ----
    pruned_nwk = os.path.join(P["output_dir"], "curated_54mags.gtdb_pruned.nwk")
    if not args.skip_prune:
        decorated = os.path.join(P["de_novo_tree"],
                                 "54_archaea.ar53.decorated.tree")
        prune_decorated_tree(decorated, keep_ids, pruned_nwk)
        rename_tree_leaves(pruned_nwk, pruned_nwk, full_rename_map)
    else:
        log.info("--skip-prune: skipping DendroPy pruning.")

    # ---- Generate iTOL annotation files ----
    # Remove any stale .txt files from previous runs so the validator
    # only sees files generated by the current configuration.
    for _old in os.listdir(P["itol_dir"]):
        if _old.endswith(".txt"):
            os.remove(os.path.join(P["itol_dir"], _old))
    log.info("Generating iTOL annotation files …")
    writers = [
        write_01_source_colorstrip(P["itol_dir"], mag_order, leaf_map),
        write_02_taxonomy_colorstrip(P["itol_dir"], mag_order, genera, leaf_map),
        # write_03_completeness_bar(P["itol_dir"], mag_order, comp, leaf_map),
        # write_04_contamination_bar(P["itol_dir"], mag_order, cont, leaf_map),
        # write_05_ani_heatmap(P["itol_dir"], mag_order, ani, leaf_map),
        # write_06_label_styles(P["itol_dir"], mag_order, leaf_map),
        write_07_dram2_binary(P["itol_dir"], mag_order, product, leaf_map),
        write_08_clade_colors(P["itol_dir"], mag_order, genera, leaf_map),
        write_09_station_colorstrip(P["itol_dir"], mag_order, stations, leaf_map),
        write_10_depth_colorstrip(P["itol_dir"], mag_order, depths, leaf_map),
    ]
    for w in writers:
        log.info("  ✔ %s", os.path.relpath(w, args.project_dir))

    # ---- Validate iTOL coverage ----
    # 01/02 must cover all leaves; 07/09-10 MAGs only; 08 is range-only
    # 03-06 are commented out (not generated)
    all_itol_leaves = itol_mag_leaves | itol_ref_leaves
    validate_itol_coverage(P["itol_dir"], itol_mag_leaves, all_itol_leaves)

    # ---- Summary ----
    log.info("══════════════════════════════════════════════════════════")
    log.info("Done.")
    log.info("  MAGs in tree         : %d", len(mag_order))
    log.info("  Reference genomes    : %d", len(SELECTED_REF_IDS))
    log.info("  Total leaves         : %d", len(mag_order) + len(SELECTED_REF_IDS))
    log.info("  Output directory     : %s", P["output_dir"])
    if not args.skip_fasttree and os.path.isfile(fasttree_nwk):
        log.info("  Re-inferred tree     : %s", fasttree_nwk)
    if not args.skip_prune and os.path.isfile(pruned_nwk):
        log.info("  GTDB-pruned tree     : %s", pruned_nwk)
    log.info("  iTOL annotations     : %s", P["itol_dir"])
    log.info("══════════════════════════════════════════════════════════")

if __name__ == "__main__":
    main()
