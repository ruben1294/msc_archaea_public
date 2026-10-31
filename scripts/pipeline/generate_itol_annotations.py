#!/usr/bin/env python3
"""Write the iTOL annotation files for the 10-MAG tree from pipeline step 18.

The tree covers Poseidoniales (Thalassarchaeaceae: MGIIb-Q1 and UBA36) and
Nitrosopumilaceae (Nitrosopelagicus). iTOL takes one file per annotation track,
so each of these lands in OUTPUT_DIR and gets dragged onto the tree separately:

    01_source_colorstrip.txt      DATASET_COLORSTRIP   MAG or reference
    02_taxonomy_colorstrip.txt    DATASET_COLORSTRIP   taxonomic group
    03_completeness_bar.txt       DATASET_SIMPLEBAR    CheckM2 completeness
    04_contamination_bar.txt      DATASET_SIMPLEBAR    CheckM2 contamination
    06_label_styles.txt           DATASET_STYLE        label formatting
    07_dram2_binary.txt           DATASET_BINARY       DRAM2 metabolic potential
    08_branch_support_colors.txt  TREE_COLORS          clade highlighting by MRCA
    09_station_colorstrip.txt     DATASET_COLORSTRIP   sampling station
    10_depth_colorstrip.txt       DATASET_COLORSTRIP   sampling depth

Leaf labels have to match the tree, so the renaming tables below are written out
by hand: they follow the mapping step 18 used to relabel the tree.
"""

import os
import sys
import csv
import re

# CONFIGURATION
# Project root: PROJECT_DIR or MSC_ARCHAEA_DIR when set, otherwise this clone.
# The metadata table is part of the repository, so it is always read from there.
REPO_ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))
PROJECT_DIR = (os.environ.get("PROJECT_DIR")
               or os.environ.get("MSC_ARCHAEA_DIR")
               or REPO_ROOT)
RESULTS_DIR = os.environ.get("RESULTS_DIR") or os.path.join(PROJECT_DIR, "results")

TAXONOMY_TSV = os.path.join(RESULTS_DIR, "03_taxonomy_assigned", "classify", "tax_.ar53.summary.tsv")
PRODUCT_TSV = os.path.join(RESULTS_DIR, "10_MAG_distill", "VISUALIZE", "product.tsv")
METADATA_TSV = os.path.join(REPO_ROOT, "data", "MAGs_metadata.tsv")
QUALITY_TSV = os.path.join(RESULTS_DIR, "01_mag_quality", "quality_report.tsv")
OUTPUT_DIR = os.path.join(RESULTS_DIR, "itol_annotations")

# Tree leaves, by their internal identifiers before relabelling
TARGET_MAGS = [
    'A4.1000.megahit.bin.28',
    'A4.BTM.megahit.bin.15',
    'A6.BTM.megahit.bin.3',
    'PO1.BTM.megahit.bin.42',
    'Y6.1000.megahit.bin.23',
    'A4.BTM.megahit.bin.25',
    'PO1.BTM.megahit.bin.7',
    'A4.BTM.megahit.bin.33',
    'A4.BTM.megahit.bin.58',
    'A6.BTM.megahit.bin.9',
]

# Reference names as they appear in 10mags_pmsf.renamed.treefile
REFERENCE_LEAVES = [
    'REF_GCA_012964335.1',  # Thalassarchaeaceae sp.
    'REF_GCA_002505495.1',  # Thalassarchaeaceae sp.
    'REF_GCA_963580235.1',  # Thalassarchaeaceae sp.
    'REF_GCA_036153975.1',  # UBA36 sp.
    'REF_GCA_036336595.1',  # UBA36 sp.
    'REF_GCA_013390375.1',  # Ca. Nitrosopelagicus sp.
    'REF_GCF_000018465.1',  # Nitrosopumilus maritimus SCM1
    'REF_GCA_000012285.1',  # Sulfolobus acidocaldarius
    'REF_GCA_000259595.1',  # Ca. Nitrososphaera viennensis
    'REF_GCA_000091665.1',  # Methanocaldococcus jannaschii
    'REF_GCF_000008545.1',  # Thermotoga maritima (outgroup)
    'REF_GCF_000008625.1',  # Aquifex aeolicus (outgroup)
    'REF_GCA_000402075.1',  # not renamed (as in the treefile)
    'REF_GCA_000008545.1',  # not renamed (as in the treefile)
    'REF_GCA_000008625.1',  # not renamed (as in the treefile)
    'REF_GCA_000009345.1',  # not renamed (as in the treefile)
    'REF_GCA_000011805.1',  # not renamed (as in the treefile)
    'REF_GCA_900232055.1',  # not renamed (as in the treefile)
]

ALL_LEAVES = TARGET_MAGS + REFERENCE_LEAVES
BACTERIAL_OUTGROUPS = {'REF_GCF_000008545.1', 'REF_GCA_000008545.1'}

# Mapping to the relabelled names, as they appear in 10mags_pmsf.renamed.treefile

# Internal name to relabelled name, for the MAGs
MAG_DISPLAY_NAMES = {
    'A4.1000.megahit.bin.28':  'MGIIb-Q1_A4_1000_28',
    'A4.BTM.megahit.bin.15':   'MGIIb-Q1_A4_BTM_15',
    'A6.BTM.megahit.bin.3':    'MGIIb-Q1_A6_BTM_3',
    'PO1.BTM.megahit.bin.42':  'MGIIb-Q1_PO1_BTM_42',
    'Y6.1000.megahit.bin.23':  'MGIIb-Q1_Y6_1000_23',
    'A4.BTM.megahit.bin.25':   'Nitrosopelagicus_A4_BTM_25',
    'PO1.BTM.megahit.bin.7':   'Nitrosopelagicus_PO1_BTM_7',
    'A4.BTM.megahit.bin.33':   'Nitrosopelagicus_A4_BTM_33',
    'A4.BTM.megahit.bin.58':   'UBA36_A4_BTM_58',
    'A6.BTM.megahit.bin.9':    'UBA36_A6_BTM_9',
}

# The same mapping for the reference genomes
REF_DISPLAY_NAMES = {
    # Thalassarchaeaceae / MGIIb-Q1
    'REF_GCA_012964335.1': 'Thalassarchaeaceae_sp._GCA_012964335.1',
    'REF_GCA_002505495.1': 'Thalassarchaeaceae_sp._GCA_002505495.1',
    'REF_GCA_963580235.1': 'Thalassarchaeaceae_sp._GCA_963580235.1',
    # UBA36
    'REF_GCA_036153975.1': 'UBA36_sp._GCA_036153975.1',
    'REF_GCA_036336595.1': 'UBA36_sp._GCA_036336595.1',
    # Nitrosopelagicus / Nitrosopumilus / Nitrososphaera
    'REF_GCA_013390375.1': 'Ca._Nitrosopelagicus_sp._GCA_013390375.1',
    'REF_GCF_000018465.1': 'Nitrosopumilus_maritimus_SCM1',
    'REF_GCA_000259595.1': 'Ca._Nitrososphaera_viennensis_EN76',
    # Sulfolobus
    'REF_GCA_000012285.1': 'Sulfolobus_acidocaldarius_DSM_639',
    # Methanocaldococcus
    'REF_GCA_000091665.1': 'Methanocaldococcus_jannaschii_DSM_2661',
    # Bacterial outgroup
    'REF_GCF_000008545.1': 'Thermotoga_maritima_MSB8',
    'REF_GCA_000008545.1': 'REF_GCA_000008545.1',
    'REF_GCF_000008625.1': 'Aquifex_aeolicus_VF5',
    'REF_GCA_000008625.1': 'REF_GCA_000008625.1',
    # Anything else keeps its name
    'REF_GCA_000402075.1': 'REF_GCA_000402075.1',
    'REF_GCA_000009345.1': 'REF_GCA_000009345.1',
    'REF_GCA_000011805.1': 'REF_GCA_000011805.1',
    'REF_GCA_900232055.1': 'REF_GCA_900232055.1',
}

# main() fills this in with:
#   MAGs: NameWithoutExtension to DisplayName, spaces become underscores
#   refs: REF_DISPLAY_NAMES
LEAF_NAME_MAP = {}


def translate_leaf(name: str) -> str:
    """Returns the renamed label that the treefile and iTOL use."""
    return LEAF_NAME_MAP.get(name, name)

# HARD-CODED VALUES (CheckM2, checked)
ANI_EXPECTED = {
    'A6.BTM.megahit.bin.3': 93.23,
    'A4.BTM.megahit.bin.15': 93.02,
    'PO1.BTM.megahit.bin.42': 93.12,
    'Y6.1000.megahit.bin.23': 94.89,
    'A4.1000.megahit.bin.28': 93.90,
    'A4.BTM.megahit.bin.33': 84.59,
    'A4.BTM.megahit.bin.25': -1,
    'PO1.BTM.megahit.bin.7': -1,
    'A4.BTM.megahit.bin.58': 93.17,
    'A6.BTM.megahit.bin.9': 93.12,
}

# Expected DRAM2 values, checked on load. The 14 fields are listed in
# DRAM2_FIELDS below.
DRAM2_EXPECTED = {
    'A4.1000.megahit.bin.28': [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    'A4.BTM.megahit.bin.15': [1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0],
    'A6.BTM.megahit.bin.3':  [1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0],
    'PO1.BTM.megahit.bin.42': [1, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0],
    'Y6.1000.megahit.bin.23': [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    'A4.BTM.megahit.bin.25': [0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0],
    'PO1.BTM.megahit.bin.7': [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0],
    'A4.BTM.megahit.bin.33': [0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 1, 0],
    'A4.BTM.megahit.bin.58': [0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
    'A6.BTM.megahit.bin.9':  [1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
}

# Taxonomic assignment per leaf. Written out rather than derived, because it
# follows the topology of this particular tree.
TAXONOMY_MAP = {
    # MGIIb-Q1 (Poseidoniales / Thalassarchaeaceae)
    'A4.1000.megahit.bin.28': ('#e76f51', 'MGIIb-Q1'),
    'A4.BTM.megahit.bin.15':  ('#e76f51', 'MGIIb-Q1'),
    'A6.BTM.megahit.bin.3':   ('#e76f51', 'MGIIb-Q1'),
    'PO1.BTM.megahit.bin.42': ('#e76f51', 'MGIIb-Q1'),
    'Y6.1000.megahit.bin.23': ('#e76f51', 'MGIIb-Q1'),
    'REF_GCA_012964335.1':    ('#e76f51', 'MGIIb-Q1'),
    'REF_GCA_002505495.1':    ('#e76f51', 'MGIIb-Q1'),
    'REF_GCA_963580235.1':    ('#e76f51', 'MGIIb-Q1'),

    # UBA36 (Poseidoniales / Thalassarchaeaceae)
    'A4.BTM.megahit.bin.58':  ('#2a9d8f', 'UBA36'),
    'A6.BTM.megahit.bin.9':   ('#2a9d8f', 'UBA36'),
    'REF_GCA_036153975.1':    ('#2a9d8f', 'UBA36'),
    'REF_GCA_036336595.1':    ('#2a9d8f', 'UBA36'),

    # Nitrososphaerales / Nitrosopumilaceae
    'A4.BTM.megahit.bin.25':  ('#264653', 'Nitrosopelagicus'),
    'PO1.BTM.megahit.bin.7':  ('#264653', 'Nitrosopelagicus'),
    'A4.BTM.megahit.bin.33':  ('#264653', 'Nitrosopelagicus'),
    'REF_GCA_013390375.1':    ('#264653', 'Nitrosopelagicus'),
    'REF_GCF_000812185.2':    ('#264653', 'Nitrosopelagicus'),  # CN25
    'REF_GCF_000018465.1':    ('#264653', 'Nitrosopumilus'),    # SCM1
    'REF_GCA_000259595.1':    ('#264653', 'Nitrososphaera'),

    # Thermoproteota (Crenarchaeota s.l. / Sulfolobales)
    'REF_GCA_000012285.1':    ('#f4a261', 'Sulfolobus'),

    # Methanobacteriota
    'REF_GCA_000091665.1':    ('#6a4c93', 'Methanocaldococcus'),

    # Bacteria (outgroup)
    'REF_GCF_000008545.1':    ('#adb5bd', 'Bacteria'),
    'REF_GCF_000008625.1':    ('#adb5bd', 'Bacteria'),
}

# DRAM2 fields for 07_dram2_binary.txt, labelled in Spanish for the figure
DRAM2_FIELDS = [
    ('Fijación de carbono (CBB)',
     ['Reductive pentose phosphate cycle (Calvin cycle)'], 0.5),
    ('Fijación de carbono (rTCA)',
     ['Reductive citrate cycle (Arnon-Buchanan cycle)'], 0.5),
    ('Fijación de carbono (WL)',
     ['Reductive acetyl-CoA pathway (Wood-Ljungdahl pathway)'], 0.5),
    ('Metanogénesis',
     ['Methanogenesis, CO2 => methane'], 0.5),
    ('Fijación de nitrógeno',
     ['Nitrogen: nitrogen => ammonia'], None),
    ('Oxidación de amoniaco',
     ['Nitrogen: Bacterial/Archaeal ammonia oxidation',
      'Nitrogen: ammonia => nitrite'], None),
    ('Reducción disimatoria de nitrato',
     ['Nitrogen: nitrate => nitrite',
      'Nitrogen: Dissimilatory nitrite reduction to ammonia (DNRA)'], None),
    ('Reducción asimilatoria de nitrato',
     ['Nitrogen: nitrate => nitrite'], None),
    ('Desnitrificación',
     ['Nitrogen: nitrite => nitric oxide',
      'Nitrogen: nitric oxide => nitrous oxide',
      'Nitrogen: nitrous oxide => nitrogen'], None),
    ('Oxidación de azufre',
     ['Sulfur: Thiosulfate oxidation by SOX complex, thiosulfate => sulfate'], None),
    ('Reducción de sulfato',
     ['Sulfur: dissimilatory sulfate reduction (and oxidation) sulfate => sulfide'], None),
    ('Respiración aeróbica',
     ['Complex IV L: Cytochrome c oxidase',
      'Complex IV L: Cytochrome c oxidase, prokaryotes',
      'Complex IV H: Cytochrome c oxidase, cbb3-type',
      'Complex IV L: Cytochrome o ubiquinol oxidase',
      'Complex IV H: Cytochrome bd ubiquinol oxidase'], 0.5),
    ('Fermentación',
     ['SCFA/alcohol conversions: Alcohol production',
      'SCFA/alcohol conversions: Butyrate, pt 1',
      'SCFA/alcohol conversions: Butyrate, pt 2',
      'SCFA/alcohol conversions: lactate D',
      'SCFA/alcohol conversions: lactate L',
      'SCFA/alcohol conversions: acetate, pt 1',
      'SCFA/alcohol conversions: acetate, pt 2'], None),
    ('Fototrofía',
     ['Photo: Photosystem I',
      'Photo: Photosystem II'], None),
]

DRAM2_FIELD_COLORS = [
    '#2a9d8f',  # CBB
    '#2a9d8f',  # rTCA
    '#2a9d8f',  # WL
    '#6a4c93',  # Methanogenesis
    '#e9c46a',  # N fixation
    '#e9c46a',  # Ammonia oxidation
    '#e9c46a',  # Dissim. NO3
    '#e9c46a',  # Assim. NO3
    '#e9c46a',  # Denitrification
    '#f4a261',  # Sulfur oxidation
    '#f4a261',  # Sulfate reduction
    '#264653',  # Aerobic respiration
    '#264653',  # Fermentation
    '#e76f51',  # Phototrophy
]

# CLADE RANGES (08_branch_support_colors.txt)
CLADE_RANGES = [
    ('A4.BTM.megahit.bin.15', 'REF_GCA_012964335.1', '#fde8e0', 'MGIIb-Q1 sublinaje A'),
    ('Y6.1000.megahit.bin.23', 'REF_GCA_963580235.1', '#fde8e0', 'MGIIb-Q1 sublinaje B'),
    ('A4.BTM.megahit.bin.58', 'REF_GCA_036336595.1', '#d4edda', 'UBA36'),
    ('A4.BTM.megahit.bin.25', 'REF_GCA_013390375.1', '#d6eaf8', 'Nitrososphaerales'),
    ('REF_GCF_000008545.1',   'REF_GCF_000008625.1', '#e9ecef', 'Bacteria (grupo externo)'),
]

# STATION AND DEPTH COLOURS
STATION_COLORS = {
    'A4':  '#e63946',
    'A6':  '#2a9d8f',
    'PO1': '#e9c46a',
    'Y6':  '#6a4c93',
}

DEPTH_COLORS = {
    'BTM':  ('#1d3557', 'Fondo'),
    '1000': ('#a8dadc', '1000 m'),
    'OMZ':  ('#e9c46a', 'ZMOR'),
}

# ══════════════════════════════════════════════════════════════════════════════
# PARSERS
# ══════════════════════════════════════════════════════════════════════════════


def parse_quality(path):
    """
    Reads the CheckM2 report (columns Name, Completeness, Contamination).
    Names written before utils/05_rename_mags.sh use underscores
    (A4_BTM_megahit_bin.58) and go through tsv_name_to_dot(); renamed ones
    already have dots and pass through unchanged.
    Returns:
      completeness : dict dot_name to float
      contamination : dict dot_name to float
    """
    completeness, contamination = {}, {}
    target_set = set(TARGET_MAGS)
    with open(path, newline='', encoding='utf-8') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            raw_name = row.get('Name', '').strip()
            dot_name = tsv_name_to_dot(raw_name)
            if dot_name not in target_set:
                continue
            try:
                completeness[dot_name] = float(row['Completeness'])
                contamination[dot_name] = float(row['Contamination'])
            except (KeyError, ValueError):
                print(f" WARNING [quality]: no se pudo parsear {raw_name}", file=sys.stderr)
    missing = target_set - set(completeness)
    if missing:
        for m in sorted(missing):
            print(f" WARNING [quality]: {m} not found en quality report", file=sys.stderr)
    return completeness, contamination


def parse_metadata(path):
    """
    Reads MAGs_metadata.tsv.
    Returns:
      order : NameWithoutExtension values in the order they appear, no header
      stations : dict NameWithoutExtension to Station
      depths   : dict NameWithoutExtension to Depth
      display_names_tree : dict NameWithoutExtension to the DisplayName as written in the tree
    """
    order, stations, depths = [], {}, {}
    display_names_tree = {}
    with open(path, newline='', encoding='utf-8') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            name = row.get('NameWithoutExtension', '').strip()
            if not name:
                continue
            order.append(name)

            # Station and Depth fall back to the filename when the column is absent
            stations[name] = row.get('Station', '').strip() or name.split('.')[0]
            depth_val = row.get('Depth', '').strip()
            if not depth_val and '.' in name:
                depth_val = name.split('.')[1]  # A4.BTM.megahit... to BTM
            depths[name] = depth_val

            # Tree labels cannot contain spaces
            disp = row.get('DisplayName', '').strip()
            if disp:
                display_names_tree[name] = disp.replace(' ', '_')

    return order, stations, depths, display_names_tree


def tsv_name_to_dot(tsv_name):
    """
    Converts a TSV-style name (underscores between the main fields) to the
    tree format, dot separated:
      A4_BTM_megahit_bin.28  to A4.BTM.megahit.bin.28
      A4_1000_megahit_bin.28 to A4.1000.megahit.bin.28
    Only the underscores between STATION, DEPTH and ASSEMBLER are replaced;
    the dot before the BinID stays.
    """
    # Pattern: (station)_(depth)_(assembler)_bin.(id)
    m = re.match(r'^([^_]+)_([^_]+)_([^_]+)_bin\.(.+)$', tsv_name)
    if m:
        return f"{m.group(1)}.{m.group(2)}.{m.group(3)}.bin.{m.group(4)}"
    return tsv_name  # no match: return it unchanged


def parse_ani(path, target_set):
    """
    Reads tax_.ar53.summary.tsv and pulls the most relevant ANI per target MAG.
    Priority: closest_genome_ani, then closest_placement_ani, then -1 (no hit).
    Returns dict: dot_name to float (ANI) or -1.
    """
    ani = {}
    with open(path, newline='', encoding='utf-8') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for row in reader:
            tsv_name = row.get('user_genome', '').strip()
            dot_name = tsv_name_to_dot(tsv_name)
            if dot_name not in target_set:
                continue
            val = -1
            for col in ('closest_genome_ani', 'closest_placement_ani'):
                raw = row.get(col, '').strip()
                if raw and raw not in ('N/A', '', 'None', 'nan'):
                    try:
                        val = float(raw)
                        break
                    except ValueError:
                        pass
            ani[dot_name] = val
    return ani


def parse_product(product_path, metadata_order):
    """
    Reads product.tsv. It has no index column, and its rows align with metadata_order.
    Returns dict: mag_name to dict{column_name to value_string}
    Only MAGs in TARGET_MAGS are included.
    """
    data = {}
    target_set = set(TARGET_MAGS)
    with open(product_path, newline='', encoding='utf-8') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for i, row in enumerate(reader):
            if i >= len(metadata_order):
                break
            mag_name = metadata_order[i]
            if mag_name in target_set:
                data[mag_name] = dict(row)
    return data


def eval_field(mag_row, columns, threshold):
    """
    Whether a metabolic field is present for one MAG.
    Logical OR across the listed columns:
      - threshold=None: a boolean column, True or '1' gives 1, False or '0' gives 0
      threshold as a float treats the column as numeric 0 to 1, and counts it
      as present at or above that value
    Returns 1 or 0.
    """
    for col in columns:
        raw = mag_row.get(col, '').strip()
        if not raw:
            continue
        if threshold is None:
            if raw.lower() in ('true', '1'):
                return 1
        else:
            try:
                if float(raw) >= threshold:
                    return 1
            except ValueError:
                pass
    return 0


def get_dram2_active_fields():
    """
    Which DRAM2 fields carry at least one positive across the MAGs.
    Returns the indices of the non-empty fields.
    """
    active = []
    n_fields = len(DRAM2_EXPECTED.get(list(DRAM2_EXPECTED.keys())[0], []))

    for field_idx in range(n_fields):
        has_any = False
        for mag, values in DRAM2_EXPECTED.items():
            if values[field_idx] == 1:
                has_any = True
                break
        if has_any:
            active.append(field_idx)

    return active

# ══════════════════════════════════════════════════════════════════════════════
# iTOL FILE WRITERS
# ══════════════════════════════════════════════════════════════════════════════


def _open(outdir, filename):
    """Opens an output file with Unix newlines."""
    return open(os.path.join(outdir, filename), 'w', newline='\n', encoding='utf-8')


def write_01_source_colorstrip(outdir):
    fname = '01_source_colorstrip.txt'
    with _open(outdir, fname) as f:
        f.write("DATASET_COLORSTRIP\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tFuente\n")
        f.write("COLOR\t#e63946\n")
        f.write("STRIP_WIDTH\t0\n")
        f.write("MARGIN\t0\n")
        f.write("BORDER_WIDTH\t0\n")
        f.write("BORDER_COLOR\t#ffffff\n")
        f.write("SHOW_STRIP_LABELS\t0\n")
        f.write("COLOR_BRANCHES\t1\n")
        f.write("LEGEND_TITLE\tFuente\n")
        f.write("LEGEND_SHAPES\t1\t1\n")
        f.write("LEGEND_COLORS\t#e63946\t#457b9d\n")
        f.write("LEGEND_LABELS\tMAG\tReferencia\n")
        f.write("DATA\n")
        for leaf in TARGET_MAGS:
            f.write(f"{translate_leaf(leaf)}\t#e63946\tMAG\n")
        for leaf in REFERENCE_LEAVES:
            f.write(f"{translate_leaf(leaf)}\t#457b9d\tReference\n")
    return fname


def write_02_taxonomy_colorstrip(outdir):
    fname = '02_taxonomy_colorstrip.txt'
    # One legend entry per distinct colour, in the order they first appear
    seen, legend_colors, legend_labels = [], [], []
    for leaf in ALL_LEAVES:
        if leaf in TAXONOMY_MAP:
            color, label = TAXONOMY_MAP[leaf]
            if label not in seen:
                seen.append(label)
                legend_colors.append(color)
                legend_labels.append(label)

    with _open(outdir, fname) as f:
        f.write("DATASET_COLORSTRIP\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tGrupo taxonómico\n")
        f.write("COLOR\t#e76f51\n")
        f.write("STRIP_WIDTH\t25\n")
        f.write("MARGIN\t5\n")
        f.write("BORDER_WIDTH\t1\n")
        f.write("BORDER_COLOR\t#ffffff\n")
        f.write("SHOW_STRIP_LABELS\t0\n")
        f.write("LEGEND_TITLE\tGrupo taxonómico\n")
        f.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(legend_colors)) + "\n")
        f.write("LEGEND_COLORS\t" + "\t".join(legend_colors) + "\n")
        f.write("LEGEND_LABELS\t" + "\t".join(legend_labels) + "\n")
        f.write("DATA\n")
        for leaf in ALL_LEAVES:
            out_leaf = translate_leaf(leaf)
            if leaf in TAXONOMY_MAP:
                color, label = TAXONOMY_MAP[leaf]
                f.write(f"{out_leaf}\t{color}\t{label}\n")
            else:
                f.write(f"{out_leaf}\t#cccccc\tUnknown\n")
    return fname


def write_03_completeness_bar(outdir, completeness, contamination):
    fname = '03_completeness_bar.txt'
    with _open(outdir, fname) as f:
        f.write("DATASET_SIMPLEBAR\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tCompletitud (%)\n")
        f.write("COLOR\t#2a9d8f\n")
        f.write("WIDTH\t150\n")
        f.write("MARGIN\t10\n")
        f.write("DATASET_SCALE\t0,20,40,60,80,100\n")
        f.write("HEIGHT_FACTOR\t0.8\n")
        f.write("BAR_SHIFT\t0\n")
        f.write("LEGEND_TITLE\tCompletitud\n")
        f.write("LEGEND_SHAPES\t1\t1\t1\n")
        f.write("LEGEND_COLORS\t#2a9d8f\t#a8dadc\t#f4a261\n")
        f.write("LEGEND_LABELS\t>=90% (calidad alta)\t50-89% (calidad media)\t<50% (calidad baja)\n")
        f.write("DATA\n")
        for mag in TARGET_MAGS:
            val = completeness.get(mag)
            if val is None:
                continue
            # MIMAG bands
            # https://www.nature.com/articles/nbt.3893
            if val >= 90:
                color = '#2a9d8f'  # High-quality
            elif val >= 50:
                color = '#a8dadc'  # Medium-quality
            else:
                color = '#f4a261'  # Low-quality
            f.write(f"{translate_leaf(mag)}\t{val}\t{color}\n")
    return fname


def write_04_contamination_bar(outdir, contamination):
    fname = '04_contamination_bar.txt'
    with _open(outdir, fname) as f:
        f.write("DATASET_SIMPLEBAR\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tContaminación (%)\n")
        f.write("COLOR\t#2a9d8f\n")
        f.write("WIDTH\t80\n")
        f.write("MARGIN\t5\n")
        f.write("DATASET_SCALE\t0,2,4,6,8,10\n")
        f.write("HEIGHT_FACTOR\t0.8\n")
        f.write("LEGEND_TITLE\tContaminación\n")
        f.write("LEGEND_SHAPES\t1\t1\t1\n")
        f.write("LEGEND_COLORS\t#2a9d8f\t#f4a261\t#e63946\n")
        f.write("LEGEND_LABELS\t<=5% (calidad alta)\t5-10% (calidad media)\t>10% (calidad baja)\n")
        f.write("DATA\n")
        for mag in TARGET_MAGS:
            val = contamination.get(mag)
            if val is None:
                continue
            # MIMAG bands (https://www.nature.com/articles/nbt.3893)
            # High-quality MAGs: ≤5%, Medium-quality MAGs: ≤10%
            if val <= 5:
                color = '#2a9d8f'  # high quality (green)
            elif val <= 10:
                color = '#f4a261'  # medium quality (orange)
            else:
                color = '#e63946'  # low quality (red)
            f.write(f"{translate_leaf(mag)}\t{val}\t{color}\n")
    return fname


def write_05_ani_heatmap(outdir, ani_data):
    fname = '05_ani_heatmap.txt'
    # Check against the expected values
    for mag, expected in ANI_EXPECTED.items():
        observed = ani_data.get(mag, None)
        if observed is None:
            print(f" WARNING [05]: {mag} not found en taxonomy TSV", file=sys.stderr)
        elif expected != -1 and observed != -1:
            if abs(observed - expected) > 0.5:
                print(
                    f" WARNING [05]: {mag} ANI observado={observed:.2f} esperado={expected:.2f}",
                    file=sys.stderr,
                )

    with _open(outdir, fname) as f:
        f.write("DATASET_HEATMAP\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tANI to closest reference\n")
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
        for mag in TARGET_MAGS:
            val = ani_data.get(mag, -1)
            if val == -1:
                f.write(f"{translate_leaf(mag)}\t\n")  # blank, which iTOL renders as COLOR_NAN
            else:
                f.write(f"{translate_leaf(mag)}\t{val:.2f}\n")
    return fname


def write_06_label_styles(outdir):
    fname = '06_label_styles.txt'
    with _open(outdir, fname) as f:
        f.write("DATASET_STYLE\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tLabel styles\n")
        f.write("COLOR\t#000000\n")
        f.write("DATA\n")
        for leaf in TARGET_MAGS:
            out_leaf = translate_leaf(leaf)
            # MAGs stand out: red, bold, slightly larger, upright
            f.write(f"{out_leaf}\tlabel\tnode\t#e63946\t1.2\tbold\t\n")
        for leaf in REFERENCE_LEAVES:
            out_leaf = translate_leaf(leaf)
            if leaf in BACTERIAL_OUTGROUPS:
                color, size = '#adb5bd', '1'
            else:
                color, size = '#457b9d', '1'
            # Reference genomes go italic
            f.write(f"{out_leaf}\tlabel\tnode\t{color}\t{size}\titalic\t\n")
    return fname


def write_07_dram2_binary(outdir, product_data):
    fname = '07_dram2_binary.txt'

    # Indices of the fields that hold data
    active_field_indices = get_dram2_active_fields()

    # Filter DRAM2_FIELDS, colours and labels
    active_dram2_fields = [DRAM2_FIELDS[i] for i in active_field_indices]
    active_dram2_colors = [DRAM2_FIELD_COLORS[i] for i in active_field_indices]
    field_labels = [f[0] for f in active_dram2_fields]
    n_fields = len(active_dram2_fields)

    # Binary values, for the active fields only
    results = {}
    missing_cols = set()
    for mag in TARGET_MAGS:
        if mag not in product_data:
            print(f" WARNING [07]: {mag} not found en product.tsv", file=sys.stderr)
            results[mag] = [0] * n_fields
            continue
        row = product_data[mag]
        vals = []
        for _, cols, threshold in active_dram2_fields:
            # Report each missing column once, not per leaf
            for c in cols:
                if c not in row and c not in missing_cols:
                    print(
                        f" WARNING [07]: columna no encontrada en product.tsv: '{c}'",
                        file=sys.stderr,
                    )
                    missing_cols.add(c)
            vals.append(eval_field(row, cols, threshold))
        results[mag] = vals

    # Check against the expected values, remapped onto the active indices
    n_mismatches = 0
    for mag, expected_full in DRAM2_EXPECTED.items():
        # Keep only the active fields
        expected = [expected_full[i] for i in active_field_indices]
        observed = results.get(mag)
        if observed is not None and observed != expected:
            n_mismatches += 1
            diff = [
                (field_labels[i], expected[i], observed[i])
                for i in range(n_fields)
                if expected[i] != observed[i]
            ]
            print(
                f" WARNING [07]: en {mag} los valores calculados difieren de los esperados:",
                file=sys.stderr,
            )
            for lbl, exp, obs in diff:
                print(f"  {lbl}: esperado={exp}, calculado={obs}", file=sys.stderr)
    if n_mismatches == 0:
        print(" [07] DRAM2 validation passed: every value matches what was expected.")

    with _open(outdir, fname) as f:
        f.write("DATASET_BINARY\n")
        f.write("SEPARATOR TAB\n")
        f.write("DATASET_LABEL\tPotencial metabólico (DRAM2)\n")
        f.write("COLOR\t#264653\n")
        f.write("FIELD_COLORS\t" + "\t".join(active_dram2_colors) + "\n")
        f.write("FIELD_LABELS\t" + "\t".join(field_labels) + "\n")
        f.write("FIELD_SHAPES\t" + "\t".join(["1"] * n_fields) + "\n")
        f.write("WIDTH\t300\n")
        f.write("MARGIN\t10\n")
        f.write("LEGEND_TITLE\tPotencial metabólico (DRAM2)\n")
        f.write("LEGEND_SHAPES\t1\t1\t1\t1\n")
        f.write("LEGEND_COLORS\t#2a9d8f\t#e9c46a\t#f4a261\t#264653\n")
        f.write("LEGEND_LABELS\tFijación de carbono\tMetabolismo del nitrógeno\tMetabolismo del azufre\tConservación de energía\n")
        f.write("DATA\n")
        for mag in TARGET_MAGS:
            vals = results.get(mag, [0] * n_fields)
            f.write(translate_leaf(mag) + "\t" + "\t".join(str(v) for v in vals) + "\n")
    return fname


def write_08_branch_support_colors(outdir):
    fname = '08_branch_support_colors.txt'
    with _open(outdir, fname) as f:
        f.write("TREE_COLORS\n")
        f.write("SEPARATOR TAB\n")
        f.write("# Resaltado de clados en iTOL mediante rangos MRCA.\n")
        f.write("# Cada línea define el ancestro más reciente común (MRCA) de dos hojas\n")
        f.write("# y colorea el fondo de ese clado en el árbol.\n")
        f.write("#\n")
        f.write("# Formato: leaf1|leaf2\trange\tCOLOR\tLABEL\n")
        f.write("# Si el par de hojas no delimita correctamente el clado deseado,\n")
        f.write("# sustituir por hojas más internas/externas del clado.\n")
        f.write("#\n")
        f.write("DATA\n")
        for leaf1, leaf2, color, label in CLADE_RANGES:
            out1 = translate_leaf(leaf1)
            out2 = translate_leaf(leaf2)
            f.write(f"{out1}|{out2}\trange\t{color}\t{label}\n")
    return fname


def write_09_station_colorstrip(outdir, stations):
    fname = '09_station_colorstrip.txt'
    # The stations the 10 MAGs come from
    present = sorted({stations.get(m, m.split('.')[0]) for m in TARGET_MAGS})
    leg_colors = [STATION_COLORS.get(s, '#999999') for s in present]

    with _open(outdir, fname) as f:
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
        for mag in TARGET_MAGS:
            station = stations.get(mag, mag.split('.')[0])
            color = STATION_COLORS.get(station, '#999999')
            f.write(f"{translate_leaf(mag)}\t{color}\t{station}\n")
    return fname


def write_10_depth_colorstrip(outdir, depths):
    fname = '10_depth_colorstrip.txt'
    # The depths the 10 MAGs come from
    present_depths = []
    seen_d = set()
    for mag in TARGET_MAGS:
        d = depths.get(mag, '')
        if d and d not in seen_d:
            present_depths.append(d)
            seen_d.add(d)

    leg_colors = [DEPTH_COLORS.get(d, ('#999999', d))[0] for d in present_depths]
    leg_labels = [DEPTH_COLORS.get(d, ('#999999', d))[1] for d in present_depths]

    with _open(outdir, fname) as f:
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
        f.write("LEGEND_SHAPES\t" + "\t".join(["1"] * len(present_depths)) + "\n")
        f.write("LEGEND_COLORS\t" + "\t".join(leg_colors) + "\n")
        f.write("LEGEND_LABELS\t" + "\t".join(leg_labels) + "\n")
        f.write("DATA\n")
        for mag in TARGET_MAGS:
            depth = depths.get(mag, '')
            color, label = DEPTH_COLORS.get(depth, ('#999999', depth))
            f.write(f"{translate_leaf(mag)}\t{color}\t{label}\n")
    return fname

# ══════════════════════════════════════════════════════════════════════════════
# POST-GENERATION CHECKS
# ══════════════════════════════════════════════════════════════════════════════


def validate_leaf_coverage(outdir):
    """
    Reads every file just written, pulls the leaf names out of the DATA section
    and compares them against the relabelled names in
    10mags_pmsf.renamed.treefile.
    """
    all_leaf_set = {translate_leaf(x) for x in ALL_LEAVES}
    print("\nLeaf coverage check")
    for fname in sorted(os.listdir(outdir)):
        if not fname.endswith('.txt'):
            continue
        referenced = set()
        in_data = False
        with open(os.path.join(outdir, fname), encoding='utf-8') as fh:
            for line in fh:
                line = line.rstrip('\n')
                if line == 'DATA':
                    in_data = True
                    continue
                if in_data and line and not line.startswith('#'):
                    # The leaf name is the first field
                    leaf = line.split('\t')[0]
                    # TREE_COLORS may address a clade as leaf1|leaf2
                    for part in leaf.split('|'):
                        referenced.add(part.strip())

        missing_from_file = all_leaf_set - referenced  # in the tree, absent from the file
        extra_in_file = referenced - all_leaf_set      # in the file, absent from the tree

        covered = len(referenced & all_leaf_set)
        total = len(all_leaf_set)

        status = "OK" if not extra_in_file and not missing_from_file else "WARN"
        print(f" {fname}: {covered}/{total} hojas cubiertas [{status}]")
        if missing_from_file:
            for m in sorted(missing_from_file):
                print(f"  MISSING, in the tree but not in this file: {m}")
        if extra_in_file:
            for e in sorted(extra_in_file):
                print(f"  EXTRA, not in the tree: {e}")
    print()

# ══════════════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════════════


def main():
    # Check the input files
    missing = [
        p for p in (TAXONOMY_TSV, PRODUCT_TSV, METADATA_TSV, QUALITY_TSV)
        if not os.path.isfile(p)
    ]
    if missing:
        for p in missing:
            print(f"ERROR: archivo de entrada not found: {p}", file=sys.stderr)
        sys.exit(1)

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    print(f"Output directory: {OUTPUT_DIR}")

    # Parse once, reuse everywhere below
    print("Parseando files de entrada...")
    metadata_order, stations, depths, display_names_tree = parse_metadata(METADATA_TSV)
    completeness, contamination = parse_quality(QUALITY_TSV)
    ani_data = parse_ani(TAXONOMY_TSV, set(TARGET_MAGS))
    product_data = parse_product(PRODUCT_TSV, metadata_order)
    print(f" Metadata: {len(metadata_order)} MAGs en orden")
    print(f" Quality: {len(completeness)} MAGs objetivo con datos de calidad")
    print(f" ANI data: {len(ani_data)} hits encontrados para {len(TARGET_MAGS)} MAGs objetivo")
    print(f" DRAM2 data: {len(product_data)} MAGs objetivo encontrados en product.tsv")

    # Build the global map of (renamed) leaf names
    global LEAF_NAME_MAP
    LEAF_NAME_MAP = {}

    # MAG labels come from MAG_DISPLAY_NAMES first, since that was read off the
    # relabelled treefile and therefore matches the tree exactly. The metadata
    # table is the fallback.
    for mag in TARGET_MAGS:
        if mag in MAG_DISPLAY_NAMES:
            LEAF_NAME_MAP[mag] = MAG_DISPLAY_NAMES[mag]
        elif mag in display_names_tree:
            LEAF_NAME_MAP[mag] = display_names_tree[mag].replace(' ', '_')
        else:
            print(
                f" WARNING [leaf_map]: {mag} sin nombre en MAG_DISPLAY_NAMES ni metadata",
                file=sys.stderr,
            )

    # References: use the hard-coded REF_DISPLAY_NAMES
    for ref in REFERENCE_LEAVES:
        new_name = REF_DISPLAY_NAMES.get(ref)
        if new_name:
            LEAF_NAME_MAP[ref] = new_name
        else:
            print(
                f" WARNING [leaf_map]: {ref} sin REF_DISPLAY_NAMES, usando nombre original",
                file=sys.stderr,
            )

    # Write the files
    print("\nWriting anotaciones iTOL...")
    generated = []
    for fn, args in [
        (write_01_source_colorstrip,      (OUTPUT_DIR,)),
        (write_02_taxonomy_colorstrip,    (OUTPUT_DIR,)),
        (write_03_completeness_bar,       (OUTPUT_DIR, completeness, contamination)),
        (write_04_contamination_bar,      (OUTPUT_DIR, contamination)),
        (write_06_label_styles,           (OUTPUT_DIR,)),
        (write_07_dram2_binary,           (OUTPUT_DIR, product_data)),
        (write_08_branch_support_colors,  (OUTPUT_DIR,)),
        (write_09_station_colorstrip,     (OUTPUT_DIR, stations)),
        (write_10_depth_colorstrip,       (OUTPUT_DIR, depths)),
    ]:
        fname = fn(*args)
        path = os.path.join(OUTPUT_DIR, fname)
        size = os.path.getsize(path)
        print(f" [OK] {fname} ({size:,} bytes)")
        generated.append(fname)

    validate_leaf_coverage(OUTPUT_DIR)

    print(f"Listo. {len(generated)} files escritos en:")
    print(f"  {OUTPUT_DIR}")
    print("\nInstrucciones iTOL:")
    print("  1. Sube el treefile (10mags_pmsf.renamed.treefile) a https://itol.embl.de")
    print("  2. Drag each .txt onto the tree to switch that annotation on.")
    print("  3. Ajusta CLADE_RANGES en 08_branch_support_colors.txt si los pares de")
    print("     leaves do not bound the clades correctly in your tree.")

if __name__ == '__main__':
    main()
