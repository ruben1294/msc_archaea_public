#!/usr/bin/env python3
"""Rewrite MAG identifiers in tree files to their display names.

IQ-TREE carries the input filenames through to its output, so the trees come out
labelled with NameWithoutExtension. This substitutes the DisplayName from
MAGs_metadata.tsv, which is what belongs on a figure.

It rewrites the files in place across every .contree, .iqtree and .treefile it
finds, so keep a copy if the originals matter.
"""

import os
import re
from pathlib import Path

import pandas as pd

# Project root, from the environment or the usual location. Every path below
# derives from it, so the script runs wherever the tree happens to live.
PROJECT_DIR = Path(os.environ.get("MSC_ARCHAEA_DIR", Path.home() / "msc_archaea"))

METADATA_FILE = str(PROJECT_DIR / "MAGs_metadata.tsv")
TREE_DIR = os.environ.get("TREE_DIR", str(PROJECT_DIR / "results" / "13_phylogeny"))
TREE_EXTENSIONS = [".contree", ".iqtree", ".treefile"]


def load_metadata(metadata_file):
    """Read MAGs_metadata.tsv and map NameWithoutExtension to DisplayName."""
    print(f"Reading metadata from: {metadata_file}")

    metadata = pd.read_csv(metadata_file, sep="\t")

    # Crear diccionario: NameWithoutExtension -> DisplayName
    name_mapping = dict(zip(metadata["NameWithoutExtension"], metadata["DisplayName"]))

    print(f"  [OK] {len(name_mapping)} MAG names loaded")
    return name_mapping


def replace_names_in_file(file_path, name_mapping):
    """Substitute every MAG name in one tree file, in place."""
    print(f"\nProcessing: {file_path.name}")

    with open(file_path, "r") as f:
        content = f.read()

    original_content = content
    replacements_made = 0

    # Reemplazar cada nombre
    for old_name, new_name in name_mapping.items():
        # Match the whole identifier only. In Newick a name is followed by a
    # colon, a closing paren, a comma, a semicolon or the end of the line, and
    # anchoring on those keeps one MAG name from matching inside another.
        pattern = re.compile(rf"\b{re.escape(old_name)}\b")

        if pattern.search(content):
            content = pattern.sub(new_name, content)
            replacements_made += 1

    # Solo escribir si hubo cambios
    if content != original_content:
        # Crear backup
        backup_path = file_path.with_suffix(file_path.suffix + ".backup")
        with open(backup_path, "w") as f:
            f.write(original_content)
        print(f"  Backup saved: {backup_path.name}")

        # Write the rewritten file back in place
        with open(file_path, "w") as f:
            f.write(content)

        print(f"  [OK] {replacements_made} names replaced")
        return replacements_made
    else:
        print("  no MAG names matched in this file")
        return 0


def process_tree_files(tree_dir, name_mapping, extensions):
    """Walk the tree directory and rewrite every file with a matching extension."""
    tree_dir = Path(tree_dir)

    if not tree_dir.exists():
        print(f"ERROR: El directorio {tree_dir} no existe")
        return

    total_files = 0
    total_replacements = 0

    for ext in extensions:
        tree_files = list(tree_dir.rglob(f"*{ext}"))

        print(f"\n{'=' * 70}")
        print(f"Looking for {ext} files...")
        print(f"  found {len(tree_files)}")
        print(f"{'=' * 70}")

        for tree_file in tree_files:
            replacements = replace_names_in_file(tree_file, name_mapping)
            total_replacements += replacements
            total_files += 1

    # Resumen final
    print(f"\n{'=' * 70}")
    print("RESUMEN")
    print(f"{'=' * 70}")
    print(f"Files rewritten: {total_files}")
    print(f"Total de reemplazos realizados: {total_replacements}")
    print("\n[OK] Finished")
    print("\nA .backup copy was written next to each file.")

# ---- EJECUCIÓN PRINCIPAL ----

if __name__ == "__main__":
    print("=" * 70)
    print("RENOMBRADO DE MAGs EN ÁRBOLES FILOGENÉTICOS")
    print("=" * 70)
    print()

    # Cargar metadatos
    name_mapping = load_metadata(METADATA_FILE)

    # Walk the tree files
    process_tree_files(TREE_DIR, name_mapping, TREE_EXTENSIONS)
