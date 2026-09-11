#!/usr/bin/env python3
"""Build MAGs_metadata.tsv from the MAG filenames and the GTDB-Tk taxonomy.

The filenames encode where each MAG came from: <station>.<depth>.megahit.bin.<n>,
so A4.BTM.megahit.bin.25 is bin 25 from station A4 at the bottom depth. This
splits that apart and joins it against the classify_wf summary.

The resulting table is the spine of the whole project. Every R script keys off
its NameWithoutExtension and DisplayName columns, and pipeline steps 16, 21 and
24 read it too.
"""

import os
from datetime import datetime
from pathlib import Path

import pandas as pd

import os
from pathlib import Path

# Project root, from the environment or the usual location. Every path below
# derives from it, so the script runs wherever the tree happens to live.
PROJECT_DIR = Path(os.environ.get("MSC_ARCHAEA_DIR", Path.home() / "msc_archaea"))

MAG_DIR = str(PROJECT_DIR / "01_raw_files" / "filtered")
TAXONOMY_DIR = str(PROJECT_DIR / "03_taxonomy_assigned" / "classify")
TAXONOMY_FILE = "tax_.ar53.summary.tsv"
OUTPUT_DIR = str(PROJECT_DIR)
EXTENSION = ".fa"
OUTPUT_TSV = "MAGs_metadata.tsv"

def extract_taxonomy_level(taxonomy_string, level_prefix):
    """
    Extract specific taxonomic level from GTDB taxonomy string
    
    Args:
        taxonomy_string: Full GTDB taxonomy (e.g., "d__Archaea;p__Thermo...;g__MGIIb-Q1;s__")
        level_prefix: Taxonomic level prefix (e.g., "g__" for genus, "s__" for species)
    
    Returns:
        Taxonomic name at specified level, or None if not found/empty
    """
    if pd.isna(taxonomy_string) or taxonomy_string == "":
        return None
    
    # Split by semicolon
    levels = taxonomy_string.split(";")
    
    # Find the level
    for level in levels:
        level = level.strip()
        if level.startswith(level_prefix):
            # Remove prefix
            taxon = level.replace(level_prefix, "").strip()
            # Return None if empty
            return taxon if taxon else None
    
    return None

def create_display_name(row):
    """
    Create display name for figures
    
    Format: Species Depth_bin#
    If no species: Genus Depth_bin#
    
    Example: 
        - "adriaticus 500m_001"
        - "MGIIb-Q1 1000m_042"
    """
    species = row.get("Species")
    genus = row.get("Genus")
    station = row.get("Station", "UNK")
    depth = row.get("Depth", "UNK")
    bin_id = row.get("BinID", "000")
    
    # Use species if available, otherwise genus
    if species and species != "N/A":
        taxon_name = species
    elif genus and genus != "N/A":
        taxon_name = genus
    else:
        taxon_name = "Unclassified"
    
    # Create display name
    display_name = f"{taxon_name} {station} {depth}_{bin_id}"
    
    return display_name

def procesar_mags(carpeta):
    """
    Process MAG files and extract metadata from filenames
    
    Args:
        carpeta: Directory containing MAG files
    
    Returns:
        pandas.DataFrame with metadata
    """
    # Get all .fa files
    archivos = [f for f in os.listdir(carpeta) if f.endswith(EXTENSION)]
    archivos.sort()
    
    if not archivos:
        print(f"Warning: No files with extension '{EXTENSION}' found in {carpeta}")
        return pd.DataFrame()
    
    datos = []
    
    for archivo in archivos:
        # Remove extension
        nombre_sin_ext = archivo.replace(EXTENSION, "")
        
        # Split by dot
        partes = nombre_sin_ext.split(".")
        
        # Parse according to expected format
        if len(partes) >= 5:
            # Formato: Station.Depth.Assembler.bin.BinNumber
            estacion = partes[0]
            profundidad = partes[1]
            ensamblador = partes[2]
#            binner = partes[3]
            bin_id_raw = ".".join(partes[4:])
            
        elif len(partes) == 5 and partes[3] == "bin":
            # Formato específico: Station.Depth.Assembler.bin.BinNumber
            estacion = partes[0]
            profundidad = partes[1]
            ensamblador = partes[2]
            binner = "bin"
            bin_id_raw = partes[4]
            
        elif len(partes) == 4:
            # Formato sin binner
            estacion = partes[0]
            profundidad = partes[1]
            ensamblador = partes[2]
            binner = "unknown"
            bin_id_raw = partes[3]
            
        else:
            # Formato desconocido
            estacion = partes[0] if len(partes) > 0 else "unknown"
            profundidad = partes[1] if len(partes) > 1 else "unknown"
            ensamblador = partes[2] if len(partes) > 2 else "unknown"
            binner = "unknown"
            bin_id_raw = ".".join(partes[3:]) if len(partes) > 3 else "unknown"
        
        # Clean bin_id (remove "bin." prefix if present)
        bin_id_clean = bin_id_raw.replace("bin.", "")
        
        # Append data
        datos.append({
            "FullName": archivo,
            "Station": estacion,
            "Depth": profundidad,
            "Assembler": ensamblador,
#            "Binner": binner,
            "BinID": bin_id_clean,
            "NameWithoutExtension": nombre_sin_ext,
            "OriginalName": nombre_sin_ext
        })
    
    return pd.DataFrame(datos)

# ==========================================
# MAIN EXECUTION
# ==========================================

def main():
    print("=" * 70)
    print("MAG Metadata Extraction with Taxonomy Integration")
    print("=" * 70)
    print()
    
    # Check if MAG directory exists
    if not os.path.isdir(MAG_DIR):
        print(f"[FAIL] Error: MAG directory not found: {MAG_DIR}")
        return 1
    
    print(f"MAG directory: {MAG_DIR}")
    print(f"Taxonomy directory: {TAXONOMY_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print()
    
    # Create output directory if it doesn't exist
    os.makedirs(OUTPUT_DIR, exist_ok=True)
    
    # Process MAGs
    print("Processing MAG files...")
    df_mags = procesar_mags(MAG_DIR)
    
    if df_mags.empty:
        print("[FAIL] Error: No MAG files processed")
        return 1
    
    print(f"[OK] Parsed {len(df_mags)} MAG files")
    print()
    
    # Load taxonomy file
    taxonomy_path = os.path.join(TAXONOMY_DIR, TAXONOMY_FILE)
    
    if os.path.isfile(taxonomy_path):
        print(f"Loading taxonomy from: {taxonomy_path}")
        
        try:
            # Read taxonomy file (assuming tab-separated)
            # First column: MAG name (user_genome)
            # Second column: classification
            df_taxonomy = pd.read_csv(taxonomy_path, sep='\t')
            
            print(f"[OK] Loaded taxonomy for {len(df_taxonomy)} MAGs")
            print(f"  Columns in taxonomy file: {list(df_taxonomy.columns)}")
            

            mag_col = df_taxonomy.columns[0]  # First column is MAG name
            tax_col = df_taxonomy.columns[1]  # Second column is classification
            
            print(f"  Using columns: '{mag_col}' (MAG) and '{tax_col}' (taxonomy)")
            print()
            
            # Rename columns for easier handling
            df_taxonomy = df_taxonomy.rename(columns={
                mag_col: 'NameWithoutExtension',
                tax_col: 'FullTaxonomy'
            })
            
            # Extract genus and species
            print("Extracting taxonomic levels...")
            df_taxonomy['Genus'] = df_taxonomy['FullTaxonomy'].apply(
                lambda x: extract_taxonomy_level(x, "g__")
            )
            df_taxonomy['Species'] = df_taxonomy['FullTaxonomy'].apply(
                lambda x: extract_taxonomy_level(x, "s__")
            )
            
            # Replace None with "N/A"
            df_taxonomy['Genus'] = df_taxonomy['Genus'].fillna("N/A")
            df_taxonomy['Species'] = df_taxonomy['Species'].fillna("N/A")
            
            # Merge with MAG metadata
            df_mags = df_mags.merge(
                df_taxonomy[['NameWithoutExtension', 'FullTaxonomy', 'Genus', 'Species']],
                on='NameWithoutExtension',
                how='left'
            )
            
            # Fill missing taxonomy with "N/A"
            df_mags['FullTaxonomy'] = df_mags['FullTaxonomy'].fillna("N/A")
            df_mags['Genus'] = df_mags['Genus'].fillna("N/A")
            df_mags['Species'] = df_mags['Species'].fillna("N/A")
            
            n_with_taxonomy = (df_mags['FullTaxonomy'] != "N/A").sum()
            n_with_genus = (df_mags['Genus'] != "N/A").sum()
            n_with_species = (df_mags['Species'] != "N/A").sum()
            
            print(f"[OK] Taxonomy merged:")
            print(f"  MAGs with taxonomy: {n_with_taxonomy}/{len(df_mags)}")
            print(f"  MAGs with genus: {n_with_genus}/{len(df_mags)}")
            print(f"  MAGs with species: {n_with_species}/{len(df_mags)}")
            print()
            
        except Exception as e:
            print(f"[FAIL] Error loading taxonomy file: {e}")
            print("  Continuing without taxonomy...")
            df_mags['FullTaxonomy'] = "N/A"
            df_mags['Genus'] = "N/A"
            df_mags['Species'] = "N/A"
    else:
        print(f"Warning: Taxonomy file not found: {taxonomy_path}")
        print("  Continuing without taxonomy...")
        df_mags['FullTaxonomy'] = "N/A"
        df_mags['Genus'] = "N/A"
        df_mags['Species'] = "N/A"
    
    # Create display name for figures
    print("Creating display names for figures...")
    df_mags['DisplayName'] = df_mags.apply(create_display_name, axis=1)
    print("[OK] Display names created")
    print()
    
    # Reorder columns
    column_order = [
        'FullName',
        'NameWithoutExtension',
        'DisplayName',
        'Station',
        'Depth',
        'Assembler',
#        'Binner',
        'BinID',
        'Genus',
        'Species',
        'FullTaxonomy',
        'OriginalName'
    ]
    
    df_mags = df_mags[column_order]
    
    # Save to TSV
    output_path = os.path.join(OUTPUT_DIR, OUTPUT_TSV)
    df_mags.to_csv(output_path, sep='\t', index=False)
    
    # Print summary
    print("=" * 70)
    print("SUMMARY")
    print("=" * 70)
    print(f"[OK] Total MAGs processed: {len(df_mags)}")
    print(f"[OK] TSV file saved to: {output_path}")
    print()
    
    # Print unique values for each column
    print("Metadata summary:")
    print(f"  Stations: {df_mags['Station'].nunique()} ({', '.join(sorted(df_mags['Station'].unique()))})")
    print(f"  Depths: {df_mags['Depth'].nunique()} ({', '.join(sorted(df_mags['Depth'].unique()))})")
    print(f"  Assemblers: {df_mags['Assembler'].nunique()} ({', '.join(sorted(df_mags['Assembler'].unique()))})")
#    print(f"  Binners: {df_mags['Binner'].nunique()} ({', '.join(sorted(df_mags['Binner'].unique()))})")
    
    if 'Genus' in df_mags.columns:
        n_genera = (df_mags['Genus'] != "N/A").sum()
        print(f"  Unique genera: {df_mags[df_mags['Genus'] != 'N/A']['Genus'].nunique()} ({n_genera} MAGs)")
    
    if 'Species' in df_mags.columns:
        n_species = (df_mags['Species'] != "N/A").sum()
        print(f"  Unique species: {df_mags[df_mags['Species'] != 'N/A']['Species'].nunique()} ({n_species} MAGs)")
    
    print()
    
    # Show first 10 records
    print("First 10 records:")
    display_cols = ['FullName', 'DisplayName', 'Station', 'Depth', 'Genus', 'Species']
    print(df_mags[display_cols].head(10).to_string(index=False))
    print()
    
    # Show examples of display names
    print("Example display names:")
    for i, name in enumerate(df_mags['DisplayName'].head(5), 1):
        print(f"  {i}. {name}")
    print()
    
    return 0

if __name__ == "__main__":
    exit(main())
