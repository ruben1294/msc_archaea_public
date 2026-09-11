# Archaeal MAGs from the deep Gulf of Mexico

Code for a master's thesis on 54 archaeal metagenome-assembled genomes recovered
from six stations and three depths in the deep Gulf of Mexico, at CICESE.

The MAGs fall into four orders: *Poseidoniales* (38), *Nitrososphaerales* (7),
MGIII (6) and *Undinarchaeales* (3). Ten of them carry no species assignment in
GTDB and are candidate new species, which is what much of the analysis here is
built to test. The questions it answers are whether those ten are genuinely
novel and how far from described relatives they sit, which metabolic traits
distinguish the orders, and whether depth or station explains the differences
between genomes of the same species.

Everything in this repository is analysis code. It assumes the MAGs already
exist; assembly and binning happened upstream and are not included.

## What's here

```
scripts/pipeline/   the numbered workflow, 01 to 26
scripts/utils/      helpers the workflow calls between steps
scripts/lib/        shared shell functions: logging, configuration, checks
scripts/R/          downstream analysis and every figure
slurm/              SLURM wrappers, written for the CICESE cluster
configs/            configuration templates and the archaea-calibrated DRAM2 rules
env/                conda environments and the R package list
data/               the small input tables the analysis reads
```

The split matters for reuse: `scripts/pipeline/` and `scripts/utils/` run
anywhere the tools are installed, while `slurm/` is site-specific and meant to be
edited. Nothing under `scripts/` needs a scheduler.

## Pipeline

Each step reads `configs/server.env` and writes into a numbered directory under
the project root. Steps 01 to 15 ran on one cluster, 16 to 26 on another, which
is why there are two configuration profiles.

| Step | What it does | Tool |
|---|---|---|
| 01 | Quality assessment, plus stripping duplicate contig IDs | CheckM2 1.0.2 |
| 02 | Taxonomic classification against GTDB | GTDB-Tk 2.6.1 |
| 03 | Gene calling, one directory per MAG | Prodigal 2.6.3 |
| 04 | Functional annotation of the predicted proteins | eggNOG-mapper 2.1.13 |
| 05 | Annotation and distillation, first pass | DRAM 1.5 |
| 06 | KO assignment. No automation: this ran on the KofamKOALA web server | KofamKOALA |
| 07 | Pathway completeness from the KO table | KEGG-Decoder, KEGG-Expander |
| 08 | Annotation and distillation, second pass. Produces what the R scripts read | DRAM2, Nextflow |
| 09 | De novo archaeal tree placing the MAGs among GTDB references | GTDB-Tk 2.6.1 |
| 10 | Marker alignment under three strategies, carried forward in parallel | MAFFT 7.520 |
| 11 | Column trimming of each alignment | trimAl 1.4.1 |
| 12 | Marker filtering and concatenation into a supermatrix | FASconCAT-G 1.06 |
| 13 | Maximum likelihood trees, three model strategies by three alignments | IQ-TREE 2.3.6 |
| 14 | Pangenome construction, one per species with at least four MAGs | anvi'o 9 |
| 15 | Opens each pangenome for manual figure layout | anvi'o 9 |
| 16 | ANI within the collection and against GTDB | skani 0.2.2 |
| 17 | Chimerism screening against two reference databases | GUNC 1.1 |
| 18-21 | A focused tree for the 10 unassigned MAGs, by two independent routes | GTDB-Tk, MAFFT, trimAl, IQ-TREE 2 |
| 22-26 | KEGG annotation, summaries, metadata layers, functional enrichment, cleanup | anvi'o 9 |

Steps 10 to 13 carry three MAFFT strategies (genafpair, localpair, globalpair)
through to three trees each, nine in total. That is deliberate: agreement across
alignment strategies is what makes the topology credible, and the per-branch
missing-data fractions reported by step 12 are the number to compare first.

## Downstream analysis

The R scripts read the pipeline output and produce the tables and figures. They
run locally rather than on a cluster.

| Script | Covers |
|---|---|
| `MAG_quality_assessment.R` | Assembly and quality metrics, CheckM2 with GUNC |
| `tax_analysis.R` | Taxonomic composition, novelty, and the skani ANI figures |
| `DRAM2_output_analysis.R` | Functional annotation: rates, traits, metabolism |
| `CNPS.cycle.R` | Carbon, nitrogen, phosphorus and sulfur cycle pathways |
| `treatment_comparison.R` | Traits aggregated by group, Fisher tests, logistic models |
| `ANI_intraspecies.R` | Within-species ANI, and what depth and station explain |
| `phylogenetic_trees_analysis.R` | Compares the nine trees from step 13 |
| `26` to `31_pangenome_*.R` | Pangenome classification, statistics, environmental RDA, drivers |
| `MDMS_figs.R` | The poster figure set, in English |
| `_pub_style.R` | Shared fonts, theme, palettes, and path resolution |

Figures are labelled in Spanish, matching the thesis. `MDMS_figs.R` is the
exception, since its poster was in English.

## Requirements

```bash
# One conda environment per tool
conda env create -f env/checkm2.yml       # and the rest, see env/README.md

# R packages, around 30 of them
Rscript env/r-packages.R
```

`env/README.md` lists which environment each step needs.

Reference databases are downloaded separately and are large: GTDB, eggNOG,
CheckM2, DRAM2, GUNC and the skani GTDB sketches. `scripts/utils/` has a
download script for each. **The release matters.** These trees were built
against GTDB r220 and the ANI computed against r226; mixing releases produces
inconsistent placements without any error to warn you.

## Running it

```bash
git clone https://github.com/ruben1294/gom-archaeal-mags.git
cd gom-archaeal-mags

cp configs/server.env.example configs/server.env
$EDITOR configs/server.env          # paths, databases, thresholds

# One step
bash scripts/pipeline/01_mag_quality.sh

# Or through SLURM, from the repository root
cp slurm/cluster.env.example slurm/cluster.env
$EDITOR slurm/cluster.env           # module and conda profile
$EDITOR slurm/01_mag_quality.slurm  # account and partition
mkdir -p logs && sbatch slurm/01_mag_quality.slurm
```

Steps 01 and 02 take `--help` and a set of path options; the rest read
everything from `configs/server.env`. All of them check their configuration and
their tools before doing any work, so a missing database or an unactivated
environment fails on the first line rather than an hour in.

The R scripts find the project root through `MSC_ARCHAEA_DIR`:

```bash
MSC_ARCHAEA_DIR=/path/to/project Rscript scripts/R/tax_analysis.R
```

## Data availability

`data/` holds the small tables the analysis reads directly:

- `MAGs_metadata.tsv`, one row per MAG with station, depth, bin identifier and
  GTDB taxonomy. Every script keys off its `NameWithoutExtension` and
  `DisplayName` columns.
- `Depths.xlsx` and `Physicochemical_properties.xlsx`, the environmental
  variables per station and depth that the RDA in `29_pangenome_env_RDA.R` uses.

The 54 MAG assemblies and the large pipeline intermediates (`raw-annotations.tsv`,
`traits.xlsx`, the anvi'o databases) are not in this repository. They will be
deposited when the manuscript is submitted, and this section will name the
accession.

## Citing

`CITATION.cff` carries the machine-readable form. GitHub renders it as a "Cite
this repository" link.

## License

MIT. See [LICENSE](LICENSE).
