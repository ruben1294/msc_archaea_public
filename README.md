# Downstream analysis of archaeal MAGs from the Gulf of Mexico

Code for my master's thesis at CICESE on 54 archaeal metagenome-assembled genomes
recovered from six stations and three depths in the deep Gulf of Mexico.

Taxonomically, the MAGs belong to four orders: *Poseidoniales* (38), *Nitrososphaerales* (7),
MGIII (6) and *Undinarchaeales* (3). Ten of them weren't assigned to any known species in
GTDB. The downstream analysis of this repo tests whether those
ten are new and how distant they are from described relatives, which metabolic
traits distinguish the four orders, and whether depth or station explains the
differences between genomes of the same species.

The code starts from the finished MAGs.

## Credits

The contig assembly and the binning that produced the MAGs analysed here are the
work of Dr. Lizt Selene Sibila Osorio Pando, in her doctoral thesis:

> Osorio Pando, L. S. S. (2025). *Estructura y función de las comunidades de
> bacterias y arqueas en las zonas mesopelágica y batipelágica del golfo de
> México* [Doctoral thesis, Doctorado en Ciencias en Ciencias de la Vida con
> orientación en Microbiología]. Centro de Investigación Científica y de
> Educación Superior de Ensenada, Baja California (CICESE), Ensenada, México.

The environmental tables under `data/` come from the same thesis:
`Depths.xlsx` is its Tabla 9 for the six stations the MAGs come from, and the
first sheet of `Physicochemical_properties.xlsx` is its Tabla 2 for the three
depths sampled here.

## What's here

```
scripts/pipeline/   the numbered workflow, 01 to 26
scripts/utils/      helpers the workflow calls between steps, and setup scripts
scripts/lib/        shared shell functions: logging, configuration, checks
scripts/R/          downstream analysis and every figure
slurm/              one SLURM wrapper per step, account and partition are given at submission
configs/            the configuration template, server.env.example
env/                conda environments and the R package list
data/               the small input tables the analysis reads
```

Nothing under `scripts/` needs a scheduler, the scripts run on any machine with
their tools installed. Cluster settings (modules, conda location, account,
partition) go in `slurm/cluster.env` and on the `sbatch` command line.

## Pipeline

Each step reads `configs/server.env` and writes into `results/NN_<name>/`, with
`NN` the step number. Steps run in dependency order: quality; taxonomy, with the
ANI that tests whether the MAGs without a species are new; gene calls and the
annotation built on them; phylogenomics, first the 54-MAG trees and then the
focused tree for the 10 unassigned MAGs; and last the pangenomes of the species
with enough (n ≥ 5) genomes.

| Step | What it does | Tool and version used in the thesis |
|---|---|---|
| 01 | Quality assessment, plus stripping duplicate contig IDs | CheckM2 |
| utils/01 | MIMAG filter: keeps the 54 of 64 bins with at least 50% completeness and under 10% contamination | `01_filter_mags.py` |
| utils/05 | Renames the MAGs, underscores to dots, so every later step sees the final names | `05_rename_mags.sh` |
| 02 | Chimerism screening against ProGenomes3 and GTDB r214 | GUNC 1.1.0 |
| 03 | Taxonomic classification | GTDB-Tk 2.5.2, GTDB R226 |
| utils/06 | Builds the MAG metadata table from the names and the GTDB-Tk output | `06_generate_metadata.py` |
| 04 | ANI within the collection and against GTDB, to test whether the 10 MAGs without a species are new | skani, GTDB R226 sketches |
| 05 | Gene calling, one directory per MAG | Prodigal 2.6.3 |
| 06 | Functional annotation of the predicted proteins | eggNOG-mapper 2.1.12 |
| 07 | Annotation and distillation, first pass | DRAM v1 |
| 08 | KO assignment. No automation: this ran on the KofamKOALA web server | KofamKOALA |
| 09 | Pathway completeness from the KO table | KEGG-Decoder, KEGG-Expander |
| 10 | Annotation and distillation, second pass. Produces what the R scripts read | DRAM2 (pinned fork), Nextflow 25.10.2 |
| 11 | De novo archaeal tree placing the MAGs among GTDB references | GTDB-Tk 2.5.2, GTDB R226 |
| 12 | Marker alignment under three strategies, carried forward in parallel | MAFFT |
| 13 | Column trimming of each alignment | trimAl 1.5.0 |
| 14 | Marker filtering and concatenation into a supermatrix | FASconCAT-G 1.06.1 |
| 15 | Maximum likelihood trees, three model strategies by three alignments | IQ-TREE 2.2.5 |
| 16-19 | A focused tree for the 10 unassigned MAGs, by two independent routes | GTDB-Tk 2.6.1, MAFFT, trimAl 1.5.0, IQ-TREE 2.2.5 |
| 20 | Pangenome construction, one per species with at least four MAGs | anvi'o 9 |
| 21-25 | KEGG and Pfam annotation, summaries, metadata layers, functional enrichment, cleanup | anvi'o 9 |
| 26 | Opens each cleaned pangenome for manual figure layout | anvi'o 9 |

GUNC flags MAGs here and
removes none. Of the 54 MAGs, 22 exceed its default threshold against GTDB r214
and 24 against ProGenomes3. Most of them have low reference representation, the
case in which GUNC itself warns that its scores are unreliable, and deep-sea
archaea are poorly represented in both references. All 54 MAGs continue to step
03, and `MAG_quality_assessment.R` plots the GUNC results in its figures 8 to 12.

Steps 12 to 15 carry three MAFFT strategies (genafpair, localpair & globalpair)
through to three trees each, nine in total, to test whether the topology holds
across alignment strategies. Start that comparison with the per-branch
missing-data fractions that step 14 reports.

Step 10 runs a fork of DRAM2, [ruben1294/DRAM](https://github.com/ruben1294/DRAM)
(branch `dev`), not the original [WrightonLabCSU/DRAM](https://github.com/WrightonLabCSU/DRAM).
The fork starts from upstream commit `90dfef6f` (2026-04-14) and adds:

- fixes to the KEGG database formatting;
- trait rules calibrated for archaea, `bin/assets/traits_rules_archaea.tsv`,
  derived from DRAM2's own and used by default. Every changed rule has an
  `ARCHAEA-MOD` note explaining it;
- the `dram_viz` container reverted to the version that still writes
  `product.tsv`;
- the run configuration of this step, `conf/gom_archaea_mags.config`, with a
  shorter test variant and a generic SLURM profile.

The fork's README lists each change with its date. `DRAM2_REVISION` in
`configs/server.env` pins commit `68a2d53`: the commit the thesis analyses ran
(`aa5889d`) plus the configuration files, with no code changes.
`slurm/DRAM2_HPC_install.slurm` clones it, and step 10 refuses any other commit
unless `DRAM2_ALLOW_OTHER_REVISION=true`.

## Downstream analysis

The R scripts read the pipeline output and produce the tables and figures. They
run locally, from any directory inside the clone: each one finds the repository
root itself, or takes it from `MSC_ARCHAEA_DIR` when that is set.

```bash
Rscript scripts/R/tax_analysis.R
```

| Script | Covers |
|---|---|
| `MAG_quality_assessment.R` | Assembly and quality metrics, CheckM2 with GUNC |
| `tax_analysis.R` | Taxonomic composition, novelty, and the skani ANI figures |
| `DRAM2_output_analysis.R` | Functional annotation: rates, traits, metabolism |
| `CNPS.cycle.R` | Carbon, nitrogen, phosphorus and sulfur cycle pathways |
| `treatment_comparison.R` | Traits aggregated by group, Fisher tests, logistic models |
| `ANI_intraspecies.R` | Within-species ANI, and what depth and station explain |
| `phylogenetic_trees_analysis.R` | Compares the nine trees from step 15 |
| `26` to `31_pangenome_*.R` | Pangenome classification, statistics, environmental RDA, drivers |
| `MDMS_figs.R` | The poster figure set, in English |
| `_pub_style.R` | Shared fonts, theme, palettes, and path resolution |

Figures are labelled in Spanish.

`MAG_quality_assessment.R` also compares the 54 MAGs with a second MAG set from
the same samples. This pipeline does not produce that set, pass its CheckM2
report with `--file1`. Without it, the script skips the comparison figures and
runs the rest.

The R outputs go under `results/` next to the pipeline's: `SO_01/` holds the
taxonomy, ANI, functional annotation and nutrient-cycle results, `SO_02/` the
phylogenomics, `SO_03/` the pangenomics, one directory per specific objective of
my thesis.

## Requirements

```bash
# One conda environment per tool
conda env create -f env/checkm2.yml       # and the rest, see env/README.md

# R packages
Rscript env/r-packages.R
```

`env/README.md` lists which environment each step needs and the version the
thesis ran. The environment files are specifications, not lockfiles, so a fresh
`conda env create` can resolve newer dependencies.

Reference databases are downloaded separately and are large. `env/README.md`
says where each comes from, and `scripts/utils/` has download or setup scripts
for GUNC, the skani sketches, DRAM2 and anvi'o. Use **GTDB R226** for both
GTDB-Tk and the skani search: mixing releases gives inconsistent placements and
raises no error.

## KEGG licensing

KEGG is not freely redistributable, its FTP data need a licence (see
[KEGG's legal notice](https://www.kegg.jp/kegg/legal.html)). The website,
the KEGG API and KofamKOALA are free for academic use only. The following steps use them:

- Step 10 annotated against a 2018 release of the KEGG database, from KEGG's
  licensed FTP distribution, enabled by `use_kegg` in the fork's
  `conf/gom_archaea_mags.config`. Without a licence, set
  `use_kegg = false`: DRAM2 still assigns KOs through KOfam, but the distillate,
  and the figures built on it, will differ from the thesis.
- Step 07 uses KEGG only if the DRAM v1 databases were prepared with a licensed
  KEGG peptide file. Otherwise, it too falls back to KOfam.
- Step 08 ran on the KofamKOALA web server, step 21 uses the KOfam and KEGG
  snapshot that anvi'o downloads, and the R scripts can look KEGG names up
  through the KEGG API with KEGGREST.

No KEGG data is redistributed in this repository.

## Running it

```bash
git clone https://github.com/ruben1294/gom-archaeal-mags.git
cd gom-archaeal-mags

cp configs/server.env.example configs/server.env
$EDITOR configs/server.env           # database paths; the project layout works as is

mkdir -p data/mags
cp /path/to/bins/*.fa data/mags/     # the input bins, one .fa per genome

# The quality block
bash scripts/pipeline/01_mag_quality.sh
python3 scripts/utils/01_filter_mags.py
bash scripts/utils/05_rename_mags.sh
bash scripts/pipeline/02_chimerism.sh

# and onward, one step at a time
bash scripts/pipeline/03_taxonomic_assignment.sh
```

Results are placed in `results/`, logs in `logs/` and scratch files in `temp/`, all
under the clone unless `PROJECT_DIR` in `configs/server.env` directs elsewhere. Each step checks its configuration and its
tools before it starts, so a missing database or an unactivated environment
stops it at once. Steps 01 and 03 take `--help`,
step 02 takes `--clean` to start over, and step 04 `--skip-gtdb` to run only the
comparisons within the collection.

If you're using the scheduler SLURM , then from the repository root:

```bash
cp slurm/cluster.env.example slurm/cluster.env
$EDITOR slurm/cluster.env            # module and conda setup

export SBATCH_ACCOUNT=<account> SBATCH_PARTITION=<partition>
mkdir -p logs
sbatch slurm/01_mag_quality.slurm
```

You need to set account, partition and node to `sbatch`, or export
`SBATCH_ACCOUNT` and `SBATCH_PARTITION` once as above.

The thesis used two clusters: steps 20 to 26 ran on the second. If you split steps across machines, keep one
`server.env` per machine and direct `CONFIG_FILE` at the one each step should
read.

## Data availability

`data/` holds the small tables the analysis reads directly:

- `MAGs_metadata.tsv`, one row per MAG with station, depth, bin identifier and
  GTDB taxonomy. The R scripts join on its `NameWithoutExtension` and
  `DisplayName` columns. `scripts/utils/06_generate_metadata.py` rebuilds it
  from the MAG names and the GTDB-Tk output.
- `Depths.xlsx` and `Physicochemical_properties.xlsx`, the environmental
  variables per station and depth that the RDA in `29_pangenome_env_RDA.R` uses
  (see [Credits](#credits)).

## Citing

`CITATION.cff` carries the machine-readable form. GitHub renders it as a "Cite
this repository" link.

## License

MIT. See [LICENSE](LICENSE).

DRAM2 is licensed under the GNU GPL v3, so the configuration of step 10 and the
archaea-calibrated trait rules are in the fork under that licence, and step 10
fetches them from there. This repository contains no code derived from DRAM2. If you use step 10, cite DRAM as well: Shaffer et al. (2020),
*Nucleic Acids Research* 48(16):8883-8900,
[doi:10.1093/nar/gkaa621](https://doi.org/10.1093/nar/gkaa621).
