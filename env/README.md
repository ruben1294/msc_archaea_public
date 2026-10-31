# Environments

Each pipeline step runs inside its own conda environment. The `.yml` files here
recreate them:

```bash
conda env create -f env/<name>.yml
```

For the R side, one script installs everything at once:

```bash
Rscript env/r-packages.R
```

| Environment | Steps | Provides | Version the thesis ran |
|---|---|---|---|
| `checkm2.yml` | 01, utils/01 | CheckM2 | not recorded in the run logs |
| `gunc.yml` | 02 | GUNC | 1.1.0 |
| `gtdbtk.yml` | 03, 11, 16, 19 | GTDB-Tk | 2.5.2 for 03 and 11, 2.6.1 for 16 and 19 |
| `skani.yml` | 04 | skani | not recorded in the run logs |
| `prodigal.yml` | 05 | Prodigal | 2.6.3 |
| `eggnogmapper.yml` | 06 | eggNOG-mapper | 2.1.12 |
| `dram.yml` | 07 | DRAM v1 | not recorded in the run logs |
| `keggdecoder.yml` | 09 | KEGG-Decoder and KEGG-Expander | not recorded in the run logs |
| `dram2.yml` | 10 | Nextflow, for the DRAM2 workflow | Nextflow 25.10.2, DRAM2 fork at `aa5889d` (configs added in `68a2d53`) |
| `phylo.yml` | 12 to 15, 17, 18 | MAFFT, trimAl, IQ-TREE 2 | trimAl 1.5.0, IQ-TREE 2.2.5; MAFFT not recorded |
| `anvio9.yml` | 20 to 26 | anvi'o 9 | 9 |

Step 08 has no environment: it ran on the KofamKOALA web server. FASconCAT-G,
for steps 14 and 17, is a single Perl script with no conda package;
`scripts/utils/install_FASconCAT-G.sh` fetches the commit the thesis used,
version 1.06.1.

The environment files are specifications, not lockfiles. Where the thesis
version is known it is pinned, but `conda env create` still resolves the
dependencies around it afresh, and where the version was not recorded the file
names the one it has always named. For a byte-identical rebuild, export the
environments you end up with (`conda env export`) and keep those.

## Reference databases

None of these ship with the environments, and most are large. The release
matters: GTDB-Tk ran against GTDB R226, and so did the skani search, and mixing
releases produces inconsistent placements without any error to warn you. Each
path goes into `configs/server.env`.

| Database | Used by | Release | How to get it | Variable |
|---|---|---|---|---|
| CheckM2 DIAMOND database | 01 | as `checkm2 database` fetches it | `checkm2 database --download --path <dir>` | `CHECKM2_DB` |
| GUNC ProGenomes3 and GTDB r214 | 02 | GUNC 1.1 | `scripts/utils/download_gunc_1.1_db.sh` | `GUNC_DB` |
| GTDB-Tk reference package | 03, 11, 16, 19 | R226 | `download-db.sh`, installed with GTDB-Tk | `GTDBTK_DB` |
| eggNOG data | 06 | as eggNOG-mapper 2.1.12 fetches it | `download_eggnog_data.py --data_dir <dir>` | `EGGNOG_DATA_DIR` |
| DRAM v1 databases | 07 | | `DRAM-setup.py prepare_databases --output_dir <dir>` | `DRAM_DB_DIR` |
| DRAM2 databases | 10 | | `scripts/utils/DRAM2_HPC_download_dbs.sh`, over Globus | `DRAM2_DB_DIR` |
| skani GTDB sketches | 04 | R226 | `scripts/utils/download_skani_gtdb_r226.sh` | `GTDB_SKETCH_DIR` |
| anvi'o COGs, KEGG and Pfam | 20, 21 | as anvi'o 9 fetches them | `scripts/utils/setup_anvio9.sh` | `KEGG_DATA_DIR`, `PFAM_DATA_DIR` |

KEGG itself is licensed; README.md, under "KEGG licensing", says which steps use
it and what to do without a licence.

## Two clusters

The work ran on two clusters. anvi'o 9 would not build on the first, over its
`datrie` dependency, so the pangenome steps, 20 to 26, ran on a
second one. Nothing in the code depends on that split: keep one copy of
`configs/server.env` per machine and set `CONFIG_FILE` to the one each step
should read.
