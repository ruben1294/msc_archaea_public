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

| Environment | Steps | Provides |
|---|---|---|
| `checkm2.yml` | 01, utils/01 | CheckM2 |
| `gtdbtk.yml` | 02, 09, 18, 21 | GTDB-Tk 2.6.1 |
| `eggnogmapper.yml` | 04 | eggNOG-mapper 2 |
| `dram.yml` | 05 | DRAM v1 |
| `dram2.yml` | 08 | Nextflow and the DRAM2 workflow |
| `keggdecoder.yml` | 07 | KEGG-Decoder and KEGG-Expander |
| `phylo.yml` | 10 to 13, 19, 20 | MAFFT, trimAl, IQ-TREE 2 |
| `prodigal.yml` | 03 | Prodigal |
| `skani.yml` | 16 | skani |
| `gunc.yml` | 17 | GUNC 1.1 |
| `anvio9.yml` | 14, 15, 22 to 26 | anvi'o 9 |

Two things these files do not carry.

**Reference databases.** GTDB, eggNOG, CheckM2, DRAM2, GUNC and the skani GTDB
sketches are downloaded separately and are large. `scripts/utils/` holds a
download script for each, and the paths they land in go into
`configs/server.env`. Which release you use matters: the trees here were built
against GTDB r220 and the ANI against r226, and mixing them silently produces
inconsistent results.

**Exact versions.** These are environment specifications, not lockfiles, so
`conda env create` resolves current versions rather than the ones these analyses
ran on. The versions actually used are listed in the Requirements table in the
top-level README. For a byte-identical rebuild, pin them.

anvi'o 9 will not build on every machine. It failed on the first cluster used
here over its `datrie` dependency, which is why the pangenome steps, 14 and 22
through 26, ran somewhere else and read `configs/server_ixachi.env` rather than
`configs/server.env`.
