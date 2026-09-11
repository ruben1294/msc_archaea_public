#!/usr/bin/env bash
set -euo pipefail

# Sets up anvi'o 9 on the IXACHI cluster.
#
# The pangenome steps moved to this machine because anvi'o 9 would not build on
# the first one, where the datrie dependency fails. env/anvio9.yml holds the
# environment specification this installs from.
module purge
module load conda/conda-2026
source /data/apps/conda/2026/miniconda3/etc/profile.d/conda.sh
conda activate anvio-9

# to setup NCBI’s COG database for quick annotation of genes with functions
# https://anvio.org/help/main/programs/anvi-setup-ncbi-cogs/
anvi-setup-ncbi-cogs

# so anvi-estimate-metabolism and anvi-reaction-network find the database of KEGG orthologs ready when you need it.
# https://anvio.org/help/main/programs/anvi-setup-kegg-data/
anvi-setup-kegg-data

# That's all for now.

