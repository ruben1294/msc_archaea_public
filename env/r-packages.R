# Installs every R package the scripts under scripts/R/ need.
#
#   Rscript env/r-packages.R
#
# Safe to rerun: anything already present is left alone.

cran <- c(
  # data handling
  "dplyr", "tidyr", "tibble", "readr", "readxl", "openxlsx", "reshape2",
  "stringr", "rlang", "xfun", "pacman",

  # figures
  "ggplot2", "ggtext", "ggalluvial", "patchwork", "scales", "RColorBrewer",
  "pheatmap", "circlize", "showtext", "sysfonts", "gridtext", "magick",

  # tables for the thesis appendix
  "flextable", "officer", "webshot2",

  # statistics
  "vegan", "lme4", "ape", "micropan",

  # the DeepL client used for the Spanish figure labels
  "deeplr"
)

# ComplexHeatmap and KEGGREST come from Bioconductor. KEGGREST is optional: the
# annotation figures fall back to bare KEGG IDs without it.
bioc <- c("ComplexHeatmap", "KEGGREST")

repos <- "https://cloud.r-project.org/"

missing_cran <- cran[!vapply(cran, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_cran)) {
  message("Installing from CRAN: ", paste(missing_cran, collapse = ", "))
  install.packages(missing_cran, repos = repos, dependencies = TRUE)
} else {
  message("All CRAN packages already installed.")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = repos)
}

missing_bioc <- bioc[!vapply(bioc, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing_bioc)) {
  message("Installing from Bioconductor: ", paste(missing_bioc, collapse = ", "))
  BiocManager::install(missing_bioc, ask = FALSE, update = FALSE)
} else {
  message("All Bioconductor packages already installed.")
}

# CNPS.cycle is not on CRAN. It provides the Ccyc/Ncyc/Pcyc/Scyc abundance
# functions that scripts/R/CNPS.cycle.R builds on; see the header of that script
# for the two bugs it works around.
if (!requireNamespace("CNPS.cycle", quietly = TRUE)) {
  message("CNPS.cycle is not installed.")
  message("  Obtain it from the authors of Yue et al. mSystems 2025,")
  message("  doi:10.1128/msystems.01021-25, then install with")
  message("  install.packages('CNPS.cycle_1.0.1.tar.gz', repos = NULL, type = 'source')")
}

message("Done.")
