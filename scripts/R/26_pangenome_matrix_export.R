# Turns the anvi-summarize output from pipeline step 23 into the flat tables the
# rest of the pangenome analysis reads. One set per species:
#
#   {Sp}_GCxMAG_binary.tsv   gene cluster by MAG, presence and absence
#   {Sp}_GC_functions.tsv    consensus function per cluster (KO, COG, Pfam, KEGG_Module)
#   {Sp}_MAG_metadata.tsv    Station, Depth and ROMZ_flag for that species' MAGs
#
# Scripts 27, 28, 29 and 31 all start from these.
#
# Only species with at least five MAGs are exported: Thalassarchaeum, MGIIb-O1
# and UBA8886. UBA57 and UBA60 have four each, enough to build a pangenome on the
# cluster but not enough for the statistics downstream, so they stop here.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(stringr)
})

# Shared style and path resolution. This block locates _pub_style.R without
# resolve_project_dir(), since that function lives inside it.
.msc_root <- Sys.getenv("MSC_ARCHAEA_DIR")
if (!nzchar(.msc_root) || !dir.exists(.msc_root)) {
  # The UNC entry is for RStudio on Windows reading a WSL tree; edit the distro
  # and user name to match yours, or just set MSC_ARCHAEA_DIR and skip it.
  .cand <- c("//wsl.localhost/Ubuntu-24.04/home/ruben/msc_archaea",
             path.expand("~/msc_archaea"))
  .hit <- .cand[dir.exists(.cand)]
  if (length(.hit) == 0) stop("Set MSC_ARCHAEA_DIR to the repository root.", call. = FALSE)
  .msc_root <- .hit[1]
}
source(file.path(.msc_root, "scripts", "R", "_pub_style.R"))

PROJECT_DIR <- resolve_project_dir()
PAN_DIR     <- file.path(PROJECT_DIR, "results", "14_pangenome")
PAN_WORK    <- file.path(PAN_DIR, "pangenome")
SO_03_DIR   <- file.path(PROJECT_DIR, "results", "SO_03")
DERIVED_DIR <- file.path(SO_03_DIR, "derived")
META_TSV    <- file.path(PROJECT_DIR, "MAGs_metadata.tsv")

dir.create(DERIVED_DIR, showWarnings = FALSE, recursive = TRUE)

# The species that carry enough MAGs to report on.
SPECIES_N5PLUS <- c("MGIIb_O1_sp002498525",
                    "Thalassarchaeum_sp002507125",
                    "UBA8886_sp003193815")

# ---- helpers ----
log_info  <- function(...) cat(sprintf("[%s] [INFO] %s\n",
                                       format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                       paste0(...)))
log_warn  <- function(...) cat(sprintf("[%s] [WARN] %s\n",
                                       format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                       paste0(...)), file = stderr())

# Must match what pipeline steps 14 and 24 do to MAG names: every character
# outside [a-zA-Z0-9_] becomes an underscore. Diverge here and the join against
# the pangenome silently matches nothing.
sanitize_genome_name <- function(x) {
  gsub("[^a-zA-Z0-9_]", "_", x)
}

# The consensus annotation for a gene cluster is whichever value occurs most
# often across its genes, ignoring blanks. All blank gives NA.
mode_str <- function(x) {
  x <- x[!is.na(x) & nchar(trimws(as.character(x))) > 0]
  if (length(x) == 0) return(NA_character_)
  tab <- sort(table(x), decreasing = TRUE)
  names(tab)[1]
}

# ---- metadata ----
log_info("Loading MAGs_metadata.tsv")
meta <- read_tsv(META_TSV, show_col_types = FALSE) %>%
  mutate(
    genome     = sanitize_genome_name(NameWithoutExtension),
    Station    = Station,
    Depth      = as.character(Depth),
    ROMZ_flag  = Depth %in% c("ROMZ", "OMZ"),
    Order      = str_extract(FullTaxonomy, "(?<=o__)[^;]+")
  ) %>%
  select(genome, Station, Depth, ROMZ_flag, Order, Genus, Species, DisplayName)

# ---- one species at a time ----
species_dirs <- list.dirs(PAN_WORK, recursive = FALSE)
if (length(species_dirs) == 0) {
  stop("No species directories found in ", PAN_WORK)
}

for (sp_dir in species_dirs) {
  sp_name <- basename(sp_dir)

  if (!(sp_name %in% SPECIES_N5PLUS)) {
    log_info(sprintf("Skipping %s (n < 5 MAGs; fuera del cohorte canónico)", sp_name))
    next
  }

  summary_dir <- file.path(sp_dir, "SUMMARY")
  gc_summary  <- file.path(summary_dir, paste0(sp_name, "_gene_clusters_summary.txt.gz"))

  if (!file.exists(gc_summary)) {
    log_warn(sprintf("Missing gene_clusters_summary for %s; run script 23 on HPC first", sp_name))
    next
  }
  log_info(sprintf("Processing species: %s", sp_name))

  # gene_clusters_summary.txt.gz is tab delimited, one row per gene call
  df <- read_tsv(gc_summary, show_col_types = FALSE, progress = FALSE)

  # The annotation sources anvi'o 9 writes. Missing ones are skipped rather than
  # treated as an error, since not every pangenome carries all of them.
  fn_cols <- intersect(
    c("KOfam_ACC", "KOfam", "KEGG_Module_ACC", "KEGG_Module",
      "KEGG_Class", "COG20_FUNCTION_ACC", "COG20_FUNCTION",
      "COG20_CATEGORY_ACC", "COG20_CATEGORY", "Pfam_ACC", "Pfam"),
    colnames(df)
  )

  # ---- gene cluster by MAG matrix ----
  bin_df <- df %>%
    distinct(gene_cluster_id, genome_name) %>%
    mutate(present = 1L) %>%
    pivot_wider(
      names_from   = genome_name,
      values_from  = present,
      values_fill  = list(present = 0L)
    ) %>%
    arrange(gene_cluster_id)

  out_bin <- file.path(DERIVED_DIR, paste0(sp_name, "_GCxMAG_binary.tsv"))
  write_tsv(bin_df, out_bin)
  log_info(sprintf("  wrote %s (%d GCs × %d MAGs)",
                   basename(out_bin), nrow(bin_df), ncol(bin_df) - 1))

  # ---- consensus function per cluster ----
  if (length(fn_cols) > 0) {
    fn_df <- df %>%
      group_by(gene_cluster_id) %>%
      summarise(
        n_genes      = n(),
        n_genomes    = n_distinct(genome_name),
        across(all_of(fn_cols), mode_str),
        .groups = "drop"
      ) %>%
      arrange(gene_cluster_id)
    out_fn <- file.path(DERIVED_DIR, paste0(sp_name, "_GC_functions.tsv"))
    write_tsv(fn_df, out_fn)
    log_info(sprintf("  wrote %s", basename(out_fn)))
  } else {
    log_warn("  no functional columns found in gene_clusters_summary; skipping GC_functions")
  }

  # ---- metadata for this species' MAGs ----
  sp_genomes <- setdiff(colnames(bin_df), "gene_cluster_id")
  meta_sp <- meta %>% filter(genome %in% sp_genomes)
  missing <- setdiff(sp_genomes, meta_sp$genome)
  if (length(missing) > 0) {
    log_warn(sprintf("  %d genome(s) not found in MAGs_metadata: %s",
                     length(missing), paste(missing, collapse = ", ")))
  }
  out_meta <- file.path(DERIVED_DIR, paste0(sp_name, "_MAG_metadata.tsv"))
  write_tsv(meta_sp, out_meta)
  log_info(sprintf("  wrote %s (%d MAGs)", basename(out_meta), nrow(meta_sp)))
}

log_info("Done. Derived TSVs in: ", DERIVED_DIR)
