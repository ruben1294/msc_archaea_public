# Asks whether where a MAG came from, its depth or station, or the water
# chemistry there, explains differences in gene cluster homogeneity and AAI, and
# in which shell and cloud clusters a genome carries. High-pressure adaptation is
# the case this was built to test.
#
# Five analyses, all per species except the ANI one, which is global:
#
#   A  Spearman, per-cluster AAI and homogeneity against depth entropy, where
#      entropy measures how confined to one depth a cluster is. Boxplots by
#      dominant-depth category.
#   B  Mantel, ANI distance between MAGs against depth (encoded ordinally),
#      station (Gower) and environment (Euclidean).
#   C  Fisher exact, cluster against station, which complements the cluster
#      against depth test in script 28. BH corrected, with the top clusters
#      ranked.
#   D  Spearman, cluster presence against Pressure_atm and O2_dissolved.
#      Logistic regression degenerates at these sample sizes, so this uses
#      rank correlation. BH corrected per variable.
#   E  Candidate list: clusters passing either test at BH < 0.10, annotated
#      against keywords for pressure, osmotic, transport, lipid and heat shock
#      adaptation.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(readxl)
  library(stringr)
  library(ggplot2)
  library(ggtext)
  library(patchwork)
  library(vegan)
})

# Null-coalescing operator, defined up front because it is used early
`%||%` <- function(a, b) if (is.null(a) || (length(a) == 1 && is.na(a))) b else a

# --- Paths --------------------------------------------------------------------
# Repository root: MSC_ARCHAEA_DIR when set, otherwise the nearest directory
# above this script (under Rscript) or above the working directory (in an
# interactive session) that holds scripts/R/_pub_style.R.
PROJECT_DIR <- local({
  root <- Sys.getenv("MSC_ARCHAEA_DIR")
  if (!nzchar(root)) {
    script <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE))
    root <- normalizePath(if (length(script) == 1) dirname(script) else getwd())
    while (!file.exists(file.path(root, "scripts", "R", "_pub_style.R")) &&
           dirname(root) != root) {
      root <- dirname(root)
    }
    if (!file.exists(file.path(root, "scripts", "R", "_pub_style.R"))) {
      stop("Repository root not found. Set MSC_ARCHAEA_DIR to it.", call. = FALSE)
    }
    Sys.setenv(MSC_ARCHAEA_DIR = root)
  }
  normalizePath(root)
})
PAN_DIR     <- file.path(PROJECT_DIR, "results", "20_pangenome")
PAN_WORK    <- file.path(PAN_DIR, "pangenome")
SO_03_DIR   <- file.path(PROJECT_DIR, "results", "SO_03")
DERIVED_DIR <- file.path(SO_03_DIR, "derived")
CLASS_DIR   <- file.path(SO_03_DIR, "classification")
OUTPUT_DIR  <- file.path(SO_03_DIR, "origin_drivers")
ENV_XLSX    <- file.path(PROJECT_DIR, "data", "Physicochemical_properties.xlsx")
ANI_TRI     <- file.path(PROJECT_DIR, "results", "04_sketching_ani",
                         "all_mags", "ani_triangle_full.tsv")
META_TSV    <- file.path(PROJECT_DIR, "data", "MAGs_metadata.tsv")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

SCRIPTS_R_DIR <- file.path(PROJECT_DIR, "scripts", "R")
source(file.path(SCRIPTS_R_DIR, "_pub_style.R"))

species_pretty <- c(
  Thalassarchaeum_sp002507125 = "*Thalassarchaeum* sp002507125",
  MGIIb_O1_sp002498525        = "MGIIb-O1 sp002498525",
  UBA8886_sp003193815         = "UBA8886 sp003193815"
)

# --- Helpers ------------------------------------------------------------------
log_info <- function(...) cat(sprintf("[%s] [INFO] %s\n",
                                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                      paste0(...)))
log_warn <- function(...) cat(sprintf("[%s] [WARN] %s\n",
                                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                      paste0(...)), file = stderr())
sanitize <- function(x) gsub("[^a-zA-Z0-9_]", "_", x)

# Shannon entropy normalised by log(k): 0 means the cluster is confined to a single
# category, 1 means it spreads evenly across all k. It separates
# depth-restricted clusters from depth-agnostic ones.
shannon_norm <- function(counts) {
  counts <- counts[counts > 0]
  if (length(counts) <= 1) return(0)
  p <- counts / sum(counts)
  -sum(p * log(p)) / log(length(counts))
}

# Case-insensitive keywords suggesting adaptation to depth, pressure or
# oligotrophic water. Curated for marine archaea rather than taken from a
# generic list, so a hit here is worth following up.
adapt_regex <- paste0(
  "(?i)",
  "pressure|barotol|piezo|",
  "trehalose|ectoine|glycine.betaine|compatible.solute|osmoprotect|osmotol|",
  "TMAO|DMSP|DMSO|",
  "membrane.fluid|fatty.acid.desat|cyclopropane|",
  "cold.shock|heat.shock|chaperon|HSP|GroEL|DnaK|",
  "oxidative.stress|peroxiredoxin|catalase|superoxide|",
  "K\\+ transport|sodium.proton|Na\\+/H\\+|",
  "ammonia|amoA|nitrite|nitrate|",
  "cytochrome.c.oxidase|cbb3|aa3"
)

# --- Load metadata and environmental data -------------------------------------
meta <- read_tsv(META_TSV, show_col_types = FALSE) %>%
  mutate(genome = sanitize(NameWithoutExtension)) %>%
  select(genome, NameWithoutExtension, Station, Depth)

env_depth <- read_xlsx(ENV_XLSX, sheet = "Physicochemical_by_depth") %>%
  rename(
    Temperature  = `mean_ConservativeTemp (°C)`,
    Salinity     = `mean_AbsoluteSalinity (g/kg)`,
    O2_dissolved = mean_DissolvedOxigen,
    AOU          = `mean_CAO (umol/kg)`,
    NO3_NO2      = mean_NO3_NO2,
    PO4          = mean_PO4
  ) %>%
  select(Depth, Temperature, Salinity, O2_dissolved, AOU, NO3_NO2, PO4)

env_station <- read_xlsx(ENV_XLSX, sheet = "Stations_depths") %>%
  mutate(Pressure_atm = Collection_depth_m / 10) %>%
  select(Station, Depth, Pressure_atm, Density_kg_m3, Collection_depth_m)

env_full <- meta %>%
  left_join(env_station, by = c("Station", "Depth")) %>%
  left_join(env_depth,   by = "Depth")

# --- Load the ANI triangle (skani, full) --------------------------------------
ani_long <- NULL
if (file.exists(ANI_TRI)) {
  ani_lines <- readLines(ANI_TRI)
  n_genomes_ani <- as.integer(ani_lines[1])
  ani_rows <- strsplit(ani_lines[2:(n_genomes_ani + 1)], "\t")
  ani_paths <- vapply(ani_rows, function(r) r[1], character(1))
  ani_keys  <- sanitize(tools::file_path_sans_ext(basename(ani_paths)))
  ani_mat   <- t(vapply(ani_rows, function(r) suppressWarnings(as.numeric(r[-1])),
                        numeric(n_genomes_ani)))
  rownames(ani_mat) <- ani_keys
  colnames(ani_mat) <- ani_keys
  # skani triangle can come back asymmetric. Taking the maximum of the two
  # directions matches how ANI is normally reported.
  ani_mat[ani_mat == 0] <- NA
  ani_sym <- pmax(ani_mat, t(ani_mat), na.rm = TRUE)
  diag(ani_sym) <- 100
  ani_long <- as_tibble(ani_sym, rownames = "g1") %>%
    pivot_longer(-g1, names_to = "g2", values_to = "ANI")
  log_info(sprintf("ANI triangle loaded: %d MAGs", n_genomes_ani))
} else {
  log_warn("ANI triangle not found: ", ANI_TRI)
}

# ---- one species at a time ----
canonical_species <- names(species_pretty)
cand_all          <- list()
ani_summary       <- list()
panel_homog       <- list()
panel_volcano_pr  <- list()

for (sp_name in canonical_species) {
  log_info(sprintf("Species: %s", sp_name))
  pretty <- species_pretty[sp_name]

  bin_file  <- file.path(DERIVED_DIR, paste0(sp_name, "_GCxMAG_binary.tsv"))
  fn_file   <- file.path(DERIVED_DIR, paste0(sp_name, "_GC_functions.tsv"))
  cls_file  <- file.path(CLASS_DIR,   paste0("pangenome_classes_", sp_name, ".tsv"))
  gcsum_gz  <- file.path(PAN_WORK, sp_name, "SUMMARY",
                         paste0(sp_name, "_gene_clusters_summary.txt.gz"))
  if (!all(file.exists(c(bin_file, cls_file, gcsum_gz)))) {
    log_warn("  inputs missing; skipping"); next
  }

  bin_df  <- read_tsv(bin_file, show_col_types = FALSE)
  cls     <- read_tsv(cls_file, show_col_types = FALSE)
  fn_df   <- if (file.exists(fn_file)) read_tsv(fn_file, show_col_types = FALSE) else NULL
  gc_sum  <- read_tsv(gcsum_gz, show_col_types = FALSE, progress = FALSE)

  bin_mat <- as.matrix(bin_df[, -1])
  rownames(bin_mat) <- bin_df$gene_cluster_id
  sp_genomes <- colnames(bin_mat)
  meta_sp    <- env_full %>% filter(genome %in% sp_genomes)

  # ---- A) per-cluster homogeneity and AAI ----
  # gc_sum holds one row per gene call, so it collapses to one row per cluster
  # by taking the first value. The indices are cluster attributes and repeat
  # identically across every gene call in the cluster.
  gc_metrics <- gc_sum %>%
    group_by(gene_cluster_id) %>%
    summarise(
      AAI_min     = first(AAI_min),
      AAI_max     = first(AAI_max),
      AAI_avg     = first(AAI_avg),
      hom_funct   = first(functional_homogeneity_index),
      hom_geom    = first(geometric_homogeneity_index),
      hom_combined = first(combined_homogeneity_index),
      .groups = "drop"
    )

  # How each cluster distributes across depths
  depth_per_mag <- meta_sp %>% select(genome, Depth)
  depth_levels_ord <- c("ROMZ", "1000", "BTM")
  gc_depth_long <- as_tibble(bin_mat, rownames = "gene_cluster_id") %>%
    pivot_longer(-gene_cluster_id, names_to = "genome", values_to = "present") %>%
    filter(present == 1) %>%
    left_join(depth_per_mag, by = "genome") %>%
    count(gene_cluster_id, Depth, name = "n")
  gc_depth_wide <- gc_depth_long %>%
    pivot_wider(names_from = Depth, values_from = n, values_fill = 0)
  for (d in depth_levels_ord) {
    if (!d %in% colnames(gc_depth_wide)) gc_depth_wide[[d]] <- 0L
  }
  depth_mat <- as.matrix(gc_depth_wide[, depth_levels_ord, drop = FALSE])
  gc_by_depth <- gc_depth_wide %>%
    select(gene_cluster_id) %>%
    mutate(
      depth_entropy = apply(depth_mat, 1, shannon_norm),
      depth_dom     = vapply(seq_len(nrow(depth_mat)), function(i) {
        v <- depth_mat[i, ]
        if (sum(v) == 0) NA_character_
        else if (sum(v > 0) == 1) paste0(depth_levels_ord[which.max(v)], "-only")
        else "mixto"
      }, character(1))
    )

  # Join against the metrics and the pangenome class
  gc_full <- cls %>%
    select(gene_cluster_id, class) %>%
    left_join(gc_metrics, by = "gene_cluster_id") %>%
    left_join(gc_by_depth %>% select(gene_cluster_id, depth_entropy, depth_dom),
              by = "gene_cluster_id")

  # Spanish labels, matching the rest of the figures
  recode_dom <- c("ROMZ-only" = "ZMOR-only", "1000-only" = "1000 m-only",
                  "BTM-only"  = "Fondo-only", "mixto" = "mixto")
  gc_full <- gc_full %>%
    mutate(depth_dom = factor(recode_dom[depth_dom],
                              levels = c("ZMOR-only", "1000 m-only",
                                         "Fondo-only", "mixto")))

  # Spearman, each metric against entropy
  metric_cols <- c("AAI_min", "AAI_max", "AAI_avg",
                   "hom_funct", "hom_geom", "hom_combined")
  spearman_homog <- bind_rows(lapply(metric_cols, function(m) {
    x <- gc_full[[m]]; y <- gc_full$depth_entropy
    keep <- !is.na(x) & !is.na(y) & is.finite(x) & is.finite(y)
    if (sum(keep) < 5) return(tibble(metric = m, n = sum(keep),
                                     rho = NA, p = NA))
    res <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman"))
    tibble(metric = m, n = sum(keep), rho = res$estimate, p = res$p.value)
  })) %>% mutate(species = sp_name, BH = p.adjust(p, method = "BH"))

  write_tsv(spearman_homog,
            file.path(OUTPUT_DIR, paste0("A_spearman_homog_", sp_name, ".tsv")))
  write_tsv(gc_full,
            file.path(OUTPUT_DIR, paste0("A_GC_full_metrics_", sp_name, ".tsv")))

  # Boxplot of hom_combined and AAI_avg by dominant depth
  long_hom <- gc_full %>%
    filter(!is.na(depth_dom)) %>%
    select(class, depth_dom, AAI_avg, hom_funct, hom_geom) %>%
    pivot_longer(c(AAI_avg, hom_funct, hom_geom),
                 names_to = "metric", values_to = "value") %>%
    mutate(metric = recode(metric,
                           AAI_avg   = "AAI promedio",
                           hom_funct = "Homog. funcional",
                           hom_geom  = "Homog. geométrica"))

  # geom_boxplot draws nothing at n = 1, since there is no IQR. The jitter keeps
  # the points visible regardless, and the crossbar puts a median line on the
  # single-observation categories, so every category is drawn as a box.
  long_hom_plot <- long_hom %>%
    filter(!is.na(value), !is.na(depth_dom)) %>%
    group_by(metric, depth_dom) %>%
    mutate(n_per_cat = dplyr::n()) %>%
    ungroup()
  fig_hom <- ggplot(long_hom_plot,
                    aes(depth_dom, value, fill = depth_dom, color = depth_dom)) +
    geom_boxplot(data = filter(long_hom_plot, n_per_cat >= 2),
                 alpha = 0.7, outlier.shape = NA, color = "grey20") +
    stat_summary(data = filter(long_hom_plot, n_per_cat == 1),
                 fun = median, geom = "crossbar",
                 width = 0.55, fatten = 1.2,
                 color = "grey20") +
    geom_jitter(width = 0.15, size = 0.9, alpha = 0.55) +
    facet_wrap(~ metric, scales = "free_y", ncol = 3) +
    scale_fill_manual(values = c("ZMOR-only"   = "#4A90D9",
                                 "1000 m-only" = "#E67E22",
                                 "Fondo-only"  = "#8E44AD",
                                 "mixto"       = "grey60"),
                      drop = FALSE) +
    scale_color_manual(values = c("ZMOR-only"   = "#4A90D9",
                                  "1000 m-only" = "#E67E22",
                                  "Fondo-only"  = "#8E44AD",
                                  "mixto"       = "grey60"),
                       drop = FALSE, guide = "none") +
    labs(title = pretty,
         x = "Distribución del clúster génico",
         y = "Valor del índice AAI (%)",
         fill = NULL) +
    theme_publication(base_size = 12) +
    theme(legend.position = "none",
          axis.text.x = element_text(angle = 30, hjust = 1),
          axis.title.x = element_markdown())
  ggsave(file.path(OUTPUT_DIR, paste0("Fig_homog_by_depthdist_", sp_name, ".png")),
         fig_hom, width = 8, height = 4.5, dpi = 300)
  panel_homog[[sp_name]] <- fig_hom

  # ---- B) per-MAG ANI against depth, station and environment (Mantel) ----
  if (!is.null(ani_long) && nrow(meta_sp) >= 4) {
    sp_keys <- meta_sp$genome
    ani_sp  <- ani_long %>%
      filter(g1 %in% sp_keys, g2 %in% sp_keys) %>%
      pivot_wider(names_from = g2, values_from = ANI) %>%
      column_to_rownames("g1")
    if (all(sp_keys %in% rownames(ani_sp))) {
      ani_sp <- as.matrix(ani_sp)[sp_keys, sp_keys]
      ani_sp[ani_sp <= 0 | is.na(ani_sp)] <- min(ani_sp[ani_sp > 0], na.rm = TRUE)
      d_ani <- as.dist(100 - ani_sp)  # distancia = 100 - ANI

      depth_ord_num <- c("ROMZ" = 1, "1000" = 2, "BTM" = 3)
      d_depth <- dist(depth_ord_num[meta_sp$Depth])
      d_stn   <- vegdist(model.matrix(~ Station - 1, data = meta_sp),
                         method = "jaccard", binary = TRUE)
      env_cols <- c("Pressure_atm", "O2_dissolved", "Temperature",
                    "Salinity", "AOU", "NO3_NO2", "PO4")
      env_cols <- intersect(env_cols, colnames(meta_sp))
      env_mat  <- as.matrix(meta_sp[, env_cols, drop = FALSE])
      d_env <- if (sum(complete.cases(env_mat)) >= 4) {
        vegdist(scale(env_mat[complete.cases(env_mat), , drop = FALSE]),
                method = "euclidean")
      } else NULL

      mant <- list(
        Depth   = tryCatch(mantel(d_ani, d_depth, permutations = 999), error = function(e) NULL),
        Station = tryCatch(mantel(d_ani, d_stn,   permutations = 999), error = function(e) NULL),
        Env     = if (!is.null(d_env) && length(attr(d_env, "Size")) > 0)
                    tryCatch(mantel(as.dist(as.matrix(d_ani)[
                      complete.cases(env_mat), complete.cases(env_mat)]),
                      d_env, permutations = 999), error = function(e) NULL)
                  else NULL
      )
      ani_summary[[sp_name]] <- bind_rows(lapply(names(mant), function(k) {
        m <- mant[[k]]
        tibble(species = sp_name, predictor = k,
               r = if (!is.null(m)) m$statistic else NA,
               p = if (!is.null(m)) m$signif    else NA,
               n = nrow(meta_sp))
      }))

      # ANI heatmap, ordered by depth
      ani_df <- as_tibble(ani_sp, rownames = "g1") %>%
        pivot_longer(-g1, names_to = "g2", values_to = "ANI") %>%
        left_join(meta_sp %>% select(genome, Depth, Station),
                  by = c("g1" = "genome")) %>%
        rename(Depth1 = Depth, Station1 = Station) %>%
        left_join(meta_sp %>% select(genome, Depth, Station),
                  by = c("g2" = "genome")) %>%
        rename(Depth2 = Depth, Station2 = Station)
      mag_order <- meta_sp %>%
        arrange(factor(Depth, levels = c("ROMZ", "1000", "BTM")), Station) %>%
        pull(genome)
      ani_df <- ani_df %>%
        mutate(g1 = factor(g1, levels = mag_order),
               g2 = factor(g2, levels = mag_order))
      fig_ani <- ggplot(ani_df, aes(g1, g2, fill = ANI)) +
        geom_tile() +
        scale_fill_viridis_c(option = "C", na.value = "white",
                             limits = c(80, 100), oob = scales::squish) +
        labs(title = pretty,
             subtitle = sprintf("ANI per-MAG | n = %d", nrow(meta_sp)),
             x = NULL, y = NULL, fill = "ANI") +
        theme_publication(base_size = 11) +
        theme(axis.text.x = element_text(angle = 60, hjust = 1, size = 8),
              axis.text.y = element_text(size = 8),
              panel.grid = element_blank())
      ggsave(file.path(OUTPUT_DIR, paste0("Fig_ANI_heatmap_", sp_name, ".png")),
             fig_ani, width = 6, height = 5, dpi = 300)
    }
  }

  # ---- C) Fisher, cluster against station (complements script 28) ----
  accessory <- cls %>% filter(class %in% c("shell", "cloud")) %>% pull(gene_cluster_id)
  fisher_stn_df <- tibble()
  if (length(accessory) > 0 && length(unique(meta_sp$Station)) >= 2) {
    rows <- lapply(accessory, function(gc) {
      if (!gc %in% rownames(bin_mat)) return(NULL)
      tab <- table(meta_sp$Station, bin_mat[gc, meta_sp$genome])
      if (any(dim(tab) < 2)) return(NULL)
      ft <- tryCatch(fisher.test(tab, simulate.p.value = TRUE, B = 5000),
                     error = function(e) NULL)
      if (is.null(ft)) return(NULL)
      tibble(gene_cluster_id = gc, p_station = ft$p.value)
    })
    fisher_stn_df <- bind_rows(rows) %>%
      mutate(BH_station = p.adjust(p_station, method = "BH"))
    if (!is.null(fn_df)) {
      fisher_stn_df <- fisher_stn_df %>%
        left_join(fn_df, by = "gene_cluster_id")
    }
    write_tsv(fisher_stn_df,
              file.path(OUTPUT_DIR, paste0("C_fisher_GCxStation_", sp_name, ".tsv")))
  }

  # --- D) Spearman GC presence vs env ---------------------------------------
  env_pred <- c("Pressure_atm", "O2_dissolved", "NO3_NO2", "Temperature")
  env_pred <- intersect(env_pred, colnames(meta_sp))
  spear_env <- list()
  for (v in env_pred) {
    x <- meta_sp[[v]]
    if (sum(!is.na(x)) < 4 || length(unique(x[!is.na(x)])) < 2) next
    rows <- lapply(accessory, function(gc) {
      if (!gc %in% rownames(bin_mat)) return(NULL)
      y <- bin_mat[gc, meta_sp$genome]
      keep <- !is.na(x) & !is.na(y)
      if (sum(keep) < 4 || length(unique(y[keep])) < 2) return(NULL)
      res <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman"))
      tibble(gene_cluster_id = gc, env_var = v,
             rho = unname(res$estimate), p = res$p.value)
    })
    df <- bind_rows(rows)
    if (nrow(df) > 0) df$BH <- p.adjust(df$p, method = "BH")
    spear_env[[v]] <- df
  }
  spear_env_df <- bind_rows(spear_env)
  if (nrow(spear_env_df) > 0) {
    if (!is.null(fn_df)) {
      spear_env_df <- spear_env_df %>%
        left_join(fn_df, by = "gene_cluster_id")
    }
    write_tsv(spear_env_df,
              file.path(OUTPUT_DIR, paste0("D_spearman_env_", sp_name, ".tsv")))

    # Volcano shape: |rho| against -log10(p) per environmental variable, with
    # the adaptive hits labelled, meaning the KOfam or COG24 annotations that
    # match adapt_regex.
    label_col <- intersect(c("KOfam", "COG24_FUNCTION", "COG20_FUNCTION", "Pfam"),
                           colnames(spear_env_df))[1]
    spear_env_df <- spear_env_df %>%
      mutate(adaptive = if (!is.na(label_col)) {
        grepl(adapt_regex, .data[[label_col]], perl = TRUE)
      } else FALSE)

    env_labels <- c(
      Pressure_atm = "Presión (atm)",
      O2_dissolved = "O<sub>2</sub> disuelto",
      NO3_NO2      = "NO<sub>3</sub><sup>−</sup>+NO<sub>2</sub><sup>−</sup>",
      Temperature  = "Temperatura"
    )
    # shape = 16, a filled circle with no border. The bordered default rendered
    # as two greys and two reds, since border and fill picked up the scale
    # separately.
    fig_volc <- ggplot(spear_env_df, aes(rho, -log10(p))) +
      geom_point(aes(color = adaptive, size = adaptive),
                 shape = 16, alpha = 0.85) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed",
                 color = "grey30") +
      facet_wrap(~ env_var, scales = "free", ncol = 2,
                 labeller = as_labeller(env_labels)) +
      scale_color_manual(values = c(`TRUE` = "#D7191C", `FALSE` = "grey45"),
                         labels = c(`TRUE` = "Candidato adaptativo",
                                    `FALSE` = "Resto"),
                         name = NULL) +
      scale_size_manual(values = c(`TRUE` = 2.0, `FALSE` = 1.0),
                        guide = "none") +
      guides(color = guide_legend(
        override.aes = list(size = 6, alpha = 1)
      )) +
      # The rho glyph does not render in this environment: the active font
      # family has no Greek glyphs, and the fallback substitutes a question mark.
      # The axis names the statistic instead; the scale is -1 to +1 either way.
      labs(title = pretty,
           x = "Coeficiente de Spearman",
           y = "-log<sub>10</sub>(*p*)") +
      theme_publication(base_size = 11) +
      theme(legend.position = "bottom",
            axis.title.y  = element_markdown(),
            strip.text    = element_markdown())
    ggsave(file.path(OUTPUT_DIR, paste0("Fig_volcano_env_", sp_name, ".png")),
           fig_volc, width = 7, height = 5.5, dpi = 300)
    panel_volcano_pr[[sp_name]] <- fig_volc
  }

  # --- E) Combined table of candidates --------------------------------------
  fisher_dep_file <- file.path(SO_03_DIR, "stats_categorical",
                               paste0("fisher_GCxDepth_", sp_name, ".tsv"))
  fisher_dep <- if (file.exists(fisher_dep_file)) {
    read_tsv(fisher_dep_file, show_col_types = FALSE) %>%
      select(gene_cluster_id, p_depth = p, BH_depth = p_BH)
  } else tibble()

  candidates <- cls %>%
    filter(class %in% c("shell", "cloud")) %>%
    select(gene_cluster_id, class, n_genomes_present, prevalence)
  if (nrow(fisher_dep) > 0)
    candidates <- candidates %>% left_join(fisher_dep, by = "gene_cluster_id")
  if (nrow(fisher_stn_df) > 0)
    candidates <- candidates %>%
      left_join(fisher_stn_df %>% select(gene_cluster_id, p_station, BH_station),
                by = "gene_cluster_id")
  if (nrow(spear_env_df) > 0) {
    spear_wide <- spear_env_df %>%
      select(gene_cluster_id, env_var, rho, BH) %>%
      pivot_wider(names_from = env_var, values_from = c(rho, BH),
                  names_glue = "{env_var}_{.value}")
    candidates <- candidates %>% left_join(spear_wide, by = "gene_cluster_id")
  }
  if (!is.null(fn_df)) candidates <- candidates %>% left_join(fn_df, by = "gene_cluster_id")

  fn_text_cols <- intersect(c("KOfam", "COG24_FUNCTION", "COG20_FUNCTION", "Pfam"),
                            colnames(candidates))
  if (length(fn_text_cols) > 0) {
    fn_concat <- apply(candidates[, fn_text_cols, drop = FALSE], 1,
                       function(r) paste(r[!is.na(r)], collapse = " | "))
    candidates$adaptive_keyword <- grepl(adapt_regex, fn_concat, perl = TRUE)
  } else {
    candidates$adaptive_keyword <- FALSE
  }
  candidates <- candidates %>%
    mutate(
      sig_origin = (!is.na(BH_depth)   & BH_depth   < 0.10) |
                   (!is.na(BH_station) & BH_station < 0.10),
      sig_env = if (any(grepl("_BH$", colnames(.))))
                  rowSums(across(matches("_BH$"), ~ !is.na(.) & . < 0.10)) > 0
                else FALSE,
      candidate_priority = sig_origin | sig_env | adaptive_keyword
    ) %>%
    arrange(desc(candidate_priority),
            pmin(BH_depth, BH_station, na.rm = TRUE))

  cand_all[[sp_name]] <- candidates %>% mutate(species = sp_name)
  write_tsv(candidates,
            file.path(OUTPUT_DIR, paste0("E_candidates_integrated_", sp_name, ".tsv")))

  log_info(sprintf("  candidates with priority flag: %d / %d accessory GCs",
                   sum(candidates$candidate_priority, na.rm = TRUE),
                   nrow(candidates)))
}

# ---- summaries and panels across species ----
if (length(ani_summary) > 0) {
  write_tsv(bind_rows(ani_summary),
            file.path(OUTPUT_DIR, "B_mantel_ANI_summary.tsv"))
}
if (length(cand_all) > 0) {
  write_tsv(bind_rows(cand_all),
            file.path(OUTPUT_DIR, "E_candidates_all_species.tsv"))
}

build_panel <- function(plot_list, title_main, ncol_max = 3) {
  if (length(plot_list) == 0) return(NULL)
  ord <- intersect(canonical_species, names(plot_list))
  plot_list <- plot_list[ord]
  ncol_grid <- if (length(plot_list) <= ncol_max) length(plot_list) else ncol_max
  patchwork::wrap_plots(plot_list, ncol = ncol_grid) +
    patchwork::plot_annotation(
      title = title_main,
      theme = ggplot2::theme(
        plot.title = ggtext::element_markdown(
          family = font_family, face = "bold", size = 18, hjust = 0.5
        )
      )
    )
}

if (length(panel_homog) > 0) {
  pp <- build_panel(panel_homog,
                    "Homogeneidad y AAI por distribución de profundidad",
                    ncol_max = 2)
  ggsave(file.path(OUTPUT_DIR, "Fig_panel_homog_by_depthdist.png"),
         pp,
         width  = 8 * min(2, length(panel_homog)),
         height = 4.5 * ceiling(length(panel_homog) / 2),
         dpi = 300, limitsize = FALSE)
  log_info("Panel A (homogeneidad) written")
}

if (length(panel_volcano_pr) > 0) {
  ord <- intersect(canonical_species, names(panel_volcano_pr))
  plot_list <- panel_volcano_pr[ord]
  n_pl      <- length(plot_list)
  ncol_grid <- if (n_pl <= 3) n_pl else 3

  # When the panels do not fill the grid, for instance three panels in a 2x2,
  # the spare cell becomes a guide_area and holds the legend at a larger size.
  # Otherwise the legend goes underneath as usual.
  n_rows         <- ceiling(n_pl / ncol_grid)
  n_cells        <- n_rows * ncol_grid
  has_guide_area <- n_cells > n_pl

  if (has_guide_area) {
    pp <- patchwork::wrap_plots(
      c(plot_list, list(patchwork::guide_area())),
      ncol = ncol_grid
    )
  } else {
    pp <- patchwork::wrap_plots(plot_list, ncol = ncol_grid)
  }

  pp <- pp +
    patchwork::plot_annotation(
      title = "Spearman: presencia del clúster génico vs variables fisicoquímicas",
      theme = ggplot2::theme(
        plot.title = ggtext::element_markdown(
          family = font_family, face = "bold", size = 18, hjust = 0.5
        )
      )
    ) +
    patchwork::plot_layout(guides = "collect")

  if (has_guide_area) {
    pp <- pp & ggplot2::theme(
      legend.position  = "right",
      legend.direction = "vertical",
      legend.text      = ggplot2::element_text(size = 17),
      legend.title     = ggplot2::element_text(size = 17, face = "bold"),
      legend.key.size  = grid::unit(2.2, "lines"),
      legend.spacing.y = grid::unit(0.6, "lines")
    )
  } else {
    pp <- pp & ggplot2::theme(
      legend.position = "bottom",
      legend.text     = ggplot2::element_text(size = 17),
      legend.title    = ggplot2::element_text(size = 17, face = "bold"),
      legend.key.size = grid::unit(2.2, "lines")
    )
  }

  ggsave(file.path(OUTPUT_DIR, "Fig_panel_volcano_env.png"),
         pp,
         width  = 7 * min(ncol_grid, n_pl),
         height = 5.5 * n_rows,
         dpi = 300, limitsize = FALSE)
  log_info("Panel D (volcano env) written")
}

log_info("Done. Outputs in: ", OUTPUT_DIR)
