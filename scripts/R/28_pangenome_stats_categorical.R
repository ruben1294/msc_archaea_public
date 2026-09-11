# Tests the gene cluster matrix against the categorical groupings: Depth,
# Station and ROMZ_flag.
#
#   PERMANOVA   vegan::adonis2 over Jaccard ~ Depth + Station, by = "margin"
#   Fisher      exact test per accessory cluster against Depth, BH corrected
#   NMDS        Jaccard ordination
#
# Count-based methods like DESeq2 or CPM do not apply to a pangenome: it records
# whether a cluster occurs in a genome, not how many reads mapped to it. These
# three work on presence and absence, which is what the data actually is.
#
# Figure labels are Spanish, matching the thesis: ROMZ becomes ZMOR, 1000 becomes
# "1000 m", BTM becomes Fondo.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggtext)
  library(stringr)
  library(vegan)
  library(patchwork)
})

PROJECT_DIR <- resolve_project_dir()
SO_03_DIR   <- file.path(PROJECT_DIR, "results", "SO_03")
DERIVED_DIR <- file.path(SO_03_DIR, "derived")
CLASS_DIR   <- file.path(SO_03_DIR, "classification")
OUTPUT_DIR  <- file.path(SO_03_DIR, "stats_categorical")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Shared fonts, theme and palettes
SCRIPTS_R_DIR <- file.path(PROJECT_DIR, "scripts", "R")
source(file.path(SCRIPTS_R_DIR, "_pub_style.R"))

# Depth labels as they appear in the figures
depth_es <- c("ROMZ" = "ZMOR", "1000" = "1000 m", "BTM" = "Fondo")
depth_levels <- c("ZMOR", "1000 m", "Fondo")

# Species name with the genus italicised
species_pretty <- c(
  Thalassarchaeum_sp002507125 = "*Thalassarchaeum* sp002507125",
  MGIIb_O1_sp002498525        = "MGIIb-O1 sp002498525",
  UBA8886_sp003193815         = "UBA8886 sp003193815"
)

# ---- helpers ----
log_info <- function(...) cat(sprintf("[%s] [INFO] %s\n",
                                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                      paste0(...)))
log_warn <- function(...) cat(sprintf("[%s] [WARN] %s\n",
                                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                      paste0(...)), file = stderr())

# ---- one species at a time ----
canonical_species <- names(species_pretty)
bin_files <- list.files(DERIVED_DIR, pattern = "_GCxMAG_binary\\.tsv$", full.names = TRUE)
if (length(bin_files) == 0) stop("No binary matrices in ", DERIVED_DIR, ". Run script 26 first.")

sp_names_all <- str_replace(basename(bin_files), "_GCxMAG_binary\\.tsv$", "")
keep_idx     <- sp_names_all %in% canonical_species
if (any(!keep_idx)) {
  log_info(sprintf("  skipping non-canonical species files: %s",
                   paste(sp_names_all[!keep_idx], collapse = ", ")))
}
bin_files <- bin_files[keep_idx]

panel_nmds <- list()  # subplots para Fig_panel_nmds

for (bin_file in bin_files) {
  sp_name <- str_replace(basename(bin_file), "_GCxMAG_binary\\.tsv$", "")
  pretty  <- species_pretty[sp_name]
  if (is.na(pretty)) pretty <- sp_name
  log_info(sprintf("Species: %s", sp_name))

  bin_df <- read_tsv(bin_file, show_col_types = FALSE)
  meta_file <- file.path(DERIVED_DIR, paste0(sp_name, "_MAG_metadata.tsv"))
  if (!file.exists(meta_file)) { log_warn("  no MAG metadata; skipping"); next }
  meta_sp <- read_tsv(meta_file, show_col_types = FALSE)

  bin_mat <- as.matrix(bin_df[, -1])
  rownames(bin_mat) <- bin_df$gene_cluster_id
  # vegan wants samples in rows, so MAGs are rows and gene clusters are columns
  mag_x_gc <- t(bin_mat)

  # Metadata has to follow the row order exactly
  meta_aligned <- meta_sp[match(rownames(mag_x_gc), meta_sp$genome), ]
  if (any(is.na(meta_aligned$Station))) {
    miss <- rownames(mag_x_gc)[is.na(meta_aligned$Station)]
    log_warn("  MAGs without metadata (skipping them): ",
             paste(miss, collapse = ", "))
    keep <- !is.na(meta_aligned$Station)
    mag_x_gc <- mag_x_gc[keep, , drop = FALSE]
    meta_aligned <- meta_aligned[keep, ]
  }

  if (nrow(mag_x_gc) < 4) {
    log_warn(sprintf("  n=%d MAGs after metadata alignment; skipping stats", nrow(mag_x_gc)))
    next
  }

  # A cluster present in every MAG, or in none, carries no information for any
  # of these tests, so it is dropped before they run.
  pres <- colSums(mag_x_gc)
  variable_gcs <- pres > 0 & pres < nrow(mag_x_gc)
  mag_x_gc_var <- mag_x_gc[, variable_gcs, drop = FALSE]
  if (ncol(mag_x_gc_var) < 5) {
    log_warn(sprintf("  only %d variable GCs; multivariate tests unreliable", ncol(mag_x_gc_var)))
  }

  # Jaccard, the appropriate distance for presence and absence
  jac <- vegdist(mag_x_gc_var, method = "jaccard", binary = TRUE)

  # ---- PERMANOVA ----
  permanova_df <- tibble()
  perm_result <- tryCatch({
    n_levels_dep  <- length(unique(meta_aligned$Depth))
    n_levels_stn  <- length(unique(meta_aligned$Station))
    formula_terms <- character(0)
    if (n_levels_dep >= 2) formula_terms <- c(formula_terms, "Depth")
    if (n_levels_stn >= 2) formula_terms <- c(formula_terms, "Station")
    if (length(formula_terms) == 0) stop("no factors with >=2 levels")
    f <- as.formula(paste("jac ~", paste(formula_terms, collapse = " + ")))
    res <- adonis2(f, data = meta_aligned, by = "margin", permutations = 999)
    as_tibble(res, rownames = "term") %>% mutate(species = sp_name)
  }, error = function(e) {
    log_warn("  PERMANOVA failed: ", e$message); tibble()
  })
  if (nrow(perm_result) > 0) {
    write_tsv(perm_result,
              file.path(OUTPUT_DIR, paste0("permanova_", sp_name, ".tsv")))
  }

  # ---- Fisher exact per accessory cluster against Depth ----
  cls_file <- file.path(CLASS_DIR, paste0("pangenome_classes_", sp_name, ".tsv"))
  if (file.exists(cls_file)) {
    cls <- read_tsv(cls_file, show_col_types = FALSE)
    accessory_gcs <- cls %>%
      filter(class %in% c("shell", "cloud")) %>%
      pull(gene_cluster_id)
    if (length(accessory_gcs) > 0 && length(unique(meta_aligned$Depth)) >= 2) {
      log_info(sprintf("  Fisher exact GC × Depth: testing %d accessory GCs",
                       length(accessory_gcs)))
      fisher_rows <- lapply(accessory_gcs, function(gc) {
        if (!gc %in% colnames(mag_x_gc)) return(NULL)
        tab <- table(meta_aligned$Depth, mag_x_gc[, gc])
        if (any(dim(tab) < 2)) return(NULL)
        ft <- tryCatch(fisher.test(tab, simulate.p.value = TRUE, B = 5000),
                       error = function(e) NULL)
        if (is.null(ft)) return(NULL)
        tibble(gene_cluster_id = gc, p = ft$p.value)
      })
      fisher_df <- bind_rows(fisher_rows) %>%
        mutate(p_BH = p.adjust(p, method = "BH")) %>%
        arrange(p_BH)
      # Attach the consensus annotation where there is one
      fn_file <- file.path(DERIVED_DIR, paste0(sp_name, "_GC_functions.tsv"))
      if (file.exists(fn_file)) {
        fisher_df <- fisher_df %>%
          left_join(read_tsv(fn_file, show_col_types = FALSE),
                    by = "gene_cluster_id")
      }
      write_tsv(fisher_df,
                file.path(OUTPUT_DIR, paste0("fisher_GCxDepth_", sp_name, ".tsv")))
    }
  } else {
    log_warn("  no classification TSV; skipping Fisher tests")
  }

  # ---- NMDS ----
  nmds_res <- tryCatch(
    metaMDS(jac, k = 2, trymax = 50, trace = FALSE),
    error = function(e) { log_warn("  NMDS failed: ", e$message); NULL }
  )
  if (!is.null(nmds_res)) {
    pts <- as_tibble(nmds_res$points, rownames = "genome") %>%
      left_join(meta_aligned, by = "genome") %>%
      mutate(Depth_es = factor(depth_es[Depth], levels = depth_levels),
             Station  = factor(Station, levels = station_levels))
    # The shape scale is fixed to the full station list rather than to the
    # stations this species happens to cover. Without that, each panel builds a
    # different scale and plot_layout(guides = "collect") cannot merge them into
    # one legend.
    station_shapes <- setNames(c(16, 17, 15, 18, 8, 4)[seq_along(station_levels)],
                               station_levels)
    fig_nmds <- ggplot(pts, aes(MDS1, MDS2, color = Depth_es, shape = Station)) +
      geom_hline(yintercept = 0, linetype = "dashed",
                 color = "black", linewidth = 0.4) +
      geom_vline(xintercept = 0, linetype = "dashed",
                 color = "black", linewidth = 0.4) +
      geom_point(size = 3, alpha = 0.85) +
      scale_color_manual(values = depth_colors_es, drop = FALSE) +
      scale_shape_manual(values = station_shapes, drop = FALSE) +
      labs(title = pretty,
           subtitle = sprintf("NMDS Jaccard | *stress* = %.3f", nmds_res$stress),
           x = "NMDS1", y = "NMDS2",
           color = "Profundidad", shape = "Estación") +
      theme_publication(base_size = 12)
    ggsave(file.path(OUTPUT_DIR, paste0("nmds_", sp_name, ".png")),
           fig_nmds, width = 6, height = 4.5, dpi = 300)
    panel_nmds[[sp_name]] <- fig_nmds
  }
}

# ---- NMDS panel across species ----
#
# plot_layout(guides = "collect") still fails to merge the legends when the
# scales observed per panel differ
# (each species covers a different subset of stations, and the internal grob
# comparison treats those guides as unlike each other even with drop = FALSE).
#
# So the legend is suppressed on every panel and one master legend is built
# separately, from a throwaway ggplot covering all six stations and all three
# depths. That legend goes in as a strip along the bottom, which is what
# guarantees exactly one box for depth and one for station.
build_nmds_panel <- function(plot_list, title_main) {
  ncol_grid <- if (length(plot_list) <= 3) length(plot_list) else 3

  station_shapes <- setNames(c(16, 17, 15, 18, 8, 4)[seq_along(station_levels)],
                             station_levels)

  legend_df <- expand.grid(
    Station  = factor(station_levels, levels = station_levels),
    Depth_es = factor(depth_levels,  levels = depth_levels),
    stringsAsFactors = FALSE
  )
  legend_p <- ggplot2::ggplot(
      legend_df,
      ggplot2::aes(x = 1, y = 1, color = Depth_es, shape = Station)
    ) +
    ggplot2::geom_point(size = 4, alpha = 0.9) +
    ggplot2::scale_color_manual(
      values = depth_colors_es, name = "Profundidad", drop = FALSE
    ) +
    ggplot2::scale_shape_manual(
      values = station_shapes, name = "Estación", drop = FALSE
    ) +
    ggplot2::guides(
      color = ggplot2::guide_legend(
        override.aes = list(shape = 16, size = 5), order = 1
      ),
      shape = ggplot2::guide_legend(
        override.aes = list(size = 4, color = "grey30"), order = 2
      )
    ) +
    ggplot2::theme_void(base_size = 12) +
    ggplot2::theme(
      legend.position  = "bottom",
      legend.direction = "horizontal",
      legend.box       = "horizontal",
      legend.title     = ggplot2::element_text(face = "bold", size = 12),
      legend.text      = ggplot2::element_text(size = 11),
      plot.margin      = grid::unit(c(0, 0, 0, 0), "lines")
    )

  panels <- patchwork::wrap_plots(plot_list, ncol = ncol_grid) &
    ggplot2::theme(legend.position = "none")

  pn <- panels / legend_p +
    patchwork::plot_layout(heights = c(1, 0.18)) +
    patchwork::plot_annotation(
      title = title_main,
      theme = ggplot2::theme(
        plot.title = ggtext::element_markdown(
          family = font_family, face = "bold",
          size = 18, hjust = 0.5
        )
      )
    )

  panel_w <- max(7, 5.6 * min(3, length(plot_list)))
  panel_h <- 4.5 * ceiling(length(plot_list) / ncol_grid) + 1.2
  list(plot = pn, w = panel_w, h = panel_h)
}

if (length(panel_nmds) > 0) {
  ord <- intersect(canonical_species, names(panel_nmds))
  panel_nmds <- panel_nmds[ord]
  pn <- build_nmds_panel(panel_nmds,
                         "NMDS Jaccard del genoma accesorio por especie")
  ggsave(file.path(OUTPUT_DIR, "Fig_panel_nmds.png"),
         pn$plot, width = pn$w, height = pn$h,
         dpi = 300, limitsize = FALSE)
  log_info("Panel NMDS written")
}

log_info("Done. Outputs in: ", OUTPUT_DIR)
