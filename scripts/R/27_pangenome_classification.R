# Sorts each species' gene clusters into core, soft-core, shell and cloud, then
# asks through Heaps' law whether the pangenome is open or closed.
#
# The bands, by fraction of MAGs carrying the cluster:
#   core       100%
#   soft-core  at least 90%, which collapses into core when n <= 7
#   shell      15% up to 90%
#   cloud      under 15%, singletons included
#
# Two results here were expected from the annotation work and serve as a check
# that the pangenome is keyed correctly:
#   cbb3 (K00404 to K00407) should land in the core of MGIIb_O1_sp002498525
#   HndC (K18006) should land in core or soft-core of UBA8886_sp003193815
#
# A cluster missing from an incomplete MAG looks the same as one that genuinely
# is not there, which biases everything toward the accessory bands. The
# sensitivity pass at the end reclassifies using only MAGs above 80% CheckM2
# completeness and reports how far the two agree.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggplot2)
  library(ggtext)
  library(stringr)
  library(patchwork)
})

# micropan is optional; Heaps' law is skipped without it
has_micropan <- requireNamespace("micropan", quietly = TRUE)
if (!has_micropan) {
  message("micropan no instalado; Heaps' law no se calculará. Instala con: install.packages('micropan')")
}

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
SO_03_DIR <- file.path(PROJECT_DIR, "results", "SO_03")
DERIVED_DIR <- file.path(SO_03_DIR, "derived")
OUTPUT_DIR <- file.path(SO_03_DIR, "classification")
QUALITY_TSV <- file.path(
  PROJECT_DIR, "results", "01_mag_quality", "quality_report.tsv"
)

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Shared fonts, theme and palettes
SCRIPTS_R_DIR <- file.path(PROJECT_DIR, "scripts", "R")
source(file.path(SCRIPTS_R_DIR, "_pub_style.R"))

# Figure titles, with the genus italicised
species_pretty <- c(
  Thalassarchaeum_sp002507125 = "*Thalassarchaeum* sp002507125",
  MGIIb_O1_sp002498525        = "MGIIb-O1 sp002498525",
  UBA8886_sp003193815         = "UBA8886 sp003193815"
)

# --- Helpers ------------------------------------------------------------------
log_info <- function(...) {
  cat(sprintf(
    "[%s] [INFO] %s\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    paste0(...)
  ))
}
log_warn <- function(...) {
  cat(sprintf(
    "[%s] [WARN] %s\n",
    format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    paste0(...)
  ), file = stderr())
}

classify_gc <- function(presence_count, n_genomes,
                        soft_core_thr = 0.90, shell_thr = 0.15) {
  prev <- presence_count / n_genomes
  case_when(
    prev == 1 ~ "core",
    prev >= soft_core_thr ~ "soft-core",
    prev >= shell_thr ~ "shell",
    TRUE ~ "cloud"
  )
}

heaps_gamma <- function(bin_mat, n_perm = 100) {
  if (!has_micropan) {
    return(list(gamma = NA_real_, intercept = NA_real_))
  }
  # micropan::heaps wants MAGs in rows
  res <- tryCatch(
    micropan::heaps(t(bin_mat), n.perm = n_perm),
    error = function(e) {
      log_warn("Heaps falló: ", e$message)
      NULL
    }
  )
  if (is.null(res)) {
    return(list(gamma = NA_real_, intercept = NA_real_))
  }
  list(gamma = unname(res["alpha"]), intercept = unname(res["intercept"]))
}

# --- Load quality (CheckM2) ---------------------------------------------------
quality <- tibble(NameWithoutExtension = character(0), Completeness = numeric(0))
if (file.exists(QUALITY_TSV)) {
  quality <- read_tsv(QUALITY_TSV, show_col_types = FALSE) %>%
    rename_with(tolower) %>%
    mutate(
      NameWithoutExtension = str_remove(name, "\\.fa$"),
      Completeness = as.numeric(completeness)
    ) %>%
    select(NameWithoutExtension, Completeness)
} else {
  log_warn(
    "Quality report not found at ", QUALITY_TSV,
    " — sensitivity analysis will be skipped"
  )
}

# Must match the sanitisation in pipeline step 20
sanitize <- function(x) gsub("[^a-zA-Z0-9_]", "_", x)

# ---- one species at a time ----
heaps_summary <- list()
class_counts <- list()
panel_classes <- list() # subplots for Fig_panel_pan_classes
panel_curves <- list() # subplots for Fig_panel_pan_curves

canonical_species <- names(species_pretty)

bin_files <- list.files(DERIVED_DIR, pattern = "_GCxMAG_binary\\.tsv$", full.names = TRUE)
if (length(bin_files) == 0) stop("No binary matrices in ", DERIVED_DIR, ". Run script 26 first.")

# Drop the near-duplicate species names that come from hyphenated artefacts
sp_names_all <- str_replace(basename(bin_files), "_GCxMAG_binary\\.tsv$", "")
keep_idx <- sp_names_all %in% canonical_species
if (any(!keep_idx)) {
  log_info(sprintf(
    "  skipping non-canonical species files: %s",
    paste(sp_names_all[!keep_idx], collapse = ", ")
  ))
}
bin_files <- bin_files[keep_idx]

for (bin_file in bin_files) {
  sp_name <- str_replace(basename(bin_file), "_GCxMAG_binary\\.tsv$", "")
  pretty <- species_pretty[sp_name] %||% sp_name
  log_info(sprintf("Species: %s", sp_name))

  bin_df <- read_tsv(bin_file, show_col_types = FALSE)
  fn_file <- file.path(DERIVED_DIR, paste0(sp_name, "_GC_functions.tsv"))
  fn_df <- if (file.exists(fn_file)) read_tsv(fn_file, show_col_types = FALSE) else NULL

  bin_mat <- as.matrix(bin_df[, -1])
  rownames(bin_mat) <- bin_df$gene_cluster_id
  n_genomes <- ncol(bin_mat)

  presence_count <- rowSums(bin_mat)
  classes <- classify_gc(presence_count, n_genomes)

  cls_df <- tibble(
    gene_cluster_id = bin_df$gene_cluster_id,
    n_genomes_present = presence_count,
    prevalence = presence_count / n_genomes,
    class = classes
  )
  if (!is.null(fn_df)) cls_df <- cls_df %>% left_join(fn_df, by = "gene_cluster_id")

  out_cls <- file.path(OUTPUT_DIR, paste0("pangenome_classes_", sp_name, ".tsv"))
  write_tsv(cls_df, out_cls)
  log_info(sprintf("  wrote %s", basename(out_cls)))

  # Cluster counts per band
  cls_counts <- cls_df %>%
    count(class, name = "n") %>%
    mutate(species = sp_name, n_genomes = n_genomes)
  class_counts[[sp_name]] <- cls_counts

  # Heaps' law
  hg <- heaps_gamma(bin_mat)
  heaps_summary[[sp_name]] <- tibble(
    species = sp_name, n_genomes = n_genomes,
    gamma = hg$gamma, intercept = hg$intercept,
    pangenome_status = case_when(
      is.na(hg$gamma) ~ NA_character_,
      hg$gamma < 1 ~ "abierto (γ<1)",
      hg$gamma > 1 ~ "cerrado (γ>1)",
      TRUE ~ "límite"
    )
  )

  # --- Figure: core/soft-core/shell/cloud bars ---------------------------------
  cls_counts_plot <- cls_counts %>%
    mutate(class = factor(class, levels = pan_class_levels)) %>%
    tidyr::complete(class, fill = list(n = 0L))
  fig_classes <- ggplot(
    cls_counts_plot,
    aes(x = class, y = n, fill = class)
  ) +
    geom_col() +
    geom_text(aes(label = n),
      vjust = -0.3, size = 4,
      family = font_family
    ) +
    scale_fill_manual(values = pan_class_colors, drop = FALSE) +
    scale_x_discrete(labels = c(
      core = "*core*", `soft-core` = "*soft-core*",
      shell = "*shell*", cloud = "*cloud*"
    )) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(
      title = pretty,
      subtitle = sprintf(
        "n = %d MAGs | %d clústeres génicos",
        n_genomes, nrow(cls_df)
      ),
      x = "Clase pangenómica", y = "Número de clústeres génicos"
    ) +
    theme_publication(base_size = 12) +
    theme(
      legend.position = "none",
      axis.text.x = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown(),
      plot.subtitle = ggtext::element_markdown()
    )
  ggsave(file.path(OUTPUT_DIR, paste0("Fig_pan_classes_", sp_name, ".png")),
    fig_classes,
    width = 5.5, height = 4.5, dpi = 300
  )
  panel_classes[[sp_name]] <- fig_classes

  # ---- accumulation curve, rarefied by hand ----
  if (n_genomes >= 3) {
    set.seed(42)
    n_iter <- 50
    accum <- expand_grid(iter = 1:n_iter, k = 1:n_genomes) %>%
      rowwise() %>%
      mutate(
        sample = list(sample.int(n_genomes, k)),
        pan = sum(rowSums(bin_mat[, sample, drop = FALSE]) >= 1),
        core = sum(rowSums(bin_mat[, sample, drop = FALSE]) == k)
      ) %>%
      ungroup() %>%
      pivot_longer(c(pan, core), names_to = "type", values_to = "count")

    fig_curves <- ggplot(accum, aes(x = k, y = count, color = type)) +
      geom_jitter(alpha = 0.3, width = 0.15, size = 1) +
      geom_smooth(se = TRUE, method = "loess", formula = y ~ x) +
      scale_color_manual(
        values = c(pan = "#7570B3", core = "#1B9E77"),
        labels = c(pan = "Pangenoma", core = "*core*")
      ) +
      labs(
        title = pretty,
        x = "Número de genomas muestreados",
        y = "Número de clústeres génicos",
        color = NULL
      ) +
      theme_publication(base_size = 12) +
      theme(
        legend.position = "bottom",
        legend.text = ggtext::element_markdown(size = 10),
        axis.title.y = ggtext::element_markdown()
      )
    ggsave(file.path(OUTPUT_DIR, paste0("Fig_pan_curves_", sp_name, ".png")),
      fig_curves,
      width = 6, height = 4.5, dpi = 300
    )
    panel_curves[[sp_name]] <- fig_curves
  }

  # --- Cross-checks (cbb3 / HndC) ---------------------------------------------
  if (!is.null(fn_df) && "KOfam" %in% colnames(fn_df)) {
    cbb3_kos <- c("K00404", "K00405", "K00406", "K00407")
    hnd_kos <- c("K18006")
    if (sp_name == "MGIIb_O1_sp002498525") {
      hits <- cls_df %>%
        filter(str_detect(KOfam, paste(cbb3_kos, collapse = "|")))
      log_info(sprintf(
        "  cbb3 hits in MGIIb-O1: %d GCs (clases: %s)",
        nrow(hits), paste(unique(hits$class), collapse = ", ")
      ))
    }
    if (sp_name == "UBA8886_sp003193815") {
      hits <- cls_df %>%
        filter(str_detect(KOfam, paste(hnd_kos, collapse = "|")))
      log_info(sprintf(
        "  HndC hits in UBA8886: %d GCs (clases: %s)",
        nrow(hits), paste(unique(hits$class), collapse = ", ")
      ))
    }
  }

  # ---- sensitivity: reclassify using only MAGs above 80% completeness ----
  if (nrow(quality) > 0) {
    keep <- quality %>%
      filter(Completeness >= 80) %>%
      pull(NameWithoutExtension)
    keep_san <- sanitize(keep)
    keep_in_sp <- intersect(keep_san, colnames(bin_mat))
    if (length(keep_in_sp) >= 3 && length(keep_in_sp) < n_genomes) {
      bin_hq <- bin_mat[, keep_in_sp, drop = FALSE]
      cls_hq <- classify_gc(rowSums(bin_hq), ncol(bin_hq))
      conc <- mean(cls_hq == classes)
      log_info(sprintf(
        "  sensibilidad (MAGs ≥80%%): n=%d, concordancia=%.2f",
        length(keep_in_sp), conc
      ))
    }
  }
}

# ---- summaries across species ----
heaps_df <- bind_rows(heaps_summary)
write_tsv(heaps_df, file.path(OUTPUT_DIR, "heaps_gamma_summary.tsv"))
log_info("Heaps summary written")

class_df <- bind_rows(class_counts)
write_tsv(class_df, file.path(OUTPUT_DIR, "class_counts_summary.tsv"))
log_info("Class counts summary written")

# ---- panels ----
build_panel <- function(plot_list, title_main) {
  if (length(plot_list) == 0) {
    return(NULL)
  }
  # Fixed order, following species_pretty
  ord <- intersect(canonical_species, names(plot_list))
  plot_list <- plot_list[ord]
  ncol_grid <- if (length(plot_list) <= 3) length(plot_list) else 3
  nrow_grid <- ceiling(length(plot_list) / ncol_grid)
  patchwork::wrap_plots(plot_list, ncol = ncol_grid) +
    patchwork::plot_annotation(
      title = title_main,
      theme = ggplot2::theme(
        plot.title = ggtext::element_markdown(
          family = font_family, face = "bold",
          size = 18, hjust = 0.5
        )
      )
    )
}

panel_w <- function(n_subplots) max(7, 5.2 * min(3, n_subplots))
panel_h <- function(n_subplots) 4.5 * ceiling(n_subplots / 3)

# Combined panel: pangenome classes along the top, accumulation curves below.
# No panel titles, since each subplot already carries its species name, and one
# tag per row rather than per subplot. wrap_elements() is what makes tag_levels
# match at the row level instead of descending into each subplot.
if (length(panel_classes) > 0 && length(panel_curves) > 0) {
  ord <- intersect(canonical_species, names(panel_classes))
  classes_list <- panel_classes[ord]
  curves_ord <- intersect(canonical_species, names(panel_curves))
  curves_list <- panel_curves[curves_ord]

  # The species name already appears above, on the class subplot, so the curve
  # subplots drop it.
  curves_list <- lapply(curves_list, function(p) p + labs(title = NULL))

  # Only the leftmost subplot in each row keeps its y axis title. The others
  # read against it, which avoids repeating the same long label three times.
  strip_y_after_first <- function(plot_list) {
    if (length(plot_list) <= 1) return(plot_list)
    for (i in seq_along(plot_list)[-1]) {
      plot_list[[i]] <- plot_list[[i]] + labs(y = NULL)
    }
    plot_list
  }
  classes_list <- strip_y_after_first(classes_list)
  curves_list <- strip_y_after_first(curves_list)

  ncol_grid <- if (length(classes_list) <= 3) length(classes_list) else 3

  row_classes <- patchwork::wrap_plots(classes_list, ncol = ncol_grid)
  row_curves <- patchwork::wrap_plots(curves_list, ncol = ncol_grid) &
    ggplot2::theme(legend.position = "bottom")
  row_curves <- row_curves + patchwork::plot_layout(guides = "collect")

  # Tag styling matches tax_analysis.R and the other panels: lowercase letters
  # with no parentheses, bold, size 20, black. wrap_elements() keeps the tags on
  # the two rows rather than on every subplot inside them.
  combined <- (patchwork::wrap_elements(full = row_classes) /
    patchwork::wrap_elements(full = row_curves)) +
    patchwork::plot_annotation(tag_levels = "a") &
    ggplot2::theme(
      plot.tag = ggplot2::element_text(
        face = "bold", size = 20, color = "black"
      )
    )

  fig_w <- panel_w(length(classes_list))
  fig_h <- panel_h(length(classes_list)) + panel_h(length(curves_list)) + 0.5
  ggsave(file.path(OUTPUT_DIR, "Fig_panel_pan_classes_curves.png"),
    combined,
    width = fig_w, height = fig_h,
    dpi = 300, limitsize = FALSE
  )
  log_info("Combined panel (classes + curves) written")
}

log_info("Done. Outputs in: ", OUTPUT_DIR)

# The tidyverse null-coalescing operator, defined locally to avoid the import
`%||%` <- function(a, b) if (is.null(a) || is.na(a)) b else a
