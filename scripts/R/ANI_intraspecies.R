#!/usr/bin/env Rscript
################################################################################
# Script: ANI_intraspecies.R
# Within-species ANI across the 54 MAGs, from the skani output of pipeline
# step 16. Several figures, each looking at the same matrix differently.
#
#   Species with an assignment:
#     Fig1   ANI distribution per species, boxplot with jitter
#     Fig2   pairwise heatmap per species
#     Fig3   ANI by depth context, same depth against different
#     Fig4   full 54 by 54 heatmap, blocked by species
#     Fig5   each MAG's contribution to its species ANI, which is what exposes outliers
#
#   Groups with no species assignment:
#     Fig1b  ANI distribution per unassigned group
#     Fig2b  pairwise heatmaps for those groups
#     Fig3b  ANI by depth for those groups
#     Fig6   novel MAGs against their genus, assigned and unassigned together
#
# Created: 2026-04-14
# Version: 2.0
################################################################################

# From RStudio: Ctrl+Shift+Enter to source.
# From the terminal:
# Rscript scripts/R/ANI_intraspecies.R

# ============================================================
# ---- 1. project root and paths ----
# ============================================================
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

PROJECT_ROOT <- resolve_project_dir()
sparse_path <- file.path(PROJECT_ROOT, "results/16_sketching_ani/all_mags/ani_triangle_sparse.tsv")
full_path <- file.path(PROJECT_ROOT, "results/16_sketching_ani/all_mags/ani_triangle_full.tsv")
novel_path <- file.path(PROJECT_ROOT, "results/16_sketching_ani/novel_mags/novel10_vs_all54.tsv")
metadata_path <- file.path(PROJECT_ROOT, "MAGs_metadata.tsv")
output_dir <- file.path(PROJECT_ROOT, "results/SO_01/ANI")

# Minimum aligned fraction, min(AF_ref, AF_query), for an ANI value from the
# novel_vs_all run (which used --min-af 0) to be trusted. Below this, the pair
# is a single homologous contig rather than a genome-wide comparison, so it does
# not go into mat_vals.
NOVEL_AF_MIN <- 5

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ============================================================
# 2. Paquetes
# ============================================================
cat("\n=== Cargando paquetes ===\n")
required_packages <- c(
  "tidyverse", # dplyr, ggplot2, readr, tidyr, stringr, forcats
  "patchwork", # composición de múltiples figuras
  "scales", # format de ejes
  "ggrepel", # etiquetas sin solapamiento
  "ggtext", # element_markdown() para texto con formato en ejes
  "RColorBrewer", # paletas de color
  "ComplexHeatmap", # heatmap 54×54 anotado (Fig4)
  "circlize", # colorRamp2 para escalas de color
  "grid", # gráficos base (requerido por ComplexHeatmap)
  "showtext", # soporte para fuentes personalizadas
  "sysfonts" # detección de fuentes del sistema
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("  Instalando:", pkg, "...\n")
    if (pkg == "ComplexHeatmap") {
      if (!requireNamespace("BiocManager", quietly = TRUE)) {
        install.packages("BiocManager", quiet = TRUE)
      }
      BiocManager::install(pkg, update = FALSE, ask = FALSE, quiet = TRUE)
    } else {
      install.packages(pkg, dependencies = TRUE, quiet = TRUE)
    }
  }
  tryCatch(
    suppressMessages(library(pkg, character.only = TRUE, quietly = TRUE)),
    error = function(e) stop("No se pudo cargar '", pkg, "': ", conditionMessage(e))
  )
  cat("  [OK]", pkg, "\n")
}
cat("\n[OK] Todos los paquetes cargados\n\n")

# ============================================================
# ---- 3. fonts, same as tax_analysis.R ----
# ============================================================
cat("=== Configurando fuentes ===\n")
showtext_auto()
showtext_opts(dpi = 300)
font_family <- "sans"
font_found <- FALSE

local_font_dir <- path.expand("~/.local/share/fonts")
local_direct <- list(
  list(
    name = "Lato", regular = "Lato-Regular.ttf", bold = "Lato-Bold.ttf",
    italic = "Lato-Italic.ttf", bolditalic = "Lato-BoldItalic.ttf"
  ),
  list(
    name = "Source Sans 3", regular = "SourceSans3-Regular.otf",
    bold = "SourceSans3-Bold.otf", italic = "SourceSans3-It.otf",
    bolditalic = "SourceSans3-BoldIt.otf"
  )
)
if (dir.exists(local_font_dir)) {
  for (lf in local_direct) {
    reg_path <- file.path(local_font_dir, lf$regular)
    if (file.exists(reg_path)) {
      tryCatch(
        {
          sysfonts::font_add(
            family = "custom_font", regular = reg_path,
            bold = if (file.exists(file.path(local_font_dir, lf$bold))) file.path(local_font_dir, lf$bold) else NULL,
            italic = if (file.exists(file.path(local_font_dir, lf$italic))) file.path(local_font_dir, lf$italic) else NULL,
            bolditalic = if (file.exists(file.path(local_font_dir, lf$bolditalic))) file.path(local_font_dir, lf$bolditalic) else NULL
          )
          font_family <- "custom_font"
          font_found <- TRUE
          cat("  [OK] Fuente:", lf$name, "(local)\n")
          break
        },
        error = function(e) {}
      )
    }
  }
}
if (!font_found) {
  for (gf in list("Lato", "Source Sans Pro", "Open Sans")) {
    tryCatch(
      {
        sysfonts::font_add_google(gf, "custom_font")
        font_family <- "custom_font"
        font_found <- TRUE
        cat("  [OK] Fuente Google:", gf, "\n")
        break
      },
      error = function(e) {}
    )
  }
}
if (!font_found) cat("  [WARN] Usando fuente sans-serif por defecto.\n")
cat("[OK] Configuración de fuentes completada\n\n")

# ============================================================
# ---- 4. publication theme, same as tax_analysis.R ----
# ============================================================
theme_publication <- function(base_size = 14) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(color = "gray88", linetype = "dotted"),
      panel.grid.major.x = element_blank(),
      axis.text          = element_text(color = "black", size = base_size - 1),
      axis.title         = element_text(face = "bold", size = base_size + 1),
      legend.position    = "right",
      legend.background  = element_rect(fill = "white", color = NA),
      legend.title       = element_text(face = "bold", size = base_size - 1),
      legend.text        = element_text(size = base_size - 3),
      strip.background   = element_rect(fill = "gray90", color = "black"),
      strip.text         = element_text(face = "bold", size = base_size - 1),
      plot.title         = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      plot.subtitle      = element_text(hjust = 0.5, size = base_size - 1, color = "gray40"),
      text               = element_text(family = font_family)
    )
}

save_fig <- function(p, path, width, height, dpi = 300) {
  ggsave(path, plot = p, width = width, height = height, dpi = dpi, bg = "white")
  cat("  [OK] Guardada:", basename(path), "\n")
}

# ============================================================
# ---- 5. project palettes ----
# ============================================================

# Profundidad — ORDEN: ZMOR to 1000 m to Fondo ( to  OMZ al final)
# Matching tax_analysis.R and the rest of the project.
depth_label_map <- c("ROMZ" = "ZMOR", "1000" = "1000 m", "BTM" = "Fondo", "OMZ" = "OMZ")
depth_factor_levels <- c("ZMOR", "1000 m", "Fondo", "OMZ")
depth_colors <- c("ZMOR" = "#4A90D9", "1000 m" = "#E67E22", "Fondo" = "#8E44AD", "OMZ" = "#27AE60")

# Estaciones
station_levels <- c("A4", "A6", "C22", "F37", "PO1", "Y6")
station_colors <- setNames(RColorBrewer::brewer.pal(6, "Set1"), station_levels)

# Species and group, assigned species plus the unassigned groups by genus
species_palette <- c(
  "Thalassarchaeum sp002507125" = "#C0392B",
  "Thalassarchaeum betae"       = "#E74C3C",
  "MGIIb-O1 sp002498525"        = "#2471A3",
  "MGIIb-Q1 sp002505495"        = "#5DADE2",
  "UBA8886 sp003193815"         = "#1E8449",
  "UBA8886 sp937893445"         = "#52BE80",
  "UBA57 sp002494565"           = "#D68910",
  "UBA60 sp002503395"           = "#CA6F1E",
  "Undinarchaeum marinum"       = "#7D3C98",
  "JAKUTA01 sp950068635"        = "#17A589",
  # Grupos N/A (posibles nuevas especies)
  "MGIIb-Q1 N/A"                = "#AED6F1",
  "Nitrosopelagicus N/A"        = "#A9DFBF",
  "UBA36 N/A"                   = "#BFC9CA"
)

# ============================================================
# ---- 6. taxonomic name formatting ----
# ============================================================
# Copied verbatim from tax_analysis.R rather than sourced, so the two scripts
# stay visually identical even if one is edited without the other.

# CDC scientific nomenclature (https://wwwnc.cdc.gov/eid/page/scientific-nomenclature):
# only genus and species (and infraspecific ranks) are italicised; ranks above
# genus — family, order, class, phylum, domain — are plain text.
ITALIC_LEVELS <- c("Genus", "Species", "Subspecies", "Strain")

should_italicize_level <- function(level) {
  !is.null(level) && length(level) == 1 && !is.na(level) && level %in% ITALIC_LEVELS
}

# Italicises a taxon name per CDC style: genus, species and below go italic
# word by word, while GTDB codes carrying digits (sp002507125, UBA60, MGIII)
# stay upright.
#   - level en {Family, Order, Class, Phylum, Domain} to texto plano completo.
# Ejemplo: "Thalassarchaeum sp002507125" to "*Thalassarchaeum* sp002507125"
taxon_italic <- function(x, level = "Species") {
  if (!should_italicize_level(level)) {
    return(ifelse(is.na(x), x, as.character(x)))
  }
  sapply(x, function(nm) {
    if (is.na(nm) || nm %in% c("Sin asignar", "N/A")) {
      return(nm)
    }
    words <- strsplit(nm, " ")[[1]]
    is_plain <- grepl("[0-9]", words) | grepl("^[A-Z]+$", words)
    r <- rle(is_plain)
    ends <- cumsum(r$lengths)
    starts <- c(1L, ends[-length(ends)] + 1L)
    chunks <- mapply(function(s, e, plain) {
      grp <- paste(words[s:e], collapse = " ")
      if (plain) grp else paste0("*", grp, "*")
    }, starts, ends, r$values, SIMPLIFY = TRUE)
    paste(chunks, collapse = " ")
  }, USE.NAMES = FALSE)
}

# Ancestor hierarchy, as in tax_analysis.R
parent_hierarchy <- list(
  Species = "Genus",
  Genus   = "Family",
  Family  = "Order",
  Order   = "Class",
  Class   = "Phylum",
  Phylum  = "Domain"
)

# Rank abbreviations
level_abbrev <- list(
  Genus   = "g.",
  Family  = "f.",
  Order   = "o.",
  Class   = "c.",
  Phylum  = "p.",
  Domain  = "d."
)

# Construye un lookup: valor_hijo to list(name, level, abbrev)
# Climbs the hierarchy until it finds an ancestor whose name differs from the
# child. At species level the comparison uses only the first word, the genus,
# since "Thalassarchaeum betae" against "Thalassarchaeum" is not a difference.
get_parent_lookup <- function(df_in, child_col) {
  ancestors <- character(0)
  cur <- child_col
  repeat {
    par <- parent_hierarchy[[cur]]
    if (is.null(par)) break
    ancestors <- c(ancestors, par)
    cur <- par
  }
  ancestors <- ancestors[ancestors %in% names(df_in)]
  if (length(ancestors) == 0) {
    return(list())
  }

  child_vals <- unique(na.omit(df_in[[child_col]]))
  if (length(child_vals) == 0) {
    return(list())
  }

  result <- list()
  for (cv in child_vals) {
    mask <- !is.na(df_in[[child_col]]) & df_in[[child_col]] == cv
    rows <- df_in[mask, , drop = FALSE]
    if (nrow(rows) == 0) next

    child_key <- if (child_col == "Species") strsplit(cv, " ")[[1]][1] else cv

    for (anc_col in ancestors) {
      anc_vals <- na.omit(rows[[anc_col]])
      if (length(anc_vals) == 0) next
      anc_val <- names(sort(table(anc_vals), decreasing = TRUE))[1]
      if (!is.na(anc_val) && anc_val != child_key && anc_val != cv) {
        result[[cv]] <- list(
          name   = anc_val,
          level  = anc_col,
          abbrev = if (!is.null(level_abbrev[[anc_col]])) level_abbrev[[anc_col]] else ""
        )
        break
      }
    }
  }
  result
}

# Formats a taxon name: CDC italics, with the parent rank in parentheses when
# the name is a code that would otherwise say nothing.
fmt_taxon <- function(nm, parent_lookup = NULL, level = "Species") {
  if (is.na(nm) || nm %in% c("Sin asignar", "N/A")) {
    return(nm)
  }
  first_word <- strsplit(nm, " ")[[1]][1]
  is_acr <- grepl("[0-9]", first_word) || grepl("^[A-Z]+$", first_word)
  if (is_acr && !is.null(parent_lookup) && length(parent_lookup) > 0 &&
    nm %in% names(parent_lookup)) {
    info <- parent_lookup[[nm]]
    if (!is.null(info) && !is.na(info$name) && nchar(trimws(info$name)) > 0) {
      abbrev_str <- if (nchar(info$abbrev) > 0) paste0(info$abbrev, " ") else ""
      return(paste0(
        taxon_italic(nm, level),
        " (", abbrev_str, taxon_italic(info$name, info$level), ")"
      ))
    }
  }
  taxon_italic(nm, level)
}

# Vectorised fmt_taxon
fmt_taxon_vec <- function(x, parent_lookup = NULL, level = "Species") {
  sapply(x, fmt_taxon,
    parent_lookup = parent_lookup, level = level, USE.NAMES = FALSE
  )
}

# Returns a labeller closed over a given parent_lookup and rank
make_label_fn <- function(parent_lookup, level = "Species") {
  function(x) fmt_taxon_vec(x, parent_lookup, level = level)
}

# Genus label for an unassigned group: italic genus followed by N/A
fmt_genus_na <- function(genus) {
  paste0(taxon_italic(genus, level = "Genus"), " N/A")
}

# Breaks a long label at the space nearest its midpoint, so ggtext axes do not
# overflow.
wrap_md_label <- function(lbl, width = 28) {
  if (nchar(lbl) <= width) {
    return(lbl)
  }
  sp <- gregexpr(" ", lbl)[[1]]
  sp <- sp[sp > 0]
  if (length(sp) == 0) {
    return(lbl)
  }
  mid <- nchar(lbl) %/% 2
  best <- sp[which.min(abs(sp - mid))]
  paste0(substr(lbl, 1, best - 1), "<br>", substr(lbl, best + 1, nchar(lbl)))
}

# ============================================================
# ---- 7. load and prepare the data ----
# ============================================================
cat("=== Cargando datos ===\n")

# — Metadata —
meta_raw <- read_tsv(metadata_path, show_col_types = FALSE)

# Split FullTaxonomy ("d__X;p__X;...;s__X") into ranks
parse_gtdb_level <- function(classification, prefix) {
  m <- regmatches(
    classification,
    regexpr(paste0("(?<=", prefix, ")[^;]+"), classification, perl = TRUE)
  )
  if (length(m) == 0 || m == "") NA_character_ else m
}

taxonomy_cols <- do.call(rbind, lapply(meta_raw$FullTaxonomy, function(tx) {
  data.frame(
    Domain = parse_gtdb_level(tx, "d__"),
    Phylum = parse_gtdb_level(tx, "p__"),
    Class = parse_gtdb_level(tx, "c__"),
    Order = parse_gtdb_level(tx, "o__"),
    Family = parse_gtdb_level(tx, "f__"),
    stringsAsFactors = FALSE
  )
}))

meta_df <- bind_cols(meta_raw, taxonomy_cols) |>
  mutate(
    # Depth as an ordered factor: ZMOR, 1000 m, Fondo, OMZ
    Depth_label = factor(
      recode(Depth, !!!depth_label_map),
      levels = depth_factor_levels
    ),
    short_label = paste0(Station, " ", Depth, "_", BinID),
    # Group: the assigned species, or "<genus> N/A" where there is none
    Species_group = if_else(Species == "N/A", paste0(Genus, " N/A"), Species)
  )
cat("  [OK] Metadata:", nrow(meta_df), "MAGs\n")

# Parent lookup at species level
# (requiere columnas Genus, Family, Order, Class, Phylum, Domain en meta_df)
species_parent_lookup <- get_parent_lookup(meta_df, "Species")
cat("  [OK] Parent lookup para", length(species_parent_lookup), "especies\n")

# ---- sparse ANI: edge list of pairs above 75% ----
# Formato actual de skani sparse: Ref_file, Query_file, ANI, Align_fraction_ref,
# Align_fraction_query, Ref_name, Query_name. Las dos primeras columnas son
# skani reports absolute .fa paths, so the basename without extension is what
# joins against meta_df$NameWithoutExtension.
sparse_raw <- read_tsv(sparse_path, show_col_types = FALSE)
sparse_raw <- sparse_raw |>
  transmute(
    Ref_MAG   = tools::file_path_sans_ext(basename(as.character(Ref_file))),
    Query_MAG = tools::file_path_sans_ext(basename(as.character(Query_file))),
    ANI       = as.numeric(ANI),
    AF_ref    = as.numeric(Align_fraction_ref),
    AF_query  = as.numeric(Align_fraction_query)
  )

# Join the metadata onto both Ref and Query
meta_sel <- meta_df |>
  select(
    NameWithoutExtension, Station, Depth, Depth_label, BinID,
    Genus, Species, Species_group, short_label
  )

sparse_df <- sparse_raw |>
  left_join(meta_sel, by = c("Ref_MAG" = "NameWithoutExtension")) |>
  rename_with(~ paste0(., "_ref"), .cols = Station:short_label) |>
  left_join(meta_sel, by = c("Query_MAG" = "NameWithoutExtension")) |>
  rename_with(~ paste0(., "_query"), .cols = Station:short_label) |>
  mutate(
    same_station = Station_ref == Station_query,
    same_depth = Depth_ref == Depth_query,
    intraspecies = (Species_ref == Species_query) & (Species_ref != "N/A"),
    intragroup_na = (Genus_ref == Genus_query) & (Species_ref == "N/A") &
      (Species_query == "N/A"),
    intragenus_mix = (Genus_ref == Genus_query) & (Species_ref != Species_query) &
      (Species_ref == "N/A" | Species_query == "N/A")
  )
cat("  [OK] Sparse ANI:", nrow(sparse_df), "pares totales\n")
cat("  [OK] Pares intra-especie (asignados):", sum(sparse_df$intraspecies, na.rm = TRUE), "\n")
cat("  [OK] Pares intra-género N/A:         ", sum(sparse_df$intragroup_na, na.rm = TRUE), "\n")

# — Matriz completa (formato phylip 54×54) —
cat("  Leyendo matriz completa 54×54...\n")
full_lines <- readLines(full_path)
n_mags <- as.integer(trimws(full_lines[1]))
mag_names <- character(n_mags)
mat_vals <- matrix(NA_real_, nrow = n_mags, ncol = n_mags)

for (i in seq_len(n_mags)) {
  parts <- strsplit(full_lines[i + 1L], "\t")[[1]]
  mag_names[i] <- tools::file_path_sans_ext(basename(parts[1]))
  vals <- as.numeric(parts[-1])
  vals[vals == 0] <- NA_real_
  vals[i] <- 100
  mat_vals[i, ] <- vals
}
rownames(mat_vals) <- colnames(mat_vals) <- mag_names
cat("  [OK] Matriz completa:", n_mags, "×", n_mags, "MAGs\n\n")

# ---- fill the gaps from the novel_vs_all run ----
# The all_mags run used skani's default --min-af of about 15%, so MAG pairs with
# a low aligned fraction never appear and sit as NA in mat_vals. The
# novel10_vs_all54.tsv run used --min-af 0 and does catch them. Only pairs at or
# above NOVEL_AF_MIN come across, which discards the artefacts where ANI reads
# 100% off a single tiny contig, and only NA cells are filled: the values
# already there came from higher aligned fractions and are the better estimate.
if (file.exists(novel_path)) {
  cat("  Leyendo corrida novel_vs_all (rescate de pares faltantes)...\n")
  novel_raw <- read_tsv(novel_path, show_col_types = FALSE)
  novel_pairs <- novel_raw |>
    transmute(
      Ref_MAG   = tools::file_path_sans_ext(basename(as.character(Ref_file))),
      Query_MAG = tools::file_path_sans_ext(basename(as.character(Query_file))),
      ANI       = as.numeric(ANI),
      AF_ref    = as.numeric(Align_fraction_ref),
      AF_query  = as.numeric(Align_fraction_query)
    ) |>
    filter(Ref_MAG != Query_MAG) |>
    mutate(AF_min = pmin(AF_ref, AF_query, na.rm = TRUE)) |>
    filter(AF_min >= NOVEL_AF_MIN) |>
    filter(Ref_MAG %in% mag_names, Query_MAG %in% mag_names)

  n_filled  <- 0L
  n_skipped <- 0L
  for (k in seq_len(nrow(novel_pairs))) {
    r <- novel_pairs$Ref_MAG[k]
    q <- novel_pairs$Query_MAG[k]
    a <- novel_pairs$ANI[k]
    if (is.na(mat_vals[r, q])) {
      mat_vals[r, q] <- a
      n_filled <- n_filled + 1L
    } else {
      n_skipped <- n_skipped + 1L
    }
    if (is.na(mat_vals[q, r])) {
      mat_vals[q, r] <- a
    }
  }
  cat(sprintf(
    "  [OK] Rescate AF >= %g%%: %d celdas rellenadas, %d preservadas (ya tenían valor)\n",
    NOVEL_AF_MIN, n_filled, n_skipped
  ))

  # sparse_df is extended too, so Fig1b and the distributions include the
  # rescued pairs with their full metadata. The pair key is unordered, which is
  # what stops a pair already present from being duplicated under
  # Ref/Query intercambiados.
  unord_key <- function(a, b) ifelse(a < b, paste(a, b, sep = "|"),
                                     paste(b, a, sep = "|"))
  existing_keys <- unord_key(sparse_raw$Ref_MAG, sparse_raw$Query_MAG)
  novel_sparse <- novel_pairs |>
    select(Ref_MAG, Query_MAG, ANI, AF_ref, AF_query) |>
    filter(!unord_key(Ref_MAG, Query_MAG) %in% existing_keys) |>
    distinct(pair = unord_key(Ref_MAG, Query_MAG),
             .keep_all = TRUE) |>
    select(-pair)
  if (nrow(novel_sparse) > 0) {
    novel_sparse_df <- novel_sparse |>
      left_join(meta_sel, by = c("Ref_MAG" = "NameWithoutExtension")) |>
      rename_with(~ paste0(., "_ref"), .cols = Station:short_label) |>
      left_join(meta_sel, by = c("Query_MAG" = "NameWithoutExtension")) |>
      rename_with(~ paste0(., "_query"), .cols = Station:short_label) |>
      mutate(
        same_station = Station_ref == Station_query,
        same_depth = Depth_ref == Depth_query,
        intraspecies = (Species_ref == Species_query) & (Species_ref != "N/A"),
        intragroup_na = (Genus_ref == Genus_query) & (Species_ref == "N/A") &
          (Species_query == "N/A"),
        intragenus_mix = (Genus_ref == Genus_query) & (Species_ref != Species_query) &
          (Species_ref == "N/A" | Species_query == "N/A")
      )
    sparse_df <- bind_rows(sparse_df, novel_sparse_df)
    cat(sprintf("  [OK] %d pares añadidos a sparse_df desde novel_vs_all\n",
                nrow(novel_sparse_df)))
  }
} else {
  cat("  [WARN] No se encontró novel_path; mat_vals queda con NAs originales\n")
}
cat("\n")

# Species with at least 2 assigned MAGs, the minimum for a single pair
species_ge2 <- meta_df |>
  filter(Species != "N/A") |>
  count(Species) |>
  filter(n >= 2) |>
  pull(Species)

# Species with at least 3, which the per-species heatmaps and Fig5 need
species_ge3 <- meta_df |>
  filter(Species != "N/A") |>
  count(Species) |>
  filter(n >= 3) |>
  pull(Species)

# Unassigned groups with at least 2 MAGs in the same genus
na_ge2 <- meta_df |>
  filter(Species == "N/A") |>
  count(Species_group) |>
  filter(n >= 2) |>
  pull(Species_group)

# Unassigned groups with at least 3, for the heatmaps
na_ge3 <- meta_df |>
  filter(Species == "N/A") |>
  count(Species_group) |>
  filter(n >= 3) |>
  pull(Species_group)

# Pull the submatrix for one Species_group out of the full matrix
get_group_submat <- function(group_val) {
  sp_mags <- meta_df |>
    filter(Species_group == group_val) |>
    pull(NameWithoutExtension) |>
    intersect(mag_names)
  if (length(sp_mags) < 2) {
    return(NULL)
  }
  mat_vals[sp_mags, sp_mags]
}

cat("[OK] Datos preparados\n\n")

# ============================================================
# ---- Fig1: within-species ANI distribution, assigned species ----
# ============================================================
# Species on the y axis with formatted names. MAGs without a species go to
# Fig1b instead, and a species holding a single MAG is dropped because it
# yields no pairs.
# ============================================================
cat("=== Fig1: Distribución ANI — especies asignadas ===\n")

fig1_data <- sparse_df |>
  filter(intraspecies, Species_ref %in% species_ge2) |>
  transmute(ANI, AF_ref, AF_query,
    group_label = Species_ref,
    Station_ref, Depth_label_ref
  )

fig1_summary <- fig1_data |>
  group_by(group_label) |>
  summarise(median_ani = median(ANI), n_pairs = n(), .groups = "drop")

mag_per_sp <- meta_df |>
  filter(Species %in% species_ge2) |>
  count(Species) |>
  rename(group_label = Species, n_mags = n)

fig1_summary <- fig1_summary |>
  left_join(mag_per_sp, by = "group_label") |>
  arrange(median_ani) |>
  mutate(
    axis_label = paste0(
      fmt_taxon_vec(group_label, species_parent_lookup),
      sprintf("<br><span style='font-size:13px;color:#555555'>n=%d | %d pares</span>",
              n_mags, n_pairs)
    )
  )

axis_labels_fig1 <- setNames(fig1_summary$axis_label, fig1_summary$group_label)
fig1_data <- fig1_data |>
  mutate(group_fct = factor(group_label, levels = fig1_summary$group_label))

p_fig1 <- ggplot(fig1_data, aes(x = ANI, y = group_fct)) +
  geom_boxplot(
    aes(fill = group_fct),
    color = "gray35",
    alpha = 0.25, outlier.shape = NA, width = 0.55, linewidth = 0.5
  ) +
  geom_jitter(
    aes(color = group_fct),
    height = 0.18, size = 2.0, alpha = 0.85, stroke = 0.3
  ) +
  # The species threshold, drawn as a bare line
  geom_vline(
    xintercept = 95, color = "#C0392B",
    linetype = "dashed", linewidth = 0.85
  ) +
  scale_x_continuous(
    limits = c(NA, 100.5),
    breaks = seq(95, 100, by = 1),
    labels = function(x) paste0(x, "%"),
    expand = expansion(mult = c(0.02, 0))
  ) +
  scale_y_discrete(labels = axis_labels_fig1) +
  scale_color_manual(values = species_palette, guide = "none") +
  scale_fill_manual(values = species_palette, guide = "none") +
  labs(x = "ANI (%)", y = NULL) +
  theme_publication(base_size = 13) +
  theme(
    axis.text.y = element_markdown(size = 12.5),
    panel.grid.major.x = element_line(color = "gray88", linetype = "dotted"),
    panel.grid.major.y = element_blank()
  )

save_fig(p_fig1, file.path(output_dir, "Fig1_ANI_distribucion_intraespecie.png"),
  width = 8.5, height = 7.5
)

# ============================================================
# ---- Fig1b: within-genus ANI distribution, unassigned groups ----
# ============================================================
cat("=== Fig1b: Distribución ANI — grupos N/A ===\n")

fig1b_data <- sparse_df |>
  filter(intragroup_na, Species_group_ref %in% na_ge2) |>
  transmute(ANI, AF_ref, AF_query,
    group_label = Species_group_ref,
    Genus_ref, Depth_label_ref
  )

if (nrow(fig1b_data) > 0) {
  fig1b_summary <- fig1b_data |>
    group_by(group_label) |>
    summarise(median_ani = median(ANI), n_pairs = n(), .groups = "drop")

  mag_per_na <- meta_df |>
    filter(Species == "N/A", Species_group %in% na_ge2) |>
    count(Species_group) |>
    rename(group_label = Species_group, n_mags = n)

  fig1b_summary <- fig1b_summary |>
    left_join(mag_per_na, by = "group_label") |>
    arrange(median_ani) |>
    mutate(
      genus_part = sub(" N/A$", "", group_label),
      axis_label = paste0(
        taxon_italic(genus_part, level = "Genus"),
        "<br><span style='font-size:13px;color:#555555'>n=",
        n_mags, " | ", n_pairs, " pares</span>"
      )
    )

  # More saturated colours here, since these panels carry fewer points
  na_palette_fig1b <- c(
    "MGIIb-Q1 N/A"        = "#2471A3",
    "Nitrosopelagicus N/A" = "#1E8449",
    "UBA36 N/A"            = "#717D7E"
  )
  axis_labels_fig1b <- setNames(fig1b_summary$axis_label, fig1b_summary$group_label)
  fig1b_data <- fig1b_data |>
    mutate(group_fct = factor(group_label, levels = fig1b_summary$group_label))

  p_fig1b <- ggplot(fig1b_data, aes(x = ANI, y = group_fct)) +
    geom_boxplot(
      aes(fill = group_fct),
      color = "gray35",
      alpha = 0.25, outlier.shape = NA, width = 0.35, linewidth = 0.5
    ) +
    geom_jitter(
      aes(color = group_fct),
      height = 0.10, size = 2.5, alpha = 0.88, stroke = 0.3
    ) +
    geom_vline(
      xintercept = 95, color = "#C0392B",
      linetype = "dashed", linewidth = 0.85
    ) +
    scale_x_continuous(
      breaks = sort(unique(c(80, seq(85, 95, by = 5), seq(95, 100, by = 1)))),
      labels = function(x) paste0(x, "%")
    ) +
    scale_y_discrete(labels = axis_labels_fig1b) +
    scale_color_manual(values = na_palette_fig1b, guide = "none") +
    scale_fill_manual(values = na_palette_fig1b, guide = "none") +
    labs(x = "ANI (%)", y = NULL) +
    theme_publication(base_size = 13) +
    theme(
      axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5),
      axis.text.y = element_markdown(size = 12.5),
      panel.grid.major.x = element_line(color = "gray88", linetype = "dotted"),
      panel.grid.major.y = element_blank()
    )

  save_fig(p_fig1b, file.path(output_dir, "Fig1b_ANI_distribucion_grupos_NA.png"),
    width = 8, height = 5
  )

  # ----------------------------------------------------------
  # Fig1c — Panel combinado: (a) Fig1 + (b) Fig1b
  # Both subplots the same size, tagged a and b
  # ----------------------------------------------------------
  cat("=== Fig1c: Panel combinado (a) Fig1 + (b) Fig1b ===\n")
  # Subplot (b) at about 55% the height of (a)
  p_fig1c <- (p_fig1 / p_fig1b) +
    plot_layout(heights = c(1, 0.45)) +
    plot_annotation(tag_levels = "a") &
    theme(
      plot.tag = element_text(face = "bold", size = 18, color = "black")
    )
  save_fig(p_fig1c,
    file.path(output_dir, "Fig1c_ANI_distribucion_panel.png"),
    width = 8.5, height = 12
  )
} else {
  cat("  [WARN] Sin pares intra-grupo N/A con n ≥ 2 MAGs. Fig1b omitida.\n")
}

# ============================================================
# ---- Fig2: pairwise heatmaps per species ----
# ============================================================
cat("=== Fig2: Heatmaps por especie ===\n")

# Build a ggplot heatmap for one group
make_species_heatmap <- function(group_val, is_na_group = FALSE) {
  submat <- get_group_submat(group_val)
  if (is.null(submat)) {
    return(NULL)
  }

  sp_mags <- rownames(submat)
  lbl_map <- setNames(meta_df$short_label, meta_df$NameWithoutExtension)
  dp_map <- setNames(as.character(meta_df$Depth_label), meta_df$NameWithoutExtension)

  long_df <- as.data.frame(as.table(submat)) |>
    rename(MAG_x = Var1, MAG_y = Var2, ANI = Freq) |>
    mutate(
      MAG_x   = factor(MAG_x, levels = sp_mags),
      MAG_y   = factor(MAG_y, levels = sp_mags),
      Lbl_x   = lbl_map[as.character(MAG_x)],
      Lbl_y   = lbl_map[as.character(MAG_y)],
      Depth_x = dp_map[as.character(MAG_x)],
      # Numeric cells print as %.1f. NA cells, which skani dropped for low
      # aligned fraction, print as "<AF" rather than blank, so the reader sees
      # a deliberate exclusion instead of missing data.
      ANI_lab = if_else(!is.na(ANI), sprintf("%.1f", ANI), "<AF")
    )

  off_diag <- long_df$ANI[as.character(long_df$MAG_x) != as.character(long_df$MAG_y)]
  ani_min <- floor(min(off_diag, na.rm = TRUE))
  col_sp <- species_palette[[group_val]]
  if (is.null(col_sp) || is.na(col_sp)) col_sp <- "#555555"

  # Assigned groups get the full fmt_taxon title; unassigned ones show only the italic genus
  title_str <- if (is_na_group) {
    taxon_italic(sub(" N/A$", "", group_val), level = "Genus")
  } else {
    fmt_taxon(group_val, species_parent_lookup)
  }

  ggplot(long_df, aes(x = Lbl_x, y = Lbl_y, fill = ANI)) +
    geom_tile(color = "white", linewidth = 0.4) +
    geom_text(
      aes(label = ANI_lab,
          color = if_else(is.na(ANI), "na", "val")),
      size = 2.2, show.legend = FALSE
    ) +
    scale_color_manual(values = c(na = "gray45", val = "black"),
                       guide = "none") +
    scale_fill_gradient(
      low = "white", high = col_sp,
      limits = c(ani_min, 100),
      na.value = "gray92",
      name = "ANI (%)"
    ) +
    scale_x_discrete(guide = guide_axis(angle = 45)) +
    scale_y_discrete(limits = rev) +
    labs(title = title_str, x = NULL, y = NULL) +
    theme_publication(base_size = 8) +
    theme(
      plot.title = element_markdown(size = 8),
      axis.text = element_text(size = 8.5),
      legend.position = "none",
      panel.border = element_rect(color = col_sp, fill = NA, linewidth = 1.2),
      aspect.ratio = 1
    )
}

# Heatmaps for assigned species with at least 3 MAGs
heatmap_plots_main <- lapply(
  species_ge3,
  make_species_heatmap,
  is_na_group = FALSE
)
names(heatmap_plots_main) <- species_ge3
heatmap_plots_main <- Filter(Negate(is.null), heatmap_plots_main)

ncols_fig2 <- 3L
nrows_fig2 <- ceiling(length(heatmap_plots_main) / ncols_fig2)

ani_caption <- paste0(
  "Celdas marcadas \"<AF\" (gris): skani no estimó ANI por baja fracción ",
  "alineada (<~13 %% del genoma). No es un dato ausente: indica que los ",
  "genomas son demasiado divergentes para una estimación fiable."
)

p_fig2 <- wrap_plots(heatmap_plots_main, ncol = ncols_fig2) +
  plot_annotation(caption = ani_caption) &
  theme(plot.caption = element_text(size = 8, color = "gray35",
                                     hjust = 0, lineheight = 1.2))

save_fig(p_fig2, file.path(output_dir, "Fig2_ANI_heatmaps_por_especie.png"),
  width = 4.2 * ncols_fig2, height = 4.0 * nrows_fig2 + 1.0
)

# ============================================================
# ---- Fig2b: heatmaps for unassigned groups with at least 3 MAGs ----
# ============================================================
cat("=== Fig2b: Heatmaps — grupos N/A ===\n")

if (length(na_ge3) > 0) {
  heatmap_plots_na <- lapply(na_ge3, make_species_heatmap, is_na_group = TRUE)
  names(heatmap_plots_na) <- na_ge3
  heatmap_plots_na <- Filter(Negate(is.null), heatmap_plots_na)

  if (length(heatmap_plots_na) > 0) {
    ncols_fig2b <- min(3L, length(heatmap_plots_na))
    p_fig2b <- wrap_plots(heatmap_plots_na, ncol = ncols_fig2b) +
      plot_annotation(caption = ani_caption) &
      theme(plot.caption = element_text(size = 8, color = "gray35",
                                         hjust = 0, lineheight = 1.2))
    save_fig(p_fig2b, file.path(output_dir, "Fig2b_ANI_heatmaps_grupos_NA.png"),
      width = 4.5 * ncols_fig2b, height = 5.5
    )
  }
} else {
  cat("  [WARN] Sin grupos N/A con ≥ 3 MAGs. Fig2b omitida.\n")
}

# ============================================================
# ---- Fig3: ANI by depth context ----
# ============================================================
# Pairs from the same depth against pairs spanning depths, for each species
# with at least 3 MAGs. The clearest case is Thalassarchaeum sp002507125, where
# the OMZ MAG sits consistently further from the rest than the 1000 m and BTM
# MAGs sit from each other.
# ============================================================
cat("=== Fig3: ANI por profundidad — especies asignadas ===\n")

fig3_data <- sparse_df |>
  filter(intraspecies, Species_ref %in% species_ge3) |>
  mutate(
    contexto_prof = if_else(same_depth, "Misma profundidad", "Profundidades distintas"),
    group_label   = Species_ref
  ) |>
  group_by(group_label) |>
  filter(n_distinct(contexto_prof) > 1) |>
  ungroup()

if (nrow(fig3_data) > 0) {
  group_delta_fig3 <- fig3_data |>
    group_by(group_label, contexto_prof) |>
    summarise(med = median(ANI), .groups = "drop") |>
    pivot_wider(names_from = contexto_prof, values_from = med) |>
    mutate(delta = abs(`Misma profundidad` - `Profundidades distintas`)) |>
    arrange(desc(delta))

  group_order_fig3 <- group_delta_fig3$group_label
  axis_labels_fig3 <- setNames(
    sapply(group_order_fig3, function(s) fmt_taxon(s, species_parent_lookup)),
    group_order_fig3
  )
  fig3_data <- fig3_data |>
    mutate(group_fct = factor(group_label, levels = group_order_fig3))

  prof_colors <- c("Misma profundidad" = "#2471A3", "Profundidades distintas" = "#CA6F1E")

  p_fig3 <- ggplot(fig3_data, aes(x = group_fct, y = ANI, fill = contexto_prof)) +
    geom_boxplot(
      alpha = 0.65, outlier.shape = NA,
      width = 0.55, position = position_dodge(0.7),
      linewidth = 0.5
    ) +
    geom_jitter(
      aes(color = contexto_prof),
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7),
      size = 1.8, alpha = 0.8, stroke = 0.3
    ) +
    geom_hline(
      yintercept = 95, color = "#C0392B",
      linetype = "dashed", linewidth = 0.75
    ) +
    scale_fill_manual(values = prof_colors, name = "Contexto\nde profundidad") +
    scale_color_manual(values = prof_colors, guide = "none") +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      breaks = c(90, 92, 94, 95, 96, 97, 98, 99, 100)
    ) +
    scale_x_discrete(labels = axis_labels_fig3) +
    labs(x = NULL, y = "ANI (%)") +
    theme_publication(base_size = 13) +
    theme(
      legend.position = "top",
      axis.text.x = element_markdown(size = 10, angle = 45, hjust = 1, vjust = 1)
    )

  save_fig(p_fig3, file.path(output_dir, "Fig3_ANI_por_profundidad.png"),
    width = 11, height = 7
  )
} else {
  cat("  [WARN] Sin grupos con ambos contextos de profundidad. Fig3 omitida.\n")
}

# ============================================================
# ---- Fig3c: ANI by station context ----
# ============================================================
cat("=== Fig3c: ANI por estación — especies asignadas ===\n")

fig3c_data <- sparse_df |>
  filter(intraspecies, Species_ref %in% species_ge3) |>
  mutate(
    contexto_est = if_else(same_station, "Misma estación", "Estaciones distintas"),
    group_label  = Species_ref
  ) |>
  group_by(group_label) |>
  filter(n_distinct(contexto_est) > 1) |>
  ungroup()

if (nrow(fig3c_data) > 0) {
  group_delta_fig3c <- fig3c_data |>
    group_by(group_label, contexto_est) |>
    summarise(med = median(ANI), .groups = "drop") |>
    pivot_wider(names_from = contexto_est, values_from = med) |>
    mutate(delta = abs(`Misma estación` - `Estaciones distintas`)) |>
    arrange(desc(delta))

  group_order_fig3c <- group_delta_fig3c$group_label
  axis_labels_fig3c <- setNames(
    sapply(group_order_fig3c, function(s) fmt_taxon(s, species_parent_lookup)),
    group_order_fig3c
  )
  fig3c_data <- fig3c_data |>
    mutate(group_fct = factor(group_label, levels = group_order_fig3c))

  est_colors <- c("Misma estación" = "#1A6634", "Estaciones distintas" = "#8E44AD")

  p_fig3c <- ggplot(fig3c_data, aes(x = group_fct, y = ANI, fill = contexto_est)) +
    geom_boxplot(
      alpha = 0.65, outlier.shape = NA,
      width = 0.55, position = position_dodge(0.7),
      linewidth = 0.5
    ) +
    geom_jitter(
      aes(color = contexto_est),
      position = position_jitterdodge(jitter.width = 0.12, dodge.width = 0.7),
      size = 1.8, alpha = 0.8, stroke = 0.3
    ) +
    geom_hline(
      yintercept = 95, color = "#C0392B",
      linetype = "dashed", linewidth = 0.75
    ) +
    scale_fill_manual(values = est_colors, name = "Contexto\nde estación") +
    scale_color_manual(values = est_colors, guide = "none") +
    scale_y_continuous(
      labels = function(x) paste0(x, "%"),
      breaks = c(90, 92, 94, 95, 96, 97, 98, 99, 100)
    ) +
    scale_x_discrete(labels = axis_labels_fig3c) +
    labs(x = NULL, y = "ANI (%)") +
    theme_publication(base_size = 13) +
    theme(
      legend.position = "top",
      axis.text.x = element_markdown(size = 10, angle = 45, hjust = 1, vjust = 1)
    )

  save_fig(p_fig3c, file.path(output_dir, "Fig3c_ANI_por_estacion.png"),
    width = 11, height = 7
  )
} else {
  cat("  [WARN] Sin grupos con ambos contextos de estación. Fig3c omitida.\n")
}

# ============================================================
# ---- Fig3b: ANI by depth, unassigned groups ----
# ============================================================
cat("=== Fig3b: ANI por profundidad — grupos N/A ===\n")

fig3b_data <- sparse_df |>
  filter(intragroup_na, Species_group_ref %in% na_ge2) |>
  mutate(
    contexto_prof = if_else(same_depth, "Misma profundidad", "Profundidades distintas"),
    group_label   = Species_group_ref
  ) |>
  group_by(group_label) |>
  filter(n_distinct(contexto_prof) > 1) |>
  ungroup()

if (nrow(fig3b_data) > 0) {
  axis_labels_fig3b <- setNames(
    sapply(unique(fig3b_data$group_label), function(g) {
      paste0(taxon_italic(sub(" N/A$", "", g), level = "Genus"), " N/A")
    }),
    unique(fig3b_data$group_label)
  )
  prof_colors <- c("Misma profundidad" = "#2471A3", "Profundidades distintas" = "#CA6F1E")

  p_fig3b <- ggplot(fig3b_data, aes(x = group_label, y = ANI, fill = contexto_prof)) +
    geom_boxplot(
      alpha = 0.65, outlier.shape = 21, outlier.size = 2.5,
      width = 0.5, position = position_dodge(0.65), linewidth = 0.5
    ) +
    geom_jitter(
      aes(color = contexto_prof),
      position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.65),
      size = 2.5, alpha = 0.85
    ) +
    geom_hline(
      yintercept = 95, color = "#C0392B",
      linetype = "dashed", linewidth = 0.75
    ) +
    scale_fill_manual(values = prof_colors, name = "Contexto") +
    scale_color_manual(values = prof_colors, guide = "none") +
    scale_x_discrete(labels = axis_labels_fig3b) +
    scale_y_continuous(labels = function(x) paste0(x, "%")) +
    labs(x = NULL, y = "ANI (%)") +
    theme_publication(base_size = 13) +
    theme(
      legend.position = "top",
      axis.text.x = element_markdown(size = 11)
    )

  save_fig(p_fig3b, file.path(output_dir, "Fig3b_ANI_por_profundidad_NA.png"),
    width = 9, height = 6
  )
} else {
  cat("  [WARN] Sin grupos N/A con ambos contextos de profundidad. Fig3b omitida.\n")
}

# ============================================================
# ---- Fig4: full 54 by 54 heatmap ----
# ============================================================
cat("=== Fig4: Heatmap completo 54×54 ===\n")

group_order_fig4 <- c(
  "Thalassarchaeum sp002507125", "Thalassarchaeum betae",
  "MGIIb-O1 sp002498525",
  "MGIIb-Q1 sp002505495",
  "UBA8886 sp003193815", "UBA8886 sp937893445",
  "UBA57 sp002494565",
  "UBA60 sp002503395",
  "Undinarchaeum marinum",
  "JAKUTA01 sp950068635"
)

meta_ordered <- meta_df |>
  filter(Species != "N/A") |>
  mutate(
    Species_group = factor(Species_group, levels = group_order_fig4),
    # Factor de profundidad ya definido correctamente (ZMOR to 1000 m to Fondo to OMZ)
    Depth_label   = factor(Depth_label, levels = depth_factor_levels)
  ) |>
  arrange(Species_group, Depth_label, Station)

ordered_mags <- meta_ordered$NameWithoutExtension
ordered_mags <- ordered_mags[ordered_mags %in% mag_names]
mat_ordered <- mat_vals[ordered_mags, ordered_mags]

short_lbls <- setNames(meta_df$short_label, meta_df$NameWithoutExtension)
row_labels <- short_lbls[ordered_mags]

row_split <- meta_ordered |>
  filter(NameWithoutExtension %in% ordered_mags) |>
  arrange(match(NameWithoutExtension, ordered_mags)) |>
  pull(Species_group)

col_fun <- colorRamp2(
  c(75, 85, 90, 95, 97, 99, 100),
  c("#FFF5F0", "#FCBBA1", "#FB6A4A", "#DE2D26", "#A50F15", "#67000D", "#330000")
)

depth_lbls_ordered <- meta_ordered |>
  filter(NameWithoutExtension %in% ordered_mags) |>
  arrange(match(NameWithoutExtension, ordered_mags)) |>
  pull(Depth_label)

station_ordered <- meta_ordered |>
  filter(NameWithoutExtension %in% ordered_mags) |>
  arrange(match(NameWithoutExtension, ordered_mags)) |>
  pull(Station)

# Colores de profundidad en orden correcto (ZMOR to 1000 m to Fondo to OMZ)
depth_colors_fig4 <- depth_colors[intersect(depth_factor_levels, unique(as.character(depth_lbls_ordered)))]

ha_row <- rowAnnotation(
  Profundidad = as.character(depth_lbls_ordered),
  Estación = station_ordered,
  col = list(
    Profundidad = depth_colors_fig4,
    Estación    = station_colors
  ),
  annotation_name_gp   = gpar(fontsize = 9, fontface = "bold"),
  annotation_name_side = "bottom",
  annotation_legend_param = list(
    Profundidad = list(
      title = "Profundidad", title_gp = gpar(fontsize = 9, fontface = "bold"),
      # Keep the levels in order
      at = intersect(depth_factor_levels, unique(as.character(depth_lbls_ordered)))
    ),
    Estación = list(title = "Estación", title_gp = gpar(fontsize = 9, fontface = "bold"))
  )
)

# Shared parameters for both Fig4 variants
fig4_common_params <- list(
  name = "ANI (%)",
  col = col_fun,
  na_col = "gray92",
  row_labels = row_labels,
  column_labels = row_labels,
  row_split = row_split,
  column_split = row_split,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  row_gap = unit(2, "mm"),
  column_gap = unit(2, "mm"),
  # Species titles at 45 degrees
  row_title_gp  = gpar(fontsize = 7.5, fontface = "italic"),
  row_title_rot = 0,       # horizontal para los bloques de fila
  column_title_gp  = gpar(fontsize = 7.5, fontface = "italic"),
  column_title_rot = 45,   # 45° para los bloques de columna (parte superior)
  # MAG labels unrotated
  row_names_gp    = gpar(fontsize = 7),
  column_names_gp = gpar(fontsize = 7),
  column_names_rot = 60,
  heatmap_legend_param = list(
    title     = "ANI (%)",
    at        = c(75, 80, 85, 90, 95, 100),
    labels    = paste0(c(75, 80, 85, 90, 95, 100), "%"),
    direction = "vertical",
    title_gp  = gpar(fontsize = 10, fontface = "bold")
  ),
  show_parent_dend_line = FALSE,
  use_raster = requireNamespace("magick", quietly = TRUE)
)

# Variant A, without metadata sidebars
ht_fig4_clean <- do.call(Heatmap, c(list(mat_ordered), fig4_common_params))

out_fig4_clean <- file.path(output_dir, "Fig4a_ANI_heatmap_sin_metadatos.png")
png(out_fig4_clean, width = 3400, height = 3400, res = 300)
draw(ht_fig4_clean, merge_legend = TRUE, padding = unit(c(5, 5, 5, 5), "mm"))
dev.off()
cat("  [OK] Guardada:", basename(out_fig4_clean), "\n")

# Variant B, with depth and station sidebars
ht_fig4_meta <- do.call(Heatmap,
  c(list(mat_ordered), fig4_common_params, list(left_annotation = ha_row)))

out_fig4_meta <- file.path(output_dir, "Fig4b_ANI_heatmap_con_metadatos.png")
png(out_fig4_meta, width = 3600, height = 3400, res = 300)
draw(ht_fig4_meta, merge_legend = TRUE, padding = unit(c(5, 5, 5, 5), "mm"))
dev.off()
cat("  [OK] Guardada:", basename(out_fig4_meta), "\n")

# ============================================================
# ---- Fig5: each MAG's contribution to its species ANI ----
# ============================================================
# Assigned species with at least 3 MAGs only.
# One boxplot per MAG, over its ANI against every conspecific.
# An outlier shows up as the lowest median in its species.
# Box colour is the focal MAG's depth, ordered ZMOR, 1000 m, Fondo, OMZ.
# ============================================================
cat("=== Fig5: Contribución por MAG ===\n")

fig5_rows <- list()
for (sp in species_ge3) {
  sp_mags <- meta_df |>
    filter(Species == sp) |>
    pull(NameWithoutExtension)
  sp_mags <- intersect(sp_mags, mag_names)
  if (length(sp_mags) < 3) next
  submat <- mat_vals[sp_mags, sp_mags]

  for (focal in sp_mags) {
    partners <- setdiff(sp_mags, focal)
    ani_vals <- submat[focal, partners]
    ani_vals <- ani_vals[!is.na(ani_vals)]
    if (length(ani_vals) == 0) next

    focal_meta <- meta_df |> filter(NameWithoutExtension == focal)
    fig5_rows[[length(fig5_rows) + 1]] <- data.frame(
      Focal_MAG = focal,
      short_label = focal_meta$short_label,
      Depth_label = as.character(focal_meta$Depth_label),
      Station = focal_meta$Station,
      Species = sp,
      ANI = ani_vals,
      stringsAsFactors = FALSE
    )
  }
}
fig5_data <- bind_rows(fig5_rows) |>
  mutate(Depth_label = factor(Depth_label, levels = depth_factor_levels))

# MAGs sort by median ANI within each species, which puts outliers on the left
mag_order_fig5 <- fig5_data |>
  group_by(Species, Focal_MAG, short_label, Depth_label) |>
  summarise(median_ani = median(ANI), .groups = "drop") |>
  arrange(Species, median_ani)

fig5_data <- fig5_data |>
  mutate(short_fct = factor(
    short_label,
    levels = mag_order_fig5$short_label[order(mag_order_fig5$Species, mag_order_fig5$median_ani)]
  )) |>
  group_by(Species, Focal_MAG) |>
  mutate(median_focal = median(ANI)) |>
  group_by(Species) |>
  mutate(is_outlier = median_focal == min(median_focal)) |>
  ungroup()

# Facet labeller: fmt_taxon, rendered as markdown
sp_labeller <- as_labeller(function(x) fmt_taxon_vec(x, species_parent_lookup))

p_fig5 <- ggplot(fig5_data, aes(x = short_fct, y = ANI)) +
  geom_boxplot(
    aes(fill = Depth_label),
    alpha = 0.6, outlier.shape = NA,
    width = 0.55, linewidth = 0.45
  ) +
  geom_jitter(
    aes(color = Depth_label),
    width = 0.15, size = 1.8, alpha = 0.8, stroke = 0.3
  ) +
  geom_hline(
    yintercept = 95, color = "#C0392B",
    linetype = "dashed", linewidth = 0.7
  ) +
  scale_fill_manual(
    values = depth_colors,
    breaks = depth_factor_levels,
    name   = "Profundidad"
  ) +
  scale_color_manual(
    values = depth_colors,
    breaks = depth_factor_levels,
    guide  = "none"
  ) +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    breaks = c(90, 92, 94, 95, 96, 97, 98, 99, 100)
  ) +
  scale_x_discrete(guide = guide_axis(angle = 50)) +
  facet_wrap(~Species, scales = "free_x", labeller = sp_labeller) +
  labs(x = NULL, y = "ANI (%)") +
  theme_publication(base_size = 12) +
  theme(
    strip.text   = element_markdown(size = 8, face = "bold"),
    axis.text.x  = element_text(size = 7.5)
  )

save_fig(p_fig5, file.path(output_dir, "Fig5_ANI_contribucion_por_MAG.png"),
  width = 14, height = 10
)

# ============================================================
# ---- Fig5b: the same, coloured by station ----
# ============================================================
cat("=== Fig5b: Contribución por MAG — por estación ===\n")

fig5_data <- fig5_data |>
  mutate(Station = factor(Station, levels = station_levels))

p_fig5b <- ggplot(fig5_data, aes(x = short_fct, y = ANI)) +
  geom_boxplot(
    aes(fill = Station),
    alpha = 0.6, outlier.shape = NA,
    width = 0.55, linewidth = 0.45
  ) +
  geom_jitter(
    aes(color = Station),
    width = 0.15, size = 1.8, alpha = 0.8, stroke = 0.3
  ) +
  geom_hline(
    yintercept = 95, color = "#C0392B",
    linetype = "dashed", linewidth = 0.7
  ) +
  scale_fill_manual(values = station_colors, name = "Estación") +
  scale_color_manual(values = station_colors, guide = "none") +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    breaks = c(90, 92, 94, 95, 96, 97, 98, 99, 100)
  ) +
  scale_x_discrete(guide = guide_axis(angle = 50)) +
  facet_wrap(~Species, scales = "free_x", labeller = sp_labeller) +
  labs(x = NULL, y = "ANI (%)") +
  theme_publication(base_size = 12) +
  theme(
    strip.text  = element_markdown(size = 8, face = "bold"),
    axis.text.x = element_text(size = 7.5)
  )

save_fig(p_fig5b, file.path(output_dir, "Fig5b_ANI_contribucion_por_MAG_estacion.png"),
  width = 14, height = 10
)

# ============================================================
# ---- statistics: does provenance shift ANI? ----
# ============================================================
# Two complementary tests, because neither alone is sufficient here:
#
#  (1) Modelo Lineal Mixto global (LMM, lme4::lmer)
#      ANI ~ ctx + (1 | Species), fitted by ML (REML = FALSE).
#      The random intercept per species absorbs the baseline ANI differences
#      between taxa, which otherwise dominate the comparison. Tested by
#      likelihood ratio against the null model without ctx. A positive beta
#      means higher ANI when the provenance is shared.
#
#  (2) Wilcoxon rank-sum per species (Mann-Whitney U)
#      Compara ANI(misma procedencia) vs. ANI(diferente procedencia)
#      run inside each species separately.
#      Nonparametric, so it assumes neither normality nor equal variance.
#      Valid at the small n these species have.
#      Benjamini-Hochberg correction across species.
#      95% CI on the median difference by bootstrap, n = 1999.
#
#  A limitation worth stating plainly:
#      MAG pairs are not independent observations. Each focal MAG appears in
#      several pairs within its species, which inflates the effective n and
#      makes the p-values anticonservative. That is why the interpretation
#      leans on the effect size, the median difference with its bootstrap
#      interval, rather than on significance.
# ============================================================
cat("=== Análisis estadístico: procedencia vs. ANI ===\n")

# ---- Paquete adicional: lme4 --------------------------------
if (!requireNamespace("lme4", quietly = TRUE))
  install.packages("lme4", dependencies = TRUE, quiet = TRUE)
suppressMessages(library(lme4))

# ---- 1. prepare the analysis data ----
# "Distintas" is the reference level in both factors, so beta measures the
# increase in ANI when the provenance is shared.
ani_stat <- sparse_df |>
  filter(intraspecies, Species_ref %in% species_ge2) |>
  transmute(
    ANI,
    Species      = Species_ref,
    same_depth   = same_depth,
    same_station = same_station,
    ctx_depth    = factor(
      if_else(same_depth,   "Misma profundidad",   "Profundidades distintas"),
      levels = c("Profundidades distintas", "Misma profundidad")
    ),
    ctx_station  = factor(
      if_else(same_station, "Misma estación",      "Estaciones distintas"),
      levels = c("Estaciones distintas", "Misma estación")
    )
  )

cat(sprintf("  Pares intra-especie: %d  |  Especies incluidas: %d\n",
            nrow(ani_stat), n_distinct(ani_stat$Species)))

# ---- 2. LMM global ------------------------------------------
fit_lmm_lrt <- function(ctx_var, data) {
  full <- lmer(as.formula(paste0("ANI ~ ", ctx_var, " + (1 | Species)")),
               data = data, REML = FALSE)
  null <- lmer(ANI ~ 1 + (1 | Species), data = data, REML = FALSE)
  lrt  <- anova(null, full)
  data.frame(
    beta  = fixef(full)[2],
    se    = sqrt(diag(vcov(full)))[2],
    chi2  = lrt$Chisq[2],
    df    = as.integer(lrt$Df[2]),
    p_LRT = lrt$`Pr(>Chisq)`[2],
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}

lmm_d <- fit_lmm_lrt("ctx_depth",   ani_stat)
lmm_s <- fit_lmm_lrt("ctx_station", ani_stat)

lmm_results <- data.frame(
  variable = c("Profundidad", "Estacion"),
  rbind(lmm_d, lmm_s),
  stringsAsFactors = FALSE
)
lmm_results$sig <- as.character(cut(lmm_results$p_LRT,
  breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
  labels = c("***", "**", "*", "·", "ns")))

for (i in seq_len(nrow(lmm_results))) {
  cat(sprintf("  LMM %-14s β=%+.3f  SE=%.3f  χ²(%d)=%.2f  p=%s  [%s]\n",
    paste0(lmm_results$variable[i], ":"),
    lmm_results$beta[i], lmm_results$se[i],
    lmm_results$df[i], lmm_results$chi2[i],
    formatC(lmm_results$p_LRT[i], format = "e", digits = 2),
    lmm_results$sig[i]))
}

# ---- 3. Wilcoxon per species, with bootstrap intervals ----
wilcox_species <- function(ctx_var, same_label, variable_label) {
  rows <- list()
  for (sp in species_ge2) {
    sub  <- ani_stat |> filter(Species == sp)
    ctxs <- unique(as.character(sub[[ctx_var]]))
    if (length(ctxs) < 2 || !same_label %in% ctxs) next

    v_same <- sub$ANI[as.character(sub[[ctx_var]]) == same_label]
    v_diff <- sub$ANI[as.character(sub[[ctx_var]]) != same_label]
    if (length(v_same) == 0 || length(v_diff) == 0) next

    wt <- tryCatch(wilcox.test(v_same, v_diff, exact = FALSE),
                   error = function(e) NULL)
    if (is.null(wt)) next

    delta <- median(v_same) - median(v_diff)

    # Bootstrap interval, only where both groups hold at least 2 observations
    if (length(v_same) >= 2 && length(v_diff) >= 2) {
      set.seed(42)
      bd <- replicate(1999,
        median(sample(v_same, replace = TRUE)) -
        median(sample(v_diff, replace = TRUE)))
      ci_lo <- quantile(bd, 0.025)
      ci_hi <- quantile(bd, 0.975)
    } else {
      ci_lo <- ci_hi <- NA_real_
    }

    rows[[sp]] <- data.frame(
      Species      = sp,
      Variable     = variable_label,
      n_same       = length(v_same),
      n_diff       = length(v_diff),
      median_same  = round(median(v_same), 4),
      median_diff  = round(median(v_diff), 4),
      delta_median = round(delta, 4),
      ci_lower     = round(as.numeric(ci_lo), 4),
      ci_upper     = round(as.numeric(ci_hi), 4),
      W_statistic  = wt$statistic,
      p_value      = round(wt$p.value, 6),
      stringsAsFactors = FALSE
    )
  }
  if (length(rows) == 0) return(NULL)
  res      <- bind_rows(rows)
  res$p_adj <- p.adjust(res$p_value, method = "BH")
  res$sig   <- as.character(cut(res$p_adj,
    breaks = c(-Inf, 0.001, 0.01, 0.05, 0.1, Inf),
    labels = c("***", "**", "*", "·", "ns")))
  res
}

res_depth   <- wilcox_species("ctx_depth",   "Misma profundidad", "Profundidad")
res_station <- wilcox_species("ctx_station", "Misma estación",    "Estación")

# ---- 4. Guardar TSVs ----------------------------------------
if (!is.null(res_depth) || !is.null(res_station)) {
  write_tsv(
    bind_rows(res_depth, res_station),
    file.path(output_dir, "stat_ANI_wilcoxon_especie.tsv")
  )
  cat("  [OK] stat_ANI_wilcoxon_especie.tsv\n")
}
write_tsv(lmm_results, file.path(output_dir, "stat_ANI_LMM_global.tsv"))
cat("  [OK] stat_ANI_LMM_global.tsv\n")

# ---- 5. Figura: forest plots + scatter ----------------------
# Panel (a): depth delta per species
# Panel (b): station delta per species
# Panel (c): depth delta against station delta, to see whether they track each other
make_forest <- function(res_df, variable_label) {
  if (is.null(res_df) || nrow(res_df) == 0) return(NULL)

  mag_per_sp_fig7 <- meta_df |>
    filter(Species %in% species_ge2) |>
    count(Species) |>
    rename(n_mags = n)

  res_df <- res_df |>
    left_join(mag_per_sp_fig7, by = "Species") |>
    mutate(
      sp_fmt   = fmt_taxon_vec(Species, species_parent_lookup),
      sp_label = paste0(
        sp_fmt,
        sprintf(
          "<br><span style='font-size:13.5px;color:#555555'>n=%d MAGs</span>",
          n_mags
        )
      ),
      sp_fct   = factor(sp_label, levels = sp_label[order(delta_median)]),
      is_sig   = p_adj < 0.05,
      pt_col   = if_else(is_sig, "#C0392B", "gray55"),
      sig_disp = if_else(is_sig & sig %in% c("*","**","***"), sig, "")
    )

  xlim <- max(abs(c(res_df$ci_lower, res_df$ci_upper,
                    res_df$delta_median, 0.4)), na.rm = TRUE) * 1.25

  ggplot(res_df, aes(x = delta_median, y = sp_fct)) +
    geom_vline(xintercept = 0, color = "gray35",
               linetype = "dashed", linewidth = 0.65) +
    geom_errorbarh(
      data = filter(res_df, !is.na(ci_lower)),
      aes(xmin = ci_lower, xmax = ci_upper, color = pt_col),
      height = 0.28, linewidth = 0.75
    ) +
    geom_point(aes(color = pt_col, shape = is_sig), size = 3.8) +
    geom_text(
      data = filter(res_df, nchar(sig_disp) > 0),
      aes(x = xlim * 0.90, label = sig_disp),
      hjust = 0.5, size = 5, color = "#C0392B", fontface = "bold"
    ) +
    scale_color_identity() +
    scale_shape_manual(values = c("TRUE" = 19, "FALSE" = 21), guide = "none") +
    scale_x_continuous(
      limits = c(-xlim, xlim),
      labels = function(x) sprintf("%.2f%%", x)
    ) +
    labs(
      x        = sprintf(
        "ANI media (misma %s vs. diferente %s)",
        tolower(variable_label), tolower(variable_label)),
      y        = NULL,
      subtitle = variable_label
    ) +
    theme_publication(base_size = 12) +
    theme(
      axis.text.y        = element_markdown(size = 10.5),
      axis.title.y       = element_text(margin = margin(r = 0)),
      panel.grid.major.y = element_line(color = "gray93"),
      panel.grid.major.x = element_line(color = "gray88", linetype = "dotted"),
      plot.subtitle      = element_text(face = "bold", size = 13, hjust = 0,
                                         color = "black"),
      plot.margin        = margin(t = 5, r = 28, b = 5, l = 5)
    )
}

p_fd <- make_forest(res_depth,   "Profundidad")
p_fs <- make_forest(res_station, "Estación")

# Panel (c): depth delta against station delta, per species
# Shows whether the same species respond to both factors.
p_scatter <- NULL
if (!is.null(res_depth) && !is.null(res_station)) {
  delta_both <- inner_join(
    res_depth   |> select(Species, delta_depth   = delta_median,
                          p_adj_depth   = p_adj),
    res_station |> select(Species, delta_station = delta_median,
                          p_adj_station = p_adj),
    by = "Species"
  ) |>
    mutate(
      sp_fmt = fmt_taxon_vec(Species, species_parent_lookup),
      sig_group = case_when(
        p_adj_depth < 0.05 & p_adj_station < 0.05 ~ "Ambas",
        p_adj_depth < 0.05   ~ "Solo profundidad",
        p_adj_station < 0.05 ~ "Solo estación",
        TRUE                 ~ "Ninguna"
      ),
      sig_group = factor(sig_group,
        levels = c("Ambas", "Solo profundidad", "Solo estación", "Ninguna"))
    )

  sig_colors_c <- c(
    "Ambas"            = "#C0392B",
    "Solo profundidad" = "#2471A3",
    "Solo estación"    = "#8E44AD",
    "Ninguna"          = "gray60"
  )

  xlim_c <- max(abs(delta_both$delta_station), 0.3, na.rm = TRUE) * 1.3
  ylim_c <- max(abs(delta_both$delta_depth),   0.3, na.rm = TRUE) * 1.3

  p_scatter <- ggplot(delta_both,
      aes(x = delta_station, y = delta_depth, color = sig_group)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               color = "gray50", linewidth = 0.55) +
    geom_vline(xintercept = 0, linetype = "dashed",
               color = "gray50", linewidth = 0.55) +
    geom_point(size = 4.5, stroke = 0.5) +
    ggrepel::geom_label_repel(
      aes(label = gsub("[*_]", "", sp_fmt)),  # richtext no soportado en ggrepel
      size = 3, fill = "white", label.size = 0.2,
      box.padding = 0.4, point.padding = 0.3,
      max.overlaps = 20, show.legend = FALSE,
      fontface = "italic"   # toda la etiqueta en cursiva como aproximación
    ) +
    scale_color_manual(values = sig_colors_c,
                       name = "Significancia\n(FDR < 0.05)") +
    scale_x_continuous(
      limits = c(-xlim_c, xlim_c),
      labels = function(x) sprintf("%.2f%%", x)
    ) +
    scale_y_continuous(
      limits = c(-ylim_c, ylim_c),
      labels = function(x) sprintf("%.2f%%", x)
    ) +
    labs(
      x = "ANI media (misma estación vs. diferente)",
      y = "ANI media (misma profundidad vs. diferente)",
      subtitle = "Profundidad y estación"
    ) +
    theme_publication(base_size = 12) +
    theme(
      legend.position = "right",
      axis.title.y = element_text(margin = margin(r = -2)),
      plot.margin = margin(t = 5, r = 5, b = 5, l = 2),
      panel.grid.major = element_line(color = "gray88", linetype = "dotted"),
      plot.subtitle = element_text(face = "bold", size = 13, hjust = 0,
                                    color = "black")
    )
}

# Assemble the statistics figure
if (!is.null(p_fd) || !is.null(p_fs)) {
  top_panels <- Filter(Negate(is.null), list(p_fd, p_fs))
  top_row    <- if (length(top_panels) == 2) top_panels[[1]] | top_panels[[2]]
                else top_panels[[1]]

  cap_txt <- sprintf(
    paste0(
      "LMM global (lme4) — ",
      "Profundidad: b=%+.3f, chi2(%d)=%.2f, p=%s %s | ",
      "Estacion: b=%+.3f, chi2(%d)=%.2f, p=%s %s\n",
      "Puntos rojos/asteriscos (*): Wilcoxon FDR-adj p < 0.05 (BH). ",
      "Barras: IC 95%% bootstrap (n = 1 999). ",
      "Limitacion: pares de MAGs no son independientes (comparten individuos)."
    ),
    lmm_results$beta[1], lmm_results$df[1], lmm_results$chi2[1],
    formatC(lmm_results$p_LRT[1], format = "e", digits = 2), lmm_results$sig[1],
    lmm_results$beta[2], lmm_results$df[2], lmm_results$chi2[2],
    formatC(lmm_results$p_LRT[2], format = "e", digits = 2), lmm_results$sig[2]
  )

  n_sp_rows <- max(
    if (!is.null(res_depth))   nrow(res_depth)   else 0,
    if (!is.null(res_station)) nrow(res_station) else 0
  )
  fig7_h_top <- max(5, n_sp_rows * 0.65 + 2)

  if (!is.null(p_scatter)) {
    # Scatter at about 70% of the top row width, centred
    scatter_row <- patchwork::wrap_plots(
      patchwork::plot_spacer(), p_scatter, patchwork::plot_spacer(),
      nrow = 1, widths = c(0.15, 0.70, 0.15)
    )
    p_stat_fig <- (top_row / scatter_row) +
      plot_layout(heights = c(fig7_h_top, 4)) +
      plot_annotation(caption = cap_txt, tag_levels = "a") &
      theme(
        plot.tag     = element_text(face = "bold", size = 20, color = "black"),
        plot.caption = element_text(size = 7.5, color = "gray40",
                                    hjust = 0, lineheight = 1.3)
      )
    save_fig(p_stat_fig,
      file.path(output_dir, "Fig7_estadistica_procedencia.png"),
      width = 15.5, height = fig7_h_top + 4.5)
  } else {
    p_stat_fig <- top_row +
      plot_annotation(caption = cap_txt, tag_levels = "a") &
      theme(
        plot.tag     = element_text(face = "bold", size = 20, color = "black"),
        plot.caption = element_text(size = 7.5, color = "gray40",
                                    hjust = 0, lineheight = 1.3)
      )
    save_fig(p_stat_fig,
      file.path(output_dir, "Fig7_estadistica_procedencia.png"),
      width = 15.5, height = fig7_h_top)
  }
  cat("  [OK] Fig7_estadistica_procedencia.png\n")
}

# ============================================================
# ---- Fig6: unassigned MAGs against their genus ----
# Disabled behind if (FALSE). Remove the guard and its closing brace to run it.
# ============================================================
if (FALSE) {
cat("=== Fig6: MAGs noveles en contexto de género ===\n")

novel_genera <- c("MGIIb-Q1", "Nitrosopelagicus", "UBA36")

fig6_pairs_list <- list()
for (gn in novel_genera) {
  gn_pairs <- sparse_df |>
    filter(Genus_ref == gn, Genus_query == gn) |>
    mutate(
      pair_cat = case_when(
        Species_ref == "N/A" & Species_query == "N/A" ~ "Novel vs. novel",
        Species_ref == "N/A" | Species_query == "N/A" ~ "Novel vs. asignado",
        TRUE ~ "Asignado vs. asignado"
      ),
      Genus_label = gn
    ) |>
    select(ANI, AF_ref, pair_cat, Genus_label)
  fig6_pairs_list[[gn]] <- gn_pairs
}

fig6_pairs <- bind_rows(fig6_pairs_list) |>
  mutate(
    pair_cat_fct = factor(
      pair_cat,
      levels = c("Asignado vs. asignado", "Novel vs. asignado", "Novel vs. novel")
    ),
    # Format the genus name in the facet strip
    Genus_label_fmt = taxon_italic(Genus_label, level = "Genus"),
    Genus_fct = factor(Genus_label_fmt, levels = unique(taxon_italic(novel_genera, level = "Genus")))
  )

pair_cat_colors <- c(
  "Asignado vs. asignado" = "#2471A3",
  "Novel vs. asignado"    = "#D68910",
  "Novel vs. novel"       = "#CA6F1E"
)

# Genus labeller; the markdown is already applied, so this just passes it through
genus_labeller <- as_labeller(function(x) x)

p_fig6a <- ggplot(fig6_pairs, aes(x = pair_cat_fct, y = ANI, fill = pair_cat_fct)) +
  geom_boxplot(
    alpha = 0.6, outlier.shape = NA,
    width = 0.5, linewidth = 0.5
  ) +
  geom_jitter(
    aes(color = pair_cat_fct),
    width = 0.15, size = 2.5, alpha = 0.85
  ) +
  geom_hline(
    yintercept = 95, color = "#C0392B",
    linetype = "dashed", linewidth = 0.8
  ) +
  scale_fill_manual(values = pair_cat_colors, guide = "none") +
  scale_color_manual(values = pair_cat_colors, guide = "none") +
  scale_y_continuous(
    labels = function(x) paste0(x, "%"),
    breaks = c(75, 80, 85, 90, 95, 97, 99, 100)
  ) +
  scale_x_discrete(guide = guide_axis(angle = 20), drop = TRUE) +
  facet_wrap(~Genus_label_fmt,
    ncol = 3, scales = "free_x",
    labeller = genus_labeller
  ) +
  labs(x = NULL, y = "ANI (%)") +
  theme_publication(base_size = 12) +
  theme(
    strip.text  = element_markdown(size = 11, face = "bold"),
    axis.text.x = element_text(size = 9)
  )

# Panel B: heatmap per genus, novel and assigned MAGs together
make_genus_heatmap <- function(genus_val) {
  gn_mags <- meta_df |>
    filter(Genus == genus_val) |>
    arrange(Depth_label, Station) |>
    pull(NameWithoutExtension) |>
    intersect(mag_names)
  if (length(gn_mags) < 2) {
    return(NULL)
  }

  submat <- mat_vals[gn_mags, gn_mags]
  lbl_map <- setNames(meta_df$short_label, meta_df$NameWithoutExtension)
  sp_map <- setNames(meta_df$Species, meta_df$NameWithoutExtension)

  long_df <- as.data.frame(as.table(submat)) |>
    rename(MAG_x = Var1, MAG_y = Var2, ANI = Freq) |>
    mutate(
      MAG_x      = factor(MAG_x, levels = gn_mags),
      MAG_y      = factor(MAG_y, levels = gn_mags),
      Lbl_x      = lbl_map[as.character(MAG_x)],
      Lbl_y      = lbl_map[as.character(MAG_y)],
      is_novel_x = sp_map[as.character(MAG_x)] == "N/A",
      is_novel_y = sp_map[as.character(MAG_y)] == "N/A",
      # "<AF" rather than a blank: skani dropped this pair for low aligned
      # fraction, which is an exclusion rather than missing data.
      ANI_lab    = if_else(!is.na(ANI), sprintf("%.1f", ANI), "<AF")
    )

  off_diag <- long_df$ANI[as.character(long_df$MAG_x) != as.character(long_df$MAG_y)]
  ani_min <- floor(min(off_diag, na.rm = TRUE) * 10) / 10

  # Text colour: grey where the cell reads <AF, black where it carries an ANI.
  # Precomputed as a vector, which avoids clashing with the colour scale used by
  # the tile borders, which map interaction(novel_x, novel_y).
  lab_col <- if_else(is.na(long_df$ANI), "gray45", "black")

  ggplot(long_df, aes(x = Lbl_x, y = Lbl_y, fill = ANI)) +
    geom_tile(aes(color = interaction(is_novel_x, is_novel_y)), linewidth = 0.8) +
    geom_text(aes(label = ANI_lab), size = 2.8, color = lab_col) +
    scale_fill_gradient2(
      low = "#FFF5F0", mid = "#FC8D59", high = "#67000D",
      midpoint = 95, limits = c(ani_min, 100),
      na.value = "gray88", guide = "none"
    ) +
    scale_color_manual(
      values = c(
        "FALSE.FALSE" = "white", "TRUE.FALSE" = "#D68910",
        "FALSE.TRUE" = "#D68910", "TRUE.TRUE" = "#CA6F1E"
      ),
      guide = "none"
    ) +
    scale_x_discrete(guide = guide_axis(angle = 40)) +
    scale_y_discrete(limits = rev) +
    labs(title = paste0(taxon_italic(genus_val, level = "Genus")), x = NULL, y = NULL) +
    theme_publication(base_size = 9) +
    theme(
      plot.title = element_markdown(size = 9, face = "bold"),
      aspect.ratio = 1, axis.text = element_text(size = 7.5)
    )
}

genus_ht_plots <- lapply(novel_genera, make_genus_heatmap)
names(genus_ht_plots) <- novel_genera
genus_ht_plots <- Filter(Negate(is.null), genus_ht_plots)

p_fig6b <- wrap_plots(genus_ht_plots, ncol = 3)

p_fig6 <- (p_fig6a / wrap_elements(p_fig6b)) +
  plot_layout(heights = c(1, 1.1)) +
  plot_annotation(tag_levels = "a", caption = ani_caption) &
  theme(
    plot.tag     = element_text(face = "bold", size = 20, color = "black"),
    plot.caption = element_text(size = 9, color = "gray35",
                                 hjust = 0, lineheight = 1.2)
  )

save_fig(p_fig6, file.path(output_dir, "Fig6_ANI_MAGs_noveles.png"),
  width = 13, height = 12.5
)
} # fin if (FALSE) — Fig6

# ============================================================
# Resumen final
# ============================================================
cat("\n=== Resumen de resultados ===\n")
cat("Directorio de salida:", output_dir, "\n\n")

output_files <- list(
  "Fig1_ANI_distribucion_intraespecie.png"       = "Distribución ANI por especie (asignadas)",
  "Fig1b_ANI_distribucion_grupos_NA.png"         = "Distribución ANI por grupo N/A",
  "Fig1c_ANI_distribucion_panel.png"             = "Panel combinado: Fig1 (a) + Fig1b (b)",
  "Fig2_ANI_heatmaps_por_especie.png"            = "Heatmaps por especie (asignadas, ≥ 3 MAGs)",
  "Fig2b_ANI_heatmaps_grupos_NA.png"             = "Heatmaps de grupos N/A (≥ 3 MAGs)",
  "Fig3_ANI_por_profundidad.png"                 = "ANI por profundidad (asignadas)",
  "Fig3b_ANI_por_profundidad_NA.png"             = "ANI por profundidad (grupos N/A)",
  "Fig3c_ANI_por_estacion.png"                   = "ANI por estación (asignadas)",
  "Fig4a_ANI_heatmap_sin_metadatos.png"          = "Heatmap 54×54 — sin barras de metadatos",
  "Fig4b_ANI_heatmap_con_metadatos.png"          = "Heatmap 54×54 — con barras de metadatos",
  "Fig5_ANI_contribucion_por_MAG.png"            = "Contribución individual — por profundidad",
  "Fig5b_ANI_contribucion_por_MAG_estacion.png"  = "Contribución individual — por estación",
  "Fig7_estadistica_procedencia.png"             = "Análisis estadístico — procedencia vs ANI",
  "stat_ANI_wilcoxon_especie.tsv"                = "Wilcoxon FDR por especie (profundidad + estación)",
  "stat_ANI_LMM_global.tsv"                      = "LMM global — beta, chi2, p-LRT"
)

for (fn in names(output_files)) {
  full_path_out <- file.path(output_dir, fn)
  marker <- if (file.exists(full_path_out)) "[OK]" else "○ (no generada)"
  cat(sprintf("  %s  %-46s  %s\n", marker, fn, output_files[[fn]]))
}

cat("\n[OK] Script completado:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")
