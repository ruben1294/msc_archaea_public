#!/usr/bin/env Rscript
# Script: tax_analysis.R
# Description: Taxonomic analysis of the archaeal MAGs classified with GTDB-Tk
#              classify_wf: bar charts and heatmaps at species, genus and
#              family level, plus novelty, ANI and biogeography figures.
# Created: 2026-03-23
# Version: 1.0

### To run in RStudio: use Source, or Ctrl+Shift+Enter.
###
# To run from the terminal:
# Rscript scripts/R/tax_analysis.R

# ---- project root ----
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

# Shared fonts, theme and palettes, and resolve_project_dir().
source(file.path(PROJECT_DIR, "scripts", "R", "_pub_style.R"))

PROJECT_ROOT <- resolve_project_dir()

# === Input and output paths ===
input_gtdb_summary <- file.path(PROJECT_ROOT, "results/03_taxonomy_assigned/classify/tax_.ar53.summary.tsv")
input_metadata <- file.path(PROJECT_ROOT, "data", "MAGs_metadata.tsv")
output_dir <- file.path(PROJECT_ROOT, "results/SO_01/tax")

# ---- packages ----
# Missing packages install themselves. If one keeps failing, install it by hand
# from the RStudio console:
#   install.packages("package_name")
# ComplexHeatmap comes from Bioconductor:
#   if (!requireNamespace("BiocManager")) install.packages("BiocManager")
#   BiocManager::install("ComplexHeatmap")

cat("\n=== Cargando paquetes ===\n")

# Version-sensitive dependencies are checked and updated before anything that
# needs them is loaded. Each block stands alone, so if RStudio restarts the
# session mid-install, rerunning the script finds that dependency already
# current and skips straight past it.

# ggtext needs xfun 0.54 or newer
if (requireNamespace("xfun", quietly = TRUE) && packageVersion("xfun") < "0.54") {
  cat(
    "[WARN] Updating xfun (installed:", as.character(packageVersion("xfun")),
    "- need >= 0.54)...\n"
  )
  install.packages("xfun", quiet = TRUE)
  cat("[OK] xfun updated. Run the script again.\n")
  stop("Restart the session or rerun the script to continue.")
}

# patchwork 1.3 needs ggplot2 3.5.2 or newer
if (requireNamespace("ggplot2", quietly = TRUE) && packageVersion("ggplot2") < "3.5.2") {
  cat(
    "[WARN] Updating ggplot2 (installed:", as.character(packageVersion("ggplot2")),
    "- need >= 3.5.2)...\n"
  )
  install.packages("ggplot2", quiet = TRUE)
  cat("[OK] ggplot2 actualizado. Vuelve a ejecutar el script.\n")
  stop("Restart the session or rerun the script to continue.")
}

required_packages <- c(
  "ggplot2", # must load, up to date, before patchwork
  "tidyverse", # data handling and plotting
  "patchwork", # combining plots
  "scales", # plot scales
  "ggrepel", # better text placement
  "RColorBrewer", # colour palettes
  "ggtext", # rich text (markdown/HTML) in ggplot2
  "ComplexHeatmap", # heatmaps (DRAM2 style)
  "circlize", # colorRamp2 for heatmap colour scales
  "grid", # grid graphics (needed by ComplexHeatmap)
  "readxl", # reading Excel files (metabolism_summary.xlsx)
  "showtext", # custom font support
  "sysfonts" # system font detection
)

for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    cat("[WARN] Instalando:", pkg, "...\n")
    tryCatch(
      install.packages(pkg, dependencies = TRUE, quiet = TRUE),
      error = function(e) {
        stop(
          "No se pudo instalar '", pkg, "'. ",
          "Instálalo manualmente desde la consola de RStudio:\n  install.packages('",
          pkg, "')\n\nError original: ", conditionMessage(e)
        )
      }
    )
  }
  tryCatch(
    suppressMessages(library(pkg, character.only = TRUE, quietly = TRUE)),
    error = function(e) {
      stop(
        "No se pudo cargar '", pkg, "'. ",
        "Intenta reinstalarlo:\n  install.packages('", pkg, "')\n\n",
        "Error original: ", conditionMessage(e)
      )
    }
  )
  cat("[OK]", pkg, "\n")
}

cat("\n[OK] Todos los paquetes cargados correctamente\n\n")

cat("=== Configurando fuentes ===\n")
showtext_auto()
showtext_opts(dpi = 300)

font_family <- "sans"
available_fonts <- sysfonts::font_files()

font_priority <- c(
  "Source Sans Pro", "SourceSansPro",
  "CMU Serif", "CMU Sans Serif", "Computer Modern", "Latin Modern",
  "Helvetica", "Helvetica Neue"
)

# Load local fonts directly (priority: Lato, then Source Sans 3)
font_found <- FALSE
local_font_dir <- path.expand("~/.local/share/fonts")
local_direct <- list(
  list(
    name = "Lato",
    regular = "Lato-Regular.ttf",
    bold = "Lato-Bold.ttf",
    italic = "Lato-Italic.ttf",
    bolditalic = "Lato-BoldItalic.ttf"
  ),
  list(
    name = "Source Sans 3",
    regular = "SourceSans3-Regular.otf",
    bold = "SourceSans3-Bold.otf",
    italic = "SourceSans3-It.otf",
    bolditalic = "SourceSans3-BoldIt.otf"
  )
)
if (dir.exists(local_font_dir)) {
  for (lf in local_direct) {
    reg_path <- file.path(local_font_dir, lf$regular)
    if (file.exists(reg_path)) {
      tryCatch(
        {
          bld_path <- file.path(local_font_dir, lf$bold)
          ita_path <- file.path(local_font_dir, lf$italic)
          bia_path <- file.path(local_font_dir, lf$bolditalic)
          sysfonts::font_add(
            family     = "custom_font",
            regular    = reg_path,
            bold       = if (file.exists(bld_path)) bld_path else NULL,
            italic     = if (file.exists(ita_path)) ita_path else NULL,
            bolditalic = if (file.exists(bia_path)) bia_path else NULL
          )
          font_family <- "custom_font"
          cat("[OK] Usando fuente:", lf$name, "(local ~/.local/share/fonts)\n")
          font_found <- TRUE
          break
        },
        error = function(e) {
          cat("[WARN] Archivos de", lf$name, "encontrados pero no se pudieron cargar:", conditionMessage(e), "\n")
        }
      )
    }
  }
}

# Fallback: look through the system font registry
for (font_name in font_priority) {
  matching_fonts <- available_fonts[grepl(font_name, available_fonts$family, ignore.case = TRUE), ]
  if (nrow(matching_fonts) > 0) {
    tryCatch(
      {
        font_path <- matching_fonts$file[1]
        sysfonts::font_add(family = "custom_font", regular = font_path)
        font_family <- "custom_font"
        cat("[OK] Usando fuente:", font_name, "\n")
        font_found <- TRUE
        break
      },
      error = function(e) {}
    )
  }
}
if (!font_found) {
  cat("[WARN] No se encontró fuente del sistema. Intentando descargar de Google Fonts...\n")
  google_priority <- list(
    list(name = "Source Sans Pro"),
    list(name = "Lato"),
    list(name = "Open Sans")
  )
  for (gf in google_priority) {
    tryCatch(
      {
        sysfonts::font_add_google(gf$name, "custom_font")
        font_family <- "custom_font"
        cat("[OK] Descargada y usando:", gf$name, "(Google Fonts)\n")
        font_found <- TRUE
        break
      },
      error = function(e) {}
    )
  }
  if (!font_found) {
    cat("[WARN] No se pudo descargar de Google Fonts. Usando fuente sans-serif por defecto.\n")
  }
}
cat("[OK] Configuración de fuentes completada\n\n")

# === Configuration ===
cat("=== Configuration ===\n")
cat("GTDB-Tk summary:", input_gtdb_summary, "\n")
cat("MAGs metadata:  ", input_metadata, "\n")
cat("Output dir:     ", output_dir, "\n\n")

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

script_start <- Sys.time()
cat(sprintf("=== Inicio del script: %s ===\n\n", format(script_start, "%Y-%m-%d %H:%M:%S")))

# === Publication theme (matches MAG_quality_assessment.R) ===
theme_publication <- function(base_size = 16) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor        = element_blank(),
      panel.grid.major.y      = element_line(color = "gray", linetype = "dotted"),
      panel.grid.major.x      = element_blank(),
      axis.text               = element_text(color = "black", size = base_size - 1),
      axis.title              = element_text(face = "bold", size = base_size + 2),
      legend.position         = "right",
      legend.background       = element_rect(fill = "white", color = NA),
      legend.title            = element_text(face = "bold", size = base_size),
      legend.text             = element_text(size = base_size - 3),
      strip.background        = element_rect(fill = "gray90", color = "black"),
      strip.text              = element_text(face = "bold", size = base_size),
      plot.title              = element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      text                    = element_text(family = font_family)
    )
}

# === Color palettes ===
# Depth
depth_levels <- c("ROMZ", "1000", "BTM")
depth_labels <- c("ROMZ" = "ZMOR", "1000" = "1000 m", "BTM" = "Fondo")
depth_factor_levels <- c("ZMOR", "1000 m", "Fondo")
depth_colors <- c("ZMOR" = "#4A90D9", "1000 m" = "#E67E22", "Fondo" = "#8E44AD")

# Stations
station_levels <- c("A4", "A6", "C22", "F37", "PO1", "Y6")
station_colors <- RColorBrewer::brewer.pal(6, "Set1")
names(station_colors) <- station_levels

# Phyla
phylum_palette <- c(
  "Thermoplasmatota"  = "#E74C3C",
  "Thermoproteota"    = "#2ECC71",
  "Undinarchaeota"    = "#3498DB",
  "Otro"              = "#C8974D"
)

# ---- === Data Loading and Processing ----
cat("=== Loading data ===\n")

# Load GTDB-Tk summary
if (!file.exists(input_gtdb_summary)) {
  stop(
    "Archivo no encontrado: ", input_gtdb_summary,
    "\nVerifica que el directorio de trabajo sea la raíz del proyecto (getwd())."
  )
}
gtdb_raw <- read_tsv(input_gtdb_summary, show_col_types = FALSE)

# Load MAGs metadata
if (!file.exists(input_metadata)) {
  stop(
    "Archivo no encontrado: ", input_metadata,
    "\nVerifica que el directorio de trabajo sea la raíz del proyecto (getwd())."
  )
}
meta_raw <- read_tsv(input_metadata, show_col_types = FALSE)

cat("[OK] Loaded GTDB-Tk summary:", nrow(gtdb_raw), "MAGs\n")
cat("[OK] Loaded metadata:", nrow(meta_raw), "MAGs\n\n")

# === Parse GTDB taxonomy string ===
# Format: d__Archaea;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species
parse_gtdb_taxonomy <- function(classification) {
  parse_level <- function(str, prefix) {
    pattern <- paste0(prefix, "([^;]*)")
    m <- regmatches(str, regexpr(pattern, str))
    if (length(m) == 0 || m == "") {
      return(NA_character_)
    }
    val <- sub(paste0("^", prefix), "", m)
    if (val == "") {
      return(NA_character_)
    }
    return(val)
  }
  tibble(
    Domain  = parse_level(classification, "d__"),
    Phylum  = parse_level(classification, "p__"),
    Class   = parse_level(classification, "c__"),
    Order   = parse_level(classification, "o__"),
    Family  = parse_level(classification, "f__"),
    Genus   = parse_level(classification, "g__"),
    Species = parse_level(classification, "s__")
  )
}

tax_parsed <- gtdb_raw %>%
  rowwise() %>%
  mutate(parse_gtdb_taxonomy(classification)) %>%
  ungroup()

# Align MAG names. GTDB-Tk reports the names it was given: the original
# underscores (A4_1000_megahit_bin.28) in the thesis run, dots when it runs after
# utils/05_rename_mags.sh as the pipeline now has it. Either form maps onto the
# metadata OriginalName through NameWithoutExtension, and the join uses that.
df <- tax_parsed %>%
  rename(OriginalName = user_genome) %>%
  mutate(OriginalName = coalesce(
    meta_raw$OriginalName[match(gsub("_", ".", OriginalName, fixed = TRUE),
                                meta_raw$NameWithoutExtension)],
    OriginalName
  )) %>%
  left_join(
    meta_raw %>% select(OriginalName, NameWithoutExtension, Station, Depth),
    by = "OriginalName"
  ) %>%
  mutate(
    # Factor ordering
    Depth = factor(Depth, levels = depth_levels),
    Station = factor(Station, levels = station_levels),
    # Novelty flags (TRUE = no assignment at that level)
    Novel_Species = is.na(Species),
    Novel_Genus = is.na(Genus),
    Novel_Family = is.na(Family),
    # ANI to closest genome (N/A to NA)
    ANI_closest = suppressWarnings(as.numeric(
      ifelse(closest_genome_ani %in% c("N/A", ""), NA, closest_genome_ani)
    )),
    AF_closest = suppressWarnings(as.numeric(
      ifelse(closest_genome_af %in% c("N/A", ""), NA, closest_genome_af)
    )),
    # Phylum for coloring
    Phylum_group = case_when(
      Phylum == "Thermoplasmatota" ~ "Thermoplasmatota",
      Phylum == "Thermoproteota" ~ "Thermoproteota",
      Phylum == "Undinarchaeota" ~ "Undinarchaeota",
      !is.na(Phylum) ~ "Otro",
      TRUE ~ "Sin asignar"
    ),
    # Depth labels for display
    Depth_label = recode(as.character(Depth), !!!depth_labels)
  )

cat("[OK] Parsed taxonomy for", nrow(df), "MAGs\n")
cat("[OK] Matched to metadata:", sum(!is.na(df$Station)), "with station info\n\n")

# ---- === Helpers ----

# CDC scientific nomenclature (https://wwwnc.cdc.gov/eid/page/scientific-nomenclature):
# only genus and species (and infraspecific ranks) are italicised; ranks above
# genus (family, order, class, phylum, domain) are plain text.
ITALIC_LEVELS <- c("Genus", "Species", "Subspecies", "Strain")

should_italicize_level <- function(level) {
  !is.null(level) && length(level) == 1 && !is.na(level) && level %in% ITALIC_LEVELS
}

# Format a taxon name for display.
#   - level in {Genus, Species, ...} to italicise per-word, with established
#     words italic and GTDB codes (digits / all-uppercase like UBA60, MGIII,
#     sp002449325) kept plain so binomials render as
#       "Thalassarchaeum betae"        to "*Thalassarchaeum betae*"
#       "Thalassarchaeum sp002507125"  to "*Thalassarchaeum* sp002507125"
#   - level in {Family, Order, Class, Phylum, Domain} to returned as plain text
#     per CDC rules (e.g. "Poseidoniales" stays plain).
#   - NA / "Sin asignar" to returned unchanged.
taxon_italic <- function(x, level = "Species") {
  if (!should_italicize_level(level)) {
    return(ifelse(is.na(x), x, as.character(x)))
  }
  sapply(x, function(nm) {
    if (is.na(nm) || nm == "Sin asignar") {
      return(nm)
    }
    words <- strsplit(nm, " ")[[1]]
    is_plain <- grepl("[0-9]", words) | grepl("^[A-Z]+$", words)
    # Group consecutive words sharing the same italic/plain flag into one chunk
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

save_fig <- function(plot, filename, width = 10, height = 7) {
  png(file.path(output_dir, filename),
    width = width, height = height, units = "in", res = 300
  )
  print(plot)
  dev.off()
  cat("[OK] Saved:", filename, "\n")
}

# Parent hierarchy: for each taxonomic level, the column one step above
parent_hierarchy <- list(
  Species = "Genus",
  Genus   = "Family",
  Family  = "Order",
  Order   = "Class",
  Class   = "Phylum",
  Phylum  = "Domain"
)

# Abbreviations shown before the parent name: "(o. Nitrososphaerales)"
level_abbrev <- list(
  Genus   = "g.",
  Family  = "f.",
  Order   = "o.",
  Class   = "c.",
  Phylum  = "p.",
  Domain  = "d."
)

# Build a lookup: child_taxon_value to list(name, level, abbrev)
#
# Climbs the hierarchy until it finds an ancestor whose name differs from the
# child's. For Species, the comparison key is the first word (genus epithet in
# binomial nomenclature), so "UBA57 sp002494565" is compared as "UBA57".
# For all other levels the full name is used.
# If the immediate parent is identical, the function keeps climbing (e.g.
# UBA57 genus to UBA57 family to Nitrososphaerales order to stops here).
get_parent_lookup <- function(df_in, child_col) {
  # Ordered list of ancestor column names
  ancestors <- character(0)
  cur <- child_col
  repeat {
    par <- parent_hierarchy[[cur]]
    if (is.null(par)) break
    ancestors <- c(ancestors, par)
    cur <- par
  }
  # Only keep ancestors that exist as columns in df_in
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

    # For species: key is the first word (genus part of binomial name)
    child_key <- if (child_col == "Species") {
      strsplit(cv, " ")[[1]][1]
    } else {
      cv
    }

    for (anc_col in ancestors) {
      anc_vals <- na.omit(rows[[anc_col]])
      if (length(anc_vals) == 0) next
      anc_val <- names(sort(table(anc_vals), decreasing = TRUE))[1]
      # Accept only if the ancestor name differs from the child
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

# Format a single taxon name with parent in parentheses if it's an acronym.
# Acronyms: first word is all-uppercase (MGIII, DPANN) or contains digits (UBA1234).
# Italicisation follows CDC rules per level (taxon and parent each evaluated
# independently). Parent is prefixed with its level abbreviation: "(o. Nitrososphaerales)".
fmt_taxon <- function(nm, parent_lookup = NULL, level = "Species") {
  if (is.na(nm) || nm == "Sin asignar") {
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

# Vectorised version suitable for scale labels= argument
fmt_taxon_vec <- function(x, parent_lookup = NULL, level = "Species") {
  sapply(x, fmt_taxon,
    parent_lookup = parent_lookup, level = level, USE.NAMES = FALSE
  )
}

# Factory: returns a label function capturing a specific parent_lookup and
# the taxonomic level of the child column being labelled.
make_label_fn <- function(parent_lookup, level = "Species") {
  function(x) fmt_taxon_vec(x, parent_lookup, level = level)
}

# ---- ---- section 1: bar charts, every combination per taxonomic rank ----
# 4 chart types × 3 taxon levels = 12 bar charts
#
# Type A: x = station, fill = taxon
# Type B: x = depth,   fill = taxon
# Type C: x = station, facet by depth, fill = taxon
# Type D: x = depth,   facet by station, fill = taxon
cat("\n=== Sección 1: Gráficos de barras por nivel taxonómico ===\n")

taxon_levels <- list(
  list(col = "Phylum", label = "Filo", prefix = "filo"),
  list(col = "Order", label = "Orden", prefix = "orden"),
  list(col = "Family", label = "Familia", prefix = "familia"),
  list(col = "Genus", label = "Género", prefix = "genero"),
  list(col = "Species", label = "Especie", prefix = "especie")
)

# Replace gray-ish colors (HSV saturation < 0.15) in a palette vector with
# saturated alternatives drawn from Set1[1:8] (all fully saturated).
# Works on named or unnamed vectors; preserves names.
drop_palette_grays <- function(pal) {
  sat <- rgb2hsv(col2rgb(pal))["s", ]
  gray_idx <- which(sat < 0.15)
  if (length(gray_idx) == 0) {
    return(pal)
  }
  backup <- RColorBrewer::brewer.pal(8, "Set1") # Set1[1:8]: no gray
  candidates <- backup[!backup %in% pal]
  for (i in seq_along(gray_idx)) {
    if (i <= length(candidates)) pal[gray_idx[i]] <- candidates[i]
  }
  pal
}

# Build a reasonably-sized discrete palette for taxa, free of gray tones.
make_taxon_palette <- function(vals) {
  vals <- sort(unique(na.omit(vals)))
  n <- length(vals)
  if (n <= 8) {
    pal <- RColorBrewer::brewer.pal(max(n, 3), "Set2")[seq_len(n)]
  } else if (n <= 12) {
    pal <- RColorBrewer::brewer.pal(n, "Set3")
  } else {
    pal <- colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n)
  }
  pal <- drop_palette_grays(pal)
  # Add amber/tan for NA/unknown (distinct from gray strip backgrounds)
  clrs <- c(setNames(pal, vals), "Sin asignar" = "#C8974D")
  return(clrs)
}

# Canonical phylum order used by the phylum-coherent palette. Stacking goes
# top to bottom in this order, so Thermoplasmatota appears at the top of each bar,
# Undinarchaeota at the bottom, the order the earlier palette gave the phylum
# chart alphabetically.
phylum_levels <- c("Thermoplasmatota", "Thermoproteota", "Undinarchaeota")

# Lightness adjustment via HSV. factor < 1 darkens; factor > 1 lightens by
# raising V toward 1 and reducing S so the result stays a recognisable shade
# of the original hue rather than washed-out gray.
adjust_lightness <- function(hex, factor) {
  hsv_val <- rgb2hsv(col2rgb(hex))
  h <- hsv_val["h", 1]; s <- hsv_val["s", 1]; v <- hsv_val["v", 1]
  if (factor < 1) {
    v <- v * factor
  } else if (factor > 1) {
    v <- min(1, v + (1 - v) * (factor - 1))
    s <- s / factor
  }
  hsv(h, max(0, min(1, s)), max(0, min(1, v)))
}

# Build a coordinated (palette, factor-levels) pair for charts at the Phylum
# or Order level so they share colors and stacking order. Children inherit
# their parent phylum's color from `base_palette`; phyla with multiple
# children get progressively lighter shades. Stacking sequence follows
# `phyla`, alphabetical within each phylum.
build_phylum_coherent <- function(df_in, child_col,
                                  base_palette = phylum_palette,
                                  phyla = phylum_levels) {
  if (child_col == "Phylum") {
    children_by_phylum <- setNames(lapply(phyla, identity), phyla)
    extra <- setdiff(unique(na.omit(df_in[[child_col]])), phyla)
  } else {
    map_df <- df_in %>%
      filter(!is.na(.data[[child_col]]), !is.na(Phylum)) %>%
      distinct(across(all_of(c(child_col, "Phylum"))))
    children_by_phylum <- lapply(split(map_df[[child_col]], map_df$Phylum), sort)
    # Reverse children within each phylum so the alphabetically last child
    # gets the base (darkest) color and earlier children get lighter shades.
    # In practice this swaps Poseidoniales (base red) with MGIII (lighter red)
    # within Thermoplasmatota; phyla with a single child are unaffected.
    children_by_phylum <- lapply(children_by_phylum, rev)
    extra <- unique(map_df[[child_col]][!map_df$Phylum %in% phyla])
  }

  pal <- character(0)
  ordered_children <- character(0)
  for (ph in phyla) {
    children <- children_by_phylum[[ph]]
    if (is.null(children) || length(children) == 0) next
    base_col <- if (ph %in% names(base_palette)) base_palette[[ph]] else "#888888"
    if (length(children) == 1) {
      pal[children[1]] <- base_col
    } else {
      shade_factors <- seq(1.0, 1.45, length.out = length(children))
      for (i in seq_along(children)) {
        pal[children[i]] <- adjust_lightness(base_col, shade_factors[i])
      }
    }
    ordered_children <- c(ordered_children, children)
  }
  if (length(extra) > 0) pal[extra] <- "#C8974D"
  pal["Sin asignar"] <- "#C8974D"
  list(
    palette = pal,
    levels  = unique(c(ordered_children, extra, "Sin asignar"))
  )
}

df_bar <- df %>%
  filter(!is.na(Station), !is.na(Depth)) %>%
  mutate(Depth_label = factor(Depth_label, levels = depth_factor_levels))

for (tl in taxon_levels) {
  col_name <- tl$col
  col_label <- tl$label
  prefix <- tl$prefix

  # Recode NA to "Sin asignar"
  df_bar_tl <- df_bar %>%
    mutate(taxon_val = ifelse(is.na(.data[[col_name]]),
      "Sin asignar", .data[[col_name]]
    ))

  # Phylum and Order share a coordinated palette/stacking so the legend and
  # vertical layout of Fig1C_filo and Fig1C_orden line up: each order takes
  # the color of its parent phylum (lighter shades when a phylum has more
  # than one order), and orders stack in the same phylum sequence.
  if (col_name %in% c("Phylum", "Order")) {
    pc <- build_phylum_coherent(df_bar_tl, col_name)
    pal <- pc$palette
    df_bar_tl <- df_bar_tl %>%
      mutate(taxon_val = factor(taxon_val, levels = pc$levels))
  } else {
    pal <- make_taxon_palette(df_bar_tl$taxon_val)
  }
  parent_lookup <- get_parent_lookup(df_bar_tl, col_name)
  label_fn <- make_label_fn(parent_lookup, level = col_name)

  # For species-level plots wrap labels that exceed ~22 chars with <br> at the
  # space nearest to the midpoint so long binomial names don't overflow the legend.
  if (col_name == "Species") {
    raw_label_fn <- label_fn
    label_fn <- function(x) {
      sapply(raw_label_fn(x), function(lbl) {
        if (nchar(lbl) > 22) {
          sp <- gregexpr(" ", lbl)[[1]]
          sp <- sp[sp > 0]
          if (length(sp) == 0) {
            return(lbl)
          }
          mid <- nchar(lbl) %/% 2
          best <- sp[which.min(abs(sp - mid))]
          paste0(substr(lbl, 1, best - 1), "<br>", substr(lbl, best + 1, nchar(lbl)))
        } else {
          lbl
        }
      }, USE.NAMES = FALSE)
    }
  }

  # Uniform bar width across all taxon levels: matches the default 0.9 spacing.
  bar_w <- 0.9

  # ---- type C: stacked by station within depth facets ----
  pC <- ggplot(df_bar_tl, aes(x = Station, fill = taxon_val)) +
    geom_bar(position = "stack", color = "black", linewidth = 0.25, width = bar_w) +
    facet_grid(Depth_label ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = pal, labels = label_fn) +
    scale_y_continuous(breaks = function(x) seq(floor(x[1]), ceiling(x[2]), by = 1)) +
    labs(
      x = "Estación de muestreo",
      y = "Número de MAGs",
      fill = col_label
    ) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.text = element_markdown(),
      strip.text.y = element_text(face = "bold", angle = 0)
    )

  n_depths_c <- length(unique(df_bar_tl$Depth_label))
  # Species canvas is only 1 in wider than others to fit the longer legend;
  # keeping it close to 10 ensures bins (and therefore bars + gaps) stay the
  # same physical size as in all other taxon-level plots.
  w_c <- if (col_name == "Species") 11 else 10
  save_fig(pC, paste0("Fig1C_barras_", prefix, "_estacion_x_profundidad.png"),
    width = w_c, height = 4 + n_depths_c * 2
  )
  if (prefix == "filo") p_panel_a <- pC
  if (prefix == "orden") p_panel_b <- pC

  # ---- type D: the transpose, stacked by depth within station facets ----
  pD <- ggplot(df_bar_tl, aes(x = Depth_label, fill = taxon_val)) +
    geom_bar(position = "stack", color = "black", linewidth = 0.25, width = bar_w) +
    facet_grid(. ~ Station, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = pal, labels = label_fn) +
    labs(
      x = "Profundidad",
      y = "Número de MAGs",
      fill = col_label
    ) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.text = element_markdown()
    )

  n_stations_d <- length(unique(df_bar_tl$Station))
  w_d <- if (col_name == "Species") 4.5 + n_stations_d * 1.9 else 4 + n_stations_d * 1.8
  save_fig(pD, paste0("Fig1D_barras_", prefix, "_profundidad_x_estacion.png"),
    width = w_d, height = 6
  )
  if (prefix == "especie") p_panel_c <- pD
}

# === SECTION 1c: combined phylum / order / species panel ===
# a (top left): Fig1C_barras_filo_estacion_x_profundidad
# b (top right): Fig1C_barras_orden_estacion_x_profundidad
# c (bottom, full width): Fig1D_barras_especie_profundidad_x_estacion
cat("\n=== Sección 1c: Panel combinado filo/orden/especie ===\n")

{
  n_depths_panel <- length(unique(df_bar$Depth_label))
  h_top <- 4 + n_depths_panel * 2 # same formula as individual pC
  h_bot <- 6 # same formula as individual pD

  # Width is driven by c's natural width (species pD), not by a+b summed.
  # This prevents patchwork from over-stretching c and pushing its legend
  # to the far right. a and b are rendered slightly narrower (w_c/2 each),
  # which is fine for 6-station bar charts.
  n_stations_panel <- length(unique(df_bar$Station))
  w_c_natural <- 4.5 + n_stations_panel * 1.9 # species pD formula
  w_panel <- max(w_c_natural, 16)

  p_combined_panel <- (p_panel_a | p_panel_b) / p_panel_c +
    plot_layout(heights = c(h_top, h_bot)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 20, color = "black"))

  save_fig(p_combined_panel, "Fig1_panel_filo_orden_especie.png",
    width = w_panel, height = h_top + h_bot
  )
}

# ---- section 1b: MAGs with no species assignment ----
# These are the candidate new archaeal species. Plotting where they sit by
# station and depth, coloured by genus, shows where the taxonomic novelty
# concentrates and which lineages account for it.
cat("\n=== Sección 1b: MAGs sin especie asignada ===\n")

df_novel_sp <- df_bar %>%
  filter(is.na(Species)) %>%
  mutate(genus_val = ifelse(is.na(Genus), "Sin asignar", Genus))

cat("  MAGs sin especie asignada:", nrow(df_novel_sp), "\n")

if (nrow(df_novel_sp) > 0) {
  pal_novel <- make_taxon_palette(df_novel_sp$genus_val)
  parent_lookup_1b <- get_parent_lookup(df_novel_sp, "Genus")
  label_fn_1b <- make_label_fn(parent_lookup_1b, level = "Genus")

  # ---- type C: facet by depth, x = station ----
  p1e <- ggplot(df_novel_sp, aes(x = Station, fill = genus_val)) +
    geom_bar(position = "stack", color = "black", linewidth = 0.25) +
    facet_grid(Depth_label ~ ., scales = "free_y", space = "free_y") +
    scale_fill_manual(values = pal_novel, labels = label_fn_1b) +
    scale_y_continuous(breaks = function(x) seq(floor(x[1]), ceiling(x[2]), by = 1)) +
    labs(
      x = "Estación de muestreo",
      y = "Número de MAGs",
      fill = "Género"
    ) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      legend.text = element_markdown(),
      strip.text.y = element_text(face = "bold", angle = 0)
    )

  n_depths_e <- length(unique(df_novel_sp$Depth_label))
  save_fig(p1e, "Fig1E_novel_sp_estacion_x_profundidad.png",
    width = 10, height = 4 + n_depths_e * 2
  )

  # ---- type D: facet by station, x = depth ----
  p1f <- ggplot(df_novel_sp, aes(x = Depth_label, fill = genus_val)) +
    geom_bar(position = "stack", color = "black", linewidth = 0.25) +
    facet_grid(. ~ Station, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = pal_novel, labels = label_fn_1b) +
    labs(
      x = "Profundidad",
      y = "Número de MAGs",
      fill = "Género"
    ) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.text = element_markdown()
    )

  n_stations_f <- length(unique(df_novel_sp$Station))
  save_fig(p1f, "Fig1F_novel_sp_profundidad_x_estacion.png",
    width = 4 + n_stations_f * 1.8, height = 6
  )
} else {
  cat("  [WARN] No hay MAGs sin especie asignada – se omiten Fig1E y Fig1F\n")
}

# ---- section 2: phylogenetic overview, phylum and class ----
# The collection falls into three phyla: Thermoplasmatota,
# Thermoproteota and Undinarchaeota. Phylum and class are the minimum context
# for describing what this set of MAGs covers.
cat("\n=== Sección 2: Distribución de clases arqueales ===\n")

df_phylo <- df %>%
  filter(!is.na(Station), !is.na(Depth)) %>%
  mutate(
    Phylum_plot = ifelse(is.na(Phylum), "Sin asignar", Phylum),
    Class_plot  = ifelse(is.na(Class), "Sin asignar", Class),
    Depth_label = factor(Depth_label, levels = depth_factor_levels)
  )

# 2C: class by station within depth facets
classes_in_data <- sort(unique(na.omit(df_phylo$Class_plot)))
n_cls <- length(classes_in_data)
if (n_cls <= 12) {
  pal_cls <- setNames(
    drop_palette_grays(
      RColorBrewer::brewer.pal(max(n_cls, 3), ifelse(n_cls <= 8, "Set2", "Set3"))[seq_len(n_cls)]
    ),
    classes_in_data
  )
} else {
  pal_cls <- setNames(
    drop_palette_grays(colorRampPalette(RColorBrewer::brewer.pal(12, "Paired"))(n_cls)),
    classes_in_data
  )
}

parent_lookup_cls2 <- get_parent_lookup(df_phylo, "Class")
label_fn_cls2 <- make_label_fn(parent_lookup_cls2, level = "Class")

p2c <- ggplot(df_phylo, aes(x = Station, fill = Class_plot)) +
  geom_bar(position = "stack", color = "black", linewidth = 0.25) +
  facet_grid(Depth_label ~ ., scales = "free_y", space = "free_y") +
  scale_fill_manual(values = pal_cls, labels = label_fn_cls2) +
  scale_y_continuous(breaks = function(x) seq(floor(x[1]), ceiling(x[2]), by = 1)) +
  labs(x = "Estación de muestreo", y = "Número de MAGs", fill = "Clase") +
  theme_publication() +
  theme(
    axis.text.x   = element_text(angle = 0, hjust = 0.5),
    legend.text   = element_markdown(),
    strip.text.y  = element_text(angle = 0, face = "bold", size = 11)
  )

p2d <- ggplot(df_phylo, aes(x = Depth_label, fill = Class_plot)) +
  geom_bar(position = "stack", color = "black", linewidth = 0.25) +
  facet_grid(. ~ Station, scales = "free_x", space = "free_x") +
  scale_fill_manual(values = pal_cls, labels = label_fn_cls2, name = "Clase") +
  scale_y_continuous(breaks = function(x) seq(floor(x[1]), ceiling(x[2]), by = 1)) +
  labs(x = "Profundidad", y = "Número de MAGs", fill = "Clase") +
  theme_publication() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.text = element_markdown()
  )

n_depths_2c <- length(unique(df_phylo$Depth_label))
save_fig(p2c, "Fig2C_clases_estacion_profundidad.png",
  width = 10, height = 4 + n_depths_2c * 2.5
)
save_fig(p2d, "Fig2D_clases_profundidad_estacion.png", width = 12, height = 6)

# ---- section 7: bubble chart of ecological niches ----
# With few categories, a bubble chart shows presence, relative abundance and
# spatial distribution at once, where a heatmap would only show the last two.
cat("\n=== Sección 7: Gráficos de burbujas por nivel taxonómico ===\n")

make_bubble_plot <- function(df_in, taxon_col, y_label, title_str) {
  df_b <- df_in %>%
    filter(!is.na(Station), !is.na(Depth)) %>%
    mutate(
      taxon_plot  = ifelse(is.na(.data[[taxon_col]]), "Sin asignar", .data[[taxon_col]]),
      Depth_label = factor(Depth_label, levels = depth_factor_levels)
    ) %>%
    count(taxon_plot, Station, Depth_label, name = "n")

  parent_lookup_b <- get_parent_lookup(df_in, taxon_col)
  label_fn_b <- make_label_fn(parent_lookup_b, level = taxon_col)

  n_taxa <- n_distinct(df_b$taxon_plot)
  h_fig <- max(5, 1.5 + n_taxa * 0.35)

  p <- ggplot(df_b, aes(x = Station, y = taxon_plot, size = n, fill = Depth_label)) +
    geom_point(shape = 21, color = "black", stroke = 0.4, alpha = 1) +
    facet_grid(. ~ Depth_label) +
    scale_size_area(
      max_size = 12,
      breaks = function(x) {
        b <- pretty(x, n = 4)
        b[b %% 1 == 0 & b >= 1]
      }
    ) +
    scale_fill_manual(values = depth_colors) +
    scale_y_discrete(labels = label_fn_b) +
    labs(
      x = "Estación de muestreo",
      y = y_label,
      size = "MAGs",
      fill = "Profundidad"
    ) +
    guides(
      fill = guide_legend(title = "Profundidad", override.aes = list(size = 6)),
      size = guide_legend(title = "MAGs")
    ) +
    theme_publication() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5),
      axis.text.y = element_markdown()
    )

  list(plot = p, height = h_fig)
}

b7 <- make_bubble_plot(
  df, "Family", "Familia",
  "Distribución de familias arqueales por estación y profundidad"
)
b7b <- make_bubble_plot(
  df, "Class", "Clase",
  "Distribución de clases arqueales por estación y profundidad"
)
b7c <- make_bubble_plot(
  df, "Order", "Orden",
  "Distribución de órdenes arqueales por estación y profundidad"
)
b7d <- make_bubble_plot(
  df, "Genus", "Género",
  "Distribución de géneros arqueales por estación y profundidad"
)

save_fig(b7$plot, "Fig7_burbujas_familia_estacion_profundidad.png", width = 13, height = b7$height)
save_fig(b7b$plot, "Fig7b_burbujas_clase_estacion_profundidad.png", width = 13, height = b7b$height)
save_fig(b7c$plot, "Fig7c_burbujas_orden_estacion_profundidad.png", width = 13, height = b7c$height)
save_fig(b7d$plot, "Fig7d_burbujas_genero_estacion_profundidad.png", width = 13, height = b7d$height)

# ---- section 8: taxonomic diversity per sample, genus richness ----
# Genus richness per station and depth is the basic alpha diversity measure
# here. A gradient across depths would bear on how the oxygen minimum zone and
# hydrostatic pressure shape the archaeal community.
cat("\n=== Sección 8: Riqueza de géneros por muestra ===\n")

df_richness <- df %>%
  filter(!is.na(Station), !is.na(Depth)) %>%
  group_by(Station, Depth, Depth_label) %>%
  summarise(
    n_MAGs = n(),
    n_genera = n_distinct(Genus, na.rm = TRUE),
    n_genera_novel = sum(is.na(Genus)),
    .groups = "drop"
  ) %>%
  mutate(Depth_label = factor(Depth_label, levels = depth_factor_levels))

# Palette for Fig8a (different from depth_colors to distinguish from Fig8b)
richness_colors <- c("ZMOR" = "#66C2A5", "1000 m" = "#FC8D62", "Fondo" = "#8DA0CB")

# 8A: genus richness, grouped bars
p8a <- ggplot(df_richness, aes(x = Station, y = n_genera, fill = Depth_label)) +
  geom_col(position = "dodge", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = richness_colors) +
  scale_y_continuous(
    breaks = function(x) seq(floor(x[1]), ceiling(x[2]), by = 1),
    expand = expansion(mult = c(0, 0.08))
  ) +
  labs(
    x = "Estación de muestreo",
    y = "Número de géneros",
    fill = "Profundidad"
  ) +
  theme_publication()

# 8B: MAGs per sample, grouped bars
p8b <- ggplot(df_richness, aes(x = Station, y = n_MAGs, fill = Depth_label)) +
  geom_col(position = "dodge", color = "black", linewidth = 0.3) +
  geom_text(aes(label = n_MAGs),
    position = position_dodge(width = 0.9), vjust = -0.4, size = 5
  ) +
  scale_fill_manual(values = richness_colors) +
  labs(
    x = "Estación de muestreo",
    y = "Número de MAGs",
    fill = "Profundidad"
  ) +
  theme_publication()

p8_combined <- (p8a / p8b) +
  plot_layout(guides = "collect") +
  plot_annotation(tag_levels = "a") &
  theme_publication() &
  theme(plot.tag = element_text(face = "bold", size = 20, color = "black"))

save_fig(p8_combined, "Fig8_riqueza_MAGs_por_muestra.png", width = 10, height = 10)

# ---- section 9: shared and exclusive genera ----
# A genus confined to one station or depth, against one found at all of them,
# bears on dispersal, on how connected the water masses are, and on how
# specific the niche is.
cat("\n=== Sección 9: Exclusividad geográfica de géneros por estación ===\n")

excl_st_colors <- c(
  "Exclusivo de una estación"         = "#E74C3C",
  "Distribución local (2-3 est.)"     = "#F39C12",
  "Ampliamente distribuido (≥4 est.)" = "#27AE60"
)

# Helper: build exclusivity-by-station table for any taxonomic level
make_excl_station_df <- function(df_in, taxon_col, n_total_stations = 6) {
  df_in %>%
    filter(!is.na(.data[[taxon_col]]), !is.na(Station)) %>%
    group_by(.data[[taxon_col]]) %>%
    summarise(n_stations = n_distinct(Station), .groups = "drop") %>%
    rename(Taxon = 1) %>%
    mutate(
      Exclusividad = case_when(
        n_stations == 1 ~ "Exclusivo de una estación",
        n_stations <= 3 ~ "Distribución local (2-3 est.)",
        TRUE ~ "Ampliamente distribuido (≥4 est.)"
      )
    )
}

# Fig9: genus
df_excl_genero <- make_excl_station_df(df, "Genus")

# ---- section 9b: station exclusivity at phylum, class, order and family ----
# The same exclusivity analysis above genus rank tests whether whole lineages
# (phyla, classes) stay at some stations or spread across all of them.
cat("\n=== Sección 9b: Exclusividad por estación (Filo, Clase, Orden, Familia) ===\n")

df_excl_filo <- make_excl_station_df(df, "Phylum")
df_excl_clase <- make_excl_station_df(df, "Class")
df_excl_orden <- make_excl_station_df(df, "Order")
df_excl_familia <- make_excl_station_df(df, "Family")

# Export detailed exclusivity tables to Excel (one sheet per taxonomic level)
cat("\n=== Exportando tabla de exclusividad a Excel ===\n")

if (!requireNamespace("openxlsx", quietly = TRUE)) {
  cat("[WARN] Instalando openxlsx...\n")
  install.packages("openxlsx", quiet = TRUE)
}
suppressMessages(library(openxlsx, quietly = TRUE))

excl_wb <- createWorkbook()

add_excl_sheet <- function(wb, sheet_name, df_excl) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, df_excl %>% arrange(Exclusividad, across(1)))
}

add_excl_sheet(excl_wb, "Genero_estacion", df_excl_genero %>% rename(Genero = Taxon))
add_excl_sheet(excl_wb, "Filo_estacion", df_excl_filo %>% rename(Filo = Taxon))
add_excl_sheet(excl_wb, "Clase_estacion", df_excl_clase %>% rename(Clase = Taxon))
add_excl_sheet(excl_wb, "Orden_estacion", df_excl_orden %>% rename(Orden = Taxon))
add_excl_sheet(excl_wb, "Familia_estacion", df_excl_familia %>% rename(Familia = Taxon))

excl_xlsx_path <- file.path(output_dir, "tax_exclusividad_resumen.xlsx")
saveWorkbook(excl_wb, excl_xlsx_path, overwrite = TRUE)
cat("[OK] Saved: tax_exclusividad_resumen.xlsx\n")

# ---- === Summary table: taxonomy per MAG ----
cat("\n=== Exportando tabla resumen ===\n")

df_summary <- df %>%
  select(
    OriginalName, NameWithoutExtension, Station, Depth,
    Domain, Phylum, Class, Order, Family, Genus, Species,
    classification_method, ANI_closest, AF_closest,
    msa_percent, red_value,
    Novel_Species, Novel_Genus, Novel_Family
  ) %>%
  arrange(Station, Depth, Phylum, Genus)

write_tsv(df_summary, file.path(output_dir, "tax_analysis_summary.tsv"))
cat("[OK] Saved: tax_analysis_summary.tsv\n")

cat("\n========================================\n")
cat("[OK] Análisis taxonómico completado\n")
cat("Resultados guardados en:", output_dir, "\n")
cat("========================================\n\n")

# ---- ---- quality and taxonomy table for the unassigned MAGs ----
cat("=== GENERANDO TABLA DE CALIDAD/TAXONOMÍA (MAGs sin especie) ===\n")

quality_report_file <- file.path(PROJECT_ROOT, "results/01_mag_quality/quality_report.tsv")
genome_stats_file <- file.path(PROJECT_ROOT, "results/10_MAG_distill/SUMMARIZE/genome_stats.tsv")
metabolism_file <- file.path(PROJECT_ROOT, "results/10_MAG_distill/SUMMARIZE/metabolism_summary.xlsx")

mags_no_species <- df %>%
  filter(Novel_Species) %>%
  pull(NameWithoutExtension)

if (length(mags_no_species) == 0) {
  cat("[WARN] No hay MAGs sin especie. Omitiendo tabla.\n\n")
} else {
  cat(sprintf("  MAGs sin especie: %d\n", length(mags_no_species)))

  tryCatch(
    {
      # ---- 1. completeness and contamination, from the quality report ----
      # Dotted names after utils/05_rename_mags.sh, underscores before it;
      # both reduce to NameWithoutExtension.
      qual_df <- read_tsv(quality_report_file, show_col_types = FALSE) %>%
        mutate(NameWithoutExtension = gsub("_", ".", Name, fixed = TRUE)) %>%
        filter(NameWithoutExtension %in% mags_no_species) %>%
        select(
          genome            = NameWithoutExtension,
          Completitud_pct   = Completeness,
          Contaminacion_pct = Contamination
        )

      # --- 2. rRNA + tRNA totales (genome_stats.tsv) ---
      gs_rrna <- read_tsv(genome_stats_file, show_col_types = FALSE) %>%
        filter(genome %in% mags_no_species) %>%
        select(genome,
          r23S       = `23S rRNA`,
          r16S       = `16S rRNA`,
          r5S        = `5S rRNA`,
          tRNA_total = `tRNA count`
        ) %>%
        mutate(
          `23S rRNA`  = ifelse(r23S >= 1, "Sí", "No"),
          `16S rRNA`  = ifelse(r16S >= 1, "Sí", "No"),
          `5S rRNA`   = ifelse(r5S >= 1, "Sí", "No"),
          tRNAs_total = as.integer(tRNA_total)
        ) %>%
        select(genome, `23S rRNA`, `16S rRNA`, `5S rRNA`, tRNAs_total)

      # ---- 3. distinct tRNAs per amino acid, from the tRNA sheet ----
      trna_raw <- read_excel(metabolism_file, sheet = "tRNA")
      mag_cols_t <- intersect(names(trna_raw)[-1], mags_no_species)
      trna_names <- as.character(trna_raw[[1]])
      valid_trna_rows <- grepl("^[A-Za-z]+ \\(", trna_names)
      amino_acids <- sub(" .*", "", trna_names)

      trna_mat <- trna_raw[, mag_cols_t, drop = FALSE] %>%
        mutate(across(everything(), ~ suppressWarnings(as.numeric(as.character(.)))))

      unique_trna <- vapply(mag_cols_t, function(mg) {
        present_idx <- which(valid_trna_rows & !is.na(trna_mat[[mg]]) & trna_mat[[mg]] > 0)
        length(unique(amino_acids[present_idx]))
      }, integer(1))

      trna_uniq_df <- tibble(
        genome       = mag_cols_t,
        tRNAs_unicos = as.integer(unique_trna)
      )

      # ---- 4. taxonomy: ANI, AF, MSA% and RED, from gtdb_raw ----
      tax_sel <- gtdb_raw %>%
        select(
          orig_name     = 1,
          genome_ani    = 6,
          genome_af     = 7,
          placement_ani = 11,
          placement_af  = 12,
          msa_pct       = 17,
          red_value     = 19
        ) %>%
        mutate(
          closest_ani = coalesce(
            suppressWarnings(as.numeric(placement_ani)),
            suppressWarnings(as.numeric(genome_ani))
          ),
          closest_af = coalesce(
            suppressWarnings(as.numeric(placement_af)),
            suppressWarnings(as.numeric(genome_af))
          )
        ) %>%
        mutate(NameWithoutExtension = gsub("_", ".", orig_name, fixed = TRUE)) %>%
        filter(NameWithoutExtension %in% mags_no_species) %>%
        select(genome = NameWithoutExtension, closest_ani, closest_af, msa_pct, red_value)

      # ---- 5. display names ----
      display_map <- setNames(
        meta_raw %>% filter(NameWithoutExtension %in% mags_no_species) %>% pull(DisplayName),
        meta_raw %>% filter(NameWithoutExtension %in% mags_no_species) %>% pull(NameWithoutExtension)
      )

      # --- 6. Assemble the table ---
      tabla_10 <- tibble(genome = mags_no_species) %>%
        mutate(`Nombre MAG` = as.character(display_map[genome])) %>%
        left_join(qual_df, by = "genome") %>%
        left_join(gs_rrna, by = "genome") %>%
        left_join(trna_uniq_df, by = "genome") %>%
        left_join(tax_sel, by = "genome") %>%
        select(
          `Nombre MAG`,
          `Completitud (%)`            = Completitud_pct,
          `Contaminación (%)`          = Contaminacion_pct,
          `23S rRNA`,
          `16S rRNA`,
          `5S rRNA`,
          `tRNAs totales`              = tRNAs_total,
          `tRNAs únicos`               = tRNAs_unicos,
          `ANI más cercano (%)`        = closest_ani,
          `AF alineamiento más cercan` = closest_af,
          `MSA (%)`                    = msa_pct,
          `RED`                        = red_value
        )

      # ---- 7. write the TSV with a UTF-8 BOM so Excel reads the accents ----
      tsv_path <- file.path(output_dir, "tabla_calidad_10MAGs.tsv")
      write_tsv(tabla_10, tsv_path)
      raw_content <- readBin(tsv_path, "raw", n = file.size(tsv_path))
      writeBin(c(as.raw(c(0xEF, 0xBB, 0xBF)), raw_content), tsv_path)
      cat("  [OK] tabla_calidad_10MAGs.tsv guardada\n")

      # --- 8. Save the .png (flextable booktabs) ---
      if (requireNamespace("flextable", quietly = TRUE) &&
        requireNamespace("officer", quietly = TRUE)) {
        library(flextable)
        library(officer)

        ft_q <- flextable(tabla_10) %>%
          colformat_double(
            j = c(
              "Completitud (%)", "Contaminación (%)",
              "ANI más cercano (%)", "AF alineamiento más cercan",
              "MSA (%)", "RED"
            ),
            digits = 2
          ) %>%
          bold(part = "header") %>%
          align(align = "center", part = "header") %>%
          align(j = "Nombre MAG", align = "left", part = "body") %>%
          align(j = -1L, align = "center", part = "body") %>%
          border_remove() %>%
          hline_top(border = fp_border(color = "black", width = 2), part = "header") %>%
          hline_bottom(border = fp_border(color = "black", width = 1.5), part = "header") %>%
          hline_bottom(border = fp_border(color = "black", width = 2), part = "body") %>%
          fontsize(size = 9, part = "all") %>%
          font(fontname = "Times New Roman", part = "all") %>%
          theme_booktabs() %>%
          autofit()

        save_as_image(ft_q,
          path = file.path(output_dir, "tabla_calidad_10MAGs.png"),
          zoom = 4, expand = 10, webshot = "webshot2"
        )
        cat("  [OK] tabla_calidad_10MAGs.png guardada\n")
      } else {
        cat("  [WARN] flextable/officer no disponibles; omitiendo PNG de tabla\n")
        cat("    Instálalos con: install.packages(c('flextable', 'officer'))\n")
      }
      cat("[OK] Tabla de calidad/taxonomía completada\n\n")
    },
    error = function(e) {
      cat("  [WARN] Error al generar tabla de calidad:", conditionMessage(e), "\n\n")
    }
  )
}

# ---- section 10: skani, novelty as a number ----
#
# GTDB-Tk places the unassigned MAGs using phylogenetic distance thresholds,
# the RED values, computed on marker trees. skani computes ANI against the
# GTDB R226 representatives instead, a continuous measure of genomic distance
# that is easier to read:
#
#   ANI >= 95% with AF >= 50%   same species (Jain et al. 2018)
#   ANI 90 to 95% with AF >= 50%  new species inside a known genus
#   AF < 50%                    no close genomic relative in GTDB R226,
#                               possibly novel at genus level or above
#
# The aligned fraction matters most. A high ANI over a low AF
# means the alignment found conserved housekeeping regions and little else,
# which is not evidence of genome-wide homology.
#
# ANI among the collection's own MAGs, from the all-vs-all matrix, tests
# whether the novel MAGs are strains of one species recovered from several
# stations and depths or independent lineages, which bears on dispersal and on
# how connected the water masses of the Gulf of Mexico are.
#
# Figures:
#   10A  lollipop, ANI to the closest GTDB representative, sized by AF_query
#   10B  scatter, every GTDB hit above 5% AF, y axis AF_query
#   11   heatmap, ANI between the 10 novel MAGs and the whole collection
cat("\n=== Sección 10: Análisis skani – novedad taxonómica cuantitativa ===\n")

skani_gtdb_file <- file.path(PROJECT_ROOT, "results/04_sketching_ani/gtdb_search/novel10_vs_gtdb.tsv")
skani_sparse_file <- file.path(PROJECT_ROOT, "results/04_sketching_ani/all_mags/ani_triangle_sparse.tsv")

if (!file.exists(skani_gtdb_file) || !file.exists(skani_sparse_file)) {
  cat("[WARN] Archivos skani no encontrados – omitiendo Sección 10\n")
} else {
  # --- Helper: strip path and extension from genome file names ---
  basename_no_ext <- function(path) {
    b <- basename(path)
    b <- sub("\\.fa\\.gz$", "", b)
    b <- sub("\\.fna\\.gz$", "", b)
    b <- sub("\\.fa$", "", b)
    b <- sub("\\.fna$", "", b)
    return(b)
  }

  # --- Helper: simplify GTDB Ref_name to a readable taxon label ---
  clean_gtdb_name <- function(nm) {
    # 1. Remove leading accession (e.g. "DAVI01000044.1 ")
    x <- sub("^[A-Z0-9_.]+\\s+", "", nm)
    # 2. Strip common prefixes
    x <- sub("^MAG TPA_asm(_mismatch)?:\\s*", "", x)
    x <- sub("^MAG:\\s*", "", x)
    x <- sub("^TPA_asm:\\s*", "", x)
    # 3. Remove trailing clauses
    x <- sub("\\s+isolate\\s+.*$", "", x)
    x <- sub("\\s+genome assembly.*$", "", x)
    x <- sub("\\s+whole genome.*$", "", x)
    x <- sub(",\\s.*$", "", x)
    trimws(x)
  }

  # --- Display name map (NameWithoutExtension to DisplayName) ---
  dmap <- setNames(meta_raw$DisplayName, meta_raw$NameWithoutExtension)
  # Italicize genus name in display labels (ggtext markdown)
  dmap <- gsub("Nitrosopelagicus", "*Nitrosopelagicus*", dmap)
  get_display <- function(ids) {
    d <- dmap[ids]
    ifelse(is.na(d), ids, d) # fall back to raw ID if not found
  }

  # --- Novel MAG IDs (those without species assignment) ---
  novel_ids <- df %>%
    filter(Novel_Species) %>%
    pull(NameWithoutExtension)
  cat(sprintf("  MAGs sin especie asignada: %d\n", length(novel_ids)))

  # ---- Load and parse novel10_vs_gtdb.tsv ----
  skani_gtdb_raw <- read_tsv(skani_gtdb_file, show_col_types = FALSE)

  skani_gtdb <- skani_gtdb_raw %>%
    rename(
      Ref_file   = 1, Query_file = 2,
      ANI        = 3, AF_ref     = 4, AF_query = 5,
      Ref_name   = 6, Query_name = 7
    ) %>%
    mutate(
      MAG_id         = basename_no_ext(Query_file),
      ref_name_clean = clean_gtdb_name(Ref_name)
    ) %>%
    filter(MAG_id %in% novel_ids)

  # Best hit per MAG: highest AF_query × ANI (penalises low-coverage hits)
  skani_best <- skani_gtdb %>%
    mutate(score = AF_query * ANI) %>%
    group_by(MAG_id) %>%
    slice_max(order_by = score, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    left_join(
      df %>% select(NameWithoutExtension, Phylum, Class, Order, Genus),
      by = c("MAG_id" = "NameWithoutExtension")
    ) %>%
    mutate(
      DisplayName = get_display(MAG_id),
      # Novelty level: AF < 50% to no full-genome match in GTDB
      Novelty = case_when(
        ANI >= 95 & AF_query >= 50 ~ "Especie conocida (ANI ≥ 95%)",
        ANI >= 90 & AF_query >= 50 ~ "Nueva especie (ANI 90–95%)",
        ANI >= 85 & AF_query >= 50 ~ "Nueva especie (ANI 85–90%)",
        TRUE ~ "Sin pariente genómico\ncompleto en GTDB"
      ),
      # Truncate label for figure
      ref_label = ifelse(nchar(ref_name_clean) > 42,
        paste0(substr(ref_name_clean, 1, 39), "\u2026"),
        ref_name_clean
      )
    ) %>%
    # Order rows by ANI ascending (lowest novelty to top of plot)
    arrange(ANI) %>%
    mutate(DisplayName = factor(DisplayName, levels = DisplayName))

  cat("  Best GTDB hits calculados\n")

  # Fig 10A: lollipop, ANI to the closest GTDB R226 representative.
  # Puts every novel MAG on the same axis against its nearest known relative.
  # Point size encodes the aligned fraction, which resolves the
  # ambiguity of an ANI computed over short conserved stretches.
  novelty_colors_10 <- c(
    "Nueva especie (ANI 90–95%)" = "#2980B9",
    "Nueva especie (ANI 85–90%)" = "#8E44AD",
    "Sin pariente genómico\ncompleto en GTDB" = "#E74C3C"
  )

  # Pre-compute label color: white text for large circles (AF ≥ 50),
  # black text for small circles so the % is always legible.
  skani_best <- skani_best %>%
    mutate(
      af_label     = paste0(round(AF_query), "%"),
      af_txt_color = ifelse(AF_query >= 50, "white", "black")
    )

  p10a <- ggplot(skani_best, aes(x = ANI, y = DisplayName)) +
    annotate("rect",
      xmin = 95, xmax = Inf, ymin = -Inf, ymax = Inf,
      fill = "#27AE60", alpha = 0.06
    ) +
    geom_point(
      aes(size = AF_query, fill = Novelty),
      shape = 21, color = "black", stroke = 0.5
    ) +
    # AF% label to the left of each bubble (x offset in data units, hjust=1 = right-align)
    geom_text(
      aes(label = af_label, x = ANI - 1.1),
      hjust = 1, size = 4, fontface = "bold", color = "gray20"
    ) +
    geom_text(
      aes(label = ref_label),
      hjust = -0.1, size = 3.5, color = "gray35", fontface = "italic"
    ) +
    geom_vline(
      xintercept = 95, linetype = "dashed",
      color = "#27AE60", linewidth = 0.8
    ) +
    annotate("text",
      x = 95.2, y = 0.5,
      label = "Umbral especie (95% ANI)",
      hjust = 0, size = 4, color = "#27AE60", fontface = "bold"
    ) +
    scale_size_continuous(
      name   = "Fracción alineada\ndel MAG *query* (%)",
      range  = c(4, 14),
      breaks = c(10, 30, 50, 70, 90),
      labels = paste0(c(10, 30, 50, 70, 90), "%")
    ) +
    scale_fill_manual(
      values = novelty_colors_10,
      name   = "Nivel de novedad",
      drop   = FALSE
    ) +
    scale_x_continuous(
      limits = c(76, 106),
      breaks = seq(80, 100, by = 5),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(
      x = "ANI (%)",
      y = NULL
    ) +
    guides(
      fill = guide_legend(override.aes = list(size = 6)),
      size = guide_legend()
    ) +
    theme_publication() +
    theme(
      legend.title     = element_markdown(),
      legend.position  = "right"
    )

  # ---- Fig 10A, two-panel form: lollipop for ANI beside bars for AF ----

  # ---- two-panel form ----
  # The ANI percentage is printed left of the lollipop, so no lookup is needed.
  p_c_ani <- ggplot(skani_best, aes(x = ANI, y = DisplayName)) +
    annotate("rect",
      xmin = 95, xmax = Inf, ymin = -Inf, ymax = Inf,
      fill = "#27AE60", alpha = 0.06
    ) +
    geom_point(aes(fill = Novelty),
      shape = 21, size = 5,
      color = "black", stroke = 0.5
    ) +
    # ANI% to the left of each bubble (x offset in data units, hjust=1 = right-align)
    geom_text(
      aes(label = paste0(round(ANI, 1), "%"), x = ANI - 0.9),
      hjust = 1, size = 4, fontface = "bold", color = "gray20"
    ) +
    geom_text(aes(label = ref_label),
      hjust = -0.1,
      size = 3.5, color = "gray35", fontface = "italic"
    ) +
    geom_vline(
      xintercept = 95, linetype = "dashed",
      color = "#27AE60", linewidth = 0.8
    ) +
    scale_fill_manual(values = novelty_colors_10, name = "Nivel de novedad") +
    scale_x_continuous(
      limits = c(76, 109), breaks = seq(80, 100, by = 5),
      expand = expansion(mult = c(0, 0))
    ) +
    labs(x = "***Average Nucleotide Identity* (ANI, %)**", y = NULL) +
    guides(fill = guide_legend(override.aes = list(size = 6))) +
    theme_publication() +
    theme(
      axis.text.y   = element_markdown(),
      axis.title.x  = element_markdown(),
      legend.title  = element_markdown()
    )

  p_c_af <- ggplot(skani_best, aes(x = AF_query, y = DisplayName, fill = Novelty)) +
    geom_col(color = "black", linewidth = 0.3, width = 0.6) +
    geom_text(aes(label = paste0(round(AF_query), "%")),
      hjust = -0.15, size = 4.5, fontface = "bold"
    ) +
    geom_vline(
      xintercept = 50, linetype = "dashed",
      color = "gray50", linewidth = 0.6
    ) +
    annotate("text",
      x = 52, y = 0.5, label = "50%",
      hjust = 0, size = 4, color = "gray45"
    ) +
    scale_fill_manual(values = novelty_colors_10, guide = "none") +
    scale_x_continuous(
      limits = c(0, 100),
      breaks = seq(0, 100, by = 25),
      expand = expansion(mult = c(0, 0.05)),
      labels = function(x) as.character(x)
    ) +
    labs(x = "Fracción alineada con *query* (%)", y = NULL) +
    theme_publication() +
    theme(
      axis.text.y   = element_blank(),
      axis.ticks.y  = element_blank(),
      axis.title.x  = element_markdown()
    )

  p10a_opC <- (p_c_ani | p_c_af) +
    plot_layout(widths = c(2.2, 1), guides = "collect") +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 20, color = "black"))
  save_fig(p10a_opC, "Fig10A_opC_dual_panel.png", width = 18, height = 6)

  # Fig 10B: scatter of ANI against AF_query, every GTDB hit above 5% AF.
  # The cloud is each novel MAG's signature in GTDB reference space. High AF
  # with high ANI marks a close relative; very low AF with high ANI
  # indicates alignment on conserved genes only (false positives). The
  # top left quadrant (AF >= 50%, ANI < 95%) is confirmed new species; the
  # bottom band (AF < 50%) is a taxon GTDB R226 does not represent well.
  skani_hits_filt <- skani_gtdb %>%
    filter(AF_query > 5) %>%
    mutate(DisplayName = get_display(MAG_id))

  # Consistent color per MAG (10 colors from Paired palette)
  novel_display_names <- get_display(novel_ids)
  pal_novel_mags <- setNames(
    colorRampPalette(RColorBrewer::brewer.pal(10, "Paired"))(length(novel_ids)),
    novel_display_names
  )

  # Quadrant shading:
  #   AF >= 50% and ANI < 95%: new species in a known genus. The alignment
  #     covers over half the genome, but identity falls short.
  #   AF < 50%: fragmentary alignment on conserved genes,
  #     only shared housekeeping regions, which is not evidence of kinship at
  #     species or genus level.
  # Not shaded:
  #   AF >= 50% and ANI >= 95%: GTDB-Tk left these MAGs without a species, so
  #     a skani hit in this corner is most likely a partial-coverage artefact.
  #     It is not labelled here.
  #   AF < 50% and ANI >= 95%: high-identity alignment with low coverage, which
  #     may reflect shared genomic islands. That reading is ambiguous,
  #     and assigning it a colour would claim more than the data supports.
  p10b <- ggplot(
    skani_hits_filt,
    aes(x = ANI, y = AF_query, fill = DisplayName)
  ) +
    # Zone 1: confirmed new species (AF >= 50%, ANI < 95%)
    annotate("rect",
      xmin = 79, xmax = 95, ymin = 50, ymax = 100,
      fill = "#AED6F1", alpha = 0.35
    ) +
    # Zone 2: fragmentary alignment, AF < 50%, not diagnostic
    annotate("rect",
      xmin = 79, xmax = 101, ymin = 0, ymax = 50,
      fill = "#F9E79F", alpha = 0.35
    ) +
    geom_jitter(
      shape = 21, size = 2.4, color = "black", stroke = 0.3,
      alpha = 0.85, width = 0.05, height = 0.3
    ) +
    geom_hline(
      yintercept = 50, linetype = "dashed",
      color = "gray45", linewidth = 0.6
    ) +
    geom_vline(
      xintercept = 95, linetype = "dashed",
      color = "#1A5276", linewidth = 0.7
    ) +
    annotate("text",
      x = 94.8, y = 53, label = "95%",
      hjust = 1, size = 4, color = "#1A5276", fontface = "bold"
    ) +
    annotate("text",
      x = 99.5, y = 51.5, label = "AF = 50%",
      hjust = 1, size = 4, color = "gray40"
    ) +
    scale_fill_manual(values = pal_novel_mags, name = "MAG") +
    scale_x_continuous(
      limits = c(79, 101),
      breaks = seq(80, 100, by = 5)
    ) +
    # Zone labels outside the plot area
    annotate("text",
      x = 79.3, y = 103,
      label = "Nueva especie (ANI < 95%, AF \u2265 50%)",
      hjust = 0, size = 3.5, color = "#1A5276", fontface = "bold"
    ) +
    annotate("text",
      x = 79.3, y = 3,
      label = "Alineamiento fragmentario (AF < 50%) \u2014 genes conservados compartidos, no parentesco gen\u00f3mico",
      hjust = 0, size = 3.5, color = "#7D6608"
    ) +
    scale_y_continuous(
      breaks = seq(0, 100, by = 25),
      labels = seq(0, 100, by = 25)
    ) +
    coord_cartesian(clip = "off", ylim = c(0, 100)) +
    labs(
      x = "***Average Nucleotide Identity* (ANI)**",
      y = "Fracción alineada \u2013 *query* (%)"
    ) +
    theme_publication() +
    theme(
      plot.margin     = margin(t = 8, r = 5, b = 8, l = 5, unit = "mm"),
      axis.title.x    = element_markdown(),
      axis.title.y    = element_markdown(),
      legend.text     = element_markdown(),
      legend.position = "right"
    )

  save_fig(p10b, "Fig10B_skani_ani_vs_af_scatter.png", width = 12, height = 7)

  # Fig 11: ANI heatmap, novel MAGs against the whole collection.
  # The all-vs-all comparison inside the collection shows whether the novel
  # MAGs are geographically spread populations of one species (same
  # species at different stations or depths, ANI >= 95%, AF >= 50%)
  # or separate lineages. The biogeography argument rests on that
  # distinction: hydrological connectivity between stations against niche
  # specificity by depth.
  skani_sparse_raw <- read_tsv(skani_sparse_file, show_col_types = FALSE)

  skani_sparse <- skani_sparse_raw %>%
    rename(
      Ref_file   = 1, Query_file  = 2,
      ANI        = 3, AF_ref      = 4, AF_query = 5,
      Ref_name   = 6, Query_name  = 7
    ) %>%
    mutate(
      Ref_id   = basename_no_ext(Ref_file),
      Query_id = basename_no_ext(Query_file)
    )

  # Extract hits where a novel MAG is either ref or query.
  # Normalise so that novel MAG = "focal" and the other = "partner".
  # AF_focal = proportion of the focal (novel) MAG that aligned.
  hits_as_ref <- skani_sparse %>%
    filter(Ref_id %in% novel_ids) %>%
    transmute(
      focal = Ref_id,
      partner = Query_id,
      ANI,
      AF_focal = AF_ref,
      AF_partner = AF_query
    )

  hits_as_query <- skani_sparse %>%
    filter(Query_id %in% novel_ids, !Ref_id %in% novel_ids) %>%
    transmute(
      focal = Query_id,
      partner = Ref_id,
      ANI,
      AF_focal = AF_query,
      AF_partner = AF_ref
    )

  hits_novel_all <- bind_rows(hits_as_ref, hits_as_query) %>%
    group_by(focal, partner) %>%
    slice_max(order_by = ANI, n = 1, with_ties = FALSE) %>%
    ungroup()

  # Determine column set: partners with AF_focal ≥ 40% + self (diagonal).
  # Only keep IDs present in meta_raw to avoid orphan columns.
  meaningful_partners <- hits_novel_all %>%
    filter(AF_focal >= 40) %>%
    pull(partner) %>%
    unique()
  col_ids <- intersect(
    union(meaningful_partners, novel_ids),
    meta_raw$NameWithoutExtension
  )

  # Order columns biologically: Station to Depth.
  col_order <- meta_raw %>%
    filter(NameWithoutExtension %in% col_ids) %>%
    mutate(
      Depth_ord = match(as.character(Depth), depth_levels),
      Depth_ord = ifelse(is.na(Depth_ord), 99L, Depth_ord)
    ) %>%
    arrange(Station, Depth_ord) %>%
    pull(NameWithoutExtension)
  # Append any IDs not found in meta_raw (safety)
  col_order <- c(col_order, setdiff(col_ids, col_order))

  # Restrict rows to novel MAGs that appear in the sparse file.
  # MAGs in novel_ids but absent from the sparse file would produce all-NA rows.
  novel_ids_plot <- novel_ids[novel_ids %in% unique(hits_novel_all$focal)]
  cat(sprintf("  MAGs novel con datos en matriz skani: %d\n", length(novel_ids_plot)))

  # Build matrix explicitly: rows = novel_ids_plot, cols = col_order.
  # All values initialised to NA (no hit = not reported by skani).
  ani_mat <- matrix(
    NA_real_,
    nrow = length(novel_ids_plot),
    ncol = length(col_order),
    dimnames = list(novel_ids_plot, col_order)
  )

  # Fill self-hits (diagonal = 100 %)
  for (id in novel_ids_plot) {
    if (id %in% col_order) ani_mat[id, id] <- 100
  }

  # Fill cross-hits from sparse matrix (AF_focal ≥ 40 %)
  hits_fill <- hits_novel_all %>%
    filter(AF_focal >= 40, focal %in% novel_ids_plot, partner %in% col_order)
  for (i in seq_len(nrow(hits_fill))) {
    ani_mat[hits_fill$focal[i], hits_fill$partner[i]] <- hits_fill$ANI[i]
  }

  # Matrix already uses NA for "no skani hit"; alias for clarity downstream
  ani_mat_na <- ani_mat

  # Build annotations BEFORE renaming axes, using original IDs as keys
  # so the annotation vectors are guaranteed to be in matrix column/row order.
  col_order_ids <- colnames(ani_mat_na) # NameWithoutExtension order
  row_order_ids <- rownames(ani_mat_na) # NameWithoutExtension order

  col_meta_11 <- tibble(NameWithoutExtension = col_order_ids) %>%
    left_join(
      meta_raw %>%
        select(NameWithoutExtension, Station, Depth) %>%
        mutate(Depth_label_ann = recode(as.character(Depth), !!!depth_labels)),
      by = "NameWithoutExtension"
    )

  col_ann_11 <- HeatmapAnnotation(
    "Estación" = col_meta_11$Station,
    "Profundidad" = col_meta_11$Depth_label_ann,
    col = list(
      "Estación"    = station_colors,
      "Profundidad" = depth_colors
    ),
    annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
    simple_anno_size = unit(4, "mm"),
    border = TRUE,
    show_legend = TRUE,
    show_annotation_name = TRUE
  )

  # Row annotation (novelty level), also built on the original IDs
  row_novelty_11 <- tibble(MAG_id = row_order_ids) %>%
    left_join(skani_best %>% select(MAG_id, Novelty), by = "MAG_id") %>%
    mutate(Novelty = ifelse(
      is.na(Novelty),
      "Sin pariente genómico\ncompleto en GTDB",
      Novelty
    ))

  # Build italic labels for Nitrosopelagicus MAGs (ComplexHeatmap uses expressions).
  # Returns a single expression() vector so ComplexHeatmap renders all labels
  # through plotmath, with no mixed list that would turn expressions into strings.
  # Plain labels use bquote(.(txt)) which plotmath renders as the literal string.
  make_hm_labels <- function(nms) {
    exprs <- lapply(nms, function(nm) {
      dn_clean <- gsub("\\*", "", get_display(nm))
      if (grepl("Nitrosopelagicus", dn_clean)) {
        rest <- sub("Nitrosopelagicus ", "", dn_clean)
        bquote(paste(italic("Nitrosopelagicus"), " ", .(rest)))
      } else {
        bquote(.(dn_clean))
      }
    })
    do.call(expression, exprs)
  }

  row_labels_11 <- make_hm_labels(rownames(ani_mat_na))
  col_labels_11 <- make_hm_labels(colnames(ani_mat_na))

  # Now rename the matrix axes (plain text without markdown, for internal use)
  rownames(ani_mat_na) <- gsub("\\*", "", get_display(rownames(ani_mat_na)))
  colnames(ani_mat_na) <- gsub("\\*", "", get_display(colnames(ani_mat_na)))

  novelty_row_cols <- c(
    "Nueva especie (ANI 90–95%)" = "#2980B9",
    "Nueva especie (ANI 85–90%)" = "#8E44AD",
    "Sin pariente genómico\ncompleto en GTDB" = "#E74C3C"
  )

  row_ann_11 <- rowAnnotation(
    "Novedad" = row_novelty_11$Novelty,
    col = list("Novedad" = novelty_row_cols),
    annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
    simple_anno_size = unit(4, "mm"),
    border = TRUE,
    show_legend = TRUE,
    show_annotation_name = TRUE
  )

  # Color scale: white (NA) to sequential blue gradient
  col_ani_11 <- colorRamp2(
    c(85, 90, 93, 95, 97.5, 100),
    c("#EFF3FF", "#BDD7E7", "#6BAED6", "#3182BD", "#08519C", "#08306B")
  )

  n_cols_11 <- ncol(ani_mat_na)
  n_rows_11 <- nrow(ani_mat_na)
  # Extra left width to accommodate long row names (display names can be ~30 chars)
  w_11 <- max(12, 4.5 + n_cols_11 * 0.55)
  h_11 <- max(5, 2.0 + n_rows_11 * 0.50)

  ht11 <- Heatmap(
    ani_mat_na,
    name = "ANI (%)",
    col = col_ani_11,
    na_col = "white",
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_dend = FALSE,
    show_column_dend = FALSE,
    row_labels = row_labels_11,
    column_labels = col_labels_11,
    row_names_gp = gpar(fontsize = 11),
    row_names_max_width = unit(7, "cm"),
    column_names_gp = gpar(fontsize = 11),
    column_names_rot = 50,
    rect_gp = gpar(col = "black", lwd = 0.4),
    border = TRUE,
    top_annotation = col_ann_11,
    right_annotation = row_ann_11,
    cell_fun = function(j, i, x, y, width, height, fill) {
      val <- ani_mat_na[i, j]
      if (!is.na(val)) {
        lbl <- if (val == 100) "100" else sprintf("%.1f", val)
        grid.text(lbl, x, y,
          gp = gpar(
            fontsize = 8,
            col = if (val >= 96) "white" else "black"
          )
        )
      }
    },
    heatmap_legend_param = list(
      title = "ANI (%)",
      at = c(85, 90, 95, 100),
      labels = c("85%", "90%", "95%", "100%"),
      title_gp = gpar(fontsize = 13, fontface = "bold"),
      labels_gp = gpar(fontsize = 11),
      border = "black",
      legend_direction = "vertical"
    )
  )

  png(file.path(output_dir, "Fig11_skani_ani_novel_vs_coleccion.png"),
    width = w_11, height = h_11, units = "in", res = 300
  )
  draw(ht11,
    merge_legend = TRUE,
    padding = unit(c(4, 25, 4, 4), "mm")
  ) # bottom, left, top, right
  dev.off()
  cat("[OK] Saved: Fig11_skani_ani_novel_vs_coleccion.png\n")

  cat("[OK] Sección 10 completada\n\n")
}

if (FALSE) { # section 12 disabled

  cat("\n=== Sección 12: Parámetros filogenómicos GTDB-Tk ===\n")

  df12 <- df %>%
    filter(!is.na(Station), !is.na(Depth)) %>%
    mutate(
      red_value_num = suppressWarnings(as.numeric(
        ifelse(red_value %in% c("N/A", "", NA_character_), NA_real_, red_value)
      )),
      msa_percent_num = suppressWarnings(as.numeric(
        ifelse(msa_percent %in% c("N/A", "", NA_character_), NA_real_, msa_percent)
      )),
      classif_method = case_when(
        grepl("ani_screen", classification_method, ignore.case = TRUE) ~ "ANI screen",
        grepl("topology", classification_method, ignore.case = TRUE) ~ "Topología filogenética",
        TRUE ~ "Otro"
      ),
      Phylum_plot = ifelse(is.na(Phylum), "Sin asignar", Phylum),
      Class_plot = ifelse(is.na(Class), "Sin asignar", Class),
      Novelty_label = ifelse(Novel_Species, "Sin especie asignada", "Especie asignada")
    )

  method_colors <- c(
    "ANI screen"               = "#27AE60",
    "Topología filogenética"   = "#E74C3C",
    "Otro"                     = "#95A5A6"
  )

  novelty_colors_12 <- c(
    "Sin especie asignada" = "#E74C3C",
    "Especie asignada"     = "#2ECC71"
  )

  # ---- Fig12A: GTDB-Tk classification method by phylum ----
  # MAGs with a close GTDB reference are classified by ANI screen
  # (fast, high confidence). MAGs without a close reference require phylogenetic
  # tree placement (topology). The ratio per phylum reveals which lineages are
  # well-characterized vs. novel in archaeal genomics databases.
  parent_lookup_phy12 <- get_parent_lookup(df12, "Phylum")
  p12a <- ggplot(df12, aes(x = Phylum_plot, fill = classif_method)) +
    geom_bar(position = "stack", color = "black", linewidth = 0.25) +
    scale_fill_manual(values = method_colors, name = "Método de\nclasificación") +
    scale_x_discrete(labels = make_label_fn(parent_lookup_phy12, level = "Phylum")) +
    scale_y_continuous(
      breaks = function(x) seq(0, ceiling(x[2]), by = 2),
      expand = expansion(mult = c(0, 0.08))
    ) +
    labs(x = "Filo", y = "Número de MAGs") +
    theme_publication() +
    theme(axis.text.x = element_markdown(angle = 30, hjust = 1, size = 10))

  save_fig(p12a, "Fig12A_metodo_clasificacion_por_filo.png", width = 8, height = 6)
  cat("[OK] Saved: Fig12A_metodo_clasificacion_por_filo.png\n")

  # ---- Fig12A_clase: GTDB-Tk classification method by class ----
  # Same as Fig12A but at the class level, which resolves within-phylum variation
  # in classification method. Classes with predominantly topology-based placement
  # represent lineages underrepresented in GTDB at species level.
  parent_lookup_cls12 <- get_parent_lookup(df12, "Class")
  p12a_clase <- ggplot(df12, aes(x = Class_plot, fill = classif_method)) +
    geom_bar(position = "stack", color = "black", linewidth = 0.25) +
    scale_fill_manual(values = method_colors, name = "Método de\nclasificación") +
    scale_x_discrete(labels = make_label_fn(parent_lookup_cls12, level = "Class")) +
    scale_y_continuous(
      breaks = function(x) seq(0, ceiling(x[2]), by = 2),
      expand = expansion(mult = c(0, 0.08))
    ) +
    labs(x = "Clase", y = "Número de MAGs") +
    theme_publication() +
    theme(axis.text.x = element_markdown(angle = 30, hjust = 1, size = 10))

  save_fig(p12a_clase, "Fig12A_clase_metodo_clasificacion_por_clase.png", width = 9, height = 6)
  cat("[OK] Saved: Fig12A_clase_metodo_clasificacion_por_clase.png\n")

  # ---- Fig12B: RED value by class, coloured by novelty ----
  # RED (Relative Evolutionary Divergence) measures a MAG's position
  # on the GTDB reference tree relative to rank-level boundaries (Parks et al. 2022).
  # Values >0.95 indicate placement near the species level; 0.85-0.95 indicates
  # within genus but outside species boundaries (consistent with novel species).
  # Only MAGs classified by tree topology have RED values; ANI-screen MAGs are
  # excluded because they are placed by identity, not tree position.
  df12_red <- df12 %>%
    filter(!is.na(red_value_num)) %>%
    mutate(Class_plot = factor(Class_plot, levels = sort(unique(Class_plot))))

  if (nrow(df12_red) >= 2) {
    p12b <- ggplot(
      df12_red,
      aes(x = Class_plot, y = red_value_num, color = Novelty_label)
    ) +
      geom_jitter(
        width = 0.18, height = 0, size = 4, alpha = 0.85, shape = 16,
        stroke = 0
      ) +
      geom_hline(
        yintercept = 0.95, linetype = "dashed",
        color = "gray40", linewidth = 0.7
      ) +
      annotate("text",
        x = 0.55, y = 0.955,
        label = "RED \u2248 0.95 (límite especie–género)", hjust = 0,
        size = 3, color = "gray35", fontface = "italic"
      ) +
      scale_color_manual(values = novelty_colors_12, guide = "none") +
      scale_x_discrete(labels = make_label_fn(parent_lookup_cls12, level = "Class")) +
      scale_y_continuous(
        limits = c(NA, 1.01),
        breaks = seq(0.88, 1.0, by = 0.02),
        expand = expansion(mult = c(0.05, 0.03))
      ) +
      labs(
        x = "Clase",
        y = "RED value (divergencia evolutiva relativa)"
      ) +
      theme_publication() +
      theme(axis.text.x = element_markdown(angle = 30, hjust = 1, size = 10))

    save_fig(p12b, "Fig12B_red_value_por_clase.png", width = 9, height = 6)
    cat("[OK] Saved: Fig12B_red_value_por_clase.png\n")
  } else {
    cat("  [WARN] Insuficientes datos de RED value – omitiendo Fig12B\n")
  }

  # ---- Fig12C: MSA%, ar53 marker coverage, by class ----
  # The GTDB ar53 marker set comprises 53 conserved archaeal genes.
  # The MSA% is the proportion of these markers identified in each MAG.
  # Low values (<30%) can indicate draft-quality assemblies or lineages whose
  # genes diverge from the reference markers, which bears on placement
  # confidence. High MSA% in topology-classified MAGs (those lacking a close
  # ANI reference) confirms good genome quality despite taxonomic novelty.
  df12_msa <- df12 %>%
    filter(!is.na(msa_percent_num)) %>%
    mutate(Class_plot = factor(Class_plot, levels = sort(unique(Class_plot))))

  if (nrow(df12_msa) >= 2) {
    p12c <- ggplot(
      df12_msa,
      aes(x = Class_plot, y = msa_percent_num, color = Novelty_label)
    ) +
      geom_jitter(
        width = 0.18, height = 0, size = 4, alpha = 0.85, shape = 16,
        stroke = 0
      ) +
      geom_hline(
        yintercept = 30, linetype = "dashed",
        color = "#E67E22", linewidth = 0.7
      ) +
      annotate("text",
        x = 0.55, y = 31.5,
        label = "Umbral m\u00ednimo recomendado (30%)", hjust = 0,
        size = 3, color = "#E67E22", fontface = "italic"
      ) +
      scale_color_manual(values = novelty_colors_12, guide = "none") +
      scale_x_discrete(labels = make_label_fn(parent_lookup_cls12, level = "Class")) +
      scale_y_continuous(
        breaks = seq(0, 100, by = 10),
        expand = expansion(mult = c(0.02, 0.05))
      ) +
      labs(
        x = "Clase",
        y = "Cobertura de marcadores ar53 (MSA, %)"
      ) +
      theme_publication() +
      theme(axis.text.x = element_markdown(angle = 30, hjust = 1, size = 10))

    save_fig(p12c, "Fig12C_msa_percent_por_clase.png", width = 9, height = 6)
    cat("[OK] Saved: Fig12C_msa_percent_por_clase.png\n")
  } else {
    cat("  [WARN] Insuficientes datos de MSA% – omitiendo Fig12C\n")
  }

  # ---- Fig12D: RED value against MSA%, phylogenetic position against quality ----
  # The joint distribution of RED value and MSA% characterises the
  # quality–placement trade-off. MAGs with high MSA% and high RED (top-right)
  # are high-quality genomes placed near the species level, the most reliable
  # classifications. MAGs with low MSA% or extreme RED may need additional
  # validation. Novel MAGs (without species) plotted here show whether their
  # classification uncertainty stems from quality (low MSA%) or novelty
  # (high MSA%, moderate RED, no ANI reference).
  df12_joint <- df12 %>% filter(!is.na(red_value_num), !is.na(msa_percent_num))

  if (nrow(df12_joint) >= 3) {
    p12d <- ggplot(
      df12_joint,
      aes(
        x = msa_percent_num, y = red_value_num,
        color = Novelty_label, shape = Phylum_plot
      )
    ) +
      geom_point(size = 4, alpha = 0.85, stroke = 0.5) +
      geom_hline(
        yintercept = 0.95, linetype = "dashed",
        color = "gray50", linewidth = 0.6
      ) +
      geom_vline(
        xintercept = 30, linetype = "dashed",
        color = "#E67E22", linewidth = 0.6
      ) +
      scale_color_manual(values = novelty_colors_12, guide = "none") +
      scale_shape_manual(
        values = c(
          "Thermoplasmatota" = 16, "Thermoproteota" = 17,
          "Undinarchaeota" = 15, "Sin asignar" = 4
        ),
        labels = make_label_fn(parent_lookup_phy12, level = "Phylum"),
        name = "Filo"
      ) +
      scale_x_continuous(
        breaks = seq(0, 100, by = 10),
        labels = function(x) paste0(x, "%")
      ) +
      scale_y_continuous(
        breaks = seq(0.85, 1.0, by = 0.02),
        expand = expansion(mult = c(0.03, 0.03))
      ) +
      labs(
        x = "Cobertura de marcadores ar53 (MSA, %)",
        y = "RED value (divergencia evolutiva relativa)"
      ) +
      theme_publication() +
      theme(
        legend.title = element_text(size = 10),
        legend.text  = element_markdown()
      )

    save_fig(p12d, "Fig12D_red_vs_msa_calidad_vs_posicion.png", width = 9, height = 6)
    cat("[OK] Saved: Fig12D_red_vs_msa_calidad_vs_posicion.png\n")
  } else {
    cat("  [WARN] Insuficientes datos para Fig12D – omitiendo\n")
  }
} # end of if (FALSE), section 12

# ---- section 13: Sankey, taxonomic flow from phylum to species ----
# The alluvial diagram shows the 54 MAGs separating as the hierarchy descends.
# One flow per MAG, band thickness by how many share a transition, colour by
# phylum, so a lineage stays traceable the whole way down.
cat("\n=== Sección 13: Sankey plot taxonómico Filo to Especie ===\n")

if (!requireNamespace("ggalluvial", quietly = TRUE)) {
  cat("[WARN] Instalando ggalluvial...\n")
  install.packages("ggalluvial", quiet = TRUE)
}
suppressMessages(library(ggalluvial, quietly = TRUE))
cat("[OK] ggalluvial\n")

# Wide form: one row per MAG, one column per taxonomic level.
df_sankey_wide <- df %>%
  filter(!is.na(Station), !is.na(Depth)) %>%
  transmute(
    Filo    = ifelse(is.na(Phylum), "Sin asignar", Phylum),
    Clase   = ifelse(is.na(Class), "Sin asignar", Class),
    Orden   = ifelse(is.na(Order), "Sin asignar", Order),
    Familia = ifelse(is.na(Family), "Sin asignar", Family),
    Genero  = ifelse(is.na(Genus), "Sin asignar", Genus),
    Especie = ifelse(is.na(Species), "Sin especie asignada", Species)
  )

# ---- Coherent vertical ordering ----
# In ggalluvial, factor level 1 = bottom, last level = top.
# Strategy: sort each level by ascending parent-factor-index (explicit match to
# avoid as.integer() failures when dplyr drops factor structure after summarise),
# then ascending n within parent. Children of the topmost parent end up at top.

order_ascending_n <- function(df, col) {
  df %>%
    count(.data[[col]]) %>%
    arrange(n) %>%
    pull(.data[[col]])
}

order_by_parent <- function(df, parent_col, child_col) {
  # Explicit loop over parents in exact factor level order (bottom to top).
  # For each parent, children are added in ascending-count order.
  # Children shared across parents are placed at the position of their first
  # (lowest) parent. Pure character arithmetic, with no factor codes.
  parent_lvls <- levels(df[[parent_col]])
  par_vec <- as.character(df[[parent_col]])
  chi_vec <- as.character(df[[child_col]])

  result <- character(0)
  for (p in parent_lvls) {
    children <- chi_vec[par_vec == p]
    child_counts <- sort(table(children)) # ascending count
    new_children <- names(child_counts)
    new_children <- new_children[!new_children %in% result]
    result <- c(result, new_children)
  }
  result
}

df_sankey_wide <- df_sankey_wide %>%
  mutate(Filo = factor(Filo, levels = order_ascending_n(., "Filo")))

df_sankey_wide <- df_sankey_wide %>%
  mutate(Clase = factor(Clase, levels = order_by_parent(., "Filo", "Clase")))

df_sankey_wide <- df_sankey_wide %>%
  mutate(Orden = factor(Orden, levels = order_by_parent(., "Clase", "Orden")))

df_sankey_wide <- df_sankey_wide %>%
  mutate(Familia = factor(Familia, levels = order_by_parent(., "Orden", "Familia")))

df_sankey_wide <- df_sankey_wide %>%
  mutate(Genero = factor(Genero, levels = order_by_parent(., "Familia", "Genero")))

# Species: ascending-count order. order_by_parent() here would place
# "Sin especie asignada" at the first genus that uses it, leaving the other
# genera with non-contiguous species positions and visible gaps in their bars.
df_sankey_wide <- df_sankey_wide %>%
  mutate(Especie = factor(Especie, levels = order_ascending_n(., "Especie")))

# ---- Color palette ----
all_taxa_wide <- sort(unique(c(
  as.character(df_sankey_wide$Filo),    as.character(df_sankey_wide$Clase),
  as.character(df_sankey_wide$Orden),   as.character(df_sankey_wide$Familia),
  as.character(df_sankey_wide$Genero),  as.character(df_sankey_wide$Especie)
)))
base_cols_s <- c(
  RColorBrewer::brewer.pal(9, "Set1"),
  RColorBrewer::brewer.pal(8, "Set2"),
  RColorBrewer::brewer.pal(12, "Set3"),
  RColorBrewer::brewer.pal(8, "Dark2"),
  RColorBrewer::brewer.pal(8, "Accent")
)
set.seed(2024)
pal_wide <- setNames(
  colorRampPalette(base_cols_s)(length(all_taxa_wide)),
  all_taxa_wide
)
pal_wide["Sin asignar"] <- "#C0C0C0"
pal_wide["Sin especie asignada"] <- "#B0B0B0"

# Aggregate to unique taxonomic paths, one row per combination, weighted by
# MAG count. This gives one ribbon per path (not one per MAG), so borders on
# the alluvium are at the taxon level, not the individual-MAG level.
df_sankey_agg <- df_sankey_wide %>%
  group_by(Filo, Clase, Orden, Familia, Genero, Especie) %>%
  summarise(n_mags = n(), .groups = "drop")

# Map every taxon name appearing in the sankey to its rank, so the per-stratum
# label callback can apply CDC italics (genus/species italic; family and above
# plain). When a name occurs at more than one rank (e.g. "UBA57" exists as a
# family and a genus in this dataset), the lower rank wins so the genus
# rendering stays italic; the family stratum then displays consistently with
# the genus stratum, which is acceptable since these GTDB code names are kept
# plain by the per-word logic anyway.
sankey_level_map <- local({
  rank_levels <- c(
    Filo = "Phylum", Clase = "Class", Orden = "Order",
    Familia = "Family", Genero = "Genus", Especie = "Species"
  )
  m <- character(0)
  for (col in names(rank_levels)) {
    vals <- unique(na.omit(as.character(df_sankey_agg[[col]])))
    for (v in vals) m[v] <- rank_levels[[col]]
  }
  m
})

# Helper: taxon_italic() applied per-stratum with rank looked up from the map
# above; non-taxonomic labels ("Sin asignar", "Sin especie asignada") pass through.
taxon_italic_sankey <- function(x) {
  out <- character(length(x))
  for (i in seq_along(x)) {
    nm <- x[i]
    if (is.na(nm) || nm %in% c("Sin asignar", "Sin especie asignada")) {
      out[i] <- nm
    } else {
      lvl <- sankey_level_map[nm]
      if (is.na(lvl)) lvl <- "Species"
      out[i] <- taxon_italic(nm, lvl)
    }
  }
  out
}

# Bar width (used also to compute label nudge)
bar_w <- 0.07

p_sankey <- ggplot(
  df_sankey_agg,
  aes(
    y = n_mags,
    axis1 = Filo, axis2 = Clase, axis3 = Orden,
    axis4 = Familia, axis5 = Genero, axis6 = Especie
  )
) +
  # Layer 1: bars with opaque fill, no border.
  geom_stratum(
    aes(fill = after_stat(stratum)),
    width = bar_w, color = NA, na.rm = TRUE
  ) +
  # Layer 2: ribbons drawn over the bars. No later opaque fill covers them,
  # so they are always visible; color=NA means no borders and no false
  # division lines anywhere in the plot.
  geom_flow(
    fill      = gray(0.65),
    color     = gray(0.65),
    linewidth = 0.45,
    width     = bar_w,
    na.rm     = TRUE
  ) +
  # Layer 3: bar borders only. A transparent fill preserves bar and ribbon
  # colors underneath; black outline marks taxon boundaries.
  geom_stratum(
    width = bar_w, fill = NA, color = "black", linewidth = 0.45, na.rm = TRUE
  ) +
  # Markdown-rendered labels to the right; n >= 1 so no taxon is suppressed
  ggtext::geom_richtext(
    stat = "stratum",
    aes(label = ifelse(
      after_stat(count) >= 1,
      paste0(
        taxon_italic_sankey(as.character(after_stat(stratum))),
        "<br><span style='font-style:normal'>(", after_stat(count), ")</span>"
      ),
      ""
    )),
    nudge_x = bar_w / 2 + 0.03,
    hjust = 0,
    size = 4.2,
    color = "black",
    fill = NA,
    label.colour = NA,
    label.padding = unit(c(0, 0, 0, 0), "pt"),
    lineheight = 0.82,
    na.rm = TRUE
  ) +
  scale_x_discrete(
    limits = c("Filo", "Clase", "Orden", "Familia", "Género", "Especie"),
    expand = expansion(add = c(0.1, 1.0))
  ) +
  scale_fill_manual(values = pal_wide, na.value = "#C0C0C0") +
  labs(x = NULL, y = NULL) +
  theme_void() +
  theme(
    axis.text.x     = element_text(face = "bold", size = 13, margin = margin(t = 2)),
    legend.position = "none"
  )

save_fig(p_sankey, "Fig13_sankey_taxonomico_filo_especie.png", width = 17, height = 11)

# === SECTION 1c (v2): combined sankey / order / species panel ===
# A variant of Fig1_panel_filo_orden_especie where subpanel (a) shows the
# Sankey instead of the stacked phylum bars. Subpanels (b) and (c) are reused
# unchanged from the original panel. It is assembled here rather than in section
# 1c because p_sankey is defined in
# section 13, and p_panel_b and p_panel_c are still in the global environment
# from section 1.
cat("\n=== Sección 1c (v2): Panel sankey + orden + especie ===\n")

{
  n_depths_panel <- length(unique(df_bar$Depth_label))
  h_top <- 4 + n_depths_panel * 2 # natural height of pC (order)
  h_bot <- 6 # natural height of pD (species)

  n_stations_panel <- length(unique(df_bar$Station))
  w_c_natural <- 4.5 + n_stations_panel * 1.9 # natural width of pD (species)
  w_panel_v2 <- max(w_c_natural, 17) # at least the natural width of the sankey

  # The Sankey takes about 63% of the top row, 17 against 10, which is its
  # natural proportion beside the order bars. Give it less and patchwork
  # compresses it until the labels overlap.
  p_combined_panel_v2 <- (p_sankey | p_panel_b) +
    plot_layout(widths = c(17, 10))
  p_combined_panel_v2 <- p_combined_panel_v2 / p_panel_c +
    plot_layout(heights = c(h_top, h_bot)) +
    plot_annotation(tag_levels = "a") &
    theme(plot.tag = element_text(face = "bold", size = 20, color = "black"))

  save_fig(p_combined_panel_v2, "Fig1_panel_sankey_orden_especie.png",
    width = w_panel_v2, height = h_top + h_bot
  )
}

cat("[OK] Sección 13 completada\n\n")

# === Figure index ===
cat("Figuras generadas:\n")
cat("  Sección 1 (Barras dobles – estación × profundidad y profundidad × estación):\n")
cat("    Fig1C_barras_{filo,orden,familia,genero,especie}_estacion_x_profundidad.png\n")
cat("    Fig1D_barras_{filo,orden,familia,genero,especie}_profundidad_x_estacion.png\n")
cat("  Sección 1b (MAGs sin especie – por género):\n")
cat("    Fig1E_novel_sp_estacion_x_profundidad.png\n")
cat("    Fig1F_novel_sp_profundidad_x_estacion.png\n")
cat("  Sección 2 (Distribución de clases):\n")
cat("    Fig2C_clases_estacion_profundidad.png\n")
cat("    Fig2D_clases_profundidad_estacion.png\n")
cat("  Sección 7 (Gráficos de burbujas – estación × profundidad):\n")
cat("    Fig7_burbujas_familia_estacion_profundidad.png\n")
cat("    Fig7b_burbujas_clase_estacion_profundidad.png\n")
cat("    Fig7c_burbujas_orden_estacion_profundidad.png\n")
cat("    Fig7d_burbujas_genero_estacion_profundidad.png\n")
cat("  Sección 8 (Riqueza por muestra):\n")
cat("    Fig8_riqueza_MAGs_por_muestra.png\n")
cat("  Sección 10 (skani – novedad taxonómica cuantitativa; sólo MAGs sin especie):\n")
cat("    Fig10A_opC_dual_panel.png           – ANI + fracción alineada (panel doble)\n")
cat("    Fig10B_skani_ani_vs_af_scatter.png  – nube ANI × AF para todos los hits GTDB\n")
cat("    Fig11_skani_ani_novel_vs_coleccion.png – heatmap ANI novel MAGs × colección\n")
cat("  Sección 12 (parámetros filogenómicos GTDB-Tk):\n")
cat("    Fig12A_metodo_clasificacion_por_filo.png\n")
cat("    Fig12A_clase_metodo_clasificacion_por_clase.png\n")
cat("    Fig12B_red_value_por_clase.png\n")
cat("    Fig12C_msa_percent_por_clase.png\n")
cat("    Fig12D_red_vs_msa_calidad_vs_posicion.png\n")
cat("  Sección 13 (Sankey taxonómico – Filo to Especie):\n")
cat("    Fig13_sankey_taxonomico_filo_especie.png\n")
cat("  Panel combinado filo/orden/especie (dos variantes):\n")
cat("    Fig1_panel_filo_orden_especie.png       (v1: barras filo + orden + especie)\n")
cat("    Fig1_panel_sankey_orden_especie.png     (v2: sankey + barras orden + especie)\n")
cat("  Tablas:\n")
cat("    tax_analysis_summary.tsv\n")
cat("    tax_exclusividad_resumen.xlsx (hojas: Genero_estacion, Filo_estacion,\n")
cat("      Clase_estacion, Orden_estacion, Familia_estacion)\n")
cat("    tabla_calidad_10MAGs.tsv / .png (MAGs sin especie)\n")

cat(sprintf("\n=== Finalización del script: %s ===\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf(
  "    Duración total: %.1f minutos\n",
  as.numeric(difftime(Sys.time(), script_start, units = "mins"))
))
