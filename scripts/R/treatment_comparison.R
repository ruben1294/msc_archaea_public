#!/usr/bin/env Rscript

# Script: treatment_comparison.R
# Aggregated trait heatmaps by group (figures 9a to 9f), Fisher tests and
# logistic models per functional trait, over both the DRAM2 traits and the
# carbon, nitrogen and sulfur cycle pathways from CNPS.cycle.
# Fecha: 2026-02-26

# === PAQUETES ===
if (!require("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(
  tidyverse, ComplexHeatmap, RColorBrewer, circlize,
  readxl, broom, grid, showtext, sysfonts, ggtext,
  install = TRUE, update = FALSE
)

# ---- fonts, same approach as tax_analysis.R ----
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

# Cargar fuentes locales directamente (prioridad: Lato > Source Sans 3)
font_found <- FALSE
local_font_dir <- path.expand("~/.local/share/fonts")
local_direct <- list(
  list(
    name    = "Lato",
    regular = "Lato-Regular.ttf",
    bold    = "Lato-Bold.ttf",
    italic  = "Lato-Italic.ttf",
    bolditalic = "Lato-BoldItalic.ttf"
  ),
  list(
    name    = "Source Sans 3",
    regular = "SourceSans3-Regular.otf",
    bold    = "SourceSans3-Bold.otf",
    italic  = "SourceSans3-It.otf",
    bolditalic = "SourceSans3-BoldIt.otf"
  )
)
if (dir.exists(local_font_dir)) {
  for (lf in local_direct) {
    reg_path <- file.path(local_font_dir, lf$regular)
    if (file.exists(reg_path)) {
      tryCatch({
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
      }, error = function(e) {
        cat("[WARN] Archivos de", lf$name, "encontrados pero no se pudieron cargar:", conditionMessage(e), "\n")
      })
    }
  }
}

# Fallback: look through the system font registry
for (font_name in font_priority) {
  matching_fonts <- available_fonts[grepl(font_name, available_fonts$family, ignore.case = TRUE), ]
  if (nrow(matching_fonts) > 0) {
    tryCatch({
      font_path <- matching_fonts$file[1]
      sysfonts::font_add(family = "custom_font", regular = font_path)
      font_family <- "custom_font"
      cat("[OK] Usando fuente:", font_name, "\n")
      font_found <- TRUE
      break
    }, error = function(e) {})
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
    tryCatch({
      sysfonts::font_add_google(gf$name, "custom_font")
      font_family <- "custom_font"
      cat("[OK] Descargada y usando:", gf$name, "(Google Fonts)\n")
      font_found <- TRUE
      break
    }, error = function(e) {})
  }
  if (!font_found) {
    cat("[WARN] No se pudo descargar de Google Fonts. Usando fuente sans-serif por defecto.\n")
  }
}
cat("[OK] Configuración de fuentes completada\n\n")

# === RUTAS ===
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

BASE_DIR    <- resolve_project_dir()
DISTILL_DIR <- file.path(BASE_DIR, "results/08_MAG_distill")
CNPS_DIR    <- file.path(BASE_DIR, "results/15_CNPS_cycle")

METADATA_FILE <- file.path(BASE_DIR, "MAGs_metadata.tsv")
TRAITS_FILE   <- file.path(DISTILL_DIR, "SUMMARIZE/traits.xlsx")

# Directorio de salida propio (subdirectorio de funct_annot)
OUTPUT_DIR  <- file.path(BASE_DIR, "results/SO_01/funct_annot/treatment_comparison")
CACHE_FILE  <- file.path(BASE_DIR, "results/SO_01/funct_annot/translation_cache_dram2_distill.csv")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("Directorio de salida:", OUTPUT_DIR, "\n\n")

# === METADATA ===
metadata <- read.delim(METADATA_FILE, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

mag_to_depth   <- setNames(metadata$Depth,   metadata$NameWithoutExtension)
mag_to_station <- setNames(metadata$Station, metadata$NameWithoutExtension)
mag_to_genus   <- setNames(as.character(metadata$Genus), metadata$NameWithoutExtension)
mag_to_genus[is.na(mag_to_genus) | trimws(mag_to_genus) %in% c("", "N/A")] <- "Desconocido"

get_depth   <- function(x) mag_to_depth[x]
get_station <- function(x) mag_to_station[x]
get_genus   <- function(x) {
  g <- mag_to_genus[x]
  g[is.na(g)] <- "Desconocido"
  as.character(g)
}

capitalize_first <- function(x) {
  ifelse(is.na(x), x, paste0(toupper(substring(x, 1, 1)), substring(x, 2)))
}

# === GLOBALES DE PROFUNDIDAD / ESTACIÓN ===
# Depth labels for the figures: ROMZ to ZMOR, 1000 to "1000 m", BTM to Fondo
depth_label_map <- c("ROMZ" = "ZMOR", "1000" = "1000 m", "BTM" = "Fondo")
translate_depth <- function(x) {
  out <- depth_label_map[as.character(x)]
  out[is.na(out)] <- as.character(x)[is.na(out)]
  unname(out)
}
depth_levels_std <- c("ZMOR", "1000 m", "Fondo")
depth_colors_std <- setNames(
  RColorBrewer::brewer.pal(3, "Set2"),
  depth_levels_std
)

station_levels <- sort(unique(as.character(metadata$Station[!is.na(metadata$Station)])))
n_st <- max(3L, min(length(station_levels), 9L))
station_colors <- setNames(
  RColorBrewer::brewer.pal(n_st, "Set1")[seq_along(station_levels)],
  station_levels
)

# ---- scientific name formatting ----
ITALIC_LEVELS <- c("Genus", "Species", "Subspecies", "Strain")

should_italicize_level <- function(level) {
  !is.null(level) && length(level) == 1 && !is.na(level) && level %in% ITALIC_LEVELS
}

# Italicises a taxon name as markdown (*name*) while leaving GTDB codes upright:
# UBA60, MGIII, sp002507125 and the like. CDC style, so genus and species are
# italic and everything above family is not.
taxon_italic <- function(x, level = "Species") {
  if (!should_italicize_level(level)) {
    return(ifelse(is.na(x), x, as.character(x)))
  }
  sapply(x, function(nm) {
    if (is.na(nm) || nm == "Sin asignar" || nm == "Desconocido") return(nm)
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

# Pull one rank out of a GTDB string ("d__Archaea;p__...;o__Foo;...")
parse_gtdb_rank <- function(taxstr, prefix) {
  m <- regmatches(taxstr, regexpr(paste0(prefix, "[^;]*"), taxstr))
  out <- rep(NA_character_, length(taxstr))
  matched <- nchar(m) > 0
  out[matched] <- sub(paste0("^", prefix), "", m[matched])
  out[is.na(out) | out == ""] <- NA_character_
  out
}

# === MANUAL TRANSLATIONS + CACHE ===
.MANUAL_TRANSLATIONS <- list(
  # DRAM2 traits
  "Aerobe"          = "Aerobio",
  "Anaerobe"        = "Anaerobio",
  "aerobe"          = "aerobio",
  "anaerobe"        = "anaerobio",
  "C1 metabolism"   = "Metabolismo C1",
  "Carbon fixation" = "Fijación de carbono",
  "Complex I"       = "Complejo I",
  "Complex II"      = "Complejo II",
  "Complex III"     = "Complejo III",
  "Complex IV"      = "Complejo IV",
  "Complex IV H"    = "Complejo IV-H",
  "Complex IV L"    = "Complejo IV-L",
  "Complex IV-H"    = "Complejo IV-H",
  "Complex IV-L"    = "Complejo IV-L",
  "Complex V"       = "Complejo V",
  "Nitrogen"        = "Nitrógeno",
  "Sulfur"          = "Azufre",
  "Photo"           = "Fotosíntesis",
  "Other"           = "Otro",
  "SCFA/alcohol"    = "AGCC/alcohol",
  "CAZy"            = "CAZy",
  # CNPS cycle pathways (Carbono)
  "Aerobic C fixation"              = "Fijación aeróbica de C",
  "Aerobic CH4 oxidation"           = "Oxidación aeróbica de CH\u2084",
  "Aerobic respiration"             = "Respiración aeróbica",
  "Anaerobic C fixation"            = "Fijación anaeróbica de C",
  "CO oxidation"                    = "Oxidación de CO",
  "Fermentation"                    = "Fermentación",
  "Methanogenesis"                  = "Metanogénesis",
    # CNPS nitrogen cycle pathways
  "Nitrification"                   = "Nitrificación",
  "Denitrification"                 = "Desnitrificación",
  "Anammox"                         = "Anammox",
  "Nitrogen fixation"               = "Fijación de N\u2082",
  "Assimilatory nitrate reduction"  = "Reducción asimilativa de NO\u2083\u207b",
  "Dissimilatory nitrate reduction" = "Reducción disimilativa de NO\u2083\u207b",
  "Nitrogen mineralization"         = "Mineralización del N",
  "Nitrogen assimilation"           = "Asimilación de N",
  "DNRA"                            = "DNRA",
  # CNPS cycle pathways (Azufre)
  "Dissimilatory sulfate reduction" = "Reducción disimilativa de SO\u2084\u00b2\u207b",
  "Assimilatory sulfate reduction"  = "Reducción asimilativa de SO\u2084\u00b2\u207b",
  "Sulfur oxidation"                = "Oxidación de S",
  "Sulfide oxidation"               = "Oxidación de sulfuro",
  "Sulfur mineralization"           = "Mineralización de S",
  "Sulfate assimilation"            = "Asimilación de SO\u2084\u00b2\u207b",
  "Sulfite assimilation"            = "Asimilación de SO\u2083\u00b2\u207b",
  "Tetrathionate cycling"           = "Ciclo del tetrationato",
  "Sulfide cycling"                 = "Ciclo del sulfuro",
  "SOX"                             = "Oxidación de tiosulfato (SOX)"
)

translate_to_spanish <- function(names_vec, translation_cache_file = CACHE_FILE) {
  translated        <- names_vec
  translation_cache <- list()
  if (file.exists(translation_cache_file)) {
    cache_df <- tryCatch(
      read.csv(translation_cache_file, stringsAsFactors = FALSE),
      error = function(e) NULL
    )
    if (!is.null(cache_df) &&
        "English" %in% names(cache_df) && "Spanish" %in% names(cache_df)) {
      translation_cache <- setNames(as.list(cache_df$Spanish), cache_df$English)
    }
  }
  for (i in seq_along(names_vec)) {
    original <- names_vec[i]
    if (original %in% names(.MANUAL_TRANSLATIONS)) {
      translated[i] <- .MANUAL_TRANSLATIONS[[original]]
    } else if (original %in% names(translation_cache)) {
      translated[i] <- as.character(translation_cache[[original]])
    }
    # no DeepL call here, so the script runs without network access
  }
  translated
}

# HELPER: cargar CSV de abundancia CNPS y convertir a data.frame de rasgos
#
# Los CSV de CNPS tienen formato:  Type | MAG1 | MAG2 | ...
# El resultado tiene formato:      genome | pathway1 | pathway2 | ...
# Los valores son binarios (presencia/ausencia: valor > 0 to 1).

load_cnps_abundance <- function(csv_path) {
  if (!file.exists(csv_path)) {
    cat("  Archivo no encontrado:", csv_path, "\n")
    return(NULL)
  }
  df <- read.csv(csv_path, stringsAsFactors = FALSE, check.names = FALSE)
  if (!"Type" %in% colnames(df)) {
    cat("  Sin columna 'Type' en:", basename(csv_path), "\n")
    return(NULL)
  }
  type_col <- df$Type
  mag_cols <- setdiff(colnames(df), "Type")
  if (length(mag_cols) == 0) return(NULL)

  mat_pw  <- as.matrix(df[, mag_cols, drop = FALSE])
  rownames(mat_pw) <- type_col

  # Transponer: filas = MAGs, columnas = pathways to binarizar
  mat_t   <- t(mat_pw)                           # MAGs × pathways
  mat_bin <- matrix(
    as.integer(mat_t > 0),
    nrow = nrow(mat_t), dimnames = dimnames(mat_t)
  )

  out <- as.data.frame(mat_bin, stringsAsFactors = FALSE, check.names = FALSE)
  out <- cbind(genome = rownames(out), out)
  rownames(out) <- NULL
  out
}

# ---- ---- aggregated trait heatmap: percent of MAGs carrying each trait ----

make_trait_aggregated_heatmap <- function(trait_df, group_by_vec, group_label_title,
                                          title_es, outfile,
                                          group_order  = NULL,
                                          cluster_rows = FALSE,
                                          cache_file   = CACHE_FILE,
                                          row_level    = NULL) {
  tryCatch({
    genome_col <- names(trait_df)[1]
    genomes    <- trait_df[[genome_col]]
    mat_raw    <- trait_df %>% select(-all_of(genome_col)) %>% as.matrix()
    mat_bin    <- matrix(as.numeric(mat_raw > 0), nrow = nrow(mat_raw),
                         dimnames = list(genomes, colnames(mat_raw)))

    keep <- colSums(mat_bin, na.rm = TRUE) > 0
    mat_bin <- mat_bin[, keep, drop = FALSE]
    if (ncol(mat_bin) == 0) {
      cat("  \u26a0 Sin traits para:", title_es, "\n")
      return(invisible(NULL))
    }

    groups <- group_by_vec[genomes]
    groups[is.na(groups)] <- "Desconocido"

    n_per_group <- table(groups)

    df_agg <- as.data.frame(mat_bin)
    df_agg$group <- groups
    agg <- df_agg %>%
      group_by(group) %>%
      summarise(across(everything(), ~ mean(. > 0, na.rm = TRUE) * 100), .groups = "drop") %>%
      column_to_rownames("group")

    keep2 <- colSums(agg, na.rm = TRUE) > 0
    agg   <- agg[, keep2, drop = FALSE]
    if (ncol(agg) == 0) return(invisible(NULL))

    if (!is.null(group_order) && !cluster_rows) {
      valid_ord <- intersect(group_order, rownames(agg))
      extra     <- setdiff(rownames(agg), group_order)
      agg       <- agg[c(valid_ord, extra), , drop = FALSE]
    }

    row_labels <- paste0(rownames(agg), " (N=", as.integer(n_per_group[rownames(agg)]), ")")
    # Row labels: only the genus goes italic
    rn_face <- if (should_italicize_level(row_level)) "italic" else "plain"

    trait_es <- capitalize_first(
      translate_to_spanish(colnames(agg), translation_cache_file = cache_file)
    )
    trait_es[is.na(trait_es)] <- colnames(agg)[is.na(trait_es)]
    colnames(agg) <- trait_es
    agg_mat        <- as.matrix(agg)
    rownames(agg_mat) <- row_labels

    col_fun <- colorRamp2(c(0, 50, 100), c("white", "#74C0FC", "#08306B"))

    ht <- Heatmap(
      agg_mat,
      name            = "% MAGs\ncon rasgo",
      col             = col_fun,
      cluster_rows    = cluster_rows,
      cluster_columns = TRUE,
      clustering_distance_rows    = "euclidean",
      clustering_distance_columns = "euclidean",
      show_row_dend   = cluster_rows,
      row_dend_side   = "left",
      show_row_names  = TRUE,
      show_column_names = TRUE,
      row_names_gp    = gpar(fontsize = 10, fontface = rn_face),
      column_names_gp = gpar(fontsize = 9),
      column_names_rot = 45,
      rect_gp  = gpar(col = "black", lwd = 0.4),
      border   = TRUE,
      row_title    = group_label_title,
      row_title_gp = gpar(fontsize = 10, fontface = "bold"),
      column_title    = title_es,
      column_title_gp = gpar(fontsize = 11, fontface = "bold"),
      heatmap_legend_param = list(
        title     = "% MAGs con\nel rasgo",
        at        = c(0, 25, 50, 75, 100),
        labels    = c("0%", "25%", "50%", "75%", "100%"),
        title_gp  = gpar(fontsize = 10, fontface = "bold"),
        labels_gp = gpar(fontsize = 9),
        border    = "black",
        legend_direction = "vertical"
      ),
      cell_fun = function(j, i, x, y, width, height, fill) {
        val <- agg_mat[i, j]
        if (!is.na(val) && val > 0)
          grid.text(sprintf("%.0f%%", val), x, y,
                    gp = gpar(fontsize = 7,
                              col = if (val > 60) "white" else "black"))
      }
    )

    n_rows <- nrow(agg_mat)
    n_cols <- ncol(agg_mat)
    w_ht   <- max(10, round(n_cols * 0.45 + 3))
    h_ht   <- max(4,  round(n_rows * 0.5  + 2))
    png(outfile, width = w_ht, height = h_ht, units = "in", res = 300)
    draw(ht)
    dev.off()
    cat("  \u2713 Guardado:", basename(outfile), "\n")
  }, error = function(e) cat("  \u26a0 Error en", title_es, ":", conditionMessage(e), "\n"))
}

# ---- HELPER: Forest plot de odds ratios ----
plot_forest <- function(results_df, title, outfile,
                        x_col     = "OR",
                        xmin_col  = "CI_low",
                        xmax_col  = "CI_high",
                        y_col     = "label",
                        color_col = "significant",
                        facet_col = NULL) {
  if (nrow(results_df) == 0) {
    cat("  \u26a0 Sin resultados para:", title, "\n")
    return(invisible(NULL))
  }

  df_plot <- results_df %>%
    filter(is.finite(.data[[x_col]]),
           .data[[x_col]] > 0,
           is.finite(.data[[xmin_col]]),
           is.finite(.data[[xmax_col]]))

  if (nrow(df_plot) == 0) {
    cat("  \u26a0 Todos los OR son infinitos o cero:", title, "\n")
    return(invisible(NULL))
  }

  p <- ggplot(df_plot, aes(x = .data[[x_col]], y = .data[[y_col]])) +
    geom_vline(xintercept = 1, linetype = "dashed", colour = "grey50") +
    geom_errorbarh(
      aes(xmin = .data[[xmin_col]], xmax = .data[[xmax_col]]),
      height    = 0.25,
      linewidth = 0.5
    ) +
    geom_point(aes(color = .data[[color_col]]), size = 2) +
    scale_x_log10() +
    scale_color_manual(
      values = c("TRUE" = "#D62728", "FALSE" = "#1F77B4"),
      labels = c("TRUE" = "p.adj < 0.05", "FALSE" = "p.adj \u2265 0.05"),
      name   = "Significancia"
    ) +
    labs(
      title = title,
      x     = "Odds Ratio (IC 95%, escala log)",
      y     = NULL
    ) +
    theme_bw(base_size = 10, base_family = font_family) +
    theme(
      axis.text.y        = ggtext::element_markdown(size = 8),
      strip.text         = ggtext::element_markdown(face = "bold", size = 9),
      panel.grid.major.y = element_line(colour = "grey90"),
      panel.grid.minor   = element_blank(),
      plot.title         = element_text(face = "bold", hjust = 0.5, size = 11)
    )

  if (!is.null(facet_col)) {
    p <- p + facet_wrap(vars(.data[[facet_col]]), scales = "free_y", ncol = 3)
  }

  n_rows_plot <- if (is.null(facet_col)) nrow(df_plot) else
    max(table(df_plot[[facet_col]]))
  n_facets    <- if (is.null(facet_col)) 1L else
    length(unique(df_plot[[facet_col]]))
  n_cols_plot <- min(3L, n_facets)
  n_rows_grid <- ceiling(n_facets / n_cols_plot)

  w_f <- max(8,  round(n_cols_plot * 4.5))
  h_f <- max(6,  round(n_rows_grid * (n_rows_plot * 0.28 + 1.5)))

  ggsave(outfile, plot = p, width = w_f, height = h_f, dpi = 300, limitsize = FALSE)
  cat("  \u2713 Guardado:", basename(outfile), "\n")
}

# ---- CARGAR DATOS ----

# ---- DRAM2 traits ----
if (!file.exists(TRAITS_FILE)) stop("traits.xlsx no encontrado: ", TRAITS_FILE)
traits_energy <- read_excel(TRAITS_FILE, sheet = "energy")
traits_carbon <- read_excel(TRAITS_FILE, sheet = "carbon")
cat("DRAM2 traits cargados\n")

# ---- CNPS cycle abundance: CSV to a binary MAG by pathway data.frame ----
cat("Cargando datos CNPS.cycle...\n")
cnps_carbon   <- load_cnps_abundance(file.path(CNPS_DIR, "carbon",   "C_abundance.csv"))
cnps_nitrogen <- load_cnps_abundance(file.path(CNPS_DIR, "nitrogen", "N_abundance.csv"))
cnps_sulfur   <- load_cnps_abundance(file.path(CNPS_DIR, "sulfur",   "S_abundance.csv"))

if (!is.null(cnps_carbon))   cat("  CNPS carbono:   ", nrow(cnps_carbon),   "MAGs,", ncol(cnps_carbon)   - 1, "pathways\n")
if (!is.null(cnps_nitrogen)) cat("  CNPS nitrógeno: ", nrow(cnps_nitrogen), "MAGs,", ncol(cnps_nitrogen) - 1, "pathways\n")
if (!is.null(cnps_sulfur))   cat("  CNPS azufre:    ", nrow(cnps_sulfur),   "MAGs,", ncol(cnps_sulfur)   - 1, "pathways\n")

# Grouping vectors, genome_id to group
depth_grp   <- setNames(translate_depth(get_depth(metadata$NameWithoutExtension)),
                        metadata$NameWithoutExtension)
station_grp <- setNames(as.character(get_station(metadata$NameWithoutExtension)),
                        metadata$NameWithoutExtension)
genus_grp   <- setNames(get_genus(metadata$NameWithoutExtension),
                        metadata$NameWithoutExtension)

# Taxonomic order, read out of FullTaxonomy in the metadata
order_raw <- parse_gtdb_rank(as.character(metadata$FullTaxonomy), "o__")
order_raw[is.na(order_raw) | trimws(order_raw) == ""] <- "Desconocido"
order_grp <- setNames(order_raw, metadata$NameWithoutExtension)
cat("Órdenes detectados:", paste(unique(order_grp), collapse = ", "), "\n")

# SECCIÓN 1: Heatmaps agregados
# Fig9a-9f:  DRAM2 traits (energy / carbon)
# Fig12a-12i: CNPS cycle (carbon / nitrogen / sulfur)
cat("=== SECCIÓN 1: HEATMAPS AGREGADOS ===\n\n")

# ---- fig9a-9f: DRAM2 traits ----
make_trait_aggregated_heatmap(
  traits_energy, depth_grp, "Profundidad",
  "Rasgos energéticos por nivel de profundidad\nDRAM2 – traits.xlsx",
  file.path(OUTPUT_DIR, "Fig9a_traits_energy_by_depth.png"),
  group_order = depth_levels_std, cluster_rows = FALSE
)
make_trait_aggregated_heatmap(
  traits_energy, station_grp, "Estación",
  "Rasgos energéticos por estación de muestreo\nDRAM2 – traits.xlsx",
  file.path(OUTPUT_DIR, "Fig9b_traits_energy_by_station.png"),
  cluster_rows = TRUE
)
make_trait_aggregated_heatmap(
  traits_energy, genus_grp, "Género",
  "Rasgos energéticos por género taxonómico\nDRAM2 – traits.xlsx",
  file.path(OUTPUT_DIR, "Fig9c_traits_energy_by_genus.png"),
  cluster_rows = TRUE, row_level = "Genus"
)
make_trait_aggregated_heatmap(
  traits_energy, order_grp, "Orden",
  "Rasgos energéticos por orden taxonómico\nDRAM2 – traits.xlsx",
  file.path(OUTPUT_DIR, "Fig9c2_traits_energy_by_order.png"),
  cluster_rows = TRUE, row_level = "Order"
)
make_trait_aggregated_heatmap(
  traits_carbon, depth_grp, "Profundidad",
  "Utilización de carbono por nivel de profundidad\nDRAM2 – traits.xlsx",
  file.path(OUTPUT_DIR, "Fig9d_traits_carbon_by_depth.png"),
  group_order = depth_levels_std, cluster_rows = FALSE
)
make_trait_aggregated_heatmap(
  traits_carbon, station_grp, "Estación",
  "Utilización de carbono por estación de muestreo\nDRAM2 – traits.xlsx",
  file.path(OUTPUT_DIR, "Fig9e_traits_carbon_by_station.png"),
  cluster_rows = TRUE
)
make_trait_aggregated_heatmap(
  traits_carbon, genus_grp, "Género",
  "Utilización de carbono por género taxonómico\nDRAM2 – traits.xlsx",
  file.path(OUTPUT_DIR, "Fig9f_traits_carbon_by_genus.png"),
  cluster_rows = TRUE, row_level = "Genus"
)
make_trait_aggregated_heatmap(
  traits_carbon, order_grp, "Orden",
  "Utilización de carbono por orden taxonómico\nDRAM2 – traits.xlsx",
  file.path(OUTPUT_DIR, "Fig9f2_traits_carbon_by_order.png"),
  cluster_rows = TRUE, row_level = "Order"
)

# ---- fig12a-12l: CNPS cycle (depth, station, genus, order) ----
cnps_heatmap_specs <- list(
  # CNPS Carbono
  list(df = cnps_carbon,   label = "Carbono (CNPS)",   fig_pfx = "Fig12a", grp = depth_grp,
       grp_title = "Profundidad", grp_lbl = "profundidad",
       ord = depth_levels_std, cluster = FALSE, row_level = NULL),
  list(df = cnps_carbon,   label = "Carbono (CNPS)",   fig_pfx = "Fig12b", grp = station_grp,
       grp_title = "Estación",    grp_lbl = "estación",      ord = NULL, cluster = TRUE,
       row_level = NULL),
  list(df = cnps_carbon,   label = "Carbono (CNPS)",   fig_pfx = "Fig12c", grp = genus_grp,
       grp_title = "Género",      grp_lbl = "género",        ord = NULL, cluster = TRUE,
       row_level = "Genus"),
  list(df = cnps_carbon,   label = "Carbono (CNPS)",   fig_pfx = "Fig12c2", grp = order_grp,
       grp_title = "Orden",       grp_lbl = "orden",         ord = NULL, cluster = TRUE,
       row_level = "Order"),
  # CNPS nitrogen
  list(df = cnps_nitrogen, label = "Nitrógeno (CNPS)", fig_pfx = "Fig12d", grp = depth_grp,
       grp_title = "Profundidad", grp_lbl = "profundidad",
       ord = depth_levels_std, cluster = FALSE, row_level = NULL),
  list(df = cnps_nitrogen, label = "Nitrógeno (CNPS)", fig_pfx = "Fig12e", grp = station_grp,
       grp_title = "Estación",    grp_lbl = "estación",      ord = NULL, cluster = TRUE,
       row_level = NULL),
  list(df = cnps_nitrogen, label = "Nitrógeno (CNPS)", fig_pfx = "Fig12f", grp = genus_grp,
       grp_title = "Género",      grp_lbl = "género",        ord = NULL, cluster = TRUE,
       row_level = "Genus"),
  list(df = cnps_nitrogen, label = "Nitrógeno (CNPS)", fig_pfx = "Fig12f2", grp = order_grp,
       grp_title = "Orden",       grp_lbl = "orden",         ord = NULL, cluster = TRUE,
       row_level = "Order"),
  # CNPS Azufre
  list(df = cnps_sulfur,   label = "Azufre (CNPS)",    fig_pfx = "Fig12g", grp = depth_grp,
       grp_title = "Profundidad", grp_lbl = "profundidad",
       ord = depth_levels_std, cluster = FALSE, row_level = NULL),
  list(df = cnps_sulfur,   label = "Azufre (CNPS)",    fig_pfx = "Fig12h", grp = station_grp,
       grp_title = "Estación",    grp_lbl = "estación",      ord = NULL, cluster = TRUE,
       row_level = NULL),
  list(df = cnps_sulfur,   label = "Azufre (CNPS)",    fig_pfx = "Fig12i", grp = genus_grp,
       grp_title = "Género",      grp_lbl = "género",        ord = NULL, cluster = TRUE,
       row_level = "Genus"),
  list(df = cnps_sulfur,   label = "Azufre (CNPS)",    fig_pfx = "Fig12i2", grp = order_grp,
       grp_title = "Orden",       grp_lbl = "orden",         ord = NULL, cluster = TRUE,
       row_level = "Order")
)

for (spec in cnps_heatmap_specs) {
  if (is.null(spec$df)) next
  make_trait_aggregated_heatmap(
    spec$df, spec$grp, spec$grp_title,
    paste0("Ciclo del ", spec$label, " por ", spec$grp_lbl,
           "\nCNPS.cycle – ", spec$fig_pfx),
    file.path(OUTPUT_DIR, paste0(spec$fig_pfx, "_cnps_",
                                 tolower(gsub(" .*", "", spec$label)), "_by_",
                                 spec$grp_lbl, ".png")),
    group_order  = spec$ord,
    cluster_rows = spec$cluster,
    row_level    = spec$row_level
  )
}

cat("\u2713 Sección 1 completada\n\n")

# ---- SECCIÓN 2: Test de Fisher (rasgo × nivel de grupo) ----
cat("=== SECCIÓN 2: TEST DE FISHER ===\n\n")

run_fisher_tests <- function(trait_df, group_vec, sheet_name, group_var_name) {
  genome_col <- names(trait_df)[1]
  genomes    <- trait_df[[genome_col]]
  mat_raw    <- trait_df %>% select(-all_of(genome_col)) %>% as.matrix()
  mat_bin    <- matrix(as.numeric(mat_raw > 0), nrow = nrow(mat_raw),
                       dimnames = list(genomes, colnames(mat_raw)))

  groups       <- group_vec[genomes]
  groups[is.na(groups)] <- "Desconocido"
  group_levels <- sort(unique(groups))

  results <- vector("list", length(colnames(mat_bin)) * length(group_levels))
  k <- 0L
  for (trait in colnames(mat_bin)) {
    y <- mat_bin[, trait]
    for (lvl in group_levels) {
      grp_bin <- as.integer(groups == lvl)
      tbl <- table(
        factor(grp_bin, levels = c(0L, 1L)),
        factor(y,       levels = c(0L, 1L))
      )
      ft <- tryCatch(fisher.test(tbl), error = function(e) NULL)
      if (is.null(ft)) next
      k <- k + 1L
      results[[k]] <- tibble(
        sheet     = sheet_name,
        group_var = group_var_name,
        level     = lvl,
        trait     = trait,
        OR        = as.numeric(ft$estimate),
        CI_low    = ft$conf.int[1],
        CI_high   = ft$conf.int[2],
        p_value   = ft$p.value
      )
    }
  }
  if (k == 0L) return(NULL)
  bind_rows(results[seq_len(k)])
}

# Lista de fuentes de datos: DRAM2 traits + CNPS cycle
sheet_list <- list(
  list(df = traits_energy, name = "energy",         label = "Energía (DRAM2)"),
  list(df = traits_carbon, name = "carbon_dram2",   label = "Carbono (DRAM2)")
)
# Add the CNPS cycles when their data is present
if (!is.null(cnps_carbon))   sheet_list[[length(sheet_list)+1]] <- list(df = cnps_carbon,   name = "carbon_cnps",   label = "Carbono (CNPS)")
if (!is.null(cnps_nitrogen)) sheet_list[[length(sheet_list)+1]] <- list(df = cnps_nitrogen, name = "nitrogen_cnps", label = "Nitrógeno (CNPS)")
if (!is.null(cnps_sulfur))   sheet_list[[length(sheet_list)+1]] <- list(df = cnps_sulfur,   name = "sulfur_cnps",   label = "Azufre (CNPS)")

grp_list <- list(
  list(vec = depth_grp,   name = "depth",   label = "Profundidad", level_rank = NULL),
  list(vec = station_grp, name = "station", label = "Estación",    level_rank = NULL),
  list(vec = genus_grp,   name = "genus",   label = "Género",      level_rank = "Genus"),
  list(vec = order_grp,   name = "order",   label = "Orden",       level_rank = "Order")
)

fisher_parts <- list()
for (sh in sheet_list) {
  for (gv in grp_list) {
    cat("  Fisher:", sh$label, "×", gv$label, "...")
    res <- run_fisher_tests(sh$df, gv$vec, sh$name, gv$name)
    if (!is.null(res)) fisher_parts[[length(fisher_parts) + 1]] <- res
    cat(" OK\n")
  }
}

if (length(fisher_parts) > 0) {
  fisher_all <- bind_rows(fisher_parts) %>%
    group_by(sheet, group_var) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    ungroup() %>%
    mutate(
      significant   = p_adj < 0.05,
      level_display = dplyr::case_when(
        group_var == "genus" ~ taxon_italic(level, level = "Genus"),
        group_var == "order" ~ taxon_italic(level, level = "Order"),
        TRUE                 ~ as.character(level)
      ),
      label = paste0(trait, " : ", level_display)
    )
  write_csv(fisher_all, file.path(OUTPUT_DIR, "fisher_results.csv"))
  cat("\n\u2713 Tabla Fisher guardada\n\n")

  # One forest plot per sheet and grouping variable
  sheet_labels <- setNames(
    sapply(sheet_list, `[[`, "label"),
    sapply(sheet_list, `[[`, "name")
  )
  grp_labels <- c("depth" = "Profundidad", "station" = "Estación", "genus" = "Género")

  fig_counter <- 10L  # Fig10a, 10b, ... continua desde las anteriores
  fig_letters <- c(letters, paste0("a", letters))  # a-z, aa-az
  fig_idx     <- 1L

  for (sh in sheet_list) {
    for (gv in grp_list) {
      sub_df <- fisher_all %>% filter(sheet == sh$name, group_var == gv$name)
      if (nrow(sub_df) == 0) { fig_idx <- fig_idx + 1L; next }
      outfile <- file.path(OUTPUT_DIR,
                           sprintf("Fig%d%s_fisher_%s_%s.png",
                                   fig_counter, fig_letters[fig_idx],
                                   sh$name, gv$name))
      plot_forest(
        sub_df,
        title    = paste0("Test de Fisher – ", sheet_labels[sh$name],
                          " por ", grp_labels[gv$name]),
        outfile  = outfile,
        facet_col = "trait"
      )
      fig_idx <- fig_idx + 1L
    }
  }
  cat("\u2713 Sección 2 completada\n\n")
} else {
  cat("\u26a0 Sin resultados del test de Fisher\n\n")
}

# ---- ---- section 3: logistic model per trait ----
cat("=== SECCIÓN 3: MODELOS LOGÍSTICOS ===\n\n")

meta_long <- tibble(
  genome_id = metadata$NameWithoutExtension,
  depth     = as.character(depth_grp[metadata$NameWithoutExtension]),
  station   = as.character(station_grp[metadata$NameWithoutExtension]),
  genus     = as.character(genus_grp[metadata$NameWithoutExtension]),
  order     = as.character(order_grp[metadata$NameWithoutExtension])
) %>%
  filter(!is.na(depth), !is.na(station),
         genus != "Desconocido", order != "Desconocido")

run_logistic_models <- function(trait_df, sheet_name) {
  genome_col <- names(trait_df)[1]
  genomes    <- trait_df[[genome_col]]
  mat_raw    <- trait_df %>% select(-all_of(genome_col)) %>% as.matrix()
  mat_bin    <- matrix(as.numeric(mat_raw > 0), nrow = nrow(mat_raw),
                       dimnames = list(genomes, colnames(mat_raw)))

  results <- vector("list", ncol(mat_bin))
  k <- 0L
  for (trait in colnames(mat_bin)) {
    trait_vec <- mat_bin[genomes, trait]
    df_model  <- meta_long %>%
      filter(genome_id %in% genomes) %>%
      mutate(trait_bin = trait_vec[genome_id])

    if (nrow(df_model) < 5 || length(unique(df_model$trait_bin)) < 2) next

    depth_ref <- if ("ZMOR" %in% df_model$depth) "ZMOR" else sort(unique(df_model$depth))[1]
    df_model$depth   <- relevel(factor(df_model$depth),   ref = depth_ref)
    df_model$station <- relevel(factor(df_model$station), ref = sort(unique(df_model$station))[1])
    df_model$genus   <- relevel(factor(df_model$genus),   ref = sort(unique(df_model$genus))[1])
    df_model$order   <- relevel(factor(df_model$order),   ref = sort(unique(df_model$order))[1])

    fit <- tryCatch(
      glm(trait_bin ~ depth + station + genus + order, family = binomial, data = df_model),
      error = function(e) NULL
    )
    if (is.null(fit)) next

    tidy_res <- tryCatch(
      broom::tidy(fit, exponentiate = TRUE, conf.int = TRUE) %>%
        filter(term != "(Intercept)") %>%
        mutate(sheet = sheet_name, trait = trait),
      error = function(e) NULL
    )
    if (is.null(tidy_res)) next
    k <- k + 1L
    results[[k]] <- tidy_res
  }
  if (k == 0L) return(NULL)
  bind_rows(results[seq_len(k)])
}

logit_parts <- list()
for (sh in sheet_list) {
  cat("  Logístico:", sh$label, "...")
  res <- run_logistic_models(sh$df, sh$name)
  if (!is.null(res)) logit_parts[[length(logit_parts) + 1]] <- res
  cat(" OK\n")
}

if (length(logit_parts) > 0) {
  logit_all <- bind_rows(logit_parts) %>%
    mutate(
      term_family = dplyr::case_when(
        startsWith(term, "depth")   ~ "depth",
        startsWith(term, "station") ~ "station",
        startsWith(term, "genus")   ~ "genus",
        startsWith(term, "order")   ~ "order",
        TRUE                        ~ "other"
      ),
      term_value = dplyr::case_when(
        term_family == "depth"   ~ sub("^depth",   "", term),
        term_family == "station" ~ sub("^station", "", term),
        term_family == "genus"   ~ sub("^genus",   "", term),
        term_family == "order"   ~ sub("^order",   "", term),
        TRUE                     ~ term
      ),
      term_value_display = dplyr::case_when(
        term_family == "genus" ~ taxon_italic(term_value, level = "Genus"),
        term_family == "order" ~ taxon_italic(term_value, level = "Order"),
        TRUE                   ~ term_value
      )
    ) %>%
    group_by(sheet, term_family) %>%
    mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
    ungroup() %>%
    mutate(
      OR          = estimate,
      CI_low      = conf.low,
      CI_high     = conf.high,
      significant = p_adj < 0.05,
      label       = paste0(term_family, ": ", term_value_display)
    )
  write_csv(logit_all, file.path(OUTPUT_DIR, "logistic_results.csv"))
  cat("\n\u2713 Tabla logística guardada\n\n")

  cat("Generando forest plots logísticos...\n")
  fig_idx2 <- 1L
  for (sh in sheet_list) {
    sub_df <- logit_all %>% filter(sheet == sh$name)
    if (nrow(sub_df) == 0) { fig_idx2 <- fig_idx2 + 1L; next }
    plot_forest(
      sub_df,
      title    = paste0("Modelo logístico – ", sheet_labels[sh$name],
                        "\nlogit(P) = \u03b20 + \u03b21·profundidad + \u03b22·estación + \u03b23·género + \u03b24·orden"),
      outfile  = file.path(OUTPUT_DIR,
                           sprintf("Fig11%s_logistic_%s.png",
                                   fig_letters[fig_idx2], sh$name)),
      facet_col = "trait"
    )
    fig_idx2 <- fig_idx2 + 1L
  }
  cat("\u2713 Sección 3 completada\n\n")
} else {
  cat("\u26a0 Sin resultados del modelo logístico\n\n")
}

cat("=== treatment_comparison.R completado ===\n")
