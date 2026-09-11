#!/usr/bin/env Rscript

# Figures for the Microbial Dark Matter symposium poster. Unlike everything else
# in this project these are labelled in English, since the poster is.
#
# Fifteen figures, each written in two formats, SVG and PNG at 600 dpi, all on a
# transparent background because the poster lays them over its own colour. The
# files keep a `_transparent` suffix from when a white-background variant was
# also emitted; that variant is disabled inside ggsave_both() and
# draw_autosize() rather than deleted, so it can be brought back.
#
# The SVG goes through the cairo device, grDevices::svg(), not svglite. svglite
# leaves background rectangles without a `fill` attribute and declares the
# default in a CSS block, and importers that drop the stylesheet, Canva among
# them, then render those rectangles black. Section 1b covers this. Formats are
# listed in OUT_FORMATS and the raster resolution in PNG_DPI.
#
#     1.  Fig_traits_combined_sidebyside_EN
#         Fig6a and 6b side by side in one horizontal panel
#     2.  Fig_traits_combined_vertical_EN
#         the same transposed and stacked, without dendrograms
#     3.  Fig_traits_combined_vertical_dend_EN
#         with dendrograms over both traits and orders
#     4.  Fig_traits_combined_vertical_dend_nobar_EN
#         the same again, without the MAG-count barplot
#     5.  Fig_annotation_by_tool_EN_axislines
#         plot_03 from the Spanish figures, laid out horizontally
#     6.  Fig_annotation_rate_per_order_EN_axislines
#         plot_03c from the Spanish figures
#     7.  Panel_Meta_Tax_Order_EN
#         the panel of the same name from CNPS.cycle.R
#     8-9.  Fig_MAGs_per_depth_by_order[_light]_EN_axislines
#         MAGs recovered per depth, stacked by order
#    10-13. Fig_ANI_unassigned_groups_withref_A_{diamond,square}[_purple]_EN
#         within-group ANI beside ANI against the GTDB reference on the same
#         row, in two shapes and two colours
#    14-15. Fig_pangenome_classes_stacked[_light]_EN_axislines
#         Fig_panel_pan_classes as a single 100% stacked bar per species
#
# The two stacked-bar figures come in two palettes: the darkened one from
# darken(), with white numerals inside the bars, and the `_light` one with the
# original colours and black numerals.
#
# The four bar figures (5, 6, 8-9 and 14-15) go through ggsave_frames(), which
# currently emits only the `_axislines` finish: a black left y axis and bottom x
# axis, no panel frame and no grid. The full-frame and `_noframe` finishes are
# commented out inside that function.
#
# Figures 6a and 6b as standalone panels are commented out in sections 6 and 7.
#
# The code is ported from four scripts in the project:
#
#   DRAM2_output_analysis.R      the trait heatmaps, make_trait_order_heatmap(),
#                                without the DeepL layer since the trait names
#                                are already English, plus the annotation bars
#                                from its sections 2 and 5
#   CNPS.cycle.R                 the combined metadata and taxonomy panel from
#                                its section G, make_combined_meta_tax_g("Order")
#   ANI_intraspecies.R           loading the skani edge list, rescuing pairs from
#                                novel_vs_all, and Fig1b. The literal translation
#                                of Fig1b, the two-track variant and the triangle
#                                variant are commented out in section 15
#   27_pangenome_classification.R  pangenome class counts, read from
#                                class_counts_summary.tsv rather than recomputed

cat("Inicio:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "") {
  user_lib <- file.path(
    Sys.getenv("HOME"), "R", "x86_64-pc-linux-gnu-library",
    paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = ".")
  )
}
if (!dir.exists(user_lib)) dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(user_lib, .libPaths()))

suppressPackageStartupMessages({
  for (pkg in c("dplyr", "tidyr", "tibble", "readr", "readxl", "ggplot2", "scales",
                "circlize", "grid", "RColorBrewer", "gridtext", "ggtext", "stringr")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, repos = "https://cloud.r-project.org/", quiet = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
  if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    if (!requireNamespace("BiocManager", quietly = TRUE)) {
      install.packages("BiocManager", repos = "https://cloud.r-project.org/")
    }
    BiocManager::install("ComplexHeatmap", ask = FALSE, update = FALSE)
  }
  library(ComplexHeatmap)
})

# =============================================================================
# 0. RUTAS
# =============================================================================
# RStudio on Windows reaches the project over a UNC path into WSL, while WSL
# itself sees a POSIX path. The UNC form is tried first and falls back, so the
# script runs unchanged in either.
BASE_DIR <- resolve_project_dir()

SCRIPTS_R_DIR <- file.path(BASE_DIR, "scripts", "R")
METADATA_FILE <- file.path(BASE_DIR, "MAGs_metadata.tsv")
TRAITS_FILE   <- file.path(BASE_DIR, "results", "08_MAG_distill", "SUMMARIZE", "traits.xlsx")
ANNOT_FILE    <- file.path(BASE_DIR, "results", "08_MAG_distill", "ANNOTATE", "raw-annotations.tsv")
# CNPS panel, section 13: per-cycle abundance and GTDB taxonomy. As in section G
# of CNPS.cycle.R, this reads the CSVs already on disk rather than the
# objetos en memoria de aquel script.
CNPS_DIR      <- file.path(BASE_DIR, "results", "15_CNPS_cycle")
GTDB_FILE     <- file.path(BASE_DIR, "results", "03_taxonomy_assigned", "classify", "tax_.ar53.summary.tsv")
# Inputs for the newer figures, sections 14 to 16.
ANI_SPARSE_FILE <- file.path(BASE_DIR, "results", "16_sketching_ani", "all_mags",
                             "ani_triangle_sparse.tsv")
ANI_NOVEL_FILE  <- file.path(BASE_DIR, "results", "16_sketching_ani", "novel_mags",
                             "novel10_vs_all54.tsv")
PAN_CLASS_FILE  <- file.path(BASE_DIR, "results", "SO_03", "classification",
                             "class_counts_summary.tsv")
# The poster material gets its own directory under results/, because
# SO_01/funct_annot/ is wiped at the start of every DRAM2_output_analysis.R run.
OUTPUT_DIR    <- file.path(BASE_DIR, "results", "MDMS_figs")

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("Directorio base  :", BASE_DIR, "\n")
cat("Directorio salida:", OUTPUT_DIR, "\n\n")

stopifnot(file.exists(METADATA_FILE), file.exists(TRAITS_FILE), file.exists(ANNOT_FILE),
          dir.exists(CNPS_DIR), file.exists(GTDB_FILE),
          file.exists(ANI_SPARSE_FILE), file.exists(PAN_CLASS_FILE))

# =============================================================================
# ---- 1. typography, matching the rest of the project ----
# =============================================================================
source(file.path(SCRIPTS_R_DIR, "_pub_style.R"))  # define font_family + showtext

# Injects the font family into any grid::gpar(), so ComplexHeatmap renders in
# the same face as the ggplots.
gpar_f <- function(...) {
  args <- list(...)
  if (is.null(args$fontfamily)) args$fontfamily <- font_family
  do.call(grid::gpar, args)
}

# _pub_style.R exports a theme_publication() too, but it scales differently from
# the annotation script: base_size 16 and element_markdown titles. It is
# overridden here with the DRAM2_output_analysis.R version so figures 4 and 5
# come out identical to their Spanish counterparts.
theme_publication <- function(base_size = 11) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.x = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "gray", linetype = "dotted"),
      axis.text          = ggplot2::element_text(color = "black"),
      axis.title         = ggplot2::element_text(face = "bold"),
      legend.position    = "right",
      legend.background  = ggplot2::element_rect(fill = "white", color = NA),
      legend.title       = ggplot2::element_text(face = "bold", size = 12),
      legend.text        = ggplot2::element_text(size = 11),
      strip.background   = ggplot2::element_rect(fill = "gray90", color = "black"),
      strip.text         = ggplot2::element_text(face = "bold", size = 13),
      plot.title         = ggplot2::element_text(face = "bold", hjust = 0.5),
      text               = ggplot2::element_text(family = font_family)
    )
}

# =============================================================================
# 1b. FORMATOS DE SALIDA — SVG (vectorial) + PNG a 600 dpi
# =============================================================================
# Every figure is written twice, same name and different extension:
#   .svg  vector, for laying out and rescaling without loss
#   .png  raster at 600 dpi, as a fallback wherever SVG is handled badly, which
#         covers Canva, Word, thumbnail previews and email
# Drop an entry from OUT_FORMATS to emit only one.
OUT_FORMATS <- c(".svg", ".png")
PNG_DPI     <- 600

# Text width is measured on the first format in the list;
# medir y pintar en el mismo tipo de dispositivo evita desajustes de ~0.5 % en
# the glyph metrics differ between the raster and vector devices.
MEASURE_EXT <- OUT_FORMATS[1]

# showtext rasterises glyphs at whatever resolution it is told, so its dpi has
# to match the raster device exactly. Leaving it at 300 while writing a 600 dpi
# PNG renders the text at half size and throws the whole layout out; measured,
# 12% of the pixels change. _pub_style.R sets 300, so it is raised here to the
# real output resolution. On the vector devices
# ajuste es indiferente.
showtext_opts(dpi = PNG_DPI)

# --- Dispositivo SVG: cairo, no svglite --------------------------------------
# svglite declares its defaults (fill: none, stroke: #000000) in a CSS <style>
# block inside <defs>, and leaves the background rectangles of the figure, the
# panel and the legend keys with no fill attribute at all. A viewer that applies
# CSS, meaning any browser, paints them transparent. One that
# ignore o descarte la hoja de estilo —el importador de Canva -- cae al valor
# falls back to the SVG initial value, which is fill: black, and the background
# rectangle, which covers the
# lienzo entero, se traga la figura en negro.
#
# The grDevices svg() device, cairo, writes the full style as attributes on each
# element and emits no CSS at all, so the file looks the same in every importer.
# That is why it is preferred here, at the cost of about 20% more file size.
SVG_CAIRO <- isTRUE(capabilities("cairo"))
if (!SVG_CAIRO) {
  message("[MDMS_figs] Sin soporte cairo: se usará svglite y se inyectarán los ",
          "estilos por omisión en el <g> raíz (ver fix_svglite_css).")
}

# ggsave() picks its device from the extension and prefers svglite for .svg
# whenever it is installed, so this wrapper has to be passed in explicitly. The
# ... swallows the arguments ggsave() adds that svg() will not take, dpi among
# them.
svg_device <- function(filename, width, height, bg = "transparent", ...) {
  if (SVG_CAIRO) {
    grDevices::svg(filename = filename, width = width, height = height, bg = bg)
  } else {
    svglite::svglite(filename, width = width, height = height, bg = bg)
  }
}

# Safety net for the svglite path: lift the defaults out of the CSS block and
# onto the root <g> as presentation attributes. Attributes inherit down, and any
# element carrying its own style still overrides what it inherits, so nothing
# changes in a viewer that applies CSS. A viewer that ignores CSS no longer
# paints the unfilled backgrounds black.
fix_svglite_css <- function(f) {
  if (SVG_CAIRO || !file.exists(f)) return(invisible(FALSE))
  txt <- readLines(f, warn = FALSE)
  i <- grep("<g class='svglite'>", txt, fixed = TRUE)[1]
  if (is.na(i)) return(invisible(FALSE))
  txt[i] <- paste0("<g class='svglite' fill='none' stroke='#000000' ",
                   "stroke-linecap='round' stroke-linejoin='round' ",
                   "stroke-miterlimit='10.00'>")
  writeLines(txt, f)
  invisible(TRUE)
}

# --- Nombres de archivo -------------------------------------------------------
# Normaliza un nombre al formato pedido e inserta un sufijo antes de la
# extension: out_file("Fig_x_EN.svg", "_axislines", ".png") gives
# "Fig_x_EN_axislines.png". Any input extension is accepted, so the literals
# scattered through the script can stay as .png or .svg without mattering.
out_file <- function(outfile, suffix = "", ext = MEASURE_EXT) {
  paste0(sub("\\.(png|svg|pdf)$", "", outfile), suffix, ext)
}

# Opens the device matching the file extension, sized in inches. Used by
# draw_autosize(), since ComplexHeatmap draws onto an already-open device rather
# than through ggsave, and by the text measurement above.
open_device <- function(f, w_in, h_in, bg_col = "transparent") {
  if (grepl("\\.svg$", f)) {
    svg_device(f, w_in, h_in, bg_col)
  } else {
    png(f, width = w_in, height = h_in, units = "in", res = PNG_DPI, bg = bg_col)
  }
}

# =============================================================================
# 2. METADATOS TAXONÓMICOS
# =============================================================================
metadata <- read.delim(METADATA_FILE, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

.extract_rank <- function(tx, prefix) {
  if (is.na(tx)) return(NA_character_)
  m <- regmatches(tx, regexpr(paste0(prefix, "[^;]+"), tx))
  if (length(m) == 0 || m == prefix) NA_character_ else sub(paste0("^", prefix), "", m)
}

# genome_id to order, from GTDB FullTaxonomy
mag_to_order <- setNames(
  sapply(metadata$FullTaxonomy, .extract_rank, prefix = "o__"),
  metadata$NameWithoutExtension
)

# Order to class, which gives a code-like order name some context.
order_to_class <- local({
  tbl <- data.frame(
    order = sapply(metadata$FullTaxonomy, .extract_rank, prefix = "o__"),
    class = sapply(metadata$FullTaxonomy, .extract_rank, prefix = "c__"),
    stringsAsFactors = FALSE
  )
  tbl <- tbl[!is.na(tbl$order) & !is.na(tbl$class) &
               nchar(trimws(tbl$order)) > 0 & nchar(trimws(tbl$class)) > 0, ]
  tbl <- tbl[!duplicated(tbl$order), ]
  setNames(tbl$class, tbl$order)
})

# Puts the class in parentheses after an order name that is really a code: first
# word all caps, MGIII or DPANN, or carrying digits, UBA57.
fmt_order_label <- function(ord_name) {
  if (is.na(ord_name)) return(ord_name)
  first_word <- strsplit(ord_name, " ")[[1]][1]
  is_acr <- grepl("[0-9]", first_word) || grepl("^[A-Z]+$", first_word)
  if (!is_acr) return(ord_name)
  cls <- if (ord_name %in% names(order_to_class)) order_to_class[[ord_name]] else NA_character_
  if (!is.na(cls) && nchar(trimws(cls)) > 0 && cls != ord_name) {
    return(paste0(ord_name, " (c. ", cls, ")"))
  }
  ord_name
}

# Agrega un data.frame (genome, rasgo1, rasgo2, ...) a nivel de Orden.
# True and False become 1 and 0 and average into a per-order proportion of MAGs.
aggregate_to_order <- function(df) {
  genome_col <- names(df)[1]
  orders <- mag_to_order[df[[genome_col]]]
  feat_df <- df %>%
    select(-all_of(genome_col)) %>%
    mutate(across(everything(), ~ {
      v <- as.character(.)
      as.numeric(ifelse(v %in% c("True", "TRUE"), 1,
                 ifelse(v %in% c("False", "FALSE"), 0, v)))
    }))
  bind_cols(tibble(order = orders), feat_df) %>%
    filter(!is.na(order)) %>%
    group_by(order) %>%
    summarise(across(everything(), ~ mean(., na.rm = TRUE)), .groups = "drop")
}

# =============================================================================
# 3. ETIQUETAS DE RASGO EN INGLÉS
# =============================================================================
# DRAM2 writes trait names with the category run onto the front with no
# separator: "Simple Sugars Glucose", "Short-Chain Fatty Acids (SCFA)
# Propionate". On an x axis those read badly. prettify_trait() only reformats:
# it separates category from substrate with a colon and normalises the
# hydrogenase notation. No trait is renamed and none are merged.
# Set PRETTIFY_LABELS to FALSE for the raw DRAM2 names.
PRETTIFY_LABELS <- TRUE

# Category prefix to short label, for the carbon sheet
CATEGORY_PREFIXES <- c(
  "Short-Chain Fatty Acids (SCFA)" = "SCFA",
  "Complex Polysaccharide"         = "Polysaccharide",
  "Microbial Compounds"            = "Microbial",
  "Aromatics/Phenolics"            = "Aromatic",
  "Organic Acids"                  = "Organic acid",
  "Simple Sugars"                  = "Sugar",
  "Amino Acids"                    = "Amino acid",
  "Lipids"                         = "Lipid"
)

# Reescrituras literales (sheet "energy"); nombre DRAM2 to etiqueta de figura
TRAIT_RENAMES <- c(
  "Aerobic Ammonia Oxidizer"                   = "Aerobic ammonia oxidizer (AOA)",
  "Aerobic CO Oxidizer"                        = "Aerobic CO oxidizer",
  "Anaerobic CO Oxidizer"                      = "Anaerobic CO oxidizer",
  "Anaerobic Methane Oxidizer ANME"            = "Anaerobic methane oxidizer (ANME)",
  "Assimilatory Sulfate Reducer"               = "Assimilatory sulfate reducer",
  "Dissimilatory Sulfate Reducer Strict"       = "Dissimilatory sulfate reducer (strict)",
  "DMSO/TMAO Reducer DmsABC"                   = "DMSO/TMAO reducer (DmsABC)",
  "Sulfide Oxidizer Sqr Fcc"                   = "Sulfide oxidizer (Sqr/Fcc)",
  "Sulfur Oxidizer SOX Complete"               = "Sulfur oxidizer (complete SOX)",
  "Thiosulfate Reducer PhsABC"                 = "Thiosulfate reducer (PhsABC)",
  "Urease Nitrogen Mineralization"             = "Urease (N mineralization)",
  "Hydrogenase Group 1 NiFe"                   = "Hydrogenase Grp1 [NiFe]",
  "Hydrogenase Group 3a NiFe F420"             = "Hydrogenase Grp3a [NiFe] F420",
  "Hydrogenase Group 3c NiFe Mvh"              = "Hydrogenase Grp3c [NiFe] Mvh",
  "Hydrogenase Group 3d NiFe Hox"              = "Hydrogenase Grp3d [NiFe] Hox",
  "Hydrogenase Group 4 NiFe Energy Converting" = "Hydrogenase Grp4 [NiFe] (Ech/Mbh)",
  "Hydrogenase FeFe Hnd"                       = "Hydrogenase [FeFe] Hnd",
  "Hydrogenase FeFe HydA"                      = "Hydrogenase [FeFe] HydA",
  "Any Hydrogenase"                            = "Any hydrogenase",
  "Hydrogen Evolver"                           = "Hydrogen evolver",
  "Iron Oxidation"                             = "Iron oxidation",
  "Iron Reduction"                             = "Iron reduction",
  "Probable Iron Reduction"                    = "Probable iron reduction",
  "Fumarate Reduction"                         = "Fumarate reduction",
  "Tetrathionate Reduction"                    = "Tetrathionate reduction",
  "annamox"                                    = "Anammox",
  "N reducer"                                  = "N reducer",
  "S reducer"                                  = "S reducer",
  "S oxidizer"                                 = "S oxidizer"
)

prettify_trait <- function(x) {
  if (!PRETTIFY_LABELS) return(x)
  vapply(x, function(nm) {
    if (is.na(nm)) return(nm)
    if (nm %in% names(TRAIT_RENAMES)) return(unname(TRAIT_RENAMES[[nm]]))
    for (pfx in names(CATEGORY_PREFIXES)) {
      if (startsWith(nm, paste0(pfx, " "))) {
        rest <- substring(nm, nchar(pfx) + 2)
        return(paste0(CATEGORY_PREFIXES[[pfx]], ": ", rest))
      }
    }
    nm
  }, character(1), USE.NAMES = FALSE)
}

# =============================================================================
# 4. CONSTRUCCIÓN DE LA MATRIZ Orden × rasgo
# =============================================================================
# Returns list(mat = proportions, order by trait; n = MAGs per order). A trait
# absent from every MAG, meaning its column sums to zero, is dropped, as in
# make_trait_order_heatmap() in the original script.
build_order_matrix <- function(trait_df, label) {
  genome_col <- names(trait_df)[1]
  agg <- aggregate_to_order(trait_df)
  if (nrow(agg) < 2) stop("Insuficientes Órdenes para ", label)

  n_per_order <- table(mag_to_order[trait_df[[genome_col]]])
  mat <- agg %>% column_to_rownames("order") %>% as.matrix()
  mat <- mat[, colSums(mat, na.rm = TRUE) > 0, drop = FALSE]
  if (ncol(mat) == 0) stop("Sin rasgos con datos para ", label)

  colnames(mat) <- prettify_trait(colnames(mat))
  cat("  ", label, ": ", nrow(mat), " Órdenes × ", ncol(mat),
      " rasgos presentes\n", sep = "")
  list(mat = mat, n = as.integer(n_per_order[rownames(mat)]))
}

# Continuous 0 to 1: the proportion of that order's MAGs carrying the trait
COL_FUN <- colorRamp2(c(0, 0.5, 1), c("white", "#74C476", "#006D2C"))
LEGEND_TITLE <- "Proportion of\nMAGs with trait"

heatmap_legend <- list(
  title     = LEGEND_TITLE,
  at        = c(0, 0.25, 0.5, 0.75, 1),
  labels    = c("0%", "25%", "50%", "75%", "100%"),
  title_gp  = gpar_f(fontsize = 10, fontface = "bold"),
  labels_gp = gpar_f(fontsize = 9),
  border    = "black"
)

# Side barplot of MAGs per order
n_mags_annotation <- function(n_vals) {
  rowAnnotation(
    `No. MAGs` = anno_barplot(n_vals,
      gp = gpar_f(fill = "#4A90D9", col = "black"),
      bar_width = 0.7, width = unit(2, "cm")),
    annotation_name_gp = gpar_f(fontsize = 9, fontface = "bold")
  )
}

# Real width in mm of the longest string in a vector, measured in the final font
# and size on the same device the figure is drawn on
# (showtext activo). Medir en un dispositivo distinto al de salida desajusta el
# the reserved label width: glyph metrics differ by about 0.5% between
# el raster y el vectorial.
text_width_mm <- function(labels, fontsize, fontface = "plain") {
  tmp <- tempfile(fileext = MEASURE_EXT)
  open_device(tmp, 8, 8, "white")
  w <- suppressWarnings(max(vapply(labels, function(s) {
    convertWidth(grobWidth(textGrob(s, gp = gpar_f(fontsize = fontsize, fontface = fontface))),
                 "mm", valueOnly = TRUE)
  }, numeric(1))))
  dev.off(); unlink(tmp)
  w
}

# Height in mm to reserve for column names rotated 45 degrees.
# ComplexHeatmap reserves the full text width even when the text is rotated,
# which leaves about 15 mm of empty canvas below the labels.
# column_names_max_height caps it at the real vertical projection, width times
# times sin 45, plus a little slack.
col_names_height_mm <- function(cn, fontsize = 9) {
  text_width_mm(cn, fontsize) * sin(pi / 4) + 2
}

# -----------------------------------------------------------------------------
# Transparent backgrounds only.
# This used to write each figure twice, X_EN.svg on white and
# X_EN_transparent.svg with an alpha channel. The poster always sits the figures
# on its own coloured background, so the white variant went unused and its
# lineas quedan comentadas -- no borradas -- en `ggsave_both()` y
# draw_autosize(), in case a white-background version is ever needed for print.
# The _transparent suffix stays in the filenames so the images already placed in
# the poster layout do not break.
# No geom uses alpha. Colours are premixed through blend_over_white() and drawn
# opaque, because against a transparent background the
# alpha deja pasar el color de atras y altera la paleta.
# -----------------------------------------------------------------------------
transparent_path <- function(outfile, ext = MEASURE_EXT) {
  out_file(outfile, "_transparent", ext)
}

# theme_bw() paints both panel and canvas white. This makes them transparent
# without touching the black panel border, which is wanted.
theme_transparent_bg <- ggplot2::theme(
  plot.background       = ggplot2::element_rect(fill = "transparent", color = NA),
  panel.background      = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.background     = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.box.background = ggplot2::element_rect(fill = "transparent", color = NA),
  legend.key            = ggplot2::element_rect(fill = "transparent", color = NA)
)

# Flattens a semitransparent colour over white and returns the solid equivalent.
# With alpha set on the geom, whatever sits behind a transparent background
# shows through the bar and the palette shifts. Premixing and drawing at
# alpha = 1 gives a result identical to the old one over white, and stable over
# any other background.
blend_over_white <- function(cols, alpha) {
  m <- grDevices::col2rgb(cols)
  grDevices::rgb(t(round(alpha * m + (1 - alpha) * 255)), maxColorValue = 255)
}

# Saves a ggplot on a transparent background. outfile is still the base name
# with no suffix; this writes only X_transparent.svg.
ggsave_both <- function(outfile, plot, width, height, dpi = PNG_DPI) {
  # [DESACTIVADO] variante de fondo blanco:
  # ggsave(outfile, plot, width = width, height = height, dpi = dpi, bg = "white")
  p_out <- plot + theme_transparent_bg
  for (ext in OUT_FORMATS) {
    f <- transparent_path(outfile, ext)
    if (identical(ext, ".svg")) {
      # Explicit device: ggsave() would pick svglite from the extension, see 1b.
      ggsave(f, p_out, width = width, height = height, bg = "transparent",
             device = svg_device)
      fix_svglite_css(f)
    } else {
      ggsave(f, p_out, width = width, height = height, dpi = dpi,
             bg = "transparent")
    }
    cat("  [OK] Guardado:", basename(f),
        sprintf("(%.1f × %.1f in)\n", width, height))
  }
}

# -----------------------------------------------------------------------------
# Shared utilities.
# -----------------------------------------------------------------------------

# Thousands separator: a space rather than a comma. The comma is the
# anglophone convention, but on a poster carrying decimals too, the ISO 31-0 and
# SI space reads without ambiguity.
comma_sp <- function(x) scales::comma(x, big.mark = " ")

# Darkens a colour by scaling its RGB channels. Used where text sits on top of
# va texto blanco: la paleta Dark2 original tiene luminancia alta en el amarillo
# the lighter greens, where white text would otherwise disappear.
darken <- function(cols, factor = 0.75) {
  m <- grDevices::col2rgb(cols)
  # Channels passed separately: rgb() given a matrix expects three rows, and
  # handing it the transpose returns a single colour instead of one per column.
  out <- grDevices::rgb(
    pmin(255, pmax(0, round(m[1, ] * factor))),
    pmin(255, pmax(0, round(m[2, ] * factor))),
    pmin(255, pmax(0, round(m[3, ] * factor))),
    maxColorValue = 255
  )
  names(out) <- names(cols)
  out
}

# --- Variantes de marco -------------------------------------------------------
# Each bar figure can be emitted in three finishes, so the poster can pick how
# much frame it needs against whatever background it lands on:
#   (base)      full panel frame with a dotted grid, from theme_publication
#   _noframe    no frame, no grid
#   _axislines  left y axis and bottom x axis only, no grid
#
# ggplot2's axis.line draws only the left and bottom axes by default, so it
# works as is even under coord_flip().
#
# panel.grid.major = element_blank() is not enough on its own. The base figures
# set panel.grid.major.x and .y explicitly, and in ggplot2 the more specific
# setting wins over the inherited one, so all four have to be blanked.
blank_grid <- ggplot2::theme(
  panel.grid.major   = ggplot2::element_blank(),
  panel.grid.major.x = ggplot2::element_blank(),
  panel.grid.major.y = ggplot2::element_blank(),
  panel.grid.minor   = ggplot2::element_blank()
)

theme_no_frame <- blank_grid + ggplot2::theme(
  panel.border = ggplot2::element_blank(),
  axis.line    = ggplot2::element_blank(),
  # With no frame and no axis line, a lone tick floats in empty space.
  axis.ticks   = ggplot2::element_blank()
)

# Axes and ticks in full black. The theme_bw() default of grey20 washes out
# against the poster's coloured background.
theme_axis_lines <- blank_grid + ggplot2::theme(
  panel.border = ggplot2::element_blank(),
  axis.line    = ggplot2::element_line(color = "black", linewidth = 0.6),
  axis.ticks   = ggplot2::element_line(color = "black", linewidth = 0.6)
)

# Writes the frame variants of one figure. The poster settled on _axislines
# alone; the other two are commented out rather than deleted.
ggsave_frames <- function(outfile, plot, width, height, dpi = PNG_DPI) {
  # Disabled: full panel frame with a dotted grid
  # ggsave_both(outfile, plot, width, height, dpi)
  # Disabled: no frame, no axes, no grid
  # ggsave_both(out_file(outfile, "_noframe"),
  #             plot + theme_no_frame, width, height, dpi)
  ggsave_both(out_file(outfile, "_axislines"),
              plot + theme_axis_lines, width, height, dpi)
}

# -----------------------------------------------------------------------------
# Exact canvas sizing.
#
# Every component of these heatmaps has a fixed size: the body in mm, the
# barplot at 2 cm, the dendrograms. That lets ComplexHeatmap report the real
# size of the drawn object. It is measured on a null device first, then the
# canvas opens at exactly that size, which removes the empty bottom margin the
# original figures carry from having their height set by eye.
#
# The "font family not found in PostScript font database" warnings from the
# measuring device are harmless; showtext resolves the font when it paints.
# -----------------------------------------------------------------------------
measure_ht <- function(ht, ...) {
  suppressWarnings({
    pdf(NULL, width = 40, height = 40)
    d <- draw(ht, ...)
    wh <- c(convertWidth(ComplexHeatmap:::width(d), "in", valueOnly = TRUE),
            convertHeight(ComplexHeatmap:::height(d), "in", valueOnly = TRUE))
    dev.off()
  })
  wh
}

draw_autosize <- function(ht, outfile, margin_in = 0.08, ...) {
  wh <- measure_ht(ht, ...)
  w_in <- wh[1] + 2 * margin_in
  h_in <- wh[2] + 2 * margin_in
  # The device bg is not enough on its own: ComplexHeatmap's draw() paints its
  # own background rectangle over the top, white by default, which covers the
  # alpha channel. That has to be cancelled too.
  # Disabled: the white-background variant, which would need only
  # `for (bg_col in c("white", "transparent"))` y elegir el nombre de archivo
  # to follow bg_col again, as the earlier version did.
  bg_col <- "transparent"
  for (ext in OUT_FORMATS) {
    f <- transparent_path(outfile, ext)
    open_device(f, w_in, h_in, bg_col)
    draw(ht, ..., background = bg_col)
    dev.off()
    fix_svglite_css(f)
    cat("  [OK] Guardado:", basename(f),
        sprintf("(%.2f × %.2f in)\n", w_in, h_in))
  }
  invisible(wh)
}

# Size and colour of the "n = X" line under each order name.
ROW_LABEL_FONTSIZE <- 10
N_LABEL_FONTSIZE   <- ROW_LABEL_FONTSIZE - 2   # dos puntos menos que el Orden
N_LABEL_COLOR      <- "#4D4D4D"                # gris oscuro

# -----------------------------------------------------------------------------
# The base order-by-trait heatmap. panel_title NULL means no title, which is
# what the standalone figures use.
#
# row_n_vals, given the MAG count per order in the row order of mat, adds an
# "n = X" line under each row label. It goes through gridtext::gt_render()
# because a single text grob cannot carry two sizes and two colours;
# ComplexHeatmap reorders the labels along with the rows.
# -----------------------------------------------------------------------------
trait_heatmap <- function(mat, name, col_width_mm = 10, row_height_mm = 10,
                          cluster_rows = TRUE, show_row_dend = TRUE,
                          show_row_names = TRUE, row_n_vals = NULL,
                          right_annotation = NULL, show_legend = TRUE,
                          panel_title = NULL) {
  fmt_rn <- vapply(rownames(mat), fmt_order_label, character(1))
  if (is.null(row_n_vals)) {
    row_lab   <- fmt_rn
    row_lab_w <- text_width_mm(fmt_rn, ROW_LABEL_FONTSIZE) + 2
  } else {
    row_lab <- gt_render(sprintf(
      "%s<br><span style='font-size:%dpt;color:%s'>n = %d</span>",
      fmt_rn, N_LABEL_FONTSIZE, N_LABEL_COLOR, as.integer(row_n_vals)
    ))
    # gridtext needs a little more slack than a plain textGrob
    row_lab_w <- text_width_mm(fmt_rn, ROW_LABEL_FONTSIZE) + 5
  }
  Heatmap(mat,
    name = name,
    col = COL_FUN,
    cluster_rows = cluster_rows,
    show_row_dend = show_row_dend,
    row_dend_side = "left",
    row_dend_width = unit(12, "mm"),
    cluster_columns = TRUE,
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "ward.D2",
    show_row_names = show_row_names,
    # Ranks above genus stay upright, per CDC style.
    row_labels = row_lab,
    row_names_side = "right",
    row_names_gp = gpar_f(fontsize = ROW_LABEL_FONTSIZE),
    row_names_max_width = unit(row_lab_w, "mm"),
    width  = unit(ncol(mat) * col_width_mm, "mm"),
    height = unit(nrow(mat) * row_height_mm, "mm"),
    column_names_gp = gpar_f(fontsize = 9),
    column_names_rot = 45,
    column_names_max_height = unit(col_names_height_mm(colnames(mat)), "mm"),
    column_title = panel_title,
    column_title_gp = gpar_f(fontsize = 12, fontface = "bold"),
    rect_gp = gpar_f(col = "black", lwd = 0.5),
    border = TRUE,
    right_annotation = right_annotation,
    show_heatmap_legend = show_legend,
    heatmap_legend_param = heatmap_legend
  )
}

# =============================================================================
# 5. CARGA DE traits.xlsx
# =============================================================================
cat("Leyendo traits.xlsx...\n")
traits_energy <- read_excel(TRAITS_FILE, sheet = "energy")
traits_carbon <- read_excel(TRAITS_FILE, sheet = "carbon")

E <- build_order_matrix(traits_energy, "energy")
C <- build_order_matrix(traits_carbon, "carbon")
cat("\n")

# =============================================================================
# ---- 6-7. disabled: figures 1 and 2, energy and carbon as standalone panels ----
# =============================================================================
# Energy and carbon traits are no longer emitted separately; the only trait
# figure that goes on the poster is the combined side-by-side one in section 9.
# These stay commented because E and C are still built in section 5, so
# uncommenting is all it takes to get the standalone panels back.
#
# cat("=== Figura 1: energy traits by order ===\n")
# {
#   ht <- trait_heatmap(E$mat, name = "energy",
#                       right_annotation = n_mags_annotation(E$n))
#   draw_autosize(ht, file.path(OUTPUT_DIR, "Fig_traits_energy_by_order_EN.svg"),
#                 merge_legend = TRUE, padding = unit(rep(2, 4), "mm"))
# }
#
# cat("=== Figura 2: carbon-substrate traits by order ===\n")
# {
#   ht <- trait_heatmap(C$mat, name = "carbon",
#                       right_annotation = n_mags_annotation(C$n))
#   draw_autosize(ht, file.path(OUTPUT_DIR, "Fig_traits_carbon_by_order_EN.svg"),
#                 merge_legend = TRUE, padding = unit(rep(2, 4), "mm"))
# }

# =============================================================================
# 8. ORDEN DE FILAS COMPARTIDO POR LOS DOS PANELES
# =============================================================================
# Both panels cover the same four orders. For the combined figure to read as one
# unit, both are ordered by the same dendrogram, computed over the concatenated
# energy and carbon matrix so every trait counts
# igual porque todos son proporciones en [0,1]).
stopifnot(identical(sort(rownames(E$mat)), sort(rownames(C$mat))))
common_orders <- rownames(E$mat)
combined_mat  <- cbind(E$mat[common_orders, , drop = FALSE],
                       C$mat[common_orders, , drop = FALSE])
row_hc <- hclust(dist(combined_mat, method = "euclidean"), method = "ward.D2")

# =============================================================================
# ---- 9. figure 3: panels side by side ----
# =============================================================================
# Un solo bloque: dendrograma de filas a la izquierda, energy y carbon lado a
# lado compartiendo filas, barplot de No. MAGs y nombres de Orden a la derecha,
# one shared legend. Each order name carries its sample size as "n = X".
cat("=== Figura 3: combinada, paneles lado a lado ===\n")
{
  mat_e <- E$mat[common_orders, , drop = FALSE]
  mat_c <- C$mat[common_orders, , drop = FALSE]

  n_vals <- E$n[match(common_orders, rownames(E$mat))]

  # Cells narrower and taller than in the standalone figures. With 29 traits on
  # one row, 10 mm per column would give an unusable 16 inch canvas.
  ht_e <- trait_heatmap(mat_e, name = "energy_panel",
                        col_width_mm = 7, row_height_mm = 12,
                        cluster_rows = row_hc, show_row_dend = TRUE,
                        show_row_names = FALSE, show_legend = TRUE,
                        panel_title = "Energy metabolism")
  ht_c <- trait_heatmap(mat_c, name = "carbon_panel",
                        col_width_mm = 7, row_height_mm = 12,
                        cluster_rows = row_hc, show_row_dend = FALSE,
                        show_row_names = TRUE, row_n_vals = n_vals,
                        right_annotation = n_mags_annotation(n_vals),
                        show_legend = FALSE,
                        panel_title = "Carbon-substrate utilization")

  draw_autosize(ht_e + ht_c,
                file.path(OUTPUT_DIR, "Fig_traits_combined_sidebyside_EN.svg"),
                ht_gap = unit(8, "mm"), merge_legend = TRUE,
                padding = unit(rep(2, 4), "mm"))
}

# =============================================================================
# 9b. FIGURA 3v — LA MISMA COMBINADA, EN VERTICAL
# =============================================================================
# The same information transposed: orders move to the x axis as four columns,
# traits to the y axis, and the two panels stack with %v% instead of sitting
# side by side. The result is tall and narrow, which is the shape a poster
# column wants.
#
# %v% requires the two heatmaps to share columns. Transposed, the columns are
# the orders, which they do share, and that is why the vertical orientation
# needs a transpose rather than just a reorder. The order dendrogram is the same
# row_hc from the horizontal figure, applied here as a column clustering, so
# both figures tell the same story.
cat("=== Figura 3v: combinada, paneles apilados (vertical) ===
")

# Heatmap transpuesto: filas = rasgos, columnas = Órdenes.
trait_heatmap_t <- function(mat_t, name, col_width_mm = 11.5, row_height_mm = 5,
                            cluster_columns = TRUE, show_column_names = FALSE,
                            column_n_vals = NULL, top_annotation = NULL,
                            show_legend = TRUE, panel_title = NULL,
                            show_row_dend = FALSE, show_column_dend = FALSE) {
  col_lab <- NULL
  if (!is.null(column_n_vals)) {
    fmt_cn  <- vapply(colnames(mat_t), fmt_order_label, character(1))
    col_lab <- gt_render(sprintf(
      "%s<br><span style='font-size:%dpt;color:%s'>n = %d</span>",
      fmt_cn, N_LABEL_FONTSIZE, N_LABEL_COLOR, as.integer(column_n_vals)
    ))
  } else if (show_column_names) {
    col_lab <- vapply(colnames(mat_t), fmt_order_label, character(1))
  }
  Heatmap(mat_t,
    name = name,
    col = COL_FUN,
    # Traits cluster by similarity, orders follow the shared dendrogram. The
    # show_*_dend flags only control whether the trees are drawn; the clustering
    # is identical either way.
    cluster_rows = TRUE,
    show_row_dend = show_row_dend,
    row_dend_side = "left",
    row_dend_width = unit(10, "mm"),
    cluster_columns = cluster_columns,
    show_column_dend = show_column_dend,
    column_dend_side = "top",
    column_dend_height = unit(10, "mm"),
    show_row_names = TRUE,
    row_names_side = "right",
    row_names_gp = gpar_f(fontsize = 9),
    show_column_names = show_column_names,
    column_labels = col_lab,
    column_names_side = "bottom",
    column_names_gp = gpar_f(fontsize = ROW_LABEL_FONTSIZE),
    column_names_rot = 45,
    width  = unit(ncol(mat_t) * col_width_mm, "mm"),
    height = unit(nrow(mat_t) * row_height_mm, "mm"),
    # The panel title goes on the left edge: stacked vertically, that is the
    # only place it does not compete with the trait labels or the orders.
    row_title = panel_title,
    row_title_gp = gpar_f(fontsize = 12, fontface = "bold"),
    rect_gp = gpar_f(col = "black", lwd = 0.5),
    border = TRUE,
    top_annotation = top_annotation,
    show_heatmap_legend = show_legend,
    heatmap_legend_param = heatmap_legend
  )
}

{
  mat_e <- E$mat[common_orders, , drop = FALSE]
  mat_c <- C$mat[common_orders, , drop = FALSE]
  n_vals <- E$n[match(common_orders, rownames(E$mat))]

  # MAG counts as a top annotation here, where the horizontal figure puts them on the side.
  top_n_mags <- HeatmapAnnotation(
    `No. MAGs` = anno_barplot(n_vals,
      gp = gpar_f(fill = "#4A90D9", col = "black"),
      bar_width = 0.7, height = unit(2, "cm")),
    annotation_name_gp = gpar_f(fontsize = 9, fontface = "bold")
  )

  # with_dend draws the trait dendrogram to the left of each panel and the order
  # dendrogram above. The order one is row_hc, the same clustering that orders
  # the rows of the horizontal figure, so both tell the same story.
  build_vertical <- function(with_dend, with_bar = TRUE) {
    ht_e_v <- trait_heatmap_t(t(mat_e), name = "energy_panel_v",
                              cluster_columns = row_hc,
                              show_column_names = FALSE,
                              top_annotation = if (with_bar) top_n_mags else NULL,
                              show_legend = TRUE,
                              panel_title = "Energy metabolism",
                              show_row_dend = with_dend,
                              # Only the top panel carries the order
                              # dendrogram. It is shared, and repeating it
                              # between panels would break the stack.
                              show_column_dend = with_dend)
    ht_c_v <- trait_heatmap_t(t(mat_c), name = "carbon_panel_v",
                              cluster_columns = row_hc,
                              show_column_names = TRUE,
                              column_n_vals = n_vals,
                              show_legend = FALSE,
                              panel_title = "Carbon-substrate utilization",
                              show_row_dend = with_dend,
                              show_column_dend = FALSE)
    ht_e_v %v% ht_c_v
  }

  # Generous margin: ComplexHeatmap reserves the height of column names rotated
  # 45 degrees but not the part that overhangs to the left of the first column,
  # and "Undinarchaeales" was being clipped at the canvas edge.
  draw_autosize(build_vertical(FALSE),
                file.path(OUTPUT_DIR, "Fig_traits_combined_vertical_EN.svg"),
                margin_in = 0.6,
                ht_gap = unit(8, "mm"), merge_legend = TRUE,
                padding = unit(rep(2, 4), "mm"))

  draw_autosize(build_vertical(TRUE),
                file.path(OUTPUT_DIR, "Fig_traits_combined_vertical_dend_EN.svg"),
                margin_in = 0.6,
                ht_gap = unit(8, "mm"), merge_legend = TRUE,
                padding = unit(rep(2, 4), "mm"))

  # The same with dendrograms but without the MAG barplot. The sample size
  # already appears as "n = X" under each order name, so the barplot repeats it,
  # which is worth dropping when poster space is tight.
  draw_autosize(build_vertical(TRUE, with_bar = FALSE),
                file.path(OUTPUT_DIR, "Fig_traits_combined_vertical_dend_nobar_EN.svg"),
                margin_in = 0.6,
                ht_gap = unit(8, "mm"), merge_legend = TRUE,
                padding = unit(rep(2, 4), "mm"))
}

# =============================================================================
# 10. CALIDAD DE ANOTACIÓN — CLASIFICACIÓN DE GENES
# =============================================================================
# Ported verbatim from section 2 of DRAM2_output_analysis.R, so figures 4 and 5
# report exactly the same numbers as their Spanish counterparts.
#
# What counts as a real functional annotation: a gene qualifies only on a hit
# against a curated or HMM functional database (KEGG,
# KOfam, dbCAN, CAMPER, sulfur, methyl, FeGenie, metals). Un hit en UniRef NO
# does not count. It is sequence homology, and the reference protein may itself
# be uncharacterised, so those are reported separately as Homology_Only.
cat("=== Anotación: leyendo raw-annotations.tsv ===\n")
dram <- read_tsv(ANNOT_FILE, show_col_types = FALSE, na = c("", "NA"))
cat("  Genes anotados:", nrow(dram), "| MAGs:", n_distinct(dram$input_fasta), "\n")

count_hits <- function(column) sum(!is.na(column) & column != "")

is_informative_uniref <- function(uniref_id, uniref_desc) {
  if (is.na(uniref_id) | uniref_id == "") return(FALSE)
  if (is.na(uniref_desc) | uniref_desc == "") return(FALSE)

  non_informative_patterns <- c(
    "uncharacterized protein",
    "hypothetical protein",
    "predicted protein",
    "putative protein$",
    "DUF[0-9]+",
    "conserved protein of unknown function",
    "protein of unknown function"
  )

  desc_lower <- tolower(uniref_desc)
  !any(sapply(non_informative_patterns, function(p) grepl(p, desc_lower)))
}

dram_classified <- dram %>%
  mutate(
    has_kegg    = !is.na(kofam_id)   & kofam_id   != "",
    has_cazy    = !is.na(dbcan_id)   & dbcan_id   != "",
    has_metals  = !is.na(metals_id)  & metals_id  != "",
    has_camper  = !is.na(camper_id)  & camper_id  != "",
    has_methyl  = !is.na(methyl_id)  & methyl_id  != "",
    has_sulfur  = !is.na(sulfur_id)  & sulfur_id  != "",
    has_fegenie = !is.na(fegenie_id) & fegenie_id != "",
    has_uniref  = !is.na(uniref_id)  & uniref_id  != "",
    has_informative_uniref = mapply(is_informative_uniref, uniref_id, uniref_description),
    has_vague_uniref = has_uniref & !has_informative_uniref,
    has_kegg_db = if ("kegg_id" %in% names(pick(everything()))) !is.na(kegg_id) & kegg_id != "" else FALSE,
    # Curated functional annotation: any HMM or curated database, UniRef excluded.
    has_curated_functional = has_kegg | has_cazy | has_metals | has_camper |
      has_methyl | has_sulfur | has_fegenie | has_kegg_db,
    Annotation_Quality = case_when(
      has_curated_functional ~ "Curated_Functional",
      has_uniref             ~ "Homology_Only",
      TRUE                   ~ "No_Annotation"
    )
  )

summary_table <- dram_classified %>%
  summarise(
    Total_Genes = n(),
    Curated_Functional_Annotation = sum(Annotation_Quality == "Curated_Functional"),
    Homology_Only = sum(Annotation_Quality == "Homology_Only"),
    No_Annotation = sum(Annotation_Quality == "No_Annotation")
  )

summary_by_mag <- dram_classified %>%
  group_by(MAG = input_fasta) %>%
  summarise(
    Total_Genes = n(),
    Curated_Functional_Annotation = sum(Annotation_Quality == "Curated_Functional"),
    .groups = "drop"
  ) %>%
  left_join(
    metadata %>% select(NameWithoutExtension, FullTaxonomy),
    by = c("MAG" = "NameWithoutExtension")
  ) %>%
  mutate(Percent_Curated_Functional = (Curated_Functional_Annotation / Total_Genes) * 100)

cat("  Genes con anotación funcional curada:",
    summary_table$Curated_Functional_Annotation,
    sprintf("(%.1f%%)\n\n", 100 * summary_table$Curated_Functional_Annotation / summary_table$Total_Genes))

# =============================================================================
# ---- 11. figure 4: genes annotated per tool, horizontal ----
# =============================================================================
# plot_03_annotation_by_tool.png under coord_flip(): the tool labels no longer
# need a 45 degree rotation, and the count moves to the right end of each bar on
# a single line.
# The KofamKOALA and eggNOG-mapper v2 counts come from separate runs over the
# same 68,539 genes; the DRAM2 count comes from summary_table.
cat("=== Figura 4: annotated genes by tool (horizontal) ===\n")
{
  datos_anot <- data.frame(
    software = c("KofamKOALA", "eggNOG-mapper v2", "DRAM2"),
    genes = c(68539, 68539, summary_table$Total_Genes),
    hits  = c(21255, 18179, summary_table$Curated_Functional_Annotation)
  ) %>%
    mutate(pct_anotacion = (hits / genes) * 100)

  # reorder ascendente: tras coord_flip() la barra mayor queda arriba.
  p_tool <- ggplot(datos_anot, aes(x = reorder(software, hits), y = hits, fill = software)) +
    geom_bar(stat = "identity", color = "black", width = 0.45) +
    geom_text(
      aes(label = paste0(comma_sp(hits), " (", round(pct_anotacion, 1), "%)")),
      hjust = -0.12, size = 4.5, family = font_family
    ) +
    # Set3 with the original 0.8 alpha already flattened; sort() reproduces the
    # alphabetical order scale_fill_brewer() assigned the colours in.
    scale_fill_manual(values = setNames(
      blend_over_white(RColorBrewer::brewer.pal(3, "Set3"), 0.8),
      sort(unique(datos_anot$software))
    )) +
    coord_flip() +
    theme_publication(base_size = 14) +
    theme(
      axis.text.y = element_text(size = 12, face = "bold"),
      legend.position = "none",
      # after coord_flip the useful grid is the one on the count axis
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray", linetype = "dotted")
    ) +
    labs(
      subtitle = paste("Total genes analyzed:", comma_sp(datos_anot$genes[1])),
      # No axis title: the tool names identify themselves.
      x = NULL,
      y = "Annotated genes"
    ) +
    # 0.32 on the right: the label now runs horizontally and needs more room
    # than the 0.15 of the vertical version, where it wrapped onto two lines.
    scale_y_continuous(labels = comma_sp, expand = expansion(mult = c(0, 0.32)))

  out <- file.path(OUTPUT_DIR, "Fig_annotation_by_tool_EN.svg")
  ggsave_frames(out, p_tool, width = 9, height = 5)
}

# =============================================================================
# 12. FIGURA 5 — TASA DE ANOTACIÓN POR ORDEN TAXONÓMICO
# =============================================================================
# plot_03c_annotation_rate_per_order.png. Bars are the per-order mean percentage
# of genes with curated functional annotation, error bars are SD, the red line
# discontinua = media global.
#
# The test runs on the raw per-MAG values, one observation per MAG.
# Kruskal-Wallis es el test omnibus apropiado: no asume normalidad ni
# equal variance, and it tolerates the badly unbalanced group sizes, which run
# from n = 3 to n = 38. It only runs where at least two groups hold two MAGs.
cat("=== Figura 5: annotation rate per order ===\n")
{
  df <- summary_by_mag %>%
    mutate(Order = vapply(FullTaxonomy, .extract_rank, character(1), prefix = "o__")) %>%
    filter(!is.na(Order), !is.na(Percent_Curated_Functional))

  agg <- df %>%
    group_by(Taxon = Order) %>%
    summarise(
      n = n(),
      mean_pct = mean(Percent_Curated_Functional, na.rm = TRUE),
      sd_pct = if (n() > 1) sd(Percent_Curated_Functional, na.rm = TRUE) else 0,
      .groups = "drop"
    ) %>%
    mutate(
      # Same label style as the pangenome class figure: the name in bold, then
      # the sample size on a second line in smaller grey. Needs
      # ggtext::element_markdown() on axis.text.y.
      TaxonLabel = vapply(Taxon, fmt_order_label, character(1)),
      # 13px is about 9.5 pt: two points more than the original 10px (7.5 pt),
      # and still under the 11 pt of the order name above it.
      TaxonLabel = paste0("**", TaxonLabel, "**",
                          "<br><span style='font-size:13px;color:#555555'>n = ",
                          n, "</span>")
    ) %>%
    arrange(mean_pct) %>%
    mutate(TaxonLabel = factor(TaxonLabel, levels = TaxonLabel))

  global_mean <- mean(df$Percent_Curated_Functional, na.rm = TRUE)

  sig_info <- NULL
  grp_sizes   <- table(df$Order)
  usable_grps <- names(grp_sizes)[grp_sizes >= 2]
  if (length(usable_grps) >= 2) {
    test_df <- df %>% filter(Order %in% usable_grps) %>% mutate(.grp = factor(Order))
    kw <- kruskal.test(Percent_Curated_Functional ~ .grp, data = test_df)
    p_val <- unname(kw$p.value)
    sig_info <- list(
      p = p_val,
      label = if (is.na(p_val)) "ns"
              else if (p_val < 0.001) "***"
              else if (p_val < 0.01)  "**"
              else if (p_val < 0.05)  "*"
              else                    "ns",
      chi = unname(kw$statistic), df = unname(kw$parameter),
      excluded = setdiff(names(grp_sizes), usable_grps)
    )
    cat(sprintf("  Kruskal-Wallis: chi2(%d) = %.3f, p = %.4g  ->  %s%s\n",
      sig_info$df, sig_info$chi, sig_info$p, sig_info$label,
      if (length(sig_info$excluded) > 0)
        paste0("  (excluidos por n<2: ", paste(sig_info$excluded, collapse = ", "), ")")
      else ""))
  }

  y_max <- max(agg$mean_pct + agg$sd_pct, na.rm = TRUE)
  # A significant result needs room on the right, past the flipped end, for the
  # bracket and its asterisk.
  y_limit_mult <- if (!is.null(sig_info) && sig_info$label != "ns") 1.32 else 1.18

  p_order <- ggplot(agg, aes(x = TaxonLabel, y = mean_pct)) +
    geom_bar(stat = "identity", fill = "steelblue", color = "black", width = 0.6) +
    geom_errorbar(aes(ymin = pmax(0, mean_pct - sd_pct), ymax = mean_pct + sd_pct),
                  width = 0.2, color = "gray30") +
    # The label sits at the far right of the error bar, clear of the SD bar.
    geom_text(aes(y = mean_pct + sd_pct, label = sprintf("%.1f%%", mean_pct)),
              hjust = -0.15, size = 4.5, family = font_family) +
    geom_hline(yintercept = global_mean, linetype = "dashed",
               color = "red", linewidth = 0.8) +
    coord_flip() +
    theme_publication(base_size = 12) +
    theme(axis.text.y = ggtext::element_markdown(size = 11)) +
    labs(
      subtitle = paste0("Global mean: ", round(global_mean, 1), "%"),
      x = "Order",
      y = "Genes with curated functional annotation (%)"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.18)),
      limits = c(0, y_max * y_limit_mult)
    )

  # Significance drawn as a bracket spanning every group with a single asterisk
  # when the omnibus is significant. It reads vertical after coord_flip.
  if (!is.null(sig_info) && sig_info$label != "ns") {
    n_taxa    <- nrow(agg)
    y_bracket <- y_max * 1.20
    tick_len  <- y_max * 0.015
    y_text    <- y_bracket + y_max * 0.035
    p_order <- p_order +
      annotate("segment", x = 1, xend = n_taxa,
               y = y_bracket, yend = y_bracket, linewidth = 0.5, color = "black") +
      annotate("segment", x = 1, xend = 1,
               y = y_bracket - tick_len, yend = y_bracket, linewidth = 0.5, color = "black") +
      annotate("segment", x = n_taxa, xend = n_taxa,
               y = y_bracket - tick_len, yend = y_bracket, linewidth = 0.5, color = "black") +
      annotate("text", x = (1 + n_taxa) / 2, y = y_text, label = sig_info$label,
               size = 6, fontface = "bold", family = font_family, vjust = 0.5)
  }

  out  <- file.path(OUTPUT_DIR, "Fig_annotation_rate_per_order_EN.svg")
  # The axis label now wraps to two lines, so each bar gets a little more height.
  h_in <- max(3, 0.70 * nrow(agg) + 2)
  ggsave_frames(out, p_order, width = 10, height = h_in)
}

# =============================================================================
# 13. FIGURA 6 — PANEL CNPS: METADATOS + TAXONOMÍA POR ORDEN
# =============================================================================
# Ported from section G of CNPS.cycle.R, make_combined_meta_tax_g("Order"),
# which produces Panel_Meta_Tax_Orden.png. One ggplot faceted by
# facet_grid(Cycle ~ Panel): metadata on the left, station and depth, taxonomic
# orders on the right, both sharing the pathway y axis.
#   bubble size   prevalence, the percentage of the group's MAGs carrying the pathway
#   colour        mean abundance across the MAGs that do carry it
#
# As in the original, this reads the abundance CSVs already written under
# results/15_CNPS_cycle/ rather than in-memory objects, so the panel rebuilds
# without rerunning CNPS.cycle.R. The pathway names in those files are already
# English; the Spanish version is what sends them through DeepL, so there is no
# translation step here.
cat("=== Figura 6: CNPS metabolic panel by order ===\n")

# Cycles in CNPS order. Phosphorus drops out only if its matrix is all zeros.
CNPS_CYCLES <- list(
  Carbon     = file.path(CNPS_DIR, "carbon",     "C_abundance.csv"),
  Nitrogen   = file.path(CNPS_DIR, "nitrogen",   "N_abundance.csv"),
  Phosphorus = file.path(CNPS_DIR, "phosphorus", "P_abundance.csv"),
  Sulfur     = file.path(CNPS_DIR, "sulfur",     "S_abundance.csv")
)

# ROMZ stays as the sampling recorded it. BTM is expanded, because it does not
# read as a word.
DEPTH_LABELS_EN <- c("ROMZ" = "ROMZ", "1000" = "1000 m", "BTM" = "Bottom")
DEPTH_LEVELS_EN <- unname(DEPTH_LABELS_EN)

capitalize_first <- function(x) {
  ifelse(is.na(x), x, paste0(toupper(substring(x, 1, 1)), substring(x, 2)))
}

# Un CSV de abundancia (rutas × MAGs) to formato largo. NULL si no hay ni un
# any value above zero, which is how the original drops an empty cycle, here phosphorus.
read_cycle_long <- function(csv_path, cycle_name) {
  if (!file.exists(csv_path)) {
    cat("  [WARN] CSV no encontrado, ciclo omitido:", basename(csv_path), "\n")
    return(NULL)
  }
  ab <- read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE)
  num_cols <- setdiff(colnames(ab), "Type")
  if (!any(ab[, num_cols] > 0, na.rm = TRUE)) {
    cat("  · ciclo sin señal, omitido:", cycle_name, "\n")
    return(NULL)
  }
  ab %>%
    tidyr::pivot_longer(cols = -Type, names_to = "MAG", values_to = "Abundance") %>%
    mutate(Cycle = cycle_name, Pathway = capitalize_first(Type)) %>%
    select(Cycle, Pathway, MAG, Abundance)
}

# y axis labels: strip the "_<Cycle>" suffix that keeps same-named pathways
# unique across cycles, turn internal commas into colons
# ("Nitrification, ammonia to hydroxylamine" to "Nitrification: ammonia to ...")
# and abbreviate the chemical species to formulas.
#
# The formulas do not use Unicode subscripts (U+2080 to U+2089). The poster font
# is Lato, whose cmap covers neither that block nor the superscript minus
# (U+207B), so they would render as empty boxes. Each label compiles to a
# plotmath expression instead, SO[4]^{"2-"}, which grid composes from glyphs
# Lato does have and which rotates like plain text.

# Chemical species to plotmath fragment. The order matters: a pattern containing
# another as a substring, "thiosulfate" over "sulfate", "(sulfide)n" over
# "sulfide", has to be consumed first.
CHEM_PLOTMATH <- c(
  "\\(sulfide\\)n"  = '(S^{"2-"})[n]',   # polisulfuro
  "thio?sulfate"    = 'S[2]*O[3]^{"2-"}', # "thisulfate" es el typo de KEGG-decoder
  "sulfate"         = 'SO[4]^{"2-"}',
  "sulfite"         = 'SO[3]^{"2-"}',
  "sulfide"         = 'S^{"2-"}',
  "nitrous oxide"   = 'N[2]*O',
  "nitric oxide"    = 'NO',
  "nitrate"         = 'NO[3]^{"-"}',
  "nitrite"         = 'NO[2]^{"-"}',
  "hydroxylamine"   = 'NH[2]*OH',
  "hydrazine"       = 'N[2]*H[4]',
  "ammonia"         = 'NH[3]',
  "CH4"             = 'CH[4]',
  "@N2@"            = 'N[2]',            # marcadores internos (ver abajo)
  "@S0@"            = 'S^0'
)

# A cleaned label to a plotmath string ready for parse(). Each species is
# replaced by a \x01<i>\x01 marker, then reassembled alternating literals
# entrecomillados y fragmentos plotmath.
#
# A space touching a formula cannot stay inside the quoted string: plotmath
# trims the trailing space of every literal, so "X: " * Y
# sale pegado ("Nitrification:NH3"). Por eso el espacio de frontera se saca de la
# and is emitted as its own token. That token is two spaces, not one: plotmath
# measures space with math-mode metrics, and a single space comes out half width
# against a formula. The ~ operator is narrower still.
pathway_plotmath <- function(lab) {
  pats <- names(CHEM_PLOTMATH)
  for (i in seq_along(pats)) {
    lab <- gsub(pats[i], paste0("\x01", i, "\x01"), lab, ignore.case = TRUE)
  }
  parts <- strsplit(lab, "\x01", fixed = TRUE)[[1]]
  out <- character(0)
  gap <- FALSE   # ¿queda un espacio pendiente antes del siguiente token?
  for (p in parts) {
    if (p == "") next
    is_chem <- grepl("^[0-9]+$", p)
    if (is_chem) {
      tok   <- CHEM_PLOTMATH[[as.integer(p)]]
      trail <- FALSE
    } else {
      if (grepl("^ ", p)) gap <- TRUE
      trail <- grepl(" $", p)
      txt   <- trimws(p)
      if (txt == "") { gap <- TRUE; next }   # literal que era solo espacio
      tok <- paste0('"', txt, '"')
    }
    if (length(out)) out <- c(out, if (gap) '* "  " *' else "*")
    out <- c(out, tok)
    gap <- trail
  }
  if (length(out) == 0) return('""')
  paste(out, collapse = " ")
}

clean_pathway_label <- function(x) {
  lab <- sub("_(Carbon|Nitrogen|Phosphorus|Sulfur)$", "", x)
  # "organ-N" es la abreviatura de KEGG-decoder; la version ES tambien la
  # expande ("N-organico" -> "nitrogeno organico").
  lab <- gsub("organ-N", "organic N", lab)
  lab <- gsub(",", ":", lab)

  # The element heading a pathway is abbreviated to its symbol; the one
  # appearing as substrate or product is marked to receive a formula. "Carbon"
  # already arrives abbreviated from the CSV, as in "Aerobic C fixation", so the
  # only thing to do is not expand it again.
  lab <- sub("^Nitrogen ", "N ", lab)
  lab <- sub("^Sulfur ",   "S ", lab)
  # With the head consumed, any remaining "nitrogen" or "sulfur" sits in the
  # substrate or product part, where they name specific species: N2 and
  # elemental S0, not the generic element.
  lab <- gsub("\\bnitrogen\\b", "@N2@", lab)
  lab <- gsub("\\bsulfur\\b",   "@S0@", lab)

  parse(text = vapply(lab, pathway_plotmath, character(1), USE.NAMES = FALSE))
}

{
  df_long <- bind_rows(lapply(names(CNPS_CYCLES), function(cy) {
    read_cycle_long(CNPS_CYCLES[[cy]], cy)
  }))
  present_cycles <- intersect(names(CNPS_CYCLES), unique(df_long$Cycle))
  df_long <- df_long %>% mutate(Cycle = factor(Cycle, levels = present_cycles))

  # ---- metadata: station and depth, in long form ----
  df_meta <- metadata %>%
    select(MAG = NameWithoutExtension, Station, Depth) %>%
    mutate(
      Depth   = factor(unname(DEPTH_LABELS_EN[as.character(Depth)]),
                       levels = DEPTH_LEVELS_EN),
      Station = factor(Station)
    )

  df_meta_long <- df_long %>%
    inner_join(df_meta, by = "MAG") %>%
    tidyr::pivot_longer(cols = c(Station, Depth),
                        names_to = "Variable", values_to = "Group") %>%
    mutate(Variable = factor(Variable,
                             levels = c("Station", "Depth"),
                             labels = c("Sampling station", "Depth")))

  # Prevalence is the percentage of the group's MAGs carrying the pathway; mean
  # abundance averages only over the MAGs that do carry it, so zeros do not
  # wash out the colour.
  df_summary_meta <- df_meta_long %>%
    filter(!is.na(Group) & Group != "") %>%
    group_by(Cycle, Pathway, Variable, Group) %>%
    summarise(
      n_mags      = n_distinct(MAG),
      Prevalence  = sum(Abundance > 0) / n_mags * 100,
      Mean_abund  = mean(Abundance[Abundance > 0], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(Prevalence > 0) %>%
    arrange(Cycle, desc(Prevalence)) %>%
    mutate(Pathway_uniq = factor(paste(Pathway, Cycle, sep = "_"),
                                 levels = unique(paste(Pathway, Cycle, sep = "_"))))

  # The canonical y order, shared by all three panels.
  y_levels_shared <- levels(df_summary_meta$Pathway_uniq)

  # ---- GTDB taxonomy at order level ----
  df_tax <- read.delim(GTDB_FILE, stringsAsFactors = FALSE) %>%
    select(user_genome, classification) %>%
    mutate(MAG = gsub("_", ".", user_genome),
           # sub() on a string with no o__ would return the whole
           # classification, so it is filtered before extracting.
           Order = ifelse(grepl(";o__[^;]", classification),
                          sub(".*;o__([^;]*).*", "\\1", classification),
                          NA_character_)) %>%
    filter(!is.na(Order)) %>%
    select(MAG, Order)

  df_tax_sum <- df_long %>%
    inner_join(df_tax, by = "MAG") %>%
    group_by(Cycle, Pathway, Order) %>%
    summarise(
      n_mags     = n_distinct(MAG),
      Prevalence = sum(Abundance > 0) / n_mags * 100,
      Mean_abund = mean(Abundance[Abundance > 0], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(Prevalence > 0) %>%
    mutate(Pathway_uniq = factor(paste(Pathway, Cycle, sep = "_"),
                                 levels = y_levels_shared))

  # --- Un solo data.frame: metadatos a la izquierda, taxones a la derecha ---
  panel_levels <- c(levels(df_summary_meta$Variable), "Taxa")

  df_meta_part <- df_summary_meta %>%
    mutate(Panel = as.character(Variable), x_label = as.character(Group)) %>%
    select(Cycle, Pathway_uniq, Panel, x_label, Prevalence, Mean_abund)

  df_tax_part <- df_tax_sum %>%
    mutate(Panel = "Taxa",
           # Code-like order names such as MGIII get the class in parentheses,
           # as in the Spanish version. Ranks above genus stay upright, per CDC.
           x_label = vapply(Order, fmt_order_label, character(1))) %>%
    select(Cycle, Pathway_uniq, Panel, x_label, Prevalence, Mean_abund)

  df_combined <- bind_rows(df_meta_part, df_tax_part) %>%
    mutate(
      Pathway_uniq = factor(Pathway_uniq, levels = y_levels_shared),
      Panel        = factor(Panel, levels = panel_levels),
      # as.character() above dropped the factor levels, so Depth is refactored
      # to come out ROMZ, 1000 m, Bottom rather than alphabetical.
      x_label = factor(x_label, levels = c(
        sort(levels(df_meta$Station)),
        DEPTH_LEVELS_EN,
        sort(unique(df_tax_part$x_label))
      ))
    )

  n_station <- length(unique(df_meta_part$x_label[df_meta_part$Panel == panel_levels[1]]))
  n_depth   <- length(unique(df_meta_part$x_label[df_meta_part$Panel == panel_levels[2]]))
  n_taxa    <- length(unique(df_tax_part$x_label))
  w_panel <- (n_station + n_depth + n_taxa) * 0.65 + 5
  h_panel <- max(7, length(y_levels_shared) * 0.22 + 2)

  p_panel <- ggplot(df_combined, aes(x = x_label, y = Pathway_uniq)) +
    geom_point(aes(size = Prevalence, fill = Mean_abund),
               shape = 21, color = "gray20", stroke = 0.5, na.rm = TRUE) +
    scale_size_continuous(range = c(2, 8), name = "Prevalence\nin MAGs (%)",
                          limits = c(0, 100)) +
    # guide_colourbar() draws the colour bar as an embedded raster, a base64
    # <image> inside the SVG. Asking for rectangles keeps the figure
    # 100 % vectorial y la barra no se pixela al ampliar el cartel. El argumento
    # was renamed in ggplot2 3.5.0, raster to display, so the argument is chosen
    # from the installed version rather than pinned to one machine's.
    scale_fill_viridis_c(option = "plasma", name = "Mean\nabundance",
                         na.value = "transparent",
                         guide = if (utils::packageVersion("ggplot2") >= "3.5.0")
                                   guide_colourbar(display = "rectangles")
                                 else guide_colourbar(raster = FALSE)) +
    scale_x_discrete(expand = expansion(add = 0.4)) +
    scale_y_discrete(labels = clean_pathway_label) +
    # scales and space set to "free": each row shows only its own cycle's
    # pathways and each column only its own groups, with no filler cells.
    facet_grid(Cycle ~ Panel, scales = "free", space = "free") +
    labs(x = NULL, y = NULL) +
    theme_publication(base_size = 13) +
    theme(
      axis.text.x  = element_text(angle = 45, hjust = 1, size = 11),
      # The same 45 degree rotation as the group labels on the x axis, so both
      # families of labels in the panel read at the same angle.
      axis.text.y  = element_text(angle = 45, hjust = 1, vjust = 0.5, size = 10),
      strip.text.x = element_text(face = "bold", size = 12),
      strip.text.y = element_text(face = "bold", size = 12),
      legend.position = "right",
      legend.box = "vertical"
    )

  ggsave_both(file.path(OUTPUT_DIR, "Panel_Meta_Tax_Order_EN.svg"),
              p_panel, width = w_panel, height = h_panel)
}

# =============================================================================
# 14. FIGURA 7 — MAGs RECUPERADOS POR PROFUNDIDAD, DESGLOSADOS POR ORDEN
# =============================================================================
# Horizontal stacked bars, one per depth (ROMZ, 1000 m, Bottom), divided by the
# order each MAG classified into. The source is
# MAGs_metadata.tsv: columna Depth y el rango o__ de FullTaxonomy, exactamente
# the same two variables that feed the section 13 panel, so the totals agree
# with the rest of the poster: 2 + 13 + 39 = 54 MAGs.
#
# The order palette comes from the project's pruned phylogenomic tree:
# Poseidoniales teal, Nitrososphaerales orange, MGIII purple, Undinarchaeales
# yellow. Poster and tree then read as one piece.
cat("=== Figura 7: MAGs per depth by taxonomic order ===\n")

# Two versions of the figure come out of this palette:
#   dark   the palette through darken(), with white counts inside the bars.
#          White would be unreadable on the original yellow, #E6AB02.
#   light  the tree palette untouched, with black counts.
ORDER_COLORS_LIGHT <- c(
  "Poseidoniales"     = "#1B9E77",
  "Nitrososphaerales" = "#D95F02",
  "MGIII"             = "#7570B3",
  "Undinarchaeales"   = "#E6AB02"
)
ORDER_COLORS <- darken(ORDER_COLORS_LIGHT, factor = 0.75)

{
  df_depth <- metadata %>%
    transmute(
      MAG   = NameWithoutExtension,
      Depth = factor(unname(DEPTH_LABELS_EN[as.character(Depth)]),
                     levels = DEPTH_LEVELS_EN),
      Order = vapply(FullTaxonomy, .extract_rank, character(1), prefix = "o__")
    ) %>%
    filter(!is.na(Depth), !is.na(Order))

  # Orders sort by total abundance. Anything absent from the palette, none at
  # present, lands at the end in grey rather than breaking the figure.
  order_levels <- df_depth %>% count(Order, sort = TRUE) %>% pull(Order)
  pick_fill <- function(pal) {
    f <- pal[order_levels]
    f[is.na(f)] <- "#BBBBBB"
    names(f) <- order_levels
    f
  }

  df_bars <- df_depth %>%
    count(Depth, Order, name = "n_mags") %>%
    mutate(Order = factor(Order, levels = order_levels)) %>%
    # `Depth` se invierte porque en el eje discreto vertical el primer nivel se
    # draws at the bottom, which puts ROMZ at the top and Bottom below, shallow to deep.
    mutate(Depth = factor(Depth, levels = rev(DEPTH_LEVELS_EN)))

  df_tot <- df_bars %>%
    group_by(Depth) %>%
    summarise(total = sum(n_mags), .groups = "drop")

  # The count goes inside the segment only when it fits. Under 3 MAGs the text
  # overflows the segment width, and ROMZ contributes two single-MAG sections.
  MIN_SEG_LABEL <- 3L

  x_max <- max(df_tot$total)

  # dark = TRUE gives the darkened palette with white numbers inside the bar;
  # dark = FALSE gives the tree palette with black numbers.
  build_depth_plot <- function(dark) {
    order_fill <- pick_fill(if (dark) ORDER_COLORS else ORDER_COLORS_LIGHT)
    seg_label_color <- if (dark) "white" else "black"

  ggplot(df_bars, aes(x = n_mags, y = Depth, fill = Order)) +
    geom_col(width = 0.62, color = "black", linewidth = 0.35) +
    geom_text(
      data = df_bars %>% filter(n_mags >= MIN_SEG_LABEL),
      aes(label = n_mags),
      position = position_stack(vjust = 0.5),
      color = seg_label_color, fontface = "bold", size = 6.4, family = font_family
    ) +
    # Row total at the right end of each bar.
    geom_text(
      data = df_tot, aes(x = total, y = Depth, label = total),
      inherit.aes = FALSE, hjust = -0.35, fontface = "bold",
      size = 6.6, family = font_family
    ) +
    scale_fill_manual(
      values = order_fill,
      labels = function(x) vapply(x, fmt_order_label, character(1)),
      name   = NULL
    ) +
    # limits already reserves room for the total label; expand at zero stops
    # ggplot adding its usual 5% per side on top and leaving half the figure
    # empty.
    scale_x_continuous(expand = expansion(mult = c(0, 0)),
                       limits = c(0, x_max * 1.10),
                       breaks = pretty(c(0, x_max))) +
    # No axis title: ROMZ, 1000 m and Bottom read as depths on their own.
    labs(x = "Recovered MAGs", y = NULL) +
    theme_publication(base_size = 14) +
    theme(
      axis.text.y        = element_text(size = 13, face = "bold"),
      # Legend underneath, as in the pangenome class figure.
      legend.position    = "bottom",
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray", linetype = "dotted")
    )
  }

  ggsave_frames(file.path(OUTPUT_DIR, "Fig_MAGs_per_depth_by_order_EN.svg"),
                build_depth_plot(TRUE), width = 9, height = 4.4)
  ggsave_frames(file.path(OUTPUT_DIR, "Fig_MAGs_per_depth_by_order_light_EN.svg"),
                build_depth_plot(FALSE), width = 9, height = 4.4)
}

# =============================================================================
# 15. FIGURAS 8-10 — ANI INTRA-GRUPO EN LOS GRUPOS SIN ESPECIE ASIGNADA (N/A)
# =============================================================================
# Portado de ANI_intraspecies.R (Fig1b_ANI_distribucion_grupos_NA.png). Se
# This reproduces that script's data loading, the sparse skani edge list plus
# the pairs rescued from the novel_vs_all run, so the sample sizes shown as
# "n = 5 | 4 pairs" match the Spanish version exactly.
#
# Se emiten tres figuras:
#   8   a literal translation of Fig1b, within-group ANI only
#   9   Fig1b plus each MAG's ANI against its closest GTDB reference genome,
#        overlaid as triangles on the same axis
#   10  two stacked tracks per group, within-group against reference, compared
#        as parallel boxplots
# Figures 9 and 10 are the two takes on showing ANI against the reference:
# 9 keeps the original geometry and adds a layer, 10 rearranges the plot so the
# two distributions sit face to face.
cat("=== Figuras 8-10: ANI in unassigned (N/A) species groups ===\n")

# theme_publication() from ANI_intraspecies.R, which scales differently from the
# one in this script: axis.text at base-1, axis.title at base+1, grid in gray88.
# Copied verbatim so these figures match their Spanish counterparts.
theme_publication_ani <- function(base_size = 13) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "gray88", linetype = "dotted"),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text          = ggplot2::element_text(color = "black", size = base_size - 1),
      axis.title         = ggplot2::element_text(face = "bold", size = base_size + 1),
      legend.position    = "right",
      legend.background  = ggplot2::element_rect(fill = "white", color = NA),
      legend.title       = ggplot2::element_text(face = "bold", size = base_size - 1),
      legend.text        = ggplot2::element_text(size = base_size - 3),
      strip.background   = ggplot2::element_rect(fill = "gray90", color = "black"),
      strip.text         = ggplot2::element_text(face = "bold", size = base_size - 1),
      plot.title         = ggplot2::element_text(face = "bold", hjust = 0.5, size = base_size + 2),
      plot.subtitle      = ggplot2::element_text(hjust = 0.5, size = base_size - 1, color = "gray40"),
      text               = ggplot2::element_text(family = font_family)
    )
}

# Minimum aligned fraction for accepting a pair from the novel_vs_all run
# (--min-af 0). Same value as NOVEL_AF_MIN in ANI_intraspecies.R.
NOVEL_AF_MIN <- 5

# Umbral de especie procariota (95 % ANI).
ANI_SPECIES_THR <- 95

# Colours for the three unassigned groups, saturated as in Fig1b.
NA_PALETTE <- c(
  "MGIIb-Q1 N/A"         = "#2471A3",
  "Nitrosopelagicus N/A" = "#1E8449",
  "UBA36 N/A"            = "#717D7E"
)

# Colours for the two comparison tracks, figures 9 and 10.
TRACK_COLORS <- c(
  "Within group (MAG-MAG)"   = "#1B9E77",
  "MAG vs. closest GTDB ref." = "#D95F02"
)

{
  # ---- metadata with the species group ----
  meta_ani <- metadata %>%
    transmute(
      MAG   = NameWithoutExtension,
      Genus, Species,
      Species_group = if_else(Species == "N/A", paste0(Genus, " N/A"), Species)
    )

  # --- Edge list sparse de skani ----------------------------------------------
  sparse_raw <- read_tsv(ANI_SPARSE_FILE, show_col_types = FALSE) %>%
    transmute(
      Ref_MAG   = tools::file_path_sans_ext(basename(as.character(Ref_file))),
      Query_MAG = tools::file_path_sans_ext(basename(as.character(Query_file))),
      ANI       = as.numeric(ANI),
      AF_ref    = as.numeric(Align_fraction_ref),
      AF_query  = as.numeric(Align_fraction_query)
    )

  # ---- rescue pairs from novel_vs_all ----
  # The all_mags run used skani's default --min-af of about 15%, so MAG pairs
  # with a low aligned fraction never appear. novel10_vs_all54.tsv ran with
  # --min-af 0 and does catch them. Filtering at AF >= NOVEL_AF_MIN discards the
  # artefacts where ANI reads 100% off one tiny contig, and only pairs absent
  # from the edge list are added.
  unord_key <- function(a, b) ifelse(a < b, paste(a, b, sep = "|"),
                                     paste(b, a, sep = "|"))
  if (file.exists(ANI_NOVEL_FILE)) {
    novel_pairs <- read_tsv(ANI_NOVEL_FILE, show_col_types = FALSE) %>%
      transmute(
        Ref_MAG   = tools::file_path_sans_ext(basename(as.character(Ref_file))),
        Query_MAG = tools::file_path_sans_ext(basename(as.character(Query_file))),
        ANI       = as.numeric(ANI),
        AF_ref    = as.numeric(Align_fraction_ref),
        AF_query  = as.numeric(Align_fraction_query)
      ) %>%
      filter(Ref_MAG != Query_MAG) %>%
      mutate(AF_min = pmin(AF_ref, AF_query, na.rm = TRUE)) %>%
      filter(AF_min >= NOVEL_AF_MIN) %>%
      filter(!unord_key(Ref_MAG, Query_MAG) %in%
               unord_key(sparse_raw$Ref_MAG, sparse_raw$Query_MAG)) %>%
      distinct(pair = unord_key(Ref_MAG, Query_MAG), .keep_all = TRUE) %>%
      select(Ref_MAG, Query_MAG, ANI, AF_ref, AF_query)
    sparse_raw <- bind_rows(sparse_raw, novel_pairs)
    cat("  ", nrow(novel_pairs), " pares añadidos desde novel_vs_all\n", sep = "")
  } else {
    cat("  [WARN] novel_vs_all no encontrado; sin rescate de pares\n")
  }

  # ---- within-genus pairs among the unassigned MAGs ----
  sparse_df <- sparse_raw %>%
    left_join(meta_ani, by = c("Ref_MAG" = "MAG")) %>%
    rename(Genus_ref = Genus, Species_ref = Species, Species_group_ref = Species_group) %>%
    left_join(meta_ani, by = c("Query_MAG" = "MAG")) %>%
    rename(Genus_query = Genus, Species_query = Species) %>%
    filter(Genus_ref == Genus_query, Species_ref == "N/A", Species_query == "N/A")

  # Unassigned groups with at least 2 MAGs, the minimum for a single pair.
  na_counts <- meta_ani %>% filter(Species == "N/A") %>% count(Species_group)
  na_ge2    <- na_counts %>% filter(n >= 2) %>% pull(Species_group)

  fig_ani_data <- sparse_df %>%
    filter(Species_group_ref %in% na_ge2) %>%
    transmute(ANI, group_label = Species_group_ref)

  stopifnot(nrow(fig_ani_data) > 0)

  # ---- each MAG against its closest GTDB reference ----
  # For a MAG with no species, GTDB-Tk leaves closest_genome_ani empty and often
  # closest_placement_ani too: no reference genome falls inside the species
  # radius. The value does survive in other_related_references, formatted as
  # "accession, s__species, radius, ANI, AF" and separated by semicolons. The
  # maximum ANI across all three sources is what gets used.
  # Advertencia: algunos de esos hits traen AF muy baja (0.15-0.17 en
  # Nitrosopelagicus), meaning the ANI rests on a small fraction of the genome.
  # The value is kept because it is what GTDB-Tk reports, but it should not
  # be read with the same confidence as an ANI backed by AF above 0.8.
  gtdb_raw <- read.delim(GTDB_FILE, stringsAsFactors = FALSE,
                         check.names = FALSE, na.strings = c("", "N/A"))
  other_col <- grep("other_related_references", colnames(gtdb_raw), value = TRUE)[1]

  best_other_ani <- function(s) {
    if (is.na(s) || !nzchar(trimws(s))) return(NA_real_)
    entries <- strsplit(s, ";", fixed = TRUE)[[1]]
    vals <- vapply(entries, function(e) {
      f <- trimws(strsplit(e, ",", fixed = TRUE)[[1]])
      if (length(f) < 4) return(NA_real_)
      suppressWarnings(as.numeric(f[4]))
    }, numeric(1), USE.NAMES = FALSE)
    if (all(is.na(vals))) NA_real_ else max(vals, na.rm = TRUE)
  }

  ref_ani <- gtdb_raw %>%
    transmute(
      OriginalName = user_genome,
      ani_genome    = suppressWarnings(as.numeric(closest_genome_ani)),
      ani_placement = suppressWarnings(as.numeric(closest_placement_ani)),
      ani_other     = vapply(.data[[other_col]], best_other_ani, numeric(1),
                             USE.NAMES = FALSE)
    ) %>%
    rowwise() %>%
    mutate(Ref_ANI = suppressWarnings(
      max(c(ani_genome, ani_placement, ani_other), na.rm = TRUE))) %>%
    ungroup() %>%
    mutate(Ref_ANI = ifelse(is.finite(Ref_ANI), Ref_ANI, NA_real_)) %>%
    select(OriginalName, Ref_ANI)

  df_ref <- metadata %>%
    filter(Species == "N/A") %>%
    transmute(MAG = NameWithoutExtension, OriginalName,
              group_label = paste0(Genus, " N/A")) %>%
    filter(group_label %in% na_ge2) %>%
    left_join(ref_ani, by = "OriginalName")

  n_no_ref <- sum(is.na(df_ref$Ref_ANI))
  cat("  MAGs sin especie asignada:", nrow(df_ref),
      "| con ANI a referencia GTDB:", sum(!is.na(df_ref$Ref_ANI)),
      "| sin referencia reportada:", n_no_ref, "\n")
  df_ref_ok <- df_ref %>% filter(!is.na(Ref_ANI))

  # ---- axis labels: italic genus plus sample size ----
  ani_summary <- fig_ani_data %>%
    group_by(group_label) %>%
    summarise(median_ani = median(ANI), n_pairs = n(), .groups = "drop") %>%
    left_join(
      meta_ani %>% filter(Species == "N/A") %>% count(Species_group, name = "n_mags"),
      by = c("group_label" = "Species_group")
    ) %>%
    arrange(median_ani) %>%
    mutate(
      genus_part = sub(" N/A$", "", group_label),
      # La version ES escribe "1 pares"; en ingles el singular si canta, asi
      # so the noun agrees with the count.
      axis_label = paste0(
        taxon_italic(genus_part, level = "Genus"),
        "<br><span style='font-size:13px;color:#555555'>n=",
        n_mags, " | ", n_pairs, ifelse(n_pairs == 1, " pair", " pairs"), "</span>"
      )
    )
  group_levels <- ani_summary$group_label
  axis_labels  <- setNames(ani_summary$axis_label, group_levels)

  fig_ani_data <- fig_ani_data %>%
    mutate(group_fct = factor(group_label, levels = group_levels))
  df_ref_ok <- df_ref_ok %>%
    mutate(group_fct = factor(group_label, levels = group_levels))

  # x axis breaks shared by all three figures: dense above 95%, where the
  # species call is decided, and sparse below it.
  ani_breaks <- sort(unique(c(80, seq(85, 95, by = 5), seq(95, 100, by = 1))))

  # ---------------------------------------------------------------------------
  # Disabled: figure 8, a literal translation of Fig1b showing within-group ANI
  # alone. The poster uses only the version below, which overlays ANI against the
  # GTDB reference. Kept commented in case the
  # figura equivalente a la de la tesis.
  # ---------------------------------------------------------------------------
  # # ---------------------------------------------------------------------------
  # # Figure 8, a literal translation of Fig1b
  # # ---------------------------------------------------------------------------
  # p_ani_base <- ggplot(fig_ani_data, aes(x = ANI, y = group_fct)) +
  #   geom_boxplot(aes(fill = group_fct), color = "gray35",
  #                alpha = 0.25, outlier.shape = NA, width = 0.35, linewidth = 0.5) +
  #   geom_jitter(aes(color = group_fct), height = 0.10, size = 2.5,
  #               alpha = 0.88, stroke = 0.3) +
  #   geom_vline(xintercept = ANI_SPECIES_THR, color = "#C0392B",
  #              linetype = "dashed", linewidth = 0.85) +
  #   scale_x_continuous(breaks = ani_breaks, labels = function(x) paste0(x, "%")) +
  #   scale_y_discrete(labels = axis_labels) +
  #   scale_color_manual(values = NA_PALETTE, guide = "none") +
  #   scale_fill_manual(values = NA_PALETTE, guide = "none") +
  #   labs(x = "ANI (%)", y = NULL) +
  #   theme_publication_ani(base_size = 13) +
  #   theme(
  #     axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5),
  #     axis.text.y = ggtext::element_markdown(size = 12.5),
  #     panel.grid.major.x = element_line(color = "gray88", linetype = "dotted"),
  #     panel.grid.major.y = element_blank()
  #   )
  #
  # ggsave_both(file.path(OUTPUT_DIR, "Fig_ANI_unassigned_groups_EN.svg"),
  #             p_ani_base, width = 8, height = 5)

  # ---------------------------------------------------------------------------
  # Figure 9, the only ANI figure that goes on the poster.
  #
  # Boxplot and points for within-group ANI, with each MAG's ANI against its
  # closest GTDB reference drawn on the same row. The reference marker always
  # falls left of the 95% line and the within-group pairs right of it, so the
  # separation reads without dodging one series against the other.
  # make_ani_overlay() takes the reference marker shape as an argument, which is
  # what lets the same plot come out in two versions.
  # ---------------------------------------------------------------------------
  # The boxplot alpha is premixed over white. On a transparent background, a
  # fill carrying alpha lets the poster colour through and shifts the palette.
  na_fill_flat <- blend_over_white(NA_PALETTE, 0.25)
  names(na_fill_flat) <- names(NA_PALETTE)

  REF_LABEL    <- "MAG vs. closest GTDB ref."
  WITHIN_LABEL <- "Within group (MAG-MAG)"

  # Alternative colour for the reference marker. The plot already carries blue
  # (#2471A3), green (#1E8449) and grey (#717D7E) on the groups, and red
  # (#C0392B) on the threshold line. Purple is the gap left in the colour wheel
  # that collides with none of them.
  REF_FILL_ALT <- "#7D3C98"

  make_ani_overlay <- function(ref_shape, ref_fill = TRACK_COLORS[[2]]) {
    ggplot(fig_ani_data, aes(x = ANI, y = group_fct)) +
      geom_boxplot(aes(fill = group_fct), color = "gray35",
                   outlier.shape = NA, width = 0.35, linewidth = 0.5) +
      geom_jitter(aes(color = group_fct, shape = WITHIN_LABEL),
                  height = 0.10, size = 4.5, stroke = 0.3) +
      # Sin position_nudge(): la serie de referencia comparte exactamente la
      # the same y position as its group's within-group pairs.
      geom_point(
        data = df_ref_ok,
        aes(x = Ref_ANI, y = group_fct, shape = REF_LABEL),
        inherit.aes = FALSE, size = 4.4, fill = ref_fill,
        color = "gray20", stroke = 0.4
      ) +
      geom_vline(xintercept = ANI_SPECIES_THR, color = "#C0392B",
                 linetype = "dashed", linewidth = 0.85) +
      scale_x_continuous(breaks = ani_breaks, labels = function(x) paste0(x, "%")) +
      scale_y_discrete(labels = axis_labels) +
      scale_color_manual(values = NA_PALETTE, guide = "none") +
      scale_fill_manual(values = na_fill_flat, guide = "none") +
      scale_shape_manual(
        values = setNames(c(16, ref_shape), c(WITHIN_LABEL, REF_LABEL)),
        breaks = c(WITHIN_LABEL, REF_LABEL),
        name = NULL
      ) +
      # No subtitle; the caption that sat above the plot was removed. The two
      # Nitrosopelagicus MAGs with no GTDB reference at all
      # reportada siguen ausentes de la capa naranja; conviene decirlo en el pie
      # belong in the poster caption, not in the figure.
      labs(x = "ANI (%)", y = NULL) +
      guides(shape = guide_legend(override.aes = list(
        size = 4.4, color = "gray20",
        fill = c("gray20", ref_fill)))) +
      theme_publication_ani(base_size = 13) +
      theme(
        axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = ggtext::element_markdown(size = 12.5),
        # Darker and heavier dotted grid than the gray88 of Fig1b: against the
        # poster background a fine dotted line disappeared.
        panel.grid.major.x = element_line(color = "gray55", linetype = "dotted",
                                          linewidth = 0.9),
        panel.grid.major.y = element_blank(),
        # Framed like the _axislines variants of the other figures, axis only,
        # Y izquierdo y eje X inferior, en negro— pero conservando la rejilla,
        # which here does help in reading each point's ANI.
        panel.border = element_blank(),
        axis.line    = element_line(color = "black", linewidth = 0.6),
        axis.ticks   = element_line(color = "black", linewidth = 0.6),
        legend.position = "bottom",
        # Legend tight against the plot; ggplot leaves about 11 pt by default.
        legend.box.spacing = unit(2, "pt"),
        legend.margin      = margin(t = 0, r = 0, b = 0, l = 0)
      )
  }

  # Four variants of the same plot: two marker shapes by two colours.
  # Shape 23 is a filled diamond, 22 a filled square. Both carry a border so
  # they do not read as the filled circle (16) used for the within-group pairs.
  #
  # Disabled: the triangle version, shape 24, replaced by the diamond
  # ggsave_both(file.path(OUTPUT_DIR, "Fig_ANI_unassigned_groups_withref_A_EN.svg"),
  #             make_ani_overlay(24), width = 9, height = 5.2)
  ani_variants <- list(
    list(file = "Fig_ANI_unassigned_groups_withref_A_diamond_EN.svg",
         shape = 23, fill = TRACK_COLORS[[2]]),
    list(file = "Fig_ANI_unassigned_groups_withref_A_square_EN.svg",
         shape = 22, fill = TRACK_COLORS[[2]]),
    list(file = "Fig_ANI_unassigned_groups_withref_A_diamond_purple_EN.svg",
         shape = 23, fill = REF_FILL_ALT),
    list(file = "Fig_ANI_unassigned_groups_withref_A_square_purple_EN.svg",
         shape = 22, fill = REF_FILL_ALT)
  )
  for (v in ani_variants) {
    ggsave_both(file.path(OUTPUT_DIR, v$file),
                make_ani_overlay(v$shape, v$fill), width = 9, height = 5.2)
  }

  # ---------------------------------------------------------------------------
  # Disabled: figure 10, two tracks per group, within-group against reference,
  # compared as parallel boxplots. Dropped in favour of the single-row overlay
  # above. Kept commented.
  # ---------------------------------------------------------------------------
  # # ---------------------------------------------------------------------------
  # # Figure 10, variant B: two tracks per group, compared face to face.
  # # Cada grupo N/A ocupa dos renglones — ANI intra-grupo y ANI contra la
  # # closest GTDB reference, in the same geometry of boxplot and points, so
  # # the separation of the two distributions around the 95% line reads
  # # directly, without one layer sitting on top of the other.
  # # ---------------------------------------------------------------------------
  # track_levels <- names(TRACK_COLORS)
  # df_tracks <- bind_rows(
  #   fig_ani_data %>% transmute(group_label, ANI, Track = track_levels[1]),
  #   df_ref_ok    %>% transmute(group_label, ANI = Ref_ANI, Track = track_levels[2])
  # ) %>%
  #   mutate(
  #     group_fct = factor(group_label, levels = group_levels),
  #     # Niveles invertidos: en un eje y discreto el dodge dibuja el primer nivel
  #     # at the bottom, and here "within group" should sit on top, matching
  #     # the order the legend lists them in, which breaks fixes.
  #     Track     = factor(Track, levels = rev(track_levels))
  #   )
  #
  # track_fill_flat <- blend_over_white(TRACK_COLORS, 0.25)
  # names(track_fill_flat) <- track_levels
  #
  # # Axis labels drop the "n = ... pairs": each group here holds two samples of
  # # different size, pairs against MAGs, so the count goes on the track itself
  # # at the end of each box.
  # axis_labels_tracks <- setNames(
  #   vapply(group_levels, function(g)
  #     taxon_italic(sub(" N/A$", "", g), level = "Genus"), character(1)),
  #   group_levels
  # )
  #
  # n_labels <- df_tracks %>%
  #   group_by(group_fct, Track) %>%
  #   summarise(n = n(), x = max(ANI), .groups = "drop") %>%
  #   mutate(label = paste0("n=", n))
  #
  # p_ani_tracks <- ggplot(df_tracks, aes(x = ANI, y = group_fct, fill = Track)) +
  #   geom_boxplot(color = "gray35", outlier.shape = NA, linewidth = 0.5,
  #                width = 0.6, position = position_dodge(width = 0.72)) +
  #   geom_point(aes(color = Track), size = 2.4, stroke = 0.3,
  #              position = position_jitterdodge(jitter.height = 0.08,
  #                                              jitter.width = 0,
  #                                              dodge.width = 0.72)) +
  #   geom_text(data = n_labels, aes(x = x, y = group_fct, label = label),
  #             inherit.aes = FALSE, size = 3.2, family = font_family,
  #             color = "gray30", hjust = -0.35,
  #             position = position_dodge2(width = 0.72)) +
  #   geom_vline(xintercept = ANI_SPECIES_THR, color = "#C0392B",
  #              linetype = "dashed", linewidth = 0.85) +
  #   scale_x_continuous(breaks = ani_breaks, labels = function(x) paste0(x, "%"),
  #                      expand = expansion(mult = c(0.04, 0.07))) +
  #   scale_y_discrete(labels = axis_labels_tracks) +
  #   scale_fill_manual(values = track_fill_flat, breaks = track_levels, name = NULL) +
  #   scale_color_manual(values = TRACK_COLORS, breaks = track_levels, name = NULL) +
  #   labs(x = "ANI (%)", y = NULL,
  #        subtitle = sprintf(
  #          "Dashed line: 95%% species threshold%s",
  #          if (n_no_ref > 0)
  #            sprintf("  ·  %d MAG%s with no GTDB reference reported are absent from the reference track",
  #                    n_no_ref, if (n_no_ref == 1) "" else "s")
  #          else "")) +
  #   theme_publication_ani(base_size = 13) +
  #   theme(
  #     axis.text.x = element_text(size = 11, angle = 90, hjust = 1, vjust = 0.5),
  #     axis.text.y = ggtext::element_markdown(size = 12.5),
  #     panel.grid.major.x = element_line(color = "gray88", linetype = "dotted"),
  #     panel.grid.major.y = element_blank(),
  #     legend.position = "bottom",
  #     plot.subtitle = element_text(hjust = 0, size = 10, color = "gray35")
  #   )
  #
  # ggsave_both(file.path(OUTPUT_DIR, "Fig_ANI_unassigned_groups_withref_B_EN.svg"),
  #             p_ani_tracks, width = 9, height = 5.6)
}

# =============================================================================
# 16. FIGURA 11 — CLASES PANGENÓMICAS POR ESPECIE (BARRAS 100 % APILADAS)
# =============================================================================
# The same information as Fig_panel_pan_classes.png from
# 27_pangenome_classification.R, in one horizontal chart: a bar per species from
# 0 to 100%, divided by pangenome class. Against the three-subplot panel, this
# makes relative composition comparable at a glance, and the cost is losing
# absolute magnitude. So each section carries its gene cluster count inside it
# and the total goes in the axis label.
#
# Source: class_counts_summary.tsv, which script 27 writes with the counts already
# clasificados (core / soft-core / shell / cloud). No se recalcula nada.
# In species with 7 MAGs or fewer, soft-core collapses into core by script 27's
# own definition, which is why only Thalassarchaeum shows a soft-core section.
cat("=== Figura 11: pangenome classes by species (100% stacked) ===\n")

# Species names as the figures set them, per CDC style: italic genus, upright
# GTDB sp... code, and the codes MGIIb-O1 and UBA8886 upright throughout.
SPECIES_PRETTY_EN <- c(
  Thalassarchaeum_sp002507125 = "*Thalassarchaeum* sp002507125",
  MGIIb_O1_sp002498525        = "MGIIb-O1 sp002498525",
  UBA8886_sp003193815         = "UBA8886 sp003193815"
)

{
  pan_raw <- read_tsv(PAN_CLASS_FILE, show_col_types = FALSE)

  pan_df <- pan_raw %>%
    filter(species %in% names(SPECIES_PRETTY_EN)) %>%
    mutate(class = factor(class, levels = pan_class_levels)) %>%
    # Classes missing from the TSV, soft-core and cloud where the count is zero,
    # are filled in so the legend and the stacking order stay identical across
    # tres barras.
    tidyr::complete(species, class, fill = list(n = 0L)) %>%
    group_by(species) %>%
    mutate(
      n_genomes = max(n_genomes, na.rm = TRUE),
      total_gc  = sum(n),
      pct       = n / total_gc * 100
    ) %>%
    ungroup()

  pan_totals <- pan_df %>% distinct(species, n_genomes, total_gc)

  # y axis follows the canonical SPECIES_PRETTY_EN order, reversed so the first
  # species sits at the top. The sample size goes under the name.
  sp_levels <- rev(intersect(names(SPECIES_PRETTY_EN), pan_df$species))
  sp_labels <- setNames(
    vapply(sp_levels, function(sp) {
      tt <- pan_totals[pan_totals$species == sp, ]
      # Wrapping a name that already carries *genus* in **...** nests the
      # markdown correctly: bold across the whole line, italic on the genus
      # alone, which is the CDC rule.
      paste0("**", SPECIES_PRETTY_EN[[sp]], "**",
             "<br><span style='font-size:12px;color:#555555'>n = ", tt$n_genomes,
             " | ", comma_sp(tt$total_gc), " gene clusters</span>")
    }, character(1)),
    sp_levels
  )

  pan_df <- pan_df %>%
    mutate(species_fct = factor(species, levels = sp_levels))

  # The count goes inside its section only where there is room, 6% of the width or more.
  MIN_PCT_LABEL <- 6

  # reverse = TRUE, because without it ggplot stacks the first factor level
  # last, putting cloud on the left and core on the right. That inverts both the
  # Spanish panel and the core-to-cloud gradient this is meant to show.
  #
  # Two versions, as in the MAGs-per-depth figure:
  #   `dark = TRUE`  -> pan_class_colors oscurecida + cifras BLANCAS dentro de la
  #                     the bar. White had no contrast against the original
  #                     yellow and light green.
  #   `dark = FALSE` -> pan_class_colors original de _pub_style.R + cifras NEGRAS.
  build_pan_plot <- function(dark) {
    pan_fill    <- if (dark) darken(pan_class_colors, factor = 0.75) else pan_class_colors
    label_color <- if (dark) "white" else "black"

  ggplot(pan_df, aes(x = pct, y = species_fct, fill = class)) +
    geom_col(width = 0.6, color = "black", linewidth = 0.35,
             position = position_stack(reverse = TRUE)) +
    geom_text(
      data = pan_df %>% filter(pct >= MIN_PCT_LABEL),
      aes(label = sprintf("%s\n%.1f%%", comma_sp(n), pct)),
      position = position_stack(vjust = 0.5, reverse = TRUE),
      color = label_color, fontface = "bold", size = 5.9,
      family = font_family, lineheight = 0.95
    ) +
    scale_fill_manual(
      values = pan_fill, drop = FALSE, name = NULL,
      labels = c(core = "Core", `soft-core` = "Soft-core",
                 shell = "Shell", cloud = "Cloud")
    ) +
    # Without a margin the "100%" at the right edge is clipped against the
    # lienzo; un 1.5 % de holgura basta y no desalinea la barra.
    scale_x_continuous(limits = c(0, 100), breaks = seq(0, 100, 20),
                       labels = function(x) paste0(x, "%"),
                       expand = expansion(mult = c(0, 0.015))) +
    scale_y_discrete(labels = sp_labels) +
    labs(x = "Gene clusters (%)", y = NULL) +
    theme_publication(base_size = 14) +
    theme(
      axis.text.y        = ggtext::element_markdown(size = 13, hjust = 1),
      legend.position    = "bottom",
      # Upright: pangenome class names are not taxon names.
      legend.text        = element_text(size = 12),
      panel.grid.major.y = element_blank(),
      panel.grid.major.x = element_line(color = "gray", linetype = "dotted")
    )
  }

  ggsave_frames(file.path(OUTPUT_DIR, "Fig_pangenome_classes_stacked_EN.svg"),
                build_pan_plot(TRUE), width = 11, height = 4.4)
  ggsave_frames(file.path(OUTPUT_DIR, "Fig_pangenome_classes_stacked_light_EN.svg"),
                build_pan_plot(FALSE), width = 11, height = 4.4)
}

cat("\n[OK] MDMS_figs.R completado — figuras de fondo transparente en", OUTPUT_DIR, "\n")
for (ext in OUT_FORMATS) {
  n <- length(list.files(OUTPUT_DIR, pattern = paste0("_transparent\\", ext, "$")))
  cat(sprintf("    %-4s %2d archivos%s\n", toupper(sub("^\\.", "", ext)), n,
              if (identical(ext, ".png")) paste0(" a ", PNG_DPI, " dpi") else ""))
}
cat("Fin:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
