# Biogeochemical cycles (C, N, P, S) across the archaeal MAGs.
#
# Data:    54 archaeal MAGs from marine metagenomes
# Groups:  station (A4, A6, C22, F37, PO1, Y6)
#          depth (ROMZ, 1000 m, BTM)
# Runs locally; the data comes from the results tree.
#
# Built on Yue et al. mSystems 2025, doi:10.1128/msystems.01021-25
#
# OUTPUT PER CYCLE (C, N, S; P when there is data):
#   *_abundance.csv                      pathway abundance table
#   *_stats_station.csv / *_stats_depth.csv – p-values (ANOVA/Kruskal-Wallis)
#   *_heatmap_grouped_station.png        abun.heatmap.g by station, package version
#   *_heatmap_grouped_depth.png          abun.heatmap.g by depth, package version
#   *_heatmap_samples.png                abun.heatmap.s per MAG, package version
#   *_heatmap_MAGs_depth_station.png     ComplexHeatmap per MAG, depth and station sidebar
#   *_heatmap_MAGs_depth.png             ComplexHeatmap per MAG, depth sidebar only
#   *_heatmap_MAGs_station.png           ComplexHeatmap per MAG, station sidebar only
#   *_heatmap_MAGs_no_sidebar.png        ComplexHeatmap per MAG, no sidebar, free clustering
#   *_heatmap_by_station.png             ComplexHeatmap, mean abundance per station
#   *_taxonomy_{phylum,family,genus}.png ComplexHeatmap, taxa by pathway
#   host_relative/{PW}_{Rank}.csv        host-relative data per pathway and rank
#   {PW}_host_ratio.png                  stacked bars, phylum by group, every pathway
#
# Under 10_MAGs/ only two ComplexHeatmap variants are produced: depth, and no sidebar.
#
# DESIGN NOTES:
#   The Xcyc.abundance() functions read the global `ko` directly rather than
#   taking it as an argument. That is how the package behaves at v1.0.1, and it
#   is why `ko` is built in the global environment below.
#   Package bug 1, in rbind: a KO absent from `ko` yields rep(0, ncol(ko)) with
#   no names, which then misaligns the bind. Worked around by padding `ko` with
#   explicit zero rows for every KO the cycle functions ask for.
#   Package bug 2, in log10: values of 1 or less give log10(1) = 0, and then
#   seq(a + 1, b, 1) fails.
#   Worked around by scaling `ko` by 100 before the cycle functions run.

script_start <- Sys.time()
cat(sprintf("=== Inicio del script: %s ===\n\n", format(script_start, "%Y-%m-%d %H:%M:%S")))

# ---- 0. Packages ----

library(CNPS.cycle)
library(dplyr)
library(tidyr)
library(ggplot2)
library(reshape2) # melt(), used internally by abun.heatmap.g
library(RColorBrewer)
library(pheatmap) # used by host.ratio()
library(grid)
library(deeplr) # DeepL translations
library(ComplexHeatmap)
library(circlize) # colorRamp2
library(showtext)
library(sysfonts)
library(ggtext)

# ---- 0a-b. fonts ----
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

# ---- translation layer (DeepL plus an on-disk cache) ----
# Same approach as DRAM2_output_analysis.R and treatment_comparison.R. The key
# is read from the environment; set DEEPL_AUTH_KEY in ~/.Renviron. Without it
# the manual table below still applies and everything else stays in English.
DEEPL_AUTH_KEY <- Sys.getenv("DEEPL_AUTH_KEY")
if (nzchar(DEEPL_AUTH_KEY)) {
  options(deeplr.auth_key = DEEPL_AUTH_KEY)
} else {
  warning("DEEPL_AUTH_KEY is not set; labels stay in English.")
}

.MANUAL_CNPS <- list(
  # Carbon cycle
  "Aerobic C fixation"              = "Fijación aeróbica de C",
  "Aerobic CH4 oxidation"           = "Oxidación aeróbica de CH\u2084",
  "Aerobic respiration"             = "Respiración aeróbica",
  "Anaerobic C fixation"            = "Fijación anaeróbica de C",
  "CO oxidation"                    = "Oxidación de CO",
  "Fermentation"                    = "Fermentación",
  "Methanogenesis"                  = "Metanogénesis",
  # Nitrogen cycle
  "Nitrification"                   = "Nitrificación",
  "Denitrification"                 = "Desnitrificación",
  "Anammox"                         = "Anammox",
  "Nitrogen fixation"               = "Fijación de N\u2082",
  "Assimilatory nitrate reduction"  = "Reducción asimilativa de NO\u2083\u207b",
  "Dissimilatory nitrate reduction" = "Reducción disimilativa de NO\u2083\u207b",
  "Nitrogen mineralization"         = "Mineralización del N",
  "Nitrogen assimilation"           = "Asimilación de N",
  "DNRA"                            = "DNRA",
  # Phosphorus cycle
  "PhoR-PhoB system"                = "Sistema PhoR-PhoB",
  "SenX3-RegX3 system"              = "Sistema SenX3-RegX3",
  # Sulfur cycle
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

capitalize_first <- function(x) {
  ifelse(is.na(x), x, paste0(toupper(substring(x, 1, 1)), substring(x, 2)))
}

translate_to_spanish <- function(names_vec, translation_cache_file = NULL) {
  if (!is.null(translation_cache_file) && file.exists(translation_cache_file)) {
    cache_df <- read.csv(translation_cache_file, stringsAsFactors = FALSE)
    cache <- as.list(setNames(cache_df$Spanish, cache_df$English))
  } else {
    cache <- list()
  }
  translated <- character(length(names_vec))
  new_trans <- 0L
  safe_str <- function(x, fb) if (length(x) == 0 || is.na(x[1]) || x[1] == "") fb else as.character(x[1])

  for (i in seq_along(names_vec)) {
    orig <- names_vec[i]
    if (orig %in% names(.MANUAL_CNPS)) {
      translated[i] <- .MANUAL_CNPS[[orig]]
      cache[[orig]] <- .MANUAL_CNPS[[orig]]
    } else if (orig %in% names(cache)) {
      translated[i] <- safe_str(cache[[orig]], orig)
    } else if (nzchar(DEEPL_AUTH_KEY)) {
      res <- tryCatch(
        deeplr::translate2(text = orig, source_lang = "EN", target_lang = "ES"),
        error = function(e) orig
      )
      res <- safe_str(res, orig)
      translated[i] <- res
      cache[[orig]] <- res
      new_trans <- new_trans + 1L
      Sys.sleep(0.1)
    } else {
      translated[i] <- orig
    }
  }
  if (!is.null(translation_cache_file) && new_trans > 0) {
    ks <- names(cache)
    vs <- vapply(cache, function(x) safe_str(x, ""), character(1))
    ok <- !is.na(ks) & ks != ""
    write.csv(data.frame(English = ks[ok], Spanish = vs[ok], stringsAsFactors = FALSE),
      translation_cache_file,
      row.names = FALSE
    )
  }
  translated
}

# ---- 0c. depths: canonical order and palette ----

DEPTH_LEVELS         <- c("ROMZ", "1000", "BTM")
DEPTH_DISPLAY_LEVELS <- c("ZMOR", "1000 m", "Fondo")   # display labels
{
  pal <- RColorBrewer::brewer.pal(3L, "Set2")
  DEPTH_COLORS <- setNames(pal, DEPTH_LEVELS)
}

# ---- 1. Paths ----

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

WSL_BASE <- PROJECT_DIR
RES_DIR <- file.path(WSL_BASE, "results")
OUT_DIR <- file.path(RES_DIR, "SO_01", "CNPS_cycle")
CNPS_CACHE_FILE <- file.path(OUT_DIR, "input_data", "translation_cache_cnps.csv")

# The KofamKOALA submission date is part of the directory name (step 08), so
# the most recent all_kofam_KOs.byMAG.tsv under it is taken, as step 09 does.
F_KO_BYMAG <- local({
  hits <- sort(list.files(file.path(RES_DIR, "08_KO_assignment"),
                          pattern = "^all_kofam_KOs\\.byMAG\\.tsv$",
                          recursive = TRUE, full.names = TRUE))
  if (length(hits) == 0) {
    stop("No all_kofam_KOs.byMAG.tsv under ", file.path(RES_DIR, "08_KO_assignment"),
         ". Run scripts/utils/03_post_KO_assign_concat_KOs.sh first.", call. = FALSE)
  }
  tail(hits, 1)
})
F_DRAM <- file.path(RES_DIR, "10_MAG_distill/ANNOTATE/raw-annotations.tsv")
F_META <- file.path(WSL_BASE, "data", "MAGs_metadata.tsv")

# ---- 2. Output directories ----

out_subdirs <- c(
  "input_data",
  "carbon", "carbon/host_relative",
  "nitrogen", "nitrogen/host_relative",
  "phosphorus", "phosphorus/host_relative",
  "sulfur", "sulfur/host_relative",
  "ordination",
  # The 10 MAGs with no species assignment
  "10_MAGs",
  "10_MAGs/carbon", "10_MAGs/carbon/host_relative",
  "10_MAGs/nitrogen", "10_MAGs/nitrogen/host_relative",
  "10_MAGs/phosphorus", "10_MAGs/phosphorus/host_relative",
  "10_MAGs/sulfur", "10_MAGs/sulfur/host_relative",
  "10_MAGs/ordination"
)
for (d in file.path(OUT_DIR, out_subdirs)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# ---- 3. Load the raw data ----

message("[1/4] Cargando metadatos...")
meta <- read.table(F_META,
  sep = "\t", header = TRUE,
  stringsAsFactors = FALSE, quote = ""
)
all_mags <- sort(meta$NameWithoutExtension)

message("[2/4] Cargando tabla KO por MAG...")
ko_long <- read.table(F_KO_BYMAG,
  sep = "\t", header = FALSE,
  col.names = c("MAG", "KO"), stringsAsFactors = FALSE
)

message("[3/4] Cargando anotaciones DRAM2 (~30 s)...")
dram_raw <- read.table(F_DRAM,
  sep = "\t", header = TRUE,
  quote = "", fill = TRUE, stringsAsFactors = FALSE,
  comment.char = ""
)
dram_ko <- dram_raw[, c("query_id", "input_fasta", "kofam_id", "kegg_id")]
rm(dram_raw)
gc()

# ---- 4. Build the `ko` table ----

message("[4/4] Preparando tablas CNPS.cycle...")

ko_wide <- ko_long |>
  dplyr::mutate(value = 1L) |>
  tidyr::pivot_wider(id_cols = KO, names_from = MAG, values_from = value, values_fill = 0L) |>
  as.data.frame()
rownames(ko_wide) <- ko_wide$KO
ko <- ko_wide[, -1]
ko <- ko[, intersect(all_mags, colnames(ko)), drop = FALSE]

# Package bug 1: pad with zero rows for the KOs that are absent
kos_required <- unique(c(
  "K00855", "K01061", "K01602", "K08684", "K02256", "K02262", "K02274", "K02276",
  "K00174", "K00175", "K00244", "K01648", "K00194", "K00197", "K03518", "K03519", "K03520",
  "K00016", "K00400", "K00401",
  "K00370", "K00371", "K00374", "K02567", "K02568", "K00362", "K00363", "K03385",
  "K15876", "K00367", "K00372", "K00366", "K00368", "K15864", "K04561", "K02305",
  "K00376", "K02588", "K02586", "K02591", "K10944", "K10945", "K10946", "K10535",
  "K05601", "K20932", "K20933", "K20934", "K20935", "K01915", "K00265", "K00266",
  "K00264", "K00284", "K00260", "K15371", "K00261", "K00262", "K01455", "K01501",
  "K01725", "K00926", "K00549", "K02575", "K15576", "K15577", "K15578", "K15579",
  "K07636", "K07657", "K07658", "K07768", "K07776",
  "K00956", "K00957", "K00955", "K00860", "K00390", "K00380", "K00381", "K00392",
  "K00958", "K00988", "K00394", "K00395", "K11180", "K11181", "K17222", "K17223",
  "K17224", "K17225", "K17226", "K17227", "K17229", "K17230", "K17218", "K08352",
  "K08354", "K04091", "K00299", "K16968", "K16969", "K15762", "K15765", "K03119",
  "K00456", "K17217", "K02045", "K02046", "K02047", "K02048", "K15551", "K15552",
  "K10831", "K15553", "K15554", "K15555", "K01739", "K10764", "K01738", "K17228"
))
kos_missing <- setdiff(kos_required, rownames(ko))
if (length(kos_missing) > 0) {
  zero_rows <- as.data.frame(matrix(
    0L,
    nrow = length(kos_missing), ncol = ncol(ko),
    dimnames = list(kos_missing, colnames(ko))
  ))
  ko <- rbind(ko, zero_rows)
  message("  ko: ", length(kos_missing), " filas de ceros añadidas")
}

# Package bug 2: scale by 100 so log10 never produces seq(1, 0, 1)
ko <- ko * 100L
message("  ko escalado ×100")

# ---- 5. Build the `Gene` table ----

gene_ko_table <- dram_ko |>
  dplyr::mutate(KO = ifelse(kofam_id != "" & !is.na(kofam_id), kofam_id, kegg_id)) |>
  dplyr::filter(!is.na(KO) & KO != "") |>
  dplyr::select(Entry = KO, GeneID = query_id) |>
  dplyr::distinct()

Gene <- gene_ko_table |>
  dplyr::group_by(Entry) |>
  dplyr::mutate(index = dplyr::row_number()) |>
  tidyr::pivot_wider(names_from = Entry, values_from = GeneID) |>
  dplyr::select(-index) |>
  as.data.frame()

write.table(gene_ko_table, file.path(OUT_DIR, "input_data", "gene_ko_entries.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

# ---- 6. Build the `tax` table ----

mag_tax_lookup <- setNames(meta$FullTaxonomy, meta$NameWithoutExtension)
annotated_genes <- dram_ko |>
  dplyr::mutate(KO = ifelse(kofam_id != "" & !is.na(kofam_id), kofam_id, kegg_id)) |>
  dplyr::filter(!is.na(KO) & KO != "") |>
  dplyr::select(GeneID = query_id, MAG = input_fasta) |>
  dplyr::distinct()

tax <- annotated_genes |>
  dplyr::mutate(V2 = mag_tax_lookup[MAG]) |>
  dplyr::select(V1 = GeneID, V2) |>
  dplyr::distinct(V1, .keep_all = TRUE) |>
  as.data.frame()

# ---- 7. Build the `abundance` table ----

abundance <- annotated_genes |>
  dplyr::mutate(present = 1L) |>
  tidyr::pivot_wider(
    id_cols = GeneID, names_from = MAG, values_from = present,
    values_fill = 0L
  ) |>
  as.data.frame()
colnames(abundance)[1] <- "V1"
for (m in setdiff(all_mags, colnames(abundance)[-1])) abundance[[m]] <- 0L
abundance <- abundance[, c("V1", all_mags)]

# ---- 8. Build the `group` tables and the 10-MAG subset ----

group_station <- data.frame(
  ID = meta$NameWithoutExtension,
  Group = factor(meta$Station),
  stringsAsFactors = FALSE
)
group_depth <- data.frame(
  ID = meta$NameWithoutExtension,
  Group = factor(meta$Depth, levels = DEPTH_LEVELS),
  stringsAsFactors = FALSE
)

# The 10 MAGs GTDB left without a species (Species empty or "s__")
mags_10 <- meta$NameWithoutExtension[
  is.na(meta$Species) |
    trimws(as.character(meta$Species)) %in% c("", "s__", "N/A", "unknown", "Unknown", "NA")
]
message(sprintf("  10-MAG subset: %d MAGs sin especie identificados", length(mags_10)))

meta_10 <- meta[meta$NameWithoutExtension %in% mags_10, ]
group_10_station <- data.frame(
  ID = meta_10$NameWithoutExtension,
  Group = factor(meta_10$Station),
  stringsAsFactors = FALSE
)
group_10_depth <- data.frame(
  ID = meta_10$NameWithoutExtension,
  Group = factor(meta_10$Depth, levels = DEPTH_LEVELS),
  stringsAsFactors = FALSE
)

# Write the input tables
write.table(ko, file.path(OUT_DIR, "input_data", "ko.tsv"),
  sep = "\t", quote = FALSE
)
write.table(group_station, file.path(OUT_DIR, "input_data", "group_station.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(group_depth, file.path(OUT_DIR, "input_data", "group_depth.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

message(sprintf(
  "\nTablas listas: ko=%dx%d  Gene=%d KOs  tax=%d genes  abundance=%dx%d\n",
  nrow(ko), ncol(ko), ncol(Gene), nrow(tax),
  nrow(abundance), ncol(abundance) - 1L
))

# ---- LOOKUP HELPERS ----

get_mag_depth <- function(mags) meta$Depth[match(mags, meta$NameWithoutExtension)]
get_mag_station <- function(mags) meta$Station[match(mags, meta$NameWithoutExtension)]
get_mag_display_name <- function(mags) meta$DisplayName[match(mags, meta$NameWithoutExtension)]

# ---- Save the host-relative results across all six taxonomic ranks ----

save_host_relative <- function(result_list, out_subdir, prefix) {
  ranks <- c("Phylum", "Class", "Order", "Family", "Genus", "Species")
  for (i in seq_along(ranks)) {
    df <- result_list[[i]]
    if (!is.null(df) && is.data.frame(df) && nrow(df) > 0) {
      write.csv(df, file.path(OUT_DIR, out_subdir, paste0(prefix, "_", ranks[i], ".csv")))
    }
  }
}

# ---- HELPER: restrict abundance_df to a subset of MAGs ----

subset_abund <- function(abundance_df, mag_subset) {
  avail <- intersect(mag_subset, setdiff(colnames(abundance_df), "Type"))
  if (length(avail) == 0) {
    return(NULL)
  }
  abundance_df[, c("Type", avail), drop = FALSE]
}

# Biogeochemical cycle heatmap, styled to match the DRAM2 figures
#
#   rows    MAGs by DisplayName, ordered ROMZ, 1000, BTM, then by station
#   cols    pathways, labelled in Spanish
#   colour  white to red, normalised completeness 0 to 1
#   ann_v   0 no sidebar, rows cluster freely
#           1 depth and station, 2 depth only, 3 station only

make_cnps_heatmap <- function(abundance_df, title_es, outfile, ann_v = 1L) {
  type_col <- as.character(abundance_df$Type)
  mat_pw <- as.matrix(abundance_df[, setdiff(colnames(abundance_df), "Type"), drop = FALSE])
  rownames(mat_pw) <- type_col

  keep <- rowSums(mat_pw, na.rm = TRUE) > 0
  if (sum(keep) == 0) {
    message("    ", title_es, " (ann_v=", ann_v, "): sin pathways – omitido")
    return(invisible(NULL))
  }
  mat_pw <- mat_pw[keep, , drop = FALSE] / 100 # scale to 0-1

  pw_es <- capitalize_first(
    translate_to_spanish(rownames(mat_pw), translation_cache_file = CNPS_CACHE_FILE)
  )
  rownames(mat_pw) <- pw_es

  mat <- t(mat_pw) # rows = MAGs
  mag_names <- rownames(mat)

  # Station palette (Set1), computed before any reordering
  all_depths <- factor(get_mag_depth(mag_names), levels = DEPTH_LEVELS)
  all_stations <- get_mag_station(mag_names)
  st_uniq <- sort(unique(na.omit(as.character(all_stations))))
  n_st <- max(3L, length(st_uniq))
  st_cols <- setNames(RColorBrewer::brewer.pal(n_st, "Set1")[seq_along(st_uniq)], st_uniq)

  ann_gp <- gpar(fontsize = 10, fontface = "bold")
  ann_sz <- unit(0.4, "cm")

  if (ann_v == 0L) {
    # No sidebar, so rows cluster freely
    display_names <- get_mag_display_name(mag_names)
    display_names[is.na(display_names)] <- mag_names[is.na(display_names)]
    rownames(mat) <- display_names
    ra <- NULL
    do_cluster <- TRUE
  } else {
    # Order ROMZ, 1000, BTM, then by station
    ord_idx <- order(all_depths, all_stations)
    mat <- mat[ord_idx, , drop = FALSE]
    depths <- all_depths[ord_idx]
    stations <- all_stations[ord_idx]
    valid_d <- as.character(na.omit(unique(depths)))

    display_names <- get_mag_display_name(mag_names[ord_idx])
    display_names[is.na(display_names)] <- mag_names[ord_idx][is.na(display_names)]
    rownames(mat) <- display_names
    do_cluster <- FALSE

    if (ann_v == 1L) {
      ra <- rowAnnotation(
        Profundidad = depths,
        `Estación` = factor(stations, levels = st_uniq),
        col = list(Profundidad = DEPTH_COLORS, `Estación` = st_cols),
        annotation_name_side = "top", simple_anno_size = ann_sz,
        annotation_name_gp = ann_gp, show_legend = TRUE,
        annotation_legend_param = list(
          Profundidad = list(
            title = "Profundidad", at = valid_d, labels = valid_d,
            title_gp = ann_gp, labels_gp = gpar(fontsize = 9),
            border = "black"
          ),
          `Estación` = list(
            title = "Estación", title_gp = ann_gp,
            labels_gp = gpar(fontsize = 9), border = "black"
          )
        )
      )
    } else if (ann_v == 2L) {
      ra <- rowAnnotation(
        Profundidad = depths,
        col = list(Profundidad = DEPTH_COLORS),
        annotation_name_side = "top", simple_anno_size = ann_sz,
        annotation_name_gp = ann_gp, show_legend = TRUE,
        annotation_legend_param = list(
          Profundidad = list(
            title = "Profundidad", at = valid_d, labels = valid_d,
            title_gp = ann_gp, labels_gp = gpar(fontsize = 9),
            border = "black"
          )
        )
      )
    } else {
      ra <- rowAnnotation(
        `Estación` = factor(stations, levels = st_uniq),
        col = list(`Estación` = st_cols),
        annotation_name_side = "top", simple_anno_size = ann_sz,
        annotation_name_gp = ann_gp, show_legend = TRUE,
        annotation_legend_param = list(
          `Estación` = list(
            title = "Estación", title_gp = ann_gp,
            labels_gp = gpar(fontsize = 9), border = "black"
          )
        )
      )
    }
  }

  col_fun <- colorRamp2(
    seq(0, 1, length.out = 6),
    c("white", "yellow", "orange", "orangered", "red", "darkred")
  )

  n_mags <- nrow(mat)
  n_pathways <- ncol(mat)
  row_fs <- if (n_mags <= 15) 9 else 7

  ht <- Heatmap(
    mat,
    name = "Abundancia",
    col = col_fun,
    cluster_rows = do_cluster,
    show_row_dend = do_cluster,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_side = "right",
    row_names_gp = gpar(fontsize = row_fs),
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 9),
    column_names_rot = 45,
    left_annotation = ra,
    rect_gp = gpar(col = "black", lwd = 0.4),
    border = TRUE,
    column_title = title_es,
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(
      title            = "Abundancia\nnormalizada",
      at               = c(0, 0.25, 0.5, 0.75, 1),
      labels           = c("0%", "25%", "50%", "75%", "100%"),
      title_gp         = gpar(fontsize = 10, fontface = "bold"),
      labels_gp        = gpar(fontsize = 9),
      border           = "black",
      legend_direction = "vertical"
    )
  )

  w <- max(8, round(n_pathways * 0.55 + 5))
  h <- max(6, round(n_mags * 0.20 + 3))
  png(outfile, width = w, height = h, units = "in", res = 300)
  draw(ht, merge_legend = TRUE)
  dev.off()
  message("    \u2713 ", basename(outfile))
}

# Mean abundance per station
#
#   rows    stations, alphabetical
#   cols    pathways, labelled in Spanish
#   colour  white to red, mean normalised abundance 0 to 1

make_cnps_station_heatmap <- function(abundance_df, title_es, outfile) {
  type_col <- as.character(abundance_df$Type)
  mat_pw <- as.matrix(abundance_df[, setdiff(colnames(abundance_df), "Type"), drop = FALSE])
  rownames(mat_pw) <- type_col

  keep <- rowSums(mat_pw, na.rm = TRUE) > 0
  if (sum(keep) == 0) {
    message("    ", title_es, " (por estación): sin pathways – omitido")
    return(invisible(NULL))
  }
  mat_pw <- mat_pw[keep, , drop = FALSE] / 100 # scale to 0-1

  pw_es <- capitalize_first(
    translate_to_spanish(rownames(mat_pw), translation_cache_file = CNPS_CACHE_FILE)
  )
  rownames(mat_pw) <- pw_es

  # MAGs × pathways
  mat <- t(mat_pw)
  mag_names <- rownames(mat)
  stations_vec <- get_mag_station(mag_names)

  station_levels <- sort(unique(na.omit(as.character(stations_vec))))
  if (length(station_levels) == 0) {
    message("    ", title_es, " (por estación): sin estaciones – omitido")
    return(invisible(NULL))
  }

  # Normalised mean per station, giving stations by pathways
  mat_st <- do.call(rbind, lapply(station_levels, function(st) {
    cols <- which(stations_vec == st & !is.na(stations_vec))
    if (length(cols) == 0) {
      return(rep(0, ncol(mat)))
    }
    if (length(cols) == 1) {
      return(mat[cols, ])
    }
    colMeans(mat[cols, , drop = FALSE], na.rm = TRUE)
  }))
  rownames(mat_st) <- station_levels

  n_st <- nrow(mat_st)
  n_pw <- ncol(mat_st)

  col_fun <- colorRamp2(
    seq(0, 1, length.out = 6),
    c("white", "yellow", "orange", "orangered", "red", "darkred")
  )

  # Rows already carry the station names, so no sidebar is needed
  ht <- Heatmap(
    mat_st,
    name = "Abundancia",
    col = col_fun,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_side = "right",
    row_names_gp = gpar(fontsize = 10),
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 9),
    column_names_rot = 45,
    rect_gp = gpar(col = "black", lwd = 0.4),
    border = TRUE,
    column_title = paste0(title_es, " - por estaci\u00f3n"),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(
      title            = "Abundancia\nmedia normalizada",
      at               = c(0, 0.25, 0.5, 0.75, 1),
      labels           = c("0%", "25%", "50%", "75%", "100%"),
      title_gp         = gpar(fontsize = 10, fontface = "bold"),
      labels_gp        = gpar(fontsize = 9),
      border           = "black",
      legend_direction = "vertical"
    )
  )

  w <- max(8, round(n_pw * 0.55 + 5))
  h <- max(4, round(n_st * 0.5 + 3))
  png(outfile, width = w, height = h, units = "in", res = 300)
  draw(ht, merge_legend = TRUE)
  dev.off()
  message("    \u2713 ", basename(outfile))
}

# Taxa by pathway heatmap
#
#   rows   taxa at genus level, or whichever rank_idx selects
#   cols   cycle pathways, labelled in Spanish
#   colour white to red, mean normalised abundance per group
#
# pw_results: named list of lists of 6 data frames (output of ACF(), Fer(), etc.)
# rank_idx  index into each list: 1 phylum, 2 class, 3 order, 4 family, 5 genus, 6 species

make_cnps_taxonomy_heatmap <- function(pw_results, rank_idx = 5L,
                                       rank_label = "G\u00e9nero",
                                       title_es, outfile) {
  mat_cols <- lapply(names(pw_results), function(pw_name) {
    res <- pw_results[[pw_name]]
    if (is.null(res) || length(res) < rank_idx) {
      return(NULL)
    }
    df <- res[[rank_idx]]
    if (is.null(df) || !is.data.frame(df) || nrow(df) == 0) {
      return(NULL)
    }
    num_cols <- which(sapply(df, is.numeric))
    if (length(num_cols) == 0) {
      return(NULL)
    }
    vals <- rowMeans(df[, num_cols, drop = FALSE], na.rm = TRUE)
    setNames(as.data.frame(vals), pw_name)
  })

  keep_idx <- !sapply(mat_cols, is.null)
  if (sum(keep_idx) == 0) {
    message("    ", title_es, " (taxonom\u00eda): sin datos – omitido")
    return(invisible(NULL))
  }

  # One column per pathway, one row per taxon. The taxa come from the union of
  # the row names, since a pathway absent from a taxon leaves it out entirely.
  dfs <- lapply(seq_along(pw_results)[keep_idx], function(i) {
    res <- pw_results[[names(pw_results)[i]]]
    df <- res[[rank_idx]]
    num_cols <- which(sapply(df, is.numeric))
    vals <- rowMeans(df[, num_cols, drop = FALSE], na.rm = TRUE)
    data.frame(
      taxon = rownames(df), value = vals,
      pathway = names(pw_results)[i],
      stringsAsFactors = FALSE
    )
  })
  combined <- do.call(rbind, dfs)

  taxa_names <- sort(unique(combined$taxon))
  pw_names <- names(pw_results)[keep_idx]

  mat <- matrix(0,
    nrow = length(taxa_names), ncol = length(pw_names),
    dimnames = list(taxa_names, pw_names)
  )
  for (i in seq_len(nrow(combined))) {
    tx <- combined$taxon[i]
    pw <- combined$pathway[i]
    if (tx %in% taxa_names && pw %in% pw_names) {
      mat[tx, pw] <- combined$value[i]
    }
  }

  # Drop rows and columns that are entirely empty
  mat <- mat[rowSums(mat, na.rm = TRUE) > 0, , drop = FALSE]
  mat <- mat[, colSums(mat, na.rm = TRUE) > 0, drop = FALSE]
  if (nrow(mat) == 0 || ncol(mat) == 0) {
    message("    ", title_es, " (taxonom\u00eda): sin datos v\u00e1lidos – omitido")
    return(invisible(NULL))
  }

  # Normalise to 0-1 against the cycle maximum
  mx <- max(mat, na.rm = TRUE)
  if (mx > 0) mat <- mat / mx

  # Pathway labels into Spanish
  col_es <- capitalize_first(
    translate_to_spanish(colnames(mat), translation_cache_file = CNPS_CACHE_FILE)
  )
  colnames(mat) <- col_es

  col_fun <- colorRamp2(
    seq(0, 1, length.out = 6),
    c("white", "yellow", "orange", "orangered", "red", "darkred")
  )

  n_taxa <- nrow(mat)
  n_pw <- ncol(mat)
  row_fs <- if (n_taxa <= 15) 9 else 7

  ht <- Heatmap(
    mat,
    name = "Abundancia",
    col = col_fun,
    cluster_rows = TRUE,
    show_row_dend = TRUE,
    cluster_columns = FALSE,
    show_row_names = TRUE,
    row_names_side = "right",
    row_names_gp = gpar(fontsize = row_fs),
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 9),
    column_names_rot = 45,
    rect_gp = gpar(col = "black", lwd = 0.4),
    border = TRUE,
    column_title = paste0(title_es, " (por ", rank_label, ")"),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    heatmap_legend_param = list(
      title            = "Abundancia\nnormalizada",
      at               = c(0, 0.25, 0.5, 0.75, 1),
      labels           = c("0%", "25%", "50%", "75%", "100%"),
      title_gp         = gpar(fontsize = 10, fontface = "bold"),
      labels_gp        = gpar(fontsize = 9),
      border           = "black",
      legend_direction = "vertical"
    )
  )

  w <- max(8, round(n_pw * 0.55 + 5))
  h <- max(5, round(n_taxa * 0.35 + 3))
  png(outfile, width = w, height = h, units = "in", res = 300)
  draw(ht, merge_legend = TRUE)
  dev.off()
  message("    \u2713 ", basename(outfile))
}

# Every figure for one cycle, into a given directory
#
# Produces:
#   1. abun.heatmap.g by station and by depth, package logic
#   2. abun.heatmap.s, package logic
#   3. make_cnps_heatmap for each annotation version, DRAM2 style, DisplayName rows
#   4. make_cnps_station_heatmap, mean abundance per station
#      only if station_heatmap = TRUE
#
# abundance_df must already be filtered to the MAG subset, if one applies.
# ann_versions: which sidebar variants to produce.
#   0 no sidebar, 1 depth and station, 2 depth, 3 station
# station_heatmap: TRUE adds the per-station heatmap, which only applies to the full 54

run_cycle_plots <- function(abundance_df, cycle_name, cycle_title_es,
                            out_dir, group_st, group_dp,
                            ann_versions = c(1L, 2L, 3L),
                            station_heatmap = TRUE) {
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  # Write the abundance CSV
  write.csv(abundance_df, file.path(out_dir, paste0(cycle_name, "_abundance.csv")),
    row.names = FALSE
  )

  n_path <- nrow(abundance_df)
  h_grp <- max(4, round(n_path * 0.6 + 2))

  # Restrict the groups to the MAGs actually present in abundance_df
  mag_cols <- setdiff(colnames(abundance_df), "Type")
  grp_st <- group_st[group_st$ID %in% mag_cols, , drop = FALSE]
  grp_st$Group <- droplevels(grp_st$Group)
  grp_dp <- group_dp[group_dp$ID %in% mag_cols, , drop = FALSE]
  grp_dp$Group <- droplevels(grp_dp$Group)

  # Spanish names for the package figures, abun.heatmap.*
  abund_es <- abundance_df
  abund_es$Type <- capitalize_first(
    translate_to_spanish(as.character(abundance_df$Type),
      translation_cache_file = CNPS_CACHE_FILE
    )
  )

  # heatmap.s takes only the rows carrying at least one non-zero value
  num_mat <- abund_es[, mag_cols, drop = FALSE]
  has_signal <- rowSums(num_mat > 0, na.rm = TRUE) > 0
  abund_es_nz <- abund_es[has_signal, , drop = FALSE]
  h_nz <- max(3, round(sum(has_signal) * 0.6 + 2))

  # One white-to-red scale throughout, overriding the package default
  .cnps_fill_scale <- function() {
    scale_fill_gradientn(
      colours = c("white", "yellow", "orange", "orangered", "red", "darkred"),
      na.value = "grey85",
      name = "Abundancia\nnormalizada"
    )
  }

  # ---- 1a. abun.heatmap.g, grouped by station ----
  tryCatch(
    {
      res <- abun.heatmap.g(abund_es, grp_st, nlevels(grp_st$Group))
      write.csv(res[[1]], file.path(out_dir, paste0(cycle_name, "_stats_station.csv")),
        row.names = FALSE
      )
      out_f <- file.path(out_dir, paste0(cycle_name, "_heatmap_grouped_station.png"))
      p <- res[[2]] + .cnps_fill_scale()
      png(out_f, width = 11, height = h_grp, units = "in", res = 300)
      suppressMessages(print(p))
      dev.off()
      message("  ", cycle_name, " grouped/station OK")
    },
    error = function(e) {
      try(dev.off(), silent = TRUE)
      message("  ", cycle_name, " grouped/station omitido: ", e$message)
    }
  )

  # ---- 1b. abun.heatmap.g, grouped by depth ----
  tryCatch(
    {
      res <- abun.heatmap.g(abund_es, grp_dp, nlevels(grp_dp$Group))
      write.csv(res[[1]], file.path(out_dir, paste0(cycle_name, "_stats_depth.csv")),
        row.names = FALSE
      )
      out_f <- file.path(out_dir, paste0(cycle_name, "_heatmap_grouped_depth.png"))
      p <- res[[2]] + .cnps_fill_scale()
      png(out_f, width = 8, height = h_grp, units = "in", res = 300)
      suppressMessages(print(p))
      dev.off()
      message("  ", cycle_name, " grouped/depth OK")
    },
    error = function(e) {
      try(dev.off(), silent = TRUE)
      message("  ", cycle_name, " grouped/depth omitido: ", e$message)
    }
  )

  # ---- 2. abun.heatmap.s, per MAG, pathways with data only ----
  if (sum(has_signal) < 1L) {
    message("  ", cycle_name, " per-sample: sin pathways con datos, omitido")
  } else {
    tryCatch(
      {
        res <- abun.heatmap.s(
          abund_es_nz, grp_st,
          paste(cycle_title_es, "- MAGs individuales")
        )
        p <- res + .cnps_fill_scale()
        out_f <- file.path(out_dir, paste0(cycle_name, "_heatmap_samples.png"))
        n_mags <- length(mag_cols)
        png(out_f,
          width = max(10, round(n_mags * 0.22 + 4)),
          height = h_nz, units = "in", res = 300
        )
        suppressMessages(print(p))
        dev.off()
        message("  ", cycle_name, " per-sample OK")
      },
      error = function(e) {
        try(dev.off(), silent = TRUE)
        message("  ", cycle_name, " per-sample omitido: ", e$message)
      }
    )
  }

  # ---- 3. ComplexHeatmap, one per requested annotation version ----
  ann_suffixes <- c(
    "0" = "no_sidebar", "1" = "depth_station",
    "2" = "depth", "3" = "station"
  )
  for (av in ann_versions) {
    make_cnps_heatmap(
      abundance_df, cycle_title_es,
      file.path(out_dir, paste0(
        cycle_name, "_heatmap_MAGs_",
        ann_suffixes[as.character(av)], ".png"
      )),
      ann_v = as.integer(av)
    )
  }

  # ---- 4. per-station heatmap, full cohort only ----
  if (station_heatmap) {
    make_cnps_station_heatmap(
      abundance_df, cycle_title_es,
      file.path(out_dir, paste0(cycle_name, "_heatmap_by_station.png"))
    )
  }
}

# ---- SECTION A: ELEMENT CYCLES (54 MAGs) ----

message("\n══ [A] Ciclos elementales – 54 MAGs ══")

# ---- carbono ----
message(" Carbono ")
C.abundance <- Ccyc.abundance(kegg.entry = ko)
run_cycle_plots(
  C.abundance, "C", "Ciclo del Carbono",
  file.path(OUT_DIR, "carbon"),
  group_station, group_depth
)

# ---- nitrogen ----
message(" Nitr\u00f3geno ")
N.abundance <- Ncyc.abundance(kegg.entry = ko)
run_cycle_plots(
  N.abundance, "N", "Ciclo del Nitr\u00f3geno",
  file.path(OUT_DIR, "nitrogen"),
  group_station, group_depth
)

# ---- phosphorus ----
message(" F\u00f3sforo ")
P.abundance <- Pcyc.abundance(kegg.entry = ko)
p_num_cols <- setdiff(colnames(P.abundance), "Type")
if (any(P.abundance[, p_num_cols] > 0)) {
  run_cycle_plots(
    P.abundance, "P", "Ciclo del F\u00f3sforo",
    file.path(OUT_DIR, "phosphorus"),
    group_station, group_depth
  )
} else {
  write.csv(P.abundance, file.path(OUT_DIR, "phosphorus", "P_abundance.csv"), row.names = FALSE)
  message("  F\u00f3sforo: KOs de PhoR-PhoB/SenX3-RegX3 ausentes (esperado en archaea)")
}

# ---- azufre ----
message(" Azufre ")
S.abundance <- Scyc.abundance(kegg.entry = ko)
run_cycle_plots(
  S.abundance, "S", "Ciclo del Azufre",
  file.path(OUT_DIR, "sulfur"),
  group_station, group_depth
)

# ---- SECTION B: ORDINATION (PCA, PCoA, NMDS), 54 MAGs ----

message("\n══ [B] Ordenaci\u00f3n – 54 MAGs ══")

run_ordination <- function(abundance_df, group_df, cycle_name, out_dir) {
  sub_mat <- abundance_df[, 2:ncol(abundance_df)]
  tryCatch(
    {
      res <- pca.arg(sub_mat, group_df)
      write.csv(res[[1]], file.path(out_dir, paste0(cycle_name, "_PCA_coords.csv")))
      ggsave(file.path(out_dir, paste0(cycle_name, "_PCA.png")), res[[2]], width = 8, height = 6, dpi = 300)
      ggsave(file.path(out_dir, paste0(cycle_name, "_PCA_ell.png")), res[[3]], width = 8, height = 6, dpi = 300)
      ggsave(file.path(out_dir, paste0(cycle_name, "_PCA_lab.png")), res[[4]], width = 8, height = 6, dpi = 300)
      message("  PCA-", cycle_name, " OK")
    },
    error = function(e) message("  PCA-", cycle_name, " omitido: ", e$message)
  )
  tryCatch(
    {
      res <- pcoa.arg(sub_mat, group_df)
      ggsave(file.path(out_dir, paste0(cycle_name, "_PCoA.png")), res[[2]], width = 8, height = 6, dpi = 300)
      ggsave(file.path(out_dir, paste0(cycle_name, "_PCoA_ell.png")), res[[3]], width = 8, height = 6, dpi = 300)
      ggsave(file.path(out_dir, paste0(cycle_name, "_PCoA_lab.png")), res[[4]], width = 8, height = 6, dpi = 300)
      message("  PCoA-", cycle_name, " OK")
    },
    error = function(e) message("  PCoA-", cycle_name, " omitido: ", e$message)
  )
  tryCatch(
    {
      res <- nmds.arg(sub_mat, group_df)
      ggsave(file.path(out_dir, paste0(cycle_name, "_NMDS.png")), res[[2]], width = 8, height = 6, dpi = 300)
      ggsave(file.path(out_dir, paste0(cycle_name, "_NMDS_ell.png")), res[[3]], width = 8, height = 6, dpi = 300)
      ggsave(file.path(out_dir, paste0(cycle_name, "_NMDS_lab.png")), res[[4]], width = 8, height = 6, dpi = 300)
      message("  NMDS-", cycle_name, " OK")
    },
    error = function(e) message("  NMDS-", cycle_name, " omitido: ", e$message)
  )
}

ord_out <- file.path(OUT_DIR, "ordination")
run_ordination(C.abundance, group_station, "C_station", ord_out)
run_ordination(N.abundance, group_station, "N_station", ord_out)
run_ordination(P.abundance, group_station, "P_station", ord_out)
run_ordination(S.abundance, group_station, "S_station", ord_out)
run_ordination(C.abundance, group_depth, "C_depth", ord_out)
run_ordination(N.abundance, group_depth, "N_depth", ord_out)
run_ordination(S.abundance, group_depth, "S_depth", ord_out)

# ---- section C: host-relative statistics, taxonomy against cycle ----
# The package functions read Gene, tax, abundance and group from the global environment.

message("\n══ [C] Host-relative stats ══")

Gene <- Gene
tax <- tax
abundance <- abundance
group <- group_station

# Binary presence
C_host <- Ccyc.host(Gene)
names(C_host) <- c("ACF", "ACH4O", "AR", "AnCF", "COo", "Fer", "Meth")
write.csv(as.data.frame(t(C_host)), file.path(OUT_DIR, "carbon", "C_host_presence.csv"))
N_host <- Ncyc.host(Gene)
write.csv(as.data.frame(t(N_host)), file.path(OUT_DIR, "nitrogen", "N_host_presence.csv"))
P_host <- Pcyc.host(Gene)
write.csv(as.data.frame(t(P_host)), file.path(OUT_DIR, "phosphorus", "P_host_presence.csv"))
S_host <- Scyc.host(Gene)
write.csv(as.data.frame(t(S_host)), file.path(OUT_DIR, "sulfur", "S_host_presence.csv"))

# Individual pathways, carbon
c_pws <- list(ACF = ACF, Fer = Fer, Meth = Meth, ACH4O = ACH4O, COo = COo, AnCF = AnCF)
c_results <- list()
for (nm in names(c_pws)) {
  r <- tryCatch(
    c_pws[[nm]](kegg.category = Gene, NR.tax = tax, profile = abundance, group = group),
    error = function(e) NULL
  )
  if (!is.null(r)) {
    save_host_relative(r, "carbon/host_relative", nm)
    c_results[[nm]] <- r
    message("    C/", nm, " OK")
  } else {
    message("    C/", nm, " omitido")
  }
}

# Individual pathways, nitrogen
n_pws <- list(
  Nitrification_NH4_NH2OH = Nitrification.NH4.NH2OH,
  Nitrification_NH2OH_NO2 = Nitrification.NH2OH.NO2,
  Nitrification_NH2OH_NH4 = Nitrification.NH2OH.NH4,
  Denitrification_NO2_NO = Denitrification.NO2.NO,
  Denitrification_NO_N2O = Denitrification.NO.N2O,
  Denitrification_N2O_N2 = Denitrification.N2O.N2,
  ANR_NO2_NH4 = ANR.NO2.NH4, ANR_NO3_NO2 = ANR.NO3.NO2,
  DNR_NO2_NH4 = DNR.NO2.NH4, DNR_NO3_NO2 = DNR.NO3.NO2,
  Anammox_NH4_N2H4 = Anammox.NH4.N2H4, Anammox_N2H4_N2 = Anammox.N2H4.N2,
  NF_N2_NH4 = NF.N2.NH4, NU_NH4 = NU.NH4, NU_NO2 = NU.NO2, NU_NO3 = NU.NO3,
  NM_ON_NH4 = NM.ON.NH4, NA_NH4_ON = NA.NH4.ON
)
n_results <- list()
for (nm in names(n_pws)) {
  r <- tryCatch(
    n_pws[[nm]](kegg.category = Gene, NR.tax = tax, profile = abundance, group = group),
    error = function(e) NULL
  )
  if (!is.null(r)) {
    save_host_relative(r, "nitrogen/host_relative", nm)
    n_results[[nm]] <- r
    message("    N/", nm, " OK")
  } else {
    message("    N/", nm, " omitido")
  }
}

# Individual pathways, sulfur
s_pws <- list(
  DRS_SO4_APS = DRS.SO4.APS, DRS_APS_SO3 = DRS.APS.SO3, DRS_SO3_H2S = DRS.SO3.H2S,
  ARS_SO4_APS = ARS.SO4.APS, ARS_APS_PAPS = ARS.APS.PAPS, ARS_PAPS_SO3 = ARS.PAPS.SO3,
  ARS_SO3_H2S = ARS.SO3.H2S, SOX = SOX, SU_SO4 = SU.SO4, SU_SO3 = SU.SO3,
  SC_T_H2S = SC.T.H2S, SC_H2S_S = SC.H2S.S, SC_H2S_Sn = SC.H2S.Sn, SM = SM
)
s_results <- list()
for (nm in names(s_pws)) {
  r <- tryCatch(
    s_pws[[nm]](kegg.category = Gene, NR.tax = tax, profile = abundance, group = group),
    error = function(e) NULL
  )
  if (!is.null(r)) {
    save_host_relative(r, "sulfur/host_relative", nm)
    s_results[[nm]] <- r
    message("    S/", nm, " OK")
  } else {
    message("    S/", nm, " omitido")
  }
}

# Individual pathways, phosphorus
p_results <- list()
for (nm_fn in list(list("RB", RB), list("X3", X3))) {
  nm <- nm_fn[[1]]
  fn <- nm_fn[[2]]
  r <- tryCatch(
    fn(kegg.category = Gene, NR.tax = tax, profile = abundance, group = group),
    error = function(e) NULL
  )
  if (!is.null(r)) {
    save_host_relative(r, "phosphorus/host_relative", nm)
    p_results[[nm]] <- r
    message("    P/", nm, " OK")
  } else {
    message("    P/", nm, " omitido")
  }
}

# Taxa by pathway heatmaps at phylum, family and genus
message("  Generando heatmaps taxonom\u00eda \u00d7 pathways...")
for (rnk in list(
  list(1L, "Filo", "phylum"), list(4L, "Familia", "family"),
  list(5L, "G\u00e9nero", "genus")
)) {
  make_cnps_taxonomy_heatmap(
    c_results, rnk[[1]], rnk[[2]],
    "Ciclo del Carbono",
    file.path(OUT_DIR, "carbon", paste0("C_taxonomy_", rnk[[3]], ".png"))
  )
  make_cnps_taxonomy_heatmap(
    n_results, rnk[[1]], rnk[[2]],
    "Ciclo del Nitr\u00f3geno",
    file.path(OUT_DIR, "nitrogen", paste0("N_taxonomy_", rnk[[3]], ".png"))
  )
  make_cnps_taxonomy_heatmap(
    s_results, rnk[[1]], rnk[[2]],
    "Ciclo del Azufre",
    file.path(OUT_DIR, "sulfur", paste0("S_taxonomy_", rnk[[3]], ".png"))
  )
  if (length(p_results) > 0) {
    make_cnps_taxonomy_heatmap(
      p_results, rnk[[1]], rnk[[2]],
      "Ciclo del F\u00f3sforo",
      file.path(OUT_DIR, "phosphorus", paste0("P_taxonomy_", rnk[[3]], ".png"))
    )
  }
}

# ---- section D: host ratio, stacked bars by phylum ----
#
# Replaces the package's host.ratio(), which fixes cellwidth and cellheight in
# pheatmap and leaves most of the canvas blank when there are only a few taxa.
# This ggplot2 version sizes itself to whatever it is given.

message("\n══ [D] Host ratio plots ══")

run_host_ratio_result <- function(result, pw_label, outfile) {
  if (is.null(result) || length(result) == 0) {
    return(invisible(NULL))
  }
  pm <- result[[1]] # phylum level of the pathway function
  if (is.null(pm) || !is.data.frame(pm) || nrow(pm) == 0) {
    return(invisible(NULL))
  }

  num_cols <- names(which(sapply(pm, is.numeric)))
  if (length(num_cols) == 0 ||
    !any(pm[, num_cols, drop = FALSE] > 0, na.rm = TRUE)) {
    message("  host.ratio: ", pw_label, " sin genes – omitido")
    return(invisible(NULL))
  }

  tryCatch(
    {
      # Convert to long format
      df <- as.data.frame(pm)
      df$Taxon <- rownames(df)
      df_long <- reshape2::melt(df,
        id.vars = "Taxon",
        variable.name = "Grupo", value.name = "Valor"
      )
      df_long$Valor[is.na(df_long$Valor)] <- 0

      # Proportions within each group, that is within each station
      df_sum <- aggregate(Valor ~ Grupo, df_long, sum)
      names(df_sum)[2] <- "Total"
      df_long <- merge(df_long, df_sum, by = "Grupo")
      df_long$Prop <- ifelse(df_long$Total > 0,
        df_long$Valor / df_long$Total, 0
      )

      n_taxa <- length(unique(df_long$Taxon))
      pal <- if (n_taxa <= 9L) {
        RColorBrewer::brewer.pal(max(3L, n_taxa), "Set1")[seq_len(n_taxa)]
      } else {
        grDevices::rainbow(n_taxa)
      }
      taxa_sorted <- sort(unique(df_long$Taxon))
      pal <- setNames(pal, taxa_sorted)

      p <- ggplot(df_long, aes(x = Grupo, y = Prop, fill = Taxon)) +
        geom_bar(
          stat = "identity", position = "stack",
          color = "black", linewidth = 0.3
        ) +
        scale_fill_manual(values = pal) +
        scale_y_continuous(
          labels = function(x) paste0(round(x * 100), "%"),
          expand = c(0, 0), limits = c(0, 1.01)
        ) +
        labs(
          title = pw_label,
          x = "Estaci\u00f3n", y = "Proporci\u00f3n relativa", fill = "Filo"
        ) +
        theme_bw(base_size = 11, base_family = font_family) +
        theme(
          panel.grid.major.x = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(face = "bold", size = 12)
        )

      w <- max(4, 1.5 + length(unique(df_long$Grupo)) * 0.7 + n_taxa * 0.15)
      png(outfile, width = w, height = 4, units = "in", res = 300)
      print(p)
      dev.off()
      message("  host.ratio: ", pw_label, " OK")
    },
    error = function(e) {
      try(dev.off(), silent = TRUE)
      message("  host.ratio: ", pw_label, " omitido: ", e$message)
    }
  )
}

for (nm in names(c_results)) {
  run_host_ratio_result(
    c_results[[nm]], paste0("C/", nm, " (Filo)"),
    file.path(OUT_DIR, "carbon", paste0(nm, "_host_ratio.png"))
  )
}
for (nm in names(n_results)) {
  run_host_ratio_result(
    n_results[[nm]], paste0("N/", nm, " (Filo)"),
    file.path(OUT_DIR, "nitrogen", paste0(nm, "_host_ratio.png"))
  )
}
for (nm in names(s_results)) {
  run_host_ratio_result(
    s_results[[nm]], paste0("S/", nm, " (Filo)"),
    file.path(OUT_DIR, "sulfur", paste0(nm, "_host_ratio.png"))
  )
}
for (nm in names(p_results)) {
  run_host_ratio_result(
    p_results[[nm]], paste0("P/", nm, " (Filo)"),
    file.path(OUT_DIR, "phosphorus", paste0(nm, "_host_ratio.png"))
  )
}

# SECTION F: THE 10 MAGs WITHOUT A SPECIES
# Same structure as sections A and B, restricted to the mags_10 subset.
# Laid out like section 16 of DRAM2_output_analysis.R.

message("\n══ [F] Subconjunto 10 MAGs sin especie ══")
message("  MAGs: ", paste(mags_10, collapse = ", "))

if (length(mags_10) == 0) {
  message("  ADVERTENCIA: ningún MAG cumple el criterio de 'sin especie'.")
  message("  Revisa la columna Species en MAGs_metadata.tsv.")
} else {
  out_10 <- file.path(OUT_DIR, "10_MAGs")

  # Subset of the abundance tables
  C10 <- subset_abund(C.abundance, mags_10)
  N10 <- subset_abund(N.abundance, mags_10)
  P10 <- subset_abund(P.abundance, mags_10)
  S10 <- subset_abund(S.abundance, mags_10)

  # ---- carbon cycle ----
  message(" Carbono (10 MAGs) ")
  if (!is.null(C10)) {
    run_cycle_plots(
      C10, "C", "Ciclo del Carbono",
      file.path(out_10, "carbon"),
      group_10_station, group_10_depth,
      ann_versions = c(2L, 0L),
      station_heatmap = FALSE
    )
  }

  # ---- nitrogen cycle ----
  message(" Nitr\u00f3geno (10 MAGs) ")
  if (!is.null(N10)) {
    run_cycle_plots(
      N10, "N", "Ciclo del Nitr\u00f3geno",
      file.path(out_10, "nitrogen"),
      group_10_station, group_10_depth,
      ann_versions = c(2L, 0L),
      station_heatmap = FALSE
    )
  }

  # ---- phosphorus cycle ----
  message(" F\u00f3sforo (10 MAGs) ")
  if (!is.null(P10)) {
    p10_num <- setdiff(colnames(P10), "Type")
    if (any(P10[, p10_num] > 0)) {
      run_cycle_plots(
        P10, "P", "Ciclo del F\u00f3sforo",
        file.path(out_10, "phosphorus"),
        group_10_station, group_10_depth,
        ann_versions = c(2L, 0L),
        station_heatmap = FALSE
      )
    } else {
      message("  F\u00f3sforo (10 MAGs): sin datos")
    }
  }

  # ---- sulfur cycle ----
  message(" Azufre (10 MAGs) ")
  if (!is.null(S10)) {
    run_cycle_plots(
      S10, "S", "Ciclo del Azufre",
      file.path(out_10, "sulfur"),
      group_10_station, group_10_depth,
      ann_versions = c(2L, 0L),
      station_heatmap = FALSE
    )
  }

  # ---- ordination, 10 MAGs ----
  message(" Ordenaci\u00f3n (10 MAGs) ")
  ord_10 <- file.path(out_10, "ordination")
  if (!is.null(C10)) {
    run_ordination(C10, group_10_station, "C_station", ord_10)
    run_ordination(C10, group_10_depth, "C_depth", ord_10)
  }
  if (!is.null(N10)) {
    run_ordination(N10, group_10_station, "N_station", ord_10)
    run_ordination(N10, group_10_depth, "N_depth", ord_10)
  }
  if (!is.null(S10)) {
    run_ordination(S10, group_10_station, "S_station", ord_10)
    run_ordination(S10, group_10_depth, "S_depth", ord_10)
  }
}

# SECTION G: COMMUNITY-LEVEL SUMMARY, TAXONOMY AND METADATA BUBBLES
#
# Absorbs what used to live in CNPS.cycle_2.R.
# Reuses: translate_to_spanish, capitalize_first, font_family,
#            WSL_BASE, OUT_DIR, CNPS_CACHE_FILE, meta
# Reads the CSVs written in section A rather than the in-memory objects.
# The _g suffix keeps these names clear of the global environment.

message("\n\u2550\u2550 [G] Resumen comunidad: gr\u00e1ficos de burbujas taxon\u00f3micas y de metadatos \u2550\u2550")

# ---- g.0  Output directory ----

dir.create(file.path(OUT_DIR, "resumen_comunidad"), recursive = TRUE, showWarnings = FALSE)

save_fig_g <- function(plot, filename, width = 10, height = 7) {
  out_path <- file.path(OUT_DIR, "resumen_comunidad", filename)
  png(out_path, width = width, height = height, units = "in", res = 300)
  print(plot)
  dev.off()
  cat("\u2713 Saved:", filename, "\n")
}

# ---- g.1 publication theme ----

theme_publication_g <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor       = element_blank(),
      panel.grid.major.x     = element_blank(),
      panel.grid.major.y     = element_line(color = "gray", linetype = "dotted"),
      axis.text              = element_text(color = "black"),
      axis.title             = element_text(face = "bold"),
      legend.position        = "right",
      legend.background      = element_rect(fill = "white", color = NA),
      legend.title           = element_text(face = "bold", size = 12),
      legend.text            = element_text(size = 11),
      strip.background       = element_rect(fill = "gray90", color = "black"),
      strip.text             = element_text(face = "bold", size = 13),
      plot.title             = element_text(face = "bold", hjust = 0.5),
      text                   = element_text(family = font_family)
    )
}

# ---- g.2 build df_long_all_g from the CSVs written in section A ----
# Reading the files rather than the in-memory objects is deliberate. Run across
# two RStudio sessions, the objects still in memory may be left over from an
# earlier run and disagree with what is on disk.

.read_abund_csv <- function(csv_path, abund_fallback) {
  if (file.exists(csv_path)) {
    # The CSVs hold the values already multiplied by 100, on a 0 to 100 scale,
    # exactly as Ccyc.abundance, Ncyc.abundance and Scyc.abundance produce them
    # from the scaled `ko`. They match the in-memory objects and need no
    # further transformation.
    return(read.csv(csv_path, check.names = FALSE, stringsAsFactors = FALSE))
  }
  message("  [G] Advertencia: CSV no encontrado, usando objeto en memoria: ", csv_path)
  abund_fallback
}

.make_long_g <- function(abund_df, nombre_ciclo) {
  if (is.null(abund_df)) {
    return(NULL)
  }
  num_cols <- setdiff(colnames(abund_df), "Type")
  if (!any(abund_df[, num_cols] > 0, na.rm = TRUE)) {
    return(NULL)
  }
  abund_df |>
    tidyr::pivot_longer(cols = -Type, names_to = "MAG", values_to = "Abundancia") |>
    dplyr::mutate(
      Ciclo = nombre_ciclo,
      Pathway_ES = capitalize_first(
        translate_to_spanish(Type, CNPS_CACHE_FILE)
      )
    )
}

C.abund_g <- .read_abund_csv(file.path(OUT_DIR, "carbon", "C_abundance.csv"), C.abundance)
N.abund_g <- .read_abund_csv(file.path(OUT_DIR, "nitrogen", "N_abundance.csv"), N.abundance)
P.abund_g <- .read_abund_csv(file.path(OUT_DIR, "phosphorus", "P_abundance.csv"), P.abundance)
S.abund_g <- .read_abund_csv(file.path(OUT_DIR, "sulfur", "S_abundance.csv"), S.abundance)

df_long_c_g <- .make_long_g(C.abund_g, "Carbono")
df_long_n_g <- .make_long_g(N.abund_g, "Nitr\u00f3geno")
df_long_p_g <- .make_long_g(P.abund_g, "F\u00f3sforo")
df_long_s_g <- .make_long_g(S.abund_g, "Azufre")

df_long_all_g <- dplyr::bind_rows(df_long_c_g, df_long_n_g, df_long_p_g, df_long_s_g)

# ---- g.3 taxonomic name formatting ----

# CDC scientific nomenclature (https://wwwnc.cdc.gov/eid/page/scientific-nomenclature):
# only genus and species (and infraspecific ranks) are italicised; ranks above
# genus — family, order, class, phylum, domain — remain plain text.
ITALIC_LEVELS_G <- c("Genus", "Species", "Subspecies", "Strain")

should_italicize_level_g <- function(level) {
  !is.null(level) && length(level) == 1 && !is.na(level) && level %in% ITALIC_LEVELS_G
}

taxon_italic_g <- function(x, level = "Species") {
  if (!should_italicize_level_g(level)) {
    return(ifelse(is.na(x), x, as.character(x)))
  }
  sapply(x, function(nm) {
    if (is.na(nm) || nm == "Sin asignar") {
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

parent_hierarchy_g <- list(
  Species = "Genus",  Genus  = "Family",
  Family  = "Order",  Order  = "Class",
  Class   = "Phylum", Phylum = "Domain"
)
level_abbrev_g <- list(
  Genus = "g.", Family = "f.", Order = "o.",
  Class = "c.", Phylum = "p.", Domain = "d."
)

get_parent_lookup_g <- function(df_in, child_col) {
  ancestors <- character(0)
  cur <- child_col
  repeat {
    par <- parent_hierarchy_g[[cur]]
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
          abbrev = if (!is.null(level_abbrev_g[[anc_col]])) level_abbrev_g[[anc_col]] else ""
        )
        break
      }
    }
  }
  result
}

fmt_taxon_g <- function(nm, parent_lookup = NULL, level = "Species") {
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
        taxon_italic_g(nm, level),
        " (", abbrev_str, taxon_italic_g(info$name, info$level), ")"
      ))
    }
  }
  taxon_italic_g(nm, level)
}

fmt_taxon_vec_g <- function(x, parent_lookup = NULL, level = "Species") {
  sapply(x, fmt_taxon_g,
    parent_lookup = parent_lookup, level = level, USE.NAMES = FALSE
  )
}

make_label_fn_g <- function(parent_lookup, level = "Species") {
  function(x) fmt_taxon_vec_g(x, parent_lookup, level = level)
}

# ---- g.4 load the GTDB taxonomy ----

input_gtdb_summary_g <- file.path(
  WSL_BASE, "results/03_taxonomy_assigned/classify/tax_.ar53.summary.tsv"
)
if (!file.exists(input_gtdb_summary_g)) {
  stop("\u26a0 [G] No se encontr\u00f3 la taxonom\u00eda GTDB: ", input_gtdb_summary_g)
}

df_gtdb_g <- read.delim(input_gtdb_summary_g, stringsAsFactors = FALSE)

df_tax_g <- df_gtdb_g |>
  dplyr::select(user_genome, classification) |>
  dplyr::mutate(MAG = gsub("_", ".", user_genome)) |>
  tidyr::separate(
    classification,
    into = c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"),
    sep = ";", fill = "right"
  ) |>
  dplyr::mutate(
    Phylum  = gsub("p__", "", Phylum),
    Class   = gsub("c__", "", Class),
    Order   = gsub("o__", "", Order),
    Family  = gsub("f__", "", Family),
    Genus   = gsub("g__", "", Genus),
    Species = gsub("s__", "", Species)
  ) |>
  dplyr::select(MAG, Domain, Phylum, Class, Order, Family, Genus, Species)

df_integrado_g <- dplyr::inner_join(df_long_all_g, df_tax_g, by = "MAG")

# ---- g.5 bubble charts per taxonomic rank ----

plot_taxonomic_bubble_g <- function(data, tax_level, filename_suffix) {
  df_summary <- data |>
    dplyr::rename(Taxon = !!dplyr::sym(tax_level)) |>
    dplyr::filter(!is.na(Taxon) & Taxon != "") |>
    dplyr::group_by(Ciclo, Pathway_ES, Taxon) |>
    dplyr::summarise(
      n_mags_in_taxon = dplyr::n_distinct(MAG),
      Prevalencia = sum(Abundancia > 0) / n_mags_in_taxon * 100,
      Abundancia_Media = mean(Abundancia[Abundancia > 0], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::filter(Prevalencia > 0) |>
    dplyr::arrange(Ciclo, dplyr::desc(Prevalencia)) |>
    dplyr::mutate(
      Pathway_ES_uniq = factor(
        paste(Pathway_ES, Ciclo, sep = "_"),
        levels = unique(paste(Pathway_ES, Ciclo, sep = "_"))
      )
    )

  n_taxa <- length(unique(df_summary$Taxon))
  n_pathways <- length(unique(df_summary$Pathway_ES_uniq))
  lookup <- get_parent_lookup_g(data, tax_level)
  lbl_fn <- make_label_fn_g(lookup, level = tax_level)

  p <- ggplot(df_summary, aes(x = Taxon, y = Pathway_ES_uniq)) +
    geom_point(
      aes(size = Prevalencia, fill = Abundancia_Media),
      shape = 21, color = "gray20", stroke = 0.5
    ) +
    scale_size_continuous(
      range = c(1.5, 8), name = "Prevalencia\nen tax\u00f3n (%)", limits = c(0, 100)
    ) +
    scale_fill_viridis_c(
      option = "plasma", name = "Abundancia\nmedia", na.value = "transparent"
    ) +
    scale_x_discrete(labels = lbl_fn, expand = expansion(add = 0.4)) +
    scale_y_discrete(labels = function(x) {
      lab <- sub("_(Carbono|Nitr\u00f3geno|Azufre|F\u00f3sforo)$", "", x)
      lab <- sub("Fijaci\u00f3n anaer\u00f3bica de C\\b", "Fijaci\u00f3n anaer\u00f3bica de carbono", lab)
      lab <- sub("Mineralizaci\u00f3n de S\\b", "Mineralizaci\u00f3n de azufre", lab)
      lab <- gsub("N-org\u00e1nico", "nitr\u00f3geno org\u00e1nico", lab)
      lab <- gsub(",", ":", lab)
      lab
    }) +
    facet_grid(Ciclo ~ ., scales = "free_y", space = "free_y") +
    labs(x = NULL, y = NULL) +
    theme_publication_g(base_size = 13) +
    theme(
      axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, face = "plain"),
      axis.text.y = element_text(size = 10),
      legend.position = "right",
      legend.box = "vertical"
    )

  h_din <- max(6, n_pathways * 0.25 + 2)
  # No large floor: 3–4 taxa would hit max(10,…)=10 and look too sparse
  w_din <- max(6, n_taxa * 0.9 + 4)

  save_fig_g(p, paste0("Resumen_Taxonomia_", filename_suffix, ".png"),
    width = w_din, height = h_din
  )
}

message("  Generando gr\u00e1ficos taxon\u00f3micos de comunidad...")
plot_taxonomic_bubble_g(df_integrado_g, "Phylum", "Filo")
plot_taxonomic_bubble_g(df_integrado_g, "Class", "Clase")
plot_taxonomic_bubble_g(df_integrado_g, "Order", "Orden")
plot_taxonomic_bubble_g(df_integrado_g, "Family", "Familia")
plot_taxonomic_bubble_g(df_integrado_g, "Genus", "G\u00e9nero")

# ---- g.6 metadata preprocessing ----

df_metadata_g <- meta |>
  dplyr::select(MAG = NameWithoutExtension, Estacion = Station, Profundidad = Depth) |>
  dplyr::mutate(
    Profundidad = dplyr::recode(Profundidad,
      "ROMZ" = "ZMOR", "1000" = "1000 m", "BTM" = "Fondo"
    ),
    Profundidad = factor(Profundidad, levels = c("ZMOR", "1000 m", "Fondo")),
    Estacion = factor(Estacion)
  )

df_integrado_meta_g <- dplyr::inner_join(df_long_all_g, df_metadata_g, by = "MAG")

df_meta_long_g <- df_integrado_meta_g |>
  tidyr::pivot_longer(
    cols = c(Estacion, Profundidad),
    names_to = "Variable", values_to = "Grupo"
  ) |>
  dplyr::mutate(
    Variable = factor(Variable,
      levels = c("Estacion", "Profundidad"),
      labels = c("Estaci\u00f3n de muestreo", "Profundidad")
    )
  )

df_summary_meta_g <- df_meta_long_g |>
  dplyr::filter(!is.na(Grupo) & Grupo != "") |>
  dplyr::group_by(Ciclo, Pathway_ES, Variable, Grupo) |>
  dplyr::summarise(
    n_mags_in_group = dplyr::n_distinct(MAG),
    Prevalencia = sum(Abundancia > 0) / n_mags_in_group * 100,
    Abundancia_Media = mean(Abundancia[Abundancia > 0], na.rm = TRUE),
    .groups = "drop"
  ) |>
  dplyr::filter(Prevalencia > 0) |>
  dplyr::arrange(Ciclo, dplyr::desc(Prevalencia)) |>
  dplyr::mutate(
    Pathway_ES_uniq = factor(
      paste(Pathway_ES, Ciclo, sep = "_"),
      levels = unique(paste(Pathway_ES, Ciclo, sep = "_"))
    ),
    # pivot_longer drops the factor levels on Profundidad, leaving Grupo as
    # character, and ggplot2 would then order it alphabetically: 1000 m, Fondo,
    # ZMOR. Refactoring inside the mutate restores the depth order, the same way
    # tax_analysis.R does: stations first (alphabetical), then ZMOR, 1000 m, Fondo.
    Grupo = factor(Grupo, levels = c(
      sort(levels(df_metadata_g$Estacion)),
      DEPTH_DISPLAY_LEVELS
    ))
  )

# ---- g.7 combined metadata figure, station and depth ----

message("  Generando gr\u00e1fico integrado de metadatos (Estaci\u00f3n + Profundidad)...")

n_grupos_g <- length(unique(df_summary_meta_g$Grupo))
h_meta_g <- max(6, length(unique(df_summary_meta_g$Pathway_ES_uniq)) * 0.25 + 2)
# Reduced: was max(9, n_grupos * 1.2 + 5) = ~15.8 in — too wide for 9 groups
w_meta_g <- max(7, n_grupos_g * 0.75 + 3)

p_meta_g <- ggplot(df_summary_meta_g, aes(x = Grupo, y = Pathway_ES_uniq)) +
  geom_point(
    aes(size = Prevalencia, fill = Abundancia_Media),
    shape = 21, color = "gray20", stroke = 0.5
  ) +
  scale_size_continuous(
    range = c(2, 9), name = "Prevalencia\nen MAGs (%)", limits = c(0, 100)
  ) +
  scale_fill_viridis_c(
    option = "mako", name = "Abundancia\nmedia", na.value = "transparent"
  ) +
  scale_x_discrete(expand = expansion(add = 0.4)) +
  scale_y_discrete(labels = function(x) {
    lab <- sub("_(Carbono|Nitr\u00f3geno|Azufre|F\u00f3sforo)$", "", x)
    lab <- sub("Fijaci\u00f3n anaer\u00f3bica de C\\b", "Fijaci\u00f3n anaer\u00f3bica de carbono", lab)
    lab <- sub("Mineralizaci\u00f3n de S\\b", "Mineralizaci\u00f3n de azufre", lab)
    lab <- gsub("N-org\u00e1nico", "nitr\u00f3geno org\u00e1nico", lab)
    lab <- gsub(",", ":", lab)
    lab
  }) +
  facet_grid(Ciclo ~ Variable, scales = "free", space = "free") +
  labs(x = NULL, y = NULL) +
  theme_publication_g(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, face = "plain", size = 11),
    axis.text.y = element_text(size = 10),
    legend.position = "right",
    legend.box = "vertical",
    strip.text.x = element_text(face = "bold", size = 12)
  )

save_fig_g(p_meta_g, "Resumen_Metadatos_Combinado.png",
  width = w_meta_g, height = h_meta_g
)

# ---- g.8 taxonomy against metadata ----

message("  Generando gr\u00e1ficos combinados Taxonom\u00eda \u00d7 Metadatos...")

df_master_g <- dplyr::inner_join(
  df_integrado_g,
  dplyr::select(df_metadata_g, MAG, Estacion, Profundidad),
  by = "MAG"
)

plot_tax_meta_combined_g <- function(data, tax_level, meta_col,
                                     file_prefix, orientation = "tax_x") {
  df_summary <- data |>
    dplyr::rename(Taxon = !!dplyr::sym(tax_level), Meta = !!dplyr::sym(meta_col)) |>
    dplyr::filter(!is.na(Taxon) & Taxon != "" & !is.na(Meta) & Meta != "") |>
    dplyr::group_by(Ciclo, Pathway_ES, Taxon, Meta) |>
    dplyr::summarise(
      n_mags_in_group = dplyr::n_distinct(MAG),
      Prevalencia = sum(Abundancia > 0) / n_mags_in_group * 100,
      Abundancia_Media = mean(Abundancia[Abundancia > 0], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::filter(Prevalencia > 0) |>
    dplyr::arrange(Ciclo, dplyr::desc(Prevalencia)) |>
    dplyr::mutate(
      Pathway_ES_uniq = factor(
        paste(Pathway_ES, Ciclo, sep = "_"),
        levels = unique(paste(Pathway_ES, Ciclo, sep = "_"))
      )
    )

  lookup <- get_parent_lookup_g(data, tax_level)
  tax_levels_ord <- unique(df_summary$Taxon)
  tax_levels_fmt <- fmt_taxon_vec_g(tax_levels_ord, lookup, level = tax_level)

  df_summary <- df_summary |>
    dplyr::mutate(
      Taxon_Fmt = fmt_taxon_vec_g(Taxon, lookup, level = tax_level),
      Taxon_Fmt = factor(Taxon_Fmt, levels = tax_levels_fmt)
    )

  # summarise drops the factor levels, so the depth order is restored here the
  # same way tax_analysis.R does it with depth_factor_levels.
  if (meta_col == "Profundidad") {
    df_summary <- df_summary |>
      dplyr::mutate(Meta = factor(Meta, levels = DEPTH_DISPLAY_LEVELS))
  }

  bubble_scale <- scale_size_continuous(
    range = c(1.5, 6.5), name = "Prevalencia\nen MAGs (%)", limits = c(0, 100)
  )

  h_din <- max(6, length(unique(df_summary$Pathway_ES_uniq)) * 0.25 + 2)
  n_taxa <- length(unique(df_summary$Taxon))

  if (orientation == "tax_x") {
    # ---- taxon on the x axis, metadata as column facets ----
    p <- ggplot(df_summary, aes(x = Taxon_Fmt, y = Pathway_ES_uniq)) +
      geom_point(
        aes(size = Prevalencia, fill = Abundancia_Media),
        shape = 21, color = "gray20", stroke = 0.5
      ) +
      bubble_scale +
      scale_fill_viridis_c(
        option = "plasma", name = "Abundancia\nmedia", na.value = "transparent"
      ) +
      scale_x_discrete(expand = expansion(add = 0.7)) +
      scale_y_discrete(labels = function(x) {
      lab <- sub("_(Carbono|Nitr\u00f3geno|Azufre|F\u00f3sforo)$", "", x)
      lab <- sub("Fijaci\u00f3n anaer\u00f3bica de C\\b", "Fijaci\u00f3n anaer\u00f3bica de carbono", lab)
      lab <- sub("Mineralizaci\u00f3n de S\\b", "Mineralizaci\u00f3n de azufre", lab)
      lab <- gsub("N-org\u00e1nico", "nitr\u00f3geno org\u00e1nico", lab)
      lab <- gsub(",", ":", lab)
      lab
    }) +
      facet_grid(Ciclo ~ Meta, scales = "free", space = "free") +
      labs(x = NULL, y = NULL) +
      theme_publication_g(base_size = 13) +
      theme(
        axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, face = "plain"),
        axis.text.y = element_text(size = 10),
        strip.text.x = element_text(face = "plain", size = 12),
        legend.position = "right",
        legend.box = "vertical",
        panel.spacing.x = unit(0.9, "lines")
      )
    # Width must account for ALL meta-level facets × n_taxa items each
    n_meta_vals_tx <- length(unique(df_summary$Meta))
    w_din <- max(8, n_meta_vals_tx * (n_taxa * 0.55 + 1.0) + 4)
  } else {
    # ---- metadata on the x axis, taxon as column facets ----
    p <- ggplot(df_summary, aes(x = Meta, y = Pathway_ES_uniq)) +
      geom_point(
        aes(size = Prevalencia, fill = Abundancia_Media),
        shape = 21, color = "gray20", stroke = 0.5
      ) +
      bubble_scale +
      scale_fill_viridis_c(
        option = "plasma", name = "Abundancia\nmedia", na.value = "transparent"
      ) +
      # Fixed: was expansion(add = 1.8) — excessive blank space at sides
      scale_x_discrete(expand = expansion(add = 0.5)) +
      scale_y_discrete(labels = function(x) {
      lab <- sub("_(Carbono|Nitr\u00f3geno|Azufre|F\u00f3sforo)$", "", x)
      lab <- sub("Fijaci\u00f3n anaer\u00f3bica de C\\b", "Fijaci\u00f3n anaer\u00f3bica de carbono", lab)
      lab <- sub("Mineralizaci\u00f3n de S\\b", "Mineralizaci\u00f3n de azufre", lab)
      lab <- gsub("N-org\u00e1nico", "nitr\u00f3geno org\u00e1nico", lab)
      lab <- gsub(",", ":", lab)
      lab
    }) +
      # space = "free_y": all column panels equal width (was "free" to unequal)
      facet_grid(Ciclo ~ Taxon_Fmt, scales = "free_y", space = "free_y") +
      labs(x = NULL, y = NULL) +
      theme_publication_g(base_size = 13) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, face = "plain", size = 11),
        axis.text.y = element_text(size = 10),
        strip.text.x = ggtext::element_markdown(
          face = "plain", size = 11, angle = 0, hjust = 0.5
        ),
        strip.background.x = element_blank(),
        strip.clip = "off",
        legend.position = "right",
        legend.box = "vertical",
        panel.spacing.x = unit(0.3, "lines")
      )
    # Width: n_taxa equal-width panels, each wide enough for bubbles + strip text.
    # panel_w = n_meta_vals * 0.45 in (bubbles) clamped ≥ 1.8 in (strip text clearance)
    n_meta_vals <- length(unique(df_summary$Meta))
    panel_w <- max(n_meta_vals * 0.45 + 0.5, 1.8)
    w_din <- max(10, n_taxa * panel_w + 4)
  }

  save_fig_g(p, paste0(file_prefix, ".png"), width = w_din, height = h_din)
}

# ---- g.9 generate every variant ----

tax_levels_g <- c("Phylum", "Class", "Order", "Family", "Genus")
tax_names_es_g <- c("Filo", "Clase", "Orden", "Familia", "Genero")

for (i in seq_along(tax_levels_g)) {
  t_level <- tax_levels_g[i]
  t_name <- tax_names_es_g[i]

  # ---- A: taxa on top, station on the x axis ----
  plot_tax_meta_combined_g(
    df_master_g, t_level, "Estacion",
    paste0("Resumen_Comb_TaxonArriba_Estacion_", t_name),
    orientation = "meta_x"
  )

  # ---- B: taxa on top, depth on the x axis ----
  plot_tax_meta_combined_g(
    df_master_g, t_level, "Profundidad",
    paste0("Resumen_Comb_TaxonArriba_Profundidad_", t_name),
    orientation = "meta_x"
  )

  # ---- C: metadata on top, station as column facets ----
  plot_tax_meta_combined_g(
    df_master_g, t_level, "Estacion",
    paste0("Resumen_Comb_MetaArriba_Estacion_", t_name),
    orientation = "tax_x"
  )

  # ---- D: left out on purpose ----
  # Putting the depths (ZMOR, 1000 m, Fondo) across the top as column headers
  # gives labels far wider than the few bubbles each panel would hold.
  # plot_tax_meta_combined_g(
  #   df_master_g, t_level, "Profundidad",
  #   paste0("Resumen_Comb_MetaArriba_Profundidad_", t_name),
  #   orientation = "tax_x"
  # )
}

# ---- g.10 metadata and taxonomy in a single figure ----
#
# One ggplot, faceted by facet_grid(Ciclo ~ Panel, scales = "free_x",
# space = "free_x"), where Panel is station, depth or taxa.
#   drop = TRUE keeps each cycle to its own pathways
#   the cycle labels sit to the right of the taxa panel
#   "Taxones" appears as the top strip of that right-hand panel

message("  Generando paneles combinados Meta + Taxonom\u00eda...")

# Canonical y levels, computed in g.7
y_levels_shared_g <- levels(df_summary_meta_g$Pathway_ES_uniq)

make_combined_meta_tax_g <- function(tax_level, tax_suffix) {
  # ---- taxonomy summary ----
  df_tax_sum <- df_integrado_g |>
    dplyr::rename(Taxon = !!dplyr::sym(tax_level)) |>
    dplyr::filter(!is.na(Taxon) & Taxon != "") |>
    dplyr::group_by(Ciclo, Pathway_ES, Taxon) |>
    dplyr::summarise(
      n_mags = dplyr::n_distinct(MAG),
      Prevalencia = sum(Abundancia > 0) / n_mags * 100,
      Abundancia_Media = mean(Abundancia[Abundancia > 0], na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::filter(Prevalencia > 0) |>
    dplyr::mutate(
      Pathway_ES_uniq = factor(
        paste(Pathway_ES, Ciclo, sep = "_"),
        levels = y_levels_shared_g
      )
    )

  lookup <- get_parent_lookup_g(df_integrado_g, tax_level)
  n_taxa <- length(unique(df_tax_sum$Taxon))

  # ---- build the combined data frame ----
  # Left panel: metadata (x = group, Panel = variable name)
  df_meta_part <- df_summary_meta_g |>
    dplyr::mutate(
      Panel   = as.character(Variable),
      x_label = as.character(Grupo)
    ) |>
    dplyr::select(
      Ciclo, Pathway_ES_uniq, Panel, x_label,
      Prevalencia, Abundancia_Media
    )

  # Right-hand panel: taxa, x is the formatted name
  df_tax_part <- df_tax_sum |>
    dplyr::mutate(
      Panel   = "Taxones",
      x_label = fmt_taxon_vec_g(Taxon, lookup, level = tax_level)
    ) |>
    dplyr::select(
      Ciclo, Pathway_ES_uniq, Panel, x_label,
      Prevalencia, Abundancia_Media
    )

  # Metadata panels on the left, taxa on the right
  panel_levels <- c(levels(df_summary_meta_g$Variable), "Taxones")

  df_combined <- dplyr::bind_rows(df_meta_part, df_tax_part) |>
    dplyr::mutate(
      Pathway_ES_uniq = factor(Pathway_ES_uniq, levels = y_levels_shared_g),
      Panel           = factor(Panel, levels = panel_levels),
      # as.character(Grupo) in df_meta_part dropped the factor levels, so
      # x_label is refactored on the canonical order. Otherwise the depth panel
      # comes out alphabetical instead of ZMOR, 1000 m, Fondo.
      x_label = factor(x_label, levels = c(
        sort(levels(df_metadata_g$Estacion)),
        DEPTH_DISPLAY_LEVELS,
        unique(df_tax_part$x_label)
      ))
    )

  # ---- dimensions ----
  n_est <- length(unique(
    df_meta_part$x_label[df_meta_part$Panel == panel_levels[1]]
  ))
  n_prof <- length(unique(
    df_meta_part$x_label[df_meta_part$Panel == panel_levels[2]]
  ))
  w_total <- (n_est + n_prof + n_taxa) * 0.65 + 5
  h_comb <- max(7, length(y_levels_shared_g) * 0.22 + 2)

  # ---- the figure ----
  p <- ggplot(df_combined, aes(x = x_label, y = Pathway_ES_uniq)) +
    geom_point(
      aes(size = Prevalencia, fill = Abundancia_Media),
      shape = 21, color = "gray20", stroke = 0.5, na.rm = TRUE
    ) +
    scale_size_continuous(
      range = c(2, 8), name = "Prevalencia\nen MAGs (%)", limits = c(0, 100)
    ) +
    scale_fill_viridis_c(
      option = "plasma", name = "Abundancia\nmedia", na.value = "transparent"
    ) +
    scale_x_discrete(expand = expansion(add = 0.4)) +
    scale_y_discrete(labels = function(x) {
      lab <- sub("_(Carbono|Nitr\u00f3geno|Azufre|F\u00f3sforo)$", "", x)
      lab <- sub("Fijaci\u00f3n anaer\u00f3bica de C\\b", "Fijaci\u00f3n anaer\u00f3bica de carbono", lab)
      lab <- sub("Mineralizaci\u00f3n de S\\b", "Mineralizaci\u00f3n de azufre", lab)
      lab <- gsub("N-org\u00e1nico", "nitr\u00f3geno org\u00e1nico", lab)
      lab <- gsub(",", ":", lab)
      lab
    }) +
    # scales = "free": x free per column, y free per row, so pathways do not repeat
    facet_grid(Ciclo ~ Panel, scales = "free", space = "free") +
    labs(x = NULL, y = NULL) +
    theme_publication_g(base_size = 13) +
    theme(
      # x axis renders markdown, which is what italicises the taxon names
      axis.text.x = ggtext::element_markdown(angle = 45, hjust = 1, size = 11),
      axis.text.y = element_text(size = 10),
      # Top strip: station, depth or taxa
      strip.text.x = element_text(face = "bold", size = 12),
      # Side strip, the cycle name, to the right of the taxa panel and vertical
      strip.text.y = element_text(face = "bold", size = 12),
      legend.position = "right",
      legend.box = "vertical"
    )

  save_fig_g(p,
    paste0("Panel_Meta_Tax_", tax_suffix, ".png"),
    width = w_total, height = h_comb
  )
}

make_combined_meta_tax_g("Phylum", "Filo")
make_combined_meta_tax_g("Class", "Clase")
make_combined_meta_tax_g("Order", "Orden")

message(
  "\u00a1Secci\u00f3n G completada! Gr\u00e1ficos en: ",
  file.path(OUT_DIR, "resumen_comunidad")
)

# ---- FIN ----

message(paste0(
  "\n---\nAn\u00e1lisis CNPS.cycle completado.\nResultados en: ", OUT_DIR
))

cat(sprintf("\n=== Finalización del script: %s ===\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf(
  "    Duración total: %.1f minutos\n",
  as.numeric(difftime(Sys.time(), script_start, units = "mins"))
))
