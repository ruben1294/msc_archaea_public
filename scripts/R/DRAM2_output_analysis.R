#!/usr/bin/env Rscript

cat("Inicio:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n\n")

# Everything the DRAM2 output feeds: annotation-rate tables, trait heatmaps,
# metabolism summaries, the KEGG-Decoder comparison, and the figures for the
# functional annotation objective.
#
# Reads the distillate from pipeline step 08 and writes into
# results/SO_01/funct_annot/. Figure labels are Spanish, matching the thesis;
# translation runs through DeepL against a cached table, so the script works
# offline once the cache exists.

user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "") {
  user_lib <- file.path(
    Sys.getenv("HOME"), "R", "x86_64-pc-linux-gnu-library",
    paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = ".")
  )
}

if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}

.libPaths(c(user_lib, .libPaths()))
fecha <- format(Sys.Date(), "%Y%m%d")

# === VERIFICACIÓN DE PAQUETES ===
cat("Verificando compatibilidad de paquetes...\n")

check_and_update_packages <- function() {
  needs_update <- FALSE

  if (!requireNamespace("rlang", quietly = TRUE) || !requireNamespace("ggplot2", quietly = TRUE)) {
    needs_update <- TRUE
  } else {
    rlang_version <- packageVersion("rlang")
    ggplot2_version <- packageVersion("ggplot2")

    min_rlang <- "1.1.0"
    min_ggplot2 <- "3.5.0"

    if (rlang_version < min_rlang || ggplot2_version < min_ggplot2) {
      cat("  [WARN] Versiones desactualizadas detectadas\n")
      needs_update <- TRUE
    }
  }

  if (needs_update) {
    cat("  to Actualizando paquetes...\n")
    install.packages("rlang", repos = "https://cloud.r-project.org/", quiet = TRUE)
    install.packages("ggplot2", repos = "https://cloud.r-project.org/", dependencies = TRUE, quiet = TRUE)
    cat("  [OK] Paquetes actualizados\n")
  } else {
    cat("  [OK] Versiones compatibles detectadas\n")
  }
}

check_and_update_packages()
cat("\n")

# === CARGAR PAQUETES ===
if (!require("pacman")) {
  install.packages("pacman")
}

required_packages <- c(
  "tidyverse", "pheatmap", "ComplexHeatmap", "RColorBrewer", "viridis",
  "scales", "showtext", "sysfonts", "grid", "circlize", "devtools", "rlang",
  "flextable", "officer", "readxl"
)
pacman::p_load(char = required_packages, install = TRUE, update = FALSE)

# KEGGREST is optional; without it the KEGG descriptions fall back to bare IDs
KEGGREST_AVAILABLE <- FALSE
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org/")
}
if (!requireNamespace("KEGGREST", quietly = TRUE)) {
  tryCatch(
    {
      BiocManager::install("KEGGREST", ask = FALSE, update = FALSE)
      library(KEGGREST)
      KEGGREST_AVAILABLE <- TRUE
      cat("[OK] KEGGREST disponible para consultar descripciones de KEGG\n")
    },
    error = function(e) {
      cat("[WARN] KEGGREST no disponible. Solo se usarán IDs de KEGG.\n")
      KEGGREST_AVAILABLE <- FALSE
    }
  )
} else {
  library(KEGGREST)
  KEGGREST_AVAILABLE <- TRUE
  cat("[OK] KEGGREST disponible\n")
}

pacman::p_load_gh("zumbov2/deeplr")

# DeepL key is read from the environment; set DEEPL_AUTH_KEY in ~/.Renviron.
# Without it translate_to_spanish() falls back to the manual table and the
# remaining labels stay in English.
DEEPL_AUTH_KEY <- Sys.getenv("DEEPL_AUTH_KEY")
if (nzchar(DEEPL_AUTH_KEY)) {
  options(deeplr.auth_key = DEEPL_AUTH_KEY)
} else {
  warning("DEEPL_AUTH_KEY is not set; labels stay in English.")
}

showtext_auto()
showtext_opts(dpi = 300)

# ---- fonts ----
font_family <- "sans"
available_fonts <- sysfonts::font_files()

font_priority <- c(
  "Source Sans Pro", "SourceSansPro",
  "CMU Serif", "CMU Sans Serif", "Computer Modern", "Latin Modern",
  "Helvetica", "Helvetica Neue"
)

# Try to load local fonts directly (highest priority: Lato > Source Sans 3)
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
        cat("[OK] Using font:", lf$name, "(local ~/.local/share/fonts)\n")
        font_found <- TRUE
        break
      }, error = function(e) {
        cat("[WARN] Found", lf$name, "files but could not load them:", conditionMessage(e), "\n")
      })
    }
  }
}

# Fallback: search system font registry
for (font_name in font_priority) {
  matching_fonts <- available_fonts[grepl(font_name, available_fonts$family, ignore.case = TRUE), ]
  if (nrow(matching_fonts) > 0) {
    tryCatch(
      {
        font_path <- matching_fonts$file[1]
        sysfonts::font_add(
          family = "custom_font",
          regular = font_path
        )
        font_family <- "custom_font"
        cat("[OK] Using font:", font_name, "\n")
        font_found <- TRUE
        break
      },
      error = function(e) {}
    )
  }
}

if (!font_found) {
  cat("[WARN] No suitable system font found. Attempting to download from Google Fonts...\n")
  google_priority <- list(
    list(name = "Source Sans Pro", family = "Source Sans Pro"),
    list(name = "Lato",            family = "Lato"),
    list(name = "Open Sans",       family = "Open Sans")
  )
  for (gf in google_priority) {
    tryCatch({
      sysfonts::font_add_google(gf$name, "custom_font")
      font_family <- "custom_font"
      cat("[OK] Downloaded and using:", gf$name, "(Google Fonts)\n")
      font_found <- TRUE
      break
    }, error = function(e) {})
  }
  if (!font_found) {
    cat("[WARN] Could not download Google Fonts. Using default sans-serif font.\n")
  }
}

cat("[OK] Font configuration complete\n\n")

# Injects the font family into any grid::gpar(), so ComplexHeatmap renders its
# titles, legends, row and column names and annotations in the same face as the
# ggplots. ComplexHeatmap draws through the base PNG device plus showtext, which
# does not inherit the ggplot theme. An explicit fontfamily is left alone.
gpar_f <- function(...) {
  args <- list(...)
  if (is.null(args$fontfamily)) args$fontfamily <- font_family
  do.call(grid::gpar, args)
}

# Italicises taxon names, following tax_analysis.R.
#
# DisplayName carries a station and bin suffix, as in "Thalassarchaeum betae A4
# BTM_21", so testing the whole string with grepl("[0-9]", ...) is always TRUE
# and nothing would ever be italicised. Only the first word, the genus, decides
# the face for the whole string:
#   upright if it carries digits (GTDB codes: UBA57, sp002507125)
#   upright if it is all caps (MGIII, DPANN)
#   upright if the name is NA or "Sin asignar"
#   italic otherwise, meaning an established genus or species name
make_row_fontface <- function(names) {
  first_word <- sub("\\s.*$", "", names)   # primera palabra del DisplayName
  ifelse(
    is.na(names) | names == "Sin asignar" |
      grepl("[0-9]", first_word) | grepl("^[A-Z]+$", first_word),
    "plain",
    "italic"
  )
}

# CDC scientific nomenclature (https://wwwnc.cdc.gov/eid/page/scientific-nomenclature):
# only genus and species (and infraspecific ranks) are italicised; ranks above
# genus — family, order, class, phylum, domain — remain plain text.
ITALIC_LEVELS <- c("Genus", "Species", "Subspecies", "Strain")

should_italicize_level <- function(level) {
  !is.null(level) && length(level) == 1 && !is.na(level) && level %in% ITALIC_LEVELS
}

# Builds plotmath expressions that italicise only the taxonomic part of a name.
# Word by word: no digits and not all caps goes italic, everything else upright.
# Only applies at genus level and below.
#   level en {Family, Order, Class, Phylum, Domain} to todo plain (CDC).
# Ejemplos (level="Species"):
#   "Thalassarchaeum betae A6 BTM_56"         to paste(italic('Thalassarchaeum betae'), ' A6 BTM_56')
#   "Thalassarchaeum sp002507125 PO1 BTM_31"  to paste(italic('Thalassarchaeum'), ' sp002507125 PO1 BTM_31')
#   "MGIIb-Q1 A4 1000_28"                     to 'MGIIb-Q1 A4 1000_28'  (todo plain)
make_taxon_label_expr <- function(names, level = "Species") {
  italic_ok <- should_italicize_level(level)
  result <- lapply(names, function(nm) {
    if (is.na(nm) || nm == "Sin asignar") return(nm)
    if (!italic_ok) {
      return(parse(text = paste0("'", gsub("'", "\\'", nm, fixed = TRUE), "'"))[[1]])
    }
    words    <- strsplit(nm, " ")[[1]]
    is_plain <- grepl("[0-9]", words) | grepl("^[A-Z]+$", words)
    if (all(is_plain)) return(nm)
    r      <- rle(is_plain)
    ends   <- cumsum(r$lengths)
    starts <- c(1L, ends[-length(ends)] + 1L)
    chunks <- character(length(r$lengths))
    for (i in seq_along(r$lengths)) {
      grp <- paste(words[starts[i]:ends[i]], collapse = " ")
      sp  <- if (starts[i] == 1L) "" else " "
      if (r$values[i]) {
        chunks[i] <- paste0("'", sp, grp, "'")
      } else {
        chunks[i] <- paste0("italic('", sp, grp, "')")
      }
    }
    if (!any(r$values)) {
      return(parse(text = paste0("italic('", nm, "')"))[[1]])
    }
    parse(text = paste0("paste(", paste(chunks, collapse = ", "), ")"))[[1]]
  })
  do.call(expression, result)
}

# Draws a list of Legend objects into the bottom right of the canvas
.draw_br_legend <- function(lgd_list, x_mm = 18, y_mm = 50) {
  lgds <- Filter(Negate(is.null), lgd_list)
  if (length(lgds) == 0) return(invisible(NULL))
  pd <- do.call(packLegend, c(lgds, list(direction = "horizontal", gap = unit(6, "mm"))))
  draw(pd,
    x    = unit(1, "npc") - unit(x_mm, "mm"),
    y    = unit(y_mm, "mm"),
    just = c("right", "bottom")
  )
}

# Builds the depth legend, plus a station legend where stations are present
.make_annotation_legends <- function(dc, stations = NULL) {
  lgds <- list()
  lgds[["depth"]] <- Legend(
    labels     = names(dc),
    legend_gp  = gpar_f(fill = unname(dc)),
    title      = "Profundidad",
    title_gp   = gpar_f(fontsize = 10, fontface = "bold"),
    labels_gp  = gpar_f(fontsize = 9)
  )
  if (!is.null(stations)) {
    unique_st  <- sort(unique(na.omit(as.character(stations))))
    n_st       <- max(length(unique_st), 3)
    st_cols    <- setNames(RColorBrewer::brewer.pal(n_st, "Set1")[seq_len(length(unique_st))], unique_st)
    lgds[["station"]] <- Legend(
      labels    = names(st_cols),
      legend_gp = gpar_f(fill = unname(st_cols)),
      title     = "Estación",
      title_gp  = gpar_f(fontsize = 10, fontface = "bold"),
      labels_gp = gpar_f(fontsize = 9)
    )
  }
  lgds
}

# ---- TEMA PERSONALIZADO ----
theme_publication <- function(base_size = 11) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray", linetype = "dotted"),
      axis.text = element_text(color = "black"),
      axis.title = element_text(face = "bold"),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = NA),
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 11),
      strip.background = element_rect(fill = "gray90", color = "black"),
      strip.text = element_text(face = "bold", size = 13),
      plot.title = element_text(face = "bold", hjust = 0.5),
      text = element_text(family = font_family)
    )
}

# ---- 1. CARGA Y PREPARACIÓN DE DATOS ----
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

OUTPUT_DIR <- file.path(PROJECT_DIR, "results", "SO_01", "funct_annot")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Clear the previous run, keeping the translation caches
cat("Limpiando archivos de corrida anterior (preservando caché)...\n")
{
  all_items <- list.files(OUTPUT_DIR, full.names = TRUE, recursive = FALSE)
  all_files <- all_items[!file.info(all_items)$isdir]
  cache_files <- all_files[grepl("cache|translation_cache|\\.rds$", basename(all_files), ignore.case = TRUE)]
  files_to_delete <- setdiff(all_files, cache_files)
  if (length(files_to_delete) > 0) {
    file.remove(files_to_delete)
    cat("  [OK] Eliminados", length(files_to_delete), "archivos de corrida anterior (caché preservado)\n")
  } else {
    cat("  [OK] No hay archivos anteriores que eliminar\n")
  }
}
cat("Directorio de salida:", OUTPUT_DIR, "\n")

METADATA_FILE <- file.path(PROJECT_DIR, "MAGs_metadata.tsv")
cat("Leyendo archivo de metadatos...\n")
metadata <- read.delim(METADATA_FILE, sep = "\t", header = TRUE, stringsAsFactors = FALSE)

cat("Leyendo archivo de anotaciones DRAM2...\n")
dram <- read_tsv(file.path(PROJECT_DIR, "results", "08_MAG_distill", "ANNOTATE", "raw-annotations.tsv"),
  show_col_types = FALSE,
  na = c("", "NA")
)

cat("Total de genes anotados:", nrow(dram), "\n")
cat("Total de MAGs:", n_distinct(dram$input_fasta), "\n")

# Verificar columnas disponibles
cat("\nColumnas de descripción disponibles:\n")
available_desc_cols <- c()
if ("uniref_description" %in% colnames(dram)) {
  cat("  [OK] uniref_description\n")
  available_desc_cols <- c(available_desc_cols, "uniref_description")
}
if ("dbcan_description" %in% colnames(dram)) {
  cat("  [OK] dbcan_description\n")
  available_desc_cols <- c(available_desc_cols, "dbcan_description")
}
if ("camper_description" %in% colnames(dram)) {
  cat("  [OK] camper_description\n")
  available_desc_cols <- c(available_desc_cols, "camper_description")
}
if ("kofam_description" %in% colnames(dram)) {
  cat("  [OK] kofam_description\n")
  available_desc_cols <- c(available_desc_cols, "kofam_description")
} else {
  cat("  [WARN] kofam_description NO disponible (se usará API o solo IDs)\n")
}
if ("kegg_id" %in% colnames(dram)) {
  cat("  [OK] kegg_id (KEGG DB – nueva base de datos)\n")
  available_desc_cols <- c(available_desc_cols, "kegg_id")
} else {
  cat("  [WARN] kegg_id NO disponible\n")
}

# ---- FUNCIONES DE MAPEO AUTOMÁTICO: ID to DESCRIPCIÓN FUNCIONAL ----

# --- 1. MAPEO DE KEGG KOs ---
get_kegg_description <- function(ko_ids, use_cache = TRUE, cache_dir = NULL) {
  if (!KEGGREST_AVAILABLE) {
    return(ko_ids) # Retornar IDs si KEGGREST no está disponible
  }

  if (is.null(cache_dir)) {
    cache_file <- ".kegg_cache.rds"
  } else {
    cache_file <- file.path(cache_dir, ".kegg_cache.rds")
  }

  if (use_cache && file.exists(cache_file)) {
    kegg_cache <- readRDS(cache_file)
  } else {
    kegg_cache <- list()
  }

  ko_ids_unique <- unique(ko_ids[!is.na(ko_ids) & ko_ids != ""])
  descriptions <- character(length(ko_ids))
  names(descriptions) <- ko_ids

  for (i in seq_along(ko_ids_unique)) {
    ko <- ko_ids_unique[i]

    if (ko %in% names(kegg_cache)) {
      descriptions[ko] <- kegg_cache[[ko]]
    } else {
      tryCatch(
        {
          result <- keggGet(ko)
          if (length(result) > 0 && !is.null(result[[1]]$NAME)) {
            desc <- result[[1]]$NAME[1]
            descriptions[ko] <- desc
            kegg_cache[[ko]] <- desc
          } else {
            descriptions[ko] <- ko
          }
          Sys.sleep(0.15)
        },
        error = function(e) {
          descriptions[ko] <- ko
        }
      )
    }

    if (i %% 10 == 0) {
      cat("    Procesados", i, "de", length(ko_ids_unique), "KOs\r")
    }
  }

  if (use_cache) {
    saveRDS(kegg_cache, cache_file)
  }

  return(descriptions[ko_ids])
}

# --- 2. MAPEO DE CAZy ---
get_cazy_description <- function(cazy_ids) {
  descriptions <- character(length(cazy_ids))

  for (i in seq_along(cazy_ids)) {
    id <- cazy_ids[i]

    if (is.na(id) || id == "") {
      descriptions[i] <- NA
      next
    }

    family_type <- gsub("([A-Z]+)[0-9_].*", "\\1", id)
    family_num <- gsub("[A-Z]+([0-9]+).*", "\\1", id)

    family_names <- c(
      "GT" = "Glycosyltransferase",
      "GH" = "Glycoside Hydrolase",
      "CE" = "Carbohydrate Esterase",
      "PL" = "Polysaccharide Lyase",
      "AA" = "Auxiliary Activity",
      "CBM" = "Carbohydrate-Binding Module"
    )

    if (family_type %in% names(family_names)) {
      descriptions[i] <- paste0(family_names[family_type], " family ", family_num)
    } else {
      descriptions[i] <- id
    }
  }

  return(descriptions)
}

# --- 3. FUNCIÓN GENERAL: EXTRAER DESCRIPCIÓN DE DRAM2 ---
get_description_from_dram <- function(ids, id_col, desc_col, data) {
  # Confirm the columns are there
  if (!id_col %in% colnames(data) || !desc_col %in% colnames(data)) {
    return(ids) # Retornar IDs si las columnas no existen
  }

  id_to_desc <- data %>%
    filter(!is.na(!!sym(id_col)) & !!sym(id_col) != "") %>%
    filter(!is.na(!!sym(desc_col)) & !!sym(desc_col) != "") %>%
    select(!!sym(id_col), !!sym(desc_col)) %>%
    distinct() %>%
    separate_rows(!!sym(id_col), !!sym(desc_col), sep = ";\\s*") %>%
    group_by(!!sym(id_col)) %>%
    summarise(description = first(!!sym(desc_col)), .groups = "drop") %>%
    deframe()

  descriptions <- id_to_desc[ids]
  descriptions[is.na(descriptions)] <- ids[is.na(descriptions)]

  return(descriptions)
}

# --- 4. FUNCIÓN PARA LIMPIAR Y TRUNCAR DESCRIPCIONES ---
clean_description <- function(desc, max_length = 60) {
  if (is.na(desc) || desc == "") {
    return(desc)
  }

  desc <- gsub("\\s+\\[.*?\\]", "", desc)
  desc <- gsub("\\s*;.*$", "", desc)
  desc <- trimws(desc)

  if (nchar(desc) > max_length) {
    desc <- paste0(substr(desc, 1, max_length - 3), "...")
  }

  return(desc)
}

# ---- 5. translate labels into Spanish through DeepL ----
translate_to_spanish <- function(names_vec, translation_cache_file = NULL) {
  cat("  Traduciendo", length(names_vec), "nombres al español...\n")

  # The cache always loads as a list, so [[ ]] never errors on a subscript
  if (!is.null(translation_cache_file) && file.exists(translation_cache_file)) {
    cache_df <- read.csv(translation_cache_file, stringsAsFactors = FALSE)
    translation_cache <- as.list(setNames(cache_df$Spanish, cache_df$English))
  } else {
    translation_cache <- list()
  }

  translated <- character(length(names_vec))
  new_translations <- 0

  # Guarantees a length-1 character, whatever comes back
  safe_string <- function(x, fallback) {
    if (length(x) == 0 || is.na(x[1]) || x[1] == "") fallback else as.character(x[1])
  }

  for (i in seq_along(names_vec)) {
    original <- names_vec[i]
    if (original %in% names(.MANUAL_TRANSLATIONS)) {
      # Manual translations win outright, bypassing both DeepL and the cache
      translated[i] <- .MANUAL_TRANSLATIONS[[original]]
      translation_cache[[original]] <- .MANUAL_TRANSLATIONS[[original]]
    } else if (original %in% names(translation_cache)) {
      # Guard against empty or NULL entries left in the cache by an earlier run
      translated[i] <- safe_string(translation_cache[[original]], original)
    } else {
      result <- tryCatch(
        deeplr::translate2(text = original, source_lang = "EN", target_lang = "ES"),
        error = function(e) original
      )
      # Guard against an empty DeepL result, which returns character(0) rather than erroring
      result <- safe_string(result, original)
      translated[i] <- result
      translation_cache[[original]] <- result
      new_translations <- new_translations + 1
      Sys.sleep(0.1)
    }
    if (i %% 10 == 0) cat("    Progreso:", i, "/", length(names_vec), "\n")
  }

  cat(
    "  [OK] Traducciones completadas:", new_translations, "nuevas,",
    length(names_vec) - new_translations, "desde caché\n"
  )

  if (!is.null(translation_cache_file)) {
    keys <- names(translation_cache)
    values <- vapply(translation_cache, function(x) safe_string(x, ""), character(1))
    keep <- !is.na(keys) & keys != ""
    cache_out <- data.frame(
      English = keys[keep], Spanish = values[keep],
      stringsAsFactors = FALSE
    )
    write.csv(cache_out, translation_cache_file, row.names = FALSE)
  }

  return(translated)
}

# Manual translations, which take priority over DeepL and the cache
.MANUAL_TRANSLATIONS <- list(
  "Aerobe"           = "Aerobio",
  "Anaerobe"         = "Anaerobio",
  "aerobe"           = "aerobio",
  "anaerobe"         = "anaerobio",
  "C1 metabolism"    = "Metabolismo C1",
  "Carbon fixation"  = "Fijación de carbono",
  "Complex I"        = "Complejo I",
  "Complex II"       = "Complejo II",
  "Complex III"      = "Complejo III",
  "Complex IV"       = "Complejo IV",
  "Complex IV H"     = "Complejo IV-H",
  "Complex IV L"     = "Complejo IV-L",
  "Complex IV-H"     = "Complejo IV-H",
  "Complex IV-L"     = "Complejo IV-L",
  "Complex V"        = "Complejo V",
  "Nitrogen"         = "Nitrógeno",
  "Sulfur"           = "Azufre",
  "Photo"            = "Fotosíntesis",
  "Other"            = "Otro",
  "SCFA/alcohol"              = "AGCC/alcohol",
  "SCFA/alcohol conversions" = "AGCC/conversión alcoholes",
  "CAZy"             = "CAZy",
  "Pyruvate"         = "Piruvato",
  "Piruvate"         = "Piruvato",
  "Acetyl"           = "acetil",
  # Rutas pyruvate/acetyl (DeepL fallaba en estas frases cortas)
  "Acetyl-CoA pathway, CO2 => acetyl-CoA"        = "Vía del acetil-CoA, CO₂ => acetil-CoA",
  "pyruvate => acetyl CoA v1"                    = "piruvato => acetil-CoA v1",
  "pyruvate => acetyl CoA v2"                    = "piruvato => acetil-CoA v2",
  "pyruvate => acetylCoA f+ formate v3"          = "piruvato => acetil-CoA + formato v3",

  # ---- energy traits, where machine translation is not reliable ----
  "Aerobic Ammonia Oxidizer"  = "Oxidador aerobio de amonio (AOA)",
  "Comammox"                  = "Comammox",
  "DNRA"                      = "DNRA",
  "Fermenter"                 = "Fermentador",
  "Fumarate Reduction"        = "Reducción de fumarato",
  "Iron Oxidation"            = "Oxidación de hierro",
  "Iron Reduction"            = "Reducción de hierro",
  "Methanogen"                = "Metanógeno",
  "Methanotroph"              = "Metanótrofo",
  "Microaerophillic"          = "Microaerofílico",
  "N reducer"                 = "Reductor de N",
  "Nitrogen fixer"            = "Fijador de N\u2082",
  "Phototroph"                = "Fotótrofo",
  "Probable Iron Reduction"   = "Posible reducción de hierro",
  "S oxidizer"                = "Oxidador de S",
  "S reducer"                 = "Reductor de S",
  "Tetrathionate Reduction"   = "Reducción de tetrationato",
  "annamox"                   = "Anammox",

  # ---- the ARCHAEA-MOD traits, from the archaea-calibrated rule set ----
  # Sulfur cycle, expanded
  "Dissimilatory Sulfate Reducer Strict" = "Reductor disimilatorio de sulfato (estricto)",
  "Assimilatory Sulfate Reducer"         = "Reductor asimilatorio de sulfato",
  "Sulfur Oxidizer SOX Complete"         = "Oxidador de azufre SOX (completo)",
  "Sulfide Oxidizer Sqr Fcc"             = "Oxidador de sulfuro (Sqr/Fcc)",
  "Thiosulfate Reducer PhsABC"           = "Reductor de tiosulfato (PhsABC)",
  "DMSO/TMAO Reducer DmsABC"             = "Reductor de DMSO/TMAO (DmsABC)",

  # Hydrogenases, classified per HydDB
  "Hydrogenase Group 1 NiFe"                  = "Hidrogenasa Grupo 1 [NiFe]",
  "Hydrogenase Group 3a NiFe F420"            = "Hidrogenasa Grupo 3a [NiFe] F420",
  "Hydrogenase Group 3c NiFe Mvh"             = "Hidrogenasa Grupo 3c [NiFe] Mvh",
  "Hydrogenase Group 3d NiFe Hox"             = "Hidrogenasa Grupo 3d [NiFe] Hox",
  "Hydrogenase Group 4 NiFe Energy Converting" = "Hidrogenasa Grupo 4 [NiFe] (EC)",
  "Hydrogenase FeFe HydA"                     = "Hidrogenasa [FeFe] HydA",
  "Hydrogenase FeFe Hnd"                      = "Hidrogenasa [FeFe] Hnd",
  "Any Hydrogenase"                           = "Cualquier hidrogenasa",
  "Hydrogenotroph"                            = "Hidrogenotrófico",
  "Hydrogen Evolver"                          = "Productor de H\u2082",

  # CO oxidation
  "Aerobic CO Oxidizer"   = "Oxidador aerobio de CO",
  "Anaerobic CO Oxidizer" = "Oxidador anaerobio de CO",

  # ANME y ureasa
  "Anaerobic Methane Oxidizer ANME" = "Oxidador anaerobio de metano (ANME)",
  "Urease Nitrogen Mineralization"  = "Ureasa (mineralización de N)"
)

# --- 6. FUNCIÓN WRAPPER PARA CREAR LABELS INFORMATIVOS ---
create_functional_labels <- function(ids, id_col, data, max_length = 50, use_kegg_api = FALSE) {
  if (id_col == "kofam_id") {
    # The KEGG columns carry no description field
    if (use_kegg_api && KEGGREST_AVAILABLE) {
      descriptions <- get_kegg_description(ids, use_cache = TRUE, cache_dir = OUTPUT_DIR)
    } else {
      # Without the API, the bare IDs are all there is
      descriptions <- ids
    }
  } else if (id_col == "dbcan_id") {
    # dbcan does carry a description
    descriptions <- get_description_from_dram(ids, "dbcan_id", "dbcan_description", data)
    missing <- is.na(descriptions) | descriptions == ids
    if (any(missing)) {
      descriptions[missing] <- get_cazy_description(ids[missing])
    }
  } else if (id_col == "uniref_id") {
    # uniref does carry a description
    descriptions <- get_description_from_dram(ids, "uniref_id", "uniref_description", data)
  } else if (id_col == "camper_id") {
    # camper does carry a description
    descriptions <- get_description_from_dram(ids, "camper_id", "camper_description", data)
  } else if (id_col == "metals_id") {
    # metals carries no description
    descriptions <- ids
  } else if (id_col == "sulfur_id") {
    # sulfur carries no description
    descriptions <- ids
  } else if (id_col == "fegenie_id") {
    # fegenie carries no description
    descriptions <- ids
  } else if (id_col == "methyl_id") {
    # methyl carries no description
    descriptions <- ids
  } else if (id_col == "kegg_id") {
    # the KEGG database does, in kegg_description
    descriptions <- get_description_from_dram(ids, "kegg_id", "kegg_description", data)
  } else {
    descriptions <- ids
  }

  descriptions <- sapply(descriptions, clean_description, max_length = max_length)
  labels <- paste0(ids, ": ", descriptions)

  return(labels)
}

# ---- HELPER FUNCTIONS ----

make_unique_names <- function(names) {
  duplicated_names <- duplicated(names) | duplicated(names, fromLast = TRUE)
  if (any(duplicated_names)) {
    names_table <- table(names)
    for (name in names(names_table[names_table > 1])) {
      indices <- which(names == name)
      names[indices] <- paste0(name, "_", seq_along(indices))
    }
  }
  return(names)
}

order_matrix_by_depth <- function(mat, annotation_df) {
  depth_order <- if (exists("depth_levels_std")) depth_levels_std else c("ZMOR", "1000 m", "Fondo")
  annotation_df$Profundidad <- factor(annotation_df$Profundidad, levels = depth_order)
  row_order <- order(annotation_df$Profundidad)
  mat_ordered <- mat[row_order, , drop = FALSE]
  annotation_ordered <- annotation_df[row_order, , drop = FALSE]
  return(list(matrix = mat_ordered, annotation = annotation_ordered))
}

# Re-sorts rows for annotation versions 4 (station) and 5 (completeness).
# Must run after order_matrix_by_depth, since it refines that ordering rather
# than replacing it.
#   dn_to_id          display_name to genome_id, for looking up the station
#   completeness_map  display_name to percent completeness, for version 5
reorder_rows_by_version <- function(mat, ann_df, annotation_version,
                                    dn_to_id = NULL, completeness_map = NULL) {
  v <- as.integer(annotation_version)
  if (v == 4L && !is.null(dn_to_id)) {
    sts <- as.character(get_station(dn_to_id[rownames(mat)]))
    sts[is.na(sts)] <- "zzz"
    ord <- order(sts)
  } else if (v == 5L && !is.null(completeness_map)) {
    compl <- completeness_map[rownames(mat)]
    compl[is.na(compl)] <- 0
    ord <- order(compl) # menor arriba to mayor abajo
  } else {
    return(list(mat = mat, ann = ann_df))
  }
  list(
    mat = mat[ord, , drop = FALSE],
    ann = if (!is.null(ann_df)) ann_df[ord, , drop = FALSE] else NULL
  )
}

# ---- 2. CLASIFICACIÓN DE CALIDAD DE ANOTACIÓN ----
# What counts as a real functional annotation: a gene qualifies only on a hit
# against a curated or HMM functional database
# (KEGG, KOfam, dbCAN, CAMPER, sulfur, methyl, FeGenie, metals). Un hit en
# UniRef does not count. It is sequence homology, and the reference protein it
# points at may itself be uncharacterised, so a UniRef hit can mean nothing more
# than "resembles something nobody has identified". Those are reported separately
# as Homology_Only, split into informative and vague by what the hit description
# actually says.
cat("\nClasificando genes por calidad de anotación...\n")

count_hits <- function(column) {
  sum(!is.na(column) & column != "")
}

is_informative_uniref <- function(uniref_id, uniref_desc) {
  if (is.na(uniref_id) | uniref_id == "") {
    return(FALSE)
  }
  if (is.na(uniref_desc) | uniref_desc == "") {
    return(FALSE)
  }

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
  is_uninformative <- any(sapply(non_informative_patterns, function(p) grepl(p, desc_lower)))

  return(!is_uninformative)
}

dram_classified <- dram %>%
  mutate(
    has_kegg = !is.na(kofam_id) & kofam_id != "",
    has_cazy = !is.na(dbcan_id) & dbcan_id != "",
    has_metals = !is.na(metals_id) & metals_id != "",
    has_camper = !is.na(camper_id) & camper_id != "",
    has_methyl = !is.na(methyl_id) & methyl_id != "",
    has_sulfur = !is.na(sulfur_id) & sulfur_id != "",
    has_fegenie = !is.na(fegenie_id) & fegenie_id != "",
    has_uniref = !is.na(uniref_id) & uniref_id != "",
    has_informative_uniref = mapply(is_informative_uniref, uniref_id, uniref_description),
    has_vague_uniref = has_uniref & !has_informative_uniref,
    has_kegg_db = if ("kegg_id" %in% names(pick(everything()))) !is.na(kegg_id) & kegg_id != "" else FALSE,
    # Curated functional annotation: a hit in any HMM or curated database, UniRef excluded.
    has_curated_functional = has_kegg | has_cazy | has_metals | has_camper |
      has_methyl | has_sulfur | has_fegenie | has_kegg_db,
    Annotation_Quality = case_when(
      has_curated_functional        ~ "Curated_Functional",
      has_uniref                    ~ "Homology_Only",
      TRUE                          ~ "No_Annotation"
    )
  )

cat("[OK] Clasificación completada\n")

# Tabla resumen general
summary_table <- dram_classified %>%
  summarise(
    Total_Genes = n(),
    dbCAZy = count_hits(dbcan_id),
    Metals = count_hits(metals_id),
    CAMPER = count_hits(camper_id),
    Methyl = count_hits(methyl_id),
    KEGG_KOfam = count_hits(kofam_id),
    KEGG_DB = if ("kegg_id" %in% names(pick(everything()))) count_hits(kegg_id) else 0L,
    Sulfur = count_hits(sulfur_id),
    FeGenie = count_hits(fegenie_id),
    UniRef_Total = count_hits(uniref_id),
    UniRef_Informative = sum(has_informative_uniref),
    UniRef_Vague = sum(has_vague_uniref),
    # The headline metric: curated functional annotation, UniRef excluded
    Curated_Functional_Annotation = sum(Annotation_Quality == "Curated_Functional"),
    # UniRef-only hits, split by how informative the description is
    Homology_Only          = sum(Annotation_Quality == "Homology_Only"),
    Homology_Only_Informative = sum(Annotation_Quality == "Homology_Only" & has_informative_uniref),
    Homology_Only_Vague       = sum(Annotation_Quality == "Homology_Only" & has_vague_uniref),
    No_Annotation = sum(Annotation_Quality == "No_Annotation")
  ) %>%
  mutate(
    Percent_Curated_Functional = (Curated_Functional_Annotation / Total_Genes) * 100,
    Percent_Homology_Only      = (Homology_Only                  / Total_Genes) * 100,
    Percent_Any_Hit            = ((Curated_Functional_Annotation + Homology_Only) / Total_Genes) * 100,
    Percent_No_Annotation      = (No_Annotation                  / Total_Genes) * 100
  )

write_csv(summary_table, file.path(OUTPUT_DIR, "01_summary_total_hits.csv"))
cat("[OK] Tabla resumen general guardada\n")

# Per-MAG summary table
summary_by_mag <- dram_classified %>%
  group_by(MAG = input_fasta) %>%
  summarise(
    Total_Genes = n(),
    dbCAZy = count_hits(dbcan_id),
    Metals = count_hits(metals_id),
    CAMPER = count_hits(camper_id),
    Methyl = count_hits(methyl_id),
    KEGG_KOfam = count_hits(kofam_id),
    KEGG_DB = if ("kegg_id" %in% names(pick(everything()))) count_hits(kegg_id) else 0L,
    Sulfur = count_hits(sulfur_id),
    FeGenie = count_hits(fegenie_id),
    UniRef_Total = count_hits(uniref_id),
    UniRef_Informative = sum(has_informative_uniref),
    UniRef_Vague = sum(has_vague_uniref),
    Curated_Functional_Annotation = sum(Annotation_Quality == "Curated_Functional"),
    Homology_Only          = sum(Annotation_Quality == "Homology_Only"),
    Homology_Only_Informative = sum(Annotation_Quality == "Homology_Only" & has_informative_uniref),
    Homology_Only_Vague       = sum(Annotation_Quality == "Homology_Only" & has_vague_uniref),
    No_Annotation = sum(Annotation_Quality == "No_Annotation"),
    Multiple_KOs = sum(grepl(";", kofam_id), na.rm = TRUE),
    # Genes supported by at least two curated databases, UniRef excluded
    Multi_DB_Hits = sum(
      (has_kegg + has_cazy + has_metals + has_camper +
        has_methyl + has_sulfur + has_fegenie + has_kegg_db) >= 2
    ),
    .groups = "drop"
  ) %>%
  left_join(
    metadata %>% select(NameWithoutExtension, Genus, Species, FullTaxonomy),
    by = c("MAG" = "NameWithoutExtension")
  ) %>%
  mutate(
    Percent_Curated_Functional = (Curated_Functional_Annotation / Total_Genes) * 100,
    Percent_Homology_Only      = (Homology_Only                  / Total_Genes) * 100,
    Percent_Any_Hit            = ((Curated_Functional_Annotation + Homology_Only) / Total_Genes) * 100
  ) %>%
  arrange(desc(Percent_Curated_Functional))

write_csv(summary_by_mag, file.path(OUTPUT_DIR, "02_summary_by_MAG.csv"))
cat("[OK] Tabla resumen por MAG guardada\n")

dram_classified %>%
  select(
    query_id, input_fasta, kofam_id, dbcan_id, uniref_id,
    uniref_description, Annotation_Quality
  ) %>%
  write_csv(file.path(OUTPUT_DIR, "02b_gene_annotation_quality.csv"))
cat("[OK] Tabla de clasificación de genes guardada\n")

# Mapeo de metadata
mag_to_display <- metadata %>%
  select(NameWithoutExtension, DisplayName) %>%
  deframe()

mag_to_depth <- metadata %>%
  select(NameWithoutExtension, Depth) %>%
  deframe()
.depth_rename <- c("ROMZ" = "ZMOR", "1000" = "1000 m", "BTM" = "Fondo")
mag_to_depth <- setNames(
  ifelse(mag_to_depth %in% names(.depth_rename), .depth_rename[mag_to_depth], mag_to_depth),
  names(mag_to_depth)
)
metadata$Depth <- ifelse(metadata$Depth %in% names(.depth_rename),
  .depth_rename[metadata$Depth], metadata$Depth
)

mag_to_station <- metadata %>%
  select(NameWithoutExtension, Station) %>%
  deframe()

convert_to_display_name <- function(mag_names) {
  display_names <- mag_to_display[mag_names]
  display_names[is.na(display_names)] <- mag_names[is.na(display_names)]
  display_names <- make_unique_names(display_names)
  return(display_names)
}

get_depth <- function(mag_names) {
  depths <- mag_to_depth[mag_names]
  return(depths)
}

get_station <- function(mag_names) {
  stations <- mag_to_station[mag_names]
  return(stations)
}

# ---- capitalize_first ----
capitalize_first <- function(x) {
  ifelse(is.na(x), x, paste0(toupper(substring(x, 1, 1)), substring(x, 2)))
}

# ---- get_genus: genome_id to genus ----
mag_to_genus <- setNames(
  as.character(metadata$Genus),
  metadata$NameWithoutExtension
)
mag_to_genus[is.na(mag_to_genus) | trimws(mag_to_genus) %in% c("", "N/A")] <- "Desconocido"

get_genus <- function(mag_names) {
  g <- mag_to_genus[mag_names]
  g[is.na(g)] <- "Desconocido"
  as.character(g)
}

# ---- mag_to_order: genome_id to order, from GTDB FullTaxonomy ----
mag_to_order <- setNames(
  sapply(metadata$FullTaxonomy, function(tx) {
    m <- regmatches(tx, regexpr("o__[^;]+", tx))
    if (length(m) == 0 || m == "o__") NA_character_ else sub("^o__", "", m)
  }),
  metadata$NameWithoutExtension
)

# Helper: agrega cualquier data.frame (genome, col1, col2, ...) a nivel de Orden.
# Boolean and 0-to-1 numeric columns both collapse to a per-order mean, which reads as a proportion.
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

# ---- order_to_class: order to class, from GTDB FullTaxonomy ----
# Gives a code-like order name some context.
order_to_class <- local({
  cls_vec <- sapply(metadata$FullTaxonomy, function(tx) {
    m <- regmatches(tx, regexpr("c__[^;]+", tx))
    if (length(m) == 0 || m == "c__") NA_character_ else sub("^c__", "", m)
  })
  ord_vec <- sapply(metadata$FullTaxonomy, function(tx) {
    m <- regmatches(tx, regexpr("o__[^;]+", tx))
    if (length(m) == 0 || m == "o__") NA_character_ else sub("^o__", "", m)
  })
  tbl <- data.frame(order = ord_vec, class = cls_vec, stringsAsFactors = FALSE)
  tbl <- tbl[!is.na(tbl$order) & !is.na(tbl$class) &
               nchar(trimws(tbl$order)) > 0 & nchar(trimws(tbl$class)) > 0, ]
  tbl <- tbl[!duplicated(tbl$order), ]
  setNames(tbl$class, tbl$order)
})

# Puts the class in parentheses after an order name that is really a code: first
# word all caps (MGIII, DPANN) or carrying digits (UBA1234). Only applies to
# labels on figures collapsed by order.
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

# ---- 3. TABLAS FORMATEADAS CON FLEXTABLE ----
cat("\n=== GENERANDO TABLAS FORMATEADAS ===\n")

.has_flextable <- requireNamespace("flextable", quietly = TRUE) &&
  requireNamespace("officer", quietly = TRUE)
if (!.has_flextable) {
  cat("[WARN] flextable/officer no disponibles; se omiten tablas formateadas (.png/.pptx)\n")
}

# Table 1: annotation software compared
cat("\nGenerando Tabla 1: Comparación de software...\n")

# The DRAM2 column counts curated functional annotation only: KEGG and KOfam,
# dbCAZy, CAMPER, sulfur, methyl, FeGenie and metals. UniRef-only hits are left
# out, since sequence homology to an uncharacterised protein is not a function.
datos_anot <- data.frame(
  software = c("KofamKOALA", "eggNOG-mapper v2", "DRAM2"),
  genes = c(68539, 68539, summary_table$Total_Genes),
  hits = c(21255, 18179, summary_table$Curated_Functional_Annotation)
)

datos_anot <- datos_anot %>%
  mutate(
    sin_anotar = genes - hits,
    pct_anotacion = (hits / genes) * 100
  )

if (.has_flextable) {
  ft1 <- flextable(datos_anot) %>%
    set_header_labels(
      software = "Software",
      genes = "Total de genes",
      hits = "Genes anotados",
      sin_anotar = "Genes sin anotar",
      pct_anotacion = "% de anotación"
    ) %>%
    bold(part = "header") %>%
    align(align = "center", part = "header") %>%
    bg(bg = "#4472C4", part = "header") %>%
    color(color = "white", part = "header") %>%
    fontsize(size = 11, part = "header") %>%
    align(j = "software", align = "left", part = "body") %>%
    align(
      j = c("genes", "hits", "sin_anotar", "pct_anotacion"),
      align = "center", part = "body"
    ) %>%
    fontsize(size = 11, part = "body") %>%
    colformat_double(j = "pct_anotacion", digits = 2, suffix = "%") %>%
    colformat_double(
      j = c("genes", "hits", "sin_anotar"),
      digits = 0, big.mark = ","
    ) %>%
    theme_booktabs() %>%
    border_outer(border = fp_border(color = "black", width = 2), part = "all") %>%
    border_inner_h(border = fp_border(color = "gray70", width = 0.5), part = "body") %>%
    border_inner_v(border = fp_border(color = "gray70", width = 0.5), part = "body") %>%
    width(j = "software", width = 2.5) %>%
    width(j = "genes", width = 2.0) %>%
    width(j = "hits", width = 2.0) %>%
    width(j = "sin_anotar", width = 2.0) %>%
    width(j = "pct_anotacion", width = 2.0) %>%
    padding(padding.left = 10, padding.right = 10, part = "all") %>%
    padding(padding.top = 6, padding.bottom = 6, part = "all") %>%
    hrule(rule = "exact", part = "body") %>%
    height(height = 0.35, part = "body") %>%
    height(height = 0.4, part = "header")
  save_as_image(ft1,
    path = file.path(OUTPUT_DIR, "tabla_01_comparacion_software.png"),
    zoom = 5, expand = 20, webshot = "webshot2"
  )
  save_as_pptx(ft1, path = file.path(OUTPUT_DIR, "tabla_01_comparacion_software.pptx"))
  cat("[OK] Tabla 1 guardada\n")
}

# Table 2: DRAM2 broken down by database
cat("\nGenerando Tabla 2: DRAM2 por base de datos...\n")

datos_dram2 <- data.frame(
  Base_de_datos = c(
    "dbCAZy", "Metals", "CAMPER", "Methyl",
    "KEGG KOfam", "KEGG DB", "Sulfur", "FeGenie", "UniRef (informativo)"
  ),
  Hits = c(
    summary_table$dbCAZy, summary_table$Metals, summary_table$CAMPER,
    summary_table$Methyl, summary_table$KEGG_KOfam, summary_table$KEGG_DB,
    summary_table$Sulfur, summary_table$FeGenie, summary_table$UniRef_Informative
  ),
  Total_genes = summary_table$Total_Genes,
  stringsAsFactors = FALSE
)

datos_dram2 <- datos_dram2 %>%
  mutate(Pct_anotacion = (Hits / Total_genes) * 100)

datos_dram2_transpuesto <- datos_dram2 %>%
  select(Base_de_datos, Hits, Pct_anotacion) %>%
  pivot_longer(
    cols = c(Hits, Pct_anotacion),
    names_to = "Métrica",
    values_to = "Valor"
  ) %>%
  pivot_wider(
    names_from = Base_de_datos,
    values_from = Valor
  )

datos_dram2_transpuesto$Métrica <- c("Hits absolutos", "% de anotación")

if (.has_flextable) {
  ft2 <- flextable(datos_dram2_transpuesto) %>%
    set_header_labels(
      Métrica = "Métrica",
      dbCAZy = "dbCAZy",
      Metals = "Metals",
      CAMPER = "CAMPER",
      Methyl = "Methyl",
      `KEGG KOfam` = "KEGG KOfam",
      Sulfur = "Sulfur",
      FeGenie = "FeGenie",
      `UniRef (informativo)` = "UniRef"
    ) %>%
    bold(part = "header") %>%
    align(align = "center", part = "header") %>%
    bg(bg = "#2E75B6", part = "header") %>%
    color(color = "white", part = "header") %>%
    fontsize(size = 10, part = "header") %>%
    align(j = "Métrica", align = "left", part = "body") %>%
    align(j = 2:9, align = "center", part = "body") %>%
    bold(j = "Métrica", part = "body") %>%
    fontsize(size = 11, part = "body") %>%
    theme_booktabs() %>%
    border_outer(border = fp_border(color = "black", width = 2), part = "all") %>%
    border_inner_h(border = fp_border(color = "gray70", width = 1.5), part = "body") %>%
    border_inner_v(border = fp_border(color = "gray70", width = 0.5), part = "body") %>%
    colformat_double(i = 1, j = 2:9, digits = 0, big.mark = ",") %>%
    colformat_double(i = 2, j = 2:9, digits = 2, suffix = "%") %>%
    width(j = "Métrica", width = 2.0) %>%
    width(j = c("dbCAZy", "Metals", "CAMPER", "Methyl", "Sulfur", "FeGenie"), width = 1.5) %>%
    width(j = "KEGG KOfam", width = 1.8) %>%
    width(j = "UniRef (informativo)", width = 1.8) %>%
    padding(padding.left = 12, padding.right = 12, part = "all") %>%
    padding(padding.top = 6, padding.bottom = 6, part = "all") %>%
    hrule(rule = "exact", part = "body") %>%
    height(height = 0.35, part = "body") %>%
    height(height = 0.4, part = "header")
  save_as_image(ft2,
    path = file.path(OUTPUT_DIR, "tabla_02_dram2_por_database.png"),
    zoom = 5, expand = 25, webshot = "webshot2"
  )
  save_as_pptx(ft2, path = file.path(OUTPUT_DIR, "tabla_02_dram2_por_database.pptx"))
  cat("[OK] Tabla 2 guardada\n")
}

# TABLA 3: DRAM2 vertical
cat("\nGenerando Tabla 3: DRAM2 vertical...\n")

if (.has_flextable) {
  ft3 <- flextable(datos_dram2) %>%
    set_header_labels(
      Base_de_datos = "Base de datos",
      Hits = "Hits absolutos",
      Total_genes = "Total de genes",
      Pct_anotacion = "% de anotación"
    ) %>%
    bold(part = "header") %>%
    align(align = "center", part = "header") %>%
    bg(bg = "#4472C4", part = "header") %>%
    color(color = "white", part = "header") %>%
    fontsize(size = 11, part = "header") %>%
    align(j = "Base_de_datos", align = "left", part = "body") %>%
    align(j = c("Hits", "Total_genes", "Pct_anotacion"), align = "center", part = "body") %>%
    fontsize(size = 11, part = "body") %>%
    colformat_double(j = "Pct_anotacion", digits = 2, suffix = "%") %>%
    colformat_double(j = c("Hits", "Total_genes"), digits = 0, big.mark = ",") %>%
    theme_booktabs() %>%
    border_outer(border = fp_border(color = "black", width = 2), part = "all") %>%
    border_inner_h(border = fp_border(color = "gray70", width = 0.5), part = "body") %>%
    border_inner_v(border = fp_border(color = "gray70", width = 0.5), part = "body") %>%
    width(j = "Base_de_datos", width = 3.0) %>%
    width(j = "Hits", width = 2.0) %>%
    width(j = "Total_genes", width = 2.0) %>%
    width(j = "Pct_anotacion", width = 2.0) %>%
    padding(padding.left = 10, padding.right = 10, part = "all") %>%
    padding(padding.top = 6, padding.bottom = 6, part = "all") %>%
    hrule(rule = "exact", part = "body") %>%
    height(height = 0.35, part = "body") %>%
    height(height = 0.4, part = "header")
  save_as_image(ft3,
    path = file.path(OUTPUT_DIR, "tabla_03_dram2_vertical.png"),
    zoom = 5, expand = 20, webshot = "webshot2"
  )
  save_as_pptx(ft3, path = file.path(OUTPUT_DIR, "tabla_03_dram2_vertical.pptx"))
  cat("[OK] Tabla 3 guardada\n")
}

# Table 4: taxonomy
cat("\nGenerando Tabla 4: Taxonomía...\n")

datos_taxonomia <- data.frame(
  Especie = c(
    "JAKUTA01 sp950068635", "MGIIb-O1 sp002498525", "N/A", "MGIIb-Q1 sp002505495",
    "N/A", "Thalassarchaeum sp002507125", "Thalassarchaeum betae", "N/A",
    "UBA57 sp002494565", "UBA60 sp002503395", "UBA8886 sp003193815",
    "UBA8886 sp937893445", "Undinarchaeum marinum"
  ),
  Frec_Especie = c(2, 7, 5, 2, 3, 13, 3, 2, 4, 4, 5, 1, 3),
  Color_display = c("", "", "", "", "", "", "", "", "", "", "", "", ""),
  Genero = c(
    "JAKUTA01", "MGIIb-O1", "MGIIb-O1", "MGIIb-O1", "Nitrosopelagicus",
    "Thalassarchaeum", "Thalassarchaeum", "UBA36", "UBA57", "UBA60",
    "UBA8886", "UBA8886", "Undinarchaeum"
  ),
  Frec_Genero = c(2, 14, 14, 14, 3, 16, 16, 2, 4, 4, 6, 6, 3),
  stringsAsFactors = FALSE
)

if (.has_flextable) {
  ft4 <- flextable(datos_taxonomia) %>%
    set_header_labels(
      Especie = "Especie", Frec_Especie = "Frecuencia", Color_display = "Color",
      Genero = "Género", Frec_Genero = "Frecuencia"
    ) %>%
    bold(part = "header") %>%
    align(align = "center", part = "header") %>%
    bg(bg = "black", part = "header") %>%
    color(color = "white", part = "header") %>%
    fontsize(size = 11, part = "header") %>%
    align(j = c("Especie"), align = "left", part = "body") %>%
    align(j = c("Frec_Especie", "Color_display", "Frec_Genero"), align = "center", part = "body") %>%
    align(j = "Genero", align = "left", part = "body") %>%
    fontsize(size = 10, part = "body") %>%
    merge_v(j = c("Genero", "Color_display", "Frec_Genero")) %>%
    valign(j = c("Genero", "Color_display", "Frec_Genero"), valign = "center", part = "body") %>%
    border_outer(border = fp_border(color = "black", width = 2), part = "all") %>%
    border_inner_h(border = fp_border(color = "gray70", width = 0.5), part = "body") %>%
    border_inner_v(border = fp_border(color = "gray70", width = 0.5), part = "all") %>%
    bg(i = 1, j = "Color_display", bg = "#8B0000", part = "body") %>%
    bg(i = 2, j = "Color_display", bg = "#228B22", part = "body") %>%
    bg(i = 3, j = "Color_display", bg = "#228B22", part = "body") %>%
    bg(i = 4, j = "Color_display", bg = "#228B22", part = "body") %>%
    bg(i = 5, j = "Color_display", bg = "#00BFFF", part = "body") %>%
    bg(i = 6, j = "Color_display", bg = "#8B008B", part = "body") %>%
    bg(i = 7, j = "Color_display", bg = "#8B008B", part = "body") %>%
    bg(i = 8, j = "Color_display", bg = "#FFA500", part = "body") %>%
    bg(i = 9, j = "Color_display", bg = "#FF0000", part = "body") %>%
    bg(i = 10, j = "Color_display", bg = "#FFD700", part = "body") %>%
    bg(i = 11, j = "Color_display", bg = "#87CEEB", part = "body") %>%
    bg(i = 12, j = "Color_display", bg = "#87CEEB", part = "body") %>%
    bg(i = 13, j = "Color_display", bg = "#DC143C", part = "body") %>%
    width(j = "Especie", width = 2.2) %>%
    width(j = "Frec_Especie", width = 1.2) %>%
    width(j = "Color_display", width = 0.8) %>%
    width(j = "Genero", width = 2.0) %>%
    width(j = "Frec_Genero", width = 1.2) %>%
    padding(padding.left = 8, padding.right = 8, part = "all") %>%
    padding(padding.top = 4, padding.bottom = 4, part = "all") %>%
    hrule(rule = "exact", part = "body") %>%
    height(height = 0.3, part = "body")
  save_as_image(ft4,
    path = file.path(OUTPUT_DIR, "tabla_04_taxonomia.png"),
    zoom = 5, expand = 20, webshot = "webshot2"
  )
  save_as_pptx(ft4, path = file.path(OUTPUT_DIR, "tabla_04_taxonomia.pptx"))
  cat("[OK] Tabla 4 guardada\n")
}

# Guardar datos CSV
write_csv(datos_anot, file.path(OUTPUT_DIR, "tabla_01_datos_comparacion_software.csv"))
write_csv(datos_dram2, file.path(OUTPUT_DIR, "tabla_02_datos_dram2_detalle.csv"))
write_csv(datos_dram2_transpuesto, file.path(OUTPUT_DIR, "tabla_02_datos_dram2_transpuesto.csv"))
write_csv(datos_taxonomia, file.path(OUTPUT_DIR, "tabla_04_datos_taxonomia.csv"))
cat("[OK] Datos de tablas guardados en CSV\n")

# ---- 4. PREPARACIÓN DE MATRICES PARA HEATMAPS ----
cat("\n=== PREPARANDO MATRICES PARA HEATMAPS ===\n")

heatmap_data <- summary_by_mag %>%
  select(MAG, dbCAZy, Metals, CAMPER, Methyl, KEGG_KOfam, KEGG_DB, Sulfur, FeGenie, UniRef_Informative) %>%
  column_to_rownames("MAG")

heatmap_percent <- summary_by_mag %>%
  mutate(
    dbCAZy_pct = (dbCAZy / Total_Genes) * 100,
    Metals_pct = (Metals / Total_Genes) * 100,
    CAMPER_pct = (CAMPER / Total_Genes) * 100,
    Methyl_pct = (Methyl / Total_Genes) * 100,
    KEGG_KOfam_pct = (KEGG_KOfam / Total_Genes) * 100,
    KEGG_DB_pct = (KEGG_DB / Total_Genes) * 100,
    Sulfur_pct = (Sulfur / Total_Genes) * 100,
    FeGenie_pct = (FeGenie / Total_Genes) * 100,
    UniRef_Informative_pct = (UniRef_Informative / Total_Genes) * 100
  ) %>%
  select(MAG, ends_with("_pct")) %>%
  column_to_rownames("MAG")

colnames(heatmap_percent) <- gsub("_pct", "", colnames(heatmap_percent))

annotation_row <- data.frame(
  Profundidad = get_depth(rownames(heatmap_data)),
  row.names = rownames(heatmap_data)
)

depth_levels <- c("ZMOR", "1000 m", "Fondo")
colors_palette <- RColorBrewer::brewer.pal(max(3, length(depth_levels)), "Set2")
depth_colors <- setNames(colors_palette[1:length(depth_levels)], depth_levels)

annotation_colors <- list(Profundidad = depth_colors)

rownames(heatmap_data) <- convert_to_display_name(rownames(heatmap_data))
rownames(heatmap_percent) <- convert_to_display_name(rownames(heatmap_percent))
rownames(annotation_row) <- convert_to_display_name(rownames(annotation_row))

ordered_data <- order_matrix_by_depth(heatmap_data, annotation_row)
heatmap_data_ordered <- ordered_data$matrix
annotation_row_ordered <- ordered_data$annotation

ordered_percent <- order_matrix_by_depth(heatmap_percent, annotation_row)
heatmap_percent_ordered <- ordered_percent$matrix
annotation_row_percent_ordered <- ordered_percent$annotation

# ---- 5. HEATMAPS GENERALES ----
cat("\n=== GENERANDO HEATMAPS GENERALES ===\n")

cat("\nGenerando heatmap general de bases de datos...\n")

db_summary_matrix <- summary_table %>%
  select(dbCAZy, Metals, CAMPER, Methyl, KEGG_KOfam, KEGG_DB, Sulfur, FeGenie, UniRef_Informative) %>%
  t() %>%
  as.data.frame()

colnames(db_summary_matrix) <- "Total_Hits"
db_summary_matrix <- as.matrix(db_summary_matrix)

col_fun_01 <- colorRamp2(
  seq(0, max(db_summary_matrix, na.rm = TRUE), length.out = 6),
  c("white", "yellow", "orange", "orangered", "red", "darkred")
)
ht_01 <- Heatmap(
  db_summary_matrix,
  name = "Hits",
  col = col_fun_01,
  cluster_rows = FALSE, cluster_columns = FALSE,
  show_row_names = TRUE, show_column_names = TRUE,
  row_names_gp = gpar_f(fontsize = 12),
  column_names_gp = gpar_f(fontsize = 14),
  column_names_rot = 0,
  rect_gp = gpar_f(col = "black", lwd = 0.5),
  border = TRUE,
  width = unit(3, "cm"),
  height = unit(nrow(db_summary_matrix) * 1, "cm"),
  column_title = "Cobertura total de hits por base de datos\n(Solo anotaciones de alta calidad)",
  column_title_gp = gpar_f(fontsize = 11, fontface = "bold"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.0f", db_summary_matrix[i, j]), x, y,
      gp = gpar_f(fontsize = 11, col = "black", fontface = "bold")
    )
  },
  heatmap_legend_param = list(
    title = "Hits",
    title_gp = gpar_f(fontsize = 10, fontface = "bold"),
    labels_gp = gpar_f(fontsize = 9),
    border = "black", legend_direction = "vertical"
  )
)

png(file.path(OUTPUT_DIR, "heatmap_01_database_hits.png"),
  width = 10, height = 6, units = "in", res = 300
)
draw(ht_01)
dev.off()
cat("[OK] Heatmap general guardado\n")

cat("\nGenerando heatmap de conteos absolutos...\n")

mat_02 <- as.matrix(heatmap_data_ordered)
col_fun_02 <- colorRamp2(
  seq(0, max(mat_02, na.rm = TRUE), length.out = 6),
  c("white", "yellow", "orange", "orangered", "red", "darkred")
)
ra_02 <- rowAnnotation(
  Profundidad = annotation_row_ordered$Profundidad,
  col = list(Profundidad = depth_colors),
  annotation_name_side = "top",
  simple_anno_size = unit(0.4, "cm"),
  annotation_name_gp = gpar_f(fontsize = 10, fontface = "bold")
)
ht_02 <- Heatmap(
  mat_02,
  name = "Hits",
  col = col_fun_02,
  cluster_rows = FALSE, cluster_columns = TRUE,
  show_row_names = TRUE, show_column_names = TRUE,
  row_names_gp = gpar_f(fontsize = 8),
  column_names_gp = gpar_f(fontsize = 10),
  column_names_rot = 45,
  left_annotation = ra_02,
  rect_gp = gpar_f(col = "black", lwd = 0.5),
  border = TRUE,
  column_title = "Hits de anotación funcional por MAG (alta calidad)",
  column_title_gp = gpar_f(fontsize = 12, fontface = "bold"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.0f", mat_02[i, j]), x, y,
      gp = gpar_f(fontsize = 7, col = "black")
    )
  },
  heatmap_legend_param = list(
    title = "Hits",
    title_gp = gpar_f(fontsize = 10, fontface = "bold"),
    labels_gp = gpar_f(fontsize = 9),
    border = "black", legend_direction = "vertical"
  )
)

png(file.path(OUTPUT_DIR, "heatmap_02_absolute_counts.png"),
  width = 12, height = 12, units = "in", res = 300
)
draw(ht_02)
dev.off()
cat("[OK] Heatmap de conteos absolutos guardado\n")

cat("\nGenerando heatmap de porcentajes...\n")

mat_03 <- as.matrix(heatmap_percent_ordered)
col_fun_03 <- colorRamp2(
  seq(0, max(mat_03, na.rm = TRUE), length.out = 6),
  c("white", "yellow", "orange", "orangered", "red", "darkred")
)
ra_03 <- rowAnnotation(
  Profundidad = annotation_row_percent_ordered$Profundidad,
  col = list(Profundidad = depth_colors),
  annotation_name_side = "top",
  simple_anno_size = unit(0.4, "cm"),
  annotation_name_gp = gpar_f(fontsize = 10, fontface = "bold")
)
ht_03 <- Heatmap(
  mat_03,
  name = "%",
  col = col_fun_03,
  cluster_rows = FALSE, cluster_columns = TRUE,
  show_row_names = TRUE, show_column_names = TRUE,
  row_names_gp = gpar_f(fontsize = 8),
  column_names_gp = gpar_f(fontsize = 10),
  column_names_rot = 45,
  left_annotation = ra_03,
  rect_gp = gpar_f(col = "black", lwd = 0.5),
  border = TRUE,
  column_title = "Cobertura de anotación funcional por MAG (% - alta calidad)",
  column_title_gp = gpar_f(fontsize = 12, fontface = "bold"),
  cell_fun = function(j, i, x, y, width, height, fill) {
    grid.text(sprintf("%.1f", mat_03[i, j]), x, y,
      gp = gpar_f(fontsize = 7, col = "black")
    )
  },
  heatmap_legend_param = list(
    title = "%",
    title_gp = gpar_f(fontsize = 10, fontface = "bold"),
    labels_gp = gpar_f(fontsize = 9),
    border = "black", legend_direction = "vertical"
  )
)

png(file.path(OUTPUT_DIR, "heatmap_03_percentages.png"),
  width = 12, height = 12, units = "in", res = 300
)
draw(ht_03)
dev.off()
cat("[OK] Heatmap de porcentajes guardado\n")

# ---- 6. HEATMAPS INDIVIDUALES CON DESCRIPCIONES FUNCIONALES ----
cat("\n=== GENERANDO HEATMAPS INDIVIDUALES CON DESCRIPCIONES ===\n")

generate_db_heatmap <- function(data, db_col, db_name, output_file, title, use_kegg_api = TRUE) {
  db_data <- data %>%
    filter(!is.na(!!sym(db_col))) %>%
    select(MAG = input_fasta, !!sym(db_col)) %>%
    separate_rows(!!sym(db_col), sep = ";\\s*") %>%
    filter(!!sym(db_col) != "") %>%
    count(MAG, !!sym(db_col))

  if (nrow(db_data) < 5) {
    cat("  [WARN] No hay suficientes datos para", db_name, "\n")
    return(NULL)
  }

  db_matrix <- db_data %>%
    pivot_wider(names_from = !!sym(db_col), values_from = n, values_fill = 0) %>%
    column_to_rownames("MAG")

  db_matrix_filtered <- db_matrix[, colSums(db_matrix > 0) >= 2]

  if (ncol(db_matrix_filtered) < 3) {
    cat("  [WARN] No hay suficientes", db_name, "presentes en múltiples MAGs\n")
    return(NULL)
  }

  if (ncol(db_matrix_filtered) > 50) {
    top_ids <- names(sort(colSums(db_matrix_filtered), decreasing = TRUE)[1:50])
    db_matrix_filtered <- db_matrix_filtered[, top_ids]
  }

  # Protein IDs are shown bare. The block below would turn them into
  # "ID: description" labels; it is left off because the descriptions are
  # uninformative for most of these assignments.
  cat("  Generando labels funcionales para", db_name, "...\n")
  original_ids <- colnames(db_matrix_filtered)
  functional_labels <- create_functional_labels(original_ids, db_col, data,
    max_length = 40, use_kegg_api = use_kegg_api
  )
  colnames(db_matrix_filtered) <- functional_labels

  annotation_row_db <- data.frame(
    Profundidad = get_depth(rownames(db_matrix_filtered)),
    row.names = rownames(db_matrix_filtered)
  )

  rownames(db_matrix_filtered) <- convert_to_display_name(rownames(db_matrix_filtered))
  rownames(annotation_row_db) <- convert_to_display_name(rownames(annotation_row_db))

  ordered_db <- order_matrix_by_depth(db_matrix_filtered, annotation_row_db)

  mat_db <- as.matrix(ordered_db$matrix)
  col_fun_db <- colorRamp2(
    seq(0, max(mat_db, na.rm = TRUE), length.out = 6),
    c("white", "yellow", "orange", "orangered", "red", "darkred")
  )
  ra_db <- rowAnnotation(
    Profundidad = ordered_db$annotation$Profundidad,
    col = list(Profundidad = depth_colors),
    annotation_name_side = "top",
    simple_anno_size = unit(0.4, "cm"),
    annotation_name_gp = gpar_f(fontsize = 8, fontface = "bold")
  )
  ht_db <- Heatmap(
    mat_db,
    name = "Hits",
    col = col_fun_db,
    cluster_rows = FALSE, cluster_columns = TRUE,
    show_row_names = TRUE, show_column_names = TRUE,
    row_names_gp = gpar_f(fontsize = 7),
    column_names_gp = gpar_f(fontsize = 6),
    column_names_rot = 90,
    left_annotation = ra_db,
    rect_gp = gpar_f(col = "black", lwd = 0.5),
    border = TRUE,
    column_title = title,
    column_title_gp = gpar_f(fontsize = 8, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Hits",
      title_gp = gpar_f(fontsize = 8, fontface = "bold"),
      labels_gp = gpar_f(fontsize = 7),
      border = "black", legend_direction = "vertical"
    )
  )

  png(output_file, width = 20, height = 12, units = "in", res = 300)
  draw(ht_db)
  dev.off()
  cat("  [OK]", db_name, "heatmap guardado\n")
}

cat("[OK] Sección de heatmaps individuales por base de datos omitida (eliminada por diseño)\n")

# ---- 7. GRÁFICOS GENERALES ----
cat("\n=== GENERANDO GRÁFICOS GENERALES ===\n")

cat("\nGenerando gráfico de calidad de anotación...\n")

quality_comparison <- summary_table %>%
  select(Curated_Functional_Annotation, Homology_Only, No_Annotation) %>%
  pivot_longer(everything(), names_to = "Category", values_to = "Count") %>%
  mutate(
    Category = factor(Category, levels = c("Curated_Functional_Annotation", "Homology_Only", "No_Annotation")),
    Category_Label = case_when(
      Category == "Curated_Functional_Annotation" ~ "Anotación funcional curada\n(KEGG/HMM)",
      Category == "Homology_Only" ~ "Sólo homología\n(UniRef, sin función curada)",
      Category == "No_Annotation" ~ "Sin anotación"
    ),
    Category_Label = factor(Category_Label, levels = c(
      "Anotación funcional curada\n(KEGG/HMM)",
      "Sólo homología\n(UniRef, sin función curada)",
      "Sin anotación"
    ))
  )

p_quality <- ggplot(quality_comparison, aes(x = Category_Label, y = Count, fill = Category_Label)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8, width = 0.45) +
  geom_text(aes(label = paste0(comma(Count), "\n(", round(Count / sum(Count) * 100, 1), "%)")),
    vjust = -0.3, size = 5, fontface = "bold"
  ) +
  scale_fill_manual(values = c(
    "Anotación funcional curada\n(KEGG/HMM)" = "#2E7D32",
    "Sólo homología\n(UniRef, sin función curada)" = "#FBC02D",
    "Sin anotación" = "#D32F2F"
  )) +
  theme_publication(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12, face = "bold"),
    legend.position = "none"
  ) +
  labs(
    subtitle = paste("Total de genes:", comma(summary_table$Total_Genes)),
    x = "Categoría de anotación",
    y = "Cantidad de genes"
  ) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15)))

ggsave(file.path(OUTPUT_DIR, "plot_01_annotation_quality_comparison.png"),
  p_quality,
  width = 12, height = 8, dpi = 300
)

cat("[OK] Gráfico de calidad guardado\n")

cat("\nGenerando gráfico de distribución de hits...\n")

distribution_data <- summary_table %>%
  select(dbCAZy, Metals, CAMPER, Methyl, KEGG_KOfam, KEGG_DB, Sulfur, FeGenie, UniRef_Informative) %>%
  pivot_longer(everything(), names_to = "Database", values_to = "Hits") %>%
  mutate(
    Database = factor(Database, levels = Database[order(-Hits)])
  )

p1 <- ggplot(distribution_data, aes(x = Database, y = Hits, fill = Database)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8) +
  geom_text(aes(label = comma(Hits)), vjust = -0.5, size = 4) +
  scale_fill_brewer(palette = "Set3") +
  theme_publication(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    title = "Hits de alta calidad por base de datos de referencia",
    subtitle = paste("Total de genes analizados:", comma(nrow(dram))),
    x = "Base de datos de referencia",
    y = "Cantidad de hits"
  ) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.1)))

ggsave(file.path(OUTPUT_DIR, "plot_02_hits_distribution.png"),
  p1,
  width = 11, height = 6, dpi = 300
)

cat("[OK] Gráfico de distribución guardado\n")

cat("\nGenerando gráfico de anotación por herramienta...\n")

p_anot_tool <- ggplot(
  datos_anot,
  aes(x = reorder(software, -hits), y = hits, fill = software)
) +
  geom_bar(stat = "identity", color = "black", alpha = 0.8, width = 0.45) +
  geom_text(
    aes(label = paste0(comma(hits), "\n(", round(pct_anotacion, 1), "%)")),
    vjust = -0.3, size = 4
  ) +
  scale_fill_brewer(palette = "Set3") +
  theme_publication(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none"
  ) +
  labs(
    subtitle = paste("Total de genes analizados:", comma(datos_anot$genes[1])),
    x = "Herramienta",
    y = "Cantidad de genes anotados"
  ) +
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.15)))

ggsave(file.path(OUTPUT_DIR, "plot_03_annotation_by_tool.png"),
  p_anot_tool,
  width = 8, height = 6, dpi = 300
)

cat("[OK] Gráfico de anotación por herramienta guardado\n")

cat("\nGenerando gráfico de tasa de anotación por MAG...\n")

summary_by_mag_display <- summary_by_mag %>%
  mutate(DisplayName = convert_to_display_name(MAG)) %>%
  arrange(Percent_Curated_Functional) %>%
  mutate(DisplayName = factor(DisplayName, levels = DisplayName))

p2 <- ggplot(summary_by_mag_display, aes(x = DisplayName, y = Percent_Curated_Functional)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black", width = 0.45) +
  geom_hline(
    yintercept = mean(summary_by_mag_display$Percent_Curated_Functional),
    linetype = "dashed", color = "red", linewidth = 1
  ) +
  coord_flip() +
  theme_publication(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    subtitle = paste0("Promedio: ", round(mean(summary_by_mag_display$Percent_Curated_Functional), 1), "%"),
    x = "MAG",
    y = "Genes con anotación funcional curada (%)"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05)))

ggsave(file.path(OUTPUT_DIR, "plot_03_annotation_rate_per_MAG.png"),
  p2,
  width = 10, height = 14, dpi = 300
)

cat("[OK] Gráfico de tasa de anotación guardado\n")

# -----------------------------------------------------------------------------
# Mean annotation rate by phylum and by order
# -----------------------------------------------------------------------------
cat("\nGenerando gráficos de tasa de anotación por filo y por orden...\n")

.extract_rank <- function(taxonomy, prefix) {
  sapply(taxonomy, function(tx) {
    if (is.na(tx)) return(NA_character_)
    m <- regmatches(tx, regexpr(paste0(prefix, "[^;]+"), tx))
    if (length(m) == 0) return(NA_character_)
    val <- sub(paste0("^", prefix), "", m)
    if (nchar(trimws(val)) == 0) NA_character_ else val
  }, USE.NAMES = FALSE)
}

summary_by_mag_tax <- summary_by_mag %>%
  mutate(
    Phylum = .extract_rank(FullTaxonomy, "p__"),
    Order  = .extract_rank(FullTaxonomy, "o__")
  )

.plot_pct_by_taxon <- function(df, taxon_col, taxon_label, outfile,
                               label_formatter = identity,
                               add_significance = FALSE) {
  agg <- df %>%
    filter(!is.na(.data[[taxon_col]])) %>%
    group_by(Taxon = .data[[taxon_col]]) %>%
    summarise(
      n = n(),
      mean_pct = mean(Percent_Curated_Functional, na.rm = TRUE),
      sd_pct = if (n() > 1) sd(Percent_Curated_Functional, na.rm = TRUE) else 0,
      .groups = "drop"
    ) %>%
    mutate(
      TaxonLabel = vapply(Taxon, label_formatter, character(1)),
      TaxonLabel = paste0(TaxonLabel, " (n=", n, ")")
    ) %>%
    arrange(mean_pct) %>%
    mutate(TaxonLabel = factor(TaxonLabel, levels = TaxonLabel))

  global_mean <- mean(df$Percent_Curated_Functional, na.rm = TRUE)

  # The test runs on the raw per-MAG values, one observation per MAG.
  # Kruskal-Wallis is the right omnibus here: it assumes neither normality nor
  # equal variance, and it tolerates the badly unbalanced group sizes, which run
  # from n = 3 to n = 38. It only runs where at least two groups hold at least
  # two MAGs.
  sig_info <- NULL
  if (add_significance) {
    test_df <- df %>%
      filter(!is.na(.data[[taxon_col]]), !is.na(Percent_Curated_Functional)) %>%
      mutate(.grp = as.factor(.data[[taxon_col]]))
    grp_sizes <- table(test_df$.grp)
    usable_grps <- names(grp_sizes)[grp_sizes >= 2]
    if (length(usable_grps) >= 2) {
      test_df <- test_df %>% filter(.grp %in% usable_grps) %>% droplevels()
      kw <- kruskal.test(Percent_Curated_Functional ~ .grp, data = test_df)
      p_val <- unname(kw$p.value)
      sig_label <- if (is.na(p_val)) "ns"
                   else if (p_val < 0.001) "***"
                   else if (p_val < 0.01)  "**"
                   else if (p_val < 0.05)  "*"
                   else                    "ns"
      sig_info <- list(
        p = p_val, label = sig_label,
        chi = unname(kw$statistic), df = unname(kw$parameter),
        n_groups = length(usable_grps),
        excluded = setdiff(names(grp_sizes), usable_grps)
      )
      cat(sprintf(
        "  [%s] Kruskal-Wallis: χ²(%d) = %.3f, p = %.4g  to  %s%s\n",
        taxon_label, sig_info$df, sig_info$chi, sig_info$p, sig_info$label,
        if (length(sig_info$excluded) > 0)
          paste0("  (excluidos por n<2: ", paste(sig_info$excluded, collapse = ", "), ")")
        else ""
      ))
    }
  }

  y_max <- max(agg$mean_pct + agg$sd_pct, na.rm = TRUE)
  # A significant result needs room on the right, past the flipped end, for the
  # bracket and its asterisk.
  y_limit_mult <- if (!is.null(sig_info) && sig_info$label != "ns") 1.32 else 1.18

  p <- ggplot(agg, aes(x = TaxonLabel, y = mean_pct)) +
    geom_bar(stat = "identity", fill = "steelblue", color = "black", width = 0.6) +
    geom_errorbar(aes(ymin = pmax(0, mean_pct - sd_pct), ymax = mean_pct + sd_pct),
                  width = 0.2, color = "gray30") +
    # The label sits at the far right of the error bar, clear of the SD bar.
    geom_text(aes(y = mean_pct + sd_pct, label = sprintf("%.1f%%", mean_pct)),
              hjust = -0.15, size = 3.5, family = font_family) +
    geom_hline(yintercept = global_mean, linetype = "dashed",
               color = "red", linewidth = 0.8) +
    coord_flip() +
    theme_publication(base_size = 12) +
    theme(axis.text.y = element_text(size = 10)) +
    labs(
      subtitle = paste0("Promedio global: ", round(global_mean, 1), "%"),
      x = taxon_label,
      y = "Genes con anotación funcional curada (%)"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.18)),
      limits = c(0, y_max * y_limit_mult)
    )

  # Significance drawn as a bracket spanning every group with one asterisk
  # above it, which is the convention for an omnibus result. It reads as
  # horizontal here and vertical after coord_flip.
  if (!is.null(sig_info) && sig_info$label != "ns") {
    n_taxa <- nrow(agg)
    y_bracket <- y_max * 1.20
    tick_len  <- y_max * 0.015
    y_text    <- y_bracket + y_max * 0.035
    p <- p +
      annotate("segment", x = 1, xend = n_taxa,
               y = y_bracket, yend = y_bracket,
               linewidth = 0.5, color = "black") +
      annotate("segment", x = 1, xend = 1,
               y = y_bracket - tick_len, yend = y_bracket,
               linewidth = 0.5, color = "black") +
      annotate("segment", x = n_taxa, xend = n_taxa,
               y = y_bracket - tick_len, yend = y_bracket,
               linewidth = 0.5, color = "black") +
      annotate("text", x = (1 + n_taxa) / 2, y = y_text,
               label = sig_info$label,
               size = 6, fontface = "bold", family = font_family,
               vjust = 0.5)
  }

  n_taxa <- nrow(agg)
  ggsave(outfile, p,
    width = 10, height = max(3, 0.55 * n_taxa + 2), dpi = 300
  )
  invisible(list(agg = agg, sig = sig_info))
}

res_phylum <- .plot_pct_by_taxon(
  summary_by_mag_tax, "Phylum", "Filo",
  file.path(OUTPUT_DIR, "plot_03b_annotation_rate_per_phylum.png")
)
agg_phylum <- res_phylum$agg
write_csv(agg_phylum,
  file.path(OUTPUT_DIR, "plot_03b_annotation_rate_per_phylum.csv")
)

res_order <- .plot_pct_by_taxon(
  summary_by_mag_tax, "Order", "Orden",
  file.path(OUTPUT_DIR, "plot_03c_annotation_rate_per_order.png"),
  label_formatter = fmt_order_label,
  add_significance = TRUE
)
agg_order <- res_order$agg
write_csv(agg_order,
  file.path(OUTPUT_DIR, "plot_03c_annotation_rate_per_order.csv")
)
# Omnibus results: Kruskal-Wallis on Percent_Curated_Functional against order
if (!is.null(res_order$sig)) {
  sig <- res_order$sig
  # Effect size for the omnibus: epsilon squared = H / (N - 1), Tomczak &
  # Tomczak 2014. It reads as the share of rank variance the grouping explains.
  # Conventionally: under 0.04 very small, 0.04 to 0.16 small,
  # 0.16-0.36 medio, ≥0.36 grande.
  test_df_omni <- summary_by_mag_tax %>%
    filter(!is.na(Order), !is.na(Percent_Curated_Functional)) %>%
    mutate(.grp = as.factor(Order))
  grp_sizes <- table(test_df_omni$.grp)
  test_df_omni <- test_df_omni %>%
    filter(.grp %in% names(grp_sizes)[grp_sizes >= 2]) %>%
    droplevels()
  N_omni <- nrow(test_df_omni)
  epsilon_sq <- sig$chi / (N_omni - 1)
  eps_interp <- {
    if (epsilon_sq < 0.04) "muy pequeño"
    else if (epsilon_sq < 0.16) "pequeño"
    else if (epsilon_sq < 0.36) "medio" else "grande"
  }
  writeLines(
    c(
      "Kruskal-Wallis omnibus test — % anotación funcional curada (KEGG/HMM, sin UniRef) por Orden",
      sprintf("  H0: la distribución del %% de anotación funcional curada es igual entre órdenes."),
      sprintf("  N total de MAGs (con n≥2 por grupo): %d", N_omni),
      sprintf("  Grupos evaluados (n≥2 MAGs): %d", sig$n_groups),
      if (length(sig$excluded) > 0)
        sprintf("  Excluidos por n<2: %s", paste(sig$excluded, collapse = ", "))
      else "  Excluidos por n<2: ninguno",
      sprintf("  Estadístico χ²(df=%d) = %.4f", sig$df, sig$chi),
      sprintf("  p-valor = %.4g", sig$p),
      sprintf("  Tamaño de efecto ε² = %.4f (%s)", epsilon_sq, eps_interp),
      sprintf("  Convención: *** p<0.001, ** p<0.01, * p<0.05, ns p≥0.05  to  %s", sig$label)
    ),
    file.path(OUTPUT_DIR, "plot_03c_annotation_rate_per_order_KWtest.txt")
  )
}

# -----------------------------------------------------------------------------
# Dunn post-hoc (1964) with the tie correction from Zar (1999), reported under
# both BH (primary) and Bonferroni (secondary). It runs only where the
# Kruskal-Wallis omnibus reached p < 0.05 and at least three groups hold two or
# more MAGs, which are the conditions under which a post-hoc means anything.
# -----------------------------------------------------------------------------
dunn_posthoc <- function(values, groups) {
  df <- data.frame(x = values, g = as.factor(groups))
  df <- df[complete.cases(df), ]
  df$g <- droplevels(df$g)

  N <- nrow(df)
  df$rk <- rank(df$x)

  # Tie correction: C = sum(t^3 - t) / [12(N - 1)]
  tab <- table(df$x)
  ties <- tab[tab > 1]
  C <- if (length(ties) > 0) sum(ties^3 - ties) / (12 * (N - 1)) else 0
  sigma_sq_base <- (N * (N + 1) / 12) - C

  grp_n  <- tapply(df$rk, df$g, length)
  grp_mr <- tapply(df$rk, df$g, mean)

  lvl <- levels(df$g)
  pairs <- t(combn(lvl, 2))
  res <- data.frame(
    Grupo_1     = pairs[, 1],
    Grupo_2     = pairs[, 2],
    n_1         = as.integer(grp_n[pairs[, 1]]),
    n_2         = as.integer(grp_n[pairs[, 2]]),
    rango_medio_1 = unname(grp_mr[pairs[, 1]]),
    rango_medio_2 = unname(grp_mr[pairs[, 2]]),
    stringsAsFactors = FALSE
  )
  res$diff_rangos <- res$rango_medio_1 - res$rango_medio_2
  res$SE <- sqrt(sigma_sq_base * (1 / res$n_1 + 1 / res$n_2))
  res$Z  <- res$diff_rangos / res$SE
  res$p_crudo       <- 2 * pnorm(-abs(res$Z))
  res$p_BH          <- p.adjust(res$p_crudo, method = "BH")
  res$p_Bonferroni  <- p.adjust(res$p_crudo, method = "bonferroni")
  # Effect size r = |Z| / sqrt(N), Rosenthal 1991, read against Cohen-like bands.
  res$r_efecto <- abs(res$Z) / sqrt(N)
  res$magnitud <- cut(res$r_efecto,
    breaks = c(-Inf, 0.1, 0.3, 0.5, Inf),
    labels = c("despreciable", "pequeño", "medio", "grande"),
    right = FALSE
  )
  res$signif_BH <- cut(res$p_BH,
    breaks = c(-Inf, 0.001, 0.01, 0.05, Inf),
    labels = c("***", "**", "*", "ns"),
    right = FALSE
  )
  res
}

if (!is.null(res_order$sig) && res_order$sig$p < 0.05 && res_order$sig$n_groups >= 3) {
  cat("\nGenerando post-hoc Dunn para % anotación funcional por Orden...\n")

  dunn_df <- summary_by_mag_tax %>%
    filter(!is.na(Order), !is.na(Percent_Curated_Functional))
  # the same groups the omnibus used, n >= 2
  ord_sizes <- table(dunn_df$Order)
  dunn_df <- dunn_df %>%
    filter(Order %in% names(ord_sizes)[ord_sizes >= 2])

  dunn_res <- dunn_posthoc(
    values = dunn_df$Percent_Curated_Functional,
    groups = dunn_df$Order
  )

  # Order labels, with the class in parentheses after a code-like name
  dunn_res$Grupo_1_fmt <- vapply(dunn_res$Grupo_1, fmt_order_label, character(1))
  dunn_res$Grupo_2_fmt <- vapply(dunn_res$Grupo_2, fmt_order_label, character(1))

  dunn_out <- dunn_res %>%
    mutate(Comparacion = paste(Grupo_1_fmt, "vs", Grupo_2_fmt)) %>%
    select(Comparacion, Grupo_1, Grupo_2, n_1, n_2,
           rango_medio_1, rango_medio_2, diff_rangos, Z,
           p_crudo, p_BH, p_Bonferroni, r_efecto, magnitud, signif_BH) %>%
    arrange(p_BH)

  write_tsv(
    dunn_out,
    file.path(OUTPUT_DIR, "plot_03d_annotation_rate_per_order_DunnPosthoc.tsv")
  )
  cat("[OK] TSV Dunn guardado: plot_03d_annotation_rate_per_order_DunnPosthoc.tsv\n")

  # flextable PNG, formatted for the appendix
  if (.has_flextable) {
    ft_dunn <- dunn_out %>%
      select(Comparacion, n_1, n_2, Z, p_crudo, p_BH, p_Bonferroni,
             r_efecto, magnitud, signif_BH) %>%
      flextable() %>%
      set_header_labels(
        Comparacion  = "Comparación",
        n_1          = "n₁",
        n_2          = "n₂",
        Z            = "Z",
        p_crudo      = "p (sin ajustar)",
        p_BH         = "p ajust. (BH/FDR)",
        p_Bonferroni = "p ajust. (Bonferroni)",
        r_efecto     = "r (efecto)",
        magnitud     = "Magnitud",
        signif_BH    = "Signif. BH"
      ) %>%
      bold(part = "header") %>%
      align(align = "center", part = "header") %>%
      bg(bg = "#4472C4", part = "header") %>%
      color(color = "white", part = "header") %>%
      fontsize(size = 10, part = "header") %>%
      align(j = "Comparacion", align = "left", part = "body") %>%
      align(j = c("n_1", "n_2", "Z", "p_crudo", "p_BH", "p_Bonferroni",
                  "r_efecto", "magnitud", "signif_BH"),
            align = "center", part = "body") %>%
      fontsize(size = 10, part = "body") %>%
      colformat_double(j = c("Z", "r_efecto"), digits = 3) %>%
      colformat_double(j = c("p_crudo", "p_BH", "p_Bonferroni"), digits = 4) %>%
      bold(j = "signif_BH", part = "body") %>%
      bg(j = "signif_BH",
         i = ~ signif_BH != "ns",
         bg = "#E2EFDA", part = "body") %>%
      theme_booktabs() %>%
      border_outer(border = fp_border(color = "black", width = 2), part = "all") %>%
      border_inner_h(border = fp_border(color = "gray70", width = 0.5), part = "body") %>%
      border_inner_v(border = fp_border(color = "gray70", width = 0.5), part = "body") %>%
      add_footer_lines(
        sprintf(
          "Omnibus Kruskal-Wallis: χ²(%d) = %.3f, p = %.4g, ε² = %.3f.  Post-hoc Dunn con corrección de empates (Zar 1999); p ajustados por BH (primaria, FDR) y Bonferroni (FWER). r de efecto = |Z|/√N (Rosenthal 1991). Variable: %% genes con anotación funcional curada (KEGG/KOfam/dbCAZy/CAMPER/sulfur/methyl/FeGenie/metals).",
          res_order$sig$df, res_order$sig$chi, res_order$sig$p,
          res_order$sig$chi / (sum(ord_sizes[ord_sizes >= 2]) - 1)
        )
      ) %>%
      fontsize(size = 8, part = "footer") %>%
      italic(part = "footer") %>%
      padding(padding.left = 8, padding.right = 8, part = "all") %>%
      padding(padding.top = 4, padding.bottom = 4, part = "all") %>%
      autofit()

    save_as_image(ft_dunn,
      path = file.path(OUTPUT_DIR, "plot_03d_annotation_rate_per_order_DunnPosthoc.png"),
      zoom = 5, expand = 20, webshot = "webshot2"
    )
    cat("[OK] PNG Dunn guardado: plot_03d_annotation_rate_per_order_DunnPosthoc.png\n")
  }
}

cat("[OK] Gráficos de tasa de anotación por filo y orden guardados\n")

cat("\nGenerando gráfico de composición funcional...\n")

top_mags <- summary_by_mag %>%
  slice_max(Curated_Functional_Annotation, n = 20) %>%
  pull(MAG)

stacked_data <- summary_by_mag %>%
  filter(MAG %in% top_mags) %>%
  mutate(DisplayName = convert_to_display_name(MAG)) %>%
  select(DisplayName, dbCAZy, Metals, CAMPER, Methyl, KEGG_KOfam, KEGG_DB, Sulfur, FeGenie, UniRef_Informative) %>%
  pivot_longer(-DisplayName, names_to = "Database", values_to = "Hits") %>%
  group_by(DisplayName) %>%
  mutate(Percentage = (Hits / sum(Hits)) * 100) %>%
  ungroup()

p3 <- ggplot(stacked_data, aes(x = DisplayName, y = Percentage, fill = Database)) +
  geom_bar(stat = "identity", color = "white", linewidth = 0.3) +
  scale_fill_brewer(palette = "Set3") +
  coord_flip() +
  theme_publication(base_size = 12) +
  theme(
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    title = "Composición de anotaciones funcionales (Top 20 MAGs)",
    x = "MAG",
    y = "Contribución relativa de la anotación total (%)",
    fill = "Base de datos"
  )

ggsave(file.path(OUTPUT_DIR, "plot_04_functional_composition_stacked.png"),
  p3,
  width = 12, height = 10, dpi = 300
)

cat("[OK] Gráfico de composición funcional guardado\n")

# ---- 8. ANÁLISIS ESPECÍFICOS CON DESCRIPCIONES FUNCIONALES ----
cat("\n=== ANÁLISIS ESPECÍFICOS CON DESCRIPCIONES ===\n")

# KEGG KOs
cat("\nAnalizando KOs de KEGG...\n")

kegg_data <- dram_classified %>%
  filter(!is.na(kofam_id)) %>%
  select(MAG = input_fasta, kofam_id) %>%
  separate_rows(kofam_id, sep = ";\\s*") %>%
  filter(kofam_id != "")

top_kos <- kegg_data %>%
  count(kofam_id, sort = TRUE) %>%
  slice_max(n, n = 30)

# KEGG gives IDs and no description column. Uncommenting the line below fetches
# the descriptions from the API, which takes a while.
top_kos$description <- get_kegg_description(top_kos$kofam_id, use_cache = TRUE, cache_dir = OUTPUT_DIR)
top_kos$label <- trimws(sub("\\s*\\[.*$", "", top_kos$description))
top_kos$label[is.na(top_kos$label) | top_kos$label == ""] <- top_kos$kofam_id[is.na(top_kos$label) | top_kos$label == ""]

p4 <- ggplot(top_kos, aes(x = reorder(label, n), y = n)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  geom_text(aes(label = n), hjust = -0.2, size = 3) +
  coord_flip() +
  theme_publication(base_size = 10) +
  theme(
    axis.text.y = element_text(size = 8)
  ) +
  labs(
    title = "Top 30 KEGG Orthologs más abundantes",
    x = "KEGG Ortholog (KO)",
    y = "Frecuencia"
  ) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.1)))

ggsave(file.path(OUTPUT_DIR, "plot_05_top30_KEGG_KOs.png"),
  p4,
  width = 10, height = 10, dpi = 300
)

write_csv(top_kos, file.path(OUTPUT_DIR, "data_01_top30_KEGG_KOs.csv"))
cat("[OK] Análisis de KOs guardado\n")

# CAZymes
cat("\nAnalizando CAZymes con descripciones...\n")

cazy_data <- dram_classified %>%
  filter(!is.na(dbcan_id)) %>%
  select(MAG = input_fasta, dbcan_id, dbcan_description) %>%
  separate_rows(dbcan_id, dbcan_description, sep = ";\\s*") %>%
  filter(dbcan_id != "")

top_cazymes <- cazy_data %>%
  count(dbcan_id, sort = TRUE) %>%
  slice_max(n, n = 30)

# The functional description is left out: it says little for most CAZy
# assignments, so the bare family ID is shown.
# top_cazymes <- top_cazymes %>%
#   left_join(...) %>%
#   mutate(description = ..., description_es = translate_to_spanish(...))
top_cazymes <- top_cazymes %>%
  mutate(label = dbcan_id)

write_csv(top_cazymes, file.path(OUTPUT_DIR, "data_02_top30_CAZymes.csv"))
cat("[OK] Datos de CAZymes guardados (plot eliminado por diseño)\n")

# UniRef informativo
cat("\nAnalizando hits informativos de UniRef...\n")

if (sum(dram_classified$has_informative_uniref) > 0) {
  uniref_informative_data <- dram_classified %>%
    filter(has_informative_uniref) %>%
    select(MAG = input_fasta, uniref_id, uniref_description) %>%
    separate_rows(uniref_id, uniref_description, sep = ";\\s*") %>%
    filter(uniref_id != "")

  top_uniref <- uniref_informative_data %>%
    count(uniref_id, sort = TRUE) %>%
    slice_max(n, n = 30)

  # The functional description is left out: it says little for most UniRef
  # assignments, so the bare ID is shown.
  # top_uniref <- top_uniref %>%
  #   left_join(...) %>%
  #   mutate(description = ..., description_es = translate_to_spanish(...))
  top_uniref <- top_uniref %>%
    mutate(label = uniref_id)

  write_csv(top_uniref, file.path(OUTPUT_DIR, "data_03_top30_UniRef_informative.csv"))
  cat("[OK] Datos de UniRef informativo guardados (plot eliminado por diseño)\n")
} else {
  cat("[WARN] No hay hits informativos de UniRef\n")
}

# CAMPER
cat("\nAnalizando genes CAMPER con descripciones...\n")

if (sum(dram_classified$has_camper) > 10) {
  camper_data <- dram_classified %>%
    filter(!is.na(camper_id)) %>%
    select(MAG = input_fasta, camper_id, camper_description) %>%
    separate_rows(camper_id, camper_description, sep = ";\\s*") %>%
    filter(camper_id != "")

  top_camper <- camper_data %>%
    count(camper_id, sort = TRUE) %>%
    slice_max(n, n = 30)

  # The functional description is left out: uninformative for most assignments.
  # The bare CAMPER ID is shown.
  # top_camper <- top_camper %>%
  #   left_join(...) %>%
  #   mutate(description = ..., description_es = translate_to_spanish(...))
  top_camper <- top_camper %>%
    mutate(label = camper_id)

  write_csv(top_camper, file.path(OUTPUT_DIR, "data_04_top30_CAMPER.csv"))
  cat("[OK] Datos de CAMPER guardados (plot eliminado por diseño)\n")
}

# Quality against annotation rate
cat("\nGenerando gráfico de calidad vs anotación...\n")

if ("Completeness" %in% colnames(summary_by_mag_display) &&
  "Contamination" %in% colnames(summary_by_mag_display)) {

  p6 <- ggplot(summary_by_mag_display, aes(x = Completeness, y = Percent_Curated_Functional)) +
    geom_point(aes(color = Contamination), alpha = 0.85, size = 4) +
    geom_smooth(method = "lm", se = TRUE, color = "red", linetype = "dashed") +
    scale_color_viridis(option = "magma", direction = -1) +
    theme_publication(base_size = 12) +
    labs(
      title = "Calidad del MAG vs. cobertura de anotación funcional curada",
      x = "Completitud (%)",
      y = "Genes con anotación funcional curada (%)",
      color = "Contaminación (%)"
    )

  ggsave(file.path(OUTPUT_DIR, "plot_09_quality_vs_annotation.png"),
    p6,
    width = 10, height = 7, dpi = 300
  )

  cat("[OK] Gráfico de calidad vs anotación guardado\n")
} else {
  cat("[WARN] Datos de calidad no disponibles\n")
}

# ---- 9. RESUMEN FINAL ----
# ---- 10. KEGG-DECODER Y KEGG-EXPANDER HEATMAPS (ComplexHeatmap) ----
cat("\n=== GENERANDO HEATMAPS DE KEGG-DECODER Y KEGG-EXPANDER ===\n")

DECODER_FILE <- file.path(PROJECT_DIR, "results", "07_metabolic_reconstruction", "keggdecoder_static", "keggdecoder_output.list")
EXPANDER_FILE <- file.path(PROJECT_DIR, "results", "07_metabolic_reconstruction", "keggexpander_output", "HMM_OUT.list")
cat("Las figuras de KEGG se guardarán en OUTPUT_DIR:", OUTPUT_DIR, "\n\n")

# --- Helper functions for KEGG matrices ---

read_kegg_list <- function(file_path) {
  df <- read.delim(file_path, sep = "\t", header = TRUE, check.names = FALSE)
  rownames(df) <- df[[1]]
  df <- df[, -1, drop = FALSE]
  mat <- as.matrix(df)
  mat <- apply(mat, 2, as.numeric)
  rownames(mat) <- rownames(df)
  return(mat)
}

filter_kegg_zero_cols <- function(mat) {
  col_sums <- colSums(mat, na.rm = TRUE)
  non_zero_cols <- col_sums > 0
  n_removed <- sum(!non_zero_cols)
  if (n_removed > 0) cat("  Eliminadas", n_removed, "rutas/funciones con todos ceros\n")
  return(mat[, non_zero_cols, drop = FALSE])
}

load_kegg_metadata <- function(mag_names) {
  metadata_matched <- metadata[match(mag_names, metadata$NameWithoutExtension), ]
  clean_col <- function(col) {
    col[is.na(col) | col == "" | col == "N/A"] <- "Unclassified"
    col
  }
  display_names <- make_unique_names(clean_col(metadata_matched$DisplayName))
  anno_df <- data.frame(
    DisplayName = display_names,
    Station = clean_col(metadata_matched$Station),
    Depth = clean_col(metadata_matched$Depth),
    Genus = clean_col(metadata_matched$Genus),
    stringsAsFactors = FALSE,
    row.names = mag_names
  )
  anno_df$Depth <- factor(anno_df$Depth, levels = c("ZMOR", "1000 m", "Fondo", "Unclassified"))
  return(anno_df)
}

kegg_order_by_depth <- function(mat, metadata_anno) {
  meta <- metadata_anno[rownames(mat), ]
  idx <- order(meta$Depth)
  list(mat = mat[idx, , drop = FALSE], metadata = meta[idx, ])
}

kegg_rename_rows <- function(mat, metadata_anno) {
  dn <- setNames(metadata_anno$DisplayName, rownames(metadata_anno))
  rownames(mat) <- dn[rownames(mat)]
  mat
}

kegg_named_colors <- function(categories, palette_name = "Set1", na_color = "gray90") {
  cats <- as.character(categories)
  cats[is.na(cats)] <- "Unclassified"
  unique_cats <- unique(cats)
  n <- length(unique_cats)
  cols <- RColorBrewer::brewer.pal(max(n, 3), palette_name)[seq_len(n)]
  named <- setNames(cols, unique_cats)
  if ("Unclassified" %in% names(named)) named["Unclassified"] <- na_color
  named
}

# Builds the cladogram filename by inserting _cluster before _ES.png. Defined
# here rather than in section 11 so the KEGG-Decoder figures can use it too.
if (!exists(".cluster_out")) {
  .cluster_out <- function(f) sub("(_ES\\.png)$", "_cluster\\1", f)
}

# ====================== KEGG-Decoder Heatmap =================================
cat("=== Generando KEGG-Decoder Heatmap ===\n")
if (file.exists(DECODER_FILE)) {
  decoder_mat <- read_kegg_list(DECODER_FILE)
  cat("  Dimensiones originales:", dim(decoder_mat), "\n")
  decoder_mat <- filter_kegg_zero_cols(decoder_mat)
  cat("  Dimensiones filtradas:", dim(decoder_mat), "\n\n")

  cache_file_decoder <- file.path(OUTPUT_DIR, "translation_cache_decoder.csv")
  colnames(decoder_mat) <- capitalize_first(translate_to_spanish(colnames(decoder_mat),
    translation_cache_file = cache_file_decoder
  ))

  mag_names_decoder <- rownames(decoder_mat)
  meta_decoder <- load_kegg_metadata(mag_names_decoder)
  decoder_mat <- kegg_rename_rows(decoder_mat, meta_decoder)
  rownames(meta_decoder) <- meta_decoder$DisplayName

  cat("  Ordenando MAGs por profundidad...\n")
  ord_dec <- kegg_order_by_depth(decoder_mat, meta_decoder)
  decoder_mat <- ord_dec$mat
  meta_decoder <- ord_dec$metadata

  depth_col_dec <- kegg_named_colors(meta_decoder$Depth, "Set2")

  col_fun_dec <- colorRamp2(
    c(0, 0.2, 0.4, 0.6, 0.8, 1),
    c("white", "yellow", "orange", "orangered", "red", "darkred")
  )

  ra_dec <- rowAnnotation(
    Profundidad = meta_decoder$Depth,
    col = list(Profundidad = depth_col_dec),
    annotation_name_side = "top",
    show_legend = FALSE,
    simple_anno_size = unit(0.4, "cm"),
    annotation_name_gp = gpar_f(fontsize = 10, fontface = "bold"),
    annotation_legend_param = list(
      Profundidad = list(
        at = na.omit(unique(as.character(meta_decoder$Depth)))
      )
    )
  )

  # Drop columns whose name came back NA from a failed translation
  decoder_mat <- decoder_mat[, !is.na(colnames(decoder_mat)), drop = FALSE]

  ht_dec <- Heatmap(
    decoder_mat,
    name = "Completitud",
    col = col_fun_dec,
    cluster_rows = FALSE, cluster_columns = TRUE,
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "ward.D2",
    show_row_names = TRUE, show_column_names = TRUE,
    show_heatmap_legend = FALSE,
    row_labels = make_taxon_label_expr(rownames(decoder_mat)),
    row_names_gp = gpar_f(fontsize = 10),
    # The longest DisplayName runs about 39 characters, roughly 10 cm at
    # fontsize 10 in italic. ComplexHeatmap defaults to 6 cm, which truncated
    # them at about 23.
    row_names_max_width = unit(11, "cm"),
    column_names_gp = gpar_f(fontsize = 8),
    column_names_rot = 45,
    left_annotation = ra_dec,
    rect_gp = gpar_f(col = "black", lwd = 0.5),
    border = TRUE,
    heatmap_legend_param = list(
      title = "Completitud\nde rutas",
      at = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0%", "25%", "50%", "75%", "100%"),
      title_gp = gpar_f(fontsize = 10, fontface = "bold"),
      labels_gp = gpar_f(fontsize = 9),
      border = "black", legend_direction = "vertical"
    )
  )

  png(
    filename = file.path(OUTPUT_DIR, "Fig1_kegg_decoder_ES.png"),
    width = 4800, height = 3000, res = 300
  )
  # padding is c(bottom, left, top, right) in mm. The right margin is set
  # explicitly so the row names do not run into the edge of the PNG.
  draw(ht_dec, padding = unit(c(2, 2, 2, 10), "mm"))
  {
    lgd_main <- Legend(
      col_fun   = col_fun_dec,
      title     = "Completitud\nde rutas",
      at        = c(0, 0.25, 0.5, 0.75, 1),
      labels    = c("0%", "25%", "50%", "75%", "100%"),
      title_gp  = gpar_f(fontsize = 10, fontface = "bold"),
      labels_gp = gpar_f(fontsize = 9), border = "black"
    )
    valid_dep_dec <- na.omit(unique(as.character(meta_decoder$Depth)))
    lgd_dep <- Legend(
      labels    = valid_dep_dec,
      legend_gp = gpar_f(fill = depth_col_dec[valid_dep_dec]),
      title     = "Profundidad",
      title_gp  = gpar_f(fontsize = 10, fontface = "bold"),
      labels_gp = gpar_f(fontsize = 9)
    )
    .draw_br_legend(list(lgd_main, lgd_dep), y_mm = 20)
  }
  dev.off()
  cat("[OK] Guardado: Fig1_kegg_decoder_ES.png\n\n")

  # ---- cladogram variant: no depth sidebar, rows cluster freely ----
  ht_dec_cluster <- Heatmap(
    decoder_mat,
    name = "Completitud",
    col = col_fun_dec,
    cluster_rows = TRUE, show_row_dend = TRUE,
    row_dend_side = "left",
    cluster_columns = TRUE,
    clustering_distance_columns = "euclidean",
    clustering_method_columns = "ward.D2",
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "ward.D2",
    show_row_names = TRUE, show_column_names = TRUE,
    show_heatmap_legend = FALSE,
    row_labels = make_taxon_label_expr(rownames(decoder_mat)),
    row_names_gp = gpar_f(fontsize = 10),
    row_names_max_width = unit(11, "cm"),
    column_names_gp = gpar_f(fontsize = 8),
    column_names_rot = 45,
    rect_gp = gpar_f(col = "black", lwd = 0.5),
    border = TRUE,
    heatmap_legend_param = list(
      title = "Completitud\nde rutas",
      at = c(0, 0.25, 0.5, 0.75, 1.0),
      labels = c("0%", "25%", "50%", "75%", "100%"),
      title_gp = gpar_f(fontsize = 10, fontface = "bold"),
      labels_gp = gpar_f(fontsize = 9),
      border = "black", legend_direction = "vertical"
    )
  )

  png(
    filename = .cluster_out(file.path(OUTPUT_DIR, "Fig1_kegg_decoder_ES.png")),
    width = 4800, height = 3000, res = 300
  )
  draw(ht_dec_cluster, padding = unit(c(2, 2, 2, 10), "mm"))
  {
    lgd_main_c <- Legend(
      col_fun   = col_fun_dec,
      title     = "Completitud\nde rutas",
      at        = c(0, 0.25, 0.5, 0.75, 1),
      labels    = c("0%", "25%", "50%", "75%", "100%"),
      title_gp  = gpar_f(fontsize = 10, fontface = "bold"),
      labels_gp = gpar_f(fontsize = 9), border = "black"
    )
    .draw_br_legend(list(lgd_main_c), y_mm = 20)
  }
  dev.off()
  cat("[OK] Guardado: Fig1_kegg_decoder_cluster_ES.png\n\n")
} else {
  cat("[WARN] Archivo KEGG-Decoder no encontrado, omitiendo...\n\n")
}

# ====================== KEGG-Expander Heatmap ================================
cat("=== Generando KEGG-Expander Heatmap ===\n")
if (file.exists(EXPANDER_FILE)) {
  expander_mat <- read_kegg_list(EXPANDER_FILE)
  cat("  Dimensiones originales:", dim(expander_mat), "\n")
  expander_mat <- filter_kegg_zero_cols(expander_mat)
  cat("  Dimensiones filtradas:", dim(expander_mat), "\n\n")

  cache_file_expander <- file.path(OUTPUT_DIR, "translation_cache_expander.csv")
  colnames(expander_mat) <- capitalize_first(translate_to_spanish(colnames(expander_mat),
    translation_cache_file = cache_file_expander
  ))

  expander_mag_names <- rownames(expander_mat)
  meta_exp <- load_kegg_metadata(expander_mag_names)
  expander_mat <- kegg_rename_rows(expander_mat, meta_exp)
  rownames(meta_exp) <- meta_exp$DisplayName

  cat("  Ordenando MAGs por profundidad...\n")
  ord_exp <- kegg_order_by_depth(expander_mat, meta_exp)
  expander_mat <- ord_exp$mat
  meta_exp <- ord_exp$metadata

  depth_col_exp <- kegg_named_colors(meta_exp$Depth, "Set2")
  max_count <- max(expander_mat, na.rm = TRUE)
  col_fun_exp <- colorRamp2(
    seq(0, max_count, length.out = 6),
    c("white", "yellow", "orange", "orangered", "red", "darkred")
  )

  ra_exp <- rowAnnotation(
    Profundidad = meta_exp$Depth,
    col = list(Profundidad = depth_col_exp),
    annotation_name_side = "top",
    show_legend = FALSE,
    simple_anno_size = unit(0.4, "cm"),
    annotation_name_gp = gpar_f(fontsize = 10, fontface = "bold"),
    annotation_legend_param = list(
      Profundidad = list(
        at = na.omit(unique(as.character(meta_exp$Depth)))
      )
    )
  )

  # Drop columns with an NA name
  expander_mat <- expander_mat[, !is.na(colnames(expander_mat)), drop = FALSE]

  ht_exp <- Heatmap(
    expander_mat,
    name = "Abundancia",
    col = col_fun_exp,
    cluster_rows = FALSE, cluster_columns = TRUE,
    show_row_names = TRUE, show_column_names = TRUE,
    show_heatmap_legend = FALSE,
    row_labels = make_taxon_label_expr(rownames(expander_mat)),
    row_names_gp = gpar_f(fontsize = 11),
    column_names_gp = gpar_f(fontsize = 6),
    column_names_rot = 45,
    left_annotation = ra_exp,
    rect_gp = gpar_f(col = "black", lwd = 0.5),
    border = TRUE,
    column_title = "Abundancia de funciones metabólicas por MAG\nKEGG-Expander",
    column_title_gp = gpar_f(fontsize = 11, fontface = "bold"),
    heatmap_legend_param = list(
      title = "Abundancia de\nfunción",
      at = round(seq(0, max_count, length.out = 5)),
      legend_direction = "vertical",
      title_gp = gpar_f(fontsize = 10, fontface = "bold"),
      labels_gp = gpar_f(fontsize = 9),
      border = "black"
    )
  )

  png(
    filename = file.path(OUTPUT_DIR, "Fig2_kegg_expander_ES.png"),
    width = 4200, height = 3000, res = 300
  )
  draw(ht_exp)
  {
    lgd_main <- Legend(
      col_fun   = col_fun_exp,
      title     = "Abundancia de\nfunción",
      at        = round(seq(0, max_count, length.out = 5)),
      title_gp  = gpar_f(fontsize = 10, fontface = "bold"),
      labels_gp = gpar_f(fontsize = 9), border = "black"
    )
    valid_dep_exp <- na.omit(unique(as.character(meta_exp$Depth)))
    lgd_dep <- Legend(
      labels    = valid_dep_exp,
      legend_gp = gpar_f(fill = depth_col_exp[valid_dep_exp]),
      title     = "Profundidad",
      title_gp  = gpar_f(fontsize = 10, fontface = "bold"),
      labels_gp = gpar_f(fontsize = 9)
    )
    .draw_br_legend(list(lgd_main, lgd_dep))
  }
  dev.off()
  cat("[OK] Guardado: Fig2_kegg_expander_ES.png\n\n")
} else {
  cat("[WARN] Archivo KEGG-Expander no encontrado, omitiendo...\n\n")
}

# ====================== Combined Decoder + Expander ==========================
cat("=== Generando Heatmap Combinado KEGG ===\n")
if (file.exists(DECODER_FILE) && file.exists(EXPANDER_FILE)) {
  dec_orig <- rownames(read_kegg_list(DECODER_FILE))
  exp_orig <- rownames(read_kegg_list(EXPANDER_FILE))
  common_mags <- intersect(dec_orig, exp_orig)
  cat("  MAGs en Decoder:", length(dec_orig), "\n")
  cat("  MAGs en Expander:", length(exp_orig), "\n")
  cat("  MAGs comunes:", length(common_mags), "\n\n")

  if (length(common_mags) >= 3) {
    dec_m <- filter_kegg_zero_cols(read_kegg_list(DECODER_FILE))
    exp_m <- filter_kegg_zero_cols(read_kegg_list(EXPANDER_FILE))

    dec_norm <- dec_m[common_mags, , drop = FALSE]
    exp_norm <- exp_m[common_mags, , drop = FALSE] /
      max(exp_m[common_mags, ], na.rm = TRUE)

    colnames(dec_norm) <- translate_to_spanish(colnames(dec_norm),
      translation_cache_file = file.path(OUTPUT_DIR, "translation_cache_decoder.csv")
    )
    colnames(exp_norm) <- translate_to_spanish(colnames(exp_norm),
      translation_cache_file = file.path(OUTPUT_DIR, "translation_cache_expander.csv")
    )

    combined_mat <- cbind(dec_norm, exp_norm)
    meta_comb <- load_kegg_metadata(common_mags)
    dn_map <- setNames(meta_comb$DisplayName, common_mags)
    rownames(combined_mat) <- dn_map[rownames(combined_mat)]
    rownames(meta_comb) <- meta_comb$DisplayName

    cat("  Ordenando MAGs por profundidad...\n")
    ord_comb <- kegg_order_by_depth(combined_mat, meta_comb)
    combined_mat <- ord_comb$mat
    meta_comb <- ord_comb$metadata

    depth_col_comb <- kegg_named_colors(meta_comb$Depth, "Set2")

    ra_comb <- rowAnnotation(
      Profundidad = meta_comb$Depth,
      col = list(Profundidad = depth_col_comb),
      annotation_name_side = "top",
      simple_anno_size = unit(0.4, "cm"),
      annotation_name_gp = gpar_f(fontsize = 10, fontface = "bold"),
      annotation_legend_param = list(
        Profundidad = list(
          at = na.omit(unique(as.character(meta_comb$Depth)))
        )
      )
    )

    # Drop columns with an NA name
    combined_mat <- combined_mat[, !is.na(colnames(combined_mat)), drop = FALSE]

    col_fun_comb <- colorRamp2(
      c(0, 0.2, 0.4, 0.6, 0.8, 1),
      c("white", "yellow", "orange", "orangered", "red", "darkred")
    )

    ht_comb <- Heatmap(
      combined_mat,
      name = "Completitud\nde rutas",
      col = col_fun_comb,
      cluster_rows = FALSE, cluster_columns = TRUE,
      show_row_names = TRUE, show_column_names = TRUE,
      row_labels = make_taxon_label_expr(rownames(combined_mat)),
      row_names_gp = gpar_f(fontsize = 11),
      column_names_gp = gpar_f(fontsize = 5),
      column_names_rot = 45,
      left_annotation = ra_comb,
      rect_gp = gpar_f(col = "black", lwd = 0.5),
      border = TRUE,
      column_title = "Completitud combinada de rutas metabólicas por MAG\nKEGG-Decoder + KEGG-Expander",
      column_title_gp = gpar_f(fontsize = 11, fontface = "bold"),
      heatmap_legend_param = list(
        title = "Completitud\nde rutas",
        at = c(0, 0.25, 0.5, 0.75, 1.0),
        labels = c("0%", "25%", "50%", "75%", "100%"),
        title_gp = gpar_f(fontsize = 10, fontface = "bold"),
        labels_gp = gpar_f(fontsize = 9),
        border = "black"
      )
    )

    png(
      filename = file.path(OUTPUT_DIR, "Fig3_kegg_combined_ES.png"),
      width = 5000, height = 3600, res = 300
    )
    draw(ht_comb)
    dev.off()
    cat("[OK] Guardado: Fig3_kegg_combined_ES.png\n\n")
  } else {
    cat("[WARN] Pocos MAGs comunes (", length(common_mags), "), omitiendo heatmap combinado\n\n")
  }
} else {
  cat("[WARN] Archivos KEGG no encontrados, omitiendo heatmap combinado\n\n")
}

cat("[OK] Sección KEGG-Decoder/Expander completada\n")
cat("  Caché de traducciones en:", OUTPUT_DIR, "\n\n")

# ---- RUTAS Y DATOS DE DRAM2 DISTILLATION (compartidos en secciones 11-15) ----
PRODUCT_FILE <- file.path(PROJECT_DIR, "results", "08_MAG_distill", "VISUALIZE", "product.tsv")
GENOME_STATS_FILE <- file.path(PROJECT_DIR, "results", "08_MAG_distill", "SUMMARIZE", "genome_stats.tsv")
QUALITY_REPORT_FILE <- file.path(PROJECT_DIR, "results", "02_mags_quality", "20251222_valentin_quality_report.tsv")
SUMMARIZED_GENOMES_FILE <- file.path(PROJECT_DIR, "results", "08_MAG_distill", "SUMMARIZE", "summarized_genomes.tsv")
METABOLISM_FILE <- file.path(PROJECT_DIR, "results", "08_MAG_distill", "SUMMARIZE", "metabolism_summary.xlsx")
TRAITS_FILE <- file.path(PROJECT_DIR, "results", "08_MAG_distill", "SUMMARIZE", "traits.xlsx")

# One translation cache shared across the sections below
DRAM_CACHE_FILE <- file.path(OUTPUT_DIR, "translation_cache_dram2_distill.csv")

# Builds the depth annotation, plus station where asked, for any vector of
# genome IDs (NameWithoutExtension, no .fa suffix).
build_row_annotation <- function(genome_ids, include_station = FALSE) {
  lv <- if (exists("depth_levels_std")) depth_levels_std else c("ZMOR", "1000 m", "Fondo")
  depths <- factor(get_depth(genome_ids), levels = lv)
  df <- data.frame(Profundidad = depths, row.names = genome_ids, stringsAsFactors = FALSE)
  if (include_station) {
    df$Estacion <- as.character(get_station(genome_ids))
  }
  df
}

# rowAnnotation covering the five annotation variants
# annotation_version: NULL (compat. anterior) | 1=prof+est+pct | 2=prof+est |
#                     3=prof | 4=est | 5=pct
# completeness: numeric 0 to 100, aligned with depths, for the percentage bar
make_depth_row_annotation <- function(depths, depth_colors, stations = NULL,
                                      annotation_version = NULL,
                                      completeness = NULL,
                                      show_legend = TRUE) {
  valid_levels <- as.character(na.omit(unique(depths)))

  # Helper interno: construir paleta de estaciones
  build_st_cols <- function(sts) {
    unique_st <- sort(unique(na.omit(as.character(sts))))
    n_st <- max(length(unique_st), 3)
    setNames(
      RColorBrewer::brewer.pal(n_st, "Set1")[seq_len(length(unique_st))],
      unique_st
    )
  }

  # Colour ramp for completeness: white through bright blue to dark blue
  col_pct <- colorRamp2(c(0, 50, 100), c("white", "#9ECAE1", "#08519C"))

  # Legend parameters for the completeness annotation
  pct_leg_params <- list(
    title      = "% Completitud\n(anotación HQ)",
    at         = c(0, 25, 50, 75, 100),
    labels     = c("0%", "25%", "50%", "75%", "100%"),
    title_gp   = gpar_f(fontsize = 10, fontface = "bold"),
    labels_gp  = gpar_f(fontsize = 9),
    border     = "black"
  )

  # Work out which variant applies
  v <- if (is.null(annotation_version)) {
    if (!is.null(stations)) 2L else 3L # comportamiento original
  } else {
    as.integer(annotation_version)
  }

  # v == 0 means no sidebar; the rows are clustered by the caller
  if (v == 0L) {
    return(NULL)
  }

  ann_gp <- gpar_f(fontsize = 10, fontface = "bold")
  ann_sz <- unit(0.4, "cm")

  if (v == 1L) {
    # depth, station and completeness
    st_cols <- build_st_cols(stations)
    rowAnnotation(
      Profundidad = depths,
      `Estación` = factor(stations, levels = names(st_cols)),
      Completitud = as.numeric(completeness),
      col = c(depth_colors, list(`Estación` = st_cols, Completitud = col_pct)),
      annotation_name_side = "top", show_legend = show_legend,
      simple_anno_size = ann_sz,
      annotation_name_gp = ann_gp,
      annotation_legend_param = list(
        Profundidad = list(at = valid_levels, labels = valid_levels),
        Completitud = pct_leg_params
      )
    )
  } else if (v == 2L) {
    # depth and station
    st_cols <- build_st_cols(stations)
    rowAnnotation(
      Profundidad = depths,
      `Estación` = factor(stations, levels = names(st_cols)),
      col = c(depth_colors, list(`Estación` = st_cols)),
      annotation_name_side = "top", show_legend = show_legend,
      simple_anno_size = ann_sz,
      annotation_name_gp = ann_gp,
      annotation_legend_param = list(
        Profundidad = list(at = valid_levels, labels = valid_levels)
      )
    )
  } else if (v == 3L) {
    # solo profundidad
    rowAnnotation(
      Profundidad = depths,
      col = depth_colors,
      annotation_name_side = "top", show_legend = show_legend,
      simple_anno_size = ann_sz,
      annotation_name_gp = ann_gp,
      annotation_legend_param = list(
        Profundidad = list(at = valid_levels, labels = valid_levels)
      )
    )
  } else if (v == 4L) {
    # station only
    st_cols <- build_st_cols(stations)
    rowAnnotation(
      `Estación` = factor(stations, levels = names(st_cols)),
      col = list(`Estación` = st_cols),
      annotation_name_side = "top", show_legend = show_legend,
      simple_anno_size = ann_sz,
      annotation_name_gp = ann_gp
    )
  } else {
    # v == 5L: solo % completitud
    rowAnnotation(
      Completitud = as.numeric(completeness),
      col = list(Completitud = col_pct),
      annotation_name_side = "top", show_legend = show_legend,
      simple_anno_size = ann_sz,
      annotation_name_gp = ann_gp,
      annotation_legend_param = list(
        Completitud = pct_leg_params
      )
    )
  }
}

# Paleta de profundidad
depth_levels_std <- c("ZMOR", "1000 m", "Fondo")
depth_colors_std <- setNames(
  RColorBrewer::brewer.pal(3, "Set2"),
  depth_levels_std
)

# ---- 11. HEATMAP DE PRODUCT.TSV – Completitud de rutas DRAM2 ----
cat("\n=== SECCIÓN 11: HEATMAP DE PRODUCT.TSV ===\n")

if (file.exists(PRODUCT_FILE) && file.exists(GENOME_STATS_FILE)) {
  genome_stats <- read_tsv(GENOME_STATS_FILE, show_col_types = FALSE)
  mag_order_gs <- genome_stats$genome # NameWithoutExtension, sin .fa

  # ---- completeness and contamination sidebars, as in product.html ----
  completeness_map_prod <- NULL
  contamination_map_prod <- NULL
  if (file.exists(QUALITY_REPORT_FILE)) {
    orig_name_col <- names(metadata)[11]  # OriginalName
    meta_orig_prod <- metadata %>%
      select(NameWithoutExtension, OriginalName = all_of(orig_name_col))
    qual_prod_df <- tryCatch(
      read_tsv(QUALITY_REPORT_FILE, show_col_types = FALSE) %>%
        inner_join(meta_orig_prod, by = c("Name" = "OriginalName")) %>%
        filter(NameWithoutExtension %in% mag_order_gs),
      error = function(e) NULL
    )
    if (!is.null(qual_prod_df) && nrow(qual_prod_df) > 0) {
      dn_keys_prod <- convert_to_display_name(qual_prod_df$NameWithoutExtension)
      completeness_map_prod <- setNames(as.numeric(qual_prod_df$Completeness), dn_keys_prod)
      contamination_map_prod <- setNames(as.numeric(qual_prod_df$Contamination), dn_keys_prod)
      cat("  [OK] Datos de calidad (product):", nrow(qual_prod_df), "MAGs\n")
    }
  }

  product_raw <- read_tsv(PRODUCT_FILE, show_col_types = FALSE)
  cat("  Dimensiones product.tsv:", nrow(product_raw), "x", ncol(product_raw), "\n")

  if (nrow(product_raw) == length(mag_order_gs)) {
    product_raw$genome <- mag_order_gs
  } else {
    cat("  [WARN] Desajuste en número de filas; asignando índices\n")
    product_raw$genome <- paste0("MAG_", seq_len(nrow(product_raw)))
  }

  # ---- assign a functional group to each product.tsv column ----
  assign_product_group <- function(col_name, is_bool) {
    if (is_bool) {
      trimws(sub(":.*$", "", col_name))
    } else if (grepl("^Complex", col_name)) {
      sub("^(Complex [IVX]+(?:\\s+[HL])?):.*$", "\\1", col_name, perl = TRUE)
    } else {
      "Carbon fixation"
    }
  }

  # --- Helper interno: construir y guardar heatmap de product.tsv ---
  # Handles both the boolean columns and the 0-to-1 numeric ones
  make_product_heatmap <- function(prod_raw, mag_order_ids,
                                   title_main, outfile,
                                   annotation_version = NULL,
                                   completeness_map = NULL,
                                   contamination_map = NULL,
                                   cache_f = DRAM_CACHE_FILE,
                                   legend_side = "right",
                                   padding = unit(c(20, 32, 5, 2), "mm"),
                                   barplot_gap = 2,
                                   extra_canvas_in = 0) {
    all_cols <- setdiff(names(prod_raw), "genome")

    # Separate boolean columns from numeric
    is_bool_col <- vapply(all_cols, function(cn) {
      vals <- prod_raw[[cn]]
      all(is.na(vals) | as.character(vals) %in% c("True", "False", "TRUE", "FALSE"))
    }, logical(1))
    bool_cols <- all_cols[is_bool_col]
    num_cols <- all_cols[!is_bool_col]

    # Columns carrying at least one value
    bool_keep <- bool_cols[vapply(bool_cols, function(cn) {
      any(as.character(prod_raw[[cn]]) %in% c("True", "TRUE"), na.rm = TRUE)
    }, logical(1))]
    num_keep <- num_cols[vapply(num_cols, function(cn) {
      vals <- suppressWarnings(as.numeric(as.character(prod_raw[[cn]])))
      any(!is.na(vals) & vals > 0)
    }, logical(1))]

    if (length(bool_keep) == 0 && length(num_keep) == 0) {
      cat("  [WARN] Sin columnas válidas en product.tsv\n")
      return(invisible(NULL))
    }

    # Matrices
    if (length(num_keep) > 0) {
      pm_num <- prod_raw %>%
        select(genome, all_of(num_keep)) %>%
        column_to_rownames("genome") %>%
        mutate(across(everything(), ~ suppressWarnings(as.numeric(as.character(.))))) %>%
        as.matrix()
    }
    if (length(bool_keep) > 0) {
      pm_bool <- prod_raw %>%
        select(genome, all_of(bool_keep)) %>%
        column_to_rownames("genome") %>%
        mutate(across(everything(), ~ as.numeric(as.character(.) %in% c("True", "TRUE")))) %>%
        as.matrix()
    }

    # Groups for column_split, translated through .MANUAL_TRANSLATIONS
    translate_group <- function(g) {
      ifelse(g %in% names(.MANUAL_TRANSLATIONS),
        as.character(.MANUAL_TRANSLATIONS[g]), g
      )
    }
    if (length(num_keep) > 0) {
      num_groups <- translate_group(
        vapply(num_keep, assign_product_group, character(1), is_bool = FALSE)
      )
      num_split <- factor(num_groups, levels = unique(num_groups))
    }
    if (length(bool_keep) > 0) {
      bool_groups <- translate_group(
        vapply(bool_keep, assign_product_group, character(1), is_bool = TRUE)
      )
      bool_split <- factor(bool_groups, levels = unique(bool_groups))
    }

    # Translations, capitalised
    # Strip the "Grupo: " prefix first, since column_split already shows the group
    strip_col_prefix <- function(x) trimws(sub("^[^:]+:\\s*", "", x))
    if (length(num_keep) > 0) {
      num_labels <- strip_col_prefix(num_keep)
      num_es <- capitalize_first(translate_to_spanish(num_labels, translation_cache_file = cache_f))
      num_es[is.na(num_es)] <- num_labels[is.na(num_es)]
      colnames(pm_num) <- num_es
      valid_num <- !is.na(colnames(pm_num))
      pm_num <- pm_num[, valid_num, drop = FALSE]
      num_split <- num_split[valid_num]
    }
    if (length(bool_keep) > 0) {
      bool_labels <- strip_col_prefix(bool_keep)
      bool_es <- capitalize_first(translate_to_spanish(bool_labels, translation_cache_file = cache_f))
      bool_es[is.na(bool_es)] <- bool_labels[is.na(bool_es)]
      colnames(pm_bool) <- bool_es
      valid_bool <- !is.na(colnames(pm_bool))
      pm_bool <- pm_bool[, valid_bool, drop = FALSE]
      bool_split <- bool_split[valid_bool]
    }

    # Display names, ordered by depth
    dn_p <- convert_to_display_name(mag_order_ids)
    dn_to_id <- setNames(mag_order_ids, dn_p) # display_name to genome_id
    anno_p <- build_row_annotation(mag_order_ids)
    rownames(anno_p) <- dn_p
    anno_p$Profundidad <- factor(anno_p$Profundidad, levels = depth_levels_std)

    ref_mat <- if (length(num_keep) > 0) pm_num else pm_bool
    rownames(ref_mat) <- dn_p
    ord_p <- order_matrix_by_depth(ref_mat, anno_p)
    row_ord <- rownames(ord_p$matrix)
    an_p <- ord_p$annotation

    if (length(num_keep) > 0) {
      rownames(pm_num) <- dn_p
      pm_num <- pm_num[row_ord, , drop = FALSE]
    }
    if (length(bool_keep) > 0) {
      rownames(pm_bool) <- dn_p
      pm_bool <- pm_bool[row_ord, , drop = FALSE]
    }

    # Re-sort for variant 4 (station) or 5 (completeness)
    v <- if (is.null(annotation_version)) 3L else as.integer(annotation_version)
    cluster_rows_p <- (v == 0L)
    if (v %in% c(4L, 5L)) {
      ref_v <- if (length(num_keep) > 0) pm_num else pm_bool
      ro <- reorder_rows_by_version(ref_v, an_p, v,
        dn_to_id = dn_to_id,
        completeness_map = completeness_map
      )
      row_ord <- rownames(ro$mat)
      an_p <- ro$ann
      if (length(num_keep) > 0) pm_num <- pm_num[row_ord, , drop = FALSE]
      if (length(bool_keep) > 0) pm_bool <- pm_bool[row_ord, , drop = FALSE]
    }

    # Sidebar
    dc_p <- depth_colors_std[levels(an_p$Profundidad)]
    dc_p <- dc_p[!is.na(names(dc_p))]
    stations_p <- if (v %in% c(1L, 2L, 4L)) {
      as.character(get_station(dn_to_id[row_ord]))
    } else {
      NULL
    }
    compl_p <- if (!is.null(completeness_map)) completeness_map[row_ord] else NULL
    contam_p <- if (!is.null(contamination_map)) contamination_map[row_ord] else NULL

    ra_p <- make_depth_row_annotation(
      an_p$Profundidad, list(Profundidad = dc_p),
      stations = stations_p, annotation_version = v, completeness = compl_p
    )

    # The quality bars go into the last heatmap's right_annotation, and the MAG
    # names ride along as show_row_names = TRUE. ComplexHeatmap always lays row
    # names out before the legend, which is what stops the two from overlapping.
    right_parts      <- list()
    right_widths     <- c()
    right_show_name  <- c()

    if (!is.null(compl_p)) {
      cv_c      <- as.numeric(compl_p)
      col_compl <- colorRamp2(c(0, 50, 100), c("white", "#6BAED6", "#08519C"))
      right_parts[["Completitud\n(%)"]] <- anno_barplot(
        cv_c, baseline = 0, axis = TRUE, which = "row",
        gp = gpar_f(fill = col_compl(cv_c), col = "grey50", lwd = 0.3),
        bar_width = 0.7,
        axis_param = list(gp = gpar_f(fontsize = 6), at = c(0, 50, 100))
      )
      right_widths    <- c(right_widths, 2.5)
      right_show_name <- c(right_show_name, TRUE)
    }
    if (!is.null(contam_p)) {
      cv_k  <- as.numeric(contam_p)
      max_k <- max(cv_k, na.rm = TRUE)
      if (is.na(max_k) || max_k == 0) max_k <- 10
      col_contam <- colorRamp2(c(0, max_k / 2, max_k), c("white", "#FC8D59", "#D73027"))
      right_parts[["Contaminación\n(%)"]] <- anno_barplot(
        cv_k, baseline = 0, axis = TRUE, which = "row",
        gp = gpar_f(fill = col_contam(cv_k), col = "grey50", lwd = 0.3),
        bar_width = 0.7,
        axis_param = list(gp = gpar_f(fontsize = 6))
      )
      right_widths    <- c(right_widths, 2.5)
      right_show_name <- c(right_show_name, TRUE)
    }
    n_qual_bars <- length(right_parts)

    # MAG names join the same annotation block
    right_parts[["   "]] <- anno_text(make_taxon_label_expr(row_ord), which = "row",
      gp = gpar_f(fontsize = 10))
    right_widths    <- c(right_widths, 6)
    right_show_name <- c(right_show_name, FALSE)

    # width must be set explicitly. ComplexHeatmap positions the legend from
    # @width, and left to recalculate it internally it can get it wrong and drop
    # the legend on top of the row names. Setting width makes width(anno) return
    # the real value, so the legend lands to the right of everything.
    total_right_cm <- sum(right_widths)
    all_right_anno <- do.call(rowAnnotation, c(
      right_parts,
      list(
        annotation_name_gp   = gpar_f(fontsize = 8, fontface = "bold"),
        annotation_width     = unit(right_widths, "cm"),
        gap                  = unit(barplot_gap, "mm"),
        show_annotation_name = right_show_name,
        width                = unit(total_right_cm, "cm")   # ← clave
      )
    ))

    # The individual heatmaps carry neither right_annotation nor row names.
    # all_right_anno is added to the HeatmapList with +, which puts the legend
    # after it.
    has_bool <- length(bool_keep) > 0
    ht_list <- NULL
    cf_cont <- colorRamp2(
      seq(0, 1, length.out = 6),
      c("white", "yellow", "orange", "orangered", "red", "darkred")
    )

    # In the clustered variant the legends are drawn by hand, bottom right
    manual_lgd <- cluster_rows_p && legend_side == "right"

    if (length(num_keep) > 0) {
      ht_list <- Heatmap(
        pm_num,
        name = "Completitud\nde ruta",
        col = cf_cont,
        cluster_rows = cluster_rows_p, cluster_columns = FALSE,
        show_row_dend = cluster_rows_p, row_dend_side = "left",
        show_row_names = FALSE,
        show_heatmap_legend = !manual_lgd,
        column_names_gp = gpar_f(fontsize = 9), column_names_rot = 45,
        left_annotation  = ra_p,
        rect_gp = gpar_f(col = "black", lwd = 0.4), border = TRUE,
        column_split = num_split,
        column_gap = unit(2, "mm"),
        column_title_gp = gpar_f(fontsize = 9, fontface = "bold"),
        heatmap_legend_param = list(
          title = "Completitud\nde ruta",
          at = c(0, 0.25, 0.5, 0.75, 1),
          labels = c("0%", "25%", "50%", "75%", "100%"),
          title_gp = gpar_f(fontsize = 10, fontface = "bold"),
          labels_gp = gpar_f(fontsize = 9),
          border = "black", legend_direction = "vertical"
        )
      )
    }
    if (has_bool) {
      ht_bool <- Heatmap(
        pm_bool,
        name = "Rasgo",
        col = c("0" = "white", "1" = "#2D6A4F"),
        cluster_rows = cluster_rows_p, cluster_columns = FALSE,
        show_row_dend = FALSE,
        show_row_names = FALSE,
        show_heatmap_legend = !manual_lgd,
        column_names_gp = gpar_f(fontsize = 9), column_names_rot = 45,
        rect_gp = gpar_f(col = "black", lwd = 0.4), border = TRUE,
        column_split = bool_split,
        column_gap = unit(2, "mm"),
        column_title_gp = gpar_f(fontsize = 9, fontface = "bold"),
        heatmap_legend_param = list(
          title = "Rasgo",
          at = c(0, 1), labels = c("Ausente", "Presente"),
          title_gp = gpar_f(fontsize = 10, fontface = "bold"),
          labels_gp = gpar_f(fontsize = 9),
          border = "black", legend_direction = "vertical"
        )
      )
      ht_list <- if (is.null(ht_list)) ht_bool else ht_list + ht_bool
    }
    # The right annotation joins the HeatmapList with an explicit width, so
    # ComplexHeatmap places the legend after it
    ht_list <- ht_list + all_right_anno

    n_rows <- if (length(num_keep) > 0) nrow(pm_num) else nrow(pm_bool)
    total_cols <- (if (length(num_keep) > 0) ncol(pm_num) else 0L) +
      (if (length(bool_keep) > 0) ncol(pm_bool) else 0L)
    left_ann_w  <- if (is.null(ra_p)) 0.0 else 4.0
    # right_widths already accounts for the 6 cm of row names
    right_ann_w <- total_right_cm * 0.3937 + 0.4
    dend_w      <- if (cluster_rows_p) 1.8 else 0
    # width=NULL en ambas versiones: ComplexHeatmap auto-size y reserva espacio de leyenda
    w_raw <- round(total_cols * 0.46 + left_ann_w + right_ann_w + dend_w + 0.5)
    w_p   <- max(24L, w_raw) + extra_canvas_in
    lgd_h_extra <- if (manual_lgd) 1.8 else 0
    h_floor_p   <- min(10L, max(8L, as.integer(n_rows * 0.40)))
    h_p <- max(h_floor_p, round(n_rows * 0.20 + 3)) + lgd_h_extra
    cat(sprintf(
      "  [DIAG %s] total_cols=%d  left=%.1f  right=%.2f  dend=%.1f  w_p=%d  cluster=%s\n",
      basename(outfile), total_cols, left_ann_w, right_ann_w, dend_w, w_p, cluster_rows_p
    ))
    png(outfile, width = w_p, height = h_p, units = "in", res = 300)
    draw(ht_list,
      column_title           = title_main,
      column_title_gp        = gpar_f(fontsize = 11, fontface = "bold"),
      merge_legend           = TRUE,
      heatmap_legend_side    = legend_side,
      annotation_legend_side = legend_side,
      padding                = padding
    )
    if (manual_lgd) {
      lgds <- list()
      if (length(num_keep) > 0)
        lgds[["compl"]] <- Legend(
          col_fun = cf_cont, title = "Completitud\nde ruta",
          at = c(0, 0.25, 0.5, 0.75, 1),
          labels = c("0%", "25%", "50%", "75%", "100%"),
          title_gp = gpar_f(fontsize = 10, fontface = "bold"),
          labels_gp = gpar_f(fontsize = 9), border = "black"
        )
      if (length(bool_keep) > 0)
        lgds[["pres"]] <- Legend(
          labels = c("Ausente", "Presente"),
          legend_gp = gpar_f(fill = c("white", "#2D6A4F")),
          title = "Rasgo",
          title_gp = gpar_f(fontsize = 10, fontface = "bold"),
          labels_gp = gpar_f(fontsize = 9), border = "black"
        )
      if (length(lgds) > 0) {
        pd <- do.call(packLegend, c(lgds, list(direction = "horizontal", gap = unit(6, "mm"))))
        draw(pd,
          x    = unit(1, "npc") - unit(18, "mm"),
          y    = unit(50, "mm"),
          just = c("right", "bottom")
        )
      }
    }
    dev.off()
    cat("[OK] Guardado:", basename(outfile), "\n")
  }

  # Cladogram filename: insert _cluster before _ES.png
  .cluster_out <- function(f) sub("(_ES\\.png)$", "_cluster\\1", f)

  # Every column, boolean and numeric
  all_prod_cols <- setdiff(names(product_raw), "genome")

  # Columns with at least one positive
  non_empty_prod <- all_prod_cols[vapply(all_prod_cols, function(cn) {
    v <- as.character(product_raw[[cn]])
    any(v %in% c("True", "TRUE") | suppressWarnings(as.numeric(v)) > 0, na.rm = TRUE)
  }, logical(1))]
  cat("  Columnas no vacías:", length(non_empty_prod), "de", length(all_prod_cols), "\n")

  # Fig4: product filtered to pathways and genes with at least one hit
  make_product_heatmap(
    product_raw %>% select(genome, all_of(non_empty_prod)),
    mag_order_gs,
    title_main = "Completitud de rutas metabólicas por MAG\nDRAM2",
    outfile = file.path(OUTPUT_DIR, "Fig4_product_completeness_ES.png"),
    completeness_map = completeness_map_prod,
    contamination_map = contamination_map_prod,
    legend_side = "right",
    padding = unit(c(20, 32, 12, 2), "mm")   # bottom/top extra para etiquetas y "Profundidad"
  )
  # Fig4, cladogram variant: no depth sidebar, rows clustered
  make_product_heatmap(
    product_raw %>% select(genome, all_of(non_empty_prod)),
    mag_order_gs,
    title_main = "Completitud de rutas metabólicas por MAG\nDRAM2",
    outfile = .cluster_out(file.path(OUTPUT_DIR, "Fig4_product_completeness_ES.png")),
    annotation_version = 0L,
    completeness_map = completeness_map_prod,
    contamination_map = contamination_map_prod,
    legend_side = "right",
    padding = unit(c(20, 30, 5, 18), "mm")
  )
  # ---- Fig4 aggregated by order ----
  cat("  Generando heatmaps de product.tsv por Orden taxonómico...\n")
  tryCatch({
    prod_sel_ord <- product_raw %>% select(genome, all_of(non_empty_prod))
    prod_agg_ord <- aggregate_to_order(prod_sel_ord)

    if (nrow(prod_agg_ord) < 2) {
      cat("  [WARN] Insuficientes Órdenes para el heatmap de product.tsv\n")
    } else {
      # Classify the original columns for column_split
      feat_cols_p   <- non_empty_prod
      is_bool_p_ord <- vapply(feat_cols_p, function(cn) {
        v <- as.character(prod_sel_ord[[cn]])
        all(is.na(v) | v %in% c("True", "False", "TRUE", "FALSE"))
      }, logical(1))
      col_split_ord <- setNames(
        vapply(feat_cols_p, function(cn) assign_product_group(cn, is_bool_p_ord[cn]), character(1)),
        feat_cols_p
      )

      col_fun_ord  <- colorRamp2(c(0, 0.5, 1), c("white", "#FFC107", "#D32F2F"))
      n_per_order  <- table(mag_to_order[mag_order_gs])

      # Traducir nombres de columnas (quitar prefijo de grupo + translate_to_spanish)
      strip_col_prefix_p <- function(x) trimws(sub("^[^:]+:\\s*", "", x))
      feat_labels_p <- strip_col_prefix_p(feat_cols_p)
      feat_es_p     <- capitalize_first(translate_to_spanish(feat_labels_p,
                         translation_cache_file = DRAM_CACHE_FILE))
      feat_es_p[is.na(feat_es_p)] <- feat_labels_p[is.na(feat_es_p)]
      names(feat_es_p) <- feat_cols_p  # orig to traducido

      # Traducir etiquetas de grupo (column_split)
      unique_grps_p <- unique(col_split_ord)
      grps_es_p     <- capitalize_first(translate_to_spanish(unique_grps_p,
                         translation_cache_file = DRAM_CACHE_FILE))
      grps_es_p[is.na(grps_es_p)] <- unique_grps_p[is.na(grps_es_p)]
      grp_map_p     <- setNames(grps_es_p, unique_grps_p)

      .make_product_ord_ht <- function(outfile_ord) {
        mat_ord_p <- prod_agg_ord %>% column_to_rownames("order") %>% as.matrix()
        keep_c    <- colSums(mat_ord_p, na.rm = TRUE) > 0
        mat_ord_p <- mat_ord_p[, keep_c, drop = FALSE]
        # Build cs_ord on the original names first, then rename
        cs_ord_raw        <- col_split_ord[colnames(mat_ord_p)]
        colnames(mat_ord_p) <- feat_es_p[colnames(mat_ord_p)]
        cs_ord            <- setNames(grp_map_p[cs_ord_raw], colnames(mat_ord_p))

        # Row labels: class in parentheses after a code-like order name
        row_lbls_p <- make_taxon_label_expr(sapply(rownames(mat_ord_p), fmt_order_label), level = "Order")

        n_vals <- as.integer(n_per_order[rownames(mat_ord_p)])
        ra_ord <- rowAnnotation(
          `No. MAGs` = anno_barplot(n_vals,
            gp = gpar_f(fill = "#4A90D9", col = "black"),
            bar_width = 0.7, width = unit(2.3, "cm")),  # +15% vs 2cm previo
          annotation_name_gp = gpar_f(fontsize = 9, fontface = "bold")
        )

        # Right margin scaled to the longest order name
        max_nm_p  <- max(nchar(sapply(rownames(mat_ord_p), fmt_order_label)))
        pad_r_p   <- max(25, ceiling(max_nm_p * 2.2))

        n_r <- nrow(mat_ord_p); n_c <- ncol(mat_ord_p)
        col_w_mm_p  <- 11L
        row_h_mm_p  <- 12L  # ~33% más bajo que los 18mm originales
        # h_p carries 50 mm more vertical slack than it used to, so the real
        # top and bottom padding fits instead of the canvas shrinking to the
        # heatmap body.
        h_p <- max(5L, ceiling((n_r * row_h_mm_p + 110L) / 25.4))
        # Horizontal canvas: body plus right_annotation, about 28 mm with the barplot
        # ampliado) + nombres (~pad_r_p) + leyenda + margen lateral
        # restaurado parcialmente (~55mm). Piso de seguridad bajado a 12 in:
        # the previous floor of 22 exceeded the real calculation and masked
        # any reduction in the buffer.
        w_p <- max(12L, ceiling((n_c * col_w_mm_p + pad_r_p + 55L) / 25.4))

        ht_ord <- Heatmap(mat_ord_p,
          name = "Proporción",
          col = col_fun_ord,
          cluster_rows = TRUE, show_row_dend = TRUE,
          row_dend_side = "left",
          cluster_columns = FALSE,
          row_labels = row_lbls_p,
          row_names_side = "right", row_names_gp = gpar_f(fontsize = 10),
          row_names_max_width = unit(pad_r_p - 5, "mm"),
          width  = unit(n_c * col_w_mm_p, "mm"),
          height = unit(n_r * row_h_mm_p, "mm"),
          column_split = cs_ord, column_gap = unit(2, "mm"),
          column_title_gp = gpar_f(fontsize = 9, fontface = "bold"),
          column_title_rot = 45,
          column_names_gp = gpar_f(fontsize = 9), column_names_rot = 45,
          rect_gp = gpar_f(col = "black", lwd = 0.4), border = TRUE,
          right_annotation = ra_ord,
          heatmap_legend_param = list(
            title = "Proporción",
            at = c(0, 0.25, 0.5, 0.75, 1),
            labels = c("0%", "25%", "50%", "75%", "100%"),
            title_gp = gpar_f(fontsize = 10, fontface = "bold"),
            labels_gp = gpar_f(fontsize = 9), border = "black"
          )
        )

        png(outfile_ord, width = w_p, height = h_p, units = "in", res = 300)
        # padding: c(bottom, left, top, right) en mm.
        # Generous vertical padding for headroom top and bottom; the side
        # padding is small but not zero, restoring about 20% of the left and
        # right margin removed in the previous iteration.
        draw(ht_ord,
          merge_legend     = TRUE,
          padding          = unit(c(25, 4, 17, 4), "mm")
        )
        dev.off()
        cat("  [OK] Guardado:", basename(outfile_ord), "\n")
      }

      .make_product_ord_ht(file.path(OUTPUT_DIR, "Fig4_product_por_orden_ES.png"))
    }
  }, error = function(e) {
    cat("  [WARN] Error en product heatmap por Orden:", conditionMessage(e), "\n")
  })
} else {
  cat("[WARN] product.tsv o genome_stats.tsv no encontrado, omitiendo sección 11\n\n")
}

# ---- 12. ESTADÍSTICAS DE CALIDAD DE GENOMAS (genome_stats.tsv) ----
cat("=== SECCIÓN 12: omitida (Fig5a y Fig5b eliminadas por diseño) ===\n\n")

# ---- 13. RASGOS METABÓLICOS PREDICHOS (traits.xlsx) ----
cat("=== SECCIÓN 13: RASGOS METABÓLICOS PREDICHOS (traits.xlsx) ===\n")

make_trait_heatmap <- function(trait_df, title_es, subtitle_tool = NULL,
                                outfile, annotation_version = NULL,
                                include_station = TRUE,
                                mag_subset = NULL,
                                completeness_map = NULL,
                                cache_file = DRAM_CACHE_FILE,
                                legend_br = TRUE,
                                legend_x_mm = 18, legend_y_mm = 20) {
  genome_col <- names(trait_df)[1]
  if (!is.null(mag_subset)) {
    trait_df <- trait_df %>% filter(!!sym(genome_col) %in% mag_subset)
  }
  genomes <- trait_df[[genome_col]]
  mat_raw <- trait_df %>%
    select(-all_of(genome_col)) %>%
    as.matrix()
  mat_bin <- matrix(as.numeric(mat_raw > 0),
    nrow = nrow(mat_raw), dimnames = list(genomes, colnames(mat_raw))
  )

  keep_cols <- colSums(mat_bin, na.rm = TRUE) > 0
  if (sum(keep_cols) == 0) {
    cat("  [WARN] Sin traits con presencia, omitiendo\n")
    return(invisible(NULL))
  }
  mat_bin <- mat_bin[, keep_cols, drop = FALSE]

  trait_names_es <- capitalize_first(
    translate_to_spanish(colnames(mat_bin), translation_cache_file = cache_file)
  )
  colnames(mat_bin) <- trait_names_es
  mat_bin <- mat_bin[, !is.na(colnames(mat_bin)), drop = FALSE]

  anno_tr <- build_row_annotation(genomes)
  dn_tr <- convert_to_display_name(genomes)
  dn_to_id_tr <- setNames(genomes, dn_tr) # display_name to genome_id
  rownames(mat_bin) <- dn_tr
  rownames(anno_tr) <- dn_tr
  anno_tr$Profundidad <- factor(anno_tr$Profundidad, levels = depth_levels_std)

  ord_tr <- order_matrix_by_depth(mat_bin, anno_tr)
  mat_ord <- as.matrix(ord_tr$matrix)
  an_tr <- ord_tr$annotation

  dc <- depth_colors_std[depth_levels_std]
  dc <- dc[!is.na(names(dc))]

  v <- if (is.null(annotation_version)) {
    if (include_station) 2L else 3L
  } else {
    as.integer(annotation_version)
  }
  cluster_rows_tr <- (v == 0L)

  # Re-sort for variant 4 (station) or 5 (completeness)
  if (v %in% c(4L, 5L)) {
    ro <- reorder_rows_by_version(mat_ord, an_tr, v,
      dn_to_id = dn_to_id_tr,
      completeness_map = completeness_map
    )
    mat_ord <- ro$mat
    an_tr <- ro$ann
  }

  # Stations come from the metadata directly, not from build_row_annotation
  stations_tr <- if (v %in% c(1L, 2L, 4L)) {
    as.character(get_station(dn_to_id_tr[rownames(mat_ord)]))
  } else {
    NULL
  }
  compl_tr <- if (!is.null(completeness_map)) completeness_map[rownames(mat_ord)] else NULL

  ra_tr <- make_depth_row_annotation(
    an_tr$Profundidad, list(Profundidad = dc),
    stations = stations_tr, annotation_version = v, completeness = compl_tr,
    show_legend = !legend_br
  )

  col_title_tr <- if (!is.null(subtitle_tool)) paste0(title_es, "\n", subtitle_tool) else title_es

  ht_tr <- Heatmap(
    mat_ord,
    name = "Presencia",
    col = c("0" = "white", "1" = "#2D6A4F"),
    cluster_rows = cluster_rows_tr, cluster_columns = TRUE,
    show_row_dend = cluster_rows_tr, row_dend_side = "left",
    clustering_distance_columns = "binary",
    clustering_method_columns = "ward.D2",
    show_row_names = TRUE, show_column_names = TRUE,
    row_labels = make_taxon_label_expr(rownames(mat_ord), level = "Order"),
    row_names_gp = gpar_f(fontsize = 10),
    column_names_gp = gpar_f(fontsize = 9), column_names_rot = 45,
    left_annotation = ra_tr,
    rect_gp = gpar_f(col = "black", lwd = 0.5), border = TRUE,
    column_title = col_title_tr,
    column_title_gp = gpar_f(fontsize = 11, fontface = "bold"),
    show_heatmap_legend = !legend_br,
    heatmap_legend_param = list(
      title = "Rasgo",
      at = c(0, 1), labels = c("Ausente", "Presente"),
      title_gp = gpar_f(fontsize = 10, fontface = "bold"),
      labels_gp = gpar_f(fontsize = 9),
      border = "black", legend_direction = "vertical"
    )
  )

  w_tr <- max(10, round(ncol(mat_ord) * 0.55 + 4))
  h_floor_tr <- min(8L, max(4L, as.integer(nrow(mat_ord) * 0.35)))
  h_tr <- max(h_floor_tr, round(nrow(mat_ord) * 0.20 + 3))
  png(outfile, width = w_tr, height = h_tr, units = "in", res = 300)
  draw(ht_tr)
  if (legend_br) {
    lgd_main <- Legend(
      labels    = c("Ausente", "Presente"),
      legend_gp = gpar_f(fill = c("white", "#2D6A4F")),
      title     = "Rasgo",
      title_gp  = gpar_f(fontsize = 10, fontface = "bold"),
      labels_gp = gpar_f(fontsize = 9), border = "black"
    )
    .draw_br_legend(c(list(lgd_main), .make_annotation_legends(dc, stations_tr)),
                    x_mm = legend_x_mm, y_mm = legend_y_mm)
  }
  dev.off()
  cat("  [OK] Guardado:", basename(outfile), "\n")
}

# Trait heatmap aggregated by order, in two variants: the normal one with a bar
# plot of MAGs per order, and the clustered one. The scale is a continuous 0 to
# 1 proportion.
make_trait_order_heatmap <- function(trait_df, title_es, outfile, subtitle_tool = NULL) {
  genome_col <- names(trait_df)[1]
  agg <- aggregate_to_order(trait_df)
  if (nrow(agg) < 2) {
    cat("  [WARN] Insuficientes Órdenes para", basename(outfile), "\n")
    return(invisible(NULL))
  }

  n_per_order <- table(mag_to_order[trait_df[[genome_col]]])
  mat <- agg %>% column_to_rownames("order") %>% as.matrix()
  keep_c <- colSums(mat, na.rm = TRUE) > 0
  mat    <- mat[, keep_c, drop = FALSE]
  if (ncol(mat) == 0) {
    cat("  [WARN] Sin columnas con datos para", basename(outfile), "\n")
    return(invisible(NULL))
  }

  cache_file_tr_ord <- file.path(OUTPUT_DIR, "translation_cache_traits.csv")
  colnames(mat) <- capitalize_first(
    translate_to_spanish(colnames(mat), translation_cache_file = cache_file_tr_ord)
  )

  col_fun_tr_ord <- colorRamp2(c(0, 0.5, 1), c("white", "#74C476", "#006D2C"))

  # Right margin scaled to the longest order name, counting the class suffix
  # that code-like names carry
  fmt_rn <- sapply(rownames(mat), fmt_order_label)
  max_name_mm  <- max(nchar(fmt_rn)) * 2.2  # ~2.2 mm por carácter a fontsize 10
  pad_right    <- max(20, ceiling(max_name_mm))     # mínimo 20 mm

  .make_ht_tr_ord <- function(out) {
    n_vals <- as.integer(n_per_order[rownames(mat)])
    ra_ord <- rowAnnotation(
      `No. MAGs` = anno_barplot(n_vals,
        gp = gpar_f(fill = "#4A90D9", col = "black"),
        bar_width = 0.7, width = unit(2, "cm")),
      annotation_name_gp = gpar_f(fontsize = 9, fontface = "bold")
    )

    n_r <- nrow(mat); n_c <- ncol(mat)
    col_w_mm_t <- 10L
    # Fig6a, energy, is tuned separately: rows about 10% thinner and narrower
    # side margins. Fig6b, carbon, and the rest keep the original values.
    is_fig6a   <- grepl("energy", out, ignore.case = TRUE)
    h_p_base   <- max(4L, round(n_r * 0.30 + 1.5))
    h_p        <- if (is_fig6a) h_p_base * 0.9 else h_p_base
    # Fig6a: canvas 6 mm narrower than before, 55 down to 49.
    w_overhead <- if (is_fig6a) 49L else 55L
    # Canvas: body, right_annotation (about 25 mm), row names (pad_right) and margins (w_overhead)
    w_p  <- max(10L, ceiling((n_c * col_w_mm_t + pad_right + w_overhead) / 25.4))

    ht <- Heatmap(mat,
      name = "Proporción de\nMAGs con rasgo",
      col = col_fun_tr_ord,
      cluster_rows = TRUE, show_row_dend = TRUE,
      row_dend_side = "left",
      cluster_columns = TRUE,
      clustering_distance_columns = "euclidean",
      clustering_method_columns = "ward.D2",
      row_labels = make_taxon_label_expr(fmt_rn, level = "Order"),
      row_names_side = "right",
      row_names_gp = gpar_f(fontsize = 10),
      row_names_max_width = unit(pad_right - 5, "mm"),
      width = unit(n_c * col_w_mm_t, "mm"),
      column_names_gp = gpar_f(fontsize = 9), column_names_rot = 45,
      rect_gp = gpar_f(col = "black", lwd = 0.5), border = TRUE,
      right_annotation = ra_ord,
      heatmap_legend_param = list(
        title = "Proporción de\nMAGs con rasgo",
        at = c(0, 0.25, 0.5, 0.75, 1),
        labels = c("0%", "25%", "50%", "75%", "100%"),
        title_gp = gpar_f(fontsize = 10, fontface = "bold"),
        labels_gp = gpar_f(fontsize = 9), border = "black"
      )
    )

    pad_lr <- if (is_fig6a) c(2L, 1L) else c(4L, 3L)  # (left, right) mm

    png(out, width = w_p, height = h_p, units = "in", res = 300)
    draw(ht, merge_legend = TRUE,
         padding = unit(c(5, pad_lr[1], 5, pad_lr[2]), "mm"))
    dev.off()
    cat("  [OK] Guardado:", basename(out), "\n")
  }

  .make_ht_tr_ord(outfile)
}

if (file.exists(TRAITS_FILE)) {
  tryCatch(
    {
      traits_energy <- read_excel(TRAITS_FILE, sheet = "energy")
      make_trait_heatmap(
        traits_energy,
        title_es = "Rasgos energéticos predichos por MAG",
        subtitle_tool = "DRAM2",
        outfile = file.path(OUTPUT_DIR, "Fig6a_traits_energy_ES.png"),
        legend_x_mm = 8, legend_y_mm = 10
      )
      make_trait_heatmap(
        traits_energy,
        title_es = "Rasgos energéticos predichos por MAG",
        subtitle_tool = "DRAM2",
        outfile = .cluster_out(file.path(OUTPUT_DIR, "Fig6a_traits_energy_ES.png")),
        annotation_version = 0L,
        legend_x_mm = 8, legend_y_mm = 10
      )
    },
    error = function(e) cat("  [WARN] Error en traits energy:", conditionMessage(e), "\n")
  )

  tryCatch(
    {
      traits_carbon <- read_excel(TRAITS_FILE, sheet = "carbon")
      make_trait_heatmap(
        traits_carbon,
        title_es = "Utilización de fuentes de carbono predicha por MAG",
        subtitle_tool = "DRAM2",
        outfile = file.path(OUTPUT_DIR, "Fig6b_traits_carbon_ES.png"),
        legend_x_mm = 8
      )
      make_trait_heatmap(
        traits_carbon,
        title_es = "Utilización de fuentes de carbono predicha por MAG",
        subtitle_tool = "DRAM2",
        outfile = .cluster_out(file.path(OUTPUT_DIR, "Fig6b_traits_carbon_ES.png")),
        annotation_version = 0L,
        legend_x_mm = 8
      )
    },
    error = function(e) cat("  [WARN] Error en traits carbon:", conditionMessage(e), "\n")
  )
  # ---- Fig6a and 6b, aggregated by order ----
  tryCatch(
    {
      traits_energy_ord <- read_excel(TRAITS_FILE, sheet = "energy")
      make_trait_order_heatmap(
        traits_energy_ord,
        title_es = "Rasgos energéticos predichos – proporción por Orden",
        subtitle_tool = "DRAM2",
        outfile = file.path(OUTPUT_DIR, "Fig6a_traits_energy_por_orden_ES.png")
      )
    },
    error = function(e) cat("  [WARN] Error en traits energy por Orden:", conditionMessage(e), "\n")
  )
  tryCatch(
    {
      traits_carbon_ord <- read_excel(TRAITS_FILE, sheet = "carbon")
      make_trait_order_heatmap(
        traits_carbon_ord,
        title_es = "Utilización de carbono predicha – proporción por Orden",
        subtitle_tool = "DRAM2",
        outfile = file.path(OUTPUT_DIR, "Fig6b_traits_carbon_por_orden_ES.png")
      )
    },
    error = function(e) cat("  [WARN] Error en traits carbon por Orden:", conditionMessage(e), "\n")
  )
  cat("[OK] Sección 13 completada\n\n")
} else {
  cat("[WARN] traits.xlsx no encontrado, omitiendo sección 13\n\n")
}

# ---- 13b. RASGOS AGREGADOS POR PROFUNDIDAD, ESTACIÓN Y GÉNERO TAXONÓMICO ----
cat("=== SECCIÓN 13b: Ver treatment_comparison.R ===\n")
# What used to live here, the aggregated heatmaps Fig9a to 9f along with the
# Fisher tests and logistic models, now lives in treatment_comparison.R.
cat("  (Gráficos 9a-9f generados por treatment_comparison.R)\n\n")

# ---- 14. RESUMEN METABÓLICO POR CATEGORÍA (metabolism_summary.xlsx) ----
cat("=== SECCIÓN 14: RESUMEN METABÓLICO (metabolism_summary.xlsx) ===\n")

# Completeness heatmap by subcategory, from one sheet of
# metabolism_summary.xlsx
make_metabolism_heatmap <- function(sheet_name, title_es, outfile,
                                    meta_cols = 5,
                                    group_col = "subcategory",
                                    make_group_label = NULL,
                                    subtitle_tool = NULL,
                                    mag_subset = NULL,
                                    include_station = FALSE,
                                    annotation_version = NULL,
                                    completeness_map = NULL,
                                    cache_file = DRAM_CACHE_FILE,
                                    max_rows = 60,
                                    legend_br = TRUE) {
  tryCatch(
    {
      df <- read_excel(METABOLISM_FILE, sheet = sheet_name)
      mag_cols <- names(df)[(meta_cols + 1):ncol(df)]
      meta_names <- names(df)[1:meta_cols]

      if (!is.null(mag_subset)) {
        mag_cols <- intersect(mag_cols, mag_subset)
        if (length(mag_cols) == 0) {
          cat("  [WARN] Sin MAGs del subconjunto en hoja", sheet_name, "\n")
          return(invisible(NULL))
        }
      }

      if (!is.null(make_group_label)) {
        df$.grp_label <- make_group_label(df)
        group_col <- ".grp_label"
      } else {
        if (!group_col %in% meta_names) {
          group_col <- intersect(c("category", "subcategory", "pathway"), meta_names)[1]
        }
        if (is.na(group_col)) {
          cat("  [WARN] No se encontró columna de agrupación en", sheet_name, "\n")
          return(invisible(NULL))
        }
      }

      df_long <- df %>%
        select(all_of(c(group_col, mag_cols))) %>%
        mutate(across(all_of(mag_cols), ~ as.numeric(. > 0))) %>%
        group_by(across(all_of(group_col))) %>%
        summarise(across(all_of(mag_cols), mean, na.rm = TRUE), .groups = "drop") %>%
        filter(rowSums(across(all_of(mag_cols))) > 0)

      if (nrow(df_long) == 0) {
        cat("  [WARN] Sin datos válidos en hoja", sheet_name, "\n")
        return(invisible(NULL))
      }

      if (nrow(df_long) > max_rows) {
        df_long <- df_long %>%
          mutate(total = rowSums(across(all_of(mag_cols)))) %>%
          slice_max(total, n = max_rows) %>%
          select(-total)
      }

      row_labels <- df_long[[group_col]]
      mat_meta <- df_long %>%
        select(all_of(mag_cols)) %>%
        as.matrix()
      rownames(mat_meta) <- row_labels

      mat_t <- t(mat_meta)
      mag_ids_meta <- rownames(mat_t)

      subcat_es <- capitalize_first(
        translate_to_spanish(colnames(mat_t), translation_cache_file = cache_file)
      )
      colnames(mat_t) <- subcat_es
      mat_t <- mat_t[, !is.na(colnames(mat_t)), drop = FALSE]

      anno_mt <- build_row_annotation(mag_ids_meta)
      dn_mt <- convert_to_display_name(mag_ids_meta)
      dn_to_id_mt <- setNames(mag_ids_meta, dn_mt)
      rownames(mat_t) <- dn_mt
      rownames(anno_mt) <- dn_mt
      anno_mt$Profundidad <- factor(anno_mt$Profundidad, levels = depth_levels_std)

      ord_mt <- order_matrix_by_depth(mat_t, anno_mt)
      mat_ord <- as.matrix(ord_mt$matrix)
      an_mt <- ord_mt$annotation

      dc_mt <- depth_colors_std[depth_levels_std]
      dc_mt <- dc_mt[!is.na(names(dc_mt))]

      v <- if (is.null(annotation_version)) {
        if (include_station) 2L else 3L
      } else {
        as.integer(annotation_version)
      }
      cluster_rows_mt <- (v == 0L)

      if (v %in% c(4L, 5L)) {
        ro <- reorder_rows_by_version(mat_ord, an_mt, v,
          dn_to_id = dn_to_id_mt,
          completeness_map = completeness_map
        )
        mat_ord <- ro$mat
        an_mt <- ro$ann
      }

      stations_mt <- if (v %in% c(1L, 2L, 4L)) {
        as.character(get_station(dn_to_id_mt[rownames(mat_ord)]))
      } else {
        NULL
      }
      compl_mt <- if (!is.null(completeness_map)) completeness_map[rownames(mat_ord)] else NULL

      ra_mt <- make_depth_row_annotation(
        an_mt$Profundidad, list(Profundidad = dc_mt),
        stations = stations_mt, annotation_version = v, completeness = compl_mt,
        show_legend = !legend_br
      )

      col_fun_mt <- colorRamp2(
        seq(0, 1, length.out = 6),
        c("white", "yellow", "orange", "orangered", "red", "darkred")
      )
      col_title_mt <- if (!is.null(subtitle_tool)) paste0(title_es, "\n", subtitle_tool) else title_es

      ht_mt <- Heatmap(
        mat_ord,
        name = "Completitud",
        col = col_fun_mt,
        cluster_rows = cluster_rows_mt, cluster_columns = TRUE,
        show_row_dend = cluster_rows_mt, row_dend_side = "left",
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "ward.D2",
        show_row_names = TRUE, show_column_names = TRUE,
        row_labels = make_taxon_label_expr(rownames(mat_ord), level = "Order"),
        row_names_gp = gpar_f(fontsize = 10),
        column_names_gp = gpar_f(fontsize = 7), column_names_rot = 45,
        left_annotation = ra_mt,
        rect_gp = gpar_f(col = "black", lwd = 0.4), border = TRUE,
        column_title = col_title_mt,
        column_title_gp = gpar_f(fontsize = 11, fontface = "bold"),
        show_heatmap_legend = !legend_br,
        heatmap_legend_param = list(
          title = "Completitud\npor categoría",
          at = c(0, 0.25, 0.5, 0.75, 1),
          labels = c("0%", "25%", "50%", "75%", "100%"),
          title_gp = gpar_f(fontsize = 10, fontface = "bold"),
          labels_gp = gpar_f(fontsize = 9),
          border = "black", legend_direction = "vertical"
        )
      )

      w_mt <- max(12, round(ncol(mat_ord) * 0.45 + 5))
      h_floor_mt <- min(8L, max(4L, as.integer(nrow(mat_ord) * 0.35)))
      h_mt <- max(h_floor_mt, round(nrow(mat_ord) * 0.20 + 3))
      png(outfile, width = w_mt, height = h_mt, units = "in", res = 300)
      draw(ht_mt)
      if (legend_br) {
        lgd_main <- Legend(
          col_fun   = col_fun_mt,
          title     = "Completitud\npor categoría",
          at        = c(0, 0.25, 0.5, 0.75, 1),
          labels    = c("0%", "25%", "50%", "75%", "100%"),
          title_gp  = gpar_f(fontsize = 10, fontface = "bold"),
          labels_gp = gpar_f(fontsize = 9), border = "black"
        )
        .draw_br_legend(c(list(lgd_main), .make_annotation_legends(dc_mt, stations_mt)), y_mm = 20)
      }
      dev.off()
      cat("  [OK] Guardado:", basename(outfile), "\n")
    },
    error = function(e) cat("  [WARN] Error en hoja", sheet_name, ":", conditionMessage(e), "\n")
  )
}

# Helper: heatmap de metales (top N genes, escala discreta 0/1/2 copias)
make_metals_heatmap <- function(title_es, outfile,
                                meta_cols = 5, top_n = 20,
                                subtitle_tool = NULL,
                                mag_subset = NULL,
                                include_station = FALSE,
                                annotation_version = NULL,
                                completeness_map = NULL,
                                cache_file = DRAM_CACHE_FILE,
                                legend_br = TRUE) {
  tryCatch(
    {
      df <- read_excel(METABOLISM_FILE, sheet = "metals")
      mag_cols <- names(df)[(meta_cols + 1):ncol(df)]

      if (!is.null(mag_subset)) {
        mag_cols <- intersect(mag_cols, mag_subset)
        if (length(mag_cols) == 0) {
          return(invisible(NULL))
        }
      }

      gene_descs <- as.character(df$gene_description)
      mat_metal <- df %>%
        select(all_of(mag_cols)) %>%
        mutate(across(everything(), ~ suppressWarnings(as.numeric(as.character(.))))) %>%
        as.matrix()
      rownames(mat_metal) <- gene_descs

      # The top N genes by total abundance
      gene_totals <- rowSums(mat_metal, na.rm = TRUE)
      top_rows <- order(gene_totals, decreasing = TRUE)[seq_len(min(top_n, nrow(mat_metal)))]
      mat_metal <- mat_metal[top_rows, , drop = FALSE]

      # Traducir etiquetas de genes
      gene_es <- capitalize_first(
        translate_to_spanish(rownames(mat_metal), translation_cache_file = cache_file)
      )
      gene_es[is.na(gene_es)] <- rownames(mat_metal)[is.na(gene_es)]
      rownames(mat_metal) <- gene_es

      # Quitar columnas todas-cero
      keep_cols <- colSums(mat_metal, na.rm = TRUE) > 0
      mat_metal <- mat_metal[, keep_cols, drop = FALSE]

      # Transponer: filas=MAGs, columnas=genes
      mat_t <- t(mat_metal)
      mag_ids_metal <- rownames(mat_t)

      anno_mt <- build_row_annotation(mag_ids_metal)
      dn_mt <- convert_to_display_name(mag_ids_metal)
      dn_to_id_mt <- setNames(mag_ids_metal, dn_mt)
      rownames(mat_t) <- dn_mt
      rownames(anno_mt) <- dn_mt
      anno_mt$Profundidad <- factor(anno_mt$Profundidad, levels = depth_levels_std)

      ord_mt <- order_matrix_by_depth(mat_t, anno_mt)
      mat_ord <- as.matrix(ord_mt$matrix)
      an_mt <- ord_mt$annotation

      dc_mt <- depth_colors_std[depth_levels_std]
      dc_mt <- dc_mt[!is.na(names(dc_mt))]

      v <- if (is.null(annotation_version)) {
        if (include_station) 2L else 3L
      } else {
        as.integer(annotation_version)
      }
      cluster_rows_me <- (v == 0L)

      if (v %in% c(4L, 5L)) {
        ro <- reorder_rows_by_version(mat_ord, an_mt, v,
          dn_to_id = dn_to_id_mt,
          completeness_map = completeness_map
        )
        mat_ord <- ro$mat
        an_mt <- ro$ann
      }

      # Integer then character, so the palette is discrete rather than continuous
      mat_int <- matrix(as.integer(pmin(mat_ord, 2, na.rm = TRUE)),
        nrow = nrow(mat_ord), dimnames = dimnames(mat_ord)
      )
      mat_int[is.na(mat_int)] <- 0L
      mat_chr <- matrix(as.character(mat_int), nrow = nrow(mat_int), dimnames = dimnames(mat_int))

      stations_mt <- if (v %in% c(1L, 2L, 4L)) {
        as.character(get_station(dn_to_id_mt[rownames(mat_ord)]))
      } else {
        NULL
      }
      compl_mt <- if (!is.null(completeness_map)) completeness_map[rownames(mat_ord)] else NULL

      ra_mt <- make_depth_row_annotation(
        an_mt$Profundidad, list(Profundidad = dc_mt),
        stations = stations_mt, annotation_version = v, completeness = compl_mt,
        show_legend = !legend_br
      )

      col_metal <- c("0" = "white", "1" = "#6BAED6", "2" = "#08519C")
      col_title_mt <- if (!is.null(subtitle_tool)) paste0(title_es, "\n", subtitle_tool) else title_es

      ht_mt <- Heatmap(
        mat_chr,
        name = "Copias",
        col = col_metal,
        cluster_rows = cluster_rows_me, cluster_columns = TRUE,
        show_row_dend = cluster_rows_me, row_dend_side = "left",
        clustering_distance_columns = "euclidean",
        clustering_method_columns = "ward.D2",
        show_row_names = TRUE, show_column_names = TRUE,
        row_labels = make_taxon_label_expr(rownames(mat_ord), level = "Order"),
        row_names_gp = gpar_f(fontsize = 10),
        column_names_gp = gpar_f(fontsize = 7), column_names_rot = 45,
        left_annotation = ra_mt,
        rect_gp = gpar_f(col = "black", lwd = 0.4), border = TRUE,
        column_title = col_title_mt,
        column_title_gp = gpar_f(fontsize = 11, fontface = "bold"),
        show_heatmap_legend = !legend_br,
        heatmap_legend_param = list(
          title = "Copias de\ngen",
          at = c("0", "1", "2"), labels = c("Ausente", "1 copia", "≥2 copias"),
          title_gp = gpar_f(fontsize = 10, fontface = "bold"),
          labels_gp = gpar_f(fontsize = 9),
          border = "black", legend_direction = "vertical"
        )
      )

      w_mt <- max(12, round(ncol(mat_ord) * 0.45 + 5))
      h_floor_mt <- min(8L, max(4L, as.integer(nrow(mat_ord) * 0.35)))
      h_mt <- max(h_floor_mt, round(nrow(mat_ord) * 0.20 + 3))
      png(outfile, width = w_mt, height = h_mt, units = "in", res = 300)
      draw(ht_mt)
      if (legend_br) {
        lgd_main <- Legend(
          labels    = c("Ausente", "1 copia", "≥2 copias"),
          legend_gp = gpar_f(fill = c("white", "#6BAED6", "#08519C")),
          title     = "Copias de\ngen",
          title_gp  = gpar_f(fontsize = 10, fontface = "bold"),
          labels_gp = gpar_f(fontsize = 9), border = "black"
        )
        .draw_br_legend(c(list(lgd_main), .make_annotation_legends(dc_mt, stations_mt)), y_mm = 20)
      }
      dev.off()
      cat("  [OK] Guardado:", basename(outfile), "\n")
    },
    error = function(e) cat("  [WARN] Error en metals heatmap:", conditionMessage(e), "\n")
  )
}

if (file.exists(METABOLISM_FILE)) {
  # Produces both the depth variant and the cladogram variant
  .make_meta_both <- function(...) {
    args <- list(...)
    # with the depth sidebar, the default
    do.call(make_metabolism_heatmap, args)
    # cladogram: no sidebar, rows clustered
    args$annotation_version <- 0L
    args$outfile <- .cluster_out(args$outfile)
    do.call(make_metabolism_heatmap, args)
  }
  .make_metals_both <- function(...) {
    args <- list(...)
    do.call(make_metals_heatmap, args)
    args$annotation_version <- 0L
    args$outfile <- .cluster_out(args$outfile)
    do.call(make_metals_heatmap, args)
  }

  .make_meta_both("Energy",
    meta_cols = 5, group_col = "pathway",
    title_es = "Completitud de genes de metabolismo energético por MAG",
    subtitle_tool = "DRAM2",
    outfile = file.path(OUTPUT_DIR, "Fig7a_metabolism_energy_ES.png")
  )
  .make_meta_both("carbon utilization",
    meta_cols = 5,
    title_es = "Completitud de genes de utilización de carbono por MAG",
    subtitle_tool = "DRAM2",
    outfile = file.path(OUTPUT_DIR, "Fig7b_metabolism_carbon_ES.png")
  )
  .make_meta_both("carbon utilization (Woodcroft)",
    meta_cols = 5,
    title_es = "Utilización de carbono (Woodcroft) por MAG",
    subtitle_tool = "DRAM2",
    outfile = file.path(OUTPUT_DIR, "Fig7c_metabolism_carbon_woodcroft_ES.png")
  )
  # Fig7d: metals – top 20 genes, escala discreta 0/1/2 copias
  .make_metals_both(
    title_es = "Genes de reducción de metales por MAG (top 20)",
    subtitle_tool = "DRAM2",
    outfile = file.path(OUTPUT_DIR, "Fig7d_metabolism_metals_ES.png")
  )
  .make_meta_both("CAMPER",
    meta_cols = 5, group_col = "pathway",
    title_es = "Genes CAMPER (síntesis de polifenoles) por MAG",
    subtitle_tool = "DRAM2",
    outfile = file.path(OUTPUT_DIR, "Fig7e_metabolism_CAMPER_ES.png")
  )
  .make_meta_both("Organic Nitrogen",
    meta_cols = 5, group_col = "pathway",
    title_es = "Rutas metabólicas asociadas a nitrógeno orgánico por MAG",
    subtitle_tool = "DRAM2",
    outfile = file.path(OUTPUT_DIR, "Fig7f_metabolism_organic_N_ES.png")
  )
  .make_meta_both("Transporters",
    meta_cols = 5, group_col = "pathway",
    title_es = "Transportadores de biomoléculas por MAG",
    subtitle_tool = "DRAM2",
    outfile = file.path(OUTPUT_DIR, "Fig7g_metabolism_transporters_ES.png")
  )
  .misc_grp_fn_54 <- function(df) {
    sub_val <- as.character(df$subcategory)
    path_val <- as.character(df$pathway)
    cat_val <- as.character(df$category)
    ifelse(!is.na(sub_val) & trimws(sub_val) != "" & trimws(sub_val) != "NA",
      paste0(cat_val, ": ", sub_val), path_val
    )
  }
  .make_meta_both("MISC",
    meta_cols = 5, make_group_label = .misc_grp_fn_54,
    title_es = "Funciones misceláneas por MAG",
    subtitle_tool = "DRAM2",
    outfile = file.path(OUTPUT_DIR, "Fig7h_metabolism_misc_ES.png")
  )
  cat("[OK] Sección 14 completada\n\n")
} else {
  cat("[WARN] metabolism_summary.xlsx no encontrado, omitiendo sección 14\n\n")
}

# ---- 15. HEATMAPS POR CATEGORÍA DE GENES (summarized_genomes.tsv) ----
cat("=== SECCIÓN 15: HEATMAPS POR CATEGORÍA DE GENES (summarized_genomes.tsv) ===\n")

if (file.exists(SUMMARIZED_GENOMES_FILE)) {
  sg <- read_tsv(SUMMARIZED_GENOMES_FILE, show_col_types = FALSE)
  # Columnas meta: gene_id, gene_description, topic_ecosystem, category, subcategory, pathway
  meta_sg <- names(sg)[1:6]
  mag_cols_sg <- names(sg)[7:ncol(sg)] # NameWithoutExtension de cada MAG
  cat("  Dimensiones summarized_genomes:", nrow(sg), "x", ncol(sg), "\n")
  cat("  MAGs en columnas:", length(mag_cols_sg), "\n")

  # Presence and absence heatmap for individual genes, by category
  make_summarized_heatmap <- function(filter_expr, title_es, outfile,
                                      max_genes = 80,
                                      subtitle_tool = NULL,
                                      mag_subset = NULL,
                                      include_station = FALSE,
                                      annotation_version = NULL,
                                      completeness_map = NULL,
                                      cache_file = DRAM_CACHE_FILE,
                                      legend_br = TRUE,
                                      legend_x_mm = 18, legend_y_mm = 20) {
    tryCatch(
      {
        sub_sg <- sg %>% filter(!!filter_expr)
        if (nrow(sub_sg) == 0) {
          cat("  [WARN] Sin genes para:", title_es, "\n")
          return(invisible(NULL))
        }

        active_mag_cols <- if (!is.null(mag_subset)) intersect(mag_cols_sg, mag_subset) else mag_cols_sg
        if (length(active_mag_cols) == 0) {
          cat("  [WARN] Sin MAGs del subconjunto para:", title_es, "\n")
          return(invisible(NULL))
        }

        mat_sg <- sub_sg %>%
          select(all_of(active_mag_cols)) %>%
          mutate(across(everything(), ~ as.numeric(. > 0))) %>%
          as.matrix()

        gene_labels <- capitalize_first(trimws(sub("\\s*\\[.*$", "", sub_sg$gene_description)))
        gene_labels[is.na(gene_labels) | gene_labels == ""] <-
          sub_sg$gene_id[is.na(gene_labels) | gene_labels == ""]
        rownames(mat_sg) <- make_unique_names(gene_labels)

        mat_sg <- mat_sg[rowSums(mat_sg, na.rm = TRUE) > 0, , drop = FALSE]
        if (nrow(mat_sg) == 0) {
          cat("  [WARN] Todos los genes ausentes\n")
          return(invisible(NULL))
        }

        mat_sg <- mat_sg[rowSums(mat_sg, na.rm = TRUE) >= 2, , drop = FALSE]
        if (nrow(mat_sg) == 0) {
          cat("  [WARN] Sin genes en ≥2 MAGs\n")
          return(invisible(NULL))
        }

        if (nrow(mat_sg) > max_genes) {
          mat_sg <- mat_sg[order(rowSums(mat_sg), decreasing = TRUE)[1:max_genes], , drop = FALSE]
        }

        mat_t <- t(mat_sg)
        mag_ids_sg <- rownames(mat_t)
        anno_sg <- build_row_annotation(mag_ids_sg)
        dn_sg <- convert_to_display_name(mag_ids_sg)
        dn_to_id_sg <- setNames(mag_ids_sg, dn_sg)
        rownames(mat_t) <- dn_sg
        rownames(anno_sg) <- dn_sg
        anno_sg$Profundidad <- factor(anno_sg$Profundidad, levels = depth_levels_std)

        ord_sg <- order_matrix_by_depth(mat_t, anno_sg)
        mat_ord_sg <- as.matrix(ord_sg$matrix)
        an_sg <- ord_sg$annotation

        dc_sg <- depth_colors_std[depth_levels_std]
        dc_sg <- dc_sg[!is.na(names(dc_sg))]

        v <- if (is.null(annotation_version)) {
          if (include_station) 2L else 3L
        } else {
          as.integer(annotation_version)
        }
        cluster_rows_sg <- (v == 0L)

        if (v %in% c(4L, 5L)) {
          ro <- reorder_rows_by_version(mat_ord_sg, an_sg, v,
            dn_to_id = dn_to_id_sg,
            completeness_map = completeness_map
          )
          mat_ord_sg <- ro$mat
          an_sg <- ro$ann
        }

        stations_sg <- if (v %in% c(1L, 2L, 4L)) {
          as.character(get_station(dn_to_id_sg[rownames(mat_ord_sg)]))
        } else {
          NULL
        }
        compl_sg <- if (!is.null(completeness_map)) completeness_map[rownames(mat_ord_sg)] else NULL

        ra_sg <- make_depth_row_annotation(
          an_sg$Profundidad, list(Profundidad = dc_sg),
          stations = stations_sg, annotation_version = v, completeness = compl_sg,
          show_legend = !legend_br
        )

        col_sg <- c("0" = "#F8F9FA", "1" = "#264653")
        col_title_sg <- if (!is.null(subtitle_tool)) paste0(title_es, "\n", subtitle_tool) else title_es

        ht_sg <- Heatmap(
          mat_ord_sg,
          name = "Presencia",
          col = col_sg,
          cluster_rows = cluster_rows_sg, cluster_columns = TRUE,
          show_row_dend = cluster_rows_sg, row_dend_side = "left",
          clustering_distance_columns = "binary",
          clustering_method_columns = "ward.D2",
          show_row_names = TRUE, show_column_names = TRUE,
          row_labels = make_taxon_label_expr(rownames(mat_ord_sg), level = "Order"),
          row_names_gp = gpar_f(fontsize = 9),
          column_names_gp = gpar_f(fontsize = 6), column_names_rot = 75,
          left_annotation = ra_sg,
          rect_gp = gpar_f(col = "black", lwd = 0.3), border = TRUE,
          column_title = col_title_sg,
          column_title_gp = gpar_f(fontsize = 10, fontface = "bold"),
          show_heatmap_legend = !legend_br,
          heatmap_legend_param = list(
            title = "Gen",
            at = c(0, 1), labels = c("Ausente", "Presente"),
            title_gp = gpar_f(fontsize = 10, fontface = "bold"),
            labels_gp = gpar_f(fontsize = 9),
            border = "black", legend_direction = "vertical"
          )
        )

        w_sg <- max(16, round(ncol(mat_ord_sg) * 0.18 + 6))
        h_floor_sg <- min(8L, max(4L, as.integer(nrow(mat_ord_sg) * 0.35)))
        h_sg <- max(h_floor_sg, round(nrow(mat_ord_sg) * 0.20 + 3))
        png(outfile, width = w_sg, height = h_sg, units = "in", res = 300)
        draw(ht_sg)
        if (legend_br) {
          lgd_main <- Legend(
            labels    = c("Ausente", "Presente"),
            legend_gp = gpar_f(fill = c("#F8F9FA", "#264653")),
            title     = "Gen",
            title_gp  = gpar_f(fontsize = 10, fontface = "bold"),
            labels_gp = gpar_f(fontsize = 9), border = "black"
          )
          .draw_br_legend(c(list(lgd_main), .make_annotation_legends(dc_sg, stations_sg)),
                          x_mm = legend_x_mm, y_mm = legend_y_mm)
        }
        dev.off()
        cat("  [OK] Guardado:", basename(outfile), "\n")
      },
      error = function(e) cat("  [WARN] Error en", title_es, ":", conditionMessage(e), "\n")
    )
  }

  .sub_tool_15 <- "DRAM2 – summarized_genomes.tsv"

  # Produces the depth and cladogram variants of make_summarized_heatmap
  .make_sg_both <- function(filter_expr, title_es, outfile, ...) {
    make_summarized_heatmap(filter_expr, title_es = title_es, outfile = outfile, ...)
    make_summarized_heatmap(filter_expr,
      title_es = title_es,
      outfile = .cluster_out(outfile),
      annotation_version = 0L, ...
    )
  }

  # --- 15.1 Ciclo TCA y carbono central ---
  .make_sg_both(expr(category == "central carbon"),
    title_es = "Genes de metabolismo central de carbono por MAG\n(ciclo TCA, glucólisis, pentosa-fosfato, etc.)",
    outfile = file.path(OUTPUT_DIR, "Fig8a_genes_central_carbon_ES.png"),
    subtitle_tool = .sub_tool_15
  )

  # ---- 15.2 electron transport chain ----
  .make_sg_both(expr(category == "Electron transport Chain"),
    title_es = "Genes de la cadena de transporte electrónico por MAG",
    outfile = file.path(OUTPUT_DIR, "Fig8b_genes_ETC_ES.png"),
    subtitle_tool = .sub_tool_15
  )

  # ---- 15.3 nitrogen cycle ----
  .make_sg_both(expr(category == "Nitrogen"),
    title_es = "Genes del ciclo del nitrógeno por MAG",
    outfile = file.path(OUTPUT_DIR, "Fig8c_genes_nitrogen_ES.png"),
    subtitle_tool = .sub_tool_15
  )

  # ---- 15.4 sulfur cycle ----
  .make_sg_both(expr(category == "Sulfur"),
    title_es = "Genes del ciclo del azufre por MAG",
    outfile = file.path(OUTPUT_DIR, "Fig8d_genes_sulfur_ES.png"),
    subtitle_tool = .sub_tool_15
  )

  # --- 15.5 Metabolismo de C1/metano ---
  .make_sg_both(expr(category %in% c("C1-methane", "C1")),
    title_es = "Genes de metabolismo de C1/metano por MAG\n(metanogénesis, metanotrofía, metabolismo de metilaminas)",
    outfile = file.path(OUTPUT_DIR, "Fig8e_genes_C1_methane_ES.png"),
    subtitle_tool = .sub_tool_15
  )

  # --- 15.6 Hidrogenasas ---
  .make_sg_both(expr(category == "Hydrogenases"),
    title_es = "Hidrogenasas por MAG",
    outfile = file.path(OUTPUT_DIR, "Fig8f_genes_hydrogenases_ES.png"),
    subtitle_tool = .sub_tool_15
  )

  # ---- 15.7 CAZymes, polysaccharide degradation ----
  .make_sg_both(expr(category == "CAZY"),
    title_es = "Enzimas activas en carbohidratos (CAZymes) por MAG\n(degradación de polisacáridos y glicoproteínas)",
    outfile = file.path(OUTPUT_DIR, "Fig8g_genes_CAZymes_ES.png"),
    subtitle_tool = .sub_tool_15, max_genes = 60
  )

  # ---- 15.8 metal reduction ----
  .make_sg_both(expr(category == "Metal Reduction"),
    title_es = "Genes de reducción de metales por MAG",
    outfile = file.path(OUTPUT_DIR, "Fig8h_genes_metal_reduction_ES.png"),
    subtitle_tool = .sub_tool_15
  )

  # ---- 15.9 pyruvate metabolism ----
  .make_sg_both(expr(category == "pyruvate metabolism"),
    title_es = "Genes del metabolismo del piruvato por MAG",
    outfile = file.path(OUTPUT_DIR, "Fig8i_genes_pyruvate_ES.png"),
    subtitle_tool = .sub_tool_15
  )

  # ---- 15.10 short-chain fatty acid conversions and fermentation ----
  .make_sg_both(expr(category == "SCFA  and alcohol conversions"),
    title_es = "Genes de conversión de ácidos grasos de cadena corta\ny fermentación por MAG",
    outfile = file.path(OUTPUT_DIR, "Fig8j_genes_SCFA_fermentation_ES.png"),
    subtitle_tool = .sub_tool_15
  )

  # ---- 15.11 photosynthesis ----
  .make_sg_both(expr(category == "Photosynthesis"),
    title_es = "Genes de fotosíntesis / fototropía por MAG",
    outfile = file.path(OUTPUT_DIR, "Fig8k_genes_photosynthesis_ES.png"),
    subtitle_tool = .sub_tool_15,
    legend_y_mm = 10
  )

  # ---- 15.12 transposition and information systems ----
  .make_sg_both(expr(category == "Information systems"),
    title_es = "Genes de sistemas de información (replicación, reparación) por MAG",
    outfile = file.path(OUTPUT_DIR, "Fig8l_genes_info_systems_ES.png"),
    subtitle_tool = .sub_tool_15, max_genes = 50
  )

  cat("[OK] Sección 15 completada\n\n")
} else {
  cat("[WARN] summarized_genomes.tsv no encontrado, omitiendo sección 15\n\n")
}

cat("\n", rep("=", 70), "\n", sep = "")
cat("RESUMEN DEL ANÁLISIS COMPLETO\n")
cat(rep("=", 70), "\n", sep = "")
cat("\nEstadísticas generales:\n")
cat("  - Total de genes:", comma(nrow(dram)), "\n")
cat("  - Total de MAGs:", n_distinct(dram$input_fasta), "\n")
cat("\nClasificación por calidad de anotación:\n")
cat(
  "  - Anotación funcional curada (KEGG/KOfam/dbCAZy/CAMPER/sulfur/methyl/FeGenie/metals):",
  comma(summary_table$Curated_Functional_Annotation),
  paste0("(", round(summary_table$Percent_Curated_Functional, 2), "%)"), "\n"
)
cat(
  "  - Sólo homología UniRef (sin función curada):", comma(summary_table$Homology_Only),
  paste0("(", round(summary_table$Percent_Homology_Only, 2), "%)"), "\n"
)
cat(
  "      · Descripción UniRef informativa:", comma(summary_table$Homology_Only_Informative), "\n"
)
cat(
  "      · Descripción UniRef vaga (uncharacterized/DUF/hypothetical):",
  comma(summary_table$Homology_Only_Vague), "\n"
)
cat(
  "  - Sin anotación:", comma(summary_table$No_Annotation),
  paste0("(", round(summary_table$Percent_No_Annotation, 2), "%)"), "\n"
)
cat(
  "  - Subtotal con cualquier hit (curado + homología):",
  paste0("(", round(summary_table$Percent_Any_Hit, 2), "%)"), "\n"
)

cat("\nDistribución de hits por base de datos (KEGG/HMM curadas + UniRef):\n")
cat("  - KEGG KOfam:", comma(summary_table$KEGG_KOfam), "\n")
cat("  - dbCAZy:", comma(summary_table$dbCAZy), "\n")
cat("  - UniRef (informativo):", comma(summary_table$UniRef_Informative), "\n")
cat("  - Metals:", comma(summary_table$Metals), "\n")
cat("  - CAMPER:", comma(summary_table$CAMPER), "\n")
cat("  - Sulfur:", comma(summary_table$Sulfur), "\n")
cat("  - FeGenie:", comma(summary_table$FeGenie), "\n")
cat("  - Methyl:", comma(summary_table$Methyl), "\n")

cat("\nArchivos generados:\n")
cat("  - 4 tablas flextable (PNG + PPTX)\n")
cat("  - Heatmaps anotaciones DRAM2 (kofam, dbCAZy, metals, CAMPER, sulfur, FeGenie, methyl)\n")
cat("  - Fig1_kegg_decoder_ES.png        : KEGG-Decoder (completitud de rutas)\n")
cat("  - Fig2_kegg_expander_ES.png       : KEGG-Expander (genes adicionales)\n")
cat("  - Fig3_kegg_combined_ES.png       : Combinación KEGG-Decoder + Expander\n")
cat("  - Fig4_product_completeness_ES.png: Completitud de rutas DRAM2 product.tsv\n")
cat("  - Fig6a_traits_energy_ES.png      : Rasgos energéticos predichos (binario)\n")
cat("  - Fig6b_traits_carbon_ES.png      : Utilización de carbono predicha (binario)\n")
cat("  - Fig7a-h_metabolism_*_ES.png     : Completitud por categoría (metabolism_summary.xlsx)\n")
cat("  - Fig8a-l_genes_*_ES.png          : Genes por categoría biológica (summarized_genomes.tsv)\n")
cat("  - 9 gráficos generales de análisis\n")
cat("  - Datos CSV y tablas resumen\n")

cat("\nTodos los resultados en:\n")
cat("  ", OUTPUT_DIR, "\n")
cat(rep("=", 70), "\n", sep = "")
cat("\n[OK] Análisis completado exitosamente!\n\n")

# ---- 16. ANÁLISIS PARA LOS 10 MAGs SIN TAXONOMÍA A NIVEL ESPECIE ----
cat("\n=== SECCIÓN 16: ANÁLISIS PARA 10 MAGs SIN TAXONOMÍA A NIVEL ESPECIE ===\n")

# The MAGs with no species-level assignment
mags_no_species <- metadata %>%
  filter(is.na(Species) | trimws(as.character(Species)) %in% c("", "N/A")) %>%
  pull(NameWithoutExtension)

cat("MAGs sin taxonomía a nivel especie:", length(mags_no_species), "\n")
if (length(mags_no_species) == 0) {
  cat("[WARN] No se encontraron MAGs sin especie. Omitiendo sección 16.\n\n")
} else {
  OUTPUT_DIR_ORIG <- OUTPUT_DIR
  OUTPUT_DIR <- file.path(OUTPUT_DIR_ORIG, "10_MAGs")
  dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
  cat("Subcarpeta 10_MAGs:", OUTPUT_DIR, "\n\n")

  # Shares the translation cache with the main run
  DRAM_CACHE_FILE_10 <- file.path(
    OUTPUT_DIR_ORIG, "translation_cache_dram2_distill.csv"
  )

  # ---- depths for the 10-MAG subset, which has no ROMZ ----
  depth_levels_std_orig <- depth_levels_std
  depth_colors_std_orig <- depth_colors_std
  depth_colors_orig <- depth_colors

  # Which depths actually hold MAGs in this subset
  depths_present_10 <- unique(na.omit(get_depth(mags_no_species)))
  depth_levels_std <- intersect(c("1000 m", "Fondo"), depths_present_10)
  depth_colors_std <- depth_colors_std_orig[depth_levels_std]
  depth_colors <- depth_colors_orig[depth_levels_std]
  cat("  Profundidades en 10 MAGs:", paste(depth_levels_std, collapse = ", "), "\n\n")

  # ---- Estaciones presentes en 10 MAGs ----
  stations_present_10 <- sort(unique(na.omit(get_station(mags_no_species))))
  n_st_10 <- max(length(stations_present_10), 3)
  station_colors_10 <- setNames(
    RColorBrewer::brewer.pal(n_st_10, "Set1")[seq_len(length(stations_present_10))],
    stations_present_10
  )

  # ---- prepare the KEGG Decoder data once ----
  kegg_dec_10_ready <- NULL
  if (file.exists(DECODER_FILE)) {
    dec_10_raw <- read_kegg_list(DECODER_FILE)
    common_dec_10 <- intersect(rownames(dec_10_raw), mags_no_species)
    if (length(common_dec_10) >= 2) {
      dec_10_raw <- filter_kegg_zero_cols(dec_10_raw[common_dec_10, , drop = FALSE])
      colnames(dec_10_raw) <- capitalize_first(translate_to_spanish(colnames(dec_10_raw),
        translation_cache_file = file.path(OUTPUT_DIR_ORIG, "translation_cache_decoder.csv")
      ))
      dec_10_raw <- dec_10_raw[, !is.na(colnames(dec_10_raw)), drop = FALSE]
      meta_dec_10 <- load_kegg_metadata(common_dec_10)
      dec_10_raw <- kegg_rename_rows(dec_10_raw, meta_dec_10)
      rownames(meta_dec_10) <- meta_dec_10$DisplayName
      ord_d10 <- kegg_order_by_depth(dec_10_raw, meta_dec_10)
      kegg_dec_10_ready <- list(mat = ord_d10$mat, meta = ord_d10$metadata)
      cat("  [OK] KEGG-Decoder pre-procesado\n")
    }
  }

  # ---- prepare the KEGG Expander data once ----
  kegg_exp_10_ready <- NULL
  if (file.exists(EXPANDER_FILE)) {
    exp_10_raw <- read_kegg_list(EXPANDER_FILE)
    common_exp_10 <- intersect(rownames(exp_10_raw), mags_no_species)
    if (length(common_exp_10) >= 2) {
      exp_10_raw <- filter_kegg_zero_cols(exp_10_raw[common_exp_10, , drop = FALSE])
      colnames(exp_10_raw) <- capitalize_first(translate_to_spanish(colnames(exp_10_raw),
        translation_cache_file = file.path(OUTPUT_DIR_ORIG, "translation_cache_expander.csv")
      ))
      meta_exp_10 <- load_kegg_metadata(common_exp_10)
      exp_10_raw <- kegg_rename_rows(exp_10_raw, meta_exp_10)
      rownames(meta_exp_10) <- meta_exp_10$DisplayName
      ord_e10 <- kegg_order_by_depth(exp_10_raw, meta_exp_10)
      kegg_exp_10_ready <- list(
        mat = ord_e10$mat, meta = ord_e10$metadata,
        max_val = max(ord_e10$mat, na.rm = TRUE)
      )
      cat("  [OK] KEGG-Expander pre-procesado\n")
    }
  }

  # ---- prepare product.tsv once ----
  prod_10_ready <- NULL
  if (file.exists(PRODUCT_FILE) && file.exists(GENOME_STATS_FILE)) {
    gs_10 <- read_tsv(GENOME_STATS_FILE, show_col_types = FALSE)
    mag_order_10 <- gs_10$genome
    prod_10_raw <- read_tsv(PRODUCT_FILE, show_col_types = FALSE)
    if (nrow(prod_10_raw) == length(mag_order_10)) {
      prod_10_raw$genome <- mag_order_10
    } else {
      prod_10_raw$genome <- paste0("MAG_", seq_len(nrow(prod_10_raw)))
    }
    prod_10_raw <- prod_10_raw %>% filter(genome %in% mags_no_species)
    prod_10_ready <- list(
      raw = prod_10_raw,
      ids = prod_10_raw$genome,
      nonempty = {
        all_c <- setdiff(names(prod_10_raw), "genome")
        all_c[vapply(all_c, function(cn) {
          v <- as.character(prod_10_raw[[cn]])
          any(v %in% c("True", "TRUE") | suppressWarnings(as.numeric(v)) > 0, na.rm = TRUE)
        }, logical(1))]
      }
    )
    cat("  [OK] product.tsv pre-procesado\n")
  }

  # ---- completeness and contamination sidebars for the 10-MAG product.tsv ----
  completeness_map_10 <- NULL
  contamination_map_10 <- NULL
  if (file.exists(QUALITY_REPORT_FILE)) {
    orig_name_col_10 <- names(metadata)[11]  # OriginalName
    meta_orig_10 <- metadata %>%
      select(NameWithoutExtension, OriginalName = all_of(orig_name_col_10))
    qual_10_df <- tryCatch(
      read_tsv(QUALITY_REPORT_FILE, show_col_types = FALSE) %>%
        inner_join(meta_orig_10, by = c("Name" = "OriginalName")) %>%
        filter(NameWithoutExtension %in% mags_no_species),
      error = function(e) NULL
    )
    if (!is.null(qual_10_df) && nrow(qual_10_df) > 0) {
      dn_keys_10 <- convert_to_display_name(qual_10_df$NameWithoutExtension)
      completeness_map_10 <- setNames(as.numeric(qual_10_df$Completeness), dn_keys_10)
      contamination_map_10 <- setNames(as.numeric(qual_10_df$Contamination), dn_keys_10)
      cat("  [OK] Datos de calidad (product 10 MAGs):", nrow(qual_10_df), "MAGs\n")
    }
  }

  # ---- Four sidebar variants ----
  ann_version_dirs <- c(
    "2" = "prof_estacion", # v2: profundidad + estación
    "3" = "profundidad", # v3: solo profundidad
    "4" = "estacion", # v4: solo estación
    "0" = "cladograma" # v0: sin barra, filas clusterizadas
  )

  for (ann_v_str in names(ann_version_dirs)) {
    ann_v <- as.integer(ann_v_str) # versión efectiva: 2L, 3L, 4L
    OUTPUT_DIR <- file.path(OUTPUT_DIR_ORIG, "10_MAGs", ann_version_dirs[ann_v_str])
    dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
    cat("\n--- Anotación versión", ann_v, "(", ann_version_dirs[ann_v_str], ") ---\n")

    # ---- KEGG Decoder ----
    if (!is.null(kegg_dec_10_ready)) {
      mat_d10 <- kegg_dec_10_ready$mat
      meta_d10 <- kegg_dec_10_ready$meta
      # Re-sort for variant 4 (station); variant 0 lets ComplexHeatmap cluster
      if (ann_v == 4L) {
        ord_d10_v <- order(meta_d10$Station)
        mat_d10 <- mat_d10[ord_d10_v, , drop = FALSE]
        meta_d10 <- meta_d10[ord_d10_v, , drop = FALSE]
      }
      dc_d10 <- kegg_named_colors(meta_d10$Depth, "Set2")
      dc_d10 <- dc_d10[names(dc_d10) != "ZMOR"]
      ra_d10 <- make_depth_row_annotation(
        meta_d10$Depth, list(Profundidad = dc_d10),
        stations = meta_d10$Station,
        annotation_version = ann_v
      )
      cf_d10 <- colorRamp2(
        c(0, 0.2, 0.4, 0.6, 0.8, 1),
        c("white", "yellow", "orange", "orangered", "red", "darkred")
      )
      ht_d10 <- Heatmap(mat_d10,
        name = "Completitud", col = cf_d10,
        cluster_rows = (ann_v == 0L), cluster_columns = TRUE,
        show_row_dend = (ann_v == 0L), row_dend_side = "left",
        clustering_distance_columns = "euclidean", clustering_method_columns = "ward.D2",
        show_row_names = TRUE, show_column_names = TRUE,
        row_labels = make_taxon_label_expr(rownames(mat_d10)),
        row_names_gp = gpar_f(fontsize = 10),
        column_names_gp = gpar_f(fontsize = 8), column_names_rot = 45,
        left_annotation = ra_d10,
        rect_gp = gpar_f(col = "black", lwd = 0.5), border = TRUE,
        column_title = "Completitud de rutas metabólicas – 10 MAGs sin especie\nKEGG-Decoder",
        column_title_gp = gpar_f(fontsize = 11, fontface = "bold"),
        heatmap_legend_param = list(
          title = "Completitud\nde rutas",
          at = c(0, 0.25, 0.5, 0.75, 1),
          labels = c("0%", "25%", "50%", "75%", "100%"),
          title_gp = gpar_f(fontsize = 10, fontface = "bold"),
          labels_gp = gpar_f(fontsize = 9), border = "black", legend_direction = "vertical"
        )
      )
      png(file.path(OUTPUT_DIR, "Fig1_kegg_decoder_ES.png"), width = 4200, height = 1500, res = 300)
      draw(ht_d10)
      dev.off()
      cat("  [OK] Fig1_kegg_decoder_ES.png\n")
    }

    # ---- KEGG Expander ----
    if (!is.null(kegg_exp_10_ready)) {
      mat_e10 <- kegg_exp_10_ready$mat
      meta_e10 <- kegg_exp_10_ready$meta
      max_e10 <- kegg_exp_10_ready$max_val
      # Re-sort for variant 4 (station); variant 0 lets ComplexHeatmap cluster
      if (ann_v == 4L) {
        ord_e10_v <- order(meta_e10$Station)
        mat_e10 <- mat_e10[ord_e10_v, , drop = FALSE]
        meta_e10 <- meta_e10[ord_e10_v, , drop = FALSE]
      }
      dc_e10 <- kegg_named_colors(meta_e10$Depth, "Set2")
      dc_e10 <- dc_e10[names(dc_e10) != "ZMOR"]
      ra_e10 <- make_depth_row_annotation(
        meta_e10$Depth, list(Profundidad = dc_e10),
        stations = meta_e10$Station,
        annotation_version = ann_v
      )
      cf_e10 <- colorRamp2(
        seq(0, max_e10, length.out = 6),
        c("white", "yellow", "orange", "orangered", "red", "darkred")
      )
      ht_e10 <- Heatmap(mat_e10,
        name = "Abundancia", col = cf_e10,
        cluster_rows = (ann_v == 0L), cluster_columns = TRUE,
        show_row_dend = (ann_v == 0L), row_dend_side = "left",
        show_row_names = TRUE, show_column_names = TRUE,
        row_labels = make_taxon_label_expr(rownames(mat_e10)),
        row_names_gp = gpar_f(fontsize = 9),
        column_names_gp = gpar_f(fontsize = 6), column_names_rot = 45,
        left_annotation = ra_e10,
        rect_gp = gpar_f(col = "black", lwd = 0.5), border = TRUE,
        column_title = "Abundancia de funciones metabólicas – 10 MAGs sin especie\nKEGG-Expander",
        column_title_gp = gpar_f(fontsize = 11, fontface = "bold"),
        heatmap_legend_param = list(
          title = "Abundancia de\nfunción",
          at = round(seq(0, max_e10, length.out = 5)),
          legend_direction = "vertical",
          title_gp = gpar_f(fontsize = 10, fontface = "bold"),
          labels_gp = gpar_f(fontsize = 9), border = "black"
        )
      )
      png(file.path(OUTPUT_DIR, "Fig2_kegg_expander_ES.png"), width = 4200, height = 1500, res = 300)
      draw(ht_e10)
      dev.off()
      cat("  [OK] Fig2_kegg_expander_ES.png\n")
    }

    # ---- product.tsv ----
    if (!is.null(prod_10_ready)) {
      cat("  Generando heatmaps product.tsv...\n")
      make_product_heatmap(
        prod_10_ready$raw %>% select(genome, all_of(prod_10_ready$nonempty)),
        prod_10_ready$ids,
        title_main = "Completitud de rutas metabólicas – 10 MAGs sin especie\nDRAM2",
        outfile = file.path(OUTPUT_DIR, "Fig4_product_completeness_ES.png"),
        annotation_version = ann_v,
        completeness_map = completeness_map_10,
        contamination_map = contamination_map_10,
        cache_f = DRAM_CACHE_FILE_10,
        barplot_gap = 5,
        extra_canvas_in = 3
      )
    }

    # ---- Traits ----
    cat("  Generando heatmaps de traits...\n")
    if (file.exists(TRAITS_FILE)) {
      tryCatch(
        {
          tr_energy_10 <- read_excel(TRAITS_FILE, sheet = "energy")
          make_trait_heatmap(tr_energy_10,
            title_es = "Rasgos energéticos – 10 MAGs sin especie",
            subtitle_tool = "DRAM2", outfile = file.path(OUTPUT_DIR, "Fig6a_traits_energy_ES.png"),
            mag_subset = mags_no_species,
            annotation_version = ann_v,
            cache_file = DRAM_CACHE_FILE_10,
            legend_br = FALSE
          )
        },
        error = function(e) cat("  [WARN] traits energy:", conditionMessage(e), "\n")
      )
      tryCatch(
        {
          tr_carbon_10 <- read_excel(TRAITS_FILE, sheet = "carbon")
          make_trait_heatmap(tr_carbon_10,
            title_es = "Utilización de carbono – 10 MAGs sin especie",
            subtitle_tool = "DRAM2", outfile = file.path(OUTPUT_DIR, "Fig6b_traits_carbon_ES.png"),
            mag_subset = mags_no_species,
            annotation_version = ann_v,
            cache_file = DRAM_CACHE_FILE_10,
            legend_br = FALSE
          )
        },
        error = function(e) cat("  [WARN] traits carbon:", conditionMessage(e), "\n")
      )
    }

    # ---- Metabolism ----
    cat("  Generando heatmaps de metabolismo...\n")
    if (file.exists(METABOLISM_FILE)) {
      .misc_grp_fn <- function(df) {
        sub_val <- as.character(df$subcategory)
        path_val <- as.character(df$pathway)
        cat_val <- as.character(df$category)
        ifelse(!is.na(sub_val) & trimws(sub_val) != "" & trimws(sub_val) != "NA",
          paste0(cat_val, ": ", sub_val), path_val
        )
      }
      make_metabolism_heatmap("Energy",
        meta_cols = 5, group_col = "pathway",
        subtitle_tool = "DRAM2", mag_subset = mags_no_species,
        annotation_version = ann_v,
        cache_file = DRAM_CACHE_FILE_10,
        legend_br = FALSE,
        title_es = "Metabolismo energético – 10 MAGs sin especie",
        outfile = file.path(OUTPUT_DIR, "Fig7a_metabolism_energy_ES.png")
      )
      make_metabolism_heatmap("carbon utilization",
        meta_cols = 5,
        subtitle_tool = "DRAM2", mag_subset = mags_no_species,
        annotation_version = ann_v,
        cache_file = DRAM_CACHE_FILE_10,
        legend_br = FALSE,
        title_es = "Utilización de carbono – 10 MAGs sin especie",
        outfile = file.path(OUTPUT_DIR, "Fig7b_metabolism_carbon_ES.png")
      )
      make_metabolism_heatmap("carbon utilization (Woodcroft)",
        meta_cols = 5,
        subtitle_tool = "DRAM2", mag_subset = mags_no_species,
        annotation_version = ann_v,
        cache_file = DRAM_CACHE_FILE_10,
        legend_br = FALSE,
        title_es = "Carbono Woodcroft – 10 MAGs sin especie",
        outfile = file.path(OUTPUT_DIR, "Fig7c_metabolism_carbon_woodcroft_ES.png")
      )
      make_metals_heatmap(
        title_es = "Genes de reducción de metales – 10 MAGs sin especie (top 20)",
        subtitle_tool = "DRAM2", mag_subset = mags_no_species,
        annotation_version = ann_v,
        cache_file = DRAM_CACHE_FILE_10,
        legend_br = FALSE,
        outfile = file.path(OUTPUT_DIR, "Fig7d_metabolism_metals_ES.png")
      )
      make_metabolism_heatmap("CAMPER",
        meta_cols = 5, group_col = "pathway",
        subtitle_tool = "DRAM2", mag_subset = mags_no_species,
        annotation_version = ann_v,
        cache_file = DRAM_CACHE_FILE_10,
        legend_br = FALSE,
        title_es = "CAMPER – 10 MAGs sin especie",
        outfile = file.path(OUTPUT_DIR, "Fig7e_metabolism_CAMPER_ES.png")
      )
      make_metabolism_heatmap("Organic Nitrogen",
        meta_cols = 5, group_col = "pathway",
        subtitle_tool = "DRAM2", mag_subset = mags_no_species,
        annotation_version = ann_v,
        cache_file = DRAM_CACHE_FILE_10,
        legend_br = FALSE,
        title_es = "Nitrógeno orgánico – 10 MAGs sin especie",
        outfile = file.path(OUTPUT_DIR, "Fig7f_metabolism_organic_N_ES.png")
      )
      make_metabolism_heatmap("Transporters",
        meta_cols = 5, group_col = "pathway",
        subtitle_tool = "DRAM2", mag_subset = mags_no_species,
        annotation_version = ann_v,
        cache_file = DRAM_CACHE_FILE_10,
        legend_br = FALSE,
        title_es = "Transportadores – 10 MAGs sin especie",
        outfile = file.path(OUTPUT_DIR, "Fig7g_metabolism_transporters_ES.png")
      )
      make_metabolism_heatmap("MISC",
        meta_cols = 5, make_group_label = .misc_grp_fn,
        subtitle_tool = "DRAM2", mag_subset = mags_no_species,
        annotation_version = ann_v,
        cache_file = DRAM_CACHE_FILE_10,
        legend_br = FALSE,
        title_es = "Funciones misceláneas – 10 MAGs sin especie",
        outfile = file.path(OUTPUT_DIR, "Fig7h_metabolism_misc_ES.png")
      )
    }

    # ---- Summarized genomes ----
    cat("  Generando heatmaps de genes por categoría...\n")
    .sub_tool_16 <- "DRAM2"
    if (file.exists(SUMMARIZED_GENOMES_FILE)) {
      make_summarized_heatmap(expr(category == "central carbon"),
        title_es = "Carbono central – 10 MAGs sin especie",
        subtitle_tool = .sub_tool_16, mag_subset = mags_no_species,
        annotation_version = ann_v,
        legend_br = FALSE,
        outfile = file.path(OUTPUT_DIR, "Fig8a_genes_central_carbon_ES.png")
      )
      make_summarized_heatmap(expr(category == "Electron transport Chain"),
        title_es = "Cadena de transporte electrónico – 10 MAGs sin especie",
        subtitle_tool = .sub_tool_16, mag_subset = mags_no_species,
        annotation_version = ann_v,
        legend_br = FALSE,
        outfile = file.path(OUTPUT_DIR, "Fig8b_genes_ETC_ES.png")
      )
      make_summarized_heatmap(expr(category == "Nitrogen"),
        title_es = "Ciclo del nitrógeno – 10 MAGs sin especie",
        subtitle_tool = .sub_tool_16, mag_subset = mags_no_species,
        annotation_version = ann_v,
        legend_br = FALSE,
        outfile = file.path(OUTPUT_DIR, "Fig8c_genes_nitrogen_ES.png")
      )
      make_summarized_heatmap(expr(category == "Sulfur"),
        title_es = "Ciclo del azufre – 10 MAGs sin especie",
        subtitle_tool = .sub_tool_16, mag_subset = mags_no_species,
        annotation_version = ann_v,
        legend_br = FALSE,
        outfile = file.path(OUTPUT_DIR, "Fig8d_genes_sulfur_ES.png")
      )
      make_summarized_heatmap(expr(category %in% c("C1-methane", "C1")),
        title_es = "Metabolismo C1/metano – 10 MAGs sin especie",
        subtitle_tool = .sub_tool_16, mag_subset = mags_no_species,
        annotation_version = ann_v,
        legend_br = FALSE,
        outfile = file.path(OUTPUT_DIR, "Fig8e_genes_C1_methane_ES.png")
      )
      make_summarized_heatmap(expr(category == "Hydrogenases"),
        title_es = "Hidrogenasas – 10 MAGs sin especie",
        subtitle_tool = .sub_tool_16, mag_subset = mags_no_species,
        annotation_version = ann_v,
        legend_br = FALSE,
        outfile = file.path(OUTPUT_DIR, "Fig8f_genes_hydrogenases_ES.png")
      )
      make_summarized_heatmap(expr(category == "CAZY"),
        title_es = "CAZymes – 10 MAGs sin especie",
        subtitle_tool = .sub_tool_16, mag_subset = mags_no_species,
        annotation_version = ann_v,
        legend_br = FALSE,
        max_genes = 60, outfile = file.path(OUTPUT_DIR, "Fig8g_genes_CAZymes_ES.png")
      )
      make_summarized_heatmap(expr(category == "Metal Reduction"),
        title_es = "Reducción de metales – 10 MAGs sin especie",
        subtitle_tool = .sub_tool_16, mag_subset = mags_no_species,
        annotation_version = ann_v,
        legend_br = FALSE,
        outfile = file.path(OUTPUT_DIR, "Fig8h_genes_metal_reduction_ES.png")
      )
      make_summarized_heatmap(expr(category == "pyruvate metabolism"),
        title_es = "Metabolismo del piruvato – 10 MAGs sin especie",
        subtitle_tool = .sub_tool_16, mag_subset = mags_no_species,
        annotation_version = ann_v,
        legend_br = FALSE,
        outfile = file.path(OUTPUT_DIR, "Fig8i_genes_pyruvate_ES.png")
      )
      make_summarized_heatmap(expr(category == "SCFA  and alcohol conversions"),
        title_es = "AGCC y fermentación – 10 MAGs sin especie",
        subtitle_tool = .sub_tool_16, mag_subset = mags_no_species,
        annotation_version = ann_v,
        legend_br = FALSE,
        outfile = file.path(OUTPUT_DIR, "Fig8j_genes_SCFA_fermentation_ES.png")
      )
      make_summarized_heatmap(expr(category == "Photosynthesis"),
        title_es = "Fotosíntesis – 10 MAGs sin especie",
        subtitle_tool = .sub_tool_16, mag_subset = mags_no_species,
        annotation_version = ann_v,
        legend_br = FALSE,
        outfile = file.path(OUTPUT_DIR, "Fig8k_genes_photosynthesis_ES.png")
      )
      make_summarized_heatmap(expr(category == "Information systems"),
        title_es = "Sistemas de información – 10 MAGs sin especie",
        subtitle_tool = .sub_tool_16, mag_subset = mags_no_species,
        annotation_version = ann_v,
        legend_br = FALSE,
        max_genes = 50, outfile = file.path(OUTPUT_DIR, "Fig8l_genes_info_systems_ES.png")
      )
    }

    cat("  [OK] Versión", ann_v, "completada en:", OUTPUT_DIR, "\n")
  } # fin loop 5 versiones de anotación

  # ---- Restaurar globales de profundidad ----
  depth_levels_std <- depth_levels_std_orig
  depth_colors_std <- depth_colors_std_orig
  depth_colors <- depth_colors_orig

  # ---- Restaurar OUTPUT_DIR ----
  OUTPUT_DIR <- OUTPUT_DIR_ORIG

  cat("\n[OK] Sección 16 completada. Resultados en:", file.path(OUTPUT_DIR, "10_MAGs"), "\n\n")
  cat("  (Tabla de calidad/taxonomía generada por tax_analysis.R)\n\n")
} # fin if length(mags_no_species) > 0

cat(rep("=", 70), "\n", sep = "")
cat("[OK] Script completo finalizado!\n\n")
cat("Fin:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
