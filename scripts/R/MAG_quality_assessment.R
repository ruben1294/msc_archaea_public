#!/usr/bin/env Rscript
# Quality analysis of the MAGs from their CheckM2 reports, following the MIMAG
# standard (Bowers et al. 2017, Nature Biotechnology), with the GUNC chimerism
# screen of pipeline step 02 alongside.
#
# Usage, from the repository root once pipeline steps 01 and 02 have run:
#
#   Rscript scripts/R/MAG_quality_assessment.R
#
# Every input has a default under results/ and data/, and the options listed by
# --help override them. Figures 1, 2, 3A and 4 compare the MAGs analysed here
# with a second MAG set from the same samples, labelled "Lizt" in the figures.
# Its CheckM2 report is not produced by this pipeline, so those figures appear
# only when it is passed with --file1.

#.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

# === Package loading ===
# Install pacman if needed and load packages
# pacman is not used: it was unavailable on the cluster where this first ran.
#if (!require("pacman")) {
#  install.packages("pacman")
#}

#required_packages <- c("tidyverse", "patchwork", "scales", "ggrepel",
#                       "viridis", "RColorBrewer", "openxlsx", "ggthemes")

# Use pacman to load and install packages if needed
#pacman::p_load(char = required_packages, install = TRUE, update = FALSE)

# Runs locally rather than on the cluster; nothing here needs a scheduler.
#
# The subplots carry no titles. Every figure below ends up inside a panel that
# has its own title, which a per-subplot title would only repeat.
#
# From the RStudio console, source() runs it with the defaults, since
# commandArgs() then carries no options.

# === Configure library paths ===
user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "") {
  user_lib <- file.path(Sys.getenv("HOME"), "R", "x86_64-pc-linux-gnu-library",
                        paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = "."))
}

# Create user library if it doesn't exist
if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}

# Set library paths with user library first
.libPaths(c(user_lib, .libPaths()))

# === Package Management Functions ===
#' Check if a package is installed
#' @param pkg Package name
#' @return TRUE if available, FALSE otherwise
check_package <- function(pkg) {
  suppressWarnings(suppressMessages(
    requireNamespace(pkg, quietly = TRUE)
  ))
}

#' Install a package (attempts only once)
#' @param pkg Package name
#' @return TRUE if successfully installed, FALSE otherwise
install_package_once <- function(pkg) {
  cat("Attempting to install package:", pkg, "...\n")

  tryCatch({
    # Try to install in user library
    install.packages(pkg, lib = user_lib, repos = "https://cloud.r-project.org/",
                     dependencies = TRUE, quiet = TRUE)

    # Verify installation succeeded
    if (check_package(pkg)) {
      cat("[OK] Successfully installed:", pkg, "\n")
      return(TRUE)
    } else {
      cat("[FAIL] Installation completed but package not available:", pkg, "\n")
      return(FALSE)
    }
  }, error = function(e) {
    cat("[FAIL] Failed to install:", pkg, "\n")
    cat("  Error:", conditionMessage(e), "\n")
    return(FALSE)
  }, warning = function(w) {
    cat("[WARN] Warning during installation of", pkg, ":", conditionMessage(w), "\n")
    # Check if package is available despite warning
    return(check_package(pkg))
  })
}

#' Load a package
#' @param pkg Package name
#' @return TRUE if successfully loaded, FALSE otherwise
load_package <- function(pkg) {
  tryCatch({
    suppressWarnings(suppressMessages(
      library(pkg, character.only = TRUE, quietly = TRUE)
    ))
    return(TRUE)
  }, error = function(e) {
    return(FALSE)
  })
}

# === Check, Install (Once), and Load Required Packages ===
cat("\n=== Checking and loading required packages ===\n")

required_packages <- c(
  "tidyverse",   # Data manipulation and visualization
  "patchwork",   # Combining plots
  "scales",      # Scale functions for visualization
  "ggrepel",     # Better text positioning in plots
  "viridis",     # Color scales
  "RColorBrewer", # Color palettes
  "openxlsx",    # Excel file writing
  "ggthemes",    # Additional ggplot themes
  "ggtext",      # Markdown/HTML text rendering (element_markdown)
  "showtext",    # Custom font support
  "sysfonts",     # System font detection
  "ggbeeswarm",
  "ggpubr",      # significance brackets and asterisks
  "rstatix"      # the Kruskal-Wallis table
)

failed_packages <- c()
installation_attempted <- FALSE

for (pkg in required_packages) {
  # Check if package is already available
  if (check_package(pkg)) {
    cat("[OK] Found:", pkg, "\n")
  } else {
    # Package not available, try to install (ONCE)
    cat("[WARN] Missing:", pkg, "\n")
    installation_attempted <- TRUE

    # Attempt installation once
    success <- install_package_once(pkg)

    if (!success) {
      failed_packages <- c(failed_packages, pkg)
    }
  }
}

# If any installations failed, exit with instructions
if (length(failed_packages) > 0) {
  cat("\n")
  cat("========================================\n")
  cat("ERROR: Failed to install required packages\n")
  cat("========================================\n\n")

  cat("Could not install the following packages:\n")
  for (pkg in failed_packages) {
    cat("  -", pkg, "\n")
  }

  cat("\n")
  cat("This is likely due to insufficient permissions or missing system dependencies.\n\n")

  cat("To fix this, install packages in your conda environment:\n\n")
  cat("conda activate R-env\n")

  # Convert R package names to conda package names
  conda_packages <- sapply(failed_packages, function(x) {
    paste0("r-", tolower(x))
  })

  cat("conda install -c conda-forge -c bioconda", paste(conda_packages, collapse = " "), "\n\n")

  cat("========================================\n")
  quit(status = 1)
}

# Now load all packages
cat("\n=== Loading packages ===\n")
for (pkg in required_packages) {
  if (load_package(pkg)) {
    cat("[OK] Loaded:", pkg, "\n")
  } else {
    cat("[FAIL] Failed to load:", pkg, "\n")
    cat("\nERROR: Package", pkg, "is installed but cannot be loaded.\n")
    cat("This might indicate a corrupted installation or dependency issue.\n")
    cat("Try reinstalling the package in your conda environment.\n\n")
    quit(status = 1)
  }
}

cat("\n[OK] All packages loaded successfully\n\n")

# === Font Configuration ===
cat("=== Configuring fonts ===\n")

# Enable showtext for rendering
showtext_auto()
showtext_opts(dpi = 300)
# Try to detect and use Aptos or similar fonts
font_family <- "sans"  # Default fallback

# Check if system fonts are available
available_fonts <- sysfonts::font_files()

# Priority list of fonts (from most to least preferred)
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
  # Check in available system fonts
  matching_fonts <- available_fonts[grepl(font_name, available_fonts$family, ignore.case = TRUE), ]

  if (nrow(matching_fonts) > 0) {
    tryCatch({
      # Build full path: sysfonts returns path (dir) + file (name) separately
      make_full_path <- function(row_idx) {
        f <- matching_fonts$file[row_idx]
        if ("path" %in% colnames(matching_fonts) && !grepl("^/", f))
          f <- file.path(matching_fonts$path[row_idx], f)
        f
      }

      # Identify regular / bold / italic / bolditalic rows via the 'style' column
      if ("style" %in% colnames(matching_fonts)) {
        st       <- tolower(trimws(matching_fonts$style))
        reg_rows <- which(st == "regular")
        bld_rows <- which(st == "bold")
        ita_rows <- which(st == "italic")
        bia_rows <- which(grepl("bold.*italic|italic.*bold", st))
        if (length(reg_rows) == 0) reg_rows <- 1L   # fallback: first file
      } else {
        reg_rows <- 1L
        bld_rows <- integer(0); ita_rows <- integer(0); bia_rows <- integer(0)
      }

      args <- list(family = "custom_font", regular = make_full_path(reg_rows[1]))
      if (length(bld_rows) > 0) args$bold      <- make_full_path(bld_rows[1])
      if (length(ita_rows) > 0) args$italic     <- make_full_path(ita_rows[1])
      if (length(bia_rows) > 0) args$bolditalic <- make_full_path(bia_rows[1])

      do.call(sysfonts::font_add, args)
      font_family <- "custom_font"
      cat("[OK] Using font:", font_name,
          if (length(bld_rows) > 0) "(bold variant loaded)\n" else "(WARNING: no bold variant found)\n")
      font_found <- TRUE
      break
    }, error = function(e) {
      # Continue to next font
    })
  }
}

# If no system font found, try Google Fonts in priority order
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

if (!font_found) {
  cat("[WARN] Using default system font (sans-serif)\n")
}

cat("[OK] Font configuration complete\n")
cat("\n")

# === CLI Parser ===
print_help <- function() {
  cat(
    "Usage: Rscript scripts/R/MAG_quality_assessment.R [OPTIONS]\n\n",
    "Every option has a default, so a run after pipeline steps 01 and 02 needs\n",
    "none of them.\n\n",
    "Options:\n",
    "  --file2 FILE            CheckM2 report of the MAGs analysed here\n",
    "                          (default: results/01_mag_quality/quality_report.tsv)\n",
    "  --file1 FILE            CheckM2 report of a second MAG set to compare against;\n",
    "                          enables Figures 1, 2, 3A and 4 (default: none)\n",
    "  --output-dir DIR        Output directory (default: results/01_mag_quality)\n",
    "  --gunc-progenomes FILE  GUNC ProGenomes3 summary from step 02 (Figs 8-12)\n",
    "  --gunc-gtdb FILE        GUNC GTDB r214 summary from step 02 (Figs 8-12)\n",
    "  --metadata FILE         MAG metadata with a DisplayName column\n",
    "                          (default: data/MAGs_metadata.tsv)\n",
    "  --help                  Show this help message and exit\n\n",
    "Description:\n",
    "  Quality analysis of MAGs from CheckM2 reports, following the MIMAG standard.\n\n",
    sep = ""
  )
}

# Starts from the defaults and overrides whatever the command line names.
parse_arguments <- function(args, defaults) {
  parsed_args <- defaults
  option_names <- c(
    "--file1" = "file1", "--file2" = "file2", "--output-dir" = "output_dir",
    "--gunc-progenomes" = "gunc_progenomes", "--gunc-gtdb" = "gunc_gtdb",
    "--metadata" = "metadata"
  )

  i <- 1
  while (i <= length(args)) {
    arg <- args[i]
    if (arg == "--help") {
      print_help()
      quit(status = 0)
    }
    if (!arg %in% names(option_names)) {
      stop("Unknown option ", arg, ". See --help.", call. = FALSE)
    }
    if (i + 1 > length(args)) {
      stop(arg, " requires an argument", call. = FALSE)
    }
    parsed_args[[option_names[[arg]]]] <- args[i + 1]
    i <- i + 2
  }

  if (!file.exists(parsed_args$file2)) {
    stop("CheckM2 report not found: ", parsed_args$file2,
         "\nRun pipeline step 01 first, or pass the report with --file2.",
         call. = FALSE)
  }
  if (!is.null(parsed_args$file1) && !file.exists(parsed_args$file1)) {
    stop("Comparison report not found: ", parsed_args$file1, call. = FALSE)
  }

  parsed_args
}

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

# Defaults for every input, so that a plain run from the terminal and one
# through source() in RStudio both read the pipeline output. The comparison set
# of --file1 has none; see the header.
cli_defaults <- list(
  file1           = NULL,
  file2           = file.path(PROJECT_DIR, "results", "01_mag_quality", "quality_report.tsv"),
  output_dir      = file.path(PROJECT_DIR, "results", "01_mag_quality"),
  gunc_progenomes = file.path(PROJECT_DIR, "results", "02_chimerism", "gunc_summary_progenomes3.tsv"),
  gunc_gtdb       = file.path(PROJECT_DIR, "results", "02_chimerism", "gunc_summary_gtdb214.tsv"),
  metadata        = file.path(PROJECT_DIR, "data", "MAGs_metadata.tsv")
)
cli_args <- parse_arguments(commandArgs(trailingOnly = TRUE), cli_defaults)

script_start <- Sys.time()
cat(sprintf("=== Inicio del script: %s ===\n\n", format(script_start, "%Y-%m-%d %H:%M:%S")))

# === Configuration ===
cat("=== Configuration ===\n")
cat("Input file 1:", if (is.null(cli_args$file1)) "none, comparison figures skipped" else cli_args$file1, "\n")
cat("Input file 2:", cli_args$file2, "\n")
cat("Output directory:", cli_args$output_dir, "\n")
if (!is.null(cli_args$gunc_progenomes))
  cat("GUNC ProGenomes3 file:", cli_args$gunc_progenomes, "\n")
if (!is.null(cli_args$gunc_gtdb))
  cat("GUNC GTDB r214 file:", cli_args$gunc_gtdb, "\n")
if (!is.null(cli_args$metadata))
  cat("Metadata file:", cli_args$metadata, "\n")
cat("\n")

input_1 <- cli_args$file1
input_2 <- cli_args$file2
output_dir <- cli_args$output_dir

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# === MIMAG quality thresholds ===
MIMAG_THRESHOLDS <- list(
  high_quality = list(completeness = 90, contamination = 5),
  medium_quality = list(completeness = 50, contamination = 10),
  near_complete = list(completeness = 90, contamination = 5)
)

# MAG names as utils/05_rename_mags.sh leaves them, dots in place of
# underscores. GUNC and the metadata's NameWithoutExtension use that form, and
# so does the CheckM2 report once renamed; one written before keeps underscores.
mag_key <- function(x) gsub("_", ".", x, fixed = TRUE)

# === Data loading and processing ===
#' Load and classify MAG quality data
#' @param file_path Path to TSV file with quality metrics
#' @param dataset_name Name to assign to this dataset
#' @return Tibble with added MIMAG quality classifications and metadata
load_and_classify <- function(file_path, dataset_name) {
  # Load file with error handling
  if (!file.exists(file_path)) {
    stop("Error: Input file not found: ", file_path)
  }

  tryCatch({
    df <- read_tsv(file_path, show_col_types = FALSE)
  }, error = function(e) {
    stop("Error reading file ", file_path, ": ", conditionMessage(e))
  })

  df <- df %>%
    mutate(
      Dataset = dataset_name,
      # Extract metadata from MAG names with validation
      Station = ifelse(is.na(str_extract(Name, "^[A-Z0-9]+")),
                       "Unknown",
                       str_extract(Name, "^[A-Z0-9]+")),
      Depth = ifelse(is.na(str_extract(Name, "(?<=[_.])[A-Z0-9]+(?=[_.]megahit)")),
                     "Unknown",
                     str_extract(Name, "(?<=[_.])[A-Z0-9]+(?=[_.]megahit)")),
      Binner = case_when(
        str_detect(Name, "metabat") ~ "MetaBAT2",
        str_detect(Name, "maxbin") ~ "MaxBin2",
        str_detect(Name, "concoct") ~ "CONCOCT",
        TRUE ~ "Unknown"
      ),
      # MIMAG quality classification
      MIMAG_Quality = case_when(
        Completeness >= 90 & Contamination <= 5 ~ "Calidad alta",
        Completeness >= 50 & Contamination < 10 ~ "Calidad media",
        Completeness < 50 & Contamination < 10 ~ "Calidad baja",
        Contamination >= 10 ~ "Contaminado"
      ),
      MIMAG_Quality = factor(MIMAG_Quality, levels = c(
        "Calidad alta",
        "Calidad media",
        "Calidad baja",
        "Contaminado"
      )),
      # Quality metric (Completeness - 5*Contamination)
      Quality_Score = Completeness - 5 * Contamination,
      # Genome size categories
      Size_Category = cut(Genome_Size / 1e6,
                          breaks = c(0, 1, 2, 3, 5, Inf),
                          labels = c("<1 Mb", "1-2 Mb", "2-3 Mb", "3-5 Mb", ">5 Mb"))
    )
  return(df)
}

# The MAGs analysed here, keyed so that the joins with GUNC and the metadata
# work whether or not the report has been through the rename.
df_2 <- load_and_classify(input_2, "Valentín") %>% mutate(Name = mag_key(Name))

# The comparison set when given; otherwise an empty table with the same
# columns, so the summaries below need no special case.
HAS_COMPARISON <- !is.null(input_1)
df_1 <- if (HAS_COMPARISON) load_and_classify(input_1, "Lizt") else df_2[0, ]
df_combined <- bind_rows(df_1, df_2)

# === Summary statistics ===
summary_stats <- df_combined %>%
  group_by(Dataset, MIMAG_Quality) %>%
  summarise(
    n = n(),
    mean_completeness = mean(Completeness),
    mean_contamination = mean(Contamination),
    mean_genome_size = mean(Genome_Size) / 1e6,
    mean_n50 = mean(Contig_N50),
    mean_gc = mean(GC_Content) * 100,
    .groups = "drop"
  )

# Export summary table
openxlsx::write.xlsx(
  list(
    "Summary_by_Quality" = summary_stats,
    "Full_Data_1" = df_1,
    "Full_Data_2" = df_2
  ),
  file = file.path(output_dir, "MAG_quality_summary.xlsx")
)

# === Theme publication ===
# 2026/02/01: Removed X-axis grid.
theme_publication <- function(base_size = 16) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.major.y = element_line(color = "gray", linetype = "dotted"),
      text         = element_text(family = font_family),
      axis.text    = element_text(color = "black", size = base_size - 1),
      axis.title   = element_text(face = "bold",   size = base_size + 2),
      axis.title.x = element_text(face = "bold",   size = base_size + 2, color = "black"),
      axis.title.y = element_text(face = "bold",   size = base_size + 2, color = "black"),
      legend.position = "right",
      legend.background = element_rect(fill = "white", color = NA),
      legend.title = element_text(face = "bold", size = base_size),
      legend.text  = element_text(size = base_size - 2),
      strip.background = element_rect(fill = "gray90", color = "black"),
      strip.text   = element_text(face = "bold", size = base_size),
      plot.title   = element_text(face = "bold", hjust = 0.5, size = base_size + 2)
    )
}

# === Taxon italic helper (same CDC logic as tax_analysis.R) ===
# CDC scientific nomenclature (https://wwwnc.cdc.gov/eid/page/scientific-nomenclature):
# only genus and species (and infraspecific ranks) are italicised; ranks above
# genus (family, order, class, phylum, domain) remain plain text.
ITALIC_LEVELS <- c("Genus", "Species", "Subspecies", "Strain")

should_italicize_level <- function(level) {
  !is.null(level) && length(level) == 1 && !is.na(level) && level %in% ITALIC_LEVELS
}

# When level qualifies for italics, words without digits and not all-uppercase
# are wrapped in *...* for ggtext; codes (UBA60, sp002507125, MGIII) stay plain
# so binomials render as "*Thalassarchaeum* sp002507125".
taxon_italic <- function(x, level = "Species") {
  if (!should_italicize_level(level)) {
    return(ifelse(is.na(x), x, as.character(x)))
  }
  sapply(x, function(nm) {
    if (is.na(nm) || nm == "Sin asignar") return(nm)
    words    <- strsplit(nm, " ")[[1]]
    is_plain <- grepl("[0-9]", words) | grepl("^[A-Z]+$", words)
    r        <- rle(is_plain)
    ends     <- cumsum(r$lengths)
    starts   <- c(1L, ends[-length(ends)] + 1L)
    chunks   <- mapply(function(s, e, plain) {
      grp <- paste(words[s:e], collapse = " ")
      if (plain) grp else paste0("*", grp, "*")
    }, starts, ends, r$values, SIMPLIFY = TRUE)
    paste(chunks, collapse = " ")
  }, USE.NAMES = FALSE)
}

# === Reference-line helpers ===
# strip_reflines(): removes geom_vline / geom_hline / geom_abline layers
# (threshold reference lines) from a ggplot or patchwork.
# Used by save_dual() to produce the _sin_lineas version of each figure.
# theme_publication() turns the vertical grid lines off
# (panel.grid.major.x = element_blank()), so neither version has them.
strip_reflines <- function(p) {
  is_refline <- function(l)
    inherits(l$geom, c("GeomVline", "GeomHline", "GeomAbline"))
  if (inherits(p, "patchwork")) {
    for (i in seq_len(length(p))) p[[i]] <- strip_reflines(p[[i]])
    return(p)
  }
  if (inherits(p, "gg")) {
    p$layers <- Filter(Negate(is_refline), p$layers)
  }
  p
}

# save_dual(): single-output wrapper kept for call-site compatibility.
save_dual <- function(plot, outfile, width, height, dpi = 300) {
  ggsave(outfile, plot, width = width, height = height, dpi = dpi)
}

# Figures 1 and 2 compare the two MAG sets, so they need --file1.
if (HAS_COMPARISON) {

  # === Figure 1: Completeness vs Contamination (MIMAG plot) ===
  cat("=== Generating Figure 1: MIMAG quality plot ===\n")
  # Alpha refers to the opacity of a geom. Values of alpha range from 0 to 1,
  # with lower values corresponding to more transparent colors.

  p1 <- ggplot(df_combined, aes(x = Contamination, y = Completeness,
                                fill = Dataset,    # fill by Dataset
                                color = Dataset)) + # outline by Dataset
    # MIMAG quality regions
    annotate("rect", xmin = 0, xmax = 5, ymin = 90, ymax = 100,
             fill = "#2ecc71", alpha = 0.1) +
    annotate("rect", xmin = 0, xmax = 10, ymin = 50, ymax = 90,
             fill = "#f39c12", alpha = 0.1) +
    annotate("rect", xmin = 5, xmax = 10, ymin = 90, ymax = 100,
             fill = "#f39c12", alpha = 0.1) +
    annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 50,
             fill = "#95a5a6", alpha = 0.1) +
    annotate("rect", xmin = 10, xmax = Inf, ymin = 0, ymax = 100,
             fill = "#e74c3c", alpha = 0.1) +
    # Reference lines
    geom_vline(xintercept = c(5, 10), linetype = "dashed", color = "gray30", linewidth = 0.5) +
    geom_hline(yintercept = c(50, 90), linetype = "dashed", color = "gray30", linewidth = 0.5) +
    # Data points, with a heavy outline
    geom_point(size = 2.7, alpha = 0.4,
               shape = 21,      # circle with an outline
               stroke = 0.8) +  # outline width
    scale_fill_manual(values = c(     # fill colours
      "Lizt" = "#FF3030",
      "Valentín" = "#006300")) +
    scale_color_manual(values = c(    # outline colours (same as the fill)
      "Lizt" = "#FF3030",
      "Valentín" = "#006300")) +
    labs(
      title = "Evaluación de calidad de los conjuntos de MAGs",
      x = "Contaminación (%)",
      y = "Completitud (%)",
      fill = "Conjunto de MAGs",    # legend title
      color = "Conjunto de MAGs"
    ) +
    guides(color = "none") +  # hide the colour legend; only fill is shown
    scale_x_continuous(limits = c(0, max(df_combined$Contamination) + 1),
                       breaks = seq(0, 15, 5)) +
    scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10))

  p1_publication <- p1 + theme_publication()
  save_dual(p1_publication, file.path(output_dir, "Fig1_pub_MIMAG_quality.png"), width = 10, height = 6)

  #p1_clean <- p1 + theme_clean() +
  #  theme(legend.background = element_rect(fill = "white", color = NA),
  #        legend.title = element_text(face = "bold"))
  #ggsave(file.path(output_dir, "Fig1_clean_MIMAG_quality.png"),
  #       p1_clean, width = 10, height = 6, dpi = 300)

  cat("[OK] Saved Figure 1\n")

  #ggsave(file.path(output_dir, "Fig1_MIMAG_quality_plot.pdf"),
   #      p1, width = 10, height = 6)

  # === Figure 2: Quality distribution by dataset ===
  cat("=== Generating Figure 2: Quality distributions ===\n")

  dataset_colors <- c("Lizt" = "#f04b82", "Valentín" = "#2899c9")

  # Build the subplots WITHOUT theme() at the end
  p2a <- ggplot(df_combined, aes(x = MIMAG_Quality, fill = Dataset)) +
    geom_bar(position = "dodge", color = "black", linewidth = 0.3) +
    geom_text(stat = "count", aes(label = after_stat(count)),
              position = position_dodge(width = 0.9), vjust = -0.5, size = 3.5) +
    scale_fill_manual(values = dataset_colors) +  # custom colours
    labs(
      title = "Distribución de calidad de MAGs por conjunto de MAGs",
      x = "Categoría de calidad MIMAG",
      y = "Número de MAGs",
      fill = "Conjunto de MAGs"
    )

  p2b <- ggplot(df_combined, aes(x = Dataset, y = Completeness, fill = Dataset)) +
    geom_violin(alpha = 0.4, color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = 21, linewidth = 0.5) +
    geom_hline(yintercept = c(50, 90), linetype = "dashed", color = "gray30", linewidth = 0.5) +
    scale_fill_manual(values = dataset_colors) +  # custom colours
    labs(
      title = "Distribución de completitud",
      x = NULL,
      y = "Completitud (%)",
      fill = "Conjunto de MAGs"
    ) +
    guides(fill = "none")

  p2c <- ggplot(df_combined, aes(x = Dataset, y = Contamination, fill = Dataset)) +
    geom_violin(alpha = 0.4, color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.2, alpha = 0.8, outlier.shape = 21, linewidth = 0.5) +
    geom_hline(yintercept = c(5, 10), linetype = "dashed", color = "gray30", linewidth = 0.5) +
    scale_fill_manual(values = dataset_colors) +  # custom colours
    labs(
      title = "Distribución de contaminación",
      x = NULL,
      y = "Contaminación (%)",
      fill = "Conjunto de MAGs"
    ) +
    guides(fill = "none")

  # Combine the plots
  p2 <- p2a / (p2b | p2c) + plot_layout(heights = c(1.2, 1))

  # Apply the themes with &
  p2_publication <- p2 & theme_publication() & theme(
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

  save_dual(p2_publication, file.path(output_dir, "Fig2_pub_quality_distributions.png"), width = 12, height = 8)

  #p2_clean <- p2 & theme_clean() & theme(
  #  axis.text.x = element_text(angle = 25, hjust = 1),
  #  legend.background = element_rect(fill = "white", color = NA),
  #  legend.title = element_text(face = "bold")
  #)

  #ggsave(file.path(output_dir, "Fig2_clean_quality_distributions.png"),
  #       p2_clean, width = 12, height = 8, dpi = 300)

  cat("[OK] Saved Figure 2\n")

}

# === Figure 3: Genome assembly metrics (SEPARATE PER DATASET) ===
cat("=== Generating Figure 3: Assembly metrics (by dataset) ===\n")

# Colour palette for the MIMAG quality bands
mimag_colors <- c(
  "Calidad alta" = "#27ae60",
  "Calidad media" = "#f2c911",
  "Calidad baja" = "#95a5a6",
  "Contaminado" = "#e74c3c"
)

# Figure 3A shows the comparison set alone, so it needs --file1 too.
if (HAS_COMPARISON) {

  # ---- FIGURE 3A: the Lizt MAG set ----

  # Lizt dataset
  df_lizt <- df_combined %>% filter(Dataset == "Lizt")

  # Subplot 3A-1: genome size against N50 (Lizt)
  p3a_lizt <- ggplot(df_lizt, aes(x = Genome_Size / 1e6, y = Contig_N50 / 1000,
                                  color = MIMAG_Quality)) +
    geom_point(size = 2.5, alpha = 0.7) +
    scale_x_log10(labels = comma) +
    scale_y_log10(labels = comma) +
    scale_color_manual(values = mimag_colors) +
    labs(
      title = "Tamaño del genoma vs N50 de contigs",
      x = "Tamaño del genoma (Mb)",
      y = "Contig N50 (kb)",
      color = "Categoría de calidad"
    )

  # Subplot 3A-2: GC vs Coding Density (Lizt)
  p3b_lizt <- ggplot(df_lizt, aes(x = GC_Content * 100, y = Coding_Density * 100,
                                  color = MIMAG_Quality)) +
    geom_point(size = 2.5, alpha = 0.7) +
    scale_color_manual(values = mimag_colors) +
    labs(
      title = "Contenido GC vs Densidad de codificación",
      x = "Contenido GC (%)",
      y = "Densidad de codificación (%)",
      color = "Categoría de calidad"
    )

  # Subplot 3A-3: quality against size (Lizt)
  p3c_lizt <- ggplot(df_lizt, aes(x = Size_Category, fill = MIMAG_Quality)) +
    geom_bar(position = "stack", color = "black", linewidth = 0.3) +
    scale_fill_manual(values = mimag_colors) +
    labs(
      title = "Calidad por tamaño de genoma",
      x = "Categoría de tamaño de genoma",
      y = "Número de MAGs",
      fill = "Categoría de calidad"
    )

  # Assemble the Lizt panel
  p3_lizt <- (p3a_lizt | p3b_lizt) / p3c_lizt +
    plot_layout(heights = c(1, 1)) +
    plot_annotation(
      title = "Métricas de ensamble de genomas - Conjunto Dra. Lizt",
      theme = theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
    )

  # Apply the themes
  p3_lizt_publication <- p3_lizt & theme_publication() & theme(
    axis.text.x = element_text(angle = 25, hjust = 1)
  )

  save_dual(p3_lizt_publication, file.path(output_dir, "Fig3A_pub_assembly_metrics_Lizt.png"), width = 12, height = 10)

  #p3_lizt_clean <- p3_lizt & theme_clean() & theme(
  #  axis.text.x = element_text(angle = 25, hjust = 1),
  #  legend.background = element_rect(fill = "white", color = NA),
  #  legend.title = element_text(face = "bold")
  #)

  #ggsave(file.path(output_dir, "Fig3A_clean_assembly_metrics_Lizt.png"),
  #       p3_lizt_clean, width = 12, height = 10, dpi = 300)

  cat("[OK] Saved Figure 3A\n")

}

# ---- FIGURE 3B: the Valentín MAG set ----

# Valentin dataset, the 54 MAGs this thesis is about
df_valentin <- df_combined %>% filter(Dataset == "Valentín")

# Subplot 3B-1: genome size against N50
p3a_valentin <- ggplot(df_valentin, aes(x = Genome_Size / 1e6, y = Contig_N50 / 1000,
                                        color = MIMAG_Quality)) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_x_log10(labels = comma) +
  scale_y_log10(labels = comma) +
  scale_color_manual(values = mimag_colors) +
  labs(
    title = "Tamaño del genoma vs N50 de contigs",
    x = "Tamaño del genoma (Mb)",
    y = "Contig N50 (kb)",
    color = "Categoría de calidad"
  )

# Subplot 3B-2: GC against coding density
p3b_valentin <- ggplot(df_valentin, aes(x = GC_Content * 100, y = Coding_Density * 100,
                                        color = MIMAG_Quality)) +
  geom_point(size = 2.5, alpha = 0.7) +
  scale_color_manual(values = mimag_colors) +
  labs(
    title = "Contenido GC vs Densidad de codificación",
    x = "Contenido GC (%)",
    y = "Densidad de codificación (%)",
    color = "Categoría de calidad"
  )

# Subplot 3B-3: quality against size
p3c_valentin <- ggplot(df_valentin, aes(x = Size_Category, fill = MIMAG_Quality)) +
  geom_bar(position = "stack", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = mimag_colors) +
  labs(
    title = "Calidad por tamaño de genoma",
    x = "Categoría de tamaño de genoma",
    y = "Número de MAGs",
    fill = "Categoría de calidad"
  )

# Assemble the panel
p3_valentin <- (p3a_valentin | p3b_valentin) / p3c_valentin +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(
    title = "Métricas de ensamble de genomas - Conjunto Dr. Valentín",
    theme = theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
  )

# Apply the themes
p3_valentin_publication <- p3_valentin & theme_publication() & theme(
  axis.text.x = element_text(angle = 25, hjust = 1)
)

save_dual(p3_valentin_publication, file.path(output_dir, "Fig3B_pub_assembly_metrics_Valentin.png"), width = 12, height = 10)

#p3_valentin_clean <- p3_valentin & theme_clean() & theme(
#  axis.text.x = element_text(angle = 25, hjust = 1),
#  legend.background = element_rect(fill = "white", color = NA),
#  legend.title = element_text(face = "bold")
#)

#ggsave(file.path(output_dir, "Fig3B_clean_assembly_metrics_Valentin.png"),
#       p3_valentin_clean, width = 12, height = 10, dpi = 300)

cat("[OK] Saved Figure 3B\n")

# === Figure 4: Binner comparison (Lizt dataset only) ===
cat("=== Generating Figure 4: Binner comparison ===\n")

if (nrow(df_1) > 0 && any(df_1$Binner != "Unknown")) {
  binner_stats <- df_1 %>%
    filter(Binner != "Unknown") %>%
    group_by(Binner, MIMAG_Quality) %>%
    summarise(n = n(), .groups = "drop")

  if (nrow(binner_stats) > 0) {
    # Subplots WITHOUT their own theme()
    p4a <- ggplot(binner_stats, aes(x = Binner, y = n, fill = MIMAG_Quality)) +
      geom_col(position = "stack", color = "black", linewidth = 0.3) +
      geom_text(aes(label = n), position = position_stack(vjust = 0.5),
                color = "black", fontface = "bold", size = 5) +
      scale_fill_manual(values = c(
        "Calidad alta" = "#27ae60",
        "Calidad media" = "#f2c911",
        "Calidad baja" = "#95a5a6",
        "Contaminado" = "#e74c3c"
      )) +
      labs(
        title = "Calidad de MAGs por herramienta de binning (conjunto Dra. Lizt)",
        x = "Herramienta de binning",
        y = "Número de MAGs",
        fill = "Categoría de calidad"
      ) +
      theme(legend.background = element_rect(fill = "white", color = NA))

    p4b <- df_1 %>%
      filter(Binner != "Unknown") %>%
      ggplot(aes(x = Binner, y = Quality_Score, fill = Binner)) +
      geom_violin(alpha = 0.4, color = "black", linewidth = 0.5) +
      geom_boxplot(width = 0.3, alpha = 0.8, outlier.shape = 21, linewidth = 0.5) +
      scale_fill_brewer(palette = "Pastel1") +
      labs(
        title = "Calidad por herramienta de binning",
        x = "Herramienta de binning",
        y = "Puntuación de calidad (Qscore)"
      ) +
      guides(fill = "none")  # hide the legend in this subplot

    p4 <- p4a / p4b

    p4_publication <- p4 & theme_publication()

    save_dual(p4_publication, file.path(output_dir, "Fig4_pub_binner_comparison.png"), width = 10, height = 8)

    #p4_clean <- p4 & theme_clean() & theme(
    #  legend.background = element_rect(fill = "white", color = NA),
    #  legend.title = element_text(face = "bold")
    #)

    #ggsave(file.path(output_dir, "Fig4_clean_binner_comparison.png"),
    #       p4_clean, width = 10, height = 8, dpi = 300)

    cat("[OK] Saved Figure 4\n")
  }
} else {
  cat("[WARN] No binning data available. Skipping Figure 4.\n")
}

# === Figure 5: Station and depth analysis ===
station_summary <- df_combined %>%
  filter(!is.na(Station), !is.na(Depth)) %>%
  group_by(Dataset, Station, Depth, MIMAG_Quality) %>%
  summarise(n = n(), .groups = "drop")

station_summary$Depth <- factor(station_summary$Depth,
                               levels = c("OMZ", "1000", "BTM"))

if (nrow(station_summary) > 0) {
  p5 <- ggplot(station_summary, aes(x = Station, y = n, fill = MIMAG_Quality)) +
    geom_col(position = "stack", color = "black", linewidth = 0.2) +
    facet_grid(Depth ~ Dataset, scales = "free_x", space = "free_x") +
    scale_fill_manual(values = c(
      "Calidad alta" = "#27ae60",
      "Calidad media" = "#f2c911",
      "Calidad baja" = "#95a5a6",
      "Contaminado" = "#e74c3c"
    )) +
  labs(
    title = "Calidad de MAGs por estación y profundidad",
    x = "Estación de muestreo",
    y = "Número de MAGs",
    fill = "Categoría de calidad"
  ) +
#  theme_publication() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.background = element_rect(fill = "white", color = NA))
  p5_publication <- p5 + theme_publication()

  #p5_clean <- p5 & theme_clean() & theme(
  #  axis.text.x = element_text(angle = 45, hjust = 1),
  #  legend.background = element_rect(fill = "white", color = NA),
  #  legend.title = element_text(face = "bold")
  #)

# Adjust figure height based on number of unique depths to improve readability
n_depths <- length(unique(station_summary$Depth))
fig_height <- max(8, 2 + n_depths * 3)

save_dual(p5_publication, file.path(output_dir, "Fig5_pub_station_depth_analysis.png"), width = 12, height = fig_height)

#ggsave(file.path(output_dir, "Fig5_clean_station_depth_analysis.png"),
#       p5_clean, width = 12, height = fig_height, dpi = 300)

#   ggsave(file.path(output_dir, "Fig5_station_depth_analysis.pdf"),
#          p5, width = 12, height = fig_height)
} else {
  warning("No station/depth data available. Skipping Figure 5.")
}

# ---- Figure 6: completeness against contamination, this study only ----
cat("=== Generating Figure 6: MIMAG quality plot (Valentín only) ===\n")

p6 <- ggplot(df_valentin, aes(x = Contamination, y = Completeness,
                              color = MIMAG_Quality)) +
  # MIMAG quality regions
  annotate("rect", xmin = 0, xmax = 5, ymin = 90, ymax = 100,
           fill = "#2ecc71", alpha = 0.1) +
  annotate("rect", xmin = 0, xmax = 10, ymin = 50, ymax = 90,
           fill = "#f39c12", alpha = 0.1) +
  annotate("rect", xmin = 5, xmax = 10, ymin = 90, ymax = 100,
           fill = "#f39c12", alpha = 0.1) +
  annotate("rect", xmin = 0, xmax = 10, ymin = 0, ymax = 50,
           fill = "#95a5a6", alpha = 0.1) +
  annotate("rect", xmin = 10, xmax = Inf, ymin = 0, ymax = 100,
           fill = "#e74c3c", alpha = 0.1) +
  # Reference lines
  geom_vline(xintercept = c(5, 10), linetype = "dashed",
             color = "gray30", linewidth = 0.5) +
  geom_hline(yintercept = c(50, 90), linetype = "dashed",
             color = "gray30", linewidth = 0.5) +
  # Data points
  geom_point(size = 3, alpha = 0.7) +
  scale_color_manual(values = mimag_colors) +
  labs(
    x = "Contaminación (%)",
    y = "Completitud (%)",
    color = "Categoría de calidad"
  ) +
  scale_x_continuous(limits = c(0, max(df_valentin$Contamination) + 1),
                     breaks = seq(0, 20, 5)) +
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, 10))

# Apply themes
p6_publication <- p6 + theme_publication()

save_dual(p6_publication, file.path(output_dir, "Fig6_pub_MIMAG_quality_Valentin.png"), width = 10, height = 6)

#p6_clean <- p6 + theme_clean() +
#  theme(legend.background = element_rect(fill = "white", color = NA),
#        legend.title = element_text(face = "bold"))

#ggsave(file.path(output_dir, "Fig6_clean_MIMAG_quality_Valentin.png"),
#       p6_clean, width = 10, height = 6, dpi = 300)

cat("[OK] Saved Figure 6\n")

# 2026/04/02:
# ---- STATISTICAL TESTS (Kruskal-Wallis) ON THE QUALITY DISTRIBUTIONS ----
cat("\n=== Realizando pruebas estadísticas (Kruskal-Wallis) ===\n")

# The test runs on the 54 high and medium quality MAGs, less ZMOR (n = 2) and
# C22 (n = 2), leaving n = 50. Two observations cannot support a group-wise
# comparison, and leaving them in would let a single MAG drive a result.
df_stats <- df_valentin %>%
  filter(MIMAG_Quality %in% c("Calidad alta", "Calidad media")) %>%
  mutate(Depth = recode(Depth, "OMZ" = "ZMOR", "1000" = "1000 m", "BTM" = "Fondo")) %>%
  filter(Depth != "ZMOR", Station != "C22")

# Long format, so every metric goes through the same test in one pass
df_long_stats <- df_stats %>%
  select(Name, Station, Depth, Completeness, Contamination, Contig_N50, Genome_Size) %>%
  pivot_longer(cols = c(Completeness, Contamination, Contig_N50, Genome_Size),
               names_to = "Metric", values_to = "Value")

# Kruskal-Wallis by depth
kw_depth <- df_long_stats %>%
  group_by(Metric) %>%
  kruskal_test(Value ~ Depth) %>%
  mutate(Variable_Agrupacion = "Profundidad")

# Kruskal-Wallis by station
kw_station <- df_long_stats %>%
  group_by(Metric) %>%
  kruskal_test(Value ~ Station) %>%
  mutate(Variable_Agrupacion = "Estacion")

# Combine and write the table
kw_results <- bind_rows(kw_depth, kw_station) %>%
  select(Variable_Agrupacion, Metric, n, statistic, df, p) %>%
  arrange(Variable_Agrupacion, Metric)

openxlsx::write.xlsx(
  list(
    "KW_Results"  = kw_results,
    "MAGs_usados" = df_stats %>%
      select(Name, Station, Depth, MIMAG_Quality, Completeness, Contamination,
             Contig_N50, Genome_Size) %>%
      arrange(Depth, Station, Name)
  ),
  file = file.path(output_dir, "Table_Kruskal_Wallis_Results.xlsx")
)
cat("[OK] Tabla de resultados estadísticos guardada (con hoja MAGs_usados).\n")

# --- COMPARISON RULES FOR THE BOXPLOTS ---
my_depth_comparisons <- list(c("1000 m", "Fondo")) # leave out ZMOR
estaciones_validas <- setdiff(unique(as.character(df_stats$Station)), "C22") # leave out C22
my_station_comparisons <- combn(estaciones_validas, 2, simplify = FALSE)

# Returns only the pairs at p < 0.05 by Wilcoxon rank-sum. The unit conversions
# applied for plotting, dividing by 1e6 or multiplying by 100, cannot change the
# outcome: Wilcoxon works on ranks, and any monotonic transformation leaves the
# ordering intact.
sig_pairs <- function(data, x_var, y_var, all_comparisons) {
  Filter(function(pair) {
    g1 <- data[[y_var]][as.character(data[[x_var]]) == pair[1]]
    g2 <- data[[y_var]][as.character(data[[x_var]]) == pair[2]]
    if (length(g1) < 2 || length(g2) < 2) return(FALSE)
    tryCatch(wilcox.test(g1, g2)$p.value < 0.05, error = function(e) FALSE)
  }, all_comparisons)
}

# Returns the stat_compare_means layer when there are significant pairs, and an
# empty list() otherwise, which ggplot2 quietly ignores when added with +.
add_sig_brackets <- function(comparisons, style = c("depth", "station")) {
  style <- match.arg(style)
  if (length(comparisons) == 0) return(list())
  if (style == "depth") {
    list(stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                            label = "p.signif", vjust = 0.5))
  } else {
    list(stat_compare_means(comparisons = comparisons, method = "wilcox.test",
                            label = "p.signif", step.increase = 0.08))
  }
}

# ---- Figure 7: distributions by depth and station, this study only ----
cat("=== Generating Figure 7: Distributions by Depth and Station (Valentín) ===\n")

# Rows with both Station and Depth recorded
df_val_filtered <- df_valentin %>%
  filter(!is.na(Station), !is.na(Depth),
         Station != "Unknown", Depth != "Unknown",
         MIMAG_Quality %in% c("Calidad alta", "Calidad media")) %>%
  mutate(
    Depth = recode(Depth, "OMZ" = "ZMOR", "1000" = "1000 m", "BTM" = "Fondo"),
    Depth = factor(Depth, levels = c("ZMOR", "1000 m", "Fondo"))
  )
cat("  Fig 7 dataset: ", nrow(df_val_filtered), "MAGs (Calidad alta + media)\n")

if (nrow(df_val_filtered) > 0) {

  # Work out the significant comparisons before plotting, so that
  # stat_compare_means only ever receives pairs at p < 0.05. Otherwise it draws
  # a bracket for every pair and most of them carry no asterisk.
  depth_cmp_completeness  <- sig_pairs(df_val_filtered, "Depth", "Completeness",   my_depth_comparisons)
  depth_cmp_contamination <- sig_pairs(df_val_filtered, "Depth", "Contamination",  my_depth_comparisons)
  depth_cmp_quality       <- sig_pairs(df_val_filtered, "Depth", "Quality_Score",  my_depth_comparisons)
  depth_cmp_genome        <- sig_pairs(df_val_filtered, "Depth", "Genome_Size",    my_depth_comparisons)
  depth_cmp_coding        <- sig_pairs(df_val_filtered, "Depth", "Coding_Density", my_depth_comparisons)
  depth_cmp_n50           <- sig_pairs(df_val_filtered, "Depth", "Contig_N50",     my_depth_comparisons)
  depth_cmp_gc            <- sig_pairs(df_val_filtered, "Depth", "GC_Content",     my_depth_comparisons)

  stn_cmp_completeness    <- sig_pairs(df_val_filtered, "Station", "Completeness",   my_station_comparisons)
  stn_cmp_contamination   <- sig_pairs(df_val_filtered, "Station", "Contamination",  my_station_comparisons)
  stn_cmp_quality         <- sig_pairs(df_val_filtered, "Station", "Quality_Score",  my_station_comparisons)
  stn_cmp_genome          <- sig_pairs(df_val_filtered, "Station", "Genome_Size",    my_station_comparisons)
  stn_cmp_coding          <- sig_pairs(df_val_filtered, "Station", "Coding_Density", my_station_comparisons)
  stn_cmp_n50             <- sig_pairs(df_val_filtered, "Station", "Contig_N50",     my_station_comparisons)
  stn_cmp_gc              <- sig_pairs(df_val_filtered, "Station", "GC_Content",     my_station_comparisons)

  # Define color palette for each metric
  color_completeness <- "#2C7BB6"
  color_contamination <- "#D7191C"
  color_quality <- "#006300"
  color_genome_size <- "#de7f2c"

  # --- Row 1: Completeness ---
  # By Depth (left)
  p7a <- ggplot(df_val_filtered, aes(x = Depth, y = Completeness)) +
#    geom_violin(alpha = 0.3, fill = color_completeness, color = "black", linewidth = 0.5) +
    #geom_dotplot(binaxis = "y", stackdir = "center", alpha = 0.9, fill = color_completeness, color = "black", dotsize = 0.5) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2, fill = color_completeness, color = "black", shape = 21, stroke = 0.3) +
#    geom_quasirandom(alpha = 0.7, size = 2.5,
#                     fill = color_completeness, color = "black",
#                     shape = 21, stroke = 0.3) +
    geom_boxplot(width = 0.2, alpha = 0.7, fill = color_completeness,
#                 outlier.shape = 21, linewidth = 0.5) +
                 linewidth = 0.5, outlier.shape = NA) +
    geom_hline(yintercept = 90, linetype = "dashed", color = "gray30", linewidth = 0.5) +
    scale_x_discrete(position = "top") +
    scale_y_continuous(limits = c(50, 100), breaks = seq(50, 100, 10),
                       expand = expansion(mult = c(0, 0.03)),
                       labels = scales::label_number(accuracy = 1)) +
    add_sig_brackets(depth_cmp_completeness, "depth") +
    labs(
      x = "Profundidad",
      y = "Completitud (%)"
    ) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 0),
          axis.title.x.top = element_text(face = "bold", color = "black"))

  # By Station (right)
  p7b <- ggplot(df_val_filtered, aes(x = Station, y = Completeness)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2, fill = color_completeness, color = "black", shape = 21, stroke = 0.3) +
    #geom_violin(alpha = 0.3, fill = color_completeness, color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.2, alpha = 0.7, fill = color_completeness,
                 #outlier.shape = 21, linewidth = 0.5) +
                 linewidth = 0.5, outlier.shape = NA) +
    geom_hline(yintercept = 90, linetype = "dashed", color = "gray30", linewidth = 0.5) +
    scale_x_discrete(position = "top") +
    scale_y_continuous(limits = c(50, 100), breaks = seq(50, 100, 10),
                       expand = expansion(mult = c(0, 0.03)),
                       labels = scales::label_number(accuracy = 1)) +
    add_sig_brackets(stn_cmp_completeness, "station") +
    labs(
      x = "Estación de muestreo",
      y = NULL
    ) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 0),
          axis.title.x.top = element_text(face = "bold", color = "black"))

  # --- Row 2: Contamination ---
  # By Depth (left)
  p7c <- ggplot(df_val_filtered, aes(x = Depth, y = Contamination)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2, fill = color_contamination, color = "black", shape = 21, stroke = 0.3) +
    #geom_violin(alpha = 0.3, fill = color_contamination, color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.2, alpha = 0.7, fill = color_contamination,
#                 outlier.shape = 21, linewidth = 0.5) +
                     linewidth = 0.5, outlier.shape = NA) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "gray30", linewidth = 0.5) +
    scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2),
                       expand = expansion(mult = c(0.03, 0)),
                       labels = scales::label_number(accuracy = 1))  +
    add_sig_brackets(depth_cmp_contamination, "depth") +
    labs(
      x = NULL,
      y = "Contaminación (%)"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # By Station (right)
  p7d <- ggplot(df_val_filtered, aes(x = Station, y = Contamination)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2, fill = color_contamination, color = "black", shape = 21, stroke = 0.3) +
#    geom_violin(alpha = 0.3, fill = color_contamination, color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.2, alpha = 0.7, fill = color_contamination,
                 #outlier.shape = 21,
                 linewidth = 0.5, outlier.shape = NA) +
    geom_hline(yintercept = 5, linetype = "dashed", color = "gray30", linewidth = 0.5) +
    scale_y_continuous(limits = c(0, NA), breaks = seq(0, 10, 2),
                       expand = expansion(mult = c(0.03, 0.25)),
                       labels = scales::label_number(accuracy = 1)) +
    add_sig_brackets(stn_cmp_contamination, "station") +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # --- Row 3: Quality Score ---
  # By Depth (left)
  p7e <- ggplot(df_val_filtered, aes(x = Depth, y = Quality_Score)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2, fill = color_quality, color = "black", shape = 21, stroke = 0.3) +
#    geom_violin(alpha = 0.3, fill = color_quality, color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.2, alpha = 0.7, fill = color_quality,
                 #outlier.shape = 21,
                 linewidth = 0.5, outlier.shape = NA) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "gray30", linewidth = 0.5)  +
    add_sig_brackets(depth_cmp_quality, "depth") +
    labs(
      x = NULL,
      y = "Puntuación de calidad"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # By Station (right)
  p7f <- ggplot(df_val_filtered, aes(x = Station, y = Quality_Score)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2, fill = color_quality, color = "black", shape = 21, stroke = 0.3) +
#    geom_violin(alpha = 0.3, fill = color_quality, color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.2, alpha = 0.7, fill = color_quality,
                 #outlier.shape = 21,
                 linewidth = 0.5, outlier.shape = NA) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "gray30", linewidth = 0.5) +
    add_sig_brackets(stn_cmp_quality, "station") +
    labs(
      x = NULL,
      y = NULL
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # --- Row 4: Genome Size  ---
  # By Depth (left)
  p7g <- ggplot(df_val_filtered, aes(x = Depth, y = Genome_Size / 1e6)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2, fill = color_genome_size, color = "black", shape = 21, stroke = 0.3) +
#    geom_violin(alpha = 0.3, fill = color_genome_size, color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.2, alpha = 0.7, fill = color_genome_size,
                 #outlier.shape = 21,
                 linewidth = 0.5, outlier.shape = NA) +
    scale_y_continuous(labels = scales::label_number(drop0trailing = TRUE))  +
    add_sig_brackets(depth_cmp_genome, "depth") +
    labs(
      x = "Profundidad",
      y = "Tamaño del genoma (Mb)"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # By Station (right)
  p7h <- ggplot(df_val_filtered, aes(x = Station, y = Genome_Size / 1e6)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2, fill = color_genome_size, color = "black", shape = 21, stroke = 0.3) +
#    geom_violin(alpha = 0.3, fill = color_genome_size, color = "black", linewidth = 0.5) +
    geom_boxplot(width = 0.2, alpha = 0.7, fill = color_genome_size,
                 #outlier.shape = 21,
                 linewidth = 0.5, outlier.shape = NA) +
    add_sig_brackets(stn_cmp_genome, "station") +
    scale_y_continuous(labels = scales::label_number(drop0trailing = TRUE)) +
    labs(
      x = "Estación de muestreo",
      y = "Tamaño del genoma (Mb)"
    ) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  # Fig 7b: additional assembly quality metrics (4 rows × 2 columns).
  # Top row uses scale_x_discrete(position = "top") for labels; lower rows
  # hide x-axis entirely.  Same wrap_plots + [[i]] pattern as Fig 7.
  color_coding_density <- "#8e44ad"
  color_n50            <- "#16a085"
  color_gc             <- "#2c3e50"

  # --- Row 1: Genome size (Mb), reusing p7g and p7h with the x axis on top ---
  p7g_top <- p7g + scale_x_discrete(position = "top") +
    labs(x = "Profundidad") +
    theme(axis.text.x.top  = element_text(angle = 45, hjust = 0,
                                          margin = margin(b = 8)),
          axis.title.x.top = element_text(face = "bold", color = "black"),
          axis.text.x      = element_blank(),
          axis.ticks.x     = element_blank())

  p7h_top <- p7h + scale_x_discrete(position = "top") +
    labs(x = "Estación de muestreo", y = NULL) +
    theme(axis.text.x.top  = element_text(angle = 45, hjust = 0,
                                          margin = margin(b = 8)),
          axis.title.x.top = element_text(face = "bold", color = "black"),
          axis.text.x      = element_blank(),
          axis.ticks.x     = element_blank())

  # --- Row 2: Coding density (%) ---
  p7b_cd_depth <- ggplot(df_val_filtered,
                          aes(x = Depth, y = Coding_Density * 100)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2,
                fill = color_coding_density, color = "black", shape = 21, stroke = 0.3) +
    geom_boxplot(width = 0.2, alpha = 0.7, fill = color_coding_density,
                 linewidth = 0.5, outlier.shape = NA) +
    scale_y_continuous(labels = scales::label_number(drop0trailing = TRUE)) +
    add_sig_brackets(depth_cmp_coding, "depth") +
    labs(x = NULL, y = "Densidad de codificación (%)") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  p7b_cd_station <- ggplot(df_val_filtered,
                            aes(x = Station, y = Coding_Density * 100)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2,
                fill = color_coding_density, color = "black", shape = 21, stroke = 0.3) +
    geom_boxplot(width = 0.2, alpha = 0.7, fill = color_coding_density,
                 linewidth = 0.5, outlier.shape = NA) +
    scale_y_continuous(labels = scales::label_number(drop0trailing = TRUE)) +
    add_sig_brackets(stn_cmp_coding, "station") +
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  # --- Row 3: Contig N50 (kb) ---
  p7b_n50_depth <- ggplot(df_val_filtered,
                           aes(x = Depth, y = Contig_N50 / 1000)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2,
                fill = color_n50, color = "black", shape = 21, stroke = 0.3) +
    geom_boxplot(width = 0.2, alpha = 0.7, fill = color_n50,
                 linewidth = 0.5, outlier.shape = NA) +
    scale_y_continuous(labels = scales::label_number(drop0trailing = TRUE)) +
    add_sig_brackets(depth_cmp_n50, "depth") +
    labs(x = NULL, y = expression(bolditalic("Contig") * bold(" N50 (kb)"))) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  p7b_n50_station <- ggplot(df_val_filtered,
                             aes(x = Station, y = Contig_N50 / 1000)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2,
                fill = color_n50, color = "black", shape = 21, stroke = 0.3) +
    geom_boxplot(width = 0.2, alpha = 0.7, fill = color_n50,
                 linewidth = 0.5, outlier.shape = NA) +
    scale_y_continuous(labels = scales::label_number(drop0trailing = TRUE)) +
    add_sig_brackets(stn_cmp_n50, "station") +
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  # --- Row 4: GC content (%) ---
  p7b_gc_depth <- ggplot(df_val_filtered,
                          aes(x = Depth, y = GC_Content * 100)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2,
                fill = color_gc, color = "black", shape = 21, stroke = 0.3) +
    geom_boxplot(width = 0.2, alpha = 0.7, fill = color_gc,
                 linewidth = 0.5, outlier.shape = NA) +
    scale_y_continuous(labels = scales::label_number(drop0trailing = TRUE)) +
    add_sig_brackets(depth_cmp_gc, "depth") +
    labs(x = NULL, y = "Contenido GC (%)") +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  p7b_gc_station <- ggplot(df_val_filtered,
                            aes(x = Station, y = GC_Content * 100)) +
    geom_jitter(width = 0.2, alpha = 0.6, size = 2,
                fill = color_gc, color = "black", shape = 21, stroke = 0.3) +
    geom_boxplot(width = 0.2, alpha = 0.7, fill = color_gc,
                 linewidth = 0.5, outlier.shape = NA) +
    scale_y_continuous(labels = scales::label_number(drop0trailing = TRUE)) +
    add_sig_brackets(stn_cmp_gc, "station") +
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

  p7b_all <- wrap_plots(p7g_top,       p7h_top,
                        p7b_cd_depth,  p7b_cd_station,
                        p7b_n50_depth, p7b_n50_station,
                        p7b_gc_depth,  p7b_gc_station,
                        ncol = 2) +
    plot_layout(guides = "collect")

  p7b_genome_pub <- p7b_all & theme_publication()
  # Top row (panels 1-2): x-axis at top. The theme is already set in the panel definition,
  # but re-apply after theme_publication() to guarantee it wins.
  for (ti in 1:2) {
    p7b_genome_pub[[ti]] <- p7b_genome_pub[[ti]] + theme(
      axis.text.x.top  = element_text(angle = 45, hjust = 0, color = "black",
                                      margin = margin(b = 8)),
      axis.title.x.top = element_text(face = "bold", color = "black")
    )
  }
  # Rows 2-4 (panels 3-8): ensure x-axis is blank
  for (li in 3:8) {
    p7b_genome_pub[[li]] <- p7b_genome_pub[[li]] + theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    )
  }

  save_dual(p7b_genome_pub, file.path(output_dir, "Fig7b_pub_assembly_metrics.png"),
            width = 14, height = 16)
  save_dual(
    (p7b_genome_pub + plot_annotation(tag_levels = "a")) &
      theme(plot.tag = element_text(face = "bold", size = 20, color = "black", family = "sans")),
    file.path(output_dir, "Fig7b_pub_assembly_metrics_labeled.png"),
    width = 14, height = 16)
  cat("[OK] Saved Figure 7b (assembly metrics)\n")

  # Fig 7d: assembly metrics by taxonomic order
  #
  # The same metrics as Fig 7b (Genome_Size, Coding_Density, Contig_N50,
  # GC_Content), grouped by Order from FullTaxonomy rather than by depth or
  # station. Each subplot carries a Kruskal-Wallis result, restricted to orders
  # with at least ORDER_MIN_N members, 3 by default.
  if (!is.null(cli_args$metadata) && file.exists(cli_args$metadata)) {
    cat("=== Generating Figure 7d: Assembly metrics by Order ===\n")

    ORDER_MIN_N <- 3L
    .order_meta_d <- read_tsv(cli_args$metadata, show_col_types = FALSE) %>%
      select(Name = NameWithoutExtension, FullTaxonomy) %>%
      mutate(Order = str_extract(FullTaxonomy, "(?<=o__)[^;]+"),
             Class = str_extract(FullTaxonomy, "(?<=c__)[^;]+")) %>%
      filter(!is.na(Order), Order != "")

    order_levels_taxa <- c("Poseidoniales", "Nitrososphaerales",
                           "MGIII", "Undinarchaeales")

    # An order name that is a code, first word all caps or carrying
    # digits, tells a reader nothing on its own, so the parent class goes in
    # front of it as "c. <Class>". Same rule as fmt_taxon() in tax_analysis.R.
    # Order names stay upright: CDC style italicises genus and below only.
    .build_order_labels <- function(order_class_df) {
      is_acronym <- function(nm) {
        fw <- strsplit(nm, " ")[[1]][1]
        grepl("[0-9]", fw) || grepl("^[A-Z]+$", fw)
      }
      ocd <- order_class_df %>% distinct(Order, Class)
      setNames(
        vapply(seq_len(nrow(ocd)), function(i) {
          o <- ocd$Order[i]; cl <- ocd$Class[i]
          if (is_acronym(o) && !is.na(cl) && nchar(cl) > 0 && cl != o) {
            paste0(o, " (c. ", cl, ")")
          } else {
            o
          }
        }, character(1)),
        ocd$Order
      )
    }
    order_text_labels <- .build_order_labels(.order_meta_d)

    df_val_order <- df_val_filtered %>%
      left_join(.order_meta_d %>% select(Name, Order), by = "Name") %>%
      filter(!is.na(Order), Order %in% order_levels_taxa) %>%
      mutate(Order = factor(Order, levels = order_levels_taxa))

    n_per_order <- df_val_order %>% count(Order)
    cat("  Fig 7d dataset:", nrow(df_val_order), "MAGs;",
        "n por orden: ",
        paste(n_per_order$Order, "=", n_per_order$n, collapse = "; "), "\n")

    # Orders large enough to compare, n >= min_n. sig_pairs() then keeps only
    # the pairs that reach p < 0.05.
    valid_orders <- as.character(n_per_order$Order[n_per_order$n >= ORDER_MIN_N])
    my_order_comparisons <- if (length(valid_orders) >= 2) {
      combn(valid_orders, 2, simplify = FALSE)
    } else {
      list()
    }

    base_order_layers <- function(fill_color) {
      list(
        geom_jitter(width = 0.18, alpha = 0.6, size = 2,
                    fill = fill_color, color = "black",
                    shape = 21, stroke = 0.3),
        geom_boxplot(width = 0.25, alpha = 0.7, fill = fill_color,
                     linewidth = 0.5, outlier.shape = NA),
        scale_x_discrete(labels = order_text_labels)
      )
    }

    # Significant pairs per metric, Wilcoxon at p < 0.05
    ord_cmp_genome <- sig_pairs(df_val_order, "Order", "Genome_Size",    my_order_comparisons)
    ord_cmp_coding <- sig_pairs(df_val_order, "Order", "Coding_Density", my_order_comparisons)
    ord_cmp_n50    <- sig_pairs(df_val_order, "Order", "Contig_N50",     my_order_comparisons)
    ord_cmp_gc     <- sig_pairs(df_val_order, "Order", "GC_Content",     my_order_comparisons)

    p7d_gen <- ggplot(df_val_order, aes(x = Order, y = Genome_Size / 1e6)) +
      base_order_layers(color_genome_size) +
      scale_y_continuous(labels = scales::label_number(drop0trailing = TRUE)) +
      add_sig_brackets(ord_cmp_genome, "station") +
      labs(x = NULL, y = "Tamaño del genoma (Mb)")

    p7d_cd <- ggplot(df_val_order, aes(x = Order, y = Coding_Density * 100)) +
      base_order_layers(color_coding_density) +
      scale_y_continuous(labels = scales::label_number(drop0trailing = TRUE)) +
      add_sig_brackets(ord_cmp_coding, "station") +
      labs(x = NULL, y = "Densidad de codificación (%)")

    p7d_n50 <- ggplot(df_val_order, aes(x = Order, y = Contig_N50 / 1000)) +
      base_order_layers(color_n50) +
      scale_y_continuous(labels = scales::label_number(drop0trailing = TRUE)) +
      add_sig_brackets(ord_cmp_n50, "station") +
      labs(x = NULL, y = expression(bolditalic("Contig") * bold(" N50 (kb)")))

    p7d_gc <- ggplot(df_val_order, aes(x = Order, y = GC_Content * 100)) +
      base_order_layers(color_gc) +
      scale_y_continuous(labels = scales::label_number(drop0trailing = TRUE)) +
      add_sig_brackets(ord_cmp_gc, "station") +
      labs(x = NULL, y = "Contenido GC (%)")

    p7d_all <- wrap_plots(p7d_gen, p7d_cd, p7d_n50, p7d_gc, ncol = 2) +
      plot_layout(guides = "collect")
    p7d_pub <- p7d_all & theme_publication() &
      theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black"))

    save_dual(p7d_pub, file.path(output_dir, "Fig7d_pub_assembly_by_order.png"),
              width = 12, height = 10)
    save_dual(
      (p7d_pub + plot_annotation(tag_levels = "a")) &
        theme(plot.tag = element_text(face = "bold", size = 20, color = "black", family = "sans")),
      file.path(output_dir, "Fig7d_pub_assembly_by_order_labeled.png"),
      width = 12, height = 10)
    cat("[OK] Saved Figure 7d (assembly metrics by Order)\n")
  } else {
    cat("[WARN] --metadata not provided; skipping Figure 7d (assembly by Order).\n")
  }

  # Fig 7: 3-row version (chimerism row added later in the GUNC block if data available).
  # wrap_plots() creates a flat patchwork where [[1]]..[[6]] address each plot
  # individually, unlike (A|B)/(C|D) which creates nested rows.
  p7 <- wrap_plots(p7a, p7b, p7c, p7d, p7e, p7f, ncol = 2) +
    plot_layout(guides = "collect") +
    plot_annotation(
      theme = theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
    )

  p7_publication <- p7 & theme_publication()
  # Panels 1-2 (top row): x-axis at top; add margin so labels don't overlap panel.
  for (top_i in 1:2) {
    p7_publication[[top_i]] <- p7_publication[[top_i]] + theme(
      axis.text.x.top  = element_text(angle = 45, hjust = 0, color = "black",
                                      margin = margin(b = 8)),
      axis.title.x.top = element_text(face = "bold", color = "black")
    )
  }
  # Panels 3-6: hide x-axis labels and ticks entirely
  for (lo_i in 3:6) {
    p7_publication[[lo_i]] <- p7_publication[[lo_i]] + theme(
      axis.text.x  = element_blank(),
      axis.ticks.x = element_blank()
    )
  }

  save_dual(p7_publication, file.path(output_dir, "Fig7_pub_distributions_depth_station_Valentin.png"), width = 14, height = 12)
  save_dual(
    (p7_publication + plot_annotation(tag_levels = "a")) &
      theme(plot.tag = element_text(face = "bold", size = 20, color = "black", family = "sans")),
    file.path(output_dir, "Fig7_pub_distributions_depth_station_Valentin_labeled.png"),
    width = 14, height = 12)

  #p7_clean <- p7 & theme_clean() & theme(
  #  axis.text.x = element_text(angle = 45, hjust = 1),
  #  legend.background = element_rect(fill = "white", color = NA),
  #  legend.title = element_text(face = "bold")
  #)

  #ggsave(file.path(output_dir, "Fig7_clean_distributions_depth_station_Valentin.png"),
  #       p7_clean, width = 14, height = 16, dpi = 300)

  cat("[OK] Saved Figure 7\n")

} else {
  cat("[WARN] No valid Station/Depth data for Valentín dataset. Skipping Figure 7.\n")
}

# GUNC chimerism, figures 8 to 12
#
# Needs --gunc-progenomes and --gunc-gtdb pointing at the gunc_summary_*.tsv
# files that pipeline/02_chimerism.sh writes.
#
# GUNC (Orakov et al. 2021, Genome Biol.) flags chimeric bins by how far the
# taxonomy of their single-copy genes spreads beyond one dominant clade. Its
# clade separation score (CSS) marks a genome chimeric above 0.45
# (pass.GUNC = False).
#
# Two databases:
#   ProGenomes3  broad prokaryotic coverage, the baseline
#   GTDB r214    finer archaeal taxonomy, preferred for these MAGs
#
# The rank where CSS peaks (maxCSS level) shows how far apart the mixed lineages
# are: phylum or class means contamination across phyla, the worst case; genus
# or species means closely related strains mixed together.
#
# Figures:
#   Fig 8   GUNC pass and fail, and the CSS distribution, per database
#   Fig 9   CSS by sampling depth and station, both databases
#   Fig 10  taxonomic rank where CSS peaks
#   Fig 11  CSS concordance, ProGenomes3 against GTDB r214
#   Fig 12  per-MAG summary: completeness, contamination and CSS

if (!is.null(cli_args$gunc_progenomes) && !is.null(cli_args$gunc_gtdb)) {

  cat("\n=== Loading GUNC chimerism data ===\n")

  # Load one GUNC summary TSV and annotate with spatial metadata
  load_gunc <- function(file_path, db_label) {
    if (!file.exists(file_path)) {
      cat("[WARN] GUNC file not found:", file_path, "\n")
      return(NULL)
    }
    df <- read_tsv(file_path, show_col_types = FALSE) %>%
      rename(CSS = clade_separation_score,
             RRS = reference_representation_score) %>%
      mutate(
        Database  = db_label,
        # Extract station (e.g. A4) and depth (e.g. BTM) from genome name
        Station   = str_extract(genome, "^[A-Z0-9]+"),
        Depth     = str_extract(genome, "(?<=\\.)[A-Z0-9]+(?=\\.megahit)"),
        Depth     = recode(Depth, "OMZ" = "ZMOR", "1000" = "1000 m", "BTM" = "Fondo"),
        Depth     = factor(Depth, levels = c("ZMOR", "1000 m", "Fondo")),
        # Normalised name for joining with CheckM2 data
        Name      = mag_key(genome),
        # read_tsv parses "True"/"False" as logical (TRUE/FALSE); convert back
        # to character "TRUE"/"FALSE" before factoring to avoid NA mismatch.
        pass.GUNC = factor(as.character(pass.GUNC), levels = c("TRUE", "FALSE")),
        taxonomic_level = factor(
          recode(as.character(taxonomic_level),
                 "kingdom" = "Reino",  "phylum"  = "Filo",
                 "class"   = "Clase",  "order"   = "Orden",
                 "family"  = "Familia","genus"   = "Género",
                 "species" = "Especie"),
          levels = c("Reino", "Filo", "Clase", "Orden", "Familia", "Género", "Especie")
        )
      )
    cat("[OK] Loaded GUNC:", db_label, "(", nrow(df), "genomes)\n")
    df
  }

  gunc_prog <- load_gunc(cli_args$gunc_progenomes, "ProGenomes3")
  gunc_gtdb <- load_gunc(cli_args$gunc_gtdb,       "GTDB r214")

  if (!is.null(gunc_prog) && !is.null(gunc_gtdb)) {

    # Combined long-format table: both databases stacked
    gunc_all <- bind_rows(gunc_prog, gunc_gtdb) %>%
      mutate(Database = factor(Database, levels = c("ProGenomes3", "GTDB r214")))

    # Shared colours and labels
    # Levels are "TRUE"/"FALSE" because read_tsv parsed the original "True"/"False"
    # strings as logical, and as.character(TRUE/FALSE) gives "TRUE"/"FALSE".
    pass_colors  <- c("TRUE" = "#27ae60", "FALSE" = "#e74c3c")
    pass_labels  <- c("TRUE" = "No quimérico", "FALSE" = "Quimérico")
    css_threshold <- 0.45

    # -------------------------------------------------------------------------
    # Load MAGs metadata for DisplayName lookup (requires --metadata argument)
    # Falls back to the raw CheckM2 name if metadata is not provided.
    # -------------------------------------------------------------------------
    if (!is.null(cli_args$metadata) && file.exists(cli_args$metadata)) {
      mag_meta <- read_tsv(cli_args$metadata, show_col_types = FALSE) %>%
        select(Name = NameWithoutExtension, DisplayName)
      gunc_all  <- gunc_all  %>% left_join(mag_meta, by = "Name") %>%
        mutate(DisplayName = if_else(is.na(DisplayName), Name, DisplayName))
      gunc_prog <- gunc_prog %>% left_join(mag_meta, by = "Name") %>%
        mutate(DisplayName = if_else(is.na(DisplayName), Name, DisplayName))
      gunc_gtdb <- gunc_gtdb %>% left_join(mag_meta, by = "Name") %>%
        mutate(DisplayName = if_else(is.na(DisplayName), Name, DisplayName))
      cat("[OK] DisplayName loaded from metadata\n")
    } else {
      gunc_all  <- gunc_all  %>% mutate(DisplayName = Name)
      gunc_prog <- gunc_prog %>% mutate(DisplayName = Name)
      gunc_gtdb <- gunc_gtdb %>% mutate(DisplayName = Name)
      if (is.null(cli_args$metadata))
        cat("[WARN] --metadata not provided; using raw MAG names in figures.\n")
    }

    # -------------------------------------------------------------------------
    # Fig 7c: quality metrics by taxonomic order, with RRS from GTDB r214.
    #
    # Fig 7 again (completeness, contamination, Quality_Score, CSS) grouped by
    # order, plus a fifth subplot for RRS. All five reference GTDB r214, and
    # each carries a Kruskal-Wallis result restricted to orders with at least
    # ORDER_MIN_N members.
    # -------------------------------------------------------------------------
    if (exists("df_val_filtered") && nrow(df_val_filtered) > 0 &&
        !is.null(cli_args$metadata) && file.exists(cli_args$metadata)) {
      cat("=== Generating Figure 7c: Quality metrics by Order (with RRS) ===\n")

      ORDER_MIN_N_C <- 3L
      .order_meta_c <- read_tsv(cli_args$metadata, show_col_types = FALSE) %>%
        select(Name = NameWithoutExtension, FullTaxonomy) %>%
        mutate(Order = str_extract(FullTaxonomy, "(?<=o__)[^;]+"),
               Class = str_extract(FullTaxonomy, "(?<=c__)[^;]+")) %>%
        filter(!is.na(Order), Order != "")

      order_levels_taxa_c <- c("Poseidoniales", "Nitrososphaerales",
                               "MGIII", "Undinarchaeales")

      # Same rules as fmt_taxon() in tax_analysis.R: the order in plain text.
      # A code-like order name gets its parent class prefixed as "c.".
      .build_order_labels_c <- function(order_class_df) {
        is_acronym <- function(nm) {
          fw <- strsplit(nm, " ")[[1]][1]
          grepl("[0-9]", fw) || grepl("^[A-Z]+$", fw)
        }
        ocd <- order_class_df %>% distinct(Order, Class)
        setNames(
          vapply(seq_len(nrow(ocd)), function(i) {
            o <- ocd$Order[i]; cl <- ocd$Class[i]
            if (is_acronym(o) && !is.na(cl) && nchar(cl) > 0 && cl != o) {
              paste0(o, " (c. ", cl, ")")
            } else {
              o
            }
          }, character(1)),
          ocd$Order
        )
      }
      order_text_labels_c <- .build_order_labels_c(.order_meta_c)

      df_val_order_c <- df_val_filtered %>%
        left_join(.order_meta_c %>% select(Name, Order), by = "Name") %>%
        left_join(gunc_gtdb %>% select(Name, CSS, RRS), by = "Name") %>%
        filter(!is.na(Order), Order %in% order_levels_taxa_c) %>%
        mutate(Order = factor(Order, levels = order_levels_taxa_c))

      cat("  Fig 7c dataset:", nrow(df_val_order_c), "MAGs;",
          sum(!is.na(df_val_order_c$RRS)), "con GUNC GTDB r214 RRS\n")

      # Valid pairs for Wilcoxon. CSS and RRS are NA for MAGs that GUNC never
      # saw, so n per order differs between panels and the pairs are recomputed
      # inside each one rather than shared.
      n_per_order_c <- df_val_order_c %>% count(Order)
      valid_orders_c <- as.character(n_per_order_c$Order[n_per_order_c$n >= ORDER_MIN_N_C])
      my_order_comparisons_c <- if (length(valid_orders_c) >= 2) {
        combn(valid_orders_c, 2, simplify = FALSE)
      } else {
        list()
      }

      base_order_layers_c <- function(fill_color) {
        list(
          geom_jitter(width = 0.18, alpha = 0.6, size = 2,
                      fill = fill_color, color = "black",
                      shape = 21, stroke = 0.3),
          geom_boxplot(width = 0.25, alpha = 0.7, fill = fill_color,
                       linewidth = 0.5, outlier.shape = NA),
          scale_x_discrete(labels = order_text_labels_c)
        )
      }

      d_css <- df_val_order_c %>% filter(!is.na(CSS))
      d_rrs <- df_val_order_c %>% filter(!is.na(RRS))

      # Significant pairs per metric, Wilcoxon at p < 0.05. For CSS and RRS the
      # valid pairs are recomputed over the non-NA rows.
      n_per_order_css <- d_css %>% count(Order)
      valid_orders_css <- as.character(
        n_per_order_css$Order[n_per_order_css$n >= ORDER_MIN_N_C])
      cmp_css_pairs <- if (length(valid_orders_css) >= 2)
        combn(valid_orders_css, 2, simplify = FALSE) else list()

      n_per_order_rrs <- d_rrs %>% count(Order)
      valid_orders_rrs <- as.character(
        n_per_order_rrs$Order[n_per_order_rrs$n >= ORDER_MIN_N_C])
      cmp_rrs_pairs <- if (length(valid_orders_rrs) >= 2)
        combn(valid_orders_rrs, 2, simplify = FALSE) else list()

      ord_cmp_completeness  <- sig_pairs(df_val_order_c, "Order", "Completeness",  my_order_comparisons_c)
      ord_cmp_contamination <- sig_pairs(df_val_order_c, "Order", "Contamination", my_order_comparisons_c)
      ord_cmp_quality       <- sig_pairs(df_val_order_c, "Order", "Quality_Score", my_order_comparisons_c)
      ord_cmp_css           <- sig_pairs(d_css,           "Order", "CSS",          cmp_css_pairs)
      ord_cmp_rrs           <- sig_pairs(d_rrs,           "Order", "RRS",          cmp_rrs_pairs)

      p7c_comp <- ggplot(df_val_order_c, aes(x = Order, y = Completeness)) +
        base_order_layers_c(color_completeness) +
        geom_hline(yintercept = 90, linetype = "dashed",
                   color = "gray30", linewidth = 0.5) +
        scale_y_continuous(limits = c(50, 100), breaks = seq(50, 100, 10),
                           labels = scales::label_number(accuracy = 1)) +
        add_sig_brackets(ord_cmp_completeness, "station") +
        labs(x = NULL, y = "Completitud (%)")

      p7c_cont <- ggplot(df_val_order_c, aes(x = Order, y = Contamination)) +
        base_order_layers_c(color_contamination) +
        geom_hline(yintercept = 5, linetype = "dashed",
                   color = "gray30", linewidth = 0.5) +
        scale_y_continuous(limits = c(0, 10), breaks = seq(0, 10, 2),
                           labels = scales::label_number(accuracy = 1)) +
        add_sig_brackets(ord_cmp_contamination, "station") +
        labs(x = NULL, y = "Contaminación (%)")

      p7c_qs <- ggplot(df_val_order_c, aes(x = Order, y = Quality_Score)) +
        base_order_layers_c(color_quality) +
        geom_hline(yintercept = 50, linetype = "dashed",
                   color = "gray30", linewidth = 0.5) +
        add_sig_brackets(ord_cmp_quality, "station") +
        labs(x = NULL, y = "Puntuación de calidad")

      p7c_css <- ggplot(d_css, aes(x = Order, y = CSS)) +
        base_order_layers_c("#9b59b6") +
        geom_hline(yintercept = css_threshold, linetype = "dashed",
                   color = "gray30", linewidth = 0.5) +
        scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
        add_sig_brackets(ord_cmp_css, "station") +
        labs(x = NULL,
             y = expression(bolditalic("Clade Separation Score") * bold(" (CSS)")))

      p7c_rrs <- ggplot(d_rrs, aes(x = Order, y = RRS)) +
        base_order_layers_c("#e67e22") +
        geom_hline(yintercept = 0.5, linetype = "dashed",
                   color = "gray30", linewidth = 0.5) +
        scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
        add_sig_brackets(ord_cmp_rrs, "station") +
        labs(x = NULL,
             y = expression(bolditalic("Reference Representation Score") * bold(" (RRS)")))

      # 5 panels: 3 cols × 2 rows; the trailing cell stays empty.
      p7c_all <- wrap_plots(p7c_comp, p7c_cont, p7c_qs, p7c_css, p7c_rrs,
                            ncol = 3) +
        plot_layout(guides = "collect")
      p7c_pub <- p7c_all & theme_publication() &
        theme(axis.text.x = element_text(angle = 30, hjust = 1, color = "black"))

      save_dual(p7c_pub, file.path(output_dir, "Fig7c_pub_quality_by_order.png"),
                width = 16, height = 10)
      save_dual(
        (p7c_pub + plot_annotation(tag_levels = "a")) &
          theme(plot.tag = element_text(face = "bold", size = 20, color = "black", family = "sans")),
        file.path(output_dir, "Fig7c_pub_quality_by_order_labeled.png"),
        width = 16, height = 10)
      cat("[OK] Saved Figure 7c (quality metrics by Order, with RRS)\n")
    } else {
      cat("[WARN] Skipping Figure 7c — needs Valentín filtered data and --metadata.\n")
    }

    # -------------------------------------------------------------------------
    # Fig 7 update: replace genome-size row with chimerism (GTDB r214) row.
    # p7a to p7f come from the block above; only regenerate them if they exist.
    # -------------------------------------------------------------------------
    if (exists("p7a")) {
      cat("=== Regenerating Figure 7 with chimerism row (GTDB r214) ===\n")

      gunc_depth_data   <- gunc_gtdb %>% filter(!is.na(Depth))
      gunc_station_data <- gunc_gtdb %>% filter(!is.na(Station))

      p7g_chim <- ggplot(gunc_depth_data,
                         aes(x = Depth, y = CSS, fill = pass.GUNC)) +
        geom_jitter(aes(color = pass.GUNC),
                    width = 0.15, alpha = 0.6, size = 2,
                    shape = 21, stroke = 0.3) +
        geom_boxplot(width = 0.2, alpha = 0.5, outlier.shape = NA,
                     linewidth = 0.5, color = "black",
                     position = "identity") +
        geom_hline(yintercept = css_threshold,
                   linetype = "dashed", color = "gray30", linewidth = 0.5) +
        scale_fill_manual(values  = pass_colors, labels = pass_labels) +
        scale_color_manual(values = pass_colors) +
        scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
        labs(x = NULL,
             y = expression(bolditalic("Clade Separation Score")*bold(" (CSS)")),
             fill = "Resultado GUNC") +
        guides(fill = "none", color = "none") +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

      p7h_chim <- ggplot(gunc_station_data,
                         aes(x = Station, y = CSS, fill = pass.GUNC)) +
        geom_jitter(aes(color = pass.GUNC),
                    width = 0.15, alpha = 0.6, size = 2,
                    shape = 21, stroke = 0.3) +
        geom_boxplot(width = 0.2, alpha = 0.5, outlier.shape = NA,
                     linewidth = 0.5, color = "black",
                     position = "identity") +
        geom_hline(yintercept = css_threshold,
                   linetype = "dashed", color = "gray30", linewidth = 0.5) +
        scale_fill_manual(values  = pass_colors, labels = pass_labels) +
        scale_color_manual(values = pass_colors) +
        scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25)) +
        labs(x = NULL, y = NULL,
             fill = "Resultado GUNC") +
        guides(fill = "none", color = "none") +
        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

      p7_gunc <- wrap_plots(p7a, p7b, p7c, p7d, p7e, p7f,
                             p7g_chim, p7h_chim, ncol = 2) +
        plot_layout(guides = "collect") +
        plot_annotation(
          theme = theme(plot.title = element_text(face = "bold", size = 18, hjust = 0.5))
        )

      p7_gunc_pub <- p7_gunc & theme_publication()
      # Panels 1-2 (top row): x-axis at top; add margin so labels clear the panel.
      for (top_i in 1:2) {
        p7_gunc_pub[[top_i]] <- p7_gunc_pub[[top_i]] + theme(
          axis.text.x.top  = element_text(angle = 45, hjust = 0, color = "black",
                                          margin = margin(b = 8)),
          axis.title.x.top = element_text(face = "bold", color = "black")
        )
      }
      # Panels 3-8: hide x-axis labels and ticks entirely
      for (lo_i in 3:8) {
        p7_gunc_pub[[lo_i]] <- p7_gunc_pub[[lo_i]] + theme(
          axis.text.x  = element_blank(),
          axis.ticks.x = element_blank()
        )
      }

      save_dual(p7_gunc_pub,
                file.path(output_dir, "Fig7_pub_distributions_depth_station_Valentin.png"),
                width = 14, height = 16)
      save_dual(
        (p7_gunc_pub + plot_annotation(tag_levels = "a")) &
          theme(plot.tag = element_text(face = "bold", size = 20, color = "black", family = "sans")),
        file.path(output_dir, "Fig7_pub_distributions_depth_station_Valentin_labeled.png"),
        width = 14, height = 16)
      cat("[OK] Saved Figure 7 (updated with chimerism row)\n")
    }

    # -------------------------------------------------------------------------
    # Fig 8: GUNC overview, pass/fail stacked bar and CSS violin
    #
    # The number of MAGs each database calls chimeric, next to the full CSS
    # distribution with the 0.45 threshold, one database beside the other.
    # -------------------------------------------------------------------------
    cat("=== Generating Figure 8: GUNC overview ===\n")

    pass_summary <- gunc_all %>%
      count(Database, pass.GUNC) %>%
      group_by(Database) %>%
      mutate(pct = n / sum(n) * 100)

    p8a <- ggplot(pass_summary,
                  aes(x = Database, y = n, fill = pass.GUNC)) +
      geom_col(color = "black", linewidth = 0.3, width = 0.28) +
      geom_text(aes(label = paste0(n, "\n(", round(pct, 1), "%)")),
                position = position_stack(vjust = 0.5),
                color = "black", fontface = "bold", size = 5) +
      scale_fill_manual(values = pass_colors, labels = pass_labels) +
      labs(x = "Base de datos de referencia", y = "Número de MAGs",
           fill = "Quimerismo") +
      guides(fill = "none")

    # position = "identity" forces both violins (TRUE / FALSE) to be centred at
    # the same x value for each Database.  Since CSS distributions are very
    # different between the two groups (non-chimeric ≈ 0, chimeric ≈ 0.5–1),
    # the violin shapes occupy different Y regions and read clearly despite the
    # overlap.  Jitter uses the same x centre without dodging.
    p8b <- ggplot(gunc_all,
                  aes(x = Database, y = CSS, fill = pass.GUNC)) +
      geom_violin(position = "identity",
                  alpha = 0.40, color = "black", linewidth = 0.4) +
      geom_jitter(aes(color = pass.GUNC),
                  width = 0.12, size = 1.8, alpha = 0.85,
                  shape = 21, stroke = 0.3) +
      geom_hline(yintercept = css_threshold,
                 linetype = "dashed", color = "black", linewidth = 0.6) +
      scale_fill_manual(values  = pass_colors, labels = pass_labels) +
      scale_color_manual(values = pass_colors) +
      scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.2)) +
      labs(x = "Base de datos de referencia",
           y = expression(bolditalic("Clade Separation Score") * bold(" (CSS)")),
           fill = "Quimerismo") +
      guides(fill = "none", color = "none")

    # Panel C: Reference Representation Score (RRS) by database.
    # Boxplot + jitter with pass.GUNC dodge so both groups are visible per DB.
    # Only fill aesthetic is mapped (no color) so no 3rd legend is generated.
    p8c <- ggplot(gunc_all,
                  aes(x = Database, y = RRS, fill = pass.GUNC)) +
      geom_boxplot(position = position_dodge(0.7), width = 0.25,
                   alpha = 0.55, outlier.shape = NA, linewidth = 0.5,
                   color = "black") +
      geom_jitter(aes(color = pass.GUNC),
                  position = position_jitterdodge(jitter.width = 0.15,
                                                   dodge.width  = 0.7),
                  size = 1.8, alpha = 0.80, shape = 21, stroke = 0.3) +
      geom_hline(yintercept = 0.5,
                 linetype = "dashed", color = "gray40", linewidth = 0.5) +
      scale_fill_manual(values  = pass_colors, labels = pass_labels) +
      scale_color_manual(values = pass_colors, guide = "none") +
      scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2)) +
      labs(x = "Base de datos de referencia",
           y = expression(bolditalic("Reference Representation Score") * bold(" (RRS)")),
           fill = "Quimerismo") +
      # No fill legend here: p8b carries the same one, and guides = "collect"
      # shows it once for the patchwork.
      guides(fill = "none")

    # ---- fig 8 panel A, grouped (dodge) variant ----
    # One bar per GUNC result (chimeric / non-chimeric) per database,
    # side by side, so the counts can be compared.
    p8a_group <- ggplot(pass_summary,
                        aes(x = Database, y = n, fill = pass.GUNC)) +
      geom_col(color = "black", linewidth = 0.3, width = 0.35,
               position = position_dodge(0.5)) +
      geom_text(aes(label = paste0(n, "\n(", round(pct, 1), "%)")),
                position = position_dodge(width = 0.5),
                vjust = -0.25, color = "black", fontface = "bold", size = 4.5) +
      scale_fill_manual(values = pass_colors, labels = pass_labels) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
      labs(x = "Base de datos de referencia", y = "Número de MAGs",
           fill = "Quimerismo") +
      guides(fill = "none")

    # Shared no-x panels for vertical variants (built once, reused by both layouts)
    p8b_nox <- p8b +
      labs(x = NULL) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

    p8c_nox <- p8c +
      labs(x = NULL) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())

    # ---- version 3: horizontal, grouped bars ----
    p8_group_pub <- (p8a_group | p8b | p8c) + plot_layout(guides = "collect")
    p8_group_pub <- p8_group_pub & theme_publication()
    save_dual(
      (p8_group_pub + plot_annotation(tag_levels = "a")) &
        theme(plot.tag = element_text(face = "bold", size = 20, color = "black", family = "sans")),
      file.path(output_dir, "Fig8_pub_agrupado_horizontal_labeled.png"),
      width = 18, height = 6)

    # ---- version 4: vertical, grouped bars ----
    p8a_group_top <- p8a_group +
      scale_x_discrete(position = "top") +
      labs(x = "Base de datos de referencia") +
      theme(axis.text.x.top  = element_text(angle = 0, hjust = 0.5, color = "black"),
            axis.title.x.top = element_text(face = "bold", color = "black"))

    p8_group_vert <- wrap_plots(p8a_group_top, p8b_nox, p8c_nox, ncol = 1) +
      plot_layout(guides = "collect")
    p8_group_vert <- p8_group_vert & theme_publication()
    p8_group_vert[[1]] <- p8_group_vert[[1]] + theme(
      axis.text.x.top  = element_text(angle = 0, hjust = 0.5, color = "black"),
      axis.title.x.top = element_text(face = "bold", color = "black")
    )
    for (vi in 2:3) {
      p8_group_vert[[vi]] <- p8_group_vert[[vi]] + theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()
      )
    }
    save_dual(
      (p8_group_vert + plot_annotation(tag_levels = "a")) &
        theme(plot.tag = element_text(face = "bold", size = 20, color = "black", family = "sans")),
      file.path(output_dir, "Fig8_pub_agrupado_vertical_labeled.png"),
      width = 8, height = 16)
    cat("[OK] Saved Figure 8 (agrupado horizontal + vertical)\n")

    # ---- fig 8 proportion variant ----
    # Panel A replaced by 100% stacked proportion bar.
    # Two layouts: horizontal (3 panels side by side) and vertical (stacked).
    p8a_prop <- ggplot(pass_summary,
                       aes(x = Database, y = pct, fill = pass.GUNC)) +
      geom_col(color = "black", linewidth = 0.3, width = 0.28) +
      geom_text(aes(label = paste0(round(pct, 0), "%")),
                position = position_stack(vjust = 0.5),
                color = "black", fontface = "bold", size = 5) +
      scale_fill_manual(values = pass_colors, labels = pass_labels) +
      scale_y_continuous(labels = function(x) paste0(x, "%")) +
      labs(x = "Base de datos de referencia", y = "MAGs (%)",
           fill = "Quimerismo") +
      guides(fill = "none")

    # Proportion horizontal
    p8_prop_pub <- (p8a_prop | p8b | p8c) + plot_layout(guides = "collect")
    p8_prop_pub <- p8_prop_pub & theme_publication()
    save_dual(
      (p8_prop_pub + plot_annotation(tag_levels = "a")) &
        theme(plot.tag = element_text(face = "bold", size = 20, color = "black", family = "sans")),
      file.path(output_dir, "Fig8_pub_proporcion_horizontal_labeled.png"),
      width = 18, height = 6)

    # Proportion vertical
    p8a_prop_top <- p8a_prop +
      scale_x_discrete(position = "top") +
      labs(x = "Base de datos de referencia") +
      theme(axis.text.x.top  = element_text(angle = 0, hjust = 0.5,
                                            color = "black"),
            axis.title.x.top = element_text(face = "bold", color = "black"))
    p8_prop_vert <- wrap_plots(p8a_prop_top, p8b_nox, p8c_nox, ncol = 1) +
      plot_layout(guides = "collect")
    p8_prop_vert <- p8_prop_vert & theme_publication()
    p8_prop_vert[[1]] <- p8_prop_vert[[1]] + theme(
      axis.text.x.top  = element_text(angle = 0, hjust = 0.5, color = "black"),
      axis.title.x.top = element_text(face = "bold", color = "black")
    )
    for (vi in 2:3) {
      p8_prop_vert[[vi]] <- p8_prop_vert[[vi]] + theme(
        axis.text.x  = element_blank(),
        axis.ticks.x = element_blank()
      )
    }
    save_dual(
      (p8_prop_vert + plot_annotation(tag_levels = "a")) &
        theme(plot.tag = element_text(face = "bold", size = 20, color = "black", family = "sans")),
      file.path(output_dir, "Fig8_pub_proporcion_vertical_labeled.png"),
      width = 8, height = 16)
    cat("[OK] Saved Figure 8 (proporcion)\n")

    # -------------------------------------------------------------------------
    # Fig 8c supplementary: proportion bar chart (RRS < 0.5 vs >= 0.5).
    # Saved separately as a complementary figure.
    # -------------------------------------------------------------------------
    cat("=== Generating Figure 8c supplementary (RRS proportions) ===\n")

    rrs_prop <- gunc_all %>%
      mutate(rrs_cat = ifelse(RRS < 0.5, "RRS < 0.5", "RRS \u2265 0.5")) %>%
      count(Database, pass.GUNC, rrs_cat) %>%
      group_by(Database, pass.GUNC) %>%
      mutate(pct = n / sum(n) * 100) %>%
      ungroup() %>%
      mutate(gunc_label = ifelse(pass.GUNC == "TRUE", "No quimérico", "Quimérico"))

    p8c_prop <- ggplot(rrs_prop,
                       aes(x = Database, y = pct, fill = rrs_cat)) +
      geom_col(color = "black", linewidth = 0.3, width = 0.6) +
      geom_text(aes(label = paste0(round(pct, 0), "%")),
                position = position_stack(vjust = 0.5),
                color = "black", fontface = "bold", size = 4.5) +
      facet_wrap(~ gunc_label, nrow = 1) +
      scale_fill_manual(values = c("RRS < 0.5" = "#f39c12", "RRS \u2265 0.5" = "#2980b9")) +
      scale_y_continuous(breaks = seq(0, 100, 20),
                         labels = function(x) paste0(x, "%")) +
      labs(x = "Base de datos de referencia", y = "MAGs (%)", fill = "Cobertura RRS") +
      theme_publication()
    save_dual(p8c_prop,
              file.path(output_dir, "Fig8c_sup_proportions_RRS.png"),
              width = 9, height = 5)
    cat("[OK] Saved Figure 8c supplementary\n")

    # -------------------------------------------------------------------------
    # Fig 9: CSS distribution by sampling depth and station
    #
    # Tests whether chimerism tracks the sampling context. A difference between
    # depths (OMZ, 1000 m, BTM) or stations could point to assembly or binning
    # artefacts tied to community complexity or read depth at those sites.
    # -------------------------------------------------------------------------
    cat("=== Generating Figure 9: CSS by depth and station ===\n")

    gunc_spatial <- gunc_all %>% filter(!is.na(Depth), !is.na(Station))

    # Helper to build CSS jitter+boxplot panel faceted by database
    make_css_panel <- function(data, x_var, x_label, title_str) {
      ggplot(data, aes(x = .data[[x_var]], y = CSS, fill = pass.GUNC)) +
        geom_jitter(aes(color = pass.GUNC),
                    width = 0.15, size = 2, alpha = 0.75,
                    shape = 21, stroke = 0.3) +
        geom_boxplot(width = 0.35, alpha = 0.45, outlier.shape = NA,
                     linewidth = 0.5, color = "black",
                     position = "identity") +
        geom_hline(yintercept = css_threshold,
                   linetype = "dashed", color = "black", linewidth = 0.5) +
        facet_wrap(~ Database, nrow = 1) +
        scale_fill_manual(values = pass_colors, labels = pass_labels) +
        scale_color_manual(values = pass_colors) +
        scale_y_continuous(limits = c(0, 1.08), breaks = seq(0, 1, 0.25),
                       labels = scales::label_number(drop0trailing = TRUE)) +
        labs(x = x_label,
             y = expression(bolditalic("Clade Separation Score") * bold(" (CSS)")),
             fill = "Quimerismo") +
        guides(fill = "none", color = "none") +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
    }

    p9a <- make_css_panel(gunc_spatial, "Depth",
                          "Profundidad", "CSS por profundidad")
    p9b <- make_css_panel(gunc_spatial, "Station",
                          "Estación de muestreo",    "CSS por estación")

    p9_pub <- (p9a / p9b + plot_layout(guides = "collect")) &
      theme_publication() &
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    save_dual(p9_pub, file.path(output_dir, "Fig9_pub_GUNC_CSS_spatial.png"),
              width = 12, height = 10)
    save_dual(
      (p9_pub + plot_annotation(tag_levels = "a")) &
        theme(plot.tag = element_text(face = "bold", size = 20, color = "black", family = "sans")),
      file.path(output_dir, "Fig9_pub_GUNC_CSS_spatial_labeled.png"),
      width = 12, height = 10)
    cat("[OK] Saved Figure 9\n")

    # -------------------------------------------------------------------------
    # Fig 9b: RRS distribution by sampling depth and station
    #
    # Mirrors Fig 9 (CSS) but shows Reference Representation Score instead.
    # Low RRS (<0.5) means the MAG's lineage is poorly covered by the reference
    # database, which makes the CSS call less reliable.
    # -------------------------------------------------------------------------
    cat("=== Generating Figure 9b: RRS by depth and station ===\n")

    make_rrs_panel <- function(data, x_var, x_label, title_str) {
      ggplot(data, aes(x = .data[[x_var]], y = RRS, fill = pass.GUNC)) +
        geom_jitter(aes(color = pass.GUNC),
                    position = position_jitterdodge(jitter.width = 0.15,
                                                    dodge.width  = 0.7),
                    size = 2, alpha = 0.75, shape = 21, stroke = 0.3) +
        geom_boxplot(width = 0.45, alpha = 0.45, outlier.shape = NA,
                     linewidth = 0.5, color = "black",
                     position = position_dodge(0.7)) +
        geom_hline(yintercept = 0.5,
                   linetype = "dashed", color = "gray40", linewidth = 0.5) +
        facet_wrap(~ Database, nrow = 1) +
        scale_fill_manual(values = pass_colors, labels = pass_labels) +
        scale_color_manual(values = pass_colors, guide = "none") +
        scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.25),
                       labels = scales::label_number(drop0trailing = TRUE)) +
        labs(x = x_label,
             y = expression(bolditalic("Reference Representation Score") * bold(" (RRS)")),
             fill = "Resultado GUNC") +
        guides(fill = "none", color = "none") +
        theme(axis.text.x = element_text(angle = 30, hjust = 1))
    }

    p9b_depth   <- make_rrs_panel(gunc_spatial, "Depth",
                                   "Profundidad", "RRS por profundidad")
    p9b_station <- make_rrs_panel(gunc_spatial, "Station",
                                   "Estación de muestreo",    "RRS por estación")

    p9b_pub <- (p9b_depth / p9b_station + plot_layout(guides = "collect")) &
      theme_publication() &
      theme(axis.text.x = element_text(angle = 30, hjust = 1))

    save_dual(p9b_pub, file.path(output_dir, "Fig9b_pub_GUNC_RRS_spatial.png"),
              width = 12, height = 10)
    save_dual(
      (p9b_pub + plot_annotation(tag_levels = "a")) &
        theme(plot.tag = element_text(face = "bold", size = 20, color = "black", family = "sans")),
      file.path(output_dir, "Fig9b_pub_GUNC_RRS_spatial_labeled.png"),
      width = 12, height = 10)
    cat("[OK] Saved Figure 9b\n")

    # -------------------------------------------------------------------------
    # Fig 10: Taxonomic rank at which CSS is maximised (maxCSS level)
    #
    # The rank where a genome reaches its maximum CSS shows how severe the
    # chimerism is. At broad ranks (kingdom, phylum, class) the contaminant is
    # phylogenetically distant, the worst case for downstream metabolic or
    # phylogenomic analyses. At genus or species level it points to closely
    # related strains assembled together.
    # -------------------------------------------------------------------------
    cat("=== Generating Figure 10: Taxonomic level of maxCSS ===\n")

    p10 <- ggplot(gunc_all, aes(x = taxonomic_level, fill = pass.GUNC)) +
      geom_bar(color = "black", linewidth = 0.3) +
      geom_text(stat = "count", aes(label = after_stat(count)),
                position = position_stack(vjust = 0.5),
                color = "black", fontface = "bold", size = 4.5) +
      facet_wrap(~ Database, nrow = 1) +
      scale_fill_manual(values = pass_colors, labels = pass_labels) +
      labs(x = "Nivel taxonómico", y = "Número de MAGs",
           fill = "Quimerismo") +
      guides(fill = "none") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    p10_pub <- p10 + theme_publication() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    save_dual(p10_pub, file.path(output_dir, "Fig10_pub_GUNC_tax_level.png"),
              width = 12, height = 6)
    cat("[OK] Saved Figure 10\n")

    # -------------------------------------------------------------------------
    # Fig 11: Database concordance scatter (ProGenomes3 CSS vs GTDB r214 CSS)
    #
    # A MAG flagged by both databases is the stronger chimerism call and the
    # first to inspect or re-bin. MAGs flagged by one database only are
    # labelled. The dotted lines mark CSS = 0.45 on each axis, and the diagonal
    # marks equal CSS in both databases.
    # -------------------------------------------------------------------------
    cat("=== Generating Figure 11: Database concordance ===\n")

    gunc_wide <- gunc_prog %>%
      select(genome, Name, DisplayName, Station, Depth,
             css_prog    = CSS,
             pass_prog   = pass.GUNC,
             rrs_prog    = RRS,
             contam_prog = contamination_portion) %>%
      inner_join(
        gunc_gtdb %>%
          select(genome,
                 css_gtdb    = CSS,
                 pass_gtdb   = pass.GUNC,
                 rrs_gtdb    = RRS,
                 contam_gtdb = contamination_portion),
        by = "genome"
      ) %>%
      mutate(
        Concordance = case_when(
          pass_prog == "TRUE"  & pass_gtdb == "TRUE"  ~ "Pasa en ambas",
          pass_prog == "FALSE" & pass_gtdb == "FALSE" ~ "Falla en ambas",
          TRUE                                         ~ "Discordante"
        ),
        Concordance = factor(Concordance,
                             levels = c("Pasa en ambas", "Falla en ambas", "Discordante"))
      )

    concordance_colors <- c(
      "Pasa en ambas"  = "#27ae60",
      "Falla en ambas" = "#e74c3c",
      "Discordante"    = "#f39c12"
    )

    n_concordant <- sum(gunc_wide$Concordance != "Discordante")
    n_total      <- nrow(gunc_wide)

    p11 <- ggplot(gunc_wide,
                  aes(x = css_prog, y = css_gtdb, color = Concordance)) +
      # Shaded quadrants: lower-left = both pass, upper-right = both fail
      annotate("rect", xmin = 0,             xmax = css_threshold,
               ymin = 0,             ymax = css_threshold,
               fill = "#27ae60", alpha = 0.05) +
      annotate("rect", xmin = css_threshold, xmax = 1.05,
               ymin = css_threshold, ymax = 1.05,
               fill = "#e74c3c", alpha = 0.05) +
      # Reference lines
      geom_abline(slope = 1, intercept = 0,
                  linetype = "dashed", color = "gray55", linewidth = 0.5) +
      geom_vline(xintercept = css_threshold,
                 linetype = "dashed", color = "black", linewidth = 0.5) +
      geom_hline(yintercept = css_threshold,
                 linetype = "dashed", color = "black", linewidth = 0.5) +
      geom_point(size = 3, alpha = 0.85) +
      scale_color_manual(values = concordance_colors) +
      scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2),
                         expand = expansion(mult = c(0.03, 0)),
                         labels = scales::label_number(drop0trailing = TRUE)) +
      scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2),
                         expand = expansion(mult = c(0.03, 0)),
                         labels = scales::label_number(drop0trailing = TRUE)) +
      labs(
        x = expression(bolditalic("Clade Separation Score") * bold(" (CSS) \u2013 ProGenomes3")),
        y = expression(bolditalic("Clade Separation Score") * bold(" (CSS) \u2013 GTDB r214")),
        color = "Concordancia"
      )

    p11_pub <- p11 + theme_publication() + theme(aspect.ratio = 1)

    save_dual(p11_pub, file.path(output_dir, "Fig11_pub_GUNC_concordance.png"),
              width = 9, height = 7.5)
    cat("[OK] Saved Figure 11\n")

    # Save concordance results as TSV (one row per MAG)
    gunc_concordance_tsv <- gunc_wide %>%
      select(
        MAG              = Name,
        DisplayName,
        Station,
        Depth,
        CSS_ProGenomes3  = css_prog,
        Pasa_ProGenomes3 = pass_prog,
        RRS_ProGenomes3  = rrs_prog,
        CSS_GTDB_r214    = css_gtdb,
        Pasa_GTDB_r214   = pass_gtdb,
        RRS_GTDB_r214    = rrs_gtdb,
        Concordancia     = Concordance
      )
    write_tsv(gunc_concordance_tsv,
              file.path(output_dir, "GUNC_concordance_per_MAG.tsv"))
    cat("[OK] Saved GUNC_concordance_per_MAG.tsv\n")

    # -------------------------------------------------------------------------
    # Fig 11b: the same concordance scatter, restricted to the 10 MAGs with no species
    # -------------------------------------------------------------------------
    cat("=== Generating Figure 11b: 10 MAGs sin especie concordance ===\n")

    gunc_wide_10 <- if (!is.null(cli_args$metadata) && file.exists(cli_args$metadata)) {
      species_meta_11 <- read_tsv(cli_args$metadata, show_col_types = FALSE) %>%
        select(Name = NameWithoutExtension, Species)
      gunc_wide %>% left_join(species_meta_11, by = "Name") %>%
        filter(is.na(Species) | Species == "N/A")
    } else {
      cat("[WARN] --metadata not provided; skipping Fig 11b.\n")
      NULL
    }

    if (!is.null(gunc_wide_10) && nrow(gunc_wide_10) > 0) {
      cat("  10 MAGs matched:", nrow(gunc_wide_10), "\n")

      n_conc_10 <- sum(gunc_wide_10$Concordance != "Discordante")
      n_tot_10  <- nrow(gunc_wide_10)

      p11b <- ggplot(gunc_wide_10,
                     aes(x = css_prog, y = css_gtdb, color = Concordance)) +
        annotate("rect", xmin = 0,             xmax = css_threshold,
                 ymin = 0,             ymax = css_threshold,
                 fill = "#27ae60", alpha = 0.05) +
        annotate("rect", xmin = css_threshold, xmax = 1.05,
                 ymin = css_threshold, ymax = 1.05,
                 fill = "#e74c3c", alpha = 0.05) +
        geom_abline(slope = 1, intercept = 0,
                    linetype = "dashed", color = "gray55", linewidth = 0.5) +
        geom_vline(xintercept = css_threshold,
                   linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_hline(yintercept = css_threshold,
                   linetype = "dashed", color = "black", linewidth = 0.5) +
        geom_point(size = 3.5, alpha = 0.9) +
        ggrepel::geom_text_repel(
          aes(label = DisplayName, color = Concordance),
          size = 4, fontface = "plain",
          min.segment.length = 0.2, segment.linewidth = 0.4,
          box.padding = 0.4, max.overlaps = Inf, show.legend = FALSE
        ) +
        scale_color_manual(values = concordance_colors) +
        scale_x_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2),
                           expand = expansion(mult = c(0.03, 0)),
                           labels = scales::label_number(drop0trailing = TRUE)) +
        scale_y_continuous(limits = c(0, 1.05), breaks = seq(0, 1, 0.2),
                           expand = expansion(mult = c(0.03, 0)),
                           labels = scales::label_number(drop0trailing = TRUE)) +
        labs(
          x = expression(bolditalic("Clade Separation Score") * bold(" (CSS) \u2013 ProGenomes3")),
          y = expression(bolditalic("Clade Separation Score") * bold(" (CSS) \u2013 GTDB r214")),
          color = "Concordancia"
        )

      p11b_pub <- p11b + theme_publication() + theme(aspect.ratio = 1)

      save_dual(p11b_pub,
                file.path(output_dir, "Fig11b_pub_GUNC_concordance_10MAGs_sin_especie.png"),
                width = 9, height = 7.5)
      cat("[OK] Saved Figure 11b (10 MAGs sin especie)\n")
    }

    # -------------------------------------------------------------------------
    # Fig 12: per-MAG quality summary
    #   Lollipop chart: completeness (CheckM2) | contamination (CheckM2) | CSS (GUNC)
    #
    # All the quality measures for the 54 MAGs in one figure. MAGs are ranked
    # top to bottom by decreasing GTDB r214 CSS, so the most chimeric come
    # first; GTDB r214 is the reference because its archaeal taxonomy is the
    # finer of the two. Dashed vertical lines mark the MIMAG thresholds for
    # completeness (50%, 90%) and contamination (5%, 10%), and the GUNC
    # threshold (CSS 0.45).
    # -------------------------------------------------------------------------
    cat("=== Generating Figure 12: Comprehensive per-MAG summary ===\n")

    checkm2_joined <- df_2 %>%
      inner_join(
        gunc_gtdb %>%
          select(Name, DisplayName, css_gtdb = CSS, pass_gtdb = pass.GUNC,
                 contam_gunc = contamination_portion, rrs_gtdb = RRS),
        by = "Name"
      ) %>%
      mutate(Quality_Score = Completeness - 5 * Contamination) %>%
      arrange(desc(Quality_Score)) %>%
      mutate(
        # Apply taxon_italic() so genus/species words render in italics via ggtext
        Label_md = taxon_italic(DisplayName),
        Label    = factor(Label_md, levels = rev(unique(Label_md))),
        gunc_cat = factor(
          ifelse(pass_gtdb == "TRUE", "No quimérico", "Quimérico"),
          levels = c("No quimérico", "Quimérico")
        )
      )

    cat("  CheckM2 × GUNC join:", nrow(checkm2_joined), "MAGs matched\n")

    if (nrow(checkm2_joined) > 0) {

      gunc_binary_colors <- c("No quimérico" = "#27ae60", "Quimérico" = "#e74c3c")

      # Build panels WITHOUT any y-axis theme yet.
      # theme_publication() (applied via &) includes axis.text = element_text(...)
      # which would override axis.text.y = element_blank() if set beforehand.
      # Instead, we apply per-panel y-axis and margin tweaks AFTER theme_publication
      # using patchwork's [[i]] indexing, so they take effect last.

      # Panel A: Quality Score (Completeness - 5*Contamination), with the MAG names
      qs_max <- max(checkm2_joined$Quality_Score, na.rm = TRUE)
      qs_lim <- ceiling(qs_max / 5) * 5  # round up to next multiple of 5
      p12a <- ggplot(checkm2_joined,
                     aes(y = Label, x = Quality_Score, color = MIMAG_Quality)) +
        geom_segment(aes(y = Label, yend = Label, x = 0, xend = Quality_Score),
                     color = "gray82", linewidth = 0.35) +
        geom_point(size = 2.2, alpha = 0.9) +
        geom_vline(xintercept = 50,
                   linetype = "dashed", color = "gray35", linewidth = 0.4) +
        scale_color_manual(values = mimag_colors) +
        scale_x_continuous(limits = c(0, qs_lim),
                           breaks = seq(0, qs_lim, 25),
                           expand = expansion(mult = c(0.02, 0.03)),
                           labels = scales::label_number(drop0trailing = TRUE)) +
        labs(x = "Puntuación de calidad", y = NULL,
             color = "Calidad MIMAG")

      # Panel B: Completeness (CheckM2), no legend, as in A
      p12b <- ggplot(checkm2_joined,
                     aes(y = Label, x = Completeness, color = MIMAG_Quality)) +
        geom_segment(aes(y = Label, yend = Label, x = 50, xend = Completeness),
                     color = "gray82", linewidth = 0.35) +
        geom_point(size = 2.2, alpha = 0.9) +
        geom_vline(xintercept = 90,
                   linetype = "dashed", color = "gray35", linewidth = 0.4) +
        scale_color_manual(values = mimag_colors) +
        scale_x_continuous(limits = c(50, 100),
                           breaks = c(50, 75, 90, 100),
                           expand = expansion(mult = c(0.02, 0.03))) +
        labs(x = "Completitud (%)", y = NULL,
             color = "Calidad MIMAG") +
        guides(color = "none")

      # Panel C: Contamination, no legend
      p12c <- ggplot(checkm2_joined,
                     aes(y = Label, x = Contamination, color = MIMAG_Quality)) +
        geom_segment(aes(y = Label, yend = Label, x = 0, xend = Contamination),
                     color = "gray82", linewidth = 0.35) +
        geom_point(size = 2.2, alpha = 0.9) +
        geom_vline(xintercept = c(5, 10),
                   linetype = "dashed", color = "gray35", linewidth = 0.4) +
        scale_color_manual(values = mimag_colors) +
        scale_x_continuous(limits = c(0, 10),
                           breaks = c(0, 2.5, 5, 7.5, 10),
                           expand = expansion(mult = c(0.02, 0.03)),
                           labels = scales::label_number(drop0trailing = TRUE)) +
        labs(x = "Contaminación (%)", y = NULL,
             color = "Calidad MIMAG") +
        guides(color = "none")

      # Panel D: CSS, GTDB r214
      p12d <- ggplot(checkm2_joined,
                     aes(y = Label, x = css_gtdb, color = gunc_cat)) +
        geom_segment(aes(y = Label, yend = Label, x = 0, xend = css_gtdb),
                     color = "gray82", linewidth = 0.35) +
        geom_point(size = 2.2, alpha = 0.9) +
        geom_vline(xintercept = css_threshold,
                   linetype = "dashed", color = "black", linewidth = 0.5) +
        scale_color_manual(values = gunc_binary_colors) +
        scale_x_continuous(limits = c(0, 1.15),
                           breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                           expand = expansion(mult = c(0, 0)),
                           labels = scales::label_number(drop0trailing = TRUE)) +
        labs(x = expression(bolditalic("Clade Separation Score") * bold(" (CSS)")),
             y = NULL,
             color = "Quimerismo")

      # Panel E: Reference Representation Score, GTDB r214
      p12e <- ggplot(checkm2_joined,
                     aes(y = Label, x = rrs_gtdb, color = gunc_cat)) +
        geom_segment(aes(y = Label, yend = Label, x = 0, xend = rrs_gtdb),
                     color = "gray82", linewidth = 0.35) +
        geom_point(size = 2.2, alpha = 0.9) +
        geom_vline(xintercept = 0.5,
                   linetype = "dashed", color = "gray40", linewidth = 0.4) +
        scale_color_manual(values = gunc_binary_colors) +
        scale_x_continuous(limits = c(0, 1.05),
                           breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                           expand = expansion(mult = c(0.02, 0.03)),
                           labels = scales::label_number(drop0trailing = TRUE)) +
        labs(x = expression(bolditalic("Reference Representation Score") * bold(" (RRS)")),
             y = NULL,
             color = "Quimerismo") +
        guides(color = "none")

      # guides = "collect" merges all legends to the right of the patchwork
      p12 <- (p12a | p12b | p12c | p12d | p12e) +
        plot_layout(widths = c(1, 1, 1, 1, 1), guides = "collect")

      # Apply global theme first, then per-panel overrides via [[i]] so that
      # theme_publication() cannot override axis.text.y = element_blank().
      p12_pub <- p12 & theme_publication()

      # Panel A: MAG names rendered via element_markdown to enable italic taxa
      p12_pub[[1]] <- p12_pub[[1]] + theme(
        axis.text.y = ggtext::element_markdown(size = 14),
        plot.margin = margin(t = 4, r = 8, b = 4, l = 5, unit = "pt")
      )
      # Panels B-E: hide y-axis entirely
      for (i in 2:5) {
        p12_pub[[i]] <- p12_pub[[i]] + theme(
          axis.text.y         = element_blank(),
          axis.ticks.y        = element_blank(),
          axis.ticks.length.y = unit(0, "pt"),
          plot.margin         = margin(t = 4, r = 6, b = 4, l = 5, unit = "pt")
        )
      }
      p12_pub[[4]] <- p12_pub[[4]] + theme(
        axis.title.x = element_text(face = "bold", color = "black")
      )

      save_dual(p12_pub,
                file.path(output_dir, "Fig12_pub_GUNC_summary_per_MAG.png"),
                width = 25, height = 13)

      # Version without RRS panel (panels A–D: Quality Score, Completeness,
      # Contamination, CSS)
      p12_no_rrs <- (p12a | p12b | p12c | p12d) +
        plot_layout(widths = c(1, 1, 1, 1), guides = "collect")
      p12_pub_no_rrs <- p12_no_rrs & theme_publication()
      p12_pub_no_rrs[[1]] <- p12_pub_no_rrs[[1]] + theme(
        axis.text.y = ggtext::element_markdown(size = 14),
        plot.margin = margin(t = 4, r = 8, b = 4, l = 5, unit = "pt")
      )
      for (i in 2:4) {
        p12_pub_no_rrs[[i]] <- p12_pub_no_rrs[[i]] + theme(
          axis.text.y         = element_blank(),
          axis.ticks.y        = element_blank(),
          axis.ticks.length.y = unit(0, "pt"),
          plot.margin         = margin(t = 4, r = 6, b = 4, l = 5, unit = "pt")
        )
      }
      p12_pub_no_rrs[[4]] <- p12_pub_no_rrs[[4]] + theme(
        axis.title.x = element_text(face = "bold", color = "black")
      )
      save_dual(p12_pub_no_rrs,
                file.path(output_dir, "Fig12_pub_GUNC_summary_per_MAG_sin_RRS.png"),
                width = 20, height = 13)
      cat("[OK] Saved Figure 12\n")

      # -----------------------------------------------------------------------
      # Fig 12b: same layout as Fig 12 but restricted to the 10 MAGs that have
      # no species-level assignment in the GTDB taxonomy (Species == "N/A").
      # These are the focal MAGs for the phylogenomic section of the thesis.
      # -----------------------------------------------------------------------
      cat("=== Generating Figure 12b: 10 MAGs sin especie ===\n")

      checkm2_10 <- if (!is.null(cli_args$metadata) && file.exists(cli_args$metadata)) {
        species_meta <- read_tsv(cli_args$metadata, show_col_types = FALSE) %>%
          select(Name = NameWithoutExtension, Species)
        checkm2_joined %>%
          left_join(species_meta, by = "Name") %>%
          filter(is.na(Species) | Species == "N/A")
      } else {
        cat("[WARN] --metadata not provided; cannot filter 10 MAGs. Skipping Fig 12b.\n")
        NULL
      }

      if (!is.null(checkm2_10) && nrow(checkm2_10) > 0) {
        cat("  10 MAGs sin especie matched:", nrow(checkm2_10), "\n")

        # Rebuild factor levels for the subset ordered by Quality_Score desc
        checkm2_10 <- checkm2_10 %>%
          arrange(desc(Quality_Score)) %>%
          mutate(Label = factor(Label_md, levels = rev(unique(Label_md))))

        qs_max_10 <- max(checkm2_10$Quality_Score, na.rm = TRUE)
        qs_lim_10 <- ceiling(qs_max_10 / 5) * 5

        p12b_a <- ggplot(checkm2_10,
                         aes(y = Label, x = Quality_Score, color = MIMAG_Quality)) +
          geom_segment(aes(y = Label, yend = Label, x = 0, xend = Quality_Score),
                       color = "gray82", linewidth = 0.35) +
          geom_point(size = 2.5, alpha = 0.9) +
          geom_vline(xintercept = 50,
                     linetype = "dashed", color = "gray35", linewidth = 0.4) +
          scale_color_manual(values = mimag_colors) +
          scale_x_continuous(limits = c(0, qs_lim_10),
                             breaks = seq(0, qs_lim_10, 25),
                             expand = expansion(mult = c(0.02, 0.03)),
                             labels = scales::label_number(drop0trailing = TRUE)) +
          labs(x = "Puntuación de calidad", y = NULL,
               color = "Calidad MIMAG")

        p12b_b <- ggplot(checkm2_10,
                         aes(y = Label, x = Completeness, color = MIMAG_Quality)) +
          geom_segment(aes(y = Label, yend = Label, x = 50, xend = Completeness),
                       color = "gray82", linewidth = 0.35) +
          geom_point(size = 2.5, alpha = 0.9) +
          geom_vline(xintercept = 90,
                     linetype = "dashed", color = "gray35", linewidth = 0.4) +
          scale_color_manual(values = mimag_colors) +
          scale_x_continuous(limits = c(50, 100), breaks = c(50, 75, 90, 100),
                             expand = expansion(mult = c(0.02, 0.03))) +
          labs(x = "Completitud (%)", y = NULL,
               color = "Calidad MIMAG") +
          guides(color = "none")

        p12b_c <- ggplot(checkm2_10,
                         aes(y = Label, x = Contamination, color = MIMAG_Quality)) +
          geom_segment(aes(y = Label, yend = Label, x = 0, xend = Contamination),
                       color = "gray82", linewidth = 0.35) +
          geom_point(size = 2.5, alpha = 0.9) +
          geom_vline(xintercept = c(5, 10),
                     linetype = "dashed", color = "gray35", linewidth = 0.4) +
          scale_color_manual(values = mimag_colors) +
          scale_x_continuous(limits = c(0, 10), breaks = c(0, 2.5, 5, 7.5, 10),
                             expand = expansion(mult = c(0.02, 0.03)),
                             labels = scales::label_number(drop0trailing = TRUE)) +
          labs(x = "Contaminación (%)", y = NULL,
               color = "Calidad MIMAG") +
          guides(color = "none")

        p12b_d <- ggplot(checkm2_10,
                         aes(y = Label, x = css_gtdb, color = gunc_cat)) +
          geom_segment(aes(y = Label, yend = Label, x = 0, xend = css_gtdb),
                       color = "gray82", linewidth = 0.35) +
          geom_point(size = 2.5, alpha = 0.9) +
          geom_vline(xintercept = css_threshold,
                     linetype = "dashed", color = "black", linewidth = 0.5) +
          scale_color_manual(values = gunc_binary_colors) +
          scale_x_continuous(limits = c(0, 1.15), breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                             expand = expansion(mult = c(0, 0)),
                             labels = scales::label_number(drop0trailing = TRUE)) +
          labs(x = expression(bolditalic("Clade Separation Score") * bold(" (CSS)")),
               y = NULL,
               color = "Quimerismo")

        p12b_e <- ggplot(checkm2_10,
                         aes(y = Label, x = rrs_gtdb, color = gunc_cat)) +
          geom_segment(aes(y = Label, yend = Label, x = 0, xend = rrs_gtdb),
                       color = "gray82", linewidth = 0.35) +
          geom_point(size = 2.5, alpha = 0.9) +
          geom_vline(xintercept = 0.5,
                     linetype = "dashed", color = "gray40", linewidth = 0.4) +
          scale_color_manual(values = gunc_binary_colors) +
          scale_x_continuous(limits = c(0, 1.05), breaks = c(0, 0.25, 0.5, 0.75, 1.0),
                             expand = expansion(mult = c(0.02, 0.03)),
                             labels = scales::label_number(drop0trailing = TRUE)) +
          labs(x = expression(bolditalic("Reference Representation Score") * bold(" (RRS)")),
               y = NULL,
               color = "Quimerismo") +
          guides(color = "none")

        p12b_patch <- (p12b_a | p12b_b | p12b_c | p12b_d | p12b_e) +
          plot_layout(widths = c(1, 1, 1, 1, 1), guides = "collect")

        p12b_pub <- p12b_patch & theme_publication()

        p12b_pub[[1]] <- p12b_pub[[1]] + theme(
          axis.text.y = ggtext::element_markdown(size = 14),
          plot.margin = margin(t = 4, r = 8, b = 4, l = 5, unit = "pt")
        )
        for (i in 2:5) {
          p12b_pub[[i]] <- p12b_pub[[i]] + theme(
            axis.text.y         = element_blank(),
            axis.ticks.y        = element_blank(),
            axis.ticks.length.y = unit(0, "pt"),
            plot.margin         = margin(t = 4, r = 6, b = 4, l = 5, unit = "pt")
          )
        }
        p12b_pub[[4]] <- p12b_pub[[4]] + theme(
          axis.title.x = element_text(face = "bold", color = "black")
        )

        save_dual(p12b_pub,
                  file.path(output_dir, "Fig12b_pub_GUNC_summary_10MAGs_sin_especie.png"),
                  width = 25, height = 7)

        # Version without RRS panel
        p12b_no_rrs <- (p12b_a | p12b_b | p12b_c | p12b_d) +
          plot_layout(widths = c(1, 1, 1, 1), guides = "collect")
        p12b_pub_no_rrs <- p12b_no_rrs & theme_publication()
        p12b_pub_no_rrs[[1]] <- p12b_pub_no_rrs[[1]] + theme(
          axis.text.y = ggtext::element_markdown(size = 14),
          plot.margin = margin(t = 4, r = 8, b = 4, l = 5, unit = "pt")
        )
        for (i in 2:4) {
          p12b_pub_no_rrs[[i]] <- p12b_pub_no_rrs[[i]] + theme(
            axis.text.y         = element_blank(),
            axis.ticks.y        = element_blank(),
            axis.ticks.length.y = unit(0, "pt"),
            plot.margin         = margin(t = 4, r = 6, b = 4, l = 5, unit = "pt")
          )
        }
        p12b_pub_no_rrs[[4]] <- p12b_pub_no_rrs[[4]] + theme(
          axis.title.x = element_text(face = "bold", color = "black")
        )
        save_dual(p12b_pub_no_rrs,
                  file.path(output_dir, "Fig12b_pub_GUNC_summary_10MAGs_sin_especie_sin_RRS.png"),
                  width = 20, height = 7)
        cat("[OK] Saved Figure 12b (10 MAGs sin especie)\n")

      } else if (!is.null(checkm2_10)) {
        cat("[WARN] No MAGs with Species == 'N/A' found after join; skipping Fig 12b.\n")
      }

    } else {
      cat("[WARN] No CheckM2–GUNC matches found; skipping Figure 12.\n")
      cat("  GUNC Name examples  :", paste(head(gunc_gtdb$Name, 3), collapse = ", "), "\n")
      cat("  CheckM2 Name examples:", paste(head(df_2$Name,     3), collapse = ", "), "\n")
    }

  } else {
    cat("[WARN] Could not load one or both GUNC files. Skipping Figures 8–12.\n")
  }

} else {
  cat("\n[WARN] GUNC files not provided (--gunc-progenomes / --gunc-gtdb).",
      "Skipping chimerism figures (8–12).\n")
}

# === Table 1: summary statistics ===
table1 <- df_combined %>%
  group_by(Dataset) %>%
  summarise(
    Total_MAGs = n(),
    High_Quality = sum(MIMAG_Quality == "Calidad alta"),
    Medium_Quality = sum(MIMAG_Quality == "Calidad media"),
    Low_Quality = sum(MIMAG_Quality == "Calidad baja"),
    Contaminated = sum(MIMAG_Quality == "Contaminado"),
    Mean_Completeness = round(mean(Completeness), 2),
    SD_Completeness = round(sd(Completeness), 2),
    Mean_Contamination = round(mean(Contamination), 2),
    SD_Contamination = round(sd(Contamination), 2),
    Mean_Genome_Size_Mb = round(mean(Genome_Size) / 1e6, 2),
    Mean_N50_kb = round(mean(Contig_N50) / 1000, 2),
    Mean_GC = round(mean(GC_Content) * 100, 2)
  )

write.csv(table1, file.path(output_dir, "Table1_summary_statistics.csv"),
          row.names = FALSE)

# CONTEXT SUMMARY OF THE FINAL 54 MAGs (STATION, DEPTH AND QUALITY)
# Faceted bar chart that opens the results section
cat("\n=== Generando gráfico de contexto de los 54 MAGs finales ===\n")

# df_valentin rather than df_combined, so this covers only the 54 MAGs of this
# study and not the comparison dataset.
df_54 <- df_valentin %>%
  filter(MIMAG_Quality %in% c("Calidad alta", "Calidad media")) %>%
  mutate(
    # Quality score, following Parks et al.
    Quality_Score = Completeness - 5 * Contamination,
    # Translate the depths
    Depth = recode(Depth, "OMZ" = "ZMOR", "1000" = "1000 m", "BTM" = "Fondo"),
    # Order the depths from surface to bottom
    Depth = factor(Depth, levels = c("ZMOR", "1000 m", "Fondo"))
  )

# 2. Summarise the data
summary_54 <- df_54 %>%
  group_by(Station, Depth) %>%
  summarise(
    n_MAGs = n(),
    Mean_Quality_Score = round(mean(Quality_Score), 2),
    .groups = "drop"
  )

# Per-depth totals, annotated down the right-hand side
depth_totals <- summary_54 %>%
  group_by(Depth) %>%
  summarise(n_total = sum(n_MAGs), .groups = "drop")

# Build the faceted bar chart
p_bar_context <- ggplot(summary_54, aes(x = Station, y = n_MAGs, fill = Mean_Quality_Score)) +
  geom_col(color = "black", width = 0.7, linewidth = 0.4) +
  geom_text(aes(label = n_MAGs), vjust = -0.5, size = 4, fontface = "bold") +
  # Depth total in the top right corner of each facet
  geom_text(
    data = depth_totals,
    aes(x = Inf, y = Inf, label = paste0("n = ", n_total)),
    hjust = 1.08, vjust = 1.6,
    size = 4.5, fontface = "bold", color = "gray20",
    inherit.aes = FALSE
  ) +
  # Legend without the heavy capitalisation
  scale_fill_viridis_c(option = "mako", name = "Puntuación de\ncalidad media") +
  # Headroom at the top, so the counts do not sit against the facet edge
  scale_y_continuous(expand = expansion(mult = c(0, 0.3))) +
  # One row per depth
  facet_grid(Depth ~ .) +
  labs(
    x = "Estación de muestreo",
    y = "Número de MAGs recuperados"
  ) +
  # Using the shared publication theme
  theme_publication(base_size = 14) +
  theme(
    strip.background = element_rect(fill = "gray90", color = "black"),
    strip.text = element_text(face = "bold", size = 12)
  )

ggsave(file.path(output_dir, "Fig_Context_Bars.png"), plot = p_bar_context, width = 8, height = 6, dpi = 300)

cat("[OK] Gráfico de contexto introductorio generado exitosamente (Fig_Context_Bars.png).\n")

# === Print summary to console ===
cat("\n=== MAG quality assessment Summary ===\n\n")
print(table1)

cat("\n=== Quality distribution ===\n")
print(summary_stats)

cat("\n=== Analysis complete! ===\n")
cat("Output directory:", output_dir, "\n")
gunc_figs <- if (!is.null(cli_args$gunc_progenomes) && !is.null(cli_args$gunc_gtdb)) 5L else 0L
cat("Figures generated:", 7L + gunc_figs,
    if (gunc_figs > 0) "(7 CheckM2 + 5 GUNC chimerism)\n" else "(CheckM2 only; use --gunc-* for chimerism figures)\n")
cat("Tables generated: 2 (CSV + XLSX)\n\n")

cat(sprintf("\n=== Finalización del script: %s ===\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S")))
cat(sprintf("    Duración total: %.1f minutos\n",
            as.numeric(difftime(Sys.time(), script_start, units = "mins"))))
