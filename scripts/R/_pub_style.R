# Shared look for every figure in the project, plus the path resolution the
# other R scripts rely on. Load it first:
#
#   source(file.path(SCRIPTS_R_DIR, "_pub_style.R"))
#
# Font selection walks four fallbacks in order: Lato or Source Sans 3 from
# ~/.local/share/fonts, then a matching system face, then Google Fonts, then
# plain "sans". Figures render on whatever machine opens them, at the cost of
# not being byte-identical across machines.
#
# The figure labels here are in Spanish because the thesis is. They are output,
# not documentation, so they stay that way.
#
# What this exports:
#   resolve_project_dir()  : project root, from MSC_ARCHAEA_DIR or the usual paths
#   font_family            : the family that was actually found
#   theme_publication()    : the ggplot2 theme every figure starts from
#   depth_levels_es        : c("ZMOR", "1000 m", "Fondo")
#   depth_colors_es        : colour per depth
#   station_levels         : c("A4", "A6", "C22", "F37", "PO1", "Y6")
#   station_colors         : Set1
#   pan_class_levels       : c("core", "soft-core", "shell", "cloud")
#   pan_class_colors       : Dark2-like
#   taxon_italic(x, level) : italicises genus and below, per CDC style

# Project root. MSC_ARCHAEA_DIR wins; otherwise the first candidate that exists.
# RStudio on Windows reaches a WSL tree through a UNC path, which is why that
# form is tried before the native one.
resolve_project_dir <- function() {
  from_env <- Sys.getenv("MSC_ARCHAEA_DIR")
  if (nzchar(from_env) && dir.exists(from_env)) {
    return(normalizePath(from_env, mustWork = FALSE))
  }

  # The UNC entry is for RStudio on Windows reading a WSL tree; edit the distro
  # and user name to match yours, or set MSC_ARCHAEA_DIR and skip it entirely.
  candidates <- c(
    "//wsl.localhost/Ubuntu-24.04/home/ruben/msc_archaea",
    path.expand("~/msc_archaea")
  )
  found <- candidates[dir.exists(candidates)]

  if (length(found) == 0) {
    stop("Project root not found. Set MSC_ARCHAEA_DIR to the repository root.",
         call. = FALSE)
  }
  found[1]
}

suppressPackageStartupMessages({
  for (pkg in c("showtext", "sysfonts", "ggplot2", "ggtext", "RColorBrewer")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg, dependencies = TRUE, quiet = TRUE)
    }
    library(pkg, character.only = TRUE)
  }
})

showtext_auto()
showtext_opts(dpi = 300)

font_family <- "sans"
available_fonts <- tryCatch(sysfonts::font_files(), error = function(e) NULL)

font_priority <- c(
  "Source Sans Pro", "SourceSansPro",
  "CMU Serif", "CMU Sans Serif", "Computer Modern", "Latin Modern",
  "Helvetica", "Helvetica Neue"
)

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
        message("[_pub_style] Font: ", lf$name, " (from ~/.local/share/fonts)")
        font_found <- TRUE
        break
      }, error = function(e) {
        message("[_pub_style] ", lf$name, " found but would not load: ", conditionMessage(e))
      })
    }
  }
}

if (!font_found && !is.null(available_fonts)) {
  for (font_name in font_priority) {
    matching_fonts <- available_fonts[grepl(font_name, available_fonts$family, ignore.case = TRUE), ]
    if (nrow(matching_fonts) > 0) {
      tryCatch({
        sysfonts::font_add(family = "custom_font", regular = matching_fonts$file[1])
        font_family <- "custom_font"
        message("[_pub_style] Font: ", font_name, " (system)")
        font_found <- TRUE
        break
      }, error = function(e) {})
    }
  }
}
if (!font_found) {
  google_priority <- list(
    list(name = "Source Sans Pro"),
    list(name = "Lato"),
    list(name = "Open Sans")
  )
  for (gf in google_priority) {
    tryCatch({
      sysfonts::font_add_google(gf$name, "custom_font")
      font_family <- "custom_font"
      message("[_pub_style] Font: ", gf$name, " (downloaded from Google Fonts)")
      font_found <- TRUE
      break
    }, error = function(e) {})
  }
}
if (!font_found) {
  message("[_pub_style] No preferred font available; falling back to 'sans'.")
}

# ---- theme_publication ----
theme_publication <- function(base_size = 16) {
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      panel.grid.minor   = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_line(color = "gray", linetype = "dotted"),
      panel.grid.major.x = ggplot2::element_blank(),
      axis.text          = ggplot2::element_text(color = "black", size = base_size - 1),
      axis.title         = ggplot2::element_text(face = "bold", size = base_size + 2),
      legend.position    = "right",
      legend.background  = ggplot2::element_rect(fill = "white", color = NA),
      legend.title       = ggplot2::element_text(face = "bold", size = base_size),
      legend.text        = ggplot2::element_text(size = base_size - 3),
      strip.background   = ggplot2::element_rect(fill = "gray90", color = "black"),
      strip.text         = ggplot2::element_text(face = "bold", size = base_size),
      plot.title         = ggtext::element_markdown(face = "bold", hjust = 0.5, size = base_size + 2),
      plot.subtitle      = ggtext::element_markdown(hjust = 0.5, size = base_size),
      text               = ggplot2::element_text(family = font_family)
    )
}

# ---- shared palettes ----
depth_levels_es  <- c("ZMOR", "1000 m", "Fondo")
depth_colors_es  <- c("ZMOR" = "#4A90D9", "1000 m" = "#E67E22", "Fondo" = "#8E44AD")
depth_recode_es  <- c("ROMZ" = "ZMOR", "1000" = "1000 m", "BTM" = "Fondo")

station_levels <- c("A4", "A6", "C22", "F37", "PO1", "Y6")
station_colors <- setNames(RColorBrewer::brewer.pal(6, "Set1"), station_levels)

pan_class_levels <- c("core", "soft-core", "shell", "cloud")
pan_class_colors <- c(
  core        = "#1B9E77",
  `soft-core` = "#66A61E",
  shell       = "#E6AB02",
  cloud       = "#D95F02"
)

# ---- taxon_italic ----
# CDC style: genus and below are italic, family and above are not. A word
# carrying a digit or written in all caps is a strain or clade code, so it
# stays upright even inside an italic name.
ITALIC_LEVELS <- c("Genus", "Species", "Subspecies", "Strain")

should_italicize_level <- function(level) {
  !is.null(level) && length(level) == 1 && !is.na(level) && level %in% ITALIC_LEVELS
}

taxon_italic <- function(x, level = "Species") {
  if (!should_italicize_level(level)) {
    return(ifelse(is.na(x), x, as.character(x)))
  }
  sapply(x, function(nm) {
    if (is.na(nm) || nm == "Sin asignar") return(nm)
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
