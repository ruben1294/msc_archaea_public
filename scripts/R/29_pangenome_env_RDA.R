# Splits the variance in accessory-genome Jaccard dissimilarity between two
# explanatory matrices, with vegan::varpart:
#
#   X1  physicochemical variables, an adaptive subset chosen against the
#       degrees of freedom available
#   X2  PVR vectors, phylogenetic eigenvector regression on the subtree for
#       that species
#
# Environment and phylogeny are confounded here: MAGs from the same depth also
# tend to be related. Fitting either alone credits it with variance the other
# could explain as well, so varpart separates the shared fraction from the two
# independent ones.
#
# Only species with at least five MAGs are processed.
#
# Physicochemical_properties.xlsx carries two sheets:
#   Physicochemical_by_depth   one row per depth (ROMZ, 1000, BTM)
#   Stations_depths            one row per station and depth
#
# Pressure is derived here rather than measured, as Collection_depth_m / 10,
# taking 1 atm per 10 m.

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(readr)
  library(readxl)
  library(ggplot2)
  library(ggtext)
  library(stringr)
  library(vegan)
  library(ape)
})

# ---- Spanish labels inside vegan's base-graphics plots ----
# plot.varpart234 and showvarparts draw English literals ("Values <0 not shown",
# "Residuals = ...") that no parameter reaches. The only way to translate them is
# to patch vegan's namespace at run time. The patch lasts for this R session and
# never touches the installed package.
local({
  ns <- asNamespace("vegan")
  if (exists("plot.varpart234", envir = ns)) {
    fn <- get("plot.varpart234", envir = ns)
    body_chr <- deparse(body(fn))
    body_chr <- gsub('"Values <"',   '"Valores <"',     body_chr, fixed = TRUE)
    body_chr <- gsub('" not shown"', '" no mostrados"', body_chr, fixed = TRUE)
    body(fn) <- parse(text = paste(body_chr, collapse = "\n"))[[1]]
    assignInNamespace("plot.varpart234", fn, ns = "vegan")
  }
  if (exists("showvarparts", envir = ns)) {
    fn <- get("showvarparts", envir = ns)
    body_chr <- deparse(body(fn))
    body_chr <- gsub('"Residuals ="', '"Residuales ="', body_chr, fixed = TRUE)
    body(fn) <- parse(text = paste(body_chr, collapse = "\n"))[[1]]
    assignInNamespace("showvarparts", fn, ns = "vegan")
  }
})

# Converts a display name carrying markdown italics (*X*) into an expression()
# that base R can render through title() and mtext(). plot.varpart uses base
# graphics, where the asterisks would otherwise print literally.
italic_title <- function(pretty) {
  if (is.null(pretty) || is.na(pretty)) return("")
  m <- regmatches(pretty, regexec("^(.*?)\\*([^*]+)\\*(.*)$", pretty))[[1]]
  if (length(m) != 4) return(pretty)
  prefix <- m[2]; italic_part <- m[3]; suffix <- m[4]
  suffix_trim <- sub("^\\s+", "", suffix)
  if (nzchar(prefix) && nzchar(suffix_trim)) {
    bquote(.(prefix) * italic(.(italic_part)) ~ .(suffix_trim))
  } else if (nzchar(suffix_trim)) {
    bquote(italic(.(italic_part)) ~ .(suffix_trim))
  } else if (nzchar(prefix)) {
    bquote(.(prefix) * italic(.(italic_part)))
  } else {
    bquote(italic(.(italic_part)))
  }
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
SO_03_DIR   <- file.path(PROJECT_DIR, "results", "SO_03")
DERIVED_DIR <- file.path(SO_03_DIR, "derived")
CLASS_DIR   <- file.path(SO_03_DIR, "classification")
OUTPUT_DIR  <- file.path(SO_03_DIR, "rda_env")
ENV_XLSX    <- file.path(PROJECT_DIR, "data", "Physicochemical_properties.xlsx")
TREE_NWK    <- file.path(PROJECT_DIR, "results", "SO_02", "phylogenomics",
                         "curated_54mags.gtdb_pruned.nwk")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# Shared fonts, theme and palettes
SCRIPTS_R_DIR <- file.path(PROJECT_DIR, "scripts", "R")
source(file.path(SCRIPTS_R_DIR, "_pub_style.R"))

species_pretty <- c(
  Thalassarchaeum_sp002507125 = "*Thalassarchaeum* sp002507125",
  MGIIb_O1_sp002498525        = "MGIIb-O1 sp002498525",
  UBA8886_sp003193815         = "UBA8886 sp003193815"
)

# --- Helpers ------------------------------------------------------------------
log_info <- function(...) cat(sprintf("[%s] [INFO] %s\n",
                                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                      paste0(...)))
log_warn <- function(...) cat(sprintf("[%s] [WARN] %s\n",
                                      format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
                                      paste0(...)), file = stderr())

# Must match the sanitisation in pipeline step 20
sanitize <- function(x) gsub("[^a-zA-Z0-9_]", "_", x)

# PVR axes from a tree: PCoA over the patristic distance matrix
phylo_pvr <- function(tree, tips_keep, n_axes = 3) {
  tree2 <- keep.tip(tree, intersect(tips_keep, tree$tip.label))
  if (length(tree2$tip.label) < 4) return(NULL)
  pdist <- cophenetic(tree2)
  pco <- cmdscale(as.dist(pdist), k = min(n_axes, length(tree2$tip.label) - 1),
                  eig = TRUE)
  axes <- as_tibble(pco$points, rownames = "tip") %>%
    rename_with(~ paste0("PVR", seq_along(.)), -tip)
  list(axes = axes, eigvals = pco$eig)
}

# --- Load environmental data -------------------------------------------------
if (!file.exists(ENV_XLSX)) {
  stop("Environmental XLSX not found: ", ENV_XLSX,
       " — run scripts/utils/09_combine_env_metadata.py first")
}
env_depth <- read_xlsx(ENV_XLSX, sheet = "Physicochemical_by_depth") %>%
  rename(
    Temperature   = `mean_ConservativeTemp (°C)`,
    Salinity      = `mean_AbsoluteSalinity (g/kg)`,
    O2_dissolved  = mean_DissolvedOxigen,
    AOU           = `mean_CAO (umol/kg)`,
    NO3_NO2       = mean_NO3_NO2,
    PO4           = mean_PO4
  ) %>%
  select(Depth, Temperature, Salinity, O2_dissolved, AOU, NO3_NO2, PO4)

env_station <- read_xlsx(ENV_XLSX, sheet = "Stations_depths") %>%
  mutate(Pressure_atm = Collection_depth_m / 10) %>%
  select(Station, Depth, Seafloor_depth_m, Density_kg_m3,
         Collection_depth_m, Pressure_atm)

log_info(sprintf("Env depth-aggregated rows: %d", nrow(env_depth)))
log_info(sprintf("Env station×depth rows:    %d", nrow(env_station)))

# ---- load the tree ----
# The pruned tree from SO_02 labels its tips with the display form, for example
# "MGIIb-O1_sp002498525_F37_BTM_20", while the binary matrices are keyed on the
# sanitised names anvi-gen-genomes-storage produced. The tips are renamed to the
# sanitised key so the phylogenetic conditioning in the db-RDA binds to
# the right genomes instead of silently matching nothing.
tree <- NULL
if (file.exists(TREE_NWK)) {
  tree <- read.tree(TREE_NWK)
  log_info(sprintf("Tree loaded: %d tips (raw)", length(tree$tip.label)))

  META_TSV <- file.path(PROJECT_DIR, "data", "MAGs_metadata.tsv")
  meta_for_tree <- read_tsv(META_TSV, show_col_types = FALSE) %>%
    mutate(parts = str_split(NameWithoutExtension, fixed("."))) %>%
    rowwise() %>%
    mutate(
      Station_raw = parts[1],
      Depth_raw   = parts[2],
      BinID_raw   = parts[length(parts)],
      tree_suffix = paste0("_", Station_raw, "_", Depth_raw, "_", BinID_raw),
      genome_key  = sanitize(NameWithoutExtension)
    ) %>%
    ungroup() %>%
    select(NameWithoutExtension, genome_key, tree_suffix)

  new_labels <- tree$tip.label
  matched <- 0L
  for (i in seq_along(tree$tip.label)) {
    tip <- tree$tip.label[i]
    hits <- meta_for_tree[endsWith(tip, meta_for_tree$tree_suffix), ]
    if (nrow(hits) == 1L) {
      new_labels[i] <- hits$genome_key
      matched <- matched + 1L
    }
  }
  tree$tip.label <- new_labels
  log_info(sprintf("Tree tips renamed: %d / %d matched MAG metadata",
                   matched, length(tree$tip.label)))
} else {
  log_warn("Tree not found at ", TREE_NWK,
           " — RDA will be unconditioned")
}

# ---- one species at a time ----
varpart_summary <- list()
panel_varpart   <- list()  # holds vp and metadata for Fig_panel_varpart_venn

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

for (bin_file in bin_files) {
  sp_name <- str_replace(basename(bin_file), "_GCxMAG_binary\\.tsv$", "")
  pretty  <- species_pretty[sp_name]; if (is.na(pretty)) pretty <- sp_name
  log_info(sprintf("Species: %s", sp_name))

  bin_df <- read_tsv(bin_file, show_col_types = FALSE)
  meta_file <- file.path(DERIVED_DIR, paste0(sp_name, "_MAG_metadata.tsv"))
  cls_file  <- file.path(CLASS_DIR,  paste0("pangenome_classes_", sp_name, ".tsv"))
  if (!file.exists(meta_file) || !file.exists(cls_file)) {
    log_warn("  missing metadata or classification; skipping"); next
  }

  meta_sp <- read_tsv(meta_file, show_col_types = FALSE)
  cls     <- read_tsv(cls_file,  show_col_types = FALSE)
  accessory_gcs <- cls %>% filter(class %in% c("shell", "cloud")) %>% pull(gene_cluster_id)
  if (length(accessory_gcs) < 5) {
    log_warn("  <5 accessory GCs; skipping RDA"); next
  }

  bin_mat <- as.matrix(bin_df[, -1])
  rownames(bin_mat) <- bin_df$gene_cluster_id
  bin_acc <- bin_mat[rownames(bin_mat) %in% accessory_gcs, , drop = FALSE]
  Y <- t(bin_acc)  # rows = MAGs, columns = GCs

  # Station and Depth join the per-station sheet; Depth alone joins the other
  env_join <- meta_sp %>%
    left_join(env_station, by = c("Station", "Depth")) %>%
    left_join(env_depth,   by = "Depth")

  # Align Y to the metadata order, then drop MAGs with incomplete environment
  env_aligned <- env_join[match(rownames(Y), env_join$genome), ]
  full_rows <- complete.cases(env_aligned %>%
    select(Temperature, Salinity, O2_dissolved, AOU, NO3_NO2, PO4,
           Pressure_atm, Density_kg_m3))
  if (sum(full_rows) < 4) {
    log_warn(sprintf("  only %d MAGs with complete env data; skipping RDA",
                     sum(full_rows)))
    next
  }
  Y <- Y[full_rows, , drop = FALSE]
  env_aligned <- env_aligned[full_rows, ]

  # PVR axes from the within-species subtree; tips are already renamed above
  # to the sanitised genome_key, so the match is direct.
  pvr <- if (!is.null(tree)) phylo_pvr(tree, rownames(Y), n_axes = 3) else NULL
  pvr_axes <- NULL
  if (!is.null(pvr)) {
    pvr_axes <- pvr$axes %>% rename(genome = tip)
    env_aligned <- env_aligned %>% left_join(pvr_axes, by = "genome")
    log_info(sprintf("  PVR axes attached: %d/%d MAGs covered",
                     sum(!is.na(env_aligned$PVR1)), nrow(env_aligned)))
    # varpart and mantel both need complete cases in the conditioning terms
    pvr_cols_present <- intersect(c("PVR1","PVR2","PVR3"), colnames(env_aligned))
    if (length(pvr_cols_present) > 0) {
      pvr_ok <- complete.cases(env_aligned[, pvr_cols_present, drop = FALSE])
      if (any(!pvr_ok)) {
        dropped <- env_aligned$genome[!pvr_ok]
        log_warn(sprintf("  %d MAG(s) without PVR (no tree match) dropped: %s",
                         length(dropped), paste(dropped, collapse = ", ")))
        Y <- Y[pvr_ok, , drop = FALSE]
        env_aligned <- env_aligned[pvr_ok, ]
      }
    }
  }

  # Numeric environmental variables for the model
  env_vars <- c("Temperature", "Salinity", "O2_dissolved", "AOU",
                "NO3_NO2", "PO4", "Pressure_atm", "Density_kg_m3")
  env_vars <- intersect(env_vars, colnames(env_aligned))

  if (nrow(Y) < 4) {
    log_warn(sprintf("  only %d MAGs after PVR filter; skipping species", nrow(Y)))
    next
  }

  # Jaccard distance, binary
  jac <- vegdist(Y, method = "jaccard", binary = TRUE)
  cond_terms <- intersect(c("PVR1", "PVR2", "PVR3"), colnames(env_aligned))

  # --- varpart, environment against phylogeny ---------------------------------
  # varpart returns NA for the adjusted R squared when the model is
  # undersaturated, which happens easily with few MAGs and many variables. To
  # get interpretable numbers back, the environmental variables and the number
  # of PVR axes are both trimmed to what the degrees of freedom can carry.
  # The constraint is residual_df = n - n_env - n_pvr - 1 >= 1.
  if (!is.null(pvr_axes) && length(cond_terms) > 0) {
    n_mags_vp     <- nrow(Y)
    priority_env  <- intersect(
      c("Pressure_atm", "O2_dissolved", "NO3_NO2", "Temperature",
        "AOU", "Salinity", "PO4", "Density_kg_m3"),
      env_vars
    )
    n_pvr_used    <- if (length(cond_terms) >= 1) 1 else 0
    max_predictors <- n_mags_vp - 2  # keeps at least 1 residual df
    n_env_used    <- max(0, max_predictors - n_pvr_used)
    n_env_used    <- min(length(priority_env), n_env_used)
    env_vars_vp   <- head(priority_env, n_env_used)
    pvr_terms_vp  <- head(cond_terms, n_pvr_used)

    log_info(sprintf("  varpart: n=%d MAGs · %d env vars (%s) · %d PVR (%s)",
                     n_mags_vp, length(env_vars_vp),
                     paste(env_vars_vp, collapse = ","),
                     length(pvr_terms_vp),
                     paste(pvr_terms_vp, collapse = ",")))

    vp <- if (length(env_vars_vp) >= 1 && length(pvr_terms_vp) >= 1) {
      tryCatch(
        varpart(jac,
                env_aligned[, env_vars_vp, drop = FALSE],
                env_aligned[, pvr_terms_vp, drop = FALSE]),
        error = function(e) { log_warn("  varpart failed: ", e$message); NULL }
      )
    } else {
      log_warn(sprintf("  varpart skipped: insufficient df (n=%d)", n_mags_vp))
      NULL
    }

    if (!is.null(vp)) {
      vp_df <- as_tibble(vp$part$indfract, rownames = "fraction") %>%
        mutate(species   = sp_name,
               env_vars  = paste(env_vars_vp,  collapse = ","),
               pvr_terms = paste(pvr_terms_vp, collapse = ","))
      write_tsv(vp_df, file.path(OUTPUT_DIR, paste0("varpart_", sp_name, ".tsv")))
      varpart_summary[[sp_name]] <- vp_df
      panel_varpart[[sp_name]] <- list(
        vp        = vp,
        n_mags    = n_mags_vp,
        env_vars  = env_vars_vp,
        pvr_terms = pvr_terms_vp,
        pretty    = pretty,
        undersat  = all(is.na(vp$part$indfract$Adj.R.squared))
      )

      # Generous canvas and margins, so the circle labels and the title are not
      # clipped. par(mar) is in line heights here, and family is what activates
      # the custom font in base graphics.
      png(file.path(OUTPUT_DIR, paste0("Fig_varpart_venn_", sp_name, ".png")),
          width = 7.5, height = 6.5, units = "in", res = 300)
      par(mar = c(3.0, 2.5, 5.0, 2.5),
          family = font_family,
          cex.main = 1.05, cex.lab = 0.95, xpd = NA)
      # Text sizes:
      #   cex      1.45, the numbers inside the Venn and "Residuales = ..."
      #   cex.main 1.5, the species name
      plot(vp,
           bg     = c("#1B9E77", "#7570B3"),
           Xnames = c("Ambiente", "Filogenia"),
           digits = 2,
           cex    = 1.45,
           id.size = 1.4)
      title(main = italic_title(pretty), line = 2.5, cex.main = 1.5)
      mtext(sprintf("n = %d MAGs | env: %s | cond: %s",
                    n_mags_vp,
                    paste(env_vars_vp, collapse = ", "),
                    paste(pvr_terms_vp, collapse = ", ")),
            side = 3, line = 0.6, cex = 0.78, col = "grey25")
      if (all(is.na(vp$part$indfract$Adj.R.squared))) {
        mtext("Modelo sub-saturado: no hay df suficientes para R² ajustado",
              side = 1, line = 1.0, cex = 0.85, col = "grey25")
      }
      dev.off()
    }
  }

}

# --- Overall summary ----------------------------------------------------------
if (length(varpart_summary) > 0) {
  bind_rows(varpart_summary) %>%
    write_tsv(file.path(OUTPUT_DIR, "varpart_summary_all_species.tsv"))
}

# ---- varpart panel across species ----
# plot.varpart draws with base graphics and does not compose with patchwork, so
# this opens one PNG device and lays the Venns out with par(mfrow), reproducing
# the title style the other panels use (bold, around 18 pt, centred).
if (length(panel_varpart) > 0) {
  ord  <- intersect(canonical_species, names(panel_varpart))
  vps  <- panel_varpart[ord]
  n_pl <- length(vps)
  ncol_grid <- if (n_pl <= 3) n_pl else 3
  n_rows    <- ceiling(n_pl / ncol_grid)

  panel_w <- 7.5 * ncol_grid
  panel_h <- 6.5 * n_rows + 0.9   # extra top margin for the overall title

  png(file.path(OUTPUT_DIR, "Fig_panel_varpart_venn.png"),
      width = panel_w, height = panel_h, units = "in", res = 300)

  # oma reserves a strip at the top for the centred panel title
  par(mfrow = c(n_rows, ncol_grid),
      oma   = c(0, 0, 3.5, 0),
      mar   = c(3.0, 2.5, 5.0, 2.5),
      family = font_family,
      cex.main = 1.05, cex.lab = 0.95, xpd = NA)

  for (item in vps) {
    # Same sizes and rules as the per-species render. italic_title() turns the
    # markdown "*X*" in the display name into expression(italic(X)), which
    # base R sets in italics.
    plot(item$vp,
         bg     = c("#1B9E77", "#7570B3"),
         Xnames = c("Ambiente", "Filogenia"),
         digits = 2,
         cex    = 1.45,
         id.size = 1.4)
    title(main = italic_title(item$pretty), line = 2.5, cex.main = 1.5)
    mtext(sprintf("n = %d MAGs | env: %s | cond: %s",
                  item$n_mags,
                  paste(item$env_vars,  collapse = ", "),
                  paste(item$pvr_terms, collapse = ", ")),
          side = 3, line = 0.6, cex = 0.78, col = "grey25")
    if (isTRUE(item$undersat)) {
      mtext("Modelo sub-saturado: no hay df suficientes para R² ajustado",
            side = 1, line = 1.0, cex = 0.85, col = "grey25")
    }
  }

  # Blank out the unused grid cells when the panels do not fill the last row
  n_cells <- n_rows * ncol_grid
  if (n_cells > n_pl) {
    for (k in seq_len(n_cells - n_pl)) plot.new()
  }

  mtext("Partición de varianza del genoma accesorio: ambiente vs filogenia",
        outer = TRUE, side = 3, line = 0.7,
        cex = 1.5, font = 2)
  dev.off()
  log_info("Panel varpart (Venn) written")
}

log_info("Done. Outputs in: ", OUTPUT_DIR)
