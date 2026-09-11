#!/usr/bin/env Rscript

# Compares the trees that pipeline step 13 produced across its nine
# combinations of alignment strategy and model.
#
# Every tree is rooted and made ultrametric first, because cophenetic
# correlation is only meaningful between trees on the same footing. Pairs are
# then ranked by that correlation, and the ten most similar get a tanglegram so
# the disagreements can be read directly.
# Fecha: 2026/02/03

user_lib <- Sys.getenv("R_LIBS_USER")
if (user_lib == "") {
  user_lib <- file.path(Sys.getenv("HOME"), "R", "x86_64-pc-linux-gnu-library",
                        paste(R.version$major, strsplit(R.version$minor, "\\.")[[1]][1], sep = "."))
}

if (!dir.exists(user_lib)) {
  dir.create(user_lib, recursive = TRUE, showWarnings = FALSE)
}

.libPaths(c(user_lib, .libPaths()))
fecha <- format(Sys.Date(), "%Y%m%d")

if (!require("pacman")) {
  install.packages("pacman")
}

required_packages <- c("dendextend", "ape", "phangorn", "corrplot")
pacman::p_load(char = required_packages, install = TRUE, update = FALSE)

# ---- fonts ----
font_family <- "sans"
available_fonts <- sysfonts::font_files()
font_priority <- c(
  "Aptos", "Aptos Display",
  "Calibri", "Carlito",
  "Arial", "Liberation Sans",
  "Open Sans", "Lato", "Roboto"
)

font_found <- FALSE
for (font_name in font_priority) {
  matching_fonts <- available_fonts[grepl(font_name, available_fonts$family, ignore.case = TRUE), ]
  if (nrow(matching_fonts) > 0) {
    tryCatch({
      font_path <- matching_fonts$file[1]
      sysfonts::font_add(
        family = "custom_font",
        regular = font_path
      )
      font_family <- "custom_font"
      cat("[OK] Using font:", font_name, "\n")
      font_found <- TRUE
      break
    }, error = function(e) {})
  }
}

if (!font_found) {
  cat("[WARN] No suitable system font found. Attempting to download from Google Fonts...\n")
  tryCatch({
    sysfonts::font_add_google("Open Sans", "custom_font")
    font_family <- "custom_font"
    cat("[OK] Downloaded and using: Open Sans (Google Fonts)\n")
    font_found <- TRUE
  }, error = function(e) {
    cat("[WARN] Could not download Google Fonts. Using default sans-serif font.\n")
  })
}

cat("[OK] Font configuration complete\n\n")

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
      plot.title = element_text(face = "bold", hjust = 0.5)
    )
}

# ---- 1. CARGA Y PREPARACIÓN DE DATOS ----

OUTPUT_DIR <- paste0(file.path(PROJECT_DIR, "results", "SO_02", "phylogeny"), "/",
                     fecha, "_analysis")
dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)
cat("Directorio de salida:", OUTPUT_DIR, "\n")

ruta <- "\\\\wsl.localhost\\Ubuntu-24.04\\home\\ruben\\msc_archaea\\results\\13_phylogeny"
archivos <- list.files(path = ruta,
                       pattern = "\\.treefile$",
                       full.names = TRUE,
                       recursive = TRUE)

# Read every tree into one list
arboles <- lapply(archivos, read.tree)

# Name each entry by its file, without the path
names(arboles) <- basename(archivos)

dends <- lapply(arboles, function(tree) {
  # Cophenetic correlation needs a rooted tree
  if(!is.rooted(tree)) {
    tree <- midpoint(tree)
  }
  # and an ultrametric one, so branch lengths are comparable between trees
  tree_ultra <- chronos(tree, quiet = TRUE)
  # Convertir a dendrograma
  as.dendrogram(as.hclust.phylo(tree_ultra))
})

# Every pair of trees
n <- length(dends)
combinaciones <- combn(n, 2)

# Results table
resultados <- data.frame(
  arbol1 = character(),
  arbol2 = character(),
  entanglement = numeric(),
  cophenetic_cor = numeric(),
  baker_gamma = numeric(),
  stringsAsFactors = FALSE
)

# Compare pairwise
for(i in 1:ncol(combinaciones)) {
  idx1 <- combinaciones[1, i]
  idx2 <- combinaciones[2, i]

  nombre1 <- names(dends)[idx1]
  nombre2 <- names(dends)[idx2]

  # dendextend compares through a dendlist
  dend_pair <- dendlist(dends[[idx1]], dends[[idx2]])

  # Comparison metrics
  entang <- entanglement(dend_pair)
  coph_cor <- cor_cophenetic(dend_pair)
  baker <- cor_bakers_gamma(dend_pair)

  # Guardar resultados
  resultados <- rbind(resultados, data.frame(
    arbol1 = nombre1,
    arbol2 = nombre2,
    entanglement = entang,
    cophenetic_cor = coph_cor,
    baker_gamma = baker
  ))

  cat(sprintf("Comparado: %s vs %s\n", nombre1, nombre2))
}

# Most similar pairs first
resultados_ordenados <- resultados[order(-resultados$cophenetic_cor), ]

# ---- 2. GUARDAR RESULTADOS ----

# Guardar tabla de resultados en formato TSV
tsv_file <- file.path(OUTPUT_DIR, "comparaciones_arboles.tsv")
write.table(resultados_ordenados, tsv_file,
            sep = "\t",
            row.names = FALSE,
            quote = FALSE)
cat("\n[OK] Tabla de resultados guardada en:", tsv_file, "\n")

# ---- ---- 3. correlation matrix across every tree ----

cat("\n--- Generando matriz de correlación ---\n")

# One dendlist holding every tree
nombres_limpios <- gsub("\\.treefile$", "", names(dends))

# do.call unpacks the list into dendlist arguments
dend_list <- do.call(dendlist, dends)
names(dend_list) <- nombres_limpios

# Correlation matrix across all trees
cat("Calculando correlaciones cophenéticas...\n")
cors <- cor.dendlist(dend_list)

# Save it as TSV
cors_file <- file.path(OUTPUT_DIR, "matriz_correlacion.tsv")
write.table(round(cors, 3), cors_file,
            sep = "\t",
            quote = FALSE,
            col.names = NA)
cat("[OK] Matriz de correlación guardada en:", cors_file, "\n")

# Correlation matrix, pie rendering
png(file.path(OUTPUT_DIR, "correlacion_matriz_pie.png"),
    width = 1600, height = 1600, res = 200)
corrplot(cors, method = "pie", type = "lower",
         tl.col = "black", tl.cex = 0.8,
         title = "Correlación entre árboles filogenéticos (cophenetic)",
         mar = c(0, 0, 2, 0))
dev.off()
cat("[OK] Gráfico de correlación (pie) guardado\n")

# Correlation matrix, colour rendering
png(file.path(OUTPUT_DIR, "correlacion_matriz_color.png"),
    width = 1600, height = 1600, res = 200)
corrplot(cors, method = "color", type = "lower",
         addCoef.col = "black", number.cex = 0.6,
         tl.col = "black", tl.cex = 0.8,
         title = "Correlación entre árboles filogenéticos (cophenetic)",
         mar = c(0, 0, 2, 0))
dev.off()
cat("[OK] Gráfico de correlación (color) guardado\n")

# Correlation matrix, circle rendering
png(file.path(OUTPUT_DIR, "correlacion_matriz_circle.png"),
    width = 1600, height = 1600, res = 200)
corrplot(cors, method = "circle", type = "lower",
         tl.col = "black", tl.cex = 0.8,
         title = "Correlación entre árboles filogenéticos (cophenetic)",
         mar = c(0, 0, 2, 0))
dev.off()
cat("[OK] Gráfico de correlación (circle) guardado\n")

# ---- 4. GENERAR Y GUARDAR TANGLEGRAMS ----

cat("\n--- Generando tanglegrams ---\n")

# Tanglegrams for the ten most similar pairs
n_tanglegrams <- min(10, nrow(resultados_ordenados))

for(i in 1:n_tanglegrams) {
  idx1 <- which(names(dends) == resultados_ordenados$arbol1[i])
  idx2 <- which(names(dends) == resultados_ordenados$arbol2[i])

  nombre1 <- resultados_ordenados$arbol1[i]
  nombre2 <- resultados_ordenados$arbol2[i]

  # Crear nombre de archivo limpio
  nombre_archivo <- paste0("tanglegram_", sprintf("%02d", i), "_",
                           gsub("\\.treefile$", "", nombre1), "_vs_",
                           gsub("\\.treefile$", "", nombre2), ".png")
  png_path <- file.path(OUTPUT_DIR, nombre_archivo)

  # Tanglegram, with the text scaled to fit
  png(png_path, width = 2000, height = 1200, res = 200)

  # Margins
  par(mar = c(1, 1, 3, 1))

  dend_pair <- dendlist(dends[[idx1]], dends[[idx2]])

  # Main title, reduced so it fits the canvas
  titulo_main <- sprintf("Cophenetic cor: %.3f | Entanglement: %.3f",
                         resultados_ordenados$cophenetic_cor[i],
                         resultados_ordenados$entanglement[i])

  # Tree names, likewise reduced
  nombre1_limpio <- gsub("\\.treefile$", "", nombre1)
  nombre2_limpio <- gsub("\\.treefile$", "", nombre2)

  tanglegram(dend_pair,
             main_left = "",  # Vacío, lo añadiremos manualmente
             main_right = "", # Vacío, lo añadiremos manualmente
             margin_inner = 10,
             lab.cex = 0.5,              # Mantiene tamaño de etiquetas MAGs
             columns_width = c(6, 1.5, 6),
             main = "",                   # Vacío, lo añadiremos manualmente
             cex_main = 0.8)             # Este parámetro no funciona en tanglegram

  # Title drawn by hand, since tanglegram does not size it well
  title(main = titulo_main, cex.main = 0.9)

  # Tree names, each aligned toward its own side: the left one right-aligned
  # (adj = 1), the right one left-aligned (adj = 0), so both sit against the
  # tanglegram rather than drifting into it.
  mtext(nombre2_limpio, side = 3, line = 0.5, at = 0.77, adj = 0, cex = 0.8)

  dev.off()

  cat(sprintf("  [%d/%d] Guardado: %s\n", i, n_tanglegrams, nombre_archivo))
}

cat("\n[OK] Análisis completo. Todos los archivos guardados en:\n")
cat("  ", OUTPUT_DIR, "\n")
