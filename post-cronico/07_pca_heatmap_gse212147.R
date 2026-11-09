# ============================================================
# 07_pca_heatmap_gse212147.R  (CORREGIDO)
# PCA y heatmap exploratorios restringidos a GSE212147 (db5 +
# db5_cchp), para evaluar si el "dataset" tecnico (posible proxy
# de secuenciador/lote) se confunde con el tejido en las muestras
# usadas para Postcronico_16w/20w/24w vs Control_24w.
# Requiere: 06_tejido_y_sexo.R ya corrido (dds_extendido.rds con
# grupo_biologico y tejido), y vsd.rds / vsd_sin_batch.rds
# generados por normalizacion.R / corregir_batch.R
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(dplyr)
  library(stringr)
})

dds_ext <- readRDS(file.path(config$carpeta_datos, "dds_extendido.rds"))
vsd <- readRDS(file.path(config$carpeta_datos, "vsd.rds"))
vsd_sin_batch <- readRDS(file.path(config$carpeta_datos, "vsd_sin_batch.rds"))

# vsd.rds NO tiene grupo_biologico/tejido (esas columnas se calculan
# recien en 01_preparacion.R + 06_tejido_y_sexo.R, sobre el objeto
# dds, no sobre vsd). Las traemos por nombre de muestra en vez de
# recalcularlas, para no depender de que 'grupo' este escrito igual
# en ambos objetos.
agregar_metadata_dds <- function(vsd_obj, dds_obj) {
  idx <- match(colnames(vsd_obj), colnames(dds_obj))
  if (any(is.na(idx))) {
    stop("Hay muestras en vsd que no aparecen en dds_extendido. Revisa que",
         " ambos objetos vengan del mismo dds_preparado.rds.")
  }
  colData(vsd_obj)$grupo_biologico <- colData(dds_obj)$grupo_biologico[idx]
  colData(vsd_obj)$tejido <- colData(dds_obj)$tejido[idx]
  vsd_obj
}

vsd <- agregar_metadata_dds(vsd, dds_ext)
vsd_sin_batch <- agregar_metadata_dds(vsd_sin_batch, dds_ext)

# ------------------------------------------------------------
# Si tienes el instrumento/plataforma real de secuenciacion (campo
# "Platform"/"Instrument" de SRA), arma datos/metadata_secuenciador.csv
# con columnas "muestra","secuenciador" y este bloque lo usa.
# Si no existe, se usa 'dataset' (db5 vs db5_cchp) como proxy del
# lote tecnico. OJO: dentro de GSE212147, db5 y db5_cchp tambien
# difieren parcialmente en tejido (db5 concentra Corteza), asi que
# 'dataset' no es un proxy 100% limpio de secuenciador -interpretar
# el resultado con esa salvedad-.
# ------------------------------------------------------------

ruta_secuenciador <- file.path(config$carpeta_datos, "metadata_secuenciador.csv")
usar_secuenciador_real <- file.exists(ruta_secuenciador)

if (usar_secuenciador_real) {
  metadata_sec <- read.csv(ruta_secuenciador, stringsAsFactors = FALSE)
  colData(vsd)$secuenciador <- metadata_sec$secuenciador[match(colnames(vsd), metadata_sec$muestra)]
  colData(vsd_sin_batch)$secuenciador <- metadata_sec$secuenciador[match(colnames(vsd_sin_batch), metadata_sec$muestra)]
  cat("\nUsando metadata_secuenciador.csv como variable de lote real.\n")
} else {
  colData(vsd)$secuenciador <- colData(vsd)$dataset
  colData(vsd_sin_batch)$secuenciador <- colData(vsd_sin_batch)$dataset
  cat("\n[AVISO] No se encontro metadata_secuenciador.csv.\n")
  cat("        Se usa 'dataset' (db5 vs db5_cchp) como proxy del lote tecnico.\n")
  cat("        Verifica en SRA (campo Platform/Instrument) si de verdad\n")
  cat("        corresponden a dos secuenciadores distintos.\n")
}

seleccionar_gse212147 <- function(vsd_obj) {
  seleccion <- colData(vsd_obj)$dataset %in% c("db5", "db5_cchp")
  vsd_sub <- vsd_obj[, seleccion]
  vsd_sub$dataset <- droplevels(factor(vsd_sub$dataset))
  vsd_sub
}

vsd_gse <- seleccionar_gse212147(vsd)
vsd_gse_sb <- seleccionar_gse212147(vsd_sin_batch)

cat("\nMuestras GSE212147 incluidas:", ncol(vsd_gse), "\n")
print(table(colData(vsd_gse)$grupo_biologico, colData(vsd_gse)$tejido, useNA = "ifany"))
print(table(colData(vsd_gse)$secuenciador, colData(vsd_gse)$tejido, useNA = "ifany"))

graficar_pca_gse <- function(vsd_obj, titulo, archivo) {
  pca_data <- plotPCA(vsd_obj, intgroup = c("tejido", "secuenciador"), returnData = TRUE)
  porcentaje_varianza <- round(100 * attr(pca_data, "percentVar"))

  p <- ggplot(pca_data, aes(x = PC1, y = PC2, color = tejido, shape = secuenciador)) +
    geom_point(size = 4, alpha = 0.85) +
    xlab(paste0("PC1: ", porcentaje_varianza[1], "% varianza")) +
    ylab(paste0("PC2: ", porcentaje_varianza[2], "% varianza")) +
    labs(
      title = titulo,
      color = "Tejido",
      shape = if (usar_secuenciador_real) "Secuenciador" else "Dataset (proxy secuenciador)"
    ) +
    scale_color_brewer(palette = "Set2") +
    theme_bw()

  ggsave(file.path(config$carpeta_figuras, archivo), p, width = 8, height = 6, dpi = 300)
  p
}

graficar_pca_gse(vsd_gse, "PCA GSE212147 - antes de remover batch", "pca_gse212147_tejido_lote.png")
graficar_pca_gse(vsd_gse_sb, "PCA GSE212147 - batch (dataset) removido", "pca_gse212147_tejido_lote_sin_batch.png")

graficar_heatmap_gse <- function(vsd_obj, titulo, archivo) {
  distancias <- dist(t(assay(vsd_obj)))
  matriz_distancias <- as.matrix(distancias)

  anotacion <- data.frame(
    tejido = colData(vsd_obj)$tejido,
    lote = colData(vsd_obj)$secuenciador,
    grupo = colData(vsd_obj)$grupo_biologico
  )
  rownames(anotacion) <- colnames(vsd_obj)

  png(file.path(config$carpeta_figuras, archivo), width = 1400, height = 1400, res = 150)
  pheatmap(
    matriz_distancias,
    clustering_distance_rows = distancias,
    clustering_distance_cols = distancias,
    annotation_row = anotacion,
    annotation_col = anotacion,
    col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
    fontsize_row = 6, fontsize_col = 6,
    main = titulo
  )
  dev.off()
}

graficar_heatmap_gse(vsd_gse, "GSE212147 - distancias (antes de batch)", "heatmap_gse212147_tejido_lote.png")
graficar_heatmap_gse(vsd_gse_sb, "GSE212147 - distancias (batch removido)", "heatmap_gse212147_tejido_lote_sin_batch.png")

cat("\nPCA y heatmaps de GSE212147 (tejido x lote) guardados en", config$carpeta_figuras, "\n")
cat("Revisa si tejido y lote se separan de forma limpia o si se confunden\n")
cat("entre si antes de decidir el diseno final de las comparaciones por tejido.\n")
