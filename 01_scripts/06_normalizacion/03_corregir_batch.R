#!/usr/bin/env Rscript

# ============================================================
# Script 3 — Remoción visual del efecto Batch para gráficos
# ============================================================
# Remueve el efecto asociado al dataset únicamente para
# visualización exploratoria.
#
# IMPORTANTE:
# Este script NO modifica el objeto utilizado para DESeq2.
# La corrección se utiliza únicamente para PCA y Heatmap.
# ============================================================

# ------------------------------------------------------------
# 1. Librerías
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
  library(limma)
})

# ------------------------------------------------------------
# 2. Directorio de trabajo
# ------------------------------------------------------------

directorio_trabajo <- "/data/Database/Normalizacion"
setwd(directorio_trabajo)

# ------------------------------------------------------------
# 3. Cargar objeto transformado por VST
# ------------------------------------------------------------

message("\n[1/4] Cargando objeto VST...")

if (!file.exists("vsd.rds")) {
  stop("No se encontró el archivo vsd.rds.")
}

vsd <- readRDS("vsd.rds")

message("      Objeto cargado correctamente.")
message("      Genes: ", nrow(vsd))
message("      Muestras: ", ncol(vsd))


# ------------------------------------------------------------
# 4. Remoción del efecto Dataset (Batch)
# ------------------------------------------------------------
# removeBatchEffect elimina matemáticamente la variación
# asociada al dataset de origen.
#
# Se conserva la información biológica de la condición
# experimental mediante la matriz de diseño.
#
# El objetivo es visualizar con mayor claridad la señal
# biológica presente en los datos.
# ------------------------------------------------------------

message("\n[2/4] Removiendo efecto Batch/Dataset...")

matriz_vst <- assay(vsd)

dataset <- colData(vsd)$dataset
condicion <- colData(vsd)$condicion

matriz_sin_batch <- removeBatchEffect(
  matriz_vst,
  batch = dataset,
  design = model.matrix(~ condicion)
)

message("      Efecto batch removido.")


# ------------------------------------------------------------
# 5. Crear nuevo objeto VST corregido
# ------------------------------------------------------------

vsd_sin_batch <- vsd

assay(vsd_sin_batch) <- matriz_sin_batch

message("      Objeto corregido generado.")

# ------------------------------------------------------------
# 6. PCA posterior a la corrección
# ------------------------------------------------------------

message("\n[3/4] Generando PCA corregido...")

pca_data <- plotPCA(
  vsd_sin_batch,
  intgroup = c("condicion", "dataset"),
  returnData = TRUE
)

porcentaje_varianza <- round(
  100 * attr(pca_data, "percentVar")
)

pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = condicion,
    shape = dataset
  )
) +
  geom_point(
    size = 4,
    alpha = 0.85
  ) +
  xlab(
    paste0(
      "PC1: ",
      porcentaje_varianza[1],
      "% varianza"
    )
  ) +
  ylab(
    paste0(
      "PC2: ",
      porcentaje_varianza[2],
      "% varianza"
    )
  ) +
  labs(
    title = "PCA Corregido",
    color = "Condición",
    shape = "Batch"
  ) +
  scale_color_manual(
    values = c(
      "control" = "#2ecc71",
      "desmielinizacion" = "#e74c3c",
      "remielinizacion" = "#3498db"
    )
  ) +
  theme_bw()

ggsave(
  "pca_sin_batch.png",
  pca_plot,
  width = 8,
  height = 6,
  dpi = 300
)

message("      PCA corregido guardado.")

# ------------------------------------------------------------
# 7. Heatmap posterior a la corrección
# ------------------------------------------------------------

message("\n[4/4] Generando Heatmap corregido...")
distancias_muestras <- dist(
  t(
    assay(vsd_sin_batch)
  )
)

matriz_distancias <- as.matrix(
  distancias_muestras
)

anotacion_muestras <- data.frame(
  condicion = colData(vsd_sin_batch)$condicion,
  batch = colData(vsd_sin_batch)$dataset,
  libreria = colData(vsd_sin_batch)$libreria
)

rownames(anotacion_muestras) <- colnames(vsd_sin_batch)

png(
  "heatmap_sin_batch.png",
  width = 1400,
  height = 1400,
  res = 150
)

pheatmap(
  matriz_distancias,
  clustering_distance_rows = distancias_muestras,
  clustering_distance_cols = distancias_muestras,
  annotation_row = anotacion_muestras,
  annotation_col = anotacion_muestras,
  col = colorRampPalette(
    rev(
      brewer.pal(9, "Blues")
    )
  )(255),
  fontsize_row = 6,
  fontsize_col = 6,
  main = "Distancia entre muestras (Batch removido)"
)

dev.off()

message("      Heatmap corregido guardado.")


# ------------------------------------------------------------
# 8. Guardar objeto corregido
# ------------------------------------------------------------

saveRDS(vsd_sin_batch, "vsd_sin_batch.rds")

message("\n============================================================")
message(" Script 3 finalizado correctamente")
message("")
message(" Archivos generados:")
message("   - vsd_sin_batch.rds")
message("   - pca_sin_batch.png")
message("   - heatmap_sin_batch.png")
message("============================================================\n")
