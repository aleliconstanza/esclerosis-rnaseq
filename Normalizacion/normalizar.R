#!/usr/bin/env Rscript

# ============================================================
# Script 2 — Normalización y análisis exploratorio
# Calcula la normalización de DESeq2 y genera gráficos
# exploratorios usando datos transformados por VST
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)        # Normalización y transformación VST
  library(ggplot2)       # Construcción del PCA
  library(pheatmap)      # Heatmap de distancias
  library(RColorBrewer)  # Paleta de colores para heatmap
})

# ------------------------------------------------------------
# 1. Rutas de entrada y salida
# ------------------------------------------------------------

directorio_trabajo <- "/data/Database/Normalizacion"
setwd(directorio_trabajo)

archivo_dds <- "dds_preparado.rds"

# ------------------------------------------------------------
# 2. Cargar objeto dds de código 1
# ------------------------------------------------------------

message("\n[1/5] Cargando objeto DESeq2 preparado...")

if (!file.exists(archivo_dds)) {
  stop("No se encontró el archivo dds_preparado.rds en la carpeta de normalización.")
}

dds <- readRDS(archivo_dds)

# ------------------------------------------------------------
# 3. Calcular factores de tamaño
# ------------------------------------------------------------
# Los size factors corrigen diferencias entre muestras 
# causadas por distinta profundidad de secuenciación.
# Esta es la normalización principal de DESeq2

message("\n[2/5] Calculando size factors...")

dds <- estimateSizeFactors(dds)

tabla_size_factors <- data.frame(
  muestra = colnames(dds),
  size_factor = sizeFactors(dds),
  condicion = colData(dds)$condicion,
  dataset = colData(dds)$dataset,
  libreria = colData(dds)$libreria
)

write.csv(tabla_size_factors,"size_factors.csv", row.names = FALSE)

message("      Size factors calculados y guardados en size_factors.csv.")

# ------------------------------------------------------------
# 4. Transformación VST
# ------------------------------------------------------------
# La transformación VST estabiliza la varianza de los conteos.
# Se usa para PCA y heatmap, no reemplaza el análisis diferencial.

message("\n[3/5] Aplicando transformación VST...")

vsd <- vst(dds, blind = TRUE)

message("      Transformación VST completada.")

# ------------------------------------------------------------
# 5. PCA con datos transformados
# ------------------------------------------------------------
# Permite observar agrupamientos generales entre muestras.
# Internamente se usa dataset, pero en el gráfico se rotula como Batch.

message("\n[4/5] Generando PCA...")
pca_data <- plotPCA(vsd,intgroup = c("condicion", "dataset"),returnData = TRUE)
porcentaje_varianza <- round(100 * attr(pca_data, "percentVar"))

pca_plot <- ggplot(
  pca_data,
  aes(x = PC1, y = PC2, color = condicion, shape = dataset)
) +
  geom_point(size = 4, alpha = 0.85) +
  xlab(paste0("PC1: ", porcentaje_varianza[1], "% varianza")) +
  ylab(paste0("PC2: ", porcentaje_varianza[2], "% varianza")) +
  labs(
    title = "PCA de muestras normalizadas",
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


ggsave("pca.png", pca_plot, width = 8, height = 6, dpi = 300)
message("      PCA guardado como pca.png.")

# ------------------------------------------------------------
# 6. Heatmap de distancias entre muestras
# ------------------------------------------------------------
# Resume la distancia global entre muestras usando VST

message("\n[5/5] Generando heatmap de distancias...")

distancias_muestras <- dist(t(assay(vsd)))
matriz_distancias <- as.matrix(distancias_muestras)

anotacion_muestras <- data.frame(
  condicion = colData(vsd)$condicion,
  batch = colData(vsd)$dataset,
  libreria = colData(vsd)$libreria
)
rownames(anotacion_muestras) <- colnames(vsd)

png("heatmap.png",width = 1400, height = 1400, res = 150)
pheatmap(
  matriz_distancias,
  clustering_distance_rows = distancias_muestras,
  clustering_distance_cols = distancias_muestras,
  annotation_row = anotacion_muestras,
  annotation_col = anotacion_muestras,
  col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
  fontsize_row = 6, fontsize_col = 6,
  main = "Distancia entre muestras normalizadas"
)

dev.off()
message("      Heatmap guardado como heatmap.png.")

# ------------------------------------------------------------
# 7. Guardar objetos
# ------------------------------------------------------------

saveRDS(dds, "dds_normalizado.rds")
saveRDS(vsd, "vsd.rds")

message("\n============================================================")
message("Script 2 finalizado correctamente.")
message("Archivos generados:")
message("  - dds_normalizado.rds")
message("  - vsd.rds")
message("  - size_factors.csv")
message("  - pca.png")
message("  - heatmap.png")
message("============================================================\n")
