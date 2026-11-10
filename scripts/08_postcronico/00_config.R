# ============================================================
# 00_config.R
# Análisis KEGG para grupo agudo/crónico/postcrónico 
# Datasets 5, 5_cchp, 9, y 11
# Configuración general del proyecto para análisis KEGG
# ============================================================

suppressPackageStartupMessages({
  library(DESeq2)
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(enrichplot)
  library(ggplot2)
  library(ggrepel)
  library(openxlsx)
})

config <- list(
  organismo = "mouse",
  padj_significativo = 0.05,
  semilla = 123,                 # semilla fija para GSEA (permutaciones reproducibles)
  usar_tejido_en_diseno = TRUE,  # agrega 'tejido' como covariable cuando el modelo lo permite
  gsea_pvalueCutoff = 1,
  minGSSize = 5,
  maxGSSize = 1000,
  excluir_datasets = c("db8"),
  # Rutas relativas a scripts/08_postcronico/ (desde donde se corre main.R)
  carpeta_datos = file.path("..", "..", "intermedios", "08_postcronico"),
  carpeta_resultados = file.path("..", "..", "resultados", "08_postcronico"),
  carpeta_excel = file.path("..", "..", "resultados", "08_postcronico", "excel"),
  carpeta_tablas = file.path("..", "..", "resultados", "08_postcronico", "tablas"),
  carpeta_figuras = file.path("..", "..", "resultados", "08_postcronico", "figuras"),
  ruta_dds = file.path("..", "..", "intermedios", "06_normalizacion", "dds_normalizado.rds"),
  ruta_vsd = file.path("..", "..", "intermedios", "06_normalizacion", "vsd.rds"),
  ruta_vsd_sin_batch = file.path("..", "..", "intermedios", "06_normalizacion", "vsd_sin_batch.rds"),
  ruta_anotacion = file.path("..", "..", "intermedios", "05_anotacion", "genes_counts_anotados.csv"),
  ruta_inflamatorios = file.path("..", "..", "intermedios", "05_anotacion", "genes_literatura_en_counts.csv"),
  ruta_sexo = file.path("..", "00_config", "metadata", "metadata_sexo.csv"),
  archivo_dds = "dds_normalizado.rds",
  archivo_anotacion = "genes_counts_anotados.csv",
  archivo_inflamatorios = "genes_literatura_en_counts.csv",
  archivo_excel = "Analisis_GSEA_KEGG.xlsx"
)

dir.create(config$carpeta_datos, showWarnings = FALSE, recursive = TRUE)
dir.create(config$carpeta_resultados, showWarnings = FALSE, recursive = TRUE)
dir.create(config$carpeta_excel, showWarnings = FALSE, recursive = TRUE)
dir.create(config$carpeta_tablas, showWarnings = FALSE, recursive = TRUE)
dir.create(config$carpeta_figuras, showWarnings = FALSE, recursive = TRUE)

comparaciones <- list(
  Agudo_vs_Control = list(grupo_ref = "Control_0w", grupo_comp = "Agudo_5w"),
  Cronico_vs_Control = list(grupo_ref = "Control_0w", grupo_comp = "Cronico_12w"),
  Cronico_vs_Agudo = list(grupo_ref = "Agudo_5w", grupo_comp = "Cronico_12w"),
  Postcronico16w_vs_Control24w = list(grupo_ref = "Control_24w", grupo_comp = "Postcronico_16w"),
  Postcronico20w_vs_Control24w = list(grupo_ref = "Control_24w", grupo_comp = "Postcronico_20w"),
  Postcronico24w_vs_Control24w = list(grupo_ref = "Control_24w", grupo_comp = "Postcronico_24w")
)

cat("\nConfiguración cargada correctamente.\n")
