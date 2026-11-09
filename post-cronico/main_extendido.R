# ============================================================
# main_extendido.R
# Orquesta el pipeline extendido (tejido + sexo) sobre el analisis
# ya existente. Debe correrse DESPUES del pipeline oficial completo
# (00_config.R ... 05_exportar.R) porque reutiliza dds_normalizado.rds,
# anotacion_genes y las funciones de 02_funciones.R.
# ============================================================

# 1) Pipeline oficial (ya existente, sin tocar) --------------
source("00_config.R")        # config, comparaciones (oficial)
source("01_preparacion.R")   # dds con grupo_biologico + anotacion_genes
source("02_funciones.R")     # funciones DESeq2/GSEA que vamos a reutilizar

# 2) Extension: tejido + sexo ---------------------------------
source("06_tejido_y_sexo.R")             # agrega tejido/grupo_ext y sexo/grupo_sexo
source("07_pca_heatmap_gse212147.R")     # PCA/heatmap tejido x lote en GSE212147
source("08_config_extendido.R")          # define comparaciones_extendidas
source("09_analisis_extendido.R")        # DESeq2 + GSEA GO/KEGG por comparacion
source("10_integracion_exportar_extendido.R")  # Excel/CSV/figuras

cat("\n============================================================\n")
cat("Pipeline extendido (tejido + sexo) finalizado.\n")
cat("Resultados en resultados/excel, resultados/tablas, resultados/figuras\n")
cat("============================================================\n")
