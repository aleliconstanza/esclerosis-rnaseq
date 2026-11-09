# ============================================================
# Archivo main.R
# Pipeline completo de análisis funcional RNA-seq
# Análisis KEGG para grupo agudo/crónico/postcrónico 
# Datasets 5, 5_cchp, 9, y 11
# ============================================================

rm(list = ls())

cat("\n=============================================\n")
cat("  Análisis funcional RNA-seq\n")
cat("  Desmielinización y Remielinización\n")
cat("=============================================\n")

cat("Ejecutando 00_config.R...\n")
source("00_config.R")
cat("Ejecutando 01_preparacion.R...\n")
source("01_preparacion.R")
cat("Ejecutando 02_funciones.R...\n")
source("02_funciones.R")
cat("Ejecutando 03_analisis.R...\n")
source("03_analisis.R")
cat("Ejecutando 04_integracion.R...\n")
source("04_integracion.R")
cat("Ejecutando 05_exportar.R...\n")
source("05_exportar.R")

cat("\n=============================================\n")
cat("Análisis finalizado correctamente.\n")
cat("Los resultados se encuentran en:\n")
cat("   resultados/\n")
cat("=============================================\n")
