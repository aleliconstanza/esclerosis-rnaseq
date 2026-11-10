# ============================================================
# main.R  (post-cronico)
# Pipeline oficial: DESeq2 + GSEA GO/KEGG para las 6 comparaciones.
# Uso desde terminal (dentro de post-cronico/):
#   Rscript main.R
# ============================================================

rm(list = ls())
inicio <- Sys.time()

source("00_config.R")
source("01_preparacion.R")
source("02_funciones.R")
source("03_analisis.R")
source("04_integracion.R")
source("05_exportar.R")

# Registro de reproducibilidad
if (exists("disenos_usados")) {
  write.csv(
    data.frame(comparacion = names(disenos_usados), diseno = unlist(disenos_usados)),
    file.path(config$carpeta_tablas, "disenos_usados.csv"), row.names = FALSE
  )
}
writeLines(capture.output(sessionInfo()), file.path(config$carpeta_resultados, "sessionInfo.txt"))

cat("\nPipeline oficial finalizado en", round(difftime(Sys.time(), inicio, units = "mins"), 1), "min.\n")
