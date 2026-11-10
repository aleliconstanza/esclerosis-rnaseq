#!/usr/bin/env Rscript
# ============================================================
# Mapeo de genes de literatura (símbolo humano) -> ENSEMBL ratón
# Entrada: genes_literatura_curado.txt (un símbolo por línea)
# Salida: genes_literatura_mapeado.csv
# ============================================================

# --- 1. Instalar dependencias (solo la primera vez) ---
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
paquetes <- c("org.Mm.eg.db", "org.Hs.eg.db", "limma")
faltantes <- paquetes[!sapply(paquetes, requireNamespace, quietly = TRUE)]
if (length(faltantes) > 0) BiocManager::install(faltantes, update = FALSE, ask = FALSE)

suppressPackageStartupMessages({
  library(org.Mm.eg.db)
  library(org.Hs.eg.db)
  library(limma)
})

# --- 2. Cargar lista de símbolos humanos (literatura) ---
genes_humano <- readLines("/data/Database/Mapeo/genes_literatura_curado.txt")
genes_humano <- toupper(trimws(genes_humano))
genes_humano <- unique(genes_humano[genes_humano != ""])

# --- 3. Símbolo humano -> símbolo de ratón ---
# La convención de ratón es Titlecase (ej. GFAP -> Gfap).
# alias2SymbolTable resuelve también sinónimos/nombres antiguos.
simbolo_raton <- limma::alias2SymbolTable(genes_humano, species = "Mm")

tabla <- data.frame(
  simbolo_humano = genes_humano,
  simbolo_raton  = simbolo_raton,
  stringsAsFactors = FALSE
)

# Para los que no resolvió alias2SymbolTable, intentar conversión simple
# de mayúsculas a Titlecase como respaldo
sin_match <- is.na(tabla$simbolo_raton)
tabla$simbolo_raton[sin_match] <- tools::toTitleCase(tolower(tabla$simbolo_humano[sin_match]))

# --- 4. Símbolo de ratón -> ENSEMBL de ratón ---
ensembl_raton <- mapIds(
  org.Mm.eg.db,
  keys = tabla$simbolo_raton,
  column = "ENSEMBL",
  keytype = "SYMBOL",
  multiVals = "first"
)

tabla$ensembl_mouse <- ensembl_raton[tabla$simbolo_raton]

# --- 5. Reporte ---
n_total <- nrow(tabla)
n_mapeados <- sum(!is.na(tabla$ensembl_mouse))

cat("============================================================\n")
cat("Genes de literatura:", n_total, "\n")
cat("Mapeados a ENSEMBL de ratón:", n_mapeados, "\n")
cat("Sin mapear:", n_total - n_mapeados, "\n")
cat("============================================================\n")

if (n_total - n_mapeados > 0) {
  cat("\nGenes SIN mapear (revisar manualmente, pueden ser sinónimos no estándar):\n")
  print(tabla[is.na(tabla$ensembl_mouse), c("simbolo_humano", "simbolo_raton")])
}

# --- 6. Guardar resultado ---
write.csv(tabla, "/data/Database/Mapeo/genes_literatura_mapeado.csv", row.names = FALSE)
cat("\nGuardado en genes_literatura_mapeado.csv\n")

