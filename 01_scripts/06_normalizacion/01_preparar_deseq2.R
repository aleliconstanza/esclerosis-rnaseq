#!/usr/bin/env Rscript

# ============================================================
# Script 1 — Preparación de conteos y metadata para DESeq2
# ============================================================
# Este script integra conteos Single-End y Paired-End generados
# por featureCounts, construye la metadata de las muestras y crea
# un objeto DESeqDataSet listo para normalización y análisis DEG.
# ============================================================

# ------------------------------------------------------------
# 1. Librerías
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)   # Objeto DESeqDataSet
  library(dplyr)    # Manipulación de tablas
  library(stringr)  # Extracción de texto desde rutas y nombres
  library(tidyr)    # Separación de rutas en columnas
})

# ------------------------------------------------------------
# 2. Rutas de entrada y salida
# ------------------------------------------------------------

ruta_single_end <- "/data/Database/Conteos/matriz_counts_single_end.txt"
ruta_paired_end <- "/data/Database/Conteos/matriz_counts_paired_end_db9.txt"

directorio_salida <- "/data/Database/Normalizacion"
dir.create(directorio_salida, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------
# 3. Leer matriz de featureCounts
# ------------------------------------------------------------
# Se eliminan las 6 primeras columnas de anotación del gen
# Las columnas posteriores corresponden a conteos por muestra

leer_conteos_featurecounts <- function(ruta_archivo) {
  
  tabla_featurecounts <- read.delim(
    ruta_archivo,
    comment.char = "#",
    check.names = FALSE
  )
  
  rownames(tabla_featurecounts) <- tabla_featurecounts$Geneid
  
  conteos <- tabla_featurecounts[, -(1:6), drop = FALSE]
  
  longitud_genes <- setNames(
    tabla_featurecounts$Length,
    tabla_featurecounts$Geneid
  )
  
  return(list(
    conteos = conteos,
    longitud_genes = longitud_genes
  ))
}

# ------------------------------------------------------------
# 4. Crear metadata desde las rutas de los BAM
# ------------------------------------------------------------
# A partir del nombre de las columnas de featureCounts se recupera:
# muestra, dataset, condición, semana y tipo de librería.

crear_metadata_desde_rutas <- function(ruta_archivo, tipo_libreria) {
  
  encabezado <- read.delim(
    ruta_archivo,
    comment.char = "#",
    check.names = FALSE,
    nrows = 0
  )
  
  rutas_bam <- colnames(encabezado)[-(1:6)]
  
  metadata <- data.frame(
    muestra = str_extract(rutas_bam, "SRR\\d+"),
    ruta_bam = rutas_bam,
    libreria = tipo_libreria,
    stringsAsFactors = FALSE
  )
  
  metadata <- metadata %>%
    mutate(ruta_relativa = str_extract(ruta_bam, "Alineamiento/.*")) %>%
    separate(
      col = ruta_relativa,
      into = c("directorio", "dataset", "grupo", "archivo_bam"),
      sep = "/",
      fill = "right"
    ) %>%
    mutate(
      condicion = case_when(
        str_detect(grupo, "^control") ~ "control",
        str_detect(grupo, "^des")     ~ "desmielinizacion",
        str_detect(grupo, "^rem")     ~ "remielinizacion",
        TRUE                          ~ "revisar"
      ),
      semana = as.numeric(str_extract(grupo, "\\d+"))
    ) %>%
    select(muestra, dataset, condicion, semana, libreria, grupo)
  
  return(metadata)
}

# ------------------------------------------------------------
# 5. Lectura de conteos
# ------------------------------------------------------------

message("\n[1/8] Leyendo matrices de conteos...")

datos_single_end <- leer_conteos_featurecounts(ruta_single_end)
datos_paired_end <- leer_conteos_featurecounts(ruta_paired_end)

message("      Single-End: ",
        nrow(datos_single_end$conteos), " genes x ",
        ncol(datos_single_end$conteos), " muestras")

message("      Paired-End: ",
        nrow(datos_paired_end$conteos), " genes x ",
        ncol(datos_paired_end$conteos), " muestras")


# ------------------------------------------------------------
# 6. Verificación de genes
# ------------------------------------------------------------
# Antes de unir las matrices, los genes deben estar en el mismo orden.

message("\n[2/8] Verificando compatibilidad entre matrices...")

if (!identical(rownames(datos_single_end$conteos),
               rownames(datos_paired_end$conteos))) {
  stop("Las matrices Single-End y Paired-End no tienen los mismos genes en el mismo orden.")
}

message("      Genes compatibles entre ambas matrices.")


# ------------------------------------------------------------
# 7. Integración Single-End + Paired-End
# ------------------------------------------------------------

message("\n[3/8] Integrando matrices Single-End y Paired-End...")

conteos_completos <- cbind(
  datos_single_end$conteos,
  datos_paired_end$conteos
)

# Se dejan los códigos SRR como nombres de columna
# La información biológica queda guardada en la metadata
colnames(conteos_completos) <- str_extract(
  colnames(conteos_completos),
  "SRR\\d+"
)

message("      Matriz integrada: ",
        nrow(conteos_completos), " genes x ",
        ncol(conteos_completos), " muestras")


# ------------------------------------------------------------
# 8. Construcción de metadata
# ------------------------------------------------------------

message("\n[4/8] Construyendo metadata de muestras...")

metadata_single_end <- crear_metadata_desde_rutas(
  ruta_archivo = ruta_single_end,
  tipo_libreria = "single_end"
)

metadata_paired_end <- crear_metadata_desde_rutas(
  ruta_archivo = ruta_paired_end,
  tipo_libreria = "paired_end"
)

metadata_muestras <- bind_rows(
  metadata_single_end,
  metadata_paired_end
)

message("      Metadata construida para ",
        nrow(metadata_muestras), " muestras")

# ------------------------------------------------------------
# 9. Sincronización matriz-metadata
# ------------------------------------------------------------
# DESeq2 requiere que las filas de metadata estén en el mismo orden
# que las columnas de la matriz de conteos.

message("\n[5/8] Sincronizando metadata con la matriz de conteos...")

rownames(metadata_muestras) <- metadata_muestras$muestra
metadata_muestras <- metadata_muestras[colnames(conteos_completos), ]

if (!identical(rownames(metadata_muestras),
               colnames(conteos_completos))) {
  stop("La metadata no coincide con las columnas de la matriz de conteos.")
}

message("      Metadata y conteos están correctamente alineados.")

# ------------------------------------------------------------
# 10. Variables del diseño experimental
# ------------------------------------------------------------

message("\n[6/8] Preparando variables del diseño experimental...")

metadata_muestras <- metadata_muestras %>%
  mutate(
    condicion = factor(
      condicion,
      levels = c("control", "desmielinizacion", "remielinizacion")
    ),
    dataset = factor(dataset),
    libreria = factor(libreria)
  )

if (any(metadata_muestras$condicion == "revisar")) {
  stop("Existen muestras con condición no reconocida. Revisar columna 'grupo'.")
}

message("      Muestras por condición:")
print(table(metadata_muestras$condicion))

message("      Muestras por dataset:")
print(table(metadata_muestras$dataset))

message("      Distribución dataset x condición:")
print(table(metadata_muestras$dataset, metadata_muestras$condicion))

message("      Distribución dataset x semana:")
print(table(metadata_muestras$dataset, metadata_muestras$semana))

# ------------------------------------------------------------
# 11. Filtro inicial de genes de baja expresión
# ------------------------------------------------------------
# Se eliminan genes con menos de 10 conteos totales en todas las muestras 
# Esto reduce ruido antes del análisis con DESeq2

message("\n[7/8] Filtrando genes con baja expresión...")

conteos_completos <- as.matrix(conteos_completos)

genes_utiles <- rowSums(conteos_completos) >= 10
conteos_filtrados <- conteos_completos[genes_utiles, ]

message("      Genes iniciales: ", nrow(conteos_completos))
message("      Genes conservados: ", nrow(conteos_filtrados))
message("      Genes eliminados: ",
        nrow(conteos_completos) - nrow(conteos_filtrados))


# ------------------------------------------------------------
# 12. Creación del objeto DESeq2
# ------------------------------------------------------------
# dataset se usa como batch

message("\n[8/8] Creando objeto DESeqDataSet...")

dds_preparado <- DESeqDataSetFromMatrix(
  countData = conteos_filtrados,
  colData = metadata_muestras,
  design = ~ dataset + condicion
)

message("      Objeto DESeqDataSet creado correctamente.")

# ------------------------------------------------------------
# 13. Guardado de resultados
# ------------------------------------------------------------

saveRDS(
  dds_preparado,
  file.path(directorio_salida, "dds_preparado.rds")
)

write.csv(
  metadata_muestras,
  file.path(directorio_salida, "metadata_muestras.csv"),
  row.names = FALSE
)

write.csv(
  conteos_filtrados,
  file.path(directorio_salida, "matriz_conteos_filtrada.csv")
)

message("\n============================================================")
message("Script 1 finalizado correctamente.")
message("Archivos generados en: ", directorio_salida)
message("  - dds_preparado.rds")
message("  - metadata_muestras.csv")
message("  - matriz_conteos_filtrada.csv")
message("============================================================\n")
