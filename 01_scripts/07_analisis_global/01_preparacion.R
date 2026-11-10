# ============================================================
# 01_preparacion.R
# Análisis KEGG para grupo agudo/crónico/postcrónico 
# Datasets 5, 5_cchp, 9, y 11
# Lectura, validación, limpieza y clasificación de grupos
# ============================================================

cat("\nCargando dds, anotación y genes curados...\n")

ruta_dds <- config$ruta_dds
ruta_anotacion <- config$ruta_anotacion
ruta_inflamatorios <- config$ruta_inflamatorios

if (!file.exists(ruta_dds)) stop("No se encontró el archivo: ", ruta_dds)
if (!file.exists(ruta_anotacion)) stop("No se encontró el archivo: ", ruta_anotacion)
if (!file.exists(ruta_inflamatorios)) stop("No se encontró el archivo: ", ruta_inflamatorios)

dds <- readRDS(ruta_dds)
genes_anotados <- read.csv(ruta_anotacion, stringsAsFactors = FALSE)
genes_literatura <- read.csv(ruta_inflamatorios, stringsAsFactors = FALSE)

if (!inherits(dds, "DESeqDataSet")) stop("El archivo RDS no corresponde a un objeto DESeqDataSet.")
if (!all(c("dataset", "grupo") %in% names(colData(dds)))) stop("La metadata debe contener las columnas 'dataset' y 'grupo'.")
if (!all(c("ensembl_id", "simbolo", "biotipo") %in% colnames(genes_anotados))) stop("La anotación debe contener: ensembl_id, simbolo y biotipo.")
if (!all(c("ensembl_id", "simbolo_humano_lit") %in% colnames(genes_literatura))) stop("La lista inflamatoria debe contener: ensembl_id y simbolo_humano_lit.")

cat("\nPreparando anotaciones...\n")

genes_anotados <- genes_anotados %>%
  mutate(ensembl_id_limpio = sub("\\..*", "", ensembl_id)) %>%
  distinct(ensembl_id_limpio, .keep_all = TRUE) %>%
  transmute(ensembl_id_limpio, simbolo_mouse = simbolo, biotipo)

genes_literatura <- genes_literatura %>%
  mutate(ensembl_id_limpio = sub("\\..*", "", ensembl_id), inflamatorio = "Sí") %>%
  distinct(ensembl_id_limpio, .keep_all = TRUE) %>%
  transmute(ensembl_id_limpio, simbolo_humano = simbolo_humano_lit, inflamatorio)

anotacion_genes <- genes_anotados %>%
  left_join(genes_literatura, by = "ensembl_id_limpio") %>%
  mutate(
    inflamatorio = ifelse(is.na(inflamatorio), "No", inflamatorio),
    simbolo_humano = ifelse(is.na(simbolo_humano), "", simbolo_humano)
  )

cat("Genes anotados preparados:", nrow(anotacion_genes), "\n")

cat("\nExcluyendo datasets indicados en config$excluir_datasets...\n")

dds <- dds[, !(colData(dds)$dataset %in% config$excluir_datasets)]
dds$dataset <- droplevels(factor(dds$dataset))

cat("\nDatasets incluidos:\n")
print(table(colData(dds)$dataset))

clasificar_grupo <- function(grupo_original) {
  grupo <- tolower(as.character(grupo_original))
  grupo <- str_replace_all(grupo, " ", "_")
  grupo <- str_replace_all(grupo, "-", "_")
  case_when(
    str_detect(grupo, "^des5w") ~ "Agudo_5w",
    str_detect(grupo, "^des12w") ~ "Cronico_12w",
    str_detect(grupo, "^rem16w") ~ "Postcronico_16w",
    str_detect(grupo, "^rem20w") ~ "Postcronico_20w",
    str_detect(grupo, "^rem24w") ~ "Postcronico_24w",
    str_detect(grupo, "^control24w|control_24w|control_24") ~ "Control_24w",
    str_detect(grupo, "control|control0w|control_0w") ~ "Control_0w",
    TRUE ~ NA_character_
  )
}

colData(dds)$grupo_original <- colData(dds)$grupo
colData(dds)$grupo_biologico <- clasificar_grupo(colData(dds)$grupo_original)
colData(dds)$grupo_biologico <- factor(colData(dds)$grupo_biologico)

cat("\nDistribución de grupos biológicos:\n")
print(table(colData(dds)$dataset, colData(dds)$grupo_biologico, useNA = "ifany"))

if (any(is.na(colData(dds)$grupo_biologico))) {
  warning("Existen muestras sin grupo biológico clasificado. Revisar la tabla anterior.")
}

cat("\nPreparación finalizada correctamente.\n")

