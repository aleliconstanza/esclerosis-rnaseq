# ============================================================
# 06_tejido_y_sexo.R
# Extension del pipeline: extraccion de tejido (Corteza/
# Hipocampo/CuerpoCalloso) y cruce opcional con sexo para las
# comparaciones Postcronico_16w/20w/24w vs Control_24w.
# Requiere haber corrido 00_config.R y 01_preparacion.R antes
# (dds con columnas 'grupo_original' y 'grupo_biologico' listas)
# ============================================================

cat("\nExtrayendo tejido a partir de la etiqueta original de grupo...\n")

extraer_tejido <- function(grupo_original) {
  grupo <- tolower(as.character(grupo_original))
  case_when(
    str_detect(grupo, "cx$") ~ "Corteza",
    str_detect(grupo, "hp$") ~ "Hipocampo",
    str_detect(grupo, "cc$") ~ "CuerpoCalloso",
    TRUE ~ NA_character_
  )
}

colData(dds)$tejido <- extraer_tejido(colData(dds)$grupo_original)

cat("\nDistribucion de tejido por grupo biologico:\n")
print(table(colData(dds)$grupo_biologico, colData(dds)$tejido, useNA = "ifany"))

# grupo_ext: combina el grupo biologico oficial con el tejido, SOLO
# para Postcronico_16w/20w/24w y Control_24w (Agudo, Cronico_12w y
# Control_0w no se tocan; si mas adelante quieres tejido tambien
# para Cronico_12w, agregalo a este vector: tiene des12wcx/cc/hp)
grupos_con_tejido <- c(
  "Postcronico_16w", "Postcronico_20w", "Postcronico_24w", "Control_24w"
)

colData(dds)$grupo_ext <- ifelse(
  colData(dds)$grupo_biologico %in% grupos_con_tejido &
    !is.na(colData(dds)$tejido),
  paste0(colData(dds)$grupo_biologico, "_", colData(dds)$tejido),
  as.character(colData(dds)$grupo_biologico)
)
colData(dds)$grupo_ext <- factor(colData(dds)$grupo_ext)

cat("\nn por grupo_ext (grupo x tejido):\n")
print(table(colData(dds)$grupo_ext))

# ------------------------------------------------------------
# Sexo (opcional): NO viene en la metadata que ya tenemos en
# datos/ (se cruzo a mano contra SRA solo para el hallazgo de
# 24w documentado en el reporte). Si quieres la comparacion
# pareada por sexo, arma datos/metadata_sexo.csv con columnas
# "muestra","sexo" (M/F) buscando el campo "sex" de BioSample
# para cada SRR de Control_24w y Postcronico_24w (25 muestras
# en total). Sin ese archivo, este bloque se salta solo.
# ------------------------------------------------------------

ruta_sexo <- config$ruta_sexo

if (file.exists(ruta_sexo)) {
  metadata_sexo <- read.csv(ruta_sexo, stringsAsFactors = FALSE)
  if (!all(c("muestra", "sexo") %in% names(metadata_sexo))) {
    stop("metadata_sexo.csv debe tener columnas 'muestra' y 'sexo' (M/F).")
  }

  colData(dds)$sexo <- metadata_sexo$sexo[match(colnames(dds), metadata_sexo$muestra)]

  cat("\nn por grupo biologico x sexo (donde este disponible):\n")
  print(table(colData(dds)$grupo_biologico, colData(dds)$sexo, useNA = "ifany"))

  # grupo_sexo: Postcronico_24w vs Control_24w restringido a
  # hembras, para evitar el confundido de sexo ya documentado
  # (Postcronico_24w = 9H/0M vs Control_24w = 10H/6M)
  colData(dds)$grupo_sexo <- ifelse(
    colData(dds)$grupo_biologico %in% c("Postcronico_24w", "Control_24w") &
      !is.na(colData(dds)$sexo) & colData(dds)$sexo == "F",
    paste0(colData(dds)$grupo_biologico, "_F"),
    as.character(colData(dds)$grupo_biologico)
  )
  colData(dds)$grupo_sexo <- factor(colData(dds)$grupo_sexo)

  cat("\nn por grupo_sexo:\n")
  print(table(colData(dds)$grupo_sexo))
} else {
  cat("\n[AVISO] No se encontro", ruta_sexo, "\n")
  cat("        La comparacion pareada por sexo (24w solo hembras) no se\n")
  cat("        ejecutara. Genera ese CSV (ver comentario arriba) antes de\n")
  cat("        correr 08_config_extendido.R si la quieres incluir.\n")
  colData(dds)$sexo <- NA_character_
  colData(dds)$grupo_sexo <- as.character(colData(dds)$grupo_biologico)
}

saveRDS(dds, file.path(config$carpeta_datos, "dds_extendido.rds"))
cat("\ndds_extendido.rds guardado con columnas nuevas: tejido, grupo_ext, sexo, grupo_sexo\n")
