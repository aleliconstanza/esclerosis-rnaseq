# ============================================================
# 08_config_extendido.R
# Nuevas comparaciones: por tejido (16w/20w/24w x Corteza/
# Hipocampo/CuerpoCalloso) y pareada por sexo (24w, solo hembras).
# Requiere dds_extendido.rds generado en 06_tejido_y_sexo.R
# ============================================================

dds <- readRDS(file.path(config$carpeta_datos, "dds_extendido.rds"))

# --- Comparaciones por tejido ---
# OJO: Postcronico_24w_Corteza tiene un solo replicado en GSE212147
# (SRR21224507) y se excluye a proposito: DESeq2 no puede estimar
# dispersion de forma confiable con un solo replicado por grupo.

comparaciones_tejido <- list(
  Postcronico16w_Corteza_vs_Control24w_Corteza = list(
    grupo_ref = "Control_24w_Corteza", grupo_comp = "Postcronico_16w_Corteza"),
  Postcronico16w_Hipocampo_vs_Control24w_Hipocampo = list(
    grupo_ref = "Control_24w_Hipocampo", grupo_comp = "Postcronico_16w_Hipocampo"),
  Postcronico16w_CuerpoCalloso_vs_Control24w_CuerpoCalloso = list(
    grupo_ref = "Control_24w_CuerpoCalloso", grupo_comp = "Postcronico_16w_CuerpoCalloso"),

  Postcronico20w_Corteza_vs_Control24w_Corteza = list(
    grupo_ref = "Control_24w_Corteza", grupo_comp = "Postcronico_20w_Corteza"),
  Postcronico20w_Hipocampo_vs_Control24w_Hipocampo = list(
    grupo_ref = "Control_24w_Hipocampo", grupo_comp = "Postcronico_20w_Hipocampo"),
  Postcronico20w_CuerpoCalloso_vs_Control24w_CuerpoCalloso = list(
    grupo_ref = "Control_24w_CuerpoCalloso", grupo_comp = "Postcronico_20w_CuerpoCalloso"),

  # Postcronico24w_Corteza excluido (n=1, ver nota arriba)
  Postcronico24w_Hipocampo_vs_Control24w_Hipocampo = list(
    grupo_ref = "Control_24w_Hipocampo", grupo_comp = "Postcronico_24w_Hipocampo"),
  Postcronico24w_CuerpoCalloso_vs_Control24w_CuerpoCalloso = list(
    grupo_ref = "Control_24w_CuerpoCalloso", grupo_comp = "Postcronico_24w_CuerpoCalloso")
)

# --- Comparacion pareada por sexo (solo si 06 generó grupo_sexo real) ---
comparaciones_sexo <- list()
if ("grupo_sexo" %in% names(colData(dds)) &&
    all(c("Postcronico_24w_F", "Control_24w_F") %in% colData(dds)$grupo_sexo)) {
  comparaciones_sexo <- list(
    Postcronico24w_hembras_vs_Control24w_hembras = list(
      grupo_ref = "Control_24w_F", grupo_comp = "Postcronico_24w_F")
  )
} else {
  cat("\n[AVISO] No hay grupo_sexo con Postcronico_24w_F/Control_24w_F.\n")
  cat("        Revisa que 01_scripts/00_config/metadata/metadata_sexo.csv exista y tenga hembras\n")
  cat("        marcadas en ambos grupos (ver 06_tejido_y_sexo.R).\n")
}

comparaciones_extendidas <- c(comparaciones_tejido, comparaciones_sexo)

cat("\nComparaciones extendidas definidas:", length(comparaciones_extendidas), "\n")
print(names(comparaciones_extendidas))
